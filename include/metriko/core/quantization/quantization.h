//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_H
#define METRIKO_QUANTIZATION_H
#include "quantization_basisloop.h"
#include "quantization_weight.h"
#include "quantization_constraint.h"
#include "quantization_evaluation.h"
#include "quantization_validation.h"

namespace metriko {
    inline void validate_quantization(const Tmesh &tmesh, const VecXd &X) {
        for (const Tquad &tquad: tmesh.tquads) {
            for (int i = 0; i < 2; i++) {
                auto thidsA = tquad.thids_by_side(i);
                auto thidsB = tquad.thids_by_side(i + 2);
                int sumA = 0;
                int sumB = 0;
                for (int thid: thidsA) { sumA += (int) X[tmesh.thalfs[thid].edge().id]; }
                for (int thid: thidsB) { sumB += (int) X[tmesh.thalfs[thid].edge().id]; }
                assert(sumA == sumB);
            }
        }
    }

    // an arc of the strip graph: from thalf `fr` (entering its tquad) to thalf
    // `to` = twin of a thalf on the opposite side. traversing the arc means the
    // loop passes the tedge of `to`.
    struct StripArc { int fr; int to; int teid; };

    inline std::vector<StripArc> build_strip_arcs(
        const std::vector<Tquad>& tquads,
        const std::vector<Thalf>& thalfs,
        const VecXi& th2quad,
        const VecXi& th2side
    ) {
        std::vector<StripArc> arcs;
        for (const Thalf& th: thalfs) {
            const auto& tq = tquads[th2quad[th.id]];
            for (const auto& d: tq.data)
                if (d.side == (th2side[th.id] + 2) % 4) {
                    const Thalf& to = thalfs[d.thid].twin();
                    arcs.push_back({th.id, to.id, to.teid});
                }
        }
        return arcs;
    }

    // pricing: find one negative-cost cycle in the strip graph (bellman-ford
    // with a virtual super source), or return empty when none exists. cycle
    // cost = sum of w[teid] over its nodes, i.e. the exact energy change of
    // applying +-1 along the loop (up to tedges hit twice, re-checked by the
    // caller).
    inline std::vector<int> find_negative_cycle(
        const int n_thalfs,
        const std::vector<StripArc>& arcs,
        const VecXd& w
    ) {
        vec<double> dist(n_thalfs, 0.);
        vec<int>    pred(n_thalfs, -1);
        int last = -1;
        for (int it = 0; it < n_thalfs; ++it) {
            last = -1;
            for (const auto& [fr, to, teid]: arcs) {
                if (dist[fr] + w[teid] < dist[to] - 1e-12) {
                    dist[to] = dist[fr] + w[teid];
                    pred[to] = fr;
                    last = to;
                }
            }
            if (last == -1) return {};   // converged: no negative cycle
        }
        // still relaxing after n passes: walk n steps back to land inside the
        // cycle, then collect it
        int x = last;
        for (int i = 0; i < n_thalfs; ++i) x = pred[x];
        std::vector<int> cyc = {x};
        for (int v = pred[x]; v != x; v = pred[v]) cyc.push_back(v);
        return cyc;
    }

    inline VecXd compute_quantization(const Tmesh& tmesh, const Mgrph& mg) {

        VecXd R(tmesh.nTE);
        for (int i = 0; i < tmesh.nTE; i++) R[i] = tmesh.tedges[i].len;

       VecXd X = VecXd::Zero(tmesh.tedges.size());
       MatXd C = compute_constraint(tmesh);
       VecXd I = VecXd::Ones(tmesh.tedges.size());
       MatXd G = construct_generating_vectors(
           tmesh,
           R,
           [](const Comparator &c1, const Comparator &c2) {
               return c1.length / std::max(c1.weight, 1e-9)
                    < c2.length / std::max(c2.weight, 1e-9);
           }
           );

       // ----- construct first step vector ----- //
       while ((X.array() == 0).any()) {
           double min = 1e+9;
           int thid = 0;
           for (auto &th: tmesh.thalfs | vw::filter([&](auto &th_) { return X[th_.edge().id] == 0; })) {
               double w = compute_weight(R[th.edge().id], X[th.edge().id], tmesh.tedges.size());
               if (w < min) {
                   min = w;
                   thid = th.id;
               }
           }
           for (auto g: G.colwise()) {
               if (g[tmesh.thalfs[thid].edge().id] != 0) {
                   X += g;
                   break;
               }
           }
           if (compute_validation(mg, tmesh, X)) {
               std::cout << "validation passed" << std::endl;
               break;
           }
       }

        assert((C * X).norm() == 0);

       // ----- construct second step vector ----- //

       double e = (X.cwiseQuotient(R) - I).norm();

        // ----- second step: exact pricing (column generation) ----- //
        // apply negative-cost loops until none exists. the pricing searches ALL
        // strip loops, not a precomputed basis, so the result is locally optimal
        // w.r.t. every single-loop +-1 move.
        auto arcs = build_strip_arcs(tmesh.tquads, tmesh.thalfs, tmesh.th2quad, tmesh.th2side);
        const int nte = (int)tmesh.tedges.size();
        auto energy = [&](const VecXd& x) { return (x.cwiseQuotient(R) - I).norm(); };

        //double e = energy(X);
        bool improved = true;
        for (int guard = 0; improved && guard < 100 * nte; ++guard) {
            improved = false;
            for (int sgn: {+1, -1}) {
                // marginal cost of changing tedge j by sgn:
                // ((X_j+sgn)/R_j - 1)^2 - (X_j/R_j - 1)^2
                VecXd w(nte);
                for (int j = 0; j < nte; ++j) {
                    w[j] = (2. * sgn * (X[j] - R[j]) + 1.) / (R[j] * R[j]);
                    if (sgn < 0 && X[j] < 1) w[j] = 1e18;   // keep X >= 0
                }
                auto cyc = find_negative_cycle((int)tmesh.thalfs.size(), arcs, w);
                if (cyc.empty()) continue;

                VecXd g = VecXd::Zero(nte);
                for (int thid: cyc) g[tmesh.thalfs[thid].teid] += 1;
                VecXd x1 = X + sgn * g;
                double e1 = energy(x1);   // exact: a loop can hit one tedge twice
                if (e1 < e && (x1.array() >= 0).all() && compute_validation(mg, tmesh, x1)) {
                    X = x1;
                    e = e1;
                    improved = true;
                }
            }
        }

       int counter = 0;
       while (counter < 30) {
           double prev_e = e;
           std::vector<std::tuple<int, int, bool>> es;
           int l = tmesh.tedges.size();
           for (int j = 0; j < l; j++) {
               es.emplace_back(compute_weight(R[j], X[j], l, false), j, false);
               es.emplace_back(compute_weight(R[j], X[j], l, true),  j, true);
           }
           rg::sort(es.begin(), es.end(), [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });

           for (int j = 0; j < l * 2; j++) {
               auto ei = std::get<1>(es[j]);
               for (auto g: G.colwise()) {
                   if (g[ei] == 0) continue;
                   VecXd x1 = X + g;
                   VecXd x2 = X - g;
                   double n1 = (x1.cwiseQuotient(R) - I).norm();
                   double n2 = (x2.cwiseQuotient(R) - I).norm();
                   if (n1 <= e && compute_validation(mg, tmesh, x1)) { X = x1; e = n1; goto exit_loops; }
                   if (n2 <= e && compute_validation(mg, tmesh, x2)) { X = x2; e = n2; goto exit_loops; }
               }
           }
           exit_loops:
           counter = prev_e == e ? counter + 1 : 0;
       }

       std::cout << "evaluation: " << e << ", norm of diff: " << (X - R).norm() << std::endl;

       assert((C * X).norm() == 0);
       return X;
   }
}

#endif
