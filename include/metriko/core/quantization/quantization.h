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

   inline VecXd compute_quantization(
       const Tmesh& tmesh,
       const VecXd& R
   ) {
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
           if (compute_validation(tmesh, X)) {
               std::cout << "validation passed" << std::endl;
               break;
           }
       }

       // ----- construct second step vector ----- //

       double e = (X.cwiseQuotient(R) - I).norm();

       int counter = 0;
       while (counter < 10) {
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
                   if (n1 <= e && compute_validation(tmesh, x1)) { X = x1; e = n1; goto exit_loops; }
                   if (n2 <= e && compute_validation(tmesh, x2)) { X = x2; e = n2; goto exit_loops; }
               }
           }
           exit_loops:
           counter = prev_e == e ? counter + 1 : 0;
       }

       std::cout << "evaluation: " << e << ", norm of diff: " << (X - R).norm() << std::endl;

       // ----- construct second step vector (trying another basis) ----- //
       MatXd G2 = construct_generating_vectors(
           tmesh,
           R,
           [](const Comparator &c1, const Comparator &c2) { return c1.length < c2.length; }
       );

       for (auto &th: tmesh.thalfs) {
           for (auto g: G2.colwise()) {
               if (g[th.edge().id] == 0) continue;

               VecXd x1 = X + g;
               VecXd x2 = X - g;
               double n1 = (x1.cwiseQuotient(R) - I).norm();
               double n2 = (x2.cwiseQuotient(R) - I).norm();

               if (n1 <= e && (x1.array() >= 0.).all() && compute_validation(tmesh, x1)) {
                   X = x1;
                   e = n1;
               }

               if (n2 <= e && (x2.array() >= 0.).all() && compute_validation(tmesh, x2)) {
                   X = x2;
                   e = n2;
               }
           }
       }

       std::cout << "evaluation: " << e << ", norm of diff: " << (X - R).norm() << std::endl;

       assert((C * X).norm() == 0);
       return X;
   }
}

#endif
