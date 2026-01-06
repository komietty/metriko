#ifndef TUTTE_H_EMESH_H
#define TUTTE_H_EMESH_H

#include "metriko/core/hmesh/hmesh.h"
#include "compute_dijkstra.h"
#include "tutte_cutting.h"

namespace metriko::tutte {
template <class T> using vec = std::vector<T>;
template <class T> using set = std::set<T>;

struct HalfData {
    Half half;
    double v0;
    double v1;
    int thid;
    int tqid;
    int twin;
    int order; // order inside thalf
    bool first = false;
    bool crash = false;

    bool operator<(const HalfData& rhs) const noexcept {
        return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order);
    }
};

struct AuxDijkData {
    Vert vert;
    int side;
    double value;
    bool operator<(const AuxDijkData& rhs) const noexcept { return value < rhs.value; }
};


struct Emesh;

struct Ehalf {
    const Emesh* em = nullptr;
    vec<Half> halfs;
    int id = -1;
    int twid = -1;
    int eqid = -1;
    bool bgn = false; // the beginning halfedge in terms of motorcycle graph
    bool end = false; // the last halfedge in terms of motorcycle graph
    double x = -1;

    Ehalf() = default;

    Ehalf(
        const Emesh* em,
        const vec<Half>& halfs,
        const int id,
        const int twin,
        const int eqid,
        const double x,
        const bool bgn = false,
        const bool end = false
    ): em(em), halfs(halfs), id(id), twid(twin), eqid(eqid), bgn(bgn), end(end), x(x) {}

    const Ehalf& twin() const;
};

struct Equad {
    const Emesh* em = nullptr;
    int id = -1;
    vec<int> ehids;
    vec<int> sides;
    Equad() = default;

    Equad(
        const Emesh* em,
        const int id,
        const vec<int>& ehids,
        const vec<int>& sides
    ): em(em), id(id), ehids(ehids), sides(sides) {}

    vec<int> ehids_by_side(int side, bool reverse = false) const {
        vec<int> res;
        for (int i = 0; i < ehids.size(); i++) {
            int j = reverse ? (int)ehids.size() - i - 1 : i;
            if (sides[j] == side) res.emplace_back(ehids[j]);
        }
        return res;
    }

    void replace_side() { }
};

struct Emesh {
    const Hmesh& hm;
    vec<Equad> equads;
    vec<Ehalf> ehalfs;

    explicit Emesh(
        const Hmesh& hm,
        const tm::Tmesh& tm,
        const set<HalfData>& hdata,
        const VecXd& X,
        const VecXd& R
    ): hm(hm) {
        equads.resize(tm.tquads.size());
        ehalfs.resize(tm.thalfs.size());

        for (const auto& tq: tm.tquads) {
            equads[tq.id] = Equad(this, tq.id, tq.thids, tq.sides);
            auto rg_tq = rg::equal_range(hdata, tq.id, {}, &HalfData::tqid);
            for (int thid: tq.thids) {
                int teid = tm.thalfs[thid].teid;
                auto rg_th = rg_tq | vw::filter([&](auto& hd) { return thid == hd.thid; }) | rg::to<vec<HalfData>>();
                auto rg_he = rg_th | vw::transform([](auto& hd) { return hd.half; }) | rg::to<vec<Half>>();
                bool bgn = rg_th.front().first;
                bool end = rg_th.front().crash;
                ehalfs[thid] = Ehalf(this, rg_he, thid, tm.thalfs[thid].twid, tq.id, X[teid], bgn, end);
            }
        }
    }

    void collapse_quad(const int qid, const double r) {
        int side = -1; // if 0 or 1, collapse
        const auto& eq = equads[qid];
        for (int i = 0; i < 2; i++) {
            int sum = rg::fold_left(eq.ehids_by_side(i), 0, [&](int acc, int ehid) { return acc + ehalfs[ehid].x; });
            if (sum == 1) side = i; // right now not 0 but 1
        }

        if (side == -1) return;
        int bgn = -1;
        int end = -1;
        int bgn_flg = false;
        int end_flg = false;
        set<AuxDijkData> aux;

        int remain_side_a = side == 0 ? 1 : 2;
        int remain_side_b = side == 0 ? 3 : 0;
        auto ehids_p = eq.ehids_by_side(side == 0 ? 0 : 1); // collapse side p
        auto ehids_q = eq.ehids_by_side(side == 0 ? 2 : 3); // collapse side q
        auto ehids_a = eq.ehids_by_side(remain_side_a);     // remain side a
        auto ehids_b = eq.ehids_by_side(remain_side_b);     // remain side b

        for (int ehid: ehids_p) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) { bgn = eh0.halfs.front().tail().id; bgn_flg = true; break; }
            if (eh1.bgn) { bgn = eh1.halfs.front().tail().id; break; }
        }

        if (bgn == -1) {
            for (int ehid: ehids_p) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.end) { bgn = eh0.halfs.back().head().id; break; }
                if (eh1.end) { bgn = eh1.halfs.back().head().id; break; }
            }
        }

        for (int ehid: ehids_q) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) { end = eh0.halfs.front().tail().id; end_flg = true; break; }
            if (eh1.bgn) { end = eh0.halfs.back().head().id; break; }
        }

        if (end == -1) {
            for (int ehid: ehids_q) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.end) { end = eh0.halfs.back().head().id; break; }
                if (eh1.end) { end = eh1.halfs.back().head().id; break; }
                //end_flg = hd0.order == 1; // needed???
            }
        }

        double sum = 0;
        for (int ehid: ehids_a) sum += ehalfs[ehid].x;

        int side_bgn = bgn_flg ? remain_side_b : remain_side_a;
        int side_end = end_flg ? remain_side_a : remain_side_b;
        aux.emplace(AuxDijkData{hm.verts[bgn], side_bgn, 0.});
        aux.emplace(AuxDijkData{hm.verts[end], side_end, sum});

        sum = 0;
        for (int ehid: ehids_a) { Ehalf& eh = ehalfs[ehid]; sum += eh.x; aux.emplace(AuxDijkData{ eh.halfs.back().head(), remain_side_a, sum }); }
        for (int ehid: ehids_b) { Ehalf& eh = ehalfs[ehid]; sum -= eh.x; aux.emplace(AuxDijkData{ eh.halfs.back().head(), remain_side_b, sum }); }

        if (aux.empty()) return;

        assert(sum == 0);
        assert(bgn != -1);
        assert(end != -1);

        std::vector<int> full_path;
        auto visit = std::vector(hm.nH, false);
        std::vector aux_sorted(aux.begin(), aux.end());

        for (int i = 0; i < aux_sorted.size() - 1; i++) {
            const auto& a0 = aux_sorted[i];
            const auto& a1 = aux_sorted[i + 1];
            const double d = std::abs(a0.value - a1.value);

            if (a0.side == remain_side_a && a1.side == remain_side_a) {
                for (int ehid: eq.ehids_by_side(remain_side_a)) {
                    Ehalf& eh = ehalfs[ehid];
                    if (eh.halfs.front().tail().id == a0.vert.id &&
                        eh.halfs.back().head().id  == a1.vert.id)
                        full_path.emplace_back(ehid);
                }
            } else if (a0.side == remain_side_b && a1.side == remain_side_b) {
                for (int ehid: eq.ehids_by_side(remain_side_b)) {
                    Ehalf& eh = ehalfs[ehid];
                    if (eh.halfs.back().tail().id == a0.vert.id &&
                        eh.halfs.front().head().id  == a1.vert.id) // should be rev. must be okay...
                        full_path.emplace_back(ehid);
                }
            } else {
                auto path = compute_dijkstra(hm, visit, a0.vert, a1.vert);
                vec<Half> path0;
                vec<Half> path1;
                for (int hid: path) {
                    visit[hid] = true;
                    path0.emplace_back(hm.halfs[hid]);
                    path1.emplace_back(hm.halfs[hid].twin());
                }

                rg::reverse(path1);

                int n = ehalfs.size();
                full_path.emplace_back(n);
                ehalfs.emplace_back(Ehalf(this, path0, n, n + 1, -1, d));
                ehalfs.emplace_back(Ehalf(this, path1, n + 1, n, -1, d));
            }
        }

        // find adjacent equads
        for (int ehid: eq.ehids_by_side(remain_side_a)) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            Equad& eq  = equads[eh1.eqid];
            // if the ehalf still alive, that is okay
            for (int ehid1: full_path) { if (ehid1 == eh0.id) break; }

            //if not, iterate over full path from the right most vert, until finding left most vert
            bool f = false;
            for (int ehid1: full_path) {
                if (eh0.halfs.front().tail() == ehalfs[ehid1].halfs.front().tail()) f = true;
                if (f) eq.replace_side();
                if (eh0.halfs.back().head() == ehalfs[ehid1].halfs.back().head()) break;
            }
        }

        for (int ehid: eq.ehids_by_side(remain_side_b)) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            Equad& eq  = equads[eh1.eqid];
            // if the ehalf still alive, that is okay
            for (int ehid1: full_path) { if (ehid1 == eh0.id) break; }

            //if not, iterate over full path from the right most vert, until finding left most vert
            bool f = false;
            for (int ehid1: full_path) {
                if (eh0.halfs.front().tail() == ehalfs[ehid1].halfs.front().tail()) f = true;
                if (f) eq.replace_side();
                if (eh0.halfs.back().head() == ehalfs[ehid1].halfs.back().head()) break;
            }
        }
    }
};

inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }
}

#endif
