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

    void debug_draw() const;
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

    void replace_ehalf(
        int ehid,
        const vec<int>& ext0,
        const vec<int>& reps,
        const vec<int>& ext1
        ) {
        auto it = rg::find(ehids, ehid);
        if (it == ehids.end()) throw std::runtime_error("side not found");
        int idx  = (int)std::distance(ehids.begin(), it);
        int side = sides[idx];
        std::cout << "ehid: " << ehid << std::endl;
        std::cout << "idx: " << idx << std::endl;
        std::cout << "side: " << side << std::endl;
        std::cout << "reps: ";
        for (int rep: reps) { std::cout << rep << ", " << std::endl; }
        std::cout << std::endl;
        ehids.erase(ehids.begin() + idx);
        sides.erase(sides.begin() + idx);
        for (int rep: reps) {
            ehids.emplace_back(rep);
            sides.emplace_back(side);
        }
    }

    void extend_ehalf() {}

    void debug_draw() const;
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

    void collapse_quad(const int qid) {
        int side = -1; // if 0 or 1, collapse
        const auto& eq = equads[qid];
        for (int i = 0; i < 2; i++) {
            int sum = (int)rg::fold_left(eq.ehids_by_side(i), 0, [&](int acc, int ehid) { return acc + ehalfs[ehid].x; });
            if (sum == 1) side = i; // right now not 0 but 1
        }

        if (side == -1) return;
        std::cout << qid << std::endl;
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
        std::cout << "bgn: " << bgn << " end: " << end << std::endl;
        std::cout << "side a: " << remain_side_a << ", side b: " << remain_side_b << std::endl;

        vec<int> path_fr;
        vec<int> path_to;
        vec<int> path_md;

        { // find path_fr (ehids)
            Vert v0 = ehalfs[ehids_p.front()].halfs.front().tail();
            Vert v1 = ehalfs[ehids_p.back()].halfs.back().head();
            if (v1.id == bgn) { for (int ehid: ehids_p) path_fr.emplace_back(ehid); }
            if (v0.id == bgn) { for (int ehid: ehids_p) path_fr.emplace_back(ehalfs[ehid].twid); }
        }

        { // find path_to (ehids)
            Vert v0 = ehalfs[ehids_q.front()].halfs.front().tail();
            Vert v1 = ehalfs[ehids_q.back()].halfs.back().head();
            if (v0.id == bgn) { for (int ehid: ehids_q) path_to.emplace_back(ehid); }
            if (v1.id == bgn) { for (int ehid: ehids_q) path_to.emplace_back(ehalfs[ehid].twid); }
        }

        auto visit = std::vector(hm.nH, false);
        vec aux_sorted(aux.begin(), aux.end());
        vec aux_visited(aux.size(), false);

        for (int i = 0; i < aux_sorted.size() - 1; i++) {
            const auto& a0 = aux_sorted[i];
            const auto& a1 = aux_sorted[i + 1];
            const Vert v0 = a0.vert;
            const Vert v1 = a1.vert;

            for (int si: {remain_side_a, remain_side_b}) {
                if (a0.side != si || a1.side != si) continue;
                for (int ehid: eq.ehids_by_side(si)) {
                    Ehalf& eh0 = ehalfs[ehid];
                    Ehalf& eh1 = ehalfs[eh0.twid];
                    Vert va = eh0.halfs.front().tail();
                    Vert vb = eh0.halfs.back().head();
                    if (va.id == v0.id && vb.id  == v1.id) { aux_visited[i] = true; path_md.emplace_back(eh0.id); for (Half h: eh0.halfs) visit[h.id] = true; }
                    if (vb.id == v0.id && va.id  == v1.id) { aux_visited[i] = true; path_md.emplace_back(eh1.id); for (Half h: eh1.halfs) visit[h.id] = true; }
                }
            }
        }

        for (int i = 0; i < aux_sorted.size() - 1; i++) {
            if (aux_visited[i]) continue;
            const auto& a0 = aux_sorted[i];
            const auto& a1 = aux_sorted[i + 1];
            const Vert v0 = a0.vert;
            const Vert v1 = a1.vert;
            const double d = std::abs(a0.value - a1.value);

            auto path = compute_dijkstra(hm, visit, v0, v1);
            vec<Half> path0;
            vec<Half> path1;
            for (int hid: path) {
                visit[hid] = true;
                path0.emplace_back(hm.halfs[hid]);
                path1.emplace_back(hm.halfs[hid].twin());
            }

            rg::reverse(path1);

            int n = ehalfs.size();
            path_md.emplace_back(n);
            ehalfs.emplace_back(Ehalf(this, path0, n, n + 1, -1, d));
            ehalfs.emplace_back(Ehalf(this, path1, n + 1, n, -1, d));
        }

        for (int ehid: path_md) {
            ehalfs[ehid].debug_draw();
        }

        // iterate over the remain side a tquads
        for (int ehid: eq.ehids_by_side(remain_side_a)) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            Equad& eq  = equads[eh1.eqid];

            // if the ehalf still alive, that is okay
            bool f0 = false;
            for (int ehid1: path_md)  if (ehid1 == eh0.id) { f0 = true; break; }
            if (f0) continue;

            //if not, iterate over full path from the right most vert, until finding left most vert
            bool f1 = false;
            vec<int> rep = {};
            // todo: if the vert equals to path_fr's tail, f1 is true
            for (int ehid1: path_md) {
                if (eh0.halfs.front().tail() == ehalfs[ehid1].halfs.front().tail()) f1 = true;
                if (f1) rep.emplace_back(ehid1);
                if (eh0.halfs.back().head() == ehalfs[ehid1].halfs.back().head()) break;
            }

            eq.replace_ehalf(eh1.id, {}, rep, {});
            eq.debug_draw();
        }

        // iterate over the remain side b tquads
        for (int ehid: eq.ehids_by_side(remain_side_b)) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            Equad& eq  = equads[eh1.eqid];
            // if the ehalf still alive, that is okay
            bool f0 = false;
            for (int ehid1: path_md)  if (ehid1 == eh0.id) {f0 = true; break;}
            if (f0) continue;

            //if not, iterate over full path from the right most vert, until finding left most vert
            bool f1 = false;
            vec<int> rep = {};
            for (int ehid1: path_md) {
                if (eh0.halfs.front().tail() == ehalfs[ehid1].halfs.front().tail()) f1 = true;
                if (f1) rep.emplace_back(ehid1);
                if (eh0.halfs.back().head() == ehalfs[ehid1].halfs.back().head()) break;
            }
        }
    }

};

inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }

inline void Ehalf::debug_draw() const {
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> val1; //
    size_t count = 0;

    for (const Half h: halfs) {
        Row3d p1 = h.tail().pos();
        Row3d p2 = h.head().pos();
        ns.emplace_back(p1.x(), p1.y(), p1.z());
        ns.emplace_back(p2.x(), p2.y(), p2.z());
        es.emplace_back(std::array{count, count + 1});
        val1.emplace_back(count / 2);
        count += 2;
    }
    auto c = polyscope::registerCurveNetwork("ehalf " + std::to_string(id), ns, es);
    c->addEdgeScalarQuantity("order", val1);
    c->resetTransform();
    c->setRadius(0.002);
}

inline void Equad::debug_draw() const {
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> val1;
    std::vector<double> val2;
    size_t count = 0;

    for (int i = 0; i < ehids.size(); i++) {
        const Ehalf& eh = em->ehalfs[ehids[i]];
        const int side = sides[i];
        for (const Half h: eh.halfs) {
            Row3d p1 = h.tail().pos();
            Row3d p2 = h.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{count, count + 1});
            val1.emplace_back(side);
            val2.emplace_back(eh.id);
            count += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork("equad " + std::to_string(id), ns, es);
    c->addEdgeScalarQuantity("order", val1);
    c->addEdgeScalarQuantity("ehids", val2);
    c->resetTransform();
    c->setRadius(0.002);
}
}

#endif
