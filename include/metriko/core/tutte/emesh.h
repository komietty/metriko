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

    Vert tail() const { return halfs.front().tail(); }
    Vert head() const { return halfs.back().head(); }

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
        const vec<int>& reps,
        const vec<int>& ext0,
        const vec<int>& ext1
        ) {
        int side;

        { // replace side
            auto it = rg::find(ehids, ehid);
            if (it == ehids.end()) throw std::runtime_error("ehid not found");
            int idx = (int)std::distance(ehids.begin(), it);
            side = sides[idx];

            std::cout << "ehid: " << ehid << std::endl;
            std::cout << "idx:  " << idx  << std::endl;
            std::cout << "side: " << side << std::endl;
            std::cout << "reps: ";
            for (int rep: reps) { std::cout << rep << ", "; }
            std::cout << std::endl;
            auto it1 = ehids.erase(ehids.begin() + idx);
            auto it2 = sides.erase(sides.begin() + idx);
            ehids.insert_range(it1, reps);
            sides.insert_range(it2, vec(reps.size(), side));
        }

        { // insert for the prev side
            int s = (side + 3) % 4;
            auto it = std::find(sides.rbegin(), sides.rend(), s);
            int idx = sides.size() - 1 - std::distance(sides.rbegin(), it) + 1;
            ehids.insert_range(ehids.begin() + idx, ext1);
            sides.insert_range(sides.begin() + idx, vec(ext1.size(), s));
        }

        { // insert for the next side
            int s = (side + 1) % 4;
            auto it = rg::find(sides, s);
            int idx = (int)std::distance(sides.begin(), it);
            ehids.insert_range(ehids.begin() + idx, ext0);
            sides.insert_range(sides.begin() + idx, vec(ext0.size(), s));
        }

    }

    void debug_draw() const;
};

struct Emesh {
    const Hmesh& hm;
    vec<Equad> equads;
    vec<Ehalf> ehalfs;
    set<int> sings;

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
                bool end = rg_th.back().crash;
                ehalfs[thid] = Ehalf(this, rg_he, thid, tm.thalfs[thid].twid, tq.id, X[teid], bgn, end);

                // add singular vertex id
                int bgn_id = rg_th.front().half.tail().id;
                if (bgn) sings.emplace(bgn_id);
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

        std::cout << "qid: " << qid << std::endl;
        std::cout << "side: " << side << std::endl;

        if (side == -1) return;
        int bgn = -1;
        int end = -1;
        int sB = -1; // side of the bgn
        int sE = -1; // side of the end
        set<AuxDijkData> aux;

        int side_a = side == 0 ? 1 : 2;                     // remain side a
        int side_b = side == 0 ? 3 : 0;                     // remain side b
        auto ehids_p = eq.ehids_by_side(side == 0 ? 0 : 1); // collapse side p
        auto ehids_q = eq.ehids_by_side(side == 0 ? 2 : 3); // collapse side q
        auto ehids_a = eq.ehids_by_side(side_a);            // remain side a
        auto ehids_b = eq.ehids_by_side(side_b);            // remain side b

        //{
        //    // get a range of ehids_p here and use then...
        //    Ehalf& eh0 = ehalfs.front();
        //    Ehalf& eh1 = ehalfs[ehalfs.back().twid];
        //    if (eh0.bgn) { bgn = eh0.tail().id; sB = side_b; }
        //    if (eh1.bgn) { bgn = eh1.tail().id; sB = side_a; }
        //    if (bgn == -1 && eh0.end) { bgn = eh0.head().id; sB = side_a; }
        //    if (bgn == -1 && eh1.end) { bgn = eh1.head().id; sB = side_b; }
        //}

        for (int ehid: ehids_p) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) { bgn = eh0.tail().id; sB = side_b; break; }
            if (eh1.bgn) { bgn = eh1.tail().id; sB = side_a; break; }
        }
        if (bgn == -1) {
            for (int ehid: ehids_p) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.end) { bgn = eh0.head().id; sB = side_a; break; }
                if (eh1.end) { bgn = eh1.head().id; sB = side_b; break; }
            }
        }

        for (int ehid: ehids_q) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) { end = eh0.tail().id; sE = side_a; break; }
            if (eh1.bgn) { end = eh1.tail().id; sE = side_b; break; }
        }
        if (end == -1) {
            for (int ehid: ehids_q) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.end) { end = eh0.head().id; sE = side_b; break; }
                if (eh1.end) { end = eh1.head().id; sE = side_a; break; }
            }
        }

        double sum = 0;
        for (int ehid: ehids_a) sum += ehalfs[ehid].x;

        aux.emplace(AuxDijkData{hm.verts[bgn], sB, 0.});
        aux.emplace(AuxDijkData{hm.verts[end], sE, sum});

        sum = 0;
        for (int ehid: ehids_a) { Ehalf& eh = ehalfs[ehid]; sum += eh.x; aux.emplace(AuxDijkData{ eh.head(), side_a, sum }); }
        for (int ehid: ehids_b) { Ehalf& eh = ehalfs[ehid]; sum -= eh.x; aux.emplace(AuxDijkData{ eh.head(), side_b, sum }); }

        if (aux.empty()) return;


        assert(bgn != -1);
        assert(end != -1);
        assert(sum == 0);
        std::cout << "bgn: " << bgn << " end: " << end << std::endl;
        std::cout << "side a: " << side_a << ", side b: " << side_b << std::endl;

        vec<int> path_mb;

        auto visit = std::vector(hm.nH, false);
        vec aux_sorted(aux.begin(), aux.end());
        vec aux_visited(aux.size(), false);
        int N = aux.size() - 1;
        path_mb.resize(N);

        for (int i = 0; i < N; i++) {
            const auto& a0 = aux_sorted[i];
            const auto& a1 = aux_sorted[i + 1];
            const Vert v0 = a0.vert;
            const Vert v1 = a1.vert;

            for (int si: {side_a, side_b}) {
                if (a0.side != si || a1.side != si) continue;
                for (int ehid: eq.ehids_by_side(si)) {
                    Ehalf& eh0 = ehalfs[ehid];
                    Ehalf& eh1 = ehalfs[eh0.twid];
                    Vert va = eh0.tail();
                    Vert vb = eh0.head();
                    if (va.id == v0.id && vb.id  == v1.id) { aux_visited[i] = true; path_mb[i] = eh0.id; for (Half h: eh0.halfs) visit[h.id] = true; }
                    if (vb.id == v0.id && va.id  == v1.id) { aux_visited[i] = true; path_mb[i] = eh1.id; for (Half h: eh1.halfs) visit[h.id] = true; }
                }
            }
        }

        for (int i = 0; i < N; i++) {
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
            path_mb[i] = n;
            bool path0_bgn_flag = this->sings.contains(path0.front().tail().id);
            bool path0_end_flag = this->sings.contains(path0.back().head().id);
            bool path1_bgn_flag = this->sings.contains(path1.front().tail().id);
            bool path1_end_flag = this->sings.contains(path1.back().head().id);
            ehalfs.emplace_back(this, path0, n, n + 1, -1, d, path0_bgn_flag, path0_end_flag); // todo check bgn and end
            ehalfs.emplace_back(this, path1, n + 1, n, -1, d, path1_bgn_flag, path1_end_flag); // todo check bgn and end
        }

        //for (int ehid: path_mb) { ehalfs[ehid].debug_draw(); }

        // iterate over remain side a tquads
        for (int ehid: eq.ehids_by_side(side_a)) {
            Ehalf& eh0 = ehalfs[ehid];     // the curr eh
            Ehalf& eh1 = ehalfs[eh0.twid]; // the twin eh
            Equad& eq1 = equads[eh1.eqid];

            bool f0 = false;
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;
            bool f4 = false;
            vec<int> rep = {};

            for (int ehid1: path_mb) if (ehid1 == eh0.id) { f0 = true; break; }
            if (f0) continue;

            if (eh0.tail() == ehalfs[ehids_p.back()].head())  f1 = true;
            if (eh0.head() == ehalfs[ehids_q.front()].tail()) f3 = true;

            for (int ehid1: path_mb) {
                if (eh0.tail() == ehalfs[ehid1].tail()) f2 = true; // not checked yet
                if (f1 || f2) rep.emplace_back(ehid1);
                if (eh0.head() == ehalfs[ehid1].head()) { f4 = true; break; }
            }

            auto rev = rep | vw::reverse | vw::transform([&](int ehid1) { return ehalfs[ehid1].twid; }) | rg::to<vec<int>>();
            eq1.replace_ehalf(
                eh1.id, rev,
                f1 && !f2 ? ehids_p : vec<int>{},
                f3 && !f4 ? ehids_q : vec<int>{}
            );
            eq1.debug_draw();
        }

        // iterate over the remain side b tquads
        for (int ehid: eq.ehids_by_side(side_b)) {
            Ehalf& eh0 = ehalfs[ehid];
            Ehalf& eh1 = ehalfs[eh0.twid];
            Equad& eq1 = equads[eh1.eqid];

            bool f0 = false;
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;
            bool f4 = false;
            vec<int> rep = {};

            for (int ehid1: path_mb) if (ehid1 == eh1.id) { f0 = true; break; }
            if (f0) continue;

            if (eh1.tail() == ehalfs[ehids_p.front()].tail()) f1 = true;
            if (eh1.head() == ehalfs[ehids_q.back()].head())  f3 = true;

            for (int ehid1: path_mb) {
                if (eh1.tail() == ehalfs[ehid1].tail()) f2 = true; // not checked yet
                if (f1 || f2) rep.emplace_back(ehid1);
                if (eh1.head() == ehalfs[ehid1].head()) { f4 = true; break; }
            }

            eq1.replace_ehalf(
                eh1.id, rep,
                f3 && !f4 ? ehids_q : vec<int>{},
                f1 && !f2 ? ehids_p : vec<int>{}
            );
            eq1.debug_draw();
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
    std::vector<double> val3;
    std::vector<double> val4; //bgn or end flag
    size_t count = 0;

    for (int i = 0; i < ehids.size(); i++) {
        const Ehalf& eh = em->ehalfs[ehids[i]];
        const int side = sides[i];
        const bool bgn = eh.bgn;
        const bool end = eh.end;
        if (end) std::cout << "end, eqid: " << id << std::endl;
        for (int j = 0; j < eh.halfs.size(); j++) {
            const Half& h = eh.halfs[j];
            Row3d p1 = h.tail().pos();
            Row3d p2 = h.head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{count, count + 1});
            val1.emplace_back(side);
            val2.emplace_back(eh.id);
            val3.emplace_back(count);
            if      (j == 0)                   { val4.emplace_back(bgn ? -1 : 0); val4.emplace_back(0); }
            else if (j == eh.halfs.size() - 1) { val4.emplace_back(0); val4.emplace_back(end ? 1 : 0); }
            else                               { val4.emplace_back(0); val4.emplace_back(0); }
            count += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork("equad " + std::to_string(id), ns, es);
    c->addEdgeScalarQuantity("side", val1);
    c->addEdgeScalarQuantity("ehids", val2);
    c->addEdgeScalarQuantity("count", val3);
    c->addNodeScalarQuantity("bgn end", val4);
    c->resetTransform();
    c->setRadius(0.002);
}
}

#endif
