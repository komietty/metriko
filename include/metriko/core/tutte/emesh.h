#ifndef TUTTE_H_EMESH_H
#define TUTTE_H_EMESH_H

#include "metriko/core/hmesh/hmesh.h"
#include "compute_dijkstra.h"

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
    double val;
    double cmp;

    // C++20のコンパイラを納得させるための「完全なルール」
    bool operator<(const AuxDijkData& rhs) const noexcept {
        if (val != rhs.val) return val < rhs.val;
        if (cmp != rhs.cmp) return cmp < rhs.cmp;
        if (vert.id != rhs.vert.id) return vert.id < rhs.vert.id;
        return side < rhs.side;
    }

    bool operator==(const AuxDijkData& rhs) const noexcept {
        return val == rhs.val && cmp == rhs.cmp && vert.id == rhs.vert.id && side == rhs.side;
    }
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

    void extend_prev(Ehalf prev) { halfs.insert(halfs.begin(), prev.halfs.begin(), prev.halfs.end()); }
    void extend_next(Ehalf next) { halfs.insert(halfs.end()  , next.halfs.begin(), next.halfs.end()); }

    void debug_draw() const;
};

struct Edata {
    int ehid;
    int side;
};

struct Equad {
    const Emesh* em = nullptr;
    int id = -1;
    vec<Edata> edata;

    Equad() = default;
    Equad(
        const Emesh* em,
        const int id,
        const vec<int>& ehids,
        const vec<int>& sides
    ): em(em), id(id) {
        edata.reserve(ehids.size());
        for (auto [ehid, side] : vw::zip(ehids, sides)) edata.emplace_back(ehid, side);
    }

    vec<int> ehids_by_side(int side, bool reverse = false) const {
        vec<int> res;
        for (int i = 0; i < edata.size(); i++) {
            int j = reverse ? (int)edata.size() - i - 1 : i;
            if (edata[j].side == side) res.emplace_back(edata[j].ehid);
        }
        return res;
    }

    vec<int> verts_inside() const;
    void replace_ehalf(int ehid, const vec<int>& reps, const vec<int>& ext0, const vec<int>& ext1);
    void debug_draw() const;
};

struct Emesh {
    const Hmesh& hm;
    vec<Equad> equads;
    vec<Ehalf> ehalfs;
    set<int> sings;

    Emesh(
        const Emesh& old_em,
        const Hmesh& new_hm
    ): hm(new_hm) {
        equads = old_em.equads;
        ehalfs = old_em.ehalfs;
        sings  = old_em.sings;
        for (auto& eq : equads) { eq.em = this; }
        for (auto& eh : ehalfs) { eh.em = this; }
    }

    explicit Emesh(
        const Hmesh& hm,
        const tm::Tmesh& tm,
        const set<HalfData>& hdata,
        const VecXd& X
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

    vec<Half> collapse_half_find_path(int ehid);
    bool collapse_half(int ehid);

    bool collapse_quad(const int qid) {
        int side = -1;
        const auto& eq = equads[qid];
        for (int i : {0, 1}) {
            auto side_ehids = eq.ehids_by_side(i);
            int sum = rg::fold_left(side_ehids | vw::transform([this](int id){ return ehalfs[id].x; }), 0., std::plus());
            if (sum == 0) { side = i; break; }
        }

        if (side == -1) return false;
        std::cout << "collapse eqid: " << qid << std::endl;

        int side_a = side == 0 ? 1 : 2;                     // remain side a
        int side_b = side == 0 ? 3 : 0;                     // remain side b
        auto ehids_p = eq.ehids_by_side(side == 0 ? 0 : 1); // collapse side p
        auto ehids_q = eq.ehids_by_side(side == 0 ? 2 : 3); // collapse side q
        auto ehids_a = eq.ehids_by_side(side_a);            // remain side a
        auto ehids_b = eq.ehids_by_side(side_b);            // remain side b

        auto find_terminal = [&](const vec<int>& ehids, int s0, int s1) -> std::pair<int, int> {
            for (int ehid: ehids) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.bgn) return {eh0.tail().id, s0};
                if (eh1.bgn) return {eh1.tail().id, s1};
            }
            for (int ehid: ehids) {
                Ehalf& eh0 = ehalfs[ehid];
                Ehalf& eh1 = ehalfs[eh0.twid];
                if (eh0.end) return {eh0.head().id, s1};
                if (eh1.end) return {eh1.head().id, s0};
            }
            throw std::runtime_error("terminal not found");
        };

        auto [bgn, sB] = find_terminal(ehids_p, side_b, side_a); // vert and side of the bgn
        auto [end, sE] = find_terminal(ehids_q, side_a, side_b); // vert and side of the end

        double sum = 0;
        double cmp = 0;
        for (int ehid: ehids_a) {
            sum += ehalfs[ehid].x;
            cmp += ehalfs[ehid].x + 1;
        }

        vec<AuxDijkData> aux;
        aux.emplace_back(AuxDijkData{hm.verts[bgn], sB, 0., 0.});
        std::cout << "aux first vert id: " << hm.verts[bgn].id << ", aux first side id: " << sB << std::endl;
        aux.emplace_back(AuxDijkData{hm.verts[end], sE, sum, cmp});
        std::cout << "aux last vert id: " << hm.verts[end].id << ", aux last side id: " << sE << std::endl;

        double sum1 = 0;
        double cmp1 = 0;
        for (int ehid: ehids_a) {
            Ehalf& eh = ehalfs[ehid];
            sum1 += eh.x;
            cmp1 += eh.x + 1;
            if (sum1 != 0 && sum1 != sum)
                aux.emplace_back(AuxDijkData{ eh.head(), side_a, sum1, cmp1});
            std::cout << "aux side a vert id: " << eh.head().id << ", aux side a side id: " << side_a << std::endl;
        }
        for (int ehid: ehids_b) {
            Ehalf& eh = ehalfs[ehid];
            sum1 -= eh.x;
            cmp1 -= eh.x + 1;
            if (sum1 != 0 && sum1 != sum)
                aux.emplace_back(AuxDijkData{ eh.head(), side_b, sum1, cmp1});
            std::cout << "aux side b vert id: " << eh.head().id << ", aux side b side id: " << side_b << ", sum: " << sum << std::endl;
        }


        if (aux.empty()) return false;
        assert(sum1 == 0);

        rg::sort(aux, [](const AuxDijkData& a, const AuxDijkData& b) {
            return std::tie(a.val, a.cmp) < std::tie(b.val, b.cmp);
        });

        // 2. val, vert.id, side が同じものを重複とみなして消す
        auto ret = rg::unique(aux, [](const AuxDijkData& a, const AuxDijkData& b) {
            return std::tie(a.val, a.vert.id, a.side) == std::tie(b.val, b.vert.id, b.side);
        });
        aux.erase(ret.begin(), ret.end());

        vec<int> path_mb;

        auto visit = vec(hm.nH, false);
        vec aux_sorted(aux.begin(), aux.end());
        vec aux_visited(aux.size(), false);
        int N = aux.size() - 1;
        path_mb.resize(N);

        std::vector<glm::vec3> verts_to_passby;
        for (const AuxDijkData &a: aux_sorted) {
            Row3d p = a.vert.pos();
            verts_to_passby.emplace_back(p.x(), p.y(), p.z());
        }
        auto vq = polyscope::registerPointCloud("verts to pass by", verts_to_passby);
        vq->setEnabled(true);
        vq->setPointRadius(0.002);
        vq->resetTransform();
        //return false;


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
                    if (va.id == v0.id && vb.id  == v1.id) {
                        std::cout << "i: passed: " << i << std::endl;
                        aux_visited[i] = true; path_mb[i] = eh0.id; for (Half h: eh0.halfs) visit[h.id] = true;
                    }
                    if (vb.id == v0.id && va.id  == v1.id) {
                        std::cout << "i: passed: " << i << std::endl;
                        aux_visited[i] = true; path_mb[i] = eh1.id; for (Half h: eh1.halfs) visit[h.id] = true;
                    }
                }
            }
        }



        for (int i = 0; i < N; i++) {
            std::cout << "i: check: " << i << std::endl;
            if (aux_visited[i]) continue;
            const auto& a0 = aux_sorted[i];
            const auto& a1 = aux_sorted[i + 1];
            const Vert v0 = a0.vert;
            const Vert v1 = a1.vert;
            const auto d = std::abs(a0.val - a1.val);

            auto path = compute_dijkstra(hm, visit, v0, v1);
            vec<Half> path0;
            vec<Half> path1;
            for (int hid: path) {
                visit[hid] = true;
                path0.emplace_back(hm.halfs[hid]);
                path1.emplace_back(hm.halfs[hid].twin());
            }
            std::cout << "path0 len: " << path0.size() << std::endl;

            rg::reverse(path1);

            int n = ehalfs.size();
            path_mb[i] = n;
            bool bgn0 = this->sings.contains(path0.front().tail().id);
            bool end0 = this->sings.contains(path0.back().head().id);
            bool bgn1 = this->sings.contains(path1.front().tail().id);
            bool end1 = this->sings.contains(path1.back().head().id);
            ehalfs.emplace_back(this, path0, n, n + 1, -1, d, bgn0, end0);
            ehalfs.emplace_back(this, path1, n + 1, n, -1, d, bgn1, end1);
        }


        // debug visit
        /*
        {
            std::vector<glm::vec3> ns;
            std::vector<std::array<size_t, 2>> es;
            std::vector<double> val1;
            size_t count = 0;

            for (int i = 0; i < hm.nH; i++) {
                if (!visit[i]) { continue; }
                Half h = hm.halfs[i];
                Row3d p1 = h.tail().pos();
                Row3d p2 = h.head().pos();
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{count, count + 1});
                val1.emplace_back(count / 2);
                count += 2;
            }
            auto c = polyscope::registerCurveNetwork("visit half for eqid: " + std::to_string(qid), ns, es);
            c->addEdgeScalarQuantity("order", val1);
            c->resetTransform();
            c->setRadius(0.002);
        }
        */

        auto path_md = path_mb | vw::transform([&](int ehid) { return ehalfs[ehid].twid; }) | rg::to<vec<int>>();
        rg::reverse(path_md);

        {
            std::vector<glm::vec3> ns;
            std::vector<std::array<size_t, 2>> es;
            std::vector<double> val1; //
            size_t count = 0;

            for (int ehid: path_md) {
                auto eh = ehalfs[ehid];
                for (Half h: eh.halfs) {
                    Row3d p1 = h.tail().pos();
                    Row3d p2 = h.head().pos();
                    ns.emplace_back(p1.x(), p1.y(), p1.z());
                    ns.emplace_back(p2.x(), p2.y(), p2.z());
                    es.emplace_back(std::array{count, count + 1});
                    val1.emplace_back(count / 2);
                    count += 2;
                }
            }

            auto c = polyscope::registerCurveNetwork("path_md of eqid" + std::to_string(qid), ns, es);
            c->addEdgeScalarQuantity("order", val1);
            c->resetTransform();
            c->setRadius(0.002);
        }
        /*
        */

        std::cout << "path_mb" << std::endl;
        replace_remain_side(eq, path_mb, ehids_p, ehids_q, side_a);
        std::cout << "path_md" << std::endl;
        replace_remain_side(eq, path_md, ehids_q, ehids_p, side_b);
        return true;
    }

    void replace_remain_side(
        const Equad &eq,
        const vec<int> &ehids_path,
        const vec<int> &ehids_prev,
        const vec<int> &ehids_next,
        int side
    ) {
        for (int ehid: eq.ehids_by_side(side)) {
            Ehalf& currEH = ehalfs[ehid];
            Ehalf& twinEH = ehalfs[currEH.twid];

            if (rg::find(ehids_path, currEH.id) != ehids_path.end()) continue;

            auto rp = vec<int>{};
            bool f1 = currEH.tail() == ehalfs[ehids_prev.back()].head();
            bool f3 = currEH.head() == ehalfs[ehids_next.front()].tail();
            bool f2 = false;
            bool f4 = false;

            for (int ehid_: ehids_path) {
                Ehalf& eh = ehalfs[ehid_];
                if (currEH.tail() == eh.tail()) f2 = true;
                if (f1 || f2) rp.emplace_back(eh.twid);
                if (currEH.head() == eh.head()) { f4 = true; break; }
            }

            rg::reverse(rp);
            std::cout << "eqid to update: " << twinEH.eqid << std::endl;

            equads[twinEH.eqid].replace_ehalf(
                twinEH.id, rp,
                f1 && !f2 ? ehids_prev : vec<int>{},
                f3 && !f4 ? ehids_next : vec<int>{}
            );
            //eq1.debug_draw();
        }
    }
};


inline const Ehalf& Ehalf::twin() const { return em->ehalfs[twid]; }

inline vec<int> Equad::verts_inside() const {
    vec<int> res;
    const Hmesh& hm = em->hm;

    vec is_boundary(hm.nH, false);
    vec visited_face(hm.nF, false);
    std::queue<int> q;

    for (auto [ehid, side] : edata) {
        const Ehalf& eh = em->ehalfs[ehid];
        for (Half h : eh.halfs) {
            is_boundary[h.id] = true;

            int fid = h.face().id;
            if (fid != -1 && !visited_face[fid]) {
                visited_face[fid] = true;
                q.push(fid);
            }
        }
    }

    vec<int> inner_faces;
    while (!q.empty()) {
        int curr_fid = q.front();
        q.pop();
        inner_faces.push_back(curr_fid);

        Face f = hm.faces[curr_fid];
        for (Half he: f.adjHalfs()) {
            if (is_boundary[he.id]) continue;

            int next_fid = he.twin().face().id;
            if (next_fid != -1 && !visited_face[next_fid]) {
                visited_face[next_fid] = true;
                q.push(next_fid);
            }
        }
    }

    vec visited_vert(hm.nV, false);
    for (int fid : inner_faces) {
        Face f = hm.faces[fid];
        for (Half he : f.adjHalfs()) {
            Vert v = he.tail();
            if (!visited_vert[v.id]) {
                visited_vert[v.id] = true;
                res.push_back(v.id);
            }
        }
    }

    return res;
}

inline void Equad::replace_ehalf(
    int ehid,             // tgt ehalf id
    const vec<int>& reps, // ehids for replacing ehalf above
    const vec<int>& ext0, // ehids for extending next side of ehalfs
    const vec<int>& ext1  // ehids for extending prev side of ehalfs
) {
    int side;

    { // replace side
        auto it = rg::find(edata, ehid, &Edata::ehid);
        assert(it != edata.end());
        side = it->side;
        auto addr = edata.erase(it);
        auto item = reps | vw::transform([side](int r) { return Edata{r, side}; });
        edata.insert_range(addr, item);
    }

    auto remove_from_twin_equad = [&](int twin_eqid, int twin_ehid) {
        auto& eq = const_cast<Equad&>(em->equads[twin_eqid]);
        auto it = rg::find(eq.edata, twin_ehid, &Edata::ehid);
        if (it != eq.edata.end()) eq.edata.erase(it);
    };

    { // insert for the prev side
        auto rev_edata = edata | vw::reverse;
        auto it = rg::find(rev_edata, (side + 3) % 4, &Edata::side);
        assert(it != rev_edata.end());
        Ehalf& eh = const_cast<Ehalf&>(em->ehalfs[it->ehid]);
        Ehalf& tw = const_cast<Ehalf&>(em->ehalfs[eh.twid]);

        for (int ehid_: ext1) {
            const Ehalf& ext_eh = em->ehalfs[ehid_];
            const Ehalf& ext_tw = em->ehalfs[ext_eh.twid];
            eh.extend_next(ext_eh);
            tw.extend_prev(ext_tw);
            remove_from_twin_equad(ext_tw.eqid, ext_tw.id);
        }
    }

    { // insert for the next side
        auto it = rg::find(edata, (side + 1) % 4, &Edata::side);
        assert(it != edata.end());
        Ehalf& eh = const_cast<Ehalf&>(em->ehalfs[it->ehid]);
        Ehalf& tw = const_cast<Ehalf&>(em->ehalfs[eh.twid]);
        for (int ehid_: std::views::reverse(ext0)) {
            const Ehalf& ext_eh = em->ehalfs[ehid_];
            const Ehalf& ext_tw = em->ehalfs[ext_eh.twid];
            eh.extend_prev(ext_eh);
            tw.extend_next(ext_tw);
            remove_from_twin_equad(ext_tw.eqid, ext_tw.id);
        }
    }
}

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
    std::vector<double> val5;
    std::vector<double> x; //x
    size_t count = 0;

    for (int i = 0; i < edata.size(); i++) {
        const Ehalf& eh = em->ehalfs[edata[i].ehid];
        const int side = edata[i].side;
        const bool bgn = eh.bgn;
        const bool end = eh.end;
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
            val5.emplace_back(id);
            x.emplace_back(eh.x);
            if      (j == 0)                   { val4.emplace_back(bgn ? -1 : 0); val4.emplace_back(0); }
            else if (j == eh.halfs.size() - 1) { val4.emplace_back(0); val4.emplace_back(end ? 1 : 0); }
            else                               { val4.emplace_back(0); val4.emplace_back(0); }
            count += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork("equad " + std::to_string(id), ns, es);
    c->addEdgeScalarQuantity("side", val1);
    c->addEdgeScalarQuantity("eqid", val5);
    c->addEdgeScalarQuantity("ehids", val2);
    c->addEdgeScalarQuantity("count", val3);
    c->addNodeScalarQuantity("bgn end", val4);
    c->addEdgeScalarQuantity("x", x)->setEnabled(true);
    c->setEnabled(false);
    c->resetTransform();
    c->setRadius(0.002);
}
}

#endif
