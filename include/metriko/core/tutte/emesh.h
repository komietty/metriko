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

    Emesh(
        const Hmesh& hm,
        const Tmesh& tm,
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

    bool collapse_half(int ehid);
    bool collapse_quad(int eqid);

    void replace_remain_side(
        const Equad &eq,
        const vec<int> &ehids_path,
        const vec<int> &ehids_prev,
        const vec<int> &ehids_next,
        int side
    );

    bool check_topology() const {
        for (const Ehalf& currEH: ehalfs) {
            auto& twinEH = ehalfs[currEH.twid];
            assert(twinEH.halfs.size() == currEH.halfs.size());
            int l = (int)currEH.halfs.size();
            for (int i = 0; i < l; i++) {
                Half currH = currEH.halfs[i];
                Half twinH = twinEH.halfs[l - i - 1];
                assert(currH.twin() == twinH);
            }
        }
        return true;
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


    std::vector<glm::vec3> pts;
    for (int vid: res) {
        Vert v = hm.verts[vid];
        pts.emplace_back(v.pos().x(), v.pos().y(), v.pos().z());
    }
    auto p = polyscope::registerPointCloud("inside pts of quad: " + std::to_string(id), pts);
    p->setEnabled(false);
    p->setPointRadius(0.0005);
    p->resetTransform();

    return res;
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

    for (auto [ehid, side] : edata) {
        const Ehalf& eh = em->ehalfs[ehid];
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
    //c->addEdgeScalarQuantity("side", val1);
    //c->addEdgeScalarQuantity("eqid", val5);
    //c->addEdgeScalarQuantity("ehids", val2);
    //c->addEdgeScalarQuantity("count", val3);
    //c->addNodeScalarQuantity("bgn end", val4);
    //c->addEdgeScalarQuantity("x", x)->setEnabled(true);
    c->setEnabled(true);
    c->setColor(glm::vec4(0, 0, 0, 1));
    c->resetTransform();
    c->setRadius(0.0003);
    c->setMaterial("flat");
}
}

#endif
