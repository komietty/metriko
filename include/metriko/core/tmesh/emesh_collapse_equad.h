#ifndef EXAMPLE_EBD_CPP_EMESH_COLLAPSE_EQUAD_H
#define EXAMPLE_EBD_CPP_EMESH_COLLAPSE_EQUAD_H
#include "metriko/core/tmesh/emesh_collapse_util.h"
#include "metriko/core/tmesh/emesh_collapse_dijkstra.h"

namespace metriko {

inline bool Emesh::collapse_equad(const int eqid) {
    auto eq = equads[eqid];
    auto side = -1;
    for (int i : {0, 1}) {
        int sum = rg::fold_left(eq.ehids(i) | vw::transform([&](int id){ return ehalfs[id].x; }) , 0., std::plus());
        if (sum == 0) { side = i; break; }
    }

    if (side == -1) return true;
    auto side_a  = side == 0 ? 1 : 2;           // remain side a
    auto side_b  = side == 0 ? 3 : 0;           // remain side b
    auto ehids_p = eq.ehids(side == 0 ? 0 : 1); // collapse side p
    auto ehids_q = eq.ehids(side == 0 ? 2 : 3); // collapse side q
    auto ehids_a = eq.ehids(side_a);            // remain side a
    auto ehids_b = eq.ehids(side_b);            // remain side b

    auto find_terminal = [&](const vec<int>& ehids, int s0, int s1) -> std::tuple<int, int, int> {
        for (int ehid: ehids) {
            auto& eh0 = ehalfs[ehid];
            auto& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) return {eh0.tail().id, eh0.head().id, s0};
            if (eh1.bgn) return {eh1.tail().id, eh1.head().id, s1};
        }
        for (int ehid: ehids) {
            auto& eh0 = ehalfs[ehid];
            auto& eh1 = ehalfs[eh0.twid];
            if (eh0.end) return {eh0.head().id, eh0.tail().id, s1};
            if (eh1.end) return {eh1.head().id, eh1.tail().id, s0};
        }
        return {-1, -1, -1};
    };

    auto [bgn, markedB, sB] = find_terminal(ehids_p, side_b, side_a); // vert, another marked vert and side of the bgn
    auto [end, markedE, sE] = find_terminal(ehids_q, side_a, side_b); // vert, another marked vert and side of the end

    if (bgn == -1 || end == -1) { std::cerr << "collapse eqid failed" << std::endl; return false; }

    double sum = 0;
    double cmp = 0;
    for (int ehid: ehids_a) {
        sum += ehalfs[ehid].x;
        cmp += ehalfs[ehid].halfs.size();
    }

    vec<AuxDijkData> aux;
    aux.push_back({hm.verts[bgn], sB, 0., -1e6});
    aux.push_back({hm.verts[end], sE, sum, 1e6});

    double sum1 = 0;
    double cmp1 = 0;
    for (int ehid: ehids_a) {
        Ehalf& eh = ehalfs[ehid];
        sum1 += eh.x;
        cmp1 += eh.halfs.size();
        if (eh.head().id == bgn)     { continue; }
        if (eh.head().id == end)     { continue; }
        if (eh.head().id == markedB) { continue; }
        if (eh.head().id == markedE) { continue; }
        aux.push_back({eh.head(), side_a, sum1, cmp1});
    }
    for (int ehid: ehids_b) {
        Ehalf& eh = ehalfs[ehid];
        sum1 -= eh.x;
        cmp1 -= eh.halfs.size();
        if (eh.head().id == bgn)     { continue; }
        if (eh.head().id == end)     { continue; }
        if (eh.head().id == markedB) { continue; }
        if (eh.head().id == markedE) { continue; }
        aux.push_back({eh.head(), side_b, sum1, cmp1});
    }


    if (aux.empty()) return true;
    assert(sum1 == 0);

    rg::sort(aux, [](const AuxDijkData& a, const AuxDijkData& b) { return std::tie(a.val, a.cmp) < std::tie(b.val, b.cmp); });

    // 2. val, vert.id, side が同じものを重複とみなして消す
    auto ret = rg::unique(aux, [](const AuxDijkData& a, const AuxDijkData& b) {
        return std::tie(a.val, a.vert.id, a.side) == std::tie(b.val, b.vert.id, b.side);
    });

    aux.erase(ret.begin(), ret.end());

    vec<int> path_mb;

    auto visit = vec(hm.nH, false);
    vec aux_sorted(aux.begin(), aux.end());
    vec aux_visited(aux.size(), false);
    int N = (int)aux.size() - 1;
    path_mb.resize(N);

    {
        std::vector<glm::vec3> verts_to_passby;
        for (const AuxDijkData &a: aux_sorted) {
            Row3d p = a.vert.pos();
            verts_to_passby.emplace_back(p.x(), p.y(), p.z());
        }
        auto vq = polyscope::registerPointCloud("verts to pass by eqid-" + std::to_string(eqid), verts_to_passby);
        vq->setEnabled(false);
        vq->setPointRadius(0.002);
        vq->resetTransform();
    }

    for (int i = 0; i < N; i++) {
        auto& a0 = aux_sorted[i];
        auto& a1 = aux_sorted[i + 1];
        Vert v0 = a0.vert;
        Vert v1 = a1.vert;

        for (int si: {side_a, side_b}) {
            if (a0.side != si || a1.side != si) continue;
            for (int ehid: eq.ehids(si)) {
                auto& eh0 = ehalfs[ehid];
                auto& eh1 = ehalfs[eh0.twid];
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
        const auto d = std::abs(a0.val - a1.val);

        vec<bool> allow_verts = vec(hm.nV, false);
        auto inside = verts_inside(eq, *this, hm);
        for (int vid: inside) { allow_verts[vid] = true; }

        // points in the collapse side a and b are not allowed to pass by
        for (auto ehids: {ehids_p, ehids_q}) {
            for (int ehid: ehids) {
                for (Half h: ehalfs[ehid].halfs) {
                    allow_verts[h.tail().id] = false;
                    allow_verts[h.head().id] = false;
                }
            }
        }
        allow_verts[v0.id] = true;
        allow_verts[v1.id] = true;

        auto path = compute_dijkstra(hm, visit, allow_verts, v0, v1);
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

        bool bgn0 = sings.contains(path0.front().tail().id);
        bool end0 = sings.contains(path0.back().head().id);
        bool bgn1 = sings.contains(path1.front().tail().id);
        bool end1 = sings.contains(path1.back().head().id);
        ehalfs.push_back({.em = this, .id = n, .twid = n + 1, .eqid = -1, .bgn = bgn0, .end = end0, .x = d, .halfs = path0});
        ehalfs.push_back({.em = this, .id = n + 1, .twid = n, .eqid = -1, .bgn = bgn1, .end = end1, .x = d, .halfs = path1});
    }

    auto path_md = path_mb | vw::transform([&](int ehid) { return ehalfs[ehid].twid; }) | rg::to<vec<int>>();
    rg::reverse(path_md);

    replace_remain_side(eq.id, side_a, path_mb, ehids_p, ehids_q);
    replace_remain_side(eq.id, side_b, path_md, ehids_q, ehids_p);
    equads[eqid].id = -1;

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

        auto c = polyscope::registerCurveNetwork("path_md of eqid" + std::to_string(eqid), ns, es);
        c->setEnabled(false);
        c->addEdgeScalarQuantity("order", val1);
        c->resetTransform();
        c->setRadius(0.001);
    }

    return true;
}

inline void Emesh::replace_remain_side(
    int eqid,
    int side,
    const vec<int> &ehids_path,
    const vec<int> &ehids_prev,
    const vec<int> &ehids_next
) {
    const auto& eq = equads[eqid];

    for (int ehid: eq.ehids(side)) {
        auto& currEH = ehalfs[ehid];
        auto& twinEH = ehalfs[currEH.twid];

        if (rg::find(ehids_path, currEH.id) != ehids_path.end()) continue;

        auto rp = vec<int>{};
        bool f1 = currEH.tail() == ehalfs[ehids_prev.back() ].head();
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

        replace_ehalf(
            twinEH.eqid,
            twinEH.id,
            rp,
            f1 && !f2 ? ehids_prev : vec<int>{},
            f3 && !f4 ? ehids_next : vec<int>{}
        );
    }
}

inline void Emesh::replace_ehalf(
    int eqid,
    int ehid,             // tgt ehalf id
    const vec<int>& reps, // ehids for replacing ehalf above
    const vec<int>& ext0, // ehids for extending next side of ehalfs
    const vec<int>& ext1  // ehids for extending prev side of ehalfs
) {
    auto& data = equads[eqid].data;
    int side;

    { // replace side
        auto it = rg::find(data, ehid, &Edata::ehid);
        assert(it != data.end());
        side = it->side;
        auto addr = data.erase(it);
        //auto item = reps | vw::transform([side](int r) { return Edata{r, side}; });
        //edata.insert_range(addr, item);

        vec<Edata> items;
        items.reserve(reps.size());
        for (int r : reps) {
            items.push_back(Edata{r, side});
            ehalfs[r].eqid = equads[eqid].id;
        }
        data.insert_range(addr, items);
    }

    auto remove_from_twin_equad = [&](int twin_eqid, int twin_ehid) {
        auto& eq = const_cast<Equad&>(equads[twin_eqid]);
        auto it = rg::find(eq.data, twin_ehid, &Edata::ehid);
        if (it != eq.data.end()) eq.data.erase(it);
    };

    { // insert for the prev side
        auto rev_edata = data | vw::reverse;
        auto it = rg::find(rev_edata, (side + 3) % 4, &Edata::side);
        assert(it != rev_edata.end());
        auto& eh = const_cast<Ehalf&>(ehalfs[it->ehid]);
        auto& tw = const_cast<Ehalf&>(ehalfs[eh.twid]);

        for (int ehid_: ext1) {
            auto& ext_eh = ehalfs[ehid_];
            auto& ext_tw = ehalfs[ext_eh.twid];
            eh.extend_next(ext_eh);
            tw.extend_prev(ext_tw);
            remove_from_twin_equad(ext_tw.eqid, ext_tw.id);
        }
    }

    { // insert for the next side
        auto it = rg::find(data, (side + 1) % 4, &Edata::side);
        assert(it != data.end());
        auto& eh = const_cast<Ehalf&>(ehalfs[it->ehid]);
        auto& tw = const_cast<Ehalf&>(ehalfs[eh.twid]);
        for (int ehid_: std::views::reverse(ext0)) {
            auto& ext_eh = ehalfs[ehid_];
            auto& ext_tw = ehalfs[ext_eh.twid];
            eh.extend_prev(ext_eh);
            tw.extend_next(ext_tw);
            remove_from_twin_equad(ext_tw.eqid, ext_tw.id);
        }
    }
}
}

#endif
