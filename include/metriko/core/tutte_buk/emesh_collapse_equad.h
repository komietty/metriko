#ifndef TMESH_H_EMESH_COLLAPSE_EQUAD_H
#define TMESH_H_EMESH_COLLAPSE_EQUAD_H
#include "./emesh.h"

namespace metriko::tutte {
inline bool Emesh::collapse_quad(const int eqid) {
    int side = -1;
    const auto& eq = equads[eqid];
    for (int i : {0, 1}) {
        auto side_ehids = eq.ehids_by_side(i);
        int sum = (int)rg::fold_left(side_ehids | vw::transform([this](int id){ return ehalfs[id].x; }), 0., std::plus());
        if (sum == 0) { side = i; break; }
    }

    if (side == -1) return true;
    std::cout << "-------- collapse eqid: " << eqid << std::endl;

    int side_a   = side == 0 ? 1 : 2;                   // remain side a
    int side_b   = side == 0 ? 3 : 0;                   // remain side b
    auto ehids_p = eq.ehids_by_side(side == 0 ? 0 : 1); // collapse side p
    auto ehids_q = eq.ehids_by_side(side == 0 ? 2 : 3); // collapse side q
    auto ehids_a = eq.ehids_by_side(side_a);            // remain side a
    auto ehids_b = eq.ehids_by_side(side_b);            // remain side b

    auto find_terminal = [&](const vec<int>& ehids, int s0, int s1) -> std::tuple<int, int, int> {
        for (int ehid: ehids) {
            const Ehalf& eh0 = ehalfs[ehid];
            const Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.bgn) return {eh0.tail().id, eh0.head().id, s0};
            if (eh1.bgn) return {eh1.tail().id, eh1.head().id, s1};
        }
        for (int ehid: ehids) {
            const Ehalf& eh0 = ehalfs[ehid];
            const Ehalf& eh1 = ehalfs[eh0.twid];
            if (eh0.end) return {eh0.head().id, eh0.tail().id, s1};
            if (eh1.end) return {eh1.head().id, eh1.tail().id, s0};
        }
        return {-1, -1, -1};
    };

    auto [bgn, markedB, sB] = find_terminal(ehids_p, side_b, side_a); // vert, other marked vert and side of the bgn
    auto [end, markedE, sE] = find_terminal(ehids_q, side_a, side_b); // vert, other marked vert and side of the end

    if (bgn == -1 || end == -1) {
        equads[eqid].debug_draw();
        std::cerr << "collapse eqid failed: " << eqid << std::endl;
        return false;
    }

    double sum = 0;
    double cmp = 0;
    for (int ehid: ehids_a) {
        sum += ehalfs[ehid].x;
        cmp += ehalfs[ehid].halfs.size();
    }

    vec<AuxDijkData> aux;
    aux.emplace_back(AuxDijkData{hm.verts[bgn], sB, 0., -1e6});
    aux.emplace_back(AuxDijkData{hm.verts[end], sE, sum, 1e6});

    double sum1 = 0;
    double cmp1 = 0;
    for (int ehid: ehids_a) {
        Ehalf& eh = ehalfs[ehid];
        sum1 += eh.x;
        cmp1 += eh.halfs.size();
        if (eh.head().id == bgn) { continue; }
        if (eh.head().id == end) { continue; }
        if (eh.head().id == markedB) { continue; }
        if (eh.head().id == markedE) { continue; }
        aux.emplace_back(AuxDijkData{ eh.head(), side_a, sum1, cmp1});
    }
    for (int ehid: ehids_b) {
        Ehalf& eh = ehalfs[ehid];
        sum1 -= eh.x;
        cmp1 -= eh.halfs.size();
        if (eh.head().id == bgn) { continue; }
        if (eh.head().id == end) { continue; }
        if (eh.head().id == markedB) { continue; }
        if (eh.head().id == markedE) { continue; }
        aux.emplace_back(AuxDijkData{ eh.head(), side_b, sum1, cmp1});
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

    /*
    std::vector<glm::vec3> verts_to_passby;
    for (const AuxDijkData &a: aux_sorted) {
        Row3d p = a.vert.pos();
        verts_to_passby.emplace_back(p.x(), p.y(), p.z());
    }
    auto vq = polyscope::registerPointCloud("verts to pass by eqid " + std::to_string(qid), verts_to_passby);
    vq->setEnabled(false);
    vq->setPointRadius(0.002);
    vq->resetTransform();
    */

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
        const auto d = std::abs(a0.val - a1.val);

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
        bool bgn0 = this->sings.contains(path0.front().tail().id);
        bool end0 = this->sings.contains(path0.back().head().id);
        bool bgn1 = this->sings.contains(path1.front().tail().id);
        bool end1 = this->sings.contains(path1.back().head().id);
        ehalfs.emplace_back(this, path0, n, n + 1, -1, d, bgn0, end0);
        ehalfs.emplace_back(this, path1, n + 1, n, -1, d, bgn1, end1);
    }

    auto path_md = path_mb | vw::transform([&](int ehid) { return ehalfs[ehid].twid; }) | rg::to<vec<int>>();
    rg::reverse(path_md);

    /*
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
        c->setEnabled(false);
        c->addEdgeScalarQuantity("order", val1);
        c->resetTransform();
        c->setRadius(0.001);
    }
    */

    replace_remain_side(eq, path_mb, ehids_p, ehids_q, side_a);
    replace_remain_side(eq, path_md, ehids_q, ehids_p, side_b);
    equads[eqid].id = -1;
    return true;
}

inline void Emesh::replace_remain_side(
    const Equad &eq,
    const vec<int> &ehids_path,
    const vec<int> &ehids_prev,
    const vec<int> &ehids_next,
    int side
) {
    for (int ehid: eq.ehids_by_side(side)) {
        const Ehalf& currEH = ehalfs[ehid];
        const Ehalf& twinEH = ehalfs[currEH.twid];

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

        equads[twinEH.eqid].replace_ehalf(
            twinEH.id, rp,
            f1 && !f2 ? ehids_prev : vec<int>{},
            f3 && !f4 ? ehids_next : vec<int>{}
        );
    }
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
        //auto item = reps | vw::transform([side](int r) { return Edata{r, side}; });
        //edata.insert_range(addr, item);

        vec<Edata> items;
        items.reserve(reps.size());
        for (int r : reps) {
            items.push_back(Edata{r, side});
            const_cast<Ehalf&>(em->ehalfs[r]).eqid = this->id;
        }
        edata.insert_range(addr, items);
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
        auto& eh = const_cast<Ehalf&>(em->ehalfs[it->ehid]);
        auto& tw = const_cast<Ehalf&>(em->ehalfs[eh.twid]);

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
        auto& eh = const_cast<Ehalf&>(em->ehalfs[it->ehid]);
        auto& tw = const_cast<Ehalf&>(em->ehalfs[eh.twid]);
        for (int ehid_: std::views::reverse(ext0)) {
            const Ehalf& ext_eh = em->ehalfs[ehid_];
            const Ehalf& ext_tw = em->ehalfs[ext_eh.twid];
            eh.extend_prev(ext_eh);
            tw.extend_next(ext_tw);
            remove_from_twin_equad(ext_tw.eqid, ext_tw.id);
        }
    }
}

}
#endif
