#include <set>
#include "./tmesh_mut.h"

namespace metriko {
bool TmeshMut::collapse_valid_snap_0(Vert v) {
    //for (int nid = 0; nid < tnodes.size(); nid++) {
    //    if (teids[nid] == -1) continue;
    //    auto* l = std::get_if<HmLocOnV>(&tnodes[nid]);
    //    if (l && l->id == v.id) return false;
    //}

    for (auto& [id, nids]: live_tedges()) {
        if (rg::any_of(nids, [&](int nid) {
            auto* l = std::get_if<HmLocOnV>(&tnodes[nid]);
            return l && l->id == v.id;
        })) return false;
    }
    return true;
}

bool TmeshMut::collapse_valid_snap_1(Vert snap_vrt, int snap_nid) {

    for (auto& tq: live_tquads()) {
    for (int side = 0; side < 4; side++) {

        std::set<int>  nids = {};
        umap<int, int> fids = {};

        for (int i: tq.thids(side)) {
        for (int j: tedges[thalfs[i].teid].nids) { nids.insert(j); }}

        for (int nid: nids) {
            if (nid == snap_nid) {
                for (Half h: hm.verts[snap_vrt.id].adjHalfs()) fids[h.face().id]++;
            } else if (auto* l = std::get_if<HmLocOnV>(&tnodes[nid])) {
                for (Half h: hm.verts[l->id].adjHalfs()) fids[h.face().id]++;
            }
        }

        for (auto c: fids | std::views::values) { if (c >= 3) return false; }

    }}

    return true;
}


void TmeshMut::collapse_tedge_snap_joint(int teid) {
    auto  last_nid = tedges[teid].nids.back();
    auto* loc = std::get_if<HmLocOnF>(&tnodes[last_nid]);
    if (!loc) return;
    vec<std::pair<int, int>> te_tails; // teid, the biggest  nid
    vec<std::pair<int, int>> te_heads; // teid, the smallest nid

    for (auto& [id, nids]: live_tedges()) {
        if (last_nid == nids.front()) te_tails.emplace_back(id, 0);
        if (last_nid == nids.back())  te_heads.emplace_back(id, nids.size());
    }

    auto p = get_ptloc_pos(hm, *loc);
    auto d_min = 1e9;
    Vert v_min = {};
    for (Vert v: hm.faces[loc->id].adjHalfs() | vw::transform(&Half::tail)) {
        auto d = (v.pos() - p).squaredNorm();
        if (d < d_min && collapse_valid_snap_0(v)) { d_min = d; v_min = v; }
    }

    for (Face f: v_min.adjHalfs() | vw::transform(&Half::face)) {
        for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(f, tnodes[nids[k]])) j = std::max(j, k); }}
        for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; for (int k = 0; k < nids.size(); k++) { if (is_in_face(f, tnodes[nids[k]])) j = std::min(j, k); }}
    }

    for (auto& [i, j]: te_tails) { auto& nids = tedges[i].nids; if (j >= 2)              nids.erase(nids.begin() + 1, nids.begin() + j);   }
    for (auto& [i, j]: te_heads) { auto& nids = tedges[i].nids; if (j + 2 < nids.size()) nids.erase(nids.begin() + j + 1, nids.end() - 1); }
    tnodes[last_nid] = HmLocOnV{.id = v_min.id};
}

void TmeshMut::collapse_tedge_snap_dedup(int teid) {
    auto& nids = tedges[teid].nids;

    // nodes shared with other chains (junctions / crossings) must stay: dropping
    // one from this chain only would let the straightened chord cross the other
    // chain inside a face without a shared vertex
    vec shared(tnodes.size(), false);
    for (auto& [id, nids_]: live_tedges()) {
        if (id == teid) continue;
        for (int nid: nids_) shared[nid] = true;
    }

    auto faces_of = [&](const HmLoc& l, vec<int>& out) {
        std::visit(overloaded{
            [&](const HmLocOnV& v) { for (Face f: hm.verts[v.id].adjHalfs() | vw::transform(&Half::face)) out.push_back(f.id); },
            [&](const HmLocOnE& e) { out.push_back(hm.edges[e.id].face0().id); out.push_back(hm.edges[e.id].face1().id); },
            [&](const HmLocOnH& h) { Half hh = hm.halfs[h.id]; out.push_back(hh.face().id); out.push_back(hh.twin().face().id); },
            [&](const HmLocOnF& f) { out.push_back(f.id); },
            [&](const auto&)       {},
        }, l);
    };

    for (int i = 1; i < nids.size() - 1;) {
        int nid_prev = nids[i - 1];
        int nid_curr = nids[i];
        int nid_next = nids[i + 1];

        if (shared[nid_curr] || nid_prev == nid_next) { ++i; continue; }

        // a face carrying all three nodes: the chord prev-next stays inside it
        vec<int> fids;
        vec<int> fids_prev;
        vec<int> fids_next;
        faces_of(tnodes[nid_curr], fids);
        faces_of(tnodes[nid_prev], fids_prev);
        faces_of(tnodes[nid_next], fids_next);
        if (nid_curr == 3120) { // 3247
            for (int fid: fids)      std::println("nids curr: {}", fid);
            for (int fid: fids_next) std::println("nids next: {}", fid);
            for (int fid: fids_prev) std::println("nids prev: {}", fid);
        }

        bool same_face = rg::any_of(fids, [&](int fid) { return rg::contains(fids_prev, fid) && rg::contains(fids_next, fid); });
        if (!same_face) { ++i; continue; }

        // drop the middle node; keep the removal only if the validation holds
        int backup = nid_curr;
        nids.erase(nids.begin() + i);
        //if (!collapse_valid_snap_1(Vert{}, nid_curr)) { nids.insert(nids.begin() + i, backup); ++i; }
    }
}

void TmeshMut::collapse_tedge_snap_inter(int teid, int nid, Vert v) {
    auto in_ring = [&](const HmLoc& l) { return rg::any_of(v.adjHalfs(), [&](Half h) { return is_in_face(h.face(), l); });};

    auto& nids = tedges[teid].nids;

    auto it = rg::find(nids, nid);
    if (it == nids.end()) return;   // trimmed away by an earlier snap
    int i_cur = rg::find(nids, nid) - nids.begin();
    int i_min = i_cur;
    int i_max = i_cur;
    while (i_min > 0               && in_ring(tnodes[nids[i_min - 1]])) --i_min;
    while (i_max < nids.size() - 1 && in_ring(tnodes[nids[i_max + 1]])) ++i_max;

    //for (int i = i_min + 1; i < i_max; i++) teids[nids[i]] = -1;
    if (i_max - i_cur >= 2) nids.erase(nids.begin() + i_cur + 1, nids.begin() + i_max);
    if (i_cur - i_min >= 2) nids.erase(nids.begin() + i_min + 1, nids.begin() + i_cur);
    tnodes[nid] = HmLocOnV{.id = v.id};
}


void TmeshMut::collapse_tedge_snap(bool flag) {
    // 1: snap joint tnodes
    for (auto& [teid, nids]: live_tedges()) { collapse_tedge_snap_joint(teid); }

    struct Cand { int nid; int eid; Vert v; double d; };
    vec candidates(tnodes.size(), vec<Cand>{});

    // 2.0: snap inter tnodes
    for (auto& [teid, nids]: live_tedges()) {
    for (int i = 1; i < nids.size() - 1; i++) {
        auto nid = nids[i];
        auto pos = get_ptloc_pos(hm, tnodes[nid]);

        std::visit(overloaded{
            [&](const HmLocOnH& l) {
                Vert v0 = hm.halfs[l.id].tail();
                Vert v1 = hm.halfs[l.id].head();
                candidates[nid].emplace_back(nid, teid, v0, (v0.pos() - pos).squaredNorm());
                candidates[nid].emplace_back(nid, teid, v1, (v1.pos() - pos).squaredNorm());
            },
            [&](const HmLocOnE& l) {
                Vert v0 = hm.edges[l.id].vert0();
                Vert v1 = hm.edges[l.id].vert1();
                candidates[nid].emplace_back(nid, teid, v0, (v0.pos() - pos).squaredNorm());
                candidates[nid].emplace_back(nid, teid, v1, (v1.pos() - pos).squaredNorm());
            },
            [&](const HmLocOnF& l) {
                Vert v0 = hm.faces[l.id].half().tail();
                Vert v1 = hm.faces[l.id].half().head();
                Vert v2 = hm.faces[l.id].half().crnr().vert();
                candidates[nid].emplace_back(nid, teid, v0, (v0.pos() - pos).squaredNorm());
                candidates[nid].emplace_back(nid, teid, v1, (v1.pos() - pos).squaredNorm());
                candidates[nid].emplace_back(nid, teid, v2, (v2.pos() - pos).squaredNorm());
            },
            [&](const auto&) {},
        }, tnodes[nid]);
    }}

    // 2.1: sort candidates in inner/outer order
    for (auto& c: candidates) rg::sort(c, {}, &Cand::d);
    rg::sort(candidates, {}, [](const vec<Cand>& c) { return c.empty() ? 1e9 : c.front().d; });

    // 2.2: snap if it's valid
    if (flag) {
        for (auto& c: candidates) {
        for (auto& [nid, eid, vrt, _]: c) {
            //if (collapse_valid_snap_0(vrt) && collapse_valid_snap_1(vrt, nid) ) {
            if (collapse_valid_snap_0(vrt)) {
                collapse_tedge_snap_inter(eid, nid, vrt);
                break;
            }
        }}
    } else {
        for (auto& c: candidates) {
            if (c.size() == 0) continue;
            auto& [nid, eid, vrt, _] = c.front();
            //if (collapse_valid_snap_0(vrt) && collapse_valid_snap_1(vrt, nid) ) { collapse_tedge_snap_inter(eid, nid, vrt); }
            if (collapse_valid_snap_0(vrt)) { collapse_tedge_snap_inter(eid, nid, vrt); }
        }
    }
}
}
