#ifndef METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#define METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#include "tmesh_mut.h"
#include <queue>

namespace metriko {

namespace {

Row3d pos(const TmeshMut& tm, int nid)      { return get_ptloc_pos(tm.hm, tm.tnodes[nid]); }
Row3d edge_dir(const TmeshMut& tm, int eid) { return tm.hm.edges[eid].half().vec(); }   // tail -> head

vec<Face> node_faces(const TmeshMut& tm, int nid) {
    return std::visit(overloaded{
        [&](const HmLocOnV& v) -> vec<Face> { return get_faces(tm.hm, tm.hm.verts[v.id]); },
        [&](const HmLocOnE& e) -> vec<Face> { return { tm.hm.edges[e.id].face0(), tm.hm.edges[e.id].face1() }; },
        [&](const HmLocOnF& f) -> vec<Face> { return { tm.hm.faces[f.id] }; },
        [&](const auto&)       -> vec<Face> { return {}; },
    }, tm.tnodes[nid]);
}

int seg_face(const TmeshMut& tm, int a, int b) {
    vec<Face> fa = node_faces(tm, a);
    vec<Face> fb = node_faces(tm, b);
    for (Face fx : fa)
    for (Face fy : fb)
        if (fx.id == fy.id) return fx.id;
    return -1;
}

vec<int> thalf_seq(const TmeshMut& tm, const ThalfMut& th) {
    vec<int> seq = tm.tedges[th.teid].nids;
    if (!th.cano) rg::reverse(seq);
    return seq;
}

}

vec<std::tuple<int, double, double>> TmeshMut::allowed_range(int tqid) const {
    const auto& data = tquads[tqid].data;
    struct L {
        const HmLoc* loc;
        Row3d dir;   // この点での境界(thalf)進行方向
    };

    vec<L> locs;
    vec<std::tuple<int, double, double>> res;
    for (auto& [thid, side] : data) {
        const auto& th = thalfs[thid];
        const auto& te = tedges[th.teid];
        int n = te.nids.size();
        for (int i = 0; i < n - 1; ++i) {
            int j  = th.cano ? i : n - 1 - i;            // 現ノードの index
            int jn = th.cano ? j + 1 : j - 1;            // 境界順での次ノード
            int nid = te.nids[j];
            Row3d dir = pos(*this, te.nids[jn]) - pos(*this, nid);
            locs.emplace_back(&tnodes[nid], dir);
        }
    }

    // 各 edge 交点を独立に左側判定し、edge ごとに内側区間 [lo,hi] を交差で絞る。
    // edge を2回横切れば lo/hi 両側が更新されて [min,max] に、1回なら片側だけ更新されて
    // [r,1] か [0,r] になる（ペア/単独を区別する必要がなく used 不要）。
    umap<int, double> lo, hi;
    for (const L& p : locs) {
        const auto* e = std::get_if<HmLocOnE>(p.loc);
        if (!e) continue;
        if (!lo.contains(e->id)) { lo[e->id] = 0.; hi[e->id] = 1.; }
        // edge の head 頂点が thalf 進行方向の左にあれば、その頂点側が内側
        Row3d n = get_ptloc_normal(hm, *p.loc);
        Row3d d = edge_dir(*this, e->id);
        if (n.dot(p.dir.cross(d)) > 0) lo[e->id] = std::max(lo[e->id], e->r);  // 内側 ⊂ [r,1]
        else                           hi[e->id] = std::min(hi[e->id], e->r);  // 内側 ⊂ [0,r]
    }
    for (auto& [eid, l] : lo) if (l < hi[eid]) res.emplace_back(eid, l, hi[eid]);

    // --- 隣接 edge を stack で flood し、内部 edge を range [0,1] で追加 ---
    // 壁 = 交差 edge（候補）と境界頂点（HmLocOnV ノード）。交差 edge の「内側端点」を種に、
    // 頂点づたいに edge を辿り、他の候補 edge か境界頂点に当たるまで内部 edge を集める。
    vec v_stop(hm.nV, false);
    vec v_seen(hm.nV, false);
    vec e_done(hm.nE, false);
    for (const L& p : locs) if (auto* v = std::get_if<HmLocOnV>(p.loc)) v_stop[v->id] = true;
    for (const auto& eid: lo | std::views::keys) e_done[eid] = true;        // 交差 edge は候補（壁）扱い

    vec<int> stack;

    auto push_v = [&](int vid) {
        if (v_stop[vid] || v_seen[vid]) return;            // 境界頂点 / 既訪問は伸ばさない
        v_seen[vid] = true;
        stack.push_back(vid);
    };

    // 種：交差 edge の内側端点（lo==0 → tail 内側 / hi==1 → head 内側）
    for (auto& [eid, l] : lo) {
        if (!(l < hi[eid])) continue;
        Edge e = hm.edges[eid];
        if (l <= 0)       push_v(e.vert0().id);
        if (hi[eid] >= 1) push_v(e.vert1().id);
    }

    while (!stack.empty()) {
        int vid = stack.back(); stack.pop_back();
        for (Half h : hm.verts[vid].adjHalfs()) {
            int eid = h.edge().id;
            if (e_done[eid]) continue;
            e_done[eid] = true;
            res.emplace_back(eid, 0., 1.);
            push_v(h.head().id);
        }
    }

    return res;
}

// Port of allowed_ranges_in_tquad for TmeshMut.
// TmeshMut has no uv (cf), so all geometry is done in 3D: tnodes (HmLoc) -> 3D via
// get_ptloc_pos, and the CCW "left = inside" test uses the face normal:
//   left(d, v)  <=>  dot(normal_f, cross(d, v)) > 0
//
// Returns: eid -> Row2d(r0, r1)  (allowed sub-range along edge.half(), tail->head).
vec<std::tuple<int, double, double>> TmeshMut::allowed_range_(int tqid) const {

    // --- 1. inside sub-range of crossed edges (CCW left = inside, via face normal) ---
    auto data = tquads[tqid].data;
    vec is_crossed(hm.nE, false);
    umap<int, double> lo, hi;

    auto handle_crossing = [&](int nid, const Row3d& d) {
        if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) {
            if (!is_crossed[e->id]) {
                is_crossed[e->id] = true;
                lo[e->id] = 0.;
                hi[e->id] = 1.;
            }
            Row3d n = get_ptloc_normal(hm, tnodes[nid]);
            bool f = n.dot(d.cross(edge_dir(*this, e->id))) > 0; // inside toward head vertex?
            if (f) lo[e->id] = std::max(lo[e->id], e->r); // inside subset of [r, 1]
            else   hi[e->id] = std::min(hi[e->id], e->r); // inside subset of [0, r]
        }
    };

    for (const auto& [thid, _] : data) {
        vec<int> seq = thalf_seq(*this, thalfs[thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(*this, fr, to);
            if (fid == -1) continue;
            Row3d d = pos(*this, to) - pos(*this, fr);
            handle_crossing(fr, d);
            handle_crossing(to, d);
        }
    }

    umap<int, Row2d> allowed;
    for (int eid = 0; eid < hm.nE; ++eid)
        if (is_crossed[eid] && lo[eid] < hi[eid])
            allowed[eid] = Row2d(lo[eid], hi[eid]);

    // --- 2. 領域の2彩色（交差パリティ）。境界曲線は閉じているので、各辺の「交差回数の偶奇」で
    //        面を BFS 彩色すると内側/外側に2分される（向き・サイズに依存しない厳密な方法）。
    //
    //   交差回数 ncross[e]:
    //     ・辺 e 上の境界ノード（HmLocOnE）の個数 = γ が e を横切る回数。
    //     ・頂点通過（HmLocOnV ノード）では γ が辺を横切らずパリティが破れる。そこで
    //       通過頂点で f_in/f_out が共有する辺に仮想交差を1つ足してピンチを閉じる。
    vec<int> ncross(hm.nE, 0);
    {
        vec<char> seen(tnodes.size(), 0);
        for (const auto& [thid, _] : data)
            for (int nid : thalf_seq(*this, thalfs[thid]))
                if (!seen[nid] && (seen[nid] = 1))
                    if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) ncross[e->id]++;
    }
    // 通過頂点の補正：境界の有向セグメントから各頂点ノードの f_in / f_out を求め、共有辺に +1。
    {
        umap<int, int> fin, fout;  // vertex-node id -> incoming / outgoing face
        for (const auto& [thid, _] : data) {
            vec<int> seq = thalf_seq(*this, thalfs[thid]);
            for (size_t i = 0; i + 1 < seq.size(); ++i) {
                int fid = seg_face(*this, seq[i], seq[i + 1]);
                if (fid == -1) continue;
                if (std::holds_alternative<HmLocOnV>(tnodes[seq[i]]))     fout[seq[i]]     = fid;
                if (std::holds_alternative<HmLocOnV>(tnodes[seq[i + 1]])) fin[seq[i + 1]]  = fid;
            }
        }
        for (auto& [nid, fi] : fin) {
            auto it = fout.find(nid);
            if (it == fout.end() || it->second == fi) continue;
            for (Half h : hm.faces[fi].adjHalfs())
                if (h.twin().face().id == it->second) { ncross[h.edge().id]++; break; }
        }
    }

    // BFS 2彩色：奇数交差の辺を渡るたびに色を反転。conflicts>0 なら彩色が破れている
    // （未処理の頂点通過 or 非分離 tquad）＝結果は信頼できない、の警告。
    vec<int> col(hm.nF, -1);
    int conflicts = 0;
    for (int s = 0; s < hm.nF; ++s) {
        if (col[s] != -1) continue;
        col[s] = 0;
        std::queue<int> q; q.push(s);
        while (!q.empty()) {
            int fid = q.front(); q.pop();
            for (Half h : hm.faces[fid].adjHalfs()) {
                int nf = h.twin().face().id;
                if (nf == -1) continue;
                int c = col[fid] ^ (ncross[h.edge().id] & 1);
                if (col[nf] == -1) { col[nf] = c; q.push(nf); }
                else if (col[nf] != c) conflicts++;
            }
        }
    }
    if (conflicts > 0)
        std::cout << "[allowed_range] tq" << tqid << ": parity 2-coloring inconsistent ("
                  << conflicts << " conflicts) — result unreliable\n";

    // 内側 = 小さい方の色クラス。2彩色はクリーン（断片化しない）ので size 判定が安定。
    // CCW 符号は曲面/非分離 tquad で反転し信頼できないため使わない。collapse 対象は十分小さい前提。
    int cnt[2] = {0, 0};
    for (int f = 0; f < hm.nF; ++f) cnt[col[f]]++;
    int inside = cnt[1] < cnt[0] ? 1 : 0;

    for (int fid = 0; fid < hm.nF; ++fid) {
        if (col[fid] != inside) continue;
        for (Half h : hm.faces[fid].adjHalfs())
            if (!is_crossed[h.edge().id]) allowed[h.edge().id] = Row2d(0.0, 1.0);
    }

    // 内部は map セマンティクス（count/上書き）が必要なので umap で組み立て、末尾でフラット化
    vec<std::tuple<int, double, double>> out;
    out.reserve(allowed.size());
    for (auto& [eid, rng] : allowed) out.emplace_back(eid, rng.x(), rng.y());
    return out;
}


// 境界点から横断トレース版。各境界 edge 交差点から、境界に直交する内向き(CCW 左)方向へ
// surface 上を直進(straightest geodesic)させ、反対側境界(交差辺)/頂点に当たるまでに
// 通過した hmesh edge を allowed にする。flood/2彩色を使わないので頂点問題を回避できる。
vec<std::tuple<int, double, double>> TmeshMut::allowed_range_trace(int tqid) const {

    auto data = tquads[tqid].data;
    vec is_crossed(hm.nE, false);
    umap<int, double> lo, hi;
    umap<int, vec<std::pair<Row3d, Row3d>>> faceseg;   // face -> その面内の境界弦（3D端点）
    struct Start { Row3d p; int fid; Row3d dir; int eid; };
    vec<Start> starts;

    auto handle_crossing = [&](int nid, const Row3d& d, int fid) {
        if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) {
            if (!is_crossed[e->id]) { is_crossed[e->id] = true; lo[e->id] = 0.; hi[e->id] = 1.; }
            Row3d n = get_ptloc_normal(hm, tnodes[nid]);
            bool f = n.dot(d.cross(edge_dir(*this, e->id))) > 0;
            if (f) lo[e->id] = std::max(lo[e->id], e->r);
            else   hi[e->id] = std::min(hi[e->id], e->r);
            Row3d dir = n.cross(d);                                  // 内向き(左) = n × d
            if (dir.norm() > 1e-12)
                starts.push_back({ get_ptloc_pos(hm, tnodes[nid]), fid, dir.normalized(), e->id });
        }
    };

    for (const auto& [thid, _] : data) {
        vec<int> seq = thalf_seq(*this, thalfs[thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(*this, fr, to);
            if (fid == -1) continue;
            Row3d pfr = get_ptloc_pos(hm, tnodes[fr]);
            Row3d pto = get_ptloc_pos(hm, tnodes[to]);
            faceseg[fid].push_back({ pfr, pto });           // 境界弦（停止判定に使う）
            Row3d d = pto - pfr;
            handle_crossing(fr, d, fid);
            handle_crossing(to, d, fid);
        }
    }

    umap<int, Row2d> allowed;
    for (int eid = 0; eid < hm.nE; ++eid)
        if (is_crossed[eid] && lo[eid] < hi[eid]) allowed[eid] = Row2d(lo[eid], hi[eid]);

    // straightest geodesic を 1 本トレース
    auto trace = [&](Row3d P, int fid, Row3d dir, int startEdge) {
        for (int step = 0; step < 1000; ++step) {
            Face f = hm.faces[fid];
            Row3d O = f.half().tail().pos(), bx = f.basisX(), by = f.basisY();
            auto to2 = [&](const Row3d& X) -> Row2d { Row3d v = X - O; return Row2d(v.dot(bx), v.dot(by)); };
            Row2d P2 = to2(P);
            Row2d D2(dir.dot(bx), dir.dot(by));
            // ray (P2 + s D2, s>0) と線分 A2-B2 の交差距離 s（無ければ inf）
            auto ray_seg = [&](const Row2d& A2, const Row2d& B2) -> double {
                double ex = B2.x() - A2.x(), ey = B2.y() - A2.y();
                double den = D2.x() * (-ey) - D2.y() * (-ex);
                if (std::abs(den) < 1e-12) return 1e30;
                double rx = A2.x() - P2.x(), ry = A2.y() - P2.y();
                double s = (rx * (-ey) - ry * (-ex)) / den;
                double u = (D2.x() * ry - D2.y() * rx) / den;
                return (s > 1e-7 && u > -1e-9 && u < 1 + 1e-9) ? s : 1e30;
            };
            double best = 1e30; int exitEid = -1; double exitU = 0; Half exitHalf = f.half();
            for (Half h : f.adjHalfs()) {
                Row2d A2 = to2(h.tail().pos()), B2 = to2(h.head().pos());
                double ex = B2.x() - A2.x(), ey = B2.y() - A2.y();
                double den = D2.x() * (-ey) - D2.y() * (-ex);
                if (std::abs(den) < 1e-12) continue;
                double rx = A2.x() - P2.x(), ry = A2.y() - P2.y();
                double s = (rx * (-ey) - ry * (-ex)) / den;
                double u = (D2.x() * ry - D2.y() * rx) / den;
                if (s > 1e-7 && u > -1e-9 && u < 1 + 1e-9 && s < best) { best = s; exitEid = h.edge().id; exitU = u; exitHalf = h; }
            }
            // 反対側境界に到達したら停止：この面内で境界弦（辺に沿わない頂点通過分も含む）を
            // 出口辺より手前で横切ったら、そこが境界。is_crossed だけでは頂点通過を捕捉できない。
            if (auto it = faceseg.find(fid); it != faceseg.end()) {
                Face ff = hm.faces[fid];
                Row3d O2 = ff.half().tail().pos(), bx2 = ff.basisX(), by2 = ff.basisY();
                auto t2 = [&](const Row3d& X) -> Row2d { Row3d v = X - O2; return Row2d(v.dot(bx2), v.dot(by2)); };
                for (auto& [A, B] : it->second)
                    if (ray_seg(t2(A), t2(B)) < best - 1e-9) return;     // 境界弦を横切る → 停止
            }
            if (exitEid == -1) break;
            if (exitU < 1e-4 || exitU > 1 - 1e-4) break;                 // 頂点近傍で停止
            if (is_crossed[exitEid] && exitEid != startEdge) break;      // 交差辺で反対側境界に到達
            if (!is_crossed[exitEid] && !allowed.count(exitEid)) allowed[exitEid] = Row2d(0.0, 1.0);
            Half tw = exitHalf.twin();
            int f2 = tw.face().id;
            if (f2 == -1) break;
            Row3d Pn  = P + dir * best;
            Row3d ev  = exitHalf.vec().normalized();
            double a  = dir.dot(ev);
            double pm = std::sqrt(std::max(0.0, 1.0 - a * a));
            Row3d perp = ev.cross(hm.faces[f2].normal());                // ⟂ ev, f2 平面内
            Row3d mid  = (exitHalf.tail().pos() + exitHalf.head().pos()) * 0.5;
            Row3d c2 = Row3d::Zero(); int cc = 0;
            for (Half h : hm.faces[f2].adjHalfs()) { c2 += h.tail().pos(); ++cc; }
            c2 /= cc;
            if (perp.dot(c2 - mid) < 0) perp = -perp;                     // f2 内向きへ
            dir = (a * ev + pm * perp.normalized()).normalized();
            P = Pn; fid = f2; startEdge = exitEid;
        }
    };

    for (auto& s : starts) trace(s.p, s.fid, s.dir, s.eid);

    vec<std::tuple<int, double, double>> out;
    out.reserve(allowed.size());
    for (auto& [eid, rng] : allowed) out.emplace_back(eid, rng.x(), rng.y());
    return out;
}

}
#endif
