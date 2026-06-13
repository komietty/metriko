#ifndef METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#define METRIKO_TMESH_MUT_COLLAPSE_RANGE_H
#include "tmesh_mut.h"
#include <queue>

namespace metriko {

// Port of allowed_ranges_in_tquad for TmeshMut.
// TmeshMut has no uv (cf), so all geometry is done in 3D: tnodes (HmLoc) -> 3D via
// get_ptloc_pos, and the CCW "left = inside" test uses the face normal:
//   left(d, v)  <=>  dot(normal_f, cross(d, v)) > 0
//
// Returns: eid -> Row2d(r0, r1)  (allowed sub-range along edge.half(), tail->head).
vec<std::tuple<int, double, double>> TmeshMut::allowed_range(int tqid) const {

    auto pos = [&](int nid) -> Row3d { return get_ptloc_pos(hm, tnodes[nid]); };

    // incident faces of a tnode
    auto node_faces = [&](int nid) -> vec<int> {
        return std::visit(overloaded{
            [&](const HmLocOnV& v) -> vec<int> {
                vec<int> fs;
                for (Half h : hm.verts[v.id].adjHalfs()) if (h.face().id != -1) fs.push_back(h.face().id);
                return fs;
            },
            [&](const HmLocOnE& e) -> vec<int> {
                vec<int> fs;
                Edge ed = hm.edges[e.id];
                if (ed.face0().id != -1) fs.push_back(ed.face0().id);
                if (ed.face1().id != -1) fs.push_back(ed.face1().id);
                return fs;
            },
            [&](const HmLocOnF& f) -> vec<int> { return { f.id }; },
            [&](const auto&)       -> vec<int> { return {}; },
        }, tnodes[nid]);
    };

    // the face a segment (between two consecutive tnodes) lies in = their common incident face
    auto seg_face = [&](int a, int b) -> int {
        vec<int> fa = node_faces(a);
        vec<int> fb = node_faces(b);
        for (int x : fa)
        for (int y : fb)
            if (x == y) return x;
        return -1;
    };

    // canonical (tail->head) 3D direction of edge eid
    auto edge_dir = [&](int eid) -> Row3d { return hm.edges[eid].half().vec(); };

    // directed node sequence of a thalf in CCW boundary order (reversed for non-canonical)
    auto thalf_seq = [&](const ThalfMut& th) -> vec<int> {
        vec<int> seq = tedges[th.teid].nids;
        if (!th.cano) rg::reverse(seq);
        return seq;
    };

    // --- 1. inside sub-range of crossed edges (CCW left = inside, via face normal) ---
    auto data = tquads[tqid].data;
    vec is_crossed(hm.nE, false);
    umap<int, double> lo, hi;

    auto handle_crossing = [&](int nid, const Row3d& d, int fid) {
        if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) {
            if (!is_crossed[e->id]) {
                is_crossed[e->id] = true;
                lo[e->id] = 0.;
                hi[e->id] = 1.;
            }
            Row3d n = hm.faces[fid].normal();
            bool f = n.dot(d.cross(edge_dir(e->id))) > 0; // inside toward head vertex?
            if (f) lo[e->id] = std::max(lo[e->id], e->r); // inside subset of [r, 1]
            else   hi[e->id] = std::min(hi[e->id], e->r); // inside subset of [0, r]
        }
    };

    for (const auto& [thid, _] : data) {
        vec<int> seq = thalf_seq(thalfs[thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(fr, to);
            if (fid == -1) continue;
            Row3d d = pos(to) - pos(fr);
            handle_crossing(fr, d, fid);
            handle_crossing(to, d, fid);
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
            for (int nid : thalf_seq(thalfs[thid]))
                if (!seen[nid] && (seen[nid] = 1))
                    if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) ncross[e->id]++;
    }
    // 通過頂点の補正：境界の有向セグメントから各頂点ノードの f_in / f_out を求め、共有辺に +1。
    {
        umap<int, int> fin, fout;  // vertex-node id -> incoming / outgoing face
        for (const auto& [thid, _] : data) {
            vec<int> seq = thalf_seq(thalfs[thid]);
            for (size_t i = 0; i + 1 < seq.size(); ++i) {
                int fid = seg_face(seq[i], seq[i + 1]);
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

    auto node_faces = [&](int nid) -> vec<int> {
        return std::visit(overloaded{
            [&](const HmLocOnV& v) -> vec<int> { vec<int> fs; for (Half h : hm.verts[v.id].adjHalfs()) if (h.face().id != -1) fs.push_back(h.face().id); return fs; },
            [&](const HmLocOnE& e) -> vec<int> { vec<int> fs; Edge ed = hm.edges[e.id]; if (ed.face0().id != -1) fs.push_back(ed.face0().id); if (ed.face1().id != -1) fs.push_back(ed.face1().id); return fs; },
            [&](const HmLocOnF& f) -> vec<int> { return { f.id }; },
            [&](const auto&)       -> vec<int> { return {}; },
        }, tnodes[nid]);
    };
    auto seg_face = [&](int a, int b) -> int {
        vec<int> fa = node_faces(a), fb = node_faces(b);
        for (int x : fa) for (int y : fb) if (x == y) return x;
        return -1;
    };
    auto thalf_seq = [&](const ThalfMut& th) -> vec<int> {
        vec<int> seq = tedges[th.teid].nids; if (!th.cano) rg::reverse(seq); return seq;
    };
    auto edge_dir = [&](int eid) -> Row3d { return hm.edges[eid].half().vec(); };

    auto data = tquads[tqid].data;
    vec is_crossed(hm.nE, false);
    umap<int, double> lo, hi;
    umap<int, vec<std::pair<Row3d, Row3d>>> faceseg;   // face -> その面内の境界弦（3D端点）
    struct Start { Row3d p; int fid; Row3d dir; int eid; };
    vec<Start> starts;

    auto handle_crossing = [&](int nid, const Row3d& d, int fid) {
        if (auto* e = std::get_if<HmLocOnE>(&tnodes[nid])) {
            if (!is_crossed[e->id]) { is_crossed[e->id] = true; lo[e->id] = 0.; hi[e->id] = 1.; }
            Row3d n = hm.faces[fid].normal();
            bool f = n.dot(d.cross(edge_dir(e->id))) > 0;
            if (f) lo[e->id] = std::max(lo[e->id], e->r);
            else   hi[e->id] = std::min(hi[e->id], e->r);
            Row3d dir = n.cross(d);                                  // 内向き(左) = n × d
            if (dir.norm() > 1e-12)
                starts.push_back({ get_ptloc_pos(hm, tnodes[nid]), fid, dir.normalized(), e->id });
        }
    };

    for (const auto& [thid, _] : data) {
        vec<int> seq = thalf_seq(thalfs[thid]);
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            int fr = seq[i], to = seq[i + 1];
            int fid = seg_face(fr, to);
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
