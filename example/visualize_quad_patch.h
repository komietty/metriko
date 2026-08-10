#ifndef TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_QUAD_PATCH_H
#define TMESH_MUT_COLLAPSE_TQUAD_CPP_VISUALIZE_QUAD_PATCH_H

namespace metriko::visualizer {
            // ------ per-quad tqid: replay the t-mesh on the quad graph ------
            // tedges anchored at singularities are walked first (their quad vertex
            // is pinned to the grid exactly); each finished walk fixes its far
            // junction combinatorially, which anchors the tedges meeting there.
            // no re-race: each walk is exactly x quad edges long.
            /*
            {
                // weld the quad soup (combinatorics only, positions unsmoothed)
                MatXd qv;
                VecXi SVI, SVJ;
                const double weps = 1e-7 * (pos.colwise().maxCoeff() - pos.colwise().minCoeff()).norm();
                igl::remove_duplicate_vertices(pos, weps, qv, SVI, SVJ);
                MatXi qidx(l, 4);
                for (int i = 0; i < l; ++i) for (int j = 0; j < 4; ++j) qidx(i, j) = SVJ(idx(i, j));

                // directed edge -> (quad, corner); vertex -> outgoing edge heads
                std::map<std::pair<int, int>, std::pair<int, int>> de;
                std::unordered_map<int, vec<int>> out;
                for (int i = 0; i < l; ++i)
                for (int j = 0; j < 4; ++j) {
                    int a = qidx(i, j), b = qidx(i, (j + 1) % 4);
                    de[{a, b}] = {i, j};
                    out[a].push_back(b);
                }
                // next outgoing edge head, rotating CCW around v from (v -> w)
                auto rot = [&](int v, int w) {
                    auto [q, j] = de.at({v, w});
                    return qidx(q, (j + 3) % 4);   // prev corner of the quad
                };
                // arrived via (t -> v): the straight continuation is two rotations
                auto straight = [&](int t, int v) { return rot(v, rot(v, t)); };

                auto find_qv = [&](const Row3d& p) {
                    for (int i = 0; i < (int)qv.rows(); ++i)
                        if ((Row3d(qv.row(i)) - p).norm() < weps * 10) return i;
                    return -1;
                };

                umap<int, int> xs;                       // teid -> steps
                umap<int, int> te2th;                    // teid -> cano thalf id
                umap<int, std::pair<int, int>> sides;    // teid -> (left, right) tqid
                for (auto& th: tmm.thalfs) {
                    if (th.id == -1 || !th.cano) continue;
                    xs[th.teid]    = (int)std::round(th.x);
                    te2th[th.teid] = th.id;
                    sides[th.teid] = {th.tqid, tmm.thalfs[th.twid].tqid};
                }

                // node -> incident tedges (with which endpoint touches the node)
                umap<int, vec<std::pair<int, bool>>> node_te;
                for (auto& [teid, nids]: tmm.live_tedges()) {
                    node_te[nids.front()].emplace_back(teid, true);
                    node_te[nids.back()].emplace_back(teid, false);
                }
                // outgoing 3d direction of a tedge at one of its endpoint nodes
                auto dir_at = [&](int teid, int nid) -> Row3d {
                    auto& nids = tmm.tedges[teid].nids;
                    int nb = nids.front() == nid ? nids[1] : nids[nids.size() - 2];
                    return (get_ptloc_pos(hm, tmm.tnodes[nb]) - get_ptloc_pos(hm, tmm.tnodes[nid])).normalized();
                };

                std::vector<double> quad_tqid(l, -1);
                std::set<std::pair<int, int>> bnd;       // traced (undirected) quad edges
                vec done(tmm.tedges.size(), false);

                struct Job { int teid; bool fwd; int s; int b; };   // walk from quad vertex s, first head b
                std::queue<Job> jobs;

                // walk x quad edges; the cano thalf's tquad lies on the LEFT of the
                // nids direction (same convention as the cut annotation). returns the
                // far vertex and the tail of the arriving edge, or {-1, -1}
                auto walk = [&](const Job& jb) -> std::pair<int, int> {
                    auto [lq, rq] = sides.at(jb.teid);
                    if (!jb.fwd) std::swap(lq, rq);
                    int steps = xs.at(jb.teid);
                    int a = jb.s, b = jb.b;
                    for (int k = 0; k < steps; ++k) {
                        quad_tqid[de.at({a, b}).first] = lq;
                        quad_tqid[de.at({b, a}).first] = rq;
                        bnd.insert(std::minmax(a, b));
                        if (k + 1 < steps) {
                            if ((int)out[b].size() != 4) return {-1, -1};   // early singular: mismatch
                            int c = straight(a, b);
                            a = b;
                            b = c;
                        }
                    }
                    return {b, a};
                };

                // seeds: endpoints at singularities, whose quad vertex is exact
                for (auto& [teid, nids]: tmm.live_tedges()) {
                    if (xs.at(teid) <= 0 || nids.size() < 2) { done[teid] = true; continue; }
                    for (bool fwd: {true, false}) {
                        int nid = fwd ? nids.front() : nids.back();
                        auto* lv = std::get_if<HmLocOnV>(&tmm.tnodes[nid]);
                        if (!lv || singular(lv->id) == 0) continue;
                        int s = find_qv(hm.verts[lv->id].pos());
                        if (s < 0) continue;
                        Row3d d3 = dir_at(teid, nid);
                        int best = -1;
                        double bd = -2;
                        for (int w: out[s]) {
                            double d = (Row3d(qv.row(w)) - Row3d(qv.row(s))).normalized().dot(d3);
                            if (d > bd) { bd = d; best = w; }
                        }
                        jobs.push({teid, fwd, s, best});
                        break;
                    }
                }

                int failed = 0;
                while (!jobs.empty()) {
                    Job jb = jobs.front();
                    jobs.pop();
                    if (done[jb.teid]) continue;
                    done[jb.teid] = true;
                    auto [vj, tail] = walk(jb);
                    if (vj < 0) { ++failed; continue; }

                    // the walk fixed the far junction: anchor the tedges meeting
                    // there. slots are assigned COMBINATORIALLY: going CCW around the
                    // node, consecutive tedges share a patch (left of the current =
                    // right of the next), and that patch spans one slot when the node
                    // is its corner, two when its boundary runs straight through
                    auto& nids = tmm.tedges[jb.teid].nids;
                    int nid = jb.fwd ? nids.back() : nids.front();
                    if ((int)out[vj].size() != 4) continue;   // singular end: those seed themselves

                    int ring[4];   // quad-edge slots, CCW from the reverse of arrival
                    ring[0] = tail;
                    for (int k = 1; k < 4; ++k) ring[k] = rot(vj, ring[k - 1]);

                    auto lr_out = [&](int te, bool afn) {     // left/right of the OUTGOING direction
                        auto [lq, rq] = sides.at(te);
                        return afn ? std::pair(lq, rq) : std::pair(rq, lq);
                    };
                    auto side_in = [&](int te, int P) {       // side index of te's thalf bounding P
                        int c = te2th.at(te);
                        int h = tmm.thalfs[c].tqid == P ? c : tmm.thalfs[c].twid;
                        return tmm.tquads[P].side_of(tmm.thalfs[h]);
                    };

                    auto ports = node_te[nid];
                    int  cur  = jb.teid;
                    bool afn  = !jb.fwd;   // the arrival tedge, oriented outgoing at this node
                    int  slot = 0;
                    for (size_t it = 1; it < ports.size(); ++it) {
                        int P = lr_out(cur, afn).first;       // the wedge CCW after cur
                        auto nx = rg::find_if(ports, [&](auto& pr) {
                            return pr.first != cur && lr_out(pr.first, pr.second).second == P;
                        });
                        if (nx == ports.end()) break;         // inconsistent incidence data
                        slot += side_in(cur, P) == side_in(nx->first, P) ? 2 : 1;
                        if (!done[nx->first]) jobs.push({nx->first, nx->second, vj, ring[slot % 4]});
                        cur = nx->first;
                        afn = nx->second;
                    }
                }
                int unreached = 0;
                for (auto& [teid, nids]: tmm.live_tedges()) if (!done[teid]) ++unreached;
                if (failed || unreached) std::println("[tqid] replay: {} failed, {} unreached", failed, unreached);

                // flood fill the patch interiors, never crossing a traced edge
                std::queue<int> que;
                for (int i = 0; i < l; ++i) if (quad_tqid[i] >= 0) que.push(i);
                while (!que.empty()) {
                    int q = que.front(); que.pop();
                    for (int j = 0; j < 4; ++j) {
                        int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
                        if (bnd.contains(std::minmax(a, b))) continue;
                        int nb = de.at({b, a}).first;
                        if (quad_tqid[nb] < 0) { quad_tqid[nb] = quad_tqid[q]; que.push(nb); }
                    }
                }
                if (int u = (int)rg::count(quad_tqid, -1.); u > 0) std::println("[tqid] {} quads unlabeled", u);
                quad->addFaceScalarQuantity("tqid", quad_tqid)->setEnabled(true);
            }
            */
}

#endif
