//
// Created by saki on 2026/07/14.
//
// demo: collapse the TmeshMut, then cut the original hmesh along the collapsed
// t-mesh (emesh_cutting) and display the resulting cut mesh with its patch
// boundaries.
#include <fstream>
#include <igl/readOBJ.h>
#include <igl/slim.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "metriko/core/tutte/tutte_cutting.h"
#include "metriko/core/tutte/tutte_params.h"
#include "metriko/core/qex/sanitization.h"
#include "metriko/core/qex/gen_q_vert.h"
#include "metriko/core/qex/gen_q_port.h"
#include "metriko/core/qex/gen_q_edge.h"
#include "metriko/core/qex/gen_q_face.h"
#include "metriko/core/qex/refinement.h"
#include "common.h"
#include "visualize_hmesh.h"
#include "visualize_qex.h"
#include "visualize_tmesh_mut.h"

using namespace metriko;

static VecXc uv2;
static MatXd V;
static MatXi F;
static VecXi matching;
static VecXi singular;
static vec<bool> seam;

int main(int argc, char** argv) {
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);
    if (!load_cache(std::format("{}.{}.cache", argv[1], argv[2]), uv2, matching, singular, seam)) throw std::runtime_error("the cache does not exist");

    ///--- gen mport, medge ---///
    auto mg = mc::Mgrph(hm, uv2, matching, singular);
    auto tm = Tmesh(mg);
    auto X  = compute_quantization(tm, mg);
    TmeshMut tmm(mg, tm, X);

    for (int i = 0; i < 10; ++i) {
        for (const ThalfMut& th0: tmm.thalfs) {
            if (th0.id == -1) continue;
            auto& th1 = tmm.thalfs[th0.twid];
            if (th1.id == -1) continue;
            if (th0.x != 0)   continue;
            auto& tq0 = tmm.tquads[th0.tqid];
            auto& tq1 = tmm.tquads[th1.tqid];
            if (tq0.thids(tq0.side_of(th0)).size() == 1) continue;
            if (tq1.thids(tq1.side_of(th1)).size() == 1) continue;
            std::cout << "th collapse: " << th0.id << std::endl;
            tmm.collapse_thalf(th0.id);
        }
        for (const TquadMut& tq: tmm.live_tquads()) {
            Tqchain chain;
            if (tmm.collapse_tquad_chain_prepare(tq.id, chain)) {
                std::cout << "tq collapse: " << tq.id << std::endl;
                tmm.collapse_tquad_chain_execute(chain);
            }
        }
    }

    auto t0 = std::chrono::steady_clock::now();
    tmm.collapse_tedge_snap(false);
    tmm.collapse_tedge_snap(true);
    for (const auto& [teid, _] : tmm.live_tedges()) { tmm.collapse_tedge_snap_dedup(teid); }
    auto t1 = std::chrono::steady_clock::now();
    std::println("[time] snap: {:.3f}s", std::chrono::duration<double>(t1 - t0).count());


    ///--- validate tedge nid chains: duplicated / backtracking nodes break the cut ---///
    for (const auto& th: tmm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tmm.tedges[th.teid].nids;
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
            Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
            if (nids[k] == nids[k + 1] || (a - b).norm() < 1e-12) std::println("[warn] tedge {} (thid {}, x {}): duplicate node at {} (nid {} / {}, locs {} / {})", th.teid, th.id, th.x, k, nids[k], nids[k + 1], loc_str(tmm.tnodes[nids[k]]), loc_str(tmm.tnodes[nids[k + 1]]));
        }
        for (size_t k = 0; k + 2 < nids.size(); ++k) {
            if (nids[k] == nids[k + 2]) std::println("[warn] tedge {} (thid {}, x {}): backtrack at {} (nid {})", th.teid, th.id, th.x, k, nids[k]);
        }
    }

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(hm, tmm, seam, matching, singular, seam1, matching1, singular1, hdata);
    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);


    // --- validate: every halfedge must have its opposite pair. ---
    {
        std::map<std::pair<int, int>, int> cnt;
        for (int i = 0; i < hm_emb->idx.rows(); ++i)
        for (int j = 0; j < 3; ++j) {
            int a = hm_emb->idx(i, j);
            int b = hm_emb->idx(i, (j + 1) % 3);
            cnt[{a, b}]++;
        }
        int unpaired = 0, duplicated = 0;
        for (const auto& [k, c]: cnt) {
            if (c > 1)                              { if (duplicated++ < 10) std::println("[error] duplicated halfedge ({} -> {}) x{}", k.first, k.second, c); }
            if (!cnt.contains({k.second, k.first})) { if (unpaired++   < 10) std::println("[error] unpaired halfedge ({} -> {})", k.first, k.second); }
        }

        if (unpaired || duplicated) throw std::runtime_error(std::format("cut mesh validation failed: {} unpaired, {} duplicated halfedges", unpaired, duplicated));
        std::println("halfedge pairing: all {} halfedges paired", cnt.size());
    }

    visualizer::visualize_init();
    auto* base = visualizer::visualize_mesh(hm.pos, hm.idx, false, "base mesh");
    auto* embd = visualizer::visualize_mesh(hm_emb->pos, hm_emb->idx, false, "base mesh");
    visualizer::visualize_tedge(tm, mg, uv2, &X);
    visualizer::visualize_seam(*hm_emb, seam1, VecXi(), "cut seam", false);
    visualizer::visualize_non_snapped_tnodes(hm, tmm, false);
    visualizer::visualize_tedge_mut_collapsed(hm, tmm, false);
    visualizer::visualize_face_collinear_error(hm, tmm, true);

    ///--- tutte parameterization (pre-SLIM initial uv) ---///
    std::sort(hdata.begin(), hdata.end());
    MatXd uv;
    if (compute_tutte_parameterization(*hm_emb, tmm, seam1, hdata, uv)) {
        embd->addParameterizationQuantity("tutte uv", uv);
        igl::SLIMData sData;

        {
            MatXd uv_init(hm_cut->nV, 2);
            for (auto v: hm_cut->verts) uv_init.row(v.id) = uv.row(v.half().next().crnr().id);

            // pin the seam (boundary) vertices softly to the tutte uv
            std::vector<int>   b_;
            std::vector<Row2d> bc_;
            for (auto v: hm_cut->verts) {
                if (!v.isBoundary()) continue;
                b_.push_back(v.id);
                bc_.emplace_back(uv_init.row(v.id));
            }
            VecXi b = Eigen::Map<VecXi>(b_.data(), b_.size());
            MatXd bc(bc_.size(), 2);
            for (int i = 0; i < bc_.size(); ++i) bc.row(i) = bc_[i];

            sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
            slim_precompute(hm_cut->pos, hm_cut->idx, uv_init, sData, sData.slim_energy, b, bc, 1e5);
            slim_solve(sData, 10);

            std::println("[slim] displacement: {}", (sData.V_o - uv_init).norm());
            auto* surf = polyscope::registerSurfaceMesh("slim result", hm_cut->pos, hm_cut->idx);
            auto* prms = surf->addVertexParameterizationQuantity("uv", sData.V_o);
            surf->setEdgeWidth(0.7);
            prms->setEnabled(true);
            prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
            prms->setCheckerSize(1);
        }

        // ------ qex on the slim result ------
        {
            // per-corner uv from the per-vertex slim result: hm_cut and hm2
            // share the face matrix, so corner (i, j) <-> vertex hm2->idx(i, j)
            VecXc cfn(hm_emb->nF * 3);
            for (int i = 0; i < hm_cut->nF; ++i) {
            for (int j = 0; j < 3; ++j) {
                int k = hm_cut->idx(i, j);
                cfn(i * 3 + j) = complex(sData.V_o(k, 0), sData.V_o(k, 1));
            }}

            qex::sanitization(*hm_emb, matching1, singular1, 4, cfn);

            std::vector<qex::Qport> q_ports;
            std::vector<qex::Qvert> vqvs, eqvs, fqvs;
            qex::generate_q_vert(*hm_emb, cfn, vqvs, eqvs, fqvs);

            visualizer::visualize_qverts(vqvs, eqvs, fqvs, false);

            qex::generate_vqvert_qport(*hm_emb, cfn, vqvs, q_ports);
            qex::generate_eqvert_qport(*hm_emb, cfn, eqvs, q_ports);
            qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);

            visualizer::visualize_qports(*hm_emb, cfn, q_ports, false);
            auto qedges = qex::generate_q_edge(*hm_emb, cfn, matching1, q_ports);
            auto qfaces = qex::generate_q_faces(q_ports, qedges);

            visualizer::visualize_qedges(qedges);
            visualizer::visualize_qfaces(qfaces);

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
    } else std::println("[tutte] compute_tutte_parameterization failed");

    { // patch boundaries (t-mesh edges on the cut mesh), colored by tquad id
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> tqids, thids, v0s;
        size_t c = 0;
        for (const auto& d: hdata) {
            Row3d a = d.half.tail().pos();
            Row3d b = d.half.head().pos();
            ns.emplace_back(a.x(), a.y(), a.z());
            ns.emplace_back(b.x(), b.y(), b.z());
            es.push_back({c, c + 1}); c += 2;
            tqids.push_back(d.tqid);
            thids.push_back(d.thid);
            v0s.push_back(d.v0);
        }
        auto* cn = polyscope::registerCurveNetwork("patch boundaries", ns, es);
        cn->addEdgeScalarQuantity("tqid", tqids)->setEnabled(true);
        cn->addEdgeScalarQuantity("thid", thids);
        cn->addEdgeScalarQuantity("v0", v0s);
        cn->setRadius(0.0015);
    }

    polyscope::show();
    return 0;
}
