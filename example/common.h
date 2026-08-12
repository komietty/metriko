#ifndef METRIKO_EXAMPLE_COMMON_H
#define METRIKO_EXAMPLE_COMMON_H
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "glm/glm.hpp"
#include "metriko/core/hmesh/utilities.h"

using namespace metriko;

namespace metriko::visualizer {
inline void visualize_tedge(
    const Tmesh& tm,
    const mc::Mgrph& mg,
    const VecXc& uv,
    const VecXd* X = nullptr,
    const std::vector<int> &selector = std::vector<int>(),
    const std::string& prefix = std::string(""),
    const bool show = true
) {
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> teids;
    vec<double> tqid1;
    vec<double> tqid2;
    vec<double> randoms;
    size_t counter = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, 30);

    vec<double> vecX;
    vec<double> vecR;
    vec<double> uvX;
    vec<double> uvY;
    vec<double> difx;
    vec<double> dify;
    vec<double> count;

    for (int i = 0; i < tm.nTE; i++) {
        const auto& te = tm.tedges[i];
        // find two thalf of tedge (created as a consecutive pair: 2*teid = cano, +1 = twin)
        const auto& th0 = tm.thalfs[te.id * 2];   // cano
        const auto& th1 = tm.thalfs[th0.twid];    // twin
        int q1 = tm.th2quad[th0.id];              // tquad on cano side
        int q2 = tm.th2quad[th1.id];              // tquad on twin side

        int random_value = distr(gen);

        if (!selector.empty() && rg::find(selector, i) == selector.end()) continue;

        for (const mc::Msgmt &ts: te.segs) {
            complex uvFr = mc::get_face_uv(mg.mnodes[ts.fr_nid], ts.face_id, mg.hm, mg.cf);
            complex uvTo = mc::get_face_uv(mg.mnodes[ts.to_nid], ts.face_id, mg.hm, mg.cf);
            Row3d p1 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvFr);
            Row3d p2 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvTo);

            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});

            teids.emplace_back(i);
            tqid1.emplace_back(q1);
            tqid2.emplace_back(q2);

            uvX.emplace_back(uvFr.real());
            uvX.emplace_back(uvTo.real());
            uvY.emplace_back(uvFr.imag());
            uvY.emplace_back(uvTo.imag());

            difx.emplace_back(uvTo.real() - uvFr.real());
            dify.emplace_back(uvTo.imag() - uvFr.imag());

            vecR.emplace_back(te.len);
            if (X != nullptr) vecX.emplace_back((*X)[i]);

            randoms.emplace_back(random_value);
            count.emplace_back(counter);
            counter += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork(prefix + "tedges", ns, es);
    c->setColor(glm::vec4(.0, .0, .0, 1.));
    c->addEdgeScalarQuantity("teid", teids);
    c->addEdgeScalarQuantity("tqid1", tqid1);
    c->addEdgeScalarQuantity("tqid2", tqid2);
    c->addEdgeScalarQuantity("R", vecR);
    if (X != nullptr) c->addEdgeScalarQuantity("X", vecX);
    c->addNodeScalarQuantity("uv x", uvX);
    c->addNodeScalarQuantity("uv y", uvY);
    //c->addEdgeScalarQuantity("dif x", difx);
    //c->addEdgeScalarQuantity("dif y", dify);
    //c->addEdgeScalarQuantity("random", randoms);
    //c->addEdgeScalarQuantity("count", count);
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.0005);
    c->setMaterial("flat");
}

inline void debug_tquad_sides(
    const Tmesh& tm,
    const mc::Mgrph& mg,
    const VecXc& uv,
    const std::string& prefix = "",
    const bool show = true
) {
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> quad_ids;
    std::vector<double> side_ids;
    std::vector<double> is_invalid; // 異常なQuadを赤く光らせる用
    size_t counter = 0;

    bool all_valid = true;

    std::cout << "========== Tquad Validation Start ==========" << std::endl;

    for (const auto& tq : tm.tquads) {
        // 1. Quadが持つ一意なsideを収集
        std::set<int> unique_sides;
        for (const auto& d : tq.data) unique_sides.insert(d.side);

        // 0, 1, 2, 3 が全て揃っているかチェック
        bool valid = (unique_sides.size() == 4) &&
                     unique_sides.count(0) && unique_sides.count(1) &&
                     unique_sides.count(2) && unique_sides.count(3);

        // 異常なQuadをコンソールに詳細出力
        if (!valid) {
            all_valid = false;
            std::cout << "[ERROR] Tquad ID: " << tq.id << " is invalid!\n";
            std::cout << "  -> Unique sides found: ";
            for (int s : unique_sides) std::cout << s << " ";
            std::cout << "\n";

            std::cout << "  -> Details (Index: thid [curv_id] = side):\n";
            for (int i = 0; i < tq.data.size(); ++i) {
                int thid = tq.data[i].thid;
                int cid = tm.thalfs[thid].edge().crv_id;
                int side = tq.data[i].side;
                std::cout << "       [" << i << "]: Thalf " << thid
                          << " [Curv " << cid << "] = Side " << side << "\n";
            }
            std::cout << "--------------------------------------------\n";
        }

        // 2. Polyscope用のジオメトリ構築 (visualize_tedge と同等)
        for (int i = 0; i < tq.data.size(); ++i) {
            int thid = tq.data[i].thid;
            int side = tq.data[i].side;
            const auto& th = tm.thalfs[thid];
            const auto& te = th.edge();

            for (const mc::Msgmt &ts: te.segs) {
                complex uvFr = mc::get_face_uv(mg.mnodes[ts.fr_nid], ts.face_id, mg.hm, mg.cf);
                complex uvTo = mc::get_face_uv(mg.mnodes[ts.to_nid], ts.face_id, mg.hm, mg.cf);

                Row3d p1 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvFr);
                Row3d p2 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvTo);

                if (!valid) {
                    std::cout << "len: " << (p2 - p1).norm() << std::endl;
                    std::cout << "p1: " << p1.transpose() << std::endl;
                    std::cout << "p2: " << p2.transpose() << std::endl;
                    std::cout << "--------------------------------------------\n";
                }

                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array<size_t, 2>{counter, counter + 1});

                quad_ids.emplace_back(tq.id);
                side_ids.emplace_back(side);
                is_invalid.emplace_back(valid ? 0.0 : 1.0); // 異常なら 1.0

                counter += 2;
            }
        }
    }

    if (all_valid) {
        std::cout << "[SUCCESS] All Tquads strictly contain sides 0, 1, 2, and 3.\n";
    }
    std::cout << "============================================" << std::endl;

    // 3. Polyscopeへの登録
    auto c = polyscope::registerCurveNetwork(prefix + "tquad_validation", ns, es);
    c->setColor(glm::vec4(0.2, 0.2, 0.2, 1.0)); // 基本は暗いグレー
    c->addEdgeScalarQuantity("quad_id", quad_ids);
    c->addEdgeScalarQuantity("side_id", side_ids);

    // 異常なQuadだけを強調表示するためのQuantity
    auto invalid_q = c->addEdgeScalarQuantity("is_invalid", is_invalid);
    invalid_q->setEnabled(true); // 自動的にこのカラーマップを表示

    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.0006); // 若干太めにして見やすく
}

inline void visualize_mapped_mnodes(
    const std::map<int, int>& mnode2dense_v,
    const Hmesh& dense_hm,
    const std::string& name = "mapped_mnodes"
) {
    std::vector<glm::vec3> pts;
    std::vector<double> nid_vals;

    for (const auto& [nid, dense_vid]: mnode2dense_v) {
        auto p = dense_hm.verts[dense_vid].pos();
        pts.emplace_back(p.x(), p.y(), p.z());
        nid_vals.push_back(nid);
    }

    auto pc = polyscope::registerPointCloud(name, pts);
    pc->addScalarQuantity("mnode_id", nid_vals);
    pc->setPointRadius(0.002);
    pc->setMaterial("flat");
}

template <typename Container>
inline void visualize_half_data(
    const Container& half_data_list,
    int target_tqid = -1, // -1 を指定するとすべての tqid を表示
    const std::string& name = "half_data"
) {
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;

    // エッジ（線分）に紐づくデータ
    std::vector<double> thids;
    std::vector<double> tqids;
    std::vector<double> orders;
    std::vector<double> is_first;
    std::vector<double> is_crash;

    // ノード（頂点）に紐づくデータ (v0, v1)
    std::vector<double> v_vals;

    size_t counter = 0;

    for (const auto& data : half_data_list) {
        // ターゲットの tqid でフィルタリング
        if (target_tqid != -1 && data.tqid != target_tqid) continue;

        // Hmeshのハーフエッジから両端の3D座標を取得
        Row3d p1 = data.half.tail().pos();
        Row3d p2 = data.half.head().pos();

        nodes.emplace_back(p1.x(), p1.y(), p1.z());
        nodes.emplace_back(p2.x(), p2.y(), p2.z());
        edges.push_back({counter, counter + 1});

        // スカラー値の記録
        thids.push_back(data.thid);
        tqids.push_back(data.tqid);
        orders.push_back(data.order);
        is_first.push_back(data.first ? 1.0 : 0.0);
        is_crash.push_back(data.crash ? 1.0 : 0.0);

        // 始点(tail)に v0、終点(head)に v1 を割り当てる
        v_vals.push_back(data.v0);
        v_vals.push_back(data.v1);

        counter += 2;
    }

    if (nodes.empty()) {
        std::cout << "[visualizer] No HalfData found for tqid: " << target_tqid << std::endl;
        return;
    }

    // Polyscope への登録
    auto* net = polyscope::registerCurveNetwork(name, nodes, edges);
    net->setRadius(0.001); // 見やすいように少し太め

    // エッジアトリビュートの追加
    net->addEdgeScalarQuantity("thid", thids);
    net->addEdgeScalarQuantity("tqid", tqids);
    net->addEdgeScalarQuantity("order", orders);
    net->addEdgeScalarQuantity("is_first", is_first);

    // crash しているエッジは赤色系のカラーマップにして目立たせる
    auto* q_crash = net->addEdgeScalarQuantity("is_crash", is_crash);
    q_crash->setColorMap("reds");

    // ノードアトリビュートの追加 (Polyscope が v0 から v1 へ自動で色を補間してくれます)
    net->addNodeScalarQuantity("v_val (v0->v1)", v_vals);

    std::cout << "[visualizer] Rendered " << (nodes.size() / 2)
              << " HalfData segments for tqid: " << target_tqid << std::endl;
}

// debug views of a collapsed TmeshMut: collapsed tedge polylines, tnodes not
// yet snapped to a vertex, and faces violating the collapse_valid_snap_1 rule
inline void visualize_tmesh_mut(
    const Hmesh& hm,
    const TmeshMut& tmm
) {
    // tnodes not snapped to a vertex, by carrier type
    {
        std::vector<glm::vec3> ps;
        std::vector<double> type, ids;
        std::set<int> seen;   // shared nodes (junctions/crossings) appear in several chains
        for (const auto& [teid, nids]: tmm.tedges) {
            if (teid == -1) continue;
            for (int nid: nids) {
                if (!seen.insert(nid).second) continue;
                double t = -1, id = -1;
                std::visit(overloaded{
                    [&](const HmLocOnE& e) { t = 0; id = e.id; },
                    [&](const HmLocOnH& h) { t = 1; id = h.id; },
                    [&](const HmLocOnF& f) { t = 2; id = f.id; },
                    [&](const auto&)       {},
                }, tmm.tnodes[nid]);
                if (t < 0) continue;   // OnV: snapped, skip
                Row3d p = get_ptloc_pos(hm, tmm.tnodes[nid]);
                ps.emplace_back(p.x(), p.y(), p.z());
                type.push_back(t);
                ids.push_back(id);
            }
        }
        auto* pc = polyscope::registerPointCloud("unsnapped tnodes", ps);
        pc->addScalarQuantity("type (0:E 1:H 2:F)", type)->setEnabled(true);
        pc->addScalarQuantity("carrier id", ids);
        pc->setPointRadius(0.002);
    }

    // collapsed tedges as polylines, colored by teid
    {
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        std::vector<double> ids;
        size_t c = 0;
        for (const auto& [id, nids]: tmm.live_tedges()) {
            for (size_t k = 0; k + 1 < nids.size(); ++k) {
                Row3d a = get_ptloc_pos(hm, tmm.tnodes[nids[k]]);
                Row3d b = get_ptloc_pos(hm, tmm.tnodes[nids[k + 1]]);
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1});
                c += 2;
                ids.push_back(id);
            }
        }
        auto* cn = polyscope::registerCurveNetwork("tedges_collapsed", ns, es);
        cn->addEdgeScalarQuantity("teid", ids)->setEnabled(true);
        cn->setRadius(0.0015);
        cn->resetTransform();
    }

    // faces holding 3+ same-side vertex-snapped nodes (collapse_valid_snap_1)
    {
        std::set<int> bad;   // hm face ids violating the snap_1 criterion
        for (auto& tq: tmm.live_tquads()) {
        for (int side = 0; side < 4; side++) {
            std::set<int>  nids;
            umap<int, int> count;
            for (int thid: tq.thids(side))
            for (int nid: tmm.tedges[tmm.thalfs[thid].teid].nids) nids.insert(nid);
            for (int nid: nids) {
                if (auto* l = std::get_if<HmLocOnV>(&tmm.tnodes[nid]))
                    for (Face f: hm.verts[l->id].adjHalfs() | vw::transform(&Half::face)) count[f.id]++;
            }
            for (auto& [fid, c]: count) if (c >= 3) bad.insert(fid);
        }}
        std::println("[collinear] {} faces violate snap_1", bad.size());

        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2>> es;
        size_t c = 0;
        for (int fid: bad) {
            for (Half h: hm.faces[fid].adjHalfs()) {
                Row3d a = h.tail().pos();
                Row3d b = h.head().pos();
                ns.emplace_back(a.x(), a.y(), a.z());
                ns.emplace_back(b.x(), b.y(), b.z());
                es.push_back({c, c + 1});
                c += 2;
            }
        }
        if (!ns.empty()) {
            auto* cn = polyscope::registerCurveNetwork("collinear faces", ns, es);
            cn->setColor({1., 0.2, 0.1});
            cn->setRadius(0.0015);
            cn->resetTransform();
        }
    }
}
}

static bool load_cache(const std::string& p, VecXc& uv2, VecXi& matching, VecXi& singular, std::vector<bool>& seam) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    int64_t nu, nm, ns, ne;
    f.read((char*)&nu, 8); f.read((char*)&nm, 8); f.read((char*)&ns, 8); f.read((char*)&ne, 8);
    uv2.resize(nu); matching.resize(nm); singular.resize(ns);
    f.read((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.read((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.read((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(ne);
    f.read(sb.data(), ne);
    seam.assign(sb.begin(), sb.end());
    return (bool)f;
}

static void save_cache(const std::string& p, const VecXc& uv2, const VecXi& matching, const VecXi& singular, const std::vector<bool>& seam) {
    std::ofstream f(p, std::ios::binary);
    int64_t nu = uv2.size(), nm = matching.size(), ns = singular.size(), ne = (int64_t)seam.size();
    f.write((char*)&nu, 8); f.write((char*)&nm, 8); f.write((char*)&ns, 8); f.write((char*)&ne, 8);
    f.write((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.write((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.write((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(seam.begin(), seam.end());
    f.write(sb.data(), ne);
}

// binary snapshot of a (collapsed) TmeshMut, so downstream demos can skip the
// motorcycle graph / quantization / collapse / snap stages
static void save_tmm(const std::string& p, const TmeshMut& tm) {
    std::ofstream f(p, std::ios::binary);
    auto wi = [&](int v)    { f.write((char*)&v, 4); };
    auto wd = [&](double v) { f.write((char*)&v, 8); };

    wi((int)tm.tnodes.size());
    for (const HmLoc& l: tm.tnodes) {
        wi((int)l.index());   // 0:OnV 1:OnE 2:OnF 3:OnH 4:OnC 5:OnP (variant order)
        std::visit(overloaded{
            [&](const HmLocOnV& v) { wi(v.id); },
            [&](const HmLocOnC& v) { wi(v.id); },
            [&](const HmLocOnE& v) { wi(v.id); wd(v.r); },
            [&](const HmLocOnH& v) { wi(v.id); wd(v.r); },
            [&](const HmLocOnF& v) { wi(v.id); wd(v.xy.real()); wd(v.xy.imag()); },
            [&](const HmLocOnP& v) { wi(v.id); wd(v.uv.real()); wd(v.uv.imag()); },
        }, l);
    }
    wi((int)tm.tedges.size());
    for (auto& te: tm.tedges) {
        wi(te.id);
        wi((int)te.nids.size());
        f.write((char*)te.nids.data(), te.nids.size() * 4);
    }
    wi((int)tm.thalfs.size());
    for (auto& th: tm.thalfs) {
        wi(th.id); wi(th.twid); wi(th.teid); wi(th.tqid);
        wi(th.cano); wi(th.bgn); wi(th.end);
        wd(th.x); wd(th.r);
    }
    wi((int)tm.tquads.size());
    for (auto& tq: tm.tquads) {
        wi(tq.id);
        wi((int)tq.data.size());
        for (auto& d: tq.data) { wi(d.thid); wi(d.side); }
    }
}

// fills a TmeshMut shell (constructed as TmeshMut(hm)). do not move the object
// afterwards: the thalf back-pointers are bound here
static bool load_tmm(const std::string& p, TmeshMut& tm) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    auto ri = [&]() { int v = 0;    f.read((char*)&v, 4); return v; };
    auto rd = [&]() { double v = 0; f.read((char*)&v, 8); return v; };

    tm.tnodes.clear(); tm.tedges.clear(); tm.thalfs.clear(); tm.tquads.clear();

    for (int n = ri(), i = 0; i < n; ++i) {
        int tag = ri(), id = ri();
        switch (tag) {
            case 0: tm.tnodes.emplace_back(HmLocOnV{id}); break;
            case 1: tm.tnodes.emplace_back(HmLocOnE{id, rd()}); break;
            case 2: { double x = rd(), y = rd(); tm.tnodes.emplace_back(HmLocOnF{id, {x, y}}); break; }
            case 3: tm.tnodes.emplace_back(HmLocOnH{id, rd()}); break;
            case 4: tm.tnodes.emplace_back(HmLocOnC{id}); break;
            case 5: { double x = rd(), y = rd(); tm.tnodes.emplace_back(HmLocOnP{id, {x, y}}); break; }
            default: return false;
        }
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        TedgeMut te{.id = ri()};
        te.nids.resize(ri());
        f.read((char*)te.nids.data(), te.nids.size() * 4);
        tm.tedges.push_back(std::move(te));
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        ThalfMut th{.tm = &tm};
        th.id   = ri(); th.twid = ri(); th.teid = ri(); th.tqid = ri();
        th.cano = ri(); th.bgn  = ri(); th.end  = ri();
        th.x    = rd(); th.r    = rd();
        tm.thalfs.push_back(th);
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        TquadMut tq{.id = ri()};
        tq.data.resize(ri());
        for (auto& d: tq.data) { d.thid = ri(); d.side = ri(); }
        tm.tquads.push_back(std::move(tq));
    }
    return (bool)f;
}
#endif

