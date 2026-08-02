#ifndef METRIKO_EXAMPLE_COMMON_H
#define METRIKO_EXAMPLE_COMMON_H
#include "metriko/core/hmesh/utilities.h"

using namespace metriko;

namespace metriko::visualizer {
inline void visualize_frosy_field(
    polyscope::SurfaceMesh* surf,
    const Hmesh& hm,
    const FaceRosyField& rawf,
    const FaceRosyField& cmbf,
    const int rosyN = 4,
    const bool show = true
) {
    MatXd rawInt(hm.nF, 2);
    MatXd cmbInt(hm.nF, 2);
    MatXd rawExt(hm.nF, 3 * rosyN);
    MatXd cmbExt(hm.nF, 3 * rosyN);
    for (Face f: hm.faces) {
        complex rc0 = rawf.field(f.id, 0);
        complex rc1 = rawf.field(f.id, 1);
        complex rc2 = rawf.field(f.id, 2);
        complex rc3 = rawf.field(f.id, 3);
        complex cc0 = cmbf.field(f.id, 0);
        complex cc1 = cmbf.field(f.id, 1);
        complex cc2 = cmbf.field(f.id, 2);
        complex cc3 = cmbf.field(f.id, 3);
        rawInt.row(f.id) = Row2d(rc0.real(), rc0.imag()).normalized();
        cmbInt.row(f.id) = Row2d(cc0.real(), cc0.imag()).normalized();

        rawExt.block(f.id, 0, 1, 3) = (rc0.real() * f.basisX() + rc0.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 3, 1, 3) = (rc1.real() * f.basisX() + rc1.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 6, 1, 3) = (rc2.real() * f.basisX() + rc2.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 9, 1, 3) = (rc3.real() * f.basisX() + rc3.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 0, 1, 3) = (cc0.real() * f.basisX() + cc0.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 3, 1, 3) = (cc1.real() * f.basisX() + cc1.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 6, 1, 3) = (cc2.real() * f.basisX() + cc2.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 9, 1, 3) = (cc3.real() * f.basisX() + cc3.imag() * f.basisY()).normalized();
    }
    auto rawFQ = surf->addFaceVectorQuantity("raw ext", rawExt.block(0, 0, rawExt.rows(), 3));
    auto cmbFQ = surf->addFaceVectorQuantity("cmb ext", cmbExt.block(0, 0, cmbExt.rows(), 3));
    rawFQ->setEnabled(show);
    cmbFQ->setEnabled(show);
    rawFQ->setVectorLengthScale(0.004);
    cmbFQ->setVectorLengthScale(0.004);
}

inline void visualize_motorcycle_graph(
    const mc::Mgrph& mg,
    const VecXc& uv,
    bool show = true
) {
    vec<glm::vec3> l0;
    vec<glm::vec3> mnodes;
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> mcid;
    vec<double> mnid;
    vec<double> mnid_flag;
    vec<double> modes_val;

    vec<int> flag_ids = vec{
        1085,
        1091,
        1110,
        1115
    };

    vec<bool> mnodes_reserved = vec(mg.mnodes.size(), false);
    size_t counter = 0;

    auto get_loc_type_value = [](const HmLoc& loc) -> double {
        if (std::holds_alternative<HmLocOnV>(loc)) return 1.;
        if (std::holds_alternative<HmLocOnE>(loc)) return 2.;
        if (std::holds_alternative<HmLocOnP>(loc)) return 3.;
        return 0.;
    };

    for (const auto& c: mg.mcurvs) {
    for (const auto& s: c.sgmts) {
        auto iFr = s.fr_nid;
        auto iTo = s.to_nid;
        auto nFr = mg.mnodes[iFr];
        auto nTo = mg.mnodes[iTo];
        auto uv1 = mc::get_face_uv(nFr, s.face_id, mg.hm, mg.cf);
        auto uv2 = mc::get_face_uv(nTo, s.face_id, mg.hm, mg.cf);
        Row3d p1 = conversion_2d_3d(mg.hm.faces[s.face_id], uv, uv1);
        Row3d p2 = conversion_2d_3d(mg.hm.faces[s.face_id], uv, uv2);
        if (abs(uv1 - uv2) < EPS) { l0.emplace_back(p1.x(), p1.y(), p1.z()); }

        if (!mnodes_reserved[iFr]) {
            mnodes.emplace_back(p1.x(), p1.y(), p1.z());
            modes_val.emplace_back(get_loc_type_value(nFr.loc));
            mnid.emplace_back(iFr);
            mnid_flag.emplace_back(rg::contains(flag_ids, iFr) ? 1. : 0.);
            mnodes_reserved[iFr] = true;
        }

        if (!mnodes_reserved[iTo]) {
            mnodes.emplace_back(p2.x(), p2.y(), p2.z());
            modes_val.emplace_back(get_loc_type_value(nTo.loc));
            mnid.emplace_back(iTo);
            mnid_flag.emplace_back(rg::contains(flag_ids, iTo) ? 1. : 0.);
            mnodes_reserved[iTo] = true;
        }

        ns.emplace_back(p1.x(), p1.y(), p1.z());
        ns.emplace_back(p2.x(), p2.y(), p2.z());
        es.emplace_back(std::array{counter, counter + 1});
        mcid.emplace_back(c.id);
        counter += 2;
    }}

    {
        auto p = polyscope::registerPointCloud("zero len edge", l0);
        p->setMaterial("flat");
        p->setPointRadius(0.003);
    }

    {
        auto p = polyscope::registerPointCloud("mnodes", mnodes);
        p->addScalarQuantity("type", modes_val);
        p->addScalarQuantity("flag", mnid_flag);
        p->addScalarQuantity("mnid", mnid);
        p->setMaterial("flat");
        p->setPointRadius(0.003);
    }

    auto c = polyscope::registerCurveNetwork("motorcycle graph", ns, es);
    c->setColor(glm::vec4(.0, .0, .0, 1.));
    auto v_mcid = c->addEdgeScalarQuantity("mcid", mcid);
    v_mcid->setEnabled(true);
    v_mcid->setColorMap("magma");
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.001);
    c->setMaterial("flat");
}

inline void visualize_node_adjacency(const mc::Mgrph& mg, const VecXc& uv, bool show = true) {
    std::vector<glm::vec3> pts;
    std::vector<double> adj_order;
    std::vector<double> nid_list;

    for (int nid = 0; nid < mg.mnodes.size(); ++nid) {
        const auto& mn = mg.mnodes[nid];
        if (mn.adj.size() < 3) continue;
        std::vector<glm::vec3> local_pts;

        for (int i = 0; i < mn.adj.size(); ++i) {
            const auto& as = mn.adj[i];
            const auto& sg = mg.mcurvs[as.x()].sgmts[as.y()];

            bool is_outgoing = sg.fr_nid == nid;
            int fid = sg.face_id;

            Row3d pA = conversion_2d_3d(mg.hm.faces[fid], uv, mc::get_face_uv(mg.mnodes[sg.fr_nid], fid, mg.hm, mg.cf));
            Row3d pB = conversion_2d_3d(mg.hm.faces[fid], uv, mc::get_face_uv(mg.mnodes[sg.to_nid], fid, mg.hm, mg.cf));
            Row3d p0 = is_outgoing ? pA : pB;
            Row3d p1 = is_outgoing ? pB : pA;
            Row3d pt = p0 * 0.85 + p1 * 0.15;
            glm::vec3 gpt(pt.x(), pt.y(), pt.z());

            pts.push_back(gpt);
            local_pts.push_back(gpt);

            adj_order.push_back(i);
            nid_list.push_back(nid);
        }
    }

    auto pc = polyscope::registerPointCloud("CCW Adjacency Points", pts);
    pc->setEnabled(show);
    pc->setPointRadius(0.002);

    auto q_order = pc->addScalarQuantity("adj_index", adj_order);
    q_order->setEnabled(true);
    q_order->setColorMap("turbo");
    pc->addScalarQuantity("node_id", nid_list);

}

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

// =======================================================================
// 1. 細分化された TrackedDenseMesh の可視化 (3D座標の再構築を含む)
// =======================================================================
//inline void visualize_tracked_mesh(
//    const TrackedDenseMesh& dmesh,
//    const Hmesh& base_hm,
//    const VecXc& base_cf,
//    const std::string& name = "subdiv_mesh",
//    bool show = true
//) {
//    // 1. 親FaceのローカルUVから、曲面上の正確な3D座標を復元する
//    std::vector pos(dmesh.num_verts, glm::vec3(0.0f));
//    std::vector visited(dmesh.num_verts, false);
//    MatXd uv(dmesh.polygons.size() * 3, 2);
//
//    for (size_t i = 0; i < dmesh.polygons.size(); ++i) {
//        const auto& poly = dmesh.polygons[i];
//        const auto& f_uvs = dmesh.uvs[i];
//        int parent_fid = dmesh.face2parent[i];
//
//        for (size_t j = 0; j < poly.size(); ++j) {
//            uv.row(i * poly.size() + j) << f_uvs[j].real(), f_uvs[j].imag();
//
//            int vid = poly[j];
//            if (!visited[vid]) {
//                // conversion_2d_3d を使ってUV平面から3D空間へ写像
//                Row3d p = conversion_2d_3d(base_hm.faces[parent_fid], base_cf, f_uvs[j]);
//                pos[vid] = glm::vec3(p.x(), p.y(), p.z());
//                visited[vid] = true;
//            }
//        }
//    }
//
//    // 2. Polyscope用のFace配列（三角形）を構築
//    std::vector<std::array<size_t, 3>> faces;
//    faces.reserve(dmesh.polygons.size());
//    for (const auto& poly : dmesh.polygons) {
//        faces.push_back({ (size_t)poly[0], (size_t)poly[1], (size_t)poly[2] });
//    }
//
//    // 3. Polyscopeに登録
//    auto surf = polyscope::registerSurfaceMesh(name, pos, faces);
//    surf->setSurfaceColor(glm::vec3(0.8f, 0.9f, 1.0f)); // 爽やかな水色
//    surf->setEdgeWidth(1.0f);
//    surf->setEdgeColor(glm::vec3(0.2f, 0.2f, 0.2f));
//    surf->setEnabled(show);
//
//    auto prms = surf->addParameterizationQuantity("uv", uv);
//    prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
//    prms->setCheckerSize(1.);
//    prms->setEnabled(false);
//}

// ---
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

template <typename PosType, typename IdxType, typename UvType>
polyscope::SurfaceMesh* visualize_mesh_with_uv(
    const PosType& pos,
    const IdxType& idx,
    const UvType& uv,
    const std::string& name = "mesh",
    const bool show = true
) {
    auto* surf = polyscope::registerSurfaceMesh(name, pos, idx);
    auto* prms = surf->addParameterizationQuantity("uv_param", uv);
    prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    prms->setCheckerSize(1.);
    prms->setEnabled(false);

    surf->setEdgeWidth(0.7);
    surf->setEnabled(show);
    surf->setMaterial("flat");
    surf->setSurfaceColor(glm::vec3(0.15, 0.15, 0.15));

    vec<glm::vec3> f_col(surf->nFaces(), glm::vec3(0.8, 0.8, 0.8));

    f_col[18]  = glm::vec3(1, 0, 0);
    f_col[27]  = glm::vec3(1, 0, 0);
    f_col[258] = glm::vec3(1, 0, 0);
    f_col[933] = glm::vec3(1, 0, 0);
    f_col[177] = glm::vec3(1, 0, 0);
    surf->addFaceColorQuantity("face color", f_col);

    return surf;
}

inline void visualize_seam(
    const Hmesh& hm,
    const vec<bool>& seam,
    const std::string& name = "seam",
    const VecXi& matching = VecXi(),
    const bool show = true
) {
    bool use_matching = matching.rows() > 0;
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> ms;
    size_t counter = 0;
    for (auto e: hm.edges) {
        if (seam[e.id]) {
            Row3d p1 = e.half().tail().pos();
            Row3d p2 = e.half().head().pos();
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});
            if (use_matching) ms.emplace_back(matching[e.id]);
            counter += 2;
        }
    }
    auto c = polyscope::registerCurveNetwork(name, ns, es);
    if(use_matching) c->addEdgeScalarQuantity("matching", ms);
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.001);
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
#endif

