#ifndef METRIKO_EXAMPLE_COMMON_H
#define METRIKO_EXAMPLE_COMMON_H
#include "metriko/core/hmesh/utilities.h"

namespace metriko::visualizer {

inline void visualize_frosy_field(
    polyscope::SurfaceMesh* surf,
    const Hmesh& hm,
    const FaceRosyField& rawf,
    const FaceRosyField& cmbf,
    const int rosyN = 4
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
    rawFQ->setEnabled(false);
    cmbFQ->setEnabled(true);
    rawFQ->setVectorLengthScale(0.004);
    cmbFQ->setVectorLengthScale(0.004);
}

inline void visualize_motorcycle_graph(
    const mc::MotorcycleGraph& graph,
    const VecXc& uv,
    bool show = true
) {
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> mcids;
    size_t counter = 0;

    for (auto& c: graph.mcurvs) {
        for (auto& s: c.sgmts) {
            Row3d p1 = conversion_2d_3d(graph.hm.faces[s.face_id], uv, mc::get_face_uv(graph.mnodes[s.fr_nid], s.face_id, graph.hm, graph.cf));
            Row3d p2 = conversion_2d_3d(graph.hm.faces[s.face_id], uv, mc::get_face_uv(graph.mnodes[s.to_nid], s.face_id, graph.hm, graph.cf));

            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});
            mcids.emplace_back(c.id);
            counter += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork("motor cycle graph", ns, es);
    c->setColor(glm::vec4(.0, .0, .0, 1.));
    auto v_mcid = c->addEdgeScalarQuantity("mcid", mcids);
    v_mcid->setEnabled(true);
    v_mcid->setColorMap("magma");
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.001);
    c->setMaterial("flat");
}

inline void visualize_node_adjacency(const mc::MotorcycleGraph& mg, const VecXc& uv, bool show = true) {
    std::vector<glm::vec3> pts;
    std::vector<double> adj_order;
    std::vector<double> nid_list;

    for (int nid = 0; nid < mg.mnodes.size(); ++nid) {
        const auto& mn = mg.mnodes[nid];
        if (mn.adj.size() < 3) continue;
        std::vector<glm::vec3> local_pts;

        for (int i = 0; i < mn.adj.size(); ++i) {
            const auto& as = mn.adj[i];
            const auto& sg = mg.mcurvs[as.curv_id].sgmts[as.sgmt_id];

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
    const mc::MotorcycleGraph& mg, // 【追加】幾何座標の参照に必須
    const VecXc& uv,               // 3D変換用のベース頂点座標など
    const VecXd* X = nullptr,
    const std::vector<int> &selector = std::vector<int>(),
    const std::string& prefix = std::string(""),
    const bool show = true
) {
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> teids;
    std::vector<double> randoms;
    size_t counter = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, 30);

    std::vector<double> vecX;
    std::vector<double> vecR;
    std::vector<double> difx;
    std::vector<double> dify;
    std::vector<double> count;

    for (int i = 0; i < tm.nTE; i++) {
        const auto& te = tm.tedges[i];
        int random_value = distr(gen);

        // セレクタによるフィルタリング
        if (!selector.empty() && rg::find(selector, i) == selector.end()) continue;

        // Tedge が保持する Msgmt の履歴を辿る
        for (const mc::Msgmt &ts: te.segs) {
            // mg を使って、ノードIDからFaceローカルなUV座標を動的に計算する
            complex uvFr = mc::get_face_uv(mg.mnodes[ts.fr_nid], ts.face_id, mg.hm, mg.cf);
            complex uvTo = mc::get_face_uv(mg.mnodes[ts.to_nid], ts.face_id, mg.hm, mg.cf);

            // ローカルUVから3D空間座標へ変換
            Row3d p1 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvFr);
            Row3d p2 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvTo);

            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});

            teids.emplace_back(i);
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
    c->addEdgeScalarQuantity("R", vecR);
    if (X != nullptr) c->addEdgeScalarQuantity("X", vecX);
    c->addEdgeScalarQuantity("difx", difx);
    c->addEdgeScalarQuantity("dify", dify);
    c->addEdgeScalarQuantity("random", randoms);
    c->addEdgeScalarQuantity("count", count);
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.0004);
    c->setMaterial("flat");
}

inline void debug_tquad_sides(
    const Tmesh& tm,
    const mc::MotorcycleGraph& mg,
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
        std::set<int> unique_sides(tq.sides.begin(), tq.sides.end());

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
            for (int i = 0; i < tq.thids.size(); ++i) {
                int thid = tq.thids[i];
                int cid = tm.thalfs[thid].edge().curv_id;
                int side = tq.sides[i];
                std::cout << "       [" << i << "]: Thalf " << thid
                          << " [Curv " << cid << "] = Side " << side << "\n";
            }
            std::cout << "--------------------------------------------\n";
        }

        // 2. Polyscope用のジオメトリ構築 (visualize_tedge と同等)
        for (int i = 0; i < tq.thids.size(); ++i) {
            int thid = tq.thids[i];
            int side = tq.sides[i];
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
}

#endif

