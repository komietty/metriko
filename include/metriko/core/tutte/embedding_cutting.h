//
// Created by saki on 2025/11/17.
//

#ifndef METRIKO_EXAMPLE_EMBEDDING_CUTTING_H
#define METRIKO_EXAMPLE_EMBEDDING_CUTTING_H
#include <set>
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh.h"

namespace metriko {
namespace impl_eb {
constexpr double EPS = 1e-10;

using Mat3x2d = Eigen::Matrix<double, 3, 2>;
using Mat2x3d = Eigen::Matrix<double, 2, 3>;

inline Mat2x3d GetAxisAlignedProjection(const Row3d &normal) {
    Row3d abs = normal.cwiseAbs();
    double max;
    Mat3x2d P;

    if (abs.z() > abs.x() && abs.z() > abs.y()) {
        // projection = mat3x2({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
        P.col(0) << 1., 0., 0.;
        P.col(1) << 0., 1., 0.;
        max = normal.z();
    } else if (abs.y() > abs.x()) {
        // projection = mat3x2({0.0, 0.0, 1.0}, {1.0, 0.0, 0.0});
        P.col(0) << 0., 0., 1.;
        P.col(1) << 1., 0., 0.;
        max = normal.y();
    } else {
        // projection = mat3x2({0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
        P.col(0) << 0., 1., 0.;
        P.col(1) << 0., 0., 1.;
        max = normal.x();
    }

    if (max < 0.) P.col(0) *= -1.; // projection[0] *= -1.0;
    return P.transpose(); // mat3x2 -> mat2x3
}

struct Polygon {
    std::vector<int> indices;
};
}

struct EmbeddingHalf {
    Half *half; // reference to the cut hmesh
    Row2d uv0; // uv for tail
    Row2d uv1; // uv for head
    int tqid; // a tracker to the tquad. Usable for visualization
    //int thid; // (maybe later use) a tracker to the thalf
};

struct AuxVertData {
    int vid; // vertex id
};

struct AuxHalfData {
    std::set<std::pair<double, int> > idcs; // the list of indices in the Half
};

struct AuxFaceData {
    int id; // face id
    std::vector<Row3i> idcs; // the triangle idcs data
};


// Beware epsilon validity must be cleared beforehand
inline void face_cutting(
    const Face &f, // face to be cut
    const VecXc &cf, // corner function of original mesh
    const std::vector<Msgmt> &sgs, // segments inside the face
    std::vector<Row3d> &vpos, // vertex position (altered)
    AuxFaceData &f_auxs, //
    std::vector<AuxHalfData> &h_auxs //
) {
    //for (AuxHalfData &ha: h_auxs) {
    //    for (int i: ha.idcs) std::cout << i << ", ";
    //    std::cout << std::endl;
    //}

    if (sgs.empty()) {
        f_auxs.id = f.id;
        f_auxs.idcs = {Row3i(
            f.half().tail().id,
            f.half().next().tail().id,
            f.half().prev().tail().id
        )};
        return;
    }

    // 1st: assigning halfedges
    std::vector<std::tuple<int, int> > halfs; // vid, vid, used flag set
    std::vector<int> idcs_f; // adjacent new vert idcs and newly added verts inside face and in halfs
    std::vector<std::tuple<Half, double, int> > idcs_h; //

    for (Half h0: f.adjHalfs()) {
        idcs_f.push_back(h0.tail().id); // verts are candidate

        //
        Half h1 = h0.twin();
        for (auto [r, i]: h_auxs[h1.id].idcs) {
            if (r > 0 && r < 1) {
                idcs_f.emplace_back(i);
                idcs_h.emplace_back(h1, r, i);
            }
        }
    }

    for (const Msgmt &s: sgs) {
        std::pair<const Mvert *, int> mvs[2] = {
            std::pair(&s.fr, -1),
            std::pair(&s.to, -1)
        };

        for (auto &[ptr, idx]: mvs) {
            const Mvert &mv = *ptr;
            Row3d p;
            std::optional<std::pair<Half, double> > val;
            if (mv.cut.has_value()) {
                auto [h, r] = mv.cut.value();
                p = h.tail().pos() * r + h.head().pos() * (1 - r);
                val = std::pair(h, r);
            } else { p = conversion_2d_3d(s.face, cf, mv.uv); }

            //std::cout << "p: " << p << std::endl;

            auto find = rg::find_if(idcs_f, [&](int i) { return (vpos[i] - p).norm() < EPS; });
            if (find != idcs_f.end()) { idx = *find; } else {
                vpos.emplace_back(p);
                int l = (int) vpos.size() - 1;
                idcs_f.emplace_back(l);
                idx = l;
                if (val.has_value()) {
                    h_auxs[val->first.id].idcs.emplace(std::pair(val->second, idx));
                    h_auxs[val->first.twin().id].idcs.emplace(std::pair(1 - val->second, idx));
                }
            }
        }

        halfs.emplace_back(mvs[0].second, mvs[1].second);
        halfs.emplace_back(mvs[1].second, mvs[0].second); // need invert
    }

    for (Half h: f.adjHalfs()) {
        const auto &idcs = h_auxs[h.id].idcs;
        if (idcs.size() > 0) {
            halfs.emplace_back(idcs.begin()->second, h.head().id);

            auto prev = idcs.begin();
            auto curr = std::next(idcs.begin());
            for (; curr != idcs.end(); ++curr, ++prev) {
                halfs.emplace_back(curr->second, prev->second);
            }
            halfs.emplace_back(h.tail().id, idcs.rbegin()->second);
        } else { halfs.emplace_back(h.tail().id, h.head().id); }
    }

    for (auto [v0, v1]: halfs) {
        std::cout << "h: " << v0 << ", " << v1 << std::endl;
    }

    std::vector<std::vector<int> > polygons;

    while (!halfs.empty()) {
        auto [v0, v1] = halfs.back();
        halfs.pop_back();

        std::vector poly = {v0, v1};
        int sta = v0;
        int cur = v1;

        while (true) {
            int prev_v = poly[poly.size() - 2];
            int curr_v = poly.back();

            Row3d p_prev = vpos[prev_v];
            Row3d p_curr = vpos[curr_v];
            Row3d d_prev = (p_curr - p_prev).normalized();

            auto best_it = halfs.end();
            double score = -2.0;

            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (prev_v == 164 && curr_v == 163 && std::get<0>(*it) == 163 && std::get<1>(*it) == 162) {
                    std::cout << "prev_v: " << prev_v << ", curr_v: " << curr_v << std::endl;
                }

                if (std::get<0>(*it) != curr_v) continue;
                if (std::get<1>(*it) == prev_v) continue;

                int cand_to = std::get<1>(*it);
                Row3d d_cand = (vpos[cand_to] - p_curr).normalized();

                // 向き判定：d_prev から d_cand への回転が CCW かどうか
                Row3d c = d_prev.cross(d_cand); // 回転軸方向
                double sign = f.normal().dot(c); // face 法線と同じ向きなら CCW
                if (c.norm() > EPS && sign <= 0) { continue; }

                double s = (-d_prev).dot(d_cand);
                if (s > score) {
                    score = s;
                    best_it = it;
                }
            }

            // CCW 候補が無ければチェーン終了
            if (best_it == halfs.end()) break;

            int next = std::get<1>(*best_it);
            poly.push_back(next);

            cur = next;
            halfs.erase(best_it);
            if (cur == sta) break;
        }

        polygons.push_back(poly);
    }
    /*

    for (int i = 0; i < halfs.size(); i++) {
        std::cout << "h: " << std::get<0>(halfs[i]) << ", " << std::get<1>(halfs[i]) << std::endl;
        auto h = halfs[i];
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        size_t counter = 0;

        Row3d p1 = vpos[std::get<0>(h)];
        Row3d p2 = vpos[std::get<1>(h)];

        ns.emplace_back(p1.x(), p1.y(), p1.z());
        ns.emplace_back(p2.x(), p2.y(), p2.z());
        es.emplace_back(std::array{counter, counter + 1});
        counter += 2;

        auto c = polyscope::registerCurveNetwork("temp-" + std::to_string(f.id) + ", " + std::to_string(i), ns, es);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.001);
    }
    */

    /*
    for (int i = 0; i < polygons.size(); i++) {
        auto poly = polygons[i];
        std::vector<glm::vec3> ns;
        std::vector<std::array<size_t, 2> > es;
        size_t counter = 0;
        for (int j = 0; j < poly.size() - 1; j++) {
            Row3d p1 = vpos[poly[j]];
            Row3d p2 = vpos[poly[j + 1]];
            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});
            counter += 2;
        }
        auto c = polyscope::registerCurveNetwork("poly-" + std::to_string(f.id) + ", " + std::to_string(i), ns, es);
        c->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.001);
    }
    */

    for (int i = 0; i < polygons.size(); i++) {
        auto poly = polygons[i];
        f_auxs.id = f.id;
        if (poly.size() == 4) {
            f_auxs.idcs.emplace_back(Row3i(poly[0], poly[1], poly[2]));
        } else if (poly.size() == 5) {
            f_auxs.idcs.emplace_back(Row3i(poly[0], poly[1], poly[2]));
            f_auxs.idcs.emplace_back(Row3i(poly[0], poly[2], poly[3]));
        } else {
            int unique_n = static_cast<int>(poly.size()) - 1;
            for (int j = 1; j + 1 < unique_n; ++j) {
                f_auxs.idcs.emplace_back(Row3i(poly[0], poly[j], poly[j + 1]));
            }
        }
    }
    /*
    */

    // newly added verts
    std::vector<glm::vec3> vis_port;
    std::vector<double> idx;
    for (int i = f.m->nV; i < vpos.size(); i++) {
        idx.emplace_back(i);
        Row3d p = vpos[i];
        vis_port.emplace_back(p.x(), p.y(), p.z());
    }


    auto vq = polyscope::registerPointCloud("newly added verts-" + std::to_string(f.id), vis_port);
    vq->setEnabled(false);
    vq->addScalarQuantity("idx", idx);
    vq->setPointRadius(0.002);
    vq->resetTransform();

    // 2nd: triangulate
}

// Cut the original mesh for tutte parameterization.
// Consider preparing a table between the original and the cut mesh for visualize
// parameterization with the original mesh (in the next step)
//
// memo: If the epsilon-snapping of segments is needed, Better done before
// creating embedded tquad as an intermediate data set because the vertex
// positions of each tedges are consistent as the whole graph.
// OR, snap a segment-edge vertex for hmesh vertex if the distance is less than epsilon (now used)
inline void embedding_cutting(
    const Hmesh &hm_in, // input hmesh
    const VecXc &cf, // input corner function of naive paramenterization
    const Tmesh &tm_in // input tmesh
    //const Hmesh &hm_out // output hmesh
) {
    using namespace impl_eb;

    std::map<int, std::vector<Msgmt> > cuts; // face id & copy of Msgmt

    std::vector<Row3d> vpos;
    for (auto p: hm_in.pos.rowwise()) { vpos.emplace_back(p); }

    // Assign segments with face index
    for (const Thalf &th: tm_in.thalfs) {
        if (!th.cano) continue;
        const Tedge &te = th.edge();
        for (const auto &s: te.segments()) {
            cuts[s.face.id].emplace_back(s);
        }
    }


    std::vector f_auxs(hm_in.nF, AuxFaceData{});
    std::vector h_auxs(hm_in.nH, AuxHalfData{});


    for (Face f: hm_in.faces) {
        //if (
        ////    f.id != 205 &&
        ////    f.id != 231 &&
        ////    f.id != 232 &&
        ////    f.id != 166 &&
        ////    f.id != 54 &&
        ////    f.id != 66 &&
        ////    f.id != 57 &&
        ////    f.id != 60 &&
        ////    f.id != 56 &&
        ////    f.id != 55 &&
        ////    f.id != 52
        //    f.id != 300
        //)
        //    continue;
        if (!cuts.contains(f.id)) { face_cutting(f, cf, {}, vpos, f_auxs[f.id], h_auxs); }
        face_cutting(f, cf, cuts[f.id], vpos, f_auxs[f.id], h_auxs);
    }

    int count = 0;
    for (AuxFaceData &fa: f_auxs) {
        count += fa.idcs.size();
    }


    MatXi face_info(count, 3);
    MatXd vert_info(vpos.size(), 3);

    int count2 = 0;
    for (AuxFaceData &fa: f_auxs) {
        for (Row3i &tri: fa.idcs) {
            face_info.row(count2) = tri;
            count2++;
        }
    }

    for (int i = 0; i < vpos.size(); i++) {
        vert_info.row(i) = vpos[i];
    }

    const auto surf = polyscope::registerSurfaceMesh("new mesh!!", vert_info, face_info);
    surf->setEdgeWidth(1);

    // debug bgn
    /*
    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    size_t counter = 0;
    std::vector<double> hasH;
    std::vector<double> sgid;
    std::vector<double> cvid;

    for (auto val: cuts[57]) {
        Row3d p1 = conversion_2d_3d(val.face, cf, val.fr.uv);
        Row3d p2 = conversion_2d_3d(val.face, cf, val.to.uv);
        ns.emplace_back(p1.x(), p1.y(), p1.z());
        ns.emplace_back(p2.x(), p2.y(), p2.z());
        hasH.emplace_back(val.fr.cut.has_value());
        hasH.emplace_back(val.to.cut.has_value());
        es.emplace_back(std::array{counter, counter + 1});
        counter += 2;
    }
    auto c = polyscope::registerCurveNetwork("segments_cut", ns, es);
    c->addNodeScalarQuantity("hasH", hasH);
    c->setEnabled(true);
    c->resetTransform();
    c->setRadius(0.0005);
    */

    // debug end
}
}

#endif
