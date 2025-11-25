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

inline Mat2x3d GetAxisAlignedProjection(const Row3d& normal) {
    Row3d abs = normal.cwiseAbs();
    double max;
    Mat3x2d P;

    if (abs.z() > abs.x() && abs.z() > abs.y()) {
        // projection = mat3x2({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
        P.col(0) << 1., 0., 0.;
        P.col(1) << 0., 1., 0.;
        max = normal.z();
    }
    else if (abs.y() > abs.x()) {
        // projection = mat3x2({0.0, 0.0, 1.0}, {1.0, 0.0, 0.0});
        P.col(0) << 0., 0., 1.;
        P.col(1) << 1., 0., 0.;
        max = normal.y();
    }
    else {
        // projection = mat3x2({0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});
        P.col(0) << 0., 1., 0.;
        P.col(1) << 0., 0., 1.;
        max = normal.x();
    }

    if (max < 0.) P.col(0) *= -1.; // projection[0] *= -1.0;
    return P.transpose(); // mat3x2 -> mat2x3
}
}

struct AuxSgmtData {
    Msgmt sg;
    Thalf th;
    Row2d v0; // val cano
    Row2d v1; // val non cano
};

struct AuxHalfData {
    std::set<std::pair<double, int>> idcs; // the list of index and ratio in the Half
};

struct MidHalfData {
    int thid;
    double v0;
    double v1;
};

struct ResHalfData {
    int i0;
    int i1;
    std::optional<MidHalfData> data;
};

struct ResFaceData {
    int id; // original face id
    std::vector<std::array<ResHalfData, 3>> halfs;
};

inline int find_tqid(
    const Tmesh& tmesh,
    const int thid_
) {
    for (auto& tq: tmesh.tquads)
        for (int thid: tq.thids)
            if (thid == thid_) return tq.id;
    throw std::runtime_error("half not found in any tquad");
}


// Beware epsilon validity must be solved beforehand
inline void face_cutting(
    const Face& f,                       // face to be cut
    const VecXc& cf,                     // corner function of original mesh
    const std::vector<AuxSgmtData>& sgs, // segments inside the face
    std::vector<Row3d>& vpos,            // vertex position (altered)
    ResFaceData& f_auxs,                 //
    std::vector<AuxHalfData>& h_auxs     //
) {
    if (sgs.empty()) {
        f_auxs.id = f.id;
        Half h0 = f.half();
        Half h1 = f.half().next();
        Half h2 = f.half().prev();
        f_auxs.halfs.emplace_back(std::array{
            ResHalfData{h0.tail().id, h0.head().id},
            ResHalfData{h1.tail().id, h1.head().id},
            ResHalfData{h2.tail().id, h2.head().id}
        });
        return;
    }

    std::vector<ResHalfData> halfs;
    std::vector<int> idcs_f; // adjacent new vert idcs and newly added verts inside face and in halfs
    std::vector<std::tuple<Half, double, int>> idcs_h; //
    std::vector<std::optional<ResHalfData>> vals; // num of halfs

    for (Half h0: f.adjHalfs()) {
        idcs_f.push_back(h0.tail().id); // verts are candidate
        Half h1 = h0.twin();
        for (auto [r, i]: h_auxs[h1.id].idcs) {
            if (r > 0 && r < 1) {
                idcs_f.emplace_back(i);
                idcs_h.emplace_back(h1, r, i);
            }
        }
    }

    for (const auto& [sg, th, v0, v1]: sgs) {
        std::pair<const Mvert*, int> mvs[2] = {
            std::pair(&sg.fr, -1),
            std::pair(&sg.to, -1)
        };

        // 1: assign inside halfs
        for (auto& [mv, idx]: mvs) {
            Row3d p;
            std::optional<std::pair<Half, double>> val;
            if (mv->cut.has_value()) {
                auto [h, r] = mv->cut.value();
                p = h.tail().pos() * r + h.head().pos() * (1 - r);
                val = std::pair(h, r);
            }
            else { p = conversion_2d_3d(sg.face, cf, mv->uv); }

            auto find = rg::find_if(idcs_f, [&](int i) { return (vpos[i] - p).norm() < EPS; });
            if (find != idcs_f.end()) { idx = *find; }
            else {
                vpos.emplace_back(p);
                int l = (int)vpos.size() - 1;
                idcs_f.emplace_back(l);
                idx = l;
                if (val.has_value()) {
                    auto [h0, r] = val.value();
                    auto h1 = h0.twin();
                    h_auxs[h0.id].idcs.emplace(r, idx);
                    h_auxs[h1.id].idcs.emplace(1 - r, idx);
                }
            }
        }

        halfs.emplace_back(ResHalfData{mvs[0].second, mvs[1].second, MidHalfData{th.id  , v0.x(), v0.y()}});
        halfs.emplace_back(ResHalfData{mvs[1].second, mvs[0].second, MidHalfData{th.twid, v1.x(), v1.y()}});
    }

    // 2: assign edge halfs
    for (Half h: f.adjHalfs()) {
        const auto& idcs = h_auxs[h.id].idcs;
        if (!idcs.empty()) {
            halfs.emplace_back(ResHalfData{idcs.begin()->second, h.head().id});

            auto prev = idcs.begin();
            auto curr = std::next(idcs.begin());
            for (; curr != idcs.end(); ++curr, ++prev) {
                halfs.emplace_back(ResHalfData{curr->second, prev->second});
            }
            halfs.emplace_back(ResHalfData{h.tail().id, idcs.rbegin()->second});
        }
        else { halfs.emplace_back(ResHalfData{h.tail().id, h.head().id}); }
    }

    while (!halfs.empty()) {
        auto h = halfs.back();
        halfs.pop_back();

        std::vector poly = {h};
        int sta = h.i0;
        int cur = h.i1;

        while (true) {
            int prev_v = poly.back().i0;
            int curr_v = poly.back().i1;

            const Row3d& p_prev = vpos[prev_v];
            const Row3d& p_curr = vpos[curr_v];
            const Row3d  d_prev = (p_curr - p_prev).normalized();

            auto best_it = halfs.end();
            double score = -2;

            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (it->i0 != curr_v) continue;
                if (it->i1 == prev_v) continue;
                int cand_to = it->i1;
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

            int next = best_it->i1;
            poly.push_back(*best_it);

            cur = next;
            halfs.erase(best_it);
            if (cur == sta) break;
        }

        f_auxs.id = f.id;
        switch (int n = poly.size()) {
        case 3:
            f_auxs.halfs.emplace_back(std::array{poly[0], poly[1], poly[2]});
            break;
        case 4:
            f_auxs.halfs.emplace_back(std::array{poly[0], poly[1], poly[2]});
            f_auxs.halfs.emplace_back(std::array{poly[0], poly[2], poly[3]});
            break;
        default:
            for (int j = 1; j < n - 1; ++j)
                f_auxs.halfs.emplace_back(std::array{poly[0], poly[j], poly[j + 1]});
            break;
        }
    }
}

// Cut the original mesh for tutte parameterization.
// Consider preparing a table between the original and the cut mesh for visualize
// parameterization with the original mesh (in the next step)
//
// memo: If the epsilon-snapping of segments is needed, Better done before
// creating embedded tquad as an intermediate data set because the vertex
// positions of each tedges are consistent as the whole graph.
// OR, snap a segment-edge vertex for hmesh vertex if the distance is less than epsilon (now used)
inline void compute_embedding_cut_hmesh(
    const Hmesh& hm, // input hmesh
    const Tmesh& tm, // input tmesh
    const VecXc& cf, // input corner function of naive parameterization
    const VecXd& R   //
) {
    using namespace impl_eb;

    std::map<int, std::vector<AuxSgmtData>> cuts; // face id & aux segment data

    for (const Thalf& th: tm.thalfs) {
        if (!th.cano) continue;
        const Tedge& te = th.edge();
        const double r = R[te.id];
        double sum = 0;
        for (const auto& s: te.segments()) {
            auto len = abs(s.diff());
            auto v0 = sum / r; sum += len;
            auto v1 = sum / r;
            auto aux = AuxSgmtData{s, th, {v0, v1}, {1 - v1, 1 - v0}};
            cuts[s.face.id].emplace_back(aux);
        }
    }

    std::vector h_aux(hm.nH, AuxHalfData{});
    std::vector f_res(hm.nF, ResFaceData{});

    std::vector<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    for (Face f: hm.faces) {
        face_cutting(f, cf, cuts[f.id], vpos, f_res[f.id], h_aux);
    }

    int count = 0;
    for (const auto& [id, hs]: f_res) { count += hs.size(); }

    MatXi face_info(count, 3);
    MatXd vert_info(vpos.size(), 3);

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> val;
    size_t count1 = 0;

    int count2 = 0;
    for (const auto& [id, halfs]: f_res) {
        for (auto& hs: halfs) {
            face_info.row(count2) = Row3i{hs[0].i0, hs[1].i0, hs[2].i0};
            count2++;

            // debug bgn
            for (const auto& h: hs) {
                if (!h.data.has_value()) continue;
                auto tqid = find_tqid(tm, h.data.value().thid);
                if (tqid != 1) continue;
                Row3d p1 = vpos[h.i0];
                Row3d p2 = vpos[h.i1];
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                val.emplace_back(h.data.value().v0);
                val.emplace_back(h.data.value().v1);
                es.emplace_back(std::array{count1, count1 + 1});
                count1 += 2;
            }
            // debug end
        }
    }

    for (int i = 0; i < vpos.size(); i++) { vert_info.row(i) = vpos[i]; }
    const auto surf = polyscope::registerSurfaceMesh("new mesh!!", vert_info, face_info);
    surf->setEdgeWidth(1);

    if (!val.empty()) {
        auto c = polyscope::registerCurveNetwork("cut half data", ns, es);
        auto v = c->addNodeScalarQuantity("val", val);
        c->setEnabled(true);
        v->setEnabled(true);
        c->resetTransform();
        c->setRadius(0.002);
    }
}
}
#endif
