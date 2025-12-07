//
// Created by saki on 2025/11/17.
//

#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <set>
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh.h"

namespace metriko::tutte {

template <class T> using vec = std::vector<T>;
template <class T> using set = std::set<T>;

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

struct HalfData {
    Half half;
    double v0;
    double v1;
    int thid;
    int tqid;
    int order;
    bool first = false;
    bool crash = false;

    bool operator<(const HalfData& rhs) const noexcept {
        return std::tie(tqid, thid, order) < std::tie(rhs.tqid, rhs.thid, rhs.order);
    }
};

struct AuxSgmt {
    mc::Msgmt sg;
    Thalf th;
    int ord0; // order cano
    int ord1; // order non cano
    Row2d v0; // val cano
    Row2d v1; // val non cano
};

struct AuxHalf2 {
    int i0;
    int i1;
    std::optional<Half> original;
    std::optional<HalfData> data;
};

using AuxHalf1 = std::set<std::pair<double, int>>; // the list of index and ratio in the Half
using AuxFace  = vec<std::array<AuxHalf2, 3>>;     // per original face, this contains halfedge data of a divided triangle

// Beware epsilon validity must be solved beforehand
inline void face_cutting(
    const Tmesh& tm,         //
    const Face& f,           // face to be cut
    const VecXc& cf,         // corner function of original mesh
    const vec<AuxSgmt>& sgs, // segments inside the face
    vec<Row3d>& vpos,        // vertex position (altered)
    AuxFace& f_auxs,         //
    vec<AuxHalf1>& h_auxs    //
) {
    if (sgs.empty()) {
        Half h0 = f.half();
        Half h1 = f.half().next();
        Half h2 = f.half().prev();
        f_auxs.emplace_back(std::array{
            AuxHalf2{h0.tail().id, h0.head().id, h0},
            AuxHalf2{h1.tail().id, h1.head().id, h1},
            AuxHalf2{h2.tail().id, h2.head().id, h2}
        });
        return;
    }

    vec<AuxHalf2> halfs;
    vec<int> idcs; // adjacent new vert idcs and newly added verts inside face and in halfs

    for (Half h0: f.adjHalfs()) {
        idcs.push_back(h0.tail().id); // verts are candidate
        Half h1 = h0.twin();
        for (auto [r, i]: h_auxs[h1.id]) {
            if (r > 0 && r < 1) { idcs.emplace_back(i); }
        }
    }

    for (const auto& [sg, th, o0, o1, v0, v1]: sgs) {
        std::pair<const mc::Mvert*, int> mvs[2] = {
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

            auto find = rg::find_if(idcs, [&](int i) { return (vpos[i] - p).norm() < EPS; });
            if (find != idcs.end()) { idx = *find; }
            else {
                vpos.emplace_back(p);
                int l = (int)vpos.size() - 1;
                idcs.emplace_back(l);
                idx = l;
                if (val.has_value()) {
                    auto [h0, r] = val.value();
                    auto h1 = h0.twin();
                    h_auxs[h0.id].emplace(r, idx);
                    h_auxs[h1.id].emplace(1 - r, idx);
                }
            }
        }

        auto [mv0, i0] = mvs[0];
        auto [mv1, i1] = mvs[1];
        bool f0 = mv0->type == mc::MvertType::First;
        bool f1 = mv1->type == mc::MvertType::First;
        bool b0 = mv0->type == mc::MvertType::HitB;
        bool b1 = mv1->type == mc::MvertType::HitB;
        halfs.emplace_back(AuxHalf2{i0, i1, std::nullopt, HalfData{Half(), v0.x(), v0.y(), th.id  , tm.th2quad(th.id)  , o0, f0, b1 }});
        halfs.emplace_back(AuxHalf2{i1, i0, std::nullopt, HalfData{Half(), v1.x(), v1.y(), th.twid, tm.th2quad(th.twid), o1, f1, b0 }});
    }

    // 2: assign edge halfs
    for (Half h: f.adjHalfs()) {
        const auto& idcs_ = h_auxs[h.id];
        if (idcs_.empty()) { halfs.emplace_back(h.tail().id, h.head().id, h); }
        else {
            halfs.emplace_back(idcs_.begin()->second, h.head().id, h);

            auto prev = idcs_.begin();
            auto curr = std::next(idcs_.begin());
            for (; curr != idcs_.end(); ++curr, ++prev) {
                halfs.emplace_back(curr->second, prev->second, h);
            }
            halfs.emplace_back(h.tail().id, idcs_.rbegin()->second, h);
        }
    }

    while (!halfs.empty()) {
        auto h = halfs.back();
        halfs.pop_back();

        vec poly = {h};
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
                Row3d c = d_prev.cross(d_cand);  // 回転軸方向
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

        switch (int n = poly.size()) {
        case 3:
            f_auxs.emplace_back(std::array{poly[0], poly[1], poly[2]});
            break;
        case 4:
            f_auxs.emplace_back(std::array{poly[0], poly[1], poly[2]});
            f_auxs.emplace_back(std::array{poly[0], poly[2], poly[3]});
            break;
        default:
            for (int j = 1; j < n - 1; ++j)
                f_auxs.emplace_back(std::array{poly[0], poly[j], poly[j + 1]});
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
inline Hmesh compute_embedding_cut_hmesh(
    const Hmesh& hm,           // input hmesh
    const Tmesh& tm,           // input tmesh
    const VecXc& cf,           // input corner function of naive parameterization
    const VecXd& R,            //
    const vec<bool>& seam0,    //
          vec<bool>& seam1,    //
    std::set<HalfData>& h_data //
) {
    std::map<int, vec<AuxSgmt>> cuts; // face id & aux segment data

    for (const Thalf& th: tm.thalfs) {
        if (!th.cano) continue;
        const Tedge& te = th.edge();
        const double r = R[te.id];
        double sum = 0;
        int order = 0;
        for (const auto& s: te.segments()) {
            auto len = abs(s.diff());
            auto v0 = sum / r; sum += len;
            auto v1 = sum / r;
            auto aux = AuxSgmt{s, th, order, te.n_segments() - order, {v0, v1}, {1 - v1, 1 - v0}};
            order++;
            cuts[s.face.id].emplace_back(aux);
        }
    }

    vec h_aux(hm.nH, AuxHalf1{});
    vec f_aux(hm.nF, AuxFace{});

    vec<Row3d> vpos;
    for (auto p: hm.pos.rowwise()) { vpos.emplace_back(p); }

    for (Face f: hm.faces) {
        face_cutting(tm, f, cf, cuts[f.id], vpos, f_aux[f.id], h_aux);
    }

    int count = 0;
    for (const auto& hs: f_aux) { count += (int)hs.size(); }

    MatXi face_info(count, 3);
    MatXd vert_info(vpos.size(), 3);

    for (int i = 0; i < vpos.size(); i++) { vert_info.row(i) = vpos[i]; }

    int count2 = 0;
    for (const auto& fd1: f_aux) {
        for (auto& fd2: fd1) {
            face_info.row(count2) << fd2[0].i0, fd2[1].i0, fd2[2].i0;
            count2++;
        }
    }

    auto hm_cut = Hmesh(vert_info, face_info);
    h_data.clear();
    seam1 = std::vector(hm_cut.nE, false);

    struct EdgeKey {
        int tail;
        int head;
        bool operator==(const EdgeKey& o) const noexcept { return tail == o.tail && head == o.head; }
    };

    struct EdgeKeyHash {
        std::size_t operator()(const EdgeKey& k) const noexcept {
            return (static_cast<std::size_t>(k.tail) << 32) ^ static_cast<std::size_t>(k.head);
        }
    };

    std::unordered_map<EdgeKey, Half, EdgeKeyHash> half_by_verts;
    half_by_verts.reserve(hm_cut.nH * 2);

    for (Half h: hm_cut.halfs) {
        half_by_verts.insert({EdgeKey{h.tail().id, h.head().id}, h});
    }

    for (const auto& fd1: f_aux) {
    for (const auto& fd2: fd1) {
    for (int i = 0; i < 3; i++) {
        const auto& [i0, i1, original, data] = fd2[i];
        auto it = half_by_verts.find({i0, i1});

        if (original.has_value() && seam0[original.value().edge().id]) { seam1[it->second.edge().id] = true; }

        if (data.has_value()) {
            const auto& v = data.value();
            h_data.emplace(HalfData{it->second, v.v0, v.v1, v.thid, v.tqid, v.order, v.first, v.crash});
        }
    }}}

    return hm_cut;
}
}
#endif
