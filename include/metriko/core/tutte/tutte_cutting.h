#ifndef METRIKO_TUTTE_CUTTING_H
#define METRIKO_TUTTE_CUTTING_H
#include <set>

#include "emesh.h"
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/tmesh/tmesh.h"

namespace metriko::tutte {

namespace impl_eb {
constexpr double EPS = 1e-10;

using Mat3x2d = Eigen::Matrix<double, 3, 2>;
using Mat2x3d = Eigen::Matrix<double, 2, 3>;

inline Mat2x3d axis_align_proj(const Row3d& n) {
    Row3d a = n.cwiseAbs();
    double max;
    Mat3x2d P;

    if (a.z() > a.x() && a.z() > a.y()) { P.col(0) << 1., 0., 0.; P.col(1) << 0., 1., 0.; max = n.z(); }
    else if (a.y() > a.x())             { P.col(0) << 0., 0., 1.; P.col(1) << 1., 0., 0.; max = n.y(); }
    else                                { P.col(0) << 0., 1., 0.; P.col(1) << 0., 0., 1.; max = n.x(); }

    if (max < 0) P.col(0) *= -1;
    return P.transpose();
}
}

struct AuxSgmt {
    const Tsgmt sg; //
    const Thalf th; //
    int ord0;       // order cano
    int ord1;       // order non cano
    Row2d v0;       // val cano
    Row2d v1;       // val non cano
    bool isBgn;     //
    bool isEnd;     //
};

struct AuxHalf2 {
    int i0;
    int i1;
    std::optional<Half> original;
    std::optional<HalfData> data;
};

using AuxHalf = std::set<std::pair<double, int>>; // the list of index and ratio in the Half
using AuxFace = vec<std::array<AuxHalf2, 3>>;     // per original face, this contains halfedge data of a divided triangle

// Beware epsilon validity must be solved beforehand
inline void face_cutting(
    const Tmesh& tm,         //
    const Hmesh& hm,         //
    const Face& f,           // face to be cut
    const VecXc& cf,         // corner function of original mesh
    const vec<AuxSgmt>& sgs, // segments inside the face
    vec<Row3d>& pos,         // vertex position (altered)
    AuxFace& f_auxs,         //
    vec<AuxHalf>& h_auxs     //
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
        for (auto [r, i]: h_auxs[h0.twin().id]) { if (r > 0 && r < 1) idcs.emplace_back(i); }
    }

    for (const auto& [sg, th, o0, o1, v0, v1, f0, b1]: sgs) {
        std::pair<const Tvert*, int> mvs[2] = { std::pair(&sg.tvFr, -1), std::pair(&sg.tvTo, -1) };

        /// 1: Assign inside halfs
        for (auto& [mv, id]: mvs) {
            Row3d p;
            std::optional<std::pair<Half, double>> val;
            if (mv->hid != -1) {
                Half h = hm.halfs[mv->hid];
                auto r = mv->rt;
                p = h.tail().pos() * r + h.head().pos() * (1 - r);
                val = std::pair(h, r);
            }
            else { p = conversion_2d_3d(sg.face, cf, mv->uv); }

            auto find = rg::find_if(idcs, [&](int i) { return (pos[i] - p).norm() < EPS; });
            if (find != idcs.end()) { id = *find; }
            else {
                pos.emplace_back(p);
                int l = (int)pos.size() - 1;
                idcs.emplace_back(l);
                id = l;
                if (val.has_value()) {
                    auto [h0, r] = val.value();
                    auto h1 = h0.twin();
                    h_auxs[h0.id].emplace(r, id);
                    h_auxs[h1.id].emplace(1 - r, id);
                }
            }
        }

        int i0 = mvs[0].second;
        int i1 = mvs[1].second;
        if (i0 == i1) continue;
        halfs.emplace_back(AuxHalf2{i0, i1, std::nullopt, HalfData{Half(), v0.x(), v0.y(), th.id  , tm.th2quad(th.id)  , -1, o0, f0, b1 }});
        halfs.emplace_back(AuxHalf2{i1, i0, std::nullopt, HalfData{Half(), v1.x(), v1.y(), th.twid, tm.th2quad(th.twid), -1, o1, false, false }});
    }

    /// 2: Assign edge halfs
    for (Half h: f.adjHalfs()) {
        const auto& idcs_ = h_auxs[h.id];
        if (idcs_.empty()) halfs.emplace_back(h.tail().id, h.head().id, h);
        else {
            halfs.emplace_back(idcs_.begin()->second, h.head().id, h);

            auto prev = idcs_.begin();
            auto curr = std::next(idcs_.begin());
            for (; curr != idcs_.end(); ++curr, ++prev) { halfs.emplace_back(curr->second, prev->second, h); }
            halfs.emplace_back(h.tail().id, idcs_.rbegin()->second, h);
        }
    }

    /// 3: Create halfedge polyline
    while (!halfs.empty()) {
        auto h = halfs.back();
        halfs.pop_back();

        vec poly = {h};
        int sta = h.i0;
        int cur = h.i1;

        while (true) {
            int prev_v = poly.back().i0;
            int curr_v = poly.back().i1;
            const Row3d& p_prev = pos[prev_v];
            const Row3d& p_curr = pos[curr_v];
            const Row3d  d_prev = (p_curr - p_prev).normalized();

            auto best_it = halfs.end();
            double score = -2;

            for (auto it = halfs.begin(); it != halfs.end(); ++it) {
                if (it->i0 != curr_v) continue;
                if (it->i1 == prev_v) continue;
                int cand_to = it->i1;
                Row3d d_cand = (pos[cand_to] - p_curr).normalized();

                // 向き判定：d_prev から d_cand への回転が CCW かどうか
                Row3d c = d_prev.cross(d_cand);  // 回転軸方向
                if (c.norm() > EPS && f.normal().dot(c) <= 0) { continue; }

                double s = (-d_prev).dot(d_cand);
                if (s > score) { score = s; best_it = it; }
            }

            // CCW 候補が無ければチェーン終了
            if (best_it == halfs.end()) break;

            int next = best_it->i1;
            poly.push_back(*best_it);

            cur = next;
            halfs.erase(best_it);
            if (cur == sta) break;
        }

        bool flag = false;
        for (auto& p: poly) if (p.i0 == 0 && p.i1 == 0) flag = true;


        if (flag) {
            std::cout << "poly: ";
            for (auto& p: poly) std::cout << p.i0 << ",  ";
            std::cout << std::endl;
        }

        /// 4: Calc simple ear clipping for convex polygon.
        ///    Find the starting vertex not to generate 0 size area
        int i = 0;
        int n = (int)poly.size();

        for (; i < n; ++i) {
            bool valid = rg::all_of(vw::iota(1, n - 1), [&](int j) {
                Row3d& p0 = pos[poly[i].i0];
                Row3d& p1 = pos[poly[(i + j    ) % n].i0];
                Row3d& p2 = pos[poly[(i + j + 1) % n].i0];
                return std::abs(f.normal().dot((p1 - p0).cross(p2 - p0))) >= EPS;
            });
            if (valid) break;
        }

        for (int j = 1; j < n - 1; ++j)
            f_auxs.emplace_back(std::array{
                poly[i],
                poly[(i + j    ) % n],
                poly[(i + j + 1) % n]
            });
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
inline std::unique_ptr<Hmesh> compute_embedding_cut_hmesh(
    const Hmesh& hm,           // input hmesh
    const Tmesh& tm,           // input tmesh
    const VecXc& cf,           // input corner function of naive parameterization
    const vec<bool>& seam0,    //
          vec<bool>& seam1,    //
    std::set<HalfData>& h_data //
) {
    h_data.clear();
    std::map<int, vec<AuxSgmt>> cuts; // face id & aux segment data

    for (const auto& th: tm.thalfs) {
        if (!th.cano) continue;
        const auto& te = th.edge();
        const auto r = te.len;
        const auto n = (int)te.segs.size() - 1;
        auto v = 0.;
        auto i = 0;
        for (const auto& s: te.segs) {
            auto l = abs(s.tvTo.uv - s.tvFr.uv);
            auto j = n - i;
            auto v0 = v / r; v += l;
            auto v1 = v / r;
            bool f0 = i == 0 && te.isBgn;
            bool b1 = i == n && te.isEnd;
            cuts[s.face.id].emplace_back(AuxSgmt{s, th, i, j, {v0, v1}, {1 - v1, 1 - v0}, f0, b1});
            i++;
        }
    }

    vec h_aux(hm.nH, AuxHalf{});
    vec f_aux(hm.nF, AuxFace{});

    vec<Row3d> pos;
    for (auto p: hm.pos.rowwise()) pos.emplace_back(p);
    for (auto f: hm.faces) face_cutting(tm, hm, f, cf, cuts[f.id], pos, f_aux[f.id], h_aux);

    MatXd V(pos.size(), 3);
    MatXi F(rg::distance(f_aux | vw::join), 3);

    for (int i = 0; i < pos.size(); i++) V.row(i) = pos[i];

    int c = 0;
    for (const auto& d : f_aux | vw::join)  F.row(c++) << d[0].i0, d[1].i0, d[2].i0;

    auto hmC = std::make_unique<Hmesh>(V, F);
    seam1 = vec(hmC->nE, false);

    auto find_half = [&](int v0, int v1) -> Half {
        for (Half h: hmC->verts[v0].adjHalfs()) if (h.head().id == v1) return h;
        throw std::runtime_error("half not found");
    };

    std::vector<glm::vec3> ns;
    std::vector<std::array<size_t, 2>> es;
    std::vector<double> val1;
    size_t count = 0;

    int counter = 0;
    for (const auto& a: f_aux | vw::join) {
    for (const auto& [i0, i1, ori, d]: a) {
        if (i0 == i1) {
            std::cout << "===========" << std::endl;
            for (int i = 0; i < 3; i++) {
                Row3d p1 = hmC->verts[a[i].i0].pos();
                Row3d p2 = hmC->verts[a[i].i1].pos();
                std::cout << "i0: " << a[i].i0  << ", i1: " << a[i].i1 << ", norm: " << (p1 - p2).norm() << std::endl;
                ns.emplace_back(p1.x(), p1.y(), p1.z());
                ns.emplace_back(p2.x(), p2.y(), p2.z());
                es.emplace_back(std::array{count, count + 1});
                val1.emplace_back(count / 2);
                count += 2;
            }
        }
        counter++;
    }}

    auto cv = polyscope::registerCurveNetwork("i0 == i1 tris", ns, es);
    cv->addEdgeScalarQuantity("order", val1);
    cv->resetTransform();
    cv->setRadius(0.002);

    //for (const auto& [i0, i1, o, d]: f_aux | vw::join | vw::join) {
    //    Half h = find_half(i0, i1);
    //    if (o.has_value() && seam0[o->edge().id]) seam1[h.edge().id] = true;
    //    if (d.has_value()) {
    //        HalfData hd = d.value();
    //        hd.half = h;
    //        h_data.insert(hd);
    //    }
    //}




    //auto temp = vec(h_data.begin(), h_data.end());
    //for (auto& hd: temp)
    //    for (int i = 0; i < temp.size(); i++)
    //        if (hd.half.twin() == temp[i].half) hd.twin = i;

    return hmC;
}
}
#endif
