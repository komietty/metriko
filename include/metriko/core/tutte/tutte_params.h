#ifndef METRIKO_TUTTE_PARAM_H
#define METRIKO_TUTTE_PARAM_H
#include <queue>
#include <stack>
#include "tutte_cutting.h"

namespace metriko::tutte {

inline Mat2d compute_rotation(int i) {
    Mat2d r0, r1, r2, r3;
    r0 << 1, 0, 0, 1;
    r1 << 0, -1, 1, 0;
    r2 << -1, 0, 0, -1;
    r3 << 0, 1, -1, 0;
    auto r = std::vector{r2, r1, r0, r3}; // need fix
    return r[i];
}

inline MatXd embedding_tutte_for_tquad(
    const int tqid,
    const vec<HalfData>& data,
    const Hmesh& hm, // the cut mesh
    const Tmesh& tm, // the tmesh of original hmesh
    const VecXd& X
) {
    const Tquad& tq = tm.tquads[tqid];
    auto data0 = vw::filter(data, [&tqid](const HalfData& d) { return d.tqid == tqid; });
    MatXd embedded_uv = MatXd::Zero(hm.nV, 2);

    auto dir = complex(1, 0);
    auto sum = complex(0, 0);
    for (int i = 0; i < 4; i++) {
        for (int thid: tq.thids_by_side(i)) {
            // 1: filter from the original container by thid and sort it by segment order
            auto data1 = vw::filter(data0, [&thid](const HalfData& d) { return d.thid == thid; });
            vec data2(rg::begin(data1), rg::end(data1));
            rg::sort(data2, {}, &HalfData::order);

            // 2: assign values
            double x = X[tm.thalfs[thid].edge().id];
            for (const HalfData& d: data2) {
                auto val = x * d.v0;
                auto vid = d.half.tail().id;
                embedded_uv(vid, 0) = val * dir.real() + sum.real();
                embedded_uv(vid, 1) = val * dir.imag() + sum.imag();
            }
            sum += x * dir;
        }
        dir *= complex(0, 1);
    }

    // find all faces
    std::queue<int> queue;
    std::set<int> visit; // unordered_set seems too wild...

    for (const HalfData& d: data0) {
        queue.emplace(d.half.face().id);
        visit.emplace(d.half.face().id);
    }

    while (!queue.empty()) {
        Face f0 = hm.faces[queue.front()];
        queue.pop();
        for (Half h0: f0.adjHalfs()) {
            Face f1 = h0.twin().face();
            if (visit.contains(f1.id) || rg::any_of(data0, [&](auto& d) { return d.half == h0; })) continue;
            queue.emplace(f1.id);
            visit.emplace(f1.id);
        }
    }

    std::set<int> verts;
    for (int fid: visit) {
        for (Half h: hm.faces[fid].adjHalfs()) verts.emplace(h.tail().id);
    }

    MatXi face_table = MatXi::Zero((int)visit.size(), hm.nF);
    MatXd vert_table = MatXd::Zero((int)verts.size(), hm.nV);
    std::unordered_map<int, int> idcs_table;

    int c0 = 0, c1 = 0;
    for (int vid: verts) { vert_table(c0, vid) = 1; idcs_table[vid] = c0; c0++; }
    for (int fid: visit) { face_table(c1, fid) = 1; c1++; }

    MatXd V  = vert_table * hm.pos;
    MatXi F  = face_table * hm.idx;
    MatXd UV = vert_table * embedded_uv;

    for (int i = 0; i < F.rows(); i++)
    for (int j = 0; j < 3; j++)
        F(i, j) = idcs_table[F(i, j)];

    // here tutte's parameterization
    auto m   = std::make_unique<Hmesh>(V, F);
    SprsD BL = boundary_snap_laplacian(*m);
    MatXd uv(m->nV, 2);
    { Eigen::SparseLU<SprsD> lu; lu.compute(BL); VecXd res = lu.solve(UV.col(0)); uv.col(0) = res; }
    { Eigen::SparseLU<SprsD> lu; lu.compute(BL); VecXd res = lu.solve(UV.col(1)); uv.col(1) = res; }

    MatXd uv_all = MatXd::Zero(hm.nC, 2);
    MatXd uv_vrt = vert_table.transpose() * uv;
    for (int fid: visit) {
        for (Half h: hm.faces[fid].adjHalfs())
            uv_all.row(h.crnr().id) = uv_vrt.row(h.next().head().id);
    }

    return uv_all;
}


struct EdgeHash { std::size_t operator()(const Edge& e) const noexcept { return std::hash<int>{}(e.id); } };
struct HalfHash { std::size_t operator()(const Half& h) const noexcept { return std::hash<int>{}(h.id); } };

// take the tutte result as the input, embed it until seam intersection.
// computes halfedges to search with in the next loop at the same time.
inline std::unordered_set<Half, HalfHash> sequential_mapping(
    const MatXd& uv_in,
    const Half half_in,
    const vec<HalfData>& boundary, // boundary of the tquad
    const vec<bool>& seam,         // need to be altered for new cut hmesh
          vec<bool>& flag,         // the flag to check a corner is already checked or not
    MatXd& uv_all
) {
    std::queue<Half> queue;
    std::unordered_set<Edge, EdgeHash> visit;
    std::unordered_set<Half, HalfHash> nextH; // the halfedges to the other tquad
    queue.push(half_in);
    visit.emplace(half_in.edge());
    if (flag[half_in.crnr().id]) return nextH;

    while (!queue.empty()) {
        auto f  = queue.front().face();
        auto c0 = f.half().crnr();
        auto c1 = f.half().next().crnr();
        auto c2 = f.half().prev().crnr();
        uv_all.row(c0.id) = uv_in.row(c0.id);
        uv_all.row(c1.id) = uv_in.row(c1.id);
        uv_all.row(c2.id) = uv_in.row(c2.id);
        flag[c0.id] = true;
        flag[c1.id] = true;
        flag[c2.id] = true;
        queue.pop();

        for (Half h: f.adjHalfs()) {
            // 1: if hit the seam, just stops
            if (seam[h.edge().id]) continue;
            // 2: if hit the visited edge, just stops
            if (visit.contains(h.edge())) continue;
            // 3: if hit boundary, puts it as a bridge to the next tquad
            if (rg::any_of(boundary, [&h](const HalfData& d){return d.half == h; })) { nextH.insert(h.twin()); continue; }
            // 4: inside tquad. add it to the queue
            visit.emplace(h.edge());
            queue.push(h.twin());
        }
    }
    return nextH;
}

// try to multiply rotation until halfedge coner values corresponds
// need to consider: is there any possibility of flip?
inline void apply_transition(
    const Half h,      // the halfedge of unfixed side
    const MatXd& mat0, // the fixed uv information
          MatXd& mat1  // the unfixed adjacent uv information
) {
    Half h0 = h.twin();
    Half h1 = h;
    Row2d uv0  = mat0.row(h0.next().crnr().id);
    Row2d uv1  = mat1.row(h1.prev().crnr().id);
    Row2d uv0a = mat0.row(h0.prev().crnr().id);
    Row2d uv1a = mat1.row(h1.next().crnr().id);

    for (int i = 0; i < 4; i++) {
        Mat2d rot = compute_rotation(i);
        Row2d res = rot * (uv1a - uv1).transpose() + uv0.transpose();
        if ((res - uv0a).norm() < 1e-6) {
            for (int j = 0; j < mat1.rows(); j++) {
                mat1.row(j) = rot * (mat1.row(j) - uv1).transpose() + uv0.transpose();
            }
            return;
        }
    }
    throw std::runtime_error("no corresponding rotation found");
}

inline MatXd compute_tutte_parameterization(
    const Hmesh& hm,           // hmesh after tutte cutting
    const Tmesh& tm,           // tmesh original
    const vec<bool>& seam,     // seam adapted to tutte cutting
    const vec<HalfData>& data, //
    const VecXd& X             //
) {
    // compute uv per tquad first...
    std::vector<MatXd> uv_per_tquad;
    for (auto& tq: tm.tquads) {
        MatXd uv = embedding_tutte_for_tquad(tq.id, data,  hm, tm, X);
        uv_per_tquad.emplace_back(uv);
    }

    MatXd uv = MatXd::Zero(hm.nC, 2);
    auto c_flag = std::vector(hm.nC, false);
    std::stack<Half> queue; // todo: when using queue, it does not work in count 519...

    { // 1: process the first tquad
        auto h = data.front().half;
        auto i = data.front().tqid;
        auto b = data | vw::filter([&i](auto& d) { return d.tqid == i; })
                      | rg::to<std::vector>();
        auto o = sequential_mapping(uv_per_tquad[i], h, b, seam, c_flag, uv);
        for (auto nh: o) queue.push(nh);
    }

    // 2: other tquads
    while (rg::any_of(c_flag, [&](auto f) { return !f; })) {
        auto h = queue.top();
        queue.pop();

        auto curr = rg::find_if(data, [&h](auto& d) { return d.half == h; });
        auto prev = rg::find_if(data, [&h](auto& d) { return d.half == h.twin(); });
        MatXd& uv_curr = uv_per_tquad[curr->tqid];
        MatXd& uv_prev = uv_per_tquad[prev->tqid];

        apply_transition(h, uv_prev, uv_curr);
        auto b = data | vw::filter([&](auto& d) { return d.tqid == curr->tqid; })
                      | rg::to<std::vector>();
        auto o = sequential_mapping(uv_curr, h, b, seam, c_flag, uv);
        for (auto nh: o) queue.push(nh);
    }

    return uv;
}
}

#endif
