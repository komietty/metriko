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
    const std::set<HalfData>& data,
    const Hmesh& hm, // the cut mesh
    const Tmesh& tm, // the tmesh of original hmesh
    const VecXd& X
) {
    MatXd embedded_uv = MatXd::Zero(hm.nV, 2);
    auto [tq_bgn, tq_end] = rg::equal_range(data, tqid, {}, &HalfData::tqid);

    auto dir = complex(1, 0);
    auto sum = complex(0, 0);

    for (int i = 0; i < 4; i++) {
        for (int thid: tm.tquads[tqid].thids_by_side(i)) {
            auto [th_bgn, th_end] = rg::equal_range(tq_bgn, tq_end, thid, {}, &HalfData::thid);
            auto x = X[tm.thalfs[thid].edge().id];
            for (auto it = th_bgn; it != th_end; ++it) {
                auto val = x * it->v0;
                auto vid = it->half.tail().id;
                embedded_uv(vid, 0) = val * dir.real() + sum.real();
                embedded_uv(vid, 1) = val * dir.imag() + sum.imag();
            }
            sum += x * dir;
        }
        dir *= complex(0, 1);
    }

    // find all faces
    std::queue<int> queue;
    auto visit = std::vector(hm.nF, false);
    auto verts = std::vector(hm.nV, false);

    for (auto it = tq_bgn; it != tq_end; ++it){
        int fid = it->half.face().id;
        queue.emplace(fid);
        visit[fid] = true;
    }

    while (!queue.empty()) {
        Face f0 = hm.faces[queue.front()];
        queue.pop();
        for (Half h0: f0.adjHalfs()) {
            Face f1 = h0.twin().face();
            if (visit[f1.id] || rg::any_of(tq_bgn, tq_end, [&](auto& d) { return d.half == h0; })) continue;
            queue.emplace(f1.id);
            visit[f1.id] = true;
        }
    }

    for (Face f: hm.faces) {
        if (!visit[f.id]) continue;
        for (Half h: f.adjHalfs()) verts[h.tail().id] = true;
    }

    MatXi face_table = MatXi::Zero(rg::count(visit, true), hm.nF);
    MatXd vert_table = MatXd::Zero(rg::count(verts, true), hm.nV);
    std::unordered_map<int, int> idcs_table;

    int c0 = 0, c1 = 0;
    for (Vert v: hm.verts) { if (verts[v.id]) { vert_table(c0, v.id) = 1; idcs_table[v.id] = c0; c0++; } }
    for (Face f: hm.faces) { if (visit[f.id]) { face_table(c1, f.id) = 1; c1++; } }

    MatXd V  = vert_table * hm.pos;
    MatXi F  = face_table * hm.idx;
    MatXd UV = vert_table * embedded_uv;

    for (int i = 0; i < F.rows(); i++)
    for (int j = 0; j < 3; j++)
        F(i, j) = idcs_table[F(i, j)];

    // here tutte's parameterization
    Eigen::SparseLU<SprsD> lu;
    lu.compute(boundary_snap_laplacian(Hmesh(V, F)));
    MatXd uv = lu.solve(UV);

    MatXd uv_all = MatXd::Zero(hm.nC, 2);
    MatXd uv_vrt = vert_table.transpose() * uv;
    for (Face f: hm.faces) {
        if (!visit[f.id]) continue;
        for (Half h: f.adjHalfs())
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
    const std::unordered_set<Half, HalfHash>& boundary, // boundary of the tquad
    const vec<bool>& seam,                              // need to be altered for new cut hmesh
          vec<bool>& flag,                              // the flag to check a face is already marked
    MatXd& uv_all
) {
    std::queue<Half> queue;
    std::unordered_set<Edge, EdgeHash> visit;
    std::unordered_set<Half, HalfHash> nextH; // the halfedges to the other tquad
    queue.push(half_in);
    visit.emplace(half_in.edge());
    if (flag[half_in.face().id]) return nextH;

    while (!queue.empty()) {
        auto f  = queue.front().face();
        auto c0 = f.half().crnr();
        auto c1 = f.half().next().crnr();
        auto c2 = f.half().prev().crnr();
        uv_all.row(c0.id) = uv_in.row(c0.id);
        uv_all.row(c1.id) = uv_in.row(c1.id);
        uv_all.row(c2.id) = uv_in.row(c2.id);
        flag[f.id] = true;
        queue.pop();

        for (Half h: f.adjHalfs()) {
            // 1: if hit the seam, just stops
            if (seam[h.edge().id]) continue;
            // 2: if hit the visited edge, just stops
            if (visit.contains(h.edge())) continue;
            // 3: if hit boundary, puts it as a bridge to the next tquad
            if (boundary.contains(h)) { nextH.insert(h.twin()); continue; }
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
            mat1 = (mat1.rowwise() - uv1) * rot.transpose();
            mat1.rowwise() += uv0;
            return;
        }
    }
    throw std::runtime_error("no corresponding rotation found");
}

inline MatXd compute_tutte_parameterization(
    const Hmesh& hm,                // hmesh after tutte cutting
    const Tmesh& tm,                // tmesh original
    const vec<bool>& seam,          // seam adapted to tutte cutting
    const std::set<HalfData>& data, //
    const VecXd& X                  //
) {
    // compute uv per tquad first...
    vec<MatXd> uv_per_tquad;
    double t0 = omp_get_wtime();
    for (auto& tq: tm.tquads) {
        MatXd uv = embedding_tutte_for_tquad(tq.id, data,  hm, tm, X);
        uv_per_tquad.emplace_back(uv);
    }
    double t1 = omp_get_wtime();
    std::cout << "[time] prepare all tquad uv: " << (t1 - t0) << " s" << std::endl;

    MatXd uv = MatXd::Zero(hm.nC, 2);
    auto flag = std::vector(hm.nF, false);
    std::stack<Half> queue; // todo: when using queue, it does not work in count 519...

    std::unordered_map<Half, const HalfData*, HalfHash> data_by_half;
    data_by_half.reserve(data.size() * 2);
    for (auto& d : data) { data_by_half.emplace(d.half, &d); }

    double t2 = omp_get_wtime();
    { // 1: process the first tquad
        auto h = data.begin()->half;
        auto i = data.begin()->tqid;
        std::unordered_set<Half, HalfHash> b;
        for (auto&d : data) { if (d.tqid == i) b.insert(d.half); }
        auto o = sequential_mapping(uv_per_tquad[i], h, b, seam, flag, uv);
        for (auto nh: o) queue.emplace(nh);
    }

    // 2: other tquads
    while (rg::any_of(flag, [&](auto f) { return !f; })) {
        auto h = queue.top();
        queue.pop();

        if (flag[h.face().id]) continue;

        auto curr = data_by_half.at(h);
        auto prev = data_by_half.at(h.twin());
        MatXd& uv_curr = uv_per_tquad[curr->tqid];
        MatXd& uv_prev = uv_per_tquad[prev->tqid];

        apply_transition(h, uv_prev, uv_curr);

        std::unordered_set<Half, HalfHash> b;
        for (auto&d : data) { if (d.tqid == curr->tqid) b.insert(d.half); }
        auto o = sequential_mapping(uv_curr, h, b, seam, flag, uv);
        for (auto nh: o) queue.emplace(nh);
    }
    double t3 = omp_get_wtime();
    std::cout << "[time] assign them to locally injective uv: " << (t3 - t2) << " s" << std::endl;

    return uv;
}
}

#endif
