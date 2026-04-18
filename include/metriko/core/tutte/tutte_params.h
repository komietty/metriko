#ifndef METRIKO_TUTTE_PARAM_H
#define METRIKO_TUTTE_PARAM_H
#include <queue>
#include <stack>
#include "tutte_cutting.h"

namespace metriko::tutte {

inline Mat2d compute_rotation(int i) {
    Mat2d r0, r1, r2, r3;
    r0 <<  1,  0,  0,  1;
    r1 <<  0, -1,  1,  0;
    r2 << -1,  0,  0, -1;
    r3 <<  0,  1, -1,  0;
    auto r = vec{r2, r1, r0, r3}; // need fix
    return r[i];
}

inline SprsD embedding_tutte_for_tquad(
    const int tqid,
    const vec<HalfData>& data,
    const Hmesh& hm, // the cut mesh
    const Emesh& tm  // the tmesh of original hmesh
) {
    vec<TripD> T;

    SprsD emb_uv(hm.nV, 2);
    auto tq_rg = rg::equal_range(data, tqid, {}, &HalfData::tqid);

    auto dir = complex(1, 0);
    auto sum = complex(0, 0);

    for (int i = 0; i < 4; i++) {
        for (int thid: tm.equads[tqid].ehids_by_side(i)) {
            //auto x = X[tm.ehalfs[thid].edge().id];
            auto x = tm.ehalfs[thid].x;
            for (auto& it: rg::equal_range(tq_rg, thid, {}, &HalfData::thid)) {
                auto val = x * it.v0;
                auto vid = it.half.tail().id;
                T.emplace_back(vid, 0, val * dir.real() + sum.real());
                T.emplace_back(vid, 1, val * dir.imag() + sum.imag());
            }
            sum += x * dir;
        }
        dir *= complex(0, 1);
    }
    emb_uv.setFromTriplets(T.begin(), T.end());

    // find all faces
    std::queue<int> queue;
    auto visit = std::vector(hm.nF, false);
    auto verts = std::vector(hm.nV, false);

    for (auto& it: tq_rg){
        int fid = it.half.face().id;
        queue.emplace(fid);
        visit[fid] = true;
    }

    while (!queue.empty()) {
        Face f0 = hm.faces[queue.front()];
        queue.pop();
        for (Half h0: f0.adjHalfs()) {
            Face f1 = h0.twin().face();
            if (visit[f1.id] || rg::contains(tq_rg, h0, &HalfData::half)) continue;
            queue.emplace(f1.id);
            visit[f1.id] = true;
        }
    }

    for (Face f: hm.faces) {
        if (!visit[f.id]) continue;
        for (Half h: f.adjHalfs()) verts[h.tail().id] = true;
    }

    int nF_sub = (int)rg::count(visit, true);
    int nV_sub = (int)rg::count(verts, true);

    SprsD f_table(nF_sub, hm.nF);
    SprsD v_table(nV_sub, hm.nV);
    std::vector<TripD> Tv, Tf;
    std::unordered_map<int, int> idcs_table;
    Tv.reserve(nV_sub);
    Tf.reserve(nF_sub);

    int c0 = 0, c1 = 0;
    for (Vert v : hm.verts) { if (verts[v.id]) { Tv.emplace_back(c0, v.id, 1.); idcs_table[v.id] = c0; ++c0; }}
    for (Face f : hm.faces) { if (visit[f.id]) { Tf.emplace_back(c1, f.id, 1.); ++c1; }}
    v_table.setFromTriplets(Tv.begin(), Tv.end());
    f_table.setFromTriplets(Tf.begin(), Tf.end());

    MatXd V  = v_table * hm.pos;
    MatXd F  = f_table * hm.idx.cast<double>();
    SprsD UV = v_table * emb_uv;

    for (int i = 0; i < F.rows(); i++)
    for (int j = 0; j < 3; j++)
        F(i, j) = idcs_table[F(i, j)];

    Eigen::SparseLU<SprsD> lu;
    lu.compute(boundary_snap_laplacian(Hmesh(V, F.cast<int>())));
    SprsD uv = lu.solve(UV);

    SprsD uv_all(hm.nC, 2);
    T.clear();
    SprsD uv_vrt = v_table.transpose() * uv;
    for (Face f: hm.faces) {
        if (!visit[f.id]) continue;
        for (Half h: f.adjHalfs()) {
            Row2d r = uv_vrt.row(h.next().head().id);
            T.emplace_back(h.crnr().id, 0, r.x());
            T.emplace_back(h.crnr().id, 1, r.y());
        }
    }
    uv_all.setFromTriplets(T.begin(), T.end());
    return uv_all;
}

// take the tutte result as the input, embed it until seam intersection.
// computes halfedges to search with in the next loop at the same time.
inline vec<int> sequential_mapping(
    const Hmesh& hm,
    const SprsD& uv_in,
    const Half half_in,
    const vec<bool>& boundary, // flag if a halfedge is the boundary of the tquad
    const vec<bool>& seam,     // need to be altered for new cut hmesh
          vec<bool>& flag,     // the flag to check a face is already marked
    MatXd& uv_all
) {
    vec<int> nextH; // the halfedges to the other tquad
    vec visit(hm.nE, false);
    std::queue<Half> Q;
    Q.push(half_in);
    visit[half_in.edge().id] = true;
    if (flag[half_in.face().id]) return nextH;

    while (!Q.empty()) {
        auto f  = Q.front().face();
        auto c0 = f.half().crnr();
        auto c1 = f.half().next().crnr();
        auto c2 = f.half().prev().crnr();
        uv_all.row(c0.id) = uv_in.row(c0.id);
        uv_all.row(c1.id) = uv_in.row(c1.id);
        uv_all.row(c2.id) = uv_in.row(c2.id);
        flag[f.id] = true;
        Q.pop();

        for (Half h: f.adjHalfs()) {
            Edge e = h.edge();
            // 1: if hit the seam, just stops
            if (seam[e.id]) continue;
            // 2: if hit the visited edge, just stops
            if (visit[e.id]) continue;
            // 3: if hit boundary, puts it as a bridge to the next tquad
            if (boundary[h.id]) { nextH.emplace_back(h.twin().id); continue; }
            // 4: inside tquad. add it to the queue
            visit[e.id] = true;
            Q.push(h.twin());
        }
    }
    return nextH;
}

// try to multiply rotation until halfedge coner values corresponds
// need to consider: is there any possibility of flip?
inline bool apply_transition(
    const bool flag,
    const Half h,    // the halfedge of unfixed side
    const SprsD& m0, // the fixed uv information
          SprsD& m1  // the unfixed adjacent uv information
) {

    Half h0 = h.twin();
    Half h1 = h;
    Vec2d uv0  = m0.row(h0.next().crnr().id).transpose();
    Vec2d uv0a = m0.row(h0.prev().crnr().id).transpose();
    Vec2d uv1  = m1.row(h1.prev().crnr().id).transpose();
    Vec2d uv1a = m1.row(h1.next().crnr().id).transpose();

    if (flag) {
        std::cout << "h0 tail: " << h0.tail().id << std::endl;
        std::cout << "h0 head: " << h0.head().id << std::endl;
        std::cout << "h1 tail: " << h1.tail().id << std::endl;
        std::cout << "h1 head: " << h1.head().id << std::endl;
        std::cout << "cid: " << h0.prev().crnr().id << std::endl;
        //uv0a = Vec2d(-8, -8.21001);
    }


    if (flag) {
        std::cout << "uv0 : " <<  uv0.transpose()  << std::endl;
        std::cout << "uv0a: " <<  uv0a.transpose() << std::endl;
        std::cout << "uv1 : " <<  uv1.transpose()  << std::endl;
        std::cout << "uv1a: " <<  uv1a.transpose() << std::endl;
    }

    for (int i = 0; i < 4; i++) {
        Mat2d rot = compute_rotation(i);
        Vec2d v1 = rot * (uv1a - uv1);
        Vec2d v2 = uv0a - uv0;

        if ((v1 - v2).norm() < 1e-9) {

            if (flag) {
                std::cout << "rot : " <<  rot  << std::endl;
                std::cout << "v1: " <<  v1.transpose() << std::endl;
                std::cout << "v2: " <<  v2.transpose()  << std::endl;
            }
            for (SprsD::InnerIterator it(m1, 0); it; ++it) {
                int ir = it.row();
                Vec2d p(it.value(), m1.coeff(ir, 1));
                p = rot * (p - uv1) + uv0;
                m1.coeffRef(ir, 0) = p.x();
                m1.coeffRef(ir, 1) = p.y();
            }
            return true;
        }
    }

    std::cout << "failed to map tutte params: " << std::endl;
    return false;
}

struct HalfHash { std::size_t operator()(const Half& h) const noexcept { return std::hash<int>{}(h.id); } };

inline bool compute_tutte_parameterization(
    const Hmesh& hm,           // hmesh after tutte cutting
    const Emesh& tm,           // tmesh original
    const vec<bool>& seam,     // seam adapted to tutte cutting
    const vec<HalfData>& data, //
    MatXd& uv
) {
    bool success = true;
    // compute uv per tquad first...
    vec<SprsD> uv_tq;
    uv_tq.resize(tm.equads.size());

    //#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < tm.equads.size(); i++) {
        if (tm.equads[i].id != -1)
            uv_tq[i] = embedding_tutte_for_tquad(i, data, hm, tm);
    }

    uv.resize(hm.nC, 2);
    uv.setZero();

    auto flag = vec(hm.nF, false);
    std::stack<int> stack;

    std::unordered_map<Half, const HalfData*, HalfHash> data_by_half;
    data_by_half.reserve(data.size());
    for (auto& d : data) { data_by_half.emplace(d.half, &d); }

    { // 1: process the first tquad
        auto h = data.begin()->half;
        auto i = data.begin()->tqid;
        vec b(hm.nH, false);
        for (auto& d: data) { if (d.tqid == i) b[d.half.id] = true; }
        for (auto nh: sequential_mapping(hm, uv_tq[i], h, b, seam, flag, uv)) stack.emplace(nh);
    }

    int count = 0;
    // 2: other tquads
    while (rg::any_of(flag, [&](auto f) { return !f; })) {
        count++;
        //std::cout << "[count] " << count << std::endl;
        auto h = hm.halfs[stack.top()];

        stack.pop();
        if (flag[h.face().id]) continue;

        auto it = data_by_half.find(h);
        if (it == data_by_half.end()) {
            std::cerr << "[Error] compute_tutte_parameterization: Halfedge " << h.id << " not found in data_by_half map." << std::endl;
            return false;
        }

        auto curr = it->second;
        SprsD& uv_curr = uv_tq[curr->tqid];
        success &= apply_transition(false, h, uv.sparseView(), uv_curr);

        vec b(hm.nH, false);
        for (auto& d: data) { if (d.tqid == curr->tqid) b[d.half.id] = true; }
        for (auto nh: sequential_mapping(hm, uv_curr, h, b, seam, flag, uv)) stack.emplace(nh);
    }
    return success;
}
}

#endif
