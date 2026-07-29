#ifndef METRIKO_QEX_REFINEMENT_H
#define METRIKO_QEX_REFINEMENT_H
#include <set>
#include <igl/AABB.h>
#include <igl/remove_duplicate_vertices.h>
#include "./common.h"

namespace metriko::qex {
    // simple laplacian smoothing of the extracted quad mesh, with vertices
    // constrained to remain on the input surface (projection after each step)
    inline void refinement_hmesh(
        const MatXd& pos,     // quad-soup vertex positions (4 per quad)
        const MatXi& idx,     // quad-soup indices (#quads x 4)
        const MatXd& V,       // input surface vertices (projection target)
        const MatXi& F,       // input surface triangles
        MatXd& pos_refined,
        MatXi& idx_refined,
        const int    iters  = 20,
        const double lambda = 0.5
    ) {
        /// 1: weld the per-quad duplicated corners into a connected quad mesh
        VecXi SVI, SVJ;
        const double eps = 1e-7 * (pos.colwise().maxCoeff() - pos.colwise().minCoeff()).norm();
        igl::remove_duplicate_vertices(pos, eps, pos_refined, SVI, SVJ);
        idx_refined.resize(idx.rows(), idx.cols());
        for (int i = 0; i < idx.rows(); ++i)
            for (int j = 0; j < idx.cols(); ++j)
                idx_refined(i, j) = SVJ(idx(i, j));

        /// 2: vertex adjacency along quad edges
        vec<std::set<int>> adj(pos_refined.rows());
        for (int i = 0; i < idx_refined.rows(); ++i)
            for (int j = 0; j < 4; ++j) {
                int a = idx_refined(i, j);
                int b = idx_refined(i, (j + 1) % 4);
                adj[a].insert(b);
                adj[b].insert(a);
            }

        /// 3: laplacian smoothing, projected back onto the input surface each step
        igl::AABB<MatXd, 3> tree;
        tree.init(V, F);
        for (int it = 0; it < iters; ++it) {
            MatXd next = pos_refined;
            for (int v = 0; v < pos_refined.rows(); ++v) {
                if (adj[v].empty()) continue;
                Row3d c = Row3d::Zero();
                for (int n: adj[v]) c += pos_refined.row(n);
                c /= (double)adj[v].size();
                next.row(v) = (1 - lambda) * pos_refined.row(v) + lambda * c;
            }
            VecXd sqrD; VecXi I; MatXd C;
            tree.squared_distance(V, F, next, sqrD, I, C);
            pos_refined = C;
        }
    }
}
#endif
