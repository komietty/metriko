//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_PARAMETERIZATION_H
#define METRIKO_PARAMETERIZATION_H
#include "../hmesh/hmesh.h"
#include "iterative_rounding/iter_rounding.h"

namespace metriko {
    class Data {
        const Hmesh &raw;
        const Hmesh &cut;
        const MatXd &ext;
        const VecXi &singlars;
        const VecXi &matching;
        const std::vector<bool> &seam;
        double gridscale;      // Global scaling of grid
        int nT;                // Number of transitions
        int nR;                // Number of vertices that are inside the raw mesh, plus nT
        int nS;                // Number of vertices that are singular but not on  boundaries
    public:
        int N;                 // Uncompressed parametric functions
        int n;                 // Independent parameteric functions
        SprsD vtrans2cut;      // Map between the whole mesh (vertex + translational jump) representation to the vertex-based representation on the cut mesh
        SprsD constraint;      // Linear constraints (resulting from non-singular nodes)
        SprsD uncompress;      // Global uncompression of n->N
        VecXi he2matching;     // Table of halfedge index to matching index
        VecXi he2transidx;     // Table of halfedge index to transiton index. Range from -nT to nT except 0. Sub 1 (or Add 1) when using as idx
        VecXi integerIdcs;     // Indices to be integer (n * transition)
        VecXi singularIdcs;    // Indices of singulars (used if singular has to be integer)
        VecXi fixedIdcs;       // Translation fixing indices
        VecXd fixedVals;       // Translation fixed values at fixed indices
        MatXi lreductor;       // Linear Reduction tying the n dofs to the full N
        VecXd fullx;           // Result: computed result
        MatXd nfn;             // Result: vertex-uv value on cut mesh
        MatXd cfn;             // Result: corner-uv value on raw mesh
        bool verbose;          // Output the integration log.
        bool seamless;         // Whether to do full translational seamless.
        bool roundSeams;       // Whether to round seams or round singularities
        bool localInjectivity; // Enforce local injectivity; might result in failure!

        explicit Data(
            const Hmesh &raw,
            const Hmesh &cut,
            const MatXd &ext,
            const VecXi &singlars,
            const VecXi &matching,
            const std::vector<bool> &seam,
            const int degree,
            const double ratio = 0.02
        ) : raw(raw),
            cut(cut),
            ext(ext),
            singlars(singlars),
            matching(matching),
            seam(seam),
            gridscale(ratio),
            seamless(false),
            roundSeams(true),
            verbose(false),
            localInjectivity(false)
        {
            N = degree;
            if (N % 2 == 0) set_sign_symmetry(N);
            else { lreductor = MatXi::Identity(N, N); n = N; }
            fixedVals = VecXd::Zero(n);
            compute_he2transidx();
            compute_he2matching();
            nR = raw.nV + nT;
            nS = std::ranges::count_if(raw.verts, [&](auto &v) { return is_inside_singular(v); });
        }

        void setup();
        bool integ();

    private:
        void compute_he2transidx();
        void compute_he2matching();
        bool is_inside_singular(Vert v) const { return !v.isBoundary() && singlars[v.id] != 0; }
        Half find_first_seam_he(Vert v) const { Half s = v.half(); for (auto h: v.adjHalfs(false)) if (seam[h.edge().id])     { s = h; break; } return s; }
        Half find_first_bndr_he(Vert v) const { Half b = v.half(); for (auto h: v.adjHalfs(false)) if (h.edge().isBoundary()) { b = h; break; } return b; }

        void set_sign_symmetry(const int N) {
            assert(N % 2 == 0);
            n = N / 2;
            lreductor.resize(N, N / 2);
            lreductor << MatXi::Identity(N / 2, N / 2),
                        -MatXi::Identity(N / 2, N / 2);
        }

        void set_tris_symmetry(const int N) {
            // the entire first N/3 lines are symmetric w.r.t. to the next two (N/3) packets,
            assert(N % 3 == 0);
            if (N % 2 == 0) {
                n = N / 3;
                lreductor.resize(N, N / 3);
                lreductor.block(0, 0, N / 2, N / 3) <<
                        MatXi::Identity(N / 3, N / 3),
                       -MatXi::Identity(N / 6, N / 6),
                        MatXi::Identity(N / 6, N / 6);
                lreductor.block(N / 2, 0, N / 2, N / 3) = -lreductor.block(0, 0, N / 2, N / 3);
            } else {
                n = 2 * N / 3;
                lreductor.resize(N, 2 * N / 3);
                lreductor <<
                    MatXi::Identity(2 * N / 3, 2 * N / 3),
                   -MatXi::Identity(N / 3, N / 3),
                   -MatXi::Identity(N / 3, N / 3);
            }
        }

        static double avg_norm(const MatXd &field, const int dim, const int N, const int nF) {
            double sum = 0.;
            for (int i = 0; i < nF; i++) {
                for (int j = 0; j < N; j++) {
                    sum += field.block(i, dim * j, 1, dim).norm();
                }
            }
            return sum / (N * nF);
        }

    };

    //----- utility function to cut mesh by seam -----//
    inline std::unique_ptr<Hmesh> compute_cut_mesh(
        const Hmesh& m,
        const std::vector<bool>& seam
    ) {
        MatXd cutV;
        MatXi cutF(m.nF, 3);
        std::vector<Row3d> cut2pos;

        for (Vert v: m.verts) {
            Half bgn = v.half();
            if (v.isBoundary()) { for (Half h : v.adjHalfs(false)) if (h.edge().isBoundary())  { bgn = h; break; } }
            else                { for (Half h : v.adjHalfs(false)) if (seam[h.edge().id])      { bgn = h; break; } }

            for (Half h : v.adjHalfs(bgn)) {
                if (h.isBoundary()) continue;
                if (seam[h.edge().id] || h.id == bgn.id) cut2pos.emplace_back(v.pos());
                for (int j = 0; j < 3; j++) {
                    const int iF = h.face().id;
                    if (m.idx(iF, j) == v.id) cutF(iF, j) = cut2pos.size() - 1;
                }
            }
        }

        cutV.resize(cut2pos.size(), 3);
        for (int iV = 0; iV < cut2pos.size(); iV++) cutV.row(iV) = cut2pos[iV];
        return std::make_unique<Hmesh>(cutV, cutF);
    }

    template <typename T1, typename T2>
    void assign_block(
        std::vector<Eigen::Triplet<T1>>& triplets,
        const Eigen::Matrix<T2, Eigen::Dynamic, Eigen::Dynamic>& mat,
        const int row_bgn,
        const int col_bgn
    ) {
        for (int i = 0; i < mat.rows(); i++) {
            for (int j = 0; j < mat.cols(); j++) {
                T2 val = mat(i, j);
                if (val != 0) // need to check it is possible to cast T1 to T2
                    triplets.emplace_back(row_bgn + i, col_bgn + j, val);
            }
        }
    }
}

#endif
