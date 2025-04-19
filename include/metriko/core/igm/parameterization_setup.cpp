#include "parameterization.h"

namespace metriko {
    void Data::compute_he2matching() {
        he2matching.resize(raw.nH);
        for (Half h: raw.halfs) {
            int m = (h.isCanonical() ? -1 : 1) * matching(h.edge().id); // better not inversed
            he2matching[h.id] = m < 0 ? (N + m % N) % N : m % N;
        }
    }

    void Data::compute_he2transidx() {
        he2transidx.resize(raw.nH);
        he2transidx.setConstant(32767);
        std::vector valence(raw.nV, 0.);
        std::vector claimed(raw.nE, false);

        for (Edge e: raw.edges)
            if (seam[e.id]) {
                valence[e.vert0().id]++;
                valence[e.vert1().id]++;
            }

        int tid = 1;

        for (Vert v: raw.verts) {
            if ((valence[v.id] == 2 && !is_inside_singular(v)) || valence[v.id] == 0) continue;
            for (Half cH: v.adjHalfs(find_first_bndr_he(v))) {
                Edge cE = cH.edge();
                if (cE.isBoundary() || !seam[cE.id] || claimed[cE.id]) continue;
                he2transidx[cH.id] = tid;
                he2transidx[cH.twin().id] = -tid;
                claimed[cE.id] = true;

                // traces seam until reaching to singular or boundary vertex
                Vert jv = cH.head();
                while (valence[jv.id] == 2 && !is_inside_singular(jv) && !jv.isBoundary()) {
                    Half nh;
                    for (Half h: jv.adjHalfs())
                        if (seam[h.edge().id] && !claimed[h.edge().id]) { nh = h; break; }

                    he2transidx[nh.id] = tid;
                    he2transidx[nh.twin().id] = -tid;
                    claimed[nh.edge().id] = true;
                    jv = nh.head();
                }
                tid++;
            }
        }
        nT = tid - 1;
    }

    void Data::setup() {
        namespace rg = std::ranges;
        namespace vw = std::views;

        // here we compute a permutation matrix
        std::vector<MatXi> constParmMats(N);
        MatXi unitPermMat = MatXi::Zero(N, N);
        for (int i = 0; i < N; i++) unitPermMat((i + 1) % N, i) = 1;

        // generate all the members of the permutation group
        constParmMats[0] = MatXi::Identity(N, N);
        for (int i = 1; i < N; i++) constParmMats[i] = unitPermMat * constParmMats[i - 1];


        std::vector<TripD> vT, cT;
        // forming the constraints and the singularity positions
        int currConstraint = 0;
        // this loop set up the transtions (vector field matching) across the cuts
        for (Vert v: raw.verts) {
            // 1: The initial corner gets the identity without any transition
            std::vector<MatXi> permMats;
            std::vector<int> permIdcs;
            permMats.emplace_back(MatXi::Identity(N, N));
            permIdcs.emplace_back(v.id);
            int iVcut = -1;

            // remenber uv = Rot * val + Transition
            for (Half h: v.adjHalfs(v.isBoundary() ? find_first_bndr_he(v) : find_first_seam_he(v))) {
                if (h.isBoundary()) break; /// last boundary
                Face f = h.face();
                int jVcut = -1;
                for (int j = 0; j < 3; j++) if (raw.idx(f.id, j) == v.id) jVcut = cut.idx(f.id, j);

                if (jVcut != iVcut) {
                    iVcut = jVcut;
                    for (int i = 0; i < permIdcs.size(); i++)
                        assign_block(vT, permMats[i], N * iVcut, N * permIdcs[i]);
                }

                Half hn = h.prev().twin();
                if (!hn.isBoundary() && seam[hn.edge().id]) {
                    auto p = constParmMats[he2matching[hn.id]];
                    auto t = he2transidx[hn.id];
                    if (t > 0) {
                        for (auto &m: permMats) m = p * m;
                        permMats.emplace_back(MatXi::Identity(N, N));
                        permIdcs.push_back(raw.nV + t - 1);
                    } else {
                        permMats.emplace_back(-MatXi::Identity(N, N));
                        permIdcs.push_back(raw.nV - t - 1);
                        for (auto &m: permMats) m = p * m;
                    }
                }
            }

            // 2: cleaning parmMats and permIdcs to see if there is a constraint or reveal singularity-from-transition
            if (!v.isBoundary()) {
                std::set temp(permIdcs.begin(), permIdcs.end());
                std::vector idcs(temp.begin(), temp.end());
                std::vector<MatXi> mats(idcs.size());

                for (int j = 0; j < idcs.size(); j++) {
                    mats[j] = MatXi::Zero(N, N);
                    for (int k = 0; k < permIdcs.size(); k++) {
                        if (idcs[j] == permIdcs[k]) mats[j] += permMats[k];
                    }
                    if (idcs[j] == v.id) mats[j] -= MatXi::Identity(N, N);
                }

                if (rg::any_of(mats, [](auto &m) { return m.cwiseAbs().maxCoeff() != 0; })) {
                    for (int j = 0; j < mats.size(); j++)
                        assign_block(cT, mats[j], N * currConstraint, N * idcs[j]);
                    currConstraint++;
                }
            }
        }

        vtrans2cut.resize(N * cut.nV, N * nR);
        vtrans2cut.setFromTriplets(vT.begin(), vT.end());
        vtrans2cut.prune(1e-3);

        constraint.resize(N * currConstraint, N * nR);
        constraint.setFromTriplets(cT.begin(), cT.end());
        constraint.prune(1e-3);


        /// filtering out barycentric symmetry, including sign symmetry.
        /// The parameterization should always only include n dof for the surface
        /// Warning: this assumes n divides N!
        /// integer variables are per single "d" packet, and the rounding is done for the N functions with projection over linRed
        std::vector<TripD> buff;
        buff.clear();
        for (int i = 0; i < N * nR; i += N) assign_block(buff, lreductor, i, i * n/N);
        uncompress.resize(N * nR, n * nR);
        uncompress.setFromTriplets(buff.begin(), buff.end());

        fixedIdcs.resize(n);
        if (nS == 0) {
            // no inner singular vertices; vertex 0 is set to (0....0)
            for (int j = 0; j < n; j++) fixedIdcs(j) = j;
        } else {
            // fixing first singularity to (0,....0)
            int iV;
            for (iV = 0; iV < raw.nV; iV++) if (is_inside_singular(raw.verts[iV])) break;
            for (int j = 0; j < n; j++) fixedIdcs(j) = n * iV + j;
        }

        //----- indices to be integer (used for rounding on seams) -----
        integerIdcs.resize(nT * n);
        integerIdcs.setZero();
        for (int i = 0; i < nT; i++) {
        for (int j = 0; j < n; j++) {
            integerIdcs(n * i + j) = n * (raw.nV + i) + j;
        }}

        //----- indices of singular (used for rounding on singulars before rounding on seams) -----
        singularIdcs.resize(n * nS);
        int c = 0;
        for (Vert v: raw.verts | vw::filter([&](auto& e) { return is_inside_singular(e); }) ) {
            for (int j = 0; j < n; j++) singularIdcs(c++) = n * v.id + j;
        }
    }
}
