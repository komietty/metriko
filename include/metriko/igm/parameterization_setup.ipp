//
// Copyright (C) 2025 Saki Komikado <komietty@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
#pragma once
namespace metriko {
    inline void RosyParameterization::compute_he2matching() {
        he2matching.resize(raw.nH);
        for (Half h: raw.halfs) {
            int m = (h.isCanonical() ? -1 : 1) * matching(h.edge().id); // better not inversed
            he2matching[h.id] = (m % N + N) % N ;
        }
    }

    inline void RosyParameterization::compute_he2transidx() {
        he2transidx.resize(raw.nH);
        he2transidx.setConstant(32767);
        vec valence(raw.nV, 0.);
        vec claimed(raw.nE, false);

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

    inline void RosyParameterization::setup() {
        // here we compute a permutation matrix
        vec<MatXi> constParmMats(N, MatXi::Zero(N, N));
        for (int k = 0; k < N; k++)
        for (int i = 0; i < N; i++)
            constParmMats[k]((i + k) % N, i) = 1;

        vec<TripD> vT, cT;
        // forming the constraints and the singularity positions
        int currConstraint = 0;
        // this loop set up the transitions (vector field matching) across the cuts
        for (Vert v: raw.verts) {
            // 1: The initial corner gets the identity without any transition
            vec<MatXi> permMats;
            vec<int> permIdcs;
            permMats.emplace_back(MatXi::Identity(N, N));
            permIdcs.emplace_back(v.id);
            int iVcut = -1;

            // remember uv = Rot * val + Transition
            for (Half h: v.adjHalfs(v.isBoundary() ? find_first_bndr_he(v) : find_first_seam_he(v))) {
                if (h.isBoundary()) break; /// last boundary
                int jVcut = cut.idx(h.face().id, h.next().crnr().id % 3);

                if (jVcut != iVcut) {
                    iVcut = jVcut;
                    for (int i = 0; i < permIdcs.size(); i++) assign_block(vT, permMats[i], N * iVcut, N * permIdcs[i]);
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
                vec idcs(temp.begin(), temp.end());
                vec<MatXi> mats(idcs.size());

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
        vec<TripD> buff;
        for (int i = 0; i < N * nR; i += N) assign_block(buff, lreductor, i, i * n/N);
        uncompress.resize(N * nR, n * nR);
        uncompress.setFromTriplets(buff.begin(), buff.end());

        //----- indices to be integer (used for rounding on seams) -----
        const auto sv = rg::find_if(raw.verts, [&](Vert v) { return is_inside_singular(v); });
        const auto iV = sv == raw.verts.end() ? 0 : sv->id;
        fixedIdcs   = VecXi::LinSpaced(n, n * iV, n * iV + n - 1);
        integerIdcs = VecXi::LinSpaced(n * nT, n * raw.nV, n * (raw.nV + nT) - 1);

        //----- indices of singular (used for rounding on singulars before rounding on seams) -----
        singularIdcs.resize(n * nS);
        int c = 0;
        for (Vert v: raw.verts | vw::filter([&](auto& e) { return is_inside_singular(e); }) ) {
            for (int j = 0; j < n; j++) singularIdcs(c++) = n * v.id + j;
        }
    }
}
