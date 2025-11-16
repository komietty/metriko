//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_EMBEDDING_TUTTE_H
#define METRIKO_EMBEDDING_TUTTE_H
#include "embedding.h"

namespace metriko {
    inline Mat2d compute_rotation(int i) {
        Mat2d r0, r1, r2, r3;
        r0 << 1, 0, 0, 1;
        r1 << 0, -1, 1, 0;
        r2 << -1, 0, 0, -1;
        r3 << 0, 1, -1, 0;
        auto r = std::vector{r2, r1, r0, r3}; // need fix
        return r[i];
    }

    inline MatXd compute_transed_uv(
        const MatXd &uv_curr,
        const MatXd &uv_twin,
        const int side_curr,
        const int side_twin,
        const Half half_curr // evaluating halfedge
    ) {
        auto dif = (side_twin - side_curr + 4) % 4;
        auto rot = compute_rotation(dif);
        MatXd uv_next = uv_twin.row(side_curr) * rot.transpose();
        Half half_twin = half_curr.twin();
        Crnr c0 = half_curr.next().crnr();
        Crnr c1 = half_twin.prev().crnr();
        Row2d offset = uv_next.row(c1.id) - uv_curr.row(c0.id);
        uv_next.rowwise() += offset;
        return uv_next;
    }

    inline complex compute_translation(
        const Tmesh &tmesh,
        const VecXd &X,
        const Tquad &curr_tq,
        const Tquad &twin_tq,
        const Thalf &flip_th,
        const complex origin
    ) {
        Thalf th = tmesh.thalfs[curr_tq.find_first_thid(0)];
        auto dir = complex(1, 0);
        auto sum = origin;

        while (true) {
            //std::cout << "add curr: " << X[th.edge().id] * dir << std::endl;
            sum += X[th.edge().id] * dir;
            if (th.id == flip_th.id) break;
            if (curr_tq.find_side(th) != curr_tq.find_side(th.next()))
                dir *= complex(0, 1);
            th = th.next();
        }

        th = th.twin();
        dir *= -1;

        while (true) {
            if (th.id == twin_tq.find_first_thid(0)) break;
            //std::cout << "add twin: " << X[th.edge().id] * dir << std::endl;
            sum += X[th.edge().id] * dir;
            if (twin_tq.find_side(th) != twin_tq.find_side(th.next()))
                dir *= complex(0, 1);
            th = th.next();
        }

        return sum;
    }


    inline MatXd embedding_tutte_for_tquad(
        const int tqid,
        const Hmesh &hmesh,
        const Tmesh &tmesh,
        const std::vector<EmbeddedTEdge> &etes,
        const VecXd &X
        //const Mat2d &rot,
        //const Row2d &oft
    ) {
        const Tquad &tq = tmesh.tquads[tqid];
        std::vector<std::pair<EmbeddedTEdge, bool> > tq_etes;
        MatXd embeded_uv = MatXd::Zero(hmesh.nV, 2);
        for (int thid: tq.thids) {
            const auto &th = tmesh.thalfs[thid];
            tq_etes.emplace_back(etes[th.edge().id], th.cannonical);

            // assign values to vertices
        }

        std::vector<glm::vec3> ps_;
        std::vector<double> vert_val_u;
        std::vector<double> vert_val_v;
        complex dir = complex(1, 0);
        complex sum = complex(0, 0);
        for (int i = 0; i < 4; i++) {
            auto thids = tq.thids_by_side(i);
            for (int thid: thids) {
                const auto &th = tmesh.thalfs[thid];
                const auto &te = th.edge();
                const auto &ete = etes[te.id];
                for (int j = 0; j < ete.vids.size(); j++) {
                    int vid = ete.vids[j];
                    double val = th.cannonical ? ete.vals[j] : X[te.id] - ete.vals[j];
                    Row3d pos = hmesh.pos.row(vid);
                    ps_.emplace_back(pos.x(), pos.y(), pos.z());
                    embeded_uv(vid, 0) = val * dir.real() + sum.real();
                    embeded_uv(vid, 1) = val * dir.imag() + sum.imag();
                    vert_val_u.emplace_back(val * dir.real() + sum.real());
                    vert_val_v.emplace_back(val * dir.imag() + sum.imag());
                }
                sum += X[te.id] * dir;
            }
            dir *= complex(0, 1);
        }

        //--- find all faces ---
        std::queue<int> queue;
        std::unordered_set<int> visit;
        for (const auto &[ete, cann]: tq_etes) {
            for (Half h: ete.halfs) {
                if (!cann) h = h.twin();
                Face f = h.face();
                queue.emplace(f.id);
                if (!visit.contains(f.id))visit.emplace(f.id);
            }
        }

        while (!queue.empty()) {
            Face ff = hmesh.faces[queue.front()];
            queue.pop();
            //assert(rg::find(visit, ff) == visit.end());
            for (Half hh: ff.adjHalfs()) {
                Face fh = hh.twin().face();
                //std::cout << "face id: " << fh.id << std::endl;
                if (visit.contains(fh.id)) continue;

                bool f2 = true;
                for (EmbeddedTEdge &ete: tq_etes | vw::keys) {
                    //std::cout << "ete flag: " << !ete.contains(hh.edge()) << std::endl;
                    f2 &= !ete.contains(hh.edge());
                }
                //std::cout << "flag: " << f2 << std::endl;
                if (f2) {
                    //std::cout << "face id: " << fh.id << std::endl;
                    queue.emplace(fh.id);
                    visit.emplace(fh.id);
                }
            }
        }

        std::vector<glm::vec3> ps;
        int counter = 0;
        std::vector<int> vidvec;
        for (int id: visit) {
            Face face = hmesh.faces[id];
            for (Half h: face.adjHalfs()) {
                vidvec.emplace_back(h.tail().id);
            }
        }

        std::set vidset(vidvec.begin(), vidvec.end());
        vidvec.clear();
        vidvec = std::vector(vidset.begin(), vidset.end());

        MatXi face_table = MatXi::Zero(visit.size(), hmesh.nF);
        MatXd vert_table = MatXd::Zero(vidset.size(), hmesh.nV);
        std::unordered_map<int, int> idcs_table;
        for (int i = 0; i < vidvec.size(); i++) {
            vert_table(i, vidvec[i]) = 1;
            idcs_table[vidvec[i]] = i;
        }

        for (int id: visit) {
            Face ff = hmesh.faces[id];
            ps.emplace_back(ff.center().x(), ff.center().y(), ff.center().z());
            face_table(counter, id) = 1;
            counter++;
        }

        MatXd V = vert_table * hmesh.pos;
        MatXi F = face_table * hmesh.idx;
        MatXd UV = vert_table * embeded_uv;

        for (int i = 0; i < F.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                F(i, j) = idcs_table[F(i, j)];
            }
        }

        // here tutte's parameterization
        auto m = std::make_unique<Hmesh>(V, F);
        SprsD L = cotan_laplacian(*m);
        SprsD M = mass_matrix(*m);
        SprsD BL = boundary_snap_laplacian(*m);
        MatXd uv(m->nV, 2);
        {
            Eigen::SparseLU<SprsD> lu;
            lu.compute(BL);
            VecXd res = lu.solve(UV.col(0));
            uv.col(0) = res;
        }
        {
            Eigen::SparseLU<SprsD> lu;
            lu.compute(BL);
            VecXd res = lu.solve(UV.col(1));
            uv.col(1) = res;
        }

        //auto s = polyscope::registerSurfaceMesh("ebd mesh " + std::to_string(tqid), V, F);
        //s->setEdgeWidth(1);
        //auto uvw = s->addVertexParameterizationQuantity("uv", uv);

        //MatXd uv_ = uv * rot.transpose();
        //uv_.rowwise() += oft;

        //auto uvw2 = s->addVertexParameterizationQuantity("uv2", uv_);
        //uvw->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        //uvw->setCheckerSize(1);
        //uvw->setEnabled(false);
        //uvw2->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        //uvw2->setCheckerSize(1);
        //uvw2->setEnabled(true);

        MatXd uv_all = MatXd::Zero(hmesh.nC, 2);
        //MatXd uv_vrt = vert_table.transpose() * uv_;
        MatXd uv_vrt = vert_table.transpose() * uv;
        for (int iF: visit) {
            Face f = hmesh.faces[iF];
            for (Half h: f.adjHalfs()) {
                uv_all.row(h.crnr().id) = uv_vrt.row(h.next().head().id);
            }
        }
        return uv_all;
    }

    inline VecXi half_to_ethalf(
        const Hmesh &hmesh,
        const std::vector<EmbeddedTHalf> &eths
    ) {
        VecXi table(hmesh.nH, -1);
        for (Half h: hmesh.halfs) {
            for (int i = 0; i < eths.size(); i++) {
                if (eths[i].contains(h)) {
                    table[h.id] = i;
                    break;
                }
            }
        }
        return table;
    }

    struct EdgeHash {
        std::size_t operator()(const Edge &e) const noexcept { return std::hash<int>{}(e.id); }
    };

    struct HalfHash {
        std::size_t operator()(const Half &h) const noexcept { return std::hash<int>{}(h.id); }
    };

    // take the tutte result as the input, embed it until seam intersection.
    // computes halfedges to search with in the next loop at the same time.
    inline std::unordered_set<Half, HalfHash> sequential_mapping(
        const MatXd &uv_in,
        const Half half_in,
        //const EmbeddedTHalf &ethf_in,
        const std::vector<bool> &seam,
        std::vector<bool> &flag, // the flag to check a corner is already checked or not
        const std::vector<Half> &tquad_boundary,
        MatXd &uv_all
    ) {
        std::queue<Half> queue;
        std::unordered_set<Edge, EdgeHash> visit;
        std::unordered_set<Half, HalfHash> nextH; // the halfedges to the other tquad
        queue.push(half_in);
        visit.emplace(half_in.edge());
        if (flag[half_in.crnr().id]) return nextH;

        while (queue.size() > 0) {
            Half hfr = queue.front();
            Face f = hfr.face();
            Crnr c0 = f.half().crnr();
            Crnr c1 = f.half().next().crnr();
            Crnr c2 = f.half().prev().crnr();

            uv_all.row(c0.id) = uv_in.row(c0.id);
            uv_all.row(c1.id) = uv_in.row(c1.id);
            uv_all.row(c2.id) = uv_in.row(c2.id);
            flag[c0.id] = true;
            flag[c1.id] = true;
            flag[c2.id] = true;
            queue.pop();
            for (Half h: f.adjHalfs()) {
                // 0: if hit seam, just stops
                if (seam[h.edge().id]) continue;
                // 1: if hit the visited edge, just stop
                if (visit.contains(h.edge())) continue;

                // 2: if hit the excluded edge, just stops
                //if (rg::find(edge_vst, h.edge()) != edge_vst.end()) continue;
                // 3: if hit ethalf where it comes from, just stops
                //if (!ethf_in.contains(h)) { edge_exclude.emplace_back(h.edge()); continue; }

                // 4: if hit boundary, puts it as a bridge to the next tquad
                if (rg::find(tquad_boundary, h) != tquad_boundary.end()) {
                    nextH.insert(h.twin());
                    continue;
                }
                // 5: inside tquad. add it to the queue
                visit.emplace(h.edge());
                queue.push(h.twin());
            }
        }
        return nextH;
    }

    // try to multiply rotation until halfedge coner values corresponds
    // need to consider: is there any possibility of flip?
    inline void apply_transition(
        const Half h, // the halfedge of unfixed side
        const MatXd &mat0, // the fixed uv information
        MatXd &mat1 // the unfixed adjacent uv information
    ) {
        Half h0 = h.twin();
        Half h1 = h;
        Row2d uv0 = mat0.row(h0.next().crnr().id);
        Row2d uv1 = mat1.row(h1.prev().crnr().id);
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
        throw new std::runtime_error("no corresponding rotation found");
    }

    inline void apply_transition_test(
        const Half h, // the halfedge of unfixed side
        const MatXd &mat0, // the fixed uv information
        MatXd &mat1 // the unfixed adjacent uv information
    ) {
        Half h0 = h.twin();
        Half h1 = h;
        Row2d uv0  = mat0.row(h0.next().crnr().id);
        Row2d uv1  = mat1.row(h1.prev().crnr().id);
        Row2d uv0a = mat0.row(h0.prev().crnr().id);
        Row2d uv1a = mat1.row(h1.next().crnr().id);
        std::cout << "corner id a: " << h0.next().crnr().id << std::endl;
        std::cout << "corner id b: " << h0.prev().crnr().id << std::endl;
        std::cout << "uv0: " << uv0 << std::endl;
        std::cout << "uv1: " << uv1 << std::endl;
        std::cout << "uv0a: " << uv0a << std::endl;
        std::cout << "uv1a: " << uv1a << std::endl;

        //for (int i = 0; i < 4; i++) {
        //    Mat2d rot = compute_rotation(i);
        //    Row2d res = rot * (uv1a - uv1).transpose() + uv0.transpose();
        //    if ((res - uv0a).norm() < 1e-6) {
        //        for (int j = 0; j < mat1.rows(); j++) {
        //            mat1.row(j) = rot * (mat1.row(j) - uv1).transpose() + uv0.transpose();
        //        }
        //        return;
        //    }
        //}
        //throw new std::runtime_error("no corresponding rotation found");
    }

    inline std::pair<int, int> find_tqid_and_thid_from_half(
        const Tmesh &tmesh,
        const std::vector<EmbeddedTHalf> &eths,
        const Half &h
    ) {
        for (auto tq: tmesh.tquads) {
            for (int thid: tq.thids) {
                if (eths[thid].contains(h)) return std::pair(tq.id, thid);
            }
        }
        throw new std::runtime_error("half not found in any tquad");
    }

    inline MatXd compute_tutte_parameterization(
        const Hmesh &hmesh,
        const Tmesh &tmesh,
        const std::vector<EmbeddedTEdge> &etes,
        const std::vector<EmbeddedTHalf> &eths,
        const std::vector<bool> &seam,
        const VecXd &X
    ) {
        // compute uv per tquad first...
        std::vector<MatXd> uv_per_tquad;
        for (const auto &tq: tmesh.tquads) {
            MatXd uv = embedding_tutte_for_tquad(tq.id, hmesh, tmesh, etes, X);
            uv_per_tquad.emplace_back(uv);
        }

        MatXd uv_all = MatXd::Zero(hmesh.nC, 2);
        auto corner_flag = std::vector(hmesh.nC, false);
        std::unordered_set<Half, HalfHash> halfs_out;
        std::unordered_set<Edge, EdgeHash> edges_visit;
        std::stack<Half> queue; // todo: when using queue, it does not work in count 519...

        // for the first tquad...
        {
            Half h = hmesh.halfs[eths[tmesh.tquads[0].thids.front()].halfs.front().id];
            auto tqid_curr = find_tqid_and_thid_from_half(tmesh, eths, h).first;
            auto boundaries = std::vector<Half>{};
            for (int thid: tmesh.tquads[tqid_curr].thids)
                for (auto hh: eths[thid].halfs) boundaries.emplace_back(hh);

            auto next_halfs = sequential_mapping(
                uv_per_tquad[tqid_curr],
                h,
                seam,
                corner_flag,
                boundaries,
                uv_all
            );
            for (auto nh: next_halfs) queue.push(nh);
        }

        int count = 0;
        while (!rg::all_of(corner_flag, [&](auto f) { return f; }) && count < 1519) {
            auto h = queue.top();
            queue.pop();
            auto tqid_curr = find_tqid_and_thid_from_half(tmesh, eths, h).first;
            auto tqid_prev = find_tqid_and_thid_from_half(tmesh, eths, h.twin()).first;

            if (count == 519) {
                std::cout << "tail: " << h.tail().pos() << ", head: " << h.head().pos() << std::endl;
                std::cout << "face id: " << h.twin().face().id << std::endl;
                std::cout << "crnr id: " << h.twin().crnr().id << std::endl;

                for (auto tq: tmesh.tquads) {
                    for (int thid: tq.thids) {
                        if (eths[thid].contains(h)) std::cout << "thid: " << thid << std::endl;
                    }
                }
                std::cout << "count: " << count << std::endl;

                apply_transition_test(
                    h,
                    uv_per_tquad[tqid_prev],
                    uv_per_tquad[tqid_curr]
                );
            } else {
                apply_transition(
                    h,
                    uv_per_tquad[tqid_prev],
                    uv_per_tquad[tqid_curr]
                );
            }

            auto boundaries = std::vector<Half>{};
            for (int thid: tmesh.tquads[tqid_curr].thids)
                for (auto hh: eths[thid].halfs) boundaries.emplace_back(hh);

            auto next_halfs = sequential_mapping(
                uv_per_tquad[tqid_curr],
                h,
                seam,
                corner_flag,
                boundaries,
                uv_all
            );
            for (auto nh: next_halfs) queue.push(nh);
            count++;
        }
        /*
        */


        /*
        std::vector<std::pair<int, int> > next_tqid_hid_pair;
        for (Half h: halfs_out) {
            auto [tqid, thid] = find_tqid_and_thid_from_half(tmesh, eths, h);
            //std::cout << "candidate: " << tqid  << std::endl;
            if (rg::none_of(next_tqid_hid_pair, [&](auto pair) { return pair.first == tqid; })) {
                next_tqid_hid_pair.emplace_back(tqid, h.id);
                //std::cout << "emplace: " << tqid << ", h.id: " << h.id << std::endl;
            }
            edges_visit.emplace(h.edge());
        }


        for (auto [tqid, hid]: next_tqid_hid_pair) {
            //std::cout << "next half: " << next_half.id << std::endl;
            //std::cout << "next tqid: " << tqid << std::endl;
            const Half next_half = hmesh.halfs[hid];
            const Tquad tq = tmesh.tquads[tqid];
            apply_transition(
                next_half,
                uv_per_tquad[0], // must be changed when for loop
                uv_per_tquad[tqid]
            );

            auto boundaries = std::vector<Half>{};
            for (int thid: tq.thids) {
                for (auto hh: eths[thid].halfs) { boundaries.emplace_back(hh); }
            }

            sequential_mapping(
                uv_per_tquad[tqid],
                next_half,
                seam,
                boundaries,
                halfs_out,
                edges_visit,
                uv_all
            );
        }
        */


        //std::vector<glm::vec3> ns;
        //std::vector<std::array<size_t, 2> > es;
        //size_t counter = 0;
        //for (auto hh: halfs_out) {
        //    Row3d p1 = hh.tail().pos();
        //    Row3d p2 = hh.head().pos();
        //    ns.emplace_back(p1.x(), p1.y(), p1.z());
        //    ns.emplace_back(p2.x(), p2.y(), p2.z());
        //    es.emplace_back(std::array{counter, counter + 1});
        //    counter += 2;
        //}
        //auto c = polyscope::registerCurveNetwork("halfs out", ns, es);
        //c->setEnabled(true);
        //c->resetTransform();
        //c->setRadius(0.0015);
        return uv_all;

        //VecXi h2eth = half_to_ethalf(hmesh, eths);
        //Half half = eths[0].halfs[0];
        //auto eth = eths[h2eth[half.id]];
        //auto tqid = tmesh.th2quad[eth.thid];
        //std::queue<std::tuple<int, int, complex>> queue;
    }

    inline void embedding_tutte(
        const Hmesh &mesh,
        const Tmesh &tmesh,
        const std::vector<EmbeddedTEdge> &etes,
        const VecXd &X
    ) {
        std::vector<int> tqids_;
        tqids_.emplace_back(12);
        std::queue<std::tuple<int, int, complex> > queue;
        std::unordered_set<int> visit;
        complex origin = complex(10, 10);
        int bgn = 12;
        queue.emplace(bgn, 2, origin);
        visit.emplace(bgn);
        int count = 0;

        Mat2d rot0;
        Mat2d rot1;
        Mat2d rot2;
        Mat2d rot3;
        rot0 << 1, 0, 0, 1;
        rot1 << 0, -1, 1, 0;
        rot2 << -1, 0, 0, -1;
        rot3 << 0, 1, -1, 0;
        auto rots = std::vector{rot2, rot1, rot0, rot3};

        MatXd uv_all(mesh.nC, 2);

        while (!queue.empty() && count++ < 10) {
            auto [tqid, rot_id, oft] = queue.front();
            const Tquad &tq = tmesh.tquads[tqid];
            queue.pop();

            //MatXd uv_for_tquad = embedding_tutte_for_tquad(
            //    mesh, tmesh, tqid, etes, X, rots[rot_id],
            //    Row2d(oft.real(), oft.imag()));
            //sequential_mapping()

            //uv_all +=

            for (int thid: tq.thids) {
                //if (tq.id != bgn) continue;
                const Thalf &curr_th = tmesh.thalfs[thid];
                const Thalf &twin_th = curr_th.twin();
                auto it = rg::find_if(tmesh.tquads, [&](const Tquad &tq_) {
                    return rg::find(tq_.thids, twin_th.id) != tq_.thids.end();
                });
                if (it != tmesh.tquads.end()) {
                    int curr_side = tq.find_side(curr_th);
                    int twin_side = it->find_side(twin_th);
                    int diff = (twin_side - curr_side + 4) % 4;
                    //std::cout << "curr side, twin side: " << curr_side << ", " << twin_side << std::endl;

                    complex res = compute_translation(tmesh, X, tq, *it, curr_th, origin);
                    if (!visit.contains(it->id)) {
                        queue.emplace(it->id, diff, res);
                        visit.emplace(it->id);
                    }
                }
            }
        }

        /// ---- visualize mesh ---- ///
        {
            const auto surf = polyscope::registerSurfaceMesh("mesh_", mesh.pos, mesh.idx);
            const auto prms = surf->addParameterizationQuantity("params", uv_all);
            surf->setEnabled(false);
            prms->setStyle(polyscope::ParamVizStyle::GRID);
            //prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
            prms->setCheckerSize(1);
        }
    }

    /*
    */
}

#endif
