//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_QUANTIZATION_RELINEARIZATION_H
#define METRIKO_QUANTIZATION_RELINEARIZATION_H

namespace metriko {
    // back up
    std::optional<Half> find_half(
        const VecXc& cfn,
        const motorcycle::Msgmt &s,
        const complex uv
    ) {
        Face f = s.face;
        for (Half h: f.adjHalfs()) {
            complex uv_tail = cfn(h.next().crnr().id);
            complex uv_head = cfn(h.prev().crnr().id);
            if (abs(uv - uv_tail) < 1e-6 || abs(uv - uv_head) < 1e-6) continue;
            double argA = std::arg(uv_head - uv_tail);
            double argB = std::arg(uv      - uv_tail);
            if (std::abs(argA - argB) < 1e-6) return std::optional(h);
        }
        return std::nullopt;
    }

    Mat2d mid2rot(const int mid) {
        Mat2d m;
        switch (mid) {
            case 0: m = Mat2d::Identity(); break;
            case 1: m << 0, 1, -1, 0;      break;  // todo: it follows he2matching order but counter intuitive...
            case 2: m << -1, 0, 0, -1;     break;
            case 3: m << 0, -1, 1, 0;      break;  // todo: it follows he2matching order but counter intuitive...
            default: throw std::exception();
        }
        return m;
    }

    struct Constraint {
        int vidFr;
        int vidTo;
        int sumX;
        int sumY;
        std::vector<tmesh::Thalf> thalfs;
        std::vector<Half> halfs;

        Constraint(
            const int vidFr_,
            const int vidTo_,
            const std::vector<tmesh::Thalf>& thalfs_
        ): vidFr(vidFr_), vidTo(vidTo_), thalfs(thalfs_) {
            for (auto& th: thalfs) {
                auto te = th.edge();
                auto ar = te.passing_halfedges();
                if (th.cannonical) {
                    for (Half &h: ar) halfs.emplace_back(h);
                } else {
                    for (Half &h: ar | std::views::reverse)
                        halfs.emplace_back(h.twin());
                }
            }
        }

        std::optional<MatXd> calc_trs(
            const Hmesh& mesh,
            const VecXi& he2matching,
            const VecXi& he2transidx,
            const VecXd& fullx,
            const VecXc& cfn,
            const int nV,
            Mat2d& first_rot
        ) {
            MatXd trs(2, fullx.size());
            trs.setZero();
            std::vector<int> ms;
            std::vector<int> ts;

            { // skip fist trans case
                Vec2d uv1 = fullx.segment(vidFr * 2, 2);
                Vec2d uv2 = Vec2d(
                    thalfs.front().uv_fr().real(),
                    thalfs.front().uv_fr().imag()
                );
                Row2d d = (uv2 - uv1).transpose();
                //if (d.norm() > 1e-6) {return std::nullopt;}
            }
            {
                Vec2d uv1 = Vec2d(
                    thalfs.back().uv_to().real(),
                    thalfs.back().uv_to().imag()
                    );
                Vec2d uv2 = fullx.segment(vidTo * 2, 2);
                Row2d d = (uv2 - uv1).transpose();
                //if (d.norm() < 1e-6) {return std::nullopt;}
            }

            get_first_trs(mesh, he2matching, he2transidx, fullx, cfn, ms, ts, first_rot);

            for (Half h: halfs) {
                int mid = he2matching(h.id);
                int tid = he2transidx(h.id);
                if (tid < 1e3) { // todo: wired number, need fix
                    ms.emplace_back(mid);
                    ts.emplace_back(tid);
                    //std::cout << "tid, mid: " << tid << ", " << mid << std::endl;
                    //if (tid > 0) return std::nullopt;
                    //std::cout << ">>> mid: " <<  mid << ", tid: " << tid << std::endl;
                    //std::cout << std::fixed << std::setprecision(3);
                    //std::cout << "h tail pos: " << h.tail().pos() << std::endl;
                    //std::cout << "h head pos: " << h.head().pos() << std::endl;
                    //std::cout.unsetf(std::ios::fixed);
                    //Row2d t = fullx.segment(nV * 2 + (abs(tid) - 1) * 2, 2);
                    //if (tid < 0) { std::cout << "  > tid < 0, trans: " << t.x() << ", " << t.y() << std::endl; }
                    //else         { std::cout << "  > tid > 0, trans: " << t.x() << ", " << t.y() << std::endl; }
                }
            }

            get_last_trs(mesh, he2matching, he2transidx, fullx, cfn, ms, ts);


            // compute vert_to rotation
            Mat2d rot_vid_to = Mat2d::Identity();
            for (int i = ts.size() - 1; i >= 0; i--) {
                rot_vid_to = mid2rot(ms[i]).inverse() * rot_vid_to;
            }

            // compute trans rotation
            for (int i = 0; i < ts.size(); i++) {
                Mat2d r = Mat2d::Identity();
                int t = ts[i];
                int m = ms[i];
                if (t > 0) { r = mid2rot(m).inverse(); }
                for (int j = i; j >= 0; j--) r = mid2rot(ms[j]) * r;
                if (t < 0) {
                    const Mat2d R = -r.inverse();
                    trs.block(0, nV * 2 + (-t - 1) * 2, 2, 2) << R + trs.block(0, nV * 2 + (-t - 1) * 2, 2, 2);
                    //std::cout << "tid : " << t << std::endl;
                    //std::cout << R << std::endl;
                    //Row2d tt = fullx.segment(nV * 2 + (-t - 1) * 2, 2);
                    //std::cout << "trs : " << complex(tt.x(), tt.y()) << std::endl;
                }
                else {
                    const Mat2d R = r.inverse();
                    trs.block(0, nV * 2 + ( t - 1) * 2, 2, 2) << R + trs.block(0, nV * 2 + ( t - 1) * 2, 2, 2);
                    //std::cout << "tid : " << t << std::endl;
                    //std::cout << R << std::endl;
                    //Row2d tt = fullx.segment(nV * 2 + (t - 1) * 2, 2);
                    //std::cout << "trs : " << complex(tt.x(), tt.y()) << std::endl;
                }
            }
            //std::cout << "vidFr: " << std::endl;
            //std::cout << -mid2rot(0) << std::endl;
            //std::cout << "vidTo: " << std::endl;
            //std::cout << rot_vid_to << std::endl;
            trs.block(0, vidFr * 2, 2, 2) << -mid2rot(0);
            trs.block(0, vidTo * 2, 2, 2) << rot_vid_to;
            return std::optional(trs);
        }

        void get_first_trs(
            const Hmesh& mesh,
            const VecXi& he2matching,
            const VecXi& he2transidx,
            const VecXd& fullx,
            const VecXc& cfn,
            std::vector<int>& mids,
            std::vector<int>& tids,
            Mat2d& first_rot
        ) {
            // get trans from uv1 to uv2
            Vec2d uv1 = fullx.segment(vidFr * 2, 2);
            Vec2d uv2 = Vec2d(
                thalfs.front().uv_fr().real(),
                thalfs.front().uv_fr().imag()
            );
            //Row2d d = (uv2 - uv1).transpose();
            //std::cout << ">>> fr uv diff: " << d << std::endl;

            Vert v = mesh.verts[vidFr];
            first_rot = Mat2d::Identity();

            // find out which half to start
            Half bgnH = v.half();
            for (const Half& h: v.adjHalfs()) {
                complex uv = cfn(h.next().crnr().id);
                if (abs(uv - complex(uv2.x(), uv2.y())) < 1e-6) { bgnH = h; break; }
            }


            for (const Half& h: v.adjHalfs(bgnH)) {
                if ((uv2 - uv1).norm() < 1e-6) break;
                int mid = he2matching[h.twin().id];
                int tid = he2transidx[h.twin().id];
                first_rot = mid2rot(mid) * first_rot;
                if (tid < 1e3) { // todo: wired number, need fix
                    mids.emplace_back(mid);
                    tids.emplace_back(tid);
                    Row2d t = fullx.segment(mesh.nV * 2 + (abs(tid) - 1) * 2, 2);
                    if (tid < 0) uv1 = mid2rot(mid) * uv1 + t.transpose();
                    else         uv1 = mid2rot(mid) * uv1 - mid2rot(mid) * t.transpose();
                    //std::cout << "  > mid: " << mid << ", tid: " << tid << ", trans : " << t.x() << ", " << t.y() << std::endl;
                }
            }
        }

        void get_last_trs(
            const Hmesh& mesh,
            const VecXi& he2matching,
            const VecXi& he2transidx,
            const VecXd& fullx,
            const VecXc& cfn,
            std::vector<int>& mids,
            std::vector<int>& tids
        ) const {
            Vec2d uv1 = Vec2d(
                thalfs.back().uv_to().real(),
                thalfs.back().uv_to().imag()
                );
            Vec2d uv2 = fullx.segment(vidTo * 2, 2);
            //Row2d d = (uv2 - uv1).transpose();
            //std::cout << ">>> to uv diff: " << d << std::endl;

            // find out which half to start
            Vert v = mesh.verts[vidTo];
            Half bgnH = v.half();
            for (const Half& h: v.adjHalfs()) {
                complex uv = cfn(h.next().crnr().id);
                if (abs(uv - complex(uv1.x(), uv1.y())) < 1e-6) { bgnH = h; break; }
            }

            for (const Half& h: v.adjHalfs(bgnH)) {
                if ((uv2 - uv1).norm() < 1e-6) break;
                int mid = he2matching[h.id];
                int tid = he2transidx[h.id];
                if (tid < 1e3) { // todo: wired number, need fix
                    mids.emplace_back(mid);
                    tids.emplace_back(tid);
                    Row2d t = fullx.segment(mesh.nV * 2 + (abs(tid) - 1) * 2, 2);
                    if (tid < 0) uv1 = mid2rot(mid) * uv1 + t.transpose();
                    else         uv1 = mid2rot(mid) * uv1 - mid2rot(mid) * t.transpose();
                    //std::cout << "  > mid: " << mid << ", tid: " << tid << ", trans : " << t.x() << ", " << t.y() << std::endl;
                }
            }
        }

        complex calc_sum(const VecXd& X, Mat2d& rot) const {
            complex sum = complex(0, 0);
            assert(thalfs.front().cannonical); // starting th is always cannonical;
            complex dif = thalfs.front().sg_fr().diff();
            complex dir = dif / abs(dif);
            Vec2d dir_ = rot.inverse() * Vec2d(dir.real(), dir.imag());
            dir = complex(dir_.x(), dir_.y());
            for (int i = 0; i < thalfs.size(); i++) {
                //std::cout<< "dir: " << dir << std::endl;
                auto& curr = thalfs[i];
                //std::cout<< "teid: " << curr.edge().id << std::endl;
                //std::cout<< "x:    " << X[curr.edge().id] << std::endl;
                sum += X[curr.edge().id] * dir;
                if (i == thalfs.size() - 1) break;
                auto& next = thalfs[i + 1];
                complex dir1 = curr.dif_to() / abs(curr.dif_to());
                complex dir2 = next.dif_fr() / abs(next.dif_fr());
                //std::cout<< "arg: " << arg(dir2 / dir1) << std::endl;
                dir = (dir2 / dir1) * dir;
            }
            return sum;
        }

    };

    inline std::optional<Constraint> generate_constraint(
        const tmesh::Tmesh &tmesh,
        const VecXd& fullx,
        const VecXd &X,
        const int thid_init,
        const int vid_to_
    ) {
        std::queue<int> q;
        std::unordered_map<int, int> m;
        std::vector<bool> visited;
        visited.assign(tmesh.nTE, false);
        q.emplace(thid_init);

        //if (tmesh.th2sing[thid_init] != 504) return std::nullopt;

        while (!q.empty()) {
            const auto &th = tmesh.thalfs[q.front()];
            const auto &te = th.edge();
            q.pop();

            int vid_fr = tmesh.th2sing[thid_init];
            int vid_to = tmesh.th2sing[th.twin().id];
            if (!th.cannonical &&
                te.seg_fr.id == 0 &&
                vid_to != vid_fr
                //&& tmesh.th2sing[th.twin().id] == 615  // todo: test
            ) {
                std::vector<tmesh::Thalf> ths;
                ths.emplace_back(th);
                //std::cout << "----- vid_fr, vid_to: " << vid_fr << ", " << vid_to << ", ";
                //std::cout << std::fixed << std::setprecision(0)
                //    << "uv fr, to: " << complex(fullx(vid_fr * 2), fullx(vid_fr * 2 + 1))
                //    << ", "          << complex(fullx(vid_to * 2), fullx(vid_to * 2 + 1))
                //    << " -------" << std::endl;
                int thid = th.id;
                while (true) {
                    thid = m[thid];
                    ths.emplace_back(tmesh.thalfs[thid]);
                    if (thid == thid_init) break;
                }
                std::ranges::reverse(ths);
                //for (auto th_: ths) {
                //    std::cout << th_.id << ", ";
                //}
                //std::cout << std::endl;
                //for (auto th_: ths) {
                //    std::cout << X[th_.edge().id] * (th_.cannonical ? 1. : -1.)  << ", ";
                //}
                //std::cout << std::endl;

                return Constraint{ vid_fr, vid_to, ths };
            }

            for (auto pair: th.adj_thalfs()) {
                int teid = pair.edge().id;
                if (!visited[teid]) {
                    visited[teid] = true;
                    m.emplace(pair.id, th.id);
                    q.emplace(pair.id);
                }
            }
        }
        return std::nullopt;
        //throw std::runtime_error("Not implemented");
    }

    inline SprsD generate_constraint_matrix(
        const tmesh::Tmesh &tmesh,
        const VecXd &xcurr
    ) {
        for (auto &th: tmesh.thalfs) {
            if (th.is_stemming_fr_singular()) {
                std::vector<tmesh::Thalf> path;
                //... find shortest way to other singular, generate std::vector<Thalf>


                // then calc diff, find vfr and vto, accumulate transitions
                // int vfr = th.sg_fr(tmesh.tedges).edge->port.vert.id;
                for (auto &th: path) {
                    //auto& fr = th.sg_fr();
                    //for (th)
                }
            }
        }
    }

    inline bool validate_constraint_matrix(tmesh::Tmesh &tmesh) {
    }
}

#endif
