#ifndef UTIL_H
#define UTIL_H
#include "metriko/core/hmesh/hmesh.h"

namespace metriko::qex {
    struct Qvert {
        Row2d uvw;
        Row3d pos;
        int idx;
        int sid; // simplex id
    };

    class Qport {
    public:
        int idx;
        int vqvert_id;
        int vid;
        int eid;
        int fid;
        Row2d uvw;
        Row2d dir;
        Row3d pos;
        bool isConnected = false;
        int next_id = -1;
        int prev_id = -1;

        Qport(
            const int idx,
            const int vqvert_id,
            const int vid,
            const int eid,
            const int fid,
            const Row2d& uvw,
            const Row2d& dir,
            const Row3d& pos
        ) : idx(idx), vqvert_id(vqvert_id), vid(vid), eid(eid), fid(fid), uvw(uvw), dir(dir), pos(pos) { }
    };

    struct Qedge {
        const Qport& port1;
        const Qport& port2;
        const Mat2d rot;
        const Row2d trs;

        Qedge(
            const Qport& port1,
            const Qport& port2,
            const Mat2d& rot,
            const Row2d& trs
        ) : port1(port1), port2(port2), rot(rot), trs(trs) { }
    };

    class Qhalf {
    public:
        const Qedge &edge;
        const int idx;
        const bool cannonical;
        const Qport& port1() const { return cannonical ? edge.port1 : edge.port2; }
        const Qport& port2() const { return cannonical ? edge.port2 : edge.port1; }

        Qhalf(const Qedge &edge, const int idx, const bool cann): edge(edge), idx(idx), cannonical(cann) { }
    };

    class Qface {
    public:
        std::vector<Qhalf> qhalfs;
    };

    inline double cross(const Row2d &a, const Row2d &b) {
        return a.x() * b.y() - a.y() * b.x();
    }

    inline double cross(const complex &a, const complex &b) {
        return a.real() * b.imag() - a.imag() * b.real();
    }

    inline bool isccw(
        const complex &a,
        const complex &b,
        const complex &c
    ) {
        return cross(b - a, c - a) > 0;
    }

    inline bool isccw(
        const Row2d &a,
        const Row2d &b,
        const Row2d &c
    ) {
        return cross(b - a, c - a) > 0;
    }

    inline Row2d uv2coef(
        const Row2d &uv1,
        const Row2d &uv2,
        const Row2d &uv3,
        const Row2d &target_uv
    ) {
        Row2d vec12 = uv2 - uv1;
        Row2d vec13 = uv3 - uv1;
        Eigen::Matrix2d mat;
        mat << vec12.x(), vec13.x(), vec12.y(), vec13.y();
        return mat.inverse() * (target_uv - uv1).transpose();
    }

    inline Row3d uv2pos(
        const Row2d &uv1,
        const Row2d &uv2,
        const Row2d &uv3,
        const Row3d &pos1,
        const Row3d &pos2,
        const Row3d &pos3,
        const Row2d &target_uv
    ) {
        Row2d coef = uv2coef(uv1, uv2, uv3, target_uv);
        return pos1 + (pos2 - pos1) * coef.x() + (pos3 - pos1) * coef.y();
    }

    inline Row2d uv2coef(
        const Face f,
        const MatXd &cfn,
        const Row2d &target_uv
    ) {
        const Crnr c1 = f.half().crnr();
        const Crnr c2 = f.half().next().crnr();
        const Crnr c3 = f.half().prev().crnr();
        return uv2coef(
            cfn.row(c1.id),
            cfn.row(c2.id),
            cfn.row(c3.id),
            target_uv);
    }

    inline Row3d uv2pos(
        const Face f,
        const MatXd &cfn,
        const Row2d &target_uv
    ) {
        const Crnr c1 = f.half().crnr();
        const Crnr c2 = f.half().next().crnr();
        const Crnr c3 = f.half().prev().crnr();
        return uv2pos(
            cfn.row(c1.id),
            cfn.row(c2.id),
            cfn.row(c3.id),
            c1.vert().pos(),
            c2.vert().pos(),
            c3.vert().pos(),
            target_uv);
    }

    inline Eigen::Matrix2d matching2rot(int matching) {
        Eigen::Matrix2d mat;
        if (matching == 0 || matching > 100) mat << 1, 0, 0, 1; // todo: tmp
        else if (matching == 1) mat << 0, -1, 1, 0;
        else if (matching == 2) mat << -1, 0, 0, -1;
        else mat << 0, 1, -1, 0;
        return mat;
    }

    inline void computeDiff(
        const Hmesh &mesh,
        const MatXd &cfn,
        const std::vector<Row2d> &trans,
        const std::vector<int> &matching,
        const Half h
    ) {
        Crnr c1 = h.next().crnr();
        Crnr c2 = h.twin().prev().crnr();
        Row2d curUV = cfn.row(c1.id);
        Row2d adjUV = cfn.row(c2.id);
        Row2d t = trans[h.id];
        int m = matching[h.id];
        if (t.norm() > 1e-5) {
            Eigen::Matrix2d rot = matching2rot(m);
            Row2d diff = curUV.transpose() - (rot * adjUV.transpose() + t.transpose());
            if (diff.norm() > 1e-10) std::cout << "diff_meshing: " << diff << std::endl;
        }
    }

    inline void compute_trs_matrix(
        const Hmesh &mesh,
        const MatXd &cfn,
        const VecXi &matching,
        const int rosyN,
        VecXi &heMatching,
        MatXd &heTranslation
    ) {
        heMatching.resize(mesh.nH);
        heTranslation.resize(mesh.nH, 2);
        for (Half h: mesh.halfs) {
            Crnr c1 = h.next().crnr();
            Crnr c2 = h.twin().prev().crnr();
            Vec2d uv1 = cfn.row(c1.id).transpose();
            Vec2d uv2 = cfn.row(c2.id).transpose();
            int m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
            m = m < 0 ? (rosyN + m % rosyN) % rosyN : m % rosyN;
            Vec2d t = uv2 - matching2rot(m) * uv1;
            heMatching(h.id) = m;
            heTranslation.row(h.id) = Row2d{std::round(t.x()), std::round(t.y())};
        }
    }
}

#endif
