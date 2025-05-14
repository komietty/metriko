#ifndef UTIL_H
#define UTIL_H
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/common/utilities.h"

namespace metriko::qex {
    struct Qvert {
        complex uv;
        Row3d pos;
        int idx;
        int sid; // simplex id

        Qvert(
            const double u,
            const double v,
            const Row3d &pos,
            const int idx,
            const int sid
        ) : uv(complex(u, v)), pos(pos), idx(idx), sid(sid) { }

        Qvert(
            const complex &uv,
            const Row3d &pos,
            const int idx,
            const int sid
        ) : uv(uv), pos(pos), idx(idx), sid(sid) { }
    };

    class Qport {
    public:
        int idx;
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
            const int vid,
            const int eid,
            const int fid,
            const Row2d &uvw,
            const Row2d &dir,
            const Row3d &pos
        ) : idx(idx), vid(vid), eid(eid), fid(fid), uvw(uvw), dir(dir), pos(pos) {
        }
    };

    struct Qedge {
        const Qport &port1;
        const Qport &port2;
        const Mat2d rot;
        const Row2d trs;

        Qedge(
            const Qport &port1,
            const Qport &port2,
            const Mat2d &rot,
            const Row2d &trs
        ) : port1(port1), port2(port2), rot(rot), trs(trs) {
        }
    };

    class Qhalf {
    public:
        const Qedge &edge;
        const int idx;
        const bool cannonical;
        const Qport &port1() const { return cannonical ? edge.port1 : edge.port2; }
        const Qport &port2() const { return cannonical ? edge.port2 : edge.port1; }

        Qhalf(const Qedge &edge, const int idx, const bool cann): edge(edge), idx(idx), cannonical(cann) {
        }
    };

    class Qface {
    public:
        std::vector<Qhalf> qhalfs;
    };

    inline Row2d convert(complex a) { return Row2d(a.real(), a.imag()); }
    inline complex convert(Row2d a) { return complex(a.x(), a.y()); }

    inline Eigen::Matrix2d matching2rot(int matching) {
        Eigen::Matrix2d mat;
        if (matching == 0 || matching > 100) mat << 1, 0, 0, 1; // todo: tmp
        else if (matching == 1) mat << 0, -1, 1, 0;
        else if (matching == 2) mat << -1, 0, 0, -1;
        else mat << 0, 1, -1, 0;
        return mat;
    }

    inline void compute_trs_matrix_tmp(
        const Hmesh &mesh,
        const VecXc &cfn,
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
            Vec2d uv1 = Vec2d(cfn(c1.id).real(), cfn(c1.id).imag());
            Vec2d uv2 = Vec2d(cfn(c2.id).real(), cfn(c2.id).imag());
            int m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
            m = m < 0 ? (rosyN + m % rosyN) % rosyN : m % rosyN;
            Vec2d t = uv2 - matching2rot(m) * uv1;
            heMatching(h.id) = m;
            heTranslation.row(h.id) = Row2d{std::round(t.x()), std::round(t.y())};
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
