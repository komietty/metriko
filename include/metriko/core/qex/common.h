#ifndef UTIL_H
#define UTIL_H
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/common/utilities.h"
#include "metriko/core/common/predicates.h"

namespace metriko::qex {
    class Qvert {
    public:
        const complex uv;
        const Row3d pos;
        const int sid;

        Qvert(
            const complex &uv,
            const Row3d &pos,
            const int sid
        ) : uv(uv), pos(pos), sid(sid) {
        }
    };

    class Qport {
    public:
        int idx;
        int vid;
        int eid;
        int fid;
        complex uv;
        complex dir;
        Row3d pos;
        bool isConnected = false;
        int next_id = -1;
        int prev_id = -1;

        Qport(
            const int idx,
            const int vid,
            const int eid,
            const int fid,
            const complex uv,
            const complex dir,
            const Row3d &pos
        ) : idx(idx), vid(vid), eid(eid), fid(fid), uv(uv), dir(dir), pos(pos) {
        }
    };

    class Qedge {
    public:
        const Qport &port1;
        const Qport &port2;

        Qedge(const Qport &port1, const Qport &port2) : port1(port1), port2(port2) {
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

    inline complex convert(Row2d a) { return {a.x(), a.y()}; }

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

    inline double orientation(complex pa, complex pb, complex pc) {
        auto pa_ = std::array{pa.real(), pa.imag()};
        auto pb_ = std::array{pb.real(), pb.imag()};
        auto pc_ = std::array{pc.real(), pc.imag()};
        return orient2dfast(pa_.data(), pb_.data(), pc_.data());
    }

    inline bool is_collinear(complex pa, complex pb, complex pc) {
        double o = orientation(pa, pb, pc);
        return abs(o) < ACCURACY;
    }

    inline bool is_inside_triangle(
        complex pa,
        complex pb,
        complex pc,
        complex uv
    ) {
        return orientation(pa, pb, uv) > 0 &&
               orientation(pb, pc, uv) > 0 &&
               orientation(pc, pa, uv) > 0;
    }

    inline bool is_inside_triangle(
        const Face f,
        const VecXc &cfn,
        const complex uv
    ) {
        auto uv1 = cfn(f.id * 3 + 0);
        auto uv2 = cfn(f.id * 3 + 1);
        auto uv3 = cfn(f.id * 3 + 2);
        return is_inside_triangle(uv1, uv2, uv3, uv);
    }
}

#endif
