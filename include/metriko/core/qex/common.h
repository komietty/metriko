#ifndef metriko_QEX_COMMON_H
#define metriko_QEX_COMMON_H
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/common/utilities.h"

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

    inline void compute_trs_matrix(
        const Hmesh &mesh,
        const VecXc &cfn,
        const VecXi &matching,
        const int rosyN,
        VecXc &heR,
        VecXc &heT
    ) {
        heR.resize(mesh.nH);
        heT.resize(mesh.nH);
        for (Half h: mesh.halfs) {
            Crnr c1 = h.next().crnr();
            Crnr c2 = h.twin().prev().crnr();
            auto uv1 = cfn(c1.id);
            auto uv2 = cfn(c2.id);
            int m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
            m = m < 0 ? (rosyN + m % rosyN) % rosyN : m % rosyN;
            auto r = get_quater_rot(m);
            auto t = uv2 - r * uv1;
            heR(h.id) = r;
            heT(h.id) = complex(std::round(t.real()), std::round(t.imag()));
        }
    }

    inline double orientation(complex pa, complex pb, complex pc) {
        double acx = pa.real() - pc.real();
        double bcx = pb.real() - pc.real();
        double acy = pa.imag() - pc.imag();
        double bcy = pb.imag() - pc.imag();
        return acx * bcy - acy * bcx;

        // equivalently...
        // auto pa_ = std::array{pa.real(), pa.imag()};
        // auto pb_ = std::array{pb.real(), pb.imag()};
        // auto pc_ = std::array{pc.real(), pc.imag()};
        // return orient2d(pa_.data(), pb_.data(), pc_.data());
    }

    inline bool is_collinear(complex pa, complex pb, complex pc) {
        double o = orientation(pa, pb, pc);
        return abs(o) < ACCURACY;
    }

    inline bool is_points_into(complex p1, complex p2, complex p3, complex uv) {
        return orientation(p1, p2, uv) > ACCURACY && orientation(p1, p3, uv) < ACCURACY;
    }

    inline bool is_inside_triangle(complex pa, complex pb, complex pc, complex uv) {
        return orientation(pa, pb, uv) > ACCURACY &&
               orientation(pb, pc, uv) > ACCURACY &&
               orientation(pc, pa, uv) > ACCURACY;
    }

    inline bool is_inside_triangle(const Face f, const VecXc &cfn, const complex uv) {
        auto uv1 = cfn(f.id * 3 + 0);
        auto uv2 = cfn(f.id * 3 + 1);
        auto uv3 = cfn(f.id * 3 + 2);
        return is_inside_triangle(uv1, uv2, uv3, uv);
    }

    inline complex nearby_grid(const complex a) {
        return {std::round(a.real()), std::round(a.imag())};
    }

    inline double nearby_grid(double x, double dir) {
        if (dir == 0) return x;
        double f = abs(fmod(x, 1.));
        if (f < ACCURACY || 1 - f < ACCURACY) x += dir * ACCURACY;
        x = dir > 0 ? std::ceil(x) : std::floor(x);
        return x;
    }

    inline complex nearby_grid(complex uv, complex dir) {
        return {
            nearby_grid(uv.real(), dir.real()),
            nearby_grid(uv.imag(), dir.imag())
        };
    }
}

#endif
