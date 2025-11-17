#ifndef METRIKO_QEX_COMMON_H
#define METRIKO_QEX_COMMON_H
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/common/utilities.h"

namespace metriko::qex {
class Qvert {
public:
    const complex uv; //
    const Row3d pos;  //
    const int sid;    // simplex id (either vid, eid or fid)

    Qvert(
        const complex uv,
        const Row3d &pos,
        const int sid
    ) : uv(uv), pos(pos), sid(sid) { }
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
    Qedge(const Qport &p1, const Qport &p2) : port1(p1), port2(p2) { }
};

class Qhalf {
public:
    const Qedge &qe; //
    const int idx;   //
    const bool cano; // cannonical flag
    const Qport &port1() const { return cano ? qe.port1 : qe.port2; }
    const Qport &port2() const { return cano ? qe.port2 : qe.port1; }
    Qhalf(const Qedge &qe, const int idx, const bool cano) : qe(qe), idx(idx), cano(cano) { }
};

class Qface {
public:
    std::vector<Qhalf> qhalfs;
};

inline void compute_trs_matrix(
    const Hmesh &hm,
    const VecXc &cf,
    const VecXi &matching,
    const int rosyN,
    VecXc &heR,
    VecXc &heT
) {
    heR.resize(hm.nH);
    heT.resize(hm.nH);
    for (Half h: hm.halfs) {
        Crnr c1 = h.next().crnr();
        Crnr c2 = h.twin().prev().crnr();
        auto uv1 = cf(c1.id);
        auto uv2 = cf(c2.id);
        int m = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
        m = m < 0 ? (rosyN + m % rosyN) % rosyN : m % rosyN;
        auto r = get_quater_rot(m);
        auto t = uv2 - r * uv1;
        heR(h.id) = r;
        heT(h.id) = complex(std::round(t.real()), std::round(t.imag()));
    }
}

inline complex nearby_grid(const complex a) {
    return {std::round(a.real()), std::round(a.imag())};
}

inline double nearby_grid(double x, double dir) {
    if (dir == 0) return x;
    double f = abs(fmod(x, 1.));
    if (f < EPS || 1 - f < EPS) x += dir * EPS;
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
