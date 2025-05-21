#ifndef SANITIZE_H
#define SANITIZE_H
#include "common.h"

namespace metriko::qex {
    inline void fix_singular_point() { }

    inline void sanitization(
        const Hmesh& mesh,
        const VecXi& matching,
        const VecXi& singular,
        const int rosyN,
        VecXc& cfn
    ) {
        VecXc cfn_ = cfn;

        //*
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(mesh, cfn, matching, rosyN, heR, heT);
        for (Vert v: mesh.verts) {
            double max = 0;
            for (Half h: v.adjHalfs()) {
                Crnr c = h.next().crnr();
                complex uv = cfn(c.id);
                max = std::max(max, abs(uv.real()));
                max = std::max(max, abs(uv.imag()));
            }

            bool init = false;
            for (Half h: v.adjHalfs()) {
                if(!init) {
                    if (singular[v.id]) {
                        fix_singular_point();
                        Crnr cCurr = h.next().crnr();
                        auto uv = cfn(cCurr.id);
                        cfn(cCurr.id) = complex(std::round(uv.real()), std::round(uv.imag()));
                    } else {
                        Crnr cCurr = h.next().crnr();
                        auto uv = cfn(cCurr.id);
                        double delta = std::pow(2, log2(max));
                        complex sign((0 < uv.real()) - (uv.real() < 0), (0 < uv.imag()) - (uv.imag() < 0));
                        cfn(cCurr.id) = (cfn(cCurr.id) + delta * sign) - delta * sign;
                    }
                } else {
                    Crnr cPrev = h.twin().prev().crnr();
                    Crnr cCurr = h.next().crnr();
                    complex r = heR(h.twin().id);
                    complex t = heT(h.twin().id);
                    cfn(cCurr.id) = r * cfn(cPrev.id) + t;
                }
                init = true;
            }
        }
        //*/
        std::cout << "diff before: " << cfn_.norm() << std::endl;
        std::cout << "diff after: " << cfn.norm() << std::endl;
        std::cout << "sanitization diff: " << (cfn - cfn_).norm() << std::endl;
    }
}

#endif
