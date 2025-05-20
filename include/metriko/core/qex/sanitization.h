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
        MatXd& cfn
    ) {
        MatXd cfn_ = cfn;

        VecXi heMatching;
        MatXd heTranslation;
        compute_trs_matrix(mesh, cfn, matching, rosyN, heMatching, heTranslation);
        for (Vert v: mesh.verts) {
            double max = 0;
            for (Half h: v.adjHalfs()) {
                Crnr c = h.next().crnr();
                Row2d uv = cfn.row(c.id);
                max = std::max(max, abs(uv.x()));
                max = std::max(max, abs(uv.y()));
            }

            int counter = 0;
            for (Half h: v.adjHalfs()) {
                if(counter == 0) {
                    if (singular[v.id]) {
                        fix_singular_point();
                        Crnr cCurr = h.next().crnr();
                        Row2d uv = cfn.row(cCurr.id);
                        cfn.row(cCurr.id) = Row2d(std::round(uv.x()), std::round(uv.y()));
                    } else {
                        Crnr cCurr = h.next().crnr();
                        Row2d uv = cfn.row(cCurr.id);
                        double delta = std::pow(2, log2(max));
                        Row2d sign((0 < uv.x()) - (uv.x() < 0), (0 < uv.y()) - (uv.y() < 0));
                        cfn.row(cCurr.id) = (cfn.row(cCurr.id) + delta * sign) - delta * sign;
                    }
                } else {
                    Crnr cPrev = h.twin().prev().crnr();
                    Crnr cCurr = h.next().crnr();
                    int m = heMatching[h.twin().id];
                    Row2d t = heTranslation.row(h.twin().id);
                    cfn.row(cCurr.id) = matching2rot(m) * cfn.row(cPrev.id).transpose() + t.transpose();
                }
                counter++;
            }
        }
        std::cout << "sanitization diff: " << (cfn - cfn_).norm() << std::endl;
    }
}

#endif
