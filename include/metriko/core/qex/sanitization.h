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
        VecXi heMatching;
        MatXd heTranslation;
        compute_trs_matrix(mesh, cfn, matching, rosyN, heMatching, heTranslation);
        for (Vert v: mesh.verts) {
            //double max = 0;
            //for (Half h: v.adjHalfs()) {
            //    Crnr c = h.next().crnr();
            //    Row2d uv = cfn.row(c.id);
            //    max = std::max(max, uv.norm()); //todo: check whether standardn norm is fine
            //}

            int counter = 0;
            for (Half h: v.adjHalfs()) {
                if(counter == 0) {
                    if (singular[v.id]) {
                        fix_singular_point();
                        Crnr cCurr = h.next().crnr();
                        Row2d uv = cfn.row(cCurr.id);
                        cfn.row(cCurr.id) = Row2d(std::round(uv.x()), std::round(uv.y()));
                    } else {
                        //Crnr cCurr = h.next().crnr();
                        //Row2d uv = cfn.row(cCurr.id);
                        //double delta = std::pow(2, log2(max));
                        //Vec2d sign(uv.x() > 0 ? 1 : -1, uv.y() > 0 ? 1 : -1);
                        //cfn.row(cCurr.id) = cfn.row(cCurr.id) + (delta * sign).transpose() - (delta * sign).transpose(); //todo: does not make sense. better checking code
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
    }
}

#endif
