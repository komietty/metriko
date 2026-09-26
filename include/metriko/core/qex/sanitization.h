#ifndef METRIKO_QEX_SANITIZATION_H
#define METRIKO_QEX_SANITIZATION_H
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
#if METRIKO_DEBUG
        VecXc buk = cfn;
#endif
        VecXc heR;
        VecXc heT;
        compute_trs_matrix(mesh, cfn, matching, rosyN, heR, heT);
        for (Vert v: mesh.verts) {
            double max = 0;
            for (Half h: v.adjHalfs()) {
                auto cc = h.next().crnr();
                auto uv = cfn(cc.id);
                max = std::max(max, abs(uv.real()));
                max = std::max(max, abs(uv.imag()));
            }

            bool init = false;
            for (Half h: v.adjHalfs()) {
                if(!init) {
                    if (singular[v.id]) {
                        auto cc = h.next().crnr();
                        auto uv = cfn(cc.id);
                        cfn(cc.id) = complex(std::round(uv.real()), std::round(uv.imag()));
                    } else {
                        auto cc = h.next().crnr();
                        auto uv = cfn(cc.id);
                        double delta = std::pow(2, log2(max));
                        complex sign((0 < uv.real()) - (uv.real() < 0), (0 < uv.imag()) - (uv.imag() < 0));

                        //cfn(cc.id) = (cfn(cc.id) + delta * sign) - delta * sign;

                        uv = (uv + delta * sign) - delta * sign;
                        // a vertex within tol of a grid POINT is that point (both coordinates): otherwise the
                        // same point is found as an edge or face q-vertex in the neighbouring faces and never
                        // pairs. snapping a single coordinate is not done here: the tracer cannot yet pass
                        // through a vertex lying exactly on an isoline
                        // todo: 根本的には、pick_next_half に「半直線が頂点に当たったら、その頂点の周りで方向が向く面へ遷移する」処理を入れる必要があります（motorcycle の update_to_oppo が corner 通過で同じことをしています）。それが入れば座標ごとの吸着に戻せます。
                        constexpr double tol = 1e-4;
                        complex g = nearby_grid(uv);
                        if (std::abs(uv.real() - g.real()) < tol && std::abs(uv.imag() - g.imag()) < tol) uv = g;
                        cfn(cc.id) = uv;
                    }
                } else {
                    Crnr cp = h.twin().prev().crnr();
                    Crnr cc = h.next().crnr();
                    complex r = heR(h.twin().id);
                    complex t = heT(h.twin().id);
                    cfn(cc.id) = r * cfn(cp.id) + t;
                }
                init = true;
            }
        }

        #if METRIKO_DEBUG
        std::cout << "Sanitization norm diff: " << (cfn - buk).norm() << std::endl;
        #endif
    }
}
#endif
