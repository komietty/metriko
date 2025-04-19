//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_GEN_QGP_EDGE_H
#define METRIKO_GEN_QGP_EDGE_H
#include "gen_qgp_vert.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko::visualizer {
    struct SplitEdge {
        Face face;
        complex uv1;
        complex uv2;
    };

    struct SplitArc {
        std::vector<SplitEdge> edges;
    };

    inline void compute_numeric_offsets(
        const Tmesh& tm,
        const Tquad& tq,
        const VecXd& R,
        const int bgnSide,
        const SplitElem& split_v_fr,
        const SplitElem& split_v_to,
        double& diff_bgn_side,
        double& diff_mid_side,
        double& diff_end_side
    ) {
        int midSide = (bgnSide + 1) % 4;
        int endSide = (bgnSide + 2) % 4;
        double dif_mid = 0;
        double dif_end = 0;

        for (int i = 0; i < tq.thids.size(); i++) {
            auto thid = tq.thids[i];
            auto side = tq.sides[i];
            if (side == midSide) dif_mid += R[tm.thalfs[thid].edge().id];
            if (side == endSide) dif_end += R[tm.thalfs[thid].edge().id];
        }
        diff_mid_side = dif_mid;
        diff_bgn_side = split_v_fr.distance;
        diff_end_side = dif_end - split_v_to.distance;
    }

    inline SplitArc gen_split_arc(
        const Tmesh& tm,
        const Tquad& tq,
        const VecXd& R,
        const VecXc& cfn_C,
        const VecXi& cmbf_matching,
        const int bgnSide,
        const SplitElem& split_v_fr,
        const SplitElem& split_v_to
    ) {
        std::vector<SplitEdge> result;

        // ----- if in the same face, just connect ----- //
        if (split_v_fr.segment.face.id == split_v_to.segment.face.id) {
            result.emplace_back(SplitEdge{split_v_fr.segment.face, split_v_fr.uv, split_v_to.uv});
            return SplitArc{ result };
        }

        double diff_bgn_side = 0;
        double diff_mid_side = 0;
        double diff_end_side = 0;

        compute_numeric_offsets(
            tm, tq, R, bgnSide, split_v_fr, split_v_to,
            diff_bgn_side,
            diff_mid_side,
            diff_end_side
        );

        double arg = -atan((diff_end_side - diff_bgn_side) / diff_mid_side);
        complex dir;
        {
            auto& th = tm.thalfs[tq.find_first_thid(bgnSide)];
            auto sign = th.cannonical ? 1. : -1.;
            dir = split_v_fr.segment.diff() * sign * 1000.;
            dir *= std::polar(1., PI / 2. + arg); // todo: need to asset param space does not flip anywhere
        }

        // ----- compute first segment ----- //
        double ratio_a2b;
        double ratio_c2d;
        bool find = false;
        Half nH;
        for (auto h: ((Face)split_v_fr.segment.face).adjHalfs()) {

            bool res = find_extended_intersection(
                split_v_fr.uv,
                split_v_fr.uv + dir,
                cfn_C(h.prev().crnr().id),
                cfn_C(h.next().crnr().id),
                ratio_a2b,
                ratio_c2d
            );

            if (res && ratio_a2b >= 0 && ratio_a2b <= 1 && ratio_c2d >= 0 && ratio_c2d <= 1) {
                nH = h;
                find = true;
                auto uv3 = lerp(
                    cfn_C(h.prev().crnr().id),
                    cfn_C(h.next().crnr().id),
                    ratio_c2d
                );
                //todo?? check whether the opposite side is in the same face.
                result.emplace_back(SplitEdge{nH.face(), split_v_fr.uv, uv3});
                break;
            }
        }
        assert(find);

        // ----- compute latter segments ----- //
        for (int i = 0; i < 200; i++) {
            nH = nH.twin();
            auto uv0 = lerp(
                cfn_C(nH.next().crnr().id),
                cfn_C(nH.prev().crnr().id),
                ratio_c2d);
            if (int m = cmbf_matching[nH.edge().id]; m != 0)
                dir *= std::polar(1., PI / 2 * (nH.isCanonical() ? 1 : -1) * m);
            nH = get_oppsite_half(cfn_C, uv0, dir, nH);
            auto uv1_ = cfn_C(nH.prev().crnr().id);
            auto uv2_ = cfn_C(nH.next().crnr().id);
            find_extended_intersection(uv0, uv0 + dir, uv1_, uv2_, ratio_a2b, ratio_c2d);
            auto uv3 = lerp(uv1_, uv2_, ratio_c2d);

            if (nH.face().id == split_v_to.segment.face.id) {
                result.emplace_back(SplitEdge{nH.face(), uv0, split_v_to.uv});
                break;
            }
            result.emplace_back(SplitEdge{nH.face(), uv0, uv3});
        }
        return SplitArc{ result };
    }

}

#endif
