//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_ITER_ROUNDING_LOOP_H
#define METRIKO_ITER_ROUNDING_LOOP_H
#include "iter_rounding_common.h"

namespace metriko {
    template <typename T>
class SeamlessIntegration: public Integration {
public:
    const SprsD& UFull;
    const double frac;
    VecXd xInit;
    VecXd xCurr;
    VecXd xPrev;
    std::vector<int> leftIdcs;
    std::vector<int> fixedIdcs;
    std::vector<int> fixedVals;
    std::vector<int> integerIdcs;  // used when rounding uv on seams
    std::vector<int> singularIdcs; // used when rounding uv on singularities
    bool roundedSingulars, roundSeams;

    void initial_solution(VecXd &x0_) const { x0_ = xInit; }

    static bool post_iteration(const VecXd &x) { return false; }

    bool initFixedIndices() {
        xPrev = xCurr;
        VecXd X = UFull * xCurr;
        double minDif = 3276700.0;
        int minIdx = -1;
        for (int i = 0; i < leftIdcs.size(); i++) {
            double v = frac * X(leftIdcs[i]);
            double d = std::fabs(v - std::round(v));
            if (d < minDif) { minIdx = i; minDif = d; }
        }

        fixedIdcs.emplace_back(leftIdcs[minIdx]);
        fixedVals.emplace_back(std::round(X(leftIdcs[minIdx])));
        leftIdcs.erase(leftIdcs.begin() + minIdx);

        // completed rounding singularities; starting to round rest of seams
        if (leftIdcs.empty() && !roundSeams && !roundedSingulars) {
            leftIdcs = integerIdcs;
            roundedSingulars = true;
        }

        //Poisson error
        gInteg = X2F * wInteg;

        //Closeness
        igl::speye(xCurr.size(), gClose);
        gClose = gClose * wClose;

        //fixedIndices constantness
        gConst.resize(fixedIdcs.size(), X.size());
        std::vector<TripD> gcT;
        for (int i = 0; i < fixedIdcs.size(); i++)
            gcT.emplace_back(i, fixedIdcs[i], 1.);

        gConst.setFromTriplets(gcT.begin(), gcT.end());
        gConst = gConst * UFull * wConst;

        if (ESize == 0) { // why needed this inside code
            SprsD J;
            VecXd E;
            objective_jacobian(VecXd::Random(UFull.cols()), E, J, false);
            ESize = E.size();
        }

        return minDif > 10e-7;
    }

    // memo: J * E is rhs of LMSolver
    void objective_jacobian(const VecXd &xCurrSmall, VecXd &E, SprsD &J, bool updateJ) {
        const VecXd X = UFull * xCurrSmall; // here UFull
        const VecXd field = X2F * xCurrSmall;
        const VecXd fInteg = field - F2;
        const VecXd fClose = xCurrSmall - xPrev;
        VecXd fConst(fixedIdcs.size());
        for (int i = 0; i < fixedIdcs.size(); i++)
            fConst(i) = X(fixedIdcs[i]) - fixedVals[i] * frac;
        compute_jacobian(field, fInteg, fClose, fConst, updateJ, E, J, true);
    }

    // assure that all values of round diff by fixed indices are less than 1e-7
    bool post_checking(const VecXd &x) {
        VecXd fullx = UFull * x;
        double max = -1;
        for (const int i : fixedIdcs) {
            const double v = frac * fullx(i);
            max = std::max(max, std::abs(v - std::round(v)));
        }
        xCurr = x;
        xInit = x;
        xPrev = x;
        return max <= 10e-7;
    }

    SeamlessIntegration(
            const T &ni,
            const VecXd& xInit_,
            const VecXd& F2_,
            const VecXi& iIdcs,
            const VecXi& sIdcs,
            const bool roundSeams_
    ): Integration(ni.G2, F2_, ni.N, ni.n, ni.localInjectivity, ni.length, 1., 10e4, 0.01, ni.wBarrier, ni.barrier),
       UFull(ni.UFull),
       xInit(xInit_),
       frac(1.),
       roundedSingulars(false),
       roundSeams(roundSeams_)
    {
        integerIdcs  = std::vector(iIdcs.data(), iIdcs.data() + iIdcs.size());
        singularIdcs = std::vector(sIdcs.data(), sIdcs.data() + sIdcs.size());
        leftIdcs = roundSeams ? integerIdcs : singularIdcs;

        xSize = UFull.cols();
        xCurr = xInit;
        xPrev = xInit;
        X2F = (G2 * UFull * length).pruned();
    }
};
}
#endif
