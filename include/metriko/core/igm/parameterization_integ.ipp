#pragma once
namespace metriko {
    inline bool RosyParameterization::integ() {
        SprsD fullx2Nfunc = vtrans2cut * uncompress;
        SprsD Cfull = constraint * uncompress;

        if (Cfull.rows() != 0) {
            // reducing the size of the constraint matrix
            Eigen::SparseQR<SprsD, Eigen::COLAMDOrdering<int> > qr;
            qr.compute(Cfull.transpose());
            int rank = qr.rank();
            VecXi idcs = qr.colsPermutation().indices(); //creating sliced permutation matrix

            std::vector<TripD> T;
            for (int k = 0; k < Cfull.outerSize(); ++k) {
            for (SprsD::InnerIterator it(Cfull, k); it; ++it) {
            for (int j = 0; j < rank; j++)
                if (it.row() == idcs(j)) T.emplace_back(j, it.col(), it.value());
            }}

            Cfull.resize(rank, Cfull.cols());
            Cfull.setFromTriplets(T.begin(), T.end());
        }

        //---- generating G: the matrix for generating a vector field from uv -----
        std::vector<TripD> eT;
        std::vector<TripD> iT;
        for (Face f: cut.faces) {
        for (Half h: f.adjHalfs()) {
            Row3d g = f.normal().cross(h.next().vec()) / (f.area() * 2);
            for (int k = 0; k < N; k++) {
                int c  = N * h.tail().id + k;
                int er = 3 * N * f.id + k * 3;
                int ir = 2 * N * f.id + k * 2;
                eT.emplace_back(er + 0, c, g.x());
                eT.emplace_back(er + 1, c, g.y());
                eT.emplace_back(er + 2, c, g.z());
                iT.emplace_back(ir + 0, c, g.dot(f.basisX()));
                iT.emplace_back(ir + 1, c, g.dot(f.basisY()));
            }
        }}
        SprsD G3(3 * N * cut.nF, N * cut.nV);
        SprsD G2(2 * N * cut.nF, N * cut.nV);
        G3.setFromTriplets(eT.begin(), eT.end());
        G2.setFromTriplets(iT.begin(), iT.end());

        VecXd rawF2Vec;
        double norm = 0;
        {
            norm = avg_norm(ext, 3, N, cut.nF);
            MatXd rawF3 = ext;
            rawF3.array() /= norm;

            MatXd rawF2(cut.nF, 2 * N);
            for (int i = 0; i < N; i++)
                rawF2.middleCols(2 * i, 2) <<
                        ext.middleCols(3 * i, 3).cwiseProduct(cut.faceBasisX).rowwise().sum(),
                        ext.middleCols(3 * i, 3).cwiseProduct(cut.faceBasisY).rowwise().sum();
            rawF2Vec = rawF2.reshaped<Eigen::RowMajor>().transpose();
        }

        //------ generating locally_injectivity_check -----
        auto locally_injectivity_check = [&](const VecXd &X) {
            if (!localInjectivity) return false;
            const VecXd NF = fullx2Nfunc * X;
            MatXd nfn_(cut.nV, N);
            for (int i = 0; i < nfn_.rows(); i++)
                nfn_.row(i) << NF.segment(N * i, N).transpose();
            for (Face f: cut.faces) {
                Half h = f.half();
                Row2d a = nfn_.block(h.tail().id, 0, 1, 2);
                Row2d b = nfn_.block(h.next().tail().id, 0, 1, 2);
                Row2d c = nfn_.block(h.prev().tail().id, 0, 1, 2);
                Row2d d1 = b - a;
                Row2d d2 = c - a;
                if (d1.x() * d2.y() - d1.y() * d2.x() < 0) return false;
            }
            if (verbose) std::cout << "now the map is locally injective" << std::endl;
            return true;
        };

        bool success = iterative_rounding(
            fixedIdcs,
            fixedVals,
            singularIdcs,
            integerIdcs,
            (cut.pos.colwise().maxCoeff() - cut.pos.colwise().minCoeff()).norm() * gridscale / norm,
            Cfull,
            G2 * fullx2Nfunc,
            N,
            n,
            cut.nF,
            seamless,
            roundSeams,
            localInjectivity,
            verbose,
            rawF2Vec,
            locally_injectivity_check,
            fullx
        );

        if (!success && verbose) std::cout << "Rounding has failed!" << std::endl;


        nfn.resize(cut.nV, N);
        cfn.resize(raw.nF, N * 3);

        VecXd NF = fullx2Nfunc * fullx;
        for (int i = 0; i < nfn.rows(); i++)
            nfn.row(i) << NF.segment(N * i, N).transpose();

        for (int i = 0; i < raw.nF; i++) {
        for (int j = 0; j < 3; j++) {
            cfn.block(i, N * j, 1, N) = nfn.row(cut.idx(i, j)).array();
        }}

        //----- check the derivative of the map is close to the original tangent field ----//
        double evaluation = 0;
        VecXd G_NF = G3 * NF;
        MatXd vf(ext.rows(), ext.cols());
        for (int i = 0; i < ext.rows(); i++) {
            vf.row(i) << G_NF.segment(3 * N * i, 3 * N).transpose();
        }
        for (int i = 0; i < ext.rows(); i++) {
            for (int j = 0; j < N; j++) {
                Row3d v1 = ext.row(i).segment(j * 3, 3).normalized();
                Row3d v2 = vf.row(i).segment(j * 3, 3).normalized();
                evaluation += (v1 - v2).norm();
            }
        }
        std::cout <<
                "Deviation of recovered tangent field from original tangent field: "
                << evaluation / (ext.rows() * N) << std::endl;

        return success;
    }
}
