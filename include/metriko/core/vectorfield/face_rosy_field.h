//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_FACE_ROSY_FIELD_H
#define METRIKO_FACE_ROSY_FIELD_H
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>
#include <igl/cut_mesh_from_singularities.h>
#include "base_field.h"
#include "../solver/sparse_linear_solver.h"

namespace metriko {
    class FaceRosyField : public BaseVectorField {
    public:
        FaceRosyField(const Hmesh &m, const int nRosy): BaseVectorField(m, nRosy) { }
        FaceRosyField(const Hmesh &m, const int nRosy, FieldType type): BaseVectorField(m, nRosy) {
            SprsC L = connectionLaplacian();
            SprsC M = galerkinMassMatrix();
            switch (type) {
                case FieldType::Smoothest: {
                    compressed = solveSmallestEig(L, M);
                    break;
                }
                case FieldType::CurvatureAligned: {
                    assert(rosyN == 2 || rosyN == 4);
                    constexpr double lambda = 0;
                    VecXc D = principalCurvatureDir();
                    if (rosyN == 4) D = D.array().square();
                    VecXc rhs = M * D / sqrt(abs((D.adjoint() * M * D)[0]));
                    SprsC lhs = L - lambda * M;
                    compressed = solveSquare(lhs, rhs);
                    break;
                }
                default: { break; }
            }
            auto r = std::polar(1., TwoPI / rosyN);
            field.resize(mesh.nF, rosyN);
            for (auto f: mesh.faces) {
                for (int i = 0; i < rosyN; i++) {
                    field(f.id, i) = pow(std::polar(1., std::arg(compressed(f.id))), 1. / rosyN) * pow(r, i);
                }
            }
        }

        SprsC connectionLaplacian() const {
            SprsC S(mesh.nF, mesh.nF);
            SprsC I(mesh.nF, mesh.nF);
            VecXc transport(mesh.nH);
            std::vector<TripC> T;

            for (Edge e: mesh.edges) {
                Half h1 = e.half();
                Half h2 = e.half().twin();
                auto v1 = std::polar(1., h1.farg());
                auto v2 = std::polar(1., h2.farg());
                auto r = -v2 / v1;
                transport[h1.id] = r;
                transport[h2.id] = complex(1, 0) / r;
            }

            for (Face f: mesh.faces) {
                for (Half h: f.adjHalfs()) {
                    if (h.twin().isBoundary()) continue;
                    auto w = 1.;
                    auto r = transport[h.twin().id];
                    Face twinF = h.twin().face();
                    T.emplace_back(f.id, f.id, w);
                    T.emplace_back(f.id, twinF.id, -w * pow(r, rosyN));
                }
            }

            S.setFromTriplets(T.begin(), T.end());
            I.setIdentity();
            return S + 1e-9 * I;
        }

        SprsC galerkinMassMatrix() const {
            SprsC S(mesh.nF, mesh.nF);
            std::vector<TripC> T;
            for (Face f: mesh.faces) { T.emplace_back(f.id, f.id, f.area()); }
            S.setFromTriplets(T.begin(), T.end());
            return S;
        }

        VecXc principalCurvatureDir() const {
            VecXc D(mesh.nF);
            for (Face f: mesh.faces) {
                complex dir{0, 0};
                for (Half h: f.adjHalfs()) {
                    auto l = h.len();
                    auto c = std::polar(1., h.farg()) * l;
                    dir += -c * c / l * h.darg();
                }
                D[f.id] = dir * 0.25;
            }
            return D;
        }

        void computeSingular(const VecXd &effort) {
            singular = VecXi::Ones(mesh.nV);
            for (auto v: mesh.verts | vw::filter([](auto _v) { return !_v.isBoundary(); })) {
                double sum = mesh.angleDefect[v.id] * rosyN;
                for (auto h: v.adjHalfs()) sum += (h.isCanonical() ? -1 : 1) * effort[h.edge().id];
                singular[v.id] = static_cast<int>(std::round(sum / TwoPI));
            }
        }

        void computeConnection() {
            connection = VecXd::Zero(mesh.nE);
            for (auto e: mesh.edges | vw::filter([](auto _e) { return !_e.isBoundary(); })) {
                Vec3d v = e.half().vec().normalized();
                complex ef(v.dot(e.face0().basisX()), v.dot(e.face0().basisY()));
                complex eg(v.dot(e.face1().basisX()), v.dot(e.face1().basisY()));
                connection(e.id) = eg / ef;
            }
        }

        void computeMatching(const MatchingType type) {
            computeConnection();
            matching = VecXi::Constant(mesh.nE, -1);
            VecXd effort = VecXd::Zero(mesh.nE);
            switch (type) {
                case MatchingType::Principal: {
                    pMatching(effort);
                    break;
                }
                case MatchingType::Curl: {
                    cMatching(effort);
                    break;
                }
                default: { break; }
            }
            computeSingular(effort);
        }

    private:
        void pMatching(VecXd &effort) {
            for (auto e: mesh.edges | vw::filter([](auto e_) { return !e_.isBoundary(); })) {
                auto conn = connection[e.id];
                auto minRot = 1000.;
                auto offset = 0;
                complex coef(1, 0);
                complex transport0 = conn * field(e.face0().id, 0);
                for (int i = 0; i < rosyN; i++) {
                    complex jf = field(e.face0().id, i);
                    complex jg = field(e.face1().id, i);
                    coef *= jg / (conn * jf);
                    double r = std::arg(jg / transport0);
                    if (abs(r) < abs(minRot)) {
                        offset = i;
                        minRot = r;
                    }
                }
                effort[e.id] = std::arg(coef);

                double residuals = 0;
                for (int i = 0; i < rosyN; i++) {
                    auto v0 = field(e.face0().id, i);
                    auto v1 = field(e.face1().id, (i + offset) % rosyN);
                    residuals += std::arg(v1 / (conn * v0));
                }
                matching[e.id] = offset - std::round((residuals - effort[e.id]) / TwoPI);
                //std::cout << residuals - effort[e.id] << std::endl;
                //matching[e.id] = offset;
            }
        }

        void cMatching(VecXd &effort) {
            for (auto e: mesh.edges | vw::filter([](auto e_) { return !e_.isBoundary(); })) {
                int iM = 0;
                double min = 32767000.;
                Vec3d v = e.half().vec().normalized();
                Face f0 = e.face0();
                Face f1 = e.face1();

                for (int i = 0; i < rosyN; i++) {
                    double curl = 0;
                    for (int j = 0; j < rosyN; j++) {
                        complex c0 = field(f0.id, j);
                        complex c1 = field(f1.id, (i + j) % rosyN);
                        Row3d v0 = c0.real() * f0.basisX() + c0.imag() * f0.basisY();
                        Row3d v1 = c1.real() * f1.basisX() + c1.imag() * f1.basisY();
                        curl += pow(v.dot(v1 - v0), 2.);
                    }
                    if (curl < min) {
                        iM = i;
                        min = curl;
                    }
                }

                matching[e.id] = iM;

                complex coef(1, 0);
                for (int i = 0; i < rosyN; i++) {
                    auto v0 = field(e.face0().id, i);
                    auto v1 = field(e.face1().id, (i + iM) % rosyN);
                    coef *= v1 / (v0 * connection[e.id]);
                }
                effort[e.id] = arg(coef);
            }
        }
    };
}

namespace metriko {
    inline std::unique_ptr<FaceRosyField> compute_combbed_field(
            const FaceRosyField& field,
            const std::vector<bool>& seam
    ) {
        auto N = field.rosyN;
        auto F = std::make_unique<FaceRosyField>(field.mesh, field.rosyN);
        F->field.resize(field.mesh.nF, N);
        F->connection = field.connection;
        F->compressed = field.compressed;
        F->singular = field.singular;
        VecXi turns = VecXi::Zero(F->mesh.nF);
        VecXi visit = VecXi::Zero(F->mesh.nF);
        std::queue<std::pair<int, int>> matchingQ;
        matchingQ.emplace(0, 0);

        do {
            std::pair<int, int> pop = matchingQ.front();
            matchingQ.pop();
            auto [i, m] = pop;
            int r = N - m;
            if (visit(i) == 1) continue;
            visit(i) = 1;
            turns(i) = m;
            F->field.block(i, 0, 1, r) = field.field.block(i, m, 1, r);
            F->field.block(i, r, 1, m) = field.field.block(i, 0, 1, m);
            Face face = F->mesh.faces[i];

            for (Half h: face.adjHalfs()) {
                Edge e = h.edge();
                bool b = e.face0().id == i;
                Face fn = b ? e.face1() : e.face0();
                int n = (b ? 1 : -1) * field.matching[e.id];
                n = (n + m + 1000 * N) % N;
                if (!e.isBoundary() && !visit(fn.id) && !seam[e.id]) matchingQ.emplace(fn.id, n);
            }
        } while (!matchingQ.empty());

        F->matching = VecXi::Constant(F->mesh.nE, -1);
        for (Edge e: F->mesh.edges | vw::filter([](auto e_) { return !e_.isBoundary();}))
            F->matching[e.id] = (turns[e.face0().id] - turns[e.face1().id] + field.matching[e.id] + 1000 * N) % N;

        return F;
    }

    inline MatXd compute_extrinsic_field(const Hmesh& mesh, const FaceRosyField& field, int rosyN) {
        MatXd ext(mesh.nF, 3 * rosyN);
        for (Face f: mesh.faces) {
            complex c0 = field.field(f.id, 0);
            complex c1 = field.field(f.id, 1);
            complex c2 = field.field(f.id, 2);
            complex c3 = field.field(f.id, 3);
            ext.block(f.id, 0, 1, 3) = (c0.real() * f.basisX() + c0.imag() * f.basisY()).normalized();
            ext.block(f.id, 3, 1, 3) = (c1.real() * f.basisX() + c1.imag() * f.basisY()).normalized();
            ext.block(f.id, 6, 1, 3) = (c2.real() * f.basisX() + c2.imag() * f.basisY()).normalized();
            ext.block(f.id, 9, 1, 3) = (c3.real() * f.basisX() + c3.imag() * f.basisY()).normalized();
        }
        return ext;
    }



    inline void cut_mesh_with_singularities(
        const MatXd &V,
        const MatXi &F,
        const std::vector<std::vector<int>> &VF,
        const std::vector<std::vector<int>> &VV,
        const MatXi &TT,
        const MatXi &TTi,
        const VecXi &singularities,
        Eigen::MatrixXi &cuts
    ) {
        //first, get a spanning tree for the mesh (no missmatch needed)
        igl::cut_mesh_from_singularities(V, F, MatXd::Zero(F.rows(), 3).eval(), cuts);

        std::set<int> vertices_in_cut;
        for (int i = 0; i < cuts.rows(); ++i) {
        for (int j = 0; j < cuts.cols(); ++j) {
            if (cuts(i, j)) vertices_in_cut.insert(F(i, j));
        }}

        //then, add all singularities one by using Dijkstra's algorithm
        for (int i = 0; i < singularities.rows(); ++i) {
            std::vector<int> path;
            VecXd min_distance;
            VecXi previous;
            int vertex_found = igl::dijkstra(singularities[i], vertices_in_cut, VV, min_distance, previous);
            if (vertex_found == -1)
                // this means that there are no cuts
                path.push_back(singularities[i]);
            else
                igl::dijkstra(vertex_found, previous, path);

            vertices_in_cut.insert(path.begin(), path.end());

            //insert to cut
            for (int ii = 0; ii < path.size() - 1; ++ii) {
                const int &v0 = path[ii];
                const int &v1 = path[ii + 1];

                std::vector<int> vf0 = VF[v0];
                rg::sort(vf0.begin(), vf0.end());
                std::vector<int> vf1 = VF[v1];
                rg::sort(vf1.begin(), vf1.end());

                std::vector<int> common_face_v(std::max(vf0.size(), vf1.size()));
                std::vector<int>::iterator it;
                it = rg::set_intersection(vf0, vf1, common_face_v.begin()).out;
                common_face_v.resize(it - common_face_v.begin());
                assert(common_face_v.size() == 2);

                const int &fi = common_face_v[0];
                int j = -1;
                for (unsigned z = 0; z < 3; ++z)
                    if ((F(fi, z) == v0 && F(fi, (z + 1) % 3) == v1) ||
                        (F(fi, z) == v1 && F(fi, (z + 1) % 3) == v0)) { j = z; }
                assert(j != -1);
                cuts(fi, j) = 1;
                cuts(TT(fi, j), TTi(fi, j)) = 1;
            }
        }
    }

    //Wrapper of the above with only vertices and faces as mesh input
    inline void cut_mesh_with_singularities(
        const MatXd &V,
        const MatXi &F,
        const VecXi &singularities,
        MatXi &cuts
    ) {
        std::vector<std::vector<int> > VF, VFi;
        std::vector<std::vector<int> > VV;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        igl::adjacency_list(F, VV);
        MatXi TT, TTi;
        igl::triangle_triangle_adjacency(F, TT, TTi);
        cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities, cuts);
    }

    inline std::vector<bool> compute_seam(const FaceRosyField &f) {
        MatXi face2cut;
        std::vector seam(f.mesh.nE, false);
        std::vector<int> svids;
        for (Vert v: f.mesh.verts)
            if (f.singular[v.id] != 0)
                svids.emplace_back(v.id);
        VecXi s = Eigen::Map<VecXi>(svids.data(), svids.size());
        cut_mesh_with_singularities(f.mesh.pos, f.mesh.idx, s, face2cut);
        for (int iF = 0; iF < f.mesh.nF; iF++) {
            for (int j = 0; j < 3; j++) {
                if (face2cut(iF, j)) seam[f.mesh.face2edge(iF, j)] = true;
            }
        }
        return seam;
    }
}

#endif
