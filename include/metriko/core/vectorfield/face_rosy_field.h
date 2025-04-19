//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_FACE_ROSY_FIELD_H
#define METRIKO_FACE_ROSY_FIELD_H
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

    std::unique_ptr<FaceRosyField> compute_combbed_field(
        const FaceRosyField &field,
        const std::vector<bool> &seam
    );

    std::vector<bool> compute_seam(const FaceRosyField &f);
}

#endif
