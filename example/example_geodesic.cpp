#include <igl/readOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/geodesic/heat_method.h"
#include "metriko/core/hmesh/hmloc.h"
#include "polyscope/curve_network.h"

using namespace metriko;

/**
 * @brief Computes the normalized negative gradient field from a scalar field (like phi).
 */
VecXc compute_descent_field(const Hmesh& hm, const VecXd& phi) {
    VecXc G = VecXc::Zero(hm.nF);

    for (Face f : hm.faces) {
        Row3d g = Row3d::Zero();
        Row3d n = f.normal();
        for (Half h : f.adjHalfs()) {
            g += n.cross(h.vec()) * phi(h.crnr().vert().id) / (2.0 * f.area());
        }
        Row3d v = -g;
        if (v.norm() > 1e-12) {
            Row3d vn = v.normalized();
            G(f.id) = complex(vn.dot(f.basisX()), vn.dot(f.basisY()));
        }
    }
    return G;
}

VecXc compute_crnr_coord_inface(const Hmesh& hm) {
    VecXc xy(hm.nC);
    for (Crnr c: hm.crnrs) {
        Face  f = c.face();
        Row3d o = f.half().tail().pos();;
        Row3d v = c.vert().pos() - o;
        xy(c.id) = complex(v.dot(f.basisX()), v.dot(f.basisY()));
    }
    return xy;
}

HmLoc compute_opposite_loc(
    const Hmesh& hm,
    const VecXc& xycf, // #nC by 1: face interior coord xy for each crnr
    const VecXc& grad, // #nF by 1: gradient vector on each face coord
    const HmLoc& loc
) {
    return std::visit(overloaded{
        [&](const HmLocOnH& l) -> HmLoc {
            Half h_ = hm.halfs[l.id].twin();
            complex o0 = xycf(h_.next().crnr().id);
            complex o1 = xycf(h_.prev().crnr().id);
            complex a0 = lerp(o0, o1, 1 - l.r);
            complex a1 = a0 + grad(h_.face().id) * 100.;
            for (Half h: {h_.next(), h_.prev()}) {
                complex b0 = xycf(h.next().crnr().id);
                complex b1 = xycf(h.prev().crnr().id);
                double r_ab, r_cd;
                if (find_strict_intersection(a0, a1, b0, b1, r_ab, r_cd)) {
                    std::cout << "hid: " << h.id << std::endl;
                    return HmLocOnH(h.id, r_cd);
                }
            }
            throw std::runtime_error("");
        },
        [&](const HmLocOnF& l) -> HmLoc {
            Face f = hm.faces[l.id];
            for (Half h: f.adjHalfs()) {
                complex a0 = l.xy;
                complex a1 = l.xy + grad(l.id) * 100.;
                complex b0 = xycf(h.next().crnr().id);
                complex b1 = xycf(h.prev().crnr().id);
                double r_ab, r_cd;
                if (find_strict_intersection(a0, a1, b0, b1, r_ab, r_cd)) {
                    std::cout << "hid: " << h.id << std::endl;
                    return HmLocOnH(h.id, r_cd);
                }
            }
            throw std::runtime_error("");
        },
        [&](const auto& l) -> HmLoc { throw std::runtime_error("not implemented"); },
    }, loc);
}

vec<HmLoc> compute_geodesic(
    const Hmesh& hm,
    const VecXc& xycf, // #nC by 1: face interior coord xy for each crnr
    const VecXc& grad, // #nF by 1: gradient vector on each face coord
    const HmLoc& loc0,
    const HmLoc& loc1
) {
    vec path = { loc0 };

    std::visit(overloaded{
        [&](const HmLocOnF& f0, const HmLocOnV& v1) {
            Face  f = hm.faces[f0.id];
            Vert  v = hm.verts[v1.id];
            HmLoc l = loc0;
            for (int i = 0; i < 7; ++i) {
                l = compute_opposite_loc(hm, xycf, grad, l);
                path.push_back(l);
            }
        },
        [&](const auto& a, const auto& b) {
            throw std::runtime_error("This combination is not supported yet.");
        }
    }, loc0, loc1);

    return path;
}



int main(int argc, char** argv) {
    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    SprsD C;
    igl::cotmatrix(V, F, C);
    SprsD L = -C;
    for (int i = 0; i < hm.nV; ++i) L.coeffRef(i, i) += 1e-10;

    SprsD M; igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    VecXd delta = VecXd::Zero(hm.nV);

    //int target_fid = 300;
    //Row3d bary(1./3., 1./3., 1./3.);
    //delta(hm.idx(target_fid, 0)) = bary.x();
    //delta(hm.idx(target_fid, 1)) = bary.y();
    //delta(hm.idx(target_fid, 2)) = bary.z();
    delta[0] = 1.;

    VecXd phi  = compute_heatflow(hm, L, M, delta);
    VecXc grad = compute_descent_field(hm, phi);
    VecXc xycf = compute_crnr_coord_inface(hm);

    // compute geodesic
    Face f = hm.faces[10];
    complex c0 = xycf(f.half().crnr().id);
    complex c1 = xycf(f.half().next().crnr().id);
    complex c2 = xycf(f.half().prev().crnr().id);
    complex center = c0 * 0.5 + c1 * 0.25 + c2 * 0.25;
    vec path = compute_geodesic(hm, xycf, grad, HmLocOnF(f.id, center), HmLocOnV(0));

    // 3. Visualization
    MatXd G(hm.nF, 3);
    for (Face f : hm.faces) {
        complex c = grad(f.id);
        G.row(f.id) = c.real() * f.basisX() + c.imag() * f.basisY();
    }

    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    auto* psMesh = polyscope::registerSurfaceMesh("mesh", V, F);
    auto* q_dst  = psMesh->addVertexScalarQuantity("geodesic distance", phi);
    auto* q_phi  = psMesh->addFaceVectorQuantity("descent direction (X_phi)", G);
    q_dst->setEnabled(true);
    q_dst->setIsolineWidth(0.05, false);
    q_dst->setColorMap("spectral");
    q_phi->setEnabled(true);
    q_phi->setVectorLengthScale(0.01, false);
    q_phi->setVectorColor({1.0, 0.5, 0.0}); // Orange

    MatXd P(path.size(), 3);
    for (size_t i = 0; i < path.size(); ++i) {
        P.row(i) = get_ptloc_pos(hm, path[i]);   // hpath.h の変換関数
    }
    auto* psPath = polyscope::registerCurveNetworkLine("geodesic path", P);
    psPath->setColor({1.0, 0.0, 0.0});   // 赤
    psPath->setRadius(0.003);


    polyscope::show();
    return 0;
}
