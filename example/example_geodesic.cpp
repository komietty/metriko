#include <igl/readOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <polyscope/surface_mesh.h>
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/geodesic/heat_method.h"

using namespace metriko;

/**
 * @brief Computes the negative gradient (descent direction) on each face.
 */
MatXd compute_all_face_gradients(const Hmesh& hm, const VecXd& phi) {
    MatXd G(hm.nF, 3);
    for (Face f : hm.faces) {
        Row3d grad = Row3d::Zero();
        Row3d n = f.normal();
        for (Half h : f.adjHalfs()) {
            // Gradient formula: sum phi_i * (n x e_i) / 2A
            // h.crnr().vert() is the vertex opposite to edge h.vec()
            grad += n.cross(h.vec()) * phi(h.crnr().vert().id) / (2.0 * f.area());
        }
        Row3d v = -grad;
        if (v.norm() > 1e-12) {
            G.row(f.id) = v.normalized();
        } else {
            G.row(f.id).setZero();
        }
    }
    return G;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: ./metriko_geodesic <mesh.obj>" << std::endl;
        return 1;
    }

    MatXd V; MatXi F;
    if (!igl::readOBJ(argv[1], V, F)) return 1;

    Hmesh hm(V, F);

    // Setup operators
    SprsD Cot; igl::cotmatrix(V, F, Cot);
    SprsD L = -Cot;
    for (int i = 0; i < hm.nV; ++i) L.coeffRef(i, i) += 1e-10;

    SprsD M; igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    VecXd delta = VecXd::Zero(hm.nV);
    delta(0) = 1.0; 

    // 1. Compute Geodesic Distance Field
    VecXd phi = compute_heatflow(hm, L, M, delta);

    // 2. Compute Descent Vector Field (Face-based)
    MatXd gradients = compute_all_face_gradients(hm, phi);

    // 3. Visualization
    polyscope::init();
    auto* psMesh = polyscope::registerSurfaceMesh("mesh", V, F);
    
    // Distance field with isolines
    auto* q_dist = psMesh->addVertexScalarQuantity("geodesic distance", phi);
    q_dist->setEnabled(true);
    q_dist->setIsolineWidth(0.05, false);
    q_dist->setColorMap("spectral");

    // Vector field (Descent directions)
    auto* q_vec = psMesh->addFaceVectorQuantity("descent direction", gradients);
    q_vec->setEnabled(true);
    q_vec->setVectorLengthScale(0.01, false);
    q_vec->setVectorRadius(0.001, false);
    q_vec->setVectorColor({1.0, 0.5, 0.0}); // Orange vectors

    polyscope::show();
    return 0;
}
