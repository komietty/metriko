#ifndef METRIKO_VISUALIZE_VECTORFIELD_H
#define METRIKO_VISUALIZE_VECTORFIELD_H
#include "./visualize_common.h"
#include "metriko/core/vectorfield/face_rosy_field.h"

namespace metriko::visualizer {
inline void visualize_frosy_field(
    polyscope::SurfaceMesh* surf,
    const Hmesh& hm,
    const FaceRosyField& rawf,
    const FaceRosyField& cmbf,
    const int rosyN = 4,
    const bool show = true
) {
    MatXd rawInt(hm.nF, 2);
    MatXd cmbInt(hm.nF, 2);
    MatXd rawExt(hm.nF, 3 * rosyN);
    MatXd cmbExt(hm.nF, 3 * rosyN);
    for (Face f: hm.faces) {
        complex rc0 = rawf.field(f.id, 0);
        complex rc1 = rawf.field(f.id, 1);
        complex rc2 = rawf.field(f.id, 2);
        complex rc3 = rawf.field(f.id, 3);
        complex cc0 = cmbf.field(f.id, 0);
        complex cc1 = cmbf.field(f.id, 1);
        complex cc2 = cmbf.field(f.id, 2);
        complex cc3 = cmbf.field(f.id, 3);
        rawInt.row(f.id) = Row2d(rc0.real(), rc0.imag()).normalized();
        cmbInt.row(f.id) = Row2d(cc0.real(), cc0.imag()).normalized();

        rawExt.block(f.id, 0, 1, 3) = (rc0.real() * f.basisX() + rc0.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 3, 1, 3) = (rc1.real() * f.basisX() + rc1.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 6, 1, 3) = (rc2.real() * f.basisX() + rc2.imag() * f.basisY()).normalized();
        rawExt.block(f.id, 9, 1, 3) = (rc3.real() * f.basisX() + rc3.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 0, 1, 3) = (cc0.real() * f.basisX() + cc0.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 3, 1, 3) = (cc1.real() * f.basisX() + cc1.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 6, 1, 3) = (cc2.real() * f.basisX() + cc2.imag() * f.basisY()).normalized();
        cmbExt.block(f.id, 9, 1, 3) = (cc3.real() * f.basisX() + cc3.imag() * f.basisY()).normalized();
    }
    auto rawFQ = surf->addFaceVectorQuantity("raw ext", rawExt.block(0, 0, rawExt.rows(), 3));
    auto cmbFQ = surf->addFaceVectorQuantity("cmb ext", cmbExt.block(0, 0, cmbExt.rows(), 3));
    rawFQ->setEnabled(show);
    cmbFQ->setEnabled(show);
    rawFQ->setVectorLengthScale(0.004);
    cmbFQ->setVectorLengthScale(0.004);
}
}
#endif
