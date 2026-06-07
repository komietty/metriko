// パイプライン結合テスト（field→param→mc→Tmesh→TmeshMut→allowed_range）。polyscope 不要。
#include <igl/readOBJ.h>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"
#include "check.h"

using namespace metriko;

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: test_tmesh <gridscale> <mesh.obj> [more.obj ...]\n"; return 2; }
    const int N = 4;
    const double scale = std::stod(argv[1]);

    for (int a = 2; a < argc; ++a) {                 // 複数メッシュをループ（Program arguments で指定）
        const char* mesh = argv[a];
        MatXd V; MatXi F;
        igl::readOBJ(mesh, V, F);
        CHECK(V.rows() > 0 && F.rows() > 0);   // 読み込み失敗(欠落/空)をクリーンに検出
        Hmesh hm(V, F);

        FaceRosyField rawf(hm, N, FieldType::Smoothest);
        rawf.computeMatching(MatchingType::Principal);
        auto seam = compute_seam(rawf);
        auto cutm = compute_cut_mesh(hm, seam);
        auto cmbf = compute_combbed_field(rawf, seam);

        MatXd ext(hm.nF, 3 * N);
        for (Face f : hm.faces)
            for (int k = 0; k < N; ++k) {
                complex c = cmbf->field(f.id, k);
                ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
            }

        RosyParameterization rp(hm, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, scale);
        rp.seamless = false;
        rp.localInjectivity = true;
        rp.verbose = false;
        rp.setup();
        rp.integ();

        VecXc uv2(hm.nF * 3);
        for (Face f : hm.faces) {
            uv2(f.id * 3 + 0) = complex{ rp.cfn(f.id, 0), rp.cfn(f.id, 1) };
            uv2(f.id * 3 + 1) = complex{ rp.cfn(f.id, 4), rp.cfn(f.id, 5) };
            uv2(f.id * 3 + 2) = complex{ rp.cfn(f.id, 8), rp.cfn(f.id, 9) };
        }

        mc::Mgrph mg(hm, uv2, cmbf->matching, cmbf->singular);
        Tmesh    tm(mg);
        VecXd    X = compute_quantization(tm, mg);
        TmeshMut tmm(mg, tm);

        // --- Tmesh / TmeshMut の不変条件 ---
        CHECK(tm.nTH == 2 * tm.nTE);
        CHECK(tmm.thalfs.size() == 2 * tmm.tedges.size());
        CHECK(tmm.tnodes.size() == mg.mnodes.size());
        CHECK(!tmm.tquads.empty());
        CHECK(tmm.tquads.size() == tm.tquads.size());

        // --- allowed_range が局所化する（空でなく、全エッジ未満）---
        auto r = tmm.allowed_range(tmm.tquads.front().id);
        CHECK(!r.empty());
        CHECK((int)r.size() < hm.nE);
        for (auto& [eid, r0, r1] : r) {       // レンジが有効
            CHECK(eid >= 0 && eid < hm.nE);
            CHECK(r0 < r1);
        }

        std::cout << "[test_tmesh] OK  " << mesh << "  nTQ=" << tm.nTQ
                  << "  allowed=" << r.size() << "/" << hm.nE << "\n";
    }
    return 0;
}
