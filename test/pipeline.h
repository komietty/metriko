#ifndef METRIKO_TEST_PIPELINE_H
#define METRIKO_TEST_PIPELINE_H
// テスト共通: OBJ から field -> param -> mc -> Tmesh -> TmeshMut までを一括構築する。
//
// 参照保持の依存を満たすため、必要なオブジェクトを本構造体が所有して生かし続ける:
//   ・mg  は hm と uv2(cf) を参照保持
//   ・tmm は hm を参照保持
//   ・matching/singular は mg 構築時のみ使用（保持しない）→ ctor 内ローカルで可
// メンバ宣言順＝破棄は逆順なので、被参照側（hm/uv2）が参照側（mg/tm/tmm）より後に消える。
// 内部に this ポインタ（Mcurv.mg 等）を持つためコピー/ムーブ不可。
#include <igl/readOBJ.h>
#include <optional>
#include <string>
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/igm/parameterization.h"
#include "metriko/core/quantization/quantization.h"
#include "metriko/core/tmesh/tmesh_mut.h"

namespace metriko {

struct TmeshPipeline {
    std::optional<Hmesh>     hm;
    VecXc                    uv2;
    std::optional<mc::Mgrph> mg;
    std::optional<Tmesh>     tm;
    VecXd                    X;
    std::optional<TmeshMut>  tmm;
    bool ok = false;

    TmeshPipeline(const std::string& mesh_path, double scale, int N = 4) {
        MatXd V; MatXi F;
        igl::readOBJ(mesh_path, V, F);
        if (V.rows() == 0 || F.rows() == 0) return;   // ok=false のまま返す
        hm.emplace(V, F);
        Hmesh& h = *hm;

        FaceRosyField rawf(h, N, FieldType::Smoothest);
        rawf.computeMatching(MatchingType::Principal);
        auto seam = compute_seam(rawf);
        auto cutm = compute_cut_mesh(h, seam);
        auto cmbf = compute_combbed_field(rawf, seam);

        MatXd ext(h.nF, 3 * N);
        for (Face f : h.faces)
            for (int k = 0; k < N; ++k) {
                complex c = cmbf->field(f.id, k);
                ext.block(f.id, 3 * k, 1, 3) = (c.real() * f.basisX() + c.imag() * f.basisY()).normalized();
            }

        RosyParameterization rp(h, *cutm, ext, cmbf->singular, cmbf->matching, seam, N, scale);
        rp.seamless = false;
        rp.localInjectivity = true;
        rp.verbose = false;
        rp.setup();
        rp.integ();

        uv2.resize(h.nF * 3);
        for (Face f : h.faces) {
            uv2(f.id * 3 + 0) = complex{ rp.cfn(f.id, 0), rp.cfn(f.id, 1) };
            uv2(f.id * 3 + 1) = complex{ rp.cfn(f.id, 4), rp.cfn(f.id, 5) };
            uv2(f.id * 3 + 2) = complex{ rp.cfn(f.id, 8), rp.cfn(f.id, 9) };
        }

        mg.emplace(h, uv2, cmbf->matching, cmbf->singular);
        tm.emplace(*mg);
        X = compute_quantization(*tm, *mg);
        tmm.emplace(*mg, *tm, X);
        ok = true;
    }

    TmeshPipeline(const TmeshPipeline&)            = delete;
    TmeshPipeline& operator=(const TmeshPipeline&) = delete;
};

} // namespace metriko

#endif
