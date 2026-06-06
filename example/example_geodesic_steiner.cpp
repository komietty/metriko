//
// Created by saki on 2026/06/06.
//
//  Steiner 点 + Dijkstra による近似測地線（実体は hpath.h の approx_shortest_path）。
//  許可領域マップ(eid -> [r0,r1]) で通れるエッジ領域を制限できる。
//  始点・終点は vert / edge / face内点 を自由に指定可能。

#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include "metriko/core/hmesh/hmesh.h"
#include "metriko/core/hmesh/hmloc.h"
#include "metriko/core/hmesh/hpath.h"

using namespace metriko;

int main(int argc, char** argv) {
    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    // 許可領域：エッジID -> [r0,r1]（マップに無いエッジ = 横切り禁止）
    umap<int, Row2d> allowed;
    for (Edge e : hm.edges) allowed[e.id] = Row2d(0.0, 1.0);
    // 制限例（各エッジ中央60%だけ通す）:
    //   for (Edge e : hm.edges) allowed[e.id] = Row2d(0.2, 0.8);

    // 面内点を重心座標で作り、面ローカル xy を持つ HmLocOnF に
    auto faceInteriorLoc = [&](int fid, double b0, double b1, double b2) -> HmLoc {
        Face  f = hm.faces[fid];
        Row3d p = b0 * f.half().tail().pos()
                + b1 * f.half().next().tail().pos()
                + b2 * f.half().prev().tail().pos();
        Row3d v = p - f.half().tail().pos();
        return HmLoc(HmLocOnF{ fid, complex(v.dot(f.basisX()), v.dot(f.basisY())) });
    };

    // 始点・終点（face / edge / vert を自由に指定できる）
    HmLoc loc_bgn = faceInteriorLoc(0, 0.6, 0.3, 0.1);
    HmLoc loc_end = faceInteriorLoc(std::min(5000, (int)hm.nF - 1), 0.2, 0.5, 0.3);
    // 例:  HmLoc loc_end = HmLoc(HmLocOnV{ 0 });
    //      HmLoc loc_end = HmLoc(HmLocOnE{ eid, 0.4 });

    const int n_div = 8;  // 単位区間あたりの分割数
    vec<HmLoc> path = approx_shortest_path(n_div, hm, loc_bgn, loc_end, allowed);

    if (path.empty()) {
        std::cerr << "empty path" << std::endl;
        return 1;
    }

    // HmLoc -> 3D 座標
    MatXd P((int)path.size(), 3);
    for (int i = 0; i < (int)path.size(); ++i) P.row(i) = get_ptloc_pos(hm, path[i]);

    // 長さ0の不正な segment が無いか確認
    int n_degenerate = 0;
    for (int i = 0; i + 1 < (int)path.size(); ++i) {
        double len = (P.row(i + 1) - P.row(i)).norm();
        if (len < 1e-9) {
            std::cerr << "  degenerate segment at [" << i << ", " << i + 1
                      << "] len=" << len << std::endl;
            ++n_degenerate;
        }
    }
    std::cout << "degenerate (zero-length) segments = " << n_degenerate << std::endl;

    MatXd ends(2, 3);
    ends.row(0) = get_ptloc_pos(hm, loc_bgn);
    ends.row(1) = get_ptloc_pos(hm, loc_end);

    // 表示
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::registerSurfaceMesh("mesh", V, F);

    auto* psPath = polyscope::registerCurveNetworkLine("steiner geodesic", P);
    psPath->setColor({1.0, 0.0, 0.0});
    psPath->setRadius(0.004);

    auto* psEnds = polyscope::registerPointCloud("endpoints", ends);
    psEnds->setPointColor({0.0, 1.0, 0.2});
    psEnds->setPointRadius(0.008);

    polyscope::show();
    return 0;
}
