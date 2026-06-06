//
// Created by saki on 2026/06/06.
//

#include <igl/readOBJ.h>
#include <igl/exact_geodesic.h>          // igl::geodesic 名前空間もこれで入る
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include "metriko/core/hmesh/hmesh.h"

using namespace metriko;
namespace geo = igl::geodesic;

int main(int argc, char** argv) {
    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);
    Hmesh hm(V, F);

    // --- geodesic::Mesh を構築（libigl ラッパーと同じ手順） ---
    std::vector<double> points(V.rows() * 3);
    std::vector<int>    faces (F.rows() * 3);
    for (int i = 0; i < (int)points.size(); ++i) points[i] = V(i / 3, i % 3);
    for (int i = 0; i < (int)faces.size();  ++i) faces[i]  = F(i / 3, i % 3);

    geo::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    geo::GeodesicAlgorithmExact algo(&mesh);

    // --- 面 fid 内の点を重心座標(b0,b1,b2)で作る（b*>0, 和=1 なら必ず内部） ---
    auto face_interior_pos = [&](int fid, double b0, double b1, double b2) -> Row3d {
        Face f = hm.faces[fid];
        Row3d p0 = f.half().tail().pos();          // v0
        Row3d p1 = f.half().next().tail().pos();   // v1
        Row3d p2 = f.half().prev().tail().pos();   // v2
        return b0 * p0 + b1 * p1 + b2 * p2;
    };
    // (面ポインタ + その面上の実座標) のペアで面内任意点を表す
    auto to_surface_point = [&](int fid, const Row3d& p) {
        return geo::SurfacePoint(&mesh.faces()[fid], p.x(), p.y(), p.z(), geo::FACE);
    };

    // --- 始点・終点：別々の面の内部点 ---
    int  srcFid = 0;
    int  dstFid = std::min(5000, (int)hm.nF - 1);
    Row3d srcPos = face_interior_pos(srcFid, 0.6, 0.3, 0.1);
    Row3d dstPos = face_interior_pos(dstFid, 0.2, 0.5, 0.3);

    geo::SurfacePoint source = to_surface_point(srcFid, srcPos);
    geo::SurfacePoint dst    = to_surface_point(dstFid, dstPos);

    // --- 厳密伝播 + 逆トレース（距離場を経由しない MMP の window 逆たどり） ---
    std::vector<geo::SurfacePoint> sources{ source };
    algo.propagate(sources);

    double dist = 0.0;
    algo.best_source(dst, dist);             // 到達距離（ガード用）
    if (dist >= 1e100) {                      // GEODESIC_INF: 非連結
        std::cerr << "target unreachable (disconnected?)" << std::endl;
        return 1;
    }

    std::vector<geo::SurfacePoint> path;
    algo.trace_back(dst, path);              // path: dst -> ... -> source
    std::reverse(path.begin(), path.end());  // source -> dst の順に
    std::cout << "exact geodesic length = " << dist
              << ",  #points = " << path.size() << std::endl;

    // SurfacePoint -> 3D 座標
    MatXd P((int)path.size(), 3);
    for (int i = 0; i < (int)path.size(); ++i)
        P.row(i) << path[i].x(), path[i].y(), path[i].z();

    // 端点（面内点）を確認用に
    MatXd ends(2, 3);
    ends.row(0) = srcPos;
    ends.row(1) = dstPos;

    // --- 表示 ---
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::registerSurfaceMesh("mesh", V, F);

    auto* psPath = polyscope::registerCurveNetworkLine("exact geodesic", P);
    psPath->setColor({1.0, 0.0, 0.0});
    psPath->setRadius(0.004);

    auto* psEnds = polyscope::registerPointCloud("endpoints", ends);
    psEnds->setPointColor({0.0, 1.0, 0.2});
    psEnds->setPointRadius(0.008);

    polyscope::show();
    return 0;
}
