#include <chrono>
#include <functional>
#include <map>
#include <set>
#include <igl/readOBJ.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include "metriko/lib.h"

using namespace metriko;

static vec<double> four_colour_patches(const MatXi& qidx, const VecXi& tqid) {
    umap<int, int> pid;
    for (int t: tqid) if (t >= 0) pid.try_emplace(t, (int)pid.size());
    const int np = pid.size();

    vec<std::set<int>> adj(np);
    std::map<std::pair<int, int>, int> edge_quad;
    for (int q = 0; q < qidx.rows(); ++q) {
        if (tqid(q) < 0) continue;
        for (int j = 0; j < qidx.cols(); ++j) {
            auto e = std::minmax(qidx(q, j), qidx(q, (j + 1) % qidx.cols()));
            auto [it, fresh] = edge_quad.try_emplace(e, q);
            if (fresh) continue;
            const int a = pid.at(tqid(q)), b = pid.at(tqid(it->second));
            if (a != b) { adj[a].insert(b); adj[b].insert(a); }
        }
    }

    vec<int> col(np, -1);
    auto uses_of = [&](int p) { std::array<int, 4> u{}; for (int r: adj[p]) if (col[r] >= 0) ++u[col[r]]; return u; };
    auto next_patch = [&] {   // dsatur: most distinct neighbour colours, then most neighbours
        int best = -1, best_sat = -1, best_deg = -1;
        for (int p = 0; p < np; ++p) {
            if (col[p] >= 0) continue;
            const auto u = uses_of(p);
            const int sat = rg::count_if(u, [](int c) { return c > 0; }), deg = adj[p].size();
            if (sat > best_sat || (sat == best_sat && deg > best_deg)) { best = p; best_sat = sat; best_deg = deg; }
        }
        return best;
    };

    long steps = 0;
    constexpr long budget = 1000000;
    std::function<bool()> search = [&] {
        const int p = next_patch();
        if (p < 0) return true;
        const auto u = uses_of(p);
        for (int c = 0; c < 4; ++c) {
            if (u[c] > 0 || ++steps > budget) continue;
            col[p] = c;
            if (search()) return true;
        }
        col[p] = -1;
        return false;
    };

    if (!search()) {
        rg::fill(col, -1);
        int conflicts = 0;
        for (int p = next_patch(); p >= 0; p = next_patch()) {
            const auto u = uses_of(p);
            col[p] = rg::min_element(u) - u.begin();
            if (u[col[p]] > 0) ++conflicts;
        }
        std::println("[quad patch] four colours: none found in {} steps, {} patches share a colour with a neighbour", budget, conflicts);
    }

    vec<double> res(tqid.size(), -1);
    for (int q = 0; q < res.size(); ++q) if (tqid(q) >= 0) res[q] = col[pid.at(tqid(q))] / 3.;
    return res;
}

static void visualize_quad_patch(const MatXd& qv, const MatXi& qidx, const VecXi& tqid) {
    auto at = [&](int v) { return Row3d(qv.row(v)); };
    struct Segments {
        vec<glm::vec3> ns;
        vec<std::array<size_t, 2>> es;

        void add(const Row3d& a, const Row3d& b) {
            es.push_back({ns.size(), ns.size() + 1});
            ns.push_back({a.x(), a.y(), a.z()});
            ns.push_back({b.x(), b.y(), b.z()});
        }
        polyscope::CurveNetwork* show(const std::string& name, double radius, bool enabled) const {
            auto* cn = polyscope::registerCurveNetwork(name, ns, es);
            cn->setMaterial("flat");
            cn->setRadius(radius);
            cn->setEnabled(enabled);
            return cn;
        }
    };

    std::map<std::pair<int, int>, int> edge_quad;
    Segments bnd;
    for (int q = 0; q < qidx.rows(); ++q)
    for (int j = 0; j < 4; ++j) {
        const int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
        auto [it, fresh] = edge_quad.try_emplace(std::minmax(a, b), q);
        if (!fresh && tqid(q) == tqid(it->second)) it->second = -1;   // inside a patch
    }
    for (auto& [e, q]: edge_quad) if (q >= 0) bnd.add(at(e.first), at(e.second));
    bnd.show("patch boundaries", 0.001, true)->setColor({0., 0., 0.});

    auto* surf = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    surf->setShadeStyle(polyscope::MeshShadeStyle::Flat);
    surf->setEdgeWidth(1.);
    surf->addFaceScalarQuantity("tqid", tqid);
    vec<double> golden(tqid.size());   // tqid spread by the golden ratio: consecutive ids land far apart on the colour map
    for (int q = 0; q < tqid.size(); ++q) golden[q] = tqid(q) < 0 ? -1 : std::fmod(tqid(q) * 0.618033988749895, 1.);
    auto* col = surf->addFaceScalarQuantity("patch colour", golden);
    col->setColorMap("coolwarm");
    col->setEnabled(true);
    auto* col4 = surf->addFaceScalarQuantity("patch colour (4)", four_colour_patches(qidx, tqid));
    col4->setColorMap("coolwarm");
    col4->setMapRange({0., 1.});
}

int main(int argc, char** argv) {
    if (argc < 3) { std::println(stderr, "usage: example <mesh.obj> <scale>"); return 2; }
    MatXd V;
    MatXi F;
    igl::readOBJ(argv[1], V, F);

    const auto t0  = std::chrono::steady_clock::now();
    const auto res = compute_quadrangulation(V, F, std::stod(argv[2]));
    if (!res) { std::println(stderr, "failed: {}", res.error()); return 1; }
    std::println("[time] total {:8.0f} ms", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());

    polyscope::options::verbosity = 0;
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    auto* base = polyscope::registerSurfaceMesh("base mesh", V, F);
    base->setEdgeWidth(0.7);
    base->setEnabled(false);
    base->setMaterial("flat");
    base->setSurfaceColor(glm::vec3(0.3, 0.3, 0.3));
    visualize_quad_patch(res->pos, res->idx, res->val.col(0));
    polyscope::show();
    return 0;
}
