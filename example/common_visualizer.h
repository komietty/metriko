#ifndef METRIKO_EXAMPLE_COMMON_VISUALIZER_H
#define METRIKO_EXAMPLE_COMMON_VISUALIZER_H
// polyscope views of every stage of the pipeline: mesh and field, motorcycle t-mesh,
// collapsed / snapped t-mesh, qex output, and the quad patches
#include <set>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/hmesh/utilities.h"
#include "metriko/nvec/face_rosy_field.h"
#include "metriko/tmesh/emesh.h"
#include "metriko/patching.h"

namespace metriko::visualizer {

inline void visualize_init() {
    //polyscope::options::buildGui = false;
    polyscope::options::verbosity = 0;
    polyscope::init();
    polyscope::view::bgColor = std::array<float, 4>{0.02, 0.02, 0.02, 1};
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
}

//------------------------------------------------------------------------------
// small builders shared by the views below
//------------------------------------------------------------------------------

inline glm::vec3 to_glm(const Row3d& p) { return {p.x(), p.y(), p.z()}; }

// a set of segments with per-segment scalars, registered as a flat curve network
struct Segments {
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    std::map<std::string, vec<double>> scalars;   // per segment

    void add(const Row3d& a, const Row3d& b) {
        es.push_back({ns.size(), ns.size() + 1});
        ns.push_back(to_glm(a));
        ns.push_back(to_glm(b));
    }
    void add(const Hmesh& hm, const vec<int>& nids, const vec<HmLoc>& tnodes) {   // a node chain of the t-mesh
        for (size_t k = 0; k + 1 < nids.size(); ++k) add(get_ptloc_pos(hm, tnodes[nids[k]]), get_ptloc_pos(hm, tnodes[nids[k + 1]]));
    }
    void scalar(const std::string& name, double v) { scalars[name].push_back(v); }

    polyscope::CurveNetwork* show(const std::string& name, double radius, bool enabled, const char* enabled_scalar = nullptr) const {
        auto* cn = polyscope::registerCurveNetwork(name, ns, es);
        cn->setMaterial("flat");
        cn->setRadius(radius);
        cn->setEnabled(enabled);
        cn->resetTransform();
        for (auto& [n, v]: scalars) cn->addEdgeScalarQuantity(n, v)->setEnabled(enabled_scalar && n == enabled_scalar);
        return cn;
    }
};

// a point set with per-point scalars
struct Points {
    vec<glm::vec3> ps;
    std::map<std::string, vec<double>> scalars;

    void add(const Row3d& p) { ps.push_back(to_glm(p)); }
    void scalar(const std::string& name, double v) { scalars[name].push_back(v); }

    polyscope::PointCloud* show(const std::string& name, double radius, bool enabled, const char* enabled_scalar = nullptr) const {
        auto* pc = polyscope::registerPointCloud(name, ps);
        pc->setPointRadius(radius);
        pc->setEnabled(enabled);
        pc->resetTransform();
        for (auto& [n, v]: scalars) pc->addScalarQuantity(n, v)->setEnabled(enabled_scalar && n == enabled_scalar);
        return pc;
    }
};

//------------------------------------------------------------------------------
// mesh, field, seam
//------------------------------------------------------------------------------

// the mesh, with its corner uv (hm.cfn) as a checker when the mesh carries one
inline polyscope::SurfaceMesh* visualize_mesh(const Hmesh& hm, const bool show = true, const std::string& name = "mesh", const bool show_uv = true) {
    auto* surf = polyscope::registerSurfaceMesh(name, hm.pos, hm.idx);
    surf->setEdgeWidth(0.7);
    surf->setEnabled(show);
    surf->setMaterial("flat");
    surf->setSurfaceColor(glm::vec3(0.3, 0.3, 0.3));
    if (hm.cfn.size() == hm.nF * 3) {
        auto* prms = surf->addParameterizationQuantity("uv_param", hm.cfn);
        prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
        prms->setCheckerSize(1.);
        prms->setEnabled(show_uv);
    }
    return surf;
}

// the first direction of the raw and the combed field, as face vectors on `surf`
inline void visualize_frosy_field(polyscope::SurfaceMesh* surf, const FaceRosyField& rawf, const FaceRosyField& cmbf, const int rosyN = 4, const bool show = true) {
    MatXd raw = compute_extrinsic_field(rawf, rosyN);
    MatXd cmb = compute_extrinsic_field(cmbf, rosyN);
    for (auto [name, ext]: {std::pair{"raw ext", &raw}, std::pair{"cmb ext", &cmb}}) {
        auto* q = surf->addFaceVectorQuantity(name, ext->leftCols(3));
        q->setEnabled(show);
        q->setVectorLengthScale(0.004);
    }
}

inline void visualize_seam(const Hmesh& hm, const vec<bool>& seam, const VecXi& matching = VecXi(), const std::string& name = "seam", const bool show = true) {
    Segments sg;
    for (auto e: hm.edges) {
        if (!seam[e.id]) continue;
        sg.add(e.half().tail().pos(), e.half().head().pos());
        if (matching.rows() > 0) sg.scalar("matching", matching[e.id]);
    }
    sg.show(name, 0.001, show);
}

// interior singulars whose uv cone angle (the sum of the uv corner angles at the vertex) differs from the one its index
// asks for, (4 - index) * pi / 2, drawn as points. on `surf`, a face scalar marks folded faces (non-positive uv area)
// with bit 1 and the one-ring faces of those singulars with bit 2. returns the number of wrongly wound singulars
inline int visualize_wrong_cones(polyscope::SurfaceMesh* surf, const Hmesh& hm, const VecXi& singular, const bool show = true) {
    auto uv = [&](int f, int j) { return hm.cfn(f * 3 + j); };
    vec<int> mark(hm.nF, 0);
    for (Face f: hm.faces) if (cross(uv(f.id, 1) - uv(f.id, 0), uv(f.id, 2) - uv(f.id, 0)) <= 0) mark[f.id] |= 1;

    Points pts;
    for (Vert v: hm.verts) {
        const int s = singular[v.id];
        if (s == 0 || v.isBoundary()) continue;
        double cone = 0;
        for (Half h: v.adjHalfs()) {
            const int f = h.face().id;
            int j = 0;
            while (hm.idx(f, j) != v.id) ++j;
            cone += std::arg((uv(f, (j + 2) % 3) - uv(f, j)) / (uv(f, (j + 1) % 3) - uv(f, j)));
        }
        const double want = (4 - s) * M_PI / 2;
        if (std::abs(cone - want) < 1e-3) continue;
        for (Half h: v.adjHalfs()) mark[h.face().id] |= 2;
        pts.add(v.pos());
        pts.scalar("index", s);
        pts.scalar("vid", v.id);
        pts.scalar("cone / pi", cone / M_PI);
        std::println("[cone] singular {} index {:+d}: uv cone {:.2f} pi, should be {:.2f} pi", v.id, s, cone / M_PI, want / M_PI);
    }
    surf->addFaceScalarQuantity("folded (1) / wrong cone ring (2)", vec<double>(mark.begin(), mark.end()))->setEnabled(false);
    pts.show("wrong cones", 0.003, show, "index");
    return (int)pts.ps.size();
}

//------------------------------------------------------------------------------
// t-mesh right after the motorcycle graph and the quantization (before the collapse)
//------------------------------------------------------------------------------

// every tedge as its node chain, with the tquads on both sides, its length R and its quantized length X (x < 0
// before the quantization)
inline void visualize_tedge(const Emesh& tm, const bool show = true) {
    Segments sg;
    for (const Ehalf& th: tm.thalfs) {
        if (th.id == -1 || !th.cano) continue;
        const auto& nids = tm.tedges[th.teid].nids;
        for (size_t k = 0; k + 1 < nids.size(); ++k) {
            sg.add(get_ptloc_pos(tm.hm, tm.tnodes[nids[k]]), get_ptloc_pos(tm.hm, tm.tnodes[nids[k + 1]]));
            sg.scalar("teid", th.teid);
            sg.scalar("tqid1", th.tqid);
            sg.scalar("tqid2", tm.thalfs[th.twid].tqid);
            sg.scalar("R", th.r);
            sg.scalar("X", th.x);
        }
    }
    sg.show("tedges", 0.0005, show)->setColor({0., 0., 0.});
}

// tquads whose boundary walk returns to a corner node it already passed. the patch is then not a
// disk, which the per-tquad tutte embedding downstream assumes it is. it happens around
// high-valence singularities: the patch leaves along one separatrix and comes back along another
inline void visualize_pinched_tquads(const Emesh& tm, const bool show = true, const double scale = 0.002) {
    Segments sg;
    Points pts;
    for (const Equad& tq: tm.live_tquads()) {
        if (tq.data.empty()) continue;
        std::map<int, int> visits;
        for (const auto& d: tq.data) visits[tm.thalfs[d.thid].nid_fr()]++;
        if (rg::none_of(visits, [](const auto& p) { return p.second > 1; })) continue;

        std::set<int> done;
        for (const auto& [thid, side]: tq.data) {
            auto nids = tm.tedges[tm.thalfs[thid].teid].nids;
            if (!tm.thalfs[thid].cano) rg::reverse(nids);
            for (size_t k = 0; k + 1 < nids.size(); ++k) {
                Row3d p = get_ptloc_pos(tm.hm, tm.tnodes[nids[k]]);
                Row3d q = get_ptloc_pos(tm.hm, tm.tnodes[nids[k + 1]]);
                sg.add(p, q);
                sg.scalar("tqid", tq.id);
                sg.scalar("side", side);
                if (auto it = visits.find(nids[k]); it != visits.end() && it->second > 1 && done.insert(nids[k]).second) {
                    pts.add(p);
                    pts.scalar("nid", nids[k]);
                    pts.scalar("tqid", tq.id);
                    pts.scalar("visits", it->second);
                }
            }
        }
    }
    if (!pts.ps.empty()) std::println("[pinched] {} pinch nodes over {} tquads", pts.ps.size(), tm.tquads.size());
    sg.show("pinched tquad boundary", scale, show, "side");
    pts.show("pinch node", scale * 2.5, show, "nid");
}

//------------------------------------------------------------------------------
// collapsed / snapped t-mesh (Emesh)
//------------------------------------------------------------------------------

// allowed corridor of a shortest-path query: the admissible sub-segment [fr, to] of every edge
inline void visualize_allowed_range(const Hmesh& hm, const vec<Erng>& allowed, const std::string& name = "allowed range", const bool show = true, const double scale = 0.0015) {
    Segments sg;
    for (const Erng& rng: allowed) {
        Edge e = hm.edges[rng.id];
        sg.add(e.lerp(rng.fr), e.lerp(rng.to));
        sg.scalar("eid", rng.id);
        sg.scalar("span", rng.span());
    }
    sg.show(name, scale, show, "eid");
}

// the live tedges as node chains on the mesh, with their quantized length and the tquads on both sides
inline void visualize_tedges(const Hmesh& hm, const Emesh& em, const std::string& name = "tedges", const bool show = true, const double scale = 0.001) {
    umap<int, const Ehalf*> cano;   // teid -> canonical thalf
    for (const auto& th: em.thalfs) if (th.id != -1 && th.cano) cano[th.teid] = &th;
    Segments sg;
    for (const auto& [teid, nids]: em.live_tedges()) {
        const size_t n0 = sg.es.size();
        sg.add(hm, nids, em.tnodes);
        const auto* th = cano.contains(teid) ? cano.at(teid) : nullptr;
        for (size_t k = n0; k < sg.es.size(); ++k) {
            sg.scalar("teid", teid);
            sg.scalar("x", th ? th->x : -1);
            sg.scalar("tqid0", th ? th->tqid : -1);
            sg.scalar("tqid1", th ? em.thalfs[th->twid].tqid : -1);
        }
    }
    sg.show(name, scale, show, "teid");
}

// one curve network per live tquad, with the side and the quantized length of every boundary thalf
inline void visualize_tquads(const Hmesh& hm, const Emesh& em, const bool show = true, const double scale = 0.001, const vec<int>& only = {}) {
    for (const auto& [id, data]: em.live_tquads()) {
        if (!only.empty() && !rg::contains(only, id)) continue;
        Segments sg;
        for (const Edata& d: data) {
            const Ehalf& th = em.thalfs[d.thid];
            const size_t n0 = sg.es.size();
            sg.add(hm, em.tedges[th.teid].nids, em.tnodes);
            for (size_t k = n0; k < sg.es.size(); ++k) {
                sg.scalar("side", d.side);
                sg.scalar("x", th.x);
                sg.scalar("r", th.r);
                sg.scalar("thid", d.thid);
            }
        }
        if (!sg.ns.empty()) sg.show(std::format("tq {:03}", id), scale, show);
    }
}

// t-nodes still lying on an edge or inside a face (OnV nodes are snapped)
inline void visualize_unsnapped_tnodes(const Hmesh& hm, const Emesh& em, const bool show = true, const double scale = 0.003) {
    Points pts;
    std::set<int> seen;   // shared nodes (junctions / crossings) appear in several chains
    for (const auto& [teid, nids]: em.live_tedges())
    for (int nid: nids) {
        if (!seen.insert(nid).second) continue;
        double type = -1, id = -1;
        std::visit(overloaded{
            [&](const HmLocOnE& e) { type = 0; id = e.id; },
            [&](const HmLocOnH& h) { type = 1; id = h.id; },
            [&](const HmLocOnF& f) { type = 2; id = f.id; },
            [&](const auto&)       {},
        }, em.tnodes[nid]);
        if (type < 0) continue;
        pts.add(get_ptloc_pos(hm, em.tnodes[nid]));
        pts.scalar("type (0:E 1:H 2:F)", type);
        pts.scalar("elem id", id);
        pts.scalar("node id", nid);
    }
    pts.show("unsnapped tnodes", scale, show, "type (0:E 1:H 2:F)");
}

//------------------------------------------------------------------------------
// qex
//------------------------------------------------------------------------------

inline void visualize_qverts(const vec<qex::Qvert>& vqvs, const vec<qex::Qvert>& eqvs, const vec<qex::Qvert>& fqvs, const double scale = 0.001, const bool show = true) {
    for (auto [name, qvs]: {std::pair{"VQV", &vqvs}, std::pair{"EQV", &eqvs}, std::pair{"FQV", &fqvs}}) {
        Points pts;
        for (const auto& q: *qvs) pts.add(q.pos);
        pts.show(name, scale, show);
    }
}

// every port slightly offset along its direction, with the ids needed to follow the tracer log
inline void visualize_qports(const Hmesh& hm, const vec<qex::Qport>& ports, const double scale = 0.001, const bool show = true) {
    auto dir_index = [](complex d) { return d.real() > 0.5 ? 0 : d.imag() > 0.5 ? 1 : d.real() < -0.5 ? 2 : 3; };
    Points pts;
    for (const auto& q: ports) {
        pts.add(hm.faces[q.fid].uv2pos(q.uv + q.dir * 0.15));
        pts.scalar("idx", q.idx);
        pts.scalar("dir", dir_index(q.dir));
        pts.scalar("vid", q.vid);
        pts.scalar("eid", q.eid);
        pts.scalar("fid", q.fid);
        pts.scalar("u", q.uv.real());
        pts.scalar("v", q.uv.imag());
        pts.scalar("next", q.next_id);
        pts.scalar("prev", q.prev_id);
    }
    pts.show("q_ports", scale, show);
}

inline void visualize_qedges(const vec<qex::Qedge>& qedges, const double scale = 0.001, const bool show = false) {
    Segments sg;
    for (auto& q: qedges) {
        sg.add(q.port1.pos, q.port2.pos);
        sg.scalar("p1 idx", q.port1.idx);
        sg.scalar("p2 idx", q.port2.idx);
    }
    sg.show("q_edges", scale, show);
}

// the patches coloured with four colours, equally spaced on the colour map (0, 1/3, 2/3, 1), so that patches sharing a
// quad edge differ: a backtracking search in dsatur order (the patch with the most distinct colours among its coloured
// neighbours first), within a step budget. four colours are not always enough off the sphere; then the greedy dsatur
// pass gives a patch with all four taken the colour fewest of its neighbours use, and the number of such conflicts is
// printed. unlabeled quads (tqid < 0) stay -1
inline vec<double> four_colour_patches(const MatXi& qidx, const VecXi& tqid) {
    umap<int, int> pid;   // tqid -> patch index
    for (int t: tqid) if (t >= 0) pid.try_emplace(t, (int)pid.size());
    const int np = pid.size();

    // adjacency: two quads on one quad edge with different patches
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

    if (!search()) {   // no proper four colouring found within the budget: greedy, with the fewest clashes
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

// the quad mesh with its patches: the patch boundaries (quad edges between two patches, or on the mesh boundary),
// the irregular vertices, and the quads coloured by patch
inline void visualize_quad_patch(const MatXd& qv, const MatXi& qidx, const VecXi& tqid) {
    auto at = [&](int v) { return Row3d(qv.row(v)); };

    std::map<std::pair<int, int>, int> edge_quad;   // quad edge -> one quad on it
    Segments bnd;
    VecXi valence = VecXi::Zero(qv.rows());
    for (int q = 0; q < qidx.rows(); ++q)
    for (int j = 0; j < 4; ++j) {
        const int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
        ++valence(a);
        auto [it, fresh] = edge_quad.try_emplace(std::minmax(a, b), q);
        if (!fresh && tqid(q) == tqid(it->second)) it->second = -1;   // inside a patch
    }
    for (auto& [e, q]: edge_quad) if (q >= 0) bnd.add(at(e.first), at(e.second));
    bnd.show("patch boundaries", 0.001, true)->setColor({0., 0., 0.});

    Points irr;
    for (int v = 0; v < qv.rows(); ++v) if (valence(v) != 4) { irr.add(at(v)); irr.scalar("valence", valence(v)); }
    irr.show("irregular vertices", 0.002, false, "valence");

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

}
#endif
