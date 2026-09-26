#ifndef METRIKO_EXAMPLE_COMMON_VISUALIZER_H
#define METRIKO_EXAMPLE_COMMON_VISUALIZER_H
// polyscope views of every stage of the pipeline: mesh and field, motorcycle t-mesh,
// collapsed / snapped t-mesh, qex output, and the quad patches
#include <set>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "metriko/core/hmesh/utilities.h"
#include "metriko/core/vectorfield/face_rosy_field.h"
#include "metriko/core/tmesh/emesh.h"
#include "metriko/patching.h"

namespace metriko::visualizer {

inline void visualize_init() {
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

template <typename PType, typename IType>
polyscope::SurfaceMesh* visualize_mesh(const PType& pos, const IType& idx, const bool show = true, const std::string& name = "mesh") {
    auto* surf = polyscope::registerSurfaceMesh(name, pos, idx);
    surf->setEdgeWidth(0.7);
    surf->setEnabled(show);
    surf->setMaterial("flat");
    surf->setSurfaceColor(glm::vec3(0.3, 0.3, 0.3));
    return surf;
}

template <typename PType, typename IType, typename UType>
polyscope::SurfaceMesh* visualize_mesh_with_uv(const PType& pos, const IType& idx, const UType& uv, const bool show_hm = true, const bool show_uv = true, const std::string& name = "mesh") {
    auto* surf = visualize_mesh(pos, idx, show_hm, name);
    auto* prms = surf->addParameterizationQuantity("uv_param", uv);
    prms->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    prms->setCheckerSize(1.);
    prms->setEnabled(show_uv);
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

//------------------------------------------------------------------------------
// motorcycle t-mesh (Tmesh over the Mgrph, before the collapse)
//------------------------------------------------------------------------------

// every tedge as its motorcycle segments, with the quantized length X when given
inline void visualize_tedge(const Tmesh& tm, const Mgrph& mg, const VecXc& uv, const VecXd* X = nullptr, const bool show = true) {
    Segments sg;
    for (const auto& te: tm.tedges) {
        const auto& th0 = tm.thalfs[te.id * 2];   // cano thalf: created as a consecutive pair (2*teid, +1 = twin)
        for (const Msgmt& s: te.segs) {
            const Face f = mg.hm.faces[s.face_id];
            complex a = get_face_uv(mg.mnodes[s.fr_nid].loc, s.face_id, mg.hm, mg.cf);
            complex b = get_face_uv(mg.mnodes[s.to_nid].loc, s.face_id, mg.hm, mg.cf);
            sg.add(conversion_2d_3d(f, uv, a), conversion_2d_3d(f, uv, b));
            sg.scalar("teid", te.id);
            sg.scalar("tqid1", tm.th2quad[th0.id]);
            sg.scalar("tqid2", tm.th2quad[th0.twid]);
            sg.scalar("R", te.len);
            if (X) sg.scalar("X", (*X)[te.id]);
        }
    }
    sg.show("tedges", 0.0005, show)->setColor({0., 0., 0.});
}

// tquads whose boundary walk returns to a corner node it already passed. the patch is then not a
// disk, which the per-tquad tutte embedding downstream assumes it is. it happens around
// high-valence singularities: the patch leaves along one separatrix and comes back along another
inline void visualize_pinched_tquads(const Tmesh& tm, const Mgrph& mg, const VecXc& uv, const bool show = true, const double scale = 0.002) {
    Segments sg;
    Points pts;
    for (const auto& tq: tm.tquads) {
        if (tq.id == -1 || tq.data.empty()) continue;
        std::map<int, int> visits;
        for (const auto& d: tq.data) visits[tm.thalfs[d.thid].nid_fr()]++;
        if (rg::none_of(visits, [](const auto& p) { return p.second > 1; })) continue;

        std::set<int> done;
        for (const auto& [thid, side]: tq.data)
        for (const Msgmt& s: tm.tedges[tm.thalfs[thid].teid].segs) {
            const Face f = mg.hm.faces[s.face_id];
            Row3d p = conversion_2d_3d(f, uv, get_face_uv(mg.mnodes[s.fr_nid].loc, s.face_id, mg.hm, mg.cf));
            Row3d q = conversion_2d_3d(f, uv, get_face_uv(mg.mnodes[s.to_nid].loc, s.face_id, mg.hm, mg.cf));
            sg.add(p, q);
            sg.scalar("tqid", tq.id);
            sg.scalar("side", side);
            if (auto it = visits.find(s.fr_nid); it != visits.end() && it->second > 1 && done.insert(s.fr_nid).second) {
                pts.add(p);
                pts.scalar("nid", s.fr_nid);
                pts.scalar("tqid", tq.id);
                pts.scalar("visits", it->second);
            }
        }
    }
    if (!pts.ps.empty()) std::println("[pinched] {} pinch nodes over {} tquads", pts.ps.size(), tm.nTQ);
    sg.show("pinched tquad boundary", scale, show, "side");
    pts.show("pinch node", scale * 2.5, show, "nid");
}

//------------------------------------------------------------------------------
// collapsed / snapped t-mesh (Emesh)
//------------------------------------------------------------------------------

// allowed corridor of a shortest-path query: the admissible sub-segment [r0, r1] of every edge
inline void visualize_allowed_range(const Hmesh& hm, const vec<std::tuple<int, double, double>>& allowed, const std::string& name = "allowed range", const bool show = true, const double scale = 0.0015) {
    Segments sg;
    for (auto& [eid, r0, r1]: allowed) {
        Edge e = hm.edges[eid];
        sg.add(e.lerp(r0), e.lerp(r1));
        sg.scalar("eid", eid);
        sg.scalar("span", r1 - r0);
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
inline void visualize_qports(const Hmesh& hm, const VecXc& uv, const vec<qex::Qport>& ports, const double scale = 0.001, const bool show = true) {
    auto dir_index = [](complex d) { return d.real() > 0.5 ? 0 : d.imag() > 0.5 ? 1 : d.real() < -0.5 ? 2 : 3; };
    Points pts;
    for (const auto& q: ports) {
        pts.add(conversion_2d_3d(hm.faces[q.fid], uv, q.uv + q.dir * 0.15));
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

// a labeled quad mesh: singular / junction nodes, tedge tracks, and the patches coloured by tquad
inline void visualize_quad_patch(const MatXd& qv, const MatXi& qidx, const QuadPatch& patch) {
    auto at = [&](int v) { return Row3d(qv.row(v)); };

    Points sing;
    for (auto& s: patch.singulars) {
        sing.add(at(s.qv));
        sing.scalar("vertex id", s.vid);
        sing.scalar("tedges", s.tedges);
        sing.scalar("tracks", s.tracks);
        sing.scalar("quad valence", s.valence);
        if (s.tracks != s.tedges) std::println("[quad patch] singular vertex {}: {} tedges but {} quad edges on a track (quad valence {})", s.vid, s.tedges, s.tracks, s.valence);
    }
    sing.show("tedge start", 0.002, false)->setPointColor({0.1, 0.8, 0.1});

    Points junc;
    for (int v: patch.junctions) junc.add(at(v));
    junc.show("tedge end", 0.0015, false)->setPointColor({0.9, 0.1, 0.1});

    Segments tracks;
    for (auto [a, b]: patch.track) tracks.add(at(a), at(b));
    tracks.show("tedge tracks on the quad mesh", 0.001, true)->setColor({0., 0., 0.});

    auto* surf = polyscope::registerSurfaceMesh("quad patch", qv, qidx);
    surf->setShadeStyle(polyscope::MeshShadeStyle::Smooth);
    surf->setEdgeWidth(1.);
    surf->addFaceScalarQuantity("tqid", patch.tqid_of_quad);
    auto* col = surf->addFaceScalarQuantity("patch colour", patch.colour);
    col->setColorMap("coolwarm");
    col->setEnabled(true);

    if (!patch.ok())
        std::println("[quad patch] anchors {} | singulars with no consistent rotation {} | with several consistent rotations {} | unreached tedges {} | junction vertices {} | track edges {} | unlabeled quads {}",
                     patch.anchors, patch.no_rotation, patch.several, patch.unreached, (int)patch.junctions.size(), (int)patch.track.size(), patch.unlabeled);
}

}
#endif
