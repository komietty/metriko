// metriko_dump <mesh> <scale> <outdir>
//
// loads the mesh as metriko_web's native cli does (quadpatch.h load_mesh), then runs lib.h's
// compute_quadrangulation_impl (calign = false, refine = true) step by step, calling the same functions in the same
// order, and writes the data of every stage to <outdir> as .npy files plus manifest.json and summary.json (FORMAT.md).
// it also runs compute_quadrangulation_impl itself and checks that both give bitwise the same quads: a drift between
// this file and lib.h makes the tool fail. the summary lines are those of metriko_web's cli.cpp
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <print>
#include <set>
#include <sstream>
#include <string>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/boundary_facets.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include "metriko/lib.h"

using namespace metriko;
namespace fs = std::filesystem;

// ------------------------------------------------------------------------------------------------ npy writer -----

template <class T> struct NpyType;
template <> struct NpyType<double>  { static constexpr const char* descr = "<f8"; };
template <> struct NpyType<int32_t> { static constexpr const char* descr = "<i4"; };
template <> struct NpyType<complex> { static constexpr const char* descr = "<c16"; };
template <> struct NpyType<uint8_t> { static constexpr const char* descr = "|u1"; };

struct Dumper {
    fs::path dir;
    std::string stage;               // prefix of the next files, e.g. "02_rawf"
    vec<std::string> manifest;       // json objects, one per file
    vec<std::pair<std::string, std::string>> meta; // json key -> json value

    template <class T>
    void write(const std::string& name, const vec<size_t>& shape, const T* data) {
        std::string file = stage + "." + name + ".npy";
        std::string sh = "(";
        size_t n = 1;
        for (size_t k = 0; k < shape.size(); ++k) { sh += std::to_string(shape[k]) + (shape.size() == 1 ? "," : (k + 1 < shape.size() ? ", " : "")); n *= shape[k]; }
        sh += ")";
        std::string header = std::format("{{'descr': '{}', 'fortran_order': False, 'shape': {}, }}", NpyType<T>::descr, sh);
        size_t total = 10 + header.size() + 1;
        header += std::string((64 - total % 64) % 64, ' ') + "\n";
        std::ofstream os(dir / file, std::ios::binary);
        os.write("\x93NUMPY\x01\x00", 8);
        const uint16_t hl = static_cast<uint16_t>(header.size());
        os.write(reinterpret_cast<const char*>(&hl), 2);
        os.write(header.data(), header.size());
        if (n) os.write(reinterpret_cast<const char*>(data), n * sizeof(T));
        std::string shj = "[";
        for (size_t k = 0; k < shape.size(); ++k) shj += std::to_string(shape[k]) + (k + 1 < shape.size() ? ", " : "");
        shj += "]";
        manifest.push_back(std::format(R"({{"file": "{}", "dtype": "{}", "shape": {}}})", file, NpyType<T>::descr, shj));
    }

    // vectors
    void put(const std::string& name, const vec<double>& v)  { write<double>(name, {v.size()}, v.data()); }
    void put(const std::string& name, const vec<int32_t>& v) { write<int32_t>(name, {v.size()}, v.data()); }
    void put(const std::string& name, const vec<complex>& v) { write<complex>(name, {v.size()}, v.data()); }
    void put(const std::string& name, const vec<bool>& v) {
        vec<uint8_t> b(v.size());
        for (size_t i = 0; i < v.size(); ++i) b[i] = v[i] ? 1 : 0;
        write<uint8_t>(name, {b.size()}, b.data());
    }
    void put(const std::string& name, const VecXd& v) { write<double>(name, {(size_t)v.size()}, v.data()); }
    void put(const std::string& name, const VecXi& v) { write<int32_t>(name, {(size_t)v.size()}, v.data()); }
    void put(const std::string& name, const VecXc& v) { write<complex>(name, {(size_t)v.size()}, v.data()); }

    // matrices, written row-major (C order)
    template <class S>
    void put_mat(const std::string& name, const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>& m) {
        const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> r = m;
        if constexpr (std::is_same_v<S, int>) write<int32_t>(name, {(size_t)r.rows(), (size_t)r.cols()}, r.data());
        else                                  write<S>(name, {(size_t)r.rows(), (size_t)r.cols()}, r.data());
    }
    void put(const std::string& name, const MatXd& m) { put_mat(name, m); }
    void put(const std::string& name, const MatXi& m) { put_mat(name, m); }
    void put(const std::string& name, const MatXc& m) { put_mat(name, m); }
    // a flat row-major table of `cols` columns
    template <class T>
    void put_rows(const std::string& name, const vec<T>& flat, size_t cols) {
        write<T>(name, {cols ? flat.size() / cols : 0, cols}, flat.data());
    }
    // a column-major sparse matrix as its nonzeros, column by column, rows ascending (compressed storage order)
    void put_sparse(const std::string& name, const SprsD& s) {
        SprsD c = s;
        c.makeCompressed();
        vec<int32_t> rc;
        vec<double> val;
        for (int k = 0; k < c.outerSize(); ++k)
            for (SprsD::InnerIterator it(c, k); it; ++it) { rc.push_back(it.row()); rc.push_back(it.col()); val.push_back(it.value()); }
        put_rows(name + "_rc", rc, 2);
        put(name + "_val", val);
        meta.emplace_back(stage + "." + name + "_shape", std::format("[{}, {}]", s.rows(), s.cols()));
    }
    void set(const std::string& key, const std::string& json) { meta.emplace_back(stage + "." + key, json); }

    void finish(const std::string& summary_json) {
        std::ofstream m(dir / "manifest.json");
        m << "{\n  \"files\": [\n";
        for (size_t i = 0; i < manifest.size(); ++i) m << "    " << manifest[i] << (i + 1 < manifest.size() ? ",\n" : "\n");
        m << "  ],\n  \"meta\": {\n";
        for (size_t i = 0; i < meta.size(); ++i) m << "    \"" << meta[i].first << "\": " << meta[i].second << (i + 1 < meta.size() ? ",\n" : "\n");
        m << "  }\n}\n";
        std::ofstream(dir / "summary.json") << summary_json;
    }
};

static std::string jnum(double x) { return std::format("{:.17g}", x); }

// ------------------------------------------------------------------------------------------------ stage dumps ----

// HmLoc: kind = variant index (0 V, 1 E, 2 F, 3 H, 4 C, 5 P), id, two parameters (E/H: r, 0; F: xy; P: uv; V/C: 0, 0)
static void loc_rows(const HmLoc& l, vec<int32_t>& ki, vec<double>& val) {
    ki.push_back(static_cast<int32_t>(l.index()));
    std::visit(overloaded{
        [&](const HmLocOnV& x) { ki.push_back(x.id); val.push_back(0); val.push_back(0); },
        [&](const HmLocOnE& x) { ki.push_back(x.id); val.push_back(x.r); val.push_back(0); },
        [&](const HmLocOnF& x) { ki.push_back(x.id); val.push_back(x.xy.real()); val.push_back(x.xy.imag()); },
        [&](const HmLocOnH& x) { ki.push_back(x.id); val.push_back(x.r); val.push_back(0); },
        [&](const HmLocOnC& x) { ki.push_back(x.id); val.push_back(0); val.push_back(0); },
        [&](const HmLocOnP& x) { ki.push_back(x.id); val.push_back(x.uv.real()); val.push_back(x.uv.imag()); },
    }, l);
}

static void dump_hmesh(Dumper& d, const Hmesh& hm, bool geometry) {
    d.put("pos", hm.pos);
    d.put("idx", hm.idx);
    d.set("counts", std::format(R"({{"nV": {}, "nF": {}, "nE": {}, "nH": {}, "nL": {}, "nC": {}}})", hm.nV, hm.nF, hm.nE, hm.nH, hm.nL, hm.nC));
    d.put("edge2vert", hm.edge2vert);
    d.put("face2edge", hm.face2edge);
    d.put("edge2face", hm.edge2face);
    d.put("vert2half", hm.vert2half);
    d.put("edge2half", hm.edge2half);
    d.put("face2half", hm.face2half);
    d.put("loop2half", hm.loop2half);
    d.put("crnr2half", hm.crnr2half);
    d.put("next", hm.next);
    d.put("prev", hm.prev);
    d.put("twin", hm.twin);
    d.put("tail", hm.tail);
    d.put("head", hm.head);
    d.put("edge", hm.edge);
    d.put("face", hm.face);
    d.put("crnr", hm.crnr);
    d.put("isBV", hm.isBV);
    if (!geometry) return;
    d.put("vertNormal", hm.vertNormal);
    d.put("vertBasisX", hm.vertBasisX);
    d.put("vertBasisY", hm.vertBasisY);
    d.put("faceNormal", hm.faceNormal);
    d.put("faceBasisX", hm.faceBasisX);
    d.put("faceBasisY", hm.faceBasisY);
    d.put("faceArea", hm.faceArea);
    d.put("baryCenter", hm.baryCenter);
    d.put("baryDualArea", hm.baryDualArea);
    d.put("circDualArea", hm.circDualArea);
    d.put("halfCotan", hm.halfCotan);
    d.put("edgeCotan", hm.edgeCotan);
    d.put("angleDefect", hm.angleDefect);
    d.put("dihedralArg", hm.dihedralArg);
    d.put("heArgOnVert", hm.heArgOnVert);
    d.put("heArgOnFace", hm.heArgOnFace);
}

static void dump_field(Dumper& d, const FaceRosyField& f) {
    d.put("field", f.field);
    d.put("connection", f.connection);
    d.put("compressed", f.compressed);
    d.put("matching", f.matching);
    d.put("singular", f.singular);
}

static void dump_mgrph(Dumper& d, const Mgrph& mg) {
    vec<complex> pdir;
    vec<int32_t> pids;
    for (const Mport& p: mg.mports) { pdir.push_back(p.dir); pids.insert(pids.end(), {p.crnr_id, p.this_id, p.next_id, p.prev_id}); }
    d.put("mports_dir", pdir);
    d.put_rows("mports_ids", pids, 4);

    vec<int32_t> jt, lki, adj_off{0}, adj;
    vec<double> lval;
    for (const Mnode& n: mg.mnodes) {
        jt.push_back(static_cast<int32_t>(n.jt));
        loc_rows(n.loc, lki, lval);
        for (const Row2i& a: n.adj) { adj.push_back(a.x()); adj.push_back(a.y()); }
        adj_off.push_back(adj.size() / 2);
    }
    d.put("mnodes_jt", jt);
    d.put_rows("mnodes_loc_ki", lki, 2);
    d.put_rows("mnodes_loc_val", lval, 2);
    d.put("mnodes_adj_offsets", adj_off);
    d.put_rows("mnodes_adj", adj, 2);

    vec<int32_t> cid, sg_off{0}, sg, bki;
    vec<bool> bend;
    vec<double> bval;
    vec<complex> bdir;
    for (const Mcurv& c: mg.mcurvs) {
        cid.push_back(c.id);
        for (const Msgmt& s: c.sgmts) sg.insert(sg.end(), {s.fr_nid, s.to_nid, s.this_id, s.curv_id, s.face_id});
        sg_off.push_back(sg.size() / 5);
        loc_rows(c.buff.loc, bki, bval);
        bdir.push_back(c.buff.dir);
        bend.push_back(c.buff.end);
    }
    d.put("mcurvs_id", cid);
    d.put("mcurvs_sgmt_offsets", sg_off);
    d.put_rows("mcurvs_sgmts", sg, 5);
    d.put_rows("mcurvs_buff_loc_ki", bki, 2);
    d.put_rows("mcurvs_buff_loc_val", bval, 2);
    d.put("mcurvs_buff_dir", bdir);
    d.put("mcurvs_buff_end", bend);
}

static void dump_emesh(Dumper& d, const Emesh& tm) {
    vec<int32_t> nki;
    vec<double> nval;
    for (const HmLoc& l: tm.tnodes) loc_rows(l, nki, nval);
    d.put_rows("tnodes_ki", nki, 2);
    d.put_rows("tnodes_val", nval, 2);
    d.put("fixed", tm.fixed);

    vec<int32_t> teid, toff{0}, tnids;
    for (const Eedge& e: tm.tedges) { teid.push_back(e.id); tnids.insert(tnids.end(), e.nids.begin(), e.nids.end()); toff.push_back(tnids.size()); }
    d.put("tedges_id", teid);
    d.put("tedges_nids_offsets", toff);
    d.put("tedges_nids", tnids);

    vec<int32_t> hi;
    vec<double> hx;
    for (const Ehalf& h: tm.thalfs) {
        hi.insert(hi.end(), {h.id, h.twid, h.teid, h.tqid, h.cano, h.bgn, h.end});
        hx.insert(hx.end(), {h.x, h.r});
    }
    d.put_rows("thalfs_ids", hi, 7);
    d.put_rows("thalfs_xr", hx, 2);

    vec<int32_t> qid, qoff{0}, qdata;
    for (const Equad& q: tm.tquads) {
        qid.push_back(q.id);
        for (const Edata& e: q.data) { qdata.push_back(e.thid); qdata.push_back(e.side); }
        qoff.push_back(qdata.size() / 2);
    }
    d.put("tquads_id", qid);
    d.put("tquads_data_offsets", qoff);
    d.put_rows("tquads_data", qdata, 2);
}

static void dump_qverts(Dumper& d, const std::string& name, const vec<qex::Qvert>& qs) {
    vec<complex> uv;
    vec<double> pos;
    vec<int32_t> sid;
    for (const auto& q: qs) { uv.push_back(q.uv); pos.insert(pos.end(), {q.pos.x(), q.pos.y(), q.pos.z()}); sid.push_back(q.sid); }
    d.put(name + "_uv", uv);
    d.put_rows(name + "_pos", pos, 3);
    d.put(name + "_sid", sid);
}

static void dump_qports(Dumper& d, const vec<qex::Qport>& qs) {
    vec<int32_t> ids;
    vec<complex> uv, dir;
    vec<double> pos;
    for (const auto& p: qs) {
        ids.insert(ids.end(), {p.idx, p.vid, p.eid, p.fid, p.isConnected, p.next_id, p.prev_id});
        uv.push_back(p.uv);
        dir.push_back(p.dir);
        pos.insert(pos.end(), {p.pos.x(), p.pos.y(), p.pos.z()});
    }
    d.put_rows("qports_ids", ids, 7);
    d.put("qports_uv", uv);
    d.put("qports_dir", dir);
    d.put_rows("qports_pos", pos, 3);
}

// ------------------------------------------------------------------------------------------- input (quadpatch.h) -

struct InputMesh { MatXd V; MatXi F; };

// metriko_web/wasm/quadpatch.h load_mesh, verbatim
static std::expected<InputMesh, std::string> load_mesh(const std::string& path) {
    InputMesh m;
    try {
        if (!igl::read_triangle_mesh(path, m.V, m.F)) return std::unexpected("failed to read the mesh");
        if (path.ends_with(".stl")) {
            MatXd V; VecXi SVI, SVJ;
            igl::remove_duplicate_vertices(m.V, 0., V, SVI, SVJ);
            m.F = m.F.unaryExpr([&](int i) { return SVJ(i); }).eval();
            m.V = V;
        }
    } catch (const std::exception& ex) { return std::unexpected(std::format("failed to read the mesh: {}", ex.what())); }

    if (m.V.rows() == 0 || m.F.rows() == 0) return std::unexpected("the mesh is empty");
    if (m.F.cols() != 3)                    return std::unexpected("the mesh is not a triangle mesh");

    vec<bool> used(m.V.rows(), false);
    for (int i = 0; i < m.F.size(); ++i) used[m.F(i)] = true;
    if (const auto n = rg::count(used, false); n > 0) return std::unexpected(std::format("the mesh has {} unreferenced vertices", n));

    if (!igl::is_edge_manifold(m.F)) return std::unexpected("the mesh is not edge-manifold");
    VecXi B;
    if (!igl::is_vertex_manifold(m.F, B)) return std::unexpected("the mesh is not vertex-manifold");
    MatXi E;
    igl::boundary_facets(m.F, E);
    if (E.rows() > 0) return std::unexpected(std::format("the mesh has {} boundary edges: metriko needs a closed mesh", E.rows()));
    return m;
}

// quadpatch.h patch_boundaries, verbatim
static vec<int> patch_boundaries(const MatXi& qidx, const VecXi& tqid) {
    std::map<std::pair<int, int>, int> edge_quad;
    for (int q = 0; q < qidx.rows(); ++q)
    for (int j = 0; j < 4; ++j) {
        const int a = qidx(q, j), b = qidx(q, (j + 1) % 4);
        auto [it, fresh] = edge_quad.try_emplace(std::minmax(a, b), q);
        if (!fresh && it->second >= 0 && tqid(q) == tqid(it->second)) it->second = -1;   // inside a patch
    }
    vec<int> bnd;
    for (auto& [e, q]: edge_quad) if (q >= 0) { bnd.push_back(e.first); bnd.push_back(e.second); }
    return bnd;
}

template <class M>
static bool same_bits(const M& a, const M& b) {
    return a.rows() == b.rows() && a.cols() == b.cols() && (a.size() == 0 || std::memcmp(a.data(), b.data(), a.size() * sizeof(typename M::Scalar)) == 0);
}

// -------------------------------------------------------------------------------------------------------- main ---

int main(int argc, char** argv) {
    if (argc < 4) { std::println(stderr, "usage: metriko_dump <mesh> <scale> <outdir>"); return 2; }
    const std::string mesh = argv[1];
    const double scale = std::stod(argv[2]);
    Dumper d{fs::path(argv[3])};
    fs::create_directories(d.dir);

    auto m = load_mesh(mesh);
    if (!m) { std::println(stderr, "failed: {}", m.error()); return 1; }
    const MatXd& V = m->V;
    const MatXi& F = m->F;
    const bool calign = false;  // compute_quadrangulation's defaults
    const bool refine = true;

    const auto t0 = std::chrono::steady_clock::now();

    // ---- lib.h compute_quadrangulation_impl, step by step ----
    d.stage = "00_input";
    auto unit = grid_unit(V, scale);
    d.put("V", V);
    d.put("F", F);
    d.set("scale", jnum(scale));
    d.set("unit", jnum(unit));
    d.put("unit", vec<double>{unit});
    d.put("scale", vec<double>{scale});

    auto hm = std::make_unique<Hmesh>(MatXd(V / unit), F);
    d.stage = "01_hm";
    dump_hmesh(d, *hm, true);

    auto rawf = FaceRosyField(*hm, 4, calign ? FieldType::CurvatureAligned : FieldType::Smoothest);
    d.stage = "02_rawf";
    dump_field(d, rawf);

    auto seam = compute_seam(rawf);
    d.stage = "03_seam";
    d.put("seam", seam);

    auto cutm = compute_cut_mesh(*hm, seam);
    d.stage = "04_cutm";
    dump_hmesh(d, *cutm, false);

    auto cmbf = compute_combbed_field(rawf, seam);
    d.stage = "05_cmbf";
    dump_field(d, *cmbf);

    auto extf = compute_extrinsic_field(*cmbf, 4);
    d.stage = "06_extf";
    d.put("extf", extf);

    auto rp = RosyParameterization(*hm, *cutm, extf, cmbf->singular, cmbf->matching, seam, 4, scale);
    rp.localInjectivity = true;
    rp.setup();
    d.stage = "07_rp_setup";
    d.set("N", std::to_string(rp.N));
    d.set("n", std::to_string(rp.n));
    d.put_sparse("vtrans2cut", rp.vtrans2cut);
    d.put_sparse("constraint", rp.constraint);
    d.put_sparse("uncompress", rp.uncompress);
    d.put("he2matching", rp.he2matching);
    d.put("he2transidx", rp.he2transidx);
    d.put("integerIdcs", rp.integerIdcs);
    d.put("singularIdcs", rp.singularIdcs);
    d.put("fixedIdcs", rp.fixedIdcs);
    d.put("fixedVals", rp.fixedVals);
    d.put("lreductor", rp.lreductor);

    const bool integ_ok = rp.integ();
    d.stage = "08_rp_integ";
    d.set("integ_returned", integ_ok ? "true" : "false");
    d.put("fullx", rp.fullx);
    d.put("nfn", rp.nfn);
    d.put("cfn", rp.cfn);

    hm->cfn = VecXc(hm->nF * 3);
    for (const Face f: hm->faces) {
        hm->cfn(f.id * 3 + 0) = complex{rp.cfn(f.id, 0), rp.cfn(f.id, 1)};
        hm->cfn(f.id * 3 + 1) = complex{rp.cfn(f.id, 4), rp.cfn(f.id, 5)};
        hm->cfn(f.id * 3 + 2) = complex{rp.cfn(f.id, 8), rp.cfn(f.id, 9)};
    }
    d.stage = "09_hm_cfn";
    d.put("cfn", hm->cfn);

    auto mg = std::make_unique<Mgrph>(*hm, cmbf->matching, cmbf->singular);
    d.stage = "10_mg";
    dump_mgrph(d, *mg);

    auto tm = std::make_unique<Emesh>(*mg);
    d.stage = "11_tm_init";
    dump_emesh(d, *tm);

    auto X  = compute_quantization(*tm, *mg);
    d.stage = "12_quantization";
    d.put("X", X);

    validate_quantization(*tm, X);
    tm->set_x(X);
    d.stage = "13_tm_setx";
    dump_emesh(d, *tm);

    int collapse_rounds = 0;
    for (int i = 0; i < 100; ++i) {
        collapse_rounds = i + 1;
        for (Ehalf th0 : tm->thalfs) {
            if (th0.id == -1 || th0.x != 0) continue;
            const auto& th1 = th0.twin();
            if (th0.tquad().thids(th0.side()).size() == 1) continue;
            if (th1.tquad().thids(th1.side()).size() == 1) continue;
            tm->collapse_thalf(th0.id);
        }

        for (const auto& [tqid, _] : tm->live_tquads())
            if (Tqaux a; tm->collapse_tquad_prepare(tqid, a)) tm->collapse_tquad_execute(a);

        if (rg::none_of(tm->thalfs, [](const Ehalf& th) { return th.id != -1 && th.x == 0; })) break;
    }
    validate_collapse_done(*tm);
    d.stage = "14_tm_collapsed";
    d.set("rounds", std::to_string(collapse_rounds));
    dump_emesh(d, *tm);

    tm->collapse_tedge_snap(false);
    d.stage = "15_tm_snap0";
    dump_emesh(d, *tm);
    tm->collapse_tedge_snap(true);
    d.stage = "16_tm_snap1";
    dump_emesh(d, *tm);
    for (const auto& [teid, _] : tm->live_tedges()) { tm->collapse_tedge_snap_dedup(teid); }
    d.stage = "17_tm_dedup";
    dump_emesh(d, *tm);

    repair_crossing_tedges(*tm);
    d.stage = "18_tm_repaired";
    dump_emesh(d, *tm);

    ///--- cut the original mesh along the collapsed t-mesh ---///
    vec<bool> seam1;
    VecXi matching1;
    VecXi singular1;
    vec<HalfData> hdata;
    auto hm_emb = compute_embedding_cut_hmesh(*hm, *tm, seam, cmbf->matching, cmbf->singular, seam1, matching1, singular1, hdata);
    d.stage = "19_emb";
    dump_hmesh(d, *hm_emb, true);
    d.put("seam1", seam1);
    d.put("matching1", matching1);
    d.put("singular1", singular1);
    {
        vec<int32_t> hi;
        vec<double> hv;
        for (const HalfData& h: hdata) { hi.insert(hi.end(), {h.half.id, h.thid, h.tqid, h.order}); hv.insert(hv.end(), {h.v0, h.v1}); }
        d.put_rows("hdata_ids", hi, 4);
        d.put_rows("hdata_vals", hv, 2);
    }

    auto hm_cut = compute_cut_mesh(*hm_emb, seam1);
    d.stage = "20_hcut";
    dump_hmesh(d, *hm_cut, false);

    MatXd uv_ebd = compute_tutte_parameterization(*hm_emb, *tm, seam1, hdata);
    MatXd uv_cut(hm_cut->nV, 2);
    for (auto v: hm_cut->verts) uv_cut.row(v.id) = uv_ebd.row(v.half().next().crnr().id);
    d.stage = "21_tutte";
    d.put("uv_ebd", uv_ebd);
    d.put("uv_cut", uv_cut);

    MatXd uv = seamless_slim(*hm_emb, *hm_cut, seam1, singular1, uv_cut, slim_max_iter, slim_rel_tol);
    d.stage = "22_slim";
    d.put("uv", uv);

    // ------ qex on the slim result ------
    hm_emb->cfn = VecXc(hm_emb->nF * 3);
    for (int i = 0; i < hm_cut->nF; ++i) {
    for (int j = 0; j < 3; ++j) {
        int k = hm_cut->idx(i, j);
        hm_emb->cfn(i * 3 + j) = complex(uv(k, 0), uv(k, 1));
    }}
    d.stage = "23_qex_cfn";
    d.put("cfn", hm_emb->cfn);

    qex::sanitization(*hm_emb, matching1, singular1, 4);
    d.stage = "24_qex_sanitized";
    d.put("cfn", hm_emb->cfn);

    vec<qex::Qport> q_ports;
    vec<qex::Qvert> vqvs, eqvs, fqvs;
    qex::generate_q_vert(*hm_emb, vqvs, eqvs, fqvs);
    d.stage = "25_qverts";
    dump_qverts(d, "v", vqvs);
    dump_qverts(d, "e", eqvs);
    dump_qverts(d, "f", fqvs);

    qex::generate_vqvert_qport(*hm_emb, vqvs, q_ports);
    qex::generate_eqvert_qport(*hm_emb, eqvs, q_ports);
    qex::generate_fqvert_qport(*hm_emb, fqvs, q_ports);
    d.stage = "26_qports";
    dump_qports(d, q_ports);

    auto q_edges = qex::generate_q_edge(*hm_emb, matching1, q_ports);
    d.stage = "27_qedges";
    {
        vec<int32_t> pe;
        for (const auto& e: q_edges) { pe.push_back(&e.port1 - q_ports.data()); pe.push_back(&e.port2 - q_ports.data()); }
        d.put_rows("qedges", pe, 2);
    }
    dump_qports(d, q_ports);   // generate_q_edge connects the ports

    auto q_faces = qex::generate_q_faces(q_ports, q_edges);
    d.stage = "28_qfaces";
    {
        vec<int32_t> off{0}, qh;
        for (const auto& f: q_faces) {
            for (const auto& h: f.qhalfs) qh.insert(qh.end(), {static_cast<int32_t>(&h.qe - q_edges.data()), h.idx, h.cano});
            off.push_back(qh.size() / 3);
        }
        d.put("qfaces_offsets", off);
        d.put_rows("qfaces_qhalfs", qh, 3);
    }

    auto [q_pos, q_idx] = extract_quad_mesh(*hm, q_faces, refine);
    auto q_val = label_quad_patches(*tm, cmbf->singular, q_faces, q_idx);
    d.stage = "29_quads";
    d.put("q_pos_unscaled", q_pos);
    q_pos *= unit;
    d.put("q_pos", q_pos);
    d.put("q_idx", q_idx);
    d.put("q_val", q_val);

    const double ms_steps = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

    // ---- the library's own path must give the same quads ----
    const auto t1 = std::chrono::steady_clock::now();
    auto ref = compute_quadrangulation_impl(V, F, scale);
    const double ms_impl = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    const bool same = same_bits(ref.q_pos, q_pos) && same_bits(ref.q_idx, q_idx) && same_bits(ref.q_val, MatXi(q_val));
    if (!same) {
        std::println(stderr, "drift: the step-by-step pipeline differs from compute_quadrangulation_impl (tools/dump/dump.cpp is out of date with lib.h)");
        return 3;
    }

    // ---- metriko_web cli.cpp summary ----
    const VecXi tqid = q_val.col(0);
    const vec<int> bnd = patch_boundaries(q_idx, tqid);
    std::set<int> ps;
    for (int t: tqid) if (t >= 0) ps.insert(t);
    double sum = 0;
    for (int i = 0; i < q_pos.rows(); ++i) for (int j = 0; j < 3; ++j) sum += q_pos(i, j);
    uint32_t h = 2166136261u;   // fnv-1a over the quad corners and the patch ids
    auto mix = [&](uint32_t x) { for (int b = 0; b < 4; ++b) { h ^= (x >> (8 * b)) & 0xff; h *= 16777619u; } };
    for (int i = 0; i < q_idx.rows(); ++i) for (int j = 0; j < 4; ++j) mix(q_idx(i, j));
    for (int t: tqid) mix(static_cast<uint32_t>(t));

    const std::string l1 = std::format("input   {} vertices, {} triangles", V.rows(), F.rows());
    const std::string l2 = std::format("quads   {} vertices, {} quads, {} patches, {} boundary edges", q_pos.rows(), q_idx.rows(), ps.size(), bnd.size() / 2);
    const std::string l3 = std::format("pos sum {:.16e}", sum);
    const std::string l4 = std::format("hash    {:08x}", h);
    std::println("{}\n{}\n{}\n{}\ntime    {:.0f} ms (step by step, with dumps) / {:.0f} ms (compute_quadrangulation_impl)", l1, l2, l3, l4, ms_steps, ms_impl);

    const std::string summary = std::format(
        "{{\n  \"mesh\": \"{}\",\n  \"scale\": {},\n  \"unit\": {},\n  \"n_input_vertices\": {},\n  \"n_input_triangles\": {},\n"
        "  \"n_quad_vertices\": {},\n  \"n_quads\": {},\n  \"n_patches\": {},\n  \"n_boundary_edges\": {},\n"
        "  \"pos_sum\": \"{:.16e}\",\n  \"hash\": \"{:08x}\",\n  \"ms_impl\": {:.0f},\n  \"ms_steps\": {:.0f},\n"
        "  \"lines\": [\"{}\", \"{}\", \"{}\", \"{}\"]\n}}\n",
        fs::path(mesh).filename().string(), jnum(scale), jnum(unit), V.rows(), F.rows(),
        q_pos.rows(), q_idx.rows(), ps.size(), bnd.size() / 2, sum, h, ms_impl, ms_steps, l1, l2, l3, l4);
    d.finish(summary);
    return 0;
}
