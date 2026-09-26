#ifndef METRIKO_EXAMPLE_COMMON_H
#define METRIKO_EXAMPLE_COMMON_H
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include "glm/glm.hpp"
#include "metriko/core/hmesh/utilities.h"

using namespace metriko;

namespace metriko::visualizer {
inline void visualize_tedge(
    const Tmesh& tm,
    const Mgrph& mg,
    const VecXc& uv,
    const VecXd* X = nullptr,
    const std::vector<int> &selector = std::vector<int>(),
    const std::string& prefix = std::string(""),
    const bool show = true
) {
    vec<glm::vec3> ns;
    vec<std::array<size_t, 2>> es;
    vec<double> teids;
    vec<double> tqid1;
    vec<double> tqid2;
    vec<double> randoms;
    size_t counter = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, 30);

    vec<double> vecX;
    vec<double> vecR;
    vec<double> uvX;
    vec<double> uvY;
    vec<double> difx;
    vec<double> dify;
    vec<double> count;

    for (int i = 0; i < tm.nTE; i++) {
        const auto& te = tm.tedges[i];
        // find two thalf of tedge (created as a consecutive pair: 2*teid = cano, +1 = twin)
        const auto& th0 = tm.thalfs[te.id * 2];   // cano
        const auto& th1 = tm.thalfs[th0.twid];    // twin
        int q1 = tm.th2quad[th0.id];              // tquad on cano side
        int q2 = tm.th2quad[th1.id];              // tquad on twin side

        int random_value = distr(gen);

        if (!selector.empty() && rg::find(selector, i) == selector.end()) continue;

        for (const Msgmt &ts: te.segs) {
            complex uvFr = get_face_uv(mg.mnodes[ts.fr_nid].loc, ts.face_id, mg.hm, mg.cf);
            complex uvTo = get_face_uv(mg.mnodes[ts.to_nid].loc, ts.face_id, mg.hm, mg.cf);
            Row3d p1 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvFr);
            Row3d p2 = conversion_2d_3d(mg.hm.faces[ts.face_id], uv, uvTo);

            ns.emplace_back(p1.x(), p1.y(), p1.z());
            ns.emplace_back(p2.x(), p2.y(), p2.z());
            es.emplace_back(std::array{counter, counter + 1});

            teids.emplace_back(i);
            tqid1.emplace_back(q1);
            tqid2.emplace_back(q2);

            uvX.emplace_back(uvFr.real());
            uvX.emplace_back(uvTo.real());
            uvY.emplace_back(uvFr.imag());
            uvY.emplace_back(uvTo.imag());

            difx.emplace_back(uvTo.real() - uvFr.real());
            dify.emplace_back(uvTo.imag() - uvFr.imag());

            vecR.emplace_back(te.len);
            if (X != nullptr) vecX.emplace_back((*X)[i]);

            randoms.emplace_back(random_value);
            count.emplace_back(counter);
            counter += 2;
        }
    }

    auto c = polyscope::registerCurveNetwork(prefix + "tedges", ns, es);
    c->setMaterial("flat");
    c->setColor(glm::vec4(.0, .0, .0, 1.));
    c->addEdgeScalarQuantity("teid", teids);
    c->addEdgeScalarQuantity("tqid1", tqid1);
    c->addEdgeScalarQuantity("tqid2", tqid2);
    c->addEdgeScalarQuantity("R", vecR);
    if (X != nullptr) c->addEdgeScalarQuantity("X", vecX);
    c->addNodeScalarQuantity("uv x", uvX);
    c->addNodeScalarQuantity("uv y", uvY);
    //c->addEdgeScalarQuantity("dif x", difx);
    //c->addEdgeScalarQuantity("dif y", dify);
    //c->addEdgeScalarQuantity("random", randoms);
    //c->addEdgeScalarQuantity("count", count);
    c->setEnabled(show);
    c->resetTransform();
    c->setRadius(0.0005);
    c->setMaterial("flat");
}

// tquads whose boundary walk returns to a corner node it already passed. the patch is then not a
// disk, which the per-tquad tutte embedding downstream assumes it is. it happens around
// high-valence singularities: the patch leaves along one separatrix and comes back along another
inline void visualize_pinched_tquads(
    const Tmesh& tm,
    const Mgrph& mg,
    const VecXc& uv,
    const bool show = true,
    const double scale = 0.002
) {
    vec<glm::vec3> ns, ps;
    vec<std::array<size_t, 2>> es;
    vec<double> eq, esd, pq, pn, pv;
    size_t c = 0;

    for (const auto& tq: tm.tquads) {
        if (tq.id == -1 || tq.data.empty()) continue;
        std::map<int, int> visits;
        for (const auto& d: tq.data) visits[tm.thalfs[d.thid].nid_fr()]++;
        if (rg::none_of(visits, [](const auto& p) { return p.second > 1; })) continue;

        std::print("[pinched] tquad {}: corner nodes", tq.id);
        for (const auto& d: tq.data) std::print(" {}", tm.thalfs[d.thid].nid_fr());
        std::println("");

        std::set<int> done;
        for (const auto& [thid, side]: tq.data) {
        for (const Msgmt& sg: tm.tedges[tm.thalfs[thid].teid].segs) {
            complex a = get_face_uv(mg.mnodes[sg.fr_nid].loc, sg.face_id, mg.hm, mg.cf);
            complex b = get_face_uv(mg.mnodes[sg.to_nid].loc, sg.face_id, mg.hm, mg.cf);
            Row3d p = conversion_2d_3d(mg.hm.faces[sg.face_id], uv, a);
            Row3d q = conversion_2d_3d(mg.hm.faces[sg.face_id], uv, b);
            ns.emplace_back(p.x(), p.y(), p.z());
            ns.emplace_back(q.x(), q.y(), q.z());
            es.emplace_back(std::array{c, c + 1});
            c += 2;
            eq.emplace_back(tq.id);
            esd.emplace_back(side);
            if (auto it = visits.find(sg.fr_nid);
                it != visits.end() && it->second > 1 && done.insert(sg.fr_nid).second) {
                ps.emplace_back(p.x(), p.y(), p.z());
                pq.emplace_back(tq.id);
                pn.emplace_back(sg.fr_nid);
                pv.emplace_back(it->second);
            }
        }}
    }
    std::println("[pinched] {} pinch nodes over {} tquads", ps.size(), tm.nTQ);

    auto* cn = polyscope::registerCurveNetwork("pinched tquad boundary", ns, es);
    cn->setMaterial("flat");
    cn->setEnabled(show);
    cn->setRadius(scale);
    cn->addEdgeScalarQuantity("side", esd)->setEnabled(true);
    cn->addEdgeScalarQuantity("tqid", eq);
    cn->resetTransform();

    auto* pc = polyscope::registerPointCloud("pinch node", ps);
    pc->setEnabled(show);
    pc->setPointRadius(scale * 2.5);
    pc->addScalarQuantity("nid", pn)->setEnabled(true);
    pc->addScalarQuantity("tqid", pq);
    pc->addScalarQuantity("visits", pv);
}

}

static bool load_cache(const std::string& p, VecXc& uv2, VecXi& matching, VecXi& singular, std::vector<bool>& seam) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    int64_t nu, nm, ns, ne;
    f.read((char*)&nu, 8); f.read((char*)&nm, 8); f.read((char*)&ns, 8); f.read((char*)&ne, 8);
    uv2.resize(nu); matching.resize(nm); singular.resize(ns);
    f.read((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.read((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.read((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(ne);
    f.read(sb.data(), ne);
    seam.assign(sb.begin(), sb.end());
    return (bool)f;
}

static void save_cache(const std::string& p, const VecXc& uv2, const VecXi& matching, const VecXi& singular, const std::vector<bool>& seam) {
    std::ofstream f(p, std::ios::binary);
    int64_t nu = uv2.size(), nm = matching.size(), ns = singular.size(), ne = (int64_t)seam.size();
    f.write((char*)&nu, 8); f.write((char*)&nm, 8); f.write((char*)&ns, 8); f.write((char*)&ne, 8);
    f.write((char*)uv2.data(),      nu * (int64_t)sizeof(complex));
    f.write((char*)matching.data(), nm * (int64_t)sizeof(int));
    f.write((char*)singular.data(), ns * (int64_t)sizeof(int));
    std::vector<char> sb(seam.begin(), seam.end());
    f.write(sb.data(), ne);
}

// binary snapshot of a (collapsed) Emesh, so downstream demos can skip the
// motorcycle graph / quantization / collapse / snap stages
static void save_emesh(const std::string& p, const Emesh& tm) {
    std::ofstream f(p, std::ios::binary);
    auto wi = [&](int v)    { f.write((char*)&v, 4); };
    auto wd = [&](double v) { f.write((char*)&v, 8); };

    wi((int)tm.tnodes.size());
    for (const HmLoc& l: tm.tnodes) {
        wi((int)l.index());   // 0:OnV 1:OnE 2:OnF 3:OnH 4:OnC 5:OnP (variant order)
        std::visit(overloaded{
            [&](const HmLocOnV& v) { wi(v.id); },
            [&](const HmLocOnC& v) { wi(v.id); },
            [&](const HmLocOnE& v) { wi(v.id); wd(v.r); },
            [&](const HmLocOnH& v) { wi(v.id); wd(v.r); },
            [&](const HmLocOnF& v) { wi(v.id); wd(v.xy.real()); wd(v.xy.imag()); },
            [&](const HmLocOnP& v) { wi(v.id); wd(v.uv.real()); wd(v.uv.imag()); },
        }, l);
    }
    wi((int)tm.tedges.size());
    for (auto& te: tm.tedges) {
        wi(te.id);
        wi((int)te.nids.size());
        f.write((char*)te.nids.data(), te.nids.size() * 4);
    }
    wi((int)tm.thalfs.size());
    for (auto& th: tm.thalfs) {
        wi(th.id); wi(th.twid); wi(th.teid); wi(th.tqid);
        wi(th.cano); wi(th.bgn); wi(th.end);
        wd(th.x); wd(th.r);
    }
    wi((int)tm.tquads.size());
    for (auto& tq: tm.tquads) {
        wi(tq.id);
        wi((int)tq.data.size());
        for (auto& d: tq.data) { wi(d.thid); wi(d.side); }
    }
}

// fills a Emesh shell (constructed as Emesh(hm)). do not move the object
// afterwards: the thalf back-pointers are bound here
static bool load_emesh(const std::string& p, Emesh& tm) {
    std::ifstream f(p, std::ios::binary);
    if (!f) return false;
    auto ri = [&]() { int v = 0;    f.read((char*)&v, 4); return v; };
    auto rd = [&]() { double v = 0; f.read((char*)&v, 8); return v; };

    tm.tnodes.clear(); tm.tedges.clear(); tm.thalfs.clear(); tm.tquads.clear();

    for (int n = ri(), i = 0; i < n; ++i) {
        int tag = ri(), id = ri();
        switch (tag) {
            case 0: tm.tnodes.emplace_back(HmLocOnV{id}); break;
            case 1: tm.tnodes.emplace_back(HmLocOnE{id, rd()}); break;
            case 2: { double x = rd(), y = rd(); tm.tnodes.emplace_back(HmLocOnF{id, {x, y}}); break; }
            case 3: tm.tnodes.emplace_back(HmLocOnH{id, rd()}); break;
            case 4: tm.tnodes.emplace_back(HmLocOnC{id}); break;
            case 5: { double x = rd(), y = rd(); tm.tnodes.emplace_back(HmLocOnP{id, {x, y}}); break; }
            default: return false;
        }
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        Eedge te{.id = ri()};
        te.nids.resize(ri());
        f.read((char*)te.nids.data(), te.nids.size() * 4);
        tm.tedges.push_back(std::move(te));
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        Ehalf th{.tm = &tm};
        th.id   = ri(); th.twid = ri(); th.teid = ri(); th.tqid = ri();
        th.cano = ri(); th.bgn  = ri(); th.end  = ri();
        th.x    = rd(); th.r    = rd();
        tm.thalfs.push_back(th);
    }
    for (int n = ri(), i = 0; i < n; ++i) {
        Equad tq{.id = ri()};
        tq.data.resize(ri());
        for (auto& d: tq.data) { d.thid = ri(); d.side = ri(); }
        tm.tquads.push_back(std::move(tq));
    }
    return (bool)f;
}
#endif

