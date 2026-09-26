#ifndef METRIKO_EXAMPLE_COMMON_IO_H
#define METRIKO_EXAMPLE_COMMON_IO_H
// binary caches of the pipeline stages, so that a demo can start from the middle:
//   .cache : per-corner uv of the integer-grid map, matching, singular, seam   (example_0 -> example_1)
//   .em    : the collapsed / snapped Emesh                                      (example_1 -> example_2)
#include <fstream>
#include "metriko/core/tmesh/emesh.h"

using namespace metriko;

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

