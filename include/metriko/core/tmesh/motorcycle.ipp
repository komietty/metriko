#pragma once
#include "../hmesh/utilities.h"

namespace metriko {
inline const Msgmt &Msgmt::next() const { return curv->sgmts[next_id]; }
inline const Msgmt &Msgmt::prev() const { return curv->sgmts[prev_id]; }

inline void MotorcycleGraph::gen_ports(
    const Hmesh &mesh,
    const VecXc &cfn,
    const VecXi &singular
) {
    std::vector<Mport> joint;
    int counter = 0;
    for (auto v: mesh.verts | vw::filter([&](auto _v){ return singular[_v.id] != 0; })) {
        std::vector<Mport> ports;
        ports.clear();
        for (Half h: v.adjHalfs()) {
            auto a = cfn(h.next().crnr().id);
            auto b = cfn(h.prev().crnr().id);
            auto c = cfn(h.crnr().id);
            auto o = orientation(a, b, c);
            std::vector<Mport> temps;
            if (o < 0) throw std::invalid_argument("Orientation is flipped! This should be ccw order");
            int r;
            for (r = 0; r < 4; r++) { if (!points_into(get_quater_rot(r), a, b, c) && r > 0) break; }
            for (int i = 0; i < 4; i++) {
                complex dir = get_quater_rot(r - i);
                if (points_into(dir, a, b, c) ||
                    std::arg(dir) == std::arg(b - a) ||
                    std::arg(dir) == std::arg(c - a)
                ) {
                    //todo last arg(c - a) is not needed?
                    temps.emplace_back(-1, v, h.face(), a, dir);
                }
            }

            rg::sort(temps, [&](const auto& p0, const auto& p1) {
                complex dir = b - a;
                double dotA = dot(p0.dir, dir);
                double dotB = dot(p1.dir, dir);
                return dotA > dotB;
            });

            for (auto &p: temps) { p.id = counter; counter++; }

            ports.insert(ports.end(), temps.begin(), temps.end());
        }

        for (int i = 0; i < ports.size(); i++) {
            int s = ports.size();
            ports[i].prev = ports[(i - 1 + s) % s].id;
            ports[i].next = ports[(i + 1 + s) % s].id;
        }

        joint.insert(joint.end(), ports.begin(), ports.end());
    }

    mports = joint;
}

//Task: need to consider matching in with fist segments(?)
inline void Mcurv::add_segment_init(
    const VecXc& cfn,
    const std::vector<Mcurv>& edges,
    const VecXi& matching
) {
    Face f = port.face;
    Half h;
    for (const Half _h: f.adjHalfs())
        if (_h.head().id != port.vert.id && _h.tail().id != port.vert.id) h = _h;
    complex uv0 = port.uv;
    complex dir = port.dir;
    add_segment(cfn, edges, *this, dir, uv0, h);
}


//Task: need to consider parallel intersection (e.g. bumpy cube case)
inline void Mcurv::add_segment_next(
    const VecXc& cfn,
    const std::vector<Mcurv>& edges,
    const VecXi& matching
) {
    if (cache.intersected) return;
    auto h   = cache.half.twin();
    auto m   = (h.isCanonical() ? 1 : -1) * matching[h.edge().id];
    auto uv0 = lerp(cfn(h.next().crnr().id), cfn(h.prev().crnr().id), cache.ratio);
    auto dir = std::polar(1., PI / 2 * m) * cache.dir;
    auto oh  = get_oppsite_half(cfn, uv0, dir, h);
    add_segment(cfn, edges, *this, dir, uv0, oh);
}

inline void Mcurv::add_segment(
    const VecXc& cfn,
    const std::vector<Mcurv> &edges,
    const Mcurv& exception,
    const complex dir,
    const complex uv0,
    const Half h
) {
    complex uv1 = cfn(h.prev().crnr().id);
    complex uv2 = cfn(h.next().crnr().id);
    double r_ab, r_cd;
    bool f = find_extended_intersection(uv0, uv0 + dir, uv1, uv2, r_ab, r_cd);
    assert(f);

    complex uv3 = lerp(uv1, uv2, r_cd);

    std::vector<std::tuple<double, double, Msgmt>> candidates;

    auto sgs = vw::all(edges) |
               vw::filter([&](auto &e) { return e != exception; }) |
               vw::transform([](auto &e) { return e.sgmts; }) |
               vw::join |
               vw::filter([&h](auto &sg) { return sg.face.id == h.face().id; });

    for (auto &sg: sgs) {
        double ab = 0;
        double cd = 0;
        if (find_strict_intersection(uv0, uv3, sg.fr.uv, sg.to.uv, ab, cd))
            candidates.emplace_back(ab, cd, sg);
    }

    if (candidates.empty()) {
        auto v1 = Mvert(graph, uv0, nullptr, None);
        auto v2 = Mvert(graph, uv3, nullptr, None, std::optional(h));
        sgmts.emplace_back(graph, this, h.face(), v1, v2);
        cache = Cache(h, r_cd, dir);
    } else {
        auto [ab, cd, sg] = rg::min(candidates, [](auto &a, auto &b) { return std::get<0>(a) < std::get<0>(b); });
        auto v1 = Mvert(graph, uv0, nullptr, None);
        auto v2 = Mvert(graph, lerp(uv0, uv3, ab), sg.curv, Crash);
        sgmts.emplace_back(graph, this, h.face(), v1, v2);
        cache = Cache(h, r_cd, dir, true);
        sg.curv->split_segment(sg, sgmts.back(), cd);
    }
}
}