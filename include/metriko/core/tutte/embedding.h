//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_EMBEDDING_H
#define METRIKO_EMBEDDING_H

namespace metriko {
    class EmbeddedTEdge {
    public:
        int teid;
        std::vector<int> vids;
        std::vector<Half> halfs;
        std::vector<double> vals;
        std::vector<std::pair<int, double>> values;

        EmbeddedTEdge(
            const Hmesh &mesh,
            int teid,
            const std::vector<int> &vids
        ): teid(teid), vids(vids) {
            vals.resize(vids.size(), 0);

            for (int i = 0; i < vids.size() - 1; i++) {
                auto it = rg::find_if(mesh.halfs, [&](Half h) {
                    return h.tail().id == vids[i] && h.head().id == vids[i + 1];
                });
                assert(it != mesh.halfs.end());
                halfs.emplace_back(*it);
            }
        }

        bool contains(Half h) {
            return rg::find(halfs, h) != halfs.end();
        }

        bool contains(Edge e) {
            return rg::find(halfs, e.half()) != halfs.end() ||
                   rg::find(halfs, e.half().twin()) != halfs.end();
        }
    };

    inline int find_closest_point(const VecXc &cfn, const complex uv, const Face f) {
        complex uv0 = cfn[f.half().crnr().id];
        complex uv1 = cfn[f.half().next().crnr().id];
        complex uv2 = cfn[f.half().prev().crnr().id];
        Vert vt0 = f.half().next().head();
        Vert vt1 = f.half().tail();
        Vert vt2 = f.half().head();

        auto pairs = std::vector{
            std::pair(uv0, vt0.id),
            std::pair(uv1, vt1.id),
            std::pair(uv2, vt2.id)
        };

        rg::sort(pairs, [&uv](auto &a, auto &b) { return norm(uv - a.first) < norm(uv - b.first); });
        return pairs[0].second;
    }

    inline void reassign_quantization_value(
        const Hmesh &mesh,
        const int x,
        EmbeddedTEdge &ete
    ) {
        std::vector<double> ls;
        for (int i = 0; i < ete.vids.size() - 1; i++) {
            int vid0 = ete.vids[i];
            int vid1 = ete.vids[i + 1];
            auto p0 = mesh.verts[vid0].pos();
            auto p1 = mesh.verts[vid1].pos();
            ls.emplace_back((p0 - p1).norm());
        }
        double sum = std::accumulate(ls.begin(), ls.end(), 0.);
        double curr = 0;

        for (int i = 1; i < ete.vids.size(); i++) {
            curr += x * ls[i - 1] / sum;
            ete.vals[i] = curr;
        }

        for (int i = 0; i < ete.vids.size(); i++) {
            ete.values.emplace_back(ete.vids[i], ete.vals[i]);
        }
    }

    inline void reassign_quantization_values(
        const Hmesh &mesh,
        const VecXd &X,
        std::vector<EmbeddedTEdge> &etes
    ) {
        for (int i = 0; i < etes.size(); i++)
            reassign_quantization_value(mesh, X[i], etes[i]);
    }

    inline EmbeddedTEdge gen_embedded_tedge(
        const Hmesh &hmesh,
        const VecXc &cfn,
        const Tedge &tedge,
        std::vector<bool> &passthroughs
    ) {
        const int fr_vid = find_closest_point(cfn, tedge.uv_fr(), tedge.seg_fr.face);
        const int to_vid = find_closest_point(cfn, tedge.uv_to(), tedge.seg_to.face);

        if (fr_vid == to_vid) { return EmbeddedTEdge(hmesh, tedge.id, {fr_vid}); }

        std::priority_queue<std::pair<int, int> > q;
        std::unordered_map<int, int> m;
        std::vector<int> visited;

        q.emplace(0, fr_vid);
        visited.emplace_back(fr_vid);

        do {
            int len = q.top().first;
            Vert v = hmesh.verts[q.top().second];
            q.pop();

            if (v.id == to_vid) {
                std::vector<int> res;
                int vid = v.id;
                res.emplace_back(vid);
                while (true) {
                    vid = m[vid];
                    res.emplace_back(vid);
                    if (vid == fr_vid) break;
                    passthroughs[vid] = true;
                }
                return EmbeddedTEdge(hmesh, tedge.id, res);
            }

            for (Half h: v.adjHalfs()) {
                int vid = h.head().id;
                if (rg::none_of(visited, [&](int id) { return id == vid; }) && !passthroughs[vid]) {
                    q.emplace(len - 1, vid);
                    m.emplace(vid, v.id);
                    visited.emplace_back(vid);
                }
            }
        } while (!q.empty());
        throw std::runtime_error("No path found!");
    }

    inline std::optional<EmbeddedTEdge> gen_embedded_tedge_easy(
        const Hmesh &hmesh,
        const VecXc &cfn,
        const Tedge &tedge,
        std::vector<bool> &passthroughs
    ) {
        std::vector<int> vids;
        for (auto &seg: tedge.segments()) {
            Face f = seg.face;
            complex fr = seg.fr.uv;
            complex to = seg.to.uv;
            complex uv0 = cfn[f.half().crnr().id];
            complex uv1 = cfn[f.half().next().crnr().id];
            complex uv2 = cfn[f.half().prev().crnr().id];
            Vert vt0 = f.half().next().head();
            Vert vt1 = f.half().tail();
            Vert vt2 = f.half().head();

            auto pairs = std::vector{
                std::pair(uv0, vt0.id),
                std::pair(uv1, vt1.id),
                std::pair(uv2, vt2.id)
            };

            rg::sort(pairs, [&fr](auto &a, auto &b) { return norm(fr - a.first) < norm(fr - b.first); });
            int vid0 = pairs[0].second;
            if (vid0 == 1852) vid0 = 1859; // todo: temp for icosphere_5!!!!
            if (vid0 == 1851) vid0 = 1858; // todo: temp for icosphere_5!!!!
            if (vid0 == 1860) vid0 = 1852; // todo: temp for icosphere_5!!!!
            if (vids.empty() || vids.back() != vid0) vids.emplace_back(vid0);

            rg::sort(pairs, [&to](auto &a, auto &b) { return norm(to - a.first) < norm(to - b.first); });
            int vid1 = pairs[0].second;
            if (vid1 == 1852) vid1 = 1859; // todo: temp for icosphere_5!!!!
            if (vid1 == 1851) vid1 = 1858; // todo: temp for icosphere_5!!!!
            if (vid1 == 1860) vid1 = 1852; // todo: temp for icosphere_5!!!!
            if (vids.empty() || vids.back() != vid1) vids.emplace_back(vid1);
        }

        if (rg::any_of(vids, [&](int vid) { return passthroughs[vid]; })) return std::nullopt;

        for (int i = 1; i < vids.size() - 1; i++) passthroughs[vids[i]] = true;
        return EmbeddedTEdge(hmesh, tedge.id, vids);
    }


    /*
    inline std::vector<EmbeddedTEdge> gen_embedded_tedges(
        const VecXc& cfn,
        const std::vector<motorcycle::Tedge>& tedges
    ) {
        std::vector<EmbeddedTEdge> etes;
        for (auto& te: tedges) etes.emplace_back(gen_embedded_tedge(cfn, te));
        return etes;
    }
    */
}

#endif
