#ifndef EXAMPLE_EBD_CPP_TMESH_MUT_H
#define EXAMPLE_EBD_CPP_TMESH_MUT_H
#include "tmesh.h"
#include "metriko/core/hmesh/hmloc.h"
#include "metriko/core/hmesh/hpath.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko {
struct TmeshMut;

struct Tqaux {
    vec<std::tuple<HmLoc, double, double, int>> checkpoints = {};  // collapse point of tquad. (loc, val, ord, side)
    std::pair<int, vec<int>> side_thids_t = {};
    std::pair<int, vec<int>> side_thids_b = {};
    std::pair<int, int> side_thid_l = {};
    std::pair<int, int> side_thid_r = {};
    bool thid_l_merge_to_ahead = false;
    bool thid_r_merge_to_ahead = false;
};

struct Tqpoint {
    HmLoc loc;
    int    val = -1;
    double ord = -1; // geometric sort key: arc-length (th.r) position along the spine, per-tquad normalized
    int    adj = -1; // adjancy count
    bool   top = false;
};

struct Tqchain {
    vec<Tqpoint> pts = {};
    vec<int> tqids   = {};
    vec<int> bounds  = {};
    vec<int> thids_z = {}; // zero length thids from left to right
    vec<int> thids_t = {};
    vec<int> thids_b = {};
    int      thid_r  = -1;
    int      thid_l  = -1;
};

struct TedgeMut {
    vec<int> nids = {};
    void insert_locs(const vec<int>& locs);
};

struct ThalfMut {
    const TmeshMut* tm = nullptr;
    int id   = -1;
    int twid = -1;
    int teid = -1;
    int tqid = -1;
    bool cano = false;
    bool bgn  = false;
    bool end  = false;
    double x  = -1;
    double r  = -1;
    const HmLoc& loc_fr() const;
    const HmLoc& loc_to() const;
};

struct TdataMut {
    int thid = -1;
    int side = -1;
};

struct TquadMut {
    int id = -1;
    vec<TdataMut> data;

    int side_of(const ThalfMut& th) const { return rg::find(data, th.id, &TdataMut::thid)->side; }

    vec<int> thids(int side) const {
        return data | vw::filter([&](auto& d) { return d.side == side; })
                    | vw::transform([](auto& d) { return d.thid; })
                    | rg::to<vec<int>>();
    }
};

struct TmeshMut {
    const Hmesh& hm;
    vec<HmLoc>    tnodes = {};
    vec<TedgeMut> tedges = {};
    vec<ThalfMut> thalfs = {};
    vec<TquadMut> tquads = {};

    explicit TmeshMut(
        const mc::Mgrph& mg,
        const Tmesh& tm,
        const VecXd& X
    ): hm(mg.hm) {
        const VecXc& cf = mg.cf;
        tnodes.reserve(mg.mnodes.size());
        tedges.reserve(tm.tedges.size());
        thalfs.reserve(tm.thalfs.size());
        tquads.reserve(tm.tquads.size());

        for (const mc::Mnode& mn : mg.mnodes) {
            tnodes.push_back(std::visit(overloaded{
                [&](const auto&     _) -> HmLoc { throw std::runtime_error("no impl"); },
                [&](const HmLocOnV& v) -> HmLoc { return HmLocOnV{v.id}; },
                [&](const HmLocOnE& e) -> HmLoc { return HmLocOnE{e.id, e.r}; },
                [&](const HmLocOnP& f) -> HmLoc {
                    Face fc = hm.faces[f.id];
                    Row3d p = conversion_2d_3d(fc, cf, f.uv);
                    Row3d v = p - fc.half().tail().pos();
                    return HmLocOnF{f.id, complex(v.dot(fc.basisX()), v.dot(fc.basisY()))};
                },
            }, mn.loc));
        }

        for (const Tedge& te : tm.tedges) {
            TedgeMut tem {};
            tem.nids.reserve(te.segs.size() + 1);
            tem.nids.push_back(te.segs.front().fr_nid);
            for (const mc::Msgmt& sg : te.segs) tem.nids.push_back(sg.to_nid);
            tedges.push_back(std::move(tem));
        }

        for (const Thalf& th : tm.thalfs) {
            thalfs.push_back({
                .tm   = this,
                .id   = th.id,
                .twid = th.twid,
                .teid = th.teid,
                .tqid = tm.th2quad[th.id],
                .cano = th.cano,
                .bgn  = th.cano && tm.tedges[th.teid].isBgn,
                .end  = th.cano && tm.tedges[th.teid].isEnd,
                .x    = X[th.teid],
                .r    = tm.tedges[th.teid].len
            });
        }

        for (const auto& [id, data] : tm.tquads) {
            TquadMut tqm;
            tqm.id = id;
            tqm.data.reserve(data.size());
            for (const auto& d : data) tqm.data.push_back({d.thid, d.side});
            tquads.push_back(std::move(tqm));
        }
    }

    int step_next(int thid) const { auto& [_, d] = tquads[thalfs[thid].tqid]; auto it = rg::find(d, thid, &TdataMut::thid); if (it == d.end()) throw std::runtime_error("step_next"); return circular_next(d, it)->thid; };
    int step_prev(int thid) const { auto& [_, d] = tquads[thalfs[thid].tqid]; auto it = rg::find(d, thid, &TdataMut::thid); if (it == d.end()) throw std::runtime_error("step_prev"); return circular_prev(d, it)->thid; };
    int count_adj_tquads(int thid0) const {
        int count = 0, thid = thid0;
        do { ++count; thid = step_next(thalfs[thid].twid); }
        while (thid != thid0 && count <= thalfs.size());
        return count;
    }

    vec<std::tuple<int, double, double>> allowed_range(int tqid) const;

    void collapse_thalf(int thid);
    bool collapse_tquad_chain_prepare(int tqid, Tqchain& chain) const;
    void collapse_tquad_chain_execute(Tqchain& chain);

    vec<int> add_new_path(const vec<HmLoc>& path, int nid0, int nid1) {
        vec<int> nids;
        for (int i = 0; i < path.size(); ++i) {
            if      (i == 0)                  nids.push_back(nid0);
            else if (i + 1 == path.size())    nids.push_back(nid1);
            else { tnodes.push_back(path[i]); nids.push_back(tnodes.size() - 1); }
        }
        return nids;
    }
};

inline const HmLoc& ThalfMut::loc_fr() const { const auto& [nids] = tm->tedges[this->teid]; return tm->tnodes[cano ? nids.front() : nids.back()]; }
inline const HmLoc& ThalfMut::loc_to() const { const auto& [nids] = tm->tedges[this->teid]; return tm->tnodes[cano ? nids.back() : nids.front()]; }

inline void TedgeMut::insert_locs(const vec<int>& locs) {
    int f = locs.front();
    int b = locs.back();
    if      (nids.front() == b) { nids.insert(nids.begin(), locs.begin(), locs.end() - 1); } // prepend [f..b-1]
    else if (nids.back()  == f) { nids.insert(nids.end(),   locs.begin() + 1, locs.end()); } // append  [f+1..b]
    else if (nids.front() == f) { vec<int> s(locs.begin() + 1, locs.end()); rg::reverse(s); nids.insert(nids.begin(), s.begin(), s.end()); } // prepend reverse([f+1..b])
    else if (nids.back()  == b) { vec<int> s(locs.begin(), locs.end() - 1); rg::reverse(s); nids.insert(nids.end(),   s.begin(), s.end()); } // append  reverse([f..b-1])
    else throw std::runtime_error("merge_into: no shared corner");
}
}
#endif
