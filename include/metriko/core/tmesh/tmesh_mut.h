#ifndef EXAMPLE_EBD_CPP_TMESH_MUT_H
#define EXAMPLE_EBD_CPP_TMESH_MUT_H
#include "tmesh.h"
#include "metriko/core/hmesh/hmloc.h"
#include "metriko/core/hmesh/hpath.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko {
struct TmeshMut;

using Terng = std::tuple<int, double, double>; // temp. to define a range for dijkstra. (eid, fr, to)

struct Tqaux {
    vec<std::tuple<HmLoc, double, int>> checkpoints;  // collapse point of tquad. (loc, val, side)
    std::pair<int, vec<int>> side_thids_t;
    std::pair<int, vec<int>> side_thids_b;
    std::pair<int, int> side_thid_l;
    std::pair<int, int> side_thid_r;
    bool thid_l_merge_to_ahead = false;
    bool thid_r_merge_to_ahead = false;
};

struct TedgeMut {
    int id = -1;
    vec<int> nids;
    void insert_locs_front(vec<int> locs) { nids.insert(nids.begin(), locs.begin(), locs.end()); }
    void insert_locs_after(vec<int> locs) { nids.insert(nids.end()  , locs.begin(), locs.end()); }
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
    int thid;
    int side;
};

struct TquadMut {
    int id = -1;
    vec<TdataMut> data;

    vec<int> thids(int side) const {
        return data | vw::filter([&](auto& d) { return d.side == side; })
                    | vw::transform([](auto& d) { return d.thid; })
                    | rg::to<vec<int>>();
    }

    int side_of(const ThalfMut& th) const {
        return rg::find(data, th.id, &TdataMut::thid)->side;
    }
};

struct TmeshMut {
    const Hmesh& hm;
    vec<HmLoc>    tnodes;
    vec<TedgeMut> tedges;
    vec<ThalfMut> thalfs;
    vec<TquadMut> tquads;

    explicit TmeshMut(
        const mc::Mgrph& mg,
        const Tmesh& tm,
        const VecXd& X
    ): hm(mg.hm) {
        const VecXc& cf = mg.cf;

        // 1. nodes: MnodeLoc -> HmLoc
        tnodes.reserve(mg.mnodes.size());
        for (const mc::Mnode& mn : mg.mnodes) {
            tnodes.push_back(std::visit(overloaded{
                [&](const auto&     _) -> HmLoc { throw std::runtime_error("no impl"); },
                [&](const HmLocOnV& v) -> HmLoc { return HmLocOnV{v.id}; },
                [&](const HmLocOnE& e) -> HmLoc { return HmLocOnE{e.id, e.r}; },
                [&](const HmLocOnP& f) -> HmLoc {
                    Face  fc = hm.faces[f.id];
                    Row3d p  = conversion_2d_3d(fc, cf, f.uv);
                    Row3d v  = p - fc.half().tail().pos();
                    return HmLocOnF{f.id, complex(v.dot(fc.basisX()), v.dot(fc.basisY()))};
                },
            }, mn.loc));
        }

        // 2. tedges: segs(Msgmt 列) -> 通過 node id の連鎖
        tedges.reserve(tm.tedges.size());
        for (const Tedge& te : tm.tedges) {
            TedgeMut tem {.id = te.id};
            tem.nids.reserve(te.segs.size() + 1);
            tem.nids.push_back(te.segs.front().fr_nid);
            for (const mc::Msgmt& sg : te.segs) tem.nids.push_back(sg.to_nid);
            tedges.push_back(std::move(tem));
        }

        // 3. thalfs: Tmesh の thalf をミラー。 quad id / bgn / end を half に載せる
        thalfs.reserve(tm.thalfs.size());
        for (const Thalf& th : tm.thalfs) {
            const Tedge& te = tm.tedges[th.teid];
            thalfs.push_back({
                .tm   = this,
                .id   = th.id,
                .twid = th.twid,
                .teid = th.teid,
                .tqid = tm.th2quad[th.id],
                .cano = th.cano,
                .bgn  = th.cano && te.isBgn,
                .end  = th.cano && te.isEnd,
                .x    = X[th.teid],
                .r    = tm.tedges[th.teid].len
            });
        }

        // 4. tquads
        tquads.reserve(tm.tquads.size());
        for (const Tquad& tq : tm.tquads) {
            TquadMut tqm;
            tqm.id = tq.id;
            tqm.data.reserve(tq.data.size());
            for (const Tdata& d : tq.data) tqm.data.push_back({d.thid, d.side});
            tquads.push_back(std::move(tqm));
        }
    }

    void collapse_thalf(int thid);
    bool collapse_tquad_prepare(int tqid, Tqaux& tqaux) const;
    void collapse_tquad_execute(int tqid, Tqaux& tqaux);
    vec<Terng> allowed_range(int tqid) const;        // 領域2彩色版
    vec<Terng> allowed_range_trace(int tqid) const;  // 境界点から横断トレース版（vert 問題回避）
};

inline const HmLoc& ThalfMut::loc_fr() const { const auto& [_, nids] = tm->tedges[this->teid]; return tm->tnodes[cano ? nids.front() : nids.back()]; }
inline const HmLoc& ThalfMut::loc_to() const { const auto& [_, nids] = tm->tedges[this->teid]; return tm->tnodes[cano ? nids.back() : nids.front()]; }
}

#endif
