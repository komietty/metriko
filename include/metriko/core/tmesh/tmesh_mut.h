#ifndef EXAMPLE_EBD_CPP_TMESH_MUT_H
#define EXAMPLE_EBD_CPP_TMESH_MUT_H
#include "tmesh.h"
#include "metriko/core/hmesh/hmloc.h"
#include "metriko/core/hmesh/utilities.h"

namespace metriko {
struct TmeshMut;

struct TedgeMut { vec<int> nids; };

struct ThalfMut {
    const TmeshMut* tm = nullptr;
    int id   = -1;
    int twid = -1;
    int teid = -1;
    int tqid = -1;
    int nxid = -1;
    int pvid = -1;
    bool cano = false;
    bool bgn  = false;
    bool end  = false;
    double x  = -1;
};

struct TdataMut {
    int thid;
    int side;
};

struct TquadMut {
    int id = -1;
    vec<TdataMut> data;

    vec<int> thids(int side) const {
        return data
            | vw::filter([&](const auto& d) { return d.side == side; })
            | vw::transform([](const auto& d) { return d.thid; })
            | rg::to<vec<int>>();
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
        const Tmesh& tm
    ): hm(mg.hm) {
        const VecXc& cf = mg.cf;

        // 1. nodes: MnodeLoc -> HmLoc (index は mnode id と一致させる)
        //    OnFace は uv 空間の座標なので 3D 経由で面ローカル xy に直す
        tnodes.reserve(mg.mnodes.size());
        for (const mc::Mnode& mn : mg.mnodes) {
            tnodes.push_back(std::visit(overloaded{
                [&](const auto&     _) -> HmLoc { throw std::runtime_error("no impl"); },
                [&](const HmLocOnV& v) -> HmLoc { return HmLocOnV{v.id}; },
                [&](const HmLocOnE& e) -> HmLoc { return HmLocOnE{e.id, e.r}; },
                [&](const HmLocOnP& f) -> HmLoc {
                    Face  fc = hm.faces[f.id];
                    Row3d p  = conversion_2d_3d(fc, cf, f.uv);           // uv -> 3D
                    Row3d v  = p - fc.half().tail().pos();               // 面ローカル原点からの差
                    return HmLocOnF{f.id, complex(v.dot(fc.basisX()), v.dot(fc.basisY()))};
                },
            }, mn.loc));
        }

        // 2. tedges: segs(Msgmt 列) -> 通過 node id の連鎖
        //    seg[i].to_nid == seg[i+1].fr_nid なので「先頭の fr + 各 seg の to」で全 node が得られる
        tedges.reserve(tm.tedges.size());
        for (const Tedge& te : tm.tedges) {
            TedgeMut tem;
            tem.nids.reserve(te.segs.size() + 1);
            tem.nids.push_back(te.segs.front().fr_nid);
            for (const mc::Msgmt& sg : te.segs) tem.nids.push_back(sg.to_nid);
            tedges.push_back(std::move(tem));
        }

        // 3. thalfs: Tmesh の thalf をミラー。quad id / bgn / end を half に載せる
        thalfs.reserve(tm.thalfs.size());
        for (const Thalf& th : tm.thalfs) {
            const Tedge& te = tm.tedges[th.teid];
            thalfs.push_back({
                .tm   = this,
                .id   = th.id,
                .twid = th.twid,
                .teid = th.teid,
                .tqid = tm.th2quad[th.id],
                .nxid = th.nxid,
                .pvid = th.pvid,
                .cano = th.cano,
                .bgn  = th.cano && te.isBgn,
                .end  = th.cano && te.isEnd,
                // .x は quantization 由来。必要なら別途設定（ctor に X が無いので既定 -1）
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
};
}

#endif
