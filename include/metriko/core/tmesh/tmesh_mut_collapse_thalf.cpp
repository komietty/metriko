#include "./tmesh_mut.h"
#include "metriko/core/hmesh/hpath.h"

namespace metriko {

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };


bool TmeshMut::collapse_thalf(int thid) {
    auto& th_crr = thalfs[thid];
    auto& th_twn = thalfs[th_crr.twid];
    auto& tq_crr = tquads[th_crr.tqid];
    auto& tq_twn = tquads[th_twn.tqid];
    auto  region = allowed_range(tq_crr.id);
    auto  it_crr = rg::find(tq_crr.data, thid     , &TdataMut::thid);
    auto  it_twn = rg::find(tq_twn.data, th_twn.id, &TdataMut::thid);
    auto  it_prv = circular_prev(tq_crr.data, it_crr);
    auto  it_nxt = circular_next(tq_crr.data, it_crr);
    auto& th_prv = thalfs[it_prv->thid];
    auto& th_nxt = thalfs[it_nxt->thid];
    auto& te_crr = tedges[th_crr.teid];
    auto& te_prv = tedges[th_prv.teid];
    auto& te_nxt = tedges[th_nxt.teid];

    assert(!(it_crr->side != it_prv->side && it_crr->side != it_nxt->side));
    assert(!(it_crr->side == it_prv->side && it_crr->side == it_nxt->side));

    bool collapse_to_prev = it_crr->side != it_prv->side;
    bool collapse_to_next = !collapse_to_prev;

    auto [p_fr, p_to] = [&]{
      if (collapse_to_prev &&  th_prv.cano) return std::pair{&th_prv.loc_fr(), &th_crr.loc_to()};
      if (collapse_to_prev && !th_prv.cano) return std::pair{&th_crr.loc_to(), &th_prv.loc_fr()};
      if (collapse_to_next &&  th_nxt.cano) return std::pair{&th_crr.loc_fr(), &th_nxt.loc_to()};
      if (collapse_to_next && !th_nxt.cano) return std::pair{&th_nxt.loc_to(), &th_crr.loc_fr()};
      throw std::runtime_error("no impl");
    }();

    auto path = approx_shortest_path(10, hm, *p_fr, *p_to, region);

    vec<int> nids = {};
    for (const HmLoc& l : path) {
        tnodes.push_back(l);               // todo: remove duplication!
        nids.push_back(tnodes.size() - 1); // todo: if exist, push that index
    }

    // 1: remove the data from tquads which have collapsed thalfs
    tq_crr.data.erase(it_crr);
    tq_twn.data.erase(it_twn);

    if (collapse_to_prev) {
        te_prv.nids = nids; // replace prev tedge data
        auto& th_twn_nxt = thalfs[th_twn.nxid];
        auto& te_twn_nxt = tedges[th_twn_nxt.teid];
        auto it_twn_nxt = rg::find(tq_twn.data, th_twn_nxt.id, &TdataMut::thid);
        assert(it_twn_nxt != tq_twn.data.end());  // mut found in tq_twin
        assert(it_twn_nxt->side == it_twn->side); // not changing side
        if (th_crr.cano) te_twn_nxt.insert_locs_after(te_crr.nids);
        else             te_twn_nxt.insert_locs_front(te_crr.nids);
    } else {
        te_nxt.nids = nids; // replace next tedge data
        auto& th_twn_prv = thalfs[th_twn.pvid];
        auto& te_twn_prv = tedges[th_twn_prv.teid];
        auto  it_twn_prv = rg::find(tq_twn.data, th_twn_prv.id, &TdataMut::thid);
        assert(it_twn_prv != tq_twn.data.end());  // mut found in tq_twin
        assert(it_twn_prv->side == it_twn->side); // not changing side
        if (th_crr.cano) te_twn_prv.insert_locs_front(te_crr.nids);
        else             te_twn_prv.insert_locs_after(te_crr.nids);

    }
    return true;
}
}
