#include "./tmesh_mut.h"
using namespace metriko;

void TmeshMut::collapse_thalf(int thid) {
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

    auto [p_fr, p_to] = [&]{
      if ( collapse_to_prev &&  th_prv.cano) return std::pair{&th_prv.loc_fr(), &th_crr.loc_to()};
      if ( collapse_to_prev && !th_prv.cano) return std::pair{&th_crr.loc_to(), &th_prv.loc_fr()};
      if (!collapse_to_prev &&  th_nxt.cano) return std::pair{&th_crr.loc_fr(), &th_nxt.loc_to()};
      if (!collapse_to_prev && !th_nxt.cano) return std::pair{&th_nxt.loc_to(), &th_crr.loc_fr()};
      throw std::runtime_error("no impl");
    }();

    auto path = approx_shortest_path(20, hm, *p_fr, *p_to, region);

    // fr/to idx already exists, otherwise all inter-points are newly added. follows this condition.
    vec<int> nids;
    int fr_idx = p_fr - tnodes.data();
    int to_idx = p_to - tnodes.data();
    for (int i = 0; i < path.size(); ++i) {
        if      (i == 0)                  nids.push_back(fr_idx);
        else if (i + 1 == path.size())    nids.push_back(to_idx);
        else { tnodes.push_back(path[i]); nids.push_back(tnodes.size() - 1); }
    }

    auto  it_twn_adj = collapse_to_prev ? circular_next(tq_twn.data, it_twn) : circular_prev(tq_twn.data, it_twn);
    auto& th_twn_adj = thalfs[it_twn_adj->thid];
    auto& te_twn_adj = tedges[th_twn_adj.teid];
    assert(it_twn_adj->side == it_twn->side);

    // 1: remove the data from tquads which have collapsed thalfs
    tq_crr.data.erase(it_crr);
    tq_twn.data.erase(it_twn);

    // 2: collapse to the prv/nxt edge
    // 3: insert missing segments to the adjacent edge
    if (collapse_to_prev) { te_prv.nids = nids; if (th_crr.cano) te_twn_adj.insert_locs_after(te_crr.nids); else te_twn_adj.insert_locs_front(te_crr.nids); }
    else                  { te_nxt.nids = nids; if (th_crr.cano) te_twn_adj.insert_locs_front(te_crr.nids); else te_twn_adj.insert_locs_after(te_crr.nids); }
}
