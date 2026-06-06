//
// Created by saki on 2026/06/06.
//

#include "./tmesh.h"

namespace metriko {

constexpr auto circular_prev = [](auto& c, auto it) { return it == c.begin() ? std::prev(c.end()) : std::prev(it); };
constexpr auto circular_next = [](auto& c, auto it) { auto n = std::next(it); return n == c.end() ? c.begin() : n; };


void Tmesh::collapse_ehalf(int thid) {
    Thalf& eh = thalfs[thid];
    Tquad& eq = tquads[th2quad[eh.id]];
    auto it0 = rg::find(eq., ehid, &Edata::ehid); // edata of ehid
    //auto it_prev = circular_prev(eq.data, it0);
    //auto it_next = circular_next(eq.data, it0);
}
}
