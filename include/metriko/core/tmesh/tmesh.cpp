#include "tmesh.h"

namespace metriko {
    const Tedge &Thalf::edge() const { return tmesh->tedges[teid]; }
    const Thalf &Thalf::twin() const { return tmesh->thalfs[twid]; }
    const Thalf &Thalf::next() const { return tmesh->thalfs[tmesh->next_thid(id)]; }
    const Thalf &Thalf::prev() const { return tmesh->thalfs[tmesh->prev_thid(id)]; }
    //todo: it was previously twin_id. even changed to id, result does not chage. This might couse something wrong...

    bool Thalf::is_stemming_fr_singular() const { return  cannonical && edge().seg_fr.id == 0; }
    bool Thalf::is_reaching_to_singular() const { return !cannonical && edge().seg_fr.id == 0; }

    complex Thalf::uv_fr() const { return cannonical ? edge().seg_fr.fr.uv : edge().seg_to.to.uv; }
    complex Thalf::uv_to() const { return cannonical ? edge().seg_to.to.uv : edge().seg_fr.fr.uv; }

    Msgmt Thalf::sg_fr() const { return cannonical ? edge().seg_fr : edge().seg_to; }
    Msgmt Thalf::sg_to() const { return cannonical ? edge().seg_to : edge().seg_fr; }

    std::vector<Thalf> Thalf::adj_thalfs() const {
        std::vector<Thalf> res;
        auto side = cannonical ? sg_to().to.side : sg_to().fr.side;
        auto &next = this->next();
        auto &twin = this->twin();
        switch (side) {
            case Right:
                res.emplace_back(next);
                res.emplace_back(next.twin().next());
                break;
            case Crash:
            case Left:
                res.emplace_back(next);
                res.emplace_back(twin.prev().twin());
                break;
            default:;
        }
        return res;
    }
}