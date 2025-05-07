#pragma once
namespace metriko {
    inline const Tedge &Thalf::edge() const { return tmesh->tedges[teid]; }
    inline const Thalf &Thalf::twin() const { return tmesh->thalfs[twid]; }
    inline const Thalf &Thalf::next() const { return tmesh->thalfs[tmesh->next_thid(id)]; }
    inline const Thalf &Thalf::prev() const { return tmesh->thalfs[tmesh->prev_thid(id)]; }
    //todo: it was previously twin_id. even changed to id, result does not chage. This might couse something wrong...

    //inline bool Thalf::is_stemming_fr_singular() const { return  cannonical && edge().seg_fr.id == 0; }
    //inline bool Thalf::is_reaching_to_singular() const { return !cannonical && edge().seg_fr.id == 0; }

    //inline complex Thalf::uv_fr() const { return cannonical ? edge().seg_fr.fr.uv : edge().seg_to.to.uv; }
    //inline complex Thalf::uv_to() const { return cannonical ? edge().seg_to.to.uv : edge().seg_fr.fr.uv; }
    inline complex Thalf::uv_fr() const { return cannonical ? edge().uv_fr() : edge().uv_to(); }
    inline complex Thalf::uv_to() const { return cannonical ? edge().uv_to() : edge().uv_fr(); }

    inline Msgmt Thalf::sg_fr() const { return cannonical ? edge().seg_fr : edge().seg_to; }
    inline Msgmt Thalf::sg_to() const { return cannonical ? edge().seg_to : edge().seg_fr; }

    inline std::vector<Thalf> Thalf::adj_thalfs() const {
        std::vector<Thalf> res;
        auto side = cannonical ? sg_to().to.type : sg_to().fr.type;
        auto &next = this->next();
        auto &twin = this->twin();
        switch (side) {
            case HitR:
                res.emplace_back(next);
                res.emplace_back(next.twin().next());
                break;
            case HitB:
            case HitL:
                res.emplace_back(next);
                res.emplace_back(twin.prev().twin());
                break;
            default:;
        }
        return res;
    }
}