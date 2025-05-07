//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_THALF_H
#define METRIKO_THALF_H
#include "tedge.h"

namespace metriko {

    class Thalf {
    public:
        Tmesh *tmesh;
        int id;
        int teid;
        int twid;
        bool cannonical;

        Thalf(
            Tmesh *tmesh,
            const int teid,
            const int id,
            const int twid,
            const bool cannonical
        ) : tmesh(tmesh),
            id(id),
            teid(teid),
            twid(twid),
            cannonical(cannonical) { }

        bool operator==(const Thalf& rhs) const {
            return id == rhs.id && twid == rhs.twid && teid == rhs.teid && cannonical == rhs.cannonical;
        }

        const Thalf& twin() const;
        const Thalf& next() const;
        const Thalf& prev() const;
        const Tedge& edge() const;

        complex uv_fr() const;
        complex uv_to() const;
        complex dif_fr() const { return sg_fr().diff() * (cannonical ? 1. : -1.); }
        complex dif_to() const { return sg_to().diff() * (cannonical ? 1. : -1.); }
        Msgmt sg_fr() const;
        Msgmt sg_to() const;

        std::vector<Thalf> adj_thalfs() const;
    };

}

#endif
