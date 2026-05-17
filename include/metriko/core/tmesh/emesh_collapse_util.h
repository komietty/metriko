#ifndef EXAMPLE_EBD_CPP_EMESH_COLLAPSE_UTIL_H
#define EXAMPLE_EBD_CPP_EMESH_COLLAPSE_UTIL_H
#include "./emesh.h"

namespace metriko {

inline vec<int> verts_inside(
    const Equad& eq,
    const Emesh& em,
    const Hmesh& hm
    ) {
    vec<int> res;
    vec is_boundary(hm.nH, false);
    vec visited_face(hm.nF, false);
    std::queue<int> q;

    for (auto [ehid, side] : eq.data) {
        const Ehalf& eh = em.ehalfs[ehid];
        for (Half h: eh.halfs) {
            is_boundary[h.id] = true;
            int fid = hm.halfs[h.id].face().id;
            if (fid != -1 && !visited_face[fid]) {
                visited_face[fid] = true;
                q.push(fid);
            }
        }
    }

    vec<int> inner_faces;
    while (!q.empty()) {
        int curr_fid = q.front();
        q.pop();
        inner_faces.push_back(curr_fid);

        Face f = hm.faces[curr_fid];
        for (Half he: f.adjHalfs()) {
            if (is_boundary[he.id]) continue;

            int next_fid = he.twin().face().id;
            if (next_fid != -1 && !visited_face[next_fid]) {
                visited_face[next_fid] = true;
                q.push(next_fid);
            }
        }
    }

    vec visited_vert(hm.nV, false);
    for (int fid : inner_faces) {
        Face f = hm.faces[fid];
        for (Half he : f.adjHalfs()) {
            Vert v = he.tail();
            if (!visited_vert[v.id]) {
                visited_vert[v.id] = true;
                res.push_back(v.id);
            }
        }
    }


    std::vector<glm::vec3> pts;
    for (int vid: res) {
        Vert v = hm.verts[vid];
        pts.emplace_back(v.pos().x(), v.pos().y(), v.pos().z());
    }
    auto p = polyscope::registerPointCloud("inside pts of quad: " + std::to_string(eq.id), pts);
    p->setEnabled(false);
    p->setPointRadius(0.002);
    p->resetTransform();

    return res;
}

struct AuxDijkData {
    Vert vert;
    int side;
    double val;
    double cmp;

    // C++20のコンパイラを納得させるための「完全なルール」
    bool operator<(const AuxDijkData& rhs) const noexcept {
        if (val != rhs.val) return val < rhs.val;
        if (cmp != rhs.cmp) return cmp < rhs.cmp;
        if (vert.id != rhs.vert.id) return vert.id < rhs.vert.id;
        return side < rhs.side;
    }

    bool operator==(const AuxDijkData& rhs) const noexcept {
        return val == rhs.val && cmp == rhs.cmp && vert.id == rhs.vert.id && side == rhs.side;
    }
};

}

#endif
