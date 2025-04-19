#include "face_rosy_field.h"
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>
#include <igl/cut_mesh_from_singularities.h>

namespace metriko {
    std::unique_ptr<FaceRosyField> compute_combbed_field(
            const FaceRosyField& field,
            const std::vector<bool>& seam
    ) {
        auto N = field.rosyN;
        auto F = std::make_unique<FaceRosyField>(field.mesh, field.rosyN);
        F->field.resize(field.mesh.nF, N);
        F->connection = field.connection;
        F->compressed = field.compressed;
        F->singular = field.singular;
        VecXi turns = VecXi::Zero(F->mesh.nF);
        VecXi visit = VecXi::Zero(F->mesh.nF);
        std::queue<std::pair<int, int>> matchingQ;
        matchingQ.emplace(0, 0);

        do {
            std::pair<int, int> pop = matchingQ.front();
            matchingQ.pop();
            auto [i, m] = pop;
            int r = N - m;
            if (visit(i) == 1) continue;
            visit(i) = 1;
            turns(i) = m;
            F->field.block(i, 0, 1, r) = field.field.block(i, m, 1, r);
            F->field.block(i, r, 1, m) = field.field.block(i, 0, 1, m);
            Face face = F->mesh.faces[i];

            for (Half h: face.adjHalfs()) {
                Edge e = h.edge();
                bool b = e.face0().id == i;
                Face fn = b ? e.face1() : e.face0();
                int n = (b ? 1 : -1) * field.matching[e.id];
                n = (n + m + 1000 * N) % N;
                if (!e.isBoundary() && !visit(fn.id) && !seam[e.id]) matchingQ.emplace(fn.id, n);
            }
        } while (!matchingQ.empty());

        F->matching = VecXi::Constant(F->mesh.nE, -1);
        for (Edge e: F->mesh.edges | std::views::filter([](auto e) { return !e.isBoundary();}))
            F->matching[e.id] = (turns[e.face0().id] - turns[e.face1().id] + field.matching[e.id] + 1000 * N) % N;

        return F;
    }

    void cut_mesh_with_singularities(
        const MatXd &V,
        const MatXi &F,
        const std::vector<std::vector<int>> &VF,
        const std::vector<std::vector<int>> &VV,
        const MatXi &TT,
        const MatXi &TTi,
        const VecXi &singularities,
        Eigen::MatrixXi &cuts
    ) {
        //first, get a spanning tree for the mesh (no missmatch needed)
        igl::cut_mesh_from_singularities(V, F, MatXd::Zero(F.rows(), 3).eval(), cuts);

        std::set<int> vertices_in_cut;
        for (int i = 0; i < cuts.rows(); ++i) {
        for (int j = 0; j < cuts.cols(); ++j) {
            if (cuts(i, j)) vertices_in_cut.insert(F(i, j));
        }}

        //then, add all singularities one by one by using Dijkstra's algorithm
        for (int i = 0; i < singularities.rows(); ++i) {
            std::vector<int> path;
            VecXd min_distance;
            VecXi previous;
            int vertex_found = igl::dijkstra(singularities[i], vertices_in_cut, VV, min_distance, previous);
            if (vertex_found == -1)
                // this means that there are no cuts
                path.push_back(singularities[i]);
            else
                igl::dijkstra(vertex_found, previous, path);

            vertices_in_cut.insert(path.begin(), path.end());

            //insert to cut
            for (int ii = 0; ii < path.size() - 1; ++ii) {
                const int &v0 = path[ii];
                const int &v1 = path[ii + 1];

                std::vector<int> vf0 = VF[v0];
                rg::sort(vf0.begin(), vf0.end());
                std::vector<int> vf1 = VF[v1];
                rg::sort(vf1.begin(), vf1.end());

                std::vector<int> common_face_v(std::max(vf0.size(), vf1.size()));
                std::vector<int>::iterator it;
                it = rg::set_intersection(vf0, vf1, common_face_v.begin()).out;
                common_face_v.resize(it - common_face_v.begin());
                assert(common_face_v.size() == 2);

                const int &fi = common_face_v[0];
                int j = -1;
                for (unsigned z = 0; z < 3; ++z)
                    if ((F(fi, z) == v0 && F(fi, (z + 1) % 3) == v1) ||
                        (F(fi, z) == v1 && F(fi, (z + 1) % 3) == v0)) { j = z; }
                assert(j != -1);
                cuts(fi, j) = 1;
                cuts(TT(fi, j), TTi(fi, j)) = 1;
            }
        }
    }

    //Wrapper of the above with only vertices and faces as mesh input
    void cut_mesh_with_singularities(
        const MatXd &V,
        const MatXi &F,
        const VecXi &singularities,
        MatXi &cuts
    ) {
        std::vector<std::vector<int> > VF, VFi;
        std::vector<std::vector<int> > VV;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        igl::adjacency_list(F, VV);
        MatXi TT, TTi;
        igl::triangle_triangle_adjacency(F, TT, TTi);
        cut_mesh_with_singularities(V, F, VF, VV, TT, TTi, singularities, cuts);
    }

    std::vector<bool> compute_seam(const FaceRosyField &f) {
        MatXi face2cut;
        std::vector seam(f.mesh.nE, false);
        std::vector<int> svids;
        for (Vert v: f.mesh.verts)
            if (f.singular[v.id] != 0)
                svids.emplace_back(v.id);
        VecXi s = Eigen::Map<VecXi>(svids.data(), svids.size());
        cut_mesh_with_singularities(f.mesh.pos, f.mesh.idx, s, face2cut);
        for (int iF = 0; iF < f.mesh.nF; iF++) {
            for (int j = 0; j < 3; j++) {
                if (face2cut(iF, j)) seam[f.mesh.face2edge(iF, j)] = true;
            }
        }
        return seam;
    }
}