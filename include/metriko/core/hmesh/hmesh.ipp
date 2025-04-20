#pragma once

namespace metriko {
    inline Half Half::next() const { return {m->next[id], m}; }
    inline Half Half::prev() const { return {m->prev[id], m}; }
    inline Half Half::twin() const { return {m->twin[id], m}; }
    inline Vert Half::tail() const { return {m->tail[id], m}; }
    inline Vert Half::head() const { return {m->head[id], m}; }
    inline Edge Half::edge() const { return {m->edge[id], m}; }
    inline Face Half::face() const { return {m->face[id], m}; }
    inline Crnr Half::crnr() const { return {m->crnr[id], m}; }

    inline Half Vert::half() const { return {m->vert2half[id], m}; }
    inline Half Edge::half() const { return {m->edge2half[id], m}; }
    inline Half Face::half() const { return {m->face2half[id], m}; }
    inline Half Loop::half() const { return {m->loop2half[id], m}; }

    inline Half Crnr::half() const { return Half{m->crnr2half[id], m}; }
    inline Vert Crnr::vert() const { return Half{m->crnr2half[id], m}.prev().tail(); }
    inline Face Crnr::face() const { return Half{m->crnr2half[id], m}.face(); }

    inline Vert Edge::vert0() const { return {m->edge2vert(id, 0), m}; }
    inline Vert Edge::vert1() const { return {m->edge2vert(id, 1), m}; }
    inline Face Edge::face0() const { return {m->edge2face(id, 0), m}; }
    inline Face Edge::face1() const { return {m->edge2face(id, 1), m}; }

    inline bool Half::isCanonical() const { return edge().half().id == id; }
    inline bool Half::isBoundary()  const { return face().id == -1; }
    inline bool Vert::isBoundary()  const { return m->isBV[id]; }
    inline bool Edge::isBoundary()  const { return half().isBoundary() || half().twin().isBoundary(); }
    inline double Edge::len() const { return Half{m->edge2half[id], m}.len(); }
    inline double Half::len() const { return vec().norm(); }
    inline double Half::cot() const { return m->halfCotan[id]; }
    inline double Edge::cot() const { return m->edgeCotan[id]; }
    inline double Half::varg() const { return m->heArgOnVert[id]; }
    inline double Half::farg() const { return m->heArgOnFace[id]; }
    inline double Half::darg() const { return m->dihedralArg[id]; }
    inline double Face::area() const { return m->faceArea[id]; }
    inline double Vert::baryArea() const { return m->baryDualArea[id]; }
    inline double Vert::circArea() const { return m->circDualArea[id]; }
    inline Row3d Half::vec() const { return head().pos() - tail().pos(); }
    inline Row3d Vert::pos() const { return m->pos.row(id); }
    inline Row3d Vert::basisX() const { return m->vertBasisX.row(id); }
    inline Row3d Vert::basisY() const { return m->vertBasisY.row(id); }
    inline Row3d Vert::normal() const { return m->vertNormal.row(id); }
    inline Row3d Face::basisX() const { return m->faceBasisX.row(id); }
    inline Row3d Face::basisY() const { return m->faceBasisY.row(id); }
    inline Row3d Face::normal() const { return m->faceNormal.row(id); }
    inline Row3d Face::center() const { return m->baryCenter.row(id); }

    inline AdjIter<AdjVH> Vert::adjHalfs(bool ccw) { return {m, m->vert2half[id], ccw}; }
    inline AdjIter<AdjFH> Face::adjHalfs(bool ccw) { return {m, m->face2half[id], ccw}; }
    inline AdjIter<AdjLH> Loop::adjHalfs(bool ccw) { return {m, m->loop2half[id], ccw}; }
    inline AdjIter<AdjVH> Vert::adjHalfs(Half h, bool ccw) { assert(h.tail().id == id); return {m, h.id, ccw}; }
    inline AdjIter<AdjFH> Face::adjHalfs(Half h, bool ccw) { assert(h.face().id == id); return {m, h.id, ccw}; }

    /*
     Created a Double-Connected Edge-List (a.k.a. "halfedge structure") from the usual
     libhedra mesh representation. This data structure is very convenient for mesh editing
     and traversing, and the data structure is again only Eigen vectors and matrices.
     Input:
      D   #F by 1 - face degrees
      F   #F by max(D) - vertex indices in face
      EV  #E by 2 - edge vertex indices
      EF  #E by 2 - edge face indices (EF(i,0) is left face, EF(i,1)=-1 if boundary
      EFi #E by 2 - position of edge in face by EF
      innerEdges vector of inner edges into EV
     Output:
     the number of halfedges can be determined by H=|HV/HE/HF|. It is 2*[Inner edges]+[Boundary Edges]
     VH   #V by 1 - Vertex to outgoing halfedge (into HE)
     EH   #E by 2 - edge to halfedge, where EH(i,0) halfedge is positively oriented, and EH(i,1)=-1 when boundary.
     FH   #F by max(D) - face to (correctly oriented) halfedge s.t. the origin vertex of FH(i,j) is F(i,j)
     HV   #H by 1 - origin vertex of the halfedge
     HE   #H by 1 - edge carrying this halfedge. It does not say which direction.
     HF   #F by 1 - face containing halfedge
     nextH, prevH, twinH - #H by 1 DCEL traversing operations. twinH(i)=-1 for boundary edges.
    */
    inline void dcel(
        const VecXi &D,
        const MatXi &F,
        const MatXi &EV,
        const MatXi &EF,
        const MatXi &EFi,
        VecXi &VH,
        MatXi &EH,
        MatXi &FH,
        VecXi &HV,
        VecXi &HE,
        VecXi &HF,
        VecXi &nextH,
        VecXi &prevH,
        VecXi &twinH
    ) {
        //doing a local halfedge structure for polygonal meshes
        EH = MatXi::Constant(EV.rows(), 2, -1);
        int numH = 0;

        for (int i = 0; i < EF.rows(); i++) {
            if (EF(i, 0) != -1)
                EH(i, 0) = numH++;
            if (EF(i, 1) != -1)
                EH(i, 1) = numH++;
        }


        //halfedges to edge
        HE.conservativeResize(numH);
        for (int i = 0; i < EH.rows(); i++) {
            if (EH(i, 0) != -1)
                HE(EH(i, 0)) = i;
            if (EH(i, 1) != -1)
                HE(EH(i, 1)) = i;
        }

        //halfedge to vertex and vice versa
        HV.conservativeResize(numH);
        VH.conservativeResize(EV.maxCoeff() + 1);
        for (int i = 0; i < EV.rows(); i++) {
            if (EH(i, 0) != -1) {
                HV(EH(i, 0)) = EV(i, 0);
                VH(EV(i, 0)) = EH(i, 0);
            }
            if (EH(i, 1) != -1) {
                HV(EH(i, 1)) = EV(i, 1);
                VH(EV(i, 1)) = EH(i, 1);
            }
        }

        //halfedge to twin
        twinH = Eigen::VectorXi::Constant(numH, -1);
        for (int i = 0; i < EH.rows(); i++)
            if (EH(i, 0) != -1 && EH(i, 1) != -1) {
                twinH(EH(i, 0)) = EH(i, 1);
                twinH(EH(i, 1)) = EH(i, 0);
            }

        //faces to halfedges and vice versa
        FH.resize(F.rows(), F.cols());
        HF.resize(numH);
        for (int i = 0; i < EF.rows(); i++) {
            if (EF(i, 0) != -1) {
                FH(EF(i, 0), EFi(i, 0)) = EH(i, 0);
                HF(EH(i, 0)) = EF(i, 0);
            }
            if (EF(i, 1) != -1) {
                FH(EF(i, 1), EFi(i, 1)) = EH(i, 1);
                HF(EH(i, 1)) = EF(i, 1);
            }
        }

        //halfedge to next and prev
        nextH.conservativeResize(HE.rows());
        prevH.conservativeResize(HE.rows());
        for (int i = 0; i < D.rows(); i++) {
            for (int j = 0; j < D(i); j++) {
                nextH(FH(i, j)) = FH(i, (j + 1) % D(i));
                prevH(FH(i, (j + 1) % D(i))) = FH(i, j);
            }
        }
    }

    inline Hmesh::Hmesh(
        const MatXd &V,
        const MatXi &F
    ) {
        igl::edge_topology(V, F, edge2vert, face2edge, edge2face);
        pos = V;
        idx = F;
        nV = static_cast<int>(pos.rows());
        nF = static_cast<int>(idx.rows());
        nE = static_cast<int>(edge2vert.rows());
        nH = static_cast<int>(edge2vert.rows() * 2);
        nC = nF * 3;
        nEularChars = nV - nE + nF;
        vert2half.assign(nV, -1);
        edge2half.assign(nE, -1);
        face2half.assign(nF, -1);
        crnr2half.assign(nC, -1);
        next.assign(nH, -1);
        prev.assign(nH, -1);
        twin.assign(nH, -1);
        head.assign(nH, -1);
        tail.assign(nH, -1);
        edge.assign(nH, -1);
        face.assign(nH, -1);
        crnr.assign(nH, -1);

        int nP = static_cast<int>(idx.cols());
        auto pair = [&](int a, int b) {
            twin[a] = b;
            twin[b] = a;
        };

        //--- setup topology ---//
        MatXi EFi, EH, FH;
        VecXi VH, HV, HE, HF, nextH, prevH, twinH;
        VecXi D = VecXi::Constant(idx.rows(), 3);
        EFi = MatXi::Constant(edge2face.rows(), 2, -1); // number of an edge inside the face
        for (int i = 0; i < edge2face.rows(); i++) {
            for (int k = 0; k < 2; k++) {
                if (edge2face(i, k) == -1) continue;
                for (int j = 0; j < 3; j++) if (face2edge(edge2face(i, k), j) == i) EFi(i, k) = j;
            }
        }
        dcel(D, idx, edge2vert, edge2face, EFi, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);

        for (int iV = 0; iV < nV; ++iV) vert2half[iV] = VH[iV];
        for (int iE = 0; iE < nE; ++iE) edge2half[iE] = EH.row(iE)[0];
        for (int iF = 0; iF < nF; ++iF) face2half[iF] = FH.row(iF)[0];
        for (int iH = 0; iH < HV.rows(); ++iH) {
            head[iH] = HV[nextH[iH]];
            tail[iH] = HV[iH];
            edge[iH] = HE[iH];
            face[iH] = HF[iH];
            next[iH] = nextH[iH];
            prev[iH] = prevH[iH];
            if (twinH[iH] != -1) twin[iH] = twinH[iH];
        }

        //--- setup boundaries ---//
        for (int iH = 0, iB = nF * nP; iH < nF * nP; ++iH) {
            if (twin[iH] != -1)
                continue;
            int jH = iH;
            int jB = iB;
            while (true) {
                tail[iB] = head[jH];
                head[iB] = tail[jH];
                edge[iB] = edge[jH];
                prev[iB] = iB - 1;
                next[iB] = iB + 1;
                pair(jH, iB);

                jH = prev[jH];
                while (twin[jH] != -1) {
                    if (jH == iH) goto loop_done;
                    jH = prev[twin[jH]];
                }
                iB++;
            }
        loop_done:
            prev[jB] = iB;
            next[iB] = jB;
            loop2half.push_back(jB);
            iB++;
        }
        nL = static_cast<int>(loop2half.size());

        //--- setup edge->half on boundary ---//
        for (int iE = 0; iE < nE; ++iE) if (edge2half[iE] == -1) edge2half[iE] = twin[EH.row(iE)[1]];

        for (int iL = 0; iL < nL; ++iL) loops.emplace_back(Loop{iL, this});
        for (int iF = 0; iF < nF; ++iF) faces.emplace_back(Face{iF, this});
        for (int iV = 0; iV < nV; ++iV) verts.emplace_back(Vert{iV, this});
        for (int iE = 0; iE < nE; ++iE) edges.emplace_back(Edge{iE, this});
        for (int iH = 0; iH < nH; ++iH) halfs.emplace_back(Half{iH, this});
        for (int iC = 0; iC < nC; ++iC) crnrs.emplace_back(Crnr{iC, this});

        //--- setup crnr_half & half_crnr ---//
        for (int iF = 0; iF < nF; ++iF) {
            for (int it = 0; it < 3; ++it) {
                int vid = F(iF, it);
                int cid = iF * 3 + it;
                Half hc = halfs[face2half[iF]];
                for (Half h: {hc, hc.next(), hc.prev()}) {
                    if (h.tail().id != vid && h.head().id != vid) {
                        crnr2half[cid] = h.id;
                        crnr[h.id] = cid;
                    }
                }
            }
        }

        vertNormal.resize(nV, 3);
        vertBasisX.resize(nV, 3);
        vertBasisY.resize(nV, 3);
        faceBasisX.resize(nF, 3);
        faceBasisY.resize(nF, 3);
        faceNormal.resize(nF, 3);
        faceArea.resize(nF);
        halfCotan.resize(nH);
        edgeCotan.resize(nE);
        angleDefect.resize(nV);
        dihedralArg.resize(nH);
        heArgOnVert.resize(nH);
        heArgOnFace.resize(nH);
        baryCenter.resize(nF, 3);
        baryDualArea.resize(nV);
        circDualArea.resize(nV);

        for (Face f: faces) {
            Row3d x = f.half().vec();
            Row3d t = f.half().prev().vec() * -1;
            Row3d n = x.cross(t);
            Row3d p = Row3d::Zero();

            //--- set up  basis of face ---//
            faceBasisX.row(f.id) = x.normalized();
            faceBasisY.row(f.id) = -x.cross(n).normalized();
            faceNormal.row(f.id) = n.normalized();

            //--- set up barycenter of face ---//
            for (Half h: f.adjHalfs()) p += h.tail().pos();
            baryCenter.row(f.id) = p / nP;

            //--- set up area of face ---//
            faceArea[f.id] = n.norm() * 0.5;
        }

        //--- set up halfedge angle on face coordinate ---//
        for (Face f: faces) {
            double sum = 0;
            for (Half h: f.adjHalfs()) {
                Row3d v1 = h.vec().normalized();
                Row3d v2 = h.prev().vec().normalized();
                if ((v1 - f.basisX()).norm() > 1e-10) sum += acos(v1.dot(v2));
                heArgOnFace[h.id] = sum;
            }
        }

        //--- set up vertex orthogonal coordinate ---//
        for (int iF = 0; iF < nF; ++iF) {
            for (int iP = 0; iP < nP; iP++) {
                vertNormal.row(idx(iF, iP)).array()
                        += faceNormal.row(iF).array() * faceArea[iF];
            }
        }
        vertNormal.rowwise().normalize();

        for (Vert v: verts) {
            Row3d v_ = pos.row(v.half().head().id) - pos.row(v.id);
            Row3d n = vertNormal.row(v.id);
            Row3d x = (v_ - v_.dot(n) * n).normalized();
            vertBasisX.row(v.id) = x;
            vertBasisY.row(v.id) = n.cross(x);
        }

        //--- set up a dihedral angle for each halfedge ---//
        for (Half h: halfs) {
            if (h.edge().isBoundary()) {
                dihedralArg[h.id] = 0;
                continue;
            }
            Row3d n1 = faceNormal.row(h.face().id);
            Row3d n2 = faceNormal.row(h.twin().face().id);
            Row3d v = h.vec() / h.len();
            Row3d c = n1.cross(n2);
            dihedralArg[h.id] = atan2(v.dot(c), n1.dot(n2));
        }

        //--- set up halfedge cotan and edge cotan ---//
        for (Half h: halfs) {
            if (h.isBoundary()) {
                halfCotan[h.id] = 0.;
                continue;
            }
            Row3d vn = h.next().vec();
            Row3d vp = h.prev().vec() * -1;
            halfCotan[h.id] = vp.dot(vn) / vp.cross(vn).norm();
        }
        for (Edge e: edges) {
            edgeCotan[e.id] = (e.half().cot() + e.half().twin().cot()) * 0.5;
        }

        //--- set up dual area of vertex ---//
        for (Vert v: verts) {
            double a1 = 0.;
            double a2 = 0.;
            for (Half h: v.adjHalfs()) {
                if (h.face().id != -1) a1 += faceArea[h.face().id];
                a2 += h.cot() * h.vec().squaredNorm();
                a2 += h.prev().cot() * h.prev().vec().squaredNorm();
            }
            baryDualArea[v.id] = a1 / 3.;
            circDualArea[v.id] = a2 * 0.125;
        }

        //--- set up a filter to check whether a vertex is on boundary or not ---//
        isBV.assign(nV, false);
        for (Edge e: edges) {
            if (e.isBoundary()) {
                isBV[e.half().tail().id] = true;
                isBV[e.half().head().id] = true;
            }
        }

        //--- set up halfedge angle on vertex coordinate, and angle defect on vertex ---//
        for (Vert v: verts) {
            double sum = 0;
            for (Half h: v.adjHalfs()) {
                heArgOnVert[h.id] = sum;
                Row3d v1 = h.vec().normalized();
                Row3d v2 = h.prev().twin().vec().normalized();
                sum += acos(v1.dot(v2));
            }
            angleDefect[v.id] = v.isBoundary() ? PI - sum : TwoPI - sum;
            for (Half h: v.adjHalfs()) heArgOnVert[h.id] *= TwoPI / sum;
        }
    }

}
