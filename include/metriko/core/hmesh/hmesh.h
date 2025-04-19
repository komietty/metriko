//
//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.

#ifndef METRIKO_HMESH_H
#define METRIKO_HMESH_H
#include <uuid/uuid.h>
#include "metriko/core/common/typedef.h"

namespace metriko {

    struct Hmesh;
    struct Half;
    struct Vert;
    struct Edge;
    struct Face;
    struct Crnr;
    struct AdjVH; // Adjacency iterator for verts and halfs
    struct AdjFH; // Adjacency iterator for faces and halfs
    struct AdjLH; // Adjacency iterator for loops and halfs
    template<typename N> struct AdjIter;

    struct Elem {
        int id;
        Hmesh* m;
        bool operator==(const Elem &rhs) const { return id == rhs.id; }
    };

    struct Face : Elem {
        Half half() const;
        Row3d basisX() const;
        Row3d basisY() const;
        Row3d normal() const;
        Row3d center() const;
        double area() const;
        AdjIter<AdjFH> adjHalfs(bool ccw = true);
        AdjIter<AdjFH> adjHalfs(Half h, bool ccw = true);
    };

    struct Edge : Elem {
        Half half() const;
        Vert vert0() const;
        Vert vert1() const;
        Face face0() const;
        Face face1() const;
        double len() const;
        double cot() const;
        bool isBoundary() const;
    };

    struct Vert : Elem {
        Half half() const;
        Row3d pos() const;
        Row3d basisX() const;
        Row3d basisY() const;
        Row3d normal() const;
        bool isBoundary() const;
        double baryArea() const;
        double circArea() const;
        AdjIter<AdjVH> adjHalfs(bool ccw = true);
        AdjIter<AdjVH> adjHalfs(Half h, bool ccw = true);
    };

    struct Crnr : Elem {
        Half half() const;
        Vert vert() const;
        Face face() const;
    };

    struct Loop : Elem {
        Half half() const;
        AdjIter<AdjLH> adjHalfs(bool ccw = true);
    };

    struct Half : Elem {
        Half next() const;
        Half prev() const;
        Half twin() const;
        Vert tail() const;
        Vert head() const;
        Edge edge() const;
        Face face() const;
        Crnr crnr() const;
        double len() const;
        double cot() const;
        double varg() const;
        double farg() const;
        double darg() const;
        bool isBoundary()  const;
        bool isCanonical() const;
        Row3d vec() const;
    };

    struct Hmesh {
        int nV;
        int nF;
        int nE;
        int nH;
        int nL;
        int nC;
        int nEularChars;

        std::vector<int> vert2half;
        std::vector<int> edge2half;
        std::vector<int> face2half;
        std::vector<int> loop2half;
        std::vector<int> crnr2half;
        std::vector<int> next;
        std::vector<int> prev;
        std::vector<int> twin;
        std::vector<int> tail;
        std::vector<int> head;
        std::vector<int> edge;
        std::vector<int> face;
        std::vector<int> crnr;
        std::vector<bool> isBV;
        std::vector<Half> halfs;
        std::vector<Vert> verts;
        std::vector<Edge> edges;
        std::vector<Face> faces;
        std::vector<Loop> loops;
        std::vector<Crnr> crnrs;
        MatXd pos;
        MatXi idx;
        MatXi edge2vert;
        MatXi edge2face;
        MatXi face2edge;
        MatXd vertBasisX;
        MatXd vertBasisY;
        MatXd vertNormal;
        MatXd faceBasisX;
        MatXd faceBasisY;
        MatXd faceNormal;
        VecXd faceArea;
        VecXd baryDualArea;
        VecXd circDualArea;
        VecXd halfCotan;
        VecXd edgeCotan;
        MatXd baryCenter;
        VecXd angleDefect;
        VecXd dihedralArg;
        VecXd heArgOnVert;
        VecXd heArgOnFace;
        Hmesh(const MatXd& V, const MatXi& F);
    };

    struct AdjBase {
        Half h;
        bool bgn;
        bool ccw;
        AdjBase(Half h, bool ccw) : h(h), bgn(false), ccw(ccw) { }
        Half operator*() const { return h; }
        bool operator!=(const AdjBase& a) const { return !bgn || h.id != a.h.id; }
    };

    struct AdjVH : AdjBase {
        AdjVH(Hmesh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
        AdjVH& operator++() { h = ccw ? h.prev().twin() : h.twin().next(); bgn = true; return *this; }
    };

    struct AdjFH: AdjBase {
        AdjFH(Hmesh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
        AdjFH& operator++() { h = ccw ? h.next() : h.prev(); bgn = true; return *this; }
    };

    struct AdjLH : AdjBase {
        AdjLH(Hmesh* m, const int hid, bool ccw): AdjBase(Half{hid, m}, ccw) { }
        AdjLH& operator++() { h = ccw ? h.next() : h.prev(); bgn = true; return *this; }
    };

    template<typename N>
    struct AdjIter {
        AdjIter(Hmesh* m, int iX, bool ccw) : bgnNav(N(m, iX, ccw)), endNav(N(m, iX, ccw)) { }
        N begin() { return bgnNav; }
        N end()   { return endNav; }
        N bgnNav, endNav;
    };
}

#endif
