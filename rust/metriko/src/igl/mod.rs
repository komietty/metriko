//--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
//--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
//--- Ports of the libigl v2.6.4 functions metriko calls, one file per libigl file (see rust/PORTING.md).

pub mod edge_topology;
pub mod cut_mesh_from_singularities;
pub mod dijkstra;
pub mod adjacency_list;
pub mod vertex_triangle_adjacency;
pub mod triangle_triangle_adjacency;
pub mod speye;
pub mod slice;
pub mod polar_svd;
pub mod local_basis;
pub mod grad;
pub mod flip_avoiding_line_search;
pub mod doublearea;
pub mod remove_duplicate_vertices;
pub mod per_face_normals;
pub mod aabb;
pub mod umeyama;
