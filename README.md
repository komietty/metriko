# metriko
![fertility](docs/gargoyle.jpg)

Metriko is a header-only C++ library designed for mesh parameterization and quad meshing. It provides a fast and robust method suitable for meshes containing over 10k vertices.

## Why Metriko?

Quad meshing is a technique used to convert triangular meshes into quadrilateral meshes. Quadrilateral meshes are highly desirable in various applications due to their topological and geometric advantages. They perform particularly well in tasks such as animation deformation, uv mapping, and CAD operations.

Many existing quad meshing algorithms rely heavily on Mixed Integer Programming (MIP) or Mixed Integer Non-Linear Programming (MINLP). These approaches typically require commercial solvers such as Gurobi, making them inaccessible to many users. While alternative algorithms exist, some of them are numerically unstable or consumes several hours of computation, and few have publicly available open-source implementations.

Metriko implements the Quantized Global Parameterization (QGP) algorithm, providing fast and robust quadrilateral parameterization. The algorithm guarantees a valid result and scales linearly with the number of vertices.

![fertility](docs/qgp_1.jpg)

## Current Features and Limitations
Metriko is currently in early development. Available features include:
- Globally optimal rotational symmetry tangent fields
- Basic integration of tangent fields
- Integer Grid Mapping (IGM) using an iterative rounding algorithm
- T-mesh implementation based on the motorcycle graph
- Quantization and mesh generation using T-mesh
- Extraction of a quad mesh(QEx) from a locally-injective uv map

Known limitations and issues:
- Quantization results can vary in quality for coarse quads
- Meshes with boundaries are currently not supported (planned for future support)
- While the quantization is robust, the meshing technique is not yet optimized for coarse quads. The reparameterization scheme described in Lyon et al. (2021) is currently under implementation.

## Usage
Metriko is a header-only library with a minimal dependency on libigl (only the core features are required). To run a demo, install viewer (polyscope), and execute the following commands:
```
git submodule update --init --recursive
cd example
mkdir build && cd build
cmake .. 
make
./metriko_example ../models/icosphere.obj 0.03
```
If SuiteSparse is installed on your machine, it will be automatically integrated via CMake, providing significantly faster matrix computations compared to Eigen.

## References
-  Bommes, D., Zimmer, H., Kobbelt, L. 2009. Mixedinteger quadrangulation. In Proc. SIGGRAPH 2009, 1–10. DOI: [10.1007/978-3-642-11620-9_5](https://doi.org/10.1007/978-3-642-11620-9_5)

- Bommes, D., Campen, M., Ebke, H.-C., Alliez, P., Kobbelt, L. 2013. Integer-grid maps for reliable quad meshing. In Proc. SIGGRAPH 2013, 98:1–98:12. DOI: [10.1145/2461912.2462005](https://doi.org/10.1145/2461912.2462005)

- Campen, M., Bommes, D., Kobbelt, L. 2015. Quantized global parametrization. ACM Trans. Graph 34(6), 192:1–192:12. DOI: [10.1145/2816795.2818140](https://doi.org/10.1145/2816795.2818140)

- Eppstein, D., GOODRICH, M. T., Kim, E., Tamstorf, R. 2008. Motorcycle Graphs: Canonical Quad Mesh Partitioning. Computer Graphics Forum 27, 5, 1477–1486. DOI: [10.1111/j.1467-8659.2008.01288.x](https://doi.org/10.1111/j.1467-8659.2008.01288.x)

- Knöppel, F., Crane, K., Pinkall, U. & Schröder, P. 2013 Globally optimal direction fields. ACM Trans. Graph. 32, 1–10. DOI: [10.1145/2461912.2462005](https://doi.org/10.1145/2461912.2462005)

- Ebke, H.-C., Bommes, D., Campen, M., Kobbelt, L. 2013. QEx: Robust quad mesh extraction. ACM Trans. Graph. 32(6), Article 168, 1–10. DOI: [10.1145/2508363.2508372](https://doi.org/10.1145/2508363.2508372) 

- Lyon, M., Campen, M., Bommes, D., Kobbelt, L. 2019. Parametrization quantization with free boundaries for trimmed quad meshing. ACM Trans. Graph. 38(4), 51:1–51:14. DOI: [10.1145/3306346.3323010](https://doi.org/10.1145/3306346.3323019)

- Lyon, M., Campen, M., Kobbelt, L. 2021. Quad layouts via constrained T-mesh quantization. Comput. Graph. Forum 40(2), 263–275. DOI: [10.1111/cgf.14263]( https://doi.org/10.1111/cgf.142634)

