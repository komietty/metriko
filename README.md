# metriko
![bimba](example/docs/bimba.jpg)

Metriko is a header-only C++ library for quad meshing. It turns a closed triangle mesh into a quad mesh, every quad labelled with the patch of the T-mesh it lies in, in a few seconds for meshes of 10k–30k vertices and well under a minute for 100k+ vertices, without any commercial solver.

## Why Metriko?

Quad meshing is a technique used to convert triangular meshes into quadrilateral meshes. Quadrilateral meshes are highly desirable in various applications due to their topological and geometric advantages. They perform particularly well in tasks such as animation deformation, uv mapping, and CAD operations.

Many existing quad meshing algorithms rely heavily on Mixed Integer Programming (MIP) or Mixed Integer Non-Linear Programming (MINLP). These approaches typically require commercial solvers such as Gurobi, making them inaccessible to many users. While alternative algorithms exist, some of them are numerically unstable or consume several hours of computation, and few have publicly available open-source implementations.

Metriko implements the Quantized Global Parameterization (QGP) algorithm, providing fast and robust quadrilateral parameterization. The algorithm guarantees a valid result and scales linearly with the number of vertices.

## Pipeline
1. Globally optimal 4-rosy tangent field, combed along a seam
2. Integer grid map (IGM) by iterative rounding, locally injective
3. Motorcycle graph of the IGM, giving the T-mesh
4. Quantization of the T-mesh: integer lengths for every arc, no degenerate patch
5. Collapse of the zero-length arcs and snapping of the T-mesh onto mesh vertices
6. Cut along the T-mesh, Tutte embedding of every patch, then SLIM relaxation
7. Quad extraction (QEx) and conformal relaxation of the quads on the input surface
8. Every quad labelled with the T-mesh patch it lies in

## Usage
Meshes with boundaries are currently not supported (planned). Metriko is header-only. It depends on [libigl](https://libigl.github.io/) (core only, which brings Eigen); if [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) and OpenMP are found by CMake they are used automatically and the linear solves get significantly faster than with Eigen alone. The usage is simple: the library offers a single function "compute_quadrangulation".

```cpp
#include "metriko/lib.h"
Eigen::MatrixXd V;   // vertices
Eigen::MatrixXi F;   // triangles
double scale = 0.01; // quad size
auto res = metriko::compute_quadrangulation(V, F, scale);
res->pos; // quad vertices
res->idx; // quad faces
res->val; // patch index of each quad
```

### Demo
The demo in `example/` shows the input and the quads coloured by patch ([polyscope](https://polyscope.run/), as a submodule).
```
git submodule update --init --recursive
cmake -S example -B example/build
cmake --build example/build --target example
./example/build/example example/models/fertility.obj 0.01
```

## Benchmark
Wall time of `compute_quadrangulation`, including the quad relaxation. Apple M5, Apple clang 17 with `-O2`, SuiteSparse and OpenMP. `scale` is the quad edge length as a fraction of the bounding box diagonal.

| mesh | vertices | triangles | scale | quads | time |
|---|---:|---:|---:|---:|---:|
| spot | 2,930 | 5,856 | 0.01 | 7,652 | 1.0 s |
| fandisk | 6,475 | 12,946 | 0.01 | 10,277 | 0.8 s |
| elephant | 12,362 | 24,732 | 0.01 | 5,796 | 3.0 s |
| fertility | 13,971 | 27,954 | 0.01 | 7,885 | 2.1 s |
| bunny | 14,290 | 28,576 | 0.01 | 8,581 | 2.2 s |
| gargoyle | 25,953 | 51,902 | 0.01 | 8,038 | 4.5 s |
| dancer | 12,482 | 24,964 | 0.005 | 6,765 | 2.2 s |
| botijo | 12,441 | 24,898 | 0.005 | 32,691 | 4.0 s |
| rolling stage | 49,988 | 100,000 | 0.005 | 41,298 | 7.4 s |
| igea | 134,345 | 268,686 | 0.005 | 35,914 | 35.1 s |
| bimba | 112,455 | 224,906 | 0.003 | 86,770 | 38.0 s |

## References
- Bommes et al. [Mixed-integer quadrangulation](https://doi.org/10.1145/1531326.1531383). SIGGRAPH 2009.
- Bommes et al. [Integer-grid maps for reliable quad meshing](https://doi.org/10.1145/2461912.2462014). SIGGRAPH 2013.
- Campen et al. [Quantized global parametrization](https://doi.org/10.1145/2816795.2818140). SIGGRAPH Asia 2015.
- Ebke et al. [QEx: Robust quad mesh extraction](https://doi.org/10.1145/2508363.2508372). SIGGRAPH Asia 2013.
- Eppstein et al. [Motorcycle graphs: canonical quad mesh partitioning](https://doi.org/10.1111/j.1467-8659.2008.01288.x). SGP 2008.
- Knöppel et al. [Globally optimal direction fields](https://doi.org/10.1145/2461912.2462005). SIGGRAPH 2013.
- Lyon et al. [Parametrization quantization with free boundaries for trimmed quad meshing](https://doi.org/10.1145/3306346.3323019). SIGGRAPH 2019.
- Lyon et al. [Quad layouts via constrained T-mesh quantization](https://doi.org/10.1111/cgf.142634). Eurographics 2021.
- Myles et al. [Robust field-aligned global parametrization](https://doi.org/10.1145/2601097.2601154). SIGGRAPH 2014.
- Rabinovich et al. [Scalable locally injective mappings](https://doi.org/10.1145/2983621). TOG 2017.
- Bouaziz et al. [Shape-Up: Shaping discrete geometry with projections](https://doi.org/10.1111/j.1467-8659.2012.03171.x). SGP 2012.
