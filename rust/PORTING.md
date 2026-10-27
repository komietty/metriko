# Porting rules: metriko C++ → Rust

The Rust crate `rust/metriko` is a faithful port of the header-only C++ library in `include/metriko`. The C++ is
the specification. `rust/PORTED_FROM` holds the C++ commit the Rust tree mirrors; `git diff $(cat rust/PORTED_FROM)
HEAD -- include/` lists what is left to port.

Faithful means: the same procedure, the same formulas, the same data structures, the same order of operations and
the same discrete decisions (sort orders, tie breaks, pivot choices, branch conditions). Never replace a step by an
approximation. When a 1:1 port is impossible, write `// DEVIATION: <what differs and why>` at the spot and list it
in the report of the work.

## 1. Files and names

- One Rust file per C++ file, same directory, same stem: `include/metriko/tmesh/emesh_collapse_th.ipp` →
  `rust/metriko/src/tmesh/emesh_collapse_th.rs`. A `.ipp` holding member functions becomes an `impl` block of the
  same type in its own file. `lib.h` → `lib.rs`. `common/loop.h` → `common/loop.rs` (`pub mod r#loop`).
- Keep the order of functions of the C++ file, and its comments (translated to Rust comment syntax).
- Names: types keep their C++ names (`Hmesh`, `Emesh`, `Ehalf`, `Mgrph`, `Qport`, ...). Functions, methods and fields
  become snake_case (`adjHalfs` → `adj_halfs`, `isBoundary` → `is_boundary`, `nV` → `n_v`, `cr_t` stays).
- Header of every file:
  ```rust
  //--- Copyright (C) 2025 Saki Komikado <komietty@gmail.com>,
  //--- This Source Code Form is subject to the terms of the Mozilla Public License v.2.0.
  //--- Port of include/metriko/<path> (see rust/PORTED_FROM).
  ```
  Third-party ports name their source instead, e.g. `//--- Port of libigl v2.6.4 include/igl/edge_topology.cpp
  (MPL-2.0).`, `//--- Port of Eigen 5.0.1 Eigen/src/LU/FullPivLU.h (MPL-2.0).`,
  `//--- Port of libc++ 20.1 <__algorithm/sort.h> (Apache-2.0 WITH LLVM-exception).`

## 2. Types

| C++ | Rust |
|---|---|
| `complex` (`std::complex<double>`) | `Complex` = `num_complex::Complex<f64>` (same type as `faer::c64`) |
| `Row3d`, `Vec3d` | `glam::DVec3` (`common::typedef::Row3d`) |
| `Row2d`, `Vec2d` | `glam::DVec2` |
| `Row2i` | `[i32; 2]` |
| `MatXd` used as an array of rows (V, positions, uv per vertex) | `Vec<DVec3>` / `Vec<DVec2>` / `Vec<[f64; N]>` |
| `MatXi` used as an array of rows (F, quads) | `Vec<[i32; 3]>` / `Vec<[i32; 4]>` |
| `MatXd` used as a matrix in linear algebra | `faer::Mat<f64>` |
| `VecXd` / `VecXc` in linear algebra | `faer::Col<f64>` / `faer::Col<Complex>` |
| `VecXd` / `VecXi` / `VecXb` as plain per-element arrays | `Vec<f64>` / `Vec<i32>` / `Vec<bool>` |
| `SprsD` / `SprsC` | `SprsD` = `SparseColMat<usize, f64>`, `SprsC` = `SparseColMat<usize, Complex>` (solver::matrix_ops) |
| `TripD` + `setFromTriplets` | `Triplet<usize, usize, f64>` + `matrix_ops::from_triplets` (duplicates are summed, as in Eigen) |
| element ids (`int`) | `i32`, `-1` keeps meaning "none / dead" exactly as in C++. Index with `as usize`. |
| `std::optional<T>` | `Option<T>`; `.value()` → `.expect("...")` (a panic is the C++ `bad_optional_access`) |
| `std::variant<...>` | `enum` with the same alternatives in the same order |
| `std::set` / `std::map` | `BTreeSet` / `BTreeMap` (same ordering for tuples / pairs) |
| `std::unordered_map` / `set` | `HashMap` / `HashSet` **only for lookups**. If the C++ iterates one and the result depends on the order, stop and report it (libc++ hash order is not reproducible by HashMap). |

Eigen-style dense vectors are not used for 2D points: parameter space points are `Complex` as in C++.

## 3. Mesh handles

`Hmesh` stores the same SoA arrays as C++. Element handles are `Copy` structs borrowing the mesh:
`Half<'a> { id: i32, m: &'a Hmesh }`, likewise `Vert`, `Edge`, `Face`, `Crnr`, `Loop`. Equality compares ids only
(C++ `Elem::operator==`). `hm.halfs()`, `hm.faces()`, ... iterate handles; `hm.half(id)` makes one. Adjacency
iterators (`adj_halfs`, `adj_crnrs`) reproduce the C++ begin/end semantics exactly (they stop when the walk returns
to its start). Mutating a mesh (e.g. `cfn`) happens through `&mut Hmesh` while no handle is alive.

`Emesh` / `Mgrph` hold `&'a Hmesh`. Their elements (`Ehalf`, `Equad`, `Mcurv`, ...) are plain data with ids; the C++
navigation methods that went through the back pointer take the container: `th.twin()` → `th.twin(&tm)`,
`th.tquad()` → `th.tquad(&tm)`, `th.side()` → `th.side(&tm)`, `th.nid_fr()` → `th.nid_fr(&tm)`, `th.loc_fr()` →
`th.loc_fr(&tm)`. `Mcurv::add_segment` is `Mgrph::add_segment(k)`. Loops that mutate the container while walking it
(`for (Ehalf th0 : tm->thalfs)`, `for (... : tm->live_tquads())`) iterate by index over the length taken at the
start, copy the element when it is reached (as the C++ range-for copies / re-reads it), and skip dead ones.

## 4. Errors

C++ throws (`METRIKO_FAIL`, `METRIKO_CHECK`, `optional::value()`, `std::invalid_argument`) and only `lib.h`'s
`compute_quadrangulation` catches. Rust mirrors this: `metriko_fail!(...)` / `metriko_check!(cond, ...)` panic with
the same message (format strings translated from `std::format` to `format!`), and the public
`compute_quadrangulation` turns the panic into `Err(String)` with `std::panic::catch_unwind`. Do not convert internal
functions to `Result`.

## 5. Numerics (bit-for-bit where the C++ is)

- Never use `mul_add` / fused operations unless the C++ does (the reference build uses `-ffp-contract=off`).
- Keep the association order of every sum and product as written in C++ (`a + b + c` is `(a + b) + c`).
- Eigen `v.normalized()` is `v / v.norm()` (componentwise division, unchanged when the norm is 0):
  use `common::typedef::normalized(v)`, never glam's `normalize()` (it multiplies by the reciprocal).
- Eigen `v.norm()` = `sqrt(squaredNorm)`, `squaredNorm` = `((x*x + y*y) + z*z)` = glam `length_squared()`;
  glam `length()` is `sqrt(dot)` and matches. `cross` / `dot` of glam match Eigen's 3D formulas.
- `std::complex` ↔ `num_complex`: `*`, `+`, `-`, complex / real, complex * real, `std::abs` (= `hypot`, `norm()`),
  `std::arg` (= `atan2`, `arg()`), `std::norm` (= `norm_sqr()`) give identical bits.
  **complex / complex is NOT num_complex's `/`**: libc++ lowers it to compiler-rt `__divdc3`, precompiled with FMA
  contraction (plain division matches only 57% of the golden connection values). Always write
  `common::stdcomplex::std_div(a, b)` (also for `a /= b`). Likewise `std::pow` / `std::exp` / `std::log` /
  `std::polar` of complex → `stdcomplex::{std_pow, std_exp, std_log, std_polar}`.
- Eigen 5.0.1 on arm64 fuses where its packet math says so, whatever `-ffp-contract`: small fixed-size products
  with a plain left-hand side go through NEON `pmadd` (a fused chain over k → `mul_add`), a transposed left-hand
  side (`A^T*A`, `A^T*b`) sums `(l0r0+l2r2)+(l1r1+l3r3)` unfused, `Matrix2d::squaredNorm` sums pairwise, and the
  `norm()` of a dynamic vector uses 2-lane packets with two accumulators (see `eigen_squared_norm` in quantization).
  Reproduce the evaluation Eigen actually performs and verify it bitwise with a C++ harness (see igl/polar_svd.rs).
- Eigen complex dot `a.dot(b)` is `sum conj(a_i) * b_i`. Eigen's vectorized reductions of long real vectors sum in a
  different order than a plain loop: this only changes last bits, keep a plain loop and note nothing.
- `std::floor`, `std::lround`, `std::acos`, `std::atan2`, ... map to the f64 methods (`round` for `lround`:
  `lround` rounds half away from zero, as `f64::round`).

## 6. Order-sensitive algorithms (`common::stdcxx`)

C++ results depend on libc++ 20.1 behavior wherever ties exist. Use:
- `stdcxx::sort_by(&mut slice, |a, b| a_less_than_b)` for every `std::sort` / `std::ranges::sort` with a custom
  comparator or projection (a port of libc++'s introsort; unstable in exactly the same way). A projection `proj`
  becomes the comparator `proj(a) < proj(b)`. Sorting ints / tuples with the default `<` may use `sort_unstable`
  (equal elements are identical).
- `std::stable_sort` → Rust `sort_by` (any stable sort gives the same result).
- `std::priority_queue` → `stdcxx::PriorityQueue` (libc++ push_heap / pop_heap port) whenever the comparator can tie.
  A `priority_queue<pair<double,int>, ..., greater<>>` pops the smallest `(d, id)` pair; ties are identical pairs, so
  a `BinaryHeap<Reverse<(OrdF64, i32)>>` is also exact there.
- `rg::min` / `rg::min_element` return the first minimum (Rust `min_by` does too). `rg::max_element` returns the
  **first** maximum but Rust `max_by` returns the last: use `stdcxx::first_max_by`.
- `std::rotate`, `std::unique`, `std::erase` map directly.

## 7. Linear algebra

- All factorizations go through `solver::sparse_solver` (`SparseLdlt`, `SparseLlt`, `SparseLu`): faer 0.24.4 behind
  an Eigen-like interface (`compute` / `analyze_pattern` / `factorize` / `info()` / `solve`). LU reports a numerically
  zero pivot like Eigen's `SparseLU` does.
- `FullPivLU`, `SparseQR` (rank revealing) and the `COLAMDOrdering` it uses are ports of Eigen 5.0.1 in
  `solver/eigen/` (faer has no rank revealing variants, and faer's COLAMD orders differently: 96/300 random
  matrices agree).
- `nonzeros(&S)` walks a matrix column by column, rows ascending, as the C++ helper does on compressed matrices.

## 8. libigl and other dependencies

Ports of the libigl v2.6.4 functions metriko calls live in `rust/metriko/src/igl/<name>.rs`, one per libigl file,
same signature shape (outputs returned instead of written to references), same output order (edge numbering,
adjacency order, tie breaks). `Eigen::umeyama` lives in `igl/umeyama.rs`.

## 9. Parallelism

`#pragma omp parallel for` loops stay sequential by default. Under the `rayon` feature they may use `par_iter` only
when iterations are independent and the result does not depend on the schedule (no reductions in different orders).

## 10. Verification

`tools/dump` (C++) writes every stage of `compute_quadrangulation_impl` to `testdata/golden/<commit>/<mesh>@<scale>/`
(see `tools/dump/FORMAT.md`). `rust/parity` compares the Rust stages against those files: discrete data must be equal,
continuous data within a tolerance stated per test.
