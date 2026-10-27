# metriko (Rust)

A pure Rust port of the C++ library in `../include/metriko`, file by file and faithful to it (`PORTING.md`).
The crate `metriko` has the same entry point:

```rust
let res = metriko::compute_quadrangulation(&v, &f, 0.01)?; // v: &[DVec3], f: &[[i32; 3]]
res.pos; // quad vertices
res.idx; // quads
res.val; // patch index of each quad
```

## Layout

| path | content |
|---|---|
| `metriko/` | the library; `src/<dir>/<file>.rs` mirrors `include/metriko/<dir>/<file>.h` |
| `metriko/src/igl/` | the libigl v2.6.4 functions metriko calls |
| `metriko/src/solver/eigen/` | Eigen 5.0.1 FullPivLU, SparseQR and COLAMD (rank decisions must match Eigen's) |
| `metriko/src/common/stdcxx.rs`, `stdcomplex.rs` | libc++ 20.1 sort / priority_queue / unordered_map order, and `std::complex` arithmetic as the reference build executes it |
| `parity/` | stage by stage comparison with the C++ (`../tools/dump/FORMAT.md`) |
| `metriko-cli/` | `metriko <mesh.obj> <scale>`: prints the summary and hash of metriko_web's cli |
| `PORTED_FROM` | the C++ commit this tree mirrors |

Linear algebra is faer 0.24.4; the sparse factorizations (LDLT, LLT, LU) differ from Eigen's in the last bits, every
discrete result (seams, singularities, quantization, T-mesh, quads, patch ids) is identical to the C++.

## Following the C++

```sh
git diff $(cat rust/PORTED_FROM) HEAD -- include/   # what is left to port
tools/golden.sh                                     # C++ stage dumps of PORTED_FROM -> testdata/golden/
cd rust && cargo test --release --workspace         # unit tests + parity against the dumps
```

After porting a C++ change, update `PORTED_FROM`, regenerate the dumps and run the parity tests.

## Known difference

On the 28 smoothest cases of `test/test_config.json`, 27 give the same quads, patch ids and cli hash as the C++.
`botijo` at 0.005 gives the same counts (32,691 quads, 226 patches) but a different hash: its `seamless_slim` does not
converge within its 50 iterations, and the energy difference between the faer and the Eigen LLT solves grows from
1e-13 (iteration 0) to 1e-3 (iteration 46), where the convergence counter takes a different value. Every stage
before SLIM matches; the C++ itself would change there with another solver backend (e.g. CHOLMOD).
