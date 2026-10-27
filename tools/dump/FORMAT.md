# Golden stage dumps of the C++ metriko

`metriko_dump <mesh> <scale> <outdir>` loads the mesh as metriko_web's native cli does (`quadpatch.h` `load_mesh`),
runs `include/metriko/lib.h` `compute_quadrangulation_impl(V, F, scale)` (defaults: `calign = false`,
`refine = true`) **step by step**, calling the same functions in the same order, and writes the data of every step.
It then runs `compute_quadrangulation_impl` itself and exits with status 3 if its quads differ bitwise from the
step-by-step ones: `dump.cpp` must be updated whenever `lib.h` changes.

`tools/golden.sh [filter]` builds the tool against the library at the commit of `rust/PORTED_FROM` (temporary
`git worktree`) and writes every case of `testdata/cases.json` to `testdata/golden/<7-char commit>/<stem>@<scale>/`.

## Build configuration (the reference)

The same as metriko_web's native cli, which reproduces the wasm module: Eigen 5.0.1 solvers (`SUITESPARSE=OFF`),
`ENABLE_OPENMP=OFF`, `IGL_PARALLEL_FOR_FORCE_SERIAL`, `-ffp-contract=off`, libigl v2.6.4, Release `-O2 -DNDEBUG`
(metriko_web uses `-O3`; the outputs are identical, checked on spot / fandisk / fertility), Apple clang 17 / libc++ 20.1.

## Files of a case directory

| file | content |
|---|---|
| `<stage>.<name>.npy` | one array, numpy format 1.0, little endian, C order (row-major) |
| `manifest.json` | `files`: every npy with dtype and shape; `meta`: scalars per stage (`"<stage>.<key>"`) |
| `summary.json` | the metriko_web cli summary: counts, `pos_sum`, `hash`, the four printed `lines`, timings |
| `dump.log` | stdout / stderr of the run (library prints included) |

dtypes: `<f8` double, `<i4` int, `<c16` `std::complex<double>` (re, im), `|u1` bool (0 / 1).

Conventions:
- Ids are the C++ `int`s, `-1` meaning none / dead, exactly as in memory (dead elements are dumped too).
- `Eigen::Matrix` with `r` rows and `c` columns → shape `(r, c)`; `VecX*` and `std::vector` → shape `(n,)`.
- Ragged lists `vec<vec<T>>` → `<name>_offsets` (`(n + 1,)`, CSR style: item `i` is rows
  `offsets[i] .. offsets[i + 1]`) and the concatenated rows.
- `HmLoc` (variant) → `<name>_ki` `(n, 2)` i4 = (kind, id) with kind = variant index **0 V, 1 E, 2 F, 3 H, 4 C, 5 P**,
  and `<name>_val` `(n, 2)` f8 = (E, H: `r`, 0), (F: `xy.re`, `xy.im`), (P: `uv.re`, `uv.im`), (V, C: 0, 0).
- Sparse matrices `SprsD` → `<name>_rc` `(nnz, 2)` i4 = (row, col) and `<name>_val` `(nnz,)` f8, in compressed
  column-major storage order (column by column, rows ascending); the shape is meta `<stage>.<name>_shape`.
- `JunctionType`: 0 None, 1 F, 2 T, 3 C.

### An `Hmesh` (`dump_hmesh`, prefix = stage)

`pos (nV, 3) f8`, `idx (nF, 3) i4`, `edge2vert (nE, 2)`, `face2edge (nF, 3)`, `edge2face (nE, 2)`, `vert2half (nV)`,
`edge2half (nE)`, `face2half (nF)`, `loop2half (nL)`, `crnr2half (nC)`, `next`, `prev`, `twin`, `tail`, `head`,
`edge`, `face`, `crnr` (all `(nH,)` i4), `isBV (nV) u1`; meta `counts` {nV, nF, nE, nH, nL, nC}. With geometry
(stages 01, 19): `vertNormal`, `vertBasisX`, `vertBasisY` `(nV, 3)`, `faceNormal`, `faceBasisX`, `faceBasisY`,
`baryCenter` `(nF, 3)`, `faceArea (nF)`, `baryDualArea`, `circDualArea`, `angleDefect` `(nV)`, `halfCotan`,
`dihedralArg`, `heArgOnVert`, `heArgOnFace` `(nH)`, `edgeCotan (nE)` (all f8). `cfn` is dumped separately.

### An `Emesh` (`dump_emesh`, prefix = stage)

| name | dtype, shape | C++ |
|---|---|---|
| `tnodes_ki`, `tnodes_val` | HmLoc encoding, n = `tnodes.size()` | `tm.tnodes` |
| `fixed` | u1 `(n,)` | `tm.fixed` |
| `tedges_id` | i4 `(nte,)` | `tedges[i].id` (-1: dead) |
| `tedges_nids_offsets`, `tedges_nids` | ragged i4 | `tedges[i].nids` |
| `thalfs_ids` | i4 `(nth, 7)` | (`id`, `twid`, `teid`, `tqid`, `cano`, `bgn`, `end`) |
| `thalfs_xr` | f8 `(nth, 2)` | (`x`, `r`) |
| `tquads_id` | i4 `(ntq,)` | `tquads[i].id` |
| `tquads_data_offsets`, `tquads_data` | ragged, rows i4 `(k, 2)` | `tquads[i].data` as (`thid`, `side`) |

## Stages

| stage | files (`<stage>.<name>`) | C++ (lib.h `compute_quadrangulation_impl`) |
|---|---|---|
| `00_input` | `V (nV,3) f8`, `F (nF,3) i4`, `unit (1,) f8`, `scale (1,) f8`; meta `scale`, `unit` | the loaded mesh; `unit = grid_unit(V, scale)` |
| `01_hm` | Hmesh with geometry | `hm = Hmesh(MatXd(V / unit), F)` |
| `02_rawf` | `field (nF,4) c16`, `connection (nE) c16`, `compressed (nF) c16`, `matching (nE) i4`, `singular (nV) i4` | `rawf = FaceRosyField(*hm, 4, Smoothest)` |
| `03_seam` | `seam (nE) u1` | `compute_seam(rawf)` |
| `04_cutm` | Hmesh, topology only | `cutm = compute_cut_mesh(*hm, seam)` |
| `05_cmbf` | as 02 | `cmbf = compute_combbed_field(rawf, seam)` |
| `06_extf` | `extf (nF,12) f8` | `compute_extrinsic_field(*cmbf, 4)` |
| `07_rp_setup` | sparse `vtrans2cut`, `constraint`, `uncompress`; `he2matching (nH)`, `he2transidx (nH)`, `integerIdcs`, `singularIdcs`, `fixedIdcs` i4; `fixedVals` f8; `lreductor (N,n) i4`; meta `N`, `n` | `rp = RosyParameterization(...)`, `rp.localInjectivity = true`, `rp.setup()`: its public members |
| `08_rp_integ` | `fullx` f8, `nfn (nV_cut, N) f8`, `cfn (nF, 12) f8`; meta `integ_returned` | `rp.integ()` |
| `09_hm_cfn` | `cfn (nC) c16` | `hm->cfn` built from `rp.cfn` columns (0,1), (4,5), (8,9) |
| `10_mg` | `mports_dir (np) c16`, `mports_ids (np,4)` = (`crnr_id`, `this_id`, `next_id`, `prev_id`); `mnodes_jt (nn)`, `mnodes_loc_ki/val`, `mnodes_adj_offsets`, `mnodes_adj (k,2)` = (`curv_id`, `sgmt_id`); `mcurvs_id (nc)`, `mcurvs_sgmt_offsets`, `mcurvs_sgmts (k,5)` = (`fr_nid`, `to_nid`, `this_id`, `curv_id`, `face_id`), `mcurvs_buff_loc_ki/val`, `mcurvs_buff_dir c16`, `mcurvs_buff_end u1` | `mg = Mgrph(*hm, cmbf->matching, cmbf->singular)` |
| `11_tm_init` | Emesh | `tm = Emesh(*mg)` |
| `12_quantization` | `X (nte) f8` | `compute_quantization(*tm, *mg)` |
| `13_tm_setx` | Emesh | after `validate_quantization`, `tm->set_x(X)` |
| `14_tm_collapsed` | Emesh; meta `rounds` (iterations of the loop run) | the collapse loop (`collapse_thalf`, `collapse_tquad_prepare/execute`), `validate_collapse_done` |
| `15_tm_snap0` | Emesh | `collapse_tedge_snap(false)` |
| `16_tm_snap1` | Emesh | `collapse_tedge_snap(true)` |
| `17_tm_dedup` | Emesh | `collapse_tedge_snap_dedup` over `live_tedges()` |
| `18_tm_repaired` | Emesh | `repair_crossing_tedges(*tm)` |
| `19_emb` | Hmesh with geometry; `seam1 (nE) u1`, `matching1 (nE) i4`, `singular1 (nV) i4`, `hdata_ids (k,4)` = (`half.id`, `thid`, `tqid`, `order`), `hdata_vals (k,2)` = (`v0`, `v1`) | `hm_emb = compute_embedding_cut_hmesh(...)` and its outputs |
| `20_hcut` | Hmesh, topology only | `hm_cut = compute_cut_mesh(*hm_emb, seam1)` |
| `21_tutte` | `uv_ebd (nC_emb, 2) f8`, `uv_cut (nV_cut, 2) f8` | `compute_tutte_parameterization(...)`, its gather to `hm_cut` vertices |
| `22_slim` | `uv (nV_cut, 2) f8` | `seamless_slim(..., slim_max_iter, slim_rel_tol)` |
| `23_qex_cfn` | `cfn (nC_emb) c16` | `hm_emb->cfn` from `uv` through `hm_cut->idx` |
| `24_qex_sanitized` | `cfn` | after `qex::sanitization(*hm_emb, matching1, singular1, 4)` |
| `25_qverts` | `{v,e,f}_uv c16`, `{v,e,f}_pos (k,3) f8`, `{v,e,f}_sid i4` | `generate_q_vert` → `vqvs`, `eqvs`, `fqvs` |
| `26_qports` | `qports_ids (k,7)` = (`idx`, `vid`, `eid`, `fid`, `isConnected`, `next_id`, `prev_id`), `qports_uv`, `qports_dir` c16, `qports_pos (k,3)` | after `generate_{v,e,f}qvert_qport` |
| `27_qedges` | `qedges (ke,2)` = positions in `q_ports` of `port1`, `port2`; `qports_*` again | `generate_q_edge` (it connects the ports) |
| `28_qfaces` | `qfaces_offsets`, `qfaces_qhalfs (k,3)` = (position of `qe` in `q_edges`, `idx`, `cano`) | `generate_q_faces` |
| `29_quads` | `q_pos_unscaled (nq_v,3) f8`, `q_pos` (times `unit`), `q_idx (nq,4) i4`, `q_val (nq) i4` | `extract_quad_mesh(*hm, q_faces, refine)`, `label_quad_patches(...)`, `q_pos *= unit` |

## What is not dumped

Only data reachable from `lib.h` is dumped; values internal to a function are not, since dumping them would mean
modifying the library: inside `rp.setup()` / `rp.integ()` (the iterative rounding: the constraint matrix after
`reduce_to_linearly_independent`, the initial solve and which branch of `iter_rounding_init` (LLT or the LU of the
KKT system) was taken, each rounding step), inside `compute_quantization` (strip graph, negative cycles), inside
`FaceRosyField` (the connection Laplacian, the inverse power iteration), the per-iteration state of
`seamless_slim`, the private members of `RosyParameterization`, and the internal state of the collapse / snap /
reroute operations between the dumped stages.

## Writing a stage test (rust/parity/tests/sNN_*.rs)

1. Iterate `parity::golden_cases()` (it panics when no dump exists, so a test cannot pass vacuously).
2. Rebuild the stage's **inputs from the dump of the previous stage**, not from the Rust output of the previous
   stage, so that each stage is checked in isolation (e.g. the Mgrph test builds `Hmesh` from `01_hm.pos/idx`,
   sets its `cfn` from `09_hm_cfn.cfn`, and takes `matching` / `singular` from `05_cmbf`). Where an input is a
   structure the dump cannot rebuild (an `Emesh` mid-collapse), chain from the last stage that matched and say so.
3. Run the Rust stage and compare with `parity::compare`:
   - discrete data (ids, topology, offsets, bools, matching, singular, junction types, quad indices, patch ids) with
     `assert_eq_i32` / `assert_eq_bool`: they must be equal;
   - floating data bitwise (`assert_bits_*`) where the computation is a straight port of plain arithmetic (mesh
     geometry, `unit`), otherwise with `assert_close_*` and a tolerance stated and justified in the test (solver
     outputs differ in the last bits between Eigen and faer).
4. Report per case (`c.name()` in every message).
