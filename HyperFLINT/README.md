# HyperFLINT

A C++17 reimplementation of **HyperIntica** (Mathematica) and
**HyperInt** (Maple), the Feynman-integral-friendly hyperlogarithm
integration packages of Panzer (arXiv:1403.3385), built on FLINT's
`fmpq_mpoly` and Antic polynomial / number-field infrastructure.

**Status: Phase 7 complete + parallel integrator + interval
rescaling + multi-variable FibrationBasis + opt-in divergence
detection.**

For performance, the correctness scope, and the algorithm, see the SubTropica paper (arXiv:2604.20954). HyperFLINT cross-validates against HyperIntica and, where available, HyperInt to 30+ digits.

### Divergence detection (opt-in)

HyperFLINT defaults to assuming the input integrand is convergent
(matching the fast HyperIntica path). Users who want boundary-
divergence scanning can set `"check_divergences": true` in the
`integration_step` JSON request. When enabled, HyperFLINT runs a
secondary accumulation pass that visits every `{logpower, power}`
pole bin and tests each non-(0, 0) bin via `test_zero_function_sym`
(built on the multi-variable `fibration_basis_sym` port). Any
non-zero residue raises `HyperFLINTDivergentIntegral` with the
specific `Log[var]^k / var^n` combination flagged.

Cost: the opt-in check roughly doubles the integration_step cost
for that variable. Convergent inputs pay nothing if the flag is
left off (the default).

**Check-pass design notes** (latest: commit `bbad9bac3`):
- Both zero- and infinity-boundary divergences are scanned
  symmetrically (caveat-1 closed `7a8b2d46e`). The exception's
  `where` field disambiguates `"var=z (zero boundary)"` vs
  `"var=z (infinity boundary)"`.
- Each bin's accumulated RegulatorSym is passed through Fragment-P2
  positive-letter closure (shared helper
  `close_positive_letters_in_regulator_sym`) before the zero test
  (caveat-2 closed `bbad9bac3`). This matches HyperIntica's
  pre-TestZeroFunction redistribution ordering
  (HyperIntica.wl:4832-4856) and prevents false positives on
  Smirnov-tst1-class integrands whose prior-step P2 closure left
  positive-integer letters in downstream RegKeys.
- When the driver passes a populated `remaining_var_indices`, the
  bin residue is tested via `fibration_basis_sym` over those
  variables, catching shuffle-identity cancellations across keys.
  Standalone `integration_step` calls (no outer schedule) fall back
  to the base-case per-term `is_zero` check and may miss such
  cancellations — so run divergence checks through the
  `hyperflint` / `hyperflint_sym` drivers when possible.

## Goal

Bit-exact parity with HyperIntica and HyperInt on every test case,
with a 3–10× end-to-end speedup on realistic Feynman-integral diagrams
(benchmarked against the five-backend `bench_hyperintica` suite that
motivated this project).

## Phase status

Each phase ships a milestone fixture set that cross-validates against
both reference backends. **Done** means all milestone fixtures land
3-way green; **partial** means some sub-phases land but follow-up
work is documented; **todo** means not started.

| Phase | Topic | Status | Notes |
|---|---|:---:|---|
| 0 | scaffolding (CMake, FLINT link, harness) | ✅ done | `compare.py` 3-way diff harness operational |
| 1 | polynomial foundation (`Poly`, `Rat`) | ✅ done | factor, gcd, derivative, eval, partial_fractions cross-ok |
| 2 | rational layer (`LinearFactors`, `PartialFractions`, `RatResidue`) | ✅ done | per-pole decomposition matches reference backends |
| 3 | word algebra (Shuffle, Hlog/Mpl, conversions, Diff) | ⚠️ partial | shuffle_product, convert_zero_one, convert_1inf_to_01, diff_hlog/diff_mpl, ConvertToHlogRegInf (Phase 3-a..d — full 7-case Hlog dispatcher, AST-level reg_tail, top-level Plus/Times/Power/Leaf driver, hyperflint expr path; 40 fixtures green Mma↔HF) done. MplAsHlog / PolyLog / full ToMZV deferred |
| 4 | series expansion (HlogSeries, MplSeries, Expand{Zero,Inf}Word) | ⚠️ partial | expand_zero_word, expand_inf_word, hlog_series (zero-limit branch), mpl_series dispatch land; full Taylor branch deferred |
| 5 | integration driver (TransformWord, IntegrateII, IntegrationStep, hyperflint) | ⚠️ partial | end-to-end Rat path lands (`hf_*` fixtures green); hyperflint now accepts `"expr"` JSON field (Phase 3-d) for Mma-style integrands. Phase 5g-i diagnostic (commit b8ce7ac6f) surfaces a multi-step laurent.cpp gap on 3+ variable integrands (Smirnov tst0 fails at step 3 with `rat_var0_coefficient: denominator is not a monomial in var`); Phase 5g-ii is the relaxation. 6d-v-vi-0 already closed `hf_two_vars_coupled` (commit e8744bb61). |
| 6 | MZV + regularized periods | ✅ done | 6a (mzv table), 6b (zero_one_period), 6c–6d-iv (zero_inf_period branches), 6e (test_zero_function stub), 6f (evaluate_periods), 6d-v-i to 6d-v-iv (SymCoef sidecar + recursive BUC + integration_step_sym + hyperflint_sym), 6d-v-vi (close hf_two_vars_coupled), 6d-v-v-i (FibrationBasis skeleton + no-vars base case, 5 fixtures Mma↔HF green) and 6d-v-v-ii (multi-variable FibrationBasis: `fibration_basis_sym` + `test_zero_function_sym`, SymCoef-valued, handles Fragment-P2 residues natively; landed 2026-04-23). |
| 7 | algebraic letters (Wm/Wp via Antic number fields) | ✅ done | 7-i (AlgebraicLetterTable + atoms), 7-ii (linear_factors deg-2 branch), 7-iii (simplify_with_vieta — Mma-faithful per-pair give-up), 7-iv (combine_wm_wp_ratios — Wm/Wp → WmOverWp atom when shrinking), 7-v (back_substitute via sqrt_disc_<i> atoms), 7-vi (algebraic_letters flag plumbed through partial_fractions/integrate_ii/transform_word/integration_step/hyperflint + factored-denom rebuild in partial_fractions + 3 PF fixtures + 3 hyperflint fixtures Mma↔HF green) + 7-vi-b/c (SubTropica wiring: STHyperFlint registers HF's letter table into HyperIntica`$HyperAlgebraicLetterTable`, multi-word regulator keys translate to shuffle products of ZeroInfPeriods, STIntegrateHF alias with auto-downgrade) all done. HF stays in symbolic Wm/Wp atoms — no Antic dependency yet. Maple skips (HyperInt.mpl has no algebraic-letter extension). |
| 8 | LibraryLink bridge + SubTropica integration | ✅ done | Phase γ.1 (find_lr_orders LibraryLink), Phase γ.2 (hyperflint_sym LibraryLink + process-global MZV cache + version guard), Phase γ.2-b (split libomp out of the dylib to avoid dual-OMP collision), Phase γ.2-c (restore OMP via dyld dynamic-lookup against Wolfram's libomp). Long-suite bench: 2.90× aggregate speedup vs HF-CLI, 5–8× on heavy cases. |
| 9 | performance tuning (memoization, threading, FF sampling) | ✅ done | Profile-driven: 13 optimization rounds cut tst0 from 248 s → 0.594 s (420× total). OpenMP parallelism on the outer `integration_step` entry loop gives additional 5–17× on tst1/tst2. See the SubTropica paper (arXiv:2604.20954) for detailed benchmarks. Apple-Silicon-specific follow-ups (QoS promotion, sharded mutexes, Step-4 closure parallelism, GCD opt-in, custom GMP build) are Phase-9b candidates. |
| 10 | release (paclet, docs, v0.1 tag) | ⏳ todo | gated on full divergence-detection (Phase 5e-iii — MVP in place, full coverage still pending) |

**Cross-validation suite**: 282/282 fixtures green as of
Phase 6d-v-v-i close (post-3-d + 5 fibration_basis fixtures =
54 new atop the Phase-7-vi baseline of 228). Mma ↔ HF is the
hard cross-check; Maple skips the algebraic-letter path
(Phase 7), the Mma-style `"expr"` input path (Phase 3), and the
fibration_basis op (Phase 6d-v-v; Maple emission not wired in
harness). The `--all` target is the regression gate before each
commit.

## Build

Requires FLINT ≥ 3.4, `libomp` for OpenMP parallelism, and a C++17
compiler.

```
brew install flint libomp              # macOS
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release \
      -DOpenMP_ROOT=/opt/homebrew/opt/libomp
cmake --build build-release
build-release/hyperflint factor 'x^2 - 1'
# (-1 + x)*(1 + x)
```

For ASAN- or TSAN-instrumented builds:
```
cmake -S . -B build-asan  -DCMAKE_BUILD_TYPE=Debug -DHF_ASAN=ON
cmake -S . -B build-tsan  -DCMAKE_BUILD_TYPE=Debug -DHF_TSAN=ON
```

### Parallelism

The top-level `hyperflint` driver parallelizes the outer entry
loop in `integration_step` via OpenMP. Control with
`OMP_NUM_THREADS`. At one-less-than-core-count the speedup over
serial is typically 5–17× on Log²-class integrands; parallel
efficiency above `cores-1` decays due to OS-scheduler contention.

```
OMP_NUM_THREADS=13 build-release/hyperflint eval-json < tst1.json
```

HF_OPENMP can be disabled at configure time with `-DHF_OPENMP=OFF`
for a strictly-serial build.

## Usage

Two front-ends:

- **Pretty CLI**: `build/hyperflint factor 'x^2 - 1'` etc.
  Pretty output for human inspection.
- **JSON CLI**: `build/hyperflint eval-json` reads a JSON request
  on stdin and writes a JSON response on stdout. This is the
  cross-validation harness's interface; the request schema is
  documented inline at `bridge/cli/main.cpp` per op.

Currently supported ops (in addition to `factor`):
`add`, `sub`, `mul`, `neg`, `pow`, `derivative`, `eval`, `subst`,
`divexact`, `gcd`, `resultant`,
`rat_*` (rational-function arithmetic),
`linear_factors`, `partial_fractions`, `pole_degree`, `rat_residue`,
`shuffle_product`, `shuffle_words`, `concat_mul`, `collect_words`,
`convert_zero_one`, `convert_1inf_to_01`, `convert_ab_to_zero_inf`,
`diff_hlog`, `diff_mpl`,
`expand_zero_word`, `expand_inf_word`, `mpl_sum`, `mpl_series`,
`hlog_zero_expand`, `hlog_series`, `series_expansion`,
`regzero_word`, `reg0`, `reg_head`, `reg_tail`,
`differentiate_wordlist`, `shuffle_symbolic`, `reglim_word`,
`transform_word`, `transform_shuffle`, `integrate_ii`,
`integration_step`, `hyperflint`,
`apply_mzv_reductions`, `zero_one_period`, `zero_inf_period`,
`break_up_contour`, `break_up_contour_sym`,
`evaluate_periods`, `test_zero_function`,
`sym_arith`, `sym_reduce`,
`algebraic_letters_clear`, `algebraic_letters_show`,
`algebraic_letters_allocate`, `simplify_with_vieta`,
`back_substitute`, `combine_wm_wp_ratios`.

## Layout

- `include/hyperflint/` — public headers
  - `core/` — `Poly`, `Rat`, `SymCoef`
  - `algebra/` — shuffle, partial-fractions, conversions, diff,
    algebraic letters (Wm/Wp)
  - `series/` — Hlog/Mpl series, Laurent expansions
  - `integrator/` — TransformWord, primitive, integration_step,
    hyperflint driver
  - `reduce/` — MZV reductions, periods, BreakUpContour
  - `symbols/` — Word, Letter
- `src/` — implementations (mirrors `include/`)
- `bridge/cli/` — command-line interface
- `test/cross/` — three-backend cross-validation harness
  (HyperIntica ↔ HyperInt ↔ HyperFLINT diff)
  - `compare.py` — diff driver
  - `run_mma.wls` / `run_maple.sh` — reference backend runners
  - `fixtures/` — JSON test cases
- `data/` — MZV reduction table (`mzv_reductions.json`)
- `docs/`
  - `function_map.md` — per-function status (1 row per HyperIntica
    public symbol)
  - `phase_6d_v_plan.md` — completed Phase 6d-v sub-phases + §9
    post-mortem
  - `phase_6d_v_vi_plan.md` — open follow-up plan to close the
    `hf_two_vars_coupled` fixture (audit-first; Tier-0 SEGV blocker)

## Cross-validation philosophy

Every primitive ports as a 3-way comparison:

1. The HF C++ implementation lands.
2. A handful of cross-validation fixtures (typically 5–10 per op)
   get added to `test/cross/fixtures/`.
3. `python3 test/cross/compare.py --all` runs each fixture
   through Mma (`run_mma.wls`), Maple (`run_maple.sh`), and HF
   (`build/hyperflint eval-json`), then diffs the canonicalized
   outputs.
4. Maple `skip` is acceptable for ops without a direct 1-to-1
   Maple equivalent (e.g., the SymCoef sidecar from Phase
   6d-v-i). HF↔Mma is the hard guarantee.
5. Adversarial review (`Agent` subagent) stress-tests any plan or
   non-trivial derivation before merge.

A function is **`cross-ok`** in `function_map.md` only after this
loop closes. **`reg-ok`** means it also rides through the
STBenchmark Long regression at the LibraryLink-bridge layer
(Phase 8 milestone).

## License

MIT — see `LICENSE`.
