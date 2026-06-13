#pragma once
// Factor-prediction table for the Fubini reduction.
// Spec: docs/superpowers/specs/2026-06-11-stfactorpredictor-design.md
//
// Single-chain replay along a fixed LR order (deliberately NOT the DP
// intersection semantics of find_lr_orders: the chain set bounds the
// letters of intermediate expressions along this specific order and is
// a superset of the DP's intersected set; entries for letters the
// integrator never materializes are wasted space, never errors).
//
// Tabulated objects, per stage (pivot = order[k]):
//   - pair entries, admissible deg-1 letters F, G (canonical reps):
//         F/lc(F) - G/lc(G) = c * Prod intern[id]^exp
//     with numerator lc(G)F - lc(F)G == -Res_pivot(F, G) (identity),
//     denominator lc(F) lc(G); exponents signed.
//   - singleton entries, every admissible letter: factor lists of each
//     coefficient w.r.t. the pivot, plus the discriminant for deg-2.
// Trial division against the stage pool (st_fubini_lr output + current
// letters); fmpq_mpoly_factor fallback on non-constant remainder sets
// the `oop` flag.  Constants are never interned; they fold into c.

#include "hyperflint/core/poly.hpp"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hyperflint {
namespace factor_table {

// c * Prod intern[id]^exp, c a rational constant (decimal string,
// fmpq_get_str format; "0" for a vanishing object).
struct FactoredObject {
    std::string c = "1";
    std::vector<std::pair<size_t, long>> factors;  // (intern id, exponent)
    bool oop = false;   // fallback factorization fired
};

struct PairEntry {
    size_t var_idx;     // wide-ctx pivot index
    size_t f_id, g_id;  // intern ids; intern_strs[f_id] < intern_strs[g_id]
    FactoredObject diff;
};

struct SingletonCoeff {
    long power;
    FactoredObject fo;
};

struct SingletonEntry {
    size_t var_idx;
    size_t id;          // intern id of the letter (canonical rep)
    long deg;
    std::vector<SingletonCoeff> coeffs;  // powers deg..0
    bool has_disc = false;
    FactoredObject disc;
};

struct StageInfo {
    size_t var_idx;
    std::vector<std::vector<size_t>> admissible;  // per group, intern ids
    std::vector<size_t> pool;                     // intern ids (union over groups)
    size_t n_pairs = 0, n_singletons = 0, n_inadmissible = 0;
    double t_build_s = 0.0;
};

// Loud guards (errors, never truncation).  max_response_mb is enforced
// at serialization time by the handler.
struct Limits {
    size_t max_pairs = 2000000;
    size_t max_singletons = 200000;
    size_t max_response_mb = 512;
};

struct Stats {
    size_t pairs_total = 0, singletons_total = 0, oop = 0;
    size_t pair_fallbacks = 0;  // must stay 0 (deg-1 pool theorem)
    double trial_s = 0.0, fallback_s = 0.0;
};

struct FactorTable {
    std::vector<std::string> intern_strs;  // canonical_prop_form strings
    std::vector<Poly> intern_polys;        // parallel to intern_strs
    std::unordered_map<std::string, size_t> intern_idx;
    std::vector<StageInfo> stages;
    std::vector<PairEntry> pairs;
    std::vector<SingletonEntry> singletons;
    Stats stats;
};

// Build the table.  `order_indices`: wide-ctx indices, a permutation of
// the integration variables (validated by the handler).  Throws
// std::runtime_error on guard violation, naming the guard and the
// count reached.  The caller is responsible for reset_lr_memos() /
// reset_lr_trace() (per-request memo hygiene).
FactorTable build(const PolyCtx& ctx,
                  const std::vector<std::vector<Poly>>& group_polys,
                  const std::vector<size_t>& order_indices,
                  bool algebraic_letters,
                  const Limits& limits);

}  // namespace factor_table
}  // namespace hyperflint
