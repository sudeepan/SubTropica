// Doppio-port phase 3: projective Cheng-Wu GAUGE SCAN + the Doppio
// keep-rules (Euler-conic strict tier, FindRoots carried-sqrt tier) as a
// VERDICT ENGINE on top of the lr_search table machinery.
//
// This is deliberately a SEPARATE entry point from find_lr_orders: the
// classic scorer stays byte-identical (its quadratic acceptance mirrors
// what HF's integrator can actually execute), while the scan implements
// the Mathematica DoppioFubini semantics (doppio_lib.wl) for
// triage/certification:
//
//   * PROJECTIVE input detection: every factor homogeneous in the
//     integration variables and  sum_i e_i deg P_i == -n  identically in
//     eps, with exponents e_i = a_i + b_i*eps supplied as integer pairs
//     (so the check is  sum a_i d_i == -n  AND  sum b_i d_i == 0).
//   * GAUGE SCAN ("ChengWu" -> "Scan"): crawl the UNGAUGED n-variable
//     system to depth n-1; the un-integrated variable IS the Cheng-Wu
//     gauge.  One subset table serves all n gauges; per gauge the
//     letters are DEHOMOGENIZED at x_gauge -> 1 before judgment
//     (dehomogenization commutes with the disc/resultant letter
//     generation and preserves irreducibility for gauge-coprime
//     letters, so this equals the explicit gauged runs -- validated
//     against the Mathematica t24 battery).
//   * KeepRule::Strict  == Doppio's strict rule: a letter may have
//     degree <= 1 in the pivot, or be a quadratic whose Euler conic
//     rationalizes (leading coefficient A and constant term C both
//     perfect squares as polynomials; FLINT fmpq_mpoly_sqrt ==
//     dpSquareQ semantics, rejecting 2x^2 and -x^2).
//   * KeepRule::FindRoots adds the relaxed tiers: TERMINAL quadratics
//     accepted (their roots are kinematic letters of the answer);
//     otherwise the ODD-MULTIPLICITY part of the pivot-discriminant is
//     examined -- empty => rational roots (zero-disc double lines
//     accepted too); pure-kinematic => kinematic sqrt letter; pending-
//     variable-dependent => the letter is accepted and the obligation
//     CARRIED: it joins the letter set of every later step of the same
//     path and is re-judged by these same tiers recursively, discharged
//     once free of pending variables.  Orders report
//     {carried_sqrts, kin_sqrts, terminal_quads} profiles.  This tier
//     is NECESSARY-only (same status as production FindRoots->True);
//     prefer the score-minimal order, treat deep-carry orders as
//     speculative.
//
// The optional chi filter (euler_filter.hpp, HF_EULER_FILTER semantics)
// can be applied to the table letters before judgment, exactly like the
// Mathematica variant C.
#pragma once
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include "hyperflint/core/poly.hpp"

namespace hyperflint {
namespace lr_scan {

enum class KeepRule { Strict, FindRoots };

// Sentinel for step_fr_judge's `gauge` argument meaning "no Cheng-Wu
// gauge" — letters are judged AS GIVEN, with no dehomogenization.  Used
// by find_lr_orders' carry-discharge path (the production per-gauge
// integrand is already gauge-fixed upstream, so there is no further
// gauge to fold), and distinguishable from any real ctx var index
// (which is always < the number of variables, far below SIZE_MAX).
inline constexpr std::size_t kNoGauge = static_cast<std::size_t>(-1);

// twist exponent  a + b*eps  (integer pairs cover the Symanzik cases:
// U^(omega - D/2) F^(-omega) at D = 4 - 2 eps has integer a, b)
struct ScanExponent {
    long a = 0;
    long b = 0;
};

struct ScanOrder {
    std::vector<size_t> order;   // the n-1 integrated vars (ctx indices)
    size_t gauge = 0;            // the un-integrated variable (ctx index)
    double score = 0.0;
    unsigned long carried_sqrts = 0;
    unsigned long kin_sqrts = 0;
    unsigned long terminal_quads = 0;
};

struct ScanResult {
    bool projective = false;     // detection verdict (scan ran iff true)
    bool truncated = false;      // order enumeration hit max_orders
    std::vector<ScanOrder> orders;   // score-ascending
};

// ---------------------------------------------------------------------
// projectivity detection (dpProjectiveQ): every factor of every group
// homogeneous in the xvar set, and sum e_i deg P_i == -n identically.
// exps shape matches group_polys (RAW groups, no boundary monomials).
// ---------------------------------------------------------------------
bool projective_input(
    const std::vector<std::vector<Poly>>& group_polys,
    const std::vector<size_t>& xvar_indices,
    const std::vector<std::vector<ScanExponent>>& exps);

// ---------------------------------------------------------------------
// the per-letter FindRoots judgment (exposed for unit tests; the port
// of dpFRJudge).  `pivot` and `pending` are ctx var indices; `pending`
// excludes the pivot and the gauge.  Returns ok + the new obligations
// (canonicalized) + counter increments.
// ---------------------------------------------------------------------
struct FrJudgment {
    bool ok = false;
    std::vector<Poly> carry;
    unsigned long kin = 0;
    unsigned long term = 0;
};
FrJudgment fr_judge(const Poly& letter, size_t pivot,
                    const std::vector<size_t>& pending,
                    const std::vector<size_t>& all_xvars);

// Euler-conic rationalizability of a deg-2-in-pivot letter (exposed for
// unit tests): leading and constant coefficients both perfect squares.
bool conic_rationalizable(const Poly& letter, size_t pivot);

// ---------------------------------------------------------------------
// Carried-obligation path state for the FindRoots tier.  `carried` holds
// the canonical (proportionality-representative) forms of the deferred
// sqrt-obligations; `carried_keys` is the dedup set; the counters
// accumulate {carried_sqrts, kin_sqrts, terminal_quads} along the path.
// One PathState per DFS branch — copied on descent so siblings do not
// share carry.
// ---------------------------------------------------------------------
struct PathState {
    std::vector<Poly>          carried;       // canonical forms
    std::set<std::string>      carried_keys;  // dedup
    unsigned long nsq = 0, nkin = 0, ntq = 0; // carried / kin / terminal
};

// ---------------------------------------------------------------------
// One carry-discharge step (the extracted core of the FindRoots tier;
// the port of the Mathematica stepScan).  Shared verbatim by BOTH
// engines so the carry semantics are identical by construction:
//   * lr_scan's per-gauge DFS (gauge = the un-integrated Cheng-Wu var),
//   * find_lr_orders' carry-discharge path (gauge = kNoGauge).
//
// Semantics (mirrors fr_judge exactly; no new carry rules):
//   1. discharge carried obligations now free of `pending` AND of the
//      pivot — i.e. their sqrt-argument no longer depends on any
//      not-yet-integrated variable;
//   2. join the (dehomogenized, when gauge != kNoGauge) `letters` with
//      the still-pending carried obligations, dedup by canonical key,
//      table-first;
//   3. fr_judge each survivor at `pivot` with `pending`; on a single
//      `!ok` the WHOLE step fails (return false);
//   4. fold the new obligations into `st.carried`, bumping nsq, and
//      accumulate kin / terminal counters.
// `pending` is the set of integration-variable ctx indices that remain
// un-integrated AFTER this step (it excludes the pivot, and — in the
// gauge scan — the gauge).  Returns true iff the step is admissible;
// `st` is updated in place.
bool step_fr_judge(const std::vector<Poly>& letters, size_t pivot,
                   size_t gauge, const std::vector<size_t>& pending,
                   PathState& st);

// ---------------------------------------------------------------------
// the gauge scan.  group_polys are the RAW per-face groups (boundary
// monomials are added internally, mirroring dpLungoCore's seed).
// max_orders caps the enumeration (truncated flag set on hit).
// ---------------------------------------------------------------------
ScanResult find_lr_orders_scan(
    const std::vector<std::vector<Poly>>& group_polys,
    const std::vector<size_t>& xvar_indices,
    const std::vector<std::vector<ScanExponent>>& exps,
    KeepRule keep_rule = KeepRule::Strict,
    bool euler_filter = false,
    size_t max_orders = 8192);

}  // namespace lr_scan
}  // namespace hyperflint
