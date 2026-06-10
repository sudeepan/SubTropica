// LR-order search.  C++ port of SubTropica's STFasterFubini2.
// See include/hyperflint/integrator/lr_search.hpp for scope.

#include "hyperflint/integrator/lr_search.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "hyperflint/algebra/euler_filter.hpp"
#include "hyperflint/integrator/lr_scan.hpp"  // carry-discharge step_fr_judge

namespace hyperflint {
namespace lr_search {

namespace {

// HF_LR_TRACE (docs/env_flags.md, Track-doppio-port section): stderr
// profiling of the find_lr_orders DP.  Level 1 = per-size summaries +
// global wall accumulators (resultant / discriminant / leading-coeff /
// factor / dedup) + the 10 slowest st_fubini_lr calls; level 2 adds one
// line per step call.  Diagnostic only — no effect on results.  The DP
// is single-threaded (handler calls it synchronously), so plain
// namespace-level state is safe.
struct LrTrace {
    int level = 0;          // 0 = off
    double t_lc = 0, t_disc = 0, t_res = 0, t_factor = 0, t_dedup = 0;
    long n_res = 0, n_factor = 0, n_steps = 0;
    size_t max_res_terms = 0;       // largest pre-factor resultant
    // top-10 slowest step calls: (wall_s, subset bits, pivot, in, out)
    struct SlowCall { double wall; uint64_t bits; size_t pivot;
                      size_t n_in, n_out; };
    std::vector<SlowCall> slow;
    // context for the step call currently executing (set by the DP loop)
    uint64_t cur_bits = 0;
    size_t cur_pivot = 0;

    void note_slow(double wall, size_t n_in, size_t n_out) {
        SlowCall c{wall, cur_bits, cur_pivot, n_in, n_out};
        slow.push_back(c);
        std::sort(slow.begin(), slow.end(),
            [](const SlowCall& a, const SlowCall& b) {
                return a.wall > b.wall; });
        if (slow.size() > 10) slow.resize(10);
    }
};
LrTrace g_lr_trace;

double now_s() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Deficiency-3 cure (2026-06-06, notes/hf_lr_search_deficiencies.md):
// the DP's inputs are massively self-similar — Standard-closure faces
// carry O(100) near-identical groups, and the same (poly set, pivot)
// step and the same resultant/discriminant factorizations recur across
// subsets, pivots, and groups.  Two value-preserving memo layers:
//   step memo:   (pivot | input poly strings) -> full st_fubini_lr output
//   factor memo: poly string -> factor bases
// Both are pure caches (results byte-identical to the uncached path),
// reset at every find_lr_orders / find_lr_orders_scan entry so memory
// stays bounded per request.  Opt-out: HF_LR_STEP_MEMO=0 /
// HF_LR_FACTOR_MEMO=0 (docs/env_flags.md).  Single-threaded by the
// same argument as LrTrace.
struct LrMemo {
    bool step_on = true;
    bool factor_on = true;
    std::unordered_map<std::string, std::vector<Poly>> step;
    std::unordered_map<std::string, std::vector<Poly>> factors;
    long step_hit = 0, step_miss = 0, fac_hit = 0, fac_miss = 0;
};
LrMemo g_lr_memo;

bool env_flag_off(const char* name) {
    const char* v = std::getenv(name);
    return v != nullptr && std::strcmp(v, "0") == 0;
}

// Poll each monomial for its exponent in `var_idx` and its coefficient.
// Not optimal but fine for MVP — FLINT's fmpq_mpoly iterators are cheap.
//
// Leaf count per monomial:
//   1 for each variable with exponent >= 1  (the variable symbol)
//   1 for each variable with exponent >= 2  (the exponent integer)
//   1 if the coefficient is not ±1         (the numeric atom)
//
// Summed over monomials of a single poly, then over polys in the list.
long leaf_count_single(const Poly& p) {
    const fmpq_mpoly_struct* raw = p.raw();
    const auto* ctx = p.ctx().raw();
    const slong nt = fmpq_mpoly_length(raw, ctx);
    const size_t nv = p.ctx().vars().size();
    if (nt == 0) return 0;

    std::vector<slong> exps(nv);
    long total = 0;
    fmpq_t coef;
    fmpq_init(coef);
    for (slong t = 0; t < nt; ++t) {
        fmpq_mpoly_get_term_coeff_fmpq(coef, raw, t, ctx);
        fmpq_mpoly_get_term_exp_si(exps.data(),
            const_cast<fmpq_mpoly_struct*>(raw), t, ctx);

        // Coefficient atom.  Mma's LeafCount counts the literal `-1`
        // atom in Times[-1, x, y] (it is an atom in the expression tree),
        // so we count any non-unit coef — whether +k, -1, -k, or 1/k.
        // Rational coefs are undercounted relative to Mma's tree
        // structure (Rational[p, q] is 2 atoms there), but Symanzik
        // inputs are integer-coefficient, so this doesn't bite.
        if (!fmpq_is_one(coef)) total += 1;

        // variable / exponent atoms
        for (size_t v = 0; v < nv; ++v) {
            if (exps[v] >= 1) total += 1;
            if (exps[v] >= 2) total += 1;
        }
    }
    fmpq_clear(coef);
    return total;
}

}  // namespace

long leaf_count_proxy(const std::vector<Poly>& polys) {
    long total = 0;
    for (const auto& p : polys) total += leaf_count_single(p);
    return total;
}

void SingCollector::observe(const Poly& b) {
    const auto* ctx = b.ctx().raw();
    // Numeric / zero factors carry no kinematic divisor.
    if (fmpq_mpoly_is_zero(b.raw(), ctx)) return;
    if (fmpq_mpoly_is_fmpq(b.raw(), ctx)) return;
    // Kinematic-only: the factor must depend on NONE of the integration
    // variables.  A factor that still carries an integration variable is
    // an integrand letter, not a kinematic divisor, and is excluded by
    // design (the order-resolved pipeline wants only the s-plane data).
    for (size_t u : b.used_var_indices()) {
        if (integration_vars.count(u)) return;
    }
    // Canonicalize with the SAME proportionality representative the
    // dedup uses (sign + leading-coefficient normalization over Q), so
    // s - 4*mm and -(s - 4*mm) collapse to one entry, matching the
    // engine's canonical form.
    std::string key = b.canonical_prop_form().to_string();
    if (seen.insert(key).second) ordered.push_back(std::move(key));
}

std::vector<Poly> dedup_proportional(const std::vector<Poly>& polys) {
    std::vector<Poly> out;
    std::unordered_set<std::string> seen;
    for (const auto& p : polys) {
        const auto* ctx = p.ctx().raw();
        if (fmpq_mpoly_is_zero(p.raw(), ctx)) continue;
        if (fmpq_mpoly_is_fmpq(p.raw(), ctx)) continue;
        std::string key = p.canonical_prop_form().to_string();
        if (seen.insert(key).second) out.push_back(p);
    }
    return out;
}

std::vector<Poly> intersect_proportional(
    const std::vector<std::vector<Poly>>& lists) {
    if (lists.empty()) return {};
    if (lists.size() == 1) return dedup_proportional(lists[0]);

    // Canonical-form sets per input list.
    std::vector<std::unordered_set<std::string>> canon_sets;
    canon_sets.reserve(lists.size());
    for (const auto& list : lists) {
        std::unordered_set<std::string> s;
        for (const auto& p : list) {
            const auto* ctx = p.ctx().raw();
            if (fmpq_mpoly_is_zero(p.raw(), ctx)) continue;
            if (fmpq_mpoly_is_fmpq(p.raw(), ctx)) continue;
            s.insert(p.canonical_prop_form().to_string());
        }
        canon_sets.push_back(std::move(s));
    }

    // Intersect all canonical-form sets.
    std::unordered_set<std::string> common = canon_sets[0];
    for (size_t i = 1; i < canon_sets.size(); ++i) {
        std::unordered_set<std::string> next;
        for (const auto& key : common) {
            if (canon_sets[i].count(key)) next.insert(key);
        }
        common = std::move(next);
        if (common.empty()) return {};
    }

    // Representative from the first list.
    std::vector<Poly> out;
    std::unordered_set<std::string> used;
    for (const auto& p : lists[0]) {
        const auto* ctx = p.ctx().raw();
        if (fmpq_mpoly_is_zero(p.raw(), ctx)) continue;
        if (fmpq_mpoly_is_fmpq(p.raw(), ctx)) continue;
        std::string key = p.canonical_prop_form().to_string();
        if (common.count(key) && used.insert(key).second) {
            out.push_back(p);
        }
    }
    return out;
}

std::vector<Poly> st_fubini_lr(const std::vector<Poly>& polys, size_t var_idx,
                               SingCollector* sings) {
    if (polys.empty()) return {};
    const PolyCtx& pctx = polys.front().ctx();
    const auto* ctx = pctx.raw();

    LrTrace& tr = g_lr_trace;
    const bool tron = tr.level > 0;
    const double t_call0 = tron ? now_s() : 0.0;

    // Step memo lookup (deficiency-3 cure; see LrMemo).
    std::string step_key;
    if (g_lr_memo.step_on) {
        step_key.reserve(polys.size() * 32 + 16);
        step_key += std::to_string(var_idx);
        for (const auto& p : polys) {
            step_key += '|';
            step_key += p.to_string();
        }
        auto it = g_lr_memo.step.find(step_key);
        if (it != g_lr_memo.step.end()) {
            ++g_lr_memo.step_hit;
            // On a step-memo hit the factor loop is skipped, so observe
            // the cached step output here.  `out` already holds the
            // irreducible factor bases of this step (post-factor,
            // post-dedup), so this is the same multiset the miss path
            // would feed the collector --- collection is independent of
            // memo state.  Cheap: at most |out| observe() calls, only
            // when collecting (sings != nullptr).
            if (sings != nullptr) {
                for (const auto& b : it->second) sings->observe(b);
            }
            return it->second;
        }
        ++g_lr_memo.step_miss;
    }

    std::vector<Poly> temp;
    temp.reserve(polys.size() * (polys.size() + 1) / 2 + polys.size() * 2);

    // Per-poly: leading coefficient + discriminant (when deg >= 1).
    for (const auto& f : polys) {
        if (fmpq_mpoly_is_zero(f.raw(), ctx)) continue;
        if (fmpq_mpoly_is_fmpq(f.raw(), ctx)) continue;
        const long n = f.degree_in_var(var_idx);
        if (n < 0) continue;
        // Leading coefficient in var_idx (returns a Poly free of var_idx).
        if (n >= 0) {
            const double t0 = tron ? now_s() : 0.0;
            Poly lc = f.coefficient_of_var(var_idx, n);
            if (tron) tr.t_lc += now_s() - t0;
            if (!fmpq_mpoly_is_zero(lc.raw(), ctx)) temp.push_back(std::move(lc));
        }
        if (n >= 1) {
            const double t0 = tron ? now_s() : 0.0;
            Poly d = f.discriminant_in_var(var_idx);
            if (tron) tr.t_disc += now_s() - t0;
            if (!fmpq_mpoly_is_zero(d.raw(), ctx)) temp.push_back(std::move(d));
        }
    }

    // Pairwise: resultant when both polys have degree >= 1 in var_idx.
    const size_t N = polys.size();
    for (size_t i = 0; i < N; ++i) {
        const Poly& fi = polys[i];
        if (fmpq_mpoly_is_zero(fi.raw(), ctx)) continue;
        if (fi.degree_in_var(var_idx) < 1) continue;
        for (size_t j = i + 1; j < N; ++j) {
            const Poly& fj = polys[j];
            if (fmpq_mpoly_is_zero(fj.raw(), ctx)) continue;
            if (fj.degree_in_var(var_idx) < 1) continue;
            const double t0 = tron ? now_s() : 0.0;
            Poly r = fi.resultant(fj, var_idx);
            if (tron) {
                tr.t_res += now_s() - t0;
                ++tr.n_res;
                tr.max_res_terms = std::max(tr.max_res_terms, r.n_terms());
            }
            if (!fmpq_mpoly_is_zero(r.raw(), ctx)) temp.push_back(std::move(r));
        }
    }

    // Factor each temp entry, flatten factor bases.  The factor memo
    // collapses recurring resultants/discriminants (the same algebraic
    // objects recur across subsets, pivots, and near-identical groups).
    std::vector<Poly> factored;
    factored.reserve(temp.size() * 2);
    for (const auto& p : temp) {
        if (fmpq_mpoly_is_zero(p.raw(), ctx)) continue;
        if (fmpq_mpoly_is_fmpq(p.raw(), ctx)) continue;
        if (g_lr_memo.factor_on) {
            std::string fkey = p.to_string();
            auto it = g_lr_memo.factors.find(fkey);
            if (it != g_lr_memo.factors.end()) {
                ++g_lr_memo.fac_hit;
                for (const auto& b : it->second) factored.push_back(b);
                continue;
            }
            ++g_lr_memo.fac_miss;
            const double t0 = tron ? now_s() : 0.0;
            auto bases = factor_bases(p);
            if (tron) { tr.t_factor += now_s() - t0; ++tr.n_factor; }
            for (const auto& b : bases) factored.push_back(b);
            g_lr_memo.factors.emplace(std::move(fkey), std::move(bases));
        } else {
            const double t0 = tron ? now_s() : 0.0;
            auto bases = factor_bases(p);
            if (tron) { tr.t_factor += now_s() - t0; ++tr.n_factor; }
            for (auto& b : bases) factored.push_back(std::move(b));
        }
    }

    const double t_d0 = tron ? now_s() : 0.0;
    std::vector<Poly> out = dedup_proportional(factored);
    // Per-face kinematic-divisor collection (order-resolved sings).
    // Observe the step's deduped factor-base set: every irreducible
    // factor that is free of all integration variables is recorded.
    // `out` is the SAME object the step memo stores and the memo-hit
    // path observes, so collection is identical whether or not the step
    // is served from cache.  observe() itself filters non-kinematic /
    // numeric factors, so this loop is a no-op for integrand letters.
    if (sings != nullptr) {
        for (const auto& b : out) sings->observe(b);
    }
    if (tron) {
        tr.t_dedup += now_s() - t_d0;
        ++tr.n_steps;
        const double wall = now_s() - t_call0;
        if (tr.slow.size() < 10 || wall > tr.slow.back().wall)
            tr.note_slow(wall, polys.size(), out.size());
        if (tr.level >= 2) {
            std::fprintf(stderr,
                "[lrtrace2] bits=%llx pivot=%zu in=%zu out=%zu wall=%.3fs\n",
                (unsigned long long) tr.cur_bits, tr.cur_pivot,
                polys.size(), out.size(), wall);
        }
    }
    if (g_lr_memo.step_on) g_lr_memo.step.emplace(std::move(step_key), out);
    return out;
}

void reset_lr_memos() {
    g_lr_memo = LrMemo{};
    g_lr_memo.step_on   = !env_flag_off("HF_LR_STEP_MEMO");
    g_lr_memo.factor_on = !env_flag_off("HF_LR_FACTOR_MEMO");
}

namespace {

// Enumerate all bitmask subsets of {0,...,n-1} of a given size.
// Uses the "next higher bitmask with same popcount" trick (Gosper's hack).
std::vector<uint64_t> subsets_of_size(size_t n, size_t k) {
    std::vector<uint64_t> out;
    if (k > n) return out;
    if (k == 0) { out.push_back(0); return out; }
    uint64_t mask = (1ull << k) - 1;
    const uint64_t limit = 1ull << n;
    while (mask < limit) {
        out.push_back(mask);
        // Gosper's hack: next bitmask with same popcount
        uint64_t c = mask & -mask;
        uint64_t r = mask + c;
        mask = (((r ^ mask) >> 2) / c) | r;
    }
    return out;
}

// ---------------------------------------------------------------------
// Carry-discharge DFS over find_lr_orders' (un-pruned) set_table.
//
// This is the gauge-FREE analog of lr_scan's per-gauge DFS: it integrates
// ALL n variables (no Cheng-Wu gauge — the production per-gauge integrand
// is already gauge-fixed upstream, so find_lr_orders sees exactly one
// gauge's system) and judges each step with the SHARED carry-discharge
// primitive lr_scan::step_fr_judge (gauge == kNoGauge).  set_table[bits]
// is the path-INDEPENDENT S-marginal letter set per group; the PathState
// is the path-DEPENDENT carried-obligation state.  That split is exactly
// why the carry verdict cannot ride the best-score subset-DP memo: two
// orders reaching the same `bits` may carry different obligations.
//
// We keep only the SCORE-MINIMAL admissible full-order (matching the DP's
// MinimalBy and lr_scan's score sort), with its carried-sqrt profile and
// the deg-2 letters encountered along it (for root_polys).  Score is the
// same accumulator as the DP / scan: sum over groups of
// leaf_count_proxy(parent-subset letters)^1.15 at each step.
struct CarryDfs {
    // set_table[bits] -> per-group letters (G lists), un-pruned.
    const std::unordered_map<uint64_t,
        std::vector<std::vector<Poly>>>& set_table;
    const std::vector<size_t>& xvar_indices;  // bit -> ctx index
    size_t n = 0;
    double INF = std::numeric_limits<double>::infinity();

    // best-so-far
    bool found = false;
    double best_score = std::numeric_limits<double>::infinity();
    std::vector<size_t> best_order;     // ctx indices
    std::vector<Poly>   best_roots;     // deg-2 letters along best path
    unsigned long best_nsq = 0, best_nkin = 0, best_ntq = 0;

    // current path's deg-2-letter accumulator (parallel to `order`)
    void dfs(uint64_t bits, std::vector<size_t>& order, double score,
             std::vector<Poly>& roots, const lr_scan::PathState& st) {
        if (order.size() == n) {
            if (!found || score < best_score) {
                found = true;
                best_score = score;
                best_order = order;
                best_roots = roots;
                best_nsq = st.nsq;
                best_nkin = st.nkin;
                best_ntq = st.ntq;
            }
            return;
        }
        // Branch-and-bound: the score is monotone non-decreasing along a
        // path (leaf_count_proxy >= 0, pow >= 0), so a partial path that
        // already exceeds the best complete score cannot improve on it.
        if (found && score >= best_score) return;

        auto it_parent = set_table.find(bits);
        if (it_parent == set_table.end()) return;  // shouldn't happen
        const std::vector<std::vector<Poly>>& parent = it_parent->second;
        const size_t G = parent.size();

        for (size_t bit = 0; bit < n; ++bit) {
            if (bits & (1ull << bit)) continue;
            const size_t var_idx = xvar_indices[bit];

            // pending AFTER this step: unset bits except the pivot (no
            // gauge to exclude — kNoGauge).
            std::vector<size_t> pending;
            pending.reserve(n);
            for (size_t b = 0; b < n; ++b)
                if (b != bit && !(bits & (1ull << b)))
                    pending.push_back(xvar_indices[b]);

            // letters of the PARENT subset, all groups flattened (same
            // source lr_scan's DFS uses), and the score extension (same
            // accumulator as the DP / scan: per-group leaf_count^1.15).
            std::vector<Poly> letters;
            double ext = score;
            std::vector<Poly> step_roots;  // deg-2-in-pivot letters here
            for (size_t g = 0; g < G; ++g) {
                const auto& gl = parent[g];
                letters.insert(letters.end(), gl.begin(), gl.end());
                ext += std::pow(
                    static_cast<double>(leaf_count_proxy(gl)), 1.15);
                for (const auto& p : gl)
                    if (p.degree_in_var(var_idx) == 2)
                        step_roots.push_back(p);
            }

            lr_scan::PathState next = st;
            if (!lr_scan::step_fr_judge(letters, var_idx,
                    lr_scan::kNoGauge, pending, next))
                continue;

            order.push_back(var_idx);
            const size_t roots_mark = roots.size();
            roots.insert(roots.end(),
                std::make_move_iterator(step_roots.begin()),
                std::make_move_iterator(step_roots.end()));
            dfs(bits | (1ull << bit), order, ext, roots, next);
            // Pop the step's roots back off (Poly is not default-
            // constructible, so erase rather than resize).
            roots.erase(roots.begin() + static_cast<std::ptrdiff_t>(roots_mark),
                        roots.end());
            order.pop_back();
        }
    }
};

}  // namespace

LrResult find_lr_orders(
    const std::vector<std::vector<Poly>>& group_polys,
    const std::vector<size_t>& xvar_indices,
    bool allow_algebraic_letters,
    SingCollector* sings,
    bool carry_discharge) {
    // Carry-discharge (Doppio FindRoots) tier is only meaningful when
    // deg-2 letters are admitted: with deg<=1 letters there is no
    // sqrt-obligation to carry, so the flag is a no-op there and the
    // classic subset-DP runs unchanged.  When active, the full
    // (un-pruned) set_table must survive the size loop so the per-path
    // DFS can read every parent subset, and the Strict early-NOLR exit
    // must be suppressed (carry can rescue a Strict-NOLR face).
    const bool do_carry = allow_algebraic_letters && carry_discharge;
    // Seed the kinematic-divisor collector (if requested) with the
    // integration-variable index set, so st_fubini_lr can decide which
    // factors are kinematic-only.  Empty/degenerate early returns below
    // leave the collector empty, which is correct (no walk => no
    // divisors).
    if (sings != nullptr) {
        sings->integration_vars.clear();
        for (size_t v : xvar_indices) sings->integration_vars.insert(v);
    }
    if (group_polys.empty()) {
        return LrResult{{}, 0.0, {}};
    }
    const size_t G = group_polys.size();
    const size_t n = xvar_indices.size();
    if (n == 0) return LrResult{{}, 0.0, {}};
    if (n > 63) {
        throw std::runtime_error(
            "find_lr_orders: > 63 integration variables (bitmask overflow)");
    }

    // Memoized state: bitmask -> per-group poly lists
    // (value.size() == G for every populated key).
    std::unordered_map<uint64_t, std::vector<std::vector<Poly>>> set_table;
    // Memoized order+score per subset.
    std::unordered_map<uint64_t, LrResult> orders_table;

    // Doppio-C Euler chi-drop filter (HF_EULER_FILTER=1, default OFF;
    // docs/env_flags.md): after the Fubini intersection at every subset,
    // each surviving letter is tested against the genuine Euler
    // discriminant of the S-marginal (chi_count_sectors via msolve) and
    // fictitious letters are dropped.  OFF-mode is byte-identical (the
    // branch below never runs).  Boundary monomials are exempt; the
    // verdict is conservative (Indeterminate/failure keep), mirroring
    // dpGenuineDKQ.  The per-call cache memoizes generic chi per
    // (group, subset) marginal.
    const char* euler_env = std::getenv("HF_EULER_FILTER");
    const bool euler_filter_on =
        euler_env != nullptr && *euler_env != '\0'
        && std::strcmp(euler_env, "0") != 0;
    std::vector<ChiFilterCache> chi_caches(euler_filter_on ? G : 0);

    // Per-request memo lifetime (bounded memory; captures all
    // within-face duplication across subsets, pivots, groups).
    reset_lr_memos();
    if (euler_filter_on) reset_chi_filter_stats();

    // HF_LR_TRACE diagnostic profiling (docs/env_flags.md).  Reset the
    // accumulators on every entry so repeated calls in one process
    // (LibraryLink transport) profile independently.
    {
        const char* tr_env = std::getenv("HF_LR_TRACE");
        g_lr_trace = LrTrace{};
        if (tr_env != nullptr && *tr_env != '\0'
            && std::strcmp(tr_env, "0") != 0) {
            g_lr_trace.level = std::atoi(tr_env) >= 2 ? 2 : 1;
            std::fprintf(stderr,
                "[lrtrace] find_lr_orders: G=%zu n=%zu algebraic=%d\n",
                G, n, allow_algebraic_letters ? 1 : 0);
        }
    }

    const double INF = std::numeric_limits<double>::infinity();

    // Seed: set[g][{}] = group_polys[g] for every g, orders[{}] = ({}, 0).
    set_table[0] = group_polys;
    orders_table[0] = LrResult{{}, 0.0, {}};

    // DP over subset size.
    for (size_t size = 1; size <= n; ++size) {
        auto subsets = subsets_of_size(n, size);
        bool any_live_at_size = false;  // any non-NOLR order at this size
        for (uint64_t bits : subsets) {
            // Step A: for each group g, build preSTable (one list per
            // pivot bit v ∈ bits), then intersect.
            std::vector<std::vector<Poly>> set_for_bits(G);
            for (size_t g = 0; g < G; ++g) {
                std::vector<std::vector<Poly>> preTable;
                preTable.reserve(size);
                for (size_t bit = 0; bit < n; ++bit) {
                    if (!(bits & (1ull << bit))) continue;
                    const uint64_t prev_bits = bits ^ (1ull << bit);
                    auto it = set_table.find(prev_bits);
                    if (it == set_table.end()) {
                        // Parent subset was dropped by memory pruning or
                        // never computed — fatal in single-group MVP.
                        throw std::runtime_error(
                            "find_lr_orders: missing parent subset state");
                    }
                    const auto& prev_polys = it->second[g];
                    if (g_lr_trace.level > 0) {
                        g_lr_trace.cur_bits = bits;
                        g_lr_trace.cur_pivot = xvar_indices[bit];
                    }
                    preTable.push_back(
                        st_fubini_lr(prev_polys, xvar_indices[bit], sings));
                }
                set_for_bits[g] = intersect_proportional(preTable);
                if (euler_filter_on && !set_for_bits[g].empty()) {
                    std::vector<size_t> subset_vars;
                    subset_vars.reserve(size);
                    for (size_t b = 0; b < n; ++b)
                        if (bits & (1ull << b))
                            subset_vars.push_back(xvar_indices[b]);
                    set_for_bits[g] = chi_filter_letters(
                        group_polys[g], subset_vars, set_for_bits[g],
                        chi_caches[g]);
                }
            }
            set_table[bits] = std::move(set_for_bits);

            // Step B: compute orders_table[bits] via DP extension.
            LrResult best{{}, INF, {}};

            // HyperIntica-parity guard.  Mirrors the runtime check in
            // linear_factors.cpp:1444-1460 (wide path) / :1212-1229
            // (narrow PIVOT path), which is gated by
            // LFForbiddenVarsScope set at integration_step.cpp:1265.
            // A deg-2 letter base polynomial whose Wm/Wp definition
            // uses any not-yet-integrated Feynman parameter is refused
            // at runtime (pushed to out.nonlinear, then partial_fractions
            // throws "nonlinear factor in denominator").  Modelling the
            // same rejection here keeps the LR-scorer in sync with the
            // integrator: gauges whose runtime guard would crash are
            // scored NOLR rather than picked optimistically.
            //
            // Forbidden set for this DP node is the set of vars that
            // are still un-integrated AFTER taking this step --- i.e.,
            // the unset bits of `bits` (the post-step bitmask), mapped
            // through xvar_indices to wide-context indices.  The set
            // is the same for every candidate `bit` at this `bits`
            // node, so it is computed once.  By construction `var_idx
            // = xvar_indices[bit]` is NEVER in forbidden_after_step
            // (its bit is set in `bits`; xvar_indices is injective by
            // handlers.cpp:367-369), so the runtime's defensive
            // `if (v == var_idx) continue` skip is structurally a
            // no-op here, but kept for line-for-line parity.
            std::vector<size_t> forbidden_after_step;
            if (allow_algebraic_letters) {
                forbidden_after_step.reserve(n);
                for (size_t b = 0; b < n; ++b) {
                    if (!(bits & (1ull << b))) {
                        forbidden_after_step.push_back(xvar_indices[b]);
                    }
                }
            }

            for (size_t bit = 0; bit < n; ++bit) {
                if (!(bits & (1ull << bit))) continue;
                const uint64_t prev_bits = bits ^ (1ull << bit);
                auto it = orders_table.find(prev_bits);
                if (it == orders_table.end()) continue;
                const LrResult& prev = it->second;
                if (prev.nolr()) continue;  // parent is NOLR; skip

                // Precondition: every poly in set[g][prev_bits] has
                // degree ≤ max_deg in x_bit.  max_deg = 1 under classic
                // FindRoots=False; = 2 when algebraic letters are allowed
                // (HF's integrator allocates Wm_i/Wp_i for deg-2 factors
                // at integration time).
                const size_t var_idx = xvar_indices[bit];
                long max_deg = allow_algebraic_letters ? 2L : 1L;
                {
                    const char* env = std::getenv("HF_LR_MAX_DEG");
                    if (env && allow_algebraic_letters) max_deg = std::atol(env);
                }
                const auto& prev_set_all = set_table[prev_bits];
                bool all_linear = true;
                double extension_score = prev.score;
                std::vector<Poly> step_root_polys;
                for (size_t g = 0; g < G; ++g) {
                    const auto& gpolys = prev_set_all[g];
                    for (const auto& p : gpolys) {
                        const long d = p.degree_in_var(var_idx);
                        if (d > max_deg) { all_linear = false; break; }
                        if (d >= 2 && allow_algebraic_letters) {
                            // Forbidden-var check (parity with the
                            // runtime guard).  If `p`'s used vars
                            // intersect forbidden_after_step, the
                            // runtime would push this factor to
                            // out.nonlinear and partial_fractions
                            // would throw -- mark this DP step NOLR
                            // so the LR-scorer doesn't score it.
                            std::vector<size_t> used = p.used_var_indices();
                            bool has_forbidden_dep = false;
                            for (size_t v : forbidden_after_step) {
                                if (v == var_idx) continue;  // see scope comment above
                                for (size_t u : used) {
                                    if (u == v) { has_forbidden_dep = true; break; }
                                }
                                if (has_forbidden_dep) break;
                            }
                            if (has_forbidden_dep) {
                                all_linear = false;
                                break;
                            }
                            // Phase 7-vii: collect the deg-2 polynomial
                            // so the caller knows which polys to turn
                            // into Wm/Wp at integration time.  Mma's
                            // STFasterFubini2 returns the same list as
                            // result[[2]] under FindRoots=True.
                            step_root_polys.push_back(p);
                        }
                    }
                    if (!all_linear) break;
                    long leaf_sum = leaf_count_proxy(gpolys);
                    extension_score +=
                        std::pow(static_cast<double>(leaf_sum), 1.15);
                }
                if (!all_linear) continue;

                if (extension_score < best.score) {
                    std::vector<size_t> new_order = prev.order;
                    new_order.push_back(xvar_indices[bit]);
                    std::vector<Poly> new_roots = prev.root_polys;
                    new_roots.insert(new_roots.end(),
                                     std::make_move_iterator(step_root_polys.begin()),
                                     std::make_move_iterator(step_root_polys.end()));
                    best = LrResult{std::move(new_order),
                                     extension_score,
                                     std::move(new_roots)};
                }
            }
            if (!best.nolr()) any_live_at_size = true;
            orders_table[bits] = best;
        }

        // Early NOLR exit (2026-06-06, value-preserving): NOLR
        // propagates — an order at size k+1 extends a non-NOLR parent
        // of size k, so if EVERY subset at this size is NOLR the final
        // verdict is NOLR with certainty.  Returning now skips the
        // remaining (and most expensive) levels of the DP; the result
        // is identical to running them.  This is what made genuinely
        // NOLR faces burn the full 2^n table before reporting.
        //
        // SUPPRESSED under carry-discharge: `any_live_at_size` is the
        // STRICT (terminal-only) verdict; a face that is Strict-NOLR at
        // some size can still be carried-LR, so exiting here would
        // under-find.  The carry DFS below needs the full set_table.
        if (!do_carry && !any_live_at_size) {
            if (g_lr_trace.level > 0) {
                std::fprintf(stderr,
                    "[lrtrace] early NOLR exit at size=%zu/%zu "
                    "(all subsets NOLR)\n", size, n);
                std::fflush(stderr);
            }
            return LrResult{{}, INF, {}};
        }

        // Memory pruning: drop size-(size-1) set entries, they're no
        // longer referenced at size-(size+1) iterations.
        //
        // SUPPRESSED under carry-discharge: the per-path DFS reads the
        // parent subset's letters at EVERY depth (size 0..n-1), so the
        // whole table must survive the size loop.  In the production
        // per-gauge regime n is the post-gauge variable count (small),
        // so retaining 2^n subset states is acceptable — and the DP
        // already materializes set_table[bits] for all bits at each
        // size before this prune, so suppressing it only changes the
        // table's LIFETIME, not its peak per-size width.
        if (!do_carry && size >= 1) {
            auto prev_subsets = subsets_of_size(n, size - 1);
            for (uint64_t k : prev_subsets) set_table.erase(k);
        }

        if (g_lr_trace.level > 0) {
            // Largest surviving letter set at this size (post-intersect).
            size_t max_set = 0, max_terms = 0;
            for (const auto& kv : set_table) {
                for (const auto& gl : kv.second) {
                    max_set = std::max(max_set, gl.size());
                    for (const auto& p : gl)
                        max_terms = std::max(max_terms, p.n_terms());
                }
            }
            std::fprintf(stderr,
                "[lrtrace] size=%zu/%zu done: steps=%ld res=%ld(%.1fs) "
                "disc=%.1fs lc=%.1fs factor=%ld(%.1fs) dedup=%.1fs "
                "maxset=%zu maxterms=%zu maxresterms=%zu "
                "stepmemo=%ld/%ld facmemo=%ld/%ld\n",
                size, n, g_lr_trace.n_steps, g_lr_trace.n_res,
                g_lr_trace.t_res, g_lr_trace.t_disc, g_lr_trace.t_lc,
                g_lr_trace.n_factor, g_lr_trace.t_factor,
                g_lr_trace.t_dedup, max_set, max_terms,
                g_lr_trace.max_res_terms,
                g_lr_memo.step_hit, g_lr_memo.step_hit + g_lr_memo.step_miss,
                g_lr_memo.fac_hit, g_lr_memo.fac_hit + g_lr_memo.fac_miss);
            if (euler_filter_on) {
                const ChiFilterStats cs = chi_filter_stats();
                std::fprintf(stderr,
                    "[lrtrace]   chi: judged=%lu dropped=%lu "
                    "boundary_exempt=%lu msolve_calls=%lu\n",
                    cs.judged, cs.dropped, cs.boundary_exempt,
                    cs.msolve_calls);
            }
            std::fflush(stderr);
        }
    }

    if (g_lr_trace.level > 0) {
        std::fprintf(stderr, "[lrtrace] slowest step calls:\n");
        for (const auto& c : g_lr_trace.slow) {
            std::fprintf(stderr,
                "[lrtrace]   wall=%.2fs bits=%llx pivot=%zu in=%zu out=%zu\n",
                c.wall, (unsigned long long) c.bits, c.pivot,
                c.n_in, c.n_out);
        }
        std::fflush(stderr);
    }

    // Carry-discharge (Doppio FindRoots) verdict: a per-path DFS over the
    // full (un-pruned) set_table, judging each step with the shared
    // step_fr_judge primitive (gauge-free).  Strictly supersedes the
    // Strict subset-DP result: every Strict-admissible order is also
    // carry-admissible (no obligation is ever generated when every letter
    // is already terminal/conic/deg<=1), so the carry DFS finds at least
    // what the DP found and possibly more.  Score-minimal order wins,
    // matching the DP's MinimalBy.  Empty set_table guards (n==0 etc.)
    // were handled by the early returns at the top of the function.
    if (do_carry) {
        // n >= 1 here (n == 0 returned early).  Seed bits=0.
        CarryDfs drv{set_table, xvar_indices, n, INF};
        std::vector<size_t> order;
        std::vector<Poly> roots;
        lr_scan::PathState st0;
        drv.dfs(0ull, order, 0.0, roots, st0);
        if (!drv.found) return LrResult{{}, INF, {}};
        LrResult out;
        out.order = std::move(drv.best_order);
        out.score = drv.best_score;
        out.root_polys = std::move(drv.best_roots);
        out.carried_sqrts = drv.best_nsq;
        out.kin_sqrts = drv.best_nkin;
        out.terminal_quads = drv.best_ntq;
        return out;
    }

    const uint64_t full = (n == 64) ? ~0ull : ((1ull << n) - 1);
    auto it = orders_table.find(full);
    if (it == orders_table.end()) return LrResult{{}, INF, {}};
    return it->second;
}

}  // namespace lr_search
}  // namespace hyperflint
