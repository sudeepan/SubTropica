// Unit tests for the factor-prediction table
// (spec docs/superpowers/specs/2026-06-11-stfactorpredictor-design.md)
// and the shared lr_letter_admissible helper extracted from the LR DP.
//
// Style matches test/unit/test_gcd_dispatch.cpp: raw main(), assert(),
// PASS lines on stdout, exit 0 on success.

#include "hyperflint/core/poly.hpp"
#include "hyperflint/integrator/factor_table.hpp"
#include "hyperflint/integrator/lr_search.hpp"

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace ft = hyperflint::factor_table;

using hyperflint::Poly;
using hyperflint::PolyCtx;

// ---------------------------------------------------------------------
// lr_letter_admissible: degree window + deg-2 forbidden-variable guard.
// The deg-2 branch is the only consumer of the forbidden set, so this
// test MUST exercise max_deg = 2 (a deg-1-only fixture would
// false-green the parity assertion; spec section 6.1c).
// ---------------------------------------------------------------------
static void test_admissible() {
    PolyCtx ctx({"x", "z", "y"});
    const size_t ix = ctx.index_of("x"), iz = ctx.index_of("z");
    Poly lin(ctx, "x+y");          // deg-1 in x
    Poly quad(ctx, "y*x^2+x+1");   // deg-2 in x, no z
    Poly quadz(ctx, "z*x^2+x+1");  // deg-2 in x, uses future var z
    Poly cubic(ctx, "x^3+y");      // deg-3 in x
    Poly free0(ctx, "y+1");        // deg-0 in x
    std::vector<size_t> forb = {iz};
    using hyperflint::lr_search::lr_letter_admissible;
    assert(lr_letter_admissible(lin, ix, forb, 1));
    assert(!lr_letter_admissible(quad, ix, forb, 1));   // deg-2 under max_deg=1
    assert(lr_letter_admissible(quad, ix, forb, 2));    // FindRoots branch
    assert(!lr_letter_admissible(quadz, ix, forb, 2));  // forbidden dep
    assert(lr_letter_admissible(quadz, ix, {}, 2));     // empty forbidden set
    assert(!lr_letter_admissible(cubic, ix, forb, 2));
    assert(!lr_letter_admissible(free0, ix, forb, 2));  // deg-0 is not a letter
    std::puts("test_admissible PASS");
}

// ---------------------------------------------------------------------
// Exactness contracts (spec 4.4), asserted on EVERY emitted entry:
//   pair:      N * den_prod == c * num_prod * lcF * lcG
//              (cross-multiplied form of F/lc - G/lc == c * Prod P^n)
//   + the deg-1 Res identity N == -Res_var(F, G);
//   singleton: coefficient * den_prod == c * num_prod.
// Unconditional, independent of the disc/resultant conjecture.
// ---------------------------------------------------------------------
static void split_products(const ft::FactorTable& t,
                           const ft::FactoredObject& fo,
                           const hyperflint::PolyCtx& ctx,
                           hyperflint::Poly& num_prod,
                           hyperflint::Poly& den_prod) {
    num_prod = hyperflint::Poly::one_of(ctx);
    den_prod = hyperflint::Poly::one_of(ctx);
    for (const auto& [id, e] : fo.factors) {
        const hyperflint::Poly& P = t.intern_polys[id];
        if (e > 0)
            for (long i = 0; i < e; ++i) num_prod = num_prod * P;
        else
            for (long i = 0; i < -e; ++i) den_prod = den_prod * P;
    }
}

static void assert_pair_contract(const ft::FactorTable& t,
                                 const ft::PairEntry& pe,
                                 const hyperflint::PolyCtx& ctx) {
    const Poly& F = t.intern_polys[pe.f_id];
    const Poly& G = t.intern_polys[pe.g_id];
    Poly lcF = F.coefficient_of_var(pe.var_idx, 1);
    Poly lcG = G.coefficient_of_var(pe.var_idx, 1);
    Poly N = lcG * F - lcF * G;
    // Res identity (sign fixed at minus).
    Poly res = F.resultant(G, pe.var_idx);
    assert((N + res).is_zero());
    Poly num_prod = Poly::one_of(ctx), den_prod = Poly::one_of(ctx);
    split_products(t, pe.diff, ctx, num_prod, den_prod);
    Poly cpoly(ctx, pe.diff.c);
    Poly lhs = N * den_prod;
    Poly rhs = cpoly * num_prod * lcF * lcG;
    assert(lhs.equal(rhs));
}

static void assert_singleton_contract(const ft::FactorTable& t,
                                      const ft::SingletonEntry& se,
                                      const hyperflint::PolyCtx& ctx) {
    const Poly L = t.intern_polys[se.id];
    for (const auto& sc : se.coeffs) {
        Poly target = L.coefficient_of_var(se.var_idx, sc.power);
        Poly num_prod = Poly::one_of(ctx), den_prod = Poly::one_of(ctx);
        split_products(t, sc.fo, ctx, num_prod, den_prod);
        Poly cpoly(ctx, sc.fo.c);
        assert((target * den_prod).equal(cpoly * num_prod));
    }
    if (se.has_disc) {
        Poly target = L.discriminant_in_var(se.var_idx);
        Poly num_prod = Poly::one_of(ctx), den_prod = Poly::one_of(ctx);
        split_products(t, se.disc, ctx, num_prod, den_prod);
        Poly cpoly(ctx, se.disc.c);
        assert((target * den_prod).equal(cpoly * num_prod));
    }
}

static void assert_all_contracts(const ft::FactorTable& t,
                                 const hyperflint::PolyCtx& ctx) {
    for (const auto& pe : t.pairs) assert_pair_contract(t, pe, ctx);
    for (const auto& se : t.singletons) assert_singleton_contract(t, se, ctx);
}

// Fixture (a): the verification example {1+xy, x+y}, order {x, y}.
// Pair entry must be -(y-1)(y+1)/y up to convention; no pair fallback.
static void test_fixture_a() {
    PolyCtx ctx({"x", "y"});
    std::vector<std::vector<Poly>> groups{
        {Poly(ctx, "1+x*y"), Poly(ctx, "x+y")}};
    ft::FactorTable t = ft::build(
        ctx, groups, {ctx.index_of("x"), ctx.index_of("y")}, false,
        ft::Limits{});
    assert(t.stats.pair_fallbacks == 0);
    assert(t.stages.size() == 2);
    assert(t.stages[0].n_pairs == 1);
    const ft::PairEntry& pe = t.pairs[0];
    assert(!pe.diff.oop);
    bool found_ym1 = false, found_yp1 = false, found_y_neg = false;
    Poly ym1 = Poly(ctx, "y-1").canonical_prop_form();
    Poly yp1 = Poly(ctx, "y+1").canonical_prop_form();
    Poly yy = Poly(ctx, "y").canonical_prop_form();
    for (const auto& [id, e] : pe.diff.factors) {
        const Poly& p = t.intern_polys[id];
        if (p.equal(ym1) && e == 1) found_ym1 = true;
        if (p.equal(yp1) && e == 1) found_yp1 = true;
        if (p.equal(yy) && e == -1) found_y_neg = true;
    }
    assert(found_ym1 && found_yp1 && found_y_neg);
    assert_all_contracts(t, ctx);
    std::puts("test_fixture_a PASS");
}

// Fixture (b): out-of-pool singleton -- trailing coefficient y^2+y+1
// is irreducible over Q and absent from the (empty-ish) pool.
static void test_fixture_b_oop_singleton() {
    PolyCtx ctx({"x", "y"});
    std::vector<std::vector<Poly>> groups{{Poly(ctx, "x + y^2 + y + 1")}};
    ft::FactorTable t = ft::build(
        ctx, groups, {ctx.index_of("x"), ctx.index_of("y")}, false,
        ft::Limits{});
    bool oop_trailing = false;
    for (const auto& se : t.singletons)
        for (const auto& sc : se.coeffs)
            if (sc.power == 0 && sc.fo.oop) oop_trailing = true;
    assert(oop_trailing);
    assert(t.stats.oop >= 1);
    assert_all_contracts(t, ctx);
    std::puts("test_fixture_b PASS");
}

// Fixture (c): deg-2 letters under algebraic_letters=true -- excluded
// from pairs, present as singletons with disc; the forbidden-variable
// guard (deg-2 branch, spec 6.1c) rejects the z-dependent quadratic.
static void test_fixture_c_deg2() {
    PolyCtx ctx({"x", "z", "y"});
    std::vector<std::vector<Poly>> groups{{Poly(ctx, "y*x^2+x+1"),
                                           Poly(ctx, "z*x^2+x+1"),
                                           Poly(ctx, "x+y")}};
    ft::FactorTable t = ft::build(
        ctx, groups,
        {ctx.index_of("x"), ctx.index_of("z"), ctx.index_of("y")}, true,
        ft::Limits{});
    assert(t.stages[0].n_inadmissible == 1);
    assert(t.stages[0].n_pairs == 0);  // only one deg-1 letter at stage x
    bool disc_found = false;
    for (const auto& se : t.singletons)
        if (se.var_idx == ctx.index_of("x") && se.deg == 2 && se.has_disc)
            disc_found = true;
    assert(disc_found);
    assert_all_contracts(t, ctx);
    std::puts("test_fixture_c PASS");
}

// Fixture (d): proportional-letter dedup within a group; identical
// groups share entries (one deduped pair total).
static void test_fixture_d_dedup() {
    PolyCtx ctx({"x", "y"});
    std::vector<std::vector<Poly>> groups{
        {Poly(ctx, "x+y"), Poly(ctx, "2*x+2*y"), Poly(ctx, "1+x*y")},
        {Poly(ctx, "x+y"), Poly(ctx, "1+x*y")}};
    ft::FactorTable t = ft::build(
        ctx, groups, {ctx.index_of("x"), ctx.index_of("y")}, false,
        ft::Limits{});
    assert(t.stages[0].n_pairs == 1);
    assert_all_contracts(t, ctx);
    std::puts("test_fixture_d PASS");
}

// Fixture (e): loud guard -- max_pairs exceeded names the guard.
static void test_fixture_e_guards() {
    PolyCtx ctx({"x", "y"});
    std::vector<std::vector<Poly>> groups{
        {Poly(ctx, "x+y"), Poly(ctx, "x+2*y"), Poly(ctx, "x+3*y")}};
    ft::Limits lim;
    lim.max_pairs = 2;  // 3 pairs at stage x > 2
    bool threw = false;
    try {
        ft::build(ctx, groups, {ctx.index_of("x"), ctx.index_of("y")},
                  false, lim);
    } catch (const std::runtime_error& e) {
        threw = true;
        assert(std::string(e.what()).find("max_pairs") != std::string::npos);
    }
    assert(threw);
    std::puts("test_fixture_e PASS");
}

// Zero monic difference: F = y*(x+y) and G = x+y have distinct
// canonical reps but proportional monic forms; the entry is exact with
// c = "0" and no factors.
static void test_zero_difference() {
    PolyCtx ctx({"x", "y"});
    std::vector<std::vector<Poly>> groups{
        {Poly(ctx, "x*y+y^2"), Poly(ctx, "x+y")}};
    ft::FactorTable t = ft::build(
        ctx, groups, {ctx.index_of("x"), ctx.index_of("y")}, false,
        ft::Limits{});
    assert(t.stages[0].n_pairs == 1);
    assert(t.pairs[0].diff.c == "0");
    assert(t.pairs[0].diff.factors.empty());
    assert_all_contracts(t, ctx);
    std::puts("test_zero_difference PASS");
}

int main() {
    test_admissible();
    test_fixture_a();
    test_fixture_b_oop_singleton();
    test_fixture_c_deg2();
    test_fixture_d_dedup();
    test_fixture_e_guards();
    test_zero_difference();
    std::puts("ALL PASS");
    return 0;
}
