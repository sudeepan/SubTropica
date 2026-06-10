// Unit test for the HF_LF_PERFPOW univariate perfect-power detector
// in linear_factors.

#include "hyperflint/algebra/linear_factors.hpp"
#include "hyperflint/core/poly.hpp"
#include "hyperflint/core/rat.hpp"
#include "hyperflint/core/zw_table.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace hyperflint;

namespace {
int g_failures = 0;

void check(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++g_failures;
    }
}
}  // namespace

int main() {
    // Fixture A: p = (x - y)^4 in ctx {x, y}.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly base = x.sub(y);
        Poly p = base.pow(4);

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.nonlinear.empty(), "A: no nonlinear factors");
        check(lf.linear.size() == 1,
              "A: exactly one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 4,
                  "A: multiplicity == 4, got " +
                  std::to_string(lf.linear[0].multiplicity));
            Poly pole_num{y};
            Rat expected_pole{std::move(pole_num)};
            Poly cross1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cross2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cross1.equal(cross2), "A: pole == y");
        }
    }

    // Fixture B: p = (x - y)^2 * (x - z)^2 -- NOT a perfect power.
    {
        PolyCtx ctx({"x", "y", "z"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly z = Poly::gen(ctx, 2);
        Poly p = x.sub(y).pow(2).mul(x.sub(z).pow(2));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.nonlinear.empty(), "B: no nonlinear factors");
        check(lf.linear.size() == 2,
              "B: exactly two linear factors, got " +
              std::to_string(lf.linear.size()));
        for (size_t i = 0; i < lf.linear.size(); ++i) {
            check(lf.linear[i].multiplicity == 2,
                  "B: factor " + std::to_string(i) +
                  " multiplicity == 2, got " +
                  std::to_string(lf.linear[i].multiplicity));
        }
    }

    // Fixture C: p = x^3 -- edge case b=0. Detector skips (c[0]=0).
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly p = x.pow(3);

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "C: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 3,
                  "C: multiplicity == 3");
        }
    }

    // --- Phase-3 fast-path fixtures (HF_LF_PERFPOW_FAST, default ON) ---
    // Fixtures A-C above already exercise the mod-p screen + the d=4
    // sqrt-chain extraction / certified-reject / detector-skip routes.

    // Fixture D: rational-scaled square, p = (3/5)*(x - y)^2. d=2 single
    // sqrt; the 3/5 lives in the fmpq content, zpoly stays (x-y)^2.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = x.sub(y).pow(2)
                     .mul(Poly::from_int(ctx, 3))
                     .divexact(Poly::from_int(ctx, 5));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "D: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 2, "D: multiplicity == 2");
            Rat expected_pole{Poly(y)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "D: pole == y");
        }
    }

    // Fixture E: poly-scaled power, p = (y + 1)*(x - z)^4. The sqrt
    // chain on c_4 = y+1 fails (not a square), so the fast route falls
    // back to the legacy pole construction; the detector must still
    // fire identically (quot = y+1 is var-free in x).
    {
        PolyCtx ctx({"x", "y", "z"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly z = Poly::gen(ctx, 2);
        Poly p = y.add(Poly::one_of(ctx)).mul(x.sub(z).pow(4));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "E: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 4, "E: multiplicity == 4");
            Rat expected_pole{Poly(z)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "E: pole == z");
        }
    }

    // Fixture F: d=8 sqrt chain (three halvings), p = (2x + 3y)^8.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = x.mul(Poly::from_int(ctx, 2))
                     .add(y.mul(Poly::from_int(ctx, 3))).pow(8);

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "F: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 8, "F: multiplicity == 8");
            // pole = -3y/2
            Rat expected_pole{
                Poly::from_int(ctx, -3).mul(y),
                Poly::from_int(ctx, 2)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "F: pole == -3y/2");
        }
    }

    // Fixture G: negative leading coefficient, p = -(x - y)^4 (lc < 0
    // exercises the sqrt-chain sign retry).
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = x.sub(y).pow(4).mul(Poly::from_int(ctx, -1));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "G: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 4, "G: multiplicity == 4");
            Rat expected_pole{Poly(y)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "G: pole == y");
        }
    }

    // Fixture H: d=6 (not in {2,4,8}): fast extraction skipped, legacy
    // route must still detect (x + y)^6.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = x.add(y).pow(6);

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "H: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 6, "H: multiplicity == 6");
        }
    }

    // --- Adversarial-review folds (2026-06-10): sqrt chain consuming a
    // nontrivial SQUARE polynomial cofactor. Only here does s carry a
    // real cofactor, exercising the divexact(s^{d-1}) cancellation.

    // Fixture I: p = (y+1)^2 * (x - y)^2 -- chain consumes (y+1).
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = y.add(Poly::one_of(ctx)).pow(2).mul(x.sub(y).pow(2));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "I: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 2, "I: multiplicity == 2");
            Rat expected_pole{Poly(y)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "I: pole == y");
        }
    }

    // Fixture J: p = y^2 * (y*x + 3)^2 -- cofactor m = y SHARES the
    // variable with A = y: s = y^2, divexact must cancel the shared
    // factor. Pole = -3/y.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = y.pow(2)
                     .mul(y.mul(x).add(Poly::from_int(ctx, 3)).pow(2));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "J: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 2, "J: multiplicity == 2");
            Rat expected_pole{Poly::from_int(ctx, -3), Poly(y)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "J: pole == -3/y");
        }
    }

    // Fixture K: p = (y+1)^4 * (2x - y)^4 -- d=4 chain consuming a
    // 4th-power cofactor. Pole = y/2.
    {
        PolyCtx ctx({"x", "y"});
        Poly x = Poly::gen(ctx, 0);
        Poly y = Poly::gen(ctx, 1);
        Poly p = y.add(Poly::one_of(ctx)).pow(4)
                     .mul(x.mul(Poly::from_int(ctx, 2)).sub(y).pow(4));

        auto zw = std::make_shared<ZWTable>(ctx);
        LinearFactorization lf = linear_factors(p, 0, zw, false);

        check(lf.linear.size() == 1,
              "K: one linear factor, got " +
              std::to_string(lf.linear.size()));
        if (!lf.linear.empty()) {
            check(lf.linear[0].multiplicity == 4, "K: multiplicity == 4");
            Rat expected_pole{Poly(y), Poly::from_int(ctx, 2)};
            Poly cr1 = lf.linear[0].pole.num().mul(expected_pole.den());
            Poly cr2 = expected_pole.num().mul(lf.linear[0].pole.den());
            check(cr1.equal(cr2), "K: pole == y/2");
        }
    }

    if (g_failures == 0) {
        std::cout << "test_lf_perfpow: all passed\n";
    }
    return g_failures > 0 ? 1 : 0;
}
