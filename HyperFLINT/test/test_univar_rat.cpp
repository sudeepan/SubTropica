// Unit tests for UnivarRat: R(y)[x] polynomial arithmetic and
// Cauchy-product power-series inversion.

#include "univar_rat.hpp"
#include "hyperflint/core/poly.hpp"
#include "hyperflint/core/rat.hpp"
#include <iostream>
#include <string>

using namespace hyperflint;

namespace {
int g_failures = 0;

void check(bool cond, const std::string& msg) {
    if (!cond) { std::cerr << "FAIL: " << msg << "\n"; ++g_failures; }
    else       { std::cerr << "  ok: " << msg << "\n"; }
}

// Semantic equality for Rats: a.num*b.den == b.num*a.den.
bool rat_eq(const Rat& a, const Rat& b) {
    return a.num().mul(b.den()).equal(b.num().mul(a.den()));
}
}  // namespace

int main() {
    PolyCtx ctx({"x", "y"});
    Poly x = Poly::gen(ctx, 0);
    Poly y = Poly::gen(ctx, 1);

    // T1: from_poly extracts coefficients in x correctly.
    //   p = y*x^2 + 3*x + 1  =>  c0=1, c1=3, c2=y
    {
        Poly p = y.mul(x.pow(2))
                  .add(Poly::from_int(ctx, 3).mul(x))
                  .add(Poly::one_of(ctx));
        UnivarRat u = UnivarRat::from_poly(p, 0);
        check(u.degree() == 2, "T1 from_poly degree");
        check(rat_eq(u.coeff(0), Rat(Poly::one_of(ctx))), "T1 c0=1");
        check(rat_eq(u.coeff(1), Rat(Poly::from_int(ctx, 3))), "T1 c1=3");
        check(rat_eq(u.coeff(2), Rat(Poly(y))), "T1 c2=y");
        check(rat_eq(u.coeff(5), Rat(Poly::zero_of(ctx))),
              "T1 out-of-range=0");
    }

    // T2: divrem: (x^2 - 1) / (x - 1) = x + 1, remainder 0.
    {
        std::vector<Rat> a_coeffs;
        a_coeffs.push_back(Rat(Poly::from_int(ctx, -1)));
        a_coeffs.push_back(Rat(Poly::zero_of(ctx)));
        a_coeffs.push_back(Rat(Poly::one_of(ctx)));
        UnivarRat a(std::move(a_coeffs), ctx);  // x^2 - 1

        std::vector<Rat> b_coeffs;
        b_coeffs.push_back(Rat(Poly::from_int(ctx, -1)));
        b_coeffs.push_back(Rat(Poly::one_of(ctx)));
        UnivarRat b(std::move(b_coeffs), ctx);  // x - 1

        auto [q, r] = a.divrem(b);
        check(q.degree() == 1, "T2 divrem q degree=1");
        check(rat_eq(q.coeff(0), Rat(Poly::one_of(ctx))), "T2 q c0=1");
        check(rat_eq(q.coeff(1), Rat(Poly::one_of(ctx))), "T2 q c1=1");
        check(r.is_zero(), "T2 rem=0");
    }

    // T3: inverse_mod: 1/(1-x) mod x^3 = 1 + x + x^2.
    {
        std::vector<Rat> f_coeffs;
        f_coeffs.push_back(Rat(Poly::one_of(ctx)));
        f_coeffs.push_back(Rat(Poly::from_int(ctx, -1)));
        UnivarRat f(std::move(f_coeffs), ctx);  // 1 - x

        UnivarRat g = univar_inverse_mod(f, 3);
        check(g.degree() == 2, "T3 inv_mod degree=2");
        check(rat_eq(g.coeff(0), Rat(Poly::one_of(ctx))), "T3 c0=1");
        check(rat_eq(g.coeff(1), Rat(Poly::one_of(ctx))), "T3 c1=1");
        check(rat_eq(g.coeff(2), Rat(Poly::one_of(ctx))), "T3 c2=1");
    }

    // T4: inverse_mod with parametric coefficients: 1/(y - x) mod x^2
    //   = 1/y + x/y^2.
    {
        std::vector<Rat> f_coeffs;
        f_coeffs.push_back(Rat(Poly(y)));
        f_coeffs.push_back(Rat(Poly::from_int(ctx, -1)));
        UnivarRat f(std::move(f_coeffs), ctx);  // y - x

        UnivarRat g = univar_inverse_mod(f, 2);
        Rat inv_y = Rat::one_of(ctx) / Rat(Poly(y));
        Rat inv_y2 = inv_y * inv_y;
        check(rat_eq(g.coeff(0), inv_y), "T4 c0=1/y");
        check(rat_eq(g.coeff(1), inv_y2), "T4 c1=1/y^2");
    }

    // T5: multiplication round-trip: (x+1)*(x-1) = x^2 - 1.
    {
        std::vector<Rat> a_coeffs;
        a_coeffs.push_back(Rat(Poly::one_of(ctx)));
        a_coeffs.push_back(Rat(Poly::one_of(ctx)));
        UnivarRat a(std::move(a_coeffs), ctx);  // x + 1

        std::vector<Rat> b_coeffs;
        b_coeffs.push_back(Rat(Poly::from_int(ctx, -1)));
        b_coeffs.push_back(Rat(Poly::one_of(ctx)));
        UnivarRat b(std::move(b_coeffs), ctx);  // x - 1

        UnivarRat p = a * b;  // x^2 - 1
        check(p.degree() == 2, "T5 mul degree=2");
        check(rat_eq(p.coeff(0), Rat(Poly::from_int(ctx, -1))), "T5 c0=-1");
        check(rat_eq(p.coeff(1), Rat(Poly::zero_of(ctx))), "T5 c1=0");
        check(rat_eq(p.coeff(2), Rat(Poly::one_of(ctx))), "T5 c2=1");
    }

    // T6: truncate.
    {
        Poly p = x.pow(3).add(x.pow(2)).add(x).add(Poly::one_of(ctx));
        UnivarRat u = UnivarRat::from_poly(p, 0);  // 1 + x + x^2 + x^3
        UnivarRat t = u.truncate(2);  // keep only 1 + x
        check(t.degree() == 1, "T6 truncate degree=1");
        check(rat_eq(t.coeff(0), Rat(Poly::one_of(ctx))), "T6 c0=1");
        check(rat_eq(t.coeff(1), Rat(Poly::one_of(ctx))), "T6 c1=1");
    }

    // T7: pow.
    {
        std::vector<Rat> f_coeffs;
        f_coeffs.push_back(Rat(Poly::one_of(ctx)));
        f_coeffs.push_back(Rat(Poly::one_of(ctx)));
        UnivarRat f(std::move(f_coeffs), ctx);  // 1 + x

        UnivarRat f3 = f.pow(3);  // 1 + 3x + 3x^2 + x^3
        check(f3.degree() == 3, "T7 pow degree=3");
        check(rat_eq(f3.coeff(0), Rat(Poly::one_of(ctx))), "T7 c0=1");
        check(rat_eq(f3.coeff(1), Rat(Poly::from_int(ctx, 3))), "T7 c1=3");
        check(rat_eq(f3.coeff(2), Rat(Poly::from_int(ctx, 3))), "T7 c2=3");
        check(rat_eq(f3.coeff(3), Rat(Poly::one_of(ctx))), "T7 c3=1");
    }

    // T8: eval via Horner: f(x) = 2x+3 evaluated at x=5 gives 13.
    {
        PolyCtx ctx1({"t"});
        std::vector<Rat> f_coeffs;
        f_coeffs.push_back(Rat(Poly::from_int(ctx1, 3)));
        f_coeffs.push_back(Rat(Poly::from_int(ctx1, 2)));
        UnivarRat f(std::move(f_coeffs), ctx1);  // 2t + 3

        Rat val = f.eval(Rat::from_int(ctx1, 5));
        check(rat_eq(val, Rat::from_int(ctx1, 13)), "T8 eval 2*5+3=13");
    }

    // T9: inverse_mod round-trip: f * inverse_mod(f, n) mod x^n = 1.
    {
        // f = 2 + 3x (constant term nonzero).
        std::vector<Rat> f_coeffs;
        f_coeffs.push_back(Rat(Poly::from_int(ctx, 2)));
        f_coeffs.push_back(Rat(Poly::from_int(ctx, 3)));
        UnivarRat f(std::move(f_coeffs), ctx);

        long n = 5;
        UnivarRat g = univar_inverse_mod(f, n);
        UnivarRat prod = (f * g).truncate(n);
        // prod should be 1 + 0*x + 0*x^2 + ...
        check(rat_eq(prod.coeff(0), Rat::one_of(ctx)),
              "T9 roundtrip c0=1");
        for (long k = 1; k < n; ++k) {
            check(prod.coeff(k).is_zero(),
                  "T9 roundtrip c" + std::to_string(k) + "=0");
        }
    }

    if (g_failures == 0)
        std::cout << "test_univar_rat: all " << 9 << " tests passed\n";
    else
        std::cerr << "test_univar_rat: " << g_failures << " FAILURES\n";
    return g_failures > 0 ? 1 : 0;
}
