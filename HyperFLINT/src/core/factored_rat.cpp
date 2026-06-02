#include "hyperflint/core/factored_rat.hpp"

#include <algorithm>
#include <stdexcept>

namespace hyperflint {

FactoredRat FactoredRat::from_poly(Poly numerator) {
    return FactoredRat(std::move(numerator));
}

void FactoredRat::push_factor(const Poly& base, long exp) {
    if (exp == 0 || base.is_one()) return;
    for (auto& f : den_factors_) {
        if (f.base.equal(base)) { f.exp += exp; return; }
    }
    den_factors_.push_back(Factor{base, exp});  // first-insertion order
}

FactoredRat FactoredRat::from_rat(const Rat& r) {
    FactoredRat fr(r.num());
    if (!r.den().is_one()) fr.push_factor(r.den(), 1);
    return fr;
}

Poly FactoredRat::expand_denominator() const {
    Poly d = Poly::one_of(numerator_.ctx());
    for (const auto& f : den_factors_) {
        d = d.mul(f.base.pow(static_cast<unsigned long>(f.exp)));
    }
    return d;
}

Rat FactoredRat::materialize_to_rat() const {
    if (den_factors_.empty()) return Rat(numerator_);
    return Rat(numerator_, expand_denominator());  // reducing ctor
}

FactoredRat FactoredRat::neg() const {
    FactoredRat r(numerator_.neg());
    r.den_factors_ = den_factors_;
    return r;
}

bool FactoredRat::is_zero() const { return numerator_.is_zero(); }

FactoredRat FactoredRat::add(const FactoredRat& b) const {
    // Common denominator = per-base MAX exponent over the union of factor
    // sets. Lift each numerator by the missing factor powers, then add.
    // Factors are matched by Poly::equal (no hash); sets are tiny.
    auto find_exp = [](const std::vector<Factor>& v, const Poly& base) -> long {
        for (const auto& f : v) if (f.base.equal(base)) return f.exp;
        return 0;
    };
    std::vector<Factor> common = den_factors_;   // first-insertion order
    for (const auto& bf : b.den_factors_) {
        bool found = false;
        for (auto& cf : common) {
            if (cf.base.equal(bf.base)) {
                cf.exp = std::max(cf.exp, bf.exp); found = true; break;
            }
        }
        if (!found) common.push_back(bf);
    }
    // Lift this->numerator_ by prod base^(common_exp - my_exp).
    Poly lifted_a = numerator_;
    for (const auto& cf : common) {
        long d = cf.exp - find_exp(den_factors_, cf.base);
        if (d > 0) lifted_a = lifted_a.mul(cf.base.pow(static_cast<unsigned long>(d)));
    }
    Poly lifted_b = b.numerator_;
    for (const auto& cf : common) {
        long d = cf.exp - find_exp(b.den_factors_, cf.base);
        if (d > 0) lifted_b = lifted_b.mul(cf.base.pow(static_cast<unsigned long>(d)));
    }
    FactoredRat r(lifted_a.add(lifted_b));
    r.den_factors_ = common;
    return r;
}

FactoredRat FactoredRat::sub(const FactoredRat& b) const {
    return add(b.neg());
}


// mul: numerator product, denominator factor-lists concatenated (dedup by Poly::equal).
FactoredRat FactoredRat::mul(const FactoredRat& b) const {
    FactoredRat r(numerator_.mul(b.numerator_));
    r.den_factors_ = den_factors_;
    // copy ours, then merge b's factors (push_factor dedups by Poly::equal).
    for (const auto& f : b.den_factors_) r.push_factor(f.base, f.exp);
    return r;
}

FactoredRat FactoredRat::reciprocal() const {
    // 1 / (num / prod f^e) = (prod f^e) / num.  Reciprocal of zero is
    // undefined; throw rather than pushing a zero denominator factor that
    // would detonate later at materialize. Mirrors Rat::div's zero contract.
    if (numerator_.is_zero()) {
        throw std::runtime_error("FactoredRat::reciprocal: reciprocal of zero");
    }
    FactoredRat r(expand_denominator());
    if (!numerator_.is_one()) r.push_factor(numerator_, 1);
    return r;
}

FactoredRat FactoredRat::div(const FactoredRat& b) const {
    return mul(b.reciprocal());
}

FactoredRat FactoredRat::derivative(size_t var_idx) const {
    // value = N / prod_i f_i^{e_i}.
    // Empty denominator: d/dx (N) = N'.
    if (den_factors_.empty()) {
        return FactoredRat(numerator_.derivative(var_idx));
    }
    // Fast path: if every factor base is independent of var_idx, then
    // d/dx (N/D) = N'/D with the SAME factors (no exponent bump, no
    // spurious factors). This is the common (and required) case.
    bool any_depends = false;
    for (const auto& f : den_factors_) {
        if (f.base.degree_in_var(var_idx) > 0) { any_depends = true; break; }
    }
    if (!any_depends) {
        FactoredRat r(numerator_.derivative(var_idx));
        r.den_factors_ = den_factors_;   // unchanged exponents
        return r;
    }
    // General path. Build the quotient-rule numerator as a single Poly:
    //   P = prod_i f_i   (each base once, exponent 1)
    //   term1 = N' * P
    //   term2 = N * sum_i ( e_i * f_i' * (P / f_i) )
    //   numerator = term1 - term2
    // Result denominator = { (f_i, e_i + 1) } for all i.
    const PolyCtx& c = numerator_.ctx();
    Poly N = numerator_;
    Poly Nprime = numerator_.derivative(var_idx);

    // P = prod_i base_i.
    Poly P = Poly::one_of(c);
    for (const auto& f : den_factors_) P = P.mul(f.base);

    Poly term1 = Nprime.mul(P);

    // sum_i e_i * f_i' * prod_{j!=i} f_j.
    // prod_{j!=i} f_j is built per i by a left/right prefix style product;
    // factor sets are tiny (letters), so the O(k^2) restart is negligible.
    Poly accum = Poly::zero_of(c);
    for (size_t i = 0; i < den_factors_.size(); ++i) {
        Poly prod_rest = Poly::one_of(c);
        for (size_t j = 0; j < den_factors_.size(); ++j) {
            if (j == i) continue;
            prod_rest = prod_rest.mul(den_factors_[j].base);
        }
        Poly fi_prime = den_factors_[i].base.derivative(var_idx);
        Poly ei = Poly::from_int(c, den_factors_[i].exp);
        accum = accum.add(ei.mul(fi_prime).mul(prod_rest));
    }
    Poly term2 = N.mul(accum);

    FactoredRat r(term1.sub(term2));
    for (const auto& f : den_factors_) r.push_factor(f.base, f.exp + 1);
    return r;
}

FactoredRat FactoredRat::pow(long n) const {
    // n == 0 -> 1 (empty denominator, numerator 1).
    // n  > 0 -> numerator^n with each denominator factor exponent scaled by n
    //           (keeps the denominator factored: no expansion).
    // n  < 0 -> reciprocal().pow(-n).
    if (n < 0) return reciprocal().pow(-n);
    if (n == 0) return FactoredRat::from_poly(Poly::one_of(numerator_.ctx()));
    FactoredRat r(numerator_.pow(static_cast<unsigned long>(n)));
    for (const auto& f : den_factors_) r.push_factor(f.base, f.exp * n);
    return r;
}
}  // namespace hyperflint
