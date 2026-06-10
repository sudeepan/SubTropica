#pragma once

// UnivarRat: univariate polynomial over the Rat field R(y)[x].
//
// Represents f(x) = c₀ + c₁x + c₂x² + ... where each cₖ ∈ R(y₁,...,yₙ),
// the field of rational functions in the remaining variables.  Built for
// the CRT partial fraction decomposition that replaces the FactoredRat
// derivative chain on multi-pole denominators (1m-tbox step 3: 826s → O(10s)).
//
// Design decisions:
//   - Dense representation (std::vector<Rat>): the polynomials we operate on
//     are low-degree in the integration variable (typically deg ≤ 20), so
//     sparsity is not a concern.
//   - The Rat type has no default constructor, so all vector operations use
//     reserve + push_back or explicit Rat(Poly::zero_of(ctx)) fill.  The
//     PolyCtx pointer is stashed at construction so we can manufacture
//     zeros on demand.
//   - Trailing zeros are stripped by normalize() after every construction.

#include "hyperflint/core/poly.hpp"
#include "hyperflint/core/rat.hpp"
#include <utility>
#include <vector>

namespace hyperflint {

class UnivarRat {
public:
    // Construct from coefficient vector (c[0] = constant term).
    // Requires ctx for manufacturing zero Rats when needed.
    UnivarRat(std::vector<Rat> coeffs, const PolyCtx& ctx);

    // Construct zero polynomial in the given ring.
    explicit UnivarRat(const PolyCtx& ctx);

    // Extract from a multivariate Poly: c[k] = coefficient_of_var(var_idx, k).
    static UnivarRat from_poly(const Poly& p, size_t var_idx);

    // Degree (-1 for zero polynomial).
    long degree() const;

    // Access coefficient.  Returns zero Rat if k is out of range.
    const Rat& coeff(long k) const;

    // Leading coefficient (undefined on zero polynomial; caller must check).
    const Rat& lead() const;

    bool is_zero() const;

    const PolyCtx& ctx() const { return *ctx_; }

    // Arithmetic.
    UnivarRat operator+(const UnivarRat& b) const;
    UnivarRat operator-(const UnivarRat& b) const;
    UnivarRat operator*(const UnivarRat& b) const;

    // Scalar operations.
    UnivarRat operator*(const Rat& s) const;
    UnivarRat operator/(const Rat& s) const;

    // Polynomial division: returns {quotient, remainder}.
    // Requires b nonzero.
    std::pair<UnivarRat, UnivarRat> divrem(const UnivarRat& b) const;

    // Remainder only.
    UnivarRat rem(const UnivarRat& b) const;

    // Power: f^n (n >= 0).
    UnivarRat pow(long n) const;

    // Truncate to degree < n (keep only coefficients 0..n-1).
    UnivarRat truncate(long n) const;

    // Evaluate at a Rat point: Horner scheme.
    Rat eval(const Rat& x) const;

    // Strip trailing zero coefficients.
    void normalize();

    // Raw coefficient access.
    const std::vector<Rat>& coefficients() const { return c_; }

private:
    std::vector<Rat> c_;       // c_[k] = coefficient of x^k
    const PolyCtx* ctx_;       // context for manufacturing zeros
    mutable Rat zero_cache_;   // lazily-constructed zero Rat for coeff() returns
};

// Cauchy-product power-series inversion:
// Given f with f.coeff(0) nonzero, compute g such that f*g ≡ 1 mod x^n.
// Returns the first n coefficients of g.
UnivarRat univar_inverse_mod(const UnivarRat& f, long n);

}  // namespace hyperflint
