// UnivarRat: univariate polynomial over R(y)[x].
// See univar_rat.hpp for the type contract and design notes.

#include "univar_rat.hpp"
#include <algorithm>
#include <stdexcept>

namespace hyperflint {

// --- Construction ---

UnivarRat::UnivarRat(std::vector<Rat> coeffs, const PolyCtx& ctx)
    : c_(std::move(coeffs)),
      ctx_(&ctx),
      zero_cache_(Poly::zero_of(ctx))
{
    normalize();
}

UnivarRat::UnivarRat(const PolyCtx& ctx)
    : ctx_(&ctx),
      zero_cache_(Poly::zero_of(ctx))
{}

UnivarRat UnivarRat::from_poly(const Poly& p, size_t var_idx) {
    const PolyCtx& ctx = p.ctx();
    long d = p.degree_in_var(var_idx);
    if (d < 0) return UnivarRat(ctx);
    std::vector<Rat> c;
    c.reserve(static_cast<size_t>(d + 1));
    for (long k = 0; k <= d; ++k)
        c.push_back(Rat(p.coefficient_of_var(var_idx, k)));
    return UnivarRat(std::move(c), ctx);
}

// --- Queries ---

long UnivarRat::degree() const {
    return static_cast<long>(c_.size()) - 1;
}

const Rat& UnivarRat::coeff(long k) const {
    if (k < 0 || static_cast<size_t>(k) >= c_.size())
        return zero_cache_;
    return c_[static_cast<size_t>(k)];
}

const Rat& UnivarRat::lead() const { return c_.back(); }

bool UnivarRat::is_zero() const { return c_.empty(); }

void UnivarRat::normalize() {
    while (!c_.empty() && c_.back().is_zero())
        c_.pop_back();
}

// --- Arithmetic ---

UnivarRat UnivarRat::operator+(const UnivarRat& b) const {
    size_t n = std::max(c_.size(), b.c_.size());
    std::vector<Rat> r;
    r.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        const Rat& ai = (i < c_.size()) ? c_[i] : zero_cache_;
        const Rat& bi = (i < b.c_.size()) ? b.c_[i] : b.zero_cache_;
        r.push_back(ai + bi);
    }
    return UnivarRat(std::move(r), *ctx_);
}

UnivarRat UnivarRat::operator-(const UnivarRat& b) const {
    size_t n = std::max(c_.size(), b.c_.size());
    std::vector<Rat> r;
    r.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        const Rat& ai = (i < c_.size()) ? c_[i] : zero_cache_;
        const Rat& bi = (i < b.c_.size()) ? b.c_[i] : b.zero_cache_;
        r.push_back(ai - bi);
    }
    return UnivarRat(std::move(r), *ctx_);
}

UnivarRat UnivarRat::operator*(const UnivarRat& b) const {
    if (is_zero() || b.is_zero()) return UnivarRat(*ctx_);
    size_t n = c_.size() + b.c_.size() - 1;
    // Build a zero-filled result vector.
    std::vector<Rat> r;
    r.reserve(n);
    for (size_t i = 0; i < n; ++i)
        r.push_back(Rat(Poly::zero_of(*ctx_)));
    for (size_t i = 0; i < c_.size(); ++i) {
        if (c_[i].is_zero()) continue;
        for (size_t j = 0; j < b.c_.size(); ++j) {
            if (b.c_[j].is_zero()) continue;
            r[i + j] = r[i + j] + c_[i] * b.c_[j];
        }
    }
    return UnivarRat(std::move(r), *ctx_);
}

UnivarRat UnivarRat::operator*(const Rat& s) const {
    if (s.is_zero()) return UnivarRat(*ctx_);
    std::vector<Rat> r;
    r.reserve(c_.size());
    for (const auto& ci : c_)
        r.push_back(ci * s);
    return UnivarRat(std::move(r), *ctx_);
}

UnivarRat UnivarRat::operator/(const Rat& s) const {
    std::vector<Rat> r;
    r.reserve(c_.size());
    for (const auto& ci : c_)
        r.push_back(ci / s);
    return UnivarRat(std::move(r), *ctx_);
}

// --- Division ---

std::pair<UnivarRat, UnivarRat>
UnivarRat::divrem(const UnivarRat& b) const {
    if (b.is_zero())
        throw std::runtime_error("UnivarRat::divrem: division by zero");
    if (is_zero() || degree() < b.degree())
        return {UnivarRat(*ctx_), UnivarRat(std::vector<Rat>(c_), *ctx_)};
    long dq = degree() - b.degree();
    // Quotient coefficients, zero-initialized.
    std::vector<Rat> q;
    q.reserve(static_cast<size_t>(dq + 1));
    for (long i = 0; i <= dq; ++i)
        q.push_back(Rat(Poly::zero_of(*ctx_)));
    // Working copy of dividend.
    std::vector<Rat> r(c_);
    const Rat& blead = b.lead();
    for (long i = dq; i >= 0; --i) {
        long ri = i + b.degree();
        if (r[static_cast<size_t>(ri)].is_zero()) continue;
        Rat qi = r[static_cast<size_t>(ri)] / blead;
        q[static_cast<size_t>(i)] = qi;
        for (long j = 0; j <= b.degree(); ++j) {
            r[static_cast<size_t>(i + j)] =
                r[static_cast<size_t>(i + j)] - qi * b.c_[static_cast<size_t>(j)];
        }
    }
    return {UnivarRat(std::move(q), *ctx_),
            UnivarRat(std::move(r), *ctx_)};
}

UnivarRat UnivarRat::rem(const UnivarRat& b) const {
    return divrem(b).second;
}

// --- Power / truncate / eval ---

UnivarRat UnivarRat::pow(long n) const {
    if (n == 0) {
        std::vector<Rat> one;
        one.push_back(Rat::one_of(*ctx_));
        return UnivarRat(std::move(one), *ctx_);
    }
    if (n == 1) return *this;
    UnivarRat half = pow(n / 2);
    UnivarRat result = half * half;
    if (n % 2 == 1) result = result * (*this);
    return result;
}

UnivarRat UnivarRat::truncate(long n) const {
    if (n <= 0) return UnivarRat(*ctx_);
    size_t m = std::min(c_.size(), static_cast<size_t>(n));
    std::vector<Rat> r(c_.begin(), c_.begin() + static_cast<long>(m));
    return UnivarRat(std::move(r), *ctx_);
}

Rat UnivarRat::eval(const Rat& x) const {
    if (is_zero()) return Rat(Poly::zero_of(*ctx_));
    Rat result = c_.back();
    for (long i = static_cast<long>(c_.size()) - 2; i >= 0; --i)
        result = result * x + c_[static_cast<size_t>(i)];
    return result;
}

// --- Cauchy-product power-series inversion mod x^n ---
//
// Given f(x) with f(0) ≠ 0, computes g(x) = f(x)^{-1} mod x^n via:
//   g[0] = 1/f[0]
//   g[j] = -(1/f[0]) * Σ_{i=1}^{j} f[i]*g[j-i]   for j = 1..n-1.
//
// Cost: O(n²) Rat multiplications, each of which is a multivariate
// rational function operation.  For the multiplicities we encounter
// (m ≤ 20), this is negligible compared to the FactoredRat derivative
// chain it replaces.

UnivarRat univar_inverse_mod(const UnivarRat& f, long n) {
    if (f.is_zero() || f.coeff(0).is_zero())
        throw std::runtime_error("univar_inverse_mod: f(0) == 0");
    const PolyCtx& ctx = f.ctx();
    Rat inv_f0 = Rat::one_of(ctx) / f.coeff(0);
    std::vector<Rat> g;
    g.reserve(static_cast<size_t>(n));
    g.push_back(inv_f0);
    for (long j = 1; j < n; ++j) {
        Rat acc(Poly::zero_of(ctx));
        for (long i = 1; i <= j; ++i) {
            if (i > f.degree()) break;
            acc = acc + f.coeff(i) * g[static_cast<size_t>(j - i)];
        }
        g.push_back(Rat(Poly::zero_of(ctx)) - inv_f0 * acc);
    }
    return UnivarRat(std::move(g), ctx);
}

}  // namespace hyperflint
