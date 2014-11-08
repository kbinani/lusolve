#pragma once

#include <vector>
#include <cmath>

namespace LUSolve
{

namespace detail
{

template<class Value>
Value abs_(Value const& v)
{
    if (v < 0) {
        return -v;
    } else {
        return v;
    }
}

template<>
float abs_(float const& v)
{
    return std::fabs(v);
}

template<>
double abs_(double const& v)
{
    return std::fabs(v);
}

template<>
long double abs_(long double const& v)
{
    return std::fabs(v);
}

template<class Value, template<class> class Matrix>
std::vector<size_t> ludecomp_(Matrix<Value> & a, size_t n)
{
    Value const zero = 0;
    Value const one = 1;

    std::vector<size_t> ps(n);
    std::vector<Value> scales(n);
    for (size_t i = 0; i < n; ++i) {  // pick up largest(abs. val.) element in each row.
        ps[i] = i;
        Value nrmrow = zero;
        size_t const ixn = i * n;
        for(size_t j = 0; j < n; ++j) {
            Value biggst = abs_(a[ixn + j]);
            if (biggst > nrmrow) {
                nrmrow = biggst;
            }
        }
        if (nrmrow > zero) {
            scales[i] = one / nrmrow;
        } else {
            throw std::runtime_error("Singular matrix");
        }
    }
    size_t const n1 = n - 1;
    for (size_t k = 0; k < n1; ++k) { // Gaussian elimination with partial pivoting.
        Value biggst = zero;
        size_t pividx = 0;
        for (int i = k; i < n; ++i) {
            Value const size = abs_(a[ps[i] * n + k]) * scales[ps[i]];
            if (size > biggst) {
                biggst = size;
                pividx = i;
            }
        }
        if (biggst <= zero) {
            throw std::runtime_error("Singular matrix");
        }
        if (pividx != k) {
            std::swap(ps[k], ps[pividx]);
        }
        Value const pivot = a[ps[k] * n + k];
        for (size_t i = k + 1; i < n; ++i) {
            size_t const psin = ps[i] * n;
            Value const mult = a[psin + k] / pivot;
            a[psin + k] = mult;
            if (mult != zero) {
                size_t const pskn = ps[k] * n;
                for (size_t j = k + 1; j < n; ++j) {
                    a[psin + j] -= mult * a[pskn + j];
                }
            }
        }
    }
    if (a[ps[n1] * n + n1] == zero) {
        throw std::runtime_error("Singular matrix");
    }

    return ps;
}

template<class Value, template<class> class Matrix>
std::vector<Value> lusolve_(Matrix<Value> const& a, std::vector<Value> const& b, std::vector<size_t> const& ps)
{
    size_t const n = ps.size();
    std::vector<Value> x;
    Value const zero = 0;
    for (size_t i = 0; i < n; ++i) {
        Value dot = zero;
        size_t const psin = ps[i] * n;
        for (size_t j = 0; j < i; ++j) {
            dot += a[psin + j] * x[j];
        }
        x.push_back(b[ps[i]] - dot);
    }
    for (size_t i = n - 1; ; --i) {
        Value dot = zero;
        size_t const psin = ps[i] * n;
        for (size_t j = i + 1; j < n; ++j) {
            dot += a[psin + j] * x[j];
        }
        x[i] = (x[i] - dot) / a[psin + i];
        if (i == 0) {
            break;
        }
    }
    return x;
}

template<class Value>
class matrix_proxy
{
public:
    explicit matrix_proxy(std::vector<std::vector<Value> > * m)
        : m_(m)
        , n_(m->size())
    {
        for (size_t i = 0; i < m->size(); ++i) {
            if ((*m)[i].size() != n_) {
                throw std::runtime_error("Matrix size mismatch");
            }
        }
    }

    Value const& operator[](size_t index) const
    {
        return get_(index);
    }

    Value & operator[](size_t index)
    {
        return get_(index);
    }

private:
    Value & get_(size_t index) const
    {
        size_t const i = index / n_;
        size_t const j = index - i * n_;
        return (*m_)[i][j];
    }

private:
    size_t const n_;
    std::vector<std::vector<Value> > * const m_;
};

} // namespace detail

template<class Value>
std::vector<Value> lusolve(std::vector<std::vector<Value> > & a, std::vector<Value> const& b)
{
    size_t const n = b.size();
    if (a.size() != n) {
        throw std::runtime_error("Matrix size mismatch");
    }

    detail::matrix_proxy<Value> prox(&a);
    std::vector<size_t> ps = detail::ludecomp_(prox, n);
    return std::move(detail::lusolve_(prox, b, ps));
}

template<class Value>
std::vector<Value> lusolve(std::vector<Value> & a, std::vector<Value> const& b)
{
    size_t const n = b.size();
    if (a.size() != n * n) {
        throw std::runtime_error("Matrix size mismatch");
    }

    std::vector<size_t> ps = detail::ludecomp_(a, n);
    return std::move(detail::lusolve_(a, b, ps));
}

} // namespace LUSolve
