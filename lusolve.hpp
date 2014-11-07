#pragma once

#include <vector>

template<class T>
std::vector<size_t> ludecomp_internal_(std::vector<T> & a, size_t n)
{
    T const zero = 0;
    T const one = 1;

    auto abs = [zero](T const& v) {
        if (v < zero) {
            return -v;
        } else {
            return v;
        }
    };

    std::vector<size_t> ps;
    std::vector<T> scales;
    for (size_t i = 0; i < n; ++i) {  // pick up largest(abs. val.) element in each row.
        ps.push_back(i);
        T nrmrow = zero;
        size_t const ixn = i * n;
        for(size_t j = 0; j < n; ++j) {
            T biggst = abs(a[ixn + j]);
            if (biggst > nrmrow) {
                nrmrow = biggst;
            }
        }
        if (nrmrow > zero) {
            scales.push_back(one / nrmrow);
        } else {
            throw std::runtime_error("Singular matrix");
        }
    }
    size_t const n1 = n - 1;
    for (size_t k = 0; k < n1; ++k) { // Gaussian elimination with partial pivoting.
        T biggst = zero;
        size_t pividx = 0;
        for (int i = k; i < n; ++i) {
            T const size = abs(a[ps[i] * n + k]) * scales[ps[i]];
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
        T const pivot = a[ps[k] * n + k];
        for (size_t i = k + 1; i < n; ++i) {
            size_t const psin = ps[i] * n;
            T const mult = a[psin + k] / pivot;
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

template<class T>
std::vector<T> lusolve_internal_(std::vector<T> const& a, std::vector<T> const& b, std::vector<size_t> const& ps)
{
    size_t const n = ps.size();
    std::vector<T> x;
    T const zero = 0;
    for (size_t i = 0; i < n; ++i) {
        T dot = zero;
        size_t const psin = ps[i] * n;
        for (size_t j = 0; j < i; ++j) {
            dot += a[psin + j] * x[j];
        }
        x.push_back(b[ps[i]] - dot);
    }
    for (size_t i = n - 1; ; --i) {
        T dot = zero;
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

template<class T>
std::vector<T> lusolve(std::vector<std::vector<T>> const& a_, std::vector<T> const& b)
{
    std::vector<T> a;
    size_t const n = a_.size();
    a.reserve(n * n);
    for (auto const& s : a_) {
        for (T const& v : s) {
            a.push_back(v);
        }
    }
    std::vector<size_t> ps = ludecomp_internal_(a, n);
    return std::move(lusolve_internal_(a, b, ps));
}
