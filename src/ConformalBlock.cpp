#include "common.h"
#include "ConformalBlock.h"
#include "gegenbauer_polynomial.hpp"
#include <iostream>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace boost::math;

float_type& RCoefficients::GetCachedRI(int type, int n) {
    return this->cacheRI[type * MaxOrder + n -1];
}

float_type& RCoefficients::GetCachedRII(int type, int l, int n) {
    boost::unordered_map<int, float_type>::iterator it = cacheR.find(type * 64 * 64 + n * 64 + l);
    if (it != cacheR.end()) {
        return it->second;
    } else {
        cacheR[type * 64 * 64 + n * 64 + l] = std::numeric_limits<float_type>::max();
        return cacheR[type * 64 * 64 + n * 64 + l];
    }
}

float_type RCoefficients::RI(int n) {
    float_type& ret = GetCachedRI(0, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RI_NonDelta(n) * rising_factorial((dlt12 + 1.0 - n) / 2, n) * rising_factorial((dlt34 + 1.0 - n) / 2, n);
    return ret;
}

// the delta independent part of RI.
float_type RCoefficients::RI_NonDelta(int n) {
    float_type fac = factorial<float_type>(n);
    i64 p = ((i64)1) << n;
    float_type ret = p / (fac * fac) * n;

    if (n % 2 == 0) return -ret;
    return ret;
}

float_type RCoefficients::RI_dDelta12(int n) {
    float_type& ret = GetCachedRI(1, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RI_dDelta(n, dlt12, dlt34);
    return ret;
}

float_type RCoefficients::RI_dDelta34(int n) {
    float_type& ret = GetCachedRI(2, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RI_dDelta(n, dlt34, dlt12);
    return ret;
}

float_type RCoefficients::RI_dDelta(int n, float_type deltaDiffA, float_type deltaDiffB) {
    if (!IsInteger((deltaDiffA + 1 - n) / 2)) {
        return digammaDiff((deltaDiffA + 1 - n) / 2, n) * RI(n) * .5;
    } else {
        int x = (int)floor((deltaDiffA + 1 - n) / 2 + EPS);
        if (x > 0 || x + n <= 0) {
            return digammaDiff((deltaDiffA + 1 - n) / 2, n) * RI(n) * .5;
        } else {                
            // -n + 1 <= x <= 0 need special handling because digammaDiff() has a pole and RI() is zero.
            // now, RI() ~ (x)_n, the (x)_n factor has a pole
            float_type ret = rising_factorial((deltaDiffB + 1.0 - n) / 2, n) * 0.5;
            ret *= factorial<float_type>(-x) * factorial<float_type>(x + n - 1);
            if (x % 2 != 0) ret = -ret;

            return ret * RI_NonDelta(n);
        }
    }
}

float_type RCoefficients::RII(int l, int n) {
    float_type& ret = GetCachedRII(0, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RII_NonDelta(l, n) * rising_factorial((dlt12 + 1.0 - n) / 2, n) * rising_factorial((dlt34 + 1.0 - n) / 2, n);
    return ret;
}

// the delta independent part of RII.
float_type RCoefficients::RII_NonDelta(int l, int n) {
    float_type ret = binomial_coefficient<float_type>(l, n) * n;
    ret /= factorial<float_type>(n);
    ret /= (((i64)1) << n);

    float_type h = d / 2;
    if (abs(h - 1) < EPS && l == n) {
        // if h == 1 and l == n, zero factors of numerator and denominator cancel out.
        ret /= rising_factorial(h + l - n, n);
        
    } else {
        ret /= rising_factorial(h + l - n, n) * rising_factorial(h + l - n - 1, n);
        ret *= rising_factorial(d + l - n - 2, n);
    }

    if (n % 2 == 0) return -ret;
    return ret;
}

float_type RCoefficients::RII_dDelta12(int l, int n) {
    float_type& ret = GetCachedRII(1, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RII_dDelta(l, n, dlt12, dlt34);
    return ret;
}

float_type RCoefficients::RII_dDelta34(int l, int n) {
    float_type& ret = GetCachedRII(2, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RII_dDelta(l, n, dlt34, dlt12);
    return ret;
}

float_type RCoefficients::RII_dDelta(int l, int n, float_type deltaDiffA, float_type deltaDiffB) {
    if (!IsInteger((deltaDiffA + 1 - n) / 2)) {
        return digammaDiff((deltaDiffA + 1 - n) / 2, n) * RII(l, n) * .5;
    } else {
        int x = (int)floor((deltaDiffA + 1 - n) / 2 + EPS);
        if (x > 0 || x + n <= 0) {
            return digammaDiff((deltaDiffA + 1 - n) / 2, n) * RII(l, n) * .5;
        } else {                
            // -n + 1 <= x <= 0 need special handling because digammaDiff() has a pole and RII() is zero.
            // now, RII() ~ (x)_n, the (x)_n factor has a pole
            float_type ret = rising_factorial((deltaDiffB + 1.0 - n) / 2, n) * 0.5;
            ret *= factorial<float_type>(-x) * factorial<float_type>(x + n - 1);
            if (x % 2 != 0) ret = -ret;

            return ret * RII_NonDelta(l, n);
        }
    }
}

float_type RCoefficients::RIII(int l, int n) {
    float_type& ret = GetCachedRII(3, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    float_type h = d / 2;
    ret = rising_factorial((dlt12 + h + l - n) / 2, n) * rising_factorial((dlt12 - h - l - n + 2) / 2, n);
    ret *= rising_factorial((dlt34 + h + l - n) / 2, n) * rising_factorial((dlt34 - h - l - n + 2) / 2, n) * RIII_NonDelta(l, n);
    return ret;
}

float_type RCoefficients::RIII_NonDelta(int l, int n) {
    float_type h = d / 2;
    float_type ret;    
    // if h is not integer, there is no poles.
    if (!IsInteger(h)) {
        ret = rising_factorial(h - n - 1, 2 * n);
        ret /= rising_factorial(h + l - n - 1, 2 * n) * rising_factorial(h + l - n, 2 * n);
    } else {
        int ih = (int)floor(h + EPS);
        if (n >= ih - 1 && n < ih + l - 1) { 
            // in this case, the numerate part is zero, and denominator is nonzero
            return .0;
        } else if (2 * n >= l + 1 && n == ih + l - 1) {
            // in this case, the zeros of numerate and denominate cancel out.
            ret = 1.0 / binomial_coefficient<float_type>(2 * n - 1, l);
            ret /= rising_factorial(h + l - n, 2 * n);
            if (l % 2 == 1) ret = -ret;
        } else {
            return 0.0; // temperary return 0 because the following workaround leads to divergent results.

            // n > h + l - 1 case, numerator is zero but denominator is float_type-zero
            // to avoid infinity, we let h -> h + EPS
/*            float_type hh = h - EPS;
            float_type inf = rising_factorial(hh - n - 1, 2 * n) / rising_factorial(hh + l - n - 1, 2 * n);
            inf /=  rising_factorial(hh + l - n, 2 * n);
            ret *= inf;*/
        }
    }

    float_type fac = factorial<float_type>(n);
    ret *= n / (fac * fac);

    if (n % 2 == 0) return -ret;
    return ret;
}

float_type RCoefficients::RIII_dDelta12(int l, int n) {
    float_type& ret = GetCachedRII(4, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RIII_dDelta(l, n, dlt12, dlt34);
    return ret;
}

float_type RCoefficients::RIII_dDelta34(int l, int n) {
    float_type& ret = GetCachedRII(5, l, n);
    if (ret != std::numeric_limits<float_type>::max()) return ret;
    ret = RIII_dDelta(l, n, dlt34, dlt12);
    return ret;
}

float_type RCoefficients::RIII_dDelta(int l, int n, float_type deltaDiffA, float_type deltaDiffB) {
    float_type h = d / 2;
    float_type ret1;
    if (!IsInteger((deltaDiffA - h - l - n + 2) / 2)) {
        ret1 = digammaDiff((deltaDiffA - h - l - n + 2) / 2, n) * RIII(l, n);
    } else {
        int x = (int)floor((deltaDiffA - h - l - n + 2) / 2 + EPS);
        if (x > 0 || x + n <= 0) {
            ret1 = digammaDiff((deltaDiffA - h - l - n + 2) / 2, n) * RIII(l, n);
        } else {                
            // -n + 1 <= x <= 0 need special handling because digammaDiff() has a pole and RII() is zero.
            // now, RII() ~ (x)_n, the (x)_n factor has a pole
            ret1 = rising_factorial((deltaDiffA + h + l - n) / 2, n);
            ret1 *= rising_factorial((deltaDiffB + h + l - n) / 2, n) * rising_factorial((deltaDiffB - h - l - n + 2) / 2, n);
            ret1 *= factorial<float_type>(-x) * factorial<float_type>(x + n - 1);
            if (x % 2 != 0) ret1 = -ret1;

            ret1 *= RIII_NonDelta(l, n);
        }
    }

    float_type ret2;
    if (!IsInteger((deltaDiffA + h + l - n) / 2)) {
        ret2 = digammaDiff((deltaDiffA + h + l - n) / 2, n) * RIII(l, n);
    } else {
        int x = (int)floor((deltaDiffA + h + l - n) / 2 + EPS);
        if (x > 0 || x + n <= 0) {
            ret2 = digammaDiff((deltaDiffA + h + l - n) / 2, n) * RIII(l, n);
        } else {                
            // -n + 1 <= x <= 0 need special handling because digammaDiff() has a pole and RII() is zero.
            // now, RII() ~ (x)_n, the (x)_n factor has a pole
            ret2 = rising_factorial((deltaDiffA - h - l - n + 2) / 2, n);
            ret2 *= rising_factorial((deltaDiffB + h + l - n) / 2, n) * rising_factorial((deltaDiffB - h - l - n + 2) / 2, n);
            ret2 *= factorial<float_type>(-x) * factorial<float_type>(x + n - 1);
            if (x % 2 != 0) ret2 = -ret2;

            ret2 *= RIII_NonDelta(l, n);
        }
    }

    return (ret1 + ret2) * 0.5;
}

float_type ConformalBlockScalars::HInfinity(int l, float_type r, float_type eta) {
    float_type ret;
    if (hCache.Get(l, &ret)) return ret;

    float_type h = d / 2;
    ret = pow(1 - r * r, 1 - h);
    ret /= pow(r * r - 2 * r * eta + 1, (1 - dlt12 + dlt34) / 2);
    ret /= pow(r * r + 2 * r * eta + 1, (1 + dlt12 - dlt34) / 2);
    ret *= factorial<float_type>(l);
    ret /= (((i64)1) << l);

    if (abs(h - 1) > EPS) {
        ret /= rising_factorial(h - 1, l);
        float_type x[1];
        x[0] = eta;
        float_type *v = gegenbauer_polynomial_value(l, 1, h - 1, x);
        ret *= v[l + 0 * (l + 1)];
        delete v;
    } else if (l != 0) {
        // for the h == 1 and l == 0 case, (h-1)_l == 1 == gegenbauer(l, h-1);
        // for the h == 1 and l > 0 case, the zeros of (h-1)_l and gegenbauer(l, h-1) cancel out
        ret *= GegenbauerDAlphaAt0(l, eta) / factorial<float_type>(l - 1);
    }

    hCache.Set(l, ret);
    return ret;
}

float_type ConformalBlockScalars::HInfinity_dDelta12(int l, float_type r, float_type eta) {
    return HInfinity_dDelta(l, r, eta);
}

float_type ConformalBlockScalars::HInfinity_dDelta34(int l, float_type r, float_type eta) {
    return -HInfinity_dDelta(l, r, eta);
}

float_type ConformalBlockScalars::HInfinity_dDelta(int l, float_type r, float_type eta) {
    return 0.5 * HInfinity(l, r, eta) * log((r * r - 2 * eta * r + 1)/(r * r + 2 * eta * r + 1));
}

float_type ConformalBlockScalars::HRecursion(float_type delta, int l, float_type r, float_type eta, int order) {
    int idelta = std::numeric_limits<int>::max();
    float_type ret;
    if (IsInteger(delta)) {
        idelta = (int)floor(delta + EPS);
        if (hCache.Get(idelta, l, order, &ret)) return ret;
    }

    ret = HInfinity(l, r, eta);
    if (order == 0) return ret;

    float_type deltaAs;
    int la, na;

    // type I
    for (int n = 1; n <= order; n++) {
        deltaAs = 1 - l - n;
        la = l + n;
        na = n;
        // this if statement is only for the following case: 
        // when the exchange operator is identity, i.e, delta==l==0,
        // then delta - deltaAs can be zero.
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) * rCoef->RI(n) / (delta - deltaAs) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    // type II
    for (int n = 1; n <= order && n <= l; n++) {
        deltaAs = l + d - 1 - n;
        la = l - n;
        na = n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) * rCoef->RII(l, n) / (delta - deltaAs) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    // type III
    for (int n = 1; n + n <= order; n++) {
        deltaAs = d / 2 - n;
        la = l;
        na = 2 * n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) * rCoef->RIII(l, n) / (delta - deltaAs) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    if (idelta != std::numeric_limits<int>::max()) {
        hCache.Set(idelta, l, order, ret);
    }

    return ret;
}

#if CONFORMAL_BLOCK_NORMALIZATION == 1

float_type ConformalBlockScalars::NormalizationFactor(int l) {
    assert(l >= 0);
    assert(l < 63);

    float_type ret = (float_type)(((i64)1)<<l);
    if (l % 2 == 1) ret = -ret;
    if (abs(d - 2) < EPS) {
        return ret;
    } else {
        return ret * rising_factorial(d/2 - 1.0, l) / rising_factorial(d - 2.0, l);
    }
}

#endif

float_type ConformalBlockScalars::evaluate(float_type delta, int l, float_type r, float_type eta, int order) {
    hCache.ClearIfNew(r, eta);
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    float_type prefactor = pow(r, delta) * NormalizationFactor(l);
    float_type hr = HRecursion(delta, l, r, eta, order);
    return  prefactor * hr;
#else
    return pow(4 * r, delta) * HRecursion(delta, l, r, eta, order);
#endif

}

float_type ConformalBlockScalars::dDelta(float_type delta, int l, float_type r, float_type eta, int order) {
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    float_type ret = log(r) * evaluate(delta, l, r, eta, order);
    ret += pow(r, delta) * NormalizationFactor(l) * dHdDelta(delta, l, r, eta, order);
#else
    float_type ret = log(4 * r) * evaluate(delta, l, r, eta, order);
    ret += pow(4 * r, delta) * dHdDelta(delta, l, r, eta, order);
#endif
    return ret;
}

float_type ConformalBlockScalars::dHdDelta(float_type delta, int l, float_type r, float_type eta, int order) {
    float_type deltaAs, ret = 0.0;
    int la, na;

    // type I
    for (int n = 1; n <= order; n++) {
        deltaAs = 1 - l - n;
        la = l + n;
        na = n;
        // this if statement is only for the following case: 
        // when the exchange operator is identity, i.e, delta==l==0,
        // then delta - deltaAs can be zero.
        if (abs(delta - deltaAs) < EPS) continue;
        ret -= pow(4 * r, na) * rCoef->RI(n) / ((delta - deltaAs) * (delta - deltaAs)) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    // type II
    for (int n = 1; n <= order && n <= l; n++) {
        deltaAs = l + d - 1 - n;
        la = l - n;
        na = n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret -= pow(4 * r, na) * rCoef->RII(l, n) / ((delta - deltaAs) * (delta - deltaAs)) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    // type III
    for (int n = 1; n + n <= order; n++) {
        deltaAs = d / 2 - n;
        la = l;
        na = 2 * n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret -= pow(4 * r, na) * rCoef->RIII(l, n) / ((delta - deltaAs) * (delta - deltaAs)) * HRecursion(deltaAs + na, la, r, eta, order - na);
    }

    return ret;
}

float_type ConformalBlockScalars::dDelta12(float_type delta, int l, float_type r, float_type eta, int order) {
    hDiff12Cache.ClearIfNew(r, eta);
    hCache.ClearIfNew(r, eta);
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    return pow(r, delta) * NormalizationFactor(l) * dHdDelta12(delta, l, r, eta, order);
#else
    return pow(4 * r, delta) * dHdDelta12(delta, l, r, eta, order);
#endif

}

float_type ConformalBlockScalars::dHdDelta12(float_type delta, int l, float_type r, float_type eta, int order) {
    int idelta = std::numeric_limits<int>::max();
    float_type ret;
    if (IsInteger(delta)) {
        idelta = (int)floor(delta + EPS);
        if (hDiff12Cache.Get(idelta, l, order, &ret)) return ret;
    }

    ret = HInfinity_dDelta12(l, r, eta);
    if (order == 0) return ret;

    float_type deltaAs;
    int la, na;

    // type I
    for (int n = 1; n <= order; n++) {
        deltaAs = 1 - l - n;
        la = l + n;
        na = n;
        // this if statement is only for the following case: 
        // when the exchange operator is identity, i.e, delta==l==0,
        // then delta - deltaAs can be zero.
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta12(deltaAs + na, la, r, eta, order - na) * rCoef->RI(n) + HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RI_dDelta12(n));
    }

    // type II
    for (int n = 1; n <= order && n <= l; n++) {
        deltaAs = l + d - 1 - n;
        la = l - n;
        na = n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta12(deltaAs + na, la, r, eta, order - na) * rCoef->RII(l, n)+ HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RII_dDelta12(l,n));
    }

    // type III
    for (int n = 1; n + n <= order; n++) {
        deltaAs = d / 2 - n;
        la = l;
        na = 2 * n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta12(deltaAs + na, la, r, eta, order - na) * rCoef->RIII(l, n)+ HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RIII_dDelta12(l,n));
    }

    if (idelta != std::numeric_limits<int>::max()) {
        hDiff12Cache.Set(idelta, l, order, ret);
    }


    return ret;
}

float_type ConformalBlockScalars::dDelta34(float_type delta, int l, float_type r, float_type eta, int order) {
    hDiff34Cache.ClearIfNew(r, eta);
    hCache.ClearIfNew(r, eta);
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    return pow(r, delta) * NormalizationFactor(l) * dHdDelta34(delta, l, r, eta, order);
#else
    return pow(4 * r, delta) * dHdDelta34(delta, l, r, eta, order);
#endif
}

float_type ConformalBlockScalars::dHdDelta34(float_type delta, int l, float_type r, float_type eta, int order) {
    int idelta = std::numeric_limits<int>::max();
    float_type ret;
    if (IsInteger(delta)) {
        idelta = (int)floor(delta + EPS);
        if (hDiff34Cache.Get(idelta, l, order, &ret)) return ret;
    }

    ret = HInfinity_dDelta34(l, r, eta);
    if (order == 0) return ret;
    float_type deltaAs;
    int la, na;

    // type I
    for (int n = 1; n <= order; n++) {
        deltaAs = 1 - l - n;
        la = l + n;
        na = n;
        // this if statement is only for the following case: 
        // when the exchange operator is identity, i.e, delta==l==0,
        // then delta - deltaAs can be zero.
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta34(deltaAs + na, la, r, eta, order - na) * rCoef->RI(n) + HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RI_dDelta34(n));
    }

    // type II
    for (int n = 1; n <= order && n <= l; n++) {
        deltaAs = l + d - 1 - n;
        la = l - n;
        na = n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta34(deltaAs + na, la, r, eta, order - na) * rCoef->RII(l, n)+ HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RII_dDelta34(l,n));
    }

    // type III
    for (int n = 1; n + n <= order; n++) {
        deltaAs = d / 2 - n;
        la = l;
        na = 2 * n;
        if (abs(delta - deltaAs) < EPS) continue;
        ret += pow(4 * r, na) / (delta - deltaAs) * (dHdDelta34(deltaAs + na, la, r, eta, order - na) * rCoef->RIII(l, n)+ HRecursion(deltaAs + na, la, r, eta, order - na) * rCoef->RIII_dDelta34(l,n));
    }

    if (idelta != std::numeric_limits<int>::max()) {
        hDiff34Cache.Set(idelta, l, order, ret);
    }

    return ret;
}


