#pragma once
#include <boost/unordered_map.hpp>
#include <assert.h>
#include <limits.h>
#include "common.h"
using namespace std;

// if CONFORMAL_BLOCK_NORMALIZATION is 1, we normalize CB 
// as [Kos, Poland, Simmons-Duffin] arxiv:1307.6865 and arxiv:1406.4858;
// otherwise, we normalize CB as [Penedones, Trevisani, Masahito, Yamazaki] arxiv:1509.00428
#define CONFORMAL_BLOCK_NORMALIZATION 1

// the precision of computation.
#define PRECISION EPS

// Holder class for the R_A coefficients.
class RCoefficients
{
private:
    static const int MaxOrder = 32;

    boost::unordered_map<int, float_type> cacheR;

    float_type* cacheRI;

    // dimension of the conformal block.
    float_type d;

    // the value of Delta1 - Delta 2.
    float_type dlt12;
    
    // the value of Delta3 - Delta 4.
    float_type dlt34;

    float_type& GetCachedRI(int type, int n);

    // the delta independent part of RI.
    float_type RI_NonDelta(int n);

    // function used by RI_dDelta12 and RI_dDelta34
    float_type RI_dDelta(int n, float_type deltaDiffA, float_type deltaDiffB);

    // the delta independent part of RII.
    float_type RII_NonDelta(int l, int n);

    // function used by RII_dDelta12 and RII_dDelta34
    float_type RII_dDelta(int l, int n, float_type deltaDiffA, float_type deltaDiffB);

    // the delta independent part of RIII.
    float_type RIII_NonDelta(int l, int n);

    // function used by RIII_dDelta12 and RIII_dDelta34
    float_type RIII_dDelta(int l, int n, float_type deltaDiffA, float_type deltaDiffB);

    // Check cache for RII and RIII,
    // type = 0, 1, ... 5. 
    float_type& GetCachedRII(int type, int l, int n);
public:
    RCoefficients(float_type d, float_type dlt12, float_type dlt34) {
        this->d = d;
        this->dlt12 = dlt12;
        this->dlt34 = dlt34;

        this->cacheRI = new float_type[3 * MaxOrder];
        for (int i = 0; i < 3 * MaxOrder; i++) {
            this->cacheRI[i] = std::numeric_limits<float_type>::max();
        }
    }

    ~RCoefficients() {
        delete[] cacheRI;
        cacheRI = 0;
    }

    // the R_{I,n} coefficients of the recursion relation.
    float_type RI(int n);

    // derivative of R_{I,n} with respect the delta12
    float_type RI_dDelta12(int n);

    // derivative of R_{I,n} with respect the delta34
    float_type RI_dDelta34(int n);

    // the R_{II,n} coefficients of the recursion relation.
    float_type RII(int l, int n);

    // derivative of R_{II,n} with respect the delta12
    float_type RII_dDelta12(int l, int n);

    // derivative of R_{II,n} with respect the delta34
    float_type RII_dDelta34(int l, int n);

    // the R_{III,n} coefficients of the recursion relation.
    float_type RIII(int l, int n);

    // derivative of R_{III,n} with respect the delta12
    float_type RIII_dDelta12(int l, int n);

    // derivative of R_{III,n} with respect the delta34
    float_type RIII_dDelta34(int l, int n);
};

// Provides cache functionalify for conformal blocks calculation.
class ComformalBlockCache
{
private:
    boost::unordered_map<int, float_type> cache;
    float_type r;
    float_type eta;
    static int ToKey(int delta, int l, int order)
    {
        return (delta * 64 * 64) + order * 64 + l;
    }

    int hits;

public:
    ComformalBlockCache() { this-> r = -1.0; this->eta = 0; this->hits = 0; };
    ~ComformalBlockCache() {};

    // clear the cache if the (r, eta) position is new.
    void ClearIfNew(float_type r, float_type eta)
    {
        if (abs(r - this->r) > EPS || abs(eta - this->eta) > EPS) {
            cache.clear();
            this->r = r;
            this->eta = eta;
            this->hits = 0;
        }
    }

    // return true if cache is hit, and false otherwise.
    // when the return value is true, the cached value is returned via the ret argument.
    bool Get(int delta, int l, int order, float_type *ret)
    {
        assert(l < 64);
        boost::unordered_map<int, float_type>::iterator it = cache.find(ToKey(delta, l, order));
        if (it == cache.end()) return false;
        hits++;
        *ret = it->second;
        return true;
    }

    int Size() { return cache.size(); }

    int Hits() { return hits; }

    // this function is only for HInfinity() and HInfinity_dDelta(), which only need on parameter l as key.
    // return true if cache is hit, and false otherwise.
    // when the return value is true, the cached value is returned via the ret argument.
    bool Get(int l, float_type *ret)
    {
        assert(l < 64);
        boost::unordered_map<int, float_type>::iterator it = cache.find(LONG_MAX - l);
        if (it == cache.end()) return false;

        *ret = it->second;
        return true;
    }

    // update cache.
    void Set(int delta, int l, int order, float_type value)
    {
        assert(l < 64);
        cache[ToKey(delta, l, order)] = value;
    }

    // update cache. for HInfinity() and HInfinity_dDelta() only.
    void Set(int l, float_type value)
    {
        assert(l < 64);
        cache[LONG_MAX - l] = value;
    }
};

// calculate the conformal block of 4-point scalars functions.
// Implements the recursion relation of the paper arxiv:1509.00428
class ConformalBlockScalars
{
private:
    // default order for calculation.
    static int const DefaultOrder = 12;

    // dimension of the conformal block.
    float_type d;

    // the value of Delta1 - Delta 2.
    float_type dlt12;
    
    // the value of Delta3 - Delta 4.
    float_type dlt34;

    // the object holds R_{A} coefficients.
    RCoefficients *rCoef;

    // caches the intermediate results.
    // cache.Get(n, l, ord, &res) returns res = HRecursion(n, l, r, eta, order[>=ord]).
    ComformalBlockCache hCache;

    // caches the intermediate results of dh/dDelta12
    ComformalBlockCache hDiff12Cache;

    // caches the intermediate results of dh/dDelta34
    ComformalBlockCache hDiff34Cache;

    // the term h_{\infinity, l} of the recursion relation.
    float_type HInfinity(int l, float_type r, float_type eta);

    // derivative the h_{\infinity, l} with respect to delta12
    float_type HInfinity_dDelta12(int l, float_type r, float_type eta);

    // derivative the h_{\infinity, l} with respect to delta34
    float_type HInfinity_dDelta34(int l, float_type r, float_type eta);

    // function used by HInfinity_dDelta12 and HInfinity_dDelta34
    float_type HInfinity_dDelta(int l, float_type r, float_type eta);

    // the recursion function h_{\Delta, l}
    float_type HRecursion(float_type delta, int l, float_type r, float_type eta, int order);

    // calculate dH/dDelta
    float_type dHdDelta(float_type delta, int l, float_type r, float_type eta, int order);

    // calculate dH/dDelta12
    float_type dHdDelta12(float_type delta, int l, float_type r, float_type eta, int order);

    // calculate dH/dDelta34
    float_type dHdDelta34(float_type delta, int l, float_type r, float_type eta, int order);
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    float_type NormalizationFactor(int l);
#endif
public:
    // constructor 
    ConformalBlockScalars(float_type d, float_type dlt12 = .0, float_type dlt34 = .0) {
        this->d = d;
        this->dlt12 = dlt12;
        this->dlt34 = dlt34;
        rCoef = new RCoefficients(this->d, this->dlt12, this->dlt34);
    }

    // constructor with integer dimension.
    ConformalBlockScalars(int d, float_type dlt12 = .0, float_type dlt34 = .0) {
        this->d = (float_type)d;
        this->dlt12 = dlt12;
        this->dlt34 = dlt34;

        rCoef = new RCoefficients(this->d, this->dlt12, this->dlt34);
    }

    ~ConformalBlockScalars() {
        delete rCoef;
    }

//    float_type evaluate(float_type delta, int l, float_type r, float_type eta);

    float_type evaluate(float_type delta, int l, float_type r, float_type eta, int order = DefaultOrder);

    float_type dDelta(float_type delta, int l, float_type r, float_type eta, int order = DefaultOrder);

    float_type dDelta12(float_type delta, int l, float_type r, float_type eta, int order = DefaultOrder);
    
    float_type dDelta34(float_type delta, int l, float_type r, float_type eta, int order = DefaultOrder);
};
