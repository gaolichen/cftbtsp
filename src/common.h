#pragma once
#include <iostream>
#include <iterator>
#include <sstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <complex>
#include <string>
#include <sstream>
#include <ctime>
#include <limits>

#define EPS 1e-10

#if WIN32
typedef __int64 i64;
#else
typedef long long i64;
#endif

typedef unsigned int uint;
typedef double float_type;
typedef std::complex<float_type> cpx_t;

extern boost::random::mt19937 rng;
extern boost::random::uniform_01<float_type> dist01;

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
  out << '[';
  if (!v.empty()) {
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b";
  }
  //out << "\b\b]";
  out << ']';
  return out;
}

template<class T> std::string ToString(T a)
{
	std::ostringstream oss;
	oss << a;
	return oss.str();
};

template <typename T>
T Sum(const std::vector<T>& v)
{
    T ret = (T)0;
    for (int i = 0; i < v.size(); i++) ret += v[i];
    return ret;
}

template<typename T>
T Z2rho(T z)
{
    return z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
}

template<typename T>
float_type Z2r(T z)
{
    T rho = z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
    return (float_type)abs(rho);
}

// return datetime as string
std::string NowToString();

// calculate digamma(z+n) - digamma(z)
float_type digammaDiff(float_type z, int n);

// returns derivative of C(n, alpha; x) with respect to alpha at alpha = 0, 
// where C(n, alpha; x) is the Gegenbauer polynomials.
float_type GegenbauerDAlphaAt0(int n, float_type x);

bool IsInteger(float_type n);

float_type random(float_type min, float_type max);

int randomint(int min, int max);

cpx_t RandomComplex(float_type maxNorm);


