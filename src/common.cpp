#include "common.h"
#include <cmath>
#include <ctime>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

boost::random::mt19937 rng(time(0));
boost::random::uniform_01<float_type> dist01;

// calculate digamma(x+n) - digamma(x) = \sum 1/(x+k), k=0,1,...n-1
float_type digammaDiff(float_type x, int n)
{
    float_type ret = 0.0;

    for (int i = n - 1; i >= 0; i--) {
        ret += 1.0/(x + i);
    }

    return ret;
}

float_type GegenbauerDAlphaAt0(int n, float_type x)
{
    if (n <= 0) {
        return 0.0;
    } else if (n == 1) {
        return 2 * x;
    } else if (n == 2) {
        return 2 * x * x - 1;
    }

    // for -1 < x < 1 it actually equals 
    // return 2 * cos(n * acos(x)) / n;

    float_type a = 2 * x, b = 2 * x * x - 1;
    float_type c;
    for (int i = 3; i <= n; i++) {
        c = (2 * (i - 1) * x * b - (i - 2) * a) / i;
        a = b; 
        b = c;
    }

    return c;
}

bool IsInteger(float_type n)
{
    return std::abs(n - floor(n + EPS)) < EPS;
}

float_type random(float_type min, float_type max)
{
    return min + (max - min) * dist01(rng);
}

int randomint(int min, int max)
{
    boost::random::uniform_int_distribution<int> dist(min,max);
    return dist(rng);
}

cpx_t RandomComplex(float_type maxNorm)
{
    float_type r = random(0.0, maxNorm);
    float_type theta = random(0.0, 2 * acos(-1.0));
    return cpx_t(r * cos(theta), r * sin(theta));
}

std::string NowToString()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%Y-%m-%d %I:%M:%S",timeinfo);
    return std::string(buffer);
}
