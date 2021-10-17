//Link to Boost
#define BOOST_TEST_DYN_LINK

//Define our Module name (prints at testing)
 #define BOOST_TEST_MODULE Conformal Bootstrap UnitTests

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cmath>
#include "test.h"
#include "../ConformalBlock.h"
#include "../BootstrapRunner.h"

BOOST_AUTO_TEST_SUITE(demo_suite)

//Name your test cases for what they test
BOOST_AUTO_TEST_CASE(demo1, * utf::label("demo"))
{
    int d = 4;
    int l = 6;
    int order = 7;
    double dlt12 = 0.6;
    double dlt34 = -0.3;
    double r = 0.124;
    double Pi = acos(-1.0);
    double eta = cos(Pi/3);
    double dlt = 8.3;
    ConformalBlockScalars cb(d, dlt12, dlt34);
    double res1 = cb.evaluate(dlt, l, r, eta, order);
    std::cout << "res1=" << res1 << std::endl;

    ConformalBlockScalars cb2(2, 0.01, 0.0);
    float_type actual = cb2.evaluate(0.0, 0, r, eta, 20);
    std::cout << "actual=" << actual << std::endl;

    std::cout << cb2.dDelta(0.0, 0, r, eta, 15) << std::endl;
    std::cout << cb2.dDelta12(0.0, 0, r, eta, 15) << std::endl;
    std::cout << cb2.dDelta34(0.0, 0, r, eta, 15) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
