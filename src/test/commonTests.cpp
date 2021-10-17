//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>

#include <cmath>
#include "test.h"
#include "../common.h"
#include "../gegenbauer_polynomial.hpp"

BOOST_DATA_TEST_CASE(GegenbauerDAlphaAt0_test, bdata::xrange(0, 20) ^ bdata::random(-10.0, 10.0), n, x)
{
    float_type actual = GegenbauerDAlphaAt0(n, x);

    float_type arr[1] = {x};
    float_type *v1 = gegenbauer_polynomial_value(n, 1, -inc, arr);
    float_type *v2 = gegenbauer_polynomial_value(n, 1, inc, arr);
    float_type ret1 = v1[n + 0 * (n + 1)];
    float_type ret2 = v2[n + 0 * (n + 1)];
    delete v1;
    delete v2;
    
    float_type expect = (ret2 - ret1) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

