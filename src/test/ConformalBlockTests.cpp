//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>

#include <sstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/timer.hpp>
#include "test.h"
#include "../ConformalBlock.h"

struct Scalar4ptTestData
{
public:
    // spacetime dimension.
    float_type d;

    // the value of Delta1 - Delta 2.
    float_type dlt12;
    
    // the value of Delta3 - Delta 4.
    float_type dlt34;
};

std::ostream& operator<< (std::ostream& out, const Scalar4ptTestData& data) {
  out << '{';
  out << "d=" << data.d << ", dlt12=" << data.dlt12 << ", dlt34=" << data.dlt34;
  out << '}';
  return out;
}

vector<Scalar4ptTestData> CreateScalar4ptTestData(int size, bool evenD = false, bool intDelta = false) {
    vector<Scalar4ptTestData> ret(size);
    for (int i = 0; i < size; i++) {
        if (evenD) {
                ret[i].d = (float_type)(2 * randomint(1, 4));
            } else {
                ret[i].d = (float_type)(2 * randomint(1, 4) + 1);
            }
        if (intDelta) {
                ret[i].dlt12 = -randomint(0, 4);
                ret[i].dlt34 = -randomint(0, 4);
            } else {
                ret[i].dlt12 = -random(0.0, 4.0);
                ret[i].dlt34 = -random(0.0, 4.0);
            }
    }
    
    return ret;
}

struct ScalarCBTestData
{
public:
    // scaling dimension of the interchange operator.
    double delta;

    // spin of the interchange opertor.
    int l;

    // radial coordinate r.
    double r;

    // eta = cos[theta]
    double eta;

    // order of r to compute
    int order;
};

std::ostream& operator<< (std::ostream& out, const ScalarCBTestData& data) {
  out << '{';
  out << "delta=" << data.delta << ", l=" << data.l << ", r=" << data.r << ", eta=" << data.eta << ", order=" << data.order;
  out << '}';
  return out;
}

vector<ScalarCBTestData> CreateScalarCBTestData(int size) {
    vector<ScalarCBTestData> ret(size);
    for (int i = 0; i < size; i++) {
        ret[i].l = randomint(0, 10);

        // TODO: apply unitary bound?
        ret[i].delta = random(0.0, 4.0) + ret[i].l - 2;
        ret[i].r = random(0.0, 0.2);
        ret[i].eta = random(0.0, 1.0);
        ret[i].order = randomint(10, 16);
    }
    
    return ret;
}

// test suite1: odd d and floating delta12 and delta34
BOOST_FIXTURE_TEST_SUITE(ScalarCB_suite1, SimpleTestFixture, * utf::label("ScalarCB"))

BOOST_DATA_TEST_CASE(RI_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    float_type actual = r.RI_dDelta12(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RI_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    float_type actual = r.RI_dDelta34(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RII_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta12(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RII_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta34(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta12(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta34(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(dDelta_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    float_type actual = cb.dDelta(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb.evaluate(sc.delta + inc, sc.l, sc.r, sc.eta, sc.order) - cb.evaluate(sc.delta - inc, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(dDelta12_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12 - inc, rc.dlt34);
    ConformalBlockScalars cb2(rc.d, rc.dlt12 + inc, rc.dlt34);
    float_type actual = cb.dDelta12(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(dDelta34_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12, rc.dlt34 - inc);
    ConformalBlockScalars cb2(rc.d, rc.dlt12, rc.dlt34 + inc);
    float_type actual = cb.dDelta34(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_AUTO_TEST_SUITE_END()

// test suite2: odd d and integer delta12 and delta34
BOOST_FIXTURE_TEST_SUITE(ScalarCB_suite2, SimpleTestFixture, * utf::label("ScalarCB"))

BOOST_DATA_TEST_CASE(RI_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    float_type actual = r.RI_dDelta12(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RI_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    float_type actual = r.RI_dDelta34(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RII_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta12(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RII_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta34(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta12_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta12(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta34_test, CreateScalar4ptTestData(TCNumber, false, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta34(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(dDelta_test,
    CreateScalar4ptTestData(TCNumber, false, true) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    float_type actual = cb.dDelta(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb.evaluate(sc.delta + inc, sc.l, sc.r, sc.eta, sc.order) - cb.evaluate(sc.delta - inc, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}


BOOST_DATA_TEST_CASE(dDelta12_test,
    CreateScalar4ptTestData(TCNumber, false, true) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12 - inc, rc.dlt34);
    ConformalBlockScalars cb2(rc.d, rc.dlt12 + inc, rc.dlt34);
    float_type actual = cb.dDelta12(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(dDelta34_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12, rc.dlt34 - inc);
    ConformalBlockScalars cb2(rc.d, rc.dlt12, rc.dlt34 + inc);
    float_type actual = cb.dDelta34(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_AUTO_TEST_SUITE_END()


// test suite3: even d and floating delta12 and delta34
BOOST_FIXTURE_TEST_SUITE(ScalarCB_suite3, SimpleTestFixture, * utf::label("ScalarCB"))

BOOST_DATA_TEST_CASE(RI_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    float_type actual = r.RI_dDelta12(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RI_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    float_type actual = r.RI_dDelta34(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RII_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta12(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RII_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta34(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta12(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, false) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta34(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(dDelta_test,
    CreateScalar4ptTestData(TCNumber, true, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    float_type actual = cb.dDelta(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb.evaluate(sc.delta + inc, sc.l, sc.r, sc.eta, sc.order) - cb.evaluate(sc.delta - inc, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}


BOOST_DATA_TEST_CASE(dDelta12_test,
    CreateScalar4ptTestData(TCNumber, true, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12 - inc, rc.dlt34);
    ConformalBlockScalars cb2(rc.d, rc.dlt12 + inc, rc.dlt34);
    float_type actual = cb.dDelta12(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(dDelta34_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12, rc.dlt34 - inc);
    ConformalBlockScalars cb2(rc.d, rc.dlt12, rc.dlt34 + inc);
    float_type actual = cb.dDelta34(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_AUTO_TEST_SUITE_END()


// test suite4: even d and integer delta12 and delta34
BOOST_FIXTURE_TEST_SUITE(ScalarCB_suite4, SimpleTestFixture, * utf::label("ScalarCB"))

BOOST_DATA_TEST_CASE(RI_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    float_type actual = r.RI_dDelta12(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RI_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, n)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    float_type actual = r.RI_dDelta34(n);
    float_type expect = (r2.RI(n) - r1.RI(n)) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(RII_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta12(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RII_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= l; n++) {
        float_type actual = r.RII_dDelta34(l, n);
        float_type expect = (r2.RII(l, n) - r1.RII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta12_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12 - inc, rc.dlt34);
    RCoefficients r2(rc.d, rc.dlt12 + inc, rc.dlt34);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta12(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(RIII_dDelta34_test, CreateScalar4ptTestData(TCNumber, true, true) ^ bdata::xrange(1, 1 + TCNumber), rc, l)
{
    RCoefficients r(rc.d, rc.dlt12, rc.dlt34);
    RCoefficients r1(rc.d, rc.dlt12, rc.dlt34 - inc);
    RCoefficients r2(rc.d, rc.dlt12, rc.dlt34 + inc);

    for (int n = 1; n <= TCNumber; n++) {
        float_type actual = r.RIII_dDelta34(l, n);
        float_type expect = (r2.RIII(l, n) - r1.RIII(l, n)) / (2 * inc);
        MY_FLOAT_EQUAL(actual, expect, tol);
    }
}

BOOST_DATA_TEST_CASE(dDelta_test,
    CreateScalar4ptTestData(TCNumber, true, true) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    float_type actual = cb.dDelta(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb.evaluate(sc.delta + inc, sc.l, sc.r, sc.eta, sc.order) - cb.evaluate(sc.delta - inc, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}


BOOST_DATA_TEST_CASE(dDelta12_test,
    CreateScalar4ptTestData(TCNumber, true, true) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12 - inc, rc.dlt34);
    ConformalBlockScalars cb2(rc.d, rc.dlt12 + inc, rc.dlt34);
    float_type actual = cb.dDelta12(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(dDelta34_test,
    CreateScalar4ptTestData(TCNumber, false, false) ^ CreateScalarCBTestData(TCNumber), rc, sc)
{
    ConformalBlockScalars cb(rc.d, rc.dlt12, rc.dlt34);
    ConformalBlockScalars cb1(rc.d, rc.dlt12, rc.dlt34 - inc);
    ConformalBlockScalars cb2(rc.d, rc.dlt12, rc.dlt34 + inc);
    float_type actual = cb.dDelta34(sc.delta, sc.l, sc.r, sc.eta, sc.order);
    float_type expect = (cb2.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order) - cb1.evaluate(sc.delta, sc.l, sc.r, sc.eta, sc.order)) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(conformalblock_suite)

// test the accuracy of the functions. The test data comes from a failure of bootstrap run.
BOOST_AUTO_TEST_CASE(Accuracy_test, * utf::label("ScalarCB.accuracy"))
{
    vector<float_type> dims = vector<float_type>({0.509292, 1.37473, 2.23385, 3.20297});
    float_type dlt12 = dims[0] - dims[3];
    float_type dlt34 = dims[0] - dims[3];
    float_type delta = dims[0];
    cpx_t rho = Z2rho(cpx_t(1-0.41671, -0.0386725));
    float_type r = abs(rho);
    float_type eta = cos(arg(rho));
    
    ConformalBlockScalars cb(3, dlt12, dlt34);
    float_type value1 = cb.evaluate(delta, 0, r, eta, 12);
    float_type value2 = cb.evaluate(delta, 0, r, eta, 13);
    MY_FLOAT_EQUAL(value1, value2, 1e-6);

    value1 = cb.dDelta(delta, 0, r, eta, 12);
    value2 = cb.dDelta(delta, 0, r, eta, 13);
    MY_FLOAT_EQUAL(value1, value2, 1e-6);

    value1 = cb.dDelta12(delta, 0, r, eta, 12);
    value2 = cb.dDelta12(delta, 0, r, eta, 13);
    MY_FLOAT_EQUAL(value1, value2, 1e-6);

    value1 = cb.dDelta34(delta, 0, r, eta, 12);
    value2 = cb.dDelta34(delta, 0, r, eta, 13);
    MY_FLOAT_EQUAL(value1, value2, 1e-6);
}

BOOST_AUTO_TEST_CASE(ConformalBlock_Performance_test, * utf::label("ScalarCB.perf"))
{
    vector<float_type> dims = vector<float_type>({0.509292, 1.37473, 2.23385, 3.20297});
    float_type dlt12 = dims[0] - dims[3];
    float_type dlt34 = dims[0] - dims[3];
    float_type delta = dims[0];
    cpx_t rho1 = Z2rho(cpx_t(1-0.41671, -0.0386725));
    float_type r1 = abs(rho1);
    float_type eta1 = cos(arg(rho1));

    cpx_t rho2 = Z2rho(cpx_t(0.598759,0.00609649));
    float_type r2 = abs(rho2);
    float_type eta2 = cos(arg(rho2));

    ConformalBlockScalars cb(3, dlt12, dlt34);
    int count = 100;

    boost::timer stopwatch;
    for (int i = 0; i < count; i++) {
        cb.evaluate(delta, 0, r1, eta1, 12);
        cb.evaluate(delta, 0, r2, eta2, 12);
    }

    std::cout << "Performance: calling ScalarCB.evaluate " << 2 * count << " times takes " << stopwatch.elapsed() << " seconds." << std::endl;

    stopwatch.restart();

    for (int i = 0; i < count; i++) {
        cb.dDelta(delta, 0, r1, eta1, 12);
        cb.dDelta(delta, 0, r2, eta2, 12);
    }

    std::cout << "Performance: calling ScalarCB.dDelta " << 2 * count << " times takes " << stopwatch.elapsed() << " seconds." << std::endl;

    stopwatch.restart();

    for (int i = 0; i < count; i++) {
        cb.dDelta12(delta, 0, r1, eta1, 12);
        cb.dDelta12(delta, 0, r2, eta2, 12);
    }

    std::cout << "Performance: calling ScalarCB.dDelta12 " << 2 * count << " times takes " << stopwatch.elapsed() << " seconds." << std::endl;

    stopwatch.restart();
}

//Name your test cases for what they test
BOOST_AUTO_TEST_CASE(evaluate_test, * boost::unit_test::tolerance(boost::test_tools::fpc::percent_tolerance(0.001)))
{
    int d = 4;
    int l = 6;
    int order = 7;
    float_type dlt12 = 0.0;
    float_type dlt34 = 0.0;
    float_type r = 0.124;
    float_type Pi = acos(-1.0);
    float_type eta = cos(Pi/3);
    float_type dlt = 8.3;
    ConformalBlockScalars cb(d, dlt12, dlt34);
    float_type res = cb.evaluate(dlt, l, r, eta, order);
    // the expacted value obtained from mathematica program
#if CONFORMAL_BLOCK_NORMALIZATION == 1
    float_type exp = 4.327452890496684 * 1e-9;
#else
    float_type exp = 0.0000470163;
#endif

    BOOST_TEST(res == exp);
}

BOOST_AUTO_TEST_SUITE_END()
