//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>

#include <vector>
#include <iostream>
#include "test.h"
#include "../BootstrapRunner.h"
using namespace std;

BOOST_FIXTURE_TEST_SUITE(CorrelatorSet_suite, SimpleTestFixture, * utf::label("correlatorset"))

BOOST_DATA_TEST_CASE(CostDerDelta_test1, CreateCfgConfigTestData(2), cfg)
{
    CftData cftData(cfg);
    int scalarRange = min(4, cftData.MaxScalarId());

    CorrelatorSet correlatorSet(&cftData, scalarRange);
    cpx_t z = RandomComplex(0.3) + 0.5;
    int op = randomint(1, scalarRange);

    float_type actual = correlatorSet.CostDerDelta(z, op);

    float_type dim = cftData.GetPrimaryDim(op);
    cftData.SetPrimaryDim(op, dim - inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value1 = correlatorSet.Cost(z);

    cftData.SetPrimaryDim(op, dim + inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value2 = correlatorSet.Cost(z);

    float_type expect = (value2 - value1) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CostDerDelta_test2, CreateCfgConfigTestData(4), cfg)
{
    CftData cftData(cfg);

    CorrelatorSet correlatorSet(&cftData, min(4, cftData.MaxScalarId()));
    cpx_t z = RandomComplex(0.3) + 0.5;

    int op = randomint(cftData.StressTensorId + 1, cftData.MaxPrimaryId());

    float_type actual = correlatorSet.CostDerDelta(z, op);

    float_type dim = cftData.GetPrimaryDim(op);
    cftData.SetPrimaryDim(op, dim - inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value1 = correlatorSet.Cost(z);

    cftData.SetPrimaryDim(op, dim + inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value2 = correlatorSet.Cost(z);

    float_type expect = (value2 - value1) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CostDerOpeCoefficient_test, CreateCfgConfigTestData(4), cfg)
{
    CftData cftData(cfg);

    int scalarRange = min(4, cftData.MaxScalarId());
    CorrelatorSet correlatorSet(&cftData, scalarRange);
    cpx_t z = RandomComplex(0.3) + 0.5;

    int op1 = randomint(1, scalarRange);
    int op2 = randomint(1, scalarRange);
    int op3 = randomint(1, cftData.MaxScalarId());
    OpeCoefficientKey key(op1, op2, op3);

    float_type actual = correlatorSet.CostDerOpeCoefficient(z, key);

    float_type coef = cftData.GetOpeCoefficient(key);
    cftData.SetOpeCoefficient(key, coef - inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value1 = correlatorSet.Cost(z);

    cftData.SetOpeCoefficient(key, coef + inc);
    correlatorSet.UpdateCftData(&cftData);
    float_type value2 = correlatorSet.Cost(z);

    float_type expect = (value2 - value1) / (2 * inc);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CostDerivatives_test, CreateCfgConfigTestData(1), cfg)
{
    CftData cftData(cfg);

    int scalarRange = min(4, cftData.MaxScalarId());
    CorrelatorSet correlatorSet(&cftData, scalarRange);
    cpx_t z = RandomComplex(0.3) + 0.5;

    int op1 = randomint(1, scalarRange);
    int op2 = randomint(1, scalarRange);
    int op3 = randomint(1, cftData.MaxScalarId());
    int op4 = randomint(1, cftData.MaxScalarId());
    OpeCoefficientKey key1(op1, op2, op3);
    OpeCoefficientKey key2(op1, op2, op4);
    vector<OpeCoefficientKey> keys = vector<OpeCoefficientKey>({key1, key2});
    vector<int> ops = vector<int>({op1, op4});
    vector<float_type> opsDerivatives;
    vector<float_type> coefsDerivatives;

    // get the derivatives results
    correlatorSet.CostDerivatives(z, ops, keys, opsDerivatives, coefsDerivatives);

    // update the cft data to be old values - inc.
    for (uint i = 0; i < ops.size(); i++) {
        float_type dim = cftData.GetPrimaryDim(ops[i]);
        cftData.SetPrimaryDim(ops[i], dim - inc);
    }

    for (uint i = 0; i < keys.size(); i++) {
        float_type coef = cftData.GetOpeCoefficient(keys[i]);
        cftData.SetOpeCoefficient(keys[i], coef - inc);
    }

    // get the cost value.
    correlatorSet.UpdateCftData(&cftData);
    float_type value1 = correlatorSet.Cost(z);

    // update the cft data to be old values - inc.
    for (uint i = 0; i < ops.size(); i++) {
        float_type dim = cftData.GetPrimaryDim(ops[i]);
        cftData.SetPrimaryDim(ops[i], dim + 2 * inc);
    }

    for (uint i = 0; i < keys.size(); i++) {
        float_type coef = cftData.GetOpeCoefficient(keys[i]);
        cftData.SetOpeCoefficient(keys[i], coef + 2 * inc);
    }

    // get the cost value.
    correlatorSet.UpdateCftData(&cftData);
    float_type value2 = correlatorSet.Cost(z);

    float_type actual = value2 - value1;
    float_type expect = .0;

    // calculate the expected value based on derivatives.
    for (uint i = 0; i < opsDerivatives.size(); i++) {
        expect += 2 * inc * opsDerivatives[i];
    }

    for (uint i = 0; i < coefsDerivatives.size(); i++) {
        expect += 2 * inc * coefsDerivatives[i];
    }

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_AUTO_TEST_CASE(Performance_test, * utf::label("perf"))
{
    CftConfig cfg = CreateCfgConfigTestData(1)[0];
    CftData cftData(cfg);

    int scalarRange = min(4, cftData.MaxScalarId());
    CorrelatorSet correlatorSet(&cftData, scalarRange);
    cpx_t z = RandomComplex(0.3) + 0.5;

    int op1 = randomint(1, scalarRange);
    int op2 = randomint(1, scalarRange);
    int op3 = randomint(1, cftData.MaxScalarId());
    int op4 = randomint(1, cftData.MaxScalarId());
    OpeCoefficientKey key1(op1, op2, op3);
    OpeCoefficientKey key2(op1, op2, op4);
    OpeCoefficientKey key3(op1, op2, op1);
    OpeCoefficientKey key4(op1, op2, op2);
    vector<OpeCoefficientKey> keys = vector<OpeCoefficientKey>({key1, key2, key3, key4});
    vector<int> ops = vector<int>({op1, op2, op3, op4});
    vector<float_type> opsDerivatives;
    vector<float_type> coefsDerivatives;

    boost::timer stopwatch;

    // get the derivatives results
    correlatorSet.CostDerivatives(z, ops, keys, opsDerivatives, coefsDerivatives);

    BOOST_TEST_MESSAGE("Time to run CostDerivatives: " << stopwatch.elapsed() << " seconds.");

    stopwatch.restart();
    CorrelatorSet correlatorSet2(&cftData, scalarRange);
    for (uint i = 0; i < ops.size(); i++) {
        float_type value = correlatorSet2.CostDerDelta(z, ops[i]);
        MY_FLOAT_EQUAL(value, opsDerivatives[i], tol);
    }

    for (uint i = 0; i < keys.size(); i++) {
        float_type value = correlatorSet2.CostDerOpeCoefficient(z, keys[i]);
        MY_FLOAT_EQUAL(value, coefsDerivatives[i], tol);
    }

    BOOST_TEST_MESSAGE("Time to run CostDerDelta and CostDerOpeCoefficient seperately: " << stopwatch.elapsed() << " seconds.");
}


// test suite end
BOOST_AUTO_TEST_SUITE_END()
