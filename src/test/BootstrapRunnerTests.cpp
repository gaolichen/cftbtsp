//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>

#include "test.h"
#include "../common.h"
#include "../BootstrapRunner.h"

BOOST_FIXTURE_TEST_SUITE(BootstrapRunner_test_suite, SimpleTestFixture, * utf::label("runner"))

BOOST_DATA_TEST_CASE(ConstraintCost_Default_test, CreateCfgConfigTestData(TCNumber), cftConfig)
{
    BoostrapConfig bootstrapConfig;
    BootstrapRunner runner(&bootstrapConfig, cftConfig);
    CftData cftData(cftConfig);

    MY_FLOAT_EQUAL(runner.ConstraintCost(&cftData), 0, tol);
}

BOOST_DATA_TEST_CASE(ConstraintCostDerivative_test1, CreateCfgConfigTestData(TCNumber), cftConfig)
{
    BoostrapConfig bootstrapConfig;
    BootstrapRunner runner(&bootstrapConfig, cftConfig);
    CftData cftData(cftConfig);

    // randomly choose two consective operator and exchange their scaling dimensions.
    int op = randomint(1, cftData.MaxPrimaryId() - 1);
    float_type dim = cftData.GetPrimaryDim(op);
    float_type dim2 = cftData.GetPrimaryDim(op + 1);
    cftData.SetPrimaryDim(op, dim2);
    cftData.SetPrimaryDim(op + 1, dim);

    vector<int> ops = vector<int>({op, op + 1});
    vector<float_type> opsDerivatives(ops.size(), .0);

    float_type cost = runner.ConstraintCost(&cftData);
    if (cftData.GetPrimarySpin(op) == cftData.GetPrimarySpin(op + 1)) {
        MY_FLOAT_EQUAL(cost, (dim2 - dim) * bootstrapConfig.ConstraintFactor, tol);
    }

    // calculate the constraint derivatives for the operators.
    runner.ConstraintCostDerivative(&cftData, ops, opsDerivatives);

    vector<float_type> dims = vector<float_type>({dim2, dim});
    // verify derivatives.
    for (int i = 0; i < ops.size(); i++) {
        cftData.SetPrimaryDim(ops[i], dims[i] - inc);
        float_type value1 = runner.ConstraintCost(&cftData);
        cftData.SetPrimaryDim(ops[i], dims[i] + inc);
        float_type value2 = runner.ConstraintCost(&cftData);

        float_type expected = (value2 - value1) / (2 * inc);

        BOOST_TEST_INFO("op=" << ops[i] << ", dim=" << dims[i]);
        MY_FLOAT_EQUAL(opsDerivatives[i], expected, tol);
    }
}

BOOST_DATA_TEST_CASE(ConstraintCostDerivative_test2, CreateCfgConfigTestData(TCNumber), cftConfig)
{
    BoostrapConfig bootstrapConfig;
    BootstrapRunner runner(&bootstrapConfig, cftConfig);
    CftData cftData(cftConfig);
    
    // randomly choose a spin.
    int l = randomint(0, cftConfig.OperatorNumbers.size() - 1);
    int opId = 1;
    for (int i = 0; i < l; i++) opId += cftConfig.OperatorNumbers[i];

    if (l == 0) {
        cftData.SetPrimaryDim(opId, cftData.D / 2.0 - 1.0 - random(0.0, .05));
    } else {
        cftData.SetPrimaryDim(opId, cftData.D + l - 2.0 - random(0.0, .05));
    }

    // verify cost should not be zero.
    float_type cost = runner.ConstraintCost(&cftData);
    BOOST_TEST(abs(cost) > EPS);

    // test ConstraintCostDerivative.
    vector<int> ops = vector<int>({opId});
    vector<float_type> opsDerivatives(ops.size(), .0);
    // calculate the constraint derivatives for the operators.
    runner.ConstraintCostDerivative(&cftData, ops, opsDerivatives);

    // verify derivatives.
    float_type dim = cftData.GetPrimaryDim(opId);
    cftData.SetPrimaryDim(opId, dim - inc);
    float_type value1 = runner.ConstraintCost(&cftData);
    cftData.SetPrimaryDim(opId, dim + inc);
    float_type value2 = runner.ConstraintCost(&cftData);

    float_type expected = (value2 - value1) / (2 * inc);
    BOOST_TEST_INFO("op=" << opId << ", dim=" << dim);
    MY_FLOAT_EQUAL(opsDerivatives[0], expected, tol);
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()
