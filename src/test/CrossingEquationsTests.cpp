//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>

#include <vector>
#include <iostream>
#include <boost/math/special_functions/binomial.hpp>
#include "test.h"
#include "../CrossingEquations.h"
#include "../EquationSolver.h"
using namespace std;

BOOST_FIXTURE_TEST_SUITE(CrossingEquations_suite, SimpleTestFixture, * utf::label("CrossingEquations"))

BOOST_DATA_TEST_CASE(Constructor_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = randomint(1, min(cftData.MaxScalarId(), 5));
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);

    int equationNumber = (int)boost::math::binomial_coefficient<double>(numberOfScalarsToBootstrap + 3, 4);
    if (numberOfScalarsToBootstrap >= 4) {
        equationNumber += (int)boost::math::binomial_coefficient<double>(numberOfScalarsToBootstrap, 4);
    }

    BOOST_TEST(equations.EquationNumber() == equationNumber);

    int parameterNumber = cftData.MaxPrimaryId() - 1;
    parameterNumber += (int)boost::math::binomial_coefficient<double>(numberOfScalarsToBootstrap + 2, 3);
    parameterNumber += (int)boost::math::binomial_coefficient<double>(numberOfScalarsToBootstrap + 1, 2) * (cftData.MaxPrimaryId() - 1 -numberOfScalarsToBootstrap);

    BOOST_TEST(equations.ParameterNumber() == parameterNumber);
}

BOOST_DATA_TEST_CASE(GetSetParameter_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = randomint(1, min(cftData.MaxScalarId(), 5));
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);

    for (uint i = 0; i < equations.ParameterNumber(); i++) {
        float_type value = random(-10.0, 10.0);
        equations.SetParameter(i, value);
        float_type actual = equations.GetParameter(i);
        BOOST_TEST_INFO("i=" << i);
        MY_FLOAT_EQUAL(actual, value, tol);
    }
}

BOOST_DATA_TEST_CASE(EquationDerivativeByParameter_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = randomint(1, min(cftData.MaxScalarId(), 4));
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);
    int equationId = randomint(0, equations.EquationNumber() - 1);
    int parameterId = randomint(0, equations.ParameterNumber() - 1);

    cpx_t input = RandomComplex(0.2) + .5;
    float_type actual = equations.EquationDerivativeByParameter(input, equationId, parameterId);

    float_type param = equations.GetParameter(parameterId);
    equations.SetParameter(parameterId, param - inc);
    equations.OnParameterUpdated();
    float_type value1 = equations.EvaluateEquation(input, equationId);
    equations.SetParameter(parameterId, param + inc);
    equations.OnParameterUpdated();
    float_type value2 = equations.EvaluateEquation(input, equationId);

    float_type expected = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("numberOfScalarsToBootstrap=" << numberOfScalarsToBootstrap << ", input=" << input << ", equationId=" << equationId << ", parameterId=" << parameterId);
    MY_FLOAT_EQUAL(actual, expected, tol);
}

BOOST_DATA_TEST_CASE(ConstraintsDerivativeByParameter_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = randomint(1, min(cftData.MaxScalarId(), 4));
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);
    int parameterId = randomint(0, equations.ParameterNumber() - 2);
    float_type param = equations.GetParameter(parameterId);
    float_type param2 = equations.GetParameter(parameterId + 1);
    equations.SetParameter(parameterId, param2);
    equations.SetParameter(parameterId + 1, param);

    float_type actual = equations.ConstraintsDerivativeByParameter(parameterId);

    equations.SetParameter(parameterId, param2 - inc);
    float_type value1 = equations.EvaluateConstraints();
    equations.SetParameter(parameterId, param2 + inc);
    float_type value2 = equations.EvaluateConstraints();

    float_type expected = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("numberOfScalarsToBootstrap=" << numberOfScalarsToBootstrap << ", parameterId=" << parameterId);
    MY_FLOAT_EQUAL(actual, expected, tol);
}

BOOST_DATA_TEST_CASE(ConstraintsDerivativeByParameter_test2, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = randomint(1, min(cftData.MaxScalarId(), 4));
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);
    int parameterId = 0;
    equations.SetParameter(parameterId, (cftData.D/2.0 - 1 - 0.3));
    float_type actual = equations.ConstraintsDerivativeByParameter(parameterId);

    float_type param = equations.GetParameter(parameterId);
    equations.SetParameter(parameterId, param - inc);
    float_type value1 = equations.EvaluateConstraints();
    equations.SetParameter(parameterId, param + inc);
    float_type value2 = equations.EvaluateConstraints();

    float_type expected = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("numberOfScalarsToBootstrap=" << numberOfScalarsToBootstrap << ", parameterId=" << parameterId << ", value=" << param);
    MY_FLOAT_EQUAL(actual, expected, tol);
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(EquationDerivativeByParameter_Debug_test, * utf::label("CrossingEquations.debug"))
{
    CftConfig cfg(2);
    cfg.OperatorNumbers=vector<int>({5, 4, 4, 5, 4, 3});
    CftData cftData(cfg);
    int numberOfScalarsToBootstrap = 3;
    CrossingEquations equations(&cftData, numberOfScalarsToBootstrap);
    int equationId = 12;
    int parameterId = 2;

    cpx_t input = cpx_t(0.499168,0.00225563);
    float_type actual = equations.EquationDerivativeByParameter(input, equationId, parameterId);
    float_type param = equations.GetParameter(parameterId);
    equations.SetParameter(parameterId, param - inc);
    equations.OnParameterUpdated();
    float_type value1 = equations.EvaluateEquation(input, equationId);
    equations.SetParameter(parameterId, param + inc);
    equations.OnParameterUpdated();
    float_type value2 = equations.EvaluateEquation(input, equationId);

    float_type expected = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("numberOfScalarsToBootstrap=" << numberOfScalarsToBootstrap << ", input=" << input << ", equationId=" << equationId << ", parameterId=" << parameterId);
    MY_FLOAT_EQUAL(actual, expected, tol);
}


