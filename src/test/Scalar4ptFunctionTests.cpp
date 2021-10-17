//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include "test.h"
#include "../Scalar4ptFunction.h"
using namespace std;

struct Scalar4ptFunctionTestData
{
    CftConfig Config;

    int Op1;
    int Op2;
    int Op3;
    int Op4;

    int OpExchange;

    cpx_t Z;
};

std::ostream& operator<< (std::ostream& out, const Scalar4ptFunctionTestData& data)
{
    out << '{';
    out << "Config=" << data.Config << ", Op1=" << data.Op1 << ", Op2=" << data.Op2;
    out << ", Op3=" << data.Op3 << ", Op4=" << data.Op4 << ", OpExchange=" << data.OpExchange;
    out << ", Z=" << data.Z;
    out << '}';
    return out;
}

Scalar4ptFunctionTestData CreateSingleScalar4ptFunctionTestData(CftConfig* config, bool useExternalAsExchangeOp)
{
    Scalar4ptFunctionTestData data;
    data.Config = *config;
    data.Op1 = randomint(1, config->OperatorNumbers[0]);
    data.Op2 = randomint(1, config->OperatorNumbers[0]);
    data.Op3 = randomint(1, config->OperatorNumbers[0]);
    data.Op4 = randomint(1, config->OperatorNumbers[0]);

    // if useExternalAsExchangeOp, OpExchange will be one of the external operators.
    if (useExternalAsExchangeOp) {
        int ops[4] = {data.Op1, data.Op2, data.Op3, data.Op4};
        data.OpExchange = ops[randomint(0,3)];
    } else {
        data.OpExchange = randomint(0, Sum(config->OperatorNumbers));
        if (data.OpExchange == data.Op1 || data.OpExchange == data.Op2
            || data.OpExchange == data.Op3 || data.OpExchange == data.Op4) {
            data.OpExchange = randomint(config->OperatorNumbers[0] + 1, Sum(config->OperatorNumbers));
        }

        // if OpExchange is stress tensor, some derivative functions will blow up.
        // so we pick another operator if OpExchange is stress tensor.
        if (data.OpExchange == config->OperatorNumbers[0] + config->OperatorNumbers[1] + 1) {
            data.OpExchange = data.OpExchange + 1;
        }
    }

    data.Z = RandomComplex(0.3) + 0.5;
    float_type r1 = Z2r(data.Z);
    float_type r2 = Z2r(-data.Z + 1.0);

    if (r1 >= 0.25 || r2 >= 0.25) {
        //std::cout << "Z=" << data.Z << ", r1=" << r1 << ", r2=" << r2 << endl;
    }

    return data;
}

vector<Scalar4ptFunctionTestData> CreateScalar4ptFunctionTestData(int size, bool useExternalAsExchangeOp)
{
    vector<CftConfig> configs = CreateCfgConfigTestData(size);
    vector<Scalar4ptFunctionTestData> ret;

    for (uint i = 0; i < configs.size(); i++) {
        ret.push_back(CreateSingleScalar4ptFunctionTestData(&configs[i], useExternalAsExchangeOp));
    }

    return ret;
}

// test suite begin
BOOST_FIXTURE_TEST_SUITE(Scalar4ptFunction_NotExternalAsExchangeOp_suite, SimpleTestFixture, * utf::label("s4pt1"))

BOOST_DATA_TEST_CASE(Evaluate_test, CreateScalar4ptFunctionTestData(TCNumber, false), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    float_type actual = s4f.CrossingDiff(sd.Z);
    float_type expect = s4f.EvaluateS(sd.Z) - s4f.EvaluateT(sd.Z);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerDelta_test, CreateScalar4ptFunctionTestData(TCNumber, false), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    float_type actual = s4f.CrossingDerDelta(sd.Z, sd.OpExchange);

    float_type dim = cftData.GetPrimaryDim(sd.OpExchange);
    cftData.SetPrimaryDim(sd.OpExchange, dim - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z);

    cftData.SetPrimaryDim(sd.OpExchange, dim + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z);
    float_type expect = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("OpExchange: dlt=" << dim << ", l=" << cftData.GetPrimarySpin(sd.OpExchange) << ", value1=" << value1 << ", value2=" << value2);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerOpeCoefficient_op_test, CreateScalar4ptFunctionTestData(TCNumber, false), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    OpeCoefficientKey key(sd.Op1, sd.Op2, sd.OpExchange);
    float_type actual = s4f.CrossingDerOpeCoefficient(sd.Z, key, sd.OpExchange);

    float_type coef = cftData.GetOpeCoefficient(key.Op1, key.Op2, key.Op3);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z, sd.OpExchange);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z, sd.OpExchange);
    float_type expect = (value2 - value1) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerOpeCoefficient_test, CreateScalar4ptFunctionTestData(TCNumber, false), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    OpeCoefficientKey key(sd.Op3, sd.Op2, sd.OpExchange);
    float_type actual = s4f.CrossingDerOpeCoefficient(sd.Z, key);

    float_type coef = cftData.GetOpeCoefficient(key.Op1, key.Op2, key.Op3);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z);
    float_type expect = (value2 - value1) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()

// test suite begin
BOOST_FIXTURE_TEST_SUITE(Scalar4ptFunction_ExternalAsExchangeOp_suite, SimpleTestFixture, * utf::label("s4pt2"))

BOOST_DATA_TEST_CASE(Evaluate_test, CreateScalar4ptFunctionTestData(TCNumber, true), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    float_type actual = s4f.CrossingDiff(sd.Z);
    float_type expect = s4f.EvaluateS(sd.Z) - s4f.EvaluateT(sd.Z);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerDelta_test, CreateScalar4ptFunctionTestData(TCNumber, true), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    float_type actual = s4f.CrossingDerDelta(sd.Z, sd.OpExchange);

    float_type dim = cftData.GetPrimaryDim(sd.OpExchange);
    cftData.SetPrimaryDim(sd.OpExchange, dim - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z);

    cftData.SetPrimaryDim(sd.OpExchange, dim + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z);
    float_type expect = (value2 - value1) / (2 * inc);

    BOOST_TEST_INFO("OpExchange: dlt=" << dim << ", l=" << cftData.GetPrimarySpin(sd.OpExchange) << ", value1=" << value1 << ", value2=" << value2);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerOpeCoefficient_op_test, CreateScalar4ptFunctionTestData(TCNumber, true), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    OpeCoefficientKey key(sd.Op1, sd.Op2, sd.OpExchange);
    float_type actual = s4f.CrossingDerOpeCoefficient(sd.Z, key, sd.OpExchange);

    float_type coef = cftData.GetOpeCoefficient(key.Op1, key.Op2, key.Op3);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z, sd.OpExchange);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z, sd.OpExchange);
    float_type expect = (value2 - value1) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}

BOOST_DATA_TEST_CASE(CrossingDerOpeCoefficient_test, CreateScalar4ptFunctionTestData(TCNumber, true), sd)
{
    CftData cftData(sd.Config);
    Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

    OpeCoefficientKey key(sd.Op3, sd.Op2, sd.OpExchange);
    float_type actual = s4f.CrossingDerOpeCoefficient(sd.Z, key);

    float_type coef = cftData.GetOpeCoefficient(key.Op1, key.Op2, key.Op3);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef - inc);
    s4f.UpdateCftData(&cftData);
    float_type value1 = s4f.CrossingDiff(sd.Z);

    cftData.SetOpeCoefficient(key.Op1, key.Op2, key.Op3, coef + inc);
    s4f.UpdateCftData(&cftData);
    float_type value2 = s4f.CrossingDiff(sd.Z);
    float_type expect = (value2 - value1) / (2 * inc);

    MY_FLOAT_EQUAL(actual, expect, tol);
}


// test suite end
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(Scalar4ptFunction_performance_test1, * utf::label("perf"))
{
    CftConfig config(3);
    config.OperatorNumbers = vector<int>({8,6,6,5,4,3,2,1});
    vector<Scalar4ptFunctionTestData> dataList;
    for (int i = 0; i < TCNumber; i++) {
        dataList.push_back(CreateSingleScalar4ptFunctionTestData(&config, true));
    }

    float_type value;
    boost::timer stopwatch;

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        value = s4f.CrossingDiff(sd.Z);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDiff " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");

    stopwatch.restart();

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        value = s4f.CrossingDerDelta(sd.Z, sd.OpExchange);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDerDelta " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");

    stopwatch.restart();

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        OpeCoefficientKey coef(sd.Op3, sd.Op2, sd.OpExchange);
        value = s4f.CrossingDerOpeCoefficient(sd.Z, coef);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDerOpeCoefficient " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");
}

BOOST_AUTO_TEST_CASE(Scalar4ptFunction_performance_test2, * utf::label("perf"))
{
    CftConfig config(3);
    config.OperatorNumbers = vector<int>({8,6,6,5,4,3,2,1});
    vector<Scalar4ptFunctionTestData> dataList;
    for (int i = 0; i < TCNumber; i++) {
        dataList.push_back(CreateSingleScalar4ptFunctionTestData(&config, false));
    }

    float_type value;
    boost::timer stopwatch;

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        value = s4f.CrossingDiff(sd.Z);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDiff " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");

    stopwatch.restart();

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        value = s4f.CrossingDerDelta(sd.Z, sd.OpExchange);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDerDelta " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");

    stopwatch.restart();

    for (uint i = 0; i < dataList.size(); i++) {
        Scalar4ptFunctionTestData sd = dataList[i];
        CftData cftData(sd.Config);
        Scalar4ptFunction s4f(&cftData, sd.Op1, sd.Op2, sd.Op3, sd.Op4);

        OpeCoefficientKey coef(sd.Op3, sd.Op2, sd.OpExchange);
        value = s4f.CrossingDerOpeCoefficient(sd.Z, coef);
    }

    BOOST_TEST_MESSAGE("Run Scalar4ptFunction.CrossingDerOpeCoefficient " << dataList.size() << " times takes " << stopwatch.elapsed() << " seconds.");
}

BOOST_AUTO_TEST_CASE(Scalar4ptFunction_Debug, * utf::label("debug"))
{
    vector<float_type> dims = vector<float_type>({0.509292, 1.37473, 2.23385, 3.20297});
    cpx_t z(0.41671,0.0386725);
//    cpx_t z(0.598759,0.00609649);

    CftConfig cftConfig(3);
    cftConfig.OperatorNumbers = vector<int>({4, 0, 1});
    CftData cftData(cftConfig);
    for (uint i = 0; i < dims.size(); i++) {
        cftData.SetPrimaryDim(i + 1, dims[i]);
    }

    Scalar4ptFunction s4f(&cftData, 1, 1, 4, 4);
    //cout << "S channel=" << s4f.EvaluateS(z) << ", T channel=" << s4f.EvaluateT(z) << endl;
    cout << "T cb=" << s4f.ConformalBlock_T(z, 1) << endl;
}

