//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include "test.h"
#include "../common.h"
#include "../CftData.h"

BOOST_FIXTURE_TEST_SUITE(CftData_suite, SimpleTestFixture, * utf::label("CftData"))

BOOST_DATA_TEST_CASE(CftData_construct_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData data(cfg);
    BOOST_TEST(data.D == cfg.D);
    BOOST_TEST(data.MaxPrimaryId() == Sum(cfg.OperatorNumbers));
    BOOST_TEST(data.MaxScalarId() == cfg.OperatorNumbers[0]);

    int stressTensorId = data.StressTensorId;
    BOOST_TEST(data.GetPrimaryInfo(stressTensorId).Spin == 2);
    BOOST_TEST(data.GetPrimaryInfo(stressTensorId).Dim == cfg.D);

    float_type Pi = acos(-1.0);
    for (int i = 1; i <= data.MaxScalarId(); i++) {
        MY_FLOAT_EQUAL(data.GetOpeCoefficient(i, i, 0), 1.0, tol);

        float_type prefactor = - data.D * data.GetPrimaryDim(i) / (data.D - 1);
        float_type sd = 2 * pow(Pi, data.D/2.0) / boost::math::tgamma<float_type>(data.D/2.0);
        MY_FLOAT_EQUAL(data.GetOpeCoefficient(i, i, data.StressTensorId), prefactor / sd, tol);
    
        for (int j = i + 1; j <= data.MaxScalarId(); j++) {
            BOOST_TEST(data.GetOpeCoefficient(i, j, 0) == 0.0);
            BOOST_TEST(data.GetOpeCoefficient(i, j, data.StressTensorId) == 0.0);
        }
    }
}

BOOST_DATA_TEST_CASE(CftData_SaveLoad_test, CreateCfgConfigTestData(TCNumber) ^ bdata::xrange(TCNumber), cfg, index)
{
    CftData data(cfg);

    string file = "CftData_SaveLoad_test" + ToString(index) + ".txt";
    data.Save(file);

    CftData data2;
    data2.LoadFromFile(file);

    BOOST_TEST(data2.D == data.D);
    BOOST_TEST(data2.StressTensorId == data.StressTensorId);
    BOOST_TEST(data2.MaxPrimaryId() == data.MaxPrimaryId());
    BOOST_TEST(data2.MaxScalarId() == data.MaxScalarId());

    for (int i = 0; i <= data.MaxPrimaryId(); i++) {
        BOOST_TEST(data2.GetPrimarySpin(i) == data.GetPrimarySpin(i));
        MY_FLOAT_EQUAL(data2.GetPrimaryDim(i), data.GetPrimaryDim(i), tol);
    }

    for (int i = 1; i <= data.MaxScalarId(); i++) {
        for (int j = i; j <= data.MaxScalarId(); j++) {
            MY_FLOAT_EQUAL(data2.GetOpeCoefficient(i, j, 0), data.GetOpeCoefficient(i, j, 0), tol);
        }
    }

    for (int i = 1; i <= data.MaxScalarId(); i++) {
        for (int j = i; j <= data.MaxScalarId(); j++) {
            for (int k = j; k <= data.MaxPrimaryId(); k++) {
                BOOST_TEST_INFO("i=" << i << ", j=" << j << ", k=" << k);
                MY_FLOAT_EQUAL(data2.GetOpeCoefficient(i, j, k), data.GetOpeCoefficient(i, j, k), tol);
            }
        }
    }

}

BOOST_DATA_TEST_CASE(CftData_PrimaryNumber_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData data(cfg);
    for (int i = 0; i < cfg.OperatorNumbers.size(); i++) {
        BOOST_TEST(data.PrimaryNumber(i) == cfg.OperatorNumbers[i]);
    }
}

BOOST_DATA_TEST_CASE(CftData_GetSetPrimaryInfo_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData data(cfg);
    int id = cfg.OperatorNumbers[0] + randomint(1, cfg.OperatorNumbers[1]);
    double dim = random(2.0, 4.0);

    data.SetPrimaryDim(id, dim);
    PrimaryInfo info = data.GetPrimaryInfo(id);

    BOOST_TEST(info.Id == id);
    BOOST_TEST(info.Spin == 1);
    MY_FLOAT_EQUAL(info.Dim, dim, tol);
}

BOOST_DATA_TEST_CASE(CftData_GetSetOpeCoefficient_test, CreateCfgConfigTestData(TCNumber), cfg)
{
    CftData data(cfg);
    int id1 = randomint(0, data.MaxPrimaryId());
    int id2 = randomint(0, data.MaxPrimaryId());
    int id3 = randomint(0, data.MaxPrimaryId());
    float_type expect = random(-2.0, 2.0);

    data.SetOpeCoefficient(id1, id2, id3, expect);
    float_type actual = data.GetOpeCoefficient(id2, id3, id1);
    MY_FLOAT_EQUAL(actual, expect, tol);
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()
