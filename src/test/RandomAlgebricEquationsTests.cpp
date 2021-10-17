//Link to Boost
 #define BOOST_TEST_DYN_LINK

//VERY IMPORTANT - include this last
#include <boost/test/unit_test.hpp>
#include <boost/timer.hpp>

#include <vector>
#include <iostream>
#include "test.h"
#include "../RandomAlgebricEquations.h"
#include "../EquationSolver.h"
using namespace std;

BOOST_FIXTURE_TEST_SUITE(RandomAlgebricEquations_suite, SimpleTestFixture, * utf::label("RandomAlgebricEquations"))

BOOST_DATA_TEST_CASE(Constructor_test, bdata::random(5, 10) ^ bdata::random(0, 2) ^ bdata::random(2, 5) ^ bdata::xrange(TCNumber), functionNumber, orderDiff, tryOrder, index)
{
    int maxOrder = orderDiff + tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    BOOST_TEST(equations.ParameterNumber(), functionNumber * (tryOrder + 1));
    BOOST_TEST(equations.SolutionNumber(), functionNumber * (maxOrder + 1));
    BOOST_TEST(equations.EquationNumber(), functionNumber);
//    MY_FLOAT_EQUAL(equations.EvaluateConstraints(), .0, tol);

    vector<float_type> inputs = equations.GenerateInputs(3);
    for (int i = 0; i < inputs.size(); i++) {
        for (int j = 0; j < equations.EquationNumber(); j++) {
            float_type res = equations.VerifySolution(inputs[i], j);
            BOOST_TEST_INFO("input=" << inputs[i] << ", equationId=" << j);
            MY_FLOAT_EQUAL(res, 0.0, tol);
        }
    }
}

BOOST_DATA_TEST_CASE(SaveLoad_test, bdata::random(5, 10) ^ bdata::random(0, 2) ^ bdata::random(2, 5) ^ bdata::xrange(TCNumber), functionNumber, orderDiff, tryOrder, index)
{
    int maxOrder = orderDiff + tryOrder;
    RandomAlgebricEquations eqs1(functionNumber, maxOrder, tryOrder);
    string file = "SaveLoad_test" + ToString(index) + ".txt";
    eqs1.SaveToFile(file);
    RandomAlgebricEquations eqs2;
    eqs2.LoadFromFile(file);

    BOOST_TEST(eqs2.ParameterNumber(), eqs1.ParameterNumber());
    BOOST_TEST(eqs2.SolutionNumber(), eqs1.SolutionNumber());
    BOOST_TEST(eqs2.EquationNumber(), eqs1.EquationNumber());

    // compare parameters and solutions.
    for (int i = 0; i < eqs2.ParameterNumber(); i++) {
        BOOST_TEST_INFO("i=" << i);
        MY_FLOAT_EQUAL(eqs2.GetParameter(i), eqs1.GetParameter(i), tol);
    }

    for (int i = 0; i < eqs2.SolutionNumber(); i++) {
        BOOST_TEST_INFO("i=" << i);
        MY_FLOAT_EQUAL(eqs2.GetSolution(i), eqs1.GetSolution(i), tol);
    }

    // compare coefficients and constTerms.
    for (int i = 0; i < eqs2.EquationNumber(); i++) {
        for (int j = 0; j < functionNumber; j++) {
            for (int k = 0; k < functionNumber; k++) {
                MY_FLOAT_EQUAL(eqs2.GetCoefficient(i, j, k), eqs1.GetCoefficient(i, j, k), tol);
            }
        }

        for (int j = 0; j <= maxOrder; j++) {
            MY_FLOAT_EQUAL(eqs2.GetConstTerm(i, j), eqs1.GetConstTerm(i, j), tol);
        }
    }
}

BOOST_DATA_TEST_CASE(GetSetParameter_test, bdata::random(5, 10) ^ bdata::random(0, 2) ^ bdata::random(2, 5) ^ bdata::xrange(TCNumber), functionNumber, orderDiff, tryOrder, index)
{
    int maxOrder = orderDiff + tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    int parameterId = randomint(0, equations.ParameterNumber() - 1);
    float_type expected = random(-10.0, 10.0);
    equations.SetParameter(parameterId, expected);
    float_type actual = equations.GetParameter(parameterId);

    MY_FLOAT_EQUAL(actual, expected, tol);
}

BOOST_DATA_TEST_CASE(EquationDerivativeByParameter_test, bdata::random(5, 10) ^ bdata::random(0, 2) ^ bdata::random(2, 5) ^ bdata::xrange(TCNumber), functionNumber, orderDiff, tryOrder, index)
{
    int maxOrder = orderDiff + tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    int parameterId = randomint(0, equations.ParameterNumber() - 1);
    int equationId = randomint(0, equations.EquationNumber() - 1);
    float_type input = random(-1.0, 1.0);

    float_type actual = equations.EquationDerivativeByParameter(input, equationId, parameterId);

    float_type param = equations.GetParameter(parameterId);
    equations.SetParameter(parameterId, param - inc);
    float_type value1 = equations.EvaluateEquation(input, equationId);

    equations.SetParameter(parameterId, param + inc);
    float_type value2 = equations.EvaluateEquation(input, equationId);

    float_type expected = (value2 - value1) / (2 * inc);
    BOOST_TEST_INFO("input=" << input << ", equationId=" << equationId << ", parameterId=" << parameterId);
    MY_FLOAT_EQUAL(actual, expected, tol);
}

BOOST_DATA_TEST_CASE(ConstraintsDerivativeByParameter_test, bdata::random(5, 10) ^ bdata::random(0, 2) ^ bdata::random(2, 5) ^ bdata::xrange(TCNumber), functionNumber, orderDiff, tryOrder, index)
{
    int maxOrder = orderDiff + tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    int parameterId = randomint(0, tryOrder);
    float_type param = random(-.1, .0);
    equations.SetParameter(parameterId, param);

    float_type actual = equations.EvaluateConstraints();
//    MY_FLOAT_EQUAL(actual, -param, tol);

    actual = equations.ConstraintsDerivativeByParameter(functionNumber + 1);
    MY_FLOAT_EQUAL(actual, .0, tol);

    actual = equations.ConstraintsDerivativeByParameter(parameterId);
    equations.SetParameter(parameterId, param - inc);
    float_type value1 = equations.EvaluateConstraints();

    equations.SetParameter(parameterId, param + inc);
    float_type value2 = equations.EvaluateConstraints();

    float_type expected = (value2 - value1) / (2 * inc);
    BOOST_TEST_INFO("tryOrder=" << tryOrder << ", parameterId=" << parameterId << ", value=" << param);
    MY_FLOAT_EQUAL(actual, expected, tol);
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(EquationSolver_suite, SimpleTestFixture, * utf::disabled())

BOOST_AUTO_TEST_CASE(EquationSolver_test1, * utf::label("EquationSolver1"))
{
    int functionNumber = randomint(3, 4);
    int tryOrder = randomint(2, 3);
    int maxOrder = tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    GradientDescentConfig config;
    config.StepsToPrintResult = 1000;
    config.InitialStepSize = 0.1;
    config.LoopNumber = 5;
    config.InitialTry = 3;
    config.SamplesEachStep = 25;
    config.RequiredAccuracy = 1.0;
    config.FilePrefix = "unittest1";
    EquationSolver<float_type> solver(&equations, &config);
    vector<float_type> inputs;
    float_type cost = solver.Run(inputs);
    equations.OutputSolutions();

    if (cost < config.RequiredAccuracy) {
        float_type diff = .0;
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            diff += (equations.GetParameter(i) - equations.GetSolution(i)) * (equations.GetParameter(i) - equations.GetSolution(i));
        }

        BOOST_TEST(diff < 1.0);
    } else {
        cout << "Warning: only local minimum is returned.";
        equations.SaveToFile("EquationSolver_test.txt");

        vector<vector<float_type> > hessian;
        solver.CalculateHessian(inputs, hessian);

        std::cout << "save hessian to file." << endl;
        ofstream out("EquationSolver_test_hessian.txt");
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            out << "{";
            for (int j = 0; j < equations.ParameterNumber(); j++) {
                if (j > 0) out << ", ";
                out << hessian[i][j];
            }
            out << "}," << endl;
        }

        out.close();
    }
}

BOOST_AUTO_TEST_CASE(EquationSolver_test2, * utf::label("EquationSolver2"))
{
    int functionNumber = randomint(3, 4);
    int tryOrder = randomint(2, 3);
    int maxOrder = randomint(1, 2) + tryOrder;
    RandomAlgebricEquations equations(functionNumber, maxOrder, tryOrder);
    GradientDescentConfig config;
    config.StepsToPrintResult = 1000;
    config.InitialStepSize = 0.1;
    config.LoopNumber = 6;
    config.InitialTry = 3;
    config.SamplesEachStep = 25;
    config.RequiredAccuracy = 1.0;
    config.FilePrefix = "unittest2";
    EquationSolver<float_type> solver(&equations, &config);
    vector<float_type> inputs;
    float_type cost = solver.Run(inputs);
    equations.OutputSolutions();

    if (cost < config.RequiredAccuracy) {
        float_type diff = .0;
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            diff += (equations.GetParameter(i) - equations.GetSolution(i)) * (equations.GetParameter(i) - equations.GetSolution(i));
        }

        BOOST_TEST(diff < 1.0);
    } else {
        cout << "Warning: only local minimum is returned.";
        equations.SaveToFile("EquationSolver_test.txt");

        vector<vector<float_type> > hessian;
        solver.CalculateHessian(inputs, hessian);

        std::cout << "save hessian to file." << endl;
        ofstream out("EquationSolver_test_hessian.txt");
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            out << "{";
            for (int j = 0; j < equations.ParameterNumber(); j++) {
                if (j > 0) out << ", ";
                out << hessian[i][j];
            }
            out << "}," << endl;
        }

        out.close();
    }
}

BOOST_AUTO_TEST_CASE(EquationSolver_test3, * utf::label("EquationSolver3"))
{
    RandomAlgebricEquations equations;
    equations.LoadFromFile("RandomAlgebricEquations_unittest3_input.txt");
    GradientDescentConfig config;
    config.StepsToPrintResult = 1000;
    config.InitialStepSize = 0.1;
    config.LoopNumber = 6;
    config.InitialTry = 3;
    config.SamplesEachStep = 20;
    config.RequiredAccuracy = 1.0;
    config.FilePrefix = "unittest3";
    EquationSolver<float_type> solver(&equations, &config);
    vector<float_type> inputs;
    float_type cost = solver.Run(inputs);
    equations.OutputSolutions();

    if (cost < config.RequiredAccuracy) {
        float_type diff = .0;
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            diff += (equations.GetParameter(i) - equations.GetSolution(i)) * (equations.GetParameter(i) - equations.GetSolution(i));
        }

        BOOST_TEST(diff < 1.0);
    } else {
        cout << "Warning: only local minimum is returned.";
        equations.SaveToFile("EquationSolver_test2.txt");

        vector<vector<float_type> > hessian;
        solver.CalculateHessian(inputs, hessian);

        std::cout << "save hessian to file." << endl;
        ofstream out("EquationSolver_test2_hessian.txt");
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            out << "{";
            for (int j = 0; j < equations.ParameterNumber(); j++) {
                if (j > 0) out << ", ";
                out << hessian[i][j];
            }
            out << "}," << endl;
        }

        out.close();
    }
}

BOOST_AUTO_TEST_CASE(EquationSolver_test4, * utf::label("EquationSolver4"))
{
    RandomAlgebricEquations equations;
    equations.LoadFromFile("RandomAlgebricEquations_unittest4_input.txt");
    GradientDescentConfig config;
    config.StepsToPrintResult = 1000;
    config.InitialStepSize = 0.1;
    config.LoopNumber = 5;
    config.InitialTry = 3;
    config.SamplesEachStep = 20;
    config.RequiredAccuracy = 1.0;
    config.FilePrefix = "unittest4";
    EquationSolver<float_type> solver(&equations, &config);
    vector<float_type> inputs;
    float_type cost = solver.Run(inputs);
    equations.OutputSolutions();

    if (cost < config.RequiredAccuracy) {
        float_type diff = .0;
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            diff += (equations.GetParameter(i) - equations.GetSolution(i)) * (equations.GetParameter(i) - equations.GetSolution(i));
        }

        BOOST_TEST(diff < 1.0);
    } else {
        cout << "Warning: only local minimum is returned.";
        equations.SaveToFile("EquationSolver_test2.txt");

        vector<vector<float_type> > hessian;
        solver.CalculateHessian(inputs, hessian);

        std::cout << "save hessian to file." << endl;
        ofstream out("EquationSolver_test2_hessian.txt");
        for (int i = 0; i < equations.ParameterNumber(); i++) {
            out << "{";
            for (int j = 0; j < equations.ParameterNumber(); j++) {
                if (j > 0) out << ", ";
                out << hessian[i][j];
            }
            out << "}," << endl;
        }

        out.close();
    }
}

// test suite end
BOOST_AUTO_TEST_SUITE_END()
