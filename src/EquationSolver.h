#pragma once
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <boost/timer.hpp>
#include "common.h"
#include "IFunctionEquations.h"
using namespace std;

struct GradientDescentConfig
{
    float_type InitialStepSize;
    int InitialTry;
    float_type ConstraintFactor;
    float_type RequiredAccuracy;

    int StepsToPrintResult;
    int SamplesEachStep;
    int LoopNumber;

    string FilePrefix;

    GradientDescentConfig()
    {
        InitialStepSize = 0.001;
        InitialTry = 5;
        ConstraintFactor = 20.0;
        StepsToPrintResult = 1;
        SamplesEachStep = 20;
        LoopNumber = 3;
        RequiredAccuracy = 1e-14;
    }
};

template<class TInput> class EquationSolver
{
private:
    IFunctionEquations<TInput> *equations;
    GradientDescentConfig* config;

    float_type GradientDescent(vector<TInput>& inputs, float_type stepSize)
    {
        boost::timer stopwatch;
        int equationNumber = equations->EquationNumber();
        int parameterNumber = equations->ParameterNumber();
        vector<float_type> bestParameters(parameterNumber, .0);
        vector<int> parameters(parameterNumber, 0);
        for (int i = 0; i < parameterNumber; i++) {
            parameters[i] = i;
        }

        float_type bestCost = std::numeric_limits<float_type>::max();
        int bestStep = -1;
        float_type cost;
        float_type constraints;
        float_type dcost = .0;

        int totSteps = 0;
        while (true) {
            cost = .0;
            vector<float_type> derivatives(parameterNumber, 0.0);
            for (uint i = 0; i < inputs.size(); i++) {
//                float_type maxEquation = .0;
//                int maxEqId = -1;
                TInput& input = inputs[i];
                vector<float_type> res(equationNumber, 0.0);
                for (int eq = 0; eq < equationNumber; eq++) {
                    res[eq] = equations->EvaluateEquation(input, eq);
//                    if (abs(res[eq]) > abs(maxEquation)) {
//                        maxEquation = res[eq];
//                        maxEqId = eq;
//                    }
                    cost += res[eq] * res[eq];
                }

//                cout << "input=" << input << ", r1=" << Z2r(input) << ", r2=" << Z2r(-input + 1.0) << ", maxEqId=" << maxEqId << ", maxEquationValue=" << maxEquation << endl;

                for (int para = 0; para < parameterNumber; para++) {
                    for (int eq = 0; eq < equationNumber; eq++) {
                        derivatives[para] += 2 * res[eq] * equations->EquationDerivativeByParameter(input, eq, para);
                    }
                }
            }

            cost /= inputs.size();
            constraints = config->ConstraintFactor * equations->EvaluateConstraints();

            if (totSteps % config->StepsToPrintResult == 0) {
                std::cout << "Time elapsed: " << stopwatch.elapsed() << " seconds." << endl;
                cout << "cost=" << setprecision(12) << cost << "(+" << constraints << "), step=" << totSteps << ", bestCost=" << setprecision(12) << bestCost << ", bestStep=" << bestStep << endl;
                equations->OutputParameters();
            }

            cost += constraints;

            if (cost < config->RequiredAccuracy || (cost >= bestCost) || bestCost - cost < dcost * stepSize / 10) {
                break;
            }

            if (cost < bestCost) {
                bestCost = cost;
                bestStep = totSteps;

                for (int para = 0; para < parameterNumber; para++) {
                    bestParameters[para] = equations->GetParameter(para);
                }
            }

            totSteps++;

            float_type gradNorm = .0;
            for (int para = 0; para < parameterNumber; para++) {
                derivatives[para] /= inputs.size();
                derivatives[para] += config->ConstraintFactor * equations->ConstraintsDerivativeByParameter(para);
                gradNorm += derivatives[para] * derivatives[para];
            }

            if (totSteps % config->StepsToPrintResult == 0) {
                cout << "gradNorm=" << sqrt(gradNorm) << ", grad="<< derivatives << endl;
            }

            float_type factor = 1.0;
            factor = sqrt(gradNorm / parameterNumber);
            dcost = sqrt(parameterNumber) * sqrt(gradNorm);

            for (int para = 0; para < parameterNumber; para++) {
                float_type value = equations->GetParameter(para);
                equations->SetParameter(para, value - derivatives[para] * stepSize / factor);
            }

            equations->OnParameterUpdated();
        }

        // set parameter back to lastParameters if cost is larger than lastCost.
        for (int para = 0; para < parameterNumber; para++) {
            equations->SetParameter(para, bestParameters[para]);
        }
        equations->OnParameterUpdated();

        return bestCost;
    }

public:
    EquationSolver(IFunctionEquations<TInput> *equations, GradientDescentConfig* config)
    {
        this->equations = equations;
        this->config = config;
    }

    float_type Run(vector<TInput>& bestInputs)
    {
        boost::timer stopwatch;

        float_type stepSize = config->InitialStepSize;
        bestInputs.resize(config->SamplesEachStep, (TInput).0);
        vector<float_type> bestParameters(equations->ParameterNumber(), .0);
        float_type bestCost = std::numeric_limits<float_type>::max();

        // find the best inputs
        for (int i = 0; i < config->InitialTry; i++) {
            cout << "==========================Try #" << i << endl;
            vector<TInput> inputs = equations->GenerateInputs(config->SamplesEachStep);
            float_type cost = GradientDescent(inputs, stepSize);
            if (cost < bestCost) {
                bestCost = cost;
                copy(inputs.begin(), inputs.end(), bestInputs.begin());

                for (int j = 0; j < equations->ParameterNumber(); j++) {
                    bestParameters[j] = equations->GetParameter(j);
                }
            }

            equations->SaveToFile(config->FilePrefix + "_try" + ToString(i) + ".txt");

            if (i + 1 != config->InitialTry) {
                equations->GenerateRandomParameters();
            }
        }

        // set the parameter to the one corresponds to best inputs.
        for (int j = 0; j < equations->ParameterNumber(); j++) {
            equations->SetParameter(j, bestParameters[j]);
        }
        equations->OnParameterUpdated();

        // perform finer gradient descent procedure.
        for (int loop = 1; loop < config->LoopNumber; loop++) {
            cout << "==========================Start loop " << loop << endl;
            stepSize /= 10.0;
            float_type cost = GradientDescent(bestInputs, stepSize);
            cout << "cost="  << setprecision(12) << cost << endl;
            cout << "Final parameter: " << endl;
            equations->OutputParameters();
        }
        equations->SaveToFile(config->FilePrefix + "_final.txt");

        cout << "Time eslapsed: " << stopwatch.elapsed() << " seconds" << endl;
    }

    void CalculateHessian(vector<TInput> &inputs, vector<vector<float_type> >& hessian) {
        int size = equations->ParameterNumber();
        hessian.resize(size, vector<float_type>(size, 0.0));

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                for (uint k = 0; k < inputs.size(); k++) {
                    for (int eq = 0; eq < equations->EquationNumber(); eq++) {
                        hessian[i][j] += 2 * (equations->EquationDerivativeByParameter(inputs[k], eq, i) * equations->EquationDerivativeByParameter(inputs[k], eq, j));
                        hessian[i][j] += 2 * equations->EvaluateEquation(inputs[k], eq) * equations->EquationDerivativeByParameter(inputs[k], eq, i, j);
                    }
                }
                hessian[i][j] /= inputs.size();
            }
        }
    }
};
