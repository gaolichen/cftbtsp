#include <vector>
#include "common.h"
#include "CftData.h"
#include "BootstrapRunner.h"
#include "EquationSolver.h"
#include "RandomAlgebricEquations.h"
#include "CrossingEquations.h"
#include <iostream>

using namespace std;

void Example()
{
    BoostrapConfig bootstrapConfig;
    bootstrapConfig.NumberOfScalarsToBootstrap = 4;
    bootstrapConfig.SamplesEachStep = 20;
    bootstrapConfig.Step = 0.001;
    bootstrapConfig.Accuracy = 0.001;
    bootstrapConfig.ConstraintFactor = 20.0;

    CftConfig cftConfig(3);
    cftConfig.OperatorNumbers = vector<int>({4, 0, 1});
    BootstrapRunner runner(&bootstrapConfig, cftConfig);
    runner.Run();
}

void RunBootstrap(string prefix)
{
    GradientDescentConfig config;
    config.StepsToPrintResult = 1;
    config.InitialStepSize = 0.01;
    config.LoopNumber = 4;
    config.InitialTry = 5;
    config.SamplesEachStep = 30;
    config.RequiredAccuracy = 1.0;
    config.ConstraintFactor = 5.0;
    config.FilePrefix = prefix;

    CftConfig cftConfig(3);
    cftConfig.OperatorNumbers = vector<int>({4, 0, 1});
    CftData cftData(cftConfig);
    CrossingEquations equations(&cftData, 4);

    EquationSolver<cpx_t> solver(&equations, &config);
    vector<cpx_t> inputs;
    float_type cost = solver.Run(inputs);
    cout << "inputs=" << inputs << endl;
    cout << "cost=" << cost << endl;
}

void BootstrapFromFile(string file)
{
    GradientDescentConfig config;
    config.StepsToPrintResult = 1;
    config.InitialStepSize = 0.00001;
    config.LoopNumber = 2;
    config.InitialTry = 1;
    config.SamplesEachStep = 20;
    config.RequiredAccuracy = .01;
    config.ConstraintFactor = 5.0;
    config.FilePrefix = "bootstrap_from_file";

    CftData cftData;
    cftData.LoadFromFile(file);
    CrossingEquations equations(&cftData, 4);

    EquationSolver<cpx_t> solver(&equations, &config);
    vector<cpx_t> inputs;
    float_type cost = solver.Run(inputs);
    cout << "inputs=" << inputs << endl;
    cout << "cost=" << cost << endl;
}


void TestRandomAlgebra(string prefix)
{
    int functionNumber = randomint(3, 4);
    int maxOrder = randomint(2, 3);
    RandomAlgebricEquations equations(functionNumber, maxOrder);
    GradientDescentConfig config;
    config.StepsToPrintResult = 1000;
    config.InitialStepSize = 0.1;
    config.LoopNumber = 5;
    config.InitialTry = 3;
    config.SamplesEachStep = 25;
    config.RequiredAccuracy = 1.0;
    config.FilePrefix = prefix;
    EquationSolver<float_type> solver(&equations, &config);
    vector<float_type> inputs;
    float_type cost = solver.Run(inputs);
    equations.OutputSolutions();
}

int main(int argc, char* argv[])
{
    string cmd = "";
    string arg = "";
	if (argc > 1) cmd = argv[1];
    if (argc > 2) arg = argv[2];
    if (cmd == "-t") {
        TestRandomAlgebra(arg);
    } else if (cmd == "-b") {
        RunBootstrap(arg);
    } else if (cmd == "-bf") {
        BootstrapFromFile(arg);
    }
    return 0;
}

