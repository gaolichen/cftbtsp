# Conformal Bootstrap Gradient Descent algorithm

A C++ implementation of Gradient Descent algorithm for [Conformal Bootstrap](https://en.wikipedia.org/wiki/Conformal_bootstrap). The idea is explained in the [file](bootstrap.pdf)

# Features
- Solve Conformal Bootstrap equations with Gradient Descent algorithm

# Build & Run
Follow the following steps in a linux system.

- Install [CMake](https://cmake.org/)
- Install C++ library [BOOST](https://www.boost.org/)
- Create a directory `build` under the project root.
- From the `build` directory, run the command line `cmake ../` then run `make`
- Execute `.src/btsptest` to run unittests.

# Examples
The following example is from `main.cpp` file.

```cpp

void RunBootstrap(string prefix) // prefix for output file name.
{
    // setup the configration for Gradient Descent algorithm.
    GradientDescentConfig config;
    config.StepsToPrintResult = 1;
    config.InitialStepSize = 0.01;
    config.LoopNumber = 4;
    config.InitialTry = 5;
    config.SamplesEachStep = 30;
    config.RequiredAccuracy = 1.0;
    config.ConstraintFactor = 5.0;
    config.FilePrefix = prefix;

    // config CFT data
    CftConfig cftConfig(3); // spacetime dimension = 3
    
    // 4 scalar operators, 0 spin-1 operator, 1 spin-2 operator (the stress tensor)
    cftConfig.OperatorNumbers = vector<int>({4, 0, 1});
    CftData cftData(cftConfig);
    
    // create crossing equations.
    // numberOfScalarsToBootstrap = 4
    CrossingEquations equations(&cftData, 4);

    // create equation solver
    EquationSolver<cpx_t> solver(&equations, &config);
    vector<cpx_t> inputs;
    
    // solve the equations
    float_type cost = solver.Run(inputs);
    cout << "inputs=" << inputs << endl;
    cout << "cost=" << cost << endl;
}
```
