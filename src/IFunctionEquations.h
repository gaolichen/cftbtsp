#pragma once
#include <vector>
#include <string>
#include "common.h"
using namespace std;

// an interface for a set of function equations.
template<class T> class IFunctionEquations
{
public:
    virtual void GenerateRandomParameters() = 0;

    virtual void SaveToFile(string file = "") = 0;

    virtual int EquationNumber() = 0;

    virtual int ParameterNumber() = 0;

    virtual void SetParameter(int parameterId, float_type value) = 0;

    virtual float_type GetParameter(int parameterId) = 0;

    virtual float_type EvaluateEquation(T input, int equationId) = 0;

    virtual float_type EquationDerivativeByParameter(T input, int equationId, int parameterId) = 0;

    virtual float_type EquationDerivativeByParameter(float_type input, int equationId, int parameterId1, int parameterId2) = 0;

    virtual float_type EvaluateConstraints() = 0;

    virtual float_type ConstraintsDerivativeByParameter(int parameterId) = 0;

    virtual vector<T> GenerateInputs(int size) = 0;

    virtual void OutputParameters() = 0;

    virtual void OnParameterUpdated() = 0;
};
