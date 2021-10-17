#pragma once
#include <vector>
#include <string>
#include "common.h"
#include "CftData.h"
#include "Scalar4ptFunction.h"
#include "IFunctionEquations.h"
using namespace std;

class CrossingEquations : public IFunctionEquations<cpx_t>
{
private:
    CftData* cftData;
    vector<Scalar4ptFunction*> correlators;
    int operatorNumber;
    float_type unitaryBoundaryMinGap;
    float_type operatorMinGap;

    // the vector is used for mapping of parameterId to OpeCoefficientKey
    vector<OpeCoefficientKey> parameter2OpeCoefficientKey;
public:
    CrossingEquations(CftData* data, int numberOfScalarsToBootstrap);
    ~CrossingEquations();

    void GenerateRandomParameters();

    void SaveToFile(string file = "");

    int EquationNumber();

    int ParameterNumber();

    float_type GetParameter(int parameterId);

    void SetParameter(int parameterId, float_type value);

    float_type EvaluateEquation(cpx_t input, int equationId);

    float_type EquationDerivativeByParameter(cpx_t input, int equationId, int parameterId);

    float_type EquationDerivativeByParameter(float_type input, int equationId, int parameterId1, int parameterId2);

    float_type EvaluateConstraints();

    float_type ConstraintsDerivativeByParameter(int parameterId);

    vector<cpx_t> GenerateInputs(int size);

    void OutputParameters();

    void OnParameterUpdated();
};
