#pragma once
#include "common.h"
#include "CftData.h"
#include "ConformalBlock.h"

// class for scalar 4-point function.
class Scalar4ptFunction
{
private:
    CftData* cftData;

    // S channel conformal block.
    ConformalBlockScalars* cbS;

    // T channel conformal block.
    ConformalBlockScalars* cbT;

    float_type ConformalBlockDecomposition(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z);

    // calculate derivative of the S-channel CB with respect to scaling dimension of the exchanged operator.
    float_type DerExchangeDelta(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z);

    // calculate derivative of the S-channel 4pt function with respect to scaling dimension of an external operator.
    float_type DerExternalDelta(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z);

    void ClearConformalBlocks();
public:
    int Op1;
    int Op2;
    int Op3;
    int Op4;

    Scalar4ptFunction(CftData* data, int op1, int op2, int op3, int op4);
    ~Scalar4ptFunction();

    void UpdateCftData(CftData* data);

    // returns S-channel conformal block via exchange of given operator.
    float_type ConformalBlock_S(cpx_t z, int op);

    // returns T-channel conformal block via exchange of given operator.
    float_type ConformalBlock_T(cpx_t z, int op);

    // returns difference of S and T channel contributions via exchange of given operator.
    float_type CrossingDiff(cpx_t z, int op);

    // returns value of the 4pt function via S channel calculation.
    float_type EvaluateS(cpx_t z);

    // returns value of the 4pt function via T channel calculation.
    float_type EvaluateT(cpx_t z);

    // returns difference of S and T channel results.
    float_type CrossingDiff(cpx_t z);

    // calculate derivative of the S-channel 4pt function with respect to scaling dimension of a given operator.
    float_type DerDeltaS(cpx_t z, int op);

    // calculate derivative of the T-channel 4pt function with respect to scaling dimension of a given operator.
    float_type DerDeltaT(cpx_t z, int op);

    // derivative of the crossing difference with respect to scaling dimension of a given operator.
    float_type CrossingDerDelta(cpx_t z, int op);

    // derivative of crossing difference associated to the given exchange operator with respect to the given OPE coefficient.
    float_type CrossingDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef, int op);

    // derivative of crossing difference with respect to the given OPE coefficient.
    float_type CrossingDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef);
};
