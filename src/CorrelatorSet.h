#pragma once
#include <vector>
#include "common.h"
#include "CftData.h"
#include "Scalar4ptFunction.h"
using namespace std;

// a set of correlators used for conformal bootstrap.
class CorrelatorSet
{
private:
    CftData* cftData;
    vector<Scalar4ptFunction*> correlators;
public:
    CorrelatorSet(CftData* data, int numberOfScalarsToBootstrap);
    ~CorrelatorSet();

    float_type Cost(cpx_t z);

    float_type CostDerDelta(cpx_t z, int op);
    
    float_type CostDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef);

    void CostDerivatives(cpx_t z, vector<int>& ops, vector<OpeCoefficientKey>& coefs, vector<float_type>& opsDerivatives, vector<float_type>& coefsDerivatives);

    void UpdateCftData(CftData* data);
};
