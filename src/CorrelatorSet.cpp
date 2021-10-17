#include "CorrelatorSet.h"

CorrelatorSet::CorrelatorSet(CftData* data, int numberOfScalarsToBootstrap)
{
    this->cftData = data;

    for (int i = 1; i <= numberOfScalarsToBootstrap; i++) {
        for (int j = i; j <= numberOfScalarsToBootstrap; j++) {
            for (int k = j; k <= numberOfScalarsToBootstrap; k++) {
                for (int l = k; l <= numberOfScalarsToBootstrap; l++) {
                    correlators.push_back(new Scalar4ptFunction(data, i, j, k, l));
                }
            }
        }
    }
}

CorrelatorSet::~CorrelatorSet()
{
    for (uint i = 0; i < correlators.size(); i++) {
        delete correlators[i];
    }

    correlators.clear();
}

float_type CorrelatorSet::Cost(cpx_t z)
{
    float_type ret = .0;

    for (uint i = 0; i < correlators.size(); i++) {
        float_type diff = correlators[i]->CrossingDiff(z);
        ret += diff * diff;
    }

    return ret;
}

float_type CorrelatorSet::CostDerDelta(cpx_t z, int op)
{
    assert(op > 0);

    float_type ret = .0;

    for (uint i = 0; i < correlators.size(); i++) {
        ret += 2 * correlators[i]->CrossingDiff(z) * correlators[i]->CrossingDerDelta(z, op);
    }

    return ret;
}

float_type CorrelatorSet::CostDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef)
{
    float_type ret = .0;

    for (uint i = 0; i < correlators.size(); i++) {
        ret += 2 * correlators[i]->CrossingDiff(z) * correlators[i]->CrossingDerOpeCoefficient(z, coef);
    }

    return ret;
}

void CorrelatorSet::CostDerivatives(cpx_t z, vector<int>& ops, vector<OpeCoefficientKey>& coefs, vector<float_type>& opsDerivatives, vector<float_type>& coefsDerivatives)
{
    vector<float_type> crossingDiffs(correlators.size());
    opsDerivatives.resize(ops.size(), 0.0);
    coefsDerivatives.resize(coefs.size(), 0.0);

    for (uint i = 0; i < correlators.size(); i++) {
        crossingDiffs[i] = correlators[i]->CrossingDiff(z);
    }

    for (uint i = 0; i < ops.size(); i++) {
        for (uint j = 0; j < correlators.size(); j++) {
            opsDerivatives[i] += 2 * crossingDiffs[j] * correlators[j]->CrossingDerDelta(z, ops[i]);
        }
    }

    for (uint i = 0; i < coefs.size(); i++) {
        for (uint j = 0; j < correlators.size(); j++) {
            coefsDerivatives[i] += 2 * crossingDiffs[j] * correlators[j]->CrossingDerOpeCoefficient(z, coefs[i]);
        }
    }
}

void CorrelatorSet::UpdateCftData(CftData* data)
{
    this->cftData = data;
    for (uint i = 0; i < correlators.size(); i++) {
        correlators[i]->UpdateCftData(data);
    }
}
