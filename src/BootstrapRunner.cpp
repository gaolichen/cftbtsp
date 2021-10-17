#include "BootstrapRunner.h"
#include <vector>
#include <boost/timer.hpp>
using namespace std;

BootstrapRunner::BootstrapRunner(BoostrapConfig* bootstrapConfig, CftConfig cftConfig)
{
    this->bootstrapConfig = bootstrapConfig;
    this->cftData = new CftData(cftConfig);
    this->correlatorSet = new CorrelatorSet(cftData, bootstrapConfig->NumberOfScalarsToBootstrap);
}

BootstrapRunner::~BootstrapRunner()
{
    if (this->cftData != 0) {
        delete cftData;
    }

    if (this->correlatorSet != 0) {
        delete this->correlatorSet;
    }
}

float_type BootstrapRunner::ConstraintCost(CftData* cftData)
{
    float_type ret = .0;
    float_type lastDim = .0, dim;
    int lastSpin = -1;

    for (uint i = 1; i <= cftData->MaxPrimaryId(); i++) {
        dim = cftData->GetPrimaryDim(i);
        if (cftData->GetPrimarySpin(i) != lastSpin) {
            lastSpin = cftData->GetPrimarySpin(i);
            if (lastSpin == 0) {
                if (dim < cftData->D / 2.0 - 1) ret += (cftData->D / 2.0 - 1 - dim);
            } else {
                if (dim < lastSpin + cftData->D - 2) ret += (lastSpin + cftData->D - 2 - dim);
            }
        } else if (dim < lastDim) {
            ret += (lastDim - dim);
        }

        lastDim = dim;
    }

    return ret * bootstrapConfig->ConstraintFactor;
}

void BootstrapRunner::ConstraintCostDerivative(CftData* cftData, vector<int>& ops, vector<float_type>& opsDerivatives)
{
    assert(ops.size() == opsDerivatives.size());

    float_type dim, lastDim, nextDim;

    for (uint i = 0; i < ops.size(); i++) {
        int l = cftData->GetPrimarySpin(ops[i]);
        dim = cftData->GetPrimaryDim(ops[i]);

        if (ops[i] > 1 && cftData->GetPrimarySpin(ops[i] - 1) == l) {
            lastDim = cftData->GetPrimaryDim(ops[i] - 1);
            if (dim < lastDim) {
                opsDerivatives[i] -= bootstrapConfig->ConstraintFactor;
            }
        }

        if (ops[i] < cftData->MaxPrimaryId() && cftData->GetPrimarySpin(ops[i] + 1) == l) {
            nextDim = cftData->GetPrimaryDim(ops[i] + 1);
            if (dim > nextDim) {
                opsDerivatives[i] += bootstrapConfig->ConstraintFactor;
            }
        }

        if (ops[i] == 1 || cftData->GetPrimarySpin(ops[i] - 1) != l) {
            if (l == 0 && dim <= cftData->D / 2.0 - 1) {
                opsDerivatives[i] -= bootstrapConfig->ConstraintFactor;
            }

            if (l > 0 && dim <= cftData->D + l - 2.0) {
                opsDerivatives[i] -= bootstrapConfig->ConstraintFactor;
            }
        }
    }
}

void BootstrapRunner::Run()
{
    vector<OpeCoefficientKey> coefs;
    vector<int> ops;

    for (uint op = 1; op <= cftData->MaxPrimaryId(); op++) {
        if (op != cftData->StressTensorId) {
            ops.push_back(op);
        }
    }

    for (uint op1 = 1; op1 <= bootstrapConfig->NumberOfScalarsToBootstrap; op1++) {
        for (uint op2 = op1; op2 <= bootstrapConfig->NumberOfScalarsToBootstrap; op2++) {
            for (uint op3 = op2; op3 <= cftData->MaxPrimaryId(); op3++) {
                if (op3 != cftData->StressTensorId) {
                    coefs.push_back(OpeCoefficientKey(op1, op2, op3));
                }
            }
        }
    }

    vector<cpx_t> zs;
    zs.reserve(bootstrapConfig->SamplesEachStep);
    for (uint i = 0; i < bootstrapConfig->SamplesEachStep; i++) {
        float_type dr = random(0.0, 0.2);
        float_type theta = random(0.0, 2 * acos(-1.0));
        zs.push_back(cpx_t(dr * cos(theta), dr * sin(theta)) + 0.5);
    }

    vector<float_type> opsDerivatives(ops.size(), .0);
    vector<float_type> opsDerivativeTot(ops.size(), .0);
    vector<float_type> coefsDerivatives(coefs.size(), .0);
    vector<float_type> coefsDerivativeTot(coefs.size(), .0);

    float_type lastCost = std::numeric_limits<float_type>::max();
    vector<int> scalarIds;
    for (uint i = 1; i <= cftData->MaxScalarId(); i++) scalarIds.push_back(i);

    boost::timer stopwatch;

    while (true) {
        std::cout << "Time elapsed: " << stopwatch.elapsed() << " seconds." << endl;
        cftData->Output(scalarIds);

        opsDerivativeTot.resize(ops.size(), .0);
        coefsDerivativeTot.resize(coefs.size(), .0);

        float_type cost = .0;

        for (uint i = 0; i < zs.size(); i++) {
            this->correlatorSet->CostDerivatives(zs[i], ops, coefs, opsDerivatives, coefsDerivatives);
            cost += this->correlatorSet->Cost(zs[i]);

            for (uint j = 0; j < opsDerivatives.size(); j++) {
                opsDerivativeTot[j] += opsDerivatives[j];
            }

            for (uint j = 0; j < coefsDerivatives.size(); j++) {
                coefsDerivativeTot[j] += coefsDerivatives[j];
            }
        }

        cost /= zs.size();
        std::cout << "pure cost = " << cost << "\t";
        cost += ConstraintCost(cftData);
        std::cout << "total cost = " << cost << std::endl;
        
        if (cost < bootstrapConfig->Accuracy) break;

        float_type gradNorm = .0;
        for (uint j = 0; j < coefsDerivativeTot.size(); j++) {
            coefsDerivativeTot[j] /= zs.size();
            gradNorm += coefsDerivativeTot[j] * coefsDerivativeTot[j];
        }

        std::cout << "Coefs Derivative: " << coefsDerivativeTot << std::endl;

        for (uint j = 0; j < opsDerivativeTot.size(); j++) {
            opsDerivativeTot[j] /= zs.size();
        }

        std::cout << "Ops Derivative: " << opsDerivativeTot << std::endl;

        ConstraintCostDerivative(cftData, ops, opsDerivativeTot);

        std::cout << "Ops Derivative with constraints: " << opsDerivativeTot << std::endl;

        for (uint j = 0; j < opsDerivativeTot.size(); j++) {
            gradNorm += opsDerivativeTot[j] * opsDerivativeTot[j];
        }

        float_type factor = gradNorm / ((ops.size() + coefs.size()));
        factor = sqrt(factor);

        for (uint j = 0; j < coefsDerivativeTot.size(); j++) {
            coefsDerivativeTot[j] /= factor;
        }

        for (uint j = 0; j < opsDerivativeTot.size(); j++) {
            opsDerivativeTot[j] /= factor;
        }

        std::cout << "Coefs Derivative normalized: " << coefsDerivativeTot << std::endl;
        std::cout << "Ops Derivative normalized: " << opsDerivativeTot << std::endl;

        for (uint j = 0; j < opsDerivativeTot.size(); j++) {
            float_type dim = cftData->GetPrimaryDim(ops[j]);
            cftData->SetPrimaryDim(ops[j], dim - opsDerivativeTot[j] * bootstrapConfig->Step);
        }

        for (uint j = 0; j < coefsDerivativeTot.size(); j++) {
            float_type coef = cftData->GetOpeCoefficient(coefs[j]);
            cftData->SetOpeCoefficient(coefs[j], coef - coefsDerivativeTot[j] * bootstrapConfig->Step);
        }

        this->correlatorSet->UpdateCftData(cftData);
        lastCost = cost;
    }
}


