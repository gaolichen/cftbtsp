#pragma once
#include "common.h"
#include "CftData.h"
#include "CorrelatorSet.h"


struct BoostrapConfig
{
    int NumberOfScalarsToBootstrap;

    int SamplesEachStep;

    float_type Step;

    float_type Accuracy;

    float_type ConstraintFactor;

    BoostrapConfig()
    {
        NumberOfScalarsToBootstrap = 5;
        SamplesEachStep = 10;
        Step = 0.001;
        Accuracy = 0.01;
        ConstraintFactor = 20.0;
    }
};

class BootstrapRunner
{
private:
    CftData* cftData;
    BoostrapConfig* bootstrapConfig;
    CorrelatorSet* correlatorSet;
public:
    BootstrapRunner(BoostrapConfig* bootstrapConfig, CftConfig cftConfig);
    ~BootstrapRunner();

    float_type ConstraintCost(CftData* cftData);

    void ConstraintCostDerivative(CftData* cftData, vector<int>& ops, vector<float_type>& opsDerivatives);

    void Run();
};
