#include "CrossingEquations.h"

CrossingEquations::CrossingEquations(CftData* data, int numberOfScalarsToBootstrap)
{
    this->cftData = data;
    unitaryBoundaryMinGap = 0.05;
    operatorMinGap = 0.1;

    // number of scale dimensions as parameter.
    this->operatorNumber = this->cftData->MaxPrimaryId() - 1;

    for (int i = 1; i <= numberOfScalarsToBootstrap; i++) {
        for (int j = i; j <= numberOfScalarsToBootstrap; j++) {
            for (int k = j; k <= this->cftData->MaxPrimaryId(); k++) {
                if (k == this->cftData->StressTensorId) continue;
                parameter2OpeCoefficientKey.push_back(OpeCoefficientKey(i, j, k));
            }
        }
    }

    // create Scalar4ptFunctions.
    for (int i = 1; i <= numberOfScalarsToBootstrap; i++) {
        for (int j = i; j <= numberOfScalarsToBootstrap; j++) {
            for (int k = j; k <= numberOfScalarsToBootstrap; k++) {
                for (int l = k; l <= numberOfScalarsToBootstrap; l++) {
                    correlators.push_back(new Scalar4ptFunction(data, i, j, k, l));
                    // when i, j, k, l are different, there is another correlator that 
                    // can produce independent crossing constraint.
                    if (i < j && j < k && k < l) {
                        correlators.push_back(new Scalar4ptFunction(data, i, j, l, k));
                    }
                }
            }
        }
    }
}

CrossingEquations::~CrossingEquations()
{
    for (uint i = 0; i < correlators.size(); i++) {
        delete correlators[i];
    }

    correlators.clear();
}

void CrossingEquations::GenerateRandomParameters()
{
    this->cftData->GenerateRandomData();
    OnParameterUpdated();
}

void CrossingEquations::SaveToFile(string file)
{
    this->cftData->Save(file);
}

int CrossingEquations::EquationNumber()
{
    return this->correlators.size();
}

int CrossingEquations::ParameterNumber()
{
    return this->operatorNumber + parameter2OpeCoefficientKey.size();
}

float_type CrossingEquations::GetParameter(int parameterId)
{
    // if the parameter is a scaling dimension.
    if (parameterId < this->operatorNumber) {
        // the scalar dimension of identity operator and stresstensor is not parameter.
        if (parameterId < this->cftData->StressTensorId) {
            return this->cftData->GetPrimaryDim(parameterId + 1);
        } else {
            return this->cftData->GetPrimaryDim(parameterId + 2);
        }
    } else {
        // in this case the parameter is an OPE coefficients of scalar operators.
        return this->cftData->GetOpeCoefficient(parameter2OpeCoefficientKey[parameterId - this->operatorNumber]);
    }
}

void CrossingEquations::SetParameter(int parameterId, float_type value)
{
    // if the parameter is a scaling dimension.
    if (parameterId < this->operatorNumber) {
        // the scalar dimension of identity operator and stresstensor is not parameter.
        if (parameterId < this->cftData->StressTensorId) {
            this->cftData->SetPrimaryDim(parameterId + 1, value);
        } else {
            this->cftData->SetPrimaryDim(parameterId + 2, value);
        }
    } else {
        // in this case the parameter is an OPE coefficients of scalar operators.
        this->cftData->SetOpeCoefficient(parameter2OpeCoefficientKey[parameterId-operatorNumber], value);
    }
}

float_type CrossingEquations::EvaluateEquation(cpx_t input, int equationId)
{
    return correlators[equationId]->CrossingDiff(input);
}

float_type CrossingEquations::EquationDerivativeByParameter(cpx_t input, int equationId, int parameterId)
{
    // if the parameter is a scaling dimension.
    if (parameterId < this->operatorNumber) {
        // the scalar dimension of identity operator and stresstensor is not parameter.
        if (parameterId < this->cftData->StressTensorId) {
            return correlators[equationId]->CrossingDerDelta(input, parameterId + 1);
        } else {
            return correlators[equationId]->CrossingDerDelta(input, parameterId + 2);
        }
    } else {
        // in this case the parameter is an OPE coefficients of scalar operators.
        return correlators[equationId]->CrossingDerOpeCoefficient(input, parameter2OpeCoefficientKey[parameterId-operatorNumber]);
    }
}

float_type CrossingEquations::EquationDerivativeByParameter(float_type input, int equationId, int parameterId1, int parameterId2)
{
    // TODO: return 0 temperarily.
    return .0;
}

float_type CrossingEquations::EvaluateConstraints()
{
    float_type ret = .0;
    float_type lastDim = .0, dim;
    int lastSpin = -1;

    for (int i = 1; i <= cftData->MaxPrimaryId(); i++) {
        dim = cftData->GetPrimaryDim(i);
        if (cftData->GetPrimarySpin(i) != lastSpin) {
            lastSpin = cftData->GetPrimarySpin(i);
            if (lastSpin == 0) {
                ret += exp((cftData->D / 2.0 - 1 - dim) / unitaryBoundaryMinGap);
            } else {
                ret += exp((lastSpin + cftData->D - 2 - dim) / unitaryBoundaryMinGap);
            }
        } else {
            ret += exp((lastDim - dim) / operatorMinGap);
        }

        lastDim = dim;
    }

    return ret;
}

float_type CrossingEquations::ConstraintsDerivativeByParameter(int parameterId)
{
    if (parameterId >= this->operatorNumber) return .0;

    int opId = parameterId + 1;
    // the scalar dimension of identity operator and stresstensor is not parameter.
    if (parameterId >= this->cftData->StressTensorId) {
        opId = parameterId + 2;
    }

    int l = cftData->GetPrimarySpin(opId);
    float_type dim = cftData->GetPrimaryDim(opId);
    float_type ret = .0;

    // the dimension of operator should greater than the dimension of next operator with same spin.
    if (opId > 1 && cftData->GetPrimarySpin(opId - 1) == l) {
        float_type lastDim = cftData->GetPrimaryDim(opId - 1);
        ret -= exp((lastDim - dim) / operatorMinGap) / operatorMinGap;
    }

    // the dimension of operator should less than the dimension of next operator with same spin.
    if (opId < cftData->MaxPrimaryId() && cftData->GetPrimarySpin(opId + 1) == l) {
        float_type nextDim = cftData->GetPrimaryDim(opId + 1);
        ret += exp((dim - nextDim) / operatorMinGap) / operatorMinGap;
    }

    // unitary boundary. 
    if (opId == 1 || cftData->GetPrimarySpin(opId - 1) != l) {
        if (l == 0) {
            ret -= exp((cftData->D / 2.0 - 1 - dim) / unitaryBoundaryMinGap) / unitaryBoundaryMinGap;
        } else {
            ret -= exp((l + cftData->D - 2 - dim) / unitaryBoundaryMinGap) / unitaryBoundaryMinGap;
        }
    }

    return ret;
}

vector<cpx_t> CrossingEquations::GenerateInputs(int size)
{
    vector<cpx_t> ret;
    for (int i = 0; i < size; i++) {
        cpx_t z;
        while (true) {
            z = RandomComplex(0.1);
            if (Z2r(z + 0.5) < 0.2 && Z2r(-z + 0.5) < 0.2) {
                break;
            }
        }
        ret.push_back(z + 0.5);
        ret.push_back(-z + 0.5);
    }

    return ret;
}

void CrossingEquations::OutputParameters()
{
    vector<int> scalarIds;
    for (int i = 1; i <= cftData->MaxScalarId(); i++) scalarIds.push_back(i);
    cftData->Output(scalarIds);
}

void CrossingEquations::OnParameterUpdated()
{
    for (uint i = 0; i < correlators.size(); i++) {
        correlators[i]->UpdateCftData(this->cftData);
    }
}

