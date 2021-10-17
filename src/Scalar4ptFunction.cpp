#define CROSSING_EQUATION_NORMALIZATION 1

#include <cmath>
#include "Scalar4ptFunction.h"
#include "ConformalBlock.h"
using namespace std;

double Z2r(cpx_t z)
{
    cpx_t rho = z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
    return abs(rho);
}

Scalar4ptFunction::Scalar4ptFunction(CftData* data, int op1, int op2, int op3, int op4) 
: Op1(op1), Op2(op2), Op3(op3), Op4(op4)
{
    cftData = 0; 
    cbS = 0; 
    cbT = 0;
    UpdateCftData(data);
}

Scalar4ptFunction::~Scalar4ptFunction()
{
    ClearConformalBlocks();
}

void Scalar4ptFunction::ClearConformalBlocks()
{
    if (cbS != 0) {
        delete cbS;
        cbS = 0;
    }

    if (cbT != 0) {
        delete cbT;
        cbT = 0;
    }
}

void Scalar4ptFunction::UpdateCftData(CftData* data) 
{
    this->cftData = data;
    float_type dlt1 = cftData->GetPrimaryDim(Op1);
    float_type dlt2 = cftData->GetPrimaryDim(Op2);
    float_type dlt3 = cftData->GetPrimaryDim(Op3);
    float_type dlt4 = cftData->GetPrimaryDim(Op4);

    ClearConformalBlocks();
    cbS = new ConformalBlockScalars(cftData->D, dlt1 - dlt2, dlt3 - dlt4);
    cbT = new ConformalBlockScalars(cftData->D, dlt3 - dlt2, dlt1 - dlt4);
}

float_type Scalar4ptFunction::ConformalBlockDecomposition(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z)
{
    cpx_t rho = z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
    float_type r = abs(rho);
    float_type eta = cos(arg(rho));

    float_type delta = cftData->GetPrimaryDim(op);
    int l = cftData->GetPrimarySpin(op);

#if CROSSING_EQUATION_NORMALIZATION == 1
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);
    return pow(u,  -(dltj + dlti) / 2) * cb->evaluate(delta, l, r, eta);
#elif CROSSING_EQUATION_NORMALIZATION == 2
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);
    float_type v = abs(-z + 1.0) * abs(-z + 1.0);

    return pow(u,  -dlti / 2) * pow(v,  dltj / 2) * cb->evaluate(delta, l, r, eta);
#else
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dltk = cftData->GetPrimaryDim(opk);

    float_type v = abs(-z + 1.0) * abs(-z + 1.0);
    return pow(v,  (dltj + dltk) / 2) * cb->evaluate(delta, l, r, eta);
#endif
}

float_type Scalar4ptFunction::ConformalBlock_S(cpx_t z, int op)
{
    return ConformalBlockDecomposition(cbS, Op1, Op2, Op3, Op4, op, z);
}

float_type Scalar4ptFunction::ConformalBlock_T(cpx_t z, int op)
{
    return ConformalBlockDecomposition(cbT, Op3, Op2, Op1, Op4, op, -z + 1.0);
}

float_type Scalar4ptFunction::CrossingDiff(cpx_t z, int op)
{
    float_type coef1 = cftData->GetOpeCoefficient(Op1, Op2, op) * cftData->GetOpeCoefficient(Op3, Op4, op);
    float_type coef2 = cftData->GetOpeCoefficient(Op3, Op2, op) * cftData->GetOpeCoefficient(Op1, Op4, op);
    return coef1 * ConformalBlock_S(z, op) - coef2 * ConformalBlock_T(z, op);
}

float_type Scalar4ptFunction::CrossingDiff(cpx_t z)
{
    float_type ret = 0.0;
    for (int op = 0; op <= cftData->MaxPrimaryId(); op++) {
        ret += CrossingDiff(z, op);
    }

    return ret;
}

float_type Scalar4ptFunction::EvaluateS(cpx_t z)
{
    float_type ret = 0.0;
    for (int op = 0; op <= cftData->MaxPrimaryId(); op++) {
        float_type coef = cftData->GetOpeCoefficient(Op1, Op2, op) * cftData->GetOpeCoefficient(Op3, Op4, op);
        ret += coef * ConformalBlock_S(z, op);
    }

    return ret;
}

float_type Scalar4ptFunction::EvaluateT(cpx_t z)
{
    float_type ret = 0.0;
    for (int op = 0; op <= cftData->MaxPrimaryId(); op++) {
        float_type coef = cftData->GetOpeCoefficient(Op3, Op2, op) * cftData->GetOpeCoefficient(Op1, Op4, op);
        ret += coef * ConformalBlock_T(z, op);
    }

    return ret;
}

float_type Scalar4ptFunction::DerExchangeDelta(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z)
{
    cpx_t rho = z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
    float_type r = abs(rho);
    float_type eta = cos(arg(rho));

    float_type delta = cftData->GetPrimaryDim(op);
    int l = cftData->GetPrimarySpin(op);

    float_type ret = cftData->GetOpeCoefficient(opi, opj, op) * cftData->GetOpeCoefficient(opk, opl, op);

#if CROSSING_EQUATION_NORMALIZATION == 1
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);
    ret *= pow(u,  -(dltj + dlti) / 2) * cb->dDelta(delta, l, r, eta);
#elif CROSSING_EQUATION_NORMALIZATION == 2
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);
    float_type v = abs(-z + 1.0) * abs(-z + 1.0);
    ret *= pow(u,  -dlti / 2) * pow(v,  dltj / 2) * cb->dDelta(delta, l, r, eta);
#else
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dltk = cftData->GetPrimaryDim(opk);

    float_type v = abs(-z + 1.0) * abs(-z + 1.0);
    ret *= pow(v,  (dltj + dltk) / 2) * cb->dDelta(delta, l, r, eta);
#endif
    return ret;
}

float_type Scalar4ptFunction::DerExternalDelta(ConformalBlockScalars* cb, int opi, int opj, int opk, int opl, int op, cpx_t z)
{
    int a = 0, b = 0, c = 0;
    int d = 0;

#if CROSSING_EQUATION_NORMALIZATION == 1
    if (opi == op) {
        a++; c--;
    }
    if (opl == op) b--;
    if (opj == op) { 
        a--; c--;
    }
    if (opk == op) b++;
#elif CROSSING_EQUATION_NORMALIZATION == 2
    if (opi == op) {
        a++; c--;
    }
    if (opj == op) { 
        a--; d++;
    }
    if (opk == op) b++;
    if (opl == op) b--;
#else
    if (opi == op) a++;
    if (opl == op) b--;
    if (opj == op) { 
        c++; a--;
    }
    if (opk == op) {
        c++; b++;
    }
#endif

    if (a == 0 && b == 0 && c == 0) return .0;

    cpx_t rho = z / (-z + 2.0 + sqrt(-z + 1.0) * 2.0);
    float_type r = abs(rho);
    float_type eta = cos(arg(rho));

#if CROSSING_EQUATION_NORMALIZATION == 1
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);

    // the common factor
    float_type factor = pow(u, -(dltj + dlti)/2);
#elif CROSSING_EQUATION_NORMALIZATION == 2
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dlti = cftData->GetPrimaryDim(opi);

    float_type u = abs(z) * abs(z);
    float_type v = abs(-z + 1.0) * abs(-z + 1.0);

    // the common factor
    float_type factor = pow(u, -dlti/2) * pow(v, dltj/2);
#else
    float_type dltj = cftData->GetPrimaryDim(opj);
    float_type dltk = cftData->GetPrimaryDim(opk);

    float_type v = abs(-z + 1.0) * abs(-z + 1.0);

    // the common factor
    float_type factor = pow(v, (dltj + dltk)/2);
#endif

    float_type ret = 0.0;

    // variable a: makes contribution v^(dltk/2 + dltj/2) * {derivative of CB with respect to dlt12}
    if (a != 0) {
        for (int id = 0; id <= cftData->MaxPrimaryId(); id++) {
            float_type delta = cftData->GetPrimaryDim(id);
            int l = cftData->GetPrimarySpin(id);
            float_type coef = cftData->GetOpeCoefficient(opi, opj, id) * cftData->GetOpeCoefficient(opk, opl, id);
            ret += a * coef * cb->dDelta12(delta, l, r, eta);
        }
    }

    // variable b: makes contribution v^(dltk/2 + dltj/2) * {derivative of CB with respect to dlt34}
    if (b != 0) {
        for (int id = 0; id <= cftData->MaxPrimaryId(); id++) {
            float_type delta = cftData->GetPrimaryDim(id);
            int l = cftData->GetPrimarySpin(id);
            float_type coef = cftData->GetOpeCoefficient(opi, opj, id) * cftData->GetOpeCoefficient(opk, opl, id);
            ret += b * coef * cb->dDelta34(delta, l, r, eta);
        }
    }

    // variable a: makes contribution v^(dltk/2 + dltj/2) * log(v) * 0.5 * {value of CB}
    if (c != 0 || d != 0) {
        float_type totCB = .0;
        for (int id = 0; id <= cftData->MaxPrimaryId(); id++) {
            float_type delta = cftData->GetPrimaryDim(id);
            int l = cftData->GetPrimarySpin(id);
            float_type coef = cftData->GetOpeCoefficient(opi, opj, id) * cftData->GetOpeCoefficient(opk, opl, id);
            totCB += coef * cb->evaluate(delta, l, r, eta);
        }
#if CROSSING_EQUATION_NORMALIZATION == 1
        ret += c * log(u) * 0.5 * totCB;
#elif CROSSING_EQUATION_NORMALIZATION == 2
        ret += (c * log(u) + d * log(v)) * 0.5 * totCB;
#else
        ret += c * log(v) * 0.5 * totCB;
#endif
    }

    return factor * ret;
}

float_type Scalar4ptFunction::DerDeltaS(cpx_t z, int op)
{
    return DerExchangeDelta(cbS, Op1, Op2, Op3, Op4, op, z) + DerExternalDelta(cbS, Op1, Op2, Op3, Op4, op, z);
}

float_type Scalar4ptFunction::DerDeltaT(cpx_t z, int op)
{
    return DerExchangeDelta(cbT, Op3, Op2, Op1, Op4, op, -z + 1.0) + DerExternalDelta(cbT, Op3, Op2, Op1, Op4, op, -z + 1.0);
}

float_type Scalar4ptFunction::CrossingDerDelta(cpx_t z, int op)
{
    return DerDeltaS(z, op) - DerDeltaT(z, op);
}

float_type Scalar4ptFunction::CrossingDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef, int op)
{
    float_type ret1 = 0.0;
    if (coef == OpeCoefficientKey(Op1, Op2, op)) {
        ret1 += cftData->GetOpeCoefficient(Op3, Op4, op);
    }

    if (coef == OpeCoefficientKey(Op3, Op4, op)) {
        ret1 += cftData->GetOpeCoefficient(Op1, Op2, op);
    }

    if (abs(ret1) > EPS) {
        ret1 *= ConformalBlock_S(z, op);
    }

    float_type ret2 = 0.0;
    if (coef == OpeCoefficientKey(Op3, Op2, op)) {
        ret2 += cftData->GetOpeCoefficient(Op1, Op4, op);
    }

    if (coef == OpeCoefficientKey(Op1, Op4, op)) {
        ret2 += cftData->GetOpeCoefficient(Op3, Op2, op);
    }

    if (abs(ret2) > EPS) {
        ret2 *= ConformalBlock_T(z, op);
    }

    return ret1 - ret2;
}

float_type Scalar4ptFunction::CrossingDerOpeCoefficient(cpx_t z, OpeCoefficientKey& coef)
{
    float_type ret = 0.0;
    for (int op = 0; op <= cftData->MaxPrimaryId(); op++) {
        ret += CrossingDerOpeCoefficient(z, coef, op);
    }

    return ret;
}
