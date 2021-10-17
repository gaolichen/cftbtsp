#include <fstream>
#include <iomanip>
#include "RandomAlgebricEquations.h"
using namespace std;


RandomAlgebricEquations::RandomAlgebricEquations(int functionNumber, int maxOrder, int tryOrder)
{
    this->functionNumber = functionNumber;
    this->equationNumber = functionNumber;
    this->maxOrder = maxOrder;
    if (tryOrder < 0) {
        this->tryOrder = maxOrder;
    } else {
        this->tryOrder = tryOrder;
    }

    AllocateData();

    for (int i = 0; i <= maxOrder; i++) {
        for (int j = 0; j < functionNumber; j++) {
            if (i == 0) {
                solutions[i][j] = random(0.0, 10.0);
            } else {
                solutions[i][j] = random(-10.0, 10.0);
            }
        }
    }

    for (int i = 0; i < equationNumber; i++) {
        for (int j = 0; j < functionNumber; j++) {
            for (int k = j; k < functionNumber; k++) {
                coefficients[i][j][k] = random(-10.0, 10.0);
                float_type factor = coefficients[i][j][k];
                if (j != k) {
                    coefficients[i][k][j] = factor;
                    factor *= 2;
                }

                for (int h = 0; h <= maxOrder; h++) {
                    for (int g = 0; g <= maxOrder; g++) {
                        constTerms[i][g + h] -= factor * solutions[h][j] * solutions[g][k];
                    }
                }
            }
        }
    }

    GenerateRandomParameters();
}

RandomAlgebricEquations::~RandomAlgebricEquations()
{
}

void RandomAlgebricEquations::GenerateRandomParameters()
{
    for (int i = 0; i <= tryOrder; i++) {
        for (int j = 0; j < functionNumber; j++) {
            if (i == 0) {
                parameters[i][j] = random(0.0, 10.0);
            } else {
                parameters[i][j] = random(-10.0, 10.0);
            }
        }
    }
}

void RandomAlgebricEquations::AllocateData()
{
    this->parameterNumber = functionNumber * (tryOrder + 1);
    this->solutionNumber = functionNumber * (maxOrder + 1);

    solutions.resize(maxOrder + 1, vector<float_type>(functionNumber, .0));
    parameters.resize(tryOrder + 1, vector<float_type>(functionNumber, .0));
    vector<vector<float_type> > coef(functionNumber, vector<float_type>(functionNumber, 0.0));
    coefficients.resize(equationNumber, coef);
    constTerms.resize(equationNumber, vector<float_type>(2 * maxOrder + 1, 0.0));
}

void RandomAlgebricEquations::SaveToFile(string file)
{
    if (file.length() == 0) {
        file = "RandomAlgebricEquations_" + ToString(equationNumber) + "_" + ToString(maxOrder) + "_" + ToString(tryOrder) + "_" + NowToString() + ".txt";
    }

    //std::cout << "Save RandomAlgebricEquations to the file " << file << std::endl;

    ofstream out(file.c_str());
    out << this->equationNumber << ' ' << this->functionNumber << ' ' << maxOrder << ' ' << tryOrder << endl;
    out << "Solutions:" << endl;

    for (int i = 0; i <= maxOrder; i++) {
        SaveVector(solutions[i], out);
    }

    out << "Parameters:" << endl;

    for (uint i = 0; i <= tryOrder; i++) {
        SaveVector(parameters[i], out);
    }

    out << "Coefficients:" << endl;

    for (int i = 0; i < equationNumber; i++) {
        for (int j = 0; j < functionNumber; j++) {
            SaveVector(coefficients[i][j], out);
        }
        out << endl;
    }

    out << "ConstTerms:" << endl;
    for (int i = 0; i < equationNumber; i++) {
        SaveVector(constTerms[i], out);
    }

    out.close();
}

void RandomAlgebricEquations::LoadFromFile(string file)
{
    string buf;
    ifstream in(file.c_str());
    in >> equationNumber >> functionNumber >> maxOrder >> tryOrder;
    AllocateData();
    in >> buf;

    for (int i = 0; i <= maxOrder; i++) {
        LoadVector(solutions[i], in);
    }
    in >> buf;

    for (int i = 0; i <= tryOrder; i++) {
        LoadVector(parameters[i], in);
    }
    in >> buf;

    for (int i = 0; i < equationNumber; i++) {
        for (int j = 0; j < functionNumber; j++) {
            LoadVector(coefficients[i][j], in);
        }
    }
    in >> buf;

    for (int i = 0; i < equationNumber; i++) {
        LoadVector(constTerms[i], in);
    }

    in.close();
}

void RandomAlgebricEquations::SaveVector(vector<float_type>& v, ofstream& out)
{
    for (uint i = 0; i < v.size(); i++) {
        if (i > 0) out << ' ';
        out << setprecision(12) << v[i];
    }
    out << endl;
}

void RandomAlgebricEquations::LoadVector(vector<float_type>& v, ifstream& in)
{
    for (uint i = 0; i < v.size(); i++) {
        in >> v[i];
    }
}

float_type RandomAlgebricEquations::GetParameter(int parameterId)
{
    int order = parameterId / functionNumber;
    return parameters[order][parameterId % functionNumber];
}

void RandomAlgebricEquations::SetParameter(int parameterId, float_type value)
{
    int order = parameterId / functionNumber;
    parameters[order][parameterId % functionNumber] = value;
}

float_type RandomAlgebricEquations::EvaluateEquation(float_type input, int equationId)
{
    return EvaluateEquation(input, equationId, parameters);
}

float_type RandomAlgebricEquations::VerifySolution(float_type input, int equationId)
{
    return EvaluateEquation(input, equationId, solutions);
}

float_type RandomAlgebricEquations::EvaluateEquation(float_type input, int equationId, vector<vector<float_type> >& polynomials)
{
    vector<vector<float_type> >& coef = coefficients[equationId];
    vector<float_type>& constTerm = constTerms[equationId];

    vector<float_type> x;
    CalculateVariables(input, polynomials, x);

    float_type res1 = .0;
    float_type res2 = .0;

    for (int i = 0; i < functionNumber; i++) {
        for (int j = 0; j < functionNumber; j++) {
            res1 += x[i] * x[j] * coef[i][j];
        }
    }

    for (int j = 2 * maxOrder; j >= 0; j--) {
        res2 = res2 * input + constTerm[j];
    }

    return res1 + res2;
}

void RandomAlgebricEquations::CalculateVariables(float_type input, vector<vector<float_type> >& polynomials, vector<float_type>& x)
{
    assert(polynomials.size() > 0);
    assert(functionNumber == polynomials[0].size());
    x.resize(functionNumber, .0);
    for (int i = 0; i < functionNumber; i++) {
        for (int j = polynomials.size() - 1; j >= 0; j--) {
            x[i] = x[i] * input + polynomials[j][i];
        }
    }
}

float_type RandomAlgebricEquations::EquationDerivativeByParameter(float_type input, int equationId, int parameterId)
{
    vector<vector<float_type> >& coef = coefficients[equationId];
    vector<float_type> x;
    CalculateVariables(input, parameters, x);
    int functionId = parameterId % functionNumber;
    int order = parameterId / functionNumber;

    float_type ret = .0;
    for (int i = 0; i < functionNumber; i++) {
        ret += coef[functionId][i] * x[i];
    }

    if (order == 0) {
        return 2 * ret;
    } else {
        return 2 * ret * pow(input, order);
    }
}

float_type RandomAlgebricEquations::EquationDerivativeByParameter(float_type input, int equationId, int parameterId1, int parameterId2)
{
    int var1 = parameterId1 % functionNumber;
    int order1 = parameterId1 / functionNumber;

    int var2 = parameterId2 % functionNumber;
    int order2 = parameterId2 / functionNumber;

    vector<vector<float_type> >& coef = coefficients[equationId];

    return 2 * coef[var1][var2] * pow(input, order1 + order2);
}

float_type RandomAlgebricEquations::EvaluateConstraints()
{
    float_type ret = .0;
    for (int i = 0; i < functionNumber; i++) {
        ret += exp(-parameters[0][i] / minGap);
//        if (parameters[0][i] < 0) ret -= parameters[0][i];
    }

    return ret;
}

float_type RandomAlgebricEquations::ConstraintsDerivativeByParameter(int parameterId)
{
    if (parameterId >= functionNumber) return .0;
    int functionId = parameterId % functionNumber;
    return -exp(-parameters[0][functionId] / minGap) / minGap;
//    if (parameters[0][functionId] < 0) return -1.0;
//    return .0;
}

vector<float_type> RandomAlgebricEquations::GenerateInputs(int size)
{
    vector<float_type> inputs(size, .0);
    for (int i = 0; i < size; i++) {
        inputs[i] = random(-.2, .2);
    }

    return inputs;
}

void RandomAlgebricEquations::OutputParameters()
{
    for (int i = 0; i <= tryOrder; i++) {
        std::cout << "function " << i << "=" << setprecision(12) << parameters[i] << std::endl;
    }
}

void RandomAlgebricEquations::OutputSolutions()
{
    std::cout << "Solution of the equations:" << std::endl;
    for (int i = 0; i <= maxOrder; i++) {
        std::cout << "function " << i << "=" << setprecision(12) << solutions[i] << std::endl;
    }
}

float_type RandomAlgebricEquations::GetSolution(int parameterId)
{
    int functionId = parameterId % functionNumber;
    return solutions[parameterId / functionNumber][functionId];
}

float_type RandomAlgebricEquations::GetCoefficient(int equationId, int var1, int var2)
{
    return coefficients[equationId][var1][var2];
}

float_type RandomAlgebricEquations::GetConstTerm(int equationId, int order)
{
    return constTerms[equationId][order];
}
