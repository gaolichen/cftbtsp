#include <iostream>
#include <fstream>
#include <cmath>
#include "CftData.h"
#include "common.h"
#include <boost/math/special_functions/gamma.hpp>

std::ostream& operator<< (std::ostream& out, const CftConfig& config)
{
    out << "{D=" << config.D << ", OperatorNumbers=" << config.OperatorNumbers << '}';
    return out;
}

CftData::CftData(CftConfig config)
{
    this->D = config.D;
    this->opsNumbers = config.OperatorNumbers;
    GenerateRandomData();
}

void CftData::GenerateRandomData()
{
    // clear the data.
    ops.clear();
    opeCoefs.clear();

    // add Identity operator.
    ops.push_back(PrimaryInfo(0, 0, 0));

    int id = 0;
    for (uint i = 0; i < this->opsNumbers.size(); i++) {
        float_type lastDim = .0;
        float_type gap;

        // the unitary bound.
        if (i == 0) {
            gap = this->D / 2.0 - 1;
        } else {
            gap = i + this->D - 2.0;
        }

        for (int j = 0; j < this->opsNumbers[i]; j++) {
            double dim;
            id++;
            // for stress tensor
            if (i == 2 && j == 0) {
                dim = (double)this->D;
            } else {
                // randomly set scaling dimsion for primaries except for that
                // scaling dimesions increase with id.
                dim = lastDim + gap + random(0, .2);
            }

            gap = .8;

            ops.push_back(PrimaryInfo(id, dim, i));
            lastDim = dim;
        }
    }

    // randomly generate OPE coefficients for three scalars operators and two scalars and one spinning operators..
    for (int i = 1; i <= this->opsNumbers[0]; i++) {
        for (int j = i; j <= this->opsNumbers[0]; j++) {
            for (uint k = j; k <= this->ops.size(); k++) {
                SetOpeCoefficient(this->ops[i].Id, this->ops[j].Id, this->ops[k].Id, random(-5.0, 5.0));
            }
        }
    }

    SetFixedCftData();
}

void CftData::SetFixedCftData()
{
    maxPrimaryId = Sum(opsNumbers);
    this->StressTensorId = opsNumbers[0] + opsNumbers[1] + 1;

    // normalize OPE coefficients of two identitical scalars to identity operator or stress tensor.
    // for identity operator the OPE coefficients are 1.
    // for stress tensor, it is -d * Dim/(d-1)  * 1/S_d, where S_d is is the volumn of the unit sphere S^(d-1).
    float_type Pi = acos(-1.0);
    for (int i = 0; i <= opsNumbers[0]; i++) {
        SetOpeCoefficient(i, i, 0, 1.0);
        if (i > 0) {
            double coef = -D * GetPrimaryDim(i) / (D - 1);
            double sd = 2 * pow(Pi, D/2.0) / boost::math::tgamma<float_type>(D/2.0);
            SetOpeCoefficient(i, i, StressTensorId, coef / sd);
        }

        for (int j = i + 1; j <= opsNumbers[0]; j++) {
            SetOpeCoefficient(i, j, 0, 0.0);
            SetOpeCoefficient(i, j, StressTensorId, 0.0);
        }
    }
}

PrimaryInfo CftData::GetPrimaryInfo(int id) const
{
    return this->ops[id];
}

void CftData::SetPrimaryDim(int id, double dim)
{
    this->ops[id].Dim = dim;
}

int CftData::PrimaryNumber(int spin) const
{
    return this->opsNumbers[spin];
}

float_type CftData::GetOpeCoefficient(int id1, int id2, int id3)
{
    return this->opeCoefs[OpeCoefficientKey(id1, id2, id3)];
}

float_type CftData::GetOpeCoefficient(OpeCoefficientKey& key)
{
    return this->opeCoefs[key];
}

void CftData::SetOpeCoefficient(int id1, int id2, int id3, float_type coef)
{
    this->opeCoefs[OpeCoefficientKey(id1, id2, id3)] = coef;
}

void CftData::SetOpeCoefficient(OpeCoefficientKey& key, float_type coef)
{
    this->opeCoefs[key] = coef;
}

void CftData::Output(vector<int> ops)
{
    std::cout << "Dimensions of operators" << std::endl;
    for (uint i = 0; i < ops.size(); i++) {
        std::cout << ops[i] << ":\t" << GetPrimaryDim(ops[i]) << std::endl;
    }

    std::cout << std::endl << "OPE coefficients: " << std::endl;

    int count = 0;

    for (uint i = 0; i < ops.size(); i++) {
        for (uint j = i; j < ops.size(); j++) {
            for (uint k = j; k < ops.size(); k++) {
                std::cout << '(' << ops[i] << ',' << ops[j] << ',' << ops[k] << "): ";
                std::cout << GetOpeCoefficient(ops[i], ops[j], ops[k]) << "\t";
                count++;
                if (count % 5 == 0) std::cout << endl;
            }
        }
    }
    std::cout << std::endl;
}

void CftData::Save(string file)
{
    ofstream out(file);
    out << this->D << ' ' << this->opsNumbers.size() << std::endl;

    out << "OperatorNubmers:" << std::endl;
    for (uint i = 0; i < this->opsNumbers.size(); i++) {
        if (i > 0) out << ' ';
        out << opsNumbers[i];
    }
    out << std::endl;
    
    out << "ScalingDimensions:" << std::endl;
    for (uint i = 1; i < ops.size(); i++) {
        if (i > 1) out << ' ';
        out << setprecision(15) << ops[i].Dim;
    }
    out << std::endl;

    out << "OpeCoefficients:" << std::endl;
    for (int i = 1; i <= this->opsNumbers[0]; i++) {
        for (int j = i; j <= this->opsNumbers[0]; j++) {
            for (uint k = j; k < this->ops.size(); k++) {
                if (k == StressTensorId) continue;
                out << i << ' ' << j << ' ' << k << ' ' << setprecision(15) << GetOpeCoefficient(i, j, k) << std::endl;
            }
        }
    }

    out.close();
}

void CftData::LoadFromFile(string file)
{
    ifstream in(file);
    int spinNumber;
    string buf;
    in >> this->D >> spinNumber;
    this->opsNumbers.resize(spinNumber, 0);

    in >> buf;
    for (int i = 0; i < spinNumber; i++) {
        in >> opsNumbers[i];
    }
    
    in >> buf;

    int id = 0;
    float_type dim;
    ops.clear();
    ops.push_back(PrimaryInfo(id, 0.0, 0));
    for (int spin = 0; spin < spinNumber; spin++) {
        for (int i = 0; i < opsNumbers[spin]; i++) {
            in >> dim;
            ops.push_back(PrimaryInfo(++id, dim, spin));
        }
    }

    in >> buf;

    int op1, op2, op3;
    float_type coef;
    while(!in.eof()) {
        in >> op1 >> op2 >> op3 >> coef;
        if (!in.good()) break;
        SetOpeCoefficient(op1, op2, op3, coef);
    }

    in.close();

    SetFixedCftData();
}

