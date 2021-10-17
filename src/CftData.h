#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "common.h"
using namespace std;


struct CftConfig
{
    // OperatorNumbers[i] = number of primary operators with spin-i.
    vector<int> OperatorNumbers;

    int D;

    CftConfig() {};

    CftConfig(int d) : D(d) {};
};

std::ostream& operator<< (std::ostream& out, const CftConfig& config);

struct PrimaryInfo
{
    int Id;

    float_type Dim;

    int Spin;

    PrimaryInfo() {};

    PrimaryInfo(int id, float_type dim, int spin) : Id(id), Dim(dim), Spin(spin) {};
};

struct OpeCoefficientKey
{
    int Op1;

    int Op2;

    int Op3;

    OpeCoefficientKey(int op1, int op2, int op3)
    {
        Op1 = op1;
        Op2 = op2;
        Op3 = op3;

        if (Op1 > Op2) swap(Op1, Op2);
        if (Op2 > Op3) swap(Op2, Op3);
        if (Op1 > Op2) swap(Op1, Op2);
    }

    bool operator<(const OpeCoefficientKey& other) const
    {
        if (Op1 != other.Op1) {
            return Op1 < other.Op1;
        } else if (Op2 != other.Op2) {
            return Op2 < other.Op2;
        } else {
            return Op3 < other.Op3;
        }
    }

    bool operator==(const OpeCoefficientKey& other) const
    {
        return Op1 == other.Op1 && Op2 == other.Op2 && Op3 == other.Op3;
    }
};

// represents the CFT data in the truncated CFT
class CftData
{
private:
    vector<PrimaryInfo> ops;
    vector<int> opsNumbers;
    map<OpeCoefficientKey, double> opeCoefs;
    int maxPrimaryId;

    void SetFixedCftData();
public:
    // spacetime dimension;
    int D;

    int StressTensorId;

    // constructor
    CftData() {};
    CftData(CftConfig config);

    // destructor
    ~CftData() {};

    // randomly generate CFT data.
    void GenerateRandomData();

    // get Operator Info for the operator with given id.
    PrimaryInfo GetPrimaryInfo(int id) const;

    float_type GetPrimaryDim(int id) const { return this->ops[id].Dim; };

    int GetPrimarySpin(int id) const { return this->ops[id].Spin; };

    // set the scaling dimension of the primary operator with given id.
    void SetPrimaryDim(int id, double dim);

    // get the number of primaries with given spin.
    int PrimaryNumber(int spin) const;

    // gets the max primary id.
    int MaxPrimaryId() { return this->maxPrimaryId; }

    // gets the max id of scalar operators.
    int MaxScalarId() { return this->opsNumbers[0]; }

    // gets OPE coefficients
    float_type GetOpeCoefficient(int id1, int id2, int id3);
    float_type GetOpeCoefficient(OpeCoefficientKey &key);

    // set the ope coefficient
    void SetOpeCoefficient(int id1, int id2, int id3, float_type coef);
    void SetOpeCoefficient(OpeCoefficientKey &key, float_type coef);

    // output the CFT data related to given operators.
    void Output(vector<int> ops);

    void Save(string file);

    void LoadFromFile(string file);
};
