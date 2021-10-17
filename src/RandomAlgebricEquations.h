#include <string>
#include "IFunctionEquations.h"
using namespace std;

class RandomAlgebricEquations : public IFunctionEquations<float_type>
{
private:
    int equationNumber;
    int parameterNumber;
    int functionNumber;
    int solutionNumber;
    vector<vector<float_type> > solutions;
    vector<vector<float_type> > parameters;
    vector<vector<vector<float_type> > > coefficients;
    vector<vector<float_type> > constTerms;
    int maxOrder;
    int tryOrder;
    static constexpr float_type minGap = 0.01;

    float_type EvaluateEquation(float_type input, int equationId, vector<vector<float_type> >& polynomials);

    void CalculateVariables(float_type input, vector<vector<float_type> >& polynomials, vector<float_type>& x);

    static void SaveVector(vector<float_type>& v, ofstream& out);
    static void LoadVector(vector<float_type>& v, ifstream& out);

    void AllocateData();
public:
    RandomAlgebricEquations() {};
    RandomAlgebricEquations(int functionNumber, int maxOrder, int tryOrder = -1);
    ~RandomAlgebricEquations();

    void GenerateRandomParameters();
    void SaveToFile(string file = "");

    void LoadFromFile(string file);

    int EquationNumber() { return this->equationNumber; }

    int ParameterNumber() { return this->parameterNumber; }

    int SolutionNumber() { return this->solutionNumber; };

    void SetParameter(int parameterId, float_type value);

    float_type GetParameter(int parameterId);

    float_type EvaluateEquation(float_type input, int equationId);

    float_type VerifySolution(float_type input, int equationId);

    float_type EquationDerivativeByParameter(float_type input, int equationId, int parameterId);

    float_type EquationDerivativeByParameter(float_type input, int equationId, int parameterId1, int parameterId2);

    float_type EvaluateConstraints();

    float_type ConstraintsDerivativeByParameter(int parameterId);

    vector<float_type> GenerateInputs(int size);

    void OutputParameters();

    void OutputSolutions();

    // do nothing.
    void OnParameterUpdated() {};

    float_type GetCoefficient(int equationId, int var1, int var2);
    float_type GetConstTerm(int equationId, int order);
    float_type GetSolution(int parameterId);
};
