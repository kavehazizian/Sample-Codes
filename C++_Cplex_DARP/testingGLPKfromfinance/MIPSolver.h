//
// MIPSolver.h
 
#ifndef __FinancialSamples__MIPSolver__
#define __FinancialSamples__MIPSolver__
 
#include "LPSolver.h"
 
class MIPSolver : public LPSolver {
public:
MIPSolver(Matrix &A, const std::vector<double> &b, const std::vector<double> &c);
MIPSolver(const MIPSolver &p);
~MIPSolver();
MIPSolver &operator=(const MIPSolver &p);
 
void setColInteger(int colNum);
void setColBinary(int colNum);
virtual ResultType solve(std::vector<double> &result, double &objValue);
};
 
#endif /* defined(__FinancialSamples__MIPSolver__) */
 
//
