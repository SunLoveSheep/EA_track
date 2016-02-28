//solution candidate definition
#include <string>
#include <vector>

using std::string;

#ifndef _SOLUTION_H
#define _SOLUTION_H

struct Solution
{
	double **S;//a set of solutions
	double **IniS;//initial solutions, more than S
	
	int EANP;//fixed size initial population size for EA
	int Pdouble;
	int IniNP;//number of solutions initialized by techniques, can be larger than P
	double *Min;
	double *Max;//min and max value of the variable
	int D;//number of variables in each solution, dimension
	double *Y;//objective function value of the solution
	double *IniY;//objective function value of the initial solutions
	int Func_num;
	double Bestever;//to record the best ever solution found
	double *Aconverge;//convergence recording
	int Tech_num;//1 for chaos, 2 for opposition
	double TargetBest;//theoritically optima point
	double Error;//if the difference between the current best and TargetBest is smaller than this Error, stop and record FE
	int *FEused;//to record how many FE used to reach the target range
	bool AdaptiveFlag;//to turn on/off the adaptive parameter updating scheme for each algorithm
	char *dir_file;//the file recording the sobol directions

	int Algorithm;//name of the optimization algorithm used
	string ConstraintHandle;//name of the method used to handle the constraints violation

	int MaxFE;//maximum FE times
	int Looptime;//loop time

	int CECorBBOB;//to decide whether CEC or BBOB benchmark functions to use
	int trialID;//to decide the rseed used for BBOB09
	int CRODecomS;
};

class SolutionOperator
{
public:
	SolutionOperator();
	virtual~SolutionOperator();

	void UpdateY();
};

#endif