#include "solution.h"
#include <iostream>
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"

using namespace std;

Solution solution;

SolutionOperator::SolutionOperator()
{
	
}

SolutionOperator::~SolutionOperator()
{
	
}

void SolutionOperator::UpdateY()
{
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	//solution.Y=new double[solution.P];
	double *temp=new double[solution.D];

	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			temp[n]=solution.S[p][n];
		}
				
		if (solution.CECorBBOB==0)//use CEC14 expensive
		{
			solution.Y[p]=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==1) //use BBOB
		{
			solution.Y[p]=bbob09.FunctionCalculation(temp,solution.Func_num);
		}
		else if (solution.CECorBBOB==2) //use CEC13
		{
			solution.Y[p]=CEC13.cec13_problems(temp,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==3) //use CEC14
		{
			solution.Y[p]=CEC14_normal.cec14_problems(temp,solution.D,solution.Func_num);
		}
	}

	delete []temp;
	//delete temp;
}