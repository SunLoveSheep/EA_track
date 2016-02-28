/*

This is a summation function that is used to control what algorithm to implement in simulation.

*/

#include "solution.h"
#include "Algorithm.h"
#include "DE.h"
#include "CRO.h"
#include "PSO.h"

Algorithm::Algorithm()
{
	
}

Algorithm::~Algorithm()
{
	
}

void Algorithm::Optimization(int loop)
{
	extern Solution solution;
	DE de;
	CRO cro;
	PSO pso;
	
	//take solution.Algorithm to decide which case (algorith) to use. 1 for DE, 2 for CRO, and 3 for PSO
	switch(solution.Algorithm)
	{
	case 1:
		de.DEprocess(loop);
		break;
	case 2:
		cro.CROprocess(loop);
		break;
	case 3:
		pso.PSOprocess(loop);
		break;
	}
}
