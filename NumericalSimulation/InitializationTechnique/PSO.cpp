#include "PSO.h"
#include "solution.h"
#include "ConstraintHandle.h"
#include "FEoutput.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "Initialization.h"
#include "bbob09functions.h"
#include <iostream>

using namespace std;

PSOparameter PSOp;

PSO::PSO()
{
	extern Solution solution;

	PSOp.c1=2;
	PSOp.c2=2;
	PSOp.w=0.8;
	PSOp.gbest=10000000000;
	
}

PSO::~PSO()
{
	
}

void PSO::Initial()
{
	extern Solution solution;

	for (int p=0;p<solution.EANP;p++)
	{
		pbest[p]=solution.Y[p];
	}

	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			pbestposition[p][n]=solution.S[p][n];
		}
	}

	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			v[p][n]=0;
		}
	}

	//vmin=new double[solution.N];
	for (int n=0;n<solution.D;n++)
	{
		vmax[n]=(solution.Max[n]-solution.Min[n]);
	}
}

void PSO::PSOparticleupdate()
{
	extern Solution solution;
	ConstraintHandle constrainthandle;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	FEoutput feoutput;
	SolutionOperator SOperator;
	Initialization initialization;
	BBOB09 bbob09;

	double *temp=new double[solution.D];
	double *vtemp=new double[solution.D];

	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			double r1,r2;
			r1=initialization.RandomGen(0.0,1.0);
			r2=initialization.RandomGen(0.0,1.0);
			vtemp[n]=v[p][n]*PSOp.w+PSOp.c1*r1*(pbestposition[p][n]-solution.S[p][n])
				+PSOp.c2*r2*(gbestposition[n]-solution.S[p][n]);
			//velocity constraints fixing required
			if (vtemp[n]>vmax[n])
			{
				vtemp[n]=vmax[n];
			}
			else if (vtemp[n]<-vmax[n])
			{
				vtemp[n]=-vmax[n];
			}
			temp[n]=solution.S[p][n]+vtemp[n];
			constrainthandle.Handleconstaint(temp);
		}		
		double tempresult=0;
		if (solution.CECorBBOB==0)
			tempresult=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempresult=bbob09.FunctionCalculation(temp,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempresult=CEC13.cec13_problems(temp,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempresult=CEC14_normal.cec14_problems(temp,solution.D,solution.Func_num);
		//elite selection:
		if (tempresult<solution.Y[p])
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=temp[n];
				v[p][n]=vtemp[n];
			}
			solution.Y[p]=tempresult;
		}
		else
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=solution.S[p][n];
				v[p][n]=v[p][n];
			}
		}
		
		FE++;
	}

	delete []temp;
	delete []vtemp;
}

void PSO::PSOpbestupdate()
{
	extern Solution solution;
	for (int p=0;p<solution.EANP;p++)
	{
		if (pbest[p]>solution.Y[p])//designed for minimization problem
		{
			pbest[p]=solution.Y[p];
			for (int n=0;n<solution.D;n++)
			{
				pbestposition[p][n]=solution.S[p][n];
			}
		}
	}
}

void PSO::PSOgbestupdate()
{
	extern Solution solution;
	
	for (int p=0;p<solution.EANP;p++)
	{
		if (PSOp.gbest>solution.Y[p])
		{
			PSOp.gbest=solution.Y[p];
			for(int n=0;n<solution.D;n++)
			{
				gbestposition[n]=solution.S[p][n];
			}
		}
	}
	
}

void PSO::PSOprocess(int loop)
{
	extern Solution solution;
	SolutionOperator SOperator;
	FEoutput feoutput;
	cec14_eotp CEC14;
	
	double bestever=1000000000;
		
	FE=0;
	SOperator.UpdateY();
	feoutput.Initial();
	pbest=new double[solution.EANP];
	pbestposition=new double*[solution.EANP];
	v=new double*[solution.EANP];
	for (int p=0;p<solution.EANP;p++)
	{
		pbestposition[p]=new double[solution.D];
		v[p]=new double[solution.D];
	}
	vmax=new double[solution.D];
	gbestposition=new double[solution.D];

	Initial();
	PSOp.gbest=1000000000;
	PSOp.w=3;

	do
	{
		PSOp.w=PSOp.w-PSOp.w*(double(FE)/double(solution.MaxFE));//less aggresive adaptation
		//if (PSOp.w<0.6)
		//{
		//	PSOp.w=0.6;
		//}
		
		PSOparticleupdate();
		PSOgbestupdate();
		PSOpbestupdate();
		
		//SOperator.UpdateY();
		//feoutput.FEOprocess(FE);
		feoutput.UpdateBest(FE);
		
		/*for (int i=0;i<solution.P;i++)
		{
			cout<<solution.Y[i]<<" ";
		}
		cout<<endl;
		getchar();*/

	}while(FE<solution.MaxFE);
	//cout<<fevariable.Besteverresult<<endl;
	//getchar();
	feoutput.Release();
	delete []gbestposition;
	delete []vmax;
	delete []pbest;
	for (int i=0;i<solution.EANP;i++)
	{
		delete []v[i];
		delete []pbestposition[i];
	}
	delete []v;
	delete []pbestposition;
}