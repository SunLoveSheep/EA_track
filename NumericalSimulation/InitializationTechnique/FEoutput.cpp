/*
This .cpp file contains functions that initiate variables & arrays used for recording convergence data from the optimization process.
*/

#include "solution.h"
#include "FEoutput.h"
#include <iostream>
#include <math.h>
#include "cec14_eotp.h"

using namespace std;

FEvariable fevariable;

FEoutput::FEoutput()
{
	extern Solution solution;
	
	//Bestever=new double[solution.D];
	//BesteverFE=0;
	//ConvergeCounter=20;
	//Convergence=new double[ConvergeCounter];

	//Loopresults=new double[solution.Looptime];
	Loopbest=0;
	Loopworst=0;
	Loopmean=0;
	Loopstd=0;
}

FEoutput::~FEoutput()
{
	//delete []Bestever;
	//delete []Convergence;
	//delete []Loopresults;
}

void FEoutput::Initial()
{
	extern Solution solution;

	Bestever=new double[solution.D];
	BesteverFE=0;
	ConvergeCounter=solution.MaxFE;
	Convergence=new double[ConvergeCounter+1];
	ConvergenceFirst10=new double[ConvergeCounter+1];
}

void FEoutput::Release()
{
	delete []Bestever;
	delete []Convergence;
	delete []ConvergenceFirst10;
}

//given the current FE, update the best ever solution. Both its objective value and variable values
void FEoutput::UpdateBest(int FE)
{
	extern Solution solution;
	
	for (int p=0;p<solution.EANP;p++)
	{
		if (fevariable.Besteverresult>solution.Y[p])
		{
			fevariable.Besteverresult=solution.Y[p];
			for (int n=0;n<solution.D;n++)
			{
				Bestever[n]=solution.S[p][n];
			}
			BesteverFE=FE;
		}
	}
}

/*void FEoutput::RecordConvergence(int FE)
{
	extern Solution solution;
	cec14_eotp CEC14;
	
	if(solution.Algorithm==2)//CRO
	{
		//record in total ConvergeCounter points
		//use FE=25 as eg. if means FE from 24->25
		if (FE%(int(solution.MaxFE/ConvergeCounter))==0)//record totally 20 points
		{
			Convergence[FE/(int(solution.MaxFE/ConvergeCounter))-1]
			=CEC14.cec14_eotp_problems(Bestever,solution.D,solution.Func_num);
		}

		//else if means FE from 24->26, FE=FE+2
		else if (((FE-1)%(int(solution.MaxFE/ConvergeCounter))==0)&&(FE!=1))
		{
			//Convergence[(FE-1)/(int(solution.MaxFE/ConvergeCounter))-1]
			//=function.FunctionEvaluationCEC13(Bestever,solution.Func_num);
			Convergence[(FE-1)/(int(solution.MaxFE/ConvergeCounter))-1]
			=CEC14.cec14_eotp_problems(Bestever,solution.D,solution.Func_num);
		}

		else if (((FE-2)%(int(solution.MaxFE/ConvergeCounter))==0)&&(FE!=2))
		{
			//Convergence[(FE-2)/(int(solution.MaxFE/ConvergeCounter))-1]
			//=function.FunctionEvaluationCEC13(Bestever,solution.Func_num);
			Convergence[(FE-2)/(int(solution.MaxFE/ConvergeCounter))-1]
			=CEC14.cec14_eotp_problems(Bestever,solution.D,solution.Func_num);
		}
	}
	else if (solution.Algorithm==1)
	{
		Convergence[FE]=CEC14.cec14_eotp_problems(Bestever,solution.D,solution.Func_num);
	}
	
}*/

//record the convergence for each loop
void FEoutput::Output(int FE)
{
	extern Solution solution;

	//the optimization has met the stop criterion
	if (FE>=solution.MaxFE)
	{
		/*
		cout<<"at FE: "<<BesteverFE<<" Best ever result found: "<<fevariable.Besteverresult<<endl;
		cout<<"The solution: "<<endl;
		for (int n=0;n<solution.D;n++)
		{
			cout<<Bestever[n]<<" ";
		}
		cout<<endl;
		*/
		//cout<<"Convergence: "<<endl;
		for (int i=0;i<ConvergeCounter;i++)//there are totally 20 points recorded
		{
			//cout<<Convergence[i]<<" ";
			//if((i+1)%5==0)
			//{
				//cout<<endl;
			//}
			
			solution.Aconverge[i]+=Convergence[i];
		}
	}
}

//record the best and worst optimal result from all performed loops
void FEoutput::Loopsort()
{	
	extern Solution solution;
	
	double max=-10000000000000, min=1000000000000000;
	
	for(int i=0;i<solution.Looptime;i++)
	{
		if(fevariable.Loopresults[i]>max)
		{
			max=fevariable.Loopresults[i];
		}
		if(fevariable.Loopresults[i]<min)
		{
			min=fevariable.Loopresults[i];
		}
	}

	Loopbest=min;
	Loopworst=max;	
}

//given the cumulated optimization results from each loop, calculate the mean and std
void FEoutput::Loopmeanstd()
{
	extern Solution solution;
	
	double sum=0;
	for (int i=0;i<solution.Looptime;i++)
	{
		sum=sum+fevariable.Loopresults[i];
	}
	Loopmean=sum/solution.Looptime;

	double sumstd=0;
	for (int i=0;i<solution.Looptime;i++)
	{
		sumstd=sumstd+(fevariable.Loopresults[i]-Loopmean)*(fevariable.Loopresults[i]-Loopmean);
	}
	Loopstd=sqrt(sumstd/solution.Looptime);	
}

void FEoutput::LOutput()
{
	cout<<endl<<"Loop best: "<<Loopbest<<endl;
	cout<<"Loop worst: "<<Loopworst<<endl;
	cout<<"Loop mean and std: "<<Loopmean<<" "<<Loopstd<<endl;
}

void FEoutput::LoopOprocess()
{
	Loopsort();
	Loopmeanstd();
	//LOutput();

	//delete []Loopresults;
}
