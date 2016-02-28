//main function, designed for minimization problem
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <windows.h>
#include <fstream>
#include "initialization.h"
#include "solution.h"
#include "Algorithm.h"
#include "FEoutput.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"
#include "bbob09supportfunctions.h"
//#include <vld.h>

using std::string;
using namespace std;

int main()
{
	srand(time(NULL));
	
	extern Solution solution;
	extern FEvariable fevariable;
	Initialization initialization;
	Algorithm algorithm;
	FEoutput feoutput;
	SolutionOperator SOperator;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	BBOB09support bbob09support;

	solution.Tech_num=1;
	solution.dir_file="new-joe-kuo-6.21201.txt";
	//*1 for random, 2 for chaos, 3 for opposition, 
	//*4 for Quasi Opposition, 5 for Quasi Interpolation, 6 for sobol
	//int Max_Tech=6;
	
	//int N=solution.Func_num*10; //dimension of each candidate
	//solution.MaxFE=500*solution.Func_num;//if =0, the result is the initialization result
	solution.Algorithm=1;
	//*1 for DE, 2 for CRO, 3 for PSO
	
	solution.ConstraintHandle="Bounceback";//Equaltobound/Bounceback
	solution.AdaptiveFlag=0;//1 turn on adaptive scheme, 0 turn off adaptive scheme

	//note if DE is used, P=5*N is suggested, revised inside function loop
	solution.TargetBest=0;
	solution.Error=0.00000001;
	solution.Looptime=10;
	
	solution.trialID=1;
	int StartFunction=3;//the first function that will be tested
	int MaxFunctionTested=1;//how many function is not being tested, 1 stands for 1 function.
	//Control variables: solution.D, MaxFE/D, RA
	solution.D=30;
	int IfUseBestP = 0;//if to use data in 2-D Array BestP, 0 is not to use, 1 is to use.
	int FixedEANPsize = 10;
	bbob09.BBOBparameterInitial();

	solution.CECorBBOB=3;//select which benchmark to use
	//0 for CEC14 expensive, 1 for BBOB09, 2 for CEC13, 3 for CEC14
	int RAorP=0;//select whether to use P or RA as x-axis
	//0 for P
	//1 for RA
	int FreeIni = 0;//decide if iniT is free or not
	//0 means non-free IniT
	//1 means free IniT, traditional assumption
	int Pini_control;
	if (solution.Algorithm==1)
	{
		Pini_control=4; //for DE
	}
	else
	{
		Pini_control=1;
	}
	//int Pstep_control=1;
	int Pstep_control=5;//for 10,000 FE/D test
	double RAini_control=0.01;//RA from 1% - 30%, while P>=1;
	double RAstep_control=0.005;//incremental of RA
	int MaxRAcount=50;//maximum number of RA tested
	int MaxPcount=100;
	int RAcount=0;
	int Pcount=0;
	solution.CRODecomS = 5;

	//int FEarray[8] = {25,50,75,100,150,200,300,500};
	int FEarray[1] = {50};
	int MaxFEperDcount = 1;
	int fecount=0;

	//simulation information:
	cout<<"Technique: "<<solution.Tech_num<<endl;
	cout<<"Algorithm: "<<solution.Algorithm<<endl;
	cout<<"Function set: "<<solution.CECorBBOB<<endl;
	cout<<"Function tested: "<<StartFunction<<"-"<<StartFunction+MaxFunctionTested-1<<endl;
	cout<<"Loop time: "<<solution.Looptime<<endl;
	cout<<"BestP? "<<IfUseBestP;
	if (IfUseBestP==0)
		cout<<" Fixed NP: "<<FixedEANPsize<<endl;
	else
		cout<<endl;
	cout<<"RAorP? "<<RAorP<<endl;
	cout<<"FE tested: ";
	for (int i=0;i<MaxFEperDcount;i++)
	{
		cout<<FEarray[i]<<" ";
	}
	cout<<endl;

	int **BestP = new int*[24];
	for (int i=0;i<24;i++)
	{
		BestP[i] = new int[MaxFEperDcount];
		
		for (int j=0;j<MaxFEperDcount;j++)
		{
			switch (solution.Algorithm)
			{
				case 1://using DE
					if (solution.Tech_num==3)
						//BestP[i][j-firstj] = BestP_DE[i][j];
						//BestP[i][j] = BestP_DE[i][j];
						BestP[i][j] = 0;
					else if (solution.Tech_num==1)
						//BestP[i][j-firstj] = BestP_DE_RNG[i][j];
						//BestP[i][j] = BestP_DE_RNG[i][j];
						BestP[i][j] = 0;
					break;
				case 2://using CRO
					if (solution.Tech_num==3)
						//BestP[i][j-firstj] = BestP_CRO[i][j];
						//BestP[i][j] = BestP_CRO[i][j];
						BestP[i][j] = 0;
					else if (solution.Tech_num==1)
						//BestP[i][j-firstj] = BestP_CRO_RNG[i][j];
						//BestP[i][j] = BestP_CRO_RNG[i][j];
						BestP[i][j] = 00;
					break;
				case 3://using PSO
					//BestP[i][j] = BestP_PSO[i][j];
					break;
			}
		}
	}
	
	solution.Min=new double[solution.D];
	solution.Max=new double[solution.D];
	solution.FEused=new int[solution.Looptime];
	//double *Loopbest=new double[solution.Looptime];
	//double *Loopworst=new double[solution.Looptime];
	fevariable.Loopresults=new double[solution.Looptime];
	
	double ***BestIniRecord=new double**[MaxFEperDcount];//to record the averaged best initial solutions
	//to record all functions on all FEperD
	double **BestResultperPorRA=new double*[MaxFEperDcount];//8 is the number of functions
	double **MeanResultperPorRA=new double*[MaxFEperDcount];//8 is the number of functions
	double **StdResultperPorRA=new double*[MaxFEperDcount];//8 is the number of functions
	//3 dimensions [FE/D][f][PorRA]
	double ***BestResult_perPorRA_allFED=new double**[MaxFEperDcount];
	double ***MeanResult_perPorRA_allFED=new double**[MaxFEperDcount];
	double ***StdResult_perPorRA_allFED=new double**[MaxFEperDcount];
	for (int k=0;k<MaxFEperDcount;k++)
	{
		BestResult_perPorRA_allFED[k]=new double*[MaxFunctionTested];
		MeanResult_perPorRA_allFED[k]=new double*[MaxFunctionTested];
		StdResult_perPorRA_allFED[k]=new double*[MaxFunctionTested];
		BestIniRecord[k]=new double*[MaxFunctionTested];
	}
	if (RAorP==0)//using P as x-axis
	{
		for (int i=0;i<MaxFEperDcount;i++)
		{
			BestResultperPorRA[i]=new double[MaxPcount+1];//to record the best result under each P value
			MeanResultperPorRA[i]=new double[MaxPcount+1];//to record the best result under each P value
			StdResultperPorRA[i]=new double[MaxPcount+1];//to record the best result under each P value
		}
		for (int k=0;k<MaxFEperDcount;k++)
		{
			for (int i=0;i<MaxFunctionTested;i++)
			{
				BestResult_perPorRA_allFED[k][i]=new double[MaxPcount+1];
				MeanResult_perPorRA_allFED[k][i]=new double[MaxPcount+1];
				StdResult_perPorRA_allFED[k][i]=new double[MaxPcount+1];
				BestIniRecord[k][i]=new double[MaxPcount+1];
			}
		}
	}
	else
	{
		for (int i=0;i<MaxFEperDcount;i++)
		{
			BestResultperPorRA[i]=new double[MaxRAcount+1];//to record the best result under each RA value
			MeanResultperPorRA[i]=new double[MaxRAcount+1];//to record the best result under each RA value
			StdResultperPorRA[i]=new double[MaxRAcount+1];//to record the best result under each RA value
		}
		for (int k=0;k<MaxFEperDcount;k++)
		{
			for (int i=0;i<MaxFunctionTested;i++)
			{
				BestResult_perPorRA_allFED[k][i]=new double[MaxRAcount+1];
				MeanResult_perPorRA_allFED[k][i]=new double[MaxRAcount+1];
				StdResult_perPorRA_allFED[k][i]=new double[MaxRAcount+1];
				BestIniRecord[k][i]=new double[MaxRAcount+1];
			}
		}
	}

	double *RAvalue;
	int *Pvalue;
	int Count=0;
	int MaxCount=0;
	if (RAorP==0)
		MaxCount=MaxPcount;
	else
		MaxCount=MaxRAcount;

	LARGE_INTEGER Frequency;//计数器频率  
	LARGE_INTEGER start_PerformanceCount;//起始计数器  
	LARGE_INTEGER end_PerformanceCount;//结束计数器 
	double run_time=0; //运行时间  
	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&start_PerformanceCount);
	int FunctionCount=0;

	int TotalProgress = solution.Looptime * MaxFunctionTested * MaxCount * MaxFEperDcount;

	while(FunctionCount<MaxFunctionTested)//test one function first
	{
		//cout<<"Current Function Num being Evaluated, "<<solution.Func_num<<endl;
		if(solution.CECorBBOB==0)
		{
			if (solution.D==10)
			{
				solution.Func_num=StartFunction*3-2;//1 as sphere, etc.
			}
			else if (solution.D==30)
			{
				solution.Func_num=StartFunction*3-1;
			}
			else if (solution.D==100)
			{
				solution.Func_num=StartFunction*3;
			}
		}
		else if (solution.CECorBBOB==1)
		{
			solution.Func_num=StartFunction+FunctionCount;
			bbob09.CalculateParameter(solution.Func_num);
		}
		else if (solution.CECorBBOB==2)
		{
			solution.Func_num=StartFunction+FunctionCount;
		}
		else if (solution.CECorBBOB==3)
		{
			solution.Func_num=StartFunction+FunctionCount;
		}
		
		if (solution.CECorBBOB==0)//initial setting for CEC14
		{
			//FileIO, read the M rotation matrix and o shift vector information
			CEC14.FileIO(solution.D,solution.Func_num);

			//different min max constraints for Function 13-18
			if((solution.Func_num>12)&&(solution.Func_num<16))
			{
				for (int n=0;n<solution.D;n++)
				{
					solution.Min[n]=-32;
					solution.Max[n]=32;
				}
			}
			else if((solution.Func_num>15)&&(solution.Func_num<19))
			{
				for (int n=0;n<solution.D;n++)
				{
					solution.Min[n]=-600;
					solution.Max[n]=600;
				}
			}
			else
			{
				for (int n=0;n<solution.D;n++)
				{
					solution.Min[n]=-20;
					solution.Max[n]=20;
				}
			}
		}
		else if (solution.CECorBBOB==1) //initial setting for BBOB09
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-5;
				solution.Max[n]=5;
			}
		}
		else if (solution.CECorBBOB==2) //initial setting for CEC13
		{
			CEC13.FileIO(solution.D,solution.Func_num);
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-100;
				solution.Max[n]=100;
			}
		}
		else if (solution.CECorBBOB==3) //initial setting for CEC13
		{
			CEC14_normal.FileIO(solution.D,solution.Func_num);
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-100;
				solution.Max[n]=100;
			}
		}
		fecount=0;
		//----------------------------------------------
		//start the FEperD loop
		while (fecount<MaxFEperDcount)
		{
			if (IfUseBestP==1)//to decide what FixedEANPsize to use
			{
				FixedEANPsize = BestP[solution.Func_num-1][fecount];
			}
			//currently D can only take value of 10, 30, and 100
			int FE_preset = FEarray[fecount]*solution.D; //500 for N=10, 1500 for N=30, ....
			int Opt_FE = 0;

			int P=0;
			int Pstep=0;
			double RA=0;//RA from 1% - 30%, while P>=1;
			//double RAmax=0;//maximum RA tested
			double RAstep=0;//incremental of RA
			RAcount=0;
			Pcount=0;
			
			if (RAorP==0)//using P as x-axis
			{
				P=Pini_control;
				if (P<FixedEANPsize)
				{
					P=FixedEANPsize;
				}
				//Pmax=Pmax_control;
				Pstep=Pstep_control;
				RAvalue=new double[MaxPcount+1];
				Pvalue=new int[MaxPcount+1];

				if ((solution.Tech_num==1)||(solution.Tech_num==2)||(solution.Tech_num==6))
				{
					//P=(int)(RA*FE_preset);
					RA = (double)(P)/(double)(FE_preset);
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}
				else //for OBL, QOBL, and QI
				{
					//P=(int)(RA*FE_preset/2);
					if ((P<4)&&(solution.Tech_num==5))
					{
						cout<<"QI must has at least 4 initial solutions!"<<endl;
						getchar();
						exit(0);
					}
					else
					{
						if (P<1)
						{
							P=1;
						}
					}
					RA = (double)(2*P)/(double)(FE_preset);
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (2*P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}
			}
			else //using RA as x-axis
			{
				RA=RAini_control;
				RAstep=RAstep_control;
				RAvalue=new double[MaxRAcount+1];
				Pvalue=new int[MaxRAcount+1];

				if ((solution.Tech_num==1)||(solution.Tech_num==2)||(solution.Tech_num==6))
				{
					P=(int)(RA*FE_preset);
					if (solution.Algorithm==1)//using DE
					{
						if (P<4)
						{
							P=4;
						}
					}
					else
					{
						if (P<1)
						{
							P=1;
						}
					}
					if (P<FixedEANPsize)
					{
						P=FixedEANPsize;
					}
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}	
				else //for OBL, QOBL, and QI
				{
					P=(int)(RA*FE_preset/2);
					if (solution.Algorithm==1)//using DE
					{
						if (P<4)
						{
							P=4;
						}
					}
					else
					{
						if ((P<4)&&(solution.Tech_num==5))
						{
							cout<<"QI must has at least 4 initial solutions!"<<endl;
							exit(0);
						}
						else
						{
							if (P<1)
							{
								P=1;
							}
						}
					}
					if (P<FixedEANPsize)
					{
						P=FixedEANPsize;
					}
					//RA = (double)(2*P)/(double)(FE_preset);
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (2*P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}
			}

			double run_time=0;
			if (RAorP==0)//P
			{
				Count=Pcount;
				MaxCount=MaxPcount;
			}
			else
			{
				Count=RAcount;
				MaxCount=MaxRAcount;
			}

			//----------------------------------------------
			//start of RA loop
			while(Count<MaxCount)
			{
				//cout<<RA<<endl;
				//same for RNG, CNG and SQNG, for OBL, QOBL and QI, = 500-100; 1500-300...
				//for PRNG, CNG and SQNG:
				if ((solution.Tech_num==1)||(solution.Tech_num==2)||(solution.Tech_num==6))
				{
					if (RAorP==0)//using P
					{
						RA = (double)(P)/(double)(FE_preset);
					}
					else //using RA
					{
						P=(int)(RA*FE_preset);
						if (P<4)
						{
							P=4;
						}
						if (P<FixedEANPsize)
						{
							P=FixedEANPsize;
						}
					}
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}
				else //for OBL, QOBL, and QI
				{
					if (RAorP==0)//using P
					{
						RA = (double)(2*P)/(double)(FE_preset);
					}
					else //using RA
					{
						P=(int)(RA*FE_preset/2);
						if (P<4)
						{
							P=4;
						}
						if (P<FixedEANPsize)
						{
							P=FixedEANPsize;
						}
					}
					if (FreeIni==0)
					{
						Opt_FE = (FEarray[fecount]*solution.D) - (2*P);
					}
					else
					{
						Opt_FE = (FEarray[fecount]*solution.D);
					}
				}
		
				solution.IniNP=P;//number of initial solutions
				solution.EANP=FixedEANPsize;//EA population size
				solution.Pdouble=solution.CRODecomS*FixedEANPsize;//upper limit of population size. only initial 2P positions
				//cout<<"RA: "<<RA<<" P: "<<P<<endl;
				solution.S=new double*[solution.CRODecomS*FixedEANPsize];
				for (int p=0;p<solution.CRODecomS*FixedEANPsize;p++)
				{
					solution.S[p]=new double[solution.D];
					for (int d=0;d<solution.D;d++)
					{
						solution.S[p][d]=0;
					}
				}
				solution.IniS=new double*[solution.IniNP];
				for (int p=0;p<solution.IniNP;p++)
				{
					solution.IniS[p]=new double[solution.D];
					for (int d=0;d<solution.D;d++)
					{
						solution.IniS[p][d]=0;
					}
				}
				solution.Y=new double[solution.CRODecomS*FixedEANPsize];
				solution.IniY=new double[solution.IniNP];
				solution.Aconverge=new double[Opt_FE+1];
				solution.MaxFE=Opt_FE;
				for (int i=0;i<solution.Looptime;i++)
				{
					solution.FEused[i]=solution.MaxFE;
				}
		
				for (int i=0;i<solution.Looptime;i++)
				{
					solution.FEused[i]=solution.MaxFE;
				}
				//to record the average convergence
				for (int i=0;i<solution.MaxFE+1;i++)
				{
					solution.Aconverge[i]=0;
				}			
			
				int loop=0;
				for (int i=0;i<solution.Looptime;i++)
				{
					fevariable.Loopresults[i]=0;
				}

				double best_ini_sum=0;
				
				while(loop<solution.Looptime)
				{
					//If in future loop test is constructed, remember to reset the parameters as a new loop
					//initial solutions regarding to the given initialization technique and problem function
					//cout<<FixedEANPsize<<" ";
					solution.EANP=FixedEANPsize;//reset P, CRO will change the population size
					//cout<<"Current Loop: "<<loop+1<<endl;
					//if ((loop+1)%10==0)
					//	cout<<"Current Loop: "<<loop+1<<endl;
					fevariable.Besteverresult=10000000000;
				
					initialization.Initial();
					best_ini_sum+=solution.Y[0];

					algorithm.Optimization(loop);
					
					fevariable.Loopresults[loop]=fevariable.Besteverresult;
			
					loop++;
					
					//progress display:
					int CurrentProgress = FunctionCount*MaxFEperDcount*MaxCount*solution.Looptime
						+ fecount*MaxCount*solution.Looptime + Count*solution.Looptime + loop;
					double Percentage = double(CurrentProgress)/double(TotalProgress);
					printf("%.2lf% % \r",Percentage*100);
				};//end of optimization loop
				
				feoutput.LoopOprocess();//get loop best, worst, mean and std;
			
				//record result to per RA vector:
				//BestResult_perPorRA_allFED[fecount][FunctionCount][Count]=feoutput.Loopbest;
				//MeanResult_perPorRA_allFED[fecount][FunctionCount][Count]=feoutput.Loopmean;
				//StdResult_perPorRA_allFED[fecount][FunctionCount][Count]=feoutput.Loopstd;
				BestResultperPorRA[fecount][Count]=feoutput.Loopbest;
				MeanResultperPorRA[fecount][Count]=feoutput.Loopmean;
				StdResultperPorRA[fecount][Count]=feoutput.Loopstd;
				BestIniRecord[fecount][FunctionCount][Count]=best_ini_sum/solution.Looptime;
				
				//record and update RA and P
				RAvalue[Count]=RA;//record RA values
				Pvalue[Count]=P;
		
				delete []solution.Y;
				delete []solution.IniY;
				for (int p=0;p<solution.CRODecomS*FixedEANPsize;p++)
				{
					delete []solution.S[p];
				}
				for (int p=0;p<solution.IniNP;p++)
				{
					delete []solution.IniS[p];
				}
				delete []solution.IniS;
				delete []solution.S;
				delete []solution.Aconverge;
		
				if (RAorP==0) //using P
				{
					if (FEarray[fecount]>5000)
					{
						if (Count<=15)
						{
							Pstep = 1;
						}
						else if (Count<=30)
						{
							Pstep = 5;
						}
						else
						{
							Pstep = 30;
						}
					}
					else
					{
						if (Count<=15)
						{
							Pstep = Pstep_control;
						}
						else
						{
							Pstep = Pstep_control;
						}
					}
					P=P+Pstep;
					//Pcount++;
				}
				else //using RA
				{
					if (Count<=15)
					{
						//RAstep = 5.0/(double)FE_preset;
						RAstep = RAstep_control;
					}
					else
					{
						RAstep = RAstep_control*2;
					}
					RA=RA+RAstep;
					//RAcount++;
				}
				Count++;
			}//end of RA or P count loop
	
			fecount++;
			//delete []Pvalue;
			//delete []RAvalue;
		}//end of FEperD loop
		//-------------------------------------------

		//---------------------------------------------------------------
		//record data to a 3-D matrix for printing all FE/D to one .csv file
		for (int i=0;i<MaxFEperDcount;i++)
		{
			for (int ra=0;ra<MaxCount;ra++)
			{
				BestResult_perPorRA_allFED[i][FunctionCount][ra]=BestResultperPorRA[i][ra];
				MeanResult_perPorRA_allFED[i][FunctionCount][ra]=MeanResultperPorRA[i][ra];
				StdResult_perPorRA_allFED[i][FunctionCount][ra]=StdResultperPorRA[i][ra];
			}
		}
		//---------------------------------------------------------------

		if (solution.CECorBBOB==0)
		{
			solution.Func_num=solution.Func_num+3;
		}
		else
		{
			solution.Func_num++;
		}
		FunctionCount++;
		//cout<<endl;
	}//end of problem loop

	QueryPerformanceCounter(&end_PerformanceCount);
	run_time = ( end_PerformanceCount.QuadPart - start_PerformanceCount.QuadPart ) / (double)Frequency.QuadPart;

	//-------------------------------------------
	//print the 3-D matrix with all FE/D data to one .csv file and record best PorRA for each F and FE/D
	int **Precord = new int*[MaxFEperDcount];
	double **RArecord = new double*[MaxFEperDcount];
	double **Resultrecord = new double*[MaxFEperDcount];
	for (int k=0;k<MaxFEperDcount;k++)
	{
		Precord[k]=new int[MaxFunctionTested];
		RArecord[k]=new double[MaxFunctionTested];
		Resultrecord[k]=new double[MaxFunctionTested];
	}
	FILE* finalfout;
	char FinalFileName[100];
	sprintf_s(FinalFileName,"ResultOutput_SumData//A%d_Tech%d_Functions%d_%d_F%d-%d_D%d_Mean_Runs%d_FreeIni%d_BestP%d.csv",solution.Algorithm,solution.Tech_num,solution.CECorBBOB,RAorP,StartFunction,StartFunction+MaxFunctionTested-1,solution.D,solution.Looptime, FreeIni,IfUseBestP);
	fstream finalfoutclear(FinalFileName,ios::out);
	finalfoutclear.close();
	errno_t err=fopen_s(&finalfout,FinalFileName,"a+");
	if (err==0)
	{
		if (IfUseBestP==1)
			fprintf(finalfout,"Algorithm %d,Tech %d,D %d, runtime, %f\n",solution.Algorithm,solution.Tech_num,solution.D,run_time);
		else
			fprintf(finalfout,"Algorithm %d,Tech %d,D %d, runtime, %f, Fixed NP, %d\n",solution.Algorithm,solution.Tech_num,solution.D,run_time,FixedEANPsize);
	}
	for (int i=0;i<MaxFunctionTested;i++)
	{
		fprintf(finalfout,"F%d",StartFunction+i);
		if (RAorP==0)
		{
			if (IfUseBestP==1)
			{
				
			}
			else
			{
				fprintf(finalfout,"\nP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout,"%d,",Pvalue[ra]);
					//fout2<<Pvalue[ra]<<",";
				}
			}
		}
		else
		{
			if (IfUseBestP==1)
			{

			}
			else
			{
				fprintf(finalfout,"\nRA,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout,"%f,",RAvalue[ra]);
				}
			}
		}
		for (int k=0;k<MaxFEperDcount;k++)
		{
			if (IfUseBestP==1)
			{
				if (RAorP==0)
					fprintf(finalfout,"\nP,BestP with not fixed NP,");
				else
					fprintf(finalfout,"\nRA,BestP with not fixed NP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout,"%d,",BestP[StartFunction+i-1][k]+(ra*Pstep_control));
					//fout2<<Pvalue[ra]<<",";
				}
				fprintf(finalfout,"\nF%d %dFE/D Mean result,%d,",StartFunction+i,FEarray[k],BestP[StartFunction+i-1][k]);
			}
			else
			{
				fprintf(finalfout,"\nF%d %dFE/D Mean result,",StartFunction+i,FEarray[k]);
			}
			double Rmin=1000000000;
			for (int ra=0;ra<MaxCount;ra++)
			{
				fprintf(finalfout,"%f,",MeanResult_perPorRA_allFED[k][i][ra]);
				if (MeanResult_perPorRA_allFED[k][i][ra]<Rmin)
				{
					Rmin=MeanResult_perPorRA_allFED[k][i][ra];
					Precord[k][i]=Pvalue[ra];
					RArecord[k][i]=RAvalue[ra];
					Resultrecord[k][i]=Rmin;
				}
			}
		}
		for (int i=0;i<20;i++)
		{
			fprintf(finalfout,"\n");
		}
	}

	fprintf(finalfout,"P/RA that achieves best Mean results:\nF and FE/D,Best P,Mean result");
	for (int i=0;i<MaxFunctionTested;i++)
	{
		for (int k=0;k<MaxFEperDcount;k++)
		{
			if (RAorP == 0)
				fprintf(finalfout,"\nF%d %dFE/D best P and Result,%d,%f",StartFunction+i,FEarray[k],Precord[k][i],Resultrecord[k][i]);
			else
				fprintf(finalfout,"\nF%d %dFE/D best RA and Result,%f,%f",StartFunction+i,FEarray[k],RArecord[k][i],Resultrecord[k][i]);
		}
		fprintf(finalfout,"\n");
	}
	fclose(finalfout);
	//-------------------------------------------

	//-------------------------------------------
	//Std result output:
	FILE* finalfout1;
	char FinalFileName1[100];
	sprintf_s(FinalFileName1,"ResultOutput_SumData//A%d_Tech%d_Functions%d_%d_F%d-%d_D%d_Std_Runs%d_FreeIni%d_BestP%d.csv",solution.Algorithm,solution.Tech_num,solution.CECorBBOB,RAorP,StartFunction,StartFunction+MaxFunctionTested-1,solution.D,solution.Looptime, FreeIni, IfUseBestP);
	fstream finalfout1clear(FinalFileName1,ios::out);
	finalfout1clear.close();
	errno_t err1=fopen_s(&finalfout1,FinalFileName1,"a+");
	if (err1==0)
	{
		if (IfUseBestP==1)
			fprintf(finalfout1,"Algorithm %d,Tech %d,D %d, runtime, %f\n",solution.Algorithm,solution.Tech_num,solution.D, run_time);
		else
			fprintf(finalfout1,"Algorithm %d,Tech %d,D %d, runtime, %f, Fixed NP, %d\n",solution.Algorithm,solution.Tech_num,solution.D, run_time, FixedEANPsize);
	}
	for (int i=0;i<MaxFunctionTested;i++)
	{
		fprintf(finalfout1,"F%d",StartFunction+i);
		if (RAorP==0)
		{
			if (IfUseBestP==1)
			{
				
			}
			else
			{
				fprintf(finalfout1,"\nP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout1,"%d,",Pvalue[ra]);
					//fout2<<Pvalue[ra]<<",";
				}
			}
		}
		else
		{
			if (IfUseBestP==1)
			{

			}
			else
			{
				fprintf(finalfout1,"\nRA,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout1,"%f,",RAvalue[ra]);
				}
			}
		}
		for (int k=0;k<MaxFEperDcount;k++)
		{
			if (IfUseBestP==1)
			{
				if (RAorP==0)
					fprintf(finalfout1,"\nP,BestP with not fixed NP,");
				else
					fprintf(finalfout1,"\nRA,BestP with not fixed NP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout1,"%d,",BestP[StartFunction+i-1][k]+(ra*Pstep_control));
					//fout2<<Pvalue[ra]<<",";
				}
				fprintf(finalfout1,"\nF%d %dFE/D Mean result,%d,",StartFunction+i,FEarray[k],BestP[StartFunction+i-1][k]);
			}
			else
			{
				fprintf(finalfout1,"\nF%d %dFE/D Mean result,",StartFunction+i,FEarray[k]);
			}
			double Rmin=1000000000;
			for (int ra=0;ra<MaxCount;ra++)
			{
				fprintf(finalfout1,"%f,",StdResult_perPorRA_allFED[k][i][ra]);
				if (StdResult_perPorRA_allFED[k][i][ra]<Rmin)
				{
					Rmin=StdResult_perPorRA_allFED[k][i][ra];
					Precord[k][i]=Pvalue[ra];
					RArecord[k][i]=RAvalue[ra];
					Resultrecord[k][i]=Rmin;
				}
			}
		}
		for (int i=0;i<20;i++)
		{
			fprintf(finalfout1,"\n");
		}
	}
	
	fprintf(finalfout1,"P that achieves best Std results:\nF and FE/D,Best P,Std result");
	for (int i=0;i<MaxFunctionTested;i++)
	{
		if (RAorP == 0)
		{
			for (int k=0;k<MaxFEperDcount;k++)
			{
				fprintf(finalfout1,"\nF%d %dFE/D best P and Result,%d,%f",StartFunction+i,FEarray[k],Precord[k][i],Resultrecord[k][i]);
			}
		}
		else
		{
			for (int k=0;k<MaxFEperDcount;k++)
			{
				fprintf(finalfout1,"\nF%d %dFE/D best RA and Result,%f,%f",StartFunction+i,FEarray[k],RArecord[k][i],Resultrecord[k][i]);
			}
		}
		fprintf(finalfout1,"\n");
	}
	fclose(finalfout1);
	//-------------------------------------------

	//-------------------------------------------
	//Best Ini Record output
	FILE* finalfout_ini;
	char FinalFileName_ini[100];
	sprintf_s(FinalFileName_ini,"ResultOutput_SumData//A%d_Tech%d_Functions%d_%d_F%d-%d_D%d_BestIni_Runs%d_FreeIni%d_BestP%d.csv",solution.Algorithm,solution.Tech_num,solution.CECorBBOB,RAorP,StartFunction,StartFunction+MaxFunctionTested-1,solution.D,solution.Looptime, FreeIni,IfUseBestP);
	fstream finalfout_iniclear(FinalFileName_ini,ios::out);
	finalfout_iniclear.close();
	errno_t err_ini=fopen_s(&finalfout_ini,FinalFileName_ini,"a+");
	if (err_ini==0)
	{
		if (IfUseBestP==1)
			fprintf(finalfout_ini,"Algorithm %d,Tech %d,D %d, runtime, %f\n",solution.Algorithm,solution.Tech_num,solution.D,run_time);
		else
			fprintf(finalfout_ini,"Algorithm %d,Tech %d,D %d, runtime, %f, Fixed NP, %d\n",solution.Algorithm,solution.Tech_num,solution.D,run_time,FixedEANPsize);
	}
	for (int i=0;i<MaxFunctionTested;i++)
	{
		fprintf(finalfout_ini,"F%d",StartFunction+i);
		if (RAorP==0)
		{
			if (IfUseBestP==1)
			{
				
			}
			else
			{
				fprintf(finalfout_ini,"\nP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout_ini,"%d,",Pvalue[ra]);
					//fout2<<Pvalue[ra]<<",";
				}
			}
		}
		else
		{
			if (IfUseBestP==1)
			{

			}
			else
			{
				fprintf(finalfout_ini,"\nRA,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout_ini,"%f,",RAvalue[ra]);
				}
			}
		}
		for (int k=0;k<MaxFEperDcount;k++)
		{
			if (IfUseBestP==1)
			{
				if (RAorP==0)
					fprintf(finalfout_ini,"\nP,BestP with not fixed NP,");
				else
					fprintf(finalfout_ini,"\nRA,BestP with not fixed NP,");
				for (int ra=0;ra<MaxCount;ra++)
				{
					fprintf(finalfout_ini,"%d,",BestP[StartFunction+i-1][k]+(ra*Pstep_control));
					//fout2<<Pvalue[ra]<<",";
				}
				fprintf(finalfout_ini,"\nF%d %dFE/D Mean result,%d,",StartFunction+i,FEarray[k],BestP[StartFunction+i-1][k]);
			}
			else
			{
				fprintf(finalfout_ini,"\nF%d %dFE/D Mean result,",StartFunction+i,FEarray[k]);
			}
			double Rmin=1000000000;
			for (int ra=0;ra<MaxCount;ra++)
			{
				fprintf(finalfout_ini,"%f,",BestIniRecord[k][i][ra]);
				if (BestIniRecord[k][i][ra]<Rmin)
				{
					Rmin=BestIniRecord[k][i][ra];
					Precord[k][i]=Pvalue[ra];
					RArecord[k][i]=RAvalue[ra];
					Resultrecord[k][i]=Rmin;
				}
			}
		}
		for (int i=0;i<20;i++)
		{
			fprintf(finalfout_ini,"\n");
		}
	}

	fprintf(finalfout_ini,"P/RA that achieves best Mean results:\nF and FE/D,Best P,Mean result");
	for (int i=0;i<MaxFunctionTested;i++)
	{
		for (int k=0;k<MaxFEperDcount;k++)
		{
			if (RAorP == 0)
				fprintf(finalfout_ini,"\nF%d %dFE/D best P and Result,%d,%f",StartFunction+i,FEarray[k],Precord[k][i],Resultrecord[k][i]);
			else
				fprintf(finalfout_ini,"\nF%d %dFE/D best RA and Result,%f,%f",StartFunction+i,FEarray[k],RArecord[k][i],Resultrecord[k][i]);
		}
		fprintf(finalfout_ini,"\n");
	}
	fclose(finalfout_ini);
	//-----------------------------------------------

	delete []solution.Min;
	delete []solution.Max;
	delete []solution.FEused;
	delete []fevariable.Loopresults;	
	bbob09.Release();
	for (int k=0;k<MaxFEperDcount;k++)
	{
		for (int i=0;i<MaxFunctionTested;i++)
		{
			delete []BestIniRecord[k][i];
			delete []BestResult_perPorRA_allFED[k][i];
			delete []MeanResult_perPorRA_allFED[k][i];
			delete []StdResult_perPorRA_allFED[k][i];
		}
		delete []BestIniRecord[k];
		delete []BestResult_perPorRA_allFED[k];
		delete []MeanResult_perPorRA_allFED[k];
		delete []StdResult_perPorRA_allFED[k];
		delete []BestResultperPorRA[k];
		delete []MeanResultperPorRA[k];
		delete []StdResultperPorRA[k];
	}
	delete []BestIniRecord;
	delete []BestResult_perPorRA_allFED;
	delete []MeanResult_perPorRA_allFED;
	delete []StdResult_perPorRA_allFED;
	delete []BestResultperPorRA;
	delete []MeanResultperPorRA;
	delete []StdResultperPorRA;
	delete []RAvalue;
	delete []Pvalue;

	cout<<"program end";
	getchar();
	return 0;
}