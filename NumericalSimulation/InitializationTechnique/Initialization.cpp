/*
This .cpp file contains functions that generate initial solutions by different techniques. The function used to read selection of technique
and call according function is also included.
*/

#include <string>
#include "solution.h"
#include "Initialization.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

using std::string;
using namespace std;

Initialization::Initialization()
{

}

Initialization::~Initialization()
{
	
}

//a common random number generator that returns a double random number between min and max. Seed is based on system time.
double Initialization::RandomGen(double min, double max)
{
	double Min = (min*1000000);
    double Max = (max*1000000);
    double Rand = rand()*rand();
    double Result;

	if (min!=max)
	{
		Result = ((int)(Rand)%(int)(Max-Min))+Min;
	}
	else
	{
		Result=min;
	}
    return Result/1000000.0;
}

//a common gaussian random number generator that returns a double random number based on miu and sigma
double Initialization::GaussRandomGen(double miu, double sigma)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
	
	X=X*sigma+miu;
    return X;
}

//This is a special function for our design. Number of initial solutions will be larger than solutions used by EA. We use this function
//to sort and select the best solutions from initial solutions and form those used for EA.
void MatchfromIniSolutionsToEASolutions(double *temparray)
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *tempX = new double[solution.D];
	double min=1000000000000000000;
	//double tempY;
	double *temparrayY = new double[solution.IniNP];
	int *tempCount = new int[solution.IniNP];

	for (int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			tempX[n]=solution.IniS[p][n];
		}
		/*if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(tempX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(tempX,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(tempX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(tempX,solution.D,solution.Func_num);
		temparrayY[p]=tempY;*/
		solution.IniY[p]=temparrayY[p];
	}

	//for (int i=0;i<solution.IniNP;i++)
	//{
	//	cout<<temparrayY[i]<<" ";
	//}
	//cout<<endl;
	sort(temparrayY,temparrayY+solution.IniNP);//sorting solutions
	//for (int i=0;i<solution.IniNP;i++)
	//{
	//	cout<<temparrayY[i]<<" ";
	//}
	//cout<<endl;
	//getchar();

	for (int p=0;p<solution.EANP;p++)
	{
		solution.Y[p]=temparrayY[p];
	}
	for (int p=0;p<solution.EANP;p++)
	{
		for (int q=0;q<solution.IniNP;q++)
		{
			if (temparrayY[p]==solution.IniY[q])
				tempCount[p]=q;
		}

		for (int d=0;d<solution.D;d++)
		{
			solution.S[p][d]=solution.IniS[tempCount[p]][d];
		}
	}

	delete []tempX;
	//delete []temparrayY;
	delete []tempCount;
}

//Random number generator
void Initialization::Random()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];
	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		for (int i=0;i<solution.D;i++)
		{
			solution.IniS[p][i]=RandomGen(solution.Min[i],solution.Max[i]);
			X[i]=solution.IniS[p][i];
		}
		if (solution.CECorBBOB==0)
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	}

	MatchfromIniSolutionsToEASolutions(Y);
	delete []X;
	delete []Y;
}

//Chaotic number generator
void Initialization::Chaos()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];
	
	//numerical precision fixed at .xxxxxxx
	cout << setiosflags(ios::fixed) << setprecision(7);

	double **R=new double*[solution.IniNP];
	for (int p=0;p<solution.IniNP;p++)
	{
		R[p]=new double[solution.D];
	}

	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		R[p][0]=RandomGen(0,1); //for each solution, the first variable is randomly initialized
		for (int i=1;i<solution.D;i++)
		{
			R[p][i]=4*R[p][i-1]*(1-R[p][i-1]);//for the rest variables of the solution, chaosly assign values
		}
	}

	//to recover R from [0,1] to the real range [min,max]
	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		for (int i=0;i<solution.D;i++)
		{
			solution.IniS[p][i]=solution.Min[i]+(solution.Max[i]-solution.Min[i])*R[p][i];
			X[i]=solution.IniS[p][i];
		}
		if (solution.CECorBBOB==0)
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	}

	MatchfromIniSolutionsToEASolutions(Y);

	//return R;
	for (int p=0;p<solution.IniNP;p++)
	{
		delete []R[p];
	}
	delete []R;
	delete []X;
	delete []Y;
}

//Opposition based initialization
void Initialization::Opposition()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];

	double tempY,tempYop; //objective function value of original and its opposition solutions
	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[0],solution.Max[0]);//randomly generate the variables
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		for (int n=0;n<solution.D;n++)
		{
			X[n]=solution.Max[n]-X[n]+solution.Min[n];//get the opposition of the solution
		}
		if (solution.CECorBBOB==0)
			tempYop=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempYop=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempYop=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempYop=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		//suppose minimization problem
		if(tempY<=tempYop)//if original better, keep the original solution
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.IniS[p][n]=solution.Max[n]-X[n]+solution.Min[n];
				//since solution.X already been opposited, opposite it back to get the original solution
			}
			Y[p]=tempY;
		}
		else// if opposition result better, use the opposition solution.
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.IniS[p][n]=X[n];
			}
			Y[p]=tempYop;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);

	delete []X;
	delete []Y;
}

//Quasi opposition initialization (QOBL). Instead of fully oppositing, QOBL take a port of the opposite position.
void Initialization::QuasiOpposition()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	/*double **R=new double*[solution.Pdouble];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<solution.Pdouble;p++)
	{
		R[p]=new double[solution.D];
	}*/
	
	double *X=new double[solution.D];
	double *QX=new double[solution.D];
	double *Y=new double[solution.IniNP];

	double Middle;//record the middle point

	double tempY,tempYop; //objective function value of original and its opposition solutions
	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[0],solution.Max[0]);//randomly generate the variables
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		for (int n=0;n<solution.D;n++)
		{
			//get the opposition of the solution
			QX[n]=solution.Max[n]-X[n]+solution.Min[n];
			Middle=(solution.Min[n]+solution.Max[n])/2;

			//quasi scheme
			if(X[n]<Middle)
			{
				QX[n]=Middle+(QX[n]-Middle)*rand();
			}
			else
			{
				QX[n]=QX[n]+(Middle-QX[n])*rand();
			}	
		}
		if (solution.CECorBBOB==0)
			tempYop=CEC14.cec14_eotp_problems(QX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempYop=bbob09.FunctionCalculation(QX,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempYop=CEC13.cec13_problems(QX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempYop=CEC14_normal.cec14_problems(QX,solution.D,solution.Func_num);

		//suppose minimization problem
		if(tempY<=tempYop)//if original better, keep the original solution
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=X[n];
				//since solution.X already been opposited, opposite it back to get the original solution
			}
			Y[p]=tempY;
		}
		else// if opposition result better, use the opposition solution.
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=QX[n];
			}
			Y[p]=tempYop;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);

	delete []X;
	delete []QX;
	delete []Y;
}

//Quasi interpolation. Utilizing a simple parabola model to find better initial solutions
void Initialization::QuasiInterpolation()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	int DoubleIniNP=2*solution.IniNP;
	double **R=new double*[DoubleIniNP];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<DoubleIniNP;p++)
	{
		R[p]=new double[solution.D];
	}
	/*double **Rfinal=new double*[solution.Pdouble];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<solution.Pdouble;p++)
	{
		Rfinal[p]=new double[solution.D];
	}*/
	
	double *X=new double[solution.D];
	double *Y=new double[DoubleIniNP];
	double *rY=new double[solution.IniNP];
	
	double tempY;
	double Min=10000000000000;
	int MinCount=0;

	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[n],solution.Max[n]);//randomly generate the variables
			R[p][n]=X[n];
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
		solution.IniY[p]=tempY;
		
		//record the best random
		if(tempY<Min)
		{
			Min=tempY;
			MinCount=p;
		}
	}

	int b,c;
	double bY,cY,MinCountY;
	for (int p=solution.IniNP;p<DoubleIniNP;p++)
	{
		//find b and c different from MinCount
		b=rand()%(solution.IniNP+1);
		while(b==MinCount)
		{
			b=rand()%(solution.IniNP+1);
		}
		c=rand()%(solution.IniNP+1);
		while((c==b)||(c==MinCount))
		{
			c=rand()%(solution.IniNP+1);
		}
		//copy the objective function value
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[b][n];	
		}
		//bY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		bY=solution.IniY[b];
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[c][n];	
		}
		//cY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		cY=solution.IniY[c];
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[MinCount][n];	
		}
		//MinCountY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		MinCountY=solution.IniY[MinCount];

		//calculate the quasi interpolation points
		for (int n=0;n<solution.D;n++)
		{
			R[p][n]=0.5*(((R[b][n]*R[b][n]-R[c][n]*R[c][n])*MinCountY+(R[c][n]*R[c][n]-R[MinCount][n]*R[MinCount][n])*bY+(R[MinCount][n]*R[MinCount][n]-R[b][n]*R[b][n])*cY)
				/((R[b][n]-R[c][n])*MinCountY+(R[c][n]-R[MinCount][n])*bY+(R[MinCount][n]-R[b][n])*cY));
		}	
	}
	
	//select the best solution.P solutions
	for (int p=0;p<solution.IniNP;p++)
	{
		Y[p]=solution.IniY[p];
	}
	for (int p=solution.IniNP;p<DoubleIniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[p][n];
		}
		//The following change can hugely reduce the RAM occupied
		if (solution.CECorBBOB==0)
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	}

	int *Ycount=new int[DoubleIniNP];//to record p after sorting Y
	for (int i=0; i<solution.Pdouble; i++)
	{
		Ycount[i]=i;
	}
	//sorting Y
	for (int i=0; i<solution.Pdouble; i++)
	{  
        for (int j=solution.Pdouble-1; j>i; j--)  
        {  
			if (Y[j]<Y[j-1])
			{
				swap(Y[j], Y[j-1]);
				swap(Ycount[j],Ycount[j-1]);
			}
		}
	}

	for (int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			solution.IniS[p][n]=R[Ycount[p]][n];
		}
		rY[p]=Y[p];
	}

	MatchfromIniSolutionsToEASolutions(rY);

	for (int p=0;p<solution.Pdouble;p++)
	{
		delete []R[p];
	}
	delete []R;
	/*for (int p=0;p<solution.Pdouble;p++)
	{
		delete []Rfinal[p];
	}
	delete []Rfinal;*/
	delete []X;
	delete []Y;
	delete []rY;
	delete []Ycount;
}

//Sobol sequence, aimming to generate evenly distributed initial solutions
void Initialization::Sobol()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	int omit_N = 5; //omit the first omit_N numbers
	int generate_P = 0;
	generate_P = omit_N + 2*solution.IniNP;
	ifstream infile(solution.dir_file,ios::in);
	if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
	getchar();
    exit(1);
	}
	char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)generate_P)/log(2.0)); 
  
  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [generate_P];
  C[0] = 1;
  for (int i=1;i<=generate_P-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  double **POINTS = new double * [generate_P];
  for (int i=0;i<=generate_P-1;i++) POINTS[i] = new double [solution.D];
  for (int j=0;j<=solution.D-1;j++) POINTS[0][j] = 0; 
  
  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [generate_P];
  X[0] = 0;
  for (int i=1;i<=generate_P-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
    
  // ----- Compute the remaining dimensions -----
  for (int j=1;j<=solution.D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
	V[i] = V[i-s] ^ (V[i-s] >> s); 
	for (unsigned k=1;k<=s-1;k++) 
	  V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [generate_P];
    X[0] = 0;
	
    for (int i=1;i<=generate_P-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
	  
      POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
   }
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
    
  double *Y=new double[solution.IniNP];
  double *rX=new double[solution.D];
  for (int p=0;p<solution.Pdouble;p++)
  {
	for (int n=0;n<solution.D;n++)
	{
		solution.IniS[p][n]=POINTS[p+omit_N][n];
		rX[n]=solution.IniS[p][n];
	}
	if (solution.CECorBBOB==0)
		Y[p]=CEC14.cec14_eotp_problems(rX,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==1)
		Y[p]=bbob09.FunctionCalculation(rX,solution.Func_num);
	else if (solution.CECorBBOB==2)
		Y[p]=CEC13.cec13_problems(rX,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==3)
		Y[p]=CEC14_normal.cec14_problems(rX,solution.D,solution.Func_num);
  }

  MatchfromIniSolutionsToEASolutions(Y);

  for (int p=0;p<generate_P;p++)
  {
	delete []POINTS[p];
  }
  delete []POINTS;
}

//general function, reads the selection and calls the corresponding initialization function
void Initialization::Initial()
{
	extern Solution solution;

	switch(solution.Tech_num)
	{
	case 1:
		Random();
		break;
	case 2:
		Chaos();
		break;
	case 3:
		Opposition();
		break;
	case 4:
		QuasiOpposition();
		break;
	case 5:
		QuasiInterpolation();
		break;
	case 6:
		Sobol();
		break;
	}
}
