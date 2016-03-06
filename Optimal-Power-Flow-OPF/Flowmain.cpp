////////////////////////////////////////////////////////////////////////////////////////
//////////////////   潮流计算主程序                /////////////////////
///////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <String.h>
#include "GeneralInfo.h"
#include "ReadData.h"
#include "NodeOptimize.h"
#include "YMatrix.h"
#include "Jacobi.h"
#include "PowerFlow.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
double M_PI=3.1415926;
#define RE 10
double Dif30=0;
double Difpf30=0;

double Modulation(double V[], int a)
{
	CReadData readdata;
	time_t start130,end130;
	time_t startpf30,endpf30;
	time(&start130);
	char filename[20]="30.txt";
	//cout<<"请输入要计算的原始数据文件名(5节点对应的是5.txt,依次类推)： "<<endl;
	//cin>>filename;

	if(!strcmp(filename,"5.txt"))
	{
       strcpy(filename,"输入数据//5.txt");
	}
    else if(!strcmp(filename,"14.txt"))
	{
	    strcpy(filename,"输入数据//14.txt");
	}
	else if(!strcmp(filename,"30.txt"))
	{
		strcpy(filename,"输入数据//30.txt");
	}
	else if(!strcmp(filename,"57.txt"))
	{
		strcpy(filename,"输入数据//57.txt");
	}
	else if(!strcmp(filename,"118.txt"))
	{
		strcpy(filename,"输入数据//118.txt");
	}
	readdata.Read(filename);
	NodeOptimize Nodeoptimize;
	Nodeoptimize.OptimizeType=readdata.BusOptimType;
	if(Nodeoptimize.OptimizeType!=0)
	{
     Nodeoptimize.AdjustNodeNo();
	}
	 extern SystemInfo m_SystemInfo;   ///系统信息

	 int *new_old=new int[m_SystemInfo.SystemTotalBusNum+RE];
     for(int i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
	 {
         new_old[i]=Nodeoptimize.new_old[i];
	 }

	
//	Jacobi jacobi;	
//	jacobi.GetNewInNum();
// 	jacobi.FormJacobi();
//  jacobi.GetNonZeroNum();

	extern BusInfo *m_BusInfo;
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路
	YMatrix Ymatrix;
	Ymatrix.SpareYMatrix();

	m_BusInfo[1].Voltage_e=V[0];
	m_BusInfo[2].Voltage_e=V[1];
	m_BusInfo[5].Voltage_e=V[2];
	m_BusInfo[8].Voltage_e=V[3];
	m_BusInfo[11].Voltage_e=V[4];
	m_BusInfo[13].Voltage_e=V[5];

	m_TransformerInfo[1].Ratio=V[6];
	m_TransformerInfo[2].Ratio=V[7];
	m_TransformerInfo[3].Ratio=V[8];
	m_TransformerInfo[4].Ratio=V[9];

	m_BusInfo[2].PG=V[11];
	m_BusInfo[5].PG=V[12];
	m_BusInfo[8].PG=V[13];
	m_BusInfo[11].PG=V[14];
	m_BusInfo[13].PG=V[15];

	Ymatrix.DY[10].B=Ymatrix.DY[10].B+V[16];
	Ymatrix.DY[12].B=Ymatrix.DY[12].B+V[17];
	Ymatrix.DY[15].B=Ymatrix.DY[15].B+V[18];
	Ymatrix.DY[17].B=Ymatrix.DY[17].B+V[19];
	Ymatrix.DY[20].B=Ymatrix.DY[20].B+V[20];
	Ymatrix.DY[21].B=Ymatrix.DY[21].B+V[21];
	Ymatrix.DY[23].B=Ymatrix.DY[23].B+V[22];
	Ymatrix.DY[24].B=Ymatrix.DY[24].B+V[23];
	Ymatrix.DY[29].B=Ymatrix.DY[29].B+V[24];
	
	PowerFlow powerflow;
	powerflow.MaxIteraNum=readdata.MaxIteraNum;   ///最大迭代数目
	powerflow.IteraError=readdata.IteraError; ///最大迭代误差
	time(&startpf30);
	powerflow.PFlow();
	time(&endpf30);
	Difpf30=Difpf30+difftime(endpf30,startpf30);
	powerflow.OutputResult(filename,new_old);
	
	//double *ReturnValue= new double[ng];
	//P[0]=m_BusInfo[1].Voltage_e;
	//P[1]=m_BusInfo[2].Voltage_e;
	//P[2]=m_BusInfo[5].Voltage_e;
	//P[3]=m_BusInfo[8].Voltage_e;
	//P[4]=m_BusInfo[11].Voltage_e;
	//P[5]=m_BusInfo[13].Voltage_e;
	//P[6]=m_TransformerInfo[1].Ratio;
	//P[7]=m_TransformerInfo[2].Ratio;
    //P[8]=m_TransformerInfo[3].Ratio;
	//P[9]=m_TransformerInfo[4].Ratio;
    
//	double P[6];
//	P[0]=m_BusInfo[1].PG;
//	P[1]=m_BusInfo[2].PG;
//	P[2]=m_BusInfo[5].PG;
//	P[3]=m_BusInfo[8].PG;
//	P[4]=m_BusInfo[11].PG;
//  P[5]=m_BusInfo[13].PG;
//	for (int x=1;x<5;x++)
//  {
//		cout<<m_TransformerInfo[x].Ratio;
//	}
//	for (int j=0;j<6;j++)
//	{
//		cout<<P[j]<<endl;
//	}
	
	//return ReturnValue;
	//delete ReturnValue;
	//delete []ReturnValue;
	double c;
	c=m_BusInfo[1].PG;
	time(&end130);
	double dif130;
	dif130=difftime(end130,start130);
	Dif30=Dif30+dif130;
	return c;
	
}


double RandomGen(double min, double max)
{
       int Min = (int)(min*1000000);
       int Max = (int)(max*1000000);
       int Rand = rand()*rand();
       
       int Result = Rand%(Max-Min)+Min;
       
       return Result/1000000.0; 
}

double AverageRandom(double min, double max)
//function used in generate Gauss variable
{
       int Rand=rand()%10001;
       return (min+Rand*(max-min)/10000);
}

double Normal(double x, double miu, double sigma)
//function used in generate Gauss variable
{
       return 1.0/sqrt(2*M_PI*sigma)*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
} 

double   GaussRandom(double miu, double sigma,double min,double max)
//generate a Gauss distribute variable 
{
         double   x; 
         double   dScope; 
         double   y; 
         do 
         {
                  x = AverageRandom(min,max);   
                  y = Normal(x,miu,sigma); 
                  dScope = AverageRandom(0,Normal(miu,miu,sigma)); 
         }while(dScope>y);
         
         return x; 
}

//CRO related functions

double GetMaxConstraints(double a[],int n)
{
    double max=a[0];
    for (int i=1;i<n;i++)
    {
        if (a[i]>max)
        {
                     max=a[i];
        }
    }
    return max;
}

double GetMinConstraints(double a[],int n)
{
    double min=a[0];
    for (int i=1;i<n;i++)
    {
        if (a[i]<min)
        {
                     min=a[i];
        }
    }
    return min;
}

double SUM(double p[], int n)
{
    double sum=0;
    for (int i=0;i<n;i++)
    {
        sum=sum+p[i];
    }
    return sum;
}

//double IncreaseOutput(double p[], double cmax[], int n, double load, double sum) //n represent the nth randomly chosen elementary in the molecule 
//{
//    double min, d;
//    min=cmax[n]-p[n];
//    d=load-sum;
//    if(min>d)
//    {
//             min=d;
//    }
//    return min+p[n];
//}

//double DecreaseOutput(double p[], double cmin[], int n, double load, double sum)
//{
//    double min, d;
//    min=p[n]-cmin[n];
//    d=sum-load;
//    if(min>d)
//    {
//             min=d;
//    }
//    return p[n]-min;
//}

void Fix_molecule(double p[], double cmax[], double cmin[], int n)
//meet the boundary constraints with mirror reflects effects
{
	for(int i=0;i<n-4;i++)
     {
         double t=RandomGen(0.9,1.1);

         if (p[i]>cmax[i])   
		 {
			 if (t<1.0)
			 {
				 p[i]=cmax[i];
			 }
			 else if (t>1.0)
			 {
				 p[i]=2*cmax[i]-p[i];
			 }
         }

		 else if (p[i]<cmin[i])   
		 {
			 if (t<1.0)
			 {
				 p[i]=cmin[i];
			 }
			 else if (t>1.0)
			 {
				 p[i]=2*cmin[i]-p[i];
			 }
         }
     }
}


//void EquilityModulation(double P[], double Cmax[], double Cmin[], int NG, double TLoad)
//{
     //fit the max and min constraints for every generator. 
//    double sum;
//    sum=SUM(P,NG-2);
    //Balance the total load and output power.
//    while(sum<TLoad)
//    {
//                    int r=rand()%(NG-2);
//                    P[r]=IncreaseOutput(P, Cmax, r, TLoad, sum);
//                    sum=SUM(P,NG-2);
//    }
//    while(sum>TLoad)
//    {
//                    int v=rand()%(NG-2);
//                    P[v]=DecreaseOutput(P, Cmin, v, TLoad, sum);
//                    sum=SUM(P,NG-2);
//    }
//}

/*double CalculatePE(double b[],double c[], int conum, double p[], int ng)
{
    extern BusInfo *m_BusInfo; 
	
	double sum=0;
    for (int i=0;i<conum;i++)
    {
         sum=sum+b[i]*p[i+10]*100+c[i]*p[i+10]*p[i+10]*10000;
    }

    double Vi[31];
	for (int m=1;m<31;m++)
	{
		Vi[m]=sqrt(m_BusInfo[m].Voltage_e*m_BusInfo[m].Voltage_e+m_BusInfo[m].Voltage_f*m_BusInfo[m].Voltage_f);
		if ((m!=1)||(m!=2)||(m!=5)||(m!=8)||(m!=11)||(m!=13))
		{
			if(Vi[m]>1.1025)
			{
				sum=sum+(Vi[m]-1.1025)*(Vi[m]-1.1025)*1000;
			}
			else if(Vi[m]<0.9025)
			{
				sum=sum+(Vi[m]-0.9025)*(Vi[m]-0.9025)*1000;
			}
		}
	}

    return sum;
}*/

double CalculatePEVD(double b[],double c[], int conum, double p[], int ng)
{
    extern BusInfo *m_BusInfo; 
	
	double sum=0;
    for (int i=0;i<conum;i++)
    {
         sum=sum+b[i]*p[i+10]*100+c[i]*p[i+10]*p[i+10]*10000;
    }

    double Vi[31];
	for (int m=1;m<31;m++)
	{
		Vi[m]=sqrt(m_BusInfo[m].Voltage_e*m_BusInfo[m].Voltage_e+
			m_BusInfo[m].Voltage_f*m_BusInfo[m].Voltage_f);
		if ((m!=1)||(m!=2)||(m!=5)||(m!=8)||(m!=11)||(m!=13))
		{
			if(Vi[m]>1.05)
			{
				sum=sum+(Vi[m]-1.05)*(Vi[m]-1.05);
			}
			else if(Vi[m]<0.95)
			{
				sum=sum+(Vi[m]-0.95)*(Vi[m]-0.95);
			}
			else
			{
				sum=sum;
			}
		}		
	}

	double Qi[14];
	Qi[0]=0;
	for (int n=1;n<14;n++)
	{
		Qi[n]=m_BusInfo[n].QG;
	}

	if (Qi[1]>1.5)
	{
		sum=sum+(Qi[1]-1.5)*(Qi[1]-1.5);//*100;
	}
	else if (Qi[1]<-0.2)
	{
		sum=sum+(Qi[1]-0.2)*(Qi[1]-0.2);//*100;
	}
	else
	{
		sum=sum;
	}
	if (Qi[2]>0.6)
	{
		sum=sum+(Qi[2]-0.6)*(Qi[2]-0.6);//*100;
	}
	else if (Qi[2]<-0.2)
	{
		sum=sum+(Qi[2]-0.2)*(Qi[2]-0.2);//*100;
	}
	else
	{
		sum=sum;
	}
	if (Qi[5]>0.63)
	{
		sum=sum+(Qi[5]-0.63)*(Qi[5]-0.63);//*100;
	}
	else if (Qi[5]<-0.15)
	{
		sum=sum+(Qi[5]-0.15)*(Qi[5]-0.15);//*100;
	}
	else
	{
		sum=sum;
	}
	if (Qi[8]>0.5)
	{
		sum=sum+(Qi[8]-0.5)*(Qi[8]-0.5);//*100;
	}
	else if (Qi[8]<-0.15)
	{
		sum=sum+(Qi[8]-0.15)*(Qi[8]-0.15);//*100;
	}
	else
	{
		sum=sum;
	}
	if (Qi[11]>0.4)
	{
		sum=sum+(Qi[11]-0.4)*(Qi[11]-0.4);//*100;
	}
	else if (Qi[11]<-0.1)
	{
		sum=sum+(Qi[11]-0.1)*(Qi[11]-0.1);//*100;
	}
	else
	{
		sum=sum;
	}
	if (Qi[13]>0.45)
	{
		sum=sum+(Qi[13]-0.45)*(Qi[13]-0.45);//*100;
	}
	else if (Qi[13]<-0.15)
	{
		sum=sum+(Qi[13]-0.15)*(Qi[13]-0.15);//*100;
	}
	else
	{
		sum=sum;
	}
	return sum;
}

double *Initialize_Revise_molecule(double buf, double ke, double Cmax[], double Cmin[], int NG, double b[], double c[],int CoNum)
{
    //randomly generat the output power for every generator.
    //double max, min;
    //max=GetMaxConstraints(Cmax,NG-2);
    //min=GetMinConstraints(Cmin,NG-2);
    
    double *P= new double[NG];
        
    for (int i=0;i<NG-4;i++)
    {
        if((i==10)||(i==16)||(i==17)||(i==18)||(i==19)||(i==20)||(i==21)||(i==22||(i==23)||(i==24)))
		{
			P[i]=0;
		}
		else
		{
			P[i]=RandomGen(Cmin[i], Cmax[i]);
		}
    } 
    Fix_molecule(P,Cmax,Cmin,NG);
	P[10]=Modulation(P,NG);
	P[NG-4]=0;
	P[NG-3]=0;//number of hits and min hits.
    P[NG-2]=CalculatePEVD(b,c,CoNum,P,NG);//should be PE, later add.
    P[NG-1]=ke; 
    
    //Fix_molecule(P, Cmax, Cmin, NG);
    //EquilityModulation(P, Cmax, Cmin, NG, TLoad);
    //After reaction, the new generated solution may have invalid output power at some individual generator.
    
    return P;
    delete P;
    delete []P;
}


//functions for reactions
//generate a near neighbor solution
double *Neighbor(double b[], double c[], int conum, double temp[], double cmax[], double cmin[], int ng)
{
       double *wnew=new double[ng];
       for (int i=0;i<ng-4;i++)
       {
           if(i<16)
		   {
			   wnew[i]=temp[i]+GaussRandom(0.0,0.003, temp[i]-cmax[i], cmax[i]-temp[i]);//0.003 is the miu and can be revised.
		   }
		   else
		   {
			   wnew[i]=temp[i]+GaussRandom(0.0,0.0005, temp[i]-cmax[i], cmax[i]-temp[i]);
		   }
       }
	   Fix_molecule(wnew,cmax,cmin,ng);
	   wnew[10]=Modulation(wnew,ng);
	   wnew[ng-4]=temp[ng-4];
	   wnew[ng-3]=temp[ng-3];
       wnew[ng-2]=CalculatePEVD(b,c,conum,wnew,ng);//update of PE to be added later 
       wnew[ng-1]=temp[ng-1];
       //fit the new solution into the constraints
       //Fix_molecule(wnew, cmax, cmin, ng);
       //EquilityModulation(wnew, cmax, cmin, ng, tload);
       
       return wnew;
       delete wnew;
       delete []wnew;
}

double *Gensyn(double b[], double c[], int conum, double w1[], double w2[], 
			   double cmax[], double cmin[], int ng)
{
       double *w=new double[ng];
       for (int i=0;i<ng-4;i++)
       {
           double t = RandomGen(0.0, 1.0);
           if (t<0.49)
           {
                     w[i]=w1[i];
           }
           else
           {
                      w[i]=w2[i];
           }
       }
	   Fix_molecule(w,cmax,cmin,ng);
	   w[10]=Modulation(w,ng);
       w[ng-4]=0;
	   w[ng-3]=0;
       w[ng-2]=CalculatePEVD(b,c,conum,w,ng);
       w[ng-1]=w1[ng-2]+w1[ng-1]+w2[ng-2]+w2[ng-1]-w[ng-2];
       
       return w;
       delete w;
       delete []w;
}

bool DecomposeJudge(double w[], double w1[], double w2[], int ng, double buffer,int alpha)
{
     bool result;
     double temp1;
	 double m1 = RandomGen(0.0, 1.0);
     double m2 = RandomGen(0.0, 1.0);
     temp1 = w[ng-2]+w[ng-1]-w1[ng-2]-w2[ng-2];
     if (((temp1>=0)||(temp1+(buffer*m1*m2)>=0))&&((w[ng-4]-w[ng-3])>alpha))
     {
                    result=1;
     }
     else
     {
         result=0;
     }
     
     return result;
}

void Decompose(double w[], double w1[], double w2[], int ng, double buffer)
{
     double temp2;
	 double m1 = RandomGen(0.0, 1.0);
     double m2 = RandomGen(0.0, 1.0);
     temp2 = w[ng-2]+w[ng-1]-w1[ng-2]-w2[ng-2];
     if (temp2>=0)
     {
                    double k = RandomGen(0.0, 1.0);
                    w1[ng-1] = temp2*k;
                    w2[ng-1] = temp2*(1-k);
                    // pe1 and pe2 are remained to be set later
     }
     else if (temp2+(buffer*m1*m2)>=0) 
     {
          
          double l = RandomGen(0.0, 1.0);
          
          w1[ng-1] = (temp2+buffer*m1*m2)*l;
          w2[ng-1] = (temp2+buffer*m1*m2)*(1-l);
          buffer=buffer+temp2-w1[ng-1]-w2[ng-1];
          // pe1 and pe2 are remained to be set later
     }
	 w1[ng-4]=0;
	 w1[ng-3]=0;
	 w2[ng-4]=0;
	 w2[ng-3]=0;
}

void Onwall(double w[], double w1[], int ng, double buffer, double kelossrate)
{
     if (w[ng-2]+w[ng-1]>=w1[ng-2])
     {
                                   double q = RandomGen(kelossrate, 1.0);
                                   w1[ng-1] = (w[ng-2]+w[ng-1]-w1[ng-2])*q;
                                   buffer=buffer+(w[ng-2]+w[ng-1]-w1[ng-2])*(1-q);
                                   
                                   for(int i=0;i<ng;i++)
                                   {
                                           w[i]=w1[i];
                                   }
								   w[ng-3]=w[ng-4];
     }
}

bool Synjudge(double w[], double w1[], double w2[], int ng, double beta)
{
     bool result=0;
     double temp3;
	 temp3=w1[ng-2]+w1[ng-1]+w2[ng-2]+w2[ng-1]-w[ng-2];//-2*w[ng-2]-200;
     if ((temp3>=0)&&(w1[ng-1]<=beta)&&(w2[ng-1]<=beta)) //RCCRO
     {
                   result=1;
     }
     else
     {
         result=0;
     }
     
     return result;
}

//void Synthesis(double w[], double w1[], int ng)
//{
//    for (int i=0;i<ng;i++)
//    {
//         w[i]=w1[i];
//     }
//	 w[ng-4]=0;
//	 w[ng-3]=0;
//}

void Ineffcoll(double b[], double c[], int conum, double w1[], 
			   double w2[], double cmax[], double cmin[], int ng)
{
     double *w01, *w02;
     w01 = Neighbor(b,c,conum,w1,cmax,cmin,ng);
     w02 = Neighbor(b,c,conum,w2,cmax,cmin,ng);
	 w01[10]=Modulation(w01,ng);
	 w02[10]=Modulation(w02,ng);

     w01[ng-2]=CalculatePEVD(b,c,conum,w01,ng);
     w02[ng-2]=CalculatePEVD(b,c,conum,w02,ng);
     double temp4;
     temp4=(w1[ng-2]+w1[ng-1]+w2[ng-2]+w2[ng-1])-(w01[ng-2]+w02[ng-2]);
     
     if (temp4>=0)
     {
                   double p = RandomGen(0.0,1.0);
                   w01[ng-1] = temp4 * p;
                   w02[ng-1] = temp4 * (1-p);
                   
                   for (int i=0;i<ng;i++)
                   {
                       w1[i]=w01[i];
                       w2[i]=w02[i];
                   }
     }
	 w1[ng-3]=w1[ng-4];
	 w2[ng-3]=w2[ng-4];
	 delete w01,w02;
}
     

int main()
{
    srand((unsigned)time(NULL));
	time_t start30,end30;
	double dif30;

    //parameter initialization
    double buffer=0; 
    double InitialKE=1000;
    double KElossrate=0.2, Molecoll=0.2;
    int alpha=500;
	double beta=0.0005;
    
    int CoNum=6, NG=29;
	//CoNum: the cost coifficients for the generators
	//NG:Number of variables.
    
    //double TLoad=2.9;
    //cout<<"The total load of the system: ";
    //cin>>TLoad;
    
    //double P[18]; //output, max and min constraints for every generator. 
	double Cmax[29]={1.1,1.1,1.1,1.1,1.1,1.1,
			1.1,1.1,1.1,1.1,
			0,0.8,0.5,0.35,0.3,0.4,
			0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
			0,0,0,0	};
	double Cmin[29]={0.95,0.95,0.95,0.95,0.95,0.95,
			0.9,0.9,0.9,0.9,
			0,0.2,0.15,0.1,0.1,0.12,
			0,0,0,0,0,0,0,0,0,
			0,0,0,0};
    
    double b[6]={2.0,1.75,1.0,3.25,3.0,3.0};
    double c[6]={0.00375,0.0175,0.0625,0.00834,0.025,0.025};
    //double Qmax[14]={0,1.5,0.6,0,0,0.63,0,0,0.5,0,0,0.4,0,0.45};
	//double Qmin[14]={0,-0.2,-0.2,0,0,-0.15,0,0,-0.15,0,0,-0.1,0,-0.15};
		
	//get the constraints of every generator:
    //cout<<"Input the max and min output power for every generator(max first): "<<endl;
    
    int Ininum=5;
    //cout<<"Please input the number of molecules to be initialized: ";
    //cin>>Ininum;

    int Reactiontimes=2500;
    //cout<<"Please set the reaction times: ";
    //cin>>Reactiontimes;
    
    double *Molecule, Solution[1200][29];

    double temp5;
	time(&start30);
    for(int h=0;h<Ininum;h++) // 5 is the number of initially generated molecules.  
    {
            Molecule=Initialize_Revise_molecule(0, InitialKE, Cmax, Cmin, NG, b, c, CoNum);
            for (int g=0;g<NG;g++)
            {
                temp5=Molecule[g];
                Solution[h][g]=temp5;
                cout<<Solution[h][g]<<"  ";
            }
            cout<<endl;
			delete Molecule;
    }
    
	int demcounter=0;
    int syncounter=0;
	int onwalcounter=0;
	int intercounter=0;
	double t;
    double minimum[2500];
	int v=0;
	//double MIN=5000;

    while(Reactiontimes>0)
    {
                       t = RandomGen(0.0, 1.0);
                       if ((t>Molecoll)||(Ininum<2)) // onwall or decomposition
                       {
                                      
									  int ranini;
									  if(Ininum!=1)
									  {
										  ranini=rand()%(Ininum-1);
									  }
									  else
									  {
										  ranini=0;
									  }
									  Solution[ranini][NG-4]=Solution[ranini][NG-4]+1;
                                      double w[29];
                                      for (int p=0;p<NG;p++)
                                	  {
                                          w[p]=Solution[ranini][p];
                                      }
                                      double *w1, *w2;
                                      w1 = Neighbor(b,c,CoNum,w,Cmax,Cmin,NG);
                                      w2 = Neighbor(b,c,CoNum,w,Cmax,Cmin,NG);
                                      //w1[NG-2]=CalculatePEVD(b,c,CoNum,w1,NG);
                                      //w2[NG-2]=CalculatePEVD(b,c,CoNum,w2,NG);
									  double *ws = Neighbor(b,c,CoNum,w,Cmax,Cmin,NG);//single on wall
                                      //ws[NG-2]=CalculatePEVD(b,c,CoNum,ws,NG);
									  //if ((w1[NG-2]<=(MIN*1.01))||(w2[NG-2]<=(MIN*1.01))||(ws[NG-2]<=(MIN*1.01)))
									  //{
                                      bool judge = DecomposeJudge(w, w1, w2, NG, buffer,alpha);
                                      if ((judge==1)||(Ininum==1))
                                      {
                                                  Decompose(w,w1,w2,NG,buffer);//void function, to update the KE, PE and buffer.
                                                  for (int p=0; p<NG;p++)
                                                  {
                                                      Solution[ranini][p]=w1[p];
                                                      Solution[Ininum][p]=w2[p];
                                                  }
                                                  Ininum=Ininum+1;
                                                  demcounter=demcounter+1;
												  Reactiontimes=Reactiontimes-1;
                                      }
                                      else
                                      {
                                          
                                          Onwall(w,ws,NG,buffer,KElossrate);
                                          for (int p=0;p<NG;p++)
                                          {
                                              Solution[ranini][p]=w[p];
                                          }
										  onwalcounter++;
										  Reactiontimes=Reactiontimes-1;
                                      }

									  
                       }
                       else if((t<=Molecoll)&&(Ininum>=2))
                       {
                                      int ranini01=rand()%(Ininum-1);
                                      int ranini02=rand()%(Ininum-1);
									  Solution[ranini01][NG-4]=Solution[ranini01][NG-4]+1;
                                      Solution[ranini02][NG-4]=Solution[ranini02][NG-4]+1;
                                      double w01[29];
                                      double w02[29];
                                      for (int p=0;p<NG;p++)
                                      {
                                          w01[p]=Solution[ranini01][p];
                                      }
                                      for (int q=0;q<NG;q++)
                                      {
                                          w02[q]=Solution[ranini02][q];
                                      }
                                      //w01[NG-2]=CalculatePEVD(b,c,CoNum,w01,NG);
                                      //w02[NG-2]=CalculatePEVD(b,c,CoNum,w02,NG);
                                      double *wsyn;
                                      wsyn=Gensyn(b,c,CoNum,w01,w02,Cmax,Cmin,NG);
                                      wsyn[NG-2]=CalculatePEVD(b,c,CoNum,wsyn,NG);
									  //if(wsyn[NG-2]<=(MIN*1.01))
									  //{
                                      bool judge1=0;
									  judge1=Synjudge(wsyn,w01,w02,NG,beta);
                                      
                                      if (judge1==1)
                                      {
                                                  //synthesis, replace one of the origin molecule by the synthesis result and delete the other origin one.
                                                  for (int i=0;i<NG;i++)
                                                  {
                                                      Solution[ranini01][i]=wsyn[i];
                                                  }
                                                  //delete the w02 from the Solution matrix.
                                                  for (int p=ranini02;p<Ininum-1;p++)
                                                  {
                                                      for (int q=0;q<NG;q++)
                                                      {
                                                          Solution[p][q]=Solution[p+1][q];
                                                      }
                                                  }
                                                  Ininum=Ininum-1;
                                                  syncounter=syncounter+1;
												  Reactiontimes=Reactiontimes-1;
                                      }
                                      
                                      else
                                      {
                                          Ineffcoll(b,c,CoNum,w01,w02,Cmax,Cmin,NG);
                                          for (int i=0;i<NG;i++)
                                          {
                                              Solution[ranini01][i]=w01[i];
                                              Solution[ranini02][i]=w02[i];
                                          }
										  intercounter++;
										  Reactiontimes=Reactiontimes-1;
                                      }
									  
									  //}
									  									  
                       }
					   /*for (int u=0;u<Ininum;u++)
					   {
							if (MIN>Solution[u][27])
							{
                                  MIN=Solution[u][27];
							}
					   }*/

					   if(Ininum==0)
					   {
						   cout<<"Error!!";
						   break;
					   }

					   double Mintemp=100000;
					   for (int k=0;k<Ininum;k++)
					   {
							if (Mintemp>Solution[k][27])
							{
                                  Mintemp=Solution[k][27];
							}
					   }
					   minimum[v]=Mintemp;
					   v++;
					   

                       
    }
    
    cout<<"Remained molecule: "<<Ininum<<endl;
	cout<<"Decom: "<<demcounter<<endl;
	cout<<"Synthe: "<<syncounter<<endl;
	cout<<"Onwall: "<<onwalcounter<<endl;
	cout<<"Intermolecule: "<<intercounter<<endl; 
    
    for (int i=0;i<Ininum;i++)
    {
        for (int j=0;j<NG;j++)
        {
			if(j==0)
			{
				cout<<"V1-V13: ";
				cout<<Solution[i][j]<<"  ";
			}
			else if(j==5)
			{
				cout<<Solution[i][j]<<"  "<<endl;
			}
			else if(j==6)
			{
				cout<<"T1-T4: ";
				cout<<Solution[i][j]<<"  ";
			}
			else if(j==9)
			{
				cout<<Solution[i][j]<<"  "<<endl;
			}
			else if(j==10)
			{
				cout<<"P1-P13: ";
				cout<<Solution[i][j]<<"  ";
			}
			else if(j==15)
			{
				cout<<Solution[i][j]<<"  "<<endl;
			}
			else if(j==16)
			{
				cout<<"Q10-Q29: ";
				cout<<Solution[i][j]<<"  ";
			}
			else if(j==24)
			{
				cout<<Solution[i][j]<<"  "<<endl;
			}
			else
			{
				cout<<Solution[i][j]<<"  ";
			}
        }
        cout<<endl;
    }

    for (int r=0;r<v+1;r=r+200)
	{
		cout<<minimum[r]<<"  ";
	}
	cout<<endl<<endl;

	double mini1[2500];
	double min1=5000;
	int r2=0;
	for (int r1=0;r1<v+1;r1++)
	{
		if (min1>minimum[r1])
		{
			min1=minimum[r1];
		}
		mini1[r2]=min1;
		r2++;
	}

	cout<<"Mini-ever result"<<endl;
	for (int r3=0;r3<v+1;r3=r3+200)
	{
		cout<<mini1[r3]<<"  ";
	}
	cout<<endl<<endl;

    double Min=1000000;
    for (int l=0;l<Ininum;l++)
    {
        if (Min>Solution[l][27])
        {
                                  Min=Solution[l][27];
        }
    }
    
    cout<<Min<<endl;
    
	time(&end30);
	dif30=difftime(end30,start30);
	printf ("It took you %.2lf seconds to run the program.\n", dif30 );
	printf ("It took you %.2lf seconds to run the whole PF program.\n", Dif30 );
	printf ("It took you %.2lf seconds to run the PF function.\n", Difpf30 );

	extern 	BusInfo *m_BusInfo;
	extern 	LineInfo *m_LineInfo;
	double I_e[39];   ///电流实部 
    double I_f[39];   ///电流虚部
	double I[39];
	for(int k=1;k<=39;k++)
	{
		  int i=m_LineInfo[k].BusINo;
          int j=m_LineInfo[k].BusJNo;
		  
		  double z=m_LineInfo[k].R*m_LineInfo[k].R+m_LineInfo[k].X*m_LineInfo[k].X;
		  double Gij=m_LineInfo[k].R/z;
		  double Bij=-m_LineInfo[k].X/z;
		  I_e[k]=(m_BusInfo[i].Voltage_e-m_BusInfo[j].Voltage_e)*Gij
			     -(m_BusInfo[i].Voltage_f-m_BusInfo[j].Voltage_f)*Bij;
          I_f[k]=(m_BusInfo[i].Voltage_e-m_BusInfo[j].Voltage_e)*Bij
			     +(m_BusInfo[i].Voltage_f-m_BusInfo[j].Voltage_f)*Gij;
		  I[k]=sqrt((I_e[k]*I_e[k])+(I_f[k]*I_f[k]));
	}
	cout<<endl;
	for (int kt=1;kt<40;kt++)
	{
		cout<<kt<<"  "<<I[kt]<<endl;
	}
	return 0;
}
