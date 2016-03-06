/*
This .cpp file contains functions that read in image grey level matrix and do smoothing & filtering to improve the SNR, signal-to-noise ratio
output is directly print in the console.
*/

#include<iostream>
#include<stdlib.h> 
#include<stdio.h> 
#include<time.h>
#include<math.h>
#include<fstream>

double PI=3.1415926;
double e=2.718282;
class Molecule
{
public:
	double **MaskMatrix;
	double PE, KE;
	double dBuffer;//difference caused to the buffer
	bool Judge;//to decide whether a certain reaction will take place or not
};

using namespace std;

//to generate a gaussian random number ~N(miu,sigma)
int gaussrand(double miu, double sigma)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
	int Y;
     
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
	Y=(double)(X);
    return Y;
}

double RandomGenerator(double min, double max)
{
       int Min = (int)(min*1000000);
       int Max = (int)(max*1000000);
       int Rand = rand()*rand();
       int Result;

	   if (min!=max)
	   {
			Result = Rand%(Max-Min)+Min;
	   }
	   else
	   {
			Result=min;
	   }
       
       return Result/1000000.0; 
}

double Absolute(double input)
{
	double output;
	if (input>=0)
	{
		output=input;
	}
	else
	{
		output=-input;
	}

	return output;
}

//testing 20*20
//read the gray level matrix from the .txt file and save it into a 2D matrix
double **ReadingFile(int Length, int Height)
{
	double **Image=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Image[i]=new double[Height];
	}
    FILE *fp;
    if ((fp = fopen("Image100100.txt", "r")) == NULL)
    {   
        printf("文件Image.txt不存在\n");
		return 0;
		for(int i=0;i<Length;i++)
		{
			delete []Image[i];
		}
		delete Image;
    }   
	for(int i=0;i<Length;i++)
	{
		for(int j=0;j<Height;j++)
		{
			fscanf(fp,"%lf \n",&Image[i][j]);//reading file, %d for int, %f for float,%f for float,%lf for double
		}
		fscanf(fp,"\n");
	}
    fclose(fp);

	return Image;

	for(int i=0;i<Length;i++)
	{
		delete []Image[i];
	}
	delete Image;
}

//file name is output.txt
void WriteToFile(double **Matrix, int Length, int Height)
{
	ofstream os("output.txt");
    if (os)
    {
        for (int i=0;i<Length;++i)
        {
            for (int j=0;j<Height;++j)
				os<< Matrix[i][j]<<"\t";
            os<<endl;
        }
    }
    else
        cerr<<"无法打开文件！"<<endl;
}

void WriteToFileMask(double **Matrix, int Length, int Height)
{
	ofstream os("outputmask.txt");
    if (os)
    {
        for (int i=0;i<Length;++i)
        {
            for (int j=0;j<Height;++j)
				os<< Matrix[i][j]<<"\t";
            os<<endl;
        }
    }
    else
        cerr<<"无法打开文件！"<<endl;
}

void Displaychecking(double **Matrix, int Length, int Height)
{
	cout<<"Display the value for checking:"<<endl;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			cout<<Matrix[c][r]<<" ";
		}
		cout<<endl;
	}
	system("pause");
}

//given a matrix of gray level values, return the mean value
double GrayMatrixMean(double **Matrix,  int Length, int Height)
{
	double Mean;
	double Sum=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Sum=Sum+Matrix[c][r];
		}
	}
	Mean=Sum/(Height*Length);

	return Mean;
}

//given a matrix of gray level values, return standard deviation value
double GrayMatrixSTD(double **Matrix, int Length, int Height, double Mean)
{
	double STD;
	double Sum=0;
	Mean=GrayMatrixMean(Matrix,Height,Length);
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Sum=Sum+(Matrix[c][r]-Mean)*(Matrix[c][r]-Mean);
		}
	}
	STD=sqrt(Sum/(Length*Height));

	return STD;
}

//fix the gray values in the resultant matrix within Min and Max values
//using the equal fixing
void GrayMatrixEqualFix(double **Result,double **Image, double **Mask, int Length, int Height, double ImageMinGray, double ImageMaxGray)
{
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			if (Result[c][r]<ImageMinGray)
			{
				Mask[c][r]=ImageMinGray-Image[c][r];
				Result[c][r]=ImageMinGray;
			}
			if (Result[c][r]>ImageMaxGray)
			{
				Mask[c][r]=ImageMaxGray-Image[c][r];
				Result[c][r]=ImageMaxGray;
			}
		}
	}
}

//fix the gray values in the resultant matrix within Min and Max values
//using the bouncing-back fixing
void GrayMatrixBounceFix(double **Image, double **Mask, int Length, int Height, double ImageMinGray, double ImageMaxGray)
{
	double **AddResult=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		AddResult[i]=new double[Height];
	}

	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			AddResult[c][r]=Image[c][r]+Mask[c][r];
		}
	}

	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			if (AddResult[c][r]<ImageMinGray)
			{
				AddResult[c][r]=2*ImageMinGray-AddResult[c][r];
				Mask[c][r]=AddResult[c][r]-Image[c][r];
			}
			if (AddResult[c][r]>ImageMaxGray)
			{
				AddResult[c][r]=2*ImageMaxGray-AddResult[c][r];
				Mask[c][r]=AddResult[c][r]-Image[c][r];
			}
		}
	}

	for(int i=0;i<Length;i++)
	{
		delete []AddResult[i];
	}
	delete AddResult;

}

//given a matrix of gray level value, return the value of objective function
//assume Objective Funtion = Mean/STD considering all elements
double GrayMatrixtoFunctionValue(double **Image, double **Mask, int Length, int Height, double ImageMinGray, double ImageMaxGray)
{
	double Mean,STD,FunctionValue;
	double MeanImage;
	double **Result=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Result[i]=new double[Height];
	}
	double **In=Image;
	double **InMask=Mask;

	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Result[c][r]=In[c][r]+InMask[c][r];
		}
	}

	//GrayMatrixBounceFix(Result,Image,Mask,Length,Height,ImageMinGray,ImageMaxGray);

	Mean=GrayMatrixMean(Result,Length,Height);
	MeanImage=GrayMatrixMean(In,Length,Height);
	STD=GrayMatrixSTD(Result,Length,Height,Mean);
	double functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Result[c][r]-MeanImage)));	
		}
	}

	FunctionValue=functionresult;

	//Displaychecking(Result,Length,Height);
	//cout<<"Mean and STD check:"<<endl;
	//cout<<Mean<<" "<<STD<<endl;

	return FunctionValue;
	for(int i=0;i<Length;i++)
	{
		delete []Result[i];
		delete []In[i];
		delete []InMask[i];
	}
	delete Result;
	delete In;
	delete InMask;
}

void CopyMolecule(Molecule *From, Molecule *To, int Length, int Height)
{
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			To->MaskMatrix[c][r]=From->MaskMatrix[c][r];
		}
	}
	To->dBuffer=From->dBuffer;
	To->Judge=From->Judge;
	To->KE=From->KE;
	To->PE=From->PE;
}

//CRO functions

double SigmoidFunction(double x)
{
	double p;
	double input=(x-25)/2.86;// x from 0 to 35, input from -6 to 6
	p=1/(1+pow(e,-input));
	//p=1-p;

	return p;
}

//Initial the Molecules for CRO
Molecule *Initialization(int InitialPopulation, int InitialPE, int InitialKE, int Length, int Height,  
						 double **Image, double MeanGray, double VarianceGray, int MinGray, int MaxGray,
						 double ImageMinGray,double ImageMaxGray)
{
	//generate number of InitialPopulation molecules
	Molecule *InitialMolecules=new Molecule[InitialPopulation];
	double sum=0;//to add up the difference of mean and imagevalue result
	//set value for the parameter and control variables 
	for (int m=0;m<InitialPopulation;m++)
	{
		InitialMolecules[m].MaskMatrix=new double*[Length];
		for (int i=0;i<Length;i++)
		{
			InitialMolecules[m].MaskMatrix[i]=new double[Height];
		}
		//sum=0;
		for (int c=0;c<Length;c++)
		{
			for (int r=0;r<Height;r++)
			{
				//normal initialization
				//InitialMolecules[m].MaskMatrix[c][r]=gaussrand(MeanGray,VarianceGray);
				
				//problem specific initialization,generate random number (0,miu-Image)
				if (MeanGray-Image[c][r]>0)//Image bit less than mean
				{
					//cout<<int(MeanGray-Image[c][r])<<" "<<MeanGray-Image[c][r]<<endl;
					InitialMolecules[m].MaskMatrix[c][r]=int(RandomGenerator(0,MeanGray-Image[c][r]));
				}
				else //Image bit bigger than mean
				{
					//cout<<int(MeanGray-Image[c][r])<<" "<<MeanGray-Image[c][r]<<endl;
					InitialMolecules[m].MaskMatrix[c][r]=int(RandomGenerator(MeanGray-Image[c][r],0));
				}
			}
			//cout<<endl;
		}
		//system("pause");
		GrayMatrixBounceFix(Image,InitialMolecules[m].MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
		InitialMolecules[m].PE=GrayMatrixtoFunctionValue(Image,InitialMolecules[m].MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
		//InitialMolecules[m].PE=1000*VarianceGray;
		//InitialMolecules[m].PE=sum;
		InitialMolecules[m].KE=InitialKE;
		InitialMolecules[m].dBuffer=0;
		InitialMolecules[m].Judge=0;
	}
	
	return InitialMolecules;
	for (int m=0;m<InitialPopulation;m++)
	{
		for(int i=0;i<Length;i++)
		{
			delete []InitialMolecules[m].MaskMatrix[i];
		}
		delete InitialMolecules[m].MaskMatrix;
	}
	delete InitialMolecules;
	delete []InitialMolecules;
	
}

//Neighbor search, given a 2D graylevel matrix, return a neighborhood matrix
double **NeighborSearch(double **Matrix, int Length, int Height, double **Image,
						double ImageMinGray, double ImageMaxGray,
						double VarianceGray, double MinGray, double MaxGray)
{
	double **Neighbor=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Neighbor[i]=new double[Height];
	}

	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			//parameters pending for changes
			//Neighbor[c][r]=NormalRandom(Matrix[c][r],VarianceGray,MinGray,MaxGray);
			Neighbor[c][r]=gaussrand(Matrix[c][r],VarianceGray);
		}
	}

	//Neighbor updated according to the Image gray level constraints
	GrayMatrixBounceFix(Image,Neighbor,Length,Height,ImageMinGray,ImageMaxGray);

	return Neighbor;
	for(int i=0;i<Length;i++)
	{
		delete []Neighbor[i];
	}
	delete Neighbor;
}

//Onwall collision, given a Molecule class, return the result of onwall collision
void OnWall(Molecule *CROmolecule, Molecule *Result, int Length, int Height, double Buffer, double KElossrate,
				double **Image, double ImageMinGray, double ImageMaxGray, double VarianceGray,double MinGray, double MaxGray)
{
	Molecule *In=CROmolecule;
	Molecule *Out=Result;
	CopyMolecule(In,Out,Length,Height);
	double TempCheck;
	double q;
	double **Result1=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Result1[i]=new double[Height];
	}
	//double sum=0;
	//Out->MaskMatrix=NeighborSearch(In->MaskMatrix,Length,Height,Image,ImageMinGray,ImageMaxGray,
		//VarianceGray,MinGray,MaxGray);
	//Out->PE=GrayMatrixtoFunctionValue(Image,Out->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	double MeanImage=GrayMatrixMean(Image,Length,Height);
	double p, randomp;//p is the probability to conduct neighbor search, randomp is to judge whether happen or not
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Result1[c][r]=Image[c][r]+In->MaskMatrix[c][r];
			//p=SigmoidFunction(Absolute(MeanImage-Image[c][r]));//deterministic
			p=SigmoidFunction(Absolute(MeanImage-Result1[c][r]));//dynamic
			randomp=RandomGenerator(0,1);
			if (randomp<=p)
			{
				Out->MaskMatrix[c][r]=gaussrand(Out->MaskMatrix[c][r],VarianceGray);
			}
			Result1[c][r]=Image[c][r]+Out->MaskMatrix[c][r];
		}
	}
	//cout<<"judge check during onwall before tempcheck: "<<Out->Judge<<endl;
	GrayMatrixBounceFix(Image,Out->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	double Mean,STD;
	Mean=GrayMatrixMean(Result1,Length,Height);
	MeanImage=GrayMatrixMean(Image,Length,Height);
	STD=GrayMatrixSTD(Result1,Length,Height,Mean);
	
	double functionresult=0;
	//distinguish formulation
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Result1[c][r]-MeanImage)));	
		}
	}

	//normal formulation
	//functionresult=-Mean/STD+100*Absolute(MeanImage-Mean);
	Out->PE=functionresult;

	//Out->PE=GrayMatrixtoFunctionValue(Image,Out->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out->PE=1000*Absolute(Mean-MeanImage)+1000*Absolute(STD);
	//Out->PE=1000*STD;
	//Out->PE=sum;
	Out->dBuffer=0;
	//check
	TempCheck=In->PE+In->KE-Out->PE;
	//cout<<"In PE: "<<In->PE<<" In KE: "<<In->KE<<endl;
	//cout<<"Out PE: "<<Out->PE<<endl;
	//cout<<"TempCheck: "<<TempCheck<<endl;
	//system("pause");
	if (TempCheck>0)
	{
		q=RandomGenerator(KElossrate,1);
		Out->KE=TempCheck*q;
		//update the molecular structure and return
		Out->dBuffer=TempCheck*(1-q);//require buffer updating in main
		Out->Judge=1;
	}
	else
	{
		Out->Judge=0;
	}
	CopyMolecule(Out,Result,Length,Height);
	//cout<<"judge check during onwall after tempcheck: "<<Out->Judge<<endl;
	for(int i=0;i<Length;i++)
	{
		delete []Result1[i];
	}
	delete Result1;
	//delete []In;
	//delete []Out;
}

//InterMolecular collision, use void to update the given two molecular structures
void InterMolecule(Molecule *CROmolecule1, Molecule *CROmolecule2, Molecule *Result1, Molecule *Result2, int Length, int Height,
				   double **Image, double ImageMinGray, double ImageMaxGray,double VarianceGray,double MinGray,double MaxGray)
{
	Molecule *In1=CROmolecule1, *In2=CROmolecule2;
	Molecule *Out1=Result1, *Out2=Result2;
	CopyMolecule(In1,Out1,Length,Height);
	CopyMolecule(In2,Out2,Length,Height);
	double TempCheck;
	//double sum=0;
	double **Resultout1=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Resultout1[i]=new double[Height];
	}
	double **Resultout2=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Resultout2[i]=new double[Height];
	}
	double MeanImage=GrayMatrixMean(Image,Length,Height);
	double p, randomp;//p is the probability to conduct neighbor search, randomp is to judge whether happen or not
	for (int c=0;c<Height;c++)
	{
		for (int r=0;r<Length;r++)
		{
			Resultout1[c][r]=Image[c][r]+In1->MaskMatrix[c][r];
			//p=SigmoidFunction(Absolute(MeanImage-Image[c][r]));//deterministic
			p=SigmoidFunction(Absolute(MeanImage-Resultout1[c][r]));//dynamic
			randomp=RandomGenerator(0,1);
			if (randomp<=p)
			{
				Out1->MaskMatrix[c][r]=gaussrand(Out1->MaskMatrix[c][r],VarianceGray);
			}
			Resultout1[c][r]=Image[c][r]+Out1->MaskMatrix[c][r];
			//sum=sum+Absolute(64-Resultout1[c][r]);
		}
	}
	GrayMatrixBounceFix(Image,Out1->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	double Mean1,STD1;
	Mean1=GrayMatrixMean(Resultout1,Length,Height);
	STD1=GrayMatrixSTD(Resultout1,Length,Height,Mean1);

	double functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Resultout1[c][r]-MeanImage)));	
		}
	}
	Out1->PE=functionresult;

	//Out1->PE=GrayMatrixtoFunctionValue(Image,Out1->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out1->PE=1000*Absolute(Mean1-MeanImage)+1000*Absolute(STD1);
	//Out1->PE=1000*STD1;
	//Out1->PE=sum;
	Out1->dBuffer=0;

	//sum=0;
	for (int c=0;c<Height;c++)
	{
		for (int r=0;r<Length;r++)
		{
			Resultout2[c][r]=Image[c][r]+In2->MaskMatrix[c][r];
			p=SigmoidFunction(Absolute(MeanImage-Image[c][r]));//deterministic
			//p=SigmoidFunction(Absolute(MeanImage-Resultout2[c][r]));//dynamic
			randomp=RandomGenerator(0,1);
			if (randomp<=p)
			{
				Out2->MaskMatrix[c][r]=gaussrand(Out2->MaskMatrix[c][r],VarianceGray);
			}
			Resultout2[c][r]=Image[c][r]+Out2->MaskMatrix[c][r];
			//sum=sum+Absolute(64-Resultout2[c][r]);
		}
	}
	GrayMatrixBounceFix(Image,Out2->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	double Mean2,STD2;
	Mean2=GrayMatrixMean(Resultout2,Length,Height);
	STD2=GrayMatrixSTD(Resultout2,Length,Height,Mean2);
	
	functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Resultout2[c][r]-MeanImage)));	
		}
	}
	Out2->PE=functionresult;
	
	//Out2->PE=GrayMatrixtoFunctionValue(Image,Out2->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out2->PE=1000*Absolute(Mean2-MeanImage)+1000*Absolute(STD2);
	//Out2->PE=1000*STD2;
	//Out2->PE=sum;
	Out2->dBuffer=0;

	TempCheck=(In1->PE+In1->KE+In2->PE+In2->KE)-Out1->PE-Out2->PE;
	double CROp;
	if (TempCheck>=0)
	{
		CROp=RandomGenerator(0,1);
		Out1->KE=TempCheck*CROp;
		Out2->KE=TempCheck*(1-CROp);
		Out1->Judge=1;
		Out2->Judge=1;
		//update molecular structure
		//CopyMolecule(Out1,Result1,Length,Height);
		//CopyMolecule(Out2,Result2,Length,Height);
	}
	else
	{
		Out1->Judge=0;
		Out2->Judge=0;
	}
	CopyMolecule(Out1,Result1,Length,Height);
	CopyMolecule(Out2,Result2,Length,Height);
	//delete In1, In2;
	//delete Out1, Out2;
	for(int i=0;i<Length;i++)
	{
		delete []Resultout1[i];
	}
	delete Resultout1;
	for(int i=0;i<Length;i++)
	{
		delete []Resultout2[i];
	}
	delete Resultout2;
}

//Decomposition, decide whether a decomposition will happen and update parameters
void Decomposition(Molecule *CROmolecule, Molecule *CROmolecule1, Molecule *CROmolecule2, 
				   int Length, int Height, double Buffer, double KElossrate,
				double **Image, double ImageMinGray, double ImageMaxGray,
				double MinGray, double MaxGray, double MeanGray, double VarianceGray)
{
	double TempCheck;
	double k,m1,m2,m3,m4;
	m1=RandomGenerator(0,1);
	m2=RandomGenerator(0,1);
	m3=RandomGenerator(0,1);
	m4=RandomGenerator(0,1);

	Molecule *In=CROmolecule;
	Molecule *Out1=CROmolecule1, *Out2=CROmolecule2;
	//CopyMolecule(In,Out1,Length,Height);
	//CopyMolecule(In,Out2,Length,Height);
	double **Result1=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Result1[i]=new double[Height];
	}
	double **Result2=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Result2[i]=new double[Height];
	}
	
	for (int c=0;c<Length/2;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Out1->MaskMatrix[c][r]=In->MaskMatrix[c][r];
			Out2->MaskMatrix[c][r]=gaussrand(0,VarianceGray);//a new gaussand random around 0
			Out1->MaskMatrix[c+Length/2][r]=gaussrand(0,VarianceGray);//a new gaussand random around 0
			Out2->MaskMatrix[c+Length/2][r]=In->MaskMatrix[c+Length/2][r];
		}
	}

	GrayMatrixBounceFix(Image,Out1->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	GrayMatrixBounceFix(Image,Out2->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	for (int c=0;c<Length/2;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Result1[c][r]=Out1->MaskMatrix[c][r]+Image[c][r];
			Result2[c][r]=Out2->MaskMatrix[c][r]+Image[c][r];
		}
	}

	double Mean1,Mean2,MeanImage,STD1,STD2;
	MeanImage=GrayMatrixMean(Image,Length,Height);
	Mean1=GrayMatrixMean(Result1,Length,Height);
	Mean2=GrayMatrixMean(Result2,Length,Height);
	STD1=GrayMatrixSTD(Result1,Length,Height,Mean1);
	STD2=GrayMatrixSTD(Result2,Length,Height,Mean2);

	double functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Result1[c][r]-MeanImage)));	
		}
	}
	Out1->PE=functionresult;

	functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Result2[c][r]-MeanImage)));	
		}
	}
	Out2->PE=functionresult;

	//Out1->PE=GrayMatrixtoFunctionValue(Image,Out1->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out2->PE=GrayMatrixtoFunctionValue(Image,Out2->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out1->PE=1000*Absolute(Mean1-MeanImage)+1000*Absolute(STD1);
	//Out2->PE=1000*Absolute(Mean2-MeanImage)+1000*Absolute(STD2);
	//Out1->PE=1000*STD1;
	//Out2->PE=1000*STD2;

	TempCheck=In->KE+In->PE-Out1->PE-Out2->PE;

	//cout<<"In PE: "<<In->PE<<" In KE: "<<In->KE<<endl;
	//cout<<"Out1 PE: "<<Out1->PE<<endl;
	//cout<<"Out2 PE: "<<Out2->PE<<endl;
	//cout<<"buffer: "<<Buffer<<endl;
	//cout<<"TempCheck: "<<TempCheck<<endl;
	//system("pause");

	if (TempCheck>0)
	{
		k=RandomGenerator(0,1);
		Out1->KE=TempCheck*k;
		Out2->KE=TempCheck*(1-k);
		Out1->dBuffer=0;
		Out1->Judge=1;
		Out2->dBuffer=0;
		Out2->Judge=1;
		//CopyMolecule(Out1,CROmolecule1,Length,Height);
		//CopyMolecule(Out2,CROmolecule2,Length,Height);
	}
	else if (TempCheck+Buffer>=0)
	{
		Out1->KE=(TempCheck+Buffer)*m1*m2;
		Out2->KE=(TempCheck+Buffer-Out1->KE)*m3*m4;
		Out1->dBuffer=TempCheck-Out1->KE-Out2->KE;//update buffer in main
		Out1->Judge=1;
		Out2->dBuffer=TempCheck-Out1->KE-Out2->KE;//update buffer in main
		Out2->Judge=1;
		//CopyMolecule(Out1,CROmolecule1,Length,Height);
		//CopyMolecule(Out2,CROmolecule2,Length,Height);
	}
	else
	{
		Out1->Judge=0;
		Out2->Judge=0;
	}
	CopyMolecule(Out1,CROmolecule1,Length,Height);
	CopyMolecule(Out2,CROmolecule2,Length,Height);
	//delete In;
	//delete Out1, Out2;
	for(int i=0;i<Length;i++)
	{
		delete []Result1[i];
	}
	delete Result1;
	for(int i=0;i<Length;i++)
	{
		delete []Result2[i];
	}
	delete Result2;
}

//Synthesis, given two molecules, return one combined result
void Synthesis(Molecule *CROmolecule1, Molecule *CROmolecule2, Molecule *Result, int Length, int Height,
				   double **Image, double ImageMinGray, double ImageMaxGray)
{
	Molecule *In1=CROmolecule1,*In2=CROmolecule2;
	Molecule *Out=Result;
	double **Resultout1=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Resultout1[i]=new double[Height];
	}
	//Generate the synthesis molecular structure
	for (int c=0;c<Length/2;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Out->MaskMatrix[c][r]=In1->MaskMatrix[c][r];
			Out->MaskMatrix[c+Length/2][r]=In2->MaskMatrix[c+Length/2][r];
			Resultout1[c][r]=Out->MaskMatrix[c][r]+Image[c][r];
		}
	}
	GrayMatrixBounceFix(Image,Out->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	double Mean1,MeanImage,STD1;
	Mean1=GrayMatrixMean(Resultout1,Length,Height);
	MeanImage=GrayMatrixMean(Image,Length,Height);
	STD1=GrayMatrixSTD(Resultout1,Length,Height,Mean1);

	double functionresult=0;
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			functionresult=functionresult+pow(e,(Absolute(Resultout1[c][r]-MeanImage)));	
		}
	}
	Out->PE=functionresult;

	//Out->PE=GrayMatrixtoFunctionValue(Image,Out->MaskMatrix,Length,Height,ImageMinGray,ImageMaxGray);
	//Out->PE=1000*Absolute(Mean1-MeanImage)+1000*Absolute(STD1);
	//Out->PE=1000*STD1;
	Out->dBuffer=0;

	double TempCheck;
	TempCheck=In1->PE+In1->KE+In2->PE+In2->KE-Out->PE;

	if (TempCheck>=0)
	{
		Out->KE=TempCheck;
		Out->Judge=1;
	}
	else
	{
		Out->Judge=0;
	}
	CopyMolecule(Out,Result,Length,Height);
	//delete In1, In2;
	//delete Out;
	for(int i=0;i<Length;i++)
	{
		delete []Resultout1[i];
	}
	delete Resultout1;
}

int main()
{
	srand(time(NULL));

	int Length=100;//Length should be a even number otherwise sth wrong in Decomposition&Synthesis
	int Height=100;//imaging matrix dimension
	double ImageMinGray=27;
	double ImageMaxGray=90;//range of gray level value
	double MinGray=0;
	double MaxGray=30;//Range of gray level value for mask matrix
	double VarianceGray=2;//Mean and Variance for generating Random Mask

	//define a dynamic 2-dimension matrix with Length*Height elements
	//Read the file value into the Image matrix
	double **Image=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Image[i]=new double[Height];
	}
	Image=ReadingFile(Length,Height);
	double MeanGray=GrayMatrixMean(Image,Length,Height);

	//CRO parameter Initialization
	int InitialPopulation=1;
	int Population=InitialPopulation;
	int MaxFE=50000;
	int FE=0;
	int MoleculePick=0, MoleculePick1=0;//for randomly picked up molecules from the population
	double InitialKE=1000;
	double InitialPE=0;
	double Buffer=0;
	double KElossrate=0.01;
	double InterMolecularCollisionRate=0,t;//t to judge whether a single or multiple molecular collision will happen
	int start=clock();
	int concounter=0;
	double Converge[500];

	int looptime=1;
	int loop=looptime;//for output
	double *Loopresult=new double[looptime];

	while (looptime>0)
	{
		FE=0;
		Population=InitialPopulation;
		Buffer=0;
	//Initial molecules
	Molecule *CROPopulation=Initialization(10,InitialPE,InitialKE,Length,Height,
		Image,MeanGray,VarianceGray,MinGray,MaxGray,ImageMinGray,ImageMaxGray);
	//to temporarily store information in the optimization process
	Molecule *Temp=Initialization(2,InitialPE,InitialKE,Length,Height,Image,MeanGray,VarianceGray,MinGray,MaxGray,ImageMinGray,ImageMaxGray);
	//CRO loop
	int OnwallCounter=0;
	int DecompositionCounter=0;
	int InterMoleculeCounter=0;
	int SynthesisCounter=0;
	int Globalcounter=0;
	double GlobalBest=100000;
	double **GlobalMatrix=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		GlobalMatrix[i]=new double[Height];
	}

	while(FE<MaxFE)
	{
		t=RandomGenerator(0,1);
		//reset Temp
		for (int i=0;i<2;i++)
		{
			Temp[i].dBuffer=0;
			Temp[i].Judge=0;
			Temp[i].KE=0;
			Temp[i].PE=0;
			for (int c=0;c<Length;c++)
			{
				for (int r=0;r<Height;r++)
				{
					Temp[i].MaskMatrix[c][r]=0;
				}
			}
		}
		
		if ((t>InterMolecularCollisionRate)||(Population==1))//single molecular reactions are considered
		{
			if (Population==1)
			{
				MoleculePick=0;
			}
			else
			{
				MoleculePick=rand()%(Population-1);
			}
			
			//see if a Decomposition can happen
			//Temp=Initialization(2,InitialPE,InitialKE,Length,Height,Image,
				//MeanGray,VarianceGray,MinGray,MaxGray,ImageMinGray,ImageMaxGray);
			Temp[0].Judge=0;
			//Decomposition(&CROPopulation[MoleculePick],&Temp[0],&Temp[1],Length,Height,Buffer,KElossrate,
				//Image,ImageMinGray,ImageMaxGray,MinGray,MaxGray,MeanGray,VarianceGray);
			
			if (Temp[0].Judge==1)//a decomposition reaction happened
			{
				//cout<<"Decomposition happens"<<endl;
				//cout<<"dbuffer: "<<Temp[0].dBuffer<<endl;
				Buffer=Buffer+Temp[0].dBuffer;
				FE=FE+2;
				Population++;

				//update solution space
				CopyMolecule(&Temp[0],&CROPopulation[MoleculePick],Length,Height);
				CopyMolecule(&Temp[1],&CROPopulation[Population-1],Length,Height);
				DecompositionCounter++;
			}
			
			else
			{				
				//Temp=Initialization(2,InitialPE,InitialKE,Length,Height,Image,
					//MeanGray,VarianceGray,MinGray,MaxGray,ImageMinGray,ImageMaxGray);
				Temp[0].Judge=0;
				OnWall(&CROPopulation[MoleculePick],&Temp[0],Length,Height,Buffer,KElossrate,
					Image,ImageMinGray,ImageMaxGray,VarianceGray,MinGray,MaxGray);
				if (Temp[0].Judge==1)//a onwall reaction happened
				{
					//cout<<"Onwall happens"<<endl;
					//cout<<"dbuffer: "<<Temp[0].dBuffer<<endl;
					Buffer=Buffer+Temp[0].dBuffer;//update buffer
					FE=FE+1;
					CopyMolecule(&Temp[0],&CROPopulation[MoleculePick],Length,Height);//replace the solution, update the solutionspace
					OnwallCounter++;
				}
				else
				{
					//cout<<"No single reaction happens"<<endl;
					FE=FE+1;
				}
				//system("pause");
			}
		}
		
		else//multiple molecular reactions are considered
		{
			MoleculePick=rand()%(Population-1);
			//do
			//{
				MoleculePick1=rand()%(Population-1);
			//}while(MoleculePick1==MoleculePick);
			Temp[0].Judge=0;
			//Synthesis(&CROPopulation[MoleculePick],&CROPopulation[MoleculePick1],&Temp[0],Length,Height,Image,ImageMinGray,ImageMaxGray);
			
			if (Temp[0].Judge==1)//a synthesis reaction happened
			{
				//cout<<"Synthesis happens"<<endl;
				FE=FE+1;
				//update solution space
				if (MoleculePick>=MoleculePick1)
				{
					CopyMolecule(&Temp[0],&CROPopulation[MoleculePick],Length,Height);
					//CROPopulation[MoleculePick]=Temp[0];
					for (int i=MoleculePick1+1;i<Population;i++)
					{
						CROPopulation[i-1]=CROPopulation[i];
					}
				}
				else
				{
					CopyMolecule(&Temp[0],&CROPopulation[MoleculePick1],Length,Height);
					//CROPopulation[MoleculePick]=Temp[0];
					for (int i=MoleculePick+1;i<Population;i++)
					{
						CROPopulation[i-1]=CROPopulation[i];
					}
				}
				Population--;
				SynthesisCounter++;
			}

			else
			{
				Temp[0].Judge=0;
				//InterMolecule(&CROPopulation[MoleculePick],&CROPopulation[MoleculePick1],&Temp[0],&Temp[1],Length,Height,Image,ImageMinGray,ImageMaxGray,
						//VarianceGray,MinGray,MaxGray);
				if (Temp[0].Judge==1)
				{
					//cout<<"InterMolecular Collision happens"<<endl;
					CopyMolecule(&Temp[0],&CROPopulation[MoleculePick],Length,Height);
					CopyMolecule(&Temp[1],&CROPopulation[MoleculePick1],Length,Height);
					FE=FE+2;
					InterMoleculeCounter++;
				}
				else
				{
					FE=FE+2;
				}
			}
		}
		
		for (int i=0;i<Population;i++)
		{
			if (CROPopulation[i].PE<GlobalBest)
			{
				GlobalBest=CROPopulation[i].PE;
				Globalcounter=i;
			}
		}

		double **Convergematrix=new double*[Length];
		for (int i=0;i<Length;i++)
		{
			Convergematrix[i]=new double[Height];
		}
		for (int c=0;c<Height;c++)
		{
			for (int r=0;r<Length;r++)
			{
				GlobalMatrix[c][r]=CROPopulation[Globalcounter].MaskMatrix[c][r];
				Convergematrix[c][r]=Image[c][r]+GlobalMatrix[c][r];
			}
		}

		double Meanconverge=GrayMatrixMean(Convergematrix,Length,Height);
		double STDconverge=GrayMatrixSTD(Convergematrix,Length,Height,Meanconverge);
		if (FE%5000==0)
		{
			cout<<"FE="<<FE<<endl;
			Converge[concounter]=STDconverge;
			concounter++;
		}
		for(int i=0;i<Length;i++)
		{
			delete []Convergematrix[i];
		}
		delete Convergematrix;

	}//end of iterations
	int end=clock();
	int Runtime=(end-start)/ CLOCKS_PER_SEC;

	cout<<"Converge: "<<endl;
	for (int i=0;i<concounter;i++)
	{
		cout<<Converge[i]<<" ";
	}
	cout<<endl;
	//output result
	
	double Min=10000000;
	int Record;
	/*for (int i=0;i<Population;i++)
	{
		if (CROPopulation[i].PE<Min)
		{
			Min=CROPopulation[i].PE;
			Record=i;
		}
	}*/

	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			//cout<<CROPopulation[Record].MaskMatrix[c][r]<<" ";
			//cout<<GlobalMatrix[c][r]<<" ";
		}
		//cout<<endl;
	}

	double **Result=new double*[Length];
	for (int i=0;i<Length;i++)
	{
		Result[i]=new double[Height];
	}
	for (int c=0;c<Length;c++)
	{
		for (int r=0;r<Height;r++)
		{
			Result[c][r]=GlobalMatrix[c][r]+Image[c][r];
			//cout<<GlobalMatrix[c][r]+Image[c][r]<<" ";
		}
		//cout<<endl;
	}

	double MeanImage, MeanResult, VarianceImage,VarianceResult;
	MeanImage=GrayMatrixMean(Image,Length,Height);
	MeanResult=GrayMatrixMean(Result,Length,Height);
	VarianceImage=GrayMatrixSTD(Image,Length,Height,MeanImage);
	VarianceResult=GrayMatrixSTD(Result,Length,Height,MeanResult);
	cout<<"Mean check:"<<endl;
	cout<<MeanImage<<" "<<MeanResult;
	cout<<endl<<"STD result:"<<endl;
	cout<<VarianceImage<<" "<<VarianceResult;
	cout<<endl<<"best PE ever:"<<endl;
	cout<<CROPopulation[Globalcounter].PE<<endl;
	
	cout<<endl;
	cout<<"Run Time: "<<Runtime<<endl;
	cout<<"Onwall Counter: "<<OnwallCounter<<endl;
	cout<<"Decomposition Counter: "<<DecompositionCounter<<endl;
	cout<<"InterMolecule Counter: "<<InterMoleculeCounter<<endl;
	cout<<"Synthesis Counter: "<<SynthesisCounter<<endl;

	//output matrix to file
	WriteToFile(Result,Length,Height);
	WriteToFileMask(CROPopulation[Globalcounter].MaskMatrix,Length,Height);

	for(int i=0;i<Length;i++)
	{
		delete []Result[i];
	}
	delete Result;

	Loopresult[looptime-1]=VarianceResult;
	cout<<"Current in the "<<looptime<<"th loop"<<endl;
	looptime--;
	}//end of the loop

	double loopmin=10000000;
	double loopmax=0;
	double loopsum=0;
	for (int i=0;i<loop;i++)
	{
		if (loopmin>Loopresult[i])
		{
			loopmin=Loopresult[i];
		}
		if (loopmax<Loopresult[i])
		{
			loopmax=Loopresult[i];
		}
		loopsum+=Loopresult[i];
	}
	cout<<"Loop Min result: "<<loopmin<<endl;
	cout<<"Loop Max result: "<<loopmax<<endl;
	cout<<"Loop Mean result: "<<loopsum/loop<<endl;

	//release the RAM
	for(int i=0;i<Length;i++)
	{
		delete []Image[i];
	}
	delete Image;
	//delete CROPopulation; //if included, will cause Heap Block Problem
	//delete Temp;
	
	system("pause");
    return 0;
}
