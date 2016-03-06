////////////////////////////////////////////////////////////////////////////////////////
//////////////////   潮流计算主程序                /////////////////////
///////////////////////////////////////////////////////////////////////////////////////
#include "iostream.h" 
#include <String.h>
#include "GeneralInfo.h"
#include "ReadData.h"
#include "NodeOptimize.h"
#include "YMatrix.h"
#include "Jacobi.h"
#include "PowerFlow.h"
#define RE 10

void Modulation(double P[], int ng)
{
	CReadData readdata;
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

//	YMatrix Ymatrix;
// 	Ymatrix.SpareYMatrix();
//	Jacobi jacobi;	
//	jacobi.GetNewInNum();
// 	jacobi.FormJacobi();
//  jacobi.GetNonZeroNum();

	extern BusInfo *m_BusInfo;
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路

	m_BusInfo[2].Voltage_e=P[1];
	m_BusInfo[5].Voltage_e=P[2];
	m_BusInfo[8].Voltage_e=P[3];
	m_BusInfo[11].Voltage_e=P[4];
	m_BusInfo[13].Voltage_e=P[5];

	m_BusInfo[2].PG=P[11];
	m_BusInfo[5].PG=P[12];
	m_BusInfo[8].PG=P[13];
	m_BusInfo[11].PG=P[14];
	m_BusInfo[13].PG=P[15];

    m_TransformerInfo[1].Ratio=P[6];
	m_TransformerInfo[2].Ratio=P[7];
	m_TransformerInfo[3].Ratio=P[8];
	m_TransformerInfo[4].Ratio=P[9];

	PowerFlow powerflow;
	powerflow.MaxIteraNum=readdata.MaxIteraNum;   ///最大迭代数目
	powerflow.IteraError=readdata.IteraError; ///最大迭代误差
	powerflow.PFlow();
	//powerflow.OutputResult(filename,new_old);
	
	//double *ReturnValue= new double[ng];
	
	P[0]=m_BusInfo[1].Voltage_e;
	P[1]=m_BusInfo[2].Voltage_e;
	P[2]=m_BusInfo[5].Voltage_e;
	P[3]=m_BusInfo[8].Voltage_e;
	P[4]=m_BusInfo[11].Voltage_e;
	P[5]=m_BusInfo[13].Voltage_e;
	P[6]=m_TransformerInfo[1].Ratio;
	P[7]=m_TransformerInfo[2].Ratio;
	P[8]=m_TransformerInfo[3].Ratio;
	P[9]=m_TransformerInfo[4].Ratio;
	P[10]=m_BusInfo[1].PG;
	P[11]=m_BusInfo[2].PG;
	P[12]=m_BusInfo[5].PG;
	P[13]=m_BusInfo[8].PG;
	P[14]=m_BusInfo[11].PG;
	P[15]=m_BusInfo[13].PG;

	cout<<P[10]<<m_BusInfo[1].PG<<"  "<<endl;

	//return ReturnValue;
	//delete ReturnValue;
	//delete []ReturnValue;

}