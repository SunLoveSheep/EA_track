/////////////////////////////////////////////////////////////////////////////
//////      Read input data from "输入数据" (inside "InputData.rar")     ////
/////                                                                    ////
/////////////////////////////////////////////////////////////////////////////

#include "ReadData.h"
#include "GeneralInfo.h"
#include <vector>
#define RE 10
using namespace std;
CReadData::CReadData()
{

}
CReadData::~CReadData()
{

}
void CReadData::Read(char *Filename)
{
	extern SystemInfo m_SystemInfo;   ///系统信息
	extern LineInfo *m_LineInfo;       ///线路信息
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路 
	extern BusInfo *m_BusInfo;        //// 节点
	/////因为每组数据的长度不知道所以使用容器读取
   vector<LineInfo> Line;
   vector<TransformerInfo> transformer;
   vector<BusInfo> generator;
   vector<BusInfo> load;
   char End[20];   ///读入结束标志
   datain.open(Filename);
   datain>>BusOptimType>>MaxIteraNum>>IteraError>>BasePower>>BaseVoltage;
   while(1) ////读入线路数据
	{ 
	   LineInfo temp;
       datain>>End;
	   if(!strcmp(End,"end")) break;
	   temp.BusINo=atoi(End);
	   datain>>temp.BusJNo>>temp.R>>temp.X>>temp.Yk;	  
	   Line.push_back(temp);
	}  
	while(1)  ///读入变压器数据
	{
      TransformerInfo temp;
       datain>>End;
	   if(!strcmp(End,"end")) break;
	   temp.BusINo=atoi(End);
	   datain>>temp.BusJNo>>temp.R>>temp.X>>temp.Ratio;
	   transformer.push_back(temp);
	}	
   while(1)   ///读入发电机节点数据
   {
      BusInfo temp;
	  datain>>End;
	    if(!strcmp(End,"end")) break;
	   temp.BusNo=atoi(End);
	   datain>>temp.BusType>>temp.PG>>temp.QG>>temp.Voltage_e;
	   generator.push_back(temp);
   }
   while(1)   ///读入负荷节点数据
   {
      BusInfo temp;
	   datain>>End;
	   if(!strcmp(End,"end")) break;
	   temp.BusNo=atoi(End);
	   datain>>temp.BusType>>temp.PD>>temp.QD>>temp.Voltage_e;
	   load.push_back(temp);
   }
   int m_num;   //统计每组数据的长度
   ////从容器中数据
   ////读入线路数据
   m_num=Line.size()+RE;
   m_LineInfo=new LineInfo[m_num];
   m_SystemInfo.GeneralLineNum=0; m_SystemInfo.LandBranchNum=0;
   for(int i=1;i<=Line.size();i++)
   {
	   m_LineInfo[i].BusINo=Line[i-1].BusINo;
	   m_LineInfo[i].BusJNo=Line[i-1].BusJNo;
       m_LineInfo[i].R=Line[i-1].R;
	   m_LineInfo[i].X=Line[i-1].X;
	   m_LineInfo[i].Yk=Line[i-1].Yk;
	   if(Line[i-1].BusINo==Line[i-1].BusJNo) m_SystemInfo.LandBranchNum+=1;
   }
   m_SystemInfo.GeneralLineNum=Line.size()-m_SystemInfo.LandBranchNum;
   /////读入变压器数据
   m_num=transformer.size()+RE;
   m_TransformerInfo=new TransformerInfo[m_num];
   for(i=1;i<=transformer.size();i++)
   {
       m_TransformerInfo[i].BusINo=transformer[i-1].BusINo;
	   m_TransformerInfo[i].BusJNo=transformer[i-1].BusJNo;
	   m_TransformerInfo[i].R=transformer[i-1].R;
	   m_TransformerInfo[i].X=transformer[i-1].X;
	   m_TransformerInfo[i].Ratio=transformer[i-1].Ratio;
   }
   m_SystemInfo.TransformerNum=transformer.size();
   ////读入节点信息
   int busnum=0;/// 节点数目
   for(i=0;i<generator.size();i++)
	   if(generator[i].BusNo>busnum) busnum=generator[i].BusNo;	   
   for(i=0;i<load.size();i++)
	  { if(load[i].BusNo>=busnum) busnum=load[i].BusNo;}                  
   m_SystemInfo.SystemTotalBusNum=busnum;
   m_BusInfo=new BusInfo[busnum+RE];
   for(i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
   {
      m_BusInfo[i].PD=0;
	  m_BusInfo[i].QD=0;
	  m_BusInfo[i].PG=0;
	  m_BusInfo[i].QG=0;
   }
    for(i=1;i<=load.size();i++)
	 {
	  m_BusInfo[load[i-1].BusNo].BusType=load[i-1].BusType;
	  m_BusInfo[load[i-1].BusNo].BusNo=load[i-1].BusNo;
	  m_BusInfo[load[i-1].BusNo].PD=load[i-1].PD;
	  m_BusInfo[load[i-1].BusNo].QD=load[i-1].QD;
	  m_BusInfo[load[i-1].BusNo].Voltage_e=load[i-1].Voltage_e;
       m_BusInfo[load[i-1].BusNo].Voltage_f=0;
	}
  for(i=1;i<=generator.size();i++)
	 {
	  m_BusInfo[generator[i-1].BusNo].BusType=generator[i-1].BusType;
	  m_BusInfo[generator[i-1].BusNo].BusNo=generator[i-1].BusNo;
	  m_BusInfo[generator[i-1].BusNo].PG=generator[i-1].PG;
	  m_BusInfo[generator[i-1].BusNo].QG=generator[i-1].QG;
	  m_BusInfo[generator[i-1].BusNo].Voltage_e=generator[i-1].Voltage_e;
      m_BusInfo[generator[i-1].BusNo].Voltage_f=0;
	  if(generator[i-1].BusType==0)  m_SystemInfo.SwingBusNo=generator[i-1].BusNo;
	}
//    for(i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
//		m_BusInfo[i].VolMod=m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+
//		                    m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f;
//	for(i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
//   cout<<m_BusInfo[i].BusNo<<"	"<<m_BusInfo[i].BusType<<"	"
//	 <<m_BusInfo[i].PG<<"	"<<m_BusInfo[i].QG<<"	"<<m_BusInfo[i].PD
//	 <<"	"<<m_BusInfo[i].QD<<"	"<<"--"<<m_BusInfo[i].Voltage_e<<endl;
}
