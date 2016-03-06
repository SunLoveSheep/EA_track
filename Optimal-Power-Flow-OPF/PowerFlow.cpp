/////////////////////////////////////////////////////////////////////////////////////////////
///////////   潮流计算  迭代，输出计算结果，初始化节点运行方式                    ///////////
//////////                                                                        //////////
//////////////////////////////////////////////////////////////////////////

#include "PowerFlow.h"
#include "Jacobi.h"
#include <math.h>
#include <iostream.h>
#include <fstream>
#include <string.h>
#include<iomanip.h>
#define  RE 100
PowerFlow::PowerFlow()
{

}

PowerFlow::~PowerFlow()
{
  delete [] delta_Votltage;
}
void PowerFlow::Initialize()
{
   extern SystemInfo m_SystemInfo;
   extern 	BusInfo *m_BusInfo;        //// 节点
for(int i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
{
	m_BusInfo[i].VolMod=m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+
		                 m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f;
	if(m_BusInfo[i].BusType==1)
	{
	m_BusInfo[i].Voltage_e=1;
	m_BusInfo[i].Voltage_f=0;
	}
}
}

void PowerFlow::PFlow()
{  
   extern SystemInfo m_SystemInfo;
   extern 	BusInfo *m_BusInfo;        //// 节点

  Initialize();
   int nodenum=m_SystemInfo.SystemTotalBusNum;
  ofstream dataout;  
	if(nodenum==5)
	{
		dataout.open("输出文件//5节点//迭代结果.txt");
	}
	else if(nodenum==14)
	{
	dataout.open("输出文件//14节点//迭代结果.txt");
	}
	else if(nodenum==30)
	{
	dataout.open("输出文件//30节点//迭代结果.txt");
	}
    else if(nodenum==57)
	{
	dataout.open("输出文件//57节点//迭代结果.txt");
	}
	else if(nodenum==118)
	{
	dataout.open("输出文件//118节点//迭代结果.txt");
	}
	dataout<<"          第"<<Itearnum<<"次迭代:"<<endl;
    dataout<<"----------------------------------------"<<endl;
	dataout<<"节点编号"<<"      "<<"电压实部"<<"       "<<"电压虚部"<<endl;
  Itearnum=0;
  double iteraerror=100;
  while(Itearnum<=MaxIteraNum&&IteraError<iteraerror)
  {

   iteraerror=GetDeltaV();
   Itearnum++;
    dataout<<"          第"<<Itearnum<<"次迭代:"<<endl;
    dataout<<"----------------------------------------"<<endl;
	dataout<<"节点编号"<<"      "<<"电压实部"<<"       "<<"电压虚部"<<endl;
   for(int k=1;k<=nodenum;k++)
   {
   dataout<<m_BusInfo[k].BusNo<<"             "<<m_BusInfo[k].Voltage_e<<"          "<<m_BusInfo[k].Voltage_f<<endl;
   }
  }
  dataout.close();
   cout<<"迭代次数"<<Itearnum<<endl;
}
double PowerFlow::GetDeltaV()
{
   extern SystemInfo m_SystemInfo;
   extern 	BusInfo  *m_BusInfo;        //// bus

   Jacobi jacobi;
   jacobi.GetNewInNum();
   jacobi.FormJacobi();
   int nodenum=m_SystemInfo.SystemTotalBusNum;
   delta_Votltage=new double[2*(nodenum-1)+RE];

 //  delta_Votltage[2*(nodenum-1)]=jacobi.delta_PQV[2*(nodenum-1)];

   for(int k=2*(nodenum-1);k>=1;k--)
   {
	   double sum=0;  ////求累加和
   for(int i=jacobi.FTIU[k];i<jacobi.FTIU[k+1];i++)
   {
	   int j=jacobi.FTJU[i];
       sum+=jacobi.FTU[i]*delta_Votltage[j];   
   }
    delta_Votltage[k]=jacobi.delta_PQV[k]-sum;
//     cout<<delta_Votltage[k]<<endl;
   }

   for(k=1;k<=nodenum;k++)
   {
	   if(m_BusInfo[k].BusNo<m_SystemInfo.SwingBusNo)
	   {
	    m_BusInfo[k].Voltage_e=m_BusInfo[k].Voltage_e-delta_Votltage[2*k-1];
        m_BusInfo[k].Voltage_f=m_BusInfo[k].Voltage_f-delta_Votltage[2*k];
	   }
	   if(m_BusInfo[k].BusNo>m_SystemInfo.SwingBusNo)
	   {
        m_BusInfo[k].Voltage_e=m_BusInfo[k].Voltage_e-delta_Votltage[2*(k-1)-1];
        m_BusInfo[k].Voltage_f=m_BusInfo[k].Voltage_f-delta_Votltage[2*(k-1)];
	   }
//	   dataout<<m_BusInfo[k].BusNo<<"      "<<m_BusInfo[k].Voltage_e<<"      "<<m_BusInfo[k].Voltage_f<<endl;
   }

   double maxerror=0;
    for(k=2*(nodenum-1);k>=1;k--)
	{
		if(fabs(jacobi.delta_PQV[k])>=maxerror)
		{
			maxerror=fabs(jacobi.delta_PQV[k]);
		}
	}	
  return maxerror;
}

void PowerFlow::OutputResult(char *filename,int *new_old)
{
      extern SystemInfo m_SystemInfo;
      extern LineInfo *m_LineInfo;       ///线路信息
	  extern TransformerInfo *m_TransformerInfo;  ///变压器支路 
	  extern BusInfo *m_BusInfo;        //// 节点 
	  
	  int generallinenum=m_SystemInfo.GeneralLineNum;
	  int transformernum=m_SystemInfo.TransformerNum;
	  int landbranchnum=m_SystemInfo.LandBranchNum;
	  int linenum=generallinenum+landbranchnum;
	  int totallinenum=m_SystemInfo.GeneralLineNum+m_SystemInfo.TransformerNum+landbranchnum;
	  int nodenum=m_SystemInfo.SystemTotalBusNum;

	  double *theta=new double[nodenum+RE];
	  double *Vi=new double[nodenum+RE];
	  double *Pl=new double[nodenum+RE];
	  double *Ql=new double[nodenum+RE];
      double *dP=new double[nodenum+RE];
      double *dQ=new double[nodenum+RE];
       
	  double Vmin=1000.0;
	  int VminNode=0;
	  for(int m=1;m<=nodenum;m++)
	  {
         theta[m]=atan((m_BusInfo[m].Voltage_f/m_BusInfo[m].Voltage_e))/3.1415926*180;
		 Vi[m]=sqrt(m_BusInfo[m].Voltage_e*m_BusInfo[m].Voltage_e+m_BusInfo[m].Voltage_f*m_BusInfo[m].Voltage_f);
	     if(Vi[m]<Vmin)
		 {
           Vmin=Vi[m];
		   VminNode=m;
		 }		 
	  }
      
	  double *I_e=new double[totallinenum+RE];   ///电流实部 
      double *I_f=new double[totallinenum+RE];   ///电流虚部
	  double *Pl_ij=new double[totallinenum+RE];  ///支路有功
      double *Pl_ji=new double[totallinenum+RE];  ///支路有功
	  double *Ql_ij=new double[totallinenum+RE];  ///支路无功
	  double *Ql_ji=new double[totallinenum+RE];  ///支路无功

	  double sumdP=0;
	  double sumdQ=0;
  
	  ////////////////计算一般支路的功率
	  for(int k=1;k<=linenum;k++)
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
		  Pl_ij[k]=I_e[k]*m_BusInfo[i].Voltage_e+I_f[k]*m_BusInfo[i].Voltage_f;	  
		  Ql_ij[k]=I_e[k]*m_BusInfo[i].Voltage_f-I_f[k]*m_BusInfo[i].Voltage_e;
          Pl_ji[k]=-I_e[k]*m_BusInfo[j].Voltage_e-I_f[k]*m_BusInfo[j].Voltage_f;	  
		  Ql_ji[k]=-I_e[k]*m_BusInfo[j].Voltage_f+I_f[k]*m_BusInfo[j].Voltage_e;

		  double Vi=m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f;
		  double Vj=m_BusInfo[j].Voltage_e*m_BusInfo[j].Voltage_e+m_BusInfo[j].Voltage_f*m_BusInfo[j].Voltage_f;
		  Ql_ij[k]=Ql_ij[k]-Vi*m_LineInfo[k].Yk;
          Ql_ji[k]=Ql_ji[k]-Vj*m_LineInfo[k].Yk;
		  
		  dP[k]=Pl_ij[k]+Pl_ji[k];
		  dQ[k]=Ql_ij[k]+Ql_ji[k];
		  sumdP+=dP[k];
          sumdQ+=dQ[k];
		  if(m_BusInfo[i].BusType==-1)
		  {
		    m_BusInfo[i].QG+=Ql_ij[k];
		  }
		  if(m_BusInfo[j].BusType==-1)
		  {
            m_BusInfo[j].QG+=Ql_ji[k];
		  }
		  if(m_BusInfo[i].BusType==0)
		  {
          m_BusInfo[i].PG+=Pl_ij[k];
            m_BusInfo[i].QG+=Ql_ij[k];
		  }
		  if(m_BusInfo[j].BusType==0)
		  {
            m_BusInfo[j].PG+=Pl_ji[k];
            m_BusInfo[j].QG+=Ql_ji[k];
 		  }
	  }

	  ////计算变压器支路的功率
	  for(k=linenum+1;k<=totallinenum;k++)
	  {
          
		  int n=k-linenum;
		  int i=m_TransformerInfo[n].BusINo;
          int j=m_TransformerInfo[n].BusJNo;
      
          double Tvoltage_e=m_BusInfo[j].Voltage_e/m_TransformerInfo[n].Ratio;
		  double Tvoltage_f=m_BusInfo[j].Voltage_f/m_TransformerInfo[n].Ratio; 

		  double z=m_TransformerInfo[n].R*m_TransformerInfo[n].R+m_TransformerInfo[n].X*m_TransformerInfo[n].X;
		  double Gij=m_TransformerInfo[n].R/z;
		  double Bij=-m_TransformerInfo[n].X/z;
		 
		  I_e[k]=(m_BusInfo[i].Voltage_e-Tvoltage_e)*Gij
			     -(m_BusInfo[i].Voltage_f-Tvoltage_f)*Bij;
          I_f[k]=(m_BusInfo[i].Voltage_e-Tvoltage_e)*Bij
			     +(m_BusInfo[i].Voltage_f-Tvoltage_f)*Gij;
		  
		  Pl_ij[k]=I_e[k]*m_BusInfo[i].Voltage_e+I_f[k]*m_BusInfo[i].Voltage_f;	  
		  Ql_ij[k]=I_e[k]*m_BusInfo[i].Voltage_f-I_f[k]*m_BusInfo[i].Voltage_e;
          Pl_ji[k]=-I_e[k]*Tvoltage_e-I_f[k]*Tvoltage_f;	  
		  Ql_ji[k]=-I_e[k]*Tvoltage_f+I_f[k]*Tvoltage_e;
		  
		  dP[k]=Pl_ij[k]+Pl_ji[k];
		  dQ[k]=Ql_ij[k]+Ql_ji[k];
		  sumdP+=dP[k];
          sumdQ+=dQ[k];
		  
		 if(m_BusInfo[i].BusType==-1)
		  {
		    m_BusInfo[i].QG+=Ql_ij[k];
		  }
		  if(m_BusInfo[j].BusType==-1)
		  {
            m_BusInfo[j].QG+=Ql_ji[k];
		  }
		  if(m_BusInfo[i].BusType==0)
		  {
          m_BusInfo[i].PG+=Pl_ij[k];
            m_BusInfo[i].QG+=Ql_ij[k];
		  }
		  if(m_BusInfo[j].BusType==0)
		  {
            m_BusInfo[j].PG+=Pl_ji[k];
            m_BusInfo[j].QG+=Ql_ji[k];
 		  }
	  }  
	  int generalnum=0;
	  for(k=1;k<=nodenum;k++)
	  {
		  if(m_BusInfo[k].BusType==-1)
		  {
			  m_BusInfo[k].QG+=m_BusInfo[k].QD;
			  generalnum++;
		  }
		  if(m_BusInfo[k].BusType==0)
		  {
			  m_BusInfo[k].PG+=m_BusInfo[k].PD;
              m_BusInfo[k].QG+=m_BusInfo[k].QD;
		  }
	  }


     ofstream dataout;  	
	 dataout<<setiosflags(ios::fixed); 

	 if(!strcmp(filename,"输入数据//5.txt"))
	 {
		 dataout.open("输出文件//5节点//计算结果.xls");
	 }
	 else if(!strcmp(filename,"输入数据//14.txt"))
	 {
         dataout.open("输出文件//14节点//计算结果.xls");
	 }
	 else if(!strcmp(filename,"输入数据//30.txt"))
	 {
         dataout.open("输出文件//30节点//计算结果.xls");
	 }
	 else if(!strcmp(filename,"输入数据//57.txt"))
	 {
         dataout.open("输出文件//57节点//计算结果.xls");
	 }
	 else if(!strcmp(filename,"输入数据//118.txt"))
	 {
		 dataout.open("输出文件//118节点//计算结果.xls");
	 }
      dataout<<"				"<<"牛顿-拉夫逊潮流计算-直角坐标"<<endl;
	  dataout<<endl;
	  dataout<<"系统基本信息:"<<endl;
	  dataout<<"节点数:"<<"	"<<nodenum<<endl;
	  dataout<<"发电机数"<<"	"<<generalnum<<endl;
      dataout<<"支路数"<<"	"<<totallinenum+m_SystemInfo.LandBranchNum<<endl;
      dataout<<"变压器输"<<"	"<<transformernum<<endl;
	  dataout<<"平衡节点编号"<<"	"<<m_SystemInfo.SwingBusNo<<endl;
      dataout<<endl<<endl;
	  dataout<<"潮流计算结果:"<<endl;
	  dataout<<"电压幅值最低的节点为:"<<"	"<<VminNode<<"	"<<"其电压幅值为:"<<"	"<<Vmin<<endl;
      dataout<<"1.线路支路功率分布:"<<endl;
	  dataout<<"节点i"<<"	"<<"节点j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=generallinenum;k++)
	  {
		  dataout<<m_LineInfo[k].BusINo<<"("<<new_old[m_LineInfo[k].BusINo]<<")"<<"	"
			  <<m_LineInfo[k].BusJNo<<"("<<new_old[m_LineInfo[k].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k]
			     <<"	"<<setprecision(6)<<Ql_ij[k]<<"	"<<setprecision(6)<<Pl_ji[k]<<"	"
				 <<setprecision(6)<<Ql_ji[k]<<"	"<<setprecision(6)<<dP[k]<<"	"
				 <<setprecision(6)<<dQ[k]<<endl;
	  }
	  dataout<<endl;
	  dataout<<"2.接地支路功率分布:"<<endl;
	    dataout<<"节点i"<<"	"<<"节点j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=landbranchnum;k++)
	  {
         dataout<<m_LineInfo[k+generallinenum].BusINo<<"("<<new_old[m_LineInfo[k+generallinenum].BusINo]<<")"
			 <<"	"<<m_LineInfo[k+generallinenum].BusJNo<<"("<<new_old[m_LineInfo[k+generallinenum].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k+generallinenum]
			     <<"	"<<setprecision(6)<<Ql_ij[k+generallinenum]<<"	"<<setprecision(6)<<Pl_ji[k+generallinenum]<<"	"
				 <<setprecision(6)<<Ql_ji[k+generallinenum]<<"	"<<setprecision(6)<<dP[k+generallinenum]<<"	"<<setprecision(6)<<dQ[k+generallinenum]<<endl;
	  }   
	   dataout<<endl;
	  dataout<<"3.变压器支路功率分布:"<<endl;
	    dataout<<"节点i"<<"	"<<"节点j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=transformernum;k++)
	  {
      dataout<<m_TransformerInfo[k].BusINo<<"("<<new_old[m_TransformerInfo[k].BusINo]<<")"
		  <<"	"<<m_TransformerInfo[k].BusJNo<<"("<<new_old[m_TransformerInfo[k].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k+linenum]
			     <<"	"<<setprecision(6)<<Ql_ij[k+linenum]<<"	"<<setprecision(6)<<Pl_ji[k+linenum]<<"	"
				 <<setprecision(6)<<Ql_ji[k+linenum]<<"	"<<setprecision(6)<<dP[k+linenum]<<"	"<<setprecision(6)<<dQ[k+linenum]<<endl;
	  }
	   dataout<<endl;
	  dataout<<"4.节点电压及相关功率:"<<endl;
	  dataout<<"节点编号"<<"	"<<"电压幅值"<<"	"<<"电压相角"<<"	"<<"负荷有功"<<"	"<<"负荷无功"
		    <<"	"<<"发电机有功"<<"	"<<"发电机无功"<<endl;
	  for(k=1;k<=nodenum;k++)
	  {
		  dataout<<m_BusInfo[k].BusNo<<"("<<new_old[m_BusInfo[k].BusNo]<<")"<<"	"<<Vi[k]<<"	"<<theta[k]<<"	"<<m_BusInfo[k].PD<<"	"
			     <<m_BusInfo[k].QD<<"	"<<m_BusInfo[k].PG<<"	"<<m_BusInfo[k].QG<<endl;
	  }
	  dataout<<"以上括号中的节点编号是原始的节点号,没有进行优化编号时为0"<<endl;
     dataout<<"系统总有功网损:"<<sumdP<<endl;
     dataout<<"系统总无功网损:"<<sumdQ<<endl;
	 
	 dataout.close();
	  delete I_e;   ///电流实部 
      delete I_f;   ///电流虚部
	  delete Pl_ij;  ///支路有功
      delete Pl_ji;  ///支路有功
	  delete Ql_ij;  ///支路无功
	  delete Ql_ji;  ///支路无功
	  delete dP;
	  delete dQ;
	  
}