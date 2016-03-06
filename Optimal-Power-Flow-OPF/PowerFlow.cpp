/////////////////////////////////////////////////////////////////////////////////////////////
///////////   ��������  �������������������ʼ���ڵ����з�ʽ                    ///////////
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
   extern 	BusInfo *m_BusInfo;        //// �ڵ�
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
   extern 	BusInfo *m_BusInfo;        //// �ڵ�

  Initialize();
   int nodenum=m_SystemInfo.SystemTotalBusNum;
  ofstream dataout;  
	if(nodenum==5)
	{
		dataout.open("����ļ�//5�ڵ�//�������.txt");
	}
	else if(nodenum==14)
	{
	dataout.open("����ļ�//14�ڵ�//�������.txt");
	}
	else if(nodenum==30)
	{
	dataout.open("����ļ�//30�ڵ�//�������.txt");
	}
    else if(nodenum==57)
	{
	dataout.open("����ļ�//57�ڵ�//�������.txt");
	}
	else if(nodenum==118)
	{
	dataout.open("����ļ�//118�ڵ�//�������.txt");
	}
	dataout<<"          ��"<<Itearnum<<"�ε���:"<<endl;
    dataout<<"----------------------------------------"<<endl;
	dataout<<"�ڵ���"<<"      "<<"��ѹʵ��"<<"       "<<"��ѹ�鲿"<<endl;
  Itearnum=0;
  double iteraerror=100;
  while(Itearnum<=MaxIteraNum&&IteraError<iteraerror)
  {

   iteraerror=GetDeltaV();
   Itearnum++;
    dataout<<"          ��"<<Itearnum<<"�ε���:"<<endl;
    dataout<<"----------------------------------------"<<endl;
	dataout<<"�ڵ���"<<"      "<<"��ѹʵ��"<<"       "<<"��ѹ�鲿"<<endl;
   for(int k=1;k<=nodenum;k++)
   {
   dataout<<m_BusInfo[k].BusNo<<"             "<<m_BusInfo[k].Voltage_e<<"          "<<m_BusInfo[k].Voltage_f<<endl;
   }
  }
  dataout.close();
   cout<<"��������"<<Itearnum<<endl;
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
	   double sum=0;  ////���ۼӺ�
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
      extern LineInfo *m_LineInfo;       ///��·��Ϣ
	  extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	  extern BusInfo *m_BusInfo;        //// �ڵ� 
	  
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
      
	  double *I_e=new double[totallinenum+RE];   ///����ʵ�� 
      double *I_f=new double[totallinenum+RE];   ///�����鲿
	  double *Pl_ij=new double[totallinenum+RE];  ///֧·�й�
      double *Pl_ji=new double[totallinenum+RE];  ///֧·�й�
	  double *Ql_ij=new double[totallinenum+RE];  ///֧·�޹�
	  double *Ql_ji=new double[totallinenum+RE];  ///֧·�޹�

	  double sumdP=0;
	  double sumdQ=0;
  
	  ////////////////����һ��֧·�Ĺ���
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

	  ////�����ѹ��֧·�Ĺ���
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

	 if(!strcmp(filename,"��������//5.txt"))
	 {
		 dataout.open("����ļ�//5�ڵ�//������.xls");
	 }
	 else if(!strcmp(filename,"��������//14.txt"))
	 {
         dataout.open("����ļ�//14�ڵ�//������.xls");
	 }
	 else if(!strcmp(filename,"��������//30.txt"))
	 {
         dataout.open("����ļ�//30�ڵ�//������.xls");
	 }
	 else if(!strcmp(filename,"��������//57.txt"))
	 {
         dataout.open("����ļ�//57�ڵ�//������.xls");
	 }
	 else if(!strcmp(filename,"��������//118.txt"))
	 {
		 dataout.open("����ļ�//118�ڵ�//������.xls");
	 }
      dataout<<"				"<<"ţ��-����ѷ��������-ֱ������"<<endl;
	  dataout<<endl;
	  dataout<<"ϵͳ������Ϣ:"<<endl;
	  dataout<<"�ڵ���:"<<"	"<<nodenum<<endl;
	  dataout<<"�������"<<"	"<<generalnum<<endl;
      dataout<<"֧·��"<<"	"<<totallinenum+m_SystemInfo.LandBranchNum<<endl;
      dataout<<"��ѹ����"<<"	"<<transformernum<<endl;
	  dataout<<"ƽ��ڵ���"<<"	"<<m_SystemInfo.SwingBusNo<<endl;
      dataout<<endl<<endl;
	  dataout<<"����������:"<<endl;
	  dataout<<"��ѹ��ֵ��͵Ľڵ�Ϊ:"<<"	"<<VminNode<<"	"<<"���ѹ��ֵΪ:"<<"	"<<Vmin<<endl;
      dataout<<"1.��·֧·���ʷֲ�:"<<endl;
	  dataout<<"�ڵ�i"<<"	"<<"�ڵ�j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=generallinenum;k++)
	  {
		  dataout<<m_LineInfo[k].BusINo<<"("<<new_old[m_LineInfo[k].BusINo]<<")"<<"	"
			  <<m_LineInfo[k].BusJNo<<"("<<new_old[m_LineInfo[k].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k]
			     <<"	"<<setprecision(6)<<Ql_ij[k]<<"	"<<setprecision(6)<<Pl_ji[k]<<"	"
				 <<setprecision(6)<<Ql_ji[k]<<"	"<<setprecision(6)<<dP[k]<<"	"
				 <<setprecision(6)<<dQ[k]<<endl;
	  }
	  dataout<<endl;
	  dataout<<"2.�ӵ�֧·���ʷֲ�:"<<endl;
	    dataout<<"�ڵ�i"<<"	"<<"�ڵ�j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=landbranchnum;k++)
	  {
         dataout<<m_LineInfo[k+generallinenum].BusINo<<"("<<new_old[m_LineInfo[k+generallinenum].BusINo]<<")"
			 <<"	"<<m_LineInfo[k+generallinenum].BusJNo<<"("<<new_old[m_LineInfo[k+generallinenum].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k+generallinenum]
			     <<"	"<<setprecision(6)<<Ql_ij[k+generallinenum]<<"	"<<setprecision(6)<<Pl_ji[k+generallinenum]<<"	"
				 <<setprecision(6)<<Ql_ji[k+generallinenum]<<"	"<<setprecision(6)<<dP[k+generallinenum]<<"	"<<setprecision(6)<<dQ[k+generallinenum]<<endl;
	  }   
	   dataout<<endl;
	  dataout<<"3.��ѹ��֧·���ʷֲ�:"<<endl;
	    dataout<<"�ڵ�i"<<"	"<<"�ڵ�j"<<"	"<<"Pij"<<"	"<<"Qij"<<"	"<<"Pji"<<"	"<<"Qji"<<"	"<<"dPij"<<"	"<<"dQij"<<endl;
	  for(k=1;k<=transformernum;k++)
	  {
      dataout<<m_TransformerInfo[k].BusINo<<"("<<new_old[m_TransformerInfo[k].BusINo]<<")"
		  <<"	"<<m_TransformerInfo[k].BusJNo<<"("<<new_old[m_TransformerInfo[k].BusJNo]<<")"<<"	"<<setprecision(6)<<Pl_ij[k+linenum]
			     <<"	"<<setprecision(6)<<Ql_ij[k+linenum]<<"	"<<setprecision(6)<<Pl_ji[k+linenum]<<"	"
				 <<setprecision(6)<<Ql_ji[k+linenum]<<"	"<<setprecision(6)<<dP[k+linenum]<<"	"<<setprecision(6)<<dQ[k+linenum]<<endl;
	  }
	   dataout<<endl;
	  dataout<<"4.�ڵ��ѹ����ع���:"<<endl;
	  dataout<<"�ڵ���"<<"	"<<"��ѹ��ֵ"<<"	"<<"��ѹ���"<<"	"<<"�����й�"<<"	"<<"�����޹�"
		    <<"	"<<"������й�"<<"	"<<"������޹�"<<endl;
	  for(k=1;k<=nodenum;k++)
	  {
		  dataout<<m_BusInfo[k].BusNo<<"("<<new_old[m_BusInfo[k].BusNo]<<")"<<"	"<<Vi[k]<<"	"<<theta[k]<<"	"<<m_BusInfo[k].PD<<"	"
			     <<m_BusInfo[k].QD<<"	"<<m_BusInfo[k].PG<<"	"<<m_BusInfo[k].QG<<endl;
	  }
	  dataout<<"���������еĽڵ�����ԭʼ�Ľڵ��,û�н����Ż����ʱΪ0"<<endl;
     dataout<<"ϵͳ���й�����:"<<sumdP<<endl;
     dataout<<"ϵͳ���޹�����:"<<sumdQ<<endl;
	 
	 dataout.close();
	  delete I_e;   ///����ʵ�� 
      delete I_f;   ///�����鲿
	  delete Pl_ij;  ///֧·�й�
      delete Pl_ji;  ///֧·�й�
	  delete Ql_ij;  ///֧·�޹�
	  delete Ql_ji;  ///֧·�޹�
	  delete dP;
	  delete dQ;
	  
}