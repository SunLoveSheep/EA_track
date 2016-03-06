//////////////////////////////////////////////////////////////////////////
/////  �ڵ��Ż����  ��̬��ţ���̬��ţ��붯̬���                   ////
//////////////////////////////////////////////////////////////////////////

#include "NodeOptimize.h"
#include <iostream.h>
#define RE 10

NodeOptimize::NodeOptimize()
{
	 extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
	int nodenum=m_SystemInfo.SystemTotalBusNum;
	old_new=new int[nodenum+RE];
	for(int k=1;k<=nodenum;k++)
		old_new[k]=0;
    new_old=new int[nodenum+RE];
	for(k=1;k<=nodenum;k++)
		new_old[k]=0;
}

NodeOptimize::~NodeOptimize()
{
    
}

void NodeOptimize::StaticOptimize()
{
    extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ�    
	int nodenum=m_SystemInfo.SystemTotalBusNum;
    int linenum=m_SystemInfo.GeneralLineNum;//+m_SystemInfo.LandBranchNum;
	int totallinenum=linenum+m_SystemInfo.TransformerNum;
	int *busIno=new int[totallinenum+RE];  ///��¼��·��Ž�С�˵�
	int *busJno=new int[totallinenum+RE];  ///��¼��·��Žϴ�˵�
    for(int i=1;i<=linenum;i++)
	{
		if(m_LineInfo[i].BusINo<=m_LineInfo[i].BusJNo)
		{
		busIno[i]=m_LineInfo[i].BusINo;
		busJno[i]=m_LineInfo[i].BusJNo;
		}
		else
        {
		busJno[i]=m_LineInfo[i].BusINo;
		busIno[i]=m_LineInfo[i].BusJNo;
		}
	}
	for(i=linenum+1;i<=totallinenum;i++)
	{
		if(m_TransformerInfo[i-linenum].BusINo<=m_TransformerInfo[i-linenum].BusJNo)
		{		
		busIno[i]=m_TransformerInfo[i-linenum].BusINo;	
		busJno[i]=m_TransformerInfo[i-linenum].BusJNo;
		}
		else
		{
		busIno[i]=m_TransformerInfo[i-linenum].BusJNo;
        busJno[i]=m_TransformerInfo[i-linenum].BusINo;
		}
	}	
	int *nodedegree=new int[nodenum+RE];   ///�ڵ��
	for(i=1;i<=nodenum;i++)
		nodedegree[i]=0;
    int newbusno=0;  ///�µĽڵ���
	int oldbusno=0;  ///�ϱ��
	for(i=1;i<=nodenum;i++)
	{
	for(int j=1;j<=totallinenum;j++)
	{
		if(busIno[j]==i) nodedegree[i]++;
		if(busJno[j]==i) nodedegree[i]++;
		/////����˫���ߵ�Ӱ��
		if(busIno[j]==busIno[j-1]&&busJno[j]==busJno[j-1]) nodedegree[i]--;
	} 

	} 
	int maxdegree=0,mindegree=10000; ////��¼���ڵ��������С�ڵ����
	for(i=1;i<=nodenum;i++)
	{
		if(nodedegree[i]>maxdegree) 
		{
			maxdegree=nodedegree[i];
		}
		if(nodedegree[i]<mindegree)
		{
			mindegree=nodedegree[i];
		}
	}
	
    for(i=mindegree;i<=maxdegree;i++)
	{
		for(int j=1;j<=nodenum;j++)
		{
			if(nodedegree[j]==i) 
			{
				newbusno++;
				old_new[j]=newbusno;				
			}
		}
	}


   for(i=1;i<=nodenum;i++)
	   cout<<i<<old_new[i]<<endl;

   	delete [] nodedegree;
	delete [] busIno;
	delete [] busJno;
}

void NodeOptimize::HalfDymOptimize()
{
    extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ�    
//    
	int nodenum=m_SystemInfo.SystemTotalBusNum;
//	old_new=new int[nodenum+RE];
//	for(int k=1;k<=nodenum;k++)
//		old_new[k]=0;
//    new_old=new int[nodenum+RE];
//	for(k=1;k<=nodenum;k++)
//		new_old[k]=0;

	int linenum=m_SystemInfo.GeneralLineNum;//+m_SystemInfo.LandBranchNum;
	int totallinenum=linenum+m_SystemInfo.TransformerNum;
	int *busIno=new int[totallinenum+RE];  ///��¼��·��Ž�С�˵�
	int *busJno=new int[totallinenum+RE];  ///��¼��·��Žϴ�˵�
    for(int i=1;i<=linenum;i++)
	{
		if(m_LineInfo[i].BusINo<=m_LineInfo[i].BusJNo)
		{
		busIno[i]=m_LineInfo[i].BusINo;
		busJno[i]=m_LineInfo[i].BusJNo;
		}
		else
        {
		busJno[i]=m_LineInfo[i].BusINo;
		busIno[i]=m_LineInfo[i].BusJNo;
		}
	}
	for(i=linenum+1;i<=totallinenum;i++)
	{
		if(m_TransformerInfo[i-linenum].BusINo<=m_TransformerInfo[i-linenum].BusJNo)
		{		
		busIno[i]=m_TransformerInfo[i-linenum].BusINo;	
		busJno[i]=m_TransformerInfo[i-linenum].BusJNo;
		}
		else
		{
		busIno[i]=m_TransformerInfo[i-linenum].BusJNo;
        busJno[i]=m_TransformerInfo[i-linenum].BusINo;
		}
	}	
	int *nodedegree=new int[nodenum+RE];   ///�ڵ��
	for(i=1;i<=nodenum;i++)
		nodedegree[i]=0;
    int newbusno=0;  ///�µĽڵ���
	int oldbusno=0;  ///�ϱ��
	while(newbusno<nodenum)
	{
		for(i=1;i<=nodenum;i++)
		{
		for(int j=1;j<=totallinenum;j++)
		{
			if(busIno[j]==i) nodedegree[i]++;
			if(busJno[j]==i) nodedegree[i]++;
			/////����˫���ߵ�Ӱ��
			if(busIno[j]==busIno[j-1]&&busJno[j]==busJno[j-1]) nodedegree[i]--;
		}
// 		cout<<nodedegree[i]<<endl;
		}
//     cout<<"------------------"<<endl;
		int max=nodedegree[1];
        for(i=1;i<=nodenum;i++)    ///�ҵ�����ѭ����,�ڵ������С�Ľڵ�
		{
			if(nodedegree[i]<=max)
			{ 
				max=nodedegree[i];  
	            oldbusno=i;   		    	
			} 
		}
		newbusno++;
//		cout<<oldbusno<<endl;
//		cout<<"-------------------"<<endl;
		old_new[oldbusno]=newbusno;
		new_old[newbusno]=oldbusno;
		/////////////////////////
//		for(int j=1;j<=totallinenum;j++)
//			cout<<busIno[j]<<"  "<<busJno[j]<<endl;
//		cout<<"-------------------"<<endl;
        //////////////////////////////
		for(i=1;i<=nodenum;i++)
	    {
			nodedegree[i]=0;	
		}
		for(i=1;i<=nodenum;i++)
		{
			nodedegree[new_old[i]]=100000;
		}
		for(int j=1;j<=totallinenum;j++)
		{
			if(busIno[j]==oldbusno) 
			{
				busIno[j]=10000+oldbusno; busJno[j]=20000+oldbusno;
			}
			if(busJno[j]==oldbusno) 
			{
				busIno[j]=10000+oldbusno; busJno[j]=20000+oldbusno;
			}
		}

	}
//   for(i=1;i<=nodenum;i++)
//	   cout<<i<<old_new[i]<<endl;
    ///��û�е�������Ľڵ���
	int no_out=0;
   for(i=1;i<=nodenum;i++)
   	   if(old_new[i]<0)  no_out++;
	   newbusno=newbusno-no_out;
	   for(i=1;i<=nodenum;i++)
		   if(old_new[i]<0)
		   {  
			   old_new[i]=newbusno; 
			   newbusno++;
		   }

   for(i=1;i<=nodenum;i++)
	   cout<<i<<old_new[i]<<endl;

	delete [] nodedegree;
	delete [] busIno;
	delete [] busJno;
    
}

void NodeOptimize::DynamicOptimize()
{
     extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ�    
//    
	int nodenum=m_SystemInfo.SystemTotalBusNum;
//	old_new=new int[nodenum+RE];
//	for(int k=1;k<=nodenum;k++)
//		old_new[k]=0;
//    new_old=new int[nodenum+RE];
//	for(k=1;k<=nodenum;k++)
//		new_old[k]=0;

	int linenum=m_SystemInfo.GeneralLineNum;//+m_SystemInfo.LandBranchNum;
	int totallinenum=linenum+m_SystemInfo.TransformerNum;
	int *busIno=new int[totallinenum+RE];  ///��¼��·��Ž�С�˵�
	int *busJno=new int[totallinenum+RE];  ///��¼��·��Žϴ�˵�
    for(int i=1;i<=linenum;i++)
	{
		if(m_LineInfo[i].BusINo<=m_LineInfo[i].BusJNo)
		{
		busIno[i]=m_LineInfo[i].BusINo;
		busJno[i]=m_LineInfo[i].BusJNo;
		}
		else
        {
		busJno[i]=m_LineInfo[i].BusINo;
		busIno[i]=m_LineInfo[i].BusJNo;
		}
	}
	for(i=linenum+1;i<=totallinenum;i++)
	{
		if(m_TransformerInfo[i-linenum].BusINo<=m_TransformerInfo[i-linenum].BusJNo)
		{		
		busIno[i]=m_TransformerInfo[i-linenum].BusINo;	
		busJno[i]=m_TransformerInfo[i-linenum].BusJNo;
		}
		else
		{
		busIno[i]=m_TransformerInfo[i-linenum].BusJNo;
        busJno[i]=m_TransformerInfo[i-linenum].BusINo;
		}
	}	
	int *nodedegree=new int[nodenum+RE];   ///�ڵ��
	for(i=1;i<=nodenum;i++)
		nodedegree[i]=0;
    int newbusno=0;  ///�µĽڵ���
	int oldbusno=0;  ///�ϱ��
	int *tempnode=new int[nodenum+RE];///������¼��ĳ��ĳ���ڵ����������нڵ�
	while(newbusno<nodenum)
	{
		for(i=1;i<=nodenum;i++)
		{
		for(int j=1;j<=nodenum;j++)
		{
           tempnode[j]=0;
		}
		int k=0;
		for( j=1;j<=totallinenum;j++)
		{
			if(busIno[j]==i) 
			{
				nodedegree[i]++;
				k++;
				tempnode[k]=busJno[j];
			}
			if(busJno[j]==i) 
			{
				nodedegree[i]++;
				k++;
               tempnode[k]=busIno[j];
			}
			/////����˫���ߵ�Ӱ��
        	if(busIno[j]==busIno[j-1]&&busJno[j]==busJno[j-1]) nodedegree[i]--;			
		}
		int Dk=0;/// ��������ڵ�֮�����е�֧·��
		for(j=1;j<=k;j++)
		{
			for(int n=j+1;n<=k;n++)
			{
				for(int m=1;m<=totallinenum;m++)
				{
					if(busIno[m]==tempnode[j]&&busJno[m]==tempnode[n])
					{
                      Dk++;
					}
				}
			}
		}
		nodedegree[i]=nodedegree[i]*(nodedegree[i]-1)/2-Dk;
// 		cout<<nodedegree[i]<<endl;
		}
//     cout<<"------------------"<<endl;
		int max=nodedegree[1];
        for(i=1;i<=nodenum;i++)    ///�ҵ�����ѭ����,�ڵ������С�Ľڵ�
		{
			if(nodedegree[i]<=max)
			{ 
				max=nodedegree[i];  
	            oldbusno=i;   		    	
			} 
		}
		newbusno++;
//		cout<<oldbusno<<endl;
//		cout<<"-------------------"<<endl;
		old_new[oldbusno]=newbusno;
		new_old[newbusno]=oldbusno;
		/////////////////////////
//		for(int j=1;j<=totallinenum;j++)
//			cout<<busIno[j]<<"  "<<busJno[j]<<endl;
//		cout<<"-------------------"<<endl;
        //////////////////////////////
		for(i=1;i<=nodenum;i++)
	    {
			nodedegree[i]=0;	
		}
		for(i=1;i<=nodenum;i++)
		{
			nodedegree[new_old[i]]=100000;
		}
		for(int j=1;j<=totallinenum;j++)
		{
			if(busIno[j]==oldbusno) 
			{
				busIno[j]=10000+oldbusno; busJno[j]=20000+oldbusno;
			}
			if(busJno[j]==oldbusno) 
			{
				busIno[j]=10000+oldbusno; busJno[j]=20000+oldbusno;
			}
		}

	}

cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	for(i=1;i<=nodenum;i++)
	  cout<<i<<old_new[i]<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	delete [] nodedegree;
	delete [] busIno;
	delete [] busJno;

}

void NodeOptimize::AdjustNodeNo()
{
    if(OptimizeType==1)  StaticOptimize();
	else if (OptimizeType==2) HalfDymOptimize();
	else if(OptimizeType==3) DynamicOptimize();

    extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ�   
	
	int linenum,transformernum,nodenum;///֧·������ѹ�������ڵ���
	linenum=m_SystemInfo.GeneralLineNum+m_SystemInfo.LandBranchNum;
    transformernum=m_SystemInfo.TransformerNum;
	nodenum=m_SystemInfo.SystemTotalBusNum;
	for(int k=1;k<=linenum;k++)
	{
		int i=m_LineInfo[k].BusINo;
		int j=m_LineInfo[k].BusJNo;
        m_LineInfo[k].BusINo=old_new[i];
        m_LineInfo[k].BusJNo=old_new[j]; 
//	  cout<<m_LineInfo[k].BusINo<<"  "<<m_LineInfo[k].BusJNo<<endl;
	}
    for(k=1;k<=transformernum;k++)
	{
		int i=m_TransformerInfo[k].BusINo;
		int j=m_TransformerInfo[k].BusJNo;
		m_TransformerInfo[k].BusINo=old_new[i];
	    m_TransformerInfo[k].BusJNo=old_new[j];
//		cout<<m_TransformerInfo[k].BusINo<<"  "<<m_TransformerInfo[k].BusJNo<<endl;
	}
	///////////�Խڵ����±�ţ����Ұ��ձ�Ŵ�С˳����������˳������
	BusInfo *businfo;    
	businfo=new BusInfo[nodenum+RE];   
	for(k=1;k<=nodenum;k++)
	{
		int i=m_BusInfo[k].BusNo;
		m_BusInfo[k].BusNo=old_new[i];	
        businfo[old_new[i]].BusNo=m_BusInfo[k].BusNo;
        businfo[old_new[i]].BusType=m_BusInfo[k].BusType;
		businfo[old_new[i]].PG=m_BusInfo[k].PG;
		businfo[old_new[i]].QG=m_BusInfo[k].QG;
		businfo[old_new[i]].PD=m_BusInfo[k].PD;
		businfo[old_new[i]].QD=m_BusInfo[k].QD;
		businfo[old_new[i]].Voltage_e=m_BusInfo[k].Voltage_e;  
		businfo[old_new[i]].Voltage_f=m_BusInfo[k].Voltage_f;  
//		cout<<k<<" "<<m_BusInfo[k].BusNo<<endl;		
	}
   	for(k=1;k<=nodenum;k++)
	{
     m_BusInfo[k].BusNo=businfo[k].BusNo;
	 m_BusInfo[k].BusType=businfo[k].BusType;
	 m_BusInfo[k].PG=businfo[k].PG;
	 m_BusInfo[k].QG=businfo[k].QG;
	 m_BusInfo[k].PD=businfo[k].PD;
	 m_BusInfo[k].QD=businfo[k].QD;
	 m_BusInfo[k].Voltage_e=businfo[k].Voltage_e;
	 m_BusInfo[k].Voltage_f=businfo[k].Voltage_f;
	 if(m_BusInfo[k].BusType==0) m_SystemInfo.SwingBusNo=m_BusInfo[k].BusNo; 
	}

//	   for(int i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
//		m_BusInfo[i].VolMod=m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+
//		                    m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f;

//for(int i=1;i<=m_SystemInfo.SystemTotalBusNum;i++)
//   cout<<m_BusInfo[i].BusNo<<"	"<<m_BusInfo[i].BusType<<"	"
//	 <<m_BusInfo[i].PG<<"	"<<m_BusInfo[i].QG<<"	"<<m_BusInfo[i].PD
// 	 <<"	"<<m_BusInfo[i].QD<<"	"<<"--"<<m_BusInfo[i].Voltage_e<<endl;
}