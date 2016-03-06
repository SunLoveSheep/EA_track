///////////////////////////////////////////////////////////////////////////////////////////
/////////////          潮流计算的修正方程的处理                    ////////////////////////
////////////  形成雅克比矩阵，形成迭代因子表，确定消去过程中形成的注入元的位置      ///////
//                      
///////////////////////////////////////////////////////////////////////////////////////////

#include "Jacobi.h"
#include <iostream.h>
#define RE 10
Jacobi::Jacobi()
{
   
}
Jacobi::~Jacobi()
{
	
	delete [] FTU;
	delete [] FTIU;
	delete [] FTJU;
	delete [] delta_PQV;
    delete [] JI;
	delete [] JJ;
	delete [] JLI;
	delete [] JLJ;
//	cout<<"2222"<<endl;
}

void Jacobi::FormJacobi()
{

	Ymatrix.SpareYMatrix();

    extern SystemInfo m_SystemInfo;   ///系统信息
    extern LineInfo *m_LineInfo;       ///线路信息
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路 
	extern BusInfo *m_BusInfo;        //// 节点 
	
	int nodenum=m_SystemInfo.SystemTotalBusNum;
    delta_PQV=new double[2*(nodenum-1)+RE];


	for(int i=1;i<=nodenum;i++)     /////逐行形成雅克比矩阵
	{	
	double ai=0,bi=0;////计算主对角元时，ai=(Gij*ej-Bij*fj)累加和，bi=(Gij*fj+Bij*ej)累加和
    int swingbusno=m_SystemInfo.SwingBusNo;
    if(i==swingbusno) i=i+1;   ///跳过平衡节点
//	cout<<Ymatrix.IY[i]<<endl;
	    int num=Ymatrix.IY[i+1]-Ymatrix.IY[i]; ////雅克比矩阵每行非零元的个数
	    int JJno=0;//// 每一行元素在JH,JN,JM,JL的下标 
		JH=new double[num+RE];   ////？？？开辟的空间不准确
		JN=new double[num+RE];
		JM=new double[num+RE];
		JL=new double[num+RE];

	////计算雅克比矩阵下三角部分
	  for(int m=1;m<Ymatrix.IY[i];m++)
	  {
		  int j=Ymatrix.JY[m];
	      double Gij=Ymatrix.UY[m].G;
		  double Bij=Ymatrix.UY[m].B;
		
		  if(j==i) 
		  {	
			  ////确定在哪一行
	 		  int I;
			  for(int x=1;x<i;x++)
			  {
				  if(m>=Ymatrix.IY[x]&&m<=Ymatrix.IY[x+1]-1) I=x;
			  }
			  ai+=Gij*m_BusInfo[I].Voltage_e-Bij*m_BusInfo[I].Voltage_f;
			  bi+=Gij*m_BusInfo[I].Voltage_f+Bij*m_BusInfo[I].Voltage_e;	  
			  ///计算下三角的非对角元
		//	if(I!=swingbusno&&j!=swingbusno+1)
			  if(I!=swingbusno)
			  {
              JJno++;
			  JH[JJno]=-(Gij*m_BusInfo[i].Voltage_e+Bij*m_BusInfo[i].Voltage_f);
              JN[JJno]=Bij*m_BusInfo[i].Voltage_e-Gij*m_BusInfo[i].Voltage_f;
			  if(m_BusInfo[i].BusType==-1)///PV节点的电压i对j节点电压偏导为0
			  {
			  JM[JJno]=0;
              JL[JJno]=0;
			  }
              if(m_BusInfo[i].BusType==1)  //PQ节点
			  {
			  JM[JJno]=JN[JJno];
              JL[JJno]=-JH[JJno];
			  }
			  }
		  }	
	  }

      //计算雅克比矩阵上三角部分	  
      int JJDno=JJno+1;  ///对角元的下标
	  JJno=JJno+1;
	  for(int k=Ymatrix.IY[i];k<Ymatrix.IY[i+1];k++)   
	  {		
		  int j=Ymatrix.JY[k];
		  double Gij=Ymatrix.UY[k].G;
		  double Bij=Ymatrix.UY[k].B;
		  ai+=Gij*m_BusInfo[j].Voltage_e-Bij*m_BusInfo[j].Voltage_f;
		  bi+=Gij*m_BusInfo[j].Voltage_f+Bij*m_BusInfo[j].Voltage_e;
		  //////计算上三角的非对角元
		  if(j!=swingbusno)
		  {
          JJno++;
		  JH[JJno]=-(Gij*m_BusInfo[i].Voltage_e+Bij*m_BusInfo[i].Voltage_f);
          JN[JJno]=Bij*m_BusInfo[i].Voltage_e-Gij*m_BusInfo[i].Voltage_f;
	      if(m_BusInfo[i].BusType==-1)   ///PV节点的电压i对j节点电压偏导为0
			  {
			  JM[JJno]=0;
              JL[JJno]=0;
			  }
          if(m_BusInfo[i].BusType==1) //PQ节点
			  {
			  JM[JJno]=JN[JJno];
              JL[JJno]=-JH[JJno];
			  }
		  }
	  } 
      ///计算雅克比矩阵对角元素
	  double Gii=Ymatrix.DY[i].G;
	  double Bii=Ymatrix.DY[i].B;
	  ai+=Gii*m_BusInfo[i].Voltage_e-Bii*m_BusInfo[i].Voltage_f;
      bi+=Gii*m_BusInfo[i].Voltage_f+Bii*m_BusInfo[i].Voltage_e;

	  JH[JJDno]=-ai-(Gii*m_BusInfo[i].Voltage_e+Bii*m_BusInfo[i].Voltage_f);
	  JN[JJDno]=-bi+(Bii*m_BusInfo[i].Voltage_e-Gii*m_BusInfo[i].Voltage_f);
	  int sno;///功率修正量的下标
	  sno=i;
	  if(i>swingbusno) sno=i-1;
	 delta_PQV[2*sno]=m_BusInfo[i].PG-m_BusInfo[i].PD-m_BusInfo[i].Voltage_e*ai-m_BusInfo[i].Voltage_f*bi;
	  
	 if(m_BusInfo[i].BusType==1)  //PQ节点
	  {
	  JM[JJDno]=bi+(Bii*m_BusInfo[i].Voltage_e-Gii*m_BusInfo[i].Voltage_f);
	  JL[JJDno]=-ai+(Gii*m_BusInfo[i].Voltage_e+Bii*m_BusInfo[i].Voltage_f);
	 delta_PQV[2*sno-1]=m_BusInfo[i].QG-m_BusInfo[i].QD-m_BusInfo[i].Voltage_f*ai+m_BusInfo[i].Voltage_e*bi;
	  }
     if(m_BusInfo[i].BusType==-1)   //pv节点
	 {
	  JM[JJDno]=-2*m_BusInfo[i].Voltage_e;
	  JL[JJDno]=-2*m_BusInfo[i].Voltage_f;
	  delta_PQV[2*sno-1]=m_BusInfo[i].VolMod-(m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f);
	 }
// 	 cout<<sno<<"PQV"<<delta_PQV[2*sno-1]<<"  "<<delta_PQV[2*sno]<<endl;
	 if(i==swingbusno-1&&swingbusno==nodenum) i=i+1; ///处理平衡节点为最后一个节点的情况
//   for(k=1;k<=JJno;k++)
//          {
//       	  cout<<i<<endl;
//            cout<<JM[k]<<" "<<JL[k]<<endl;
//       	 cout<<JH[k]<<" "<<JN[k]<<endl;
//       	 cout<<"++++"<<endl;
//           }
  //            cout<<"-------------------"<<endl;
    FormFactorTable(i);
	
	 delete [] JH;
	 delete [] JM;
	 delete [] JN;
	 delete [] JL;
}
}

//void Jacobi::FormFactorTable(int i,double *f_JH,double *f_JM,double *f_JN,double *f_JL)  ///消去
void Jacobi::FormFactorTable(int i)
{ 

    extern SystemInfo m_SystemInfo;   ///系统信息
   
    int nodenum=m_SystemInfo.SystemTotalBusNum;
	int Lineno=i;   ///形成的是雅克比的哪一行
    if(i>=m_SystemInfo.SwingBusNo)
		Lineno=i-1; 

	FTIU[1]=1;
	double *temp=new double[2*nodenum+RE];/////用来存储将要消去行的元素
	
	for(int s=1;s<=2;s++)/////////////因为每次消去的的雅克比矩阵都有两行，所以需要两次循环
	{
	for(int k=1;k<=2*nodenum;k++)
		temp[k]=0;
	int jno=0;  /////雅克比矩阵元素的位置
	/////先读入下三角的元素
	for(int m=1;m<Lineno;m++)
	{
	for(k=JLJ[m];k<JLJ[m+1];k++)
	{
      if(JLI[k]==Lineno)  
	  { 
		 jno++;
		 if(1==s)
		 {	
	     temp[2*m-1]=JM[jno];
	     temp[2*m]=JL[jno];
		 }
		 if(2==s)
		 {
         temp[2*m-1]=JH[jno];
	     temp[2*m]=JN[jno];
		 }
	  }
	}
	}
	///////////////读入对角元素
	jno=jno+1;
	if(1==s)
	{
	temp[2*Lineno-1]=JM[jno];
	temp[2*Lineno]=JL[jno];
	}
    if(2==s)
	{
	temp[2*Lineno-1]=JH[jno];
	temp[2*Lineno]=JN[jno];
	}

	///////////读入上三角元素
	for(k=JI[Lineno];k<JI[Lineno+1];k++)
	{
      jno++;
	  int j=JJ[k];
	  if(1==s)
	  {	 
	  temp[2*j-1]=JM[jno];
	  temp[2*j]=JL[jno];
	  }
      if(2==s)
	  {
	  temp[2*j-1]=JH[jno];
	  temp[2*j]=JN[jno];
	  }
	}
	jno=0;
	int k2=0;
	if(1==s) k2=2*Lineno-1;
	if(2==s) k2=2*Lineno;

//	for(m=1;m<=2*(nodenum-1);m++)
//		if(temp[m]!=0) cout<<k2<<" "<<m<<" "<<temp[m]<<endl;
// 	cout<<"---------------"<<endl;

//  cout<<k2-1<<" "<<FTIU[k2-1]<<endl;

	//////消去运算

    for(k=1;k<k2;k++)
	{		
	  double a=0;
      if(temp[k]!=0) 
	  {
	   a=temp[k];   ////
	   temp[k]=0;

	   for(int n=k+1;n<=2*(nodenum-1);n++)
	   {
		//   if(n==k) 
		//	  {
		//		  temp[k]=0;
		//	  }
		   for(int x=FTIU[k];x<FTIU[k+1];x++)
		   {
			  int j=FTJU[x];
			  if(n==j)
			  {
				  temp[n]=temp[n]-a*FTU[x];			
			  }
			  else 
			  {
				  temp[n]=temp[n];
			  }
		   }
	     
	   } 
	   delta_PQV[k2]=delta_PQV[k2]-a*delta_PQV[k];		 
	  }
	}  

//////////////////delta_PQV还有问题	v 
//	cout<<delta_PQV[k2]<<endl;
 //     cout<<temp[k2]<<endl;
	 delta_PQV[k2]=delta_PQV[k2]/temp[k2];
//	cout<<k2<<" "<<delta_PQV[k2]<<endl;
 


	///////////////将消去的结果读入因子表中

    int num=0;  ///每行因子表非零元的个数
	if(1==s) 
	{
		num=FTIU[2*Lineno-1];
	}
	if(2==s) 
	{
		num=FTIU[2*Lineno];
	}
	int k1=0;
	if(1==s) k1=2*Lineno-1;
	if(2==s) k1=2*Lineno;
	for(k=k1+1;k<=2*(nodenum-1);k++)
	{
		if(temp[k]!=0)   
		{				
			FTJU[num]=k;
			FTU[num]=temp[k]/temp[k1];		
//			cout<<"(((((((((((((("<<num<<" "<<m<<endl;
//			cout<<FTIU[num]<<endl;
		    num++;
		}
	}

 //    cout<<"***"<<num<<endl;
		FTIU[k2+1]=num;
 //      cout<<FTIU[k2+1]<<endl;
//		cout<<endl;
////    cout<<2*Lineno-1<<"^^^"<<FTIU[2*Lineno-1]<<endl;
//		cout<<endl;
////    cout<<2*Lineno<<"^^^"<<FTIU[2*Lineno]<<endl;
//
//	   if(1==s)
//	   {
//		   cout<<k1<<endl;
//		for(int k=FTIU[2*Lineno-1];k<FTIU[2*Lineno];k++)
//			cout<<FTJU[k]<<" "<<FTU[k]<<endl;
//		cout<<"******************************"<<endl;
//	   }
//	   if(2==s)
//	   {
//		  cout<<k1<<endl;
//		  for(int k=FTIU[2*Lineno];k<FTIU[2*Lineno+1];k++)
//			   cout<<FTJU[k]<<" "<<FTU[k]<<endl;
//		  	cout<<"******************************"<<endl;
//  	   }
//	   
//		cout<<"--------------"<<endl;
	
 	}
	delete [] temp;
//	for(int u=1;u<=nodenum-1;u++)
//     cout<<JLI[u]<<endl;
//	cout<<"================"<<endl;

}

void Jacobi::GetNonZeroNum()
{
	extern SystemInfo m_SystemInfo;   ///系统信息
    extern LineInfo *m_LineInfo;       ///线路信息
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路 
	extern BusInfo *m_BusInfo;        //// 节点 

}
int Jacobi::GetNewInNum()
{
	Ymatrix.SpareYMatrix();
	extern SystemInfo m_SystemInfo;   ///系统信息
	int nodenum=m_SystemInfo.SystemTotalBusNum;

	/////利用导纳阵的稀疏存贮结构得到雅克比矩阵的存贮结构
	 // 求出JI，JJ
	 int swingbusno=m_SystemInfo.SwingBusNo;    ///平衡节点编号
	 int swingnum=Ymatrix.IY[swingbusno+1]-Ymatrix.IY[swingbusno];   //////平衡节点所在行的非零元个数
	 int lineno=m_SystemInfo.GeneralLineNum+m_SystemInfo.TransformerNum;  /////线路数目，非零元的总数目
     int UJnum=lineno-m_SystemInfo.MultiLineNum; ///非零元的总数目，总支路数-多回线的数目
     JI=new int[nodenum+RE];    ///按行存储
     JJ=new int[UJnum-swingnum+RE];    
     JLI=new int[UJnum-swingnum+RE];///将下三角元素按列存储
	 JLJ=new int[nodenum+RE]; 
	 
	 int swingJnum=0; ///平衡节点所在列的上三角元素个数
	 for(int i=1;i<=nodenum;i++)
	 {
		 int k=0;
		 if(i>=swingbusno) 
		 {
			 k=i+1;
			 JLJ[i]=JI[i]=Ymatrix.IY[k]-swingnum-swingJnum;
	//		 cout<<i<<" "<<JI[i]<<endl;
		 }
		 else
		 {
	     JLJ[i]=JI[i]=Ymatrix.IY[i]-swingJnum;		
	//	 cout<<i<<" "<<JI[i]<<endl;
		 for(int j=Ymatrix.IY[i];j<Ymatrix.IY[i+1];j++)
			 {
				 if(Ymatrix.JY[j]==swingbusno) 
				 {
                       swingJnum++;
				 }
			 }
		 }
   //    cout<<endl;
	 }

	  swingJnum=0;
	  int J=0;
	for(i=1;i<=nodenum;i++)
	{
		if(i==68) 
		{
	//		cout<<"ererere"<<endl;
		}
		if(i!=swingbusno)
		{		
		if(i==swingbusno+1) 
			{
				swingJnum=swingJnum+swingnum;
			}
        for(int k=Ymatrix.IY[i];k<Ymatrix.IY[i+1];k++)
		{
			int j=Ymatrix.JY[k];	
            if(j>swingbusno)
			{	
             JLI[k-swingJnum]=JJ[k-swingJnum]=j-1;
//	    	cout<<i<<" "<<JJ[k-swingJnum]<<" "<<k-swingJnum<<endl;
			}
			if(j<swingbusno) 
			{
               JLI[k-swingJnum]=JJ[k-swingJnum]=j;
//			   cout<<i<<" "<<JJ[k-swingJnum]<<" "<<k-swingJnum<<endl;
			}
	      if(j==swingbusno) 
			{
                    swingJnum++;
		  }
		}
//		cout<<"iiiiiiiiiiiiiiiiiiiiiii"<<endl;
		}
	}



//	 for(int s1=1;s1<nodenum;s1++)
//	 {
//		 cout<<JLJ[s1]<<"----"<<endl;
//	 for(int s2=JLJ[s1];s2<JLJ[s1+1];s2++)
//	  {
//         cout<<s1<<" "<<JLI[s2]<<endl;
//	  }
//	cout<<"------------------------"<<endl;
//	 }

    int *mark=new int[nodenum+RE]; /////消去过程中产生非零员的标志，0为不存在，1为存在
	int addnum=0;////新增非零元的个数
//	int *I_addnum;/////每一行新增元素的个数
//	int *J_addnum;///每行新增元素的行号
 //   I_addnum=new int[nodenum-1+RE];
	for(i=1;i<=nodenum;i++)
		mark[i]=0;
	/////用第一行非零员的列号初始化mark矩阵
    for(i=Ymatrix.IY[1];i<Ymatrix.IY[2];i++)
	{
		int j=Ymatrix.JY[i];
		mark[j]=1;
	}
	///////如果前i-1行中任意一个有第j列元素，第i行没有第j列元素，就增加一个 
	for(i=2;i<=nodenum;i++)
	{   
		bool IJ=0;   ///第i行第i-1列是否有非零元
		for(int m=1;m<JI[i];m++)
		{
			if(i==JJ[m]) IJ=1;  
		}
		if(1==IJ)
		{
		for(int j=i+1;j<=nodenum;j++)
		{
			int exit=0;
			if(1==mark[j]) 
			{
				for(int k=JI[i];k<JI[i+1];k++)
				{
					int m=JJ[k];
				    if(j==m)  exit++;      	
				}
			if(exit==0) addnum++;
 			}		
		}
      	for( j=JI[i];j<JI[i+1];j++)
		{
			int k=JJ[j];
            if(mark[k]==0)
			 {
			   mark[k]=1;  
			 }
		}
		}   
	}	
//	cout<<addnum<<"%%"<<endl;
	
	int Unum=4*(Ymatrix.IY[nodenum+1]-1);	
	int a=Unum+2*addnum+RE;
    FTU=new double[Unum+4*addnum+RE];  ////因子表的上三角元素 
    FTJU=new int[Unum+4*addnum+RE]; ///因子表上三角元素的列号
	FTIU=new int[2*nodenum+RE]; /////因子表的每行第一个非零元地址			

	 for(i=0;i<=2*nodenum;i++)
		 FTIU[i]=0;	 
	 delete [] mark;
     return addnum;

}

