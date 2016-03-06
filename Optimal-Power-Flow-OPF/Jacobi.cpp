///////////////////////////////////////////////////////////////////////////////////////////
/////////////          ����������������̵Ĵ���                    ////////////////////////
////////////  �γ��ſ˱Ⱦ����γɵ������ӱ�ȷ����ȥ�������γɵ�ע��Ԫ��λ��      ///////
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

    extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ� 
	
	int nodenum=m_SystemInfo.SystemTotalBusNum;
    delta_PQV=new double[2*(nodenum-1)+RE];


	for(int i=1;i<=nodenum;i++)     /////�����γ��ſ˱Ⱦ���
	{	
	double ai=0,bi=0;////�������Խ�Ԫʱ��ai=(Gij*ej-Bij*fj)�ۼӺͣ�bi=(Gij*fj+Bij*ej)�ۼӺ�
    int swingbusno=m_SystemInfo.SwingBusNo;
    if(i==swingbusno) i=i+1;   ///����ƽ��ڵ�
//	cout<<Ymatrix.IY[i]<<endl;
	    int num=Ymatrix.IY[i+1]-Ymatrix.IY[i]; ////�ſ˱Ⱦ���ÿ�з���Ԫ�ĸ���
	    int JJno=0;//// ÿһ��Ԫ����JH,JN,JM,JL���±� 
		JH=new double[num+RE];   ////���������ٵĿռ䲻׼ȷ
		JN=new double[num+RE];
		JM=new double[num+RE];
		JL=new double[num+RE];

	////�����ſ˱Ⱦ��������ǲ���
	  for(int m=1;m<Ymatrix.IY[i];m++)
	  {
		  int j=Ymatrix.JY[m];
	      double Gij=Ymatrix.UY[m].G;
		  double Bij=Ymatrix.UY[m].B;
		
		  if(j==i) 
		  {	
			  ////ȷ������һ��
	 		  int I;
			  for(int x=1;x<i;x++)
			  {
				  if(m>=Ymatrix.IY[x]&&m<=Ymatrix.IY[x+1]-1) I=x;
			  }
			  ai+=Gij*m_BusInfo[I].Voltage_e-Bij*m_BusInfo[I].Voltage_f;
			  bi+=Gij*m_BusInfo[I].Voltage_f+Bij*m_BusInfo[I].Voltage_e;	  
			  ///���������ǵķǶԽ�Ԫ
		//	if(I!=swingbusno&&j!=swingbusno+1)
			  if(I!=swingbusno)
			  {
              JJno++;
			  JH[JJno]=-(Gij*m_BusInfo[i].Voltage_e+Bij*m_BusInfo[i].Voltage_f);
              JN[JJno]=Bij*m_BusInfo[i].Voltage_e-Gij*m_BusInfo[i].Voltage_f;
			  if(m_BusInfo[i].BusType==-1)///PV�ڵ�ĵ�ѹi��j�ڵ��ѹƫ��Ϊ0
			  {
			  JM[JJno]=0;
              JL[JJno]=0;
			  }
              if(m_BusInfo[i].BusType==1)  //PQ�ڵ�
			  {
			  JM[JJno]=JN[JJno];
              JL[JJno]=-JH[JJno];
			  }
			  }
		  }	
	  }

      //�����ſ˱Ⱦ��������ǲ���	  
      int JJDno=JJno+1;  ///�Խ�Ԫ���±�
	  JJno=JJno+1;
	  for(int k=Ymatrix.IY[i];k<Ymatrix.IY[i+1];k++)   
	  {		
		  int j=Ymatrix.JY[k];
		  double Gij=Ymatrix.UY[k].G;
		  double Bij=Ymatrix.UY[k].B;
		  ai+=Gij*m_BusInfo[j].Voltage_e-Bij*m_BusInfo[j].Voltage_f;
		  bi+=Gij*m_BusInfo[j].Voltage_f+Bij*m_BusInfo[j].Voltage_e;
		  //////���������ǵķǶԽ�Ԫ
		  if(j!=swingbusno)
		  {
          JJno++;
		  JH[JJno]=-(Gij*m_BusInfo[i].Voltage_e+Bij*m_BusInfo[i].Voltage_f);
          JN[JJno]=Bij*m_BusInfo[i].Voltage_e-Gij*m_BusInfo[i].Voltage_f;
	      if(m_BusInfo[i].BusType==-1)   ///PV�ڵ�ĵ�ѹi��j�ڵ��ѹƫ��Ϊ0
			  {
			  JM[JJno]=0;
              JL[JJno]=0;
			  }
          if(m_BusInfo[i].BusType==1) //PQ�ڵ�
			  {
			  JM[JJno]=JN[JJno];
              JL[JJno]=-JH[JJno];
			  }
		  }
	  } 
      ///�����ſ˱Ⱦ���Խ�Ԫ��
	  double Gii=Ymatrix.DY[i].G;
	  double Bii=Ymatrix.DY[i].B;
	  ai+=Gii*m_BusInfo[i].Voltage_e-Bii*m_BusInfo[i].Voltage_f;
      bi+=Gii*m_BusInfo[i].Voltage_f+Bii*m_BusInfo[i].Voltage_e;

	  JH[JJDno]=-ai-(Gii*m_BusInfo[i].Voltage_e+Bii*m_BusInfo[i].Voltage_f);
	  JN[JJDno]=-bi+(Bii*m_BusInfo[i].Voltage_e-Gii*m_BusInfo[i].Voltage_f);
	  int sno;///�������������±�
	  sno=i;
	  if(i>swingbusno) sno=i-1;
	 delta_PQV[2*sno]=m_BusInfo[i].PG-m_BusInfo[i].PD-m_BusInfo[i].Voltage_e*ai-m_BusInfo[i].Voltage_f*bi;
	  
	 if(m_BusInfo[i].BusType==1)  //PQ�ڵ�
	  {
	  JM[JJDno]=bi+(Bii*m_BusInfo[i].Voltage_e-Gii*m_BusInfo[i].Voltage_f);
	  JL[JJDno]=-ai+(Gii*m_BusInfo[i].Voltage_e+Bii*m_BusInfo[i].Voltage_f);
	 delta_PQV[2*sno-1]=m_BusInfo[i].QG-m_BusInfo[i].QD-m_BusInfo[i].Voltage_f*ai+m_BusInfo[i].Voltage_e*bi;
	  }
     if(m_BusInfo[i].BusType==-1)   //pv�ڵ�
	 {
	  JM[JJDno]=-2*m_BusInfo[i].Voltage_e;
	  JL[JJDno]=-2*m_BusInfo[i].Voltage_f;
	  delta_PQV[2*sno-1]=m_BusInfo[i].VolMod-(m_BusInfo[i].Voltage_e*m_BusInfo[i].Voltage_e+m_BusInfo[i].Voltage_f*m_BusInfo[i].Voltage_f);
	 }
// 	 cout<<sno<<"PQV"<<delta_PQV[2*sno-1]<<"  "<<delta_PQV[2*sno]<<endl;
	 if(i==swingbusno-1&&swingbusno==nodenum) i=i+1; ///����ƽ��ڵ�Ϊ���һ���ڵ�����
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

//void Jacobi::FormFactorTable(int i,double *f_JH,double *f_JM,double *f_JN,double *f_JL)  ///��ȥ
void Jacobi::FormFactorTable(int i)
{ 

    extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
   
    int nodenum=m_SystemInfo.SystemTotalBusNum;
	int Lineno=i;   ///�γɵ����ſ˱ȵ���һ��
    if(i>=m_SystemInfo.SwingBusNo)
		Lineno=i-1; 

	FTIU[1]=1;
	double *temp=new double[2*nodenum+RE];/////�����洢��Ҫ��ȥ�е�Ԫ��
	
	for(int s=1;s<=2;s++)/////////////��Ϊÿ����ȥ�ĵ��ſ˱Ⱦ��������У�������Ҫ����ѭ��
	{
	for(int k=1;k<=2*nodenum;k++)
		temp[k]=0;
	int jno=0;  /////�ſ˱Ⱦ���Ԫ�ص�λ��
	/////�ȶ��������ǵ�Ԫ��
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
	///////////////����Խ�Ԫ��
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

	///////////����������Ԫ��
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

	//////��ȥ����

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

//////////////////delta_PQV��������	v 
//	cout<<delta_PQV[k2]<<endl;
 //     cout<<temp[k2]<<endl;
	 delta_PQV[k2]=delta_PQV[k2]/temp[k2];
//	cout<<k2<<" "<<delta_PQV[k2]<<endl;
 


	///////////////����ȥ�Ľ���������ӱ���

    int num=0;  ///ÿ�����ӱ����Ԫ�ĸ���
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
	extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ� 

}
int Jacobi::GetNewInNum()
{
	Ymatrix.SpareYMatrix();
	extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
	int nodenum=m_SystemInfo.SystemTotalBusNum;

	/////���õ������ϡ������ṹ�õ��ſ˱Ⱦ���Ĵ����ṹ
	 // ���JI��JJ
	 int swingbusno=m_SystemInfo.SwingBusNo;    ///ƽ��ڵ���
	 int swingnum=Ymatrix.IY[swingbusno+1]-Ymatrix.IY[swingbusno];   //////ƽ��ڵ������еķ���Ԫ����
	 int lineno=m_SystemInfo.GeneralLineNum+m_SystemInfo.TransformerNum;  /////��·��Ŀ������Ԫ������Ŀ
     int UJnum=lineno-m_SystemInfo.MultiLineNum; ///����Ԫ������Ŀ����֧·��-����ߵ���Ŀ
     JI=new int[nodenum+RE];    ///���д洢
     JJ=new int[UJnum-swingnum+RE];    
     JLI=new int[UJnum-swingnum+RE];///��������Ԫ�ذ��д洢
	 JLJ=new int[nodenum+RE]; 
	 
	 int swingJnum=0; ///ƽ��ڵ������е�������Ԫ�ظ���
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

    int *mark=new int[nodenum+RE]; /////��ȥ�����в�������Ա�ı�־��0Ϊ�����ڣ�1Ϊ����
	int addnum=0;////��������Ԫ�ĸ���
//	int *I_addnum;/////ÿһ������Ԫ�صĸ���
//	int *J_addnum;///ÿ������Ԫ�ص��к�
 //   I_addnum=new int[nodenum-1+RE];
	for(i=1;i<=nodenum;i++)
		mark[i]=0;
	/////�õ�һ�з���Ա���кų�ʼ��mark����
    for(i=Ymatrix.IY[1];i<Ymatrix.IY[2];i++)
	{
		int j=Ymatrix.JY[i];
		mark[j]=1;
	}
	///////���ǰi-1��������һ���е�j��Ԫ�أ���i��û�е�j��Ԫ�أ�������һ�� 
	for(i=2;i<=nodenum;i++)
	{   
		bool IJ=0;   ///��i�е�i-1���Ƿ��з���Ԫ
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
    FTU=new double[Unum+4*addnum+RE];  ////���ӱ��������Ԫ�� 
    FTJU=new int[Unum+4*addnum+RE]; ///���ӱ�������Ԫ�ص��к�
	FTIU=new int[2*nodenum+RE]; /////���ӱ��ÿ�е�һ������Ԫ��ַ			

	 for(i=0;i<=2*nodenum;i++)
		 FTIU[i]=0;	 
	 delete [] mark;
     return addnum;

}

