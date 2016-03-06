////////////////////////////////////////////////////////////////////////////////////////
///////   ���ɾ��� �γ�ϡ�赼����                               //////////
//////////////////////////////////////////////////////////////////////////

#include "YMatrix.h"
#include <iostream.h>
#include <fstream.h>
#define RE 10
YMatrix::YMatrix()
{
	
}
YMatrix::~YMatrix()
{
	
}
void YMatrix::SpareYMatrix()
{
	extern SystemInfo m_SystemInfo;   ///ϵͳ��Ϣ
    extern LineInfo *m_LineInfo;       ///��·��Ϣ
	extern TransformerInfo *m_TransformerInfo;  ///��ѹ��֧· 
	extern BusInfo *m_BusInfo;        //// �ڵ� 
	
	int nodenum=m_SystemInfo.SystemTotalBusNum;
	int linenum=m_SystemInfo.GeneralLineNum;   ////һ��֧·��Ŀ
	int totallinenum=m_SystemInfo.GeneralLineNum+m_SystemInfo.LandBranchNum;  ///��֧·��(һ��+�ӵ�)
	int transformernum=m_SystemInfo.TransformerNum;  ///��ѹ����
	int UYnum=m_SystemInfo.GeneralLineNum+m_SystemInfo.TransformerNum;///Ҳ��������·����.������Ԫ�صĸ���
    
	////����˫���ߵ�Ӱ��
   	int *busIno=new int[UYnum+RE];  //////��¼��·��Ž�С�˵�
	int *busJno=new int[UYnum+RE];  //////��¼��·��Žϴ�˵�
	int mullinenun=0;  ////����ߵ���Ŀ
    for(int i=1;i<=linenum;i++)
	{
		busIno[i]=m_LineInfo[i].BusINo;
		busJno[i]=m_LineInfo[i].BusJNo;
	}
	for(i=linenum+1;i<=UYnum;i++)
	{	
		
		busIno[i]=m_TransformerInfo[i-linenum].BusINo;	
		busJno[i]=m_TransformerInfo[i-linenum].BusJNo;
	}
	for(i=1;i<=UYnum;i++)
	{
		for(int j=1;j<=UYnum;j++)
			if(i!=j&&busIno[i]==busIno[j]&&busJno[i]==busJno[j]) 
				mullinenun++;
	}	
    mullinenun=mullinenun/2;    
	m_SystemInfo.MultiLineNum=mullinenun;
	UYnum=UYnum-mullinenun;  //�õ�ʵ�ʵ�������Ԫ�ظ���	
	
	UY=new CYMatrix[UYnum+RE];   ////���ɾ���������Ԫ��
	DY=new CYMatrix[nodenum+RE]; ////���ɾ���Խ�Ԫ
    JY=new int[UYnum+RE];        ///������Ԫ�ص��к�
	IY=new int[nodenum+1+RE];      ///������Ԫ��ÿ�е�һ������Ԫ����UY�е�λ��(�׵�ַ)
	for(i=1;i<=nodenum+1;i++)
		DY[i].G=DY[i].B=IY[i]=0;
	IY[1]=1;
	///////�γɵ�����
	int JYno=0;///��������Ԫ�ص����
    int *JY_lineno=new int[UYnum+RE];  ///��¼������Ԫ�ص��кţ������ж��Ƿ��У��޳�˫����
	
	for(int k=1;k<=nodenum;k++)
	{	
		int I,J; //���������С��ΪI�����ΪJ
		for(int m=1;m<=linenum;m++)   ////��ͨ��·��Ӧ�Ľڵ�
		{
			int i=m_LineInfo[m].BusINo;
			int j=m_LineInfo[m].BusJNo;
			
			if(i<j) { I=i; J=j;}
			else    { I=j; J=i;} 
			if(k==I)
			{
			    double mul_G=0,mul_B=0;
				if(JY[JYno]==J&&JY_lineno[JYno]==I)
				{
				 mul_G=UY[JYno].G;
                 mul_B=UY[JYno].B;
				}

				if(JY[JYno]!=J||JY_lineno[JYno]!=I) JYno++;
				JY[JYno]=J;
				JY_lineno[JYno]=I;
				double r=m_LineInfo[m].R;
				double x=m_LineInfo[m].X;
				double z=r*r+x*x;       
				UY[JYno].G=-r/z;
				UY[JYno].B=x/z;
				
				DY[i].G+=-UY[JYno].G;
				DY[i].B+=-UY[JYno].B+m_LineInfo[m].Yk;
				
				DY[j].G+=-UY[JYno].G;
				DY[j].B+=-UY[JYno].B+m_LineInfo[m].Yk;
				
				if(JY[JYno]==J&&JY_lineno[JYno]==I)
				{
	              UY[JYno].G+=mul_G;
				  UY[JYno].B+=mul_B;
				}
			}		
		}
		for(m=1;m<=transformernum;m++)   ///��ѹ����Ӧ�Ľڵ�
		{
			int i=m_TransformerInfo[m].BusINo;
			int j=m_TransformerInfo[m].BusJNo;
			if(i<j) { I=i; J=j;}
			else    { I=j; J=i;} 
			
			if(k==I)
			{       
				double mul_G=0,mul_B=0;
				if(JY[JYno]==J&&JY_lineno[JYno]==I)
				{
				 mul_G=UY[JYno].G;
                 mul_B=UY[JYno].B;
				}

				if(JY[JYno]!=J||JY_lineno[JYno]!=I) JYno++;
				JY[JYno]=J;	
				JY_lineno[JYno]=I;
				double r=m_TransformerInfo[m].R;
				double x=m_TransformerInfo[m].X;
				double z=r*r+x*x;    
				double TK=m_TransformerInfo[m].Ratio;
				//UY[JYno].G=(-r/z)/TK;
				//UY[JYno].B=(x/z)/TK;		

				UY[JYno].G=(-r/z)/TK;
				UY[JYno].B=(x/z)/TK;


				DY[i].G+=-UY[JYno].G+((TK-1)/TK)*(r/z);
				DY[i].B+=-UY[JYno].B+((TK-1)/TK)*(-x/z);
				
				DY[j].G+=-UY[JYno].G+(1-TK)/(TK*TK)*(r/z);
				DY[j].B+=-UY[JYno].B+(1-TK)/(TK*TK)*(-x/z);
				if(JY[JYno]==J&&JY_lineno[JYno]==I)
				{
	              UY[JYno].G+=mul_G;
				  UY[JYno].B+=mul_B;
				}
			}		
		}
        IY[k+1]=1+JYno;

//			cout<<k<<" "<<IY[k]<<"  " <<JY[JYno]<<endl;	
    }
//	for(k=1;k<=nodenum;k++)
//	{
//		for(int i=IY[k];i<IY[k+1];i++)
//		{
//			int j=JY[i];
//			cout<<k<<" "<<j<<" "<<i<<endl;
//		}
//		cout<<endl;
//	}
	////����ӵ�֧·
	for(k=linenum+1;k<=totallinenum;k++)
	{
		int i=m_LineInfo[k].BusINo;
		int j=m_LineInfo[k].BusJNo;
		DY[i].G+=m_LineInfo[k].R;
		DY[j].B+=m_LineInfo[k].X;
	}
//for(k=1;k<=nodenum;k++)
// cout<<k<<"  "<<DY[k].G<<"+"<<DY[k].B<<"j"<<endl;
	//////////////////�ԷǶԽ�Ԫÿ�а����кŴ�С��������
	for(k=1;k<=nodenum;k++)
	{
		for(int j=IY[k];j<IY[k+1];j++)
		{
			for(i=IY[k];i<IY[k+1]-1;i++)
			{
				if(JY[i]>JY[i+1])
				{
					int temp1=JY[i];
					JY[i]=JY[i+1];
					JY[i+1]=temp1;
					double temp2=UY[i].G;
					UY[i].G=UY[i+1].G;
					UY[i+1].G=temp2;
					temp2=UY[i].B;
					UY[i].B=UY[i+1].B;
					UY[i+1].B=temp2;
				}
			}
			
		}
	}
//     for(k=1;k<=nodenum;k++)
//	{
//		for(i=IY[k];i<IY[k+1];i++)
//			cout<<k<<JY[i]<<" "<<UY[i].G<<"+"<<UY[i].B<<"j"<<"      ";
//		cout<<endl;
 //  	}
	ofstream dataout;
	if(nodenum==5)
	{
		dataout.open("����ļ�//5�ڵ�//���ɾ���.txt");
	}
	else if(nodenum==14)
	{
	dataout.open("����ļ�//14�ڵ�//���ɾ���.txt");
	}
	else if(nodenum==30)
	{
	dataout.open("����ļ�//30�ڵ�//���ɾ���.txt");
	}
    else if(nodenum==57)
	{
	dataout.open("����ļ�//57�ڵ�//���ɾ���.txt");
	}
	else if(nodenum==118)
	{
	dataout.open("����ļ�//118�ڵ�//���ɾ���.txt");
	}

    dataout<<"���ɾ���Խ�Ԫ:"<<endl;
for(k=1;k<=nodenum;k++)
{
	dataout<<k<<"  "<<DY[k].G<<"+"<<DY[k].B<<"j"<<endl;
}	
dataout<<endl;
dataout<<"���ɾ���ǶԽ�Ԫ:"<<endl;
     for(k=1;k<=nodenum;k++)
 	{
 		for(i=IY[k];i<IY[k+1];i++)
 			dataout<<k<<JY[i]<<" "<<UY[i].G<<"+"<<UY[i].B<<"j"<<"      ";
 		dataout<<endl;
   	}
	dataout<<"����Ӧ�Ľڵ��ǽ��б��֮��Ľڵ�";
	dataout.close();
	delete [] busIno;
	delete [] busJno;
	delete [] JY_lineno;
}