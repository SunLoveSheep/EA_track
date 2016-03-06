////////////////////////////////////////////////////////////////////////////////////////
///////   导纳矩阵 形成稀疏导纳阵                               //////////
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
	extern SystemInfo m_SystemInfo;   ///系统信息
    extern LineInfo *m_LineInfo;       ///线路信息
	extern TransformerInfo *m_TransformerInfo;  ///变压器支路 
	extern BusInfo *m_BusInfo;        //// 节点 
	
	int nodenum=m_SystemInfo.SystemTotalBusNum;
	int linenum=m_SystemInfo.GeneralLineNum;   ////一般支路数目
	int totallinenum=m_SystemInfo.GeneralLineNum+m_SystemInfo.LandBranchNum;  ///总支路数(一般+接地)
	int transformernum=m_SystemInfo.TransformerNum;  ///变压器数
	int UYnum=m_SystemInfo.GeneralLineNum+m_SystemInfo.TransformerNum;///也就是总线路数据.上三角元素的个数
    
	////消除双回线的影响
   	int *busIno=new int[UYnum+RE];  //////记录线路编号较小端点
	int *busJno=new int[UYnum+RE];  //////记录线路编号较大端点
	int mullinenun=0;  ////多回线的数目
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
	UYnum=UYnum-mullinenun;  //得到实际的上三角元素个数	
	
	UY=new CYMatrix[UYnum+RE];   ////导纳矩阵上三角元素
	DY=new CYMatrix[nodenum+RE]; ////导纳矩阵对角元
    JY=new int[UYnum+RE];        ///上三角元素的列号
	IY=new int[nodenum+1+RE];      ///上三角元素每行第一个非零元素在UY中的位置(首地址)
	for(i=1;i<=nodenum+1;i++)
		DY[i].G=DY[i].B=IY[i]=0;
	IY[1]=1;
	///////形成导纳阵
	int JYno=0;///上三角列元素的序号
    int *JY_lineno=new int[UYnum+RE];  ///记录上三角元素的行号，用来判断是否换行，剔除双回线
	
	for(int k=1;k<=nodenum;k++)
	{	
		int I,J; //两个编号中小的为I，大的为J
		for(int m=1;m<=linenum;m++)   ////普通线路对应的节点
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
		for(m=1;m<=transformernum;m++)   ///变压器对应的节点
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
	////补充接地支路
	for(k=linenum+1;k<=totallinenum;k++)
	{
		int i=m_LineInfo[k].BusINo;
		int j=m_LineInfo[k].BusJNo;
		DY[i].G+=m_LineInfo[k].R;
		DY[j].B+=m_LineInfo[k].X;
	}
//for(k=1;k<=nodenum;k++)
// cout<<k<<"  "<<DY[k].G<<"+"<<DY[k].B<<"j"<<endl;
	//////////////////对非对角元每行按照行号从小到大排序
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
		dataout.open("输出文件//5节点//导纳矩阵.txt");
	}
	else if(nodenum==14)
	{
	dataout.open("输出文件//14节点//导纳矩阵.txt");
	}
	else if(nodenum==30)
	{
	dataout.open("输出文件//30节点//导纳矩阵.txt");
	}
    else if(nodenum==57)
	{
	dataout.open("输出文件//57节点//导纳矩阵.txt");
	}
	else if(nodenum==118)
	{
	dataout.open("输出文件//118节点//导纳矩阵.txt");
	}

    dataout<<"导纳矩阵对角元:"<<endl;
for(k=1;k<=nodenum;k++)
{
	dataout<<k<<"  "<<DY[k].G<<"+"<<DY[k].B<<"j"<<endl;
}	
dataout<<endl;
dataout<<"导纳矩阵非对角元:"<<endl;
     for(k=1;k<=nodenum;k++)
 	{
 		for(i=IY[k];i<IY[k+1];i++)
 			dataout<<k<<JY[i]<<" "<<UY[i].G<<"+"<<UY[i].B<<"j"<<"      ";
 		dataout<<endl;
   	}
	dataout<<"所对应的节点是进行编号之后的节点";
	dataout.close();
	delete [] busIno;
	delete [] busJno;
	delete [] JY_lineno;
}