//////读入数据类  

#ifndef _READDATA_H_
#define _READDATA_H_ 
#include <fstream.h>
class CReadData
{
public: 
    CReadData();
	virtual~CReadData();
	ifstream datain;
	/////计算控制变量
	int BusOptimType;  ////节点优化编号方式
	int MaxIteraNum;   ///最大迭代数目
	double IteraError; ///最大迭代误差
    double BasePower;  ///基准容量
	double BaseVoltage;  ///基准电压
  
	void Read(char *);   //读入数据
};
#endif