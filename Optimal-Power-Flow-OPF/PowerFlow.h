#ifndef _POWERFLOW_H
#define _POWERFLOW_H

class PowerFlow
{
public:
	PowerFlow();
	virtual~PowerFlow();
	void PFlow();       //结算潮流
	double GetDeltaV(); ///修正方程求解
	void OutputResult(char *,int *);  ///输出结果
	void Initialize();     ///初始化运行方式
	double *delta_Votltage;   ///电压修正值
    int MaxIteraNum;   ///最大迭代数目
	double IteraError; ///最大迭代误差
	int Itearnum;    //迭代次数
};
#endif