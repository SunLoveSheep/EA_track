#ifndef _POWERFLOW_H
#define _POWERFLOW_H

class PowerFlow
{
public:
	PowerFlow();
	virtual~PowerFlow();
	void PFlow();       //���㳱��
	double GetDeltaV(); ///�����������
	void OutputResult(char *,int *);  ///������
	void Initialize();     ///��ʼ�����з�ʽ
	double *delta_Votltage;   ///��ѹ����ֵ
    int MaxIteraNum;   ///��������Ŀ
	double IteraError; ///���������
	int Itearnum;    //��������
};
#endif