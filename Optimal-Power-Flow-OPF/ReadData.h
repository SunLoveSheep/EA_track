//////����������  

#ifndef _READDATA_H_
#define _READDATA_H_ 
#include <fstream.h>
class CReadData
{
public: 
    CReadData();
	virtual~CReadData();
	ifstream datain;
	/////������Ʊ���
	int BusOptimType;  ////�ڵ��Ż���ŷ�ʽ
	int MaxIteraNum;   ///��������Ŀ
	double IteraError; ///���������
    double BasePower;  ///��׼����
	double BaseVoltage;  ///��׼��ѹ
  
	void Read(char *);   //��������
};
#endif