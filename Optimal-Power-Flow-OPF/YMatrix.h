/////���ɾ����� 
#ifndef _YMATRIX_H
#define _YMATRIX_H
#include "GeneralInfo.h"

struct CYMatrix
{
	double G;
	double B;
};
class YMatrix
{
public:
	YMatrix();
	virtual~YMatrix();
	void SpareYMatrix();   ///�γ�ϡ�赼�ɾ���
	void Output(int *);
	CYMatrix *DY; ////���ɾ���Խ�Ԫ
	CYMatrix *UY; ////���ɾ���������Ԫ��
	int *JY;   ///������Ԫ�ص��к�
	int *IY;   ///������Ԫ��ÿ�е�һ������Ԫ����UY�е�λ��(�׵�ַ) 

};
#endif

