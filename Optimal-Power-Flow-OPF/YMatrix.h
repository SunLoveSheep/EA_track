/////导纳矩阵类 
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
	void SpareYMatrix();   ///形成稀疏导纳矩阵
	void Output(int *);
	CYMatrix *DY; ////导纳矩阵对角元
	CYMatrix *UY; ////导纳矩阵上三角元素
	int *JY;   ///上三角元素的列号
	int *IY;   ///上三角元素每行第一个非零元素在UY中的位置(首地址) 

};
#endif

