//////////////////////////////////////////////////////////////////////////////////////////
////////   雅克比类            //////
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef  _JACOBI_H
#define _JACOBI_H
#include "YMatrix.h"
class  Jacobi
{
public:
     Jacobi();
	 virtual~Jacobi();
     double *JH;   ///有功对电压实部求导
	 double *JN;   ///有功对电压虚部求导
	 double *JM;   ///无功对电压虚部求导
	 double *JL;   ///无功对电压虚部求导	 
	 
	 void FormJacobi();
	 void FormFactorTable(int);
	 void GetNonZeroNum();////得到消去过程中新增加的非零元，从而预留空间
	 int GetNewInNum();  ////关于消去过程中注入元,i=1表示仅计算个数，i=2表示将注入元插入

	 int *JI;    ///雅克比矩阵行首地址,按行存储
	 int *JJ;    ///雅克比矩阵元素列号，按行存储
	 int *JLI;   ///雅克比矩阵元素的行号，按列存储
     int *JLJ;   ///雅克比矩阵元素的列首地址，按列存储
	 	 
	 YMatrix Ymatrix;
	 int *FTIU; /////因子表的每行第一个非零元地址
	 int *FTJU; ///因子表上三角元素的列号
     double *FTU;  ////因子表的上三角元素
     double *delta_PQV;


};
#endif