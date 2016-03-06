//////////////////////////////////////////////////////////////////////////////////////////
////////   �ſ˱���            //////
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef  _JACOBI_H
#define _JACOBI_H
#include "YMatrix.h"
class  Jacobi
{
public:
     Jacobi();
	 virtual~Jacobi();
     double *JH;   ///�й��Ե�ѹʵ����
	 double *JN;   ///�й��Ե�ѹ�鲿��
	 double *JM;   ///�޹��Ե�ѹ�鲿��
	 double *JL;   ///�޹��Ե�ѹ�鲿��	 
	 
	 void FormJacobi();
	 void FormFactorTable(int);
	 void GetNonZeroNum();////�õ���ȥ�����������ӵķ���Ԫ���Ӷ�Ԥ���ռ�
	 int GetNewInNum();  ////������ȥ������ע��Ԫ,i=1��ʾ�����������i=2��ʾ��ע��Ԫ����

	 int *JI;    ///�ſ˱Ⱦ������׵�ַ,���д洢
	 int *JJ;    ///�ſ˱Ⱦ���Ԫ���кţ����д洢
	 int *JLI;   ///�ſ˱Ⱦ���Ԫ�ص��кţ����д洢
     int *JLJ;   ///�ſ˱Ⱦ���Ԫ�ص����׵�ַ�����д洢
	 	 
	 YMatrix Ymatrix;
	 int *FTIU; /////���ӱ��ÿ�е�һ������Ԫ��ַ
	 int *FTJU; ///���ӱ�������Ԫ�ص��к�
     double *FTU;  ////���ӱ��������Ԫ��
     double *delta_PQV;


};
#endif