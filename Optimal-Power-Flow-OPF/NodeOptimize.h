//////////�ڵ��Ż���     ////////

#ifndef _NODEOPTIMIZE_H
#define _NODEOPTIMIZE_H

#include "GeneralInfo.h"

class NodeOptimize
{
public:
	NodeOptimize();
	virtual~NodeOptimize();
	void StaticOptimize();   ///��̬�Ż���� 
	void HalfDymOptimize();  ///�붯̬�Ż����
	void DynamicOptimize();  ///��̬�Ż����
	void AdjustNodeNo();     ///�����ڵ���
	int OptimizeType;        ///�ڵ��Ż�����
	int *old_new;   ////�ϱ�Ŷ�Ӧ���±�ţ���old_new[i]=j��ʾԭ�ڵ�����i���½ڵ�����j
	int *new_old;   ////�±�Ŷ�Ӧ���ϱ�ţ���old_new[i]=j��ʾ�½ڵ�����i��ԭ�ڵ�����j
};
#endif