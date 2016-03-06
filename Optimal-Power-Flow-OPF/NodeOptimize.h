//////////节点优化类     ////////

#ifndef _NODEOPTIMIZE_H
#define _NODEOPTIMIZE_H

#include "GeneralInfo.h"

class NodeOptimize
{
public:
	NodeOptimize();
	virtual~NodeOptimize();
	void StaticOptimize();   ///静态优化编号 
	void HalfDymOptimize();  ///半动态优化编号
	void DynamicOptimize();  ///动态优化编号
	void AdjustNodeNo();     ///调整节点编号
	int OptimizeType;        ///节点优化类型
	int *old_new;   ////老编号对应的新编号，即old_new[i]=j表示原节点编号是i，新节点编号是j
	int *new_old;   ////新编号对应的老编号，即old_new[i]=j表示新节点编号是i，原节点编号是j
};
#endif