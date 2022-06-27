#pragma once
#include"include/armadillo"
#include "factorial.h"
#include"material.h"
#include"CFDInfo.h"
#include"myNastran.h"
using namespace std;
using namespace arma;
using namespace arith;

//单个梁结构
class SingleBone
{
public:
	SingleBone() { BoneConnP_ID = 0; MainConnP_ID = 0; }

	mat BoneP;						//梁骨架节点
	mat BonePForce;					//骨架梁节点的力
	mat GetConnBone(mat MainBone);	//生成并返回与主梁连接的分支骨架梁
	mat BoneE(int fstID = 0);		//梁骨架单元(fstID=第一个参数的id)
	mat StructBoneP;				//记录有连接结构的梁分支结构 at GetConnBone
	mat StructBoneE;				//at BoneE
	mat StructBoneF;				//at BoneE
	mat ConnBoneP;

	int BoneConnP_ID;				//粱身上要与主梁连接的点编号
	int MainConnP_ID;				//主骨架上要与梁连接的点编号
};

//整体梁结构
class Bone
{
public:
	SingleBone MainBone;			//身体骨架梁
	vector<SingleBone> BrachBone;	//翼，尾梁

	mat Bone_Pf;					//节点气动力
	mat Bone_P;
	mat Bone_E;
	int CalcBone();					//计算生成总体梁节点 单元 节点气动力
	int SaveAsTcp(string DATname);	//将骨架结构打印为文本格式
	int SaveAsAbaqus(string DATname);//将骨架结构打印为abaqus格式
};

//-----梁壳组合单元构建--------
//单个梁壳组合部件
class ShellandBeam
{
public:
	string name;	//梁壳结构名称
	mat node;		//节点
	mat node_f;		//节点力
	mat elem_aero;	//气动壳单元
	mat elem_strc;	//结构梁单元

	mat beamID;		//记录梁结构
	mat elem_h;		//气动壳单元对应的厚度
	mat Xbar;		//x方向梁的位置
	mat Zbar;		//z方向梁的位置

	mat LinkPtId;	//连接点的编号,节点变更时需要更新
	mat ConnPtId;	//运行GetConnPtId后方可使用
	mat GetConnPtId(const mat& MainBone);

	int write_abaqus(ofstream& ofs);
private:

};

//整个梁壳组合结构
class Struct1D
{
public:
	SingleBone MainBone;			//身体骨架梁
	vector<ShellandBeam> BarchStruc;//分支梁壳结构

	void getConnPoint();
	
	int SaveAsAbaqus(string DATname);//将骨架结构打印为abaqus格式
	int SaveAsTecplt(string DatName);//将东西输出为Tecplot格式
	
	int SaveAsNastran(string DatName, int SOL = 101);//输出为Nastran格式 默认格式为静力学分析
};

