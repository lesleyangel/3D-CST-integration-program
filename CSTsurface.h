#pragma once
#include <iostream>
#include "include/armadillo"
#include "Bone.h"
#include "Curve.h"
//using namespace arma;
//using namespace std;
//using namespace arith;



//网格划分的设定
class MeshDefine
{
public:
	enum Style
	{
		normal,	 //正规的
		oblique, //倾斜的
	};
	int NFaiU;				//上表面横向网格数
	int NFaiL;				//下表面横向网格数
	int NEta;				//轴向网格数
	int NHeight;			//高度方向网格数
	mat MeshRefineRatio;	//网格密度修正因子
};

//单个cst曲面类
class CSTsurface
{
public:
	//------------------形状参数--------------------
	mat Origin;				//初始位置（X偏移，Y偏移，Z偏移）
	mat Rotation;			//旋转因子（X 轴 ，Y 轴 ，Z 轴 ）
	mat Slope;				//倾斜因子（X正向，Y正向，Z正向）
	mat Scale;				//缩放因子(X比例,YUPP比例,YLOW比例,Z比例)
	mat Length;				//缩放因子(X比例,YUPP比例,YLOW比例,Z比例)
	ClassFunc class_upp;
	ClassFunc class_low;
	mat NS;					//头部截面形状因子（上表面N1/N2,下表面N1/N2）
	mat NE;					//尾部截面形状因子（上表面N1/N2,下表面N1/N2）
	mat M;					//侧面导引形状因子（左导引M1/M2,右导引M1/M2）
	mat T;					//中面导引形状因子（上导引T1/T2,下导引T1/T2）
	double Ratio;			//形状因子过渡系数
	//------------------网格参数--------------------
	int NFaiU;				//上表面横向网格数
	int NFaiL;				//下表面横向网格数
	int NEta;				//轴向网格数
	int NHeight;			//高度方向网格数

	bool Is_MeshFront;		//是否缝合头部
	bool Is_MeshBack;		//是否缝合尾部
	mat BUPP;				//上表面控制因子
	mat BLOW;				//下表面控制因子
	mat DUPP;				//侧向表面控制因子
	mat DLOW;				//侧向表面控制因子
	int PriFuncType;		//基函数类型（=0 Bernstein基函数，=i B样条基函数）

	mat MeshRefineRatio;	//网格密度修正因子
	//------------------网格单元参数--------------------
	Grid GridUpp;	Grid GridLow;	MeshInfo Info;
	mat N1mUPP;//过渡截面参数 上表面左
	mat N2mUPP;//过渡截面参数 上表面右
	mat N1mLOW;//过渡截面参数 下表面左
	mat N2mLOW;//过渡截面参数 下表面右
	mat FaiU;	mat FaiL;		mat EtaU;	mat EtaL;mat Eta;
	//------------------网格生成--------------------
	CSTsurface();			//初始化参数
	void CST3D();			//计算生成CST曲面网格信息	
	void EvalCurvature();	//曲率评估 计算得出曲面曲率值	
	void GridRefine();		//网格修正（根据曲率修正网格）
	//------------------一维结构--------------------
	mat MakeBone(double bias1 = 0,double bias2 = 0);//(偏置量[-1,1])生成该曲面的梁结构点*****4.12更新
	mat ClacBonePForce();	//计算骨架梁节点上的力和力矩
	mat ClacBoneS();		//计算骨架梁节点对应的壳面积(还没写)
	//------------------二维结构--------------------
	int MakeStruct2D(mat X = mat(0, 0), mat Z = mat(0,0));//生成二维结构单元
	Grid StructPE;
	mat siteZend;//单个曲面末端与下一个曲面的连接交界面的坐标
	//------------------二维结构--------------------
	//新的生成方式：主要定义结构梁不同，定义均在气动网格基础上
	int Struct2dXnum, Struct2dZnum;
	int MakeStruct2D_byNum();//生成二维结构单元
	
	vector<PointSet> node_2D;//存储结构壳单元所需要的节点，最外层vector序号是一层节点

	
	//------------------单面壳&加强梁结构--------------------
	int MakeStruct1D(mat X = mat(0, 0), mat Z = mat(0, 0));//
	//mat point1D;
	mat elem_aero;
	mat elem_strc;
	ShellandBeam aero_strc;
	//------------------网格修正新方法--------------------
	void RefineMesh(const int refinetime = 1,bool i_x_free = false);//2021/7/29更新 计算整体正则坐标 参数代表调用修正算法次数
	//网格参数
	PointSet PUPP;				//CST上表面节点参数
	PointSet PLOW;				//CST下表面节点参数
private:
	
	//------------------网格修正新方法--------------------
	void GetLengthMesh(string uporlow = "upp", bool ifpow = false, bool ifrefine = false);//2021/7/19更新 计算全部坐标都自由的正则坐标
	void GetLengthMesh_z(string uporlow = "upp", bool ifpow = false);//2021/7/19更新 计算z方向正则坐标
	void GetLengthMash_x(bool ifpow = false);//2021/8/5更新 修正x方向（轴向）坐标
	
	void MeshNumCorrection();//网格数修正
	mat MakeS(int PriFuncType, int K, mat B, mat fai, mat eta);//生成形状函数S
	PointSet SketchCST(mat fai, mat eta,ClassFunc shape_func, mat B, mat D, double Ratio,mat&N1m,mat&N2m);//CST基本数学模型
	PointSet SketchCST2(mat fai, mat eta, double N1, double N2, double N3, double N4,
		double M1, double M2, double T1, double T2, mat B, mat D, double Ratio,mat&N1m,mat&N2m);//CST基本数学模型
	void EvalGrid();		//评估网格信息 计算得出曲面面积等信息	
	void GetNet(bool ifrefresh = false);			//生成默认网格节点列表：FaiU、FaiL、Eta	
	//void MeshPlane() {}	//调用GMesh生成前后面网格

	mat BonePoint;			//MakeBone计算出的返回值BonePoint
	
};