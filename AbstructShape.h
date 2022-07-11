#pragma once
#include<iostream>
#include<vector>
#include"include/armadillo"
#include"factorial.h"
#include"CSTsurface.h"
#include"Bone.h"
using namespace arma;
using namespace std;



class GridStruct
{
public:
	GridStruct() { NFaiU = 40; NFaiL = 40; NEta = 20; }
	Grid GridUpp;	Grid GridLow;
	int NFaiU;	int NFaiL;	int NEta;
	MeshInfo Info;
	Grid GridFront;	Grid GridBack;
};

//抽象飞行器部件类
class AbstructShape
{
public:
	vector<CSTsurface> CSTsf;	//曲面参数
	vector<GridStruct> GdStruct;//返回值参数

	int ShapeType;				//曲面类型（0 头身部；1 有连接的机翼；2 无中间连接的机翼和垂尾）
	int GridRefineType;			//网格加密参照(=-1不加密，=0完全修正，=i部分修正)

	AbstructShape();			//默认加密类型为不加密
	void BacsShapeCSTs();		//多体连接网格生成

	virtual void BuildShape() = 0;					//纯虚函数，生成部件
	virtual vector<SingleBone> BuildBone() = 0;		//纯虚函数，生成部件梁结构
	virtual vector<ShellandBeam> BuildStrc1D() = 0;	//纯虚函数，生成部件壳梁结合结构
	void SaveToAero(string name);					//生成用于aero_calc计算气动的网格单元文件**4.23
	void MakeAeroInfo(string aeropath,int num);		//生成命令行调用aero_calc的文件
	static void copyM(vector<vector<double>>& v, mat& m);	//vector向mat类型数据转换

	void SaveTecplotMesh();		//将计算数据导出为Tecplot形式

private:
	void NotGridRefine();		//不进行加密
	void AllGridRefine();		//采用整体曲率进行加密
	void PrtGridRefine();		//采用部分part曲率进行加密
};