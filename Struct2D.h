#pragma once
#include"include/armadillo"
#include "factorial.h"
#include"material.h"
#include"CFDInfo.h"
#include"myNastran.h"
using namespace std;
using namespace arma;
using namespace arith;

//-----全壳单元结构模型构建
//单个部件2d结构类
class StructPart
{
public:
	StructPart(vector<PointSet> node_2D);
	void setSite(mat siteX, mat siteZ) { this->siteX_2D = siteX;	this->siteZ_2D = siteZ; }
	void setProperty(Property& p) { this->p = p; }
	void setIsFixedMass(bool isFixedMass) { this->isFixedMass = isFixedMass; }
	int calcStructPart(mat Origin = mat(0, 0), mat Rotation = mat(0, 0));
	void SaveAsTecplot(ofstream& ofs);
	void printAnalysisRes(string respath, bool isERROR = false);
	//int AeroelasticAnalysis(string fileName, Property p);//输出Nastran格式

	int calcAeroForce(string aeroPath);
	int calcAeroForce_AVL(string aeroPath, string exepath);
	double calcMass();
	void SaveAsNastran(string fileName, int SOL);
	static int ijk2ID(int i, int j, int k, int numi, int numj, int numk);//用于计算三阶张量拉直后的id
	int updateNode(const map<int, Point&> Displacements);
	string name;	//部件结构名
	mat node;		//节点
	mat disp_node;
	int updateDispNode(const map<int, Point>& Displacements);
	//面元法计算节点力
	map<int, Point> node_f;//作用在表面的力（静气弹两次差值）
	map<int, Point> node_f_initial;//直接作用力
	//涡格法计算结果
	map<int, AVLres_surface> res;
	map<int, Point> elem_p;
	//nastran计算结果
	pair<int, double> MaxElemStress;//最大应力结果
	pair<int, double> MaxElemStrain;//最大应变结果
	pair<int, Point> MaxNodeDisp;//最大节点位移
	double MaxTwisting;

	mat elem_aero;	//气动壳单元
	mat elem_strcX;	//结构壳单元 垂直于X方向的
	mat elem_strcZ;	//结构壳单元 垂直于Z方向的
	//
	int getnumi() { return numi; }
	int getnumj() { return numj; }
	int getnumk() { return numk; }
private:
	int checkSite();
	void calcNode();
	void calcElem();
	void clearUniqueNode();
	int getPointID(int i, int j, int k);//得到节点id（清除重复节点之后的）
	//double calcElemArea(int elemID);//获得气动单元的单元面积
	void AddOriginRotation(mat Origin, mat Rotation);//为部件添加旋转和平移
	int printAeroCalcMesh(string filepath);
	int subduction_node_f(map<int, Point>& nf);
	int calcElemPress_AVL(map<int, Point>& nf);
	int calcElemPress_AVL_refine(map<int, Point>& nf,const vector<vector<Field<double>>> &xfoilRes);
	bool isFixedMass;//是否固定质量
	Property p;//材料属性列表
	vector<PointSet> node_2D;
	mat siteX_2D;
	mat siteZ_2D;
	int numi;
	int numj;
	int numk;
};
//整个飞行器2d结构类
class Struct2D
{
public:
	Struct2D() { PartList.clear(); }
	vector<StructPart> PartList;
	//
	int SaveAsTecplot(string dataname);//输出Tecplot格式
	int SaveAsNastran(string fileName, int SOL);//输出Nastran格式
	int AeroAnalysis(string fileName, string exepath);
	int AeroelasticAnalysis(string fileName, string exepath);//静气动弹性迭代分析计算
	double CalcTotalMass(string outputFilePath = "");
};