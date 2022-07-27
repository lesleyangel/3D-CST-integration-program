#pragma once
#include"CST_Instantiation.h"
// #include<boost/any.hpp>
#include"Bone.h"
#include"Struct2D.h"
#include"material.h"

// using namespace boost;
using namespace std;


class Aircraft
{
public:
	Aircraft();							//初始化飞行器
	~Aircraft();						//提供析构函数，释放部件

	void AddShape(AbstructShape* Shape);	//为飞行器添加部件
	void CalculateAircraft();				//01.计算生成飞行器模型
	void SaveTecPlotAll(string DATAname);	//05.将气动网格数据导出为TecPlot（需要01）
	void GetBone();							//04.得到结构梁模型（需要01）
	void calcAeroForce(string path);		//调用aero_calc计算气动力
	double calcK(string savepath = "none");	//计算飞行器的整体升阻比 21/7/22添加
	double getVol(string saepath = "none");	//获取飞行器的体积
	Bone bone;

	//------------结构生成函数-------------------
	void CalcStruct2D();
	void FindOptimalStruct2D();
	//void SaveStruct2DAsTecplot(string savepath);
	Struct2D struct2d;
	//
	void CalcStruct1D();
	Struct1D struct1d;
	//
	//------------数据读取接口-------------------
	void setEXEworkPath(string path) { exePath = path; }
	int RunFromFile(string filename);
	
private:

	//void AssemblePart();//装配部件,清理重复网格
	string Name;
	string exePath;
	string docPath;
	bool isTecPlotFile;		//生成几何网格的Tecplot文件
	bool isAeroForce;		//调用aerocalc计算气动力
	bool isLiftDragRatio;	//面元法计算气动升阻比
	bool isVolume;			//计算几何外形的容积
	bool isBoneInp;			//生成骨架模型的Abaqus分析文件
	bool isBoneTec;			//生成骨架模型的Tecplot文件
	bool isStructInp;		//生成梁壳单元模型的Abaqus分析文件
	bool isStructBDF;		//生成梁壳单元模型的Nastran分析文件
	bool isStructTec;		//生成梁壳单元模型的Tecplot模型文件
	bool isStruct2D_Tec;	//生成壳单元模型的Tecplot模型文件
	bool isStruct2D_Nas;	//生成壳单元模型的Nastran分析文件
	bool isStruct2D_AeroStruct;	//分析壳单元模型的静气动弹性并输出结果
	bool isStruct2D_Struct101;	//分析壳单元模型的静力学并输出结果
	bool isStruct2D_LDratio;	//分析壳单元模型的气动特性并输出结果
	bool isStruct2D_TotalMass;	//分析壳单元模型的结构总质量并输出结果
	bool isStruct2D_FixedMass;	//设置壳单元模型结构截面属性约束为质量固定
	vector<AbstructShape*> CSTall;
	void compareKeyWords(stringstream& ss, string str, bool& keywords);
	//结构材料属性参数
	Property m_property;

	//------气动力计算函数
	void SaveToClac(string path);		//02.存储外形网格用于气动计算（需要01）
	void ReadAeroFromFile(string path);	//03.把气动节点力信息存入内存（需要06）
	void CalcAero(string path);			//06.根据导出的气动网格 调用aero_calc计算气动力（需要02）

	int AeroStrcStruct2D(int XSnum, int ZSnum, CSTsurface& surf, double& maxStress, double& maxDisp, double& maxTwist);//Struct2D静气弹分析
};