#pragma once
#include "Bone.h"
class myNastran
{
public:
	myNastran() 
	{ 
		NasPath.clear(); 
		BDFPath.clear(); 
		ResPath.clear(); 
		KillNastranBatPath.clear(); 
		logPath.clear(); 
		xdbPath.clear(); 
		pchPath.clear();
		Displacements.clear();
		ElemStress.clear();
		MaxElemStress = pair<int, double>();
	}
	string NasPath;
	string BDFPath;
	string ResPath;
	string KillNastranBatPath;
	void CalcFilePath();
	int NastranCalc(); //调用nastran计算有限元模型
	int ReadResPCH(bool isOnlyReadDisp);
	map<int, Point> GetDisplacements() { return Displacements; }
	map<int, double> GetElemStress() { return ElemStress; }
	pair<int, Point> GetMaxDisp() { return MaxNodeDisp; }
	pair<int, double> GetMaxElemStress() { return MaxElemStress; }
	pair<int, double> GetMaxElemStrain() { return MaxElemStrain; }
	//
	int ReadResPCH(vector<int> CQUAD4id, vector<int>CBARid);	//读取有限元模型中的pch文件
	vector<double> GetMaxStress() { return maxstress; }
	vector<int> GetMaxStressID() { return maxstressID; }
	vector<double> GetMaxStrain() { return maxstrain; }
	vector<int> GetMaxStrainID() { return maxstrainID; }
	vector<double> GetMaxStressBar() { return maxstress2; }
	vector<int> GetMaxStressIDBar() { return maxstressID2; }
	vector<double> GetMaxStrainBar() { return maxstrain2; }
	vector<int> GetMaxStrainIDBar() { return maxstrainID2; }
private:
	string logPath;
	string xdbPath;
	string pchPath;
	//
	map<int, Point> Displacements;//节点位移结果记录 <节点编号，位移值>
	map<int, double> ElemStress;//单元应力结果记录 <单元编号，单元应变>
	map<int, double> ElemStrain;//单元应变结果记录 <单元编号，单元应变>
	pair<int, Point>MaxNodeDisp;//最大节点位移记录（这里记录最大）
	pair<int, double> MaxElemStress;//最大单元应力记录
	pair<int, double> MaxElemStrain;//最大单元应力记录
	//
	vector<double> maxstress;//存储最大应力值
	vector<int> maxstressID;//存储最大应力单元编号
	vector<double> maxstrain;//存储最大应变值
	vector<int> maxstrainID;//存储最大应变单元编号

	vector<double> maxstress2;//存储最大应力值
	vector<int> maxstressID2;//存储最大应力单元编号
	vector<double> maxstrain2;//存储最大应变值
	vector<int> maxstrainID2;//存储最大应变单元编号
};