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
	int NastranCalc(); //����nastran��������Ԫģ��
	int ReadResPCH(bool isOnlyReadDisp);
	map<int, Point> GetDisplacements() { return Displacements; }
	map<int, double> GetElemStress() { return ElemStress; }
	pair<int, Point> GetMaxDisp() { return MaxNodeDisp; }
	pair<int, double> GetMaxElemStress() { return MaxElemStress; }
	pair<int, double> GetMaxElemStrain() { return MaxElemStrain; }
	//
	int ReadResPCH(vector<int> CQUAD4id, vector<int>CBARid);	//��ȡ����Ԫģ���е�pch�ļ�
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
	map<int, Point> Displacements;//�ڵ�λ�ƽ����¼ <�ڵ��ţ�λ��ֵ>
	map<int, double> ElemStress;//��ԪӦ�������¼ <��Ԫ��ţ���ԪӦ��>
	map<int, double> ElemStrain;//��ԪӦ������¼ <��Ԫ��ţ���ԪӦ��>
	pair<int, Point>MaxNodeDisp;//���ڵ�λ�Ƽ�¼�������¼���
	pair<int, double> MaxElemStress;//���ԪӦ����¼
	pair<int, double> MaxElemStrain;//���ԪӦ����¼
	//
	vector<double> maxstress;//�洢���Ӧ��ֵ
	vector<int> maxstressID;//�洢���Ӧ����Ԫ���
	vector<double> maxstrain;//�洢���Ӧ��ֵ
	vector<int> maxstrainID;//�洢���Ӧ�䵥Ԫ���

	vector<double> maxstress2;//�洢���Ӧ��ֵ
	vector<int> maxstressID2;//�洢���Ӧ����Ԫ���
	vector<double> maxstrain2;//�洢���Ӧ��ֵ
	vector<int> maxstrainID2;//�洢���Ӧ�䵥Ԫ���
};