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
	Aircraft();							//��ʼ��������
	~Aircraft();						//�ṩ�����������ͷŲ���

	void AddShape(AbstructShape* Shape);	//Ϊ��������Ӳ���
	void CalculateAircraft();				//01.�������ɷ�����ģ��
	void SaveTecPlotAll(string DATAname);	//05.�������������ݵ���ΪTecPlot����Ҫ01��
	void GetBone();							//04.�õ��ṹ��ģ�ͣ���Ҫ01��
	void calcAeroForce(string path);		//����aero_calc����������
	double calcK(string savepath = "none");	//�������������������� 21/7/22���
	double getVol(string saepath = "none");	//��ȡ�����������
	Bone bone;

	//------------�ṹ���ɺ���-------------------
	void CalcStruct2D();
	void FindOptimalStruct2D();
	//void SaveStruct2DAsTecplot(string savepath);
	Struct2D struct2d;
	//
	void CalcStruct1D();
	Struct1D struct1d;
	//
	//------------���ݶ�ȡ�ӿ�-------------------
	void setEXEworkPath(string path) { exePath = path; }
	int RunFromFile(string filename);
	
private:

	//void AssemblePart();//װ�䲿��,�����ظ�����
	string Name;
	string exePath;
	string docPath;
	bool isTecPlotFile;		//���ɼ��������Tecplot�ļ�
	bool isAeroForce;		//����aerocalc����������
	bool isLiftDragRatio;	//��Ԫ���������������
	bool isVolume;			//���㼸�����ε��ݻ�
	bool isBoneInp;			//���ɹǼ�ģ�͵�Abaqus�����ļ�
	bool isBoneTec;			//���ɹǼ�ģ�͵�Tecplot�ļ�
	bool isStructInp;		//�������ǵ�Ԫģ�͵�Abaqus�����ļ�
	bool isStructBDF;		//�������ǵ�Ԫģ�͵�Nastran�����ļ�
	bool isStructTec;		//�������ǵ�Ԫģ�͵�Tecplotģ���ļ�
	bool isStruct2D_Tec;	//���ɿǵ�Ԫģ�͵�Tecplotģ���ļ�
	bool isStruct2D_Nas;	//���ɿǵ�Ԫģ�͵�Nastran�����ļ�
	bool isStruct2D_AeroStruct;	//�����ǵ�Ԫģ�͵ľ��������Բ�������
	bool isStruct2D_Struct101;	//�����ǵ�Ԫģ�͵ľ���ѧ��������
	bool isStruct2D_LDratio;	//�����ǵ�Ԫģ�͵��������Բ�������
	bool isStruct2D_TotalMass;	//�����ǵ�Ԫģ�͵Ľṹ��������������
	bool isStruct2D_FixedMass;	//���ÿǵ�Ԫģ�ͽṹ��������Լ��Ϊ�����̶�
	vector<AbstructShape*> CSTall;
	void compareKeyWords(stringstream& ss, string str, bool& keywords);
	//�ṹ�������Բ���
	Property m_property;

	//------���������㺯��
	void SaveToClac(string path);		//02.�洢�������������������㣨��Ҫ01��
	void ReadAeroFromFile(string path);	//03.�������ڵ�����Ϣ�����ڴ棨��Ҫ06��
	void CalcAero(string path);			//06.���ݵ������������� ����aero_calc��������������Ҫ02��

	int AeroStrcStruct2D(int XSnum, int ZSnum, CSTsurface& surf, double& maxStress, double& maxDisp, double& maxTwist);//Struct2D����������
};