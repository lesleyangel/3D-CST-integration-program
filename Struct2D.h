#pragma once
#include"include/armadillo"
#include "factorial.h"
#include"material.h"
#include"CFDInfo.h"
#include"myNastran.h"
using namespace std;
using namespace arma;
using namespace arith;

//-----ȫ�ǵ�Ԫ�ṹģ�͹���
//��������2d�ṹ��
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
	//int AeroelasticAnalysis(string fileName, Property p);//���Nastran��ʽ

	int calcAeroForce(string aeroPath);
	int calcAeroForce_AVL(string aeroPath, string exepath);
	double calcMass();
	void SaveAsNastran(string fileName, int SOL);
	static int ijk2ID(int i, int j, int k, int numi, int numj, int numk);//���ڼ�������������ֱ���id
	int updateNode(const map<int, Point&> Displacements);
	string name;	//�����ṹ��
	mat node;		//�ڵ�
	mat disp_node;
	int updateDispNode(const map<int, Point>& Displacements);
	//��Ԫ������ڵ���
	map<int, Point> node_f;//�����ڱ�����������������β�ֵ��
	map<int, Point> node_f_initial;//ֱ��������
	//�и񷨼�����
	map<int, AVLres_surface> res;
	map<int, Point> elem_p;
	//nastran������
	pair<int, double> MaxElemStress;//���Ӧ�����
	pair<int, double> MaxElemStrain;//���Ӧ����
	pair<int, Point> MaxNodeDisp;//���ڵ�λ��
	double MaxTwisting;

	mat elem_aero;	//�����ǵ�Ԫ
	mat elem_strcX;	//�ṹ�ǵ�Ԫ ��ֱ��X�����
	mat elem_strcZ;	//�ṹ�ǵ�Ԫ ��ֱ��Z�����
	//
	int getnumi() { return numi; }
	int getnumj() { return numj; }
	int getnumk() { return numk; }
private:
	int checkSite();
	void calcNode();
	void calcElem();
	void clearUniqueNode();
	int getPointID(int i, int j, int k);//�õ��ڵ�id������ظ��ڵ�֮��ģ�
	//double calcElemArea(int elemID);//���������Ԫ�ĵ�Ԫ���
	void AddOriginRotation(mat Origin, mat Rotation);//Ϊ���������ת��ƽ��
	int printAeroCalcMesh(string filepath);
	int subduction_node_f(map<int, Point>& nf);
	int calcElemPress_AVL(map<int, Point>& nf);
	int calcElemPress_AVL_refine(map<int, Point>& nf,const vector<vector<Field<double>>> &xfoilRes);
	bool isFixedMass;//�Ƿ�̶�����
	Property p;//���������б�
	vector<PointSet> node_2D;
	mat siteX_2D;
	mat siteZ_2D;
	int numi;
	int numj;
	int numk;
};
//����������2d�ṹ��
class Struct2D
{
public:
	Struct2D() { PartList.clear(); }
	vector<StructPart> PartList;
	//
	int SaveAsTecplot(string dataname);//���Tecplot��ʽ
	int SaveAsNastran(string fileName, int SOL);//���Nastran��ʽ
	int AeroAnalysis(string fileName, string exepath);
	int AeroelasticAnalysis(string fileName, string exepath);//���������Ե�����������
	double CalcTotalMass(string outputFilePath = "");
};