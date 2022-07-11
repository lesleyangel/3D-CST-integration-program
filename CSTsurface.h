#pragma once
#include <iostream>
#include "include/armadillo"
#include "Bone.h"
#include "Curve.h"
//using namespace arma;
//using namespace std;
//using namespace arith;



//���񻮷ֵ��趨
class MeshDefine
{
public:
	enum Style
	{
		normal,	 //�����
		oblique, //��б��
	};
	int NFaiU;				//�ϱ������������
	int NFaiL;				//�±������������
	int NEta;				//����������
	int NHeight;			//�߶ȷ���������
	mat MeshRefineRatio;	//�����ܶ���������
};

//����cst������
class CSTsurface
{
public:
	//------------------��״����--------------------
	mat Origin;				//��ʼλ�ã�Xƫ�ƣ�Yƫ�ƣ�Zƫ�ƣ�
	mat Rotation;			//��ת���ӣ�X �� ��Y �� ��Z �� ��
	mat Slope;				//��б���ӣ�X����Y����Z����
	mat Scale;				//��������(X����,YUPP����,YLOW����,Z����)
	mat Length;				//��������(X����,YUPP����,YLOW����,Z����)
	ClassFunc class_upp;
	ClassFunc class_low;
	mat NS;					//ͷ��������״���ӣ��ϱ���N1/N2,�±���N1/N2��
	mat NE;					//β��������״���ӣ��ϱ���N1/N2,�±���N1/N2��
	mat M;					//���浼����״���ӣ�����M1/M2,�ҵ���M1/M2��
	mat T;					//���浼����״���ӣ��ϵ���T1/T2,�µ���T1/T2��
	double Ratio;			//��״���ӹ���ϵ��
	//------------------�������--------------------
	int NFaiU;				//�ϱ������������
	int NFaiL;				//�±������������
	int NEta;				//����������
	int NHeight;			//�߶ȷ���������

	bool Is_MeshFront;		//�Ƿ���ͷ��
	bool Is_MeshBack;		//�Ƿ���β��
	mat BUPP;				//�ϱ����������
	mat BLOW;				//�±����������
	mat DUPP;				//��������������
	mat DLOW;				//��������������
	int PriFuncType;		//���������ͣ�=0 Bernstein��������=i B������������

	mat MeshRefineRatio;	//�����ܶ���������
	//------------------����Ԫ����--------------------
	Grid GridUpp;	Grid GridLow;	MeshInfo Info;
	mat N1mUPP;//���ɽ������ �ϱ�����
	mat N2mUPP;//���ɽ������ �ϱ�����
	mat N1mLOW;//���ɽ������ �±�����
	mat N2mLOW;//���ɽ������ �±�����
	mat FaiU;	mat FaiL;		mat EtaU;	mat EtaL;mat Eta;
	//------------------��������--------------------
	CSTsurface();			//��ʼ������
	void CST3D();			//��������CST����������Ϣ	
	void EvalCurvature();	//�������� ����ó���������ֵ	
	void GridRefine();		//��������������������������
	//------------------һά�ṹ--------------------
	mat MakeBone(double bias1 = 0,double bias2 = 0);//(ƫ����[-1,1])���ɸ���������ṹ��*****4.12����
	mat ClacBonePForce();	//����Ǽ����ڵ��ϵ���������
	mat ClacBoneS();		//����Ǽ����ڵ��Ӧ�Ŀ����(��ûд)
	//------------------��ά�ṹ--------------------
	int MakeStruct2D(mat X = mat(0, 0), mat Z = mat(0,0));//���ɶ�ά�ṹ��Ԫ
	Grid StructPE;
	mat siteZend;//��������ĩ������һ����������ӽ����������
	//------------------��ά�ṹ--------------------
	//�µ����ɷ�ʽ����Ҫ����ṹ����ͬ����������������������
	int Struct2dXnum, Struct2dZnum;
	int MakeStruct2D_byNum();//���ɶ�ά�ṹ��Ԫ
	
	vector<PointSet> node_2D;//�洢�ṹ�ǵ�Ԫ����Ҫ�Ľڵ㣬�����vector�����һ��ڵ�

	
	//------------------�����&��ǿ���ṹ--------------------
	int MakeStruct1D(mat X = mat(0, 0), mat Z = mat(0, 0));//
	//mat point1D;
	mat elem_aero;
	mat elem_strc;
	ShellandBeam aero_strc;
	//------------------���������·���--------------------
	void RefineMesh(const int refinetime = 1,bool i_x_free = false);//2021/7/29���� ���������������� ����������������㷨����
	//�������
	PointSet PUPP;				//CST�ϱ���ڵ����
	PointSet PLOW;				//CST�±���ڵ����
private:
	
	//------------------���������·���--------------------
	void GetLengthMesh(string uporlow = "upp", bool ifpow = false, bool ifrefine = false);//2021/7/19���� ����ȫ�����궼���ɵ���������
	void GetLengthMesh_z(string uporlow = "upp", bool ifpow = false);//2021/7/19���� ����z������������
	void GetLengthMash_x(bool ifpow = false);//2021/8/5���� ����x������������
	
	void MeshNumCorrection();//����������
	mat MakeS(int PriFuncType, int K, mat B, mat fai, mat eta);//������״����S
	PointSet SketchCST(mat fai, mat eta,ClassFunc shape_func, mat B, mat D, double Ratio,mat&N1m,mat&N2m);//CST������ѧģ��
	PointSet SketchCST2(mat fai, mat eta, double N1, double N2, double N3, double N4,
		double M1, double M2, double T1, double T2, mat B, mat D, double Ratio,mat&N1m,mat&N2m);//CST������ѧģ��
	void EvalGrid();		//����������Ϣ ����ó������������Ϣ	
	void GetNet(bool ifrefresh = false);			//����Ĭ������ڵ��б�FaiU��FaiL��Eta	
	//void MeshPlane() {}	//����GMesh����ǰ��������

	mat BonePoint;			//MakeBone������ķ���ֵBonePoint
	
};