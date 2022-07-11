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

//���������������
class AbstructShape
{
public:
	vector<CSTsurface> CSTsf;	//�������
	vector<GridStruct> GdStruct;//����ֵ����

	int ShapeType;				//�������ͣ�0 ͷ����1 �����ӵĻ���2 ���м����ӵĻ���ʹ�β��
	int GridRefineType;			//������ܲ���(=-1�����ܣ�=0��ȫ������=i��������)

	AbstructShape();			//Ĭ�ϼ�������Ϊ������
	void BacsShapeCSTs();		//����������������

	virtual void BuildShape() = 0;					//���麯�������ɲ���
	virtual vector<SingleBone> BuildBone() = 0;		//���麯�������ɲ������ṹ
	virtual vector<ShellandBeam> BuildStrc1D() = 0;	//���麯�������ɲ���������Ͻṹ
	void SaveToAero(string name);					//��������aero_calc��������������Ԫ�ļ�**4.23
	void MakeAeroInfo(string aeropath,int num);		//���������е���aero_calc���ļ�
	static void copyM(vector<vector<double>>& v, mat& m);	//vector��mat��������ת��

	void SaveTecplotMesh();		//���������ݵ���ΪTecplot��ʽ

private:
	void NotGridRefine();		//�����м���
	void AllGridRefine();		//�����������ʽ��м���
	void PrtGridRefine();		//���ò���part���ʽ��м���
};