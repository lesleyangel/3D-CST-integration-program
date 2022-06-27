#pragma once
#include"include/armadillo"
#include "factorial.h"
#include"material.h"
#include"CFDInfo.h"
#include"myNastran.h"
using namespace std;
using namespace arma;
using namespace arith;

//�������ṹ
class SingleBone
{
public:
	SingleBone() { BoneConnP_ID = 0; MainConnP_ID = 0; }

	mat BoneP;						//���Ǽܽڵ�
	mat BonePForce;					//�Ǽ����ڵ����
	mat GetConnBone(mat MainBone);	//���ɲ��������������ӵķ�֧�Ǽ���
	mat BoneE(int fstID = 0);		//���Ǽܵ�Ԫ(fstID=��һ��������id)
	mat StructBoneP;				//��¼�����ӽṹ������֧�ṹ at GetConnBone
	mat StructBoneE;				//at BoneE
	mat StructBoneF;				//at BoneE
	mat ConnBoneP;

	int BoneConnP_ID;				//������Ҫ���������ӵĵ���
	int MainConnP_ID;				//���Ǽ���Ҫ�������ӵĵ���
};

//�������ṹ
class Bone
{
public:
	SingleBone MainBone;			//����Ǽ���
	vector<SingleBone> BrachBone;	//��β��

	mat Bone_Pf;					//�ڵ�������
	mat Bone_P;
	mat Bone_E;
	int CalcBone();					//���������������ڵ� ��Ԫ �ڵ�������
	int SaveAsTcp(string DATname);	//���Ǽܽṹ��ӡΪ�ı���ʽ
	int SaveAsAbaqus(string DATname);//���Ǽܽṹ��ӡΪabaqus��ʽ
};

//-----������ϵ�Ԫ����--------
//����������ϲ���
class ShellandBeam
{
public:
	string name;	//���ǽṹ����
	mat node;		//�ڵ�
	mat node_f;		//�ڵ���
	mat elem_aero;	//�����ǵ�Ԫ
	mat elem_strc;	//�ṹ����Ԫ

	mat beamID;		//��¼���ṹ
	mat elem_h;		//�����ǵ�Ԫ��Ӧ�ĺ��
	mat Xbar;		//x��������λ��
	mat Zbar;		//z��������λ��

	mat LinkPtId;	//���ӵ�ı��,�ڵ���ʱ��Ҫ����
	mat ConnPtId;	//����GetConnPtId�󷽿�ʹ��
	mat GetConnPtId(const mat& MainBone);

	int write_abaqus(ofstream& ofs);
private:

};

//����������Ͻṹ
class Struct1D
{
public:
	SingleBone MainBone;			//����Ǽ���
	vector<ShellandBeam> BarchStruc;//��֧���ǽṹ

	void getConnPoint();
	
	int SaveAsAbaqus(string DATname);//���Ǽܽṹ��ӡΪabaqus��ʽ
	int SaveAsTecplt(string DatName);//���������ΪTecplot��ʽ
	
	int SaveAsNastran(string DatName, int SOL = 101);//���ΪNastran��ʽ Ĭ�ϸ�ʽΪ����ѧ����
};

