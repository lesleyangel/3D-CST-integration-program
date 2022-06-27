#include <iomanip>
#include <direct.h>
#include "Bone.h"
#include "material.h"
#include "CSTsurface.h"
#include "myNastran.h"
#include "interpolation.h"
#define pi 3.1415926


//-------------------�������ṹ-------------------
// 
//���ɲ��������������ӵķ�֧�Ǽ���
mat SingleBone::GetConnBone(mat MainBone)
{
	cout << "������������������" << endl;

	//int maxID = 0;
	double Len0 = sqrt(
		pow(MainBone(0, 0) - BoneP(BoneConnP_ID, 0), 2) +
		pow(MainBone(0, 1) - BoneP(BoneConnP_ID, 1), 2) +
		pow(MainBone(0, 2) - BoneP(BoneConnP_ID, 2), 2));
	double minLen = Len0;
	//������������Ǽ�����Ѱ�ҵ�ǰ���������ϵ����ӵ��֮��ľ��루��̾��룩
	for (size_t i = 1; i < MainBone.n_rows; i++)
	{
		double len = sqrt(
			pow(MainBone(i, 0) - BoneP(BoneConnP_ID, 0), 2) +
			pow(MainBone(i, 1) - BoneP(BoneConnP_ID, 1), 2) +
			pow(MainBone(i, 2) - BoneP(BoneConnP_ID, 2), 2));
		if (minLen > len)
		{
			MainConnP_ID = i;
			minLen = len;
		}
	}
	double delta = sqrt(
		pow(BoneP(BoneConnP_ID + 1, 0) - BoneP(BoneConnP_ID, 0), 2) +
		pow(BoneP(BoneConnP_ID + 1, 1) - BoneP(BoneConnP_ID, 1), 2) +
		pow(BoneP(BoneConnP_ID + 1, 2) - BoneP(BoneConnP_ID, 2), 2));
	//����������Ϣ���������
	int ConnENum = (int)(minLen / delta);//�������ĵ�Ԫ����(�Ƚڵ������1)
	if (ConnENum < 1) ConnENum = 1;
	//ConnENum = 10;
	ConnBoneP = zeros(ConnENum + 1, 3);
	for (int i = 0; i < ConnENum; i++)
	{
		ConnBoneP(i, 0) = MainBone(MainConnP_ID, 0) + (BoneP(BoneConnP_ID, 0) - MainBone(MainConnP_ID, 0)) * i / ConnENum;
		ConnBoneP(i, 1) = MainBone(MainConnP_ID, 1) + (BoneP(BoneConnP_ID, 1) - MainBone(MainConnP_ID, 1)) * i / ConnENum;
		ConnBoneP(i, 2) = MainBone(MainConnP_ID, 2) + (BoneP(BoneConnP_ID, 2) - MainBone(MainConnP_ID, 2)) * i / ConnENum;
	}
	ConnBoneP(ConnENum, 0) = BoneP(BoneConnP_ID, 0);
	ConnBoneP(ConnENum, 1) = BoneP(BoneConnP_ID, 1);
	ConnBoneP(ConnENum, 2) = BoneP(BoneConnP_ID, 2);

	mat connBoneF = zeros(ConnBoneP.n_rows, 6);

	if (BoneConnP_ID != 0)
	{
		StructBoneP = join_cols(ConnBoneP.rows(1, ConnBoneP.n_rows - 1), BoneP.rows(0, BoneConnP_ID - 1));//���ɴ������ӽṹ��������
		StructBoneP = join_cols(StructBoneP, BoneP.rows(BoneConnP_ID + 1, BoneP.n_rows - 1));

		
		if (BonePForce.n_rows!=0)
		{
			connBoneF.row(connBoneF.n_rows - 1) = BonePForce.row(BoneConnP_ID);
			StructBoneF = join_cols(connBoneF.rows(1, connBoneF.n_rows - 1), BonePForce.rows(0, BoneConnP_ID - 1));//��װ���ڵ����
			StructBoneF = join_cols(StructBoneF, BonePForce.rows(BoneConnP_ID + 1, BonePForce.n_rows - 1));
		}
	}
	else if (ConnBoneP.n_rows < 3)
	{
		StructBoneP = join_cols(ConnBoneP.row(1), BoneP.rows(BoneConnP_ID + 0, BoneP.n_rows - 1));
		if (BonePForce.n_rows != 0)
		{
			StructBoneF = join_cols(connBoneF.row(1), BonePForce.rows(BoneConnP_ID + 0, BoneP.n_rows - 1));//��װ���ڵ����
		}
	}
	else
	{
		StructBoneP = join_cols(ConnBoneP.rows(1, ConnBoneP.n_rows - 2), BoneP.rows(BoneConnP_ID + 0, BoneP.n_rows - 1));
		if (BonePForce.n_rows != 0)
		{
			StructBoneF = join_cols(connBoneF.rows(1, connBoneF.n_rows - 2), BonePForce.rows(BoneConnP_ID + 0, BoneP.n_rows - 1));//��װ���ڵ����
		}
	}

	return ConnBoneP;
}

//���Ǽܵ�Ԫ(fstID=��һ��������id
mat SingleBone::BoneE(int fstID/* = 0*/)
{
	if (ConnBoneP.n_rows > 1)//���ڸ�֦
	{
		//���ӽṹ����Ԫ
		mat connBoneE = zeros(ConnBoneP.n_rows - 1, 2);
		for (size_t i = 0; i < ConnBoneP.n_rows - 1; i++)
		{
			connBoneE(i, 0) = i;
			connBoneE(i, 1) = i + 1.0;
		}
		//connBoneE(0, 0) = MainConnP_ID;

		//ԭ���ṹ����Ԫ
		mat brchBoneE = zeros(BoneP.n_rows - 1, 2);

		for (size_t i = 0; i < BoneP.n_rows - 1; i++)
		{
			brchBoneE(i, 0) = i;
			brchBoneE(i, 1) = i + 1.0;
		}
		for (size_t i = 0; i < brchBoneE.n_rows; i++)
		{
			if (brchBoneE(i, 0) == BoneConnP_ID) brchBoneE(i, 0) = -1;//���ӵ㴦ID��Ϊ-1
			if (brchBoneE(i, 1) == BoneConnP_ID) brchBoneE(i, 1) = -1;
			if (brchBoneE(i, 0) > BoneConnP_ID) brchBoneE(i, 0) -= 1;//���ӵ�֮���ID��1
			if (brchBoneE(i, 1) > BoneConnP_ID) brchBoneE(i, 1) -= 1;
		}
		brchBoneE += connBoneE.n_rows + 1.0;//�ṹ�ڵ��� ����ID ���Ų��ô��
		StructBoneE = join_cols(connBoneE, brchBoneE) + fstID - 1;
		StructBoneE(0, 0) = MainConnP_ID;
		return StructBoneE;
	}
	else//����û�����Ӳ���
	{
		if (BoneP.n_rows == 0)//��ֹû�в�������������Ϊ-1
		{
			return zeros(0, 2);
		}
		mat brchBoneE = zeros(BoneP.n_rows - 1, 2);

		for (size_t i = 0; i < BoneP.n_rows - 1; i++)
		{
			brchBoneE(i, 0) = i;
			brchBoneE(i, 1) = i + 1.0;
		}
		return brchBoneE;
	}
}


//-------------------�������ṹ-------------------
// 
//���������������ڵ� ��Ԫ �ڵ�������
int Bone::CalcBone()
{
	Bone_P = MainBone.BoneP;
	Bone_E = MainBone.BoneE();
	if (Bone_E.n_rows == 0)
	{
		cout << "���������ڣ�����Ҫ�����������ӣ���" << endl;
		return 0;
	}
	Bone_Pf = MainBone.BonePForce;

	int fstID = MainBone.BoneP.n_rows;
	for (size_t i = 0; i < BrachBone.size(); i++)
	{
		//cout << "��" << i << "�������Ľڵ�͵�Ԫ" << endl;
		BrachBone[i].GetConnBone(MainBone.BoneP);//���������͸�֦����������
		mat strcE = BrachBone[i].BoneE(fstID);
		mat strcP = BrachBone[i].StructBoneP;

		fstID += strcE.n_rows;

		Bone_E = join_cols(Bone_E, strcE);//����֦��Ԫ�ӵ���������Ԫ��
		Bone_P = join_cols(Bone_P, strcP);
		Bone_Pf = join_cols(Bone_Pf, BrachBone[i].StructBoneF);
	}
	Bone_E += 1;
	if (Bone_P.n_rows == Bone_Pf.n_rows)
	{
		return 1;
	}
	else
	{
		cout << "�ڵ����ͽڵ�δ��ȷ��Ӧ�����ڵ��� = " << Bone_P.n_rows << " �ڵ����� = " << Bone_Pf.n_rows << "��" << endl;
		return 0;
	}
}

//���Ǽܽṹ��ӡΪ�ı���ʽ
int Bone::SaveAsTcp(string name)
{
	ofstream ofs;
	ofs.open(name, ios::trunc);

	ofs << "Title = \"Structal BONE\"" << endl;
	ofs << "variables = \"x\",\"y\",\"z\"" << endl;
	ofs << "zone N=" << Bone_P.n_rows << ",E=" << Bone_E.n_rows << ",datapacking=block" << endl;
	ofs << "zonetype = FELINESEG" << endl;
	ofs << trans(Bone_P) << endl;
	ofs << Bone_E << endl;

	ofs.close();
	cout << "-------------------" << endl;
	cout << "����TecPlot��ʽ���ṹ����ɹ���" << endl;
	cout << "-------------------" << endl;
	return 0;
}

//���Ǽܽṹ��ӡΪabaqus��ʽ
int Bone::SaveAsAbaqus(string name)
{
	ofstream ofs;
	ofs.open(name, ios::trunc);
	//--------------------------.inp�ļ�--------------------------
	ofs << "*Heading" << endl;
	ofs << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;
	ofs << "*Part, name=beam" << endl;
	//--------------------------д��ڵ�--------------------------
	ofs << "*Node" << endl;
	for (size_t i = 0; i < Bone_P.n_rows; i++)
	{
		ofs << setw(8) << i + 1
			<< "," << setw(13) << Bone_P(i, 0)
			<< "," << setw(13) << Bone_P(i, 1)
			<< "," << setw(13) << Bone_P(i, 2) << endl;
	}
	//--------------------------д�뵥Ԫ--------------------------
	ofs << "*Element, type = B31" << endl;
	for (size_t i = 0; i < Bone_E.n_rows; i++)
	{
		ofs << setw(8) << i + 1 << ","
			<< setw(8) << Bone_E(i, 0) << ","
			<< setw(8) << Bone_E(i, 1) << endl;
	}
	ofs << "*Elset, elset=_I1, internal, generate" << endl;
	ofs << setw(5) << 1 << "," << setw(5) << Bone_E.n_rows << ",    1" << endl;
	//�������
	ofs << "*Beam Section, elset=_I1, material=JIEGOUGANG, temperature=GRADIENTS, section=CIRC" << endl;
	ofs << 5 << "." << endl;//���ý���ֱ��Ϊ5
	ofs << "1., 1., 1." << endl;//���ý��淽��
	ofs << "*End Part" << endl;
	//--------------------------����װ��--------------------------
	ofs << "*Assembly, name=Assembly" << endl;
	ofs << "*Instance, name=BEAM-1, part=BEAM" << endl;//
	ofs << "*End Instance" << endl;
	//�������е㼯��
	for (size_t i = 0; i < Bone_P.n_rows; i++)
	{
		ofs << "*Nset, nset=_PickedSet" << i + 1 << ", internal, instance=BEAM-1" << endl;
		ofs << " " << i + 1 << "," << endl;
	}
	ofs << "*End Assembly" << endl;
	//--------------------------�������--------------------------
	ofs << "*Material, name = JIEGOUGANG" << endl;
	ofs << "* Density" << endl;//�ܶ�
	ofs << "7280.," << endl;
	ofs << "* Elastic" << endl;//����ģ��
	ofs << "2.1e+11, 0.28" << endl;
	//-------------------------- ������ --------------------------
	ofs << "*Step, name = Step - 1, nlgeom = NO, inc = 1000" << endl;
	ofs << "* Static" << endl;
	ofs << "0.01, 1., 1e-08, 1." << endl;
	//--------------------------�غ�����--------------------------
	for (size_t i = 0; i < Bone_P.n_rows; i++)
	{
		ofs << "*Cload" << endl;//�ڵ㼯����
		if (fabs(Bone_Pf(i, 0)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 1, " << Bone_Pf(i, 0) << endl;//force_x
		if (fabs(Bone_Pf(i, 1)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 2, " << Bone_Pf(i, 1) << endl;//force_y
		if (fabs(Bone_Pf(i, 2)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 3, " << Bone_Pf(i, 2) << endl;//force_z
		if (fabs(Bone_Pf(i, 3)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 4, " << Bone_Pf(i, 3) << endl;//force_z
		if (fabs(Bone_Pf(i, 4)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 5, " << Bone_Pf(i, 4) << endl;//force_z
		if (fabs(Bone_Pf(i, 5)) > 0)
			ofs << "_PickedSet" << i + 1 << ", 6, " << Bone_Pf(i, 5) << endl;//force_z
	}
	ofs << "*Inertia Relief" << endl;//����Ϊ�����ͷ�
	//--------------------------�������--------------------------
	ofs << "*Restart, write, frequency=0" << endl;
	ofs << "*Output, field, variable=PRESELECT" << endl;
	ofs << "*Output, history, variable=PRESELECT" << endl;
	ofs << "*End Step" << endl;
	//--------------------------д�����--------------------------

	ofs.close();
	cout << "-------------------" << endl;
	//cout << "SaveAbaqus Structal BONE Completed!" << endl;
	cout << "���ɽṹ����ԪABAQUS .inp�ļ��ɹ���" << endl;
	cout << "-------------------" << endl;
	return 0;
}


//--------------������ϵ�Ԫ����------------------
mat ShellandBeam::GetConnPtId(const mat& MainBone)
{
	ConnPtId = zeros(LinkPtId.n_rows, 1);
	for (size_t id = 0; id < LinkPtId.n_rows; id++)
	{
		int BoneConnP_ID = (int)LinkPtId(id);
		int MainConnP_ID = 0;
		//================
		double Len0 = sqrt(
			pow(MainBone(0, 0) - node(BoneConnP_ID, 0), 2) +
			pow(MainBone(0, 1) - node(BoneConnP_ID, 1), 2) +
			pow(MainBone(0, 2) - node(BoneConnP_ID, 2), 2));
		double minLen = Len0;
		//������������Ǽ�����Ѱ�ҵ�ǰ���������ϵ����ӵ��֮��ľ��루��̾��룩
		for (size_t i = 1; i < MainBone.n_rows; i++)
		{
			double len = sqrt(
				pow(MainBone(i, 0) - node(BoneConnP_ID, 0), 2) +
				pow(MainBone(i, 1) - node(BoneConnP_ID, 1), 2) +
				pow(MainBone(i, 2) - node(BoneConnP_ID, 2), 2));
			if (minLen > len)
			{
				MainConnP_ID = i;
				minLen = len;
			}
		}
		ConnPtId(id) = MainConnP_ID;
	}
	return ConnPtId;
}

int ShellandBeam::write_abaqus(ofstream& ofs)
{
	ofs << "*Part, name=" << name << endl;
	//--------------------------д��ڵ�--------------------------
	ofs << "*Node" << endl;
	for (size_t i = 0; i < node.n_rows; i++)
	{
		ofs << setw(8) << i + 1
			<< "," << setw(13) << node(i, 0)
			<< "," << setw(13) << node(i, 1)
			<< "," << setw(13) << node(i, 2) << endl;
	}

	//--------------------------д������Ԫ--------------------------
	ofs << "*Element, type = B31" << endl;
	for (size_t i = 0; i < elem_strc.n_rows; i++)
	{
		ofs << setw(8) <<  i + 1 << ","
			<< setw(8) << elem_strc(i, 0)+1 << ","
			<< setw(8) << elem_strc(i, 1)+1 << endl;
	}
	//-------------------����Elset------------------
	int begin = 0, end = 0, numID = 0;
	for (size_t i = 1; i < elem_strc.n_rows; i++)
	{
		if (elem_strc(i, 3) != elem_strc(i-1, 3))
		{
			end = i - 1;
			string elsetName = name + "_b" + to_string(int(elem_strc(i - 1, 2))) + to_string(int(elem_strc(i - 1, 3))) + to_string(numID);
			ofs << "*Elset, elset=" << elsetName <<", internal, generate" << endl;
			ofs << setw(5) << begin+1 << "," << setw(5) << end+1 << ",    1" << endl;
			//-----
			ofs << "*Beam Section, elset="<< elsetName <<", material=Material-1, temperature=GRADIENTS, section=RECT" << endl;
			ofs << "150, 150" << endl;
			ofs << "0,1,0" << endl;
			begin = i;
			numID++;
		}
	}
	string elsetName = name + "_b" + to_string(int(elem_strc(elem_strc.n_rows - 1, 2))) + to_string(int(elem_strc(elem_strc.n_rows - 1, 3))) + to_string(numID);
	ofs << "*Elset, elset=" << elsetName << ", internal, generate" << endl;
	ofs << setw(5) << begin + 1 << "," << setw(5) << elem_strc.n_rows << ",    1" << endl;
	//-----
	ofs << "*Beam Section, elset=" << elsetName << ", material=Material-1, temperature=GRADIENTS, section=RECT" << endl;
	ofs << "150, 150" << endl;
	ofs << "0,1,0" << endl;

	//--------------------------д��ǵ�Ԫ--------------------------
	ofs << "*Element, type = S4R" << endl;
	
	for (size_t i = 0; i < elem_aero.n_rows; i++)
	{
		ofs << setw(8) << elem_strc.n_rows + i + 1 << ","
			<< setw(8) << elem_aero(i, 0)+1 << ","
			<< setw(8) << elem_aero(i, 1)+1 << ","
			<< setw(8) << elem_aero(i, 2)+1 << ","
			<< setw(8) << elem_aero(i, 3)+1 << endl;
	}
	//----------------�ǵ�ԪElset----------------
	ofs << "*Elset, elset="<< name <<"_s, internal, generate" << endl;
	ofs << setw(5) << elem_strc.n_rows + 1 << "," << setw(5) << elem_strc.n_rows + elem_aero.n_rows << ",    1" << endl;
	//-----
	ofs << "*Shell Section, elset=" << name << "_s, material=Material-1" << endl;
	ofs << "5, 5" << endl;//0.1, 5


	ofs << "*End Part" << endl;
	return 0;
}


void Struct1D::getConnPoint()
{
	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		BarchStruc[i].GetConnPtId(MainBone.BoneP);
	}
}

int Struct1D::SaveAsAbaqus(string DATname)
{
	ofstream ofs;
	ofs.open(DATname, ios::trunc);

	//------------------���ⶨ��-------------------
	ofs << "*Heading" << endl;
	ofs << "*Preprint, echo=NO, model=NO, history=NO, contact=NO" << endl;

	//---
	ofs << "*Part, name=MainBone" << endl;
	ofs << "*Node" << endl;
	for (size_t i = 0; i < MainBone.BoneP.n_rows; i++)
	{
		ofs << setw(8) << i + 1
			<< "," << setw(13) << MainBone.BoneP(i, 0)
			<< "," << setw(13) << MainBone.BoneP(i, 1)
			<< "," << setw(13) << MainBone.BoneP(i, 2) << endl;
	}
	ofs << "*Element, type = B31" << endl;
	for (size_t i = 0; i < MainBone.BoneE().n_rows; i++)
	{
		ofs << setw(8) << i + 1 << ","
			<< setw(8) << MainBone.BoneE()(i, 0)+1 << ","
			<< setw(8) << MainBone.BoneE()(i, 1)+1 << endl;
	}
	ofs << "*Elset, elset=MainBoneE, internal, generate" << endl;
	ofs << setw(5) << 1 << "," << setw(5) << MainBone.BoneE().n_rows << ",    1" << endl;
	//�������
	ofs << "*Beam Section, elset=MainBoneE, material=JIEGOUGANG, temperature=GRADIENTS, section=CIRC" << endl;
	ofs << 1000 << "." << endl;//���ý���ֱ��Ϊ5
	ofs << "0., 1., 0." << endl;//���ý��淽��
	ofs << "*End Part" << endl;
	//--------------����Partģ�鶨��---------------
	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		BarchStruc[i].write_abaqus(ofs);
	}

	//--------------ASSEMBLYģ�鶨��---------------
	ofs << "*Assembly, name=Assembly" << endl;

	ofs << "*Instance, name=MainBone-1, part=MainBone" << endl;//װ��ͷ������
	ofs << "*End Instance" << endl;
	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		ofs << "*Instance, name=" << BarchStruc[i].name << "-1, part=" << BarchStruc[i].name << endl;//װ����������
		ofs << "*End Instance" << endl;
	}
	for (size_t i = 0; i < BarchStruc.size(); i++)//����Լ��
	{
		for (size_t j = 0; j < BarchStruc[i].LinkPtId.n_rows; j++)
		{
			ofs << "*Nset, nset=LinkP-" << i << "-" << j << ", internal, instance=" << BarchStruc[i].name << "-1" << endl;
			ofs << BarchStruc[i].LinkPtId(j)+1 << "," << endl;
			ofs << "*Nset, nset=ConnP-" << i << "-" << j << ", internal, instance=MAINBONE-1" << endl;
			ofs << BarchStruc[i].ConnPtId(j)+1 << "," << endl;
			ofs << "*MPC" << endl;//����MPCԼ��
			ofs << "BEAM,  LinkP-" << i << "-" << j << ", ConnP-" << i << "-" << j << endl;
		}
	}
	for (size_t i = 0; i < MainBone.BoneP.n_rows; i++)
	{
		ofs << "*Nset, nset=MB_P_" << i << ", internal, instance=MainBone-1" << endl;
		ofs << i + 1 << "," << endl;
	}
	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		for (size_t j = 0; j < BarchStruc[i].node.n_rows; j++)
		{
			ofs << "*Nset, nset=BarchP-" << i << "-" << j << ", internal, instance=" << BarchStruc[i].name << "-1" << endl;
			ofs << j + 1 << "," << endl;
		}
	}

	ofs << "*End Assembly" << endl;

	//--------------Materialģ�鶨��---------------
	ofs << "*Material, name=Material-1" << endl;
	ofs << "*Density" << endl;
	ofs << "2760.," << endl;
	ofs << "*Elastic" << endl;
	ofs << " 73e+06, 0.324" << endl;

	ofs << "*Material, name = JIEGOUGANG" << endl;
	ofs << "* Density" << endl;//�ܶ�
	ofs << "2280.," << endl;
	ofs << "* Elastic" << endl;//����ģ��
	ofs << "2.1e+7, 0.28" << endl;
	ofs << "** ----------------------------------------------------------------" << endl;
	//** ----------------------------------------------------------------
	//--------------STEP������ģ�鶨��---------------
	ofs << "*Step, name=Step-1, nlgeom=NO" << endl;
	ofs << "*Static" << endl;
	ofs << "1., 1., 1e-05, 1." << endl;

	//-------------------�غɶ���--------------------
	ofs << "*Inertia Relief" << endl;//�����ͷ�
	for (size_t i = 0; i < MainBone.BoneP.n_rows; i++)
	{
		ofs << "*Cload" << endl;
		for (int z = 0; z < 6; z++)
		{
			if (MainBone.BonePForce(i, z) != 0)
			{
				ofs << "MB_P_" << i << ", " << z + 1 << ", " << MainBone.BonePForce(i, z) << endl;
			}
		}
	}
	
	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		for (size_t j = 0; j < BarchStruc[i].node.n_rows; j++)
		{
			ofs << "*Cload" << endl;
			for (int z = 0; z < 6; z++)
			{
				if (BarchStruc[i].node_f(j, z) != 0)
				{
					ofs << "BarchP-" << i << "-" << j << ", " << z + 1 << ", " << BarchStruc[i].node_f(j, z) << endl;
				}
			}
		}
	}

	//-------------------�������--------------------
	ofs << "*Restart, write, frequency=0" << endl;
	ofs << "*Output, field, variable=PRESELECT" << endl;
	ofs << "*Output, history, variable=PRESELECT" << endl;
	ofs << "*End Step" << endl;
	ofs.close();
	return 0;
}

int Struct1D::SaveAsTecplt(string DatName)
{
	ofstream ofs, ofs1;
	ofs.open(DatName + "_strc1D_strc.dat", ios::trunc);
	ofs1.open(DatName + "_strc1D_aero.dat", ios::trunc);

	ofs << "Title = \"Structal MainBody\"" << endl;
	ofs << "variables = \"x\",\"y\",\"z\"" << endl;
	ofs << "zone N=" << MainBone.BoneP.n_rows << ",E=" << MainBone.BoneE().n_rows << ",datapacking=block" << endl;
	ofs << "zonetype = FELINESEG" << endl;
	ofs << trans(MainBone.BoneP) << endl;
	ofs << MainBone.BoneE() + 1 << endl;

	ofs1 << "Title = \"Structal MainBody\"" << endl;
	ofs1 << "variables = \"x\",\"y\",\"z\"" << endl;
	ofs1 << "zone N=" << MainBone.BoneP.n_rows << ",E=" << MainBone.BoneE().n_rows << ",datapacking=block" << endl;
	ofs1 << "zonetype = FELINESEG" << endl;
	ofs1 << trans(MainBone.BoneP) << endl;
	ofs1 << MainBone.BoneE() + 1 << endl;

	for (size_t i = 0; i < BarchStruc.size(); i++)
	{
		ofs << "Title = \"Structal BEAM"<<i<<"\"" << endl;
		ofs << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs << "zone N=" << BarchStruc[i].node.n_rows << ",E=" << BarchStruc[i].elem_strc.n_rows << ",datapacking=block" << endl;
		ofs << "zonetype = FELINESEG" << endl;
		ofs << trans(BarchStruc[i].node) << endl;
		ofs << BarchStruc[i].elem_strc.cols(0,1) + 1 << endl;
	
		ofs1 << "Title = \"Structal AERO" << i << "\"" << endl;
		ofs1 << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs1 << "zone N=" << BarchStruc[i].node.n_rows << ",E=" << BarchStruc[i].elem_aero.n_rows << ",datapacking=block" << endl;
		ofs1 << "zonetype = fequadrilateral" << endl;
		ofs1 << trans(BarchStruc[i].node) << endl;
		ofs1 << BarchStruc[i].elem_aero + 1 << endl;
	}
	ofs.close();
	ofs1.close();

	return 0;
}

int Struct1D::SaveAsNastran(string DatName, int SOL/* = 101*/)
{
	std::cout << "$----------------------------------------$" << endl;
	std::cout << "$        MSC.Nastran�����ʽ���ڴ�ӡ     $" << endl;
	std::cout << "$----------------------------------------$" << endl;

	ofstream ofs, ofs1, ofs2;
	string fileName = DatName + "_mesh.bdf";
	string fileName1 = DatName + "_property.bdf";
	string fileName2 = DatName + "_force.bdf";
	ofs.open(fileName, ios::trunc);//�ڵ������ļ�
	ofs1.open(fileName1, ios::trunc);//������ǵ�Ԫ����
	ofs2.open(fileName2, ios::trunc);//�ڵ�����Լ��

	ofs << fixed << setprecision(8);//�����������Ϊ8
	ofs1 << showpoint << setprecision(3);//ǿ�����С����
	ofs2 << fixed << setprecision(5);//�����������Ϊ5

	int ID = 0;//�ڵ���
	int EID = 0;//��Ԫ���
	int PID = 0;//����or�ǵ�Ԫ���Ա��
	int SID = 0;//�غ�orԼ���������

	//-----------------------------��Ԫ�ڵ���Ϣ----------------------------------
	ofs << "$----------------------------------------" << endl;
	ofs << "$NAME   MainBone" << endl;
	ofs << "$----------------------------------------" << endl;
	
	//----------------node----------------
	for (size_t i = 0; i < MainBone.BoneP.n_rows; i++)
	{
		ID++;
		
		ofs << "GRID*   " << setw(16) << ID << setw(16) << 0 << setw(16) <<MainBone.BoneP(i, 0)
			<< setw(16) << MainBone.BoneP(i, 1) << setw(8) << ID << endl;
		ofs << "*" << setw(7) << ID << setw(16) << MainBone.BoneP(i, 2) << endl;

		//----------------node force----------------
		if (abs(MainBone.BonePForce(i, 0)) > 1e-5 || abs(MainBone.BonePForce(i, 1)) > 1e-5 || abs(MainBone.BonePForce(i, 2)) > 1e-5)
		{
			ofs2 << "FORCE*  " << setw(16) << ++SID << setw(16) << ID << setw(16) << 0 << setw(16) << 1.0 << "*" << endl;
			ofs2 << "*       "
				<< setw(16) << MainBone.BonePForce(i, 0)
				<< setw(16) << MainBone.BonePForce(i, 1)
				<< setw(16) << MainBone.BonePForce(i, 2) << endl;
		}
		if (abs(MainBone.BonePForce(i, 3)) > 1e-5 || abs(MainBone.BonePForce(i, 4)) > 1e-5 || abs(MainBone.BonePForce(i, 5)) > 1e-5)
		{
			ofs2 << "MOMENT* " << setw(16) << ++SID << setw(16) << ID << setw(16) << 0 << setw(16) << 1.0 << "*" << endl;
			ofs2 << "*       "
				<< setw(16) << MainBone.BonePForce(i, 3)
				<< setw(16) << MainBone.BonePForce(i, 4)
				<< setw(16) << MainBone.BonePForce(i, 5) << endl;
		}
	}
	//----------------property----------------
	//PID++;
	PBARL mb;
	mb.name = "MainBone";
	mb.PID = ++PID;
	mb.TYPE = "TUBE";//��������-Բ��
	mb.DIM.push_back(100);//�⾶
	mb.DIM.push_back(50);//�ھ�
	mb.printPBARL(ofs1);
	MAT1 mt; mt.MID = PID;
	mt.printMAT1(ofs1);
	//----------------element----------------
	mat MainElem = MainBone.BoneE();
	ofs << setprecision(2);//�����������Ϊ2
	for (size_t i = 0; i < MainElem.n_rows; i++)
	{
		//EID++;
		mat orntBar = MainBone.BoneP.row((size_t)MainElem(i, 1)) - MainBone.BoneP.row((size_t)MainElem(i, 0));
		mat orntX = {-orntBar(2), 0, orntBar(0)};//���ﶨ����y���˵ķ���Ϊ����ָ�ɷ���
		ofs << "CBAR    " << setw(8) << ++EID << setw(8) << PID
			<< setw(8) << (int)MainElem(i, 0) + 1
			<< setw(8) << (int)MainElem(i, 1) + 1
			<< setw(8) << orntX(0) 
			<< setw(8) << orntX(1) 
			<< setw(8) << orntX(2) << endl;
	}

	//ÿ�����������
	for (size_t n = 0; n < BarchStruc.size(); n++)
	{
		ofs << "$----------------------------------------" << endl;
		ofs << "$NAME   " << BarchStruc[n].name << "-" << n << endl;
		ofs << "$----------------------------------------" << endl;
		ofs1 << "$----------------------------------------" << endl;
		ofs1 << "$NAME   " << BarchStruc[n].name << "-" << n << endl;
		ofs1 << "$----------------------------------------" << endl;
		
		int IDplus = ID+1;
		//----------------node----------------
		ofs << setprecision(8);//�����������Ϊ8
		for (size_t i = 0; i < BarchStruc[n].node.n_rows; i++)
		{
			ID++;
			ofs << "GRID*   " << setw(16) << ID << setw(16) << 0 << setw(16) << BarchStruc[n].node(i, 0)
				<< setw(16) << BarchStruc[n].node(i, 1) << setw(8) << ID << endl;
			ofs << "*" << setw(7) << ID << setw(16) << BarchStruc[n].node(i, 2) << endl;

			//----------------node force----------------
			if (abs(BarchStruc[n].node_f(i, 0)) > 1e5 || abs(BarchStruc[n].node_f(i, 1)) > 1e5 || abs(BarchStruc[n].node_f(i, 2)) > 1e5)
			{
				ofs2 << "FORCE*  " << setw(16) << ++SID << setw(16) << ID << setw(16) << 0 << setw(16) << 1.0 << "*" << endl;
				ofs2 << "*       "
					<< setw(16) << BarchStruc[n].node_f(i, 0)
					<< setw(16) << BarchStruc[n].node_f(i, 1)
					<< setw(16) << BarchStruc[n].node_f(i, 2) << endl;
			}
			if (abs(BarchStruc[n].node_f(i, 3)) > 1e5 || abs(BarchStruc[n].node_f(i, 4)) > 1e5 || abs(BarchStruc[n].node_f(i, 5)) > 1e5)
			{
				ofs2 << "MOMENT* " << setw(16) << ++SID << setw(16) << ID << setw(16) << 0 << setw(16) << 1.0 << "*" << endl;
				ofs2 << "*       "
					<< setw(16) << BarchStruc[n].node_f(i, 3)
					<< setw(16) << BarchStruc[n].node_f(i, 4)
					<< setw(16) << BarchStruc[n].node_f(i, 5) << endl;
			}
			
			
		}
		//----------------constraint----------------
		ofs << "$----------------------------------------$" << endl;
		ofs << "$               MPC Setting              $" << endl;
		ofs << "$----------------------------------------$" << endl;
		for (size_t i = 0; i < BarchStruc[n].LinkPtId.n_rows; i++)
		{
			//RBE2
			ofs << "RBE2    " << setw(8) << ++EID
				<< setw(8) << (int)BarchStruc[n].ConnPtId(i) + 1 << setw(8) << "123456"
				<< setw(8) << (int)BarchStruc[n].LinkPtId(i) + IDplus << endl;
			
			//MPC ����ʧ��
			//ofs2 << "MPC     " << setw(8) << ++SID
			//	<< setw(8) << BarchStruc[n].ConnPtId(i) + 1 << setw(8) << 123456 << setw(8) << 1
			//	<< setw(8) << BarchStruc[n].LinkPtId(i) + IDplus << setw(8) << 123456 << setw(8) << -1 << endl;
		}

		/*PID++;*/
		for (size_t i = 0; i < BarchStruc[n].elem_strc.n_rows; i++)
		{
			if (i == 0)
			{
				//----------------property----------------
				PBARL p;
				p.PID = ++PID;
				p.name = BarchStruc[n].name + to_string(PID);
				p.TYPE = "I";//��������-���ָ�
				p.DIM.resize(6);//�⾶
				p.DIM[0] = 10; p.DIM[1] = 8; p.DIM[2] = 8;
				p.DIM[3] = 3 ; p.DIM[4] = 3; p.DIM[5] = 3;
				p.printPBARL(ofs1);
				MAT1 m; m.MID = PID;
				m.printMAT1(ofs1);
				
			}
			else
			{ 
				if (BarchStruc[n].elem_strc(i, 3) != BarchStruc[n].elem_strc(i - 1, 3))
				{
					//----------------property----------------
					PBARL p;
					p.PID = ++PID;
					p.name = BarchStruc[n].name + to_string(PID);
					p.TYPE = "I";//��������-���ָ�
					p.DIM.resize(6);//�⾶
					p.DIM[0] = 10; p.DIM[1] = 8; p.DIM[2] = 8;
					p.DIM[3] = 3; p.DIM[4] = 3; p.DIM[5] = 3;
					p.printPBARL(ofs1);
					MAT1 m; m.MID = PID;
					m.printMAT1(ofs1);
				}
			}
			//----------------element----------------
			
			ofs << setprecision(2);//�����������Ϊ2
			mat orntBar = BarchStruc[n].node.row((size_t)BarchStruc[n].elem_strc(i, 1)) - BarchStruc[n].node.row((size_t)BarchStruc[n].elem_strc(i, 0));
			mat orntX = { -orntBar(2), 0, orntBar(0) };//���ﶨ����y���˵ķ���Ϊ����ָ�ɷ���
			ofs << "CBAR    " << setw(8) << ++EID << setw(8) << PID
				<< setw(8) << (int)BarchStruc[n].elem_strc(i, 0) + IDplus
				<< setw(8) << (int)BarchStruc[n].elem_strc(i, 1) + IDplus
				<< setw(8) << orntX(0)
				<< setw(8) << orntX(1)
				<< setw(8) << orntX(2) << endl;
		}

		//----------------property----------------
		PSHELL ps;
		ps.PID = ++PID;
		ps.name = BarchStruc[n].name + to_string(PID);
		ps.T = 0.5;
		ps.printPSHELL(ofs1);
		MAT1 m; m.MID = PID;
		m.printMAT1(ofs1);

		for (size_t i = 0; i < BarchStruc[n].elem_aero.n_rows; i++)//�������������Ԫ
		{
			//EID++;
			ofs << "CQUAD4  " << setw(8) << ++EID << setw(8) << PID 
				<< setw(8) << (int)BarchStruc[n].elem_aero(i, 0) + IDplus
				<< setw(8) << (int)BarchStruc[n].elem_aero(i, 1) + IDplus
				<< setw(8) << (int)BarchStruc[n].elem_aero(i, 2) + IDplus
				<< setw(8) << (int)BarchStruc[n].elem_aero(i, 3) + IDplus << endl;
		}
	}


	//--------------------------------------------------
	ofstream ofs3;
	switch (SOL)
	{
	case 101://��������
		ofs3.open(DatName + "_Header.bdf", ios::trunc);
		//ִ�п��Ʋ���
		ofs3 << "SOL 101" << endl;
		//������������
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Aircraft" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SUBCASE 1" << endl;//��һ������
		ofs3 << "   SUBTITLE=Default" << endl;
		ofs3 << "   LOAD = " << ++SID << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;
		ofs3 << "   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL" << endl;
		//BEGIN BULK
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		//ofs3 << "MDLPRM   HDF5    0" << endl;//����'.h5'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		ofs3 << "PARAM,INREL,-2" << endl;//��������ͷ�

		ofs3 << "include '" << fileName.substr(fileName.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << fileName1.substr(fileName1.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << fileName2.substr(fileName2.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		//����������������غ�
		ofs3 << "LOAD    " << setw(8) << SID << setw(8) << "1.";
		for (int i = 1; i < SID; i++)
		{
			ofs3 << setw(8) << "1." << setw(8) << i;
			if (i % 4 == 3)	ofs3 << endl << "        ";
		}
		ofs3 << endl;
		ofs3 << "ENDDATA" << endl;
		break;

	case 200://�Ż�����
		ofs3.open(DatName + "_Header.bdf", ios::trunc);
		//ִ�п��Ʋ���
		ofs3 << "SOL 200" << endl;
		//������������
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Aircraft" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SUBCASE 1" << endl;//��һ������
		ofs3 << "   SUBTITLE=Default" << endl;
		ofs3 << "   LOAD = " << ++SID << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;
		ofs3 << "   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL" << endl;
		//BEGIN BULK
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		//ofs3 << "MDLPRM   HDF5    0" << endl;//����'.h5'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		ofs3 << "PARAM,INREL,-2" << endl;//��������ͷ�

		ofs3 << "include '" << fileName.substr(fileName.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << fileName1.substr(fileName1.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << fileName2.substr(fileName2.find_last_of("/") + 1) << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		//����������������غ�
		ofs3 << "LOAD    " << setw(8) << SID << setw(8) << "1.";
		for (int i = 1; i < SID; i++)
		{
			ofs3 << setw(8) << "1." << setw(8) << i;
			if (i % 4 == 3)	ofs3 << endl << "        ";
		}
		ofs3 << endl;
		ofs3 << "ENDDATA" << endl;
		break;

	default:
		cout << "û���������Ӧ��SOL������ͣ����������룡" << endl;
		break;
	}
	

	ofs.close();
	ofs1.close();
	ofs2.close();
	ofs3.close();

	cout << "$----------------------------------------$" << endl;
	cout << "$        MSC.Nastran�����ʽ��ӡ���     $" << endl;
	cout << "$----------------------------------------$" << endl;
	return 0;
}

