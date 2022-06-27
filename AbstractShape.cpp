#include "CST_Instantiation.h"
#include <fstream>
#include <iomanip>
// #include<boost/assign.hpp>
// #include<boost/foreach.hpp>
#include<direct.h>
// using namespace boost::assign;

//��vector<double>����ת����mat
void AbstructShape::copyM(vector<vector<double>>& v, mat& m)
{
	vector<vector<double>>::iterator s = v.begin();
	//cout << v.size()<<" "<<(*s).size() << endl;//for test
	m = zeros(v.size(), (*s).size());
	//cout << "��ʼ�� m��" << endl;//for test
	//cout << m << endl;//for test
	int i = 0;
	for (vector<vector<double>>::iterator it = v.begin(); it != v.end(); it++)
	{
		int j = 0;
		//(*it)---���� vector<int>
		for (vector<double>::iterator vit = (*it).begin(); vit != (*it).end(); vit++)
		{
			//cout << *vit << " ";
			m(i, j) = *vit;
			j++;
		}
		i++;
	}
	//cout << "������ m��" << endl;//for test
	//cout << m << endl;//for test
}

void AbstructShape::SaveTecplotMesh()
{
	ofstream ofs;
	ofs.open("test.dat", ios::app);
	for(CSTsurface i: CSTsf)
	{
		ofs << "Title = \"example: 2d finite - element data\"" << endl;
		ofs << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs << "zone N=" << i.GridUpp.P.n_rows << ",E=" << i.GridUpp.E.n_rows << ",datapacking=block" << endl;
		ofs << "zonetype=fequadrilateral" << endl;
		ofs << trans(i.GridUpp.P) << endl;
		ofs << i.GridUpp.E << endl;

		ofs << "Title = \"example: 2d finite - element data\"" << endl;
		ofs << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs << "zone N=" << i.GridLow.P.n_rows << ",E=" << i.GridLow.E.n_rows << ",datapacking=block" << endl;
		ofs << "zonetype=fequadrilateral" << endl;
		ofs << trans(i.GridLow.P) << endl;
		ofs << i.GridLow.E << endl;
	}
	ofs.close();
}

//////////////////////////////////////////////
//Ĭ�ϼ�������Ϊ������
AbstructShape::AbstructShape()
{
	GridRefineType = 0;
	ShapeType = 0;
}

//���������������������
void AbstructShape::BacsShapeCSTs()
{
	int PartNum = CSTsf.size();
	GdStruct.reserve(CSTsf.size());

	if (GridRefineType > PartNum || GridRefineType < -1)
	{
		cout << "ERROR: Woring GridRefineType (value range: [-1, PartNum])" << endl;
		system("pause");
		exit(1);//������󣬹ҵ�����
	}
	else if (GridRefineType == -1)//�����м���
	{
		cout << "�������У������м��ܼ��㡣����" << endl;
		NotGridRefine();
	}
	else if (GridRefineType == 0)//������ȫ����
	{
		cout << "�������У���ȫ�������㡣����" << endl;
		AllGridRefine();
	}
	else//���ò�������
	{
		cout << "�������У������������㡣����" << endl;
		PrtGridRefine();
	}
}

//�����м���
void AbstructShape::NotGridRefine()
{
	const size_t PartNum = CSTsf.size();
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].CST3D();
		//CSTsf[iPart].RefineMesh(2, true);//
		CSTsf[iPart].RefineMesh(2, false);//
		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}

//������ȫ����
void AbstructShape::AllGridRefine()
{
	const size_t PartNum = CSTsf.size();

	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		CSTsf[iPart].EvalCurvature();
	}
	//�Ƚ�ÿ��CSTsf���������ֵ
	mat CURyzUpp = CSTsf[0].Info.CURyzUpp;
	mat CURyzLow = CSTsf[0].Info.CURyzLow;

	for (size_t iPart = 1; iPart < PartNum; iPart++)
	{
		//CURyzUpp = max(CURyzUpp, CSTsf[iPart].Info.CURyzUpp);
		//CURyzLow = max(CURyzLow, CSTsf[iPart].Info.CURyzLow);
		for (arma::uword i = 0; i < CURyzUpp.n_elem; i++)
		{
			CURyzUpp(i) = max(CURyzUpp(i), CSTsf[iPart].Info.CURyzUpp(i));
		}
		for (arma::uword i = 0; i < CURyzLow.n_elem; i++)
		{
			CURyzLow(i) = max(CURyzLow(i), CSTsf[iPart].Info.CURyzLow(i));
		}
	}
	cout << "�������У����ʷ�������ɣ�" << endl;
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].Info.CURyzUpp = CURyzUpp;
		CSTsf[iPart].Info.CURyzLow = CURyzLow;
		CSTsf[iPart].GridRefine();
		CSTsf[iPart].CST3D();

		//BoneLine.push_back(CSTsf[iPart].MakeBone());//******4.12

		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}

//���ò�������
void AbstructShape::PrtGridRefine()
{
	const size_t PartNum = CSTsf.size();
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 1;
		CSTsf[iPart].EvalCurvature();
	}
	MeshInfo Info0 = CSTsf[GridRefineType - 1].Info;
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].Info.CURyzUpp = Info0.CURyzUpp;
		CSTsf[iPart].Info.CURyzLow = Info0.CURyzLow;
		CSTsf[iPart].GridRefine();
		CSTsf[iPart].CST3D();

		//BoneLine.push_back(CSTsf[iPart].MakeBone());//******4.12

		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}
/////////////////////////////////////////////



void AbstructShape::SaveToAero(string name)
{
	ofstream ofs;
	ofs.open(name, ios::trunc);

	ofs << "$..~~~~~~GFFF����~~~~~" << endl;
	ofs << "$..~~~~~~ֻ��������������ڵ������κ��Ľڵ��ı��Σ�~~~~~" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$..         NODES" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$..CAA" << endl;

	int num = 0;
	for (size_t i = 0; i < CSTsf.size(); i++)//�ϱ���
	{
		for (arma::uword j = 0; j < CSTsf[i].GridUpp.P.n_rows; j++)
		{
			num++;

			ofs << "GRID*" << setw(19) << num << setw(32) << CSTsf[i].GridUpp.P(j, 0)
				<< setw(16) << CSTsf[i].GridUpp.P(j, 1) << setw(8) << num << endl;
			ofs << "*" << setw(7) << num << setw(16) << CSTsf[i].GridUpp.P(j, 2) << endl;
		}
	}
	for (size_t i = 0; i < CSTsf.size(); i++)//�±���
	{
		for (size_t j = 0; j < CSTsf[i].GridLow.P.n_rows; j++)
		{
			num++;

			ofs << "GRID*" << setw(19) << num << setw(32) << CSTsf[i].GridLow.P(j, 0)
				<< setw(16) << CSTsf[i].GridLow.P(j, 1) << setw(8) << num << endl;
			ofs << "*" << setw(7) << num << setw(16) << CSTsf[i].GridLow.P(j, 2) << endl;
		}
	}

	ofs << "$..~~~~~~����������GFFF����~~~~~" << endl;
	ofs << "$----------------------------------------" << endl;
	ofs << "$.." << endl;
	ofs << "$..          ELEMENTS" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$----------------------------------------" << endl;
	ofs << "$.." << endl;
	ofs << "$..          MESH PART: MSHPartSmartSurf.1(CAA)" << endl;
	ofs << "$.." << endl;

	num = 0;
	int pointNum = 0;
	for (size_t i = 0; i < CSTsf.size(); i++)
	{
		
		for (size_t j = 0; j < CSTsf[i].GridUpp.E.n_rows; j++)
		{
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1
				<< setw(8) << CSTsf[i].GridUpp.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 1) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 2) + pointNum << endl;
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1
				<< setw(8) << CSTsf[i].GridUpp.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 2) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 3) + pointNum << endl;
		}
		pointNum += CSTsf[i].GridUpp.P.n_rows;
	}
	for (size_t i = 0; i < CSTsf.size(); i++)
	{

		for (size_t j = 0; j < CSTsf[i].GridLow.E.n_rows; j++)
		{
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1 + CSTsf.size()
				<< setw(8) << CSTsf[i].GridLow.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 1) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 2) + pointNum << endl;
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1 + CSTsf.size()
				<< setw(8) << CSTsf[i].GridLow.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 2) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 3) + pointNum << endl;
		}
		pointNum += CSTsf[i].GridLow.P.n_rows;
	}

	ofs << "$.." << endl;
	ofs << "ENDDATA" << endl;
	ofs << endl;
	ofs << "��������GFFF����--����ѧԺ ������ҵ��ѧ" << endl;
}

void AbstructShape::MakeAeroInfo(string aeropath,int num)
{
	//string name = "mesh//mesh" + to_string(num) + "ָ��.txt";
	string datpath = aeropath + "mesh" + to_string(num) + ".dat";//���������ļ�
	string name = aeropath + "mesh" + to_string(num) + "_cmd.txt";//����ָ���ļ�
	string outputpath = aeropath + "mesh" + to_string(num);//����ļ���ַ
	int state = _mkdir(outputpath.c_str());
	int size = CSTsf.size() * 2;
	ofstream ofs;
	ofs.open(name, ios::trunc);
	ofs << "//BEGIN//" << endl;
	ofs << endl;
	ofs << "//////////////////����״̬(���飩///////////////////////" << endl;
	ofs << "mach	nAlpha   nBeta  nH x�᷽��(1Ϊͷ��ָ��β����-1Ϊ������)" << endl;
	ofs << "1       1        1      1  1" << endl;
	ofs << "mach    //�����" << endl;
	ofs << "10 " << endl;
	ofs << "alpha	//����(deg)" << endl;
	ofs << "30 " << endl;
	ofs << "beta 	//�໬��(deg)" << endl;
	ofs << "0 " << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////////���и߶�(m)////////////////////////" << endl;
	ofs << "20000" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////////��������///////////////////////////" << endl;
	ofs << "��������� " << size << endl;
	ofs << "//////�������ID" << endl;
	for (int i = 0; i < size; i++)	ofs << i + 1 << " ";
	ofs << endl;
	ofs << "//////�����������" << endl;
	for (int i = 0; i < size; i++)	ofs << 2 << " ";
	ofs << endl;
	ofs << "//////�����������" << endl;
	for (int i = 0; i < size; i++)	ofs << 1 << " ";
	ofs << endl;
	ofs << "///////�����ⷨ�߷���" << endl;
	for (int i = 0; i < size; i++)
	{
		int mid = size / 2;
		int direction = -1;
		if (i < mid)
		{
			direction = 1;
		}
		ofs << direction << " ";
	}
	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////������////////////////////////////" << endl;
	ofs << "���������� 0" << endl;
	ofs << "/////ÿ������������ռһ�У���������������ö�Ӧ���̶���ת��(0�ޣ�1X,2Y,3Z)������֮�䲻Ҫ���ո��ɶ��Ÿ�������ת����P2->P1����ʱ����תΪ������λm" << endl;
	ofs << "^^^��Ӧ����^^^^^^^��ת��(��1����)^^^^^^^��ת��(��2����)^^^^^^^�̶���ת��^^^��ƫ����^^^��ƫ������^^^^^^^^^^^" << endl;
	ofs << "" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////sRef�ο����//cRef�ο�����//�ο��ߴ�(m)////////" << endl;
	ofs << "1 1 1 0 0" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////�Ƿ���ճ��(��Ϊ1����Ϊ0)/////////////" << endl;
	ofs << "����ճ�Լ���(0,1)    0" << endl;
	ofs << "ͷ���뾶(m)           0.03" << endl;
	ofs << "�����¶�(K)           300" << endl;
	ofs << "���������            0.7" << endl;
	ofs << "����ģ��(0,1)        1" << endl;
	ofs << "����ϵ��              0" << endl;
	ofs << "�����¶�(K)           280" << endl;
	ofs << "�ڱ����¶�(K)         350" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "////////////////�Ƕ�������ѡ��(��Ϊ1����Ϊ0)////////////" << endl;
	ofs << "�Ƕ�������   �������:ʱ����(t)    �ڵ��ٶ�rpt�ĵ� " << endl;
	ofs << "0            0.0005                  velocity.rpt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////����ļ�//////////////////////////" << endl;
	ofs << "//�����ļ���λ(m/mm):		mm" << endl;
	ofs << "//��������                  CATIA" << endl;
	ofs << "�����ļ���                  " << datpath << endl;//C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh"
	ofs << "���·����                  " << outputpath << endl;//C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh\\mesh"
	ofs << "�������ļ�����            CalResult - ȫ����״̬.txt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "�Ƿ����ѹ���ֲ��ļ�plot.plt(plot3D��ʽ):     YES" << endl;
	ofs << "�Ƿ������Ԫ���ֲ��ļ�face.dat:               NO" << endl;
	ofs << "�Ƿ�����ڵ����ֲ��ļ�node.dat:               YES" << endl;
	ofs << "�Ƿ�����¶ȷֲ��ļ�temperature.plt(��Ҫճ�Լ���) NO" << endl;
	ofs << "�Ƿ�������طֲ��ļ�qload.plt(��Ҫճ�Լ���)       NO" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << endl;
	ofs << endl;
	ofs <<  endl;
	ofs << "//END//" << endl;
}
