#include "NasPrinter.h"
#include <direct.h>
#include <iomanip>
using namespace std;
#define C_INT const int&
#define C_DOUBEL const double&
#define C_STR const std::string&
#define C_VECTOR(type) const std::vector<type>&

void Printer::init(int size,int precision)
{
	m_precision = precision;
	str_list.clear();
	str_list.reserve(size);
}

NasPrinter::NasPrinter() 
{ 
	ssHeader.clear();	ssHeader.str("");
	ssForce.clear();	ssForce.str("");
	ssForce.clear();	ssMesh.str("");
	ssProperty.clear();	ssProperty.str("");
}
int NasPrinter::PrintBDF(const std::string& path, const std::string & name, OutputStyle type)
{	
	int state = _mkdir(path.c_str());
	switch (type)
	{
	case OutputStyle::onefile:{
		ofstream ofs;
		ofs.open(path + name + "_Header.bdf", ios::trunc);
		ofs << ssHeader.str();
		ofs << ssMesh.str();
		ofs << ssForce.str();
		ofs << ssProperty.str();
		ofs << "ENDDATA" << endl;
		ofs.close();
		return 0;
	}	break;
	case OutputStyle::multifile: {
		ofstream ofsHeader,ofsForce,ofsMesh,ofsProperty;
		//
		ofsHeader.open(path + name + "_Header.bdf", ios::trunc);
		ofsHeader << ssHeader.str();
		ofsHeader << "include '" << (name + "_mesh.bdf'") << endl;		//ͷ�ļ��ڲ�ʹ�����·��
		ofsHeader << "include '" << (name + "_property.bdf'") << endl;	//ͷ�ļ��ڲ�ʹ�����·��
		ofsHeader << "include '" << (name + "_force.bdf'") << endl;		//ͷ�ļ��ڲ�ʹ�����·��
		ofsHeader << "ENDDATA" << endl;
		ofsHeader.close();
		//
		ofsForce.open(path + name + "_force.bdf", ios::trunc);
		ofsForce << ssForce.str();
		ofsForce.close();
		//
		ofsMesh.open(path + name + "_mesh.bdf", ios::trunc);
		ofsMesh << ssMesh.str();
		ofsMesh.close();
		//
		ofsProperty.open(path + name + "_property.bdf", ios::trunc);
		ofsProperty << ssProperty.str();
		ofsProperty.close();
		return 0;
	}	break;
	default:
		return -1;
		break;
	}
	
	return state;
}

void NasPrinter::AddCardShort(std::ofstream& ofs, const std::vector<std::string>& words)
{
	int listNum = 0;
	ofs << std::left;
	ofs << setw(8) << words[0].substr(0, 8);//��һ���ؼ���
	ofs << std::right;
	for (size_t i = 1; i < words.size() - 1; i++)
	{
		listNum++;
		ofs << setw(8) << words[i].substr(0, 8);
		if (listNum % 8 == 0)	ofs << "\n        ";
	}
	ofs << setw(8) << words[words.size() - 1];

	ofs << endl;
}

void NasPrinter::AddCardLong(std::ofstream& ofs, const std::vector<std::string>& words)
{
	int listNum = 0;
	ofs << std::left;

	ofs << setw(8) << words[0].substr(0,8) + "*";//��һ���ؼ���
	ofs << std::right;
	for (size_t i = 1; i < words.size() - 1; i++)
	{
		listNum++;
		ofs << setw(16) << words[i].substr(0, 16);
		if (listNum % 4 == 0)	ofs << "*\n*       ";
	}
	ofs << setw(16) << words[words.size() - 1];

	ofs << endl;
}

void NasPrinter::AddCardShort(std::stringstream& ofs, const std::vector<std::string>& words)
{
	//ofs.str("");//��ջ�����
	int listNum = 0;
	ofs << std::left;
	ofs << setw(8) << words[0].substr(0, 8);;//��һ���ؼ���
	ofs << std::right;
	for (size_t i = 1; i < words.size()-1; i++)
	{
		listNum++;
		ofs << setw(8) << words[i].substr(0, 8);;
		if (listNum % 8 == 0)	ofs << "\n        ";
	}
	ofs << setw(8) << words[words.size() - 1];

	ofs << endl;
}

void NasPrinter::AddCardLong(std::stringstream& ofs, const std::vector<std::string>& words)
{
	//ofs.str("");//��ջ�����
	int listNum = 0;
	ofs << std::left;

	ofs << setw(8) << words[0].substr(0, 8) + "*";//��һ���ؼ���
	ofs << std::right;
	for (size_t i = 1; i < words.size() - 1; i++)
	{
		listNum++;
		ofs << setw(16) << words[i].substr(0, 16);
		if (listNum % 4 == 0)	ofs << "*\n*       ";
	}
	ofs << setw(16) << words[words.size() - 1];

	ofs << endl;
}

string Printer::Printer::double2str(const double& db, int size)
{
	stringstream ss;
	const double fabsdb = abs(db);
	int sign = 1;
	if (db < 0)
	{
		size--;
		sign = -1;
	}
	//����
	const int xx = (int)fabsdb;
	if ((fabsdb - xx) < 0.4 && abs(xx - fabsdb) < 1e-5)//����
	{
		return to_string((int)db) + ".";
	}
	if ((fabsdb - xx >= 0.4 && abs(xx - fabsdb + 1) < 1e-5))//����
	{
		return to_string((int)db + sign) + ".";
	}
	//������
	if (fabsdb > std::pow(10.0, size))
	{
		ss << setprecision(size - 6) << db;
	}
	else if(fabsdb >= 1)
	{
		ss << setprecision(size - 1) << db;
	}
	else if (fabsdb > 1e-1)
	{
		ss << setprecision(size - 2) << db;
	}
	else if (fabsdb > 1e-2)
	{
		ss << setprecision(size - 3) << db;
	}
	else if (fabsdb > 1e-3)
	{
		ss << setprecision(size - 4) << db;
	}
	else if (fabsdb > 1e-4)
	{
		ss << setprecision(size - 5) << db;
	}
	else
	{
		ss << setprecision(size - 6) << db;
	}

	string res = ss.str();
	if (res.find(".") == string::npos)
	{
		res += ".";
	}
	return res;
	
}

void Printer::add_usize(const int &usz)
{
	if (usz < 0)
		add("");
	else
		add(usz);
}

void NasPrinter::addGRID(C_INT CID, C_INT CP, const Point & pt)
{
	m_print.init(6, 16);
	m_print.add("GRID");	//��Ƭ��
	m_print.add(CID);		//����ϵ���
	m_print.add(CP);		//x����
	m_print.add(pt);		//z����
	AddCardLong(ssMesh, m_print.str_list);
	// AddCardLong(ssMesh, {
	// 	"GRID",						//��Ƭ��
	// 	Printer::int2str(CID),				//�ڵ��� CID
	// 	Printer::int2str(CP),				//����ϵ���
	// 	Printer::double2str(pt.getX(),16),	//x����
	// 	Printer::double2str(pt.getY(),16),	//y����
	// 	Printer::double2str(pt.getZ(),16),	//z����
	// });
}

void NasPrinter::addCMASS2(C_INT EID, C_DOUBEL M, C_INT G1, C_INT C1)
{
	m_print.init(5, 8);
	m_print.add("CMASS2");	//��Ƭ��
	m_print.add(EID);		//��Ԫ��� EID
	m_print.add(M);		//��Ԫ������M
	m_print.add(G1);		//�ڵ�� G1
	m_print.add(C1);		//������ C1
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CMASS2",					//��Ƭ��
	// 	Printer::int2str(EID),				//��Ԫ��� EID
	// 	Printer::double2str(M, 8),			//��Ԫ������M
	// 	Printer::int2str(G1),				//�ڵ�� G1
	// 	Printer::int2str(C1),				//������ C1
	// });
}

void NasPrinter::addCBAR(C_INT EID, C_INT PID, C_INT GA, C_INT GB, C_DOUBEL X1, C_DOUBEL X2, C_DOUBEL X3)
{
	m_print.init(8, 8);
	m_print.add("CBAR");	//��Ƭ��
	m_print.add(EID);		//��Ԫ��� EID
	m_print.add(PID);		//��Ԫ���Ա�� PID
	m_print.add(GA);		//�ڵ�1
	m_print.add(GB);		//�ڵ�2
	m_print.add(X1);		//һά��Ԫ����ָ��
	m_print.add(X2);		//һά��Ԫ����ָ��
	m_print.add(X3);		//һά��Ԫ����ָ��
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CBAR",									//��Ƭ��
	// 	Printer::int2str(EID),					//��Ԫ��� EID
	// 	Printer::int2str(PID),					//��Ԫ���Ա�� PID
	// 	Printer::int2str(GA),					//�ڵ�1
	// 	Printer::int2str(GB),					//�ڵ�2
	// 	Printer::double2str(X1, 8),				//һά��Ԫ����ָ��
	// 	Printer::double2str(X2, 8),				//һά��Ԫ����ָ��
	// 	Printer::double2str(X3, 8),				//һά��Ԫ����ָ��
	// });
}

void NasPrinter::addCTRIA3(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3)
{
	m_print.init(6, 8);
	m_print.add("CTRIA3");	//��Ƭ��
	m_print.add(EID);		//��Ԫ��� EID
	m_print.add(PID);		//��Ԫ���Ա�� PID
	m_print.add(G1);		//�ڵ�1 CID
	m_print.add(G2);		//�ڵ�2 CID
	m_print.add(G3);		//�ڵ�3 CID
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CTRIA3",				//��Ƭ��
	// 	Printer::int2str(EID),	//��Ԫ��� EID
	// 	Printer::int2str(PID),	//��Ԫ���Ա�� PID
	// 	Printer::int2str(G1),	//�ڵ�1 CID
	// 	Printer::int2str(G2),	//�ڵ�2 CID
	// 	Printer::int2str(G3),	//�ڵ�3 CID
	// });
}

void NasPrinter::addCQUAD4(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3, C_INT G4)
{
	m_print.init(7, 8);
	m_print.add("CQUAD4");	//��Ƭ��
	m_print.add(EID);		//��Ԫ��� EID
	m_print.add(PID);		//��Ԫ���Ա�� PID
	m_print.add(G1);		//�ڵ�1 CID
	m_print.add(G2);		//�ڵ�2 CID
	m_print.add(G3);		//�ڵ�3 CID
	m_print.add(G4);		//�ڵ�4 CID
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CQUAD4",				//��Ƭ��
	// 	Printer::int2str(EID),	//��Ԫ��� EID
	// 	Printer::int2str(PID),	//��Ԫ���Ա�� PID
	// 	Printer::int2str(G1),	//�ڵ�1 CID
	// 	Printer::int2str(G2),	//�ڵ�2 CID
	// 	Printer::int2str(G3),	//�ڵ�3 CID
	// 	Printer::int2str(G4),	//�ڵ�4 CID
	// });
}

void NasPrinter::addRBE3(C_INT EID, C_INT REDGRID, C_INT REFC, C_DOUBEL WT1, C_INT C1, C_VECTOR(int) G1j)
{
	m_print.init(7 + G1j.size(), 8);
	m_print.add("RBE3");	//��Ƭ��
	m_print.add(EID);		//���ݿ����
	m_print.add("");		//
	m_print.add(REDGRID);	//�ο��ڵ���	
	m_print.add(REFC);		//�ο��ڵ����ɶȷ�����
	m_print.add(WT1);		//�ڵ�G1j��λ�����ļ�Ȩֵ������ϵ����
	m_print.add(C1);		//�ڵ�G1j��λ�Ʒ�����
	m_print.add(G1j);		//����ڵ�
	AddCardShort(ssMesh, m_print.str_list);
	// vector<string> strList = {
	// 	"RBE3",							//��Ƭ��
	// 	Printer::int2str(EID),			//���ݿ����
	// 	"",								//
	// 	Printer::int2str(REDGRID),		//�ο��ڵ���
	// 	Printer::int2str(REFC),			//�ο��ڵ����ɶȷ�����
	// 	Printer::double2str(WT1, 8),	//�ڵ�G1j��λ�����ļ�Ȩֵ������ϵ����
	// 	Printer::int2str(C1),			//�ڵ�G1j��λ�Ʒ�����
	// };
	// strList.reserve(strList.size() + G1j.size());
	// for (size_t j = 0; j < G1j.size(); j++)
	// {
	// 	strList.push_back(Printer::int2str(G1j[j]));//����ڵ�
	// }
	// AddCardShort(ssMesh, strList);
}

void NasPrinter::addPLOAD(C_INT SID, C_DOUBEL P, C_INT G1, C_INT G2, C_INT G3, C_INT G4)
{
	m_print.init(7, 8);
	m_print.add("PLOAD");	//��Ƭ��
	m_print.add(SID);		//�غɱ�� SID
	m_print.add(P);		//ѹ��ֵ
	m_print.add(G1);		//�ڵ�1
	m_print.add(G2);		//�ڵ�2
	m_print.add(G3);		//�ڵ�3
	m_print.add_usize(G4);	//�ڵ�4 ȡĬ��ֵ-1ʱĬ�����Ϊ��
	AddCardShort(ssForce, m_print.str_list);
	// AddCardShort(ssForce, {
	// 	"PLOAD",					//��Ƭ��
	// 	Printer::int2str(SID),		//�غɱ�� SID
	// 	Printer::double2str(P, 8),	//ѹ��ֵ
	// 	Printer::int2str(G1),		//�ڵ�1
	// 	Printer::int2str(G2),		//�ڵ�2
	// 	Printer::int2str(G3),		//�ڵ�3
	// 	Printer::int2str(G4),		//�ڵ�4 ȡĬ��ֵ-1ʱĬ�����Ϊ��
	// });
}

void NasPrinter::addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3)
{
	m_print.init(8, 16);
	m_print.add("FORCE");	//��Ƭ��
	m_print.add(SID);		//SID
	m_print.add(G);		//�ڵ���
	m_print.add(CID);		//�ֲ�����ϵ��� Ĭ�Ͽɲ���=0��ȫ������ϵ��
	m_print.add(F);		//�������ķ�ֵ
	m_print.add(N1);		//������������ϵ�е�ʸ������
	m_print.add(N2);		//������������ϵ�е�ʸ������
	m_print.add(N3);		//������������ϵ�е�ʸ������
	AddCardLong(ssForce, m_print.str_list);
	// AddCardLong(ssForce, {
	// 	"FORCE",						//��Ƭ��
	// 	Printer::int2str(SID),			//SID
	// 	Printer::int2str(G),			//�ڵ���
	// 	Printer::int2str(CID),			//�ֲ�����ϵ��� Ĭ�Ͽɲ���=0��ȫ������ϵ��
	// 	Printer::double2str(F , 16),	//�������ķ�ֵ
	// 	Printer::double2str(N1,16),		//������������ϵ�е�ʸ������
	// 	Printer::double2str(N2,16),		//������������ϵ�е�ʸ������
	// 	Printer::double2str(N3,16),		//������������ϵ�е�ʸ������
	// });
}

void NasPrinter::addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, const Point & pt)
{
	addFORCE(SID, G, CID, F, pt.getX(), pt.getY(), pt.getZ());
}

void NasPrinter::addGRAV(C_INT SID, C_INT CID, C_DOUBEL A, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3)
{
	m_print.init(7, 16);
	m_print.add("GRAV");	//��Ƭ��
	m_print.add(SID);		//�غɱ�� SID
	m_print.add(CID);		//����ϵ���
	m_print.add(A);		//���ٶȷ�ֵ
	m_print.add(N1);		//���ٶȷ�������1
	m_print.add(N2);		//���ٶȷ�������2
	m_print.add(N3);		//���ٶȷ�������3
	AddCardLong(ssForce, m_print.str_list);
	// AddCardLong(ssForce, {
	// 	"GRAV",							//��Ƭ��
	// 	Printer::int2str(SID),			//�غɱ�� SID
	// 	Printer::int2str(CID),			//����ϵ���
	// 	Printer::double2str(A, 8),		//���ٶȷ�ֵ
	// 	Printer::double2str(N1, 8),		//���ٶȷ�������1
	// 	Printer::double2str(N2, 8),		//���ٶȷ�������2
	// 	Printer::double2str(N3, 8),		//���ٶȷ�������3
	// });
}

void NasPrinter::addLOAD(C_INT SID, C_DOUBEL S, C_VECTOR(DoubleIntPair) LOADi)
{
	m_print.init(3 + LOADi.size() * 2, 8);
	m_print.add("LOAD");	//��Ƭ��
	m_print.add(SID);		//���ݿ����
	m_print.add(S);		//��ϵ��
	m_print.add(LOADi);	//�����غɵ�ϵ��, �����غ����ݿ����
	AddCardShort(ssForce, m_print.str_list);
	// vector<string> strList = {
	// 	"LOAD",						//��Ƭ��
	// 	Printer::int2str(SID),		//���ݿ����
	// 	Printer::double2str(S, 8),	//��ϵ��
	// };
	// strList.reserve(strList.size() + LOADi.size() * 2);
	// for (size_t i = 0; i < LOADi.size(); i++)
	// {
	// 	strList.push_back(Printer::double2str(LOADi[i].first, 8));	//�����غɵ�ϵ��
	// 	strList.push_back(Printer::int2str(LOADi[i].second));			//�����غ����ݿ����
	// }
	// AddCardShort(ssForce, strList);
}

void NasPrinter::addSPC1(C_INT SID, C_INT C, C_VECTOR(int) Gi)
{
	m_print.init(3 + Gi.size(), 8);
	m_print.add("SPC1");	//��Ƭ��
	m_print.add(SID);		//���ݿ����
	m_print.add(C);		//��ϵ��
	m_print.add(Gi);		//�����غɵ�ϵ��
	AddCardShort(ssForce, m_print.str_list);
	// vector<string> strList = {
	// 	"SPC1",						//��Ƭ��
	// 	Printer::int2str(SID),		//���ݿ����
	// 	Printer::int2str(C),		//��ϵ��
	// };
	// strList.reserve(strList.size() + Gi.size());
	// for (size_t i = 0; i < Gi.size(); i++)
	// {
	// 	strList.push_back(Printer::int2str(Gi[i]));	//�����غɵ�ϵ��
	// }
	// AddCardShort(ssForce, strList);
}

void NasPrinter::addPSHELL(const PSHELL & ps)
{
	ssProperty << "$-- Property " << ps.name << " --$" << endl;
	m_print.init(4, 8);
	m_print.add("PSHELL");	//��Ƭ��
	m_print.add(ps.PID);	//��Ԫ���Ա�� PID
	m_print.add(ps.MID);	//�������Ա�� MID
	m_print.add(ps.T);		//���
	AddCardShort(ssProperty, m_print.str_list);
	// AddCardShort(ssProperty, {
	// 	"PSHELL",						//��Ƭ��
	// 	Printer::int2str(ps.PID),		//��Ԫ���Ա�� PID
	// 	Printer::int2str(ps.MID),		//�������Ա�� MID
	// 	Printer::double2str(ps.T, 8),	//���
	// });
}

void NasPrinter::addPBARL(const PBARL & pb)
{
	ssProperty << "$-- Property " << pb.name << " --$" << endl;
	m_print.init(9 + pb.DIM.size(), 8);
	m_print.add("PBARL");	//��Ƭ��
	m_print.add(pb.PID);	//��Ԫ���Ա�� PID
	m_print.add(pb.MID);	//�������Ա�� MID
	m_print.add("");		//
	m_print.add(pb.TYPE);	//��Ԫ������������
	m_print.add("");		//
	m_print.add("");		//
	m_print.add("");		//
	m_print.add("");		//
	m_print.add(pb.DIM);	//�������Բ���
	AddCardShort(ssProperty, m_print.str_list);
	// vector<string> strList = {
	// 	"PBARL",			//��Ƭ��
	// 	Printer::int2str(pb.PID),	//�������Ա�� PID
	// 	Printer::int2str(pb.MID),
	// 	"",
	// 	pb.TYPE,
	// 	"","","",""	
	// };
	// //const size_t orginSize = strList.size();
	// strList.reserve(strList.size() + pb.DIM.size());
	// for (size_t i = 0; i < pb.DIM.size(); i++)
	// {
	// 	strList.push_back(Printer::double2str(pb.DIM[i], 8));	//
	// }
	// AddCardShort(ssProperty, strList);
}

void NasPrinter::addMAT8(const MAT8 & mt)
{
	ssProperty << "$-- Material " << mt.name << " --$" << endl;
	m_print.init(9, 8);
	m_print.add("MAT8");
	m_print.add(mt.MID);
	m_print.add(mt.E1);
	m_print.add(mt.E2);
	m_print.add(mt.NU12);
	m_print.add(mt.G12);
	m_print.add(mt.G1Z);
	m_print.add(mt.G2Z);
	m_print.add(mt.RHO);
	AddCardShort(ssProperty, m_print.str_list);
	// AddCardShort(ssProperty, {
	// 	"MAT8",
	// 	Printer::int2str(mt.MID),
	// 	Printer::double2str(mt.E1  , 8),
	// 	Printer::double2str(mt.E2  , 8),
	// 	Printer::double2str(mt.NU12, 8),
	// 	Printer::double2str(mt.G12 , 8),
	// 	Printer::double2str(mt.G1Z , 8),
	// 	Printer::double2str(mt.G2Z , 8),
	// 	Printer::double2str(mt.RHO , 8),
	// });
}

void NasPrinter::addMAT1(const MAT1 & mt)
{
	ssProperty << "$-- Material " << mt.name << " --$" << endl;
	m_print.init(6, 8);
	m_print.add("MAT1");
	m_print.add(mt.MID);
	m_print.add(mt.E);
	m_print.add("");
	m_print.add(mt.Nu);
	m_print.add(mt.RHO);
	AddCardShort(ssProperty, m_print.str_list);
	// AddCardShort(ssProperty, {
	// 	"MAT1",
	// 	Printer::int2str(mt.MID),
	// 	Printer::double2str(mt.E, 8),
	// 	"",
	// 	Printer::double2str(mt.Nu, 8),
	// 	Printer::double2str(mt.RHO, 8),
	// });
}

void NasPrinter::addEIGRL(C_INT SID, C_INT ND, C_INT MSGLVL)
{
	//Lanczos����ȡģ̬
	m_print.init(5, 8);
	m_print.add("EIGRL");	//��Ƭ��
	m_print.add(SID);		//����ţ�����>0��
	m_print.add("");		//
	m_print.add(ND);		//Number of Desired Roots��������������>0,��ո�
	m_print.add(MSGLVL);	//Diagnostic Output Level�Ի�����Ϣ������0����Ϣ��1ÿ��ƽ�ơ���ӡ��õ�����ֵ��>1���ӶԻ������
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	Printer::int2str(SID),		//����ţ�����>0��
	// 	"","",
	// 	Printer::int2str(ND),		//Number of Desired Roots��������������>0,��ո�
	// 	Printer::int2str(MSGLVL),	//Diagnostic Output Level�Ի�����Ϣ������0����Ϣ��1ÿ��ƽ�ơ���ӡ��õ�����ֵ��>1���ӶԻ������
	// });
}

void NasPrinter::addDOPTPRM(C_VECTOR(std::string) command)
{
	m_print.init(1 + command.size(), 8);
	m_print.add("DOPTPRM");
	m_print.add(command);
	AddCardShort(ssHeader, m_print.str_list);
	// vector<string> strList = {
	// 	"DOPTPRM",
	// };
	// strList.reserve(strList.size() + command.size());
	// for (size_t i = 0; i < command.size(); i++)
	// {
	// 	strList.push_back(command[i]);
	// }
	// AddCardShort(ssHeader, strList);
}

void NasPrinter::addDESVAR(C_INT ID, C_STR LABEL, C_DOUBEL XINIT, C_DOUBEL XLB, C_DOUBEL XUB, C_DOUBEL DELXV)
{
	//���������ṩһ����ֵ���ļ���
	//������Ʊ����Ŀ�Ƭ
	m_print.init(7, 8);
	m_print.add("DESVAR");	//��Ƭ��
	m_print.add(ID);		//��Ʊ�����ID�ţ�int>0��
	m_print.add(LABEL);	//�û�ָ������Ʊ�������
	m_print.add(XINIT);	//��Ʊ����ĳ�ʼֵ��double XLB<=XINIT<=XUB��
	m_print.add(XLB);		//��Ʊ��������ޣ�double nasĬ��ֵ-1.0e20��
	m_print.add(XUB);		//��Ʊ��������ޣ�double nasĬ��ֵ+1.0e20��
	m_print.add(DELXV);	//��Ʊ����ĸñ仯�����ʼֵ�ı�ֵ�����ֵ��double>0��
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	"DESVAR",						//��Ƭ��
	// 	Printer::int2str(ID),			//��Ʊ�����ID�ţ�int>0��
	// 	LABEL,							//�û�ָ������Ʊ�������
	// 	Printer::double2str(XINIT,8),	//��Ʊ����ĳ�ʼֵ��double XLB<=XINIT<=XUB��
	// 	Printer::double2str(XLB, 8),	//��Ʊ��������ޣ�double nasĬ��ֵ-1.0e20��
	// 	Printer::double2str(XUB, 8),	//��Ʊ��������ޣ�double nasĬ��ֵ+1.0e20��
	// 	Printer::double2str(DELXV,8),	//��Ʊ����ĸñ仯�����ʼֵ�ı�ֵ�����ֵ��double>0��
	// });
}

void NasPrinter::addDVPREL1(C_INT ID, C_STR TYPE, C_INT PID, C_STR PNAME, C_VECTOR(DoubleIntPair) DVID_COEFi, C_DOUBEL C0, C_DOUBEL PMIN, C_DOUBEL PMAX)
{
	m_print.init(9 + DVID_COEFi.size() * 2, 8);
	m_print.add("DVPREL1");	//��Ƭ��
	m_print.add(ID);			//����Ʊ������ӵ�ID�ţ�int>0��
	m_print.add(TYPE);			//�������͵����֣��硰PBAR������PBEAM���ȣ�string��
	m_print.add(PID);			//���Ե�ID�ţ�int>0��
	m_print.add(PNAME);		//�������֣��硰T������A���������Կ�Ƭ����Field����δ֪���ߵ�Ԫ���Ա����ַ���Word����λ�ã�sring��int>0��
	m_print.add(PMIN);			//�������������Сֵ��double��nasĬ��ֵ-1.0e20��
	m_print.add(PMAX);			//��������������ֵ��double��nasĬ��ֵ+1.0e20��
	m_print.add(C0);			//���Թ�ϵ�ĳ����double��nasĬ��ֵΪ0.0��
	m_print.add("");			//
	m_print.add(DVID_COEFi);	//DVIDi DESVAR��Ƭ��ID�ţ�int>0��, COEFi ���Թ�ϵ��һ����ϵ����double��
	AddCardShort(ssHeader, m_print.str_list);
	// vector<string> strList = {
	// 	"DVPREL1",
	// 	Printer::int2str(ID),			//����Ʊ������ӵ�ID�ţ�int>0��
	// 	TYPE,							//�������͵����֣��硰PBAR������PBEAM���ȣ�string��
	// 	Printer::int2str(PID),			//���Ե�ID�ţ�int>0��
	// 	PNAME,							//�������֣��硰T������A���������Կ�Ƭ����Field����δ֪���ߵ�Ԫ���Ա����ַ���Word����λ�ã�sring��int>0��
	// 	Printer::double2str(PMIN, 8),	//�������������Сֵ��double��nasĬ��ֵ-1.0e20��
	// 	Printer::double2str(PMAX, 8),	//��������������ֵ��double��nasĬ��ֵ+1.0e20��
	// 	Printer::double2str(C0, 8),		//���Թ�ϵ�ĳ����double��nasĬ��ֵΪ0.0��
	// 	"",
	// };
	// strList.reserve(strList.size() + DVID_COEFi.size());
	// for (size_t i = 0; i < DVID_COEFi.size(); i++)
	// {
	// 	strList.push_back(Printer::int2str(DVID_COEFi[i].second));//DVIDi DESVAR��Ƭ��ID�ţ�int>0��
	// 	strList.push_back(Printer::double2str(DVID_COEFi[i].first, 8));//COEFi ���Թ�ϵ��һ����ϵ����double��
	// }
	// AddCardShort(ssHeader, strList);
}

void NasPrinter::addDRESP1(C_INT ID, C_STR LABLE, C_STR RTYPE, C_STR PTYPE, C_INT REGION, C_INT ATTA, C_INT ATTB, C_VECTOR(int) ATTi)
{
	m_print.init(8 + ATTi.size(), 8);
	m_print.add("DRESP1");
	m_print.add(ID);
	m_print.add(LABLE);
	m_print.add(RTYPE);
	m_print.add(PTYPE);
	m_print.add_usize(REGION);
	m_print.add_usize(ATTA);
	m_print.add_usize(ATTB);
	m_print.add(ATTi);
	AddCardShort(ssHeader, m_print.str_list);
	// vector<string> strList = {
	// 	"DRESP1",
	// 	Printer::int2str(ID),
	// 	LABLE,
	// 	RTYPE,
	// 	PTYPE,
	// 	Printer::int2str(REGION),
	// 	Printer::int2str(ATTA),
	// 	Printer::int2str(ATTB),
	// };
	// strList.reserve(strList.size() + ATTi.size());
	// for (size_t i = 0; i < ATTi.size(); i++)
	// {
	// 	strList.push_back(Printer::int2str(ATTi[i]));
	// }
	// AddCardShort(ssHeader, strList);
}

void NasPrinter::addDRESP1_WEIGHT(C_INT ID, C_STR LABEL)
{
	addDRESP1(ID, LABEL, "WEIGHT", "", -1, -1, -1, vector<int>(0));
}

void NasPrinter::addDRESP1_FREQ(C_INT ID, C_STR LABEL, C_INT ATTA)
{
	addDRESP1(ID, LABEL, "FREQ","", -1, ATTA, -1, vector<int>(0));
}

void NasPrinter::addDRESP1_STRESS(C_INT ID, C_STR LABEL, C_STR PTYPE, C_INT ATTA, C_VECTOR(int) ATTi)
{
	addDRESP1(ID, LABEL, "STRESS", PTYPE, -1, ATTA, -1, ATTi);
}

void NasPrinter::addDCONSTR(C_INT DCID, C_INT RID, C_DOUBEL LALLOW, C_DOUBEL UALLOW)
{
	//����Լ��
	m_print.init(5, 8);
	m_print.add("DCONSTR");//��Ƭ��
	m_print.add(DCID);		//���Լ����ID�ţ�int>0��
	m_print.add(RID);		//�����Ӧ��DRESPi��Ƭ�������ID�ţ�int>0��
	m_print.add(LALLOW);	//��Ӧ�������ޣ�double nasĬ��ֵ+1.0e20��
	m_print.add(UALLOW);	//��Ӧ�������ޣ�double nasĬ��ֵ+1.0e20��
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	"DCONSTR",
	// 	Printer::int2str(DCID),			//���Լ����ID�ţ�int>0��
	// 	Printer::int2str(RID),			//�����Ӧ��DRESPi��Ƭ�������ID�ţ�int>0��
	// 	Printer::double2str(LALLOW, 8),	//��Ӧ�������ޣ�double nasĬ��ֵ+1.0e20��
	// 	Printer::double2str(UALLOW, 8),	//��Ӧ�������ޣ�double nasĬ��ֵ+1.0e20��
	// });
}

