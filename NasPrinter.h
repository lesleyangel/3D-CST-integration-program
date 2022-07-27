#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <math.h>
#include "material.h"
#include "factorial.h"

#define C_INT const int&
#define C_DOUBEL const double&
#define C_STR const std::string&
#define C_VECTOR(type) const std::vector<type>&


class Printer
{
public:
	void init(int size, int precision);
	// template <class T>
	// std::string to_str(const T &sth);
	template <class T>
	void add(const T &sth);
	void add_usize(const int &usz);
	template <class T>
	void add(const std::vector<T> &vec);
	template <class T1,class T2>
	void add(const std::vector<pair<T1,T2>> &vec_pair);
	static std::string double2str(const double& sth, int precision);
	static std::string int2str(const int& n) { return std::to_string(n); }

	std::vector<std::string> str_list;
	int m_precision;
};

template <class T>
void Printer::add(const T &chr)
{
	// cout << "Printer::to_str(const T &sth) �����˴�������ͣ�" << endl;
	str_list.push_back(chr);
}
template <>
inline void Printer::add(const int &inn)
{
	str_list.push_back(std::to_string(inn));
}
template <>
inline void Printer::add(const double &dou)
{
	str_list.push_back(double2str(dou, m_precision));
}
template <>
inline void Printer::add(const std::string &str)
{
	str_list.push_back(str);
}
template <>
inline void Printer::add(const Point &pot)
{
	add(pot.getX());
	add(pot.getY());
	add(pot.getZ());
}
template <class T>
inline void Printer::add(const std::vector<T> &vec)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		add(vec[i]);
	}
}
template <class T1,class T2>
inline void Printer::add(const std::vector<pair<T1,T2>> &par)
{
	for (size_t i = 0; i < par.size(); i++)
	{
		add(par[i].first);
		add(par[i].second);
	}
}
//
class NasPrinter
{
	typedef std::pair<double, int> DoubleIntPair;
public:
	enum OutputStyle
	{//�����ļ���style
		onefile,   //���ļ���дbdf
		multifile, //���ļ���дbdf
	};
	NasPrinter();
	int PrintBDF(const std::string& path, const std::string& name, OutputStyle type = OutputStyle::onefile);

	void AddCardShort(std::ofstream& ofs, const std::vector<std::string>& words);
	void AddCardLong(std::ofstream& ofs, const std::vector<std::string>& words);
	void AddCardShort(std::stringstream& ss, const std::vector<std::string>& words);//��ֱ��ofstream��
	void AddCardLong(std::stringstream& ss, const std::vector<std::string>& words);//��ֱ��ofstream��

	
	std::stringstream ssHeader;
	std::stringstream ssForce;
	std::stringstream ssMesh;
	std::stringstream ssProperty;
	//����ģ��
	void addGRID(C_INT CID, C_INT CP, const Point& pt);//��ӽڵ㿨Ƭ
	void addCMASS2(C_INT EID, C_DOUBEL M, C_INT G1, C_INT C1);//��Ӽ��������� ��ά��Ԫ
	void addCBAR(C_INT EID, C_INT PID, C_INT GA, C_INT GB, C_DOUBEL X1, C_DOUBEL X2, C_DOUBEL X3);//��Ӹ˵�Ԫ һά��Ԫ
	void addCTRIA3(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3);//��������ε�Ԫ ��ά��Ԫ
	void addCQUAD4(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3, C_INT G4);//����ı��ε�Ԫ ��ά��Ԫ
	void addRBE3(C_INT EID, C_INT REDGRID, C_INT REFC, C_DOUBEL WT1, C_INT C1, C_VECTOR(int) G1j);//R��Ԫ
	//Լ�����غ�
	void addPLOAD(C_INT SID, C_DOUBEL P, C_INT G1, C_INT G2, C_INT G3, C_INT G4 = -1);//���ѹ��
	void addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3);//��ӽڵ㼯����
	void addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, const Point& pt);//��ӽڵ㼯����
	void addGRAV(C_INT SID, C_INT CID, C_DOUBEL A, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3);//��ӹ����غ�
	void addLOAD(C_INT SID, C_DOUBEL S, C_VECTOR(DoubleIntPair) LOADi);//��Ӿ�̬����غ�

	void addSPC1(C_INT SID, C_INT C, C_VECTOR(int) Gi); // ��ӵ���Լ��
	//���ϼ�����
	void addPSHELL(const PSHELL& ps);
	void addPBARL(const PBARL& pb);
	void addMAT8(const MAT8& mt);
	void addMAT1(const MAT1& mt);
	//ģ̬����
	void addEIGRL(C_INT SID, C_INT ND, C_INT MSGLVL);//���ģ̬��ȡ��Ƭ

	//�Ż�����
	void addDOPTPRM(C_VECTOR(std::string)command);
	//��Ʊ�����Ƭ
	void addDESVAR(C_INT ID, C_STR LABEL, C_DOUBEL XINIT, C_DOUBEL XLB, C_DOUBEL XUB, C_DOUBEL DELXV);//�����Ʊ�����Ƭ
	void addDVPREL1(C_INT ID, C_STR TYPE, C_INT PID, C_STR PNAME, C_VECTOR(DoubleIntPair) DVID_COEFi, C_DOUBEL C0=0, C_DOUBEL PMIN=-1e20, C_DOUBEL PMAX=1e20);
	//��Ӧ��Ƭ
	void addDRESP1(C_INT ID, C_STR LABLE, C_STR RTYPE, C_STR PTYPE, C_INT REGION, C_INT ATTA, C_INT ATTB, C_VECTOR(int)ATTi);//��������Ӧ
	void addDRESP1_WEIGHT(C_INT ID, C_STR LABEL);//��������Ӧ ������Ӧ
	void addDRESP1_FREQ(C_INT ID, C_STR LABEL, C_INT ATTA);//��������Ӧ ����ģ̬��Ӧ ����ATTA=ģ̬��Ӧ����
	void addDRESP1_STRESS(C_INT ID, C_STR LABEL, C_STR PTYPE, C_INT ATTA, C_VECTOR(int)ATTi);//��������Ӧ Ӧ����Ӧ ����ATTA=Ӧ��Ӧ����룻ATTi=������Ԫ���
	void addDCONSTR(C_INT DCID, C_INT RID, C_DOUBEL LALLOW, C_DOUBEL UALLOW);//Լ����Ƭ
private:
	std::string m_path;
	std::string m_Nmae;
	std::string positiveInt2str(const int& n) { return (n < 0) ? "" : std::to_string(n); }//�������Ϊ������Ĭ�������
	Printer m_print;
};
#undef C_INT
#undef C_DOUBEL
#undef C_STR
#undef C_VECTOR
