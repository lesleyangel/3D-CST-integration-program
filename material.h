#pragma once
#include"include/armadillo"
using namespace arma;
using namespace std;

//����Ԫ�������Կ�Ƭ
class PBARL
{
public:
	PBARL() { PID = 0; MID = -1; }

	string name;
	int PID;
	int MID;//MID���û�и�ֵ���ӡʱĬ�����PID;
	string TYPE;
	vector<double> DIM;
	void printPBARL(ofstream& ofs);
};

//�ǵ�Ԫ���Կ�Ƭ
class PSHELL
{
public:
	PSHELL() { PID = 0; MID = -1; T = 0.5; }

	string name;
	int PID;
	int MID;//MID���û�и�ֵ���ӡʱĬ�����PID;
	double T;//���
	void printPSHELL(ofstream& ofs);
};

//�������Կ�Ƭ
class MAT1
{
public:
	MAT1() { MID = 0; E = 2.1e11; Nu = 0.3; RHO = 2700; }

	string name;
	int MID;
	double E;
	double Nu;
	double RHO;
	void printMAT1(ofstream& ofs);
};

//��ά�����������Բ���
class MAT8
{
public:
	MAT8() { MID = 0; E1 = E2 = 2.1e+011; NU12 = 0.3; G12 = G1Z = G2Z = 1.5e6; RHO = 2.7e+003; }
	std::string name;
	int MID;
	double E1;
	double E2;
	double NU12;
	double G12;
	double G1Z;
	double G2Z;
	double RHO;
	//void printMAT8(std::ofstream& ofs);
};

class Property
{
public:
	Property() { pshell_list.clear(); pbarl_list.clear(); mat1_list.clear(); }
	void readFromFile(ifstream& ifs);//��CST�����ļ��ڶ�ȡProperty����
	bool isNotEmpty();
	void SetPSHELLList(map<int, PSHELL> pshell_list);// { this->pshell_list = pshell_list; }
	void SetPBARLList(map<int, PBARL> pbarl_list);// { this->pbarl_list = pbarl_list; }
	void SetMAT1List(map<int, MAT1> mat1_list);// { this->mat1_list = mat1_list; }
	int add_PSHRLL(PSHELL &&ps);
	const map<int, PSHELL> &getPSHELLlist() { return pshell_list; }
	const map<int, MAT1> &getMAT1List() { return mat1_list; }
	PSHELL getPSHELL(int PID);
	PBARL getPBARL(int PID);
	MAT1 getMAT1(int MID);
	void printPSHELL_all(ofstream& ofs);
	void printPBARL_all(ofstream& ofs);
	void printMAT1_all(ofstream& ofs);
private:
	map<int, PSHELL> pshell_list;
	map<int, PBARL> pbarl_list;
	map<int, MAT1> mat1_list;
};
inline bool Property::isNotEmpty()
{
	if (pshell_list.size() == 0 && pbarl_list.size() == 0 && mat1_list.size() == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}