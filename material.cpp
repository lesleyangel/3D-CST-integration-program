#include"material.h"
#include <iomanip>

void PBARL::printPBARL(ofstream& ofs)
{
	MID = (MID < 0) ? PID : MID;
	//ofs <<  setprecision(3);//设置输出精度为8
	ofs << "$-- Property " << name << " --$" << endl;
	ofs << "PBARL   " << setw(8) << PID << setw(8) << PID << setw(16) << TYPE << endl;
	ofs << "        ";
	for (size_t i = 0; i < DIM.size(); i++)
	{
		ofs << setw(8) << DIM[i];
	} 
	ofs << endl;
}

void PSHELL::printPSHELL(ofstream& ofs)
{
	MID = (MID < 0) ? PID : MID;
	//ofs << setprecision(3);//设置输出精度为8
	ofs << "$-- Property " << name << " --$" << endl;
	ofs << "PSHELL  " << setw(8) << PID << setw(8) << MID << setw(8) << T
		<< setw(8) << 1 << setw(16) << 1 << endl;
}

void MAT1::printMAT1(ofstream& ofs)
{
	//ofs <<  setprecision(3);//设置输出精度为8
	//ofs << "$-- Material " << name << " --$" << endl;
	ofs << "MAT1    " << setw(8) << MID << setw(8) << E << setw(16) << Nu
		<< setw(8) << RHO << endl;
}

//PSHELL Property::getPSHELL(int PID)
//{
//	if (0 <= PID && PID < pshell_List.size())
//	{
//		return pshell_List[PID];
//	}
//	else
//	{
//		cout << "梁截面属性编号已超出边界！" << endl;
//		return PSHELL();
//	}
//}
//
//PBARL Property::getPBARL(int PID)
//{
//	if (0 <= PID && PID < pbarl_List.size())
//	{
//		return pbarl_List[PID];
//	}
//	else
//	{
//		cout << "壳单元属性编号已超出边界！" << endl;
//		return PBARL();
//	}
//}
//
//MAT1 Property::getMAT1(int MID)
//{
//	if (0 <= MID && MID < mat1_List.size())
//	{
//		return mat1_List[MID];
//	}
//	else
//	{
//		cout << "梁截面属性编号已超出边界！" << endl;
//		return MAT1();
//	}
//}

//
void Property::SetMAT1List(map<int, MAT1> m_list)
{
	mat1_list.clear();
	for (map<int,MAT1>::iterator it = m_list.begin(); it != m_list.end(); it++)
	{
		it->second.MID = it->first;
	}
	this->mat1_list = m_list;
}

void Property::SetPBARLList(map<int, PBARL> p_list)
{
	pbarl_list.clear();
	for (map<int, PBARL>::iterator it = p_list.begin(); it != p_list.end(); it++)
	{
		it->second.PID = it->first;
	}
	this->pbarl_list = p_list;
}

int_fast16_t Property::add_PSHRLL(PSHELL&& ps)
{
	map<int, PSHELL>::iterator it = pshell_list.find(ps.PID);
	if (it != pshell_list.end())
	{
		cout << "Property::add_PSHRLL(PSHELL&& ps) 输入的id属性已经存在！" << endl;
		return -1;
	}
	pshell_list.insert(pair<int, PSHELL>(ps.PID, ps));
	return 0;
}

void Property::readFromFile(ifstream& ifs)
{
	map<int, PSHELL> ps_list;
	map<int, MAT1> m1_list;
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//名称清除空格行
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据

		getline(ss, str_line, ',');//获取首个string
		if (str_line[0] == '$') {}//跳过
		else if (!str_line.compare("PSHELL"))
		{
			pair<int, PSHELL> temp;
			stringstream  iss;
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.first; iss.clear();
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.second.MID; iss.clear();
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.second.T; iss.clear();
			ps_list.insert(temp);
		}
		else if (!str_line.compare("MAT1"))
		{
			pair<int, MAT1> temp;
			stringstream  iss;
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.first; iss.clear();
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.second.E; iss.clear();
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.second.Nu; iss.clear();
			getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> temp.second.RHO; iss.clear();
			m1_list.insert(temp);
		}
	}
	SetPSHELLList(ps_list);
	SetMAT1List(m1_list);
}

//从CST输入文件内读取Property属性



void Property::SetPSHELLList(map<int, PSHELL> p_list)
{ 
	pshell_list.clear();
	for (map<int, PSHELL>::iterator it = p_list.begin(); it != p_list.end(); it++)
	{
		it->second.PID = it->first;
	}
	this->pshell_list = p_list;
}

MAT1 Property::getMAT1(int MID)
{
	map<int, MAT1>::iterator it = mat1_list.find(MID);
	if (it != mat1_list.end())
	{
		return it->second;
	}
	else
	{
		cout << "MAT1已越界 已自动补齐新的MID！请检查mat1_list" << endl;
		MAT1 temp;
		temp.MID = MID;
		mat1_list.insert(pair<int, MAT1>(MID, temp));
		return temp;
	}
	//return MAT1();
}

void Property::printPSHELL_all(ofstream& ofs)
{
	for (map<int, PSHELL>::iterator it = pshell_list.begin(); it != pshell_list.end(); it++)
	{
		it->second.printPSHELL(ofs);
	}
}
void Property::printPBARL_all(ofstream& ofs)
{
	for (map<int, PBARL>::iterator it = pbarl_list.begin(); it != pbarl_list.end(); it++)
	{
		it->second.printPBARL(ofs);
	}
}
void Property::printMAT1_all(ofstream& ofs)
{
	for (map<int, MAT1>::iterator it = mat1_list.begin(); it != mat1_list.end(); it++)
	{
		it->second.printMAT1(ofs);
	}
}

PBARL Property::getPBARL(int PID)
{
	map<int, PBARL>::iterator it = pbarl_list.find(PID);
	if (it != pbarl_list.end())
	{
		return it->second;
	}
	else
	{
		cout << "PBARL已越界 已自动补齐新的PID！请检查pbarl_list" << endl;
		PBARL temp;
		temp.PID = PID;
		temp.MID = 1;
		pbarl_list.insert(pair<int, PBARL>(PID, temp));
		return temp;
	}
	return PBARL();
}

PSHELL Property::getPSHELL(int PID)
{
	map<int, PSHELL>::iterator it = pshell_list.find(PID);
	if (it != pshell_list.end())
	{
		return it->second;
	}
	else
	{
		cout << "PSHELL已越界 已自动补齐新的PID = " << PID << "！请检查pshell_list" << endl;
		PSHELL temp;
		temp.PID = PID;
		temp.MID = 1;
		pshell_list.insert(pair<int, PSHELL>(PID, temp));
		return temp;
	}
	return PSHELL();
}