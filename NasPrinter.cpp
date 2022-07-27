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
		ofsHeader << "include '" << (name + "_mesh.bdf'") << endl;		//头文件内部使用相对路径
		ofsHeader << "include '" << (name + "_property.bdf'") << endl;	//头文件内部使用相对路径
		ofsHeader << "include '" << (name + "_force.bdf'") << endl;		//头文件内部使用相对路径
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
	ofs << setw(8) << words[0].substr(0, 8);//第一个关键字
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

	ofs << setw(8) << words[0].substr(0,8) + "*";//第一个关键字
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
	//ofs.str("");//清空缓冲区
	int listNum = 0;
	ofs << std::left;
	ofs << setw(8) << words[0].substr(0, 8);;//第一个关键字
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
	//ofs.str("");//清空缓冲区
	int listNum = 0;
	ofs << std::left;

	ofs << setw(8) << words[0].substr(0, 8) + "*";//第一个关键字
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
	//整型
	const int xx = (int)fabsdb;
	if ((fabsdb - xx) < 0.4 && abs(xx - fabsdb) < 1e-5)//四舍
	{
		return to_string((int)db) + ".";
	}
	if ((fabsdb - xx >= 0.4 && abs(xx - fabsdb + 1) < 1e-5))//五入
	{
		return to_string((int)db + sign) + ".";
	}
	//浮点型
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
	m_print.add("GRID");	//卡片名
	m_print.add(CID);		//坐标系编号
	m_print.add(CP);		//x坐标
	m_print.add(pt);		//z坐标
	AddCardLong(ssMesh, m_print.str_list);
	// AddCardLong(ssMesh, {
	// 	"GRID",						//卡片名
	// 	Printer::int2str(CID),				//节点编号 CID
	// 	Printer::int2str(CP),				//坐标系编号
	// 	Printer::double2str(pt.getX(),16),	//x坐标
	// 	Printer::double2str(pt.getY(),16),	//y坐标
	// 	Printer::double2str(pt.getZ(),16),	//z坐标
	// });
}

void NasPrinter::addCMASS2(C_INT EID, C_DOUBEL M, C_INT G1, C_INT C1)
{
	m_print.init(5, 8);
	m_print.add("CMASS2");	//卡片名
	m_print.add(EID);		//单元编号 EID
	m_print.add(M);		//单元的质量M
	m_print.add(G1);		//节点号 G1
	m_print.add(C1);		//分量号 C1
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CMASS2",					//卡片名
	// 	Printer::int2str(EID),				//单元编号 EID
	// 	Printer::double2str(M, 8),			//单元的质量M
	// 	Printer::int2str(G1),				//节点号 G1
	// 	Printer::int2str(C1),				//分量号 C1
	// });
}

void NasPrinter::addCBAR(C_INT EID, C_INT PID, C_INT GA, C_INT GB, C_DOUBEL X1, C_DOUBEL X2, C_DOUBEL X3)
{
	m_print.init(8, 8);
	m_print.add("CBAR");	//卡片名
	m_print.add(EID);		//单元编号 EID
	m_print.add(PID);		//单元属性编号 PID
	m_print.add(GA);		//节点1
	m_print.add(GB);		//节点2
	m_print.add(X1);		//一维单元方向指定
	m_print.add(X2);		//一维单元方向指定
	m_print.add(X3);		//一维单元方向指定
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CBAR",									//卡片名
	// 	Printer::int2str(EID),					//单元编号 EID
	// 	Printer::int2str(PID),					//单元属性编号 PID
	// 	Printer::int2str(GA),					//节点1
	// 	Printer::int2str(GB),					//节点2
	// 	Printer::double2str(X1, 8),				//一维单元方向指定
	// 	Printer::double2str(X2, 8),				//一维单元方向指定
	// 	Printer::double2str(X3, 8),				//一维单元方向指定
	// });
}

void NasPrinter::addCTRIA3(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3)
{
	m_print.init(6, 8);
	m_print.add("CTRIA3");	//卡片名
	m_print.add(EID);		//单元编号 EID
	m_print.add(PID);		//单元属性编号 PID
	m_print.add(G1);		//节点1 CID
	m_print.add(G2);		//节点2 CID
	m_print.add(G3);		//节点3 CID
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CTRIA3",				//卡片名
	// 	Printer::int2str(EID),	//单元编号 EID
	// 	Printer::int2str(PID),	//单元属性编号 PID
	// 	Printer::int2str(G1),	//节点1 CID
	// 	Printer::int2str(G2),	//节点2 CID
	// 	Printer::int2str(G3),	//节点3 CID
	// });
}

void NasPrinter::addCQUAD4(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3, C_INT G4)
{
	m_print.init(7, 8);
	m_print.add("CQUAD4");	//卡片名
	m_print.add(EID);		//单元编号 EID
	m_print.add(PID);		//单元属性编号 PID
	m_print.add(G1);		//节点1 CID
	m_print.add(G2);		//节点2 CID
	m_print.add(G3);		//节点3 CID
	m_print.add(G4);		//节点4 CID
	AddCardShort(ssMesh, m_print.str_list);
	// AddCardShort(ssMesh, {
	// 	"CQUAD4",				//卡片名
	// 	Printer::int2str(EID),	//单元编号 EID
	// 	Printer::int2str(PID),	//单元属性编号 PID
	// 	Printer::int2str(G1),	//节点1 CID
	// 	Printer::int2str(G2),	//节点2 CID
	// 	Printer::int2str(G3),	//节点3 CID
	// 	Printer::int2str(G4),	//节点4 CID
	// });
}

void NasPrinter::addRBE3(C_INT EID, C_INT REDGRID, C_INT REFC, C_DOUBEL WT1, C_INT C1, C_VECTOR(int) G1j)
{
	m_print.init(7 + G1j.size(), 8);
	m_print.add("RBE3");	//卡片名
	m_print.add(EID);		//数据卡编号
	m_print.add("");		//
	m_print.add(REDGRID);	//参考节点编号	
	m_print.add(REFC);		//参考节点自由度分量号
	m_print.add(WT1);		//节点G1j的位移量的加权值（比例系数）
	m_print.add(C1);		//节点G1j的位移分量号
	m_print.add(G1j);		//多个节点
	AddCardShort(ssMesh, m_print.str_list);
	// vector<string> strList = {
	// 	"RBE3",							//卡片名
	// 	Printer::int2str(EID),			//数据卡编号
	// 	"",								//
	// 	Printer::int2str(REDGRID),		//参考节点编号
	// 	Printer::int2str(REFC),			//参考节点自由度分量号
	// 	Printer::double2str(WT1, 8),	//节点G1j的位移量的加权值（比例系数）
	// 	Printer::int2str(C1),			//节点G1j的位移分量号
	// };
	// strList.reserve(strList.size() + G1j.size());
	// for (size_t j = 0; j < G1j.size(); j++)
	// {
	// 	strList.push_back(Printer::int2str(G1j[j]));//多个节点
	// }
	// AddCardShort(ssMesh, strList);
}

void NasPrinter::addPLOAD(C_INT SID, C_DOUBEL P, C_INT G1, C_INT G2, C_INT G3, C_INT G4)
{
	m_print.init(7, 8);
	m_print.add("PLOAD");	//卡片名
	m_print.add(SID);		//载荷编号 SID
	m_print.add(P);		//压力值
	m_print.add(G1);		//节点1
	m_print.add(G2);		//节点2
	m_print.add(G3);		//节点3
	m_print.add_usize(G4);	//节点4 取默认值-1时默认输出为空
	AddCardShort(ssForce, m_print.str_list);
	// AddCardShort(ssForce, {
	// 	"PLOAD",					//卡片名
	// 	Printer::int2str(SID),		//载荷编号 SID
	// 	Printer::double2str(P, 8),	//压力值
	// 	Printer::int2str(G1),		//节点1
	// 	Printer::int2str(G2),		//节点2
	// 	Printer::int2str(G3),		//节点3
	// 	Printer::int2str(G4),		//节点4 取默认值-1时默认输出为空
	// });
}

void NasPrinter::addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3)
{
	m_print.init(8, 16);
	m_print.add("FORCE");	//卡片名
	m_print.add(SID);		//SID
	m_print.add(G);		//节点编号
	m_print.add(CID);		//局部坐标系编号 默认可不填=0（全局坐标系）
	m_print.add(F);		//集中力的幅值
	m_print.add(N1);		//集中力在坐标系中的矢量方向
	m_print.add(N2);		//集中力在坐标系中的矢量方向
	m_print.add(N3);		//集中力在坐标系中的矢量方向
	AddCardLong(ssForce, m_print.str_list);
	// AddCardLong(ssForce, {
	// 	"FORCE",						//卡片名
	// 	Printer::int2str(SID),			//SID
	// 	Printer::int2str(G),			//节点编号
	// 	Printer::int2str(CID),			//局部坐标系编号 默认可不填=0（全局坐标系）
	// 	Printer::double2str(F , 16),	//集中力的幅值
	// 	Printer::double2str(N1,16),		//集中力在坐标系中的矢量方向
	// 	Printer::double2str(N2,16),		//集中力在坐标系中的矢量方向
	// 	Printer::double2str(N3,16),		//集中力在坐标系中的矢量方向
	// });
}

void NasPrinter::addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, const Point & pt)
{
	addFORCE(SID, G, CID, F, pt.getX(), pt.getY(), pt.getZ());
}

void NasPrinter::addGRAV(C_INT SID, C_INT CID, C_DOUBEL A, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3)
{
	m_print.init(7, 16);
	m_print.add("GRAV");	//卡片名
	m_print.add(SID);		//载荷编号 SID
	m_print.add(CID);		//坐标系编号
	m_print.add(A);		//加速度幅值
	m_print.add(N1);		//加速度方向坐标1
	m_print.add(N2);		//加速度方向坐标2
	m_print.add(N3);		//加速度方向坐标3
	AddCardLong(ssForce, m_print.str_list);
	// AddCardLong(ssForce, {
	// 	"GRAV",							//卡片名
	// 	Printer::int2str(SID),			//载荷编号 SID
	// 	Printer::int2str(CID),			//坐标系编号
	// 	Printer::double2str(A, 8),		//加速度幅值
	// 	Printer::double2str(N1, 8),		//加速度方向坐标1
	// 	Printer::double2str(N2, 8),		//加速度方向坐标2
	// 	Printer::double2str(N3, 8),		//加速度方向坐标3
	// });
}

void NasPrinter::addLOAD(C_INT SID, C_DOUBEL S, C_VECTOR(DoubleIntPair) LOADi)
{
	m_print.init(3 + LOADi.size() * 2, 8);
	m_print.add("LOAD");	//卡片名
	m_print.add(SID);		//数据卡编号
	m_print.add(S);		//总系数
	m_print.add(LOADi);	//单个载荷的系数, 单个载荷数据卡编号
	AddCardShort(ssForce, m_print.str_list);
	// vector<string> strList = {
	// 	"LOAD",						//卡片名
	// 	Printer::int2str(SID),		//数据卡编号
	// 	Printer::double2str(S, 8),	//总系数
	// };
	// strList.reserve(strList.size() + LOADi.size() * 2);
	// for (size_t i = 0; i < LOADi.size(); i++)
	// {
	// 	strList.push_back(Printer::double2str(LOADi[i].first, 8));	//单个载荷的系数
	// 	strList.push_back(Printer::int2str(LOADi[i].second));			//单个载荷数据卡编号
	// }
	// AddCardShort(ssForce, strList);
}

void NasPrinter::addSPC1(C_INT SID, C_INT C, C_VECTOR(int) Gi)
{
	m_print.init(3 + Gi.size(), 8);
	m_print.add("SPC1");	//卡片名
	m_print.add(SID);		//数据卡编号
	m_print.add(C);		//总系数
	m_print.add(Gi);		//单个载荷的系数
	AddCardShort(ssForce, m_print.str_list);
	// vector<string> strList = {
	// 	"SPC1",						//卡片名
	// 	Printer::int2str(SID),		//数据卡编号
	// 	Printer::int2str(C),		//总系数
	// };
	// strList.reserve(strList.size() + Gi.size());
	// for (size_t i = 0; i < Gi.size(); i++)
	// {
	// 	strList.push_back(Printer::int2str(Gi[i]));	//单个载荷的系数
	// }
	// AddCardShort(ssForce, strList);
}

void NasPrinter::addPSHELL(const PSHELL & ps)
{
	ssProperty << "$-- Property " << ps.name << " --$" << endl;
	m_print.init(4, 8);
	m_print.add("PSHELL");	//卡片名
	m_print.add(ps.PID);	//单元属性编号 PID
	m_print.add(ps.MID);	//材料属性编号 MID
	m_print.add(ps.T);		//厚度
	AddCardShort(ssProperty, m_print.str_list);
	// AddCardShort(ssProperty, {
	// 	"PSHELL",						//卡片名
	// 	Printer::int2str(ps.PID),		//单元属性编号 PID
	// 	Printer::int2str(ps.MID),		//材料属性编号 MID
	// 	Printer::double2str(ps.T, 8),	//厚度
	// });
}

void NasPrinter::addPBARL(const PBARL & pb)
{
	ssProperty << "$-- Property " << pb.name << " --$" << endl;
	m_print.init(9 + pb.DIM.size(), 8);
	m_print.add("PBARL");	//卡片名
	m_print.add(pb.PID);	//单元属性编号 PID
	m_print.add(pb.MID);	//材料属性编号 MID
	m_print.add("");		//
	m_print.add(pb.TYPE);	//单元截面属性类型
	m_print.add("");		//
	m_print.add("");		//
	m_print.add("");		//
	m_print.add("");		//
	m_print.add(pb.DIM);	//截面属性参数
	AddCardShort(ssProperty, m_print.str_list);
	// vector<string> strList = {
	// 	"PBARL",			//卡片名
	// 	Printer::int2str(pb.PID),	//材料属性编号 PID
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
	//Lanczos法提取模态
	m_print.init(5, 8);
	m_print.add("EIGRL");	//卡片名
	m_print.add(SID);		//集编号（整数>0）
	m_print.add("");		//
	m_print.add(ND);		//Number of Desired Roots期望根数（整数>0,或空格）
	m_print.add(MSGLVL);	//Diagnostic Output Level对话层信息，其中0无信息，1每次平移。打印获得的特征值，>1增加对话输出层
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	Printer::int2str(SID),		//集编号（整数>0）
	// 	"","",
	// 	Printer::int2str(ND),		//Number of Desired Roots期望根数（整数>0,或空格）
	// 	Printer::int2str(MSGLVL),	//Diagnostic Output Level对话层信息，其中0无信息，1每次平移。打印获得的特征值，>1增加对话输出层
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
	//后续可以提供一个数值检测的集制
	//定义设计变量的卡片
	m_print.init(7, 8);
	m_print.add("DESVAR");	//卡片名
	m_print.add(ID);		//设计变量的ID号（int>0）
	m_print.add(LABEL);	//用户指定的设计变量名字
	m_print.add(XINIT);	//设计变量的初始值（double XLB<=XINIT<=XUB）
	m_print.add(XLB);		//设计变量的下限（double nas默认值-1.0e20）
	m_print.add(XUB);		//设计变量的上限（double nas默认值+1.0e20）
	m_print.add(DELXV);	//设计变量的该变化量与初始值的比值的最大值（double>0）
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	"DESVAR",						//卡片名
	// 	Printer::int2str(ID),			//设计变量的ID号（int>0）
	// 	LABEL,							//用户指定的设计变量名字
	// 	Printer::double2str(XINIT,8),	//设计变量的初始值（double XLB<=XINIT<=XUB）
	// 	Printer::double2str(XLB, 8),	//设计变量的下限（double nas默认值-1.0e20）
	// 	Printer::double2str(XUB, 8),	//设计变量的上限（double nas默认值+1.0e20）
	// 	Printer::double2str(DELXV,8),	//设计变量的该变化量与初始值的比值的最大值（double>0）
	// });
}

void NasPrinter::addDVPREL1(C_INT ID, C_STR TYPE, C_INT PID, C_STR PNAME, C_VECTOR(DoubleIntPair) DVID_COEFi, C_DOUBEL C0, C_DOUBEL PMIN, C_DOUBEL PMAX)
{
	m_print.init(9 + DVID_COEFi.size() * 2, 8);
	m_print.add("DVPREL1");	//卡片名
	m_print.add(ID);			//该设计变量连接的ID号（int>0）
	m_print.add(TYPE);			//属性类型的名字，如“PBAR”，“PBEAM”等（string）
	m_print.add(PID);			//属性的ID号（int>0）
	m_print.add(PNAME);		//属性名字，如“T”、“A”，或属性卡片字域（Field）的未知或者单元属性表中字符（Word）的位置（sring或int>0）
	m_print.add(PMIN);			//该属性允许的最小值（double，nas默认值-1.0e20）
	m_print.add(PMAX);			//该属性允许的最大值（double，nas默认值+1.0e20）
	m_print.add(C0);			//线性关系的常数项（double，nas默认值为0.0）
	m_print.add("");			//
	m_print.add(DVID_COEFi);	//DVIDi DESVAR卡片的ID号（int>0）, COEFi 线性关系的一次项系数（double）
	AddCardShort(ssHeader, m_print.str_list);
	// vector<string> strList = {
	// 	"DVPREL1",
	// 	Printer::int2str(ID),			//该设计变量连接的ID号（int>0）
	// 	TYPE,							//属性类型的名字，如“PBAR”，“PBEAM”等（string）
	// 	Printer::int2str(PID),			//属性的ID号（int>0）
	// 	PNAME,							//属性名字，如“T”、“A”，或属性卡片字域（Field）的未知或者单元属性表中字符（Word）的位置（sring或int>0）
	// 	Printer::double2str(PMIN, 8),	//该属性允许的最小值（double，nas默认值-1.0e20）
	// 	Printer::double2str(PMAX, 8),	//该属性允许的最大值（double，nas默认值+1.0e20）
	// 	Printer::double2str(C0, 8),		//线性关系的常数项（double，nas默认值为0.0）
	// 	"",
	// };
	// strList.reserve(strList.size() + DVID_COEFi.size());
	// for (size_t i = 0; i < DVID_COEFi.size(); i++)
	// {
	// 	strList.push_back(Printer::int2str(DVID_COEFi[i].second));//DVIDi DESVAR卡片的ID号（int>0）
	// 	strList.push_back(Printer::double2str(DVID_COEFi[i].first, 8));//COEFi 线性关系的一次项系数（double）
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
	//定义约束
	m_print.init(5, 8);
	m_print.add("DCONSTR");//卡片名
	m_print.add(DCID);		//设计约束的ID号（int>0）
	m_print.add(RID);		//设计响应（DRESPi卡片）定义的ID号（int>0）
	m_print.add(LALLOW);	//响应量的下限（double nas默认值+1.0e20）
	m_print.add(UALLOW);	//响应量的上限（double nas默认值+1.0e20）
	AddCardShort(ssHeader, m_print.str_list);
	// AddCardShort(ssHeader, {
	// 	"DCONSTR",
	// 	Printer::int2str(DCID),			//设计约束的ID号（int>0）
	// 	Printer::int2str(RID),			//设计响应（DRESPi卡片）定义的ID号（int>0）
	// 	Printer::double2str(LALLOW, 8),	//响应量的下限（double nas默认值+1.0e20）
	// 	Printer::double2str(UALLOW, 8),	//响应量的上限（double nas默认值+1.0e20）
	// });
}

