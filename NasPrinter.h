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

enum class OutputStyle
{
	onefile, multifile
};
class NasPrinter
{
	typedef std::pair<double, int> DoubleIntPair;
public:
	NasPrinter();
	int PrintBDF(const std::string& path, const std::string& name, OutputStyle type = OutputStyle::onefile);

	void AddCardShort(std::ofstream& ofs, const std::vector<std::string>& words);
	void AddCardLong(std::ofstream& ofs, const std::vector<std::string>& words);
	void AddCardShort(std::stringstream& ss, const std::vector<std::string>& words);//比直接ofstream快
	void AddCardLong(std::stringstream& ss, const std::vector<std::string>& words);//比直接ofstream快
	static std::string double2str(const double& db, int size);
	static std::string int2str(const int& n) { return std::to_string(n); }
	
	std::stringstream ssHeader;
	std::stringstream ssForce;
	std::stringstream ssMesh;
	std::stringstream ssProperty;
	//几何模型
	void addGRID(C_INT CID, C_INT CP, const arith::Point& pt);//添加节点卡片
	void addCMASS2(C_INT EID, C_DOUBEL M, C_INT G1, C_INT C1);//添加集中质量点 零维单元
	void addCBAR(C_INT EID, C_INT PID, C_INT GA, C_INT GB, C_DOUBEL X1, C_DOUBEL X2, C_DOUBEL X3);//添加杆单元 一维单元
	void addCTRIA3(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3);//添加三角形单元 二维单元
	void addCQUAD4(C_INT EID, C_INT PID, C_INT G1, C_INT G2, C_INT G3, C_INT G4);//添加四边形单元 二维单元
	void addRBE3(C_INT EID, C_INT REDGRID, C_INT REFC, C_DOUBEL WT1, C_INT C1, C_VECTOR(int) G1j);//R单元
	//约束和载荷
	void addPLOAD(C_INT SID, C_DOUBEL P, C_INT G1, C_INT G2, C_INT G3, C_INT G4 = -1);//添加压力
	void addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3);//添加节点集中力
	void addFORCE(C_INT SID, C_INT G, C_INT CID, C_DOUBEL F, const arith::Point& pt);//添加节点集中力
	void addGRAV(C_INT SID, C_INT CID, C_DOUBEL A, C_DOUBEL N1, C_DOUBEL N2, C_DOUBEL N3);//添加惯性载荷
	void addLOAD(C_INT SID, C_DOUBEL S, C_VECTOR(DoubleIntPair) LOADi);//添加静态组合载荷
	//材料及属性
	void addPSHELL(const PSHELL& ps);
	void addPBARL(const PBARL& pb);
	void addMAT8(const MAT8& mt);
	void addMAT1(const MAT1& mt);
	//模态分析
	void addEIGRL(C_INT SID, C_INT ND, C_INT MSGLVL);//添加模态提取卡片

	//优化分析
	void addDOPTPRM(C_VECTOR(std::string)command);
	//设计变量卡片
	void addDESVAR(C_INT ID, C_STR LABEL, C_DOUBEL XINIT, C_DOUBEL XLB, C_DOUBEL XUB, C_DOUBEL DELXV);//添加设计变量卡片
	void addDVPREL1(C_INT ID, C_STR TYPE, C_INT PID, C_STR PNAME, C_VECTOR(DoubleIntPair) DVID_COEFi, C_DOUBEL C0=0, C_DOUBEL PMIN=-1e20, C_DOUBEL PMAX=1e20);
	//响应卡片
	void addDRESP1(C_INT ID, C_STR LABLE, C_STR RTYPE, C_STR PTYPE, C_INT REGION, C_INT ATTA, C_INT ATTB, C_VECTOR(int)ATTi);//定义Ⅰ阶相应
	void addDRESP1_WEIGHT(C_INT ID, C_STR LABEL);//定义Ⅰ阶相应 重量响应
	void addDRESP1_FREQ(C_INT ID, C_STR LABEL, C_INT ATTA);//定义Ⅰ阶相应 正交模态响应 输入ATTA=模态响应阶数
	void addDRESP1_STRESS(C_INT ID, C_STR LABEL, C_STR PTYPE, C_INT ATTA, C_VECTOR(int)ATTi);//定义Ⅰ阶相应 应力响应 输入ATTA=应力应变代码；ATTi=包含单元编号
	void addDCONSTR(C_INT DCID, C_INT RID, C_DOUBEL LALLOW, C_DOUBEL UALLOW);//约束卡片
private:
	std::string m_path;
	std::string m_Nmae;
	std::string positiveInt2str(const int& n) { return (n < 0) ? "" : std::to_string(n); }//如果输入为负数则默认输出空

};
#undef C_INT
#undef C_DOUBEL
#undef C_STR
#undef C_VECTOR
