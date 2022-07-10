#pragma once
#include <vector>
using namespace std;
//#include"armadillo"
#include"factorial.h"
//#include"CSTsurface.h"
#include<string>
//using namespace arma;

//面元法输入信息参数
class AeroInfo
{
public:
	AeroInfo();
	vector<int> CalcBlockID;
	vector<int> CalcBlockArea;
	vector<int> CalcBlockType;
	vector<int> CalcBlockOrient;
	vector<double> Ma;
	vector<double> Alpha;
	vector<double> Beta;
	vector<double> Hm;
	double sRef;//参考面积 m2
	double cRef;//参考长度 m
	string meshPath;//网格文件路径
	string outputPath; //输出路径（文件夹）
	string outputName;//输出结果文件名（包含后缀）

	int CalcBlockIDNum() { return CalcBlockID.size(); }
	int MaNum() { return Ma.size(); }
	int AlphaNum() { return Alpha.size(); }
	int BetaNum() { return Beta.size(); }
	int HmNum() { return Hm.size(); }

	int	ifViscous;					//0   是否启用粘性计算
	double HeadRadius;				//0.03头部半径
	double WallTemperature;			//300 壁面温度
	double WallEmissivity;			//0.7 壁面辐射率
	int AirType;					//1   气体模型（0，1）
	double HeatTransferCoefficient;	//0   换热系数
	double MediumTemperature;		//280 介质系数
	double InnerWallTemperature;	//350 内壁面温度
};
//面元法求解器
class AeroCalc
{
public:
	AeroCalc() { ai = AeroInfo(); CommandfileName.clear(); }
	int printCommand(string CommandfileName);//打印输入文件
	int CalcAero();//通过命令行调用气动求解程序
	void setAeroInfo(AeroInfo ai) { this->ai = ai; }//设置气动输入参数
	string getCommandfileName() { return CommandfileName; }
	int getAeroForce(map<int, Point>& nodeForce);

private:
	int checkData();
	AeroInfo ai;
	string CommandfileName;
};


//AVL计算输入参数
class AVLInfo
{
public:
	//AVL截面类
	class AVLSectionInfo
	{
	public:
		string sectionName{ "" };
		Point Orgin{ Point() };
		//double Orgin_x, Orgin_y, Orgin_z;
		double Scale{ 0 };
		double Angle{ 0 };
		vector<Point> pointList;//存储翼型文件的节点序列
		void printSectionFile(string path);
	};
	AVLInfo();
	string name;
	double Ma;
	double Alpha;
	double V;
	double rho;
	double _q;
	//double q() { return rho * V * V / 2.0; }
	double grav;
	double IYsym, IZsym, Zsym; //是否有对称面 0代表没有 1 代表有
	double Sref, Cref, Bref;   //参考面积 弦长 展长 翼型等线显示为白色虚线
	//double Xref, Yref, Zref;
	int Nchordwise, Nspanwise; //马蹄窝在展长和弦长上的数量（一侧机翼）
	double Angle;			   //所有截面角度的改变
	Point Scale;			   //尺寸的缩放
	Point Orgin;			   //位置的平移
	vector<AVLSectionInfo> sectionInfo;
};


//AVL求解翼面输出信息
class AVLres_surface
{
public:
	//AVL求解单元条信息
	struct AVLres_strip
	{
		//AVL求解单元信息
		struct AVLres_elem
		{
			int elemID{ 0 };

			Point site;
			double Dx{ 0 };
			double Slope{ 0 };
			double dCp{ 0 };
			double dCl() { return 1 / sqrt(1 + Slope * Slope) * dCp; }
			double dCd() { return abs(Slope) / sqrt(1 + Slope * Slope) * dCp; }
		};
		int StripID{ 0 };
		int Chordwise{ 0 };
		int FirstVortex{ 0 };

		Point Orgin;
		double Chord{ 0 };
		double Incidence{ 0 };
		double StripArea{ 0 };
		double StripWidth{ 0 };
		double cl{ 0 }, cd{ 0 };
		map<int,AVLres_elem> elemRes;
	};
	AVLres_surface() { CLsurf = CDsurf = Surface_area = 1e-6; SurfaceID = Chordwise = Spanwise = FirstStrip = 0; }
	int SurfaceID;
	string SurfName;
	int Chordwise;
	int Spanwise;
	int FirstStrip;

	double Surface_area;
	double CLsurf;
	double CDsurf;
	map<int,AVLres_strip> stripRes;
	AVLInfo ai;
};

class AVLServer
{
public:
	AVLServer() { CommandfileName = "none"; partname = "testPart"; };
	int printCommand(string CommandfileName);//打印输入文件
	int CalcAero();//通过命令行调用气动求解程序
	void setAVLexePath(string exepath = "");
	void setAVLInfo(const AVLInfo& ai) { this->ai = ai; }//设置气动输入参数
	//string getCommandfileName() { return CommandfileName; }
	int getAeroForce(map<int, AVLres_surface>& res);

private:
	//int checkData();
	string AVLexePath;
	string workpath;
	string partname;
	string CommandfileName;
	string ResFilePath;
	AVLInfo ai;

};

//XFoil求解参数
//求解参数文件
class XFoilInfo
{
	typedef tuple<double, double, double> tp;

public:
	XFoilInfo() = default;
	
	string foilName{"undefine"};//翼型名称
	vector<Point> foilList{};	//翼型散点形状参数
	//
	string workPath{""};		//工作文件夹路径
	string shapePath() { return workPath + "\\" + foilName + ".dat"; } //输入文件翼型参数的路径
	string inputPath() { return workPath + "\\" + foilName + "_input.dat"; }//输入文件翼型参数的路径
	string resultPath() { return workPath + "\\" + foilName + ".cp"; }//输入文件翼型参数的路径
	//
	double VISC{0};				//雷诺数定义
	double Ma{0};				//马赫数
	int ITER{50};				//迭代的步数
	int outputNum{160};			//输出的节点个数 默认160、最大300节点
	// int setAlpha(double first, double secend, int gap);
	// const tp &Alpha() { return alpha; }


	tp alpha{0,0,0};//记录攻角的起始范围和间隔

};

class XFoilSolver
{
public:
	XFoilSolver(XFoilInfo &info) : m_info(info){}
	int solve();
	string ExePath{"xfoil.exe"};
	const vector<Field<double>> getRes() { return m_res; }

private:
	int printShape();//打印翼型文件
	int printInput();//打印翼型文件
	int readRes();
	vector<Field<double>> m_res{};
	XFoilInfo &m_info;//输入参数
};