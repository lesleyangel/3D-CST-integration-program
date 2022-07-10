#pragma once
#include <vector>
using namespace std;
//#include"armadillo"
#include"factorial.h"
//#include"CSTsurface.h"
#include<string>
//using namespace arma;

//��Ԫ��������Ϣ����
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
	double sRef;//�ο���� m2
	double cRef;//�ο����� m
	string meshPath;//�����ļ�·��
	string outputPath; //���·�����ļ��У�
	string outputName;//�������ļ�����������׺��

	int CalcBlockIDNum() { return CalcBlockID.size(); }
	int MaNum() { return Ma.size(); }
	int AlphaNum() { return Alpha.size(); }
	int BetaNum() { return Beta.size(); }
	int HmNum() { return Hm.size(); }

	int	ifViscous;					//0   �Ƿ�����ճ�Լ���
	double HeadRadius;				//0.03ͷ���뾶
	double WallTemperature;			//300 �����¶�
	double WallEmissivity;			//0.7 ���������
	int AirType;					//1   ����ģ�ͣ�0��1��
	double HeatTransferCoefficient;	//0   ����ϵ��
	double MediumTemperature;		//280 ����ϵ��
	double InnerWallTemperature;	//350 �ڱ����¶�
};
//��Ԫ�������
class AeroCalc
{
public:
	AeroCalc() { ai = AeroInfo(); CommandfileName.clear(); }
	int printCommand(string CommandfileName);//��ӡ�����ļ�
	int CalcAero();//ͨ�������е�������������
	void setAeroInfo(AeroInfo ai) { this->ai = ai; }//���������������
	string getCommandfileName() { return CommandfileName; }
	int getAeroForce(map<int, Point>& nodeForce);

private:
	int checkData();
	AeroInfo ai;
	string CommandfileName;
};


//AVL�����������
class AVLInfo
{
public:
	//AVL������
	class AVLSectionInfo
	{
	public:
		string sectionName{ "" };
		Point Orgin{ Point() };
		//double Orgin_x, Orgin_y, Orgin_z;
		double Scale{ 0 };
		double Angle{ 0 };
		vector<Point> pointList;//�洢�����ļ��Ľڵ�����
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
	double IYsym, IZsym, Zsym; //�Ƿ��жԳ��� 0����û�� 1 ������
	double Sref, Cref, Bref;   //�ο���� �ҳ� չ�� ���͵�����ʾΪ��ɫ����
	//double Xref, Yref, Zref;
	int Nchordwise, Nspanwise; //��������չ�����ҳ��ϵ�������һ�����
	double Angle;			   //���н���Ƕȵĸı�
	Point Scale;			   //�ߴ������
	Point Orgin;			   //λ�õ�ƽ��
	vector<AVLSectionInfo> sectionInfo;
};


//AVL������������Ϣ
class AVLres_surface
{
public:
	//AVL��ⵥԪ����Ϣ
	struct AVLres_strip
	{
		//AVL��ⵥԪ��Ϣ
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
	int printCommand(string CommandfileName);//��ӡ�����ļ�
	int CalcAero();//ͨ�������е�������������
	void setAVLexePath(string exepath = "");
	void setAVLInfo(const AVLInfo& ai) { this->ai = ai; }//���������������
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

//XFoil������
//�������ļ�
class XFoilInfo
{
	typedef tuple<double, double, double> tp;

public:
	XFoilInfo() = default;
	
	string foilName{"undefine"};//��������
	vector<Point> foilList{};	//����ɢ����״����
	//
	string workPath{""};		//�����ļ���·��
	string shapePath() { return workPath + "\\" + foilName + ".dat"; } //�����ļ����Ͳ�����·��
	string inputPath() { return workPath + "\\" + foilName + "_input.dat"; }//�����ļ����Ͳ�����·��
	string resultPath() { return workPath + "\\" + foilName + ".cp"; }//�����ļ����Ͳ�����·��
	//
	double VISC{0};				//��ŵ������
	double Ma{0};				//�����
	int ITER{50};				//�����Ĳ���
	int outputNum{160};			//����Ľڵ���� Ĭ��160�����300�ڵ�
	// int setAlpha(double first, double secend, int gap);
	// const tp &Alpha() { return alpha; }


	tp alpha{0,0,0};//��¼���ǵ���ʼ��Χ�ͼ��

};

class XFoilSolver
{
public:
	XFoilSolver(XFoilInfo &info) : m_info(info){}
	int solve();
	string ExePath{"xfoil.exe"};
	const vector<Field<double>> getRes() { return m_res; }

private:
	int printShape();//��ӡ�����ļ�
	int printInput();//��ӡ�����ļ�
	int readRes();
	vector<Field<double>> m_res{};
	XFoilInfo &m_info;//�������
};