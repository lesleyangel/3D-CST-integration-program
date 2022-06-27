#pragma once
#include"AbstructShape.h"

class SMB1
{
public:
	vector<double> Origin;				//Origin X/Y/Z
	vector<double> Rotation;			//Rotation X/Y/Z
	vector<double> Length;				//XBody YuppBody YlowBody ZBody
	vector<double> N1, N2, M1, T1;		//形状因子
	double Ratio1;						//形状因子过渡系数
	int NFaiU, NFaiL, NEta, NHeight;	//网格数
	vector<vector<double>> BUPP, BLOW, DUPP, DLOW;//侧向表面控制因子
	int GridRefineType;

	SMB1();
	SMB1(vector<double>& Origin,						//Origin X/Y/Z
		vector<double>& Rotation,						//Rotation X/Y/Z
		vector<double>& Length,							//XBody YuppBody YlowBody ZBody
		vector<double>& N1, vector<double>& N2,
		vector<double>& M1, vector<double>& T1,			//形状因子
		double& Ratio1,									//形状因子过渡系数
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,//网格数
		vector<vector<double>>& BUPP,		vector<vector<double>>& BLOW,
		vector<vector<double>>& DUPP,		vector<vector<double>>& DLOW,//表面控制因子
		int& GirdRefineType);					
	void readFromFile(ifstream& ifs);
};

class ShapeMainBody1 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { mat A; return A; };
	//void SaveToAero(string name);//生成用于aero_calc计算气动的网格单元文件**4.23
	vector<ShellandBeam> BuildStrc1D();

	ShapeMainBody1();
	ShapeMainBody1(const SMB1& smb1);

	SMB1 m_Body1;
};

class SMB2
{
public:
	vector<double> Origin;
	vector<double> Rotation;
	vector<double> LHead;//XHead, YuppHead, YlowHead, ZHead
	vector<double> LBody;//XBody, YuppBody, YlowBody, ZBody
	vector<double> N1, N2, N3;
	vector<double> M1, T1;
	double Ratio1, Ratio2;
	int NFaiU, NFaiL, NEta, NHeight;
	vector<vector<double>> BUPP1,  BLOW1,  DUPP1,  DLOW1;
	vector<vector<double>> BUPP2,  BLOW2,  DUPP2,  DLOW2;
	int GridRefineType;

	SMB2();
	SMB2(
		vector<double>& Origin,
		vector<double>& Rotation,
		vector<double>& LHead,//XHead, YuppHead, YlowHead, ZHead
		vector<double>& LBody,//XBody, YuppBody, YlowBody, ZBody
		vector<double>& N1, vector<double>& N2, vector<double>& N3,
		vector<double>& M1, vector<double>& T1,
		double& Ratio1, double& Ratio2,
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
		vector<vector<double>>& BUPP1, vector<vector<double>>& BLOW1,
		vector<vector<double>>& DUPP1, vector<vector<double>>& DLOW1,
		vector<vector<double>>& BUPP2, vector<vector<double>>& BLOW2,
		vector<vector<double>>& DUPP2, vector<vector<double>>& DLOW2,
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeMainBody2 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();
	vector<ShellandBeam> BuildStrc1D();

	ShapeMainBody2();
	ShapeMainBody2(SMB2& smb2);

	SMB2 m_Body2;
};

class SMB3
{
public:
	vector<double> Origin;
	vector<double> Rotation;
	vector<double> Length;//XHead, XBody, XTail, YuppBody, YlowBody, ZBody
	vector<double> N1, N2, N3, N4;
	vector<double> M1,M3, T1, T3;
	double Ratio1, Ratio2, Ratio3;
	int NFaiU, NFaiL, NEta, NHeight;
	vector<vector<double>> BUPP1, BLOW1, DUPP1, DLOW1;
	vector<vector<double>> BUPP2, BLOW2, DUPP2, DLOW2;
	vector<vector<double>> BUPP3, BLOW3, DUPP3, DLOW3;
	int GridRefineType;

	SMB3();
	SMB3(
		vector<double>& Origin,
		vector<double>& Rotation,
		vector<double>& Length,//XHead, XBody, XTail, YuppBody, YlowBody, ZBody
		vector<double>& N1, vector<double>& N2, vector<double>& N3, vector<double>& N4,
		vector<double>& M1, vector<double>& M3, vector<double>& T1, vector<double>& T3,
		double& Ratio1, double& Ratio2, double& Ratio3,
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
		vector<vector<double>>& BUPP1, vector<vector<double>>& BLOW1,
		vector<vector<double>>& DUPP1, vector<vector<double>>& DLOW1,
		vector<vector<double>>& BUPP2, vector<vector<double>>& BLOW2,
		vector<vector<double>>& DUPP2, vector<vector<double>>& DLOW2,
		vector<vector<double>>& BUPP3, vector<vector<double>>& BLOW3,
		vector<vector<double>>& DUPP3, vector<vector<double>>& DLOW3, 
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeMainBody3 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { mat A; return A; };
	vector<ShellandBeam> BuildStrc1D();

	ShapeMainBody3();
	ShapeMainBody3(SMB3& smb3);

	SMB3 m_Body3;
};

class SW1
{
public:
	vector<double> Origin;				//Origin X/Y/Z
	vector<double> Ninner, Nouter;		//形状因子
	double Ratio1;						//形状因子过渡系数
	double SpanL, RootChordL, TipRootRatio;
	double SweepBackAngle, Thickness;
	int NFaiU, NFaiL, NEta, NHeight;	//网格数
	vector<vector<double>> BUPP, BLOW, DUPP, DLOW;//侧向表面控制因子
	int GridRefineType;
	int Struct2dXnum, Struct2dZnum;
	bool ifUseMid;//是否生成中间的连接网格(默认为1存在)

	SW1();
	SW1(
		vector<double>& Origin,				//Origin X/Y/Z
		vector<double>& Ninner, vector<double>& Nouter,	//形状因子
		double& Ratiol,					//形状因子过渡系数
		double& SpanL, double& RootChordL, double& TipRootRatio,
		double& SweepBackAngle, double& Thickness,
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,	//网格数
		vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
		vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,//侧向表面控制因子
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeWing1 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { mat A; return A; };
	vector<ShellandBeam> BuildStrc1D();

	ShapeWing1();
	ShapeWing1(SW1& sw1);

	SW1 m_Wing1;
};

class SW2 
{
public:
	vector<double> Origin;				//Origin X/Y/Z
	vector<double> Ninner, Nmiddile, Nouter;		//形状因子
	double Ratio_inn,Ratio_out;						//形状因子过渡系数
	double SpanL_inn, SpanL_out, RootChordL_inn, TipRootRatio_inn, TipRootRatio_out;
	double SweepBackAngle_inn, SweepBackAngle_out, Thickness_inn;
	int NFaiU, NFaiL, NEta_inn,NEta_out, NHeight;	//网格数
	vector<vector<double>> BUPP1, BLOW1, DUPP1, DLOW1;//侧向表面控制因子
	vector<vector<double>> BUPP2, BLOW2, DUPP2, DLOW2;//侧向表面控制因子
	int GridRefineType;
	//int Struct2dXnum, Struct2dZnum;

	bool ifUseMid;//是否生成中间的连接网格(默认为1存在)

	SW2();
	SW2(vector<double>& Origin,				//Origin X/Y/Z
		vector<double>& Ninner, vector<double>& Nmiddile, vector<double>& Nouter,		//形状因子
		double& Ratio_inn, double& Ratio_out,						//形状因子过渡系数
		double& SpanL_inn, double& SpanL_out, double& RootChordL_inn,
		double& TipRootRatio_inn, double& TipRootRatio_out,
		double& SweepBackAngle_inn, double& SweepBackAngle_out, double& Thickness_inn,
		int& NFaiU, int& NFaiL, int& NEta_inn, int& NEta_out, int& NHeight,	//网格数
		vector<vector<double>>& BUPP1, vector<vector<double>>& BLOW1,
		vector<vector<double>>& DUPP1, vector<vector<double>>& DLOW1,//侧向表面控制因子
		vector<vector<double>> BUPP2, vector<vector<double>>& BLOW2,
		vector<vector<double>>& DUPP2, vector<vector<double>>& DLOW2,//侧向表面控制因子
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeWing2 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { mat A; return A; };
	vector<ShellandBeam> BuildStrc1D();

	ShapeWing2();
	ShapeWing2(SW2& sw2);

	SW2 m_Wing2;
};

class ST1
{
public:
	vector<double> Origin;
	vector<double> Ninner, Nouter;
	double Ratio1;
	double SpanL, RootChordL, TipRootRatio, SweepBackAngle, Thickness;
	int NFaiU, NFaiL, NEta, NHeight;
	int GridRefineType;
	vector<vector<double>> BUPP, BLOW, DUPP, DLOW;

	ST1();
	ST1(
		vector<double>& Origin,
		vector<double>& Ninner, vector<double>& Nouter,
		double& Ratio1,
		double& SpanL, double& RootChordL, double& TipRootRatio, double& SweepBackAngle, double& Thickness,
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
		vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
		vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeTail1 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { vector<SingleBone> a; return a; };
	vector<ShellandBeam> BuildStrc1D();

	ShapeTail1();
	ShapeTail1(ST1& st1);

	ST1 m_Tail1;
};

class ST2
{
public:
	vector<double> Origin;
	vector<double> Ninner, Nouter;
	double Ratio1;
	double SpanL, RootChordL, TipRootRatio, SweepBackAngle, Thickness,SideAngle;
	int NFaiU, NFaiL, NEta, NHeight;
	int GridRefineType;
	vector<vector<double>> BUPP, BLOW, DUPP, DLOW;

	ST2();
	ST2(
		vector<double>& Origin,
		vector<double>& Ninner, vector<double>& Nouter,
		double& Ratio1,
		double& SpanL, double& RootChordL, double& TipRootRatio, 
		double& SweepBackAngle, double& Thickness,double& SideAngel,
		int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
		vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
		vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,
		int& GridRefineType);
	void readFromFile(ifstream& ifs);
};

class ShapeTail2 :public AbstructShape
{
public:
	void BuildShape();
	vector<SingleBone> BuildBone();// { vector<SingleBone> a; return a; };
	vector<ShellandBeam> BuildStrc1D();

	ShapeTail2();
	ShapeTail2(ST2& st2);

	ST2 m_Tail2;
};