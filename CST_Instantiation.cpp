#include"CST_Instantiation.h"
#include<fstream>
// #include<boost/assign.hpp>
// #include<boost/foreach.hpp>
#include<cmath>
// using namespace boost::assign;
//using namespace arma;
#define pi 3.1415926
#define ELSE_IF_STR_TO(name, res)    \
	else if (!str_tmp.compare(name)) \
		STR_TO(res)
#define ELSE_IF_STR_TO_VECTOR(name, res) \
	else if (!str_tmp.compare(name))     \
		STR_TO_VECTOR(res)
#define ELSE_IF_STR_TO_VVECTOR(name, res) \
	else if (!str_tmp.compare(name))      \
		STR_TO_VVECTOR(res)
#define STR_TO(res)                \
	{                              \
		stringstream iss;          \
		getline(ss, str_tmp, ','); \
		iss << str_tmp;            \
		iss >> res;                \
	}
#define STR_TO_VECTOR(res)                      \
	{                                           \
		stringstream iss;                       \
		for (size_t i = 0; i < res.size(); i++) \
		{                                       \
			getline(ss, str_tmp, ',');          \
			iss << str_tmp;                     \
			iss >> res[i];                      \
			iss.clear();                        \
		}                                       \
	}
#define STR_TO_VVECTOR(res)                         \
	{                                               \
		res.clear();                                \
		while (getline(ss, str_tmp, ';'))           \
		{                                           \
			vector<double> row;                     \
			stringstream ss_first(str_tmp);         \
			while (getline(ss_first, str_tmp, ',')) \
			{                                       \
				double temp;                        \
				stringstream iss;                   \
				iss << str_tmp;                     \
				iss >> temp;                        \
				row.push_back(temp);                \
			}                                       \
			res.push_back(row);                     \
		}                                           \
	}
//---------------------------------------------------------------
//------------------------ShapeMainBody1-------------------------
//初始化ShapeMainBody1参数
SMB1::SMB1()
{
	Origin = {0, 0, 0};
	Rotation = {0, 0, 0};
	Length = {3000, 400, 150, 1000};
	N1 = {0.5, 0.5, 5, 5};
	N2 = {0.5, 0.5, 5, 5};
	M1 = {0.5, 0, 0.5, 0};
	T1 = {0.5, 0, 0.5, 0};
	Ratio1 = 1;
	NFaiU = 10;
	NFaiL = 10;
	NEta = 10;
	NHeight = 20;
	vector<double> v1;
	v1 = {1};
	BUPP = {v1};
	BLOW = {v1};
	DUPP = {v1};
	DLOW = {v1};
	GridRefineType = -1;
}

SMB1::SMB1(
	vector<double>& Origin,							//Origin X/Y/Z
	vector<double>& Rotation,						//Rotation X/Y/Z
	vector<double>& Length,							//XBody YuppBody YlowBody ZBody
	vector<double>& N1, vector<double>& N2,
	vector<double>& M1, vector<double>& T1,			//形状因子
	double& Ratio1,									//形状因子过渡系数
	int& NFaiU, int& NFaiL, int& NEta, int& NHeight,//网格数
	vector<vector<double>>& BUPP,	vector<vector<double>>& BLOW,
	vector<vector<double>>& DUPP,	vector<vector<double>>& DLOW,//表面控制因子
	int& GirdRefineType)					
{
	this->Origin = Origin;
	this->Rotation = Rotation;
	this->Length = Length;
	this->N1 = N1;
	this->N2 = N2;
	this->M1 = M1;
	this->T1 = T1;
	this->Ratio1 = Ratio1;
	this->NFaiU = NFaiU;
	this->NFaiL = NFaiL;
	this->NEta = NEta;
	this->NHeight = NHeight;
	
	this->BUPP = BUPP;
	this->BLOW = BLOW;
	this->DUPP = DUPP;
	this->DLOW = DLOW;
	this->GridRefineType = GirdRefineType;
}

void SMB1::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"			,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Rotation"		,this-> Rotation)
		ELSE_IF_STR_TO_VECTOR ("Length"			,this-> Length)
		ELSE_IF_STR_TO_VECTOR ("N1"				,this-> N1)
		ELSE_IF_STR_TO_VECTOR ("N2"				,this-> N2)
		ELSE_IF_STR_TO_VECTOR ("M1"				,this-> M1)
		ELSE_IF_STR_TO_VECTOR ("T1"				,this-> T1)
		ELSE_IF_STR_TO        ("Ratio1"			,this->Ratio1)
		ELSE_IF_STR_TO        ("NFaiU"			,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"			,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"			,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"		,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP"			,this->BUPP)
		ELSE_IF_STR_TO_VVECTOR("BLOW"			,this->BLOW)
		ELSE_IF_STR_TO_VVECTOR("DUPP"			,this->DUPP)
		ELSE_IF_STR_TO_VVECTOR("DLOW"			,this->DLOW)
		ELSE_IF_STR_TO        ("GridRefineType"	,this->GridRefineType)
	}
}

ShapeMainBody1::ShapeMainBody1() { ShapeType = 0; }

ShapeMainBody1::ShapeMainBody1(const SMB1& smb1)
{
	m_Body1 = smb1;
	ShapeType = 0;//曲面类型（0 头身部；1 有连接的机翼；2 无中间连接的机翼和垂尾）
}

//ShapeMainBody1中重写BuildShape
void ShapeMainBody1::BuildShape()
{
	cout << "Shape Main Body1: 参数正在传入。。。" << endl;
	CSTsurface Head;
	
	Head.Origin = m_Body1.Origin;
	Head.Rotation = m_Body1.Rotation;
	Head.Slope = { 0,0,0 };
	Head.Scale = { 1,1,1,1 };
	Head.Length = m_Body1.Length;
	Head.NS = m_Body1.N1;
	Head.NE = m_Body1.N2;
	Head.M = m_Body1.M1;
	Head.T = m_Body1.T1;
	Head.Ratio = m_Body1.Ratio1;
	copyM(m_Body1.BUPP, Head.BUPP);
	copyM(m_Body1.BLOW, Head.BLOW);
	copyM(m_Body1.DUPP, Head.DUPP);
	copyM(m_Body1.DLOW, Head.DLOW);
	Head.NFaiU = m_Body1.NFaiU;
	Head.NFaiL = m_Body1.NFaiL;
	Head.NEta = m_Body1.NEta;
	Head.NHeight = m_Body1.NHeight;
	Head.Is_MeshFront = 0;
	Head.Is_MeshBack = 0;

	CSTsf.push_back(Head);
	CSTsf.reserve(1);
	GridRefineType = m_Body1.GridRefineType;

	cout << "Shape Main Body1: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Main Body1: 计算结束！" << endl;
	/*Head.MakeBone();*/
}

//ShapeMainBody1中重写BuildBone
vector<SingleBone> ShapeMainBody1::BuildBone()
{
	mat bone = CSTsf[0].MakeBone();
	mat boneF = CSTsf[0].ClacBonePForce();//节点力大小

	vector<SingleBone> SB;
	SingleBone sb;
	sb.BoneP = bone;
	sb.BonePForce = boneF;
	SB.push_back(sb);
	return SB;
}

vector<ShellandBeam> ShapeMainBody1::BuildStrc1D()
{
	return vector<ShellandBeam>();
}



//---------------------------------------------------------------
//------------------------ShapeMainBody2-------------------------
//初始化ShapeMainBody2参数
SMB2::SMB2()
{
	Origin = {0,0,0};
	Rotation = {0,0,0};
	LHead = {(5000.0),(3350.0 * 0.52),(3350.0 * 0.48),(3350.0)};//XHead, YuppHead, YlowHead, ZHead
	LBody = {(18000),(3350 * 0.52),(3350 * 0.48),(3350)};//XBody, YuppBody, YlowBody, ZBody
	N1 = {0.4, 0.4, 0.4, 0.4};
	N2 = {0.5, 0.5, 0.2, 0.2};
	N3 = N2;
	M1 ={ 0.5, 0, 0.5, 0};
	T1 ={ 0.5, 0, 0.5, 0};
	Ratio1 = 1, Ratio2 = 1;
	NFaiU = 30, NFaiL = 40, NEta = 50, NHeight = 20;
	BUPP1 = {{1}}, BLOW1 = {{1}},
	DUPP1 = {{1}}, DLOW1 = {{1}},
	BUPP2 = {{1}}, BLOW2 = {{1}},
	DUPP2 = {{1}}, DLOW2 = {{1}};
	GridRefineType = 0;
}

SMB2::SMB2(
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
	int& GirdRefineType)
{
	this->Origin = Origin;
	this->Rotation = Rotation;
	this->LHead = LHead;//XHead, YuppHead, YlowHead, ZHead
	this->LBody = LBody;//XBody, YuppBody, YlowBody, ZBody
	this->N1 = N1; this->N2 = N2; this->N3 = N3;
	this->M1 = M1; this->T1 = T1;
	this->Ratio1 = Ratio1; this->Ratio2 = Ratio2;
	this->NFaiU = NFaiU; this->NFaiL = NFaiL; this->NEta = NEta; this->NHeight = NHeight;
	this->BUPP1 = BUPP1; this->BLOW1 = BLOW1; this->DUPP1 = DUPP1; this->DLOW1 = DLOW1;
	this->BUPP2 = BUPP2; this->BLOW2 = BLOW2; this->DUPP2 = DUPP2; this->DLOW2 = DLOW2;
	this->GridRefineType = GirdRefineType;
}

void SMB2::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"			,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Rotation"		,this-> Rotation)
		ELSE_IF_STR_TO_VECTOR ("LHead"			,this-> LHead)
		ELSE_IF_STR_TO_VECTOR ("LBody"			,this-> LBody)
		ELSE_IF_STR_TO_VECTOR ("N1"				,this-> N1)
		ELSE_IF_STR_TO_VECTOR ("N2"				,this-> N2)
		ELSE_IF_STR_TO_VECTOR ("N3"				,this-> N3)
		ELSE_IF_STR_TO_VECTOR ("M1"				,this-> M1)
		ELSE_IF_STR_TO_VECTOR ("T1"				,this-> T1)
		ELSE_IF_STR_TO        ("Ratio1"			,this->Ratio1)
		ELSE_IF_STR_TO        ("Ratio2"			,this->Ratio2)
		ELSE_IF_STR_TO        ("NFaiU"			,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"			,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"			,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"		,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP1"			,this->BUPP1)
		ELSE_IF_STR_TO_VVECTOR("BLOW1"			,this->BLOW1)
		ELSE_IF_STR_TO_VVECTOR("DUPP1"			,this->DUPP1)
		ELSE_IF_STR_TO_VVECTOR("DLOW1"			,this->DLOW1)
		ELSE_IF_STR_TO_VVECTOR("BUPP2"			,this->BUPP2)
		ELSE_IF_STR_TO_VVECTOR("BLOW2"			,this->BLOW2)
		ELSE_IF_STR_TO_VVECTOR("DUPP2"			,this->DUPP2)
		ELSE_IF_STR_TO_VVECTOR("DLOW2"			,this->DLOW2)
		ELSE_IF_STR_TO        ("GridRefineType"	,this->GridRefineType)
	}
}

ShapeMainBody2::ShapeMainBody2() { ShapeType = 0; }

ShapeMainBody2::ShapeMainBody2(SMB2& smb2)
{
	m_Body2 = smb2;
	ShapeType = 0;
}

void ShapeMainBody2::BuildShape()
{
	cout << "Shape Main Body2: 参数正在传入。。。" << endl;
	CSTsurface Head;
	CSTsurface Body;

	Head.Origin = m_Body2.Origin;
	Head.Rotation = m_Body2.Rotation;
	Head.Slope = { 0,0,0 };
	Head.Scale = { 1,1,1,1 };
	Head.Length = m_Body2.LHead;
	Head.NS = m_Body2.N1;
	Head.NE = m_Body2.N2;
	Head.Ratio = m_Body2.Ratio1;
	Head.M = m_Body2.M1;
	Head.T = m_Body2.T1;
	Head.NFaiU = m_Body2.NFaiU;
	Head.NFaiL = m_Body2.NFaiL;
	Head.NEta = m_Body2.NEta;
	Head.NHeight = m_Body2.NHeight;
	copyM(m_Body2.BUPP1, Head.BUPP);
	copyM(m_Body2.BLOW1, Head.BLOW);
	copyM(m_Body2.DUPP1, Head.DUPP);
	copyM(m_Body2.DLOW1, Head.DLOW);
	Head.Is_MeshFront = 0;
	Head.Is_MeshBack = 0;

	Body.Origin = m_Body2.Origin;
	Body.Origin(0) += m_Body2.LHead[0];
	Body.Rotation = m_Body2.Rotation;
	Body.Slope = { 0,0,0 };
	Body.Scale = { 1, (m_Body2.LBody[1] / m_Body2.LHead[1]), (m_Body2.LBody[2] / m_Body2.LHead[2]), (m_Body2.LBody[3] / m_Body2.LHead[3]) };
	Body.Length = m_Body2.LHead;
	Body.Length(0) = m_Body2.LBody[0];
	Body.NS = m_Body2.N2;
	Body.NE = m_Body2.N3;
	Body.Ratio = m_Body2.Ratio2;
	Body.M = { 0,0,0,0 };
	Body.T = { 0,0,0,0 };
	Body.NFaiU = m_Body2.NFaiU;
	Body.NFaiL = m_Body2.NFaiL;
	Body.NEta = m_Body2.NEta;
	Body.NHeight = m_Body2.NHeight;
	copyM(m_Body2.BUPP2, Body.BUPP);
	copyM(m_Body2.BLOW2, Body.BLOW);
	copyM(m_Body2.DUPP2, Body.DUPP);
	copyM(m_Body2.DLOW2, Body.DLOW);
	Body.Is_MeshFront = 0;
	Body.Is_MeshBack = 0;//后面再打开

	CSTsf.push_back(Head);
	CSTsf.push_back(Body);
	GridRefineType = m_Body2.GridRefineType;

	cout << "Shape Main Body2: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Main Body2: 计算结束！" << endl;

	//Head.MakeBone();
}

vector<SingleBone> ShapeMainBody2::BuildBone()
{
	mat bone = CSTsf[0].MakeBone();//第一段的梁节点
	mat boneF = CSTsf[0].ClacBonePForce();//节点力大小
	for (size_t i = 1; i < CSTsf.size(); i++)//后面的梁节点
	{
		mat pList = CSTsf[i].MakeBone();
		bone = join_cols(bone, pList.rows(1, pList.n_rows-1));//依次添加到总部件梁节点内

		mat fList = CSTsf[i].ClacBonePForce();
		if (fList.n_rows != 0)
		{
			boneF.row(boneF.n_rows - 1) += fList.row(0);
			boneF = join_cols(boneF, fList.rows(1, fList.n_rows - 1));
		}
	}
	vector<SingleBone> SB;
	SingleBone sb;
	sb.BoneP = bone;
	sb.BonePForce = boneF;
	SB.push_back(sb);

	return SB;
}

vector<ShellandBeam> ShapeMainBody2::BuildStrc1D()
{
	vector<ShellandBeam> A(0);
	return A;
}
//---------------------------------------------------------------
//------------------------ShapeMainBody3-------------------------
//初始化ShapeMainBody3参数
SMB3::SMB3()
{
	this->Origin = { 0,0,0 };	this->Rotation = { 0,0,0 };	
	this->Length = { 1700,3800,2500,850,650,1600 };
	for (size_t i = 0; i < this->Length.size(); i++)
	{
		this->Length[i] *= 2.4;
	}
	this->N1 = {0.5,0.5,0.2,0.2};	this->N2 = N1;	
	this->N3 = N2;	this->N4 = { 0.1,0.1,0.1,0.1 };
	this->M1 = { 0.5,0.0,0.25,0.0 };	this->M3 = { 0.0,0.0,0.0,0.0 };	
	this->T1 = { 0.3,0.0,0.3 ,0.0 };	this->T3 = { 0.0,0.0,0.0,0.0 };
	this->Ratio1 = 1;	this->Ratio2 = 1;	this->Ratio3 = 1;
	this->NFaiU = 40;	this->NFaiL = 40;
	this->NEta = 40;		this->NHeight = 20;
	this->BUPP1 = { {1} };	this->BLOW1 = { {1} };
	this->DUPP1 = { {1} };	this->DLOW1 = { {1} };
	this->BUPP2 = { {1} };	this->BLOW2 = { {1} };
	this->DUPP2 = { {1} };	this->DLOW2 = { {1} };
	this->BUPP3 = { {1} };	this->BLOW3 = { {1} };
	this->DUPP3 = { {1,1.1,1.2,1.2} };	this->DLOW3 = { {1,1.1,1.2,1.2} };
	this->GridRefineType = 0;
}

SMB3::SMB3(
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
	int& GridRefineType)
{
	this->Origin = Origin;	this->Rotation = Rotation;	this->Length = Length;
	this->N1 = N1;	this->N2 = N2;	this->N3 = N3;	this->N4 = N4;
	this->M1 = M1;	this->M3 = M3;	this->T1 = T1;	this->T3 = T3;
	this->Ratio1 = Ratio1;	this->Ratio2 = Ratio2;	this->Ratio3 = Ratio3;
	this->NFaiU = NFaiU;	this->NFaiL = NFaiL;
	this->NEta = NEta;		this->NHeight = NHeight;
	this->BUPP1 = BUPP1;	this->BLOW1 = BLOW1;
	this->DUPP1 = DUPP1;	this->DLOW1 = DLOW1;
	this->BUPP2 = BUPP2;	this->BLOW2 = BLOW2;
	this->DUPP2 = DUPP2;	this->DLOW2 = DLOW2;
	this->BUPP3 = BUPP3;	this->BLOW3 = BLOW3;
	this->DUPP3 = DUPP3;	this->DLOW3 = DLOW3;
}

void SMB3::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"			,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Rotation"		,this-> Rotation)
		ELSE_IF_STR_TO_VECTOR ("Length"			,this-> Length)
		ELSE_IF_STR_TO_VECTOR ("N1"				,this-> N1)
		ELSE_IF_STR_TO_VECTOR ("N2"				,this-> N2)
		ELSE_IF_STR_TO_VECTOR ("N3"				,this-> N3)
		ELSE_IF_STR_TO_VECTOR ("N4"				,this-> N4)
		ELSE_IF_STR_TO_VECTOR ("M1"				,this-> M1)
		ELSE_IF_STR_TO_VECTOR ("M3"				,this-> M3)
		ELSE_IF_STR_TO_VECTOR ("T1"				,this-> T1)
		ELSE_IF_STR_TO_VECTOR ("T3"				,this-> T3)
		ELSE_IF_STR_TO        ("Ratio1"			,this->Ratio1)
		ELSE_IF_STR_TO        ("Ratio2"			,this->Ratio2)
		ELSE_IF_STR_TO        ("Ratio3"			,this->Ratio3)
		ELSE_IF_STR_TO        ("NFaiU"			,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"			,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"			,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"		,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP1"			,this->BUPP1)
		ELSE_IF_STR_TO_VVECTOR("BLOW1"			,this->BLOW1)
		ELSE_IF_STR_TO_VVECTOR("DUPP1"			,this->DUPP1)
		ELSE_IF_STR_TO_VVECTOR("DLOW1"			,this->DLOW1)
		ELSE_IF_STR_TO_VVECTOR("BUPP2"			,this->BUPP2)
		ELSE_IF_STR_TO_VVECTOR("BLOW2"			,this->BLOW2)
		ELSE_IF_STR_TO_VVECTOR("DUPP2"			,this->DUPP2)
		ELSE_IF_STR_TO_VVECTOR("DLOW2"			,this->DLOW2)
		ELSE_IF_STR_TO_VVECTOR("BUPP3"			,this->BUPP3)
		ELSE_IF_STR_TO_VVECTOR("BLOW3"			,this->BLOW3)
		ELSE_IF_STR_TO_VVECTOR("DUPP3"			,this->DUPP3)
		ELSE_IF_STR_TO_VVECTOR("DLOW3"			,this->DLOW3)
		ELSE_IF_STR_TO        ("GridRefineType"	,this->GridRefineType)
	}
}

ShapeMainBody3::ShapeMainBody3() { ShapeType = 0; }

vector<ShellandBeam> ShapeMainBody3::BuildStrc1D()
{
	return vector<ShellandBeam>();
}

ShapeMainBody3::ShapeMainBody3(SMB3& smb3)
{
	m_Body3 = smb3;
	ShapeType = 0;
}

void ShapeMainBody3::BuildShape()
{
	cout << "Shape Main Body3: 参数正在传入。。。" << endl;
	shared_ptr<CSTsurface> Head(new CSTsurface);
	shared_ptr<CSTsurface> Body(new CSTsurface);
	shared_ptr<CSTsurface> Tail(new CSTsurface);

	(*Head).Origin = m_Body3.Origin;
	(*Head).Rotation = m_Body3.Rotation;
	(*Head).Slope = { 0,0,0 };
	(*Head).Scale = { 1,1,1,1 };
	(*Head).Length = { m_Body3.Length[0],m_Body3.Length[3] ,m_Body3.Length[4] ,m_Body3.Length[5] };
	(*Head).NS = m_Body3.N1;
	(*Head).NE = m_Body3.N2;
	(*Head).Ratio = m_Body3.Ratio1;
	(*Head).M = m_Body3.M1;
	(*Head).T = m_Body3.T1;
	(*Head).NFaiU = m_Body3.NFaiU;
	(*Head).NFaiL = m_Body3.NFaiL;
	(*Head).NEta = m_Body3.NEta;
	(*Head).NHeight = m_Body3.NHeight;
	copyM(m_Body3.BUPP1, (*Head).BUPP);
	copyM(m_Body3.BLOW1, (*Head).BLOW);
	copyM(m_Body3.DUPP1, (*Head).DUPP);
	copyM(m_Body3.DLOW1, (*Head).DLOW);
	(*Head).Is_MeshFront = 0;
	(*Head).Is_MeshBack = 0;

	(*Body).Origin = m_Body3.Origin;
	//(*Body).Origin(0) += m_Body3.Length[0];
	(*Body).Rotation = m_Body3.Rotation;
	(*Body).Slope = { 0,0,0 };
	(*Body).Scale = { 1,1,1,1 };
	(*Body).Length = { m_Body3.Length[1],m_Body3.Length[3],m_Body3.Length[4],m_Body3.Length[5] };
	(*Body).NS = m_Body3.N2;
	(*Body).NE = m_Body3.N3;
	(*Body).Ratio = m_Body3.Ratio2;
	(*Body).M = { 0,0,0,0 };
	(*Body).T = { 0,0,0,0 };
	(*Body).NFaiU = m_Body3.NFaiU;
	(*Body).NFaiL = m_Body3.NFaiL;
	(*Body).NEta = m_Body3.NEta;
	(*Body).NHeight = m_Body3.NHeight;
	copyM(m_Body3.BUPP2, (*Body).BUPP);
	copyM(m_Body3.BLOW2, (*Body).BLOW);
	copyM(m_Body3.DUPP2, (*Body).DUPP);
	copyM(m_Body3.DLOW2, (*Body).DLOW);
	(*Body).Is_MeshFront = 0;
	(*Body).Is_MeshBack = 0;

	(*Tail).Origin = m_Body3.Origin;
	//(*Tail).Origin(0) += m_Body3.Length[0]+ m_Body3.Length[1];
	(*Tail).Rotation = m_Body3.Rotation;
	(*Tail).Slope = { 0,0,0 };
	(*Tail).Scale = { 1,1,1,1 };
	(*Tail).Length = { m_Body3.Length[2],m_Body3.Length[3],m_Body3.Length[4],m_Body3.Length[5] };
	(*Tail).NS = m_Body3.N3;
	(*Tail).NE = m_Body3.N4;
	(*Tail).Ratio = m_Body3.Ratio3;
	(*Tail).M = m_Body3.M3;
	(*Tail).T = m_Body3.T3;
	(*Tail).NFaiU = m_Body3.NFaiU;
	(*Tail).NFaiL = m_Body3.NFaiL;
	(*Tail).NEta = m_Body3.NEta;
	(*Tail).NHeight = m_Body3.NHeight;
	copyM(m_Body3.BUPP3, (*Tail).BUPP);
	copyM(m_Body3.BLOW3, (*Tail).BLOW);
	copyM(m_Body3.DUPP3, (*Tail).DUPP);
	copyM(m_Body3.DLOW3, (*Tail).DLOW);
	(*Tail).Is_MeshFront = 0;
	(*Tail).Is_MeshBack = 0;//后面再打开

	CSTsf.push_back(*Head);
	CSTsf.push_back(*Body);
	CSTsf.push_back(*Tail);
	GridRefineType = m_Body3.GridRefineType;

	cout << "Shape Main Body3: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Main Body3: 计算结束！" << endl;
}

vector<SingleBone> ShapeMainBody3::BuildBone()
{
	mat bone = CSTsf[0].MakeBone();
	mat boneF = CSTsf[0].ClacBonePForce();//节点力大小
	for (size_t i = 1; i < CSTsf.size(); i++)
	{
		mat pList = CSTsf[i].MakeBone();
		bone = join_cols(bone, pList.rows(1, pList.n_rows - 1));

		mat fList = CSTsf[i].ClacBonePForce();
		if (fList.n_rows != 0)
		{
			boneF.row(boneF.n_rows - 1) += fList.row(0);
			boneF = join_cols(boneF, fList.rows(1, fList.n_rows - 1));
		}
	}

	vector<SingleBone> SB;
	SingleBone sb;
	sb.BoneP = bone;
	sb.BonePForce = boneF;
	SB.push_back(sb);
	return SB;
}

//---------------------------------------------------------------
//--------------------------ShapeWing1---------------------------
//初始化ShapeWing1参数
SW1::SW1()
{
	//M6 for test
	this->Origin = { 0,0,500 };
	this->Ninner = { 0.5,1.2,0.5,1.2 };
	this->Nouter = { 0.5,1.2,0.5,1.2 };
	this->Ratio1 = 1;
	this->SpanL = 1196.3;
	this->RootChordL = 1000;
	this->TipRootRatio = 0.5;
	this->SweepBackAngle = 45;
	this->Thickness = 100;
	this->NFaiU = 20;
	this->NFaiL = 20;
	this->NEta = 20;
	this->NHeight = 3;
	this->BUPP = { {1} };
	this->BLOW = { {1} };
	this->DUPP = { {1} };
	this->DLOW = { {1} };
	this->GridRefineType = -1;
	this->ifUseMid = true;
	Struct2dXnum = 0;
	Struct2dZnum = 0;
}

SW1::SW1(
	vector<double>& Origin,				//Origin X/Y/Z
	vector<double>& Ninner, vector<double>& Nouter,	//形状因子
	double& Ratiol,					//形状因子过渡系数
	double& SpanL, double& RootChordL, double& TipRootRatio,
	double& SweepBackAngle, double& Thickness,
	int& NFaiU, int& NFaiL, int& NEta, int& NHeight,	//网格数
	vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
	vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,//侧向表面控制因子
	int& GridRefineType)
{
	this->Origin = Origin;
	this->Ninner = Ninner;
	this->Nouter = Nouter;
	this->Ratio1 = Ratiol;
	this->SpanL = SpanL;
	this->RootChordL = RootChordL;
	this->TipRootRatio = TipRootRatio;
	this->SweepBackAngle = SweepBackAngle;
	this->Thickness = Thickness;
	this->NFaiU = NFaiU;
	this->NFaiL = NFaiL;
	this->NEta = NEta;
	this->NHeight = NHeight;
	this->BUPP = BUPP;
	this->BLOW = BLOW;
	this->DUPP = DUPP;
	this->DLOW = DLOW;
	this->GridRefineType = GridRefineType;
}

void SW1::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"			,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Ninner"			,this-> Ninner)
		ELSE_IF_STR_TO_VECTOR ("Nouter"			,this-> Nouter)
		ELSE_IF_STR_TO        ("Ratio1"			,this->Ratio1)
		ELSE_IF_STR_TO        ("SpanL"			,this->SpanL)
		ELSE_IF_STR_TO        ("RootChordL"		,this->RootChordL)
		ELSE_IF_STR_TO        ("TipRootRatio"	,this->TipRootRatio)
		ELSE_IF_STR_TO        ("SweepBackAngle"	,this->SweepBackAngle)
		ELSE_IF_STR_TO        ("Thickness"		,this->Thickness)
		ELSE_IF_STR_TO        ("NFaiU"			,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"			,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"			,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"		,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP"			,this->BUPP)
		ELSE_IF_STR_TO_VVECTOR("BLOW"			,this->BLOW)
		ELSE_IF_STR_TO_VVECTOR("DUPP"			,this->DUPP)
		ELSE_IF_STR_TO_VVECTOR("DLOW"			,this->DLOW)
		ELSE_IF_STR_TO        ("GridRefineType"	,this->GridRefineType)
		ELSE_IF_STR_TO        ("Struct2dXnum"	,this->Struct2dXnum)
		ELSE_IF_STR_TO        ("Struct2dZnum"	,this->Struct2dZnum)
		ELSE_IF_STR_TO        ("ifUseMid"		,this->ifUseMid)
	}
}

ShapeWing1::ShapeWing1() { ShapeType = 2; }

ShapeWing1::ShapeWing1(SW1& sw1)
{
	m_Wing1 = sw1;
	ShapeType = 2;
}

void ShapeWing1::BuildShape()
{
	cout << "Shape Wing1: 参数正在传入。。。" << endl;
	shared_ptr<CSTsurface> Left(new CSTsurface);
	shared_ptr<CSTsurface> Righ(new CSTsurface);
	shared_ptr<CSTsurface> Midd(new CSTsurface);

	(*Left).Origin = m_Wing1.Origin;
	(*Left).Rotation = { 0,90,180 };
	//(*Left).Rotation = { 0,0,0 };
	//(*Left).Slope = { 0,0,m_Wing1.SweepBackAngle };
	(*Left).Slope = { 0,0,atan(tan(m_Wing1.SweepBackAngle * pi / 180) - (1 - m_Wing1.TipRootRatio) * m_Wing1.RootChordL / (2 * m_Wing1.SpanL)) * 180 / pi };
	(*Left).Scale = { 1,m_Wing1.TipRootRatio,m_Wing1.TipRootRatio,m_Wing1.TipRootRatio };
	(*Left).Length = { m_Wing1.SpanL, m_Wing1.Thickness/2.0, m_Wing1.Thickness/2.0, m_Wing1.RootChordL };
	(*Left).NS = m_Wing1.Ninner;
	(*Left).NE = m_Wing1.Nouter;
	(*Left).NFaiU = m_Wing1.NFaiU;
	(*Left).NFaiL = m_Wing1.NFaiL;
	(*Left).NEta = m_Wing1.NEta;
	(*Left).NHeight = m_Wing1.NHeight;
	(*Left).Ratio = m_Wing1.Ratio1;
	(*Left).M = { 0,0,0,0 };
	(*Left).T = { 0,0,0,0 };
	(*Left).Is_MeshFront = 0;
	(*Left).Is_MeshBack = 0;// = 1;
	(*Left).Struct2dXnum = m_Wing1.Struct2dXnum;
	(*Left).Struct2dZnum = m_Wing1.Struct2dZnum;
	copyM(m_Wing1.BUPP, (*Left).BUPP);
	copyM(m_Wing1.BLOW, (*Left).BLOW);
	copyM(m_Wing1.DUPP, (*Left).DUPP);
	copyM(m_Wing1.DLOW, (*Left).DLOW);

	(*Righ).Origin = m_Wing1.Origin;
	(*Righ).Origin[2] = -m_Wing1.Origin[2];
	(*Righ).Rotation = { 0,90,0 };
	//(*Righ).Slope = { 0,0,m_Wing1.SweepBackAngle };
	//(*Righ).Slope = { 0,0,10 };
	(*Righ).Slope = { 0,0,atan(tan(m_Wing1.SweepBackAngle * pi / 180) - (1 - m_Wing1.TipRootRatio) * m_Wing1.RootChordL / (2 * m_Wing1.SpanL)) * 180 / pi };
	(*Righ).Scale = { 1,m_Wing1.TipRootRatio,m_Wing1.TipRootRatio,m_Wing1.TipRootRatio };
	(*Righ).Length = { m_Wing1.SpanL, m_Wing1.Thickness / 2.0, m_Wing1.Thickness / 2.0, m_Wing1.RootChordL };
	(*Righ).NS = m_Wing1.Ninner;
	(*Righ).NE = m_Wing1.Nouter;
	(*Righ).NFaiU = m_Wing1.NFaiU;
	(*Righ).NFaiL = m_Wing1.NFaiL;
	(*Righ).NEta = m_Wing1.NEta;
	(*Righ).NHeight = m_Wing1.NHeight;
	(*Righ).Ratio = m_Wing1.Ratio1;
	(*Righ).M = { 0,0,0,0 };
	(*Righ).T = { 0,0,0,0 };
	(*Righ).Is_MeshFront = 0;
	(*Righ).Is_MeshBack = 0;// = 1;
	(*Righ).Struct2dXnum = m_Wing1.Struct2dXnum;
	(*Righ).Struct2dZnum = m_Wing1.Struct2dZnum;
	copyM(m_Wing1.BUPP, (*Righ).BUPP);
	copyM(m_Wing1.BLOW, (*Righ).BLOW);
	copyM(m_Wing1.DUPP, (*Righ).DUPP);
	copyM(m_Wing1.DLOW, (*Righ).DLOW);

	if (m_Wing1.ifUseMid)//确认生成中间网格
	{
		(*Midd).Origin = m_Wing1.Origin;
		//(*Midd).Origin[0] = m_Wing1.Origin[0] + m_Wing1.RootChordL ;
		//(*Midd).Origin[2] = -m_Wing1.Origin[2];
		(*Midd).Rotation = { 0,90,0 };
		(*Midd).Slope = { 0,0,0 };
		(*Midd).Scale = { 1,1,1,1 };
		(*Midd).Length = { m_Wing1.Origin[2] * 2, m_Wing1.Thickness / 2, m_Wing1.Thickness / 2, m_Wing1.RootChordL };
		(*Midd).NS = m_Wing1.Ninner;
		(*Midd).NE = m_Wing1.Ninner;
		(*Midd).NFaiU = m_Wing1.NFaiU;
		(*Midd).NFaiL = m_Wing1.NFaiL;
		(*Midd).NEta = m_Wing1.NEta;
		(*Midd).NHeight = m_Wing1.NHeight;
		(*Midd).Ratio = m_Wing1.Ratio1;;
		(*Midd).M = { 0,0,0,0 };
		(*Midd).T = { 0,0,0,0 };
		(*Midd).Is_MeshFront = 0;
		(*Midd).Is_MeshBack = 0;
		vector<vector<double>> Temp = {{1}};
		copyM(Temp, (*Midd).BUPP);
		copyM(Temp, (*Midd).BLOW);
		copyM(Temp, (*Midd).DUPP);
		copyM(Temp, (*Midd).DLOW);
	}
	

	CSTsf.push_back(*Left);
	CSTsf.push_back(*Righ);
	if (m_Wing1.ifUseMid)
	{
		CSTsf.push_back(*Midd);
	}
	

	GridRefineType = m_Wing1.GridRefineType;

	cout << "Shape Wing1: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Wing1: 计算结束！" << endl;
}

vector<SingleBone> ShapeWing1::BuildBone()
{
	if (m_Wing1.ifUseMid)//（已经加力21.5.8）
	{
		mat bone = flipud(CSTsf[0].MakeBone());
		mat pList2 = CSTsf[2].MakeBone();
		bone = join_cols(bone, pList2.rows(1, pList2.n_rows - 1));
		mat pList1 = CSTsf[1].MakeBone();
		bone = join_cols(bone, pList1.rows(1, pList1.n_rows - 1));

		mat boneF = flipud(CSTsf[0].ClacBonePForce());//节点力大小
		mat fList2 = CSTsf[2].ClacBonePForce();//force
		if (fList2.n_rows != 0)
		{
			boneF.row(boneF.n_rows - 1) += fList2.row(0);//force
			boneF = join_cols(boneF, fList2.rows(1, fList2.n_rows - 1));//force
		}
		mat fList1 = CSTsf[1].ClacBonePForce();//force
		if (fList1.n_rows != 0)
		{
			boneF.row(boneF.n_rows - 1) += fList1.row(0);//force
			boneF = join_cols(boneF, fList1.rows(1, fList1.n_rows - 1));//force
		}
		vector<SingleBone> SB;
		SingleBone sb;
		sb.BoneP = bone;
		sb.BonePForce = boneF;//force
		sb.BoneConnP_ID = (bone.n_rows) / 2;
		SB.push_back(sb);
		return SB;
	}
	else//（已经加力21.5.8）
	{
		vector<SingleBone> SB;
		SingleBone sb;
		sb.BoneP = CSTsf[0].MakeBone();
		sb.BonePForce = CSTsf[0].ClacBonePForce();
		sb.BoneConnP_ID = 0;
		SB.push_back(sb);

		SingleBone sb2;
		sb2.BoneP = CSTsf[1].MakeBone();
		sb2.BonePForce = CSTsf[1].ClacBonePForce();
		sb2.BoneConnP_ID = 0;
		SB.push_back(sb2);

		return SB;
	}

	//return CSTsf[2].MakeBone();
}

vector<ShellandBeam> ShapeWing1::BuildStrc1D()
{
	vector<ShellandBeam> A(2);

	mat X = linspace(0, 1, 5);
	mat Z = linspace(0, 1, 5);
	CSTsf[0].MakeStruct1D(X, Z);
	CSTsf[1].MakeStruct1D(X, Z);

	A[0] = CSTsf[0].aero_strc;
	A[0].name = "SW1_left";
	A[1] = CSTsf[1].aero_strc;
	A[1].name = "SW1_righ";

	return A;
}

//---------------------------------------------------------------
//--------------------------ShapeWing2---------------------------
//初始化ShapeWing2参数
SW2::SW2()
{
	this->Origin = { 8000,-900,1730 };
	this->Ninner = { 0.8,1.2,0.8,1.2 };
	this->Nmiddile = { 0.9,1,0.9,1 };
	this->Nouter = { 0.9,1,0.9,1 };
	this->Ratio_inn = 1;
	this->Ratio_out = 1;
	this->SpanL_inn = 520;
	//this->SpanL_inn = 1200;//for test
	this->SpanL_out = 4500;
	this->RootChordL_inn = 12000;
	this->TipRootRatio_inn = (double)5475/(double)12000;
	this->TipRootRatio_out = (double)1469/(double)5475;
	this->SweepBackAngle_inn = atan((12000.0-5475.0)/520.0) *180 / pi;
	this->SweepBackAngle_out = atan2(5475.0 - 1469, 4500.0) *180 / pi;
	this->Thickness_inn = 600;
	this->NFaiU = 30;
	this->NFaiL = 30;
	this->NEta_inn = 10;
	this->NEta_out = 10;
	this->NHeight = 3;
	this->BUPP1 = { {1} };	this->BLOW1 = { {1} };
	this->DUPP1 = { {1} };	this->DLOW1 = { {1} };
	this->BUPP2 = { {1} };	this->BLOW2 = { {1} };
	this->DUPP2 = { {1} };	this->DLOW2 = { {1} };
	this->GridRefineType = 0;
	this->ifUseMid = true;
}

SW2::SW2(
	vector<double>& Origin,				//Origin X/Y/Z
	vector<double>& Ninner, vector<double>& Nmiddile, vector<double>& Nouter,//形状因子
	double& Ratio_inn, double& Ratio_out,						//形状因子过渡系数
	double& SpanL_inn, double& SpanL_out, double& RootChordL_inn,
	double& TipRootRatio_inn, double& TipRootRatio_out,
	double& SweepBackAngle_inn, double& SweepBackAngle_out, double& Thickness_inn,
	int& NFaiU, int& NFaiL, int& NEta_inn, int& NEta_out, int& NHeight,	//网格数
	vector<vector<double>>& BUPP1, vector<vector<double>>& BLOW1,
	vector<vector<double>>& DUPP1, vector<vector<double>>& DLOW1,//侧向表面控制因子
	vector<vector<double>> BUPP2, vector<vector<double>>& BLOW2,
	vector<vector<double>>& DUPP2, vector<vector<double>>& DLOW2,//侧向表面控制因子
	int& GridRefineType)
{
	this->Origin = Origin;
	this->Ninner = Ninner;
	this->Nmiddile = Nmiddile;
	this->Nouter = Nouter;
	this->Ratio_inn = Ratio_inn;
	this->Ratio_out = Ratio_out;
	this->SpanL_inn = SpanL_inn;
	this->SpanL_out = SpanL_out;
	this->RootChordL_inn = RootChordL_inn;
	this->TipRootRatio_inn = TipRootRatio_inn;
	this->TipRootRatio_out = TipRootRatio_out;
	this->SweepBackAngle_inn = SweepBackAngle_inn;
	this->SweepBackAngle_out = SweepBackAngle_out;
	this->Thickness_inn = Thickness_inn;
	this->NFaiU = NFaiU;
	this->NFaiL = NFaiL;
	this->NEta_inn = NEta_inn;
	this->NEta_out = NEta_out;
	this->NHeight = NHeight;
	this->BUPP1 = BUPP1;	this->BLOW1 = BLOW1;
	this->DUPP1 = DUPP1;	this->DLOW1 = DLOW1;
	this->BUPP2 = BUPP2;	this->BLOW2 = BLOW2;
	this->DUPP2 = DUPP2;	this->DLOW2 = DLOW2;
	this->GridRefineType = GridRefineType;
}

void SW2::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"				,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Ninner"				,this-> Ninner)
		ELSE_IF_STR_TO_VECTOR ("Nmiddile"			,this-> Nmiddile)
		ELSE_IF_STR_TO_VECTOR ("Nouter"				,this-> Nouter)
		ELSE_IF_STR_TO        ("Ratio_inn"			,this->Ratio_inn)
		ELSE_IF_STR_TO        ("Ratio_out"			,this->Ratio_out)
		ELSE_IF_STR_TO        ("SpanL_inn"			,this->SpanL_inn)
		ELSE_IF_STR_TO        ("SpanL_out"			,this->SpanL_out)
		ELSE_IF_STR_TO        ("RootChordL_inn"		,this->RootChordL_inn)
		ELSE_IF_STR_TO        ("TipRootRatio_inn"	,this->TipRootRatio_inn)
		ELSE_IF_STR_TO        ("TipRootRatio_out"	,this->TipRootRatio_out)
		ELSE_IF_STR_TO        ("SweepBackAngle_inn"	,this->SweepBackAngle_inn)
		ELSE_IF_STR_TO        ("SweepBackAngle_out"	,this->SweepBackAngle_out)
		ELSE_IF_STR_TO        ("Thickness_inn"		,this->Thickness_inn)
		ELSE_IF_STR_TO        ("NFaiU"				,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"				,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta_inn"			,this->NEta_inn)
		ELSE_IF_STR_TO        ("NEta_out"			,this->NEta_out)
		ELSE_IF_STR_TO        ("NHeight"			,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP1"				,this->BUPP1)
		ELSE_IF_STR_TO_VVECTOR("BLOW1"				,this->BLOW1)
		ELSE_IF_STR_TO_VVECTOR("DUPP1"				,this->DUPP1)
		ELSE_IF_STR_TO_VVECTOR("DLOW1"				,this->DLOW1)
		ELSE_IF_STR_TO_VVECTOR("BUPP2"				,this->BUPP2)
		ELSE_IF_STR_TO_VVECTOR("BLOW2"				,this->BLOW2)
		ELSE_IF_STR_TO_VVECTOR("DUPP2"				,this->DUPP2)
		ELSE_IF_STR_TO_VVECTOR("DLOW2"				,this->DLOW2)
		ELSE_IF_STR_TO        ("GridRefineType"		,this->GridRefineType)
		ELSE_IF_STR_TO        ("ifUseMid"			,this->ifUseMid)
	}
}

ShapeWing2::ShapeWing2() { ShapeType = 2; }

ShapeWing2::ShapeWing2(SW2& sw2)
{
	m_Wing2 = sw2;
	ShapeType = 2;
}

void ShapeWing2::BuildShape()
{
	cout << "Shape Wing2: 参数正在传入。。。" << endl;
	shared_ptr<CSTsurface> LeftInn(new CSTsurface);
	shared_ptr<CSTsurface> LeftOut(new CSTsurface);
	shared_ptr<CSTsurface> RighInn(new CSTsurface);
	shared_ptr<CSTsurface> RighOut(new CSTsurface);
	shared_ptr<CSTsurface> Middle_(new CSTsurface);
	
	(*LeftInn).Origin = m_Wing2.Origin;
	(*LeftInn).Origin(0) = m_Wing2.Origin[0] + m_Wing2.RootChordL_inn / 2.0;
	(*LeftInn).Rotation = { 0,90,180 };
	//(*LeftInn).Slope = { 0,0,m_Wing2.SweepBackAngle_inn };//原始数据
	(*LeftInn).Scale = { 1,0.8,0.8,m_Wing2.TipRootRatio_inn };//原始数据
	(*LeftInn).Slope = { 0,0,atan(tan(m_Wing2.SweepBackAngle_inn*pi/180)-(1-m_Wing2.TipRootRatio_inn)*m_Wing2.RootChordL_inn/(2*m_Wing2.SpanL_inn))*180/pi };//后掠角转换为物理意义
	(*LeftInn).Length = { m_Wing2.SpanL_inn,m_Wing2.Thickness_inn / 2.0,m_Wing2.Thickness_inn / 2.0,m_Wing2.RootChordL_inn };
	(*LeftInn).NS = m_Wing2.Ninner;
	(*LeftInn).NE = m_Wing2.Nmiddile;
	(*LeftInn).NFaiU = m_Wing2.NFaiU;
	(*LeftInn).NFaiL = m_Wing2.NFaiL;
	(*LeftInn).NEta = m_Wing2.NEta_inn;
	(*LeftInn).NHeight = m_Wing2.NHeight;
	(*LeftInn).Ratio = m_Wing2.Ratio_inn;
	(*LeftInn).M = { 0,0,0,0 };
	(*LeftInn).T = { 0,0,0,0 };
	(*LeftInn).Is_MeshFront = 0;
	(*LeftInn).Is_MeshBack = 0;
	copyM(m_Wing2.BUPP1, (*LeftInn).BUPP);
	copyM(m_Wing2.BLOW1, (*LeftInn).BLOW);
	copyM(m_Wing2.DUPP1, (*LeftInn).DUPP);
	copyM(m_Wing2.DLOW1, (*LeftInn).DLOW);

	(*LeftOut).Origin = { m_Wing2.Origin[0] + m_Wing2.SpanL_inn * tan(m_Wing2.SweepBackAngle_inn*pi/180) + m_Wing2.RootChordL_inn * m_Wing2.TipRootRatio_inn/2.0,
						   m_Wing2.Origin[1] + m_Wing2.SpanL_inn * sin(0) ,
						   m_Wing2.Origin[2] + m_Wing2.SpanL_inn * cos(0) };

	(*LeftOut).Rotation = { 0,90,180 };
	(*LeftOut).Slope = { 0,0,m_Wing2.SweepBackAngle_out };
	(*LeftOut).Slope = { 0,0,atan(tan(m_Wing2.SweepBackAngle_out*pi/180) - (1 - m_Wing2.TipRootRatio_out) * m_Wing2.RootChordL_inn*m_Wing2.TipRootRatio_inn / (2 * m_Wing2.SpanL_out))*180/pi };//测试数据
	(*LeftOut).Scale = { 1.0,0.6,0.6,m_Wing2.TipRootRatio_out };
	(*LeftOut).Length = { m_Wing2.SpanL_out,m_Wing2.Thickness_inn / 2.0 * 0.8, m_Wing2.Thickness_inn / 2.0 * 0.8, m_Wing2.RootChordL_inn * m_Wing2.TipRootRatio_inn };
	(*LeftOut).NS = m_Wing2.Nmiddile;
	(*LeftOut).NE = m_Wing2.Nouter;
	(*LeftOut).NFaiU = m_Wing2.NFaiU;
	(*LeftOut).NFaiL = m_Wing2.NFaiL;
	(*LeftOut).NEta = m_Wing2.NEta_out;
	(*LeftOut).NHeight = m_Wing2.NHeight;
	(*LeftOut).Ratio = m_Wing2.Ratio_out;
	(*LeftOut).M = { 0,0,0,0 };
	(*LeftOut).T = { 0,0,0,0 };
	(*LeftOut).Is_MeshFront = 0;
	(*LeftOut).Is_MeshBack = 0;
	copyM(m_Wing2.BUPP2, (*LeftOut).BUPP);
	copyM(m_Wing2.BLOW2, (*LeftOut).BLOW);
	copyM(m_Wing2.DUPP2, (*LeftOut).DUPP);
	copyM(m_Wing2.DLOW2, (*LeftOut).DLOW);

	(*RighInn) = (*LeftInn);
	(*RighInn).Origin[2] = -(*LeftInn).Origin[2];
	(*RighInn).Rotation = { 0,90,0 };

	(*RighOut) = (*LeftOut);
	(*RighOut).Origin[2] = -(*LeftOut).Origin[2];
	(*RighOut).Rotation = { 0,90,0 };

	if (m_Wing2.ifUseMid)
	{
		(*Middle_).Origin = m_Wing2.Origin;
		(*Middle_).Origin[0] = m_Wing2.Origin[0] + m_Wing2.RootChordL_inn / 2;
		//(*Middle_).Origin[2] = -(*Middle_).Origin[2];
		(*Middle_).Rotation = { 0,90,0 };
		(*Middle_).Slope = { 0,0,0 };
		(*Middle_).Scale = { 1,1,1,1 };
		(*Middle_).Length = { m_Wing2.Origin[2] * 2,m_Wing2.Thickness_inn / 2,m_Wing2.Thickness_inn / 2,m_Wing2.RootChordL_inn };
		(*Middle_).NS = m_Wing2.Ninner;
		(*Middle_).NE = m_Wing2.Ninner;
		(*Middle_).NFaiU = m_Wing2.NFaiU;
		(*Middle_).NFaiL = m_Wing2.NFaiL;
		(*Middle_).NEta = m_Wing2.NEta_inn;
		(*Middle_).NHeight = m_Wing2.NHeight;
		(*Middle_).Ratio = 1;
		(*Middle_).M = { 0,0,0,0 };
		(*Middle_).T = { 0,0,0,0 };
		(*Middle_).Is_MeshFront = 0;
		(*Middle_).Is_MeshBack = 0;
		vector<vector<double>> Temp = {{1}};
		copyM(Temp, (*Middle_).BUPP);
		copyM(Temp, (*Middle_).BLOW);
		copyM(Temp, (*Middle_).DUPP);
		copyM(Temp, (*Middle_).DLOW);
	}
	

	//CSTsf += LeftInn, LeftOut, RighInn, RighOut, Middle_;
	CSTsf.push_back(*LeftInn);
	CSTsf.push_back(*LeftOut);
	CSTsf.push_back(*RighInn);
	CSTsf.push_back(*RighOut);
	if (m_Wing2.ifUseMid)
	{
		CSTsf.push_back(*Middle_);
	}

	GridRefineType = m_Wing2.GridRefineType;

	cout << "Shape Wing2: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Wing2: 计算结束！" << endl;
}

vector<SingleBone> ShapeWing2::BuildBone()
{
	if (m_Wing2.ifUseMid)//两段连接为一体的（还没有加力！）
	{
		mat bone = flipud(CSTsf[1].MakeBone());

		mat pList0 = flipud(CSTsf[0].MakeBone());
		bone = join_cols(bone, pList0.rows(1, pList0.n_rows - 1));

		pList0 = CSTsf[4].MakeBone();
		bone = join_cols(bone, pList0.rows(1, pList0.n_rows - 1));

		pList0 = CSTsf[2].MakeBone();
		bone = join_cols(bone, pList0.rows(1, pList0.n_rows - 1));

		pList0 = CSTsf[3].MakeBone();
		bone = join_cols(bone, pList0.rows(1, pList0.n_rows - 1));

		//cout << bone << endl;
		vector<SingleBone> SB;
		SingleBone sb;
		sb.BoneP = bone;
		sb.BoneConnP_ID = bone.n_rows / 2;
		SB.push_back(sb);
		return SB;
	}
	else//单独两个机翼（已经加力21.5.8）
	{
		vector<SingleBone> SB;
		SingleBone sb;
		sb.BoneP = CSTsf[0].MakeBone(0.54);
		
		mat pList0 = (CSTsf[1].MakeBone());
		sb.BoneP = join_cols(sb.BoneP, pList0.rows(1, pList0.n_rows - 1));
		sb.BoneConnP_ID = 0;
		
		sb.BonePForce = CSTsf[0].ClacBonePForce();//force
		mat fList0 = CSTsf[1].ClacBonePForce();//force
		if (fList0.n_rows != 0)
		{
			sb.BonePForce.row(sb.BonePForce.n_rows - 1) += fList0.row(0);//force
			sb.BonePForce = join_cols(sb.BonePForce, fList0.rows(1, fList0.n_rows - 1));//force
		}

		SB.push_back(sb);
		//-------------------------

		SingleBone sb2;
		sb2.BoneP = CSTsf[2].MakeBone(0.54);
		
		mat pList1 = (CSTsf[3].MakeBone());
		sb2.BoneP = join_cols(sb2.BoneP, pList1.rows(1, pList0.n_rows - 1));
		sb2.BoneConnP_ID = 0;

		sb2.BonePForce = CSTsf[2].ClacBonePForce();//force
		mat fList1 = CSTsf[3].ClacBonePForce();//force
		if (fList1.n_rows != 0)
		{
			sb2.BonePForce.row(sb2.BonePForce.n_rows - 1) += fList1.row(0);//force
			sb2.BonePForce = join_cols(sb2.BonePForce, fList1.rows(1, fList1.n_rows - 1));//force
		}

		SB.push_back(sb2);

		return SB;
	}
}

vector<ShellandBeam> ShapeWing2::BuildStrc1D()
{
	vector<ShellandBeam> A(2);

	mat X = linspace(0, 1, 5);
	mat Z = linspace(0, 1, 5);
	CSTsf[0].MakeStruct1D(X, Z);
	CSTsf[1].MakeStruct1D(X, Z);
	CSTsf[2].MakeStruct1D(X, Z);
	CSTsf[3].MakeStruct1D(X, Z);
	//----------左机翼----------
	uniqueMat left(join_cols(CSTsf[0].aero_strc.node, CSTsf[1].aero_strc.node));
	mat e_aero = join_cols(CSTsf[0].aero_strc.elem_aero, CSTsf[1].aero_strc.elem_aero + CSTsf[0].aero_strc.node.n_rows);
	CSTsf[1].aero_strc.elem_strc.cols(0, 1) += CSTsf[0].aero_strc.node.n_rows;
	mat e_strc = join_cols(CSTsf[0].aero_strc.elem_strc, CSTsf[1].aero_strc.elem_strc);
	e_strc = left.E_new(e_strc);

	//后半段的垂直于z的第一个梁与前面重复，去掉
	int deleID0 = 0,deleID1 = 0;
	for (size_t i = e_strc.n_rows / 2 - 1; i < e_strc.n_rows; i++)
	{
		if (e_strc(i, 2) == 1)//xBeam
		{
			if (e_strc(i - 1, 2) == 0)
			{
				deleID0 = i;
			}
			if (e_strc(i + 1, 3) == 1)
			{
				deleID1 = i;
				break;
			}
		}
	}


	e_strc = join_cols(e_strc.rows(0, deleID0 - 1), e_strc.rows(deleID1 + 1, e_strc.n_rows - 1));
	A[0].name = "SW2_left";
	A[0].node = left.P_new();
	A[0].elem_aero = left.E_new(e_aero);
	A[0].elem_strc = e_strc;
	A[0].elem_h = join_cols(CSTsf[0].aero_strc.elem_h, CSTsf[1].aero_strc.elem_h);
	A[0].LinkPtId = CSTsf[0].aero_strc.LinkPtId;
	A[0].node_f = left.renewF(join_cols(CSTsf[0].aero_strc.node_f, CSTsf[1].aero_strc.node_f));
	//----------右机翼----------
	uniqueMat righ(join_cols(CSTsf[2].aero_strc.node, CSTsf[3].aero_strc.node));
	mat e_aero2 = join_cols(CSTsf[2].aero_strc.elem_aero, CSTsf[3].aero_strc.elem_aero + CSTsf[2].aero_strc.node.n_rows);
	CSTsf[3].aero_strc.elem_strc.cols(0, 1) += CSTsf[2].aero_strc.node.n_rows;
	mat e_strc2 = join_cols(CSTsf[2].aero_strc.elem_strc, CSTsf[3].aero_strc.elem_strc);
	e_strc2 = righ.E_new(e_strc2);

	//后半段的垂直于z的第一个梁与前面重复，去掉
	deleID0 = deleID1 = 0;
	for (size_t i = e_strc2.n_rows / 2 - 1; i < e_strc2.n_rows; i++)
	{
		if (e_strc2(i, 2) == 1)//xBeam
		{
			if (e_strc2(i - 1, 2) == 0)
			{
				deleID0 = i;
			}
			if (e_strc2(i + 1, 3) == 1)
			{
				deleID1 = i;
				break;
			}
		}
	}
	e_strc2 = join_cols(e_strc2.rows(0, deleID0 - 1), e_strc2.rows(deleID1 + 1, e_strc2.n_rows - 1));
	A[1].name = "SW2_righ";
	A[1].node = righ.P_new();
	A[1].elem_aero = righ.E_new(e_aero2);
	A[1].elem_strc = e_strc2;
	A[1].elem_h = join_cols(CSTsf[2].aero_strc.elem_h, CSTsf[3].aero_strc.elem_h);
	A[1].LinkPtId = CSTsf[2].aero_strc.LinkPtId;
	A[1].node_f = righ.renewF(join_cols(CSTsf[2].aero_strc.node_f, CSTsf[3].aero_strc.node_f));
	
	return A;
}

//---------------------------------------------------------------
//--------------------------ShapeTail1---------------------------
//初始化ShapeTail1参数
ST1::ST1()
{
	this->Origin = { 15000, 1574, 0 };
	this->Ninner = { 0.9, 1, 0.9, 1 };
	this->Nouter = { 0.9, 1, 0.9, 1 };
	this->Ratio1 = 1;
	this->SpanL = 5000;
	this->RootChordL = 7200;
	this->TipRootRatio = 0.48;
	this->SweepBackAngle = 40;
	this->Thickness = 500;
	//this->SideAngle = 0;
	this->NFaiU = 20;
	this->NFaiL = 20;
	this->NEta = 20;
	this->NHeight = 3;
	this->BUPP = { {1} };
	this->BLOW = { {1} };
	this->DUPP = { {1} };
	this->DLOW = { {1} };
	this->GridRefineType = -1;
}

ST1::ST1(
	vector<double>& Origin,
	vector<double>& Ninner, vector<double>& Nouter,
	double& Ratio1,
	double& SpanL, double& RootChordL, double& TipRootRatio, double& SweepBackAngle, double& Thickness,
	int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
	vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
	vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,
	int& GridRefineType)
{
	this->Origin = Origin;
	this->Ninner = Ninner;
	this->Nouter = Nouter;
	this->Ratio1 = Ratio1;
	this->SpanL = SpanL;
	this->RootChordL = RootChordL;
	this->TipRootRatio = TipRootRatio;
	this->SweepBackAngle = SweepBackAngle;
	this->Thickness = Thickness;
	//this->SideAngle = 0;
	this->NFaiU = NFaiU;
	this->NFaiL = NFaiL;
	this->NEta = NEta;
	this->NHeight = NHeight;
	this->BUPP = BUPP;
	this->BLOW = BLOW;
	this->DUPP = DUPP;
	this->DLOW = DLOW;
	this->GridRefineType = GridRefineType;
}

void ST1::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"        ,this->Origin)
		ELSE_IF_STR_TO_VECTOR ("Ninner"        ,this->Ninner)
		ELSE_IF_STR_TO_VECTOR ("Nouter"        ,this->Nouter)
		ELSE_IF_STR_TO        ("Ratio1"        ,this->Ratio1)
		ELSE_IF_STR_TO        ("SpanL"         ,this->SpanL)
		ELSE_IF_STR_TO        ("RootChordL"    ,this->RootChordL)
		ELSE_IF_STR_TO        ("TipRootRatio"  ,this->TipRootRatio)
		ELSE_IF_STR_TO        ("SweepBackAngle",this->SweepBackAngle)
		ELSE_IF_STR_TO        ("Thickness"     ,this->Thickness)
		ELSE_IF_STR_TO        ("NFaiU"         ,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"         ,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"          ,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"       ,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP"          ,this->BUPP)
		ELSE_IF_STR_TO_VVECTOR("BLOW"          ,this->BLOW)
		ELSE_IF_STR_TO_VVECTOR("DUPP"          ,this->DUPP)
		ELSE_IF_STR_TO_VVECTOR("DLOW"          ,this->DLOW)
		ELSE_IF_STR_TO        ("GridRefineType",this->GridRefineType)
	}
}

ShapeTail1::ShapeTail1() { ShapeType = 2; }
vector<ShellandBeam> ShapeTail1::BuildStrc1D()
{
	vector<ShellandBeam> A(1);

	mat X = linspace(0, 1, 5);
	mat Z = linspace(0, 1, 5);
	CSTsf[0].MakeStruct1D(X, Z);
	
	A[0] = CSTsf[0].aero_strc;
	A[0].name = "ST1";

	return A;
}
ShapeTail1::ShapeTail1(ST1& st1)
{
	m_Tail1 = st1;
	ShapeType = 2;
}


void ShapeTail1::BuildShape()
{
	cout << "Shape Tail1: 参数正在传入。。。" << endl;
	shared_ptr<CSTsurface> Tail(new CSTsurface);

	(*Tail).Origin = m_Tail1.Origin;
	(*Tail).Origin[0] += m_Tail1.RootChordL/2.0;
	(*Tail).Rotation = { 0,90,90 };
	(*Tail).Slope = { 0,0,atan(tan(m_Tail1.SweepBackAngle * pi / 180) - (1 - m_Tail1.TipRootRatio) * m_Tail1.RootChordL / (2 * m_Tail1.SpanL)) * 180 / pi };
	(*Tail).Scale = { 1,m_Tail1.TipRootRatio,m_Tail1.TipRootRatio,m_Tail1.TipRootRatio };
	(*Tail).Length = { m_Tail1.SpanL, m_Tail1.Thickness / 2.0, m_Tail1.Thickness / 2.0, m_Tail1.RootChordL };
	(*Tail).NS = m_Tail1.Ninner;
	(*Tail).NE = m_Tail1.Nouter;
	(*Tail).NFaiU = m_Tail1.NFaiU;
	(*Tail).NFaiL = m_Tail1.NFaiL;
	(*Tail).NEta = m_Tail1.NEta;
	(*Tail).NHeight = m_Tail1.NHeight;
	(*Tail).Ratio = m_Tail1.Ratio1;
	(*Tail).M = { 0,0,0,0 };
	(*Tail).T = { 0,0,0,0 };
	(*Tail).Is_MeshFront = 0;
	(*Tail).Is_MeshBack = 0;// = 1;
	copyM(m_Tail1.BUPP, (*Tail).BUPP);
	copyM(m_Tail1.BLOW, (*Tail).BLOW);
	copyM(m_Tail1.DUPP, (*Tail).DUPP);
	copyM(m_Tail1.DLOW, (*Tail).DLOW);

	CSTsf.push_back(*Tail);
	GridRefineType = m_Tail1.GridRefineType;

	cout << "Shape Tail1: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Tail1: 计算结束！" << endl;
}

vector<SingleBone> ShapeTail1::BuildBone()//（已经加力2021.5.8）
{
	vector<SingleBone> SB;
	SingleBone sb;
	sb.BoneP = CSTsf[0].MakeBone();
	sb.BonePForce = CSTsf[0].ClacBonePForce();
	sb.BoneConnP_ID = 0;
	SB.push_back(sb);
	return SB;
}

//---------------------------------------------------------------
//--------------------------ShapeTail2---------------------------
//初始化ShapeTail2参数
ST2::ST2()
{
	this->Origin = { 15000, 1574, 500 };
	this->Ninner = { 0.9, 1, 0.9, 1 };
	this->Nouter = { 0.9, 1, 0.9, 1 };
	this->Ratio1 = 1;
	this->SpanL = 5000;
	this->RootChordL = 7200;
	this->TipRootRatio = 0.48;
	this->SweepBackAngle = 40;
	this->Thickness = 500;
	this->SideAngle = 70;
	this->NFaiU = 20;
	this->NFaiL = 20;
	this->NEta = 20;
	this->NHeight = 3;
	this->BUPP = { {1} };
	this->BLOW = { {1} };
	this->DUPP = { {1} };
	this->DLOW = { {1} };
	this->GridRefineType = -1;
}

ST2::ST2(
	vector<double>& Origin,
	vector<double>& Ninner, vector<double>& Nouter,
	double& Ratio1,
	double& SpanL, double& RootChordL, double& TipRootRatio,
	double& SweepBackAngle, double& Thickness, double& SideAngel,
	int& NFaiU, int& NFaiL, int& NEta, int& NHeight,
	vector<vector<double>>& BUPP, vector<vector<double>>& BLOW,
	vector<vector<double>>& DUPP, vector<vector<double>>& DLOW,
	int& GridRefineType)
{
	this->Origin = Origin;
	this->Ninner = Ninner;
	this->Nouter = Nouter;
	this->Ratio1 = Ratio1;
	this->SpanL = SpanL;
	this->RootChordL = RootChordL;
	this->TipRootRatio = TipRootRatio;
	this->SweepBackAngle = SweepBackAngle;
	this->Thickness = Thickness;
	this->SideAngle = SideAngel;
	this->NFaiU = NFaiU;
	this->NFaiL = NFaiL;
	this->NEta = NEta;
	this->NHeight = NHeight;
	this->BUPP = BUPP;
	this->BLOW = BLOW;
	this->DUPP = DUPP;
	this->DLOW = DLOW;
	this->GridRefineType = GridRefineType;
}

void ST2::readFromFile(ifstream& ifs)
{
	string str_line;
	while (getline(ifs, str_line) && str_line.compare("*End"))//读到"*End"时停止该部分
	{
		str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//清除空格
		stringstream ss(str_line);
		string str_tmp;//存储下方判断所使用的一行数据
		getline(ss, str_tmp, ',');//获取首个string
		if (str_tmp[0] == '$') {}//跳过
		ELSE_IF_STR_TO_VECTOR ("Origin"			,this-> Origin)
		ELSE_IF_STR_TO_VECTOR ("Ninner"			,this-> Ninner)
		ELSE_IF_STR_TO_VECTOR ("Nouter"			,this-> Nouter)
		ELSE_IF_STR_TO        ("Ratio1"			,this->Ratio1)
		ELSE_IF_STR_TO        ("SpanL"			,this->SpanL)
		ELSE_IF_STR_TO        ("RootChordL"		,this->RootChordL)
		ELSE_IF_STR_TO        ("TipRootRatio"	,this->TipRootRatio)
		ELSE_IF_STR_TO        ("SweepBackAngle"	,this->SweepBackAngle)
		ELSE_IF_STR_TO        ("SideAngle"		,this->SideAngle)
		ELSE_IF_STR_TO        ("Thickness"		,this->Thickness)
		ELSE_IF_STR_TO        ("NFaiU"			,this->NFaiU)
		ELSE_IF_STR_TO        ("NFaiL"			,this->NFaiL)
		ELSE_IF_STR_TO        ("NEta"			,this->NEta)
		ELSE_IF_STR_TO        ("NHeight"		,this->NHeight)
		ELSE_IF_STR_TO_VVECTOR("BUPP"			,this->BUPP)
		ELSE_IF_STR_TO_VVECTOR("BLOW"			,this->BLOW)
		ELSE_IF_STR_TO_VVECTOR("DUPP"			,this->DUPP)
		ELSE_IF_STR_TO_VVECTOR("DLOW"			,this->DLOW)
		ELSE_IF_STR_TO        ("GridRefineType"	,this->GridRefineType)
	}
}

ShapeTail2::ShapeTail2() { ShapeType = 2; }
ShapeTail2::ShapeTail2(ST2& st2)
{
	m_Tail2 = st2;
	ShapeType = 2;
}

vector<ShellandBeam> ShapeTail2::BuildStrc1D()
{
	vector<ShellandBeam> A(2);

	mat X = linspace(0, 1, 5);
	mat Z = linspace(0, 1, 5);
	CSTsf[0].MakeStruct1D(X, Z);
	CSTsf[1].MakeStruct1D(X, Z);

	A[0] = CSTsf[0].aero_strc;
	A[0].name = "ST1_left";
	A[1] = CSTsf[1].aero_strc;
	A[1].name = "ST1_righ";

	return A;
}

void ShapeTail2::BuildShape()
{
	cout << "Shape Tail2: 参数正在传入。。。" << endl;
	shared_ptr<CSTsurface> Tail(new CSTsurface);
	shared_ptr<CSTsurface> Tail2(new CSTsurface);

	(*Tail).Origin = m_Tail2.Origin;
	(*Tail).Origin[0] += m_Tail2.RootChordL / 2.0;
	(*Tail).Origin[2] = -m_Tail2.Origin[2];
	(*Tail).Rotation = { 0,90,m_Tail2.SideAngle };
	(*Tail).Slope = { 0,0,atan(tan(m_Tail2.SweepBackAngle * pi / 180) - (1 - m_Tail2.TipRootRatio) * m_Tail2.RootChordL / (2 * m_Tail2.SpanL)) * 180 / pi };
	(*Tail).Scale = { 1,m_Tail2.TipRootRatio,m_Tail2.TipRootRatio,m_Tail2.TipRootRatio };
	(*Tail).Length = { m_Tail2.SpanL, m_Tail2.Thickness / 2.0, m_Tail2.Thickness / 2.0, m_Tail2.RootChordL };
	(*Tail).NS = m_Tail2.Ninner;
	(*Tail).NE = m_Tail2.Nouter;
	(*Tail).NFaiU = m_Tail2.NFaiU;
	(*Tail).NFaiL = m_Tail2.NFaiL;
	(*Tail).NEta = m_Tail2.NEta;
	(*Tail).NHeight = m_Tail2.NHeight;
	(*Tail).Ratio = m_Tail2.Ratio1;
	(*Tail).M = { 0,0,0,0 };
	(*Tail).T = { 0,0,0,0 };
	(*Tail).Is_MeshFront = 0;
	(*Tail).Is_MeshBack = 0;// = 1;
	copyM(m_Tail2.BUPP, (*Tail).BUPP);
	copyM(m_Tail2.BLOW, (*Tail).BLOW);
	copyM(m_Tail2.DUPP, (*Tail).DUPP);
	copyM(m_Tail2.DLOW, (*Tail).DLOW);

	(*Tail2) = (*Tail);
	(*Tail2).Origin[2] = m_Tail2.Origin[2];
	(*Tail2).Rotation = { 0,-90,m_Tail2.SideAngle };
	(*Tail2).Slope[2] = -(*Tail).Slope[2];
	
	CSTsf.push_back(*Tail);
	CSTsf.push_back(*Tail2);

	GridRefineType = m_Tail2.GridRefineType;

	cout << "Shape Tail2: 开始运行！" << endl;
	BacsShapeCSTs();
	cout << "Shape Tail2: 计算结束！" << endl;
	
}

vector<SingleBone> ShapeTail2::BuildBone()//（已经加力2021.5.8）
{
	vector<SingleBone> SB;
	SingleBone sb;
	sb.BoneP = CSTsf[0].MakeBone();
	sb.BonePForce = CSTsf[0].ClacBonePForce();
	sb.BoneConnP_ID = 0;
	SB.push_back(sb);

	
	SingleBone sb2;
	sb2.BoneP = CSTsf[1].MakeBone();
	sb2.BonePForce = CSTsf[1].ClacBonePForce();
	sb2.BoneConnP_ID = 0;
	//cout << sb2.StructBone(sb.BoneP) << endl;
	SB.push_back(sb2);
	
	return SB;
}