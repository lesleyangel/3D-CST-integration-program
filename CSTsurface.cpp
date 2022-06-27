#include <iostream>
#include "include/armadillo"
#include "factorial.h"
#include "CSTsurface.h"
#include "Struct2D.h"
using namespace arma;
using namespace std;

#define pi 3.1415926

//网格其他信息类
MeshInfo::MeshInfo()
{
	SurfaceAera = 0;
	Volume = 0;
	planforarea = 0;
}

//单个cst曲面类
//初始化参数
CSTsurface::CSTsurface()
{
	Origin = zeros(1, 3);
	Rotation = zeros(1, 3);
	Length = { 3000,400,150,1000 };
	Slope = zeros(1, 3);
	Scale = { 1,1,1,1 };
	//class
	NS = { 0.5,0.5,5,5 };
	NE = { 0.5,0.5,5,5 };
	M = { 0.5,0,0.5,0 };
	T = { 0.5,0,0.5,0 };
	Ratio = 1;

	BUPP = { 1 };	BLOW = { 1 };
	DUPP = { 1 };	DLOW = { 1 };

	BUPP = { {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7},
			 {1,2,3,4,5,6,7} };

	MeshRefineRatio = { 
		{81,81,9,9, 9,45,9, 9,9,81},
		{ 9, 9,1,1, 1, 5,1, 1,1, 9},
		{ 9, 9,1,1, 1, 5,1, 1,1, 9},
		{ 9, 9,1,1, 1, 5,1, 1,1, 9},
		{ 9, 9,1,1, 1, 5,1, 1,1, 9},
		{81,81,9,9, 9,45,9, 9,9,81} };
	//MeshRefineRatio = {
	//	{81,81,9,9, 9,45,9, 9,9,81},
	//	{ 9, 9,1,1, 1, 5,1, 1,1, 9},
	//	{ 9, 9,1,1, 9, 5,1, 1,1, 9},
	//	{ 9, 9,1,1, 9, 5,1, 1,1, 9},
	//	{ 9, 9,1,1, 1, 5,1, 1,1, 9},
	//	{81,81,9,9, 9,45,9, 9,9,81} };
	
	MeshRefineRatio = {
	{81,36,9,9,9,81},
	{ 9, 4,1,1,1, 9},
	{ 9, 4,1,1,1, 9},
	{ 9, 4,1,1,1, 9},
	{ 9, 4,1,1,1, 9},
	{81,36,9,9,9,81} };
	MeshRefineRatio = {
	{ 13, 5,1,1,1, 3},
	{ 13, 5,1,1,1, 3},
	{ 13, 5,1,1,1, 3},
	{ 13, 5,1,1,1, 3},};
	MeshRefineRatio = 1;
	PriFuncType = 0;// 1;

	Is_MeshFront = 0;	Is_MeshBack = 0;
	NFaiU = 40;	NFaiL = 40;	NEta = 60;	NHeight = 20;

	//网格数修正
	MeshNumCorrection();

	//生成网格列表：FaiU、FaiL、Eta
	GetNet();

	BonePoint.clear();
}

//计算生成CST曲面网格信息
void CSTsurface::CST3D()
{
	//参数传入
	double OriginX = Origin(0);
	double OriginY = Origin(1);
	double OriginZ = Origin(2);

	double RotationX = Rotation(0) / 180.0 * pi;
	double RotationY = Rotation(1) / 180.0 * pi;
	double RotationZ = Rotation(2) / 180.0 * pi;

	double Xlength = Length(0);
	double Ylength_upp = Length(1);
	double Ylength_low = Length(2);
	double Zlength = Length(3);

	double SlopeX = Slope(0) / 180.0 * pi;
	double SlopeY = Slope(1) / 180.0 * pi;
	double SlopeZ = Slope(2) / 180.0 * pi;

	double ScaleX = Scale(0);
	double ScaleYUPP = Scale(1);
	double ScaleYLOW = Scale(2);
	double ScaleZ = Scale(3);

	//类别函数控制参数
	double N1s = NS(0);	double N2s = NS(1);	double N3s = NS(2);	double N4s = NS(3);	//截面
	double N1e = NE(0);	double N2e = NE(1);	double N3e = NE(2);	double N4e = NE(3);	//截面
	double M1 = M(0);	double M2 = M(1);	double M3 = M(2);	double M4 = M(3);	//轴线
	double T1 = T(0);	double T2 = T(1);	double T3 = T(2);	double T4 = T(3);	//侧向

	//生成网格列表：FaiU、FaiL、Eta
	GetNet();

	N1mUPP = zeros(1, EtaU.n_cols);
	N2mUPP = zeros(1, EtaU.n_cols);
	N1mLOW = zeros(1, EtaU.n_cols);
	N2mLOW = zeros(1, EtaU.n_cols);
	//CST3D
	PUPP = SketchCST(FaiU, EtaU, N1s, N2s, N1e, N2e, M1, M2, T1, T2, BUPP, DUPP, Ratio, N1mUPP, N2mUPP);
	PLOW = SketchCST(FaiL, EtaL, N3s, N4s, N3e, N4e, M3, M4, T3, T4, BLOW, DLOW, Ratio, N1mLOW, N2mLOW);
	
	
	//网格修正
	//GridRefine();

	//Scale X Y Z
	//还原为标准尺寸 normal (0)
	PUPP.X = PUPP.X * Xlength;
	PUPP.Y = PUPP.Y * Ylength_upp;
	PUPP.Z = PUPP.Z * Zlength / 2.0;
	PLOW.X = PLOW.X * Xlength;
	PLOW.Y = PLOW.Y * Ylength_low;
	PLOW.Z = PLOW.Z * Zlength / 2.0;
	
	//还原比例刻度 （用于翼的跟梢比）scale (1)
	Eta = trans(linspace(0, 1, NEta));
	for (int j = 0; j < NEta; j++)
	{
		for (int i = 0; i < NFaiU; i++)
		{
			PUPP.X(i, j) = PUPP.X(i, j) * ((ScaleX - 1) * EtaU(i, j) + 1);
			PUPP.Y(i, j) = PUPP.Y(i, j) * ((ScaleYUPP - 1) * EtaU(i, j) + 1);
			PUPP.Z(i, j) = PUPP.Z(i, j) * ((ScaleZ - 1) * EtaU(i, j) + 1);
		}
		//PUPP.X.col(i) = PUPP.X.col(i) * ((ScaleX    - 1) * Eta(i) + 1);
		//PUPP.Y.col(i) = PUPP.Y.col(i) * ((ScaleYUPP - 1) * Eta(i) + 1);
		//PUPP.Z.col(i) = PUPP.Z.col(i) * ((ScaleZ    - 1) * Eta(i) + 1);
		for (int i = 0; i < NFaiL; i++)
		{
			PLOW.X(i, j) = PLOW.X(i, j) * ((ScaleX - 1) * EtaL(i, j) + 1);
			PLOW.Y(i, j) = PLOW.Y(i, j) * ((ScaleYLOW - 1) * EtaL(i, j) + 1);
			PLOW.Z(i, j) = PLOW.Z(i, j) * ((ScaleZ - 1) * EtaL(i, j) + 1);
		}
		//PLOW.X.col(j) = PLOW.X.col(j) * ((ScaleX    - 1) * Eta(j) + 1);
		//PLOW.Y.col(j) = PLOW.Y.col(j) * ((ScaleYLOW - 1) * Eta(j) + 1);
		//PLOW.Z.col(j) = PLOW.Z.col(j) * ((ScaleZ    - 1) * Eta(j) + 1);
	}
	
	//还原图形斜率（用于后掠角） slope (2)
	if (SlopeY != 0)
	{
		for (size_t i = 0; i < PUPP.Y.n_rows; i++)
		{
			for (size_t j = 0; j < PUPP.Y.n_cols; j++)
			{
				PUPP.Y(i, j) = PUPP.Y(i, j) + PUPP.X(i, j) * tan(SlopeY);
			}
		}
		for (size_t i = 0; i < PLOW.Y.n_rows; i++)
		{
			for (size_t j = 0; j < PLOW.Y.n_cols; j++)
			{
				PLOW.Y(i, j) = PLOW.Y(i, j) + PLOW.X(i, j) * tan(SlopeY);
			}
		}
	}
	if (SlopeZ != 0)
	{
		for (size_t i = 0; i < PUPP.Z.n_rows; i++)
		{
			for (size_t j = 0; j < PUPP.Z.n_cols; j++)
			{
				//mat Temp = repmat(PUPP.Z.row(1), NFaiU, 1);//??
				PUPP.Z(i, j) = PUPP.Z(i, j) + PUPP.X(i, j) * tan(SlopeZ);// -Temp(j) - Zlength / 2;
			}
		}
	
		for (size_t i = 0; i < PLOW.Z.n_rows; i++)
		{
			for (size_t j = 0; j < PLOW.Z.n_cols; j++)
			{
				//mat Temp = repmat(PLOW.Z.row(1), NFaiL, 1);//??
				PLOW.Z(i, j) = PLOW.Z(i, j) + PLOW.X(i, j) * tan(SlopeZ);// -Temp(j) - Zlength / 2;
			}
		}
	}

	//评估网格信息
	EvalGrid();

	//参数点数据设置
	mat Coord1 = join_rows(reshape(trans(PUPP.X), NFaiU * NEta, 1), reshape(trans(PUPP.Y), NFaiU * NEta, 1));
	Coord1 = join_rows(Coord1, reshape(trans(PUPP.Z), NFaiU * NEta, 1));
	mat Coord2 = join_rows(reshape(trans(PLOW.X), NFaiL * NEta, 1), reshape(trans(PLOW.Y), NFaiL * NEta, 1));
	Coord2 = join_rows(Coord2, reshape(trans(PLOW.Z), NFaiL * NEta, 1));

	//还原角度 通过Rotation Matrix
	mat Lrotation = { { cos(RotationZ) * cos(RotationY) , sin(RotationZ) , -cos(RotationZ) * sin(RotationY) },
					 {-sin(RotationZ) * cos(RotationY) * cos(RotationX) + sin(RotationY) * sin(RotationX), cos(RotationZ) * cos(RotationX)	, sin(RotationZ) * sin(RotationY) * cos(RotationX) + cos(RotationY) * sin(RotationX)},
					 { sin(RotationZ) * cos(RotationY) * sin(RotationX) + sin(RotationY) * cos(RotationX),-cos(RotationZ) * sin(RotationX)	,-sin(RotationZ) * sin(RotationY) * sin(RotationX) + cos(RotationY) * cos(RotationX)} };

	mat RCoord1 = Coord1 * Lrotation;
	RCoord1.col(0) = RCoord1.col(0) + OriginX;
	RCoord1.col(1) = RCoord1.col(1) + OriginY;
	RCoord1.col(2) = RCoord1.col(2) + OriginZ;

	Coord2.col(1) = -Coord2.col(1);//下侧加负号->反方向
	mat RCoord2 = Coord2 * Lrotation;
	RCoord2.col(0) = RCoord2.col(0) + OriginX;
	RCoord2.col(1) = RCoord2.col(1) + OriginY;
	RCoord2.col(2) = RCoord2.col(2) + OriginZ;

	//计算返回节点和单元参数
	//上表面节点P、单元E
	GridUpp.P = RCoord1;
	GridUpp.E = zeros((NFaiU - 1) * (NEta - 1), 4);
	int id = 0;
	for (int i = 1; i < NFaiU; i++)
	{
		for (int j = 1; j < NEta; j++)
		{
			GridUpp.E(id, 3) = (double)NEta * ((double)i - 1) + j;
			GridUpp.E(id, 2) = (double)NEta * ((double)i - 1) + j + 1;
			GridUpp.E(id, 1) = (double)NEta * i + j + 1;
			GridUpp.E(id, 0) = (double)NEta * i + j;
			id++;
		}
	}
	//下表面节点P、单元E
	GridLow.P = RCoord2;
	GridLow.E = zeros((NFaiL - 1) * (NEta - 1), 4);
	id = 0;
	for (int i = 1; i < NFaiL; i++)
	{
		for (int j = 1; j < NEta; j++)
		{
			GridLow.E(id, 3) = (double)NEta * ((double)i - 1) + j;
			GridLow.E(id, 2) = (double)NEta * ((double)i - 1) + j + 1;
			GridLow.E(id, 1) = (double)NEta * i + j + 1;
			GridLow.E(id, 0) = (double)NEta * i + j;
			id++;
		}
	}
}

//计算整体正则坐标 参数代表调用修正算法次数
void CSTsurface::RefineMesh(const int refinetime, bool i_x_free/* = false*/)
{
	if (refinetime < 1)//如果设置的网格优化次数小于1 那么直接跳出该函数
	{
		return;
	}
	int row_upp = PUPP.X.n_rows;
	int col_upp = PUPP.X.n_cols;
	
	if (i_x_free)
	{
		//第一次使用带幂次的网格修正算法
		if (refinetime > 0)
		{
			GetLengthMesh("upp", true);
			GetLengthMesh("low", true);
			//GetLengthMash_x(true);
			CST3D();
		}

		double vol = Info.Volume;
		//后几次使用不带幂次的网格修正算法
		if (refinetime > 1)
		{
			for (int i = 1; i < refinetime; i++)
			{
				GetLengthMesh("upp", false, true);
				GetLengthMesh("low", false, true);
				//GetLengthMash_x(false);
				CST3D();
				vol = Info.Volume;
			}
		}
	}
	else
	{
		//第一次使用带幂次的网格修正算法
		if (refinetime > 0)
		{
			GetLengthMesh_z("upp", true);
			GetLengthMesh_z("low", true);
			GetLengthMash_x(true);
			CST3D();
		}

		double vol = Info.Volume;
		//后几次使用不带幂次的网格修正算法
		if (refinetime > 1)
		{
			for (int i = 1; i < refinetime; i++)
			{
				GetLengthMesh_z("upp", false);
				GetLengthMesh_z("low", false);
				GetLengthMash_x(false);
				CST3D();
				vol = Info.Volume;
			}
		}
	}

}

//计算正则坐标
void CSTsurface::GetLengthMesh(string uporlow, bool ifpow, bool ifrefine)
{
	//PUPP.printP();
	PointSet PUPP; mat FaiU; mat EtaU; double T; mat N1m; mat N2m;
	//修正上表面
	if (!uporlow.compare("upp"))
	{
		PUPP = this->PUPP;
		FaiU = this->FaiU;
		EtaU = this->EtaU;
		T = this->T[0];
		N1m = this->N1mUPP;
		N2m = this->N2mUPP;
	}
	//修正下表面
	else
	{
		PUPP = this->PLOW;
		FaiU = this->FaiL;
		EtaU = this->EtaL;
		T = this->T[2];
		N1m = this->N1mLOW;
		N2m = this->N2mLOW;
	}
	int row_upp = PUPP.X.n_rows;
	int col_upp = PUPP.X.n_cols;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	mat Supp_x = zeros(row_upp, col_upp);//沿着x方向的距离矩阵
	mat Supp_x_total = zeros(row_upp, 1);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	mat Supp_z = zeros(row_upp, col_upp);//沿着z方向的距离矩阵
	mat Supp_z_total = zeros(1, col_upp);

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	mat LinealMeasure = MakeS(1, 3, trans(MeshRefineRatio), FaiU, EtaU);
	// 
	//计算每两个点x方向之间的距离
	for (int i = 0; i < row_upp; i++)
	{
		for (int j = 1; j < col_upp; j++)
		{
			Supp_x(i, j) = sqrt(
				pow(PUPP.X(i, j - 1) - PUPP.X(i, j), 2) +
				pow(PUPP.Y(i, j - 1) - PUPP.Y(i, j), 2) +
				pow(PUPP.Z(i, j - 1) - PUPP.Z(i, j), 2));
			if (ifrefine)			
				Supp_x(i, j) *= LinealMeasure(i, j);
			Supp_x_total(i, 0) += Supp_x(i, j);
		}
	}

	for (int j = 1; j < col_upp; j++)//距离求出正则坐标
	{
		Supp_x.col(j) += Supp_x.col(j - 1);
	}
	for (int i = 0; i < row_upp; i++)//正则坐标归一化
	{
		if (Supp_x(i, col_upp - 1) != 0)
		{
			Supp_x.row(i) /= Supp_x(i, col_upp - 1);
		}
	}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//计算每两个点z方向之间的距离
	for (int i = 1; i < row_upp; i++)
	{
		for (int j = 0; j < col_upp; j++)
		{
			Supp_z(i, j) = sqrt(
				pow(PUPP.X(i - 1, j) - PUPP.X(i, j), 2) +
				pow(PUPP.Y(i - 1, j) - PUPP.Y(i, j), 2) +
				pow(PUPP.Z(i - 1, j) - PUPP.Z(i, j), 2));
			if (ifrefine)			
				Supp_z(i, j) *= LinealMeasure(i, j);
			Supp_z_total(0, j) += Supp_z(i, j);
		}
	}
	for (int i = 1; i < row_upp; i++)//距离求出正则坐标
	{
		Supp_z.row(i) += Supp_z.row(i - 1);
	}
	for (int j = 0; j < col_upp; j++)//正则坐标归一化
	{
		if (Supp_z(row_upp - 1, j) != 0)
		{
			Supp_z.col(j) /= Supp_z(row_upp - 1, j);
		}
	}

	//初始化EtaU FaiU
	mat FaiU_new = zeros(row_upp, col_upp);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	mat EtaU_new = zeros(row_upp, col_upp);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	for (int i = 0; i < row_upp; i++)
	{
		for (int j = 0; j < col_upp; j++)
		{
			FaiU_new(i, j) = (double)i / ((double)row_upp - 1);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			EtaU_new(i, j) = (double)j / ((double)col_upp - 1);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}
	}

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//修正EtaU FaiU 每个网格单独计算
	for (int i = 0; i < row_upp; i++)
	{
		for (int n = 0; n < col_upp; n++)
		{
			for (int j = 1; j < col_upp; j++)
			{
				if (EtaU_new(i, n) <= Supp_x(i, j))
				{
					double d_eta = EtaU(i, j) - EtaU(i, j - 1);
					double d_s = Supp_x(i, j) - Supp_x(i, j - 1);
					if (d_s == 0)
					{
						d_s = 1;
					}
					double powNum = 1;
					if (ifpow && T!=0)
					{
						powNum = 1 / T;
					}
					EtaU_new(i, n) = EtaU(i, j - 1) + d_eta / pow(d_s, powNum) * pow(EtaU_new(i,n) - Supp_x(i, j - 1),powNum);
					//if (n == 2)
					//{
					//	cout << EtaU_new(ii, n);
					//	cout << "--- " << " " << j << d_eta << " " << d_s << " " << powNum << endl;
					//}
					break;
				}
			}
		}
	}
	////cout << EtaU_new << endl;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//修正 每个网格非独立，整体加密
	for (int j = 0; j < col_upp; j++)
	{
		for (int n = 0; n < row_upp; n++)
		{
			for (int i = 1; i < row_upp; i++)
			{
				if (FaiU_new(n, j) <= Supp_z(i, j))
				{
					double d_fai = FaiU(i, j) - FaiU(i - 1, j);
					double d_s = Supp_z(i, j) - Supp_z(i - 1, j);
					if (d_s == 0)
					{
						d_s = 1;
					}
					double powNum = 1;
					if (i < 2)
					{
						if (ifpow)
						{
							powNum = 1 / N1m(j);
						}
						FaiU_new(n, j) = FaiU(i - 1, j) + d_fai / pow(d_s, powNum) * pow(FaiU_new(n, j) - Supp_z(i - 1, j), powNum);
					}
					else if (i > row_upp - 3)
					{
						if (ifpow)
						{
							powNum = 1 / N2m(j);
						}
						FaiU_new(n, j) = FaiU(i, j) - d_fai / pow(d_s, powNum) * pow(Supp_z(i, j) - FaiU_new(n, j), powNum);
					}
					else
					{
						FaiU_new(n, j) = FaiU(i - 1, j) + d_fai / pow(d_s, 1) * pow(FaiU_new(n, j) - Supp_z(i - 1, j), 1);
					}
					break;
				}
			}
		}
	}

	//FaiU(1, FaiU.n_cols-1) *= 100;
	//FaiU.row(1) *= 0.001;
	//cout << 1-FaiU_new << endl;
	if (!uporlow.compare("upp"))
	{
		this->FaiU = FaiU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		this->EtaU = EtaU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	else
	{
		this->FaiL = FaiU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		this->EtaL = EtaU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	//cout << "$----------------------------------------------$" << endl;
	//cout << FaiU_new << endl;
	//cout << 1 - FaiU_new << endl;
}
//计算正则坐标
void CSTsurface::GetLengthMesh_z(string uporlow, bool ifpow)
{
	//PUPP.printP();
	PointSet PUPP; mat FaiU; mat EtaU; double T; mat N1m; mat N2m;
	//修正上表面
	if (!uporlow.compare("upp"))
	{
		PUPP = this->PUPP;
		FaiU = this->FaiU;
		EtaU = this->EtaU;
		T = this->T[0];
		N1m = this->N1mUPP;
		N2m = this->N2mUPP;
	}
	//修正下表面
	else
	{
		PUPP = this->PLOW;
		FaiU = this->FaiL;
		EtaU = this->EtaL;
		T = this->T[2];
		N1m = this->N1mLOW;
		N2m = this->N2mLOW;
	}
	int row_upp = PUPP.X.n_rows;
	int col_upp = PUPP.X.n_cols;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//mat Supp_x = zeros(row_upp, col_upp);//沿着x方向的距离矩阵
	//mat Supp_x_total = zeros(row_upp, 1);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	mat Supp_z = zeros(row_upp, col_upp);//沿着z方向的距离矩阵
	mat Supp_z_total = zeros(1, col_upp);
	
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	////计算每两个点x方向之间的距离
	//for (int i = 0; i < row_upp; i++)
	//{
	//	for (int j = 1; j < col_upp; j++)
	//	{
	//		Supp_x(i, j) = sqrt(
	//			pow(PUPP.X(i, j-1) - PUPP.X(i, j), 2) +
	//			pow(PUPP.Y(i, j-1) - PUPP.Y(i, j), 2) +
	//			pow(PUPP.Z(i, j-1) - PUPP.Z(i, j), 2));
	//		
	//		Supp_x_total(i, 0) += Supp_x(i, j);
	//	}
	//}
	//
	//for (int j = 1; j < col_upp; j++)//距离求出正则坐标
	//{
	//	Supp_x.col(j) += Supp_x.col(j - 1);
	//}
	//for (int i = 0; i < row_upp; i++)//正则坐标归一化
	//{
	//	if (Supp_x(i, col_upp - 1) != 0)
	//	{
	//		Supp_x.row(i) /= Supp_x(i, col_upp - 1);
	//	}
	//}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	

	mat LinealMeasure = MakeS(1, 3, trans(MeshRefineRatio), FaiU, EtaU);
	//cout << LinealMeasure << endl;


	//计算每两个点z方向之间的距离
	for (int i = 1; i < row_upp; i++)//沿着纵向
	{
		for (int j = 0; j < col_upp; j++)//沿着轴向
		{
			Supp_z(i, j) = sqrt(
				pow(PUPP.X(i - 1, j) - PUPP.X(i, j), 2) +
				pow(PUPP.Y(i - 1, j) - PUPP.Y(i, j), 2) +
				pow(PUPP.Z(i - 1, j) - PUPP.Z(i, j), 2));
			Supp_z(i, j) *= LinealMeasure(i, j);//曲面扭曲
			Supp_z_total(0, j) += Supp_z(i, j);
		}
	}
	for (int i = 1; i < row_upp; i++)//距离求出正则坐标
	{
		Supp_z.row(i) += Supp_z.row(i - 1);
	}
	for (int j = 0; j < col_upp; j++)//正则坐标归一化
	{
		if (Supp_z(row_upp - 1, j) != 0)
		{
			Supp_z.col(j) /= Supp_z(row_upp - 1, j);
		}
	}
	
	//初始化EtaU FaiU
	mat FaiU_new = zeros(row_upp, col_upp);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//mat EtaU_new = zeros(row_upp, col_upp);
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	for (int i = 0; i < row_upp; i++)
	{
		for (int j = 0; j < col_upp; j++)
		{
			FaiU_new(i, j) = (double)i / ((double)row_upp - 1);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//EtaU_new(i, j) = (double)j / ((double)col_upp - 1);
			//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}
	}

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	////修正EtaU FaiU 每个网格单独计算
	//for (int i = 0; i < row_upp; i++)
	//{
	//	for (int n = 0; n < col_upp; n++)
	//	{
	//		for (int j = 1; j < col_upp; j++)
	//		{
	//			if (EtaU_new(i, n) <= Supp_x(i, j))
	//			{
	//				double d_eta = EtaU(i, j) - EtaU(i, j - 1);
	//				double d_s = Supp_x(i, j) - Supp_x(i, j - 1);
	//				if (d_s == 0)
	//				{
	//					d_s = 1;
	//				}
	//				double powNum = 1;
	//				if (ifpow && T!=0)
	//				{
	//					powNum = 1 / T;
	//				}
	//				EtaU_new(i, n) = EtaU(i, j - 1) + d_eta / pow(d_s, powNum) * pow(EtaU_new(i,n) - Supp_x(i, j - 1),powNum);
	//				//if (n == 2)
	//				//{
	//				//	cout << EtaU_new(ii, n);
	//				//	cout << "--- " << " " << j << d_eta << " " << d_s << " " << powNum << endl;
	//				//}
	//				break;
	//			}
	//		}
	//	}
	//}
	////cout << EtaU_new << endl;
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//修正 每个网格非独立，整体加密
	for (int j = 0; j < col_upp; j++)
	{
		for (int n = 0; n < row_upp; n++)
		{
			for (int i = 1; i < row_upp; i++)
			{
				if (FaiU_new(n, j) <= Supp_z(i, j))
				{
					double d_fai = FaiU(i, j) - FaiU(i-1, j);
					double d_s = Supp_z(i, j) - Supp_z(i-1, j);
					if (d_s == 0)
					{
						d_s = 1;
					}
					double powNum = 1;
					if (i < 2)
					{
						if (ifpow)
						{
							powNum = 1 / N1m(j);
						}
						FaiU_new(n, j) = FaiU(i - 1, j) + d_fai / pow(d_s, powNum) * pow(FaiU_new(n, j) - Supp_z(i - 1, j), powNum);
					}
					else if (i > row_upp - 3)
					{
						if (ifpow)
						{
							powNum = 1 / N2m(j);
						}
						FaiU_new(n, j) = FaiU(i, j) - d_fai / pow(d_s, powNum) * pow(Supp_z(i, j) - FaiU_new(n, j), powNum);
					}
					else
					{
						FaiU_new(n, j) = FaiU(i - 1, j) + d_fai / pow(d_s, 1) * pow(FaiU_new(n, j) - Supp_z(i - 1, j), 1);
					}
					break;
				}
			}
		}
	}

	//FaiU(1, FaiU.n_cols-1) *= 100;
	//FaiU.row(1) *= 0.001;
	//cout << 1-FaiU_new << endl;
	if (!uporlow.compare("upp"))
	{
		this->FaiU = FaiU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//this->EtaU = EtaU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	else
	{
		this->FaiL = FaiU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//this->EtaL = EtaU_new;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	//cout << "$----------------------------------------------$" << endl;
	//cout << FaiU_new << endl;
	//cout << 1 - FaiU_new << endl;
}

void CSTsurface::GetLengthMash_x(bool ifpow)
{
	int row_upp = PUPP.X.n_rows;
	int col_upp = PUPP.X.n_cols;
	mat Supp_x = zeros(row_upp, col_upp);//沿着x方向的距离矩阵
	mat Supp_x_total = zeros(row_upp, 1);


	mat LinealMeasureUpp = MakeS(1, 3, trans(MeshRefineRatio), FaiU, EtaU);
	mat LinealMeasureLow = MakeS(1, 3, trans(MeshRefineRatio), FaiL, EtaL);

	//计算每两个点x方向之间的距离
	for (int i = 0; i < row_upp; i++)
	{
		for (int j = 1; j < col_upp; j++)
		{
			Supp_x(i, j) = sqrt(
				pow(PUPP.X(i, j - 1) - PUPP.X(i, j), 2) +
				pow(PUPP.Y(i, j - 1) - PUPP.Y(i, j), 2) +
				pow(PUPP.Z(i, j - 1) - PUPP.Z(i, j), 2));
			Supp_x(i, j) *= LinealMeasureUpp(i, j);//曲面扭曲
			Supp_x_total(i, 0) += Supp_x(i, j);
		}
	}
	//mat lengthUpp_x = Supp_x;
	for (int j = 1; j < col_upp; j++)//距离求出正则坐标
	{
		Supp_x.col(j) += Supp_x.col(j - 1);
	}
	for (int i = 0; i < row_upp; i++)//正则坐标归一化
	{
		if (Supp_x(i, col_upp - 1) != 0)
		{
			Supp_x.row(i) /= Supp_x(i, col_upp - 1);
		}
	}

	int row_low = PLOW.X.n_rows;
	int col_low = PLOW.X.n_cols;
	mat Slow_x = zeros(row_low, col_low);//沿着x方向的距离矩阵
	mat Slow_x_total = zeros(row_low, 1);
	//计算每两个点x方向之间的距离
	for (int i = 0; i < row_low; i++)
	{
		for (int j = 1; j < col_low; j++)
		{
			Slow_x(i, j) = sqrt(
				pow(PLOW.X(i, j - 1) - PLOW.X(i, j), 2) +
				pow(PLOW.Y(i, j - 1) - PLOW.Y(i, j), 2) +
				pow(PLOW.Z(i, j - 1) - PLOW.Z(i, j), 2));
			Slow_x(i, j) *= LinealMeasureLow(i, j);//曲面扭曲
			Slow_x_total(i, 0) += Slow_x(i, j);
		}
	}
	mat lengthLow_x = Slow_x;
	for (int j = 1; j < col_low; j++)//距离求出正则坐标
	{
		Slow_x.col(j) += Slow_x.col(j - 1);
	}
	for (int i = 0; i < row_low; i++)//正则坐标归一化
	{
		if (Slow_x(i, col_low - 1) != 0)
		{
			Slow_x.row(i) /= Slow_x(i, col_low - 1);
		}
	}
	//-----------------------------
	//记录最大值   
	int maxID = 0; double max = 0; bool ifUpp = true;
	for (int i = 0; i < row_upp; i++)
	{
		if (max < Supp_x_total(i, 0))
		{
			max = Supp_x_total(i, 0);
			maxID = i;
			ifUpp = true;
		}
	}
	for (int i = 0; i < row_low; i++)
	{
		if (max < Slow_x_total(i,0))
		{
			max = Slow_x_total(i, 0);
			maxID = i;
			ifUpp = false;
		}
	}
	//初始化EtaU FaiU
	mat EtaL_new = zeros(row_low, col_low);
	mat EtaU_new = zeros(row_upp, col_upp);
	for (int i = 0; i < row_upp; i++)
	{
		for (int j = 0; j < col_upp; j++)
		{
			EtaU_new(i, j) = (double)j / ((double)col_upp - 1);
		}
	}
	for (int i = 0; i < row_low; i++)
	{
		for (int j = 0; j < col_low; j++)
		{
			EtaL_new(i, j) = (double)j / ((double)col_low - 1);
		}
	}
	//修正EtaU FaiU
	int n_row, n_col;
	double T;
	mat S_x, Eta, Eta_new = EtaL_new.row(0);
	if (ifUpp)
	{
		n_row = row_upp;
		n_col = col_upp;
		T = this->T[0];
		S_x = Supp_x;
		Eta = EtaU;
	}
	else
	{
		n_row = row_low;
		n_col = col_low;
		T = this->T[2];
		S_x = Slow_x;
		Eta = EtaL;
	}

	for (int n = 0; n < n_col; n++)
	{
		for (int j = 1; j < n_col; j++)
		{
			int i = 1;
			if (Eta_new(0, n) <= S_x(i, j))
			{
				double d_eta = Eta(i, j) - Eta(i, j - 1);
				double d_s = S_x(i, j) - S_x(i, j - 1);
				if (d_s == 0)
				{
					d_s = 1;
				}
				double powNum = 1;
				if (ifpow && T != 0)
				{
					powNum = 1 / T;
				}
				Eta_new(0, n) = Eta(i, j - 1) + d_eta / pow(d_s, powNum) * pow(Eta_new(0, n) - S_x(i, j - 1), powNum);
				break;
			}
		}
	}
	
	for (int i = 0; i < row_upp; i++)
	{
		EtaU_new.row(i) = Eta_new;
	}
	for (int i = 0; i < row_low; i++)
	{
		EtaL_new.row(i) = Eta_new;
	}
	
	this->EtaU = EtaU_new;
	this->EtaL = EtaL_new;
	
}

//网格数修正
void CSTsurface::MeshNumCorrection()
{
	if (NFaiU % 2 == 0)
	{
		NFaiU = NFaiU + 1;
	}
	if (NFaiL % 2 == 0)
	{
		NFaiL = NFaiL + 1;
	}
	if (NEta % 2 == 0)
	{
		NEta = NEta + 1;
	}
	FaiU = linspace(0, 1, NFaiU);
	FaiL = linspace(0, 1, NFaiL);
	Eta = linspace(0, 1, NEta);
}

//生成形状函数S
mat CSTsurface::MakeS(int PriFuncType, int K, mat B, mat fai, mat eta)
{
	int fainums = fai.n_rows;//向量fai的参数个数
	int etanums = eta.n_cols;//向量eta的参数个数
	mat S = zeros(fainums, etanums);
	
	int NDIM = B.n_rows - 1;//B 权重因子
	int MDIM = B.n_cols - 1;
	
	K = min(NDIM, K);
	K = min(MDIM, K);
	
	mat NBF = zeros(NDIM + K + 1, K + 1);
	mat MBF = zeros(MDIM + K + 1, K + 1);
	//cout <<"K = "<< K << endl;
	for (int i = 0; i < fainums; i++)
	{
		//cout << "============== i = " << i << "==============" << endl;
		for (int j = 0; j < etanums; j++)
		{
			//cout <<" ............. j = " << j << ".............." << endl;
			//方式1：使用bernstrin多项式为基函数
			if (PriFuncType == 0)
			{
				double Kn;
				double Km;
				for (int n = 0; n < NDIM + 1; n++)
				{
					Kn = (double)factorial(NDIM) / ((double)factorial(n) * factorial(NDIM - n));
					for (int m = 0; m < MDIM + 1; m++)
					{
						Km = (double)factorial(MDIM) / ((double)factorial(m) * factorial(MDIM - m));
						S(i, j) = S(i, j) + B(n, m) * Kn * pow(fai(i,j), n) * pow((1 - fai(i,j)), (NDIM - n)) * Km * pow(eta(i,j), m) * pow((1 - eta(i,j)), (MDIM - m));
					}
				}
			}
			//方式2：使用B样条基函数
			if (PriFuncType != 0)
			{
				for (int k = 0; k < K + 1; k++)
				{
					//cout << "k = " << k << endl;
					//生成准均匀B样条的节点参数tn
					mat tnleft = zeros(1, K);
					mat tnright = ones(1, K);
					int tnNum = NDIM + K + 2 - 2 * K;
					mat tn = zeros(1, tnNum);
					for (int i = 0; i < tnNum; i++)
					{
						tn(0, i) = i / ((double)tnNum - 1);
					}
					tn = join_rows(tnleft, tn);
					tn = join_rows(tn, tnright);
					if (k == 0)
					{
						for (int n = 0; n < NDIM + K + 1; n++)
						{
							//int kID = k + 1;
							//int nID = n + 1;
							if ((tn(n) <= fai(i,j)) && (fai(i,j) < tn(n + 1)))//=====================
								NBF(n, k) = 1;
							else if( (tn(n) < fai(i,j)) && (fai(i,j) <= tn(n + 1)))//(i == fainums - 1) &&********7.21修改
								NBF(n, k) = 1;
							else
								NBF(n, k) = 0;
						}
						//cout << NBF << endl;
					}
					else
					{
						for (int n = 0; n < NDIM + K - k + 1; n++)
						{
							double alpha = (fai(i,j) - tn(n)) / (tn(n + k) - tn(n));
							double beta = (tn(n + k + 1) - fai(i,j)) / (tn(n + k + 1) - tn(n + 1));
							if (isnan(alpha) || isinf(alpha)) alpha = 0;
							if (isnan(beta) || isinf(beta))   beta = 0;
							NBF(n, k) = alpha * NBF(n, k - 1) + beta * NBF(n + 1, k - 1);
						}
					}
					//生成准均匀B样条的节点参数tm
					mat tmleft = zeros(1, K);
					mat tmright = ones(1, K);
					int tmNum = MDIM + K + 2 - 2 * K;
					mat tm = zeros(1, tmNum);
					for (int i = 0; i < tmNum; i++)
						tm(0, i) = i / ((double)tmNum - 1);
					tm = join_rows(tmleft, tm);
					tm = join_rows(tm, tmright);
					//cout << "tm = " << tm << endl;
					if (k == 0)
					{
						for (int n = 0; n < MDIM + K + 1; n++)
						{
							//int kID = k + 1;
							//int nID = n + 1;
							if (tm(n) <= eta(i,j) && eta(i,j) < tm(n + 1))//==================
								MBF(n, k) = 1;
							else
								MBF(n, k) = 0;
							if ((j == etanums - 1) && (tm(n) < eta(i,j)) && (eta(i,j) <= tm(n + 1)))
								MBF(n, k) = 1;
						}
					}
					else
					{
						for (int m = 0; m < MDIM + K - k + 1; m++)
						{
							double alpha = (eta(i,j) - tm(m)) / (tm(m + k) - tm(m));
							double beta = (tm(m + k + 1) - eta(i,j)) / (tm(m + k + 1) - tm(m + 1));
							if (isnan(alpha) || isinf(alpha)) alpha = 0;
							if (isnan(beta) || isinf(beta))   beta = 0;
							MBF(m, k) = alpha * MBF(m, k - 1) + beta * MBF(m + 1, k - 1);
							//MBF(m, k) = 1 * MBF(m, k - 1) + 2 * MBF(m + 1, k - 1);
						}
					}
				}
				//cout << NBF << endl;
				for (int n = 0; n < NDIM + 1; n++)
					for (int m = 0; m < MDIM + 1; m++)
						S(i, j) = S(i, j) + B(n, m) * NBF(n, K) * MBF(m, K);
		

			}//cout << MBF << endl;
		}
	}
	return S;
}

//CST基本数学模型
PointSet CSTsurface::SketchCST(mat fai, mat eta, double N1, double N2, double N3, double N4,
	double M1, double M2, double T1, double T2, mat B, mat D, double Ratio, mat&N1m,mat&N2m)//Ratio(int2double)
{
	//初始配置
	int fainums = fai.n_rows;//向量fai的参数个数
	int etanums = eta.n_cols;//向量eta的参数个数
	int K = 3;//形状函数S的次数（阶数-1）。默认为3
	PointSet p;

	//纵向修正
	int NDIM = B.n_rows - 1;//权重因子B 行数-1
	int MDIM = B.n_cols - 1;//权重因子B 列数-1

	if (NDIM == 0 && MDIM != 0)
	{
		B = repmat(B, K + 1, 1);
	}
	else if (MDIM == 0 && NDIM != 0)
	{
		B = repmat(B, 1, K + 1);
	}

	mat S = MakeS(PriFuncType, K, B, fai, eta);
	//cout << "SketchCST中的S = "<<endl << S << endl;//测试使用

	//侧向修正
	NDIM = D.n_rows - 1;//侧向权重因子D 行数-1
	MDIM = D.n_cols - 1;//侧向权重因子D 列数-1

	if (NDIM == 0 && MDIM != 0)
	{
		D = repmat(D, K + 1, 1);
	}
	else if (MDIM == 0 && NDIM != 0)
	{
		D = repmat(D, 1, K + 1);
	}
	mat St = MakeS(PriFuncType, K, D, fai, eta);

//	cout << fai << endl;
	//fai = linspace(0, 1, fai.size());
	// cout << fai << endl;
//	eta = trans(linspace(0, 1, eta.size()));

	//区域限制
	double N_max = N1 / (N1 + N2);
	double M_max = M1 / (M1 + M2);
	double T_max = T1 / (T1 + T2);

	double eta_max = pow(M_max, M1) * pow(1 - M_max, M2);
	if (isnan(eta_max))
	{
		eta_max = 1;
	}

	double tao_max = pow(T_max, T1) * pow(1 - T_max, T2);
	if (isnan(tao_max))
	{
		tao_max = 1;
	}

	//修正分析
	p.X = eta;// repmat(eta, fainums, 1);
	p.Z = zeros(fainums, etanums);

	mat Zmax = zeros(fainums, etanums);
	for (int i = 0; i < fainums; i++)
	{
		for (int n = 0; n < etanums; n++)
		{
			//cout << 1 - eta(i, n) << endl;
			Zmax(i, n) = pow(eta(i, n), T1) * pow(1 - eta(i, n), T2) / tao_max;//归一化
		}
	}
	//mat Zmax = pow(eta, T1) * pow(1 - eta, T2) / tao_max;
	mat Zmin = -Zmax;
	for (int i = 0; i < fainums; i++)
	{
		for (int j = 0; j < etanums; j++)
		{
			p.Z(i, j) = Zmin(i,j) + (Zmax(i,j) - Zmin(i,j)) * fai(i,j);
		}
		//p.Z.row(i) = Zmin + (Zmax - Zmin) * fai(i);
	}
	for (int i = 0; i < fainums; i++)
	{
		for (int j = 0; j < etanums; j++)
		{
			p.Z(i, j) = p.Z(i, j) * St(i, j);

		}
	}
	//cout << p.Z << endl;
	//p.Z = p.Z * St;

	mat ETA = zeros(fainums, etanums);
	for (int i = 0; i < fainums; i++)
	{
		for (int n = 0; n < etanums; n++)
		{
			ETA(i, n) = pow(eta(i, n), M1) * pow(1 - eta(i, n), M2) / eta_max;
		}
	}
	//mat ETA = pow(eta, M1) * pow(1 - eta, M2) / eta_max;
	mat FAI = zeros(fainums, etanums);


	for (int j = 0; j < etanums; j++)
	{
		N1m(j) = N1 + pow((double)j / ((double)etanums - 1), Ratio) * (N3 - N1);
		N2m(j) = N2 + pow((double)j / ((double)etanums - 1), Ratio) * (N4 - N2);
		double N_max = N1m(j) / (N1m(j) + N2m(j));
		double fai_max = pow(N_max, N1m(j)) * pow(1 - N_max, N2m(j));
		
		if (isnan(fai_max))
		{
			fai_max = 1;
		}
		for (int n = 0; n < fainums; n++)
		{
			FAI(n, j) = pow(fai(n,j), N1m(j)) * pow(1 - fai(n,j), N2m(j)) * ETA(n, j) / fai_max;
		}
	}
	//cout << N1m << endl << N2m << endl;
	p.Y = zeros(fainums, etanums);
	for (int i = 0; i < fainums; i++)
	{
		for (int j = 0; j < etanums; j++)
		{
			p.Y(i, j) = FAI(i, j) * S(i, j);
		}
	}
	return p;
}

//曲率评估 计算得出曲面曲率值->Info
void CSTsurface::EvalCurvature()
{
	
	//参数传入
	double OriginX = Origin(0);
	double OriginY = Origin(1);
	double OriginZ = Origin(2);
	double RotationX = Rotation(0) / 180.0 * pi;
	double RotationY = Rotation(1) / 180.0 * pi;
	double RotationZ = Rotation(2) / 180.0 * pi;
	double Xlength = Length(0);
	double Ylength_upp = Length(1);
	double Ylength_low = Length(2);
	double Zlength = Length(3);
	double SlopeX = Slope(0) / 180 * pi;
	double SlopeY = Slope(1) / 180 * pi;
	double SlopeZ = Slope(2) / 180 * pi;
	double ScaleX = Scale(0);
	double ScaleYUPP = Scale(1);
	double ScaleYLOW = Scale(2);
	double ScaleZ = Scale(3);
	//类别函数控制参数
	double N1s = NS(0);	double N2s = NS(1);	double N3s = NS(2);	double N4s = NS(3);	//截面
	double N1e = NE(0);	double N2e = NE(1);	double N3e = NE(2);	double N4e = NE(3);	//截面
	double M1 = M(0);	double M2 = M(1);	double M3 = M(2);	double M4 = M(3);	//轴线
	double T1 = T(0);	double T2 = T(1);	double T3 = T(2);	double T4 = T(3);	//侧向

	//网格数修正
	MeshNumCorrection();
	//生成网格列表：FaiU、FaiL、Eta
	GetNet();
	Info.FaiU = FaiU;
	Info.FaiL = FaiL;
	Info.EtaU = EtaU;
	Info.EtaL = EtaL;

	N1mUPP = zeros(1, EtaU.n_cols);
	N2mUPP = zeros(1, EtaU.n_cols);
	N1mLOW = zeros(1, EtaU.n_cols);
	N2mLOW = zeros(1, EtaU.n_cols);
	//初始分布计算
	PointSet PUPPcs = SketchCST(FaiU, EtaU, N1s, N2s, N1e, N2e, M1, M2, T1, T2, BUPP, DUPP, Ratio, N1mUPP, N2mUPP);
	PointSet PLOWcs = SketchCST(FaiL, EtaL, N3s, N3s, N3e, N4e, M3, M4, T3, T4, BLOW, DLOW, Ratio, N1mLOW, N2mLOW);

	//距离修正
	mat XYDispUpp = zeros(NFaiU, NEta - 1);
	mat XYDispLow = zeros(NFaiL, NEta - 1);
	mat YZDispUpp = zeros(NFaiU - 1, NEta);
	mat YZDispLow = zeros(NFaiL - 1, NEta);

	for (int i = 0; i < NFaiU; i++)
	{
		for (int j = 0; j < NEta - 1; j++)
		{
			XYDispUpp(i, j) = sqrt(pow(PUPPcs.X(i, j + 1) - PUPPcs.X(i, j), 2) + pow(PUPPcs.Y(i, j + 1) - PUPPcs.Y(i, j), 2));
		}
	}

	for (int i = 0; i < NFaiL; i++)
	{
		for (int j = 0; j < NEta - 1; j++)
		{
			XYDispLow(i, j) = sqrt(pow(PLOWcs.X(i, j + 1) - PLOWcs.X(i, j), 2) + pow(PLOWcs.Y(i, j + 1) - PLOWcs.Y(i, j), 2));
		}
	}

	for (int j = 0; j < NEta; j++)
	{
		for (int i = 0; i < NFaiU - 1; i++)
		{
			YZDispUpp(i, j) = sqrt(pow(PUPPcs.Z(i + 1, j) - PUPPcs.Z(i, j), 2) + pow(PUPPcs.Y(i + 1, j) - PUPPcs.Y(i, j), 2));
		}
		for (int i = 0; i < NFaiL - 1; i++)
		{
			YZDispLow(i, j) = sqrt(pow(PLOWcs.Z(i + 1, j) - PLOWcs.Z(i, j), 2) + pow(PLOWcs.Y(i + 1, j) - PLOWcs.Y(i, j), 2));
		}

	}

	mat CURyx = join_cols(XYDispUpp, XYDispLow);
	CURyx = max(CURyx);//寻找每列的最大值，存成一行
	mat CURyzUpp = max(YZDispUpp, 1);//寻找行的最大值，存成一列
	mat CURyzLow = max(YZDispLow, 1);//寻找行的最大值，存成一列

	//曲率计算
	for (size_t i = 0; i < CURyx.n_elem; i++)
	{
		if (CURyx(i) == 0)
		{
			CURyx(i) = 1;
		}
	}
	for (size_t i = 0; i < CURyzUpp.n_elem; i++)
	{
		if (CURyzUpp(i) == 0)
		{
			CURyzUpp(i) = 1;
		}
	}
	for (size_t i = 0; i < CURyzLow.n_elem; i++)
	{
		if (CURyzLow(i) == 0)
		{
			CURyzLow(i) = 1;
		}
	}

	mat MaxMin;//计算下列向量的最大最小值，MaxMin为中间变量
	MaxMin = max(CURyx, 1);
	double CURyxMax = MaxMin(0, 0);
	MaxMin = min(CURyx, 1);
	double CURyxMin = MaxMin(0, 0);
	MaxMin = max(CURyzUpp);
	double CURyzUppMax = MaxMin(0, 0);
	MaxMin = min(CURyzUpp);
	double CURyzUppMin = MaxMin(0, 0);
	MaxMin = max(CURyzLow);
	double CURyzLowMax = MaxMin(0, 0);
	MaxMin = min(CURyzLow);
	double CURyzLowMin = MaxMin(0, 0);

	double DallowMax = min(CURyxMax, 10.0 / (double)NEta);
	double DallowMin = max(CURyxMin, 0.1 / (double)NEta);
	Info.CURyx = zeros(1, CURyx.n_elem);
	for (size_t i = 0; i < CURyx.n_elem; i++)
	{
		Info.CURyx(i) = (CURyx(i) - CURyxMin) / (CURyxMax - CURyxMin) * (DallowMax - DallowMin) + DallowMin;
	}

	DallowMax = min(CURyzUppMax, 10.0 / (double)NFaiU);
	DallowMin = max(CURyzUppMin, 0.1 / (double)NFaiU);
	Info.CURyzUpp = zeros(CURyzUpp.n_elem, 1);
	for (size_t i = 0; i < CURyzUpp.n_elem; i++)
	{
		Info.CURyzUpp(i) = (CURyzUpp(i) - CURyzUppMin) / (CURyzUppMax - CURyzUppMin) * (DallowMax - DallowMin) + DallowMin;
	}

	DallowMax = min(CURyzLowMax, 10.0 / (double)NFaiL);
	DallowMin = max(CURyzLowMin, 0.1 / (double)NFaiL);
	Info.CURyzLow = zeros(CURyzLow.n_elem, 1);
	for (size_t i = 0; i < CURyzLow.n_elem; i++)
	{
		Info.CURyzLow(i) = (CURyzLow(i) - CURyzLowMin) / (CURyzLowMax - CURyzLowMin) * (DallowMax - DallowMin) + DallowMin;
	}

}

//网格修正（根据曲率修正网格）
void CSTsurface::GridRefine()
{
	Eta = zeros(1, Info.CURyx.n_elem + 1);
	Eta(0) = 0;
	for (size_t i = 1; i < Info.CURyx.n_elem + 1; i++)
	{
		Eta(i) = 1 / Info.CURyx(i - 1);
	}
	for (size_t i = 1; i < Eta.n_elem; i++)
	{
		mat Temp = max(Eta, 1);
		double temp = Temp(0, 0);
		Eta(i) = Eta(i) / temp;
	}
	//cout << Eta << endl;
	FaiU = zeros(1, Info.CURyzUpp.n_elem + 1);
	FaiU(0) = 0;
	for (size_t i = 1; i < Info.CURyzUpp.n_elem + 1; i++)
	{
		FaiU(i) = 1 / Info.CURyzUpp(i - 1);
	}
	//cout << FaiU << endl;
	for (size_t i = 1; i < FaiU.n_elem; i++)
	{
		mat Temp = max(FaiU, 1);
		double temp = Temp(0, 0);
		FaiU(i) = FaiU(i) / temp;
	}
	//cout << this->FaiU << endl;
	FaiL = zeros(1, Info.CURyzLow.n_elem + 1);
	FaiL(0) = 0;
	for (size_t i = 1; i < Info.CURyzLow.n_elem + 1; i++)
	{
		FaiL(i) = 1 / Info.CURyzLow(i - 1);
	}
	for (size_t i = 1; i < FaiL.n_elem; i++)
	{
		mat Temp = max(FaiL, 1);
		double temp = Temp(0, 0);
		FaiL(i) = FaiL(i) / temp;
	}
}

//评估网格信息 计算得出曲面面积等信息->Info
void CSTsurface::EvalGrid()
{
	double VolumeBody = 0;//体积
	double SurfaceAeraBody = 0;//zy平面截面面积
	double planformarea = 0;//xz平面面积

	for (int ix = 0; ix < NEta - 1; ix++)
	{
		double area = 0;
		for (int iu = 0; iu < NFaiU - 1; iu++)
		{
			area += (-1) * PUPP.Z(iu, ix) * PUPP.Y(iu + 1, ix) - (-1) * PUPP.Z(iu + 1, ix) * PUPP.Y(iu, ix);
			SurfaceAeraBody = SurfaceAeraBody + (sqrt(pow(PUPP.Z(iu, ix) - PUPP.Z(iu + 1, ix), 1) + pow(PUPP.Y(iu, ix) - PUPP.Y(iu + 1, ix), 1))
				+ sqrt(pow(PUPP.Z(iu, ix + 1) - PUPP.Z(iu + 1, ix + 1), 1) + pow(PUPP.Y(iu, ix + 1) - PUPP.Y(iu + 1, ix + 1), 1))) / 2
				* (PUPP.X(0, ix + 1) - PUPP.X(0, ix));
			planformarea = planformarea + (PUPP.X(0, ix + 1) - PUPP.X(0, ix)) * abs(PUPP.Z(iu + 1, ix) - PUPP.Z(iu, ix));
		}
		for (int il = NFaiL - 1; il > 0; il--)
		{
			area += PLOW.Z(il, ix) * PLOW.Y(il - 1, ix) - PLOW.Z(il - 1, ix) * PLOW.Y(il, ix);
			SurfaceAeraBody = SurfaceAeraBody + 
				(sqrt(pow(PLOW.Z(il, ix) - PLOW.Z(il - 1, ix), 1) + pow(PLOW.Y(il, ix) - PLOW.Y(il - 1, ix), 1))+ 
				sqrt(pow(PLOW.Z(il, ix + 1) - PLOW.Z(il - 1, ix + 1), 1) + pow(PLOW.Y(il, ix + 1) - PLOW.Y(il - 1, ix + 1), 1))) / 2
				* (PLOW.X(0, ix + 1) - PLOW.X(0, ix));
		}
		VolumeBody = VolumeBody + (PUPP.X(0, ix + 1) - PUPP.X(0, ix)) * area / 2;
	}
	Info.SurfaceAera = abs(SurfaceAeraBody);
	Info.Volume = abs(VolumeBody);
	Info.planforarea = abs(planformarea);

	
}

//生成默认网格节点列表：FaiU、FaiL、Eta
void CSTsurface::GetNet(bool ifrefresh)
{
	if (ifrefresh)//强制更新
	{
		FaiU = zeros(NFaiU, NEta);
		EtaU = zeros(NFaiU, NEta);
		for (int i = 0; i < NFaiU; i++)
		{
			for (int j = 0; j < NEta; j++)
			{
				FaiU(i, j) = (double)i / ((double)NFaiU - 1);
				EtaU(i, j) = (double)j / ((double)NEta - 1);
			}
		}
		//cout << FaiU << endl;
	
		FaiL = zeros(NFaiL, NEta);
		EtaL = zeros(NFaiL, NEta);
		for (int i = 0; i < NFaiL; i++)
		{
			for (int j = 0; j < NEta; j++)
			{
				FaiL(i, j) = (double)i / ((double)NFaiL - 1);
				EtaL(i, j) = (double)j / ((double)NEta - 1);
			}
		}
	}
	else//如果原先已有就不更新
	{
		if (FaiU.n_rows != NFaiU || EtaU.n_cols != NEta)
		{
			FaiU = zeros(NFaiU, NEta);
			EtaU = zeros(NFaiU, NEta);
			for (int i = 0; i < NFaiU; i++)
			{
				for (int j = 0; j < NEta; j++)
				{
					FaiU(i, j) = (double)i / ((double)NFaiU - 1);
					EtaU(i, j) = (double)j / ((double)NEta - 1);
				}
			}
		}


		//cout << FaiU << endl;
		if (FaiL.n_rows != NFaiL || EtaL.n_cols != NEta)
		{
			FaiL = zeros(NFaiL, NEta);
			EtaL = zeros(NFaiL, NEta);
			for (int i = 0; i < NFaiL; i++)
			{
				for (int j = 0; j < NEta; j++)
				{
					FaiL(i, j) = (double)i / ((double)NFaiL - 1);
					EtaL(i, j) = (double)j / ((double)NEta - 1);
				}
			}
		}
	}

	
	//cout << FaiU << endl; cout << EtaU << endl;
	//Eta = zeros(1, NEta);
	//for (int i = 0; i < NEta; i++)
	//{
	//	Eta(0, i) = (double)i / ((double)NEta - 1);
	//}

}

//调用GMesh生成前后面网格
//void CSTsurface::MeshPlane() {}



//生成该曲面的梁结构点*****4.12更新
mat CSTsurface::MakeBone(double bias1,double bias2)//(偏置量bias取 - 1to1)
{
	if (bias1 < -1 || bias1 > 1)
	{
		cout << "CSTsurface::GetNet--bias取值错误！！" << endl;
	}
	if (bias2 < -1 || bias2 > 1)
	{
		cout << "CSTsurface::GetNet--bias取值错误！！" << endl;
	}
	int numX = PUPP.X.n_cols;
	int numZupp = PUPP.X.n_rows;
	int numZlow = PLOW.X.n_rows;
	int midZupp = (numZupp ) / 2 ;
	int midZlow = numZlow / 2;
	BonePoint = zeros(numX, 3);//上表面中心线

	double delta1 = 0; double delta2 = 0;

	for (int i = 0; i < numX; i++)//单个切片内简化到一个点上
	{
		double lenAll = 0;//lenAll：切片的周长
		for (int j = 0; j < numZupp - 1; j++)
		{
			lenAll += sqrt(
				pow(PUPP.X(j + 1, i) - PUPP.X(j, i), 2) +
				pow(PUPP.Y(j + 1, i) - PUPP.Y(j, i), 2) +
				pow(PUPP.Z(j + 1, i) - PUPP.Z(j, i), 2));
		}
		for (int j = 0; j < numZlow - 1; j++)
		{
			lenAll += sqrt(
				pow(PLOW.X(j + 1, i) - PLOW.X(j, i), 2) +
				pow(PLOW.Y(j + 1, i) - PLOW.Y(j, i), 2) +
				pow(PLOW.Z(j + 1, i) - PLOW.Z(j, i), 2));
		}

		//计算形心
		double x_star = 0; double y_star = 0; double z_star = 0;
		for (int j = 0; j < numZupp - 1; j++)
		{
			double len = sqrt(
				pow(PUPP.X(j + 1, i) - PUPP.X(j, i), 2) +
				pow(PUPP.Y(j + 1, i) - PUPP.Y(j, i), 2) +
				pow(PUPP.Z(j + 1, i) - PUPP.Z(j, i), 2));
			double x = (PUPP.X(j + 1, i) + PUPP.X(j, i)) / 2;
			double y = (PUPP.Y(j + 1, i) + PUPP.Y(j, i)) / 2;
			double z = (PUPP.Z(j + 1, i) + PUPP.Z(j, i)) / 2;
			x_star += len * x;
			y_star += len * y;
			z_star += len * z;
		}
		for (int j = 0; j < numZlow - 1; j++)
		{
			double len = sqrt(
				pow(PLOW.X(j + 1, i) - PLOW.X(j, i), 2) +
				pow(PLOW.Y(j + 1, i) - PLOW.Y(j, i), 2) +
				pow(PLOW.Z(j + 1, i) - PLOW.Z(j, i), 2));

			double x = (PLOW.X(j + 1, i) + PLOW.X(j, i)) / 2;
			double y = (PLOW.Y(j + 1, i) + PLOW.Y(j, i)) / 2;
			double z = (PLOW.Z(j + 1, i) + PLOW.Z(j, i)) / 2;
			x_star += len * x;
			y_star -= len * y;
			z_star += len * z;
		}

		if (lenAll == 0) lenAll = 1;
		double xc = x_star / lenAll;
		double yc = y_star / lenAll;
		double zc = z_star / lenAll;
		BonePoint(i, 0) = PUPP.X(midZupp, i);
		BonePoint(i, 1) = yc;
		//double bias = bias1 + (bias2 - bias1) * (PUPP.X(0, i) - PUPP.X(0, 0)) / (PUPP.X(0, numX - 1) - PUPP.X(0, 0));
		
		if (i ==0)
		{
			if (bias1 < 0) delta1 = (zc - PLOW.Z(0, i)) * bias1;
			else delta1 = (PLOW.Z((numZlow - 1), i) - zc) * bias1;
		}
		if (i == numX-1)
		{
			if (bias2 < 0) delta2 = (zc - PLOW.Z(0, i)) * bias2;
			else delta2 = (PLOW.Z((numZlow - 1), i) - zc) * bias2;
		}
		
		BonePoint(i, 2) = zc;// +delta1 + (delta2 - delta1) * (PUPP.X(0, i) - PUPP.X(0, 0)) / (PUPP.X(0, numX - 1) - PUPP.X(0, 0));// PUPP.Z(midZupp, i);
	}
	for (int i = 0; i < numX; i++)//单个切片内简化到一个点上
		BonePoint(i,2) += delta1 + (delta2 - delta1) * (PUPP.X(0, i) - PUPP.X(0, 0)) / (PUPP.X(0, numX - 1) - PUPP.X(0, 0));

	//cout << "未平移和旋转的骨架节点坐标" << endl;
	//cout << LineMid << endl;

	double OriginX = Origin(0);
	double OriginY = Origin(1);
	double OriginZ = Origin(2);

	double RotationX = Rotation(0) / 180.0 * pi;
	double RotationY = Rotation(1) / 180.0 * pi;
	double RotationZ = Rotation(2) / 180.0 * pi;

	//还原角度 通过Rotation Matrix
	mat Lrotation = { { cos(RotationZ) * cos(RotationY) , sin(RotationZ) , -cos(RotationZ) * sin(RotationY) },
					 {-sin(RotationZ) * cos(RotationY) * cos(RotationX) + sin(RotationY) * sin(RotationX), cos(RotationZ) * cos(RotationX)	, sin(RotationZ) * sin(RotationY) * cos(RotationX) + cos(RotationY) * sin(RotationX)},
					 { sin(RotationZ) * cos(RotationY) * sin(RotationX) + sin(RotationY) * cos(RotationX),-cos(RotationZ) * sin(RotationX)	,-sin(RotationZ) * sin(RotationY) * sin(RotationX) + cos(RotationY) * cos(RotationX)} };

	BonePoint = BonePoint * Lrotation;
	BonePoint.col(0) = BonePoint.col(0) + OriginX;
	BonePoint.col(1) = BonePoint.col(1) + OriginY;
	BonePoint.col(2) = BonePoint.col(2) + OriginZ;
	//LineMid = flipud(LineMid);
	//cout << "这就是骨架节点坐标" << endl;
	//cout << LineMid << endl;
	return BonePoint;
}

//计算骨架梁节点上的力和力矩
mat CSTsurface::ClacBonePForce()
{
	
	mat BonePForce = zeros(BonePoint.n_rows,6);
	if (GridUpp.nodeForce.n_rows == 0)
	{
		return mat();
	}
	int numX = PUPP.X.n_cols;
	int numZupp = PUPP.X.n_rows;
	int numZlow = PLOW.X.n_rows;
	for (int i = 0; i < numX; i++)//每一个环
	{
		for (int j = 0; j < numZupp; j++)//每一个环内的点
		{
			int id = j * numX + i;
			mat force = { GridUpp.nodeForce(id, 0), GridUpp.nodeForce(id, 1), GridUpp.nodeForce(id, 2) };//该点的力
			mat siteF = { GridUpp.P(id, 0), GridUpp.P(id, 1), GridUpp.P(id, 2) };//该点的坐标
			mat siteP = { BonePoint(i, 0), BonePoint(i, 1), BonePoint(i, 2) };//结构梁的坐标
			mat moment = cross(siteF - siteP, force);

			//cout << "site F :" << siteF << endl;
			//cout << "siteF - siteP :" << siteF - siteP ;
			//cout << "   moment     :" << moment << endl;
			
			BonePForce(i, 0) += force(0, 0);//fx
			BonePForce(i, 1) += force(0, 1);//fy
			BonePForce(i, 2) += force(0, 2);//fz
			BonePForce(i, 3) += moment(0, 0);//mx
			BonePForce(i, 4) += moment(0, 1);//my
			BonePForce(i, 5) += moment(0, 2);//mz
		}
		for (int j = 0; j < numZlow; j++)//每一个环内的点
		{
			int id = j * numX + i;
			mat force = { GridLow.nodeForce(id, 0), GridLow.nodeForce(id, 1), GridLow.nodeForce(id, 2) };//该点的力
			mat siteF = { GridLow.P(id, 0), GridLow.P(id, 1), GridLow.P(id, 2) };//该点的坐标
			mat siteP = { BonePoint(i, 0), BonePoint(i, 1), BonePoint(i, 2) };//结构梁的坐标
			mat moment = cross(siteF - siteP, force);

			BonePForce(i, 0) += force(0, 0);//fx
			BonePForce(i, 1) += force(0, 1);//fy
			BonePForce(i, 2) += force(0, 2);//fz
			BonePForce(i, 3) += moment(0, 0);//mx
			BonePForce(i, 4) += moment(0, 1);//my
			BonePForce(i, 5) += moment(0, 2);//mz
		}
	} 
	return BonePForce;
}

mat CSTsurface::ClacBoneS()
{
	mat BoneS = zeros(BonePoint.n_rows, 6);


	return BoneS;
}

int CSTsurface::MakeStruct2D(mat X /*= mat(0, 0)*/, mat Z /*= mat(0, 0)*/)
{
	//----------1、数据处理/检查----------
	//让高度方向网格至少有2
	if (NHeight < 2)
	{
		NHeight = 2;//（上下边界表皮2个）
	}
	int numi = PUPP.X.n_rows;
	int numj = PUPP.X.n_cols;
	int numk = NHeight;
	
	if (PLOW.X.n_elem != PUPP.X.n_elem )
	{
		cout << "MakcStruct错误！！上下网格数不一致！" << endl;
		return -1;
	}

	if (X.n_rows > 1 && X.n_cols > 1)
	{
		cout << "MakcStruct错误！！ X必须是一维向量！" << endl;
		return -1;
	}
	if (Z.n_rows > 1 && Z.n_cols > 1)
	{
		cout << "MakcStruct错误！！ Z必须是一维向量！" << endl;
		return -1;
	}
	
	//----------2、转换为正则坐标&使用参数 X Z改变结构梁位置----------5.30 add 
	
	mat PLOWX, PUPPX, PLOWY, PUPPY, PLOWZ, PUPPZ;
	if (X.n_elem != 0 && Z.n_elem != 0)
	{
		PLOWX = PUPPX = PLOWY = PUPPY = PLOWZ = PUPPZ = zeros(Z.n_elem, X.n_elem);
		mat PLOWXm, PUPPXm, PLOWYm, PUPPYm, PLOWZm, PUPPZm;
		double lenX = 0;//lenX：第一个周长非零的切片周长(x方向)
		mat siteX = zeros(numj);
		for (int i = 0; i < numi; i++)
		{
			for (int j = 0; j < numj - 1; j++)
			{
				lenX += sqrt(
					pow(PUPP.X(i, j + 1) - PUPP.X(i, j), 2) +
					pow(PUPP.Y(i, j + 1) - PUPP.Y(i, j), 2) +
					pow(PUPP.Z(i, j + 1) - PUPP.Z(i, j), 2));
				siteX(j + 1) = lenX;
			}
			if (lenX > 0)	break;
		}
		siteX = siteX / lenX;

		mat siteZ = linspace(0, 1, numi);
		interp2(siteX, siteZ, PLOW.X, X, siteZ, PLOWXm);
		interp2(siteX, siteZ, PUPP.X, X, siteZ, PUPPXm);
		interp2(siteX, siteZ, PLOW.Y, X, siteZ, PLOWYm);
		interp2(siteX, siteZ, PUPP.Y, X, siteZ, PUPPYm);
		interp2(siteX, siteZ, PLOW.Z, X, siteZ, PLOWZm);
		interp2(siteX, siteZ, PUPP.Z, X, siteZ, PUPPZm);
		numi = PUPPXm.n_rows; numj = PUPPXm.n_cols;
		
		siteZ = zeros(numi);
		for (int j = 0; j < numj; j++)
		{
			double lenZ = 0;//lenZ：第一个周长非零的切片周长(z方向)
			for (int i = 0; i < numi - 1; i++)
			{
				lenZ += sqrt(
					pow(PUPPXm(i + 1, j) - PUPPXm(i, j), 2) +
					pow(PUPPYm(i + 1, j) - PUPPYm(i, j), 2) +
					pow(PUPPZm(i + 1, j) - PUPPZm(i, j), 2));
				siteZ(i + 1) = lenZ;
			}
			if (lenZ > 0)
			{
				siteZ = siteZ / lenZ;

				mat temp;
				interp1(siteZ, PLOWXm.col(j), Z, temp); PLOWX.col(j) = temp;
				interp1(siteZ, PUPPXm.col(j), Z, temp); PUPPX.col(j) = temp;
				interp1(siteZ, PLOWYm.col(j), Z, temp); PLOWY.col(j) = temp;
				interp1(siteZ, PUPPYm.col(j), Z, temp); PUPPY.col(j) = temp;
				interp1(siteZ, PLOWZm.col(j), Z, temp); PLOWZ.col(j) = temp;
				interp1(siteZ, PUPPZm.col(j), Z, temp); PUPPZ.col(j) = temp;
			}
			else
			{
				PLOWX.col(j) = ones(PUPPX.n_rows, 1) * PLOWXm(0, j);
				PUPPX.col(j) = ones(PUPPX.n_rows, 1) * PUPPXm(0, j);
				PLOWY.col(j) = ones(PUPPX.n_rows, 1) * PLOWYm(0, j);
				PUPPY.col(j) = ones(PUPPX.n_rows, 1) * PUPPYm(0, j);
				PLOWZ.col(j) = ones(PUPPX.n_rows, 1) * PLOWZm(0, j);
				PUPPZ.col(j) = ones(PUPPX.n_rows, 1) * PUPPZm(0, j);
			}
		}
		numi = PUPPX.n_rows; numj = PUPPX.n_cols;
	}
	else
	{
		PLOWX = PLOW.X; PUPPX = PUPP.X;
		PLOWY = PLOW.Y; PUPPY = PUPP.Y;
		PLOWZ = PLOW.Z; PUPPZ = PUPP.Z;
	}
	//----------------------------------------


	//----------3、点和单元的生成----------
	//生成结构所需的点集合，每个mat是一个高度的点集
	vector<PointSet> points(NHeight);
	for (int i = 0; i < NHeight; i++)
	{
		points[i].X = PLOWX + (PUPPX - PLOWX) * i / (NHeight - 1.0);
		points[i].Y = -PLOWY + (PUPPY + PLOWY) * i / (NHeight - 1.0);//y方向PLOW缺少负号，故取负
		points[i].Z = PLOWZ + (PUPPZ - PLOWZ) * i / (NHeight - 1.0);
	}
	
	int elemNum = (numi - 1) * (numk - 1) * numj + (numj - 1) * (numk - 1) * numi;
	mat Elem = zeros(elemNum, 4);
	//生成结构单元 两个循环一纵一横
	int elemID = 0;
	for (int j = 0; j < numj; j++)
	{
		for (int i = 0; i < numi - 1; i++)
		{
			for (int k = 0; k < numk-1; k++)
			{
				Elem(elemID, 0) = StructPart::ijk2ID(i, j, k, numi, numj, numk);
				Elem(elemID, 1) = StructPart::ijk2ID(i + 1, j, k, numi, numj, numk);
				Elem(elemID, 2) = StructPart::ijk2ID(i + 1, j, k + 1, numi, numj, numk);
				Elem(elemID, 3) = StructPart::ijk2ID(i, j, k + 1, numi, numj, numk);
				elemID++;
			}
		}
	}
	for (int i = 0; i < numi; i++)
	{
		for (int j = 0; j < numj - 1; j++)
		{
			for (int k = 0; k < numk - 1; k++)
			{
				Elem(elemID, 0) = StructPart::ijk2ID(i, j, k, numi, numj, numk);
				Elem(elemID, 1) = StructPart::ijk2ID(i, j + 1, k, numi, numj, numk);
				Elem(elemID, 2) = StructPart::ijk2ID(i, j + 1, k + 1, numi, numj, numk);
				Elem(elemID, 3) = StructPart::ijk2ID(i, j, k + 1, numi, numj, numk);
				elemID++;
			}
		}
	}
	StructPE.E = Elem;
	StructPE.P = zeros(numi * numj * numk, 3);
	//节点按列排序
	int id = 0;
	for (int k = 0; k < numk; k++)
	{
		for (int j = 0; j < numj; j++)
		{
			for (int i = 0; i < numi; i++)
			{
				StructPE.P(id, 0) = points[k].X(i, j);
				StructPE.P(id, 1) = points[k].Y(i, j);
				StructPE.P(id, 2) = points[k].Z(i, j);
				id++;
			}
		}
	}

	//----------4、点的旋转和平移----------
	double OriginX = Origin(0);
	double OriginY = Origin(1);
	double OriginZ = Origin(2);

	double RotationX = Rotation(0) / 180.0 * pi;
	double RotationY = Rotation(1) / 180.0 * pi;
	double RotationZ = Rotation(2) / 180.0 * pi;

	//还原角度 通过Rotation Matrix
	mat Lrotation = { { cos(RotationZ) * cos(RotationY) , sin(RotationZ) , -cos(RotationZ) * sin(RotationY) },
					 {-sin(RotationZ) * cos(RotationY) * cos(RotationX) + sin(RotationY) * sin(RotationX), cos(RotationZ) * cos(RotationX)	, sin(RotationZ) * sin(RotationY) * cos(RotationX) + cos(RotationY) * sin(RotationX)},
					 { sin(RotationZ) * cos(RotationY) * sin(RotationX) + sin(RotationY) * cos(RotationX),-cos(RotationZ) * sin(RotationX)	,-sin(RotationZ) * sin(RotationY) * sin(RotationX) + cos(RotationY) * cos(RotationX)} };

	StructPE.P = StructPE.P * Lrotation;
	StructPE.P.col(0) = StructPE.P.col(0) + OriginX;
	StructPE.P.col(1) = StructPE.P.col(1) + OriginY;
	StructPE.P.col(2) = StructPE.P.col(2) + OriginZ;

	//----------5、清除重复网格----------
	cout << "开始清除重复网格......" ;
	uniqueMat do_it(StructPE);//.P, StructPE.E
	StructPE.P = do_it.P_new();
	StructPE.E = do_it.E_new();
	cout << "清除重复网格成功！！！" << endl;
	return 0;
}

int CSTsurface::MakeStruct2D_byNum()
{
	//1、输入数据处理&检测
	NHeight = (NHeight < 2) ? 2 : NHeight; //高度方向至少网格数是2 
	int numi = PUPP.X.n_rows;//基本曲面轴向网格点数量
	int numj = PUPP.X.n_cols;//基本曲面横向网格点数量
	int numk = NHeight;//高度方向网格点数量
	//
	if (PLOW.X.n_elem != PUPP.X.n_elem)
	{
		cout << "MakeStruct2D_byNum错误！！上下网格数不一致！" << endl;
		return -1;
	}
	//2、结构单元点的生成
	node_2D.resize(NHeight);
	for (int i = 0; i < numk; i++)
	{
		node_2D[i].X = PLOW.X + (PUPP.X - PLOW.X) * i / (NHeight - 1.0);
		node_2D[i].Y = -PLOW.Y + (PUPP.Y + PLOW.Y) * i / (NHeight - 1.0);//y方向PLOW缺少负号，故取负
		node_2D[i].Z = PLOW.Z + (PUPP.Z - PLOW.Z) * i / (NHeight - 1.0);
		node_2D[i].ifuse = zeros(PLOW.X.n_rows, PLOW.X.n_cols);//默认取值为0 代表节点不会被使用
	}

	return 0;
}

int CSTsurface::MakeStruct1D(mat X, mat Z)
{
	//----------1、数据处理/检查----------
	int numi = PUPP.X.n_rows;
	int numj = PUPP.X.n_cols;

	if (PLOW.X.n_elem != PUPP.X.n_elem)
	{
		cout << "MakcStruct错误！！上下网格数不一致！" << endl;
		return 0;
	}
	//处理数据X
	double delta_x = 1.0 / (NEta - 1.0);
	mat X_uniq = unique(abs(X)); 
	if (X_uniq(X_uniq.n_rows - 1) > 1)//如果有大于一的数则归一化
	{
		X_uniq = X_uniq / X_uniq(X_uniq.n_rows - 1);
	}
	for (size_t i = 0; i < X_uniq.n_rows; i++)//以delta来过滤向量
	{
		X_uniq(i) = int(X_uniq(i) / delta_x);
	}
	X_uniq = unique(X_uniq);
	//处理数据Z
	double delta_z = 1.0 / (NFaiL - 1.0);
	mat Z_uniq = unique(abs(Z));
	if (Z_uniq(Z_uniq.n_rows - 1) > 1)//如果有大于一的数则归一化
	{
		Z_uniq = Z_uniq / Z_uniq(Z_uniq.n_rows - 1);
	}
	for (size_t i = 0; i < Z_uniq.n_rows; i++)//以delta来过滤向量
	{
		Z_uniq(i) = int(Z_uniq(i) / delta_z);
	}
	Z_uniq = unique(Z_uniq);

	//----------2、气动网格生成------------
	//point1D = (GridUpp.P + GridLow.P) / 2;//节点坐标序列
	mat pt_h = abs(GridUpp.P.col(1) - GridLow.P.col(1));//节点对应的高度
	elem_aero = GridLow.E-1;//气动单元序列
	mat elem_h = zeros(elem_aero.n_rows, 1);//气动单元对应高度
	//cout << elem_aero(4999, 0) << "  " << elem_aero(4999, 1) << "   " << elem_aero(4999, 2) << "   " << elem_aero(4999, 3) << "   " << endl;
	for (size_t i = 0; i < elem_aero.n_rows; i++)
	{
		elem_h(i) = (
			pt_h((size_t)elem_aero(i, 0),0) + 
			pt_h((size_t)elem_aero(i, 1),0) + 
			pt_h((size_t)elem_aero(i, 2),0) + 
			pt_h((size_t)elem_aero(i, 3),0)) / 4;
	}

	//----------3、结构单元生成------------
	int Enum_beam0 = Z_uniq.n_rows * (numj - 1);//旋转前垂直于z方向梁单元数
	int Enum_beam1 = X_uniq.n_rows * (numi - 1);//旋转前垂直于x方向梁单元数
	elem_strc = zeros(Enum_beam0 + Enum_beam1, 4);//规定 一二列是节点编号，第三列是梁的纵横，第四列是id
	
	aero_strc.LinkPtId = zeros(Z_uniq.n_rows, 1);
	for (size_t i = 0; i < Z_uniq.n_rows; i++)//连接梁
	{
		aero_strc.LinkPtId(i) = NEta * Z_uniq(i);//连接点的编号
		for (int j = 0; j < numj - 1; j++)
		{
			elem_strc(i * (numj - 1) + j, 0) = NEta * Z_uniq(i) + j;
			elem_strc(i * (numj - 1) + j, 1) = NEta * Z_uniq(i) + j + 1;
			elem_strc(i * (numj - 1) + j, 2) = 0;
			elem_strc(i * (numj - 1) + j, 3) = i;
		}
	}
	for (size_t j = 0; j < X_uniq.n_rows; j++)//平行梁
	{
		for (int i = 0; i < numi - 1; i++)
		{
			elem_strc(Enum_beam0 + j * (numi - 1) + i, 0) = (double)NEta * i + X_uniq(j);
			elem_strc(Enum_beam0 + j * (numi - 1) + i, 1) = (double)NEta * (i + 1.0) + X_uniq(j);
			elem_strc(Enum_beam0 + j * (numi - 1) + i, 2) = 1;
			elem_strc(Enum_beam0 + j * (numi - 1) + i, 3) = j;
		}
	}

	//----------4、参数传递与输出------------
	// 	   暂时没有重复节点 先不做清除
	//uniqueMat peA(point1D, elem_aero);
	//uniqueMat peS(point1D, elem_strc);
	aero_strc.node = (GridUpp.P + GridLow.P) / 2;// peA.P_new();
	aero_strc.elem_aero = elem_aero;// peA.E_new();
	aero_strc.elem_strc = elem_strc;// peS.E_new(2);
	
	aero_strc.node_f = zeros(aero_strc.node.n_rows, 6);
	mat moment = zeros(aero_strc.node.n_rows, 3);
	if (GridLow.nodeForce.n_rows != 0 && GridUpp.nodeForce.n_rows != 0)
	{
		aero_strc.node_f.cols(0, 2) = (GridLow.nodeForce + GridUpp.nodeForce) / 2;
		for (size_t i = 0; i < aero_strc.node.n_rows; i++)
		{
			moment.row(i) += cross((GridUpp.P.row(i) - aero_strc.node.row(i)), GridUpp.nodeForce.row(i));
			moment.row(i) += cross((GridLow.P.row(i) - aero_strc.node.row(i)), GridLow.nodeForce.row(i));
		}
	}
	
	
	
	aero_strc.node_f.cols(3, 5) = moment;

	aero_strc.elem_h = elem_h;
	aero_strc.Xbar = X_uniq;
	aero_strc.Zbar = Z_uniq;

	return 1;
}


