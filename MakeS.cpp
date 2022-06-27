//#include<iostream>
//#include"armadillo"
//#include"factorial.h"
//using namespace arma;
//using namespace std;
//
//#define pi 3.1415926
//
//////包含CST模型信息的类SETUP
////class SETUP
////{
////public:
////	mat Origin;				//初始位置（X偏移，Y偏移，Z偏移）
////	mat Rotation;			//旋转因子（X 轴 ，Y 轴 ，Z 轴 ）
////	mat Slope;				//倾斜因子（X正向，Y正向，Z正向）
////	mat Scale;				//缩放因子(X比例,YUPP比例,YLOW比例,Z比例)
////	mat Length;				//缩放因子(X比例,YUPP比例,YLOW比例,Z比例)
////	mat NS;					//头部截面形状因子（上表面N1/N2,下表面N1/N2）
////	mat NE;					//尾部截面形状因子（上表面N1/N2,下表面N1/N2）
////	mat M;					//侧面导引形状因子（左导引M1/M2,右导引M1/M2）
////	mat T;					//中面导引形状因子（上导引T1/T2,下导引T1/T2）
////	double Ratio;			//形状因子过渡系数
////	int NFaiU;				//上表面横向网格数
////	int NFaiL;				//下表面横向网格数
////	int NEta;				//轴向网格数
////	int NHeight;			//高度方向网格数
////	int GridRefineType;		//网格加密参照part(-1=不进行加密，0采用完全修正，i采用部分修正)
////	//int Max_GradsRatio;		//最大网格加密梯度
////	bool Is_MeshFront;		//是否缝合头部
////	bool Is_MeshBack;		//是否缝合尾部
////	mat BUPP;				//上表面控制因子
////	mat BLOW;				//下表面控制因子
////	mat DUPP;				//侧向表面控制因子
////	mat DLOW;				//侧向表面控制因子
////
////	//曲率评估参数
////	bool IsEvalCurvature;	//是否进行曲率评估
////	mat CURyx;			//轴向曲率
////	mat CURyzUpp;		//上表面侧向曲率
////	mat CURyzLow;		//下表面侧向曲率
////};
////
//////点类
////class Point
////{
////public:
////	mat X;
////	mat Y;
////	mat Z;
////};
////
//////网格其他信息类
////class MeshInfo
////{
////public:
////	mat NFaiU;
////	mat NFaiL;
////	mat NEta;
////
////	mat CURyx;
////	mat CURyzUpp;
////	mat CURyzLow;
////	//评估网格信息
////	double SurfaceAera;
////	double Volume;
////	double planforarea;
////};
////
//////网格类
////class Grid
////{
////public:
////	MeshInfo Info;			//网格其他信息
////	mat Upp;				//上表面网格
////	mat Low;				//下表面网格
////	mat Front;				//前表面网格
////	mat Back;				//后表面网格
////	
////};
//
////生成形状函数S
//mat MakeS(int method, int K, mat B, mat fai, mat eta)
//{
//	int fainums = fai.n_cols;//向量fai的参数个数
//	int etanums = eta.n_cols;//向量eta的参数个数
//	mat S = zeros(fainums, etanums);
//
//	int NDIM = B.n_rows - 1;//B 权重因子
//	int MDIM = B.n_cols - 1;
//
//	K = min(NDIM, K);
//	K = min(MDIM, K);
//	mat NBF = zeros(NDIM + K + 1, K + 1);
//	mat MBF = zeros(MDIM + K + 1, K + 1);
//
//	for (int i = 0; i < fainums ; i++)
//	{
//		for (int j = 0; j < etanums ; j++)
//		{
//			cout << "==========================================" << endl;//for test
//			cout << "fai i = " << i ;//for test
//			cout << " and eta j = " << j << endl;//for test
//			//方式1：使用bernstrin多项式为基函数
//			if (method == 1)
//			{
//				float Kn;
//				float Km;
//				for (int n = 0; n < NDIM + 1; n++)
//				{
//					Kn = (double)factorial(NDIM) / ((double)factorial(n) * factorial(NDIM - n));
//					for (int m = 0; m < MDIM + 1; m++)
//					{
//						Km = (double)factorial(MDIM) / ((double)factorial(m) * factorial(MDIM - m));
//						S(i, j) = S(i, j) + B(n , m ) * Kn * pow(fai(i), n) * pow((1 - fai(i)), (NDIM - n)) * Km * pow(eta(j), m) * pow((1 - eta(j)), (MDIM - m));
//					}
//				}
//			}
//			//方式2：使用B样条基函数
//			if (method == 2)
//			{
//				for (int k = 0; k < K + 1; k++)
//				{
//					cout << "------------k = " << k << "-------------" << endl;//for test
//					//生成准均匀B样条的节点参数tn
//					mat tnleft = zeros(1, K);
//					mat tnright = ones(1, K);
//					int tnNum = NDIM + K + 2 - 2 * K;
//					mat tn = zeros(1, tnNum);
//					for (int i = 0; i < tnNum; i++)
//					{
//						tn(0, i) = i / ((double)tnNum - 1);
//					}
//					tn = join_rows(tnleft, tn);
//					tn = join_rows(tn, tnright);
//					if (k == 0)
//					{
//						for (int n = 0; n < NDIM + K + 1; n++)
//						{
//							//int kID = k + 1;
//							//int nID = n + 1;
//							if ((tn(n) <= fai(i)) && (fai(i) < tn(n + 1)))
//							{
//								NBF(n, k) = 1;
//							}
//							else
//							{
//								NBF(n, k) = 0;
//							}
//						}
//					}
//					else
//					{
//						for (int n = 0; n < NDIM + K - k + 1; n++)
//						{
//							double alpha = (fai(i) - tn(n)) / (tn(n + k) - tn(n));
//							double beta = (tn(n + k + 1) - fai(i)) / (tn(n + k + 1) - tn(n + 1));
//							//double alpha = (fai(i) - tn(n)) / (tn(n + k -1) - tn(n));
//							//double beta = (tn(n + k) - fai(i)) / (tn(n + k) - tn(n + 1));
//							if (isnan(alpha) || isinf(alpha))
//							{
//								alpha = 0;
//							}
//							if (isnan(beta) || isinf(beta))
//							{
//								beta = 0;
//							}
//							NBF(n, k) = alpha * NBF(n, k - 1) + beta * NBF(n + 1, k - 1);
//							//NBF(n, k) = 0.07*n;//for test
//							cout << "tn	";//for test
//							cout << "alpha = " << alpha << "	beta = " << beta << endl;//for test
//							//cout << "NBF(n ,k) = " << endl << NBF << endl;
//						}
//						
//					}
//					//生成准均匀B样条的节点参数tm
//					mat tmleft = zeros(1, K);
//					mat tmright = ones(1, K);
//					int tmNum = MDIM + K + 2 - 2 * K;
//					mat tm = zeros(1, tmNum);
//					for (int i = 0; i < tmNum; i++)
//					{
//						tm(0, i) = i / ((double)tmNum - 1);
//					}
//					tm = join_rows(tmleft, tm);
//					tm = join_rows(tm, tmright);
//					if (k == 0)
//					{
//						for (int n = 0; n < MDIM + K + 1; n++)
//						{
//							//int kID = k + 1;
//							//int nID = n + 1;
//							if (tm(n) <= eta(j) && eta(j) <= tm(n + 1))
//							{
//								MBF(n, k) = 1;
//							}
//							else
//							{
//								MBF(n, k) = 0;
//							}
//						}
//					}
//					else
//					{
//						for (int m = 0; m < MDIM + K - k + 1; m++)
//						{
//							double alpha = (eta(j) - tm(m)) / (tm(m + k) - tm(m));
//							double beta = (tm(m + k + 1) - eta(j)) / (tm(m + k + 1) - tm(m + 1));
//							if (isnan(alpha) || isinf(alpha))
//							{
//								alpha = 0;
//							}
//							if (isnan(beta) || isinf(beta))
//							{
//								beta = 0;
//							}
//							MBF(m, k) = alpha * MBF(m, k - 1) + beta * MBF(m + 1, k - 1);
//							
//							cout << "	tm	";//for test
//							cout << "alpha = " << alpha << "	beta = " << beta ;//for test
//							cout << "	MBF(m ,k) = " << MBF(m, k) << endl;
//							//MBF(m, k) = 1 * MBF(m, k - 1) + 2 * MBF(m + 1, k - 1);//for test
//						}
//					}
//				}
//				cout << "--------------------" << endl;
//				cout << NBF << endl;
//				cout << "--------------------" << endl;
//				cout << MBF << endl;
//				for (int n = 0; n < NDIM + 1; n++)
//				{
//					for (int m = 0; m < MDIM + 1; m++)
//					{
//						S(i, j) = S(i, j) + B(n, m) * NBF(n, K) * MBF(m, K);
//						
//					}
//				}
//				//cout << "S = " << S << endl;
//			}
//		}
//	}
//	return S;
//}
//
////CST基本数学模型
////Point SketchCST(mat fai, mat eta, double N1, double N2, double N3, double N4, double M1, double M2, double T1, double T2, mat B, mat D, int Ratio)
////{
////	//初始配置
////	int fainums = fai.n_cols;//向量fai的参数个数
////	int etanums = eta.n_cols;//向量eta的参数个数
////	int method = 2;//插值方法：1、Bernstein多项式；2、B样条
////	int K = 2;//形状函数S的次数（阶数-1）。默认为3
////	Point p;
////
////	//纵向修正
////	int NDIM = B.n_rows - 1;//权重因子B 行数-1
////	int MDIM = B.n_cols - 1;//权重因子B 列数-1
////	
////	if (NDIM == 0 && MDIM != 0)
////	{
////		B = repmat(B, K + 1, 1);
////	}
////	else if (MDIM == 0 && NDIM != 0)
////	{
////		B = repmat(B, 1, K + 1);
////	}
////
////	mat S = MakeS(method, K, B, fai, eta);
////
////	//侧向修正
////	NDIM = D.n_rows - 1;//侧向权重因子D 行数-1
////	MDIM = D.n_cols - 1;//侧向权重因子D 列数-1
////
////	if (NDIM == 0 && MDIM != 0)
////	{
////		D = repmat(D, K + 1, 1);
////	}
////	else if (MDIM == 0 && NDIM != 0)
////	{
////		D = repmat(D, 1, K + 1);
////	}
////	mat St = MakeS(method, K, D, fai, eta);
////
////	//区域限制
////	double N_max = N1 / (N1 + N2);
////	double M_max = M1 / (M1 + M2);
////	double T_max = T1 / (T1 + T2);
////	
////	double eta_max = pow(M_max, M1) * pow(1 - M_max, M2);
////	if (isnan(eta_max))
////	{
////		eta_max = 1;
////	}
////
////	double tao_max = pow(T_max, T1) * pow(1 - T_max, T2);
////	if (isnan(tao_max))
////	{
////		tao_max = 1;
////	}
////
////	//修正分析
////	p.X = repmat(eta, fainums, 1);
////	p.Z = zeros(fainums, etanums);
////
////	mat Zmax = zeros(1,etanums);
////	for (int n = 0; n < etanums; n++)
////	{
////		Zmax(0,n) = pow(eta(n), T1) * pow(1 - eta(n), T2) / tao_max;
////	}
////	//mat Zmax = pow(eta, T1) * pow(1 - eta, T2) / tao_max;
////	mat Zmin = -Zmax;
////	for (int i = 0; i < fainums; i++)
////	{
////		for (int j = 0; j < etanums; j++)
////		{
////			p.Z(i, j) = Zmin(j) + (Zmax(j) - Zmin(j)) * fai(i);
////		}
////		//p.Z.row(i) = Zmin + (Zmax - Zmin) * fai(i);
////	}
////	for (int i = 0; i < fainums; i++)
////	{
////		for (int j = 0; j < etanums; j++)
////		{
////			p.Z(i, j) = p.Z(i, j) * St(i, j);
////		}
////	}
////	//p.Z = p.Z * St;
////
////	mat ETA = zeros(1,etanums);
////	for (int n = 0; n < etanums; n++)
////	{
////		ETA(0,n) = pow(eta(n), M1) * pow(1 - eta(n), M2) / eta_max;
////	}
////	//mat ETA = pow(eta, M1) * pow(1 - eta, M2) / eta_max;
////	mat FAI = zeros(fainums, etanums);
////	for (int j = 0; j < etanums; j++)
////	{
////		double N1m = N1 + pow(j / (etanums - 1), Ratio) * ((double)N3 - N1);
////		double N2m = N2 + pow(j / (etanums - 1), Ratio) * ((double)N4 - N2);
////		double N_max = N1m / (N1m + N2m);
////		double fai_max = pow(N_max, N1m) * pow(1 - N_max, N2m);
////		if (isnan(fai_max))
////		{
////			fai_max = 1;
////		}
////		for (int n = 0; n < fainums; n++)
////		{
////			FAI(n,j) = pow(fai(n), N1m) * pow(1 - fai(n), N2m) * ETA(0,j) / fai_max;
////		}
////		//FAI.col(j) = pow(fai, N1m) * pow(1 - fai, N2m) * ETA(j) / fai_max;
////	}
////	p.Y = zeros(fainums, etanums);
////	for (int i = 0; i < fainums; i++)
////	{
////		for (int j = 0; j < etanums; j++)
////		{
////			p.Y(i,j) = FAI(i,j) * S(i,j);
////		}
////	}
////	//p.Y = FAI * S;
////
////	return p;
////}
////
//////CST实际建模模型
////Grid BaceShapeCST(SETUP setup)
////{
////	Grid grid;
////	
////	//参数传入
////	double OriginX = setup.Origin(1);
////	double OriginY = setup.Origin(2);
////	double OriginZ = setup.Origin(3);
////
////	double RotationX = setup.Rotation(1) / 180 * pi;
////	double RotationY = setup.Rotation(2) / 180 * pi;
////	double RotationZ = setup.Rotation(3) / 180 * pi;
////	
////	double Xlength = setup.Length(1);
////	double Ylength_upp = setup.Length(2);
////	double Ylength_low = setup.Length(3);
////	double Zlength = setup.Length(4);
////
////	double SlopeX = setup.Slope(1) / 180 * pi;
////	double SlopeY = setup.Slope(2) / 180 * pi;
////	double SlopeZ = setup.Slope(3) / 180 * pi;
////
////	double ScaleX = setup.Scale(1);
////	double ScaleYUPP = setup.Scale(2);
////	double ScaleYLOW = setup.Scale(3);
////	double ScaleZ = setup.Scale(4);
////
////	//类别函数控制参数
////	double N1s = setup.NS(1);	double N2s = setup.NS(2);	double N3s = setup.NS(3);	double N4s = setup.NS(4);	//截面
////	double N1e = setup.NE(1);	double N2e = setup.NE(2);	double N3e = setup.NE(3);	double N4e = setup.NE(4);	//截面
////	double M1 = setup.M(1);		double M2 = setup.M(2);		double M3 = setup.M(3);		double M4 = setup.M(4);		//轴线
////	double T1 = setup.T(1);		double T2 = setup.T(2);		double T3 = setup.T(3);		double T4 = setup.T(4);		//侧向
////
////	mat BUPP = setup.BUPP;    mat BLOW = setup.BLOW;//上下表面形状控制因子
////	mat DUPP = setup.DUPP;    mat DLOW = setup.DLOW;//左右表面形状控制因子
////	double Ratio = setup.Ratio;
////
////	//int GridRefinrType = setup.GridRefineType;
////	bool Is_MeshFront = setup.Is_MeshFront;
////	bool Is_MeshBack = setup.Is_MeshBack;
////
////	//节点个数
////	int NFaiU = setup.NFaiU;	int NFaiL = setup.NFaiL;	int NEta = setup.NEta;	int NHeight = setup.NHeight;
////	
////	
////	//网格数修正
////	if (NFaiU % 2 == 0)
////	{
////		NFaiU = NFaiU + 1;
////	}
////	if (NFaiL % 2 == 0)
////	{
////		NFaiL = NFaiL + 1;
////	}
////	if (NEta % 2 == 0)
////	{
////		NEta = NEta + 1;
////	}
////	mat FaiU = linspace(0, 1, NFaiU);
////	mat FaiL = linspace(0, 1, NFaiL);
////	mat Eta = linspace(0, 1, NEta);
////
////	grid.Info.NFaiU = NFaiU;
////	grid.Info.NFaiL = NFaiL;
////	grid.Info.NEta = NEta;
////
////	//曲率评估
////	if (setup.IsEvalCurvature)
////	{
////		//初始分布计算
////		Point PUPP = SketchCST(FaiU, Eta, N1s, N2s, N1e, N2e, M1, M2, T1, T2, BUPP, DUPP, Ratio);
////		Point PLOW = SketchCST(FaiL, Eta, N3s, N3s, N3e, N4e, M3, M4, T3, T4, BLOW, DLOW, Ratio);
////
////		//距离修正
////		mat XYDispUpp = zeros(NFaiU, NEta - 1);
////		mat XYDispLow = zeros(NFaiL, NEta - 1);
////		mat YZDispUpp = zeros(NFaiU - 1, NEta);
////		mat YZDispLow = zeros(NFaiU - 1, NEta);
////		
////		for (int i = 0; i < NFaiU; i++)
////		{
////			for (int j = 0; j < NEta - 1; j++)
////			{
////				XYDispUpp(i, j) = sqrt(pow(PUPP.X(i, j + 1) - PUPP.X(i, j), 2) + pow(PUPP.Y(i, j + 1) - PUPP.Y(i, j), 2));
////			}
////		}
////		for (int i = 0; i < NFaiL; i++)
////		{
////			for (int j = 0; j < NEta - 1; j++)
////			{
////				XYDispLow(i, j) = sqrt(pow(PLOW.X(i, j + 1) - PLOW.X(i, j), 2) + pow(PLOW.Y(i, j + 1) - PLOW.Y(i, j), 2));
////			}
////		}
////		for (int j = 0; j < NEta; j++)
////		{
////			for (int i = 0; i < NFaiU - 1; i++)
////			{
////				YZDispUpp(i, j) = sqrt(pow(PUPP.Z(i + 1, j) - PUPP.Z(i, j), 2) + pow(PUPP.Y(i + 1, j) - PUPP.Y(i, j), 2));
////			}
////			for (int i = 0; i < NFaiL - 1; i++)
////			{
////				YZDispLow(i, j) = sqrt(pow(PLOW.Z(i + 1, j) - PLOW.Z(i, j), 2) + pow(PLOW.Y(i + 1, j) - PLOW.Y(i, j), 2));
////			}
////		}
////		mat CURyx = join_cols(XYDispUpp, XYDispLow);
////		CURyx = max(CURyx);//寻找每列的最大值，存成一行
////		mat CURyzUpp = max(YZDispUpp, 1);//寻找行的最大值，存成一列
////		mat CURyzLow = max(YZDispLow, 1);//寻找行的最大值，存成一列    %%这里行列有些混乱，出错需检查
////
////		//曲率计算
////		for (int i = 0; i < CURyx.n_elem; i++)
////		{
////			if (CURyx(i) == 0)
////			{
////				CURyx(i) = 1;
////			}
////		}
////		for (int i = 0; i < CURyzUpp.n_elem; i++)
////		{
////			if (CURyzUpp(i) == 0)
////			{
////				CURyzUpp(i) = 1;
////			}
////		}
////		for (int i = 0; i < CURyzLow.n_elem; i++)
////		{
////			if (CURyzLow(i) == 0)
////			{
////				CURyzLow(i) = 1;
////			}
////		}
////
////		mat MaxMin;//计算下列向量的最大最小值，MaxMin为中间变量
////		MaxMin = max(CURyx, 1);
////		double CURyxMax = MaxMin(0, 0);
////		MaxMin = min(CURyx, 1);
////		double CURyxMin = MaxMin(0, 0);
////		MaxMin = max(CURyzUpp);
////		double CURyzUppMax = MaxMin(0, 0);
////		MaxMin = min(CURyzUpp);
////		double CURyzUppMin = MaxMin(0, 0);
////		MaxMin = max(CURyzLow);
////		double CURyzLowMax = MaxMin(0, 0);
////		MaxMin = min(CURyzLow);
////		double CURyzLowMin = MaxMin(0, 0);
////								
////		double DallowMax = min(CURyxMax, 10.0 / (double)NEta);
////		double DallowMin = max(CURyxMin, 0.1 / (double)NEta);
////		grid.Info.CURyx = zeros(1, CURyx.n_elem);
////		for (int i = 0; i < CURyx.n_elem; i++)
////		{
////			grid.Info.CURyx(i) = (CURyx(i) - CURyxMin) / (CURyxMax - CURyxMin) * (DallowMax - DallowMin) + DallowMin;
////		}
////
////		DallowMax = min(CURyzUppMax, 10.0 / (double)NFaiU);
////		DallowMin = max(CURyzUppMin, 0.1 / (double)NFaiU);
////		grid.Info.CURyzUpp = zeros(CURyzUpp.n_elem, 1);
////		for (int i = 0; i < CURyzUpp.n_elem; i++)
////		{
////			grid.Info.CURyzUpp(i) = (CURyzUpp(i) - CURyzUppMin) / (CURyzUppMax - CURyzUppMin) * (DallowMax - DallowMin) + DallowMin;
////		}
////
////		DallowMax = min(CURyzLowMax, 10.0 / (double)NFaiL);
////		DallowMin = max(CURyzLowMin, 0.1 / (double)NFaiL);
////		grid.Info.CURyzLow = zeros(CURyzLow.n_elem, 1);
////		for (int i = 0; i < CURyzLow.n_elem; i++)
////		{
////			grid.Info.CURyzLow(i) = (CURyzLow(i) - CURyzLowMin) / (CURyzLowMax - CURyzLowMin) * (DallowMax - DallowMin) + DallowMin;
////		}
////		return grid;
////	}
////	else if (setup.GridRefineType >= 0)//是否采用网格修正
////	{
////		Eta = zeros(1, setup.CURyx.n_elem+1);
////		Eta(0) = 0;
////		for (int i = 1; i < setup.CURyx.n_elem+1; i++)
////		{
////			Eta(i) = 1 / setup.CURyx(i);
////		}
////		for (int i = 1; i < Eta.n_elem; i++)
////		{
////			mat Temp = max(Eta, 2);
////			double temp = Temp(0, 0);
////			Eta(i) = Eta(i) / temp;
////		}
////
////		FaiU = zeros(1, setup.CURyzUpp.n_elem + 1);
////		FaiU(0) = 0;
////		for (int i = 1; i < setup.CURyzUpp.n_elem + 1; i++)
////		{
////			FaiU(i) = 1 / setup.CURyzUpp(i);
////		}
////		for (int i = 1; i < FaiU.n_elem; i++)
////		{
////			mat Temp = max(FaiU, 2);
////			double temp = Temp(0, 0);
////			FaiU(i) = FaiU(i) / temp;
////		}
////
////		FaiL = zeros(1, setup.CURyzLow.n_elem + 1);
////		FaiL(0) = 0;
////		for (int i = 1; i < setup.CURyzLow.n_elem + 1; i++)
////		{
////			FaiL(i) = 1 / setup.CURyzLow(i);
////		}
////		for (int i = 1; i < FaiL.n_elem; i++)
////		{
////			mat Temp = max(FaiL, 2);
////			double temp = Temp(0, 0);
////			FaiL(i) = FaiL(i) / temp;
////		}
////	}
////
////	//CST 3D
////	Point PUPP = SketchCST(FaiU, Eta, N1s, N2s, N1e, N2e, M1, M2, T1, T2, BUPP, DUPP, Ratio);
////	Point PLOW = SketchCST(FaiL, Eta, N3s, N4s, N3e, N4e, M3, M4, T3, T4, BLOW, DLOW, Ratio);
////
////	//Scale X Y Z
////	//还原为标准尺寸 normal
////	PUPP.X = PUPP.X * Xlength;
////	PUPP.Y = PUPP.Y * Ylength_upp;
////	PUPP.Z = PUPP.Z * Zlength / 2;
////	PLOW.X = PLOW.X * Xlength;
////	PLOW.Y = PLOW.Y * Ylength_low;
////	PLOW.Z = PLOW.Z * Zlength / 2;
////
////	//还原比例刻度 （用于翼的跟梢比）scale
////	for (int i = 0; i < NEta; i++)
////	{
////		PUPP.X.col(i) = PUPP.X.col(i) * ((ScaleX - 1) * Eta(i) + 1);
////		PUPP.Y.col(i) = PUPP.Y.col(i) * ((ScaleYUPP - 1) * Eta(i) + 1);
////		PUPP.Z.col(i) = PUPP.Z.col(i) * ((ScaleZ - 1) * Eta(i) + 1);
////		PLOW.X.col(i) = PLOW.X.col(i) * ((ScaleX - 1) * Eta(i) + 1);
////		PLOW.Y.col(i) = PLOW.Y.col(i) * ((ScaleYLOW - 1) * Eta(i) + 1);
////		PLOW.Z.col(i) = PLOW.Z.col(i) * ((ScaleZ - 1) * Eta(i) + 1);
////	}
////
////	//还原图形斜率（用于后掠角） slope
////	if (SlopeY != 0)
////	{
////		for (int i = 0; i < PUPP.Y.n_rows; i++)
////		{
////			for (int j = 0; j < PUPP.Y.n_cols; j++)
////			{
////				PUPP.Y(i, j) = PUPP.Y(i, j) + PUPP.X(i, j) * tan(SlopeY);
////			}
////		}
////		for (int i = 0; i < PLOW.Y.n_rows; i++)
////		{
////			for (int j = 0; j < PLOW.Y.n_cols; j++)
////			{
////				PLOW.Y(i, j) = PLOW.Y(i, j) + PLOW.X(i, j) * tan(SlopeY);
////			}
////		}
////	}
////
////	if (SlopeZ != 0)
////	{
////		for (int i = 0; i < PUPP.Z.n_rows; i++)
////		{
////			for (int j = 0; j < PUPP.Z.n_cols; j++)
////			{
////				mat Temp = repmat(PUPP.Z.row(1), NFaiU, 1);
////				PUPP.Z(i, j) = PUPP.Z(i, j) + PUPP.X(i, j) * tan(SlopeZ) - Temp(j) - Zlength / 2;
////			}
////		}
////		for (int i = 0; i < PLOW.Z.n_rows; i++)
////		{
////			for (int j = 0; j < PLOW.Z.n_cols; j++)
////			{
////				mat Temp = repmat(PLOW.Z.row(1), NFaiL, 1);
////				PLOW.Z(i, j) = PLOW.Z(i, j) + PLOW.X(i, j) * tan(SlopeZ) - Temp(j) - Zlength / 2;
////			}
////		}
////	}
////
////	//评估网格信息
////	double VolumeBody = 0;//体积
////	double SurfaceAeraBody = 0;//zy平面截面面积
////	double planformarea = 0;//xz平面面积
////
////	for (int ix = 0; ix < NEta - 1; ix++)
////	{
////		double area = 0;
////		for (int iu = 0; iu < NFaiU - 1; iu++)
////		{
////			area = area + (-1) * PUPP.Z(iu, ix) * PUPP.Y(iu + 1, ix) - (-1) * PUPP.Z(iu + 1, ix) * PUPP.Y(iu, ix);
////			SurfaceAeraBody = SurfaceAeraBody + (sqrt(pow(PUPP.Z(iu, ix) - PUPP.Z(iu + 1, ix), 2) + pow(PUPP.Y(iu, ix) - PUPP.Y(iu + 1, ix), 2))
////				+ sqrt(pow(PUPP.Z(iu, ix + 1) - PUPP.Z(iu + 1, ix + 1), 2) + pow(PUPP.Y(iu, ix + 1) - PUPP.Y(iu + 1, ix + 1), 2))) / 2
////				* (PUPP.X(1, ix + 1) - PUPP.X(1, ix));
////			planformarea = planformarea + (PUPP.X(1, ix + 1) - PUPP.X(1, ix)) * abs(PUPP.Z(iu + 1, ix) - PUPP.Z(iu, ix));
////		}
////		for (int il = NFaiL; il > 1; il--)
////		{
////			area = area + (-1) * PLOW.Z(il, ix) * PLOW.Y(il - 1, ix) - (-1) * PLOW.Z(il - 1, ix) * PLOW.Y(il, ix);
////			SurfaceAeraBody = SurfaceAeraBody + (sqrt(pow(PLOW.Z(il, ix) - PLOW.Z(il - 1, ix), 2) + pow(PLOW.Y(il, ix) - PLOW.Y(il - 1, ix), 2))
////				+ sqrt(pow(PLOW.Z(il, ix + 1) - PLOW.Z(il - 1, ix + 1), 2) + pow(PLOW.Y(il, ix + 1) - PLOW.Y(il - 1, ix + 1), 2))) / 2
////				* (PLOW.X(1, ix + 1) - PLOW.X(1, ix));
////		}
////		VolumeBody = VolumeBody + (PUPP.X(1, ix + 1) - PUPP.X(1, ix)) * area / 2;
////	}
////	grid.Info.SurfaceAera = abs(SurfaceAeraBody);
////	grid.Info.Volume = abs(VolumeBody);
////	grid.Info.planforarea = abs(planformarea);
////
////
////	
////	return grid;
////}
//
////测试函数MakeS
//void testMakeS()
//{
//	//mat S = MakeS(2, 0,
//	//	{ {1},
//	//	  
//	//	
//	//	  {1} },
//	//	{ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1},
//	//	{ 0});
//
//	mat S = MakeS(2,2, 
//		{ {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7},
//		  {1,2,3,4,5,6,7} },
//		{ 0,0.1,0.2,0.3,0.4,0.5,0.6 },
//		{ 0,0.1,0.2,0.3,0.4,0.5,0.6 });
//
//	cout << "testMakeS:" << endl;
//	cout << S << endl;
//}
//
//double  N(double t, int i, int k, mat T)
//{
//	double Nr;
//	if (k == 1)
//	{
//		if (T(i) <= t || t <= T(i + 1))
//		{
//			Nr = 1;
//		}
//		else
//		{
//			Nr = 0;
//		}
//	}
//	else
//	{
//		double alpha = (t - T(i)) / (T(i + k - 1) - T(i));
//		double beta = (T(i + k) - t) / (T(i + k) - T(i + 1));
//		if (isnan(alpha) || isinf(alpha))
//		{
//			alpha = 0;
//		}
//		if (isnan(beta) || isinf(beta))
//		{
//			beta = 0;
//		}
//		Nr = alpha * N(t, i, k - 1, T) + beta * N(t, i + 1, k - 1, T);
//	}
//	//cout << "N(" << i << ", " << k << ") = " << Nr << endl;
//	return Nr;
//}
//
//mat MakeSself(int K, mat B, mat fai, mat eta)
//{
//	int fainums = fai.n_cols;//向量fai的参数个数
//	int etanums = eta.n_cols;//向量eta的参数个数
//	mat S = zeros(fainums, etanums);
//
//	int NDIM = B.n_rows - 1;//B 权重因子
//	int MDIM = B.n_cols - 1;
//
//	K = min(NDIM, K);
//	K = min(MDIM, K);
//	mat NBF = zeros(NDIM + K + 1, K + 1);
//	mat MBF = zeros(MDIM + K + 1, K + 1);
//
//
//	//生成准均匀B样条的节点参数Tn
//	mat tnleft = zeros(1, K);
//	mat tnright = ones(1, K);
//	int tnNum = NDIM + K + 2 - 2 * K;
//	mat Tn = zeros(1, tnNum);
//	for (int i = 0; i < tnNum; i++)
//	{
//		Tn(0, i) = i / ((double)tnNum - 1);
//	}
//	Tn = join_rows(tnleft, Tn);
//	Tn = join_rows(Tn, tnright);
//
//	//生成准均匀B样条的节点参数Tm
//	mat tmleft = zeros(1, K);
//	mat tmright = ones(1, K);
//	int tmNum = MDIM + K + 2 - 2 * K;
//	mat Tm = zeros(1, tmNum);
//	for (int i = 0; i < tmNum; i++)
//	{
//		Tm(0, i) = i / ((double)tmNum - 1);
//	}
//	Tm = join_rows(tmleft, Tm);
//	Tm = join_rows(Tm, tmright);
//
//	for (int i = 0; i < fainums; i++)
//	{
//		for (int j = 0; j < etanums; j++)
//		{
//			cout << "---------------------------------------" << endl;
//			cout << "(i,j) = " << i << "," << j << endl;
//			for (int n = 0; n < NDIM + 1; n++)
//			{
//				for (int m = 0; m < MDIM + 1; m++)
//				{
//					
//					cout << B(n, m) * N(fai(i), n, K, Tn) * N(eta(j), m, K, Tm) << endl;
//					S(i, j) = S(i, j) + B(n, m) * N(fai(i), n, K, Tn) * N(eta(j), m, K, Tm);
//
//				}
//			}
//		}
//	}
//	return S;
//}
//
////测试函数MakeS
//void testMakeSelf()
//{
//	//mat S = MakeSself(3, 
//	//	{ {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7},
//	//	  {1,2,3,4,5,6,7} },
//	//	{ 0,0.1,0.2,0.3,0.4,0.5,0.6 },
//	//	{ 0,0.1,0.2,0.3,0.4,0.5,0.6 });
//
//	mat S = MakeSself(3,
//		{ {1,2,3,4,5},
//		  {1,2,3,4,5},
//		  {1,2,3,4,5},
//		  {1,2,3,4,5},
//		  {1,2,3,4,5} },
//		{ 0,0.1,0.2,0.3,0.4,0.5,0.6 },
//		{ 0,0.1,0.2,0.3,0.4,0.5,0.6 });
//
//	cout << "testMakeSself:" << endl;
//	cout << S << endl;
//}
//
//
//
//////测试函数SketchCST
////void testSketchCST()
////{
////	Point p = SketchCST
////	(
////		{ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 },
////		{ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 },
////		0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, { 1 }, { 1 }, 1
////	);
////	cout << "test SketchCST" << endl;
////	cout << "this is X:" << endl;
////	cout << p.X << endl;
////	cout << "this is Y:" << endl;
////	cout << p.Y << endl;
////	cout << "this is Z:" << endl;
////	cout << p.Z << endl;
////}
//
////int main()
////{
////	testMakeS();
////
////	//testMakeSelf();
////
////	int K = 0;
////
////	mat tnleft = zeros(1, K);
////	mat tnright = ones(1, K);
////	int tnNum = 1 + K + 2 - 2 * K;
////	mat tn = zeros(1, tnNum);
////	for (int i = 0; i < tnNum; i++)
////	{
////		tn(0, i) = i / ((double)tnNum - 1);
////	}
////	tn = join_rows(tnleft, tn);
////	tn = join_rows(tn, tnright);
////
////	cout << tn << endl;
////
////	system("pause");
////
////	return 0;
////}