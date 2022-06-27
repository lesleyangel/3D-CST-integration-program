#include "CST_Instantiation.h"
#include <fstream>
#include <iomanip>
// #include<boost/assign.hpp>
// #include<boost/foreach.hpp>
#include<direct.h>
// using namespace boost::assign;

//将vector<double>参数转换到mat
void AbstructShape::copyM(vector<vector<double>>& v, mat& m)
{
	vector<vector<double>>::iterator s = v.begin();
	//cout << v.size()<<" "<<(*s).size() << endl;//for test
	m = zeros(v.size(), (*s).size());
	//cout << "初始化 m：" << endl;//for test
	//cout << m << endl;//for test
	int i = 0;
	for (vector<vector<double>>::iterator it = v.begin(); it != v.end(); it++)
	{
		int j = 0;
		//(*it)---容器 vector<int>
		for (vector<double>::iterator vit = (*it).begin(); vit != (*it).end(); vit++)
		{
			//cout << *vit << " ";
			m(i, j) = *vit;
			j++;
		}
		i++;
	}
	//cout << "结束的 m：" << endl;//for test
	//cout << m << endl;//for test
}

void AbstructShape::SaveTecplotMesh()
{
	ofstream ofs;
	ofs.open("test.dat", ios::app);
	for(CSTsurface i: CSTsf)
	{
		ofs << "Title = \"example: 2d finite - element data\"" << endl;
		ofs << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs << "zone N=" << i.GridUpp.P.n_rows << ",E=" << i.GridUpp.E.n_rows << ",datapacking=block" << endl;
		ofs << "zonetype=fequadrilateral" << endl;
		ofs << trans(i.GridUpp.P) << endl;
		ofs << i.GridUpp.E << endl;

		ofs << "Title = \"example: 2d finite - element data\"" << endl;
		ofs << "variables = \"x\",\"y\",\"z\"" << endl;
		ofs << "zone N=" << i.GridLow.P.n_rows << ",E=" << i.GridLow.E.n_rows << ",datapacking=block" << endl;
		ofs << "zonetype=fequadrilateral" << endl;
		ofs << trans(i.GridLow.P) << endl;
		ofs << i.GridLow.E << endl;
	}
	ofs.close();
}

//////////////////////////////////////////////
//默认加密类型为不加密
AbstructShape::AbstructShape()
{
	GridRefineType = 0;
	ShapeType = 0;
}

//对曲面网格进行曲率修正
void AbstructShape::BacsShapeCSTs()
{
	int PartNum = CSTsf.size();
	GdStruct.reserve(CSTsf.size());

	if (GridRefineType > PartNum || GridRefineType < -1)
	{
		cout << "ERROR: Woring GridRefineType (value range: [-1, PartNum])" << endl;
		system("pause");
		exit(1);//输入错误，挂掉程序
	}
	else if (GridRefineType == -1)//不进行加密
	{
		cout << "正在运行：不进行加密计算。。。" << endl;
		NotGridRefine();
	}
	else if (GridRefineType == 0)//采用完全修正
	{
		cout << "正在运行：完全修正计算。。。" << endl;
		AllGridRefine();
	}
	else//采用部分修正
	{
		cout << "正在运行：部分修正计算。。。" << endl;
		PrtGridRefine();
	}
}

//不进行加密
void AbstructShape::NotGridRefine()
{
	const size_t PartNum = CSTsf.size();
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].CST3D();
		//CSTsf[iPart].RefineMesh(2, true);//
		CSTsf[iPart].RefineMesh(2, false);//
		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}

//采用完全修正
void AbstructShape::AllGridRefine()
{
	const size_t PartNum = CSTsf.size();

	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		CSTsf[iPart].EvalCurvature();
	}
	//比较每个CSTsf的曲率最大值
	mat CURyzUpp = CSTsf[0].Info.CURyzUpp;
	mat CURyzLow = CSTsf[0].Info.CURyzLow;

	for (size_t iPart = 1; iPart < PartNum; iPart++)
	{
		//CURyzUpp = max(CURyzUpp, CSTsf[iPart].Info.CURyzUpp);
		//CURyzLow = max(CURyzLow, CSTsf[iPart].Info.CURyzLow);
		for (arma::uword i = 0; i < CURyzUpp.n_elem; i++)
		{
			CURyzUpp(i) = max(CURyzUpp(i), CSTsf[iPart].Info.CURyzUpp(i));
		}
		for (arma::uword i = 0; i < CURyzLow.n_elem; i++)
		{
			CURyzLow(i) = max(CURyzLow(i), CSTsf[iPart].Info.CURyzLow(i));
		}
	}
	cout << "正在运行：曲率分析已完成！" << endl;
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].Info.CURyzUpp = CURyzUpp;
		CSTsf[iPart].Info.CURyzLow = CURyzLow;
		CSTsf[iPart].GridRefine();
		CSTsf[iPart].CST3D();

		//BoneLine.push_back(CSTsf[iPart].MakeBone());//******4.12

		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}

//采用部分修正
void AbstructShape::PrtGridRefine()
{
	const size_t PartNum = CSTsf.size();
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 1;
		CSTsf[iPart].EvalCurvature();
	}
	MeshInfo Info0 = CSTsf[GridRefineType - 1].Info;
	for (size_t iPart = 0; iPart < PartNum; iPart++)
	{
		//CSTsf[iPart].IsEvalCurvature = 0;
		CSTsf[iPart].Info.CURyzUpp = Info0.CURyzUpp;
		CSTsf[iPart].Info.CURyzLow = Info0.CURyzLow;
		CSTsf[iPart].GridRefine();
		CSTsf[iPart].CST3D();

		//BoneLine.push_back(CSTsf[iPart].MakeBone());//******4.12

		GridStruct temp;
		temp.GridUpp = CSTsf[iPart].GridUpp;
		temp.GridLow = CSTsf[iPart].GridLow;
		temp.NFaiU = CSTsf[iPart].NFaiU;
		temp.NFaiL = CSTsf[iPart].NFaiL;
		temp.NEta = CSTsf[iPart].NEta;
		temp.Info = CSTsf[iPart].Info;
		GdStruct.push_back(temp);
	}
}
/////////////////////////////////////////////



void AbstructShape::SaveToAero(string name)
{
	ofstream ofs;
	ofs.open(name, ios::trunc);

	ofs << "$..~~~~~~GFFF制作~~~~~" << endl;
	ofs << "$..~~~~~~只能输出面网格（三节点三角形和四节点四边形）~~~~~" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$..         NODES" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$..CAA" << endl;

	int num = 0;
	for (size_t i = 0; i < CSTsf.size(); i++)//上表面
	{
		for (arma::uword j = 0; j < CSTsf[i].GridUpp.P.n_rows; j++)
		{
			num++;

			ofs << "GRID*" << setw(19) << num << setw(32) << CSTsf[i].GridUpp.P(j, 0)
				<< setw(16) << CSTsf[i].GridUpp.P(j, 1) << setw(8) << num << endl;
			ofs << "*" << setw(7) << num << setw(16) << CSTsf[i].GridUpp.P(j, 2) << endl;
		}
	}
	for (size_t i = 0; i < CSTsf.size(); i++)//下表面
	{
		for (size_t j = 0; j < CSTsf[i].GridLow.P.n_rows; j++)
		{
			num++;

			ofs << "GRID*" << setw(19) << num << setw(32) << CSTsf[i].GridLow.P(j, 0)
				<< setw(16) << CSTsf[i].GridLow.P(j, 1) << setw(8) << num << endl;
			ofs << "*" << setw(7) << num << setw(16) << CSTsf[i].GridLow.P(j, 2) << endl;
		}
	}

	ofs << "$..~~~~~~网格生成由GFFF制作~~~~~" << endl;
	ofs << "$----------------------------------------" << endl;
	ofs << "$.." << endl;
	ofs << "$..          ELEMENTS" << endl;
	ofs << "$.." << endl;
	ofs << "$.." << endl;
	ofs << "$----------------------------------------" << endl;
	ofs << "$.." << endl;
	ofs << "$..          MESH PART: MSHPartSmartSurf.1(CAA)" << endl;
	ofs << "$.." << endl;

	num = 0;
	int pointNum = 0;
	for (size_t i = 0; i < CSTsf.size(); i++)
	{
		
		for (size_t j = 0; j < CSTsf[i].GridUpp.E.n_rows; j++)
		{
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1
				<< setw(8) << CSTsf[i].GridUpp.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 1) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 2) + pointNum << endl;
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1
				<< setw(8) << CSTsf[i].GridUpp.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 2) + pointNum
				<< setw(8) << CSTsf[i].GridUpp.E(j, 3) + pointNum << endl;
		}
		pointNum += CSTsf[i].GridUpp.P.n_rows;
	}
	for (size_t i = 0; i < CSTsf.size(); i++)
	{

		for (size_t j = 0; j < CSTsf[i].GridLow.E.n_rows; j++)
		{
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1 + CSTsf.size()
				<< setw(8) << CSTsf[i].GridLow.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 1) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 2) + pointNum << endl;
			num++;
			ofs << "CTRIA3" << setw(10) << num << setw(8) << i + 1 + CSTsf.size()
				<< setw(8) << CSTsf[i].GridLow.E(j, 0) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 2) + pointNum
				<< setw(8) << CSTsf[i].GridLow.E(j, 3) + pointNum << endl;
		}
		pointNum += CSTsf[i].GridLow.P.n_rows;
	}

	ofs << "$.." << endl;
	ofs << "ENDDATA" << endl;
	ofs << endl;
	ofs << "本程序由GFFF制作--航天学院 西北工业大学" << endl;
}

void AbstructShape::MakeAeroInfo(string aeropath,int num)
{
	//string name = "mesh//mesh" + to_string(num) + "指令.txt";
	string datpath = aeropath + "mesh" + to_string(num) + ".dat";//输入网格文件
	string name = aeropath + "mesh" + to_string(num) + "_cmd.txt";//输入指令文件
	string outputpath = aeropath + "mesh" + to_string(num);//输出文件地址
	int state = _mkdir(outputpath.c_str());
	int size = CSTsf.size() * 2;
	ofstream ofs;
	ofs.open(name, ios::trunc);
	ofs << "//BEGIN//" << endl;
	ofs << endl;
	ofs << "//////////////////计算状态(数组）///////////////////////" << endl;
	ofs << "mach	nAlpha   nBeta  nH x轴方向(1为头部指向尾部，-1为反方向)" << endl;
	ofs << "1       1        1      1  1" << endl;
	ofs << "mach    //马赫数" << endl;
	ofs << "10 " << endl;
	ofs << "alpha	//攻角(deg)" << endl;
	ofs << "30 " << endl;
	ofs << "beta 	//侧滑角(deg)" << endl;
	ofs << "0 " << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////////飞行高度(m)////////////////////////" << endl;
	ofs << "20000" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////////分区类型///////////////////////////" << endl;
	ofs << "计算分区数 " << size << endl;
	ofs << "//////计算分区ID" << endl;
	for (int i = 0; i < size; i++)	ofs << i + 1 << " ";
	ofs << endl;
	ofs << "//////计算分区区域" << endl;
	for (int i = 0; i < size; i++)	ofs << 2 << " ";
	ofs << endl;
	ofs << "//////计算分区类型" << endl;
	for (int i = 0; i < size; i++)	ofs << 1 << " ";
	ofs << endl;
	ofs << "///////分区外法线方向" << endl;
	for (int i = 0; i < size; i++)
	{
		int mid = size / 2;
		int direction = -1;
		if (i < mid)
		{
			direction = 1;
		}
		ofs << direction << " ";
	}
	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////控制面////////////////////////////" << endl;
	ofs << "控制面总数 0" << endl;
	ofs << "/////每个控制面设置占一行，分区请与分区设置对应，固定旋转轴(0无，1X,2Y,3Z)，坐标之间不要留空格，由逗号隔开，旋转轴由P2->P1，逆时针旋转为正，单位m" << endl;
	ofs << "^^^对应分区^^^^^^^旋转轴(点1坐标)^^^^^^^旋转轴(点2坐标)^^^^^^^固定旋转轴^^^舵偏角数^^^舵偏角数组^^^^^^^^^^^" << endl;
	ofs << "" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////sRef参考面积//cRef参考长度//参考尺寸(m)////////" << endl;
	ofs << "1 1 1 0 0" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "/////////////////是否考虑粘性(是为1，否为0)/////////////" << endl;
	ofs << "启用粘性计算(0,1)    0" << endl;
	ofs << "头部半径(m)           0.03" << endl;
	ofs << "壁面温度(K)           300" << endl;
	ofs << "壁面辐射率            0.7" << endl;
	ofs << "气体模型(0,1)        1" << endl;
	ofs << "换热系数              0" << endl;
	ofs << "介质温度(K)           280" << endl;
	ofs << "内壁面温度(K)         350" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "////////////////非定常计算选项(是为1，否为0)////////////" << endl;
	ofs << "非定常开启   计算参数:时间间隔(t)    节点速度rpt文档 " << endl;
	ofs << "0            0.0005                  velocity.rpt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////相关文件//////////////////////////" << endl;
	ofs << "//网格文件单位(m/mm):		mm" << endl;
	ofs << "//网格类型                  CATIA" << endl;
	ofs << "网格文件：                  " << datpath << endl;//C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh"
	ofs << "输出路径：                  " << outputpath << endl;//C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh\\mesh"
	ofs << "输出结果文件名：            CalResult - 全弹箭状态.txt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "是否输出压力分布文件plot.plt(plot3D格式):     YES" << endl;
	ofs << "是否输出面元力分布文件face.dat:               NO" << endl;
	ofs << "是否输出节点力分布文件node.dat:               YES" << endl;
	ofs << "是否输出温度分布文件temperature.plt(需要粘性计算) NO" << endl;
	ofs << "是否输出热载分布文件qload.plt(需要粘性计算)       NO" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << endl;
	ofs << endl;
	ofs <<  endl;
	ofs << "//END//" << endl;
}
