#include <iomanip>
#include <direct.h>
#include"Struct2D.h"
#include "myNastran.h"
#include "interpolation.h"
#include "NasPrinter.h"
#define pi 3.1415926

StructPart::StructPart(const vector<PointSet>& node_2d):node_2D(node_2d)
{
	this->isFixedMass = false;
	// this->node_2D = node_2D;
	this->name = "testpart";
	numi = node_2D[0].X.n_rows;
	numj = node_2D[0].X.n_cols;
	numk = node_2D.size();
	node_f.clear();
	node_f_initial.clear();
	res.clear();
	MaxElemStress = pair<int, double>(0, 4e8);
	MaxNodeDisp = pair<int, double>(0, 10.0);
	MaxTwisting = 10.0;
}

int StructPart::calcStructPart(mat Origin /*= mat(0, 0)*/, mat Rotation/* = mat(0, 0)*/)
{
	int state = checkSite();
	if (state < 0)		return -1;
	calcNode();
	clearUniqueNode();
	calcElem();
	AddOriginRotation(Origin, Rotation);

	return 0;
}

int StructPart::checkSite()
{
	if (siteX_2D.n_rows > 1 && siteX_2D.n_cols > 1)
	{
		cout << "MakeStruct2D_byNum错误！！ X必须是一维向量！" << endl;
		return -1;
	}
	if (siteZ_2D.n_rows > 1 && siteZ_2D.n_cols > 1)
	{
		cout << "MakeStruct2D_byNum错误！！ Z必须是一维向量！" << endl;
		return -1;
	}
	return 0;
}

void StructPart::calcNode()
{
	node = zeros(numi * numj * numk, 5);
	//节点按列排序
	int id = 0;
	for (int k = 0; k < numk; k++)
	{
		for (int j = 0; j < numj; j++)
		{
			for (int i = 0; i < numi; i++)
			{
				node(id, 0) = node_2D[k].X(i, j);
				node(id, 1) = node_2D[k].Y(i, j);
				node(id, 2) = node_2D[k].Z(i, j);
				//node(id, 3) = node_2D[k].ifuse(i, j);//记录节点是否使用
				node(id, 4) = id;//记录节点编号，重复的节点编号会和重复
				id++;
			}
		}
	}
}

void StructPart::clearUniqueNode()
{
	//针对翼型 通常只用 检查j前后方向面多少个点是重复的
	int i;
	for (int k = 0; k < numk; k++)
	{
		for (int j = 0; j < numj; j++)
		{
			i = 0;
			{
				int id = ijk2ID(i, j, k, numi, numj, numk);
				node(id, 4) = ijk2ID(i, j, 0, numi, numj, numk);
			}
		}
		for (int j = 0; j < numj; j++)
		{
			i = numi - 1;
			{
				int id = ijk2ID(i, j, k, numi, numj, numk);
				node(id, 4) = ijk2ID(i, j, 0, numi, numj, numk);
			}
		}
	}
}

void StructPart::calcElem()
{
	//3、整理输入的隔框位置矩阵
	siteX_2D = unique(siteX_2D);
	siteZ_2D = unique(siteZ_2D);
	//数据修正 清除不正确的数据点 X方向
	//没写 默认输入都是正确的吧先
	//4、输出单元
	int PID = 0;
	//气动单元
	int PIDnum = (numj - 1) / 5;
	int aeroElemNum = 2 * (numi - 1) * (numj - 1);
	elem_aero = zeros(aeroElemNum, 6);
	//下表面
	int elemID = 0;
	for (int i = 0; i < numi - 1; i++)
	{
		for (int j = 0; j < numj - 1; j++)
		{
			int k = 0;
			elem_aero(elemID, 0) = getPointID(i, j, k);			node((size_t)elem_aero(elemID, 0), 3) = 1;
			elem_aero(elemID, 1) = getPointID(i + 1, j, k);		node((size_t)elem_aero(elemID, 1), 3) = 1;
			elem_aero(elemID, 2) = getPointID(i + 1, j + 1, k);	node((size_t)elem_aero(elemID, 2), 3) = 1;
			elem_aero(elemID, 3) = getPointID(i, j + 1, k);		node((size_t)elem_aero(elemID, 3), 3) = 1;
			if (PIDnum == 0)
			{
				elem_aero(elemID, 4) = 0;
			}
			else
			{
				elem_aero(elemID, 4) = j / PIDnum;//第五列代表材料属性编号PID
			}

			mat area =
				cross(node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2)); //在四边形四个点弯曲程度不高时可以使用此公式
			double S = Point(area).Norm2() / 2.0;
			elem_aero(elemID, 5) = S;//第六列代表单元面积
			elemID++;
		}
	}
	//上表面
	//PID++;

	for (int i = 0; i < numi - 1; i++)
	{
		for (int j = 0; j < numj - 1; j++)
		{
			int k = numk - 1;
			elem_aero((size_t)elemID, 0) = getPointID(i, j, k);			node((size_t)elem_aero(elemID, 0), 3) = 1;
			elem_aero((size_t)elemID, 1) = getPointID(i + 1, j, k);		node((size_t)elem_aero(elemID, 1), 3) = 1;
			elem_aero((size_t)elemID, 2) = getPointID(i + 1, j + 1, k);	node((size_t)elem_aero(elemID, 2), 3) = 1;
			elem_aero((size_t)elemID, 3) = getPointID(i, j + 1, k);		node((size_t)elem_aero(elemID, 3), 3) = 1;
			if (PIDnum == 0)
			{
				elem_aero(elemID, 4) = 0;
			}
			else
			{
				elem_aero(elemID, 4) = j / PIDnum;//第五列代表材料属性编号PID
			}
			mat area = cross(node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2));//在四边形四个点弯曲程度不高时可以使用此公式
			double S = Point(area).Norm2() / 2.0;
			elem_aero(elemID, 5) = S;//第六列代表单元面积
			elemID++;
		}
	}
	//
	map<int, PSHELL> ps_list = p.getPSHELLlist();
	// 
	//垂直于轴向分布的结构单元 siteX_2D表征其分布位置x坐标(翼肋)
	int strcXElemNum = (numi - 1) * siteX_2D.n_elem * (numk - 1);
	elem_strcX = zeros(strcXElemNum, 6);
	elemID = 0;
	PID = 6;
	if (isFixedMass)//如果isFixedMass为是，则更新材料属性
	{
		ps_list[PID + 1].T = ps_list[PID + 1].T * (double)20 / (double)siteX_2D.n_elem;
	}
	//
	const int rib_end_id = 99;
	PSHELL p_rib_end = p.getPSHELL(PID+1);
	p_rib_end.T *= x_T_ratio;
	p_rib_end.PID = rib_end_id + 1;
	p.add_PSHRLL(move(p_rib_end));

	
	//
	for (size_t jid = 0; jid < siteX_2D.n_elem; jid++)//每一根结构梁
	{
		PID = (jid == siteX_2D.n_elem - 2) ? rib_end_id : 6;
		int j = (int)siteX_2D(jid);
		if (j == numj)
			j = numj - 1;
		if (0 <= j && j < numj)
		{
			for (int i = 0; i < numi - 1; i++)
			{
				for (int k = 0; k < numk - 1; k++)
				{
					elem_strcX((size_t)elemID, 0) = getPointID(i, j, k);			node((size_t)elem_strcX(elemID, 0), 3) = 1;
					elem_strcX((size_t)elemID, 1) = getPointID(i + 1, j, k);		node((size_t)elem_strcX(elemID, 1), 3) = 1;
					elem_strcX((size_t)elemID, 2) = getPointID(i + 1, j, k + 1);	node((size_t)elem_strcX(elemID, 2), 3) = 1;
					elem_strcX((size_t)elemID, 3) = getPointID(i, j, k + 1);		node((size_t)elem_strcX(elemID, 3), 3) = 1;
					elem_strcX((size_t)elemID, 4) = PID;//第五列代表材料属性编号
					mat area =
						cross(node.row(((size_t)elem_strcX(elemID, 0))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 1))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 1))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 2))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 2))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 3))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 3))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 0))).cols(0, 2)); //在四边形四个点弯曲程度不高时可以使用此公式
					double S = Point(area).Norm2() / 2.0;
					elem_strcX(elemID, 5) = S;//第六列代表单元面积
					elemID++;
				}
			}
		}
	}
	//垂直于横向分布的结构单元 siteZ_2D表征其分布位置z坐标(翼梁)
	int strcZElemNum = siteZ_2D.n_elem * (numj - 1) * (numk - 1);
	elem_strcZ = zeros(strcZElemNum, 6);
	elemID = 0;
	PID = 7;
	if (isFixedMass)//如果isFixedMass为是，则更新材料属性
	{
		ps_list[PID + 1].T = ps_list[PID + 1].T * (double)2 / (double)siteZ_2D.n_elem;
	}
	const int spar_end_id = 100;
	PSHELL p_spar_end = p.getPSHELL(PID+1);
	p_spar_end.T *= z_T_ratio;
	p_spar_end.PID = spar_end_id + 1;
	p.add_PSHRLL(move(p_spar_end));
	for (size_t iid = 0; iid < siteZ_2D.n_elem; iid++)//每一根结构梁
	{
		PID = (iid == siteZ_2D.n_elem - 1) ? spar_end_id : 7;
		int i = (int)siteZ_2D(iid);
		if (0 <= i && i <= numi)
		{
			for (int j = 0; j < numj - 1; j++)
			{
				for (int k = 0; k < numk - 1; k++)
				{
					elem_strcZ((size_t)elemID, 0) = getPointID(i, j, k);			node((size_t)elem_strcZ(elemID, 0), 3) = 1;
					elem_strcZ((size_t)elemID, 1) = getPointID(i, j + 1, k);		node((size_t)elem_strcZ(elemID, 1), 3) = 1;
					elem_strcZ((size_t)elemID, 2) = getPointID(i, j + 1, k + 1);	node((size_t)elem_strcZ(elemID, 2), 3) = 1;
					elem_strcZ((size_t)elemID, 3) = getPointID(i, j, k + 1);		node((size_t)elem_strcZ(elemID, 3), 3) = 1;
					elem_strcZ((size_t)elemID, 4) = PID;//第五列代表材料属性编号
					mat area =
						cross(node.row(((size_t)elem_strcZ(elemID, 0))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 1))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 1))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 2))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 2))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 3))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 3))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 0))).cols(0, 2)); //在四边形四个点弯曲程度不高时可以使用此公式
					double S = Point(area).Norm2() / 2.0;
					elem_strcZ(elemID, 5) = S;//第六列代表单元面积
					elemID++;
				}
			}
		}
	}
	if (isFixedMass)//如果isFixedMass为是，则更新材料属性
	{
		p.SetPSHELLList(ps_list);
	}
}

void StructPart::AddOriginRotation(mat Origin, mat Rotation)
{
	//----------4、点的旋转和平移----------
	if (Rotation.n_elem == 3)
	{
		double RotationX = Rotation(0) / 180.0 * pi;
		double RotationY = Rotation(1) / 180.0 * pi;
		double RotationZ = Rotation(2) / 180.0 * pi;

		//还原角度 通过Rotation Matrix
		mat Lrotation = { { cos(RotationZ) * cos(RotationY) , sin(RotationZ) , -cos(RotationZ) * sin(RotationY) },
						 {-sin(RotationZ) * cos(RotationY) * cos(RotationX) + sin(RotationY) * sin(RotationX), cos(RotationZ) * cos(RotationX)	, sin(RotationZ) * sin(RotationY) * cos(RotationX) + cos(RotationY) * sin(RotationX)},
						 { sin(RotationZ) * cos(RotationY) * sin(RotationX) + sin(RotationY) * cos(RotationX),-cos(RotationZ) * sin(RotationX)	,-sin(RotationZ) * sin(RotationY) * sin(RotationX) + cos(RotationY) * cos(RotationX)} };

		node.cols(0, 2) = node.cols(0, 2) * Lrotation;
	}
	if (Origin.n_elem == 3)
	{
		double OriginX = Origin(0);
		double OriginY = Origin(1);
		double OriginZ = Origin(2);
		node.col(0) = node.col(0) + OriginX;
		node.col(1) = node.col(1) + OriginY;
		node.col(2) = node.col(2) + OriginZ;
	}
}

inline int StructPart::getPointID(int i, int j, int k)
{
	int id = ijk2ID(i, j, k, numi, numj, numk);
	return (int)node(id, 4);
}

void StructPart::SaveAsTecplot(ofstream& ofs)
{
	mat elem = join_cols(elem_strcZ.cols(0, 3), elem_strcX.cols(0, 3));
	elem = join_cols(elem, elem_aero.cols(0, 3));

	toTecplot(node, elem + 1, ofs, name);
}

void StructPart::printAnalysisRes(string respath, bool isERROR/* = false*/)
{
	//打印数据
	ofstream ofs;
	ofs.open(respath, ios::trunc);
	if (isERROR == false)
	{
		ofs << "CLsurf     = " << res[1].CLsurf << endl;
		ofs << "CDsurf     = " << res[1].CDsurf << endl;
		ofs << "CL/CD      = " << res[1].CLsurf / res[1].CDsurf << endl;
		ofs << "MaxStress  = " << MaxElemStress.second << endl;
		ofs << "MaxStrain  = " << MaxElemStrain.second << endl;
		ofs << "MaxDisp    = " << MaxNodeDisp.second.getY() << endl;
		ofs << "MaxTwisting= " << MaxTwisting << endl;
	}
	else
	{
		ofs << "CLsurf     = " << 0.0 << endl;
		ofs << "CDsurf     = " << 1.0 << endl;
		ofs << "CL/CD      = " << 0.0 << endl;
		ofs << "MaxStress  = " << 4e8 << endl;
		ofs << "MaxStrain  = " << 10.0 << endl;
		ofs << "MaxDisp    = " << 10.0 << endl;
		ofs << "MaxTwisting= " << 10.0 << endl;
	}
	ofs.close();
}

int StructPart::calcAeroForce(string aeroPath)
{
	//气动输入信息定义
	AeroInfo ai;
	ai.CalcBlockID = { 1,2 };
	ai.CalcBlockArea = { 2,2 };
	ai.CalcBlockType = { 1,1 };
	ai.CalcBlockOrient = { -1,1 };
	ai.Ma = { 10 };
	ai.Alpha = { 30 };
	ai.Beta = { 0 };
	ai.Hm = { 20000 };
	ai.sRef = 240;
	ai.cRef = 3;
	ai.meshPath = aeroPath + "/aerocalc_mesh.bdf";
	ai.outputPath = aeroPath;
	ai.outputName = "struct2d_res.txt";
	printAeroCalcMesh(ai.meshPath);
	//气动计算
	AeroCalc ac;
	ac.setAeroInfo(ai);
	ac.printCommand(aeroPath + "/aero_input.txt");
	ac.CalcAero();

	if (node_f.size() == 0)
	{
		ac.getAeroForce(node_f_initial);
		node_f = node_f_initial;
	}
	else
	{
		map<int, Point> nf;
		nf = node_f_initial;//暂存上一步的气动力
		int state = ac.getAeroForce(node_f_initial);//更新气动力
		subduction_node_f(nf);
	}

	return 0;
}

int StructPart::calcAeroForce_AVL(string aeroPath, string exepath)
{
	
	mat aero_node = (disp_node.n_rows == 0) ? node : disp_node;
	//先创建截面
	AVLInfo ai;
	ai.sectionInfo.resize(numj);
	// ofstream ofs;
    // ofs.open("log_max_y.txt", ios::app);
	// ofs << "---------------\n";
	for (int j = 0; j < numj; j++)
	{
		AVLInfo::AVLSectionInfo &section = ai.sectionInfo[j];
		const int firstNodeID = getPointID(0, j, 0);//翼型顶点的坐标
		const int endNodeID = getPointID(numi - 1, j, 0);
		Point firstNode = Point(aero_node(firstNodeID, 0), aero_node(firstNodeID, 1), aero_node(firstNodeID, 2));
		Point endNode = Point(aero_node(endNodeID, 0), aero_node(endNodeID, 1), aero_node(endNodeID, 2));
		section.Orgin = Point(aero_node(firstNodeID, 0), aero_node(firstNodeID, 2), aero_node(firstNodeID, 1));
		section.Scale = (endNode - firstNode).getX();
		section.Angle = 0;
		section.sectionName = "S" + to_string(j);
		//
		section.pointList.clear();
		double maxDanweiY = 0;
		//下表面节点
		int k = 0;
		for (int i = numi - 1; i > 0; i--)
		{
			int nodeID = getPointID(i, j, k);
			Point temp = Point(aero_node(nodeID, 0), aero_node(nodeID, 1), aero_node(nodeID, 2)) - firstNode;

			maxDanweiY = (maxDanweiY < abs(temp.getY() / section.Scale)) ? abs(temp.getY() / section.Scale) : maxDanweiY;
			// section.pointList.push_back(temp / section.Scale);
			section.pointList.push_back(temp);
		}
		//上表面节点
		k = numk - 1;
		for (int i = 0; i < numi; i++)//i末值取1是因为避免重复书写（0，0）点
		{
			int nodeID = getPointID(i, j, k);
			Point temp = Point(aero_node(nodeID, 0), aero_node(nodeID, 1), aero_node(nodeID, 2)) - firstNode;

			maxDanweiY = (maxDanweiY < abs(temp.getY() / section.Scale)) ? abs(temp.getY() / section.Scale) : maxDanweiY;
			// section.pointList.push_back(temp / section.Scale);
			section.pointList.push_back(temp);
        }
        // ofs << maxDanweiY << "\n";       
		if (maxDanweiY > 0.3)
		{
			return -1;
		}
    }
    // ofs.close();
	//输入参数设置
	
	ai.Nchordwise = 20;
	ai.Nspanwise = numj;// numj;//121;//
	ai.Sref = 231.75;
	ai.Cref = 15;
	ai.Angle = 4;//4是文章里算例用的；4.5是气动验证用的
	ai.Alpha = 0.5;//;
	ai.Ma = 0.76;
	ai._q = 5500;//23895//13895
	//ai.Angle = 0;//4
	ai.rho = 10000;
	ai.Scale = Point(1, 1, 1);
	ai.Orgin = Point(0, 0, 0);

	AVLServer as;
	as.setAVLInfo(ai);
	as.setAVLexePath(exepath);
	as.printCommand(aeroPath);
	as.CalcAero();
	int state = as.getAeroForce(res);
	int num = 1;
	while (state < 0)
	{
		num++;
		ai.Nspanwise++;
		ai.sectionInfo.erase(ai.sectionInfo.end() - 2);
		// vector<AVLInfo::AVLSectionInfo> new_section_info(ai.sectionInfo.size() / 2);
		// for (size_t n = 0; n < ai.sectionInfo.size() / 2; n++)
		// {
		// 	new_section_info[n] = ai.sectionInfo[n * 2];
		// }
		// if((ai.sectionInfo.size() - 1) / 2 * 2 != ai.sectionInfo.size() - 1)
		// 	new_section_info.push_back(ai.sectionInfo[ai.sectionInfo.size() - 1]);
		as.setAVLInfo(ai);
		as.printCommand(aeroPath);
		as.CalcAero();
		state = as.getAeroForce(res);
		if (num > 20)//更正二十步之后还无法计算
		{
			return -2;//气动力无法计算
		}
	}
	vector<XFoilInfo> xfoilList(numj);//XFoil求解参数设置
	vector<vector<Field<double>>> Xfoil_res(numj);
	// for (size_t i = 0; i < (size_t)1; i++)//(size_t)numj; i++)
	// {
	// 	xfoilList[i].foilName = ai.sectionInfo[i].sectionName+"_sfoil";
	// 	xfoilList[i].foilList = ai.sectionInfo[i].pointList;
	// 	xfoilList[i].alpha = {ai.Alpha, ai.Alpha+3, 1};
	// 	xfoilList[i].Ma = ai.Ma;
	// 	xfoilList[i].VISC = 500000;//随便给了一个雷诺数
	// 	xfoilList[i].outputNum = 200;//xfoilList[i].foilList.size() < 300 ? xfoilList[i].foilList.size() : 299;
	// 	xfoilList[i].workPath = aeroPath;
	// 	XFoilSolver xfoilSolver(xfoilList[i]);
	// 	xfoilSolver.ExePath = exepath + "/xfoil.exe";
	// 	if (xfoilSolver.solve() < 0)
	// 	{
	// 		return -3;
	// 	}
	// 	Xfoil_res[i] = xfoilSolver.getRes();
	// }

	calcElemPress_AVL_refine(node_f, Xfoil_res);
	//calcElemPress_AVL(node_f);

	return 0;
}

double StructPart::calcMass()
{
	if (elem_aero.n_cols < 6 || elem_strcX.n_cols < 6 || elem_strcZ.n_cols < 6)//如果没有第六列
	{
		cout << "当前没有计算单元面积，无法计算质量 返回值：-1.0" << endl;
		return -1.0;
	}
	if (!p.isNotEmpty())
	{
		cout << "当前没有设置材料属性，无法计算质量 返回值：-1.0" << endl;
		return -1.0;
	}
	int PID = 0;
	int MID = 0;
	double MassTotal = 0;
	for (size_t i = 0; i < elem_aero.n_rows; i++)
	{
		PID = (int)elem_aero(i, 4) + 1;
		MID = p.getPSHELL(PID).MID;
		MassTotal += elem_aero(i, 5) * p.getPSHELL(PID).T * p.getMAT1(MID).RHO;
	}
	for (size_t i = 0; i < elem_strcX.n_rows; i++)
	{
		PID = (int)elem_strcX(i, 4) + 1;
		MID = p.getPSHELL(PID).MID;
		MassTotal += elem_strcX(i, 5) * p.getPSHELL(PID).T * p.getMAT1(MID).RHO;
	}
	for (size_t i = 0; i < elem_strcZ.n_rows; i++)
	{
		PID = (int)elem_strcZ(i, 4) + 1;
		MID = p.getPSHELL(PID).MID;
		MassTotal += elem_strcZ(i, 5) * p.getPSHELL(PID).T * p.getMAT1(MID).RHO;
	}
	return MassTotal;
}

int StructPart::calcElemPress_AVL(map<int, Point>& nf)
{
	double Xle = res[1].stripRes[1].Orgin.getX();
	double Chord = res[1].stripRes[1].Chord;

	vector<double> ChordList(res[1].stripRes[1].Chordwise);
	for (int i = 0; i < res[1].stripRes[1].Chordwise; i++)
	{
		int elemID = res[1].stripRes.begin()->second.FirstVortex + i;
		double siteX = res[1].stripRes[1].elemRes[elemID].site.getX();
		ChordList[i] = (siteX - Xle) / Chord;
	}
	vector<double> StripList(res[1].Spanwise);
	vector<double> StripXleList(res[1].Spanwise);
	vector<double> StripZleList(res[1].Spanwise);
	vector<double> StripChordList(res[1].Spanwise);
	for (int i = 0; i < res[1].Spanwise; i++)
	{
		double siteY = res[1].stripRes[res[1].FirstStrip + i].elemRes.begin()->second.site.getY();
		StripList[i] = siteY;
		StripXleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getX();
		StripZleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getZ();
		StripChordList[i] = res[1].stripRes[res[1].FirstStrip + i].Chord;
	}
	vector<double> CdList(res[1].Spanwise * res[1].Chordwise);
	vector<double> ClList(res[1].Spanwise * res[1].Chordwise);
	for (int i = 0; i < res[1].Spanwise; i++)
	{
		for (int j = 0; j < res[1].Chordwise; j++)
		{
			int stripID = res[1].FirstStrip + i;
			int elemID = res[1].stripRes[stripID].FirstVortex + j;
			CdList[res[1].Chordwise * i + j] = res[1].stripRes[stripID].elemRes[elemID].dCd();
			ClList[res[1].Chordwise * i + j] = res[1].stripRes[stripID].elemRes[elemID].dCl();
		}
	}

	double areas = 0;
	double Fx = 0; double Fy = 0;
	mat info = zeros(elem_aero.n_rows / 2, 3);//用于存储单元：面积、fx、fy
	//注意 AVL坐标系和CST使用的坐标系y、z轴是反的
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//只要上表面
	{
		//计算面积
		mat node0 = (node.row((size_t)elem_aero(i, 0)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 0)).cols(0, 2)) / 2.0;
		node0(1) = Interpolation1(StripList, StripZleList, node0(2));
		mat node1 = (node.row((size_t)elem_aero(i, 1)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 1)).cols(0, 2)) / 2.0;
		node1(1) = Interpolation1(StripList, StripZleList, node1(2));
		mat node2 = (node.row((size_t)elem_aero(i, 2)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 2)).cols(0, 2)) / 2.0;
		node2(1) = Interpolation1(StripList, StripZleList, node2(2));
		mat node3 = (node.row((size_t)elem_aero(i, 3)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 3)).cols(0, 2)) / 2.0;
		node3(1) = Interpolation1(StripList, StripZleList, node3(2));
		mat area = cross(node0, node1) + cross(node1, node2) + cross(node2, node3) + cross(node3, node0);//在四边形四个点弯曲程度不高时可以使用此公式

		double S = Point(area).Norm2() / 2.0;

		//计算参考点压力系数
		mat nodeMid = (node0 + node3) * 3.0 / 8.0 + (node1 + node2) * 1.0 / 8.0;//取单元靠前1/4位置处的点为压力参考点
		double x = nodeMid(0);
		double z = nodeMid(2);
		double px = Interpolation2(ChordList, StripList, CdList, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z), z);
		double py = Interpolation2(ChordList, StripList, ClList, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z), z);
		areas += S;
		Fx += S * px;
		Fy += S * py;
		info(i, 0) = S;
		info(i, 1) = S * px;//阻力
		info(i, 2) = S * py;//升力
	}
	double ratioFx = res[1].CDsurf * res[1].Surface_area / Fx;
	double ratioFy = res[1].CLsurf * res[1].Surface_area / Fy;
	info.col(1) *= ratioFx;
	info.col(2) *= ratioFy;
	nf.clear();
	double _q = res[1].ai._q;
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//只要上表面
	{
		nf[1 + (size_t)elem_aero(i, 0)] += Field<double>(info(i, 1), info(i, 2), 0) * 3.0 / 16.0 * _q * 2.0;/// 8.0
		nf[1 + (size_t)elem_aero(i, 1)] += Field<double>(info(i, 1), info(i, 2), 0) * 3.0 / 16.0 * _q * 2.0;/// 8.0
		nf[1 + (size_t)elem_aero(i, 2)] += Field<double>(info(i, 1), info(i, 2), 0) * 1.0 / 16.0 * _q * 2.0;/// 8.0
		nf[1 + (size_t)elem_aero(i, 3)] += Field<double>(info(i, 1), info(i, 2), 0) * 1.0 / 16.0 * _q * 2.0;/// 8.0
		nf[1 + (size_t)elem_aero(i + elem_aero.n_rows / 2, 0)] += Field<double>(info(i, 1), info(i, 2), 0) * 3.0 / 16.0 * _q * -1.0;
		nf[1 + (size_t)elem_aero(i + elem_aero.n_rows / 2, 1)] += Field<double>(info(i, 1), info(i, 2), 0) * 3.0 / 16.0 * _q * -1.0;
		nf[1 + (size_t)elem_aero(i + elem_aero.n_rows / 2, 2)] += Field<double>(info(i, 1), info(i, 2), 0) * 1.0 / 16.0 * _q * -1.0;
		nf[1 + (size_t)elem_aero(i + elem_aero.n_rows / 2, 3)] += Field<double>(info(i, 1), info(i, 2), 0) * 1.0 / 16.0 * _q * -1.0;
		//nf[1 + elem_aero(i, 0)].setP(1);
		//nf[1 + elem_aero(i, 1)].setP(1);
		//nf[1 + elem_aero(i, 2)].setP(1);
		//nf[1 + elem_aero(i, 3)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 0)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 1)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 2)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 3)].setP(1);
	}
	return 0;
}

int StructPart::calcElemPress_AVL_refine(map<int, Point>& nf, const vector<vector<Field<double>>> &xfoilRes)
{
	//上表面机翼形状，横坐标取值在0-1之间 P值为cp
	vector<Field<double>> foilUpp = getPlistFromfile("E:\\Matlab-Xfoil\\DWPupp.txt");//这是两组归一化的节点坐标参数和对应的压力系数
	vector<Field<double>> foilLow = getPlistFromfile("E:\\Matlab-Xfoil\\DWPlow.txt");//这是两组归一化的节点坐标参数和对应的压力系数
	//读取上下表面坐标和真实CP分布值
	auto clac = [](const vector<Field<double>> &foil, vector<double> &xsite, vector<double> &cp)
	{
		xsite.resize(foil.size());
		cp.resize(foil.size());
		for (size_t i = 0; i < foil.size(); i++)
		{
			xsite[i] = foil[i].getX() / foil[foil.size()-1].getX();//获得归一化的x坐标列表
			//ysite.push_back(foil[i].getY());
			cp[i] = foil[i].getP();//获得压力系数列表
		}
	};

	vector<double> xsiteupp, /*ysiteupp,*/ cpupp;
	clac(foilUpp, xsiteupp, cpupp);//读取上下表面坐标和真实CP分布值
	vector<double> xsitelow, /*ysitelow,*/ cplow;
	clac(foilLow, xsitelow, cplow);//读取上下表面坐标和真实CP分布值
	//计算xfoil计算出的总升阻力系数
	double cl_S = 0, cd_S = 0;
	auto calcS = [&cl_S, &cd_S](const vector<Field<double>> &foil) mutable
	{
		for (size_t i = 0; i < foil.size() - 1; i++)
		{
			double dx = foil[i + 1].getX() - foil[i].getX();
			double dy = foil[i + 1].getY() - foil[i].getY();
			double cp = foil[i + 1].getP() * 0.75 + foil[i].getP() * 0.25;
			cl_S += cp * dx;
			cd_S += cp * dy;
		}
	};
	calcS(foilUpp);//计算xfoil计算出的总升阻力系数
	calcS(foilLow);//计算xfoil计算出的总升阻力系数
	//AVL计算结果
	vector<double> StripList(res[1].Spanwise);
	vector<double> StripXleList(res[1].Spanwise);
	vector<double> StripZleList(res[1].Spanwise);
	vector<double> StripChordList(res[1].Spanwise);

	vector<double> StripClList(res[1].Spanwise);
	vector<double> StripArea(res[1].Spanwise);
	vector<double> StripRatioList(res[1].Spanwise);
	for (int i = 0; i < res[1].Spanwise; i++)//计算每一条截面
	{
		double siteY = res[1].stripRes[res[1].FirstStrip + i].elemRes.begin()->second.site.getY();
		StripList[i] = siteY;//每一条截面的起始位置
		StripXleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getX();
		StripZleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getZ();
		StripChordList[i] = res[1].stripRes[res[1].FirstStrip + i].Chord;//每一条截面的宽度

		StripClList[i] = res[1].stripRes[res[1].FirstStrip + i].cl;//每一个截面的cl
		StripArea[i] = res[1].stripRes[res[1].FirstStrip + 1].StripArea;//每一个截面的参考面积
		StripRatioList[i] = StripClList[i] / cl_S;//每一条截面升力系数与XFoil算出的cl的比值系数
	}

	double areas = 0;
	double Fx = 0;
	double Fy = 0;
	mat infoUpp = zeros(elem_aero.n_rows / 2, 3);//用于存储单元：面积、fx、fy
	mat infoLow = zeros(elem_aero.n_rows / 2, 3);//用于存储单元：面积、fx、fy
	enum class Side{
		upp,
		low,
	};
	using vec_d = const vector<double> &;
	//注意 AVL坐标系和CST使用的坐标系y、z轴是反的
	auto calcElemForce = [&StripList,&StripXleList,&StripChordList,&StripRatioList,this](Side side,vec_d xsite,vec_d cp,mat&info,double&areas,double&Fx,double&Fy) {	 
		for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//只要上表面 每个单元
		{
			const int id = (side == Side::upp) ? i : (i + elem_aero.n_rows / 2);
			const int sign = (side == Side::upp) ? 1 : -1; //下表面取y的负号
			//计算面积
			const mat node0 = node.row((size_t)elem_aero(id, 0)).cols(0, 2);//获得单元的四个节点坐标
			const mat node1 = node.row((size_t)elem_aero(id, 1)).cols(0, 2);
			const mat node2 = node.row((size_t)elem_aero(id, 2)).cols(0, 2);
			const mat node3 = node.row((size_t)elem_aero(id, 3)).cols(0, 2);
			const double dz = (node3(2) - node0(2)) / 2.0 + (node2(2) - node1(2)) / 2.0;
			const double dx = (node1(0) - node0(0)) / 2.0 + (node2(0) - node3(0)) / 2.0;
			const double dy = ((node1(1) - node0(1)) / 2.0 + (node2(1) - node3(1)) / 2.0) * sign;

			//计算参考点压力系数
			const mat nodeMid = (node0 + node3) * 3.0 / 8.0 + (node1 + node2) * 1.0 / 8.0;//取单元靠前1/4位置处的点为压力参考点
			const double x = nodeMid(0);
			const double z = nodeMid(2);
			const double cp_temp = Interpolation1(xsite, cp, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z));
			const double ratio_temp = Interpolation1(StripList, StripRatioList, z);
			areas += dz * dx;
			Fx += cp_temp * ratio_temp * dz * dy;
			Fy += cp_temp * ratio_temp * dz * dx;
			info(i, 0) = ratio_temp * dz * dx;
			info(i, 1) = cp_temp * ratio_temp * dz * dy;//阻力
			info(i, 2) = cp_temp * ratio_temp * dz * dx;//升力
		}
	};
	calcElemForce(Side::upp, xsiteupp, cpupp, infoUpp, areas, Fx, Fy);
	calcElemForce(Side::low, xsitelow, cplow, infoLow, areas, Fx, Fy);
	

	double ratioAe = res[1].Surface_area / areas;
	double ratioFx = res[1].CDsurf * res[1].Surface_area / Fx;
	double ratioFy = res[1].CLsurf * res[1].Surface_area / Fy;
	//info.col(1) *= ratioFy;
	//info.col(2) *= ratioFy;
	nf.clear();
	double _q = res[1].ai._q;
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//只要上表面
	{
		//上表面
		nf[1 + (int)elem_aero(i, 0)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 3.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 1)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 3.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 2)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 1.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 3)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 1.0 / 8.0 * _q;/// 8.0
		//下表面
		nf[1 + (int)elem_aero(i + elem_aero.n_rows / 2, 0)] += Field<double>(infoLow(i, 1), infoLow(i, 2), 0) * 3.0 / 8.0 * _q;
		nf[1 + (int)elem_aero(i + elem_aero.n_rows / 2, 1)] += Field<double>(infoLow(i, 1), infoLow(i, 2), 0) * 3.0 / 8.0 * _q;
		nf[1 + (int)elem_aero(i + elem_aero.n_rows / 2, 2)] += Field<double>(infoLow(i, 1), infoLow(i, 2), 0) * 1.0 / 8.0 * _q;
		nf[1 + (int)elem_aero(i + elem_aero.n_rows / 2, 3)] += Field<double>(infoLow(i, 1), infoLow(i, 2), 0) * 1.0 / 8.0 * _q;
		//nf[1 + elem_aero(i, 0)].setP(1);
		//nf[1 + elem_aero(i, 1)].setP(1);
		//nf[1 + elem_aero(i, 2)].setP(1);
		//nf[1 + elem_aero(i, 3)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 0)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 1)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 2)].setP(1);
		//nf[1 + elem_aero(i + elem_aero.n_rows / 2, 3)].setP(1);
	}
	return 0;
}

void StructPart::SaveAsNastran1(string fileName, int SOL)
{
	ofstream ofs, ofs1, ofs2;
	ofs.open(fileName + "_mesh.bdf", ios::trunc);
	ofs1.open(fileName + "_property.bdf", ios::trunc);
	ofs2.open(fileName + "_force.bdf", ios::out);
	string names = fileName.substr(fileName.find_last_of("/") + 1);
	ofs << fixed << setprecision(8);//设置输出精度为8
	ofs1 << showpoint << setprecision(2);//强制输出小数点
	ofs2 << fixed << setprecision(5);//设置输出精度为5
	int ID = 0;//节点编号
	int EID = 0;//单元编号
	int PID = 0;//截面or壳单元属性编号
	int SID = 1;//载荷or约束工况编号

	//-----------------------------单元节点信息----------------------------------
	//----------------node----------------
	for (size_t i = 0; i < node.n_rows; i++)
	{
		ID++;
		if (node(i, 3) == 1)
		{
			ofs << "GRID*   "
				<< setw(16) << ID
				<< setw(32) << node(i, 0)
				<< setw(16) << node(i, 1)
				<< setw(8) << ID << endl;
			ofs << "*" << setw(7) << ID << setw(16) << node(i, 2) << endl;

			//----------------node force----------------
			map<int, Point>::iterator it = node_f.find(ID);
			if (it != node_f.end() && it->second.Norm1() > 1e-4)
			{
				//return it->second;
				ofs2 << "FORCE*  "
					 << setw(16) << SID
					 << setw(16) << ID
					 << setw(16) << 0
					 << setw(16) << 1.0 << "*" << endl;
				ofs2 << "*       "
					 << setw(16) << it->second.getX()
					 << setw(16) << it->second.getY()
					 << setw(16) << it->second.getZ() << endl;
			}
		}
	}
	int SID_F_all = 1;
	//位移约束添加
	int SID_SPC = 2;
	int SPC1num = numi * 2 - 2;
	ofs2 << "SPC1    " << setw(8) << SID_SPC << setw(8) << "123456";
	for (int n = 0; n < SPC1num; n++)//翼根处蒙皮的约束
	{
		int i, j, k;
		if (n < numi)//shangbiaomian
		{
			i = n; k = 0; j = 0;
		}
		else
		{
			i = n - numi + 1; k = numk - 1; j = 0;
		}
		ofs2 << setw(8) << ijk2ID(i, j, k, numi, numj, numk) + 1;
		if (n % 8 == 5)	ofs2 << endl << "        ";
	}

	for (size_t m = 0; m < siteZ_2D.n_elem; m++)///*翼梁的个数*/
	{
		int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//第n根梁的第一个单元的单元id
		for (int k = 0; k < numk - 1; k++)
		{
			int n = SPC1num + m * (numk - 1) + k;//书接上文循环中的n
			ofs2 << setw(8) << (int)elem_strcZ(firstID + k, 3) + 1;
			if (n % 8 == 5)	ofs2 << endl << "        ";
		}
	}

	//----------------element----------------
	for (size_t i = 0; i < elem_aero.n_rows; i++)//气动壳单元
	{
		//EID++;
		PID = (int)elem_aero(i, 4) + 1;
		ofs << "CQUAD4  "
			<< setw(8) << ++EID
			<< setw(8) << PID
			<< setw(8) << (int)elem_aero(i, 0) + 1
			<< setw(8) << (int)elem_aero(i, 1) + 1
			<< setw(8) << (int)elem_aero(i, 2) + 1
			<< setw(8) << (int)elem_aero(i, 3) + 1 << endl;
	}

	for (size_t i = 0; i < elem_strcX.n_rows; i++)//结构壳单元 翼肋
	{
		//EID++;
		PID = (int)elem_strcX(i, 4) + 1;
		mat nodeID = unique(elem_strcX.row(i).cols(0, 3));
		if (nodeID.n_elem == 3)
		{
			mat equallist = uniqueMat::isPointEqual(trans(elem_strcX.row(i).cols(0, 3)));
			ofs << "CTRIA3  "
				<< setw(8) << ++EID
				<< setw(8) << PID;
			for (size_t id = 0; id < 4; id++)
			{
				if (equallist(id) == id)
				{
					ofs << setw(8) << (int)elem_strcX(i, id) + 1;
				}
			}
			ofs << endl;
		}
		else
		{
			ofs << "CQUAD4  "
				<< setw(8) << ++EID
				<< setw(8) << PID
				<< setw(8) << (int)elem_strcX(i, 0) + 1
				<< setw(8) << (int)elem_strcX(i, 1) + 1
				<< setw(8) << (int)elem_strcX(i, 2) + 1
				<< setw(8) << (int)elem_strcX(i, 3) + 1 << endl;
		}

	}
	for (size_t i = 0; i < elem_strcZ.n_rows; i++)//结构壳单元 翼梁
	{
		PID = (int)elem_strcZ(i, 4) + 1;
		ofs << "CQUAD4  "
			<< setw(8) << ++EID
			<< setw(8) << PID
			<< setw(8) << (int)elem_strcZ(i, 0) + 1
			<< setw(8) << (int)elem_strcZ(i, 1) + 1
			<< setw(8) << (int)elem_strcZ(i, 2) + 1
			<< setw(8) << (int)elem_strcZ(i, 3) + 1 << endl;
	}

	//----------------property----------------
	p.printPSHELL_all(ofs1);
	p.printMAT1_all(ofs1);

	//--------------------------------------------------
	ofstream ofs3;
	switch (SOL)
	{
	case 101:
		ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//执行控制部分
		//ofs3 << "NASTRAN parallel=8" << endl;//对应组合工况的SID
		ofs3 << "SOL 101" << endl;
		//工况控制命令
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SUBCASE 1" << endl;//第一个工况
		ofs3 << "   SUBTITLE=Default" << endl;
		ofs3 << "   LOAD = " << SID_F_all << endl;//对应组合工况的SID
		ofs3 << "   SPC  = " << SID_SPC << endl;//对应组合工况的SID

		ofs3 << "   STRESS(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//输出单元应力 ,PUNCH
		//ofs3 << "   STRAIN(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//
		ofs3 << "   DISPLACEMENT(SORT1,PUNCH,REAL)=ALL" << endl;//
		 //ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;//
		//BEGIN BULK
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		ofs3 << "PARAM   POST     0" << endl;//生成'.XBD'类型可视化文件
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		//ofs3 << "PARAM,INREL,-2" << endl;//定义惯性释放

		ofs3 << "include '" << (names + "_mesh.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "include '" << (names + "_property.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "include '" << (names + "_force.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "ENDDATA" << endl;
		break;

	case 105://屈曲分析
		ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//执行控制部分
		ofs3 << "SOL 105" << endl;
		ofs3 << "TIME 10000" << endl;
		//工况控制命令
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SPC  = " << SID_SPC << endl;//对应组合工况的SID
		ofs3 << "SUBCASE 1" << endl;//第一个工况
		ofs3 << "   LOAD = " << SID_F_all << endl;//对应组合工况的SID
		ofs3 << "SUBCASE 2" << endl;//第二个工况
		ofs3 << "   METHOD = 10" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		ofs3 << "PARAM   POST     0" << endl;//生成'.XBD'类型可视化文件
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		ofs3 << "include '" << (names + "_mesh.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "include '" << (names + "_property.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "include '" << (names + "_force.bdf") << "'" << endl;//头文件内部使用相对路径
		ofs3 << "EIGRL   10                      5" << endl;//模态输出控制关键字 输出前五阶模态
		ofs3 << "ENDDATA" << endl;
		break;

	default:
		cout << "当前nastran求解类型为：" << SOL << "，输入错误，位于void StructPart::SaveAsNastran(string fileName)" << endl;
		break;
	}
	ofs.close();
	ofs1.close();
	ofs2.close();
	ofs3.close();
	cout << "struct-2D网格(*.BDF)导出完成！(SOL = " << SOL << ")" << endl;
}

void StructPart::SaveAsNastran(string fileName, int SOL)
{
	NasPrinter np;

	// ofstream ofs, ofs1, ofs2;
	// ofs.open(fileName + "_mesh.bdf", ios::trunc);
	// ofs1.open(fileName + "_property.bdf", ios::trunc);
	// ofs2.open(fileName + "_force.bdf", ios::out);
	// string names = fileName.substr(fileName.find_last_of("/") + 1);
	// ofs << fixed << setprecision(8);//设置输出精度为8
	// ofs1 << showpoint << setprecision(2);//强制输出小数点
	// ofs2 << fixed << setprecision(5);//设置输出精度为5
	int ID = 0;//节点编号
	int EID = 0;//单元编号
	int PID = 0;//截面or壳单元属性编号
	int SID = 1;//载荷or约束工况编号

	//-----------------------------单元节点信息----------------------------------
	//----------------node----------------
	for (size_t i = 0; i < node.n_rows; i++)
	{
		ID++;
		if (node(i, 3) == 1)
		{
			// ofs << "GRID*   "
			// 	<< setw(16) << ID
			// 	<< setw(32) << node(i, 0)
			// 	<< setw(16) << node(i, 1)
			// 	<< setw(8) << ID << endl;
			// ofs << "*" << setw(7) << ID << setw(16) << node(i, 2) << endl;
			np.addGRID(ID, 0, {node(i, 0), node(i, 1), node(i, 2)});

			//----------------node force----------------
			map<int, Point>::iterator it = node_f.find(ID);
			if (it != node_f.end() && it->second.Norm1() > 1e-4)
			{
				//return it->second;
				// ofs2 << "FORCE*  "
				// 	 << setw(16) << SID
				// 	 << setw(16) << ID
				// 	 << setw(16) << 0
				// 	 << setw(16) << 1.0 << "*" << endl;
				// ofs2 << "*       "
				// 	 << setw(16) << it->second.getX()
				// 	 << setw(16) << it->second.getY()
				// 	 << setw(16) << it->second.getZ() << endl;
				np.addFORCE(SID, ID, 0, 1.0, it->second);
			}
		}
	}
	int SID_F_all = 1;
	//位移约束添加
	int SID_SPC = 2;
	int SPC1num = numi * 2 - 2;
	vector<int> spc_list(0);
	spc_list.reserve(SPC1num + siteZ_2D.n_elem * (numk - 1));
	for (int n = 0; n < SPC1num; n++)//翼根处蒙皮的约束
	{
		int i, j, k;
		if (n < numi)//shangbiaomian
		{
			i = n; k = 0; j = 0;
		}
		else
		{
			i = n - numi + 1; k = numk - 1; j = 0;
		}
		spc_list.push_back(ijk2ID(i, j, k, numi, numj, numk) + 1);
	}

	for (size_t m = 0; m < siteZ_2D.n_elem; m++)//翼梁的个数
	{
		int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//第n根梁的第一个单元的单元id
		for (int k = 0; k < numk - 1; k++)
		{
			int n = SPC1num + m * (numk - 1) + k;//书接上文循环中的n
			spc_list.push_back((int)elem_strcZ(firstID + k, 3) + 1);
		}
	}
	np.addSPC1(SID_SPC, 123456, spc_list);
	// ofs2 << "SPC1    " << setw(8) << SID_SPC << setw(8) << "123456";
	
	// for (int n = 0; n < SPC1num; n++)//翼根处蒙皮的约束
	// {
	// 	int i, j, k;
	// 	if (n < numi)//shangbiaomian
	// 	{
	// 		i = n; k = 0; j = 0;
	// 	}
	// 	else
	// 	{
	// 		i = n - numi + 1; k = numk - 1; j = 0;
	// 	}
	// 	ofs2 << setw(8) << ijk2ID(i, j, k, numi, numj, numk) + 1;
	// 	if (n % 8 == 5)	ofs2 << endl << "        ";
	// }

	// for (size_t m = 0; m < siteZ_2D.n_elem; m++)//翼梁的个数
	// {
	// 	int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//第n根梁的第一个单元的单元id
	// 	for (int k = 0; k < numk - 1; k++)
	// 	{
	// 		int n = SPC1num + m * (numk - 1) + k;//书接上文循环中的n
	// 		ofs2 << setw(8) << (int)elem_strcZ(firstID + k, 3) + 1;
	// 		if (n % 8 == 5)	ofs2 << endl << "        ";
	// 	}
	// }

	//----------------element----------------
	for (size_t i = 0; i < elem_aero.n_rows; i++)//气动壳单元
	{
		//EID++;
		PID = (int)elem_aero(i, 4) + 1;
		np.addCQUAD4(
			++EID,
			PID,
			(int)elem_aero(i, 0) + 1,
			(int)elem_aero(i, 1) + 1,
			(int)elem_aero(i, 2) + 1,
			(int)elem_aero(i, 3) + 1);
		// ofs << "CQUAD4  "
		// 	<< setw(8) << ++EID
		// 	<< setw(8) << PID
		// 	<< setw(8) << (int)elem_aero(i, 0) + 1
		// 	<< setw(8) << (int)elem_aero(i, 1) + 1
		// 	<< setw(8) << (int)elem_aero(i, 2) + 1
		// 	<< setw(8) << (int)elem_aero(i, 3) + 1 << endl;
	}

	for (size_t i = 0; i < elem_strcX.n_rows; i++)//结构壳单元 翼肋
	{
		//EID++;
		PID = (int)elem_strcX(i, 4) + 1;
		mat nodeID = unique(elem_strcX.row(i).cols(0, 3));
		if (nodeID.n_elem == 3)
		{
			mat equallist = uniqueMat::isPointEqual(trans(elem_strcX.row(i).cols(0, 3)));
			np.ssMesh << "CTRIA3  "
				<< setw(8) << ++EID
				<< setw(8) << PID;
			for (size_t id = 0; id < 4; id++)
			{
				if (equallist(id) == id)
				{
					np.ssMesh << setw(8) << (int)elem_strcX(i, id) + 1;
				}
			}
			np.ssMesh << endl;
		}
		else
		{
			np.addCQUAD4(
				++EID,
				PID,
				(int)elem_strcX(i, 0) + 1,
				(int)elem_strcX(i, 1) + 1,
				(int)elem_strcX(i, 2) + 1,
				(int)elem_strcX(i, 3) + 1);
		}

	}
	for (size_t i = 0; i < elem_strcZ.n_rows; i++)//结构壳单元 翼梁
	{
		PID = (int)elem_strcZ(i, 4) + 1;
		np.addCQUAD4(
				++EID,
				PID,
				(int)elem_strcZ(i, 0) + 1,
				(int)elem_strcZ(i, 1) + 1,
				(int)elem_strcZ(i, 2) + 1,
				(int)elem_strcZ(i, 3) + 1);
		// ofs << "CQUAD4  "
		// 	<< setw(8) << ++EID
		// 	<< setw(8) << PID
		// 	<< setw(8) << (int)elem_strcZ(i, 0) + 1
		// 	<< setw(8) << (int)elem_strcZ(i, 1) + 1
		// 	<< setw(8) << (int)elem_strcZ(i, 2) + 1
		// 	<< setw(8) << (int)elem_strcZ(i, 3) + 1 << endl;
	}

    //----------------property----------------
    map<int,PSHELL> p_list = p.getPSHELLlist();
	for (auto it = p_list.begin(); it != p_list.end(); it++)
	{
        (*it).second.MID2 = 1; // 缺了这两行会出问题
        (*it).second.MID3 = 1; // 之前的计算都使用了这两行 应该是添加了额外的刚度
        np.addPSHELL((*it).second);
	}
	for (auto it = p.getMAT1List().begin(); it != p.getMAT1List().end(); it++)
	{
		np.addMAT1((*it).second);
	}
	// p.printPSHELL_all(ofs1);
	// p.printMAT1_all(ofs1);

	//--------------------------------------------------
	ofstream ofs3;
	switch (SOL)
	{
	case 101:
		// ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//执行控制部分
		//ofs3 << "NASTRAN parallel=8" << endl;//对应组合工况的SID
		np.ssHeader << "SOL 101" << endl;
		//工况控制命令
		np.ssHeader << "CEND" << endl;
		np.ssHeader << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		np.ssHeader << "ECHO = NONE" << endl;
		np.ssHeader << "SUBCASE 1" << endl;//第一个工况
		np.ssHeader << "   SUBTITLE=Default" << endl;
		np.ssHeader << "   LOAD = " << SID_F_all << endl;//对应组合工况的SID
		np.ssHeader << "   SPC  = " << SID_SPC << endl;//对应组合工况的SID

		np.ssHeader << "   STRESS(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//输出单元应力 ,PUNCH
		//ofs3 << "   STRAIN(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//
		np.ssHeader << "   DISPLACEMENT(SORT1,PUNCH,REAL)=ALL" << endl;//
		 //ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;//
		//BEGIN BULK
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "$$                                Bulk Data Cards                               $" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "BEGIN BULK" << endl;
		np.ssHeader << "PARAM   POST     0" << endl;//生成'.XBD'类型可视化文件
		np.ssHeader << "PARAM   PRTMAXIM YES" << endl;
		//ofs3 << "PARAM,INREL,-2" << endl;//定义惯性释放

		// np.ssHeader << "include '" << (names + "_mesh.bdf") << "'" << endl;//头文件内部使用相对路径
		// np.ssHeader << "include '" << (names + "_property.bdf") << "'" << endl;//头文件内部使用相对路径
		// np.ssHeader << "include '" << (names + "_force.bdf") << "'" << endl;//头文件内部使用相对路径
		// np.ssHeader << "ENDDATA" << endl;
		break;

	case 105://屈曲分析
		// ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//执行控制部分
		np.ssHeader << "SOL 105" << endl;
		np.ssHeader << "TIME 10000" << endl;
		//工况控制命令
		np.ssHeader << "CEND" << endl;
		np.ssHeader << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		np.ssHeader << "ECHO = NONE" << endl;
		np.ssHeader << "SPC  = " << SID_SPC << endl;//对应组合工况的SID
		np.ssHeader << "SUBCASE 1" << endl;//第一个工况
		np.ssHeader << "   LOAD = " << SID_F_all << endl;//对应组合工况的SID
		np.ssHeader << "SUBCASE 2" << endl;//第二个工况
		np.ssHeader << "   METHOD = 10" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "$$                                Bulk Data Cards                               $" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "BEGIN BULK" << endl;
		np.ssHeader << "PARAM   POST     0" << endl;//生成'.XBD'类型可视化文件
		np.ssHeader << "PARAM   PRTMAXIM YES" << endl;
		// np.ssHeader << "include '" << (names + "_mesh.bdf") << "'" << endl;//头文件内部使用相对路径
		// np.ssHeader << "include '" << (names + "_property.bdf") << "'" << endl;//头文件内部使用相对路径
		// np.ssHeader << "include '" << (names + "_force.bdf") << "'" << endl;//头文件内部使用相对路径
		np.ssHeader << "EIGRL   10                      5" << endl;//模态输出控制关键字 输出前五阶模态
		// np.ssHeader << "ENDDATA" << endl;
		break;

	default:
		cout << "当前nastran求解类型为：" << SOL << "，输入错误，位于void StructPart::SaveAsNastran(string fileName)" << endl;
		break;
	}
	// ofs.close();
	// ofs1.close();
	// ofs2.close();
	// ofs3.close();
	string::size_type iPos = (fileName.find_last_of('\\') + 1) == 0 ? fileName.find_last_of('/') + 1 : fileName.find_last_of('\\') + 1;
	string datapath = fileName.substr(0, iPos);//获取文件路径
    string names = fileName.substr(iPos);
    iPos = (names.find_last_of('\\') + 1) == 0 ? names.find_last_of('/') + 1
                                               : names.find_last_of('\\') + 1;
    datapath += names.substr(0, iPos);
    names = names.substr(iPos);
	np.PrintBDF(datapath, names, NasPrinter::onefile);
	cout << "struct-2D网格(*.BDF)导出完成！(SOL = " << SOL << ")" << endl;
}

int StructPart::printAeroCalcMesh(string filepath)
{
	ofstream ofs;
	ofs.open(filepath, ios::trunc);
	ofs << fixed << setprecision(8);//设置输出精度为8
	int ID = 0;//节点编号
	int EID = 0;//单元编号
	int PID = 0;//截面or壳单元属性编号

	//-----------------------------单元节点信息----------------------------------
	//----------------node----------------
	for (size_t i = 0; i < node.n_rows; i++)
	{
		ID++;
		if (node(i, 3) == 1)
		{
			ofs << "GRID*   " << setw(16) << ID << setw(32) << node(i, 0)
				<< setw(16) << node(i, 1) << setw(8) << ID << endl;
			ofs << "*" << setw(7) << ID << setw(16) << node(i, 2) << endl;
		}
	}
	//----------------element----------------
	for (size_t i = 0; i < elem_aero.n_rows; i++)//气动壳单元
	{
		PID = (int)elem_aero(i, 4) + 1;
		ofs << "CTRIA3  " << setw(8) << ++EID << setw(8) << PID
			<< setw(8) << (int)elem_aero(i, 0) + 1
			<< setw(8) << (int)elem_aero(i, 1) + 1
			<< setw(8) << (int)elem_aero(i, 2) + 1 << endl;
		ofs << "CTRIA3  " << setw(8) << ++EID << setw(8) << PID
			<< setw(8) << (int)elem_aero(i, 0) + 1
			<< setw(8) << (int)elem_aero(i, 2) + 1
			<< setw(8) << (int)elem_aero(i, 3) + 1 << endl;
	}
	ofs.close();
	return 0;
}

int StructPart::subduction_node_f(map<int, Point>& nf)
{
	if (nf.size() != node_f_initial.size())
	{
		cout << "函数 StructPart::subduction_node_f 出错！气动力列表无法对应" << endl;
		return -1;
	}
	else
	{
		map<int, Point>::iterator it_nf = nf.begin();
		map<int, Point>::iterator it_n_f_i = node_f_initial.begin();
		for (map<int, Point>::iterator it_n_f = node_f.begin(); it_n_f != node_f.end(); it_n_f++, it_nf++, it_n_f_i++)
		{
			//it_0->second = (it_0->second) * (it_0->second.getP());
			it_n_f->second = (it_n_f_i->second) - (it_nf->second);
			/*it_n_f->second.setP(1);*/
		}
	}
	return 0;
}

int StructPart::ijk2ID(int i, int j, int k, int numi, int numj, int numk)//用于计算三阶张量拉直后的id
{
	return k * numi * numj + j * numi + i;
}

int StructPart::updateNode(const map<int, Point&> Displacements)
{
	int state = 0;//若返回0说明收敛
	double deltaDisp = 0;
	for (size_t i = 0; i < node.size(); i++)
	{
		int ID = i + 1;
		map<int, Point&>::const_iterator it = Displacements.find(ID);
		if (it != Displacements.end())//如果节点位移存在
		{
			//Point<int> node_0(node(i, 0), node(i, 1), node(i, 2));
			double delta = it->second.Norm1();// / node_0.Norm1();//位移量
			if (delta > DISP_CONDITIONAL_CONVERGENCE)//1e-3
			{
				state = 1;//说明位移明显 需要继续迭代
			}
			node(i, 0) += it->second.getX();
			node(i, 1) += it->second.getY();
			node(i, 2) += it->second.getZ();
		}
	}
	return state;
}

int StructPart::updateDispNode(const map<int, Point>& Displacements)
{
	int state = 0;//若返回0说明收敛
	double deltaDisp = 0;
	mat disp_node_initial;//存储上一步的受力后的节点坐标
	if (disp_node.n_rows == 0)
	{
		disp_node_initial = node;
	}
	else
	{
		disp_node_initial = disp_node;
	}
	disp_node = node;
	for (size_t i = 0; i < node.size(); i++)
	{
		int ID = i + 1;
		map<int, Point>::const_iterator it = Displacements.find(ID);
		if (it != Displacements.end())//如果节点位移存在
		{
			disp_node(i, 0) += it->second.getX();
			disp_node(i, 1) += it->second.getY();
			disp_node(i, 2) += it->second.getZ();

			double delta = Point(disp_node.row(i) - disp_node_initial.row(i)).Norm1();
			if (delta > DISP_CONDITIONAL_CONVERGENCE)//1e-3
			{
				state = 1;//说明位移明显 需要继续迭代
			}
		}
	}
	return state;
	//return 0;
}

int Struct2D::SaveAsTecplot(string dataname)
{
	ofstream ofs;
	ofs.open(dataname + "_strc2D.dat", ios::trunc);
	for (size_t i = 0; i < PartList.size(); i++)
	{
		PartList[i].SaveAsTecplot(ofs);
	}
	ofs.close();
	return 0;
}

int Struct2D::SaveAsNastran(string fileName, int SOL)
{
	//for (int i = 0; i < 1; i++)
	int i = 0;
	{
		PartList[i].SaveAsNastran(fileName + "_" + to_string(i), SOL);
	}
	return 0;
}

int Struct2D::AeroAnalysis(string fileName, string exepath)
{
	int i = 0;
	{
		// 
		//AVL计算气动力
		string aeropath = fileName + "_avl";
		int a = _mkdir(aeropath.c_str());
		int state = PartList[i].calcAeroForce_AVL(aeropath, exepath);
		if (state < 0)
		{
			cout << "气动力计算失败！请检查AVL相关设置" << endl;
			return -2;
		}

		PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_AvlRes.txt");
	}
	return 0;
}

int Struct2D::AeroelasticAnalysis(string fileName, string exepath, const int iteration_num)
{
	//const int maxIterationNum = iteration_num;//最大静气弹迭代次数
	//for (int i = 0; i < 1; i++)
	int i = 0;
	{
		int iterationNum = 0;
		int updatestate = 1;
		PartList[i].disp_node.clear();
		while (updatestate == 1)
		{
			iterationNum++;

			// 
			//AVL计算气动力
			string aeropath = fileName + "_avl";
			int a = _mkdir(aeropath.c_str());
			int state = PartList[i].calcAeroForce_AVL(aeropath, exepath);//
			if (state < 0)
			{
				cout << "气动力计算失败！请检查AVL相关设置" << endl;
				//返回错误结果的离谱值
				PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_A&S_Res.txt", true);
				return -2;
			}
			//面元法计算气动力
			//PartList[i].calcAeroForce(aeropath);
			//
			//PartList[i].setProperty(p);
			PartList[i].SaveAsNastran(fileName + "_" + to_string(i), 101);
			myNastran nas;
			nas.NasPath = "E:\\Nastran2020\\bin\\nastranw.exe";
			nas.BDFPath = fileName + "_" + to_string(i) + "_Header.bdf";
			nas.CalcFilePath();
			nas.NastranCalc();
			bool isOnlyReadDisp = true;//只读取节点位移数据
			nas.ReadResPCH(isOnlyReadDisp);


			int updatestate = PartList[i].updateDispNode(nas.GetDisplacements());
			bool if_convergent = false;
			if (iteration_num == 1)
			{// 如果最大迭代次数为1 则一步直接输出结果
				if_convergent = true;
			}
			else if (updatestate == 1)
			{// 否则判断收敛性
				cout << "迭代未收敛；已迭代次数 = " << iterationNum << endl;
			}
			else if (iterationNum > iteration_num)
			{// 若迭代次数超过最大迭代次数限制，输出假结果
				cout << "迭代已超过最大迭代次数；已迭代次数 = " << iterationNum << endl;
				//返回错误结果的离谱值
				PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_A&S_Res.txt", true);
				break;
			}
			else
			{
				if_convergent = true;
			}
			// 是否输出结果并退出分析
			if (if_convergent)
			{
				////多余的添加 屈曲分析
				//PartList[i].SaveAsNastran(fileName + "_" + to_string(i) + "_105", 105);
				//myNastran nas105;
				//nas105.NasPath = "E:\\Nastran2020\\bin\\nastranw.exe";
				//nas105.BDFPath = fileName + "_" + to_string(i) + "_105" + "_Header.bdf";
				//nas105.CalcFilePath();
				//nas105.NastranCalc();
				////return -2;//计算完成一轮之后暂时退出，不进行循环
				////
				if (iteration_num == 1)
					cout << "静力分析完成！ 迭代次数 = 1" << endl;
				else
					cout << "迭代已收敛！迭代次数 = " << iterationNum << endl;

				//读取数据
				bool isOnlyReadDisp = false;//读取节点位移、单元应变和应力数据
				nas.ReadResPCH(isOnlyReadDisp);
				PartList[i].MaxElemStress = nas.GetMaxElemStress();
				PartList[i].MaxElemStrain = nas.GetMaxElemStrain();
				PartList[i].MaxNodeDisp = nas.GetMaxDisp();
				map<int, Point> disp = nas.GetDisplacements();
				int numi = PartList[i].getnumi();
				int numj = PartList[i].getnumj();
				int numk = PartList[i].getnumk();
				int pointID1 = StructPart::ijk2ID(0, numj - 1, 0, numi, numj, numk) + 1;
				int pointID2 = StructPart::ijk2ID(numi - 1, numj - 1, 0, numi, numj, numk) + 1;
				PartList[i].MaxTwisting = abs(disp[pointID1].getY() - disp[pointID2].getY());
				PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_A&S_Res.txt");
				break;
			}
		}
	}
	return 0;
}

double Struct2D::CalcTotalMass(string outputFilePath /*= ""*/)
{
	double totalmass = 0;
	int i = 0;
	{
		double partmass = PartList[i].calcMass();
		if (partmass < 0)
		{
			return -1;
		}
		else
		{
			totalmass += partmass;
		}
	}
	if (outputFilePath.size() != 0)
	{
		ofstream ofs;
		ofs.open(outputFilePath, ios::trunc);
		ofs << "TotalMass = " << totalmass << endl;
		ofs.close();
	}
	return totalmass;
}

