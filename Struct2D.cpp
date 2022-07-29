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
		cout << "MakeStruct2D_byNum���󣡣� X������һά������" << endl;
		return -1;
	}
	if (siteZ_2D.n_rows > 1 && siteZ_2D.n_cols > 1)
	{
		cout << "MakeStruct2D_byNum���󣡣� Z������һά������" << endl;
		return -1;
	}
	return 0;
}

void StructPart::calcNode()
{
	node = zeros(numi * numj * numk, 5);
	//�ڵ㰴������
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
				//node(id, 3) = node_2D[k].ifuse(i, j);//��¼�ڵ��Ƿ�ʹ��
				node(id, 4) = id;//��¼�ڵ��ţ��ظ��Ľڵ��Ż���ظ�
				id++;
			}
		}
	}
}

void StructPart::clearUniqueNode()
{
	//������� ͨ��ֻ�� ���jǰ��������ٸ������ظ���
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
	//3����������ĸ���λ�þ���
	siteX_2D = unique(siteX_2D);
	siteZ_2D = unique(siteZ_2D);
	//�������� �������ȷ�����ݵ� X����
	//ûд Ĭ�����붼����ȷ�İ���
	//4�������Ԫ
	int PID = 0;
	//������Ԫ
	int PIDnum = (numj - 1) / 5;
	int aeroElemNum = 2 * (numi - 1) * (numj - 1);
	elem_aero = zeros(aeroElemNum, 6);
	//�±���
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
				elem_aero(elemID, 4) = j / PIDnum;//�����д���������Ա��PID
			}

			mat area =
				cross(node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2)) +
				cross(node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2)); //���ı����ĸ��������̶Ȳ���ʱ����ʹ�ô˹�ʽ
			double S = Point(area).Norm2() / 2.0;
			elem_aero(elemID, 5) = S;//�����д���Ԫ���
			elemID++;
		}
	}
	//�ϱ���
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
				elem_aero(elemID, 4) = j / PIDnum;//�����д���������Ա��PID
			}
			mat area = cross(node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 1))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 2))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2))
				+ cross(node.row(((size_t)elem_aero(elemID, 3))).cols(0, 2), node.row(((size_t)elem_aero(elemID, 0))).cols(0, 2));//���ı����ĸ��������̶Ȳ���ʱ����ʹ�ô˹�ʽ
			double S = Point(area).Norm2() / 2.0;
			elem_aero(elemID, 5) = S;//�����д���Ԫ���
			elemID++;
		}
	}
	//
	map<int, PSHELL> ps_list = p.getPSHELLlist();
	// 
	//��ֱ������ֲ��Ľṹ��Ԫ siteX_2D������ֲ�λ��x����(����)
	int strcXElemNum = (numi - 1) * siteX_2D.n_elem * (numk - 1);
	elem_strcX = zeros(strcXElemNum, 6);
	elemID = 0;
	PID = 6;
	if (isFixedMass)//���isFixedMassΪ�ǣ�����²�������
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
	for (size_t jid = 0; jid < siteX_2D.n_elem; jid++)//ÿһ���ṹ��
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
					elem_strcX((size_t)elemID, 4) = PID;//�����д���������Ա��
					mat area =
						cross(node.row(((size_t)elem_strcX(elemID, 0))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 1))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 1))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 2))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 2))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 3))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcX(elemID, 3))).cols(0, 2), node.row(((size_t)elem_strcX(elemID, 0))).cols(0, 2)); //���ı����ĸ��������̶Ȳ���ʱ����ʹ�ô˹�ʽ
					double S = Point(area).Norm2() / 2.0;
					elem_strcX(elemID, 5) = S;//�����д���Ԫ���
					elemID++;
				}
			}
		}
	}
	//��ֱ�ں���ֲ��Ľṹ��Ԫ siteZ_2D������ֲ�λ��z����(����)
	int strcZElemNum = siteZ_2D.n_elem * (numj - 1) * (numk - 1);
	elem_strcZ = zeros(strcZElemNum, 6);
	elemID = 0;
	PID = 7;
	if (isFixedMass)//���isFixedMassΪ�ǣ�����²�������
	{
		ps_list[PID + 1].T = ps_list[PID + 1].T * (double)2 / (double)siteZ_2D.n_elem;
	}
	const int spar_end_id = 100;
	PSHELL p_spar_end = p.getPSHELL(PID+1);
	p_spar_end.T *= z_T_ratio;
	p_spar_end.PID = spar_end_id + 1;
	p.add_PSHRLL(move(p_spar_end));
	for (size_t iid = 0; iid < siteZ_2D.n_elem; iid++)//ÿһ���ṹ��
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
					elem_strcZ((size_t)elemID, 4) = PID;//�����д���������Ա��
					mat area =
						cross(node.row(((size_t)elem_strcZ(elemID, 0))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 1))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 1))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 2))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 2))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 3))).cols(0, 2)) +
						cross(node.row(((size_t)elem_strcZ(elemID, 3))).cols(0, 2), node.row(((size_t)elem_strcZ(elemID, 0))).cols(0, 2)); //���ı����ĸ��������̶Ȳ���ʱ����ʹ�ô˹�ʽ
					double S = Point(area).Norm2() / 2.0;
					elem_strcZ(elemID, 5) = S;//�����д���Ԫ���
					elemID++;
				}
			}
		}
	}
	if (isFixedMass)//���isFixedMassΪ�ǣ�����²�������
	{
		p.SetPSHELLList(ps_list);
	}
}

void StructPart::AddOriginRotation(mat Origin, mat Rotation)
{
	//----------4�������ת��ƽ��----------
	if (Rotation.n_elem == 3)
	{
		double RotationX = Rotation(0) / 180.0 * pi;
		double RotationY = Rotation(1) / 180.0 * pi;
		double RotationZ = Rotation(2) / 180.0 * pi;

		//��ԭ�Ƕ� ͨ��Rotation Matrix
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
	//��ӡ����
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
	//����������Ϣ����
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
	//��������
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
		nf = node_f_initial;//�ݴ���һ����������
		int state = ac.getAeroForce(node_f_initial);//����������
		subduction_node_f(nf);
	}

	return 0;
}

int StructPart::calcAeroForce_AVL(string aeroPath, string exepath)
{
	
	mat aero_node = (disp_node.n_rows == 0) ? node : disp_node;
	//�ȴ�������
	AVLInfo ai;
	ai.sectionInfo.resize(numj);
	// ofstream ofs;
    // ofs.open("log_max_y.txt", ios::app);
	// ofs << "---------------\n";
	for (int j = 0; j < numj; j++)
	{
		AVLInfo::AVLSectionInfo &section = ai.sectionInfo[j];
		const int firstNodeID = getPointID(0, j, 0);//���Ͷ��������
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
		//�±���ڵ�
		int k = 0;
		for (int i = numi - 1; i > 0; i--)
		{
			int nodeID = getPointID(i, j, k);
			Point temp = Point(aero_node(nodeID, 0), aero_node(nodeID, 1), aero_node(nodeID, 2)) - firstNode;

			maxDanweiY = (maxDanweiY < abs(temp.getY() / section.Scale)) ? abs(temp.getY() / section.Scale) : maxDanweiY;
			// section.pointList.push_back(temp / section.Scale);
			section.pointList.push_back(temp);
		}
		//�ϱ���ڵ�
		k = numk - 1;
		for (int i = 0; i < numi; i++)//iĩֵȡ1����Ϊ�����ظ���д��0��0����
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
	//�����������
	
	ai.Nchordwise = 20;
	ai.Nspanwise = numj;// numj;//121;//
	ai.Sref = 231.75;
	ai.Cref = 15;
	ai.Angle = 4;//4�������������õģ�4.5��������֤�õ�
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
		if (num > 20)//������ʮ��֮���޷�����
		{
			return -2;//�������޷�����
		}
	}
	vector<XFoilInfo> xfoilList(numj);//XFoil����������
	vector<vector<Field<double>>> Xfoil_res(numj);
	// for (size_t i = 0; i < (size_t)1; i++)//(size_t)numj; i++)
	// {
	// 	xfoilList[i].foilName = ai.sectionInfo[i].sectionName+"_sfoil";
	// 	xfoilList[i].foilList = ai.sectionInfo[i].pointList;
	// 	xfoilList[i].alpha = {ai.Alpha, ai.Alpha+3, 1};
	// 	xfoilList[i].Ma = ai.Ma;
	// 	xfoilList[i].VISC = 500000;//������һ����ŵ��
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
	if (elem_aero.n_cols < 6 || elem_strcX.n_cols < 6 || elem_strcZ.n_cols < 6)//���û�е�����
	{
		cout << "��ǰû�м��㵥Ԫ������޷��������� ����ֵ��-1.0" << endl;
		return -1.0;
	}
	if (!p.isNotEmpty())
	{
		cout << "��ǰû�����ò������ԣ��޷��������� ����ֵ��-1.0" << endl;
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
	mat info = zeros(elem_aero.n_rows / 2, 3);//���ڴ洢��Ԫ�������fx��fy
	//ע�� AVL����ϵ��CSTʹ�õ�����ϵy��z���Ƿ���
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//ֻҪ�ϱ���
	{
		//�������
		mat node0 = (node.row((size_t)elem_aero(i, 0)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 0)).cols(0, 2)) / 2.0;
		node0(1) = Interpolation1(StripList, StripZleList, node0(2));
		mat node1 = (node.row((size_t)elem_aero(i, 1)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 1)).cols(0, 2)) / 2.0;
		node1(1) = Interpolation1(StripList, StripZleList, node1(2));
		mat node2 = (node.row((size_t)elem_aero(i, 2)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 2)).cols(0, 2)) / 2.0;
		node2(1) = Interpolation1(StripList, StripZleList, node2(2));
		mat node3 = (node.row((size_t)elem_aero(i, 3)).cols(0, 2) + node.row((size_t)elem_aero(i + elem_aero.n_rows / 2, 3)).cols(0, 2)) / 2.0;
		node3(1) = Interpolation1(StripList, StripZleList, node3(2));
		mat area = cross(node0, node1) + cross(node1, node2) + cross(node2, node3) + cross(node3, node0);//���ı����ĸ��������̶Ȳ���ʱ����ʹ�ô˹�ʽ

		double S = Point(area).Norm2() / 2.0;

		//����ο���ѹ��ϵ��
		mat nodeMid = (node0 + node3) * 3.0 / 8.0 + (node1 + node2) * 1.0 / 8.0;//ȡ��Ԫ��ǰ1/4λ�ô��ĵ�Ϊѹ���ο���
		double x = nodeMid(0);
		double z = nodeMid(2);
		double px = Interpolation2(ChordList, StripList, CdList, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z), z);
		double py = Interpolation2(ChordList, StripList, ClList, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z), z);
		areas += S;
		Fx += S * px;
		Fy += S * py;
		info(i, 0) = S;
		info(i, 1) = S * px;//����
		info(i, 2) = S * py;//����
	}
	double ratioFx = res[1].CDsurf * res[1].Surface_area / Fx;
	double ratioFy = res[1].CLsurf * res[1].Surface_area / Fy;
	info.col(1) *= ratioFx;
	info.col(2) *= ratioFy;
	nf.clear();
	double _q = res[1].ai._q;
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//ֻҪ�ϱ���
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
	//�ϱ��������״��������ȡֵ��0-1֮�� PֵΪcp
	vector<Field<double>> foilUpp = getPlistFromfile("E:\\Matlab-Xfoil\\DWPupp.txt");//���������һ���Ľڵ���������Ͷ�Ӧ��ѹ��ϵ��
	vector<Field<double>> foilLow = getPlistFromfile("E:\\Matlab-Xfoil\\DWPlow.txt");//���������һ���Ľڵ���������Ͷ�Ӧ��ѹ��ϵ��
	//��ȡ���±����������ʵCP�ֲ�ֵ
	auto clac = [](const vector<Field<double>> &foil, vector<double> &xsite, vector<double> &cp)
	{
		xsite.resize(foil.size());
		cp.resize(foil.size());
		for (size_t i = 0; i < foil.size(); i++)
		{
			xsite[i] = foil[i].getX() / foil[foil.size()-1].getX();//��ù�һ����x�����б�
			//ysite.push_back(foil[i].getY());
			cp[i] = foil[i].getP();//���ѹ��ϵ���б�
		}
	};

	vector<double> xsiteupp, /*ysiteupp,*/ cpupp;
	clac(foilUpp, xsiteupp, cpupp);//��ȡ���±����������ʵCP�ֲ�ֵ
	vector<double> xsitelow, /*ysitelow,*/ cplow;
	clac(foilLow, xsitelow, cplow);//��ȡ���±����������ʵCP�ֲ�ֵ
	//����xfoil���������������ϵ��
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
	calcS(foilUpp);//����xfoil���������������ϵ��
	calcS(foilLow);//����xfoil���������������ϵ��
	//AVL������
	vector<double> StripList(res[1].Spanwise);
	vector<double> StripXleList(res[1].Spanwise);
	vector<double> StripZleList(res[1].Spanwise);
	vector<double> StripChordList(res[1].Spanwise);

	vector<double> StripClList(res[1].Spanwise);
	vector<double> StripArea(res[1].Spanwise);
	vector<double> StripRatioList(res[1].Spanwise);
	for (int i = 0; i < res[1].Spanwise; i++)//����ÿһ������
	{
		double siteY = res[1].stripRes[res[1].FirstStrip + i].elemRes.begin()->second.site.getY();
		StripList[i] = siteY;//ÿһ���������ʼλ��
		StripXleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getX();
		StripZleList[i] = res[1].stripRes[res[1].FirstStrip + i].Orgin.getZ();
		StripChordList[i] = res[1].stripRes[res[1].FirstStrip + i].Chord;//ÿһ������Ŀ��

		StripClList[i] = res[1].stripRes[res[1].FirstStrip + i].cl;//ÿһ�������cl
		StripArea[i] = res[1].stripRes[res[1].FirstStrip + 1].StripArea;//ÿһ������Ĳο����
		StripRatioList[i] = StripClList[i] / cl_S;//ÿһ����������ϵ����XFoil�����cl�ı�ֵϵ��
	}

	double areas = 0;
	double Fx = 0;
	double Fy = 0;
	mat infoUpp = zeros(elem_aero.n_rows / 2, 3);//���ڴ洢��Ԫ�������fx��fy
	mat infoLow = zeros(elem_aero.n_rows / 2, 3);//���ڴ洢��Ԫ�������fx��fy
	enum class Side{
		upp,
		low,
	};
	using vec_d = const vector<double> &;
	//ע�� AVL����ϵ��CSTʹ�õ�����ϵy��z���Ƿ���
	auto calcElemForce = [&StripList,&StripXleList,&StripChordList,&StripRatioList,this](Side side,vec_d xsite,vec_d cp,mat&info,double&areas,double&Fx,double&Fy) {	 
		for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//ֻҪ�ϱ��� ÿ����Ԫ
		{
			const int id = (side == Side::upp) ? i : (i + elem_aero.n_rows / 2);
			const int sign = (side == Side::upp) ? 1 : -1; //�±���ȡy�ĸ���
			//�������
			const mat node0 = node.row((size_t)elem_aero(id, 0)).cols(0, 2);//��õ�Ԫ���ĸ��ڵ�����
			const mat node1 = node.row((size_t)elem_aero(id, 1)).cols(0, 2);
			const mat node2 = node.row((size_t)elem_aero(id, 2)).cols(0, 2);
			const mat node3 = node.row((size_t)elem_aero(id, 3)).cols(0, 2);
			const double dz = (node3(2) - node0(2)) / 2.0 + (node2(2) - node1(2)) / 2.0;
			const double dx = (node1(0) - node0(0)) / 2.0 + (node2(0) - node3(0)) / 2.0;
			const double dy = ((node1(1) - node0(1)) / 2.0 + (node2(1) - node3(1)) / 2.0) * sign;

			//����ο���ѹ��ϵ��
			const mat nodeMid = (node0 + node3) * 3.0 / 8.0 + (node1 + node2) * 1.0 / 8.0;//ȡ��Ԫ��ǰ1/4λ�ô��ĵ�Ϊѹ���ο���
			const double x = nodeMid(0);
			const double z = nodeMid(2);
			const double cp_temp = Interpolation1(xsite, cp, (x - Interpolation1(StripList, StripXleList, z)) / Interpolation1(StripList, StripChordList, z));
			const double ratio_temp = Interpolation1(StripList, StripRatioList, z);
			areas += dz * dx;
			Fx += cp_temp * ratio_temp * dz * dy;
			Fy += cp_temp * ratio_temp * dz * dx;
			info(i, 0) = ratio_temp * dz * dx;
			info(i, 1) = cp_temp * ratio_temp * dz * dy;//����
			info(i, 2) = cp_temp * ratio_temp * dz * dx;//����
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
	for (size_t i = 0; i < elem_aero.n_rows / 2; i++)//ֻҪ�ϱ���
	{
		//�ϱ���
		nf[1 + (int)elem_aero(i, 0)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 3.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 1)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 3.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 2)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 1.0 / 8.0 * _q;/// 8.0
		nf[1 + (int)elem_aero(i, 3)] += Field<double>(infoUpp(i, 1), infoUpp(i, 2), 0) * 1.0 / 8.0 * _q;/// 8.0
		//�±���
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
	ofs << fixed << setprecision(8);//�����������Ϊ8
	ofs1 << showpoint << setprecision(2);//ǿ�����С����
	ofs2 << fixed << setprecision(5);//�����������Ϊ5
	int ID = 0;//�ڵ���
	int EID = 0;//��Ԫ���
	int PID = 0;//����or�ǵ�Ԫ���Ա��
	int SID = 1;//�غ�orԼ���������

	//-----------------------------��Ԫ�ڵ���Ϣ----------------------------------
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
	//λ��Լ�����
	int SID_SPC = 2;
	int SPC1num = numi * 2 - 2;
	ofs2 << "SPC1    " << setw(8) << SID_SPC << setw(8) << "123456";
	for (int n = 0; n < SPC1num; n++)//�������Ƥ��Լ��
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

	for (size_t m = 0; m < siteZ_2D.n_elem; m++)///*�����ĸ���*/
	{
		int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//��n�����ĵ�һ����Ԫ�ĵ�Ԫid
		for (int k = 0; k < numk - 1; k++)
		{
			int n = SPC1num + m * (numk - 1) + k;//�������ѭ���е�n
			ofs2 << setw(8) << (int)elem_strcZ(firstID + k, 3) + 1;
			if (n % 8 == 5)	ofs2 << endl << "        ";
		}
	}

	//----------------element----------------
	for (size_t i = 0; i < elem_aero.n_rows; i++)//�����ǵ�Ԫ
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

	for (size_t i = 0; i < elem_strcX.n_rows; i++)//�ṹ�ǵ�Ԫ ����
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
	for (size_t i = 0; i < elem_strcZ.n_rows; i++)//�ṹ�ǵ�Ԫ ����
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
		//ִ�п��Ʋ���
		//ofs3 << "NASTRAN parallel=8" << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "SOL 101" << endl;
		//������������
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SUBCASE 1" << endl;//��һ������
		ofs3 << "   SUBTITLE=Default" << endl;
		ofs3 << "   LOAD = " << SID_F_all << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "   SPC  = " << SID_SPC << endl;//��Ӧ��Ϲ�����SID

		ofs3 << "   STRESS(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//�����ԪӦ�� ,PUNCH
		//ofs3 << "   STRAIN(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//
		ofs3 << "   DISPLACEMENT(SORT1,PUNCH,REAL)=ALL" << endl;//
		 //ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;//
		//BEGIN BULK
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		ofs3 << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		//ofs3 << "PARAM,INREL,-2" << endl;//��������ͷ�

		ofs3 << "include '" << (names + "_mesh.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << (names + "_property.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << (names + "_force.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "ENDDATA" << endl;
		break;

	case 105://��������
		ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//ִ�п��Ʋ���
		ofs3 << "SOL 105" << endl;
		ofs3 << "TIME 10000" << endl;
		//������������
		ofs3 << "CEND" << endl;
		ofs3 << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		ofs3 << "ECHO = NONE" << endl;
		ofs3 << "SPC  = " << SID_SPC << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "SUBCASE 1" << endl;//��һ������
		ofs3 << "   LOAD = " << SID_F_all << endl;//��Ӧ��Ϲ�����SID
		ofs3 << "SUBCASE 2" << endl;//�ڶ�������
		ofs3 << "   METHOD = 10" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "$$                                Bulk Data Cards                               $" << endl;
		ofs3 << "$$------------------------------------------------------------------------------$" << endl;
		ofs3 << "BEGIN BULK" << endl;
		ofs3 << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		ofs3 << "PARAM   PRTMAXIM YES" << endl;
		ofs3 << "include '" << (names + "_mesh.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << (names + "_property.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "include '" << (names + "_force.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		ofs3 << "EIGRL   10                      5" << endl;//ģ̬������ƹؼ��� ���ǰ���ģ̬
		ofs3 << "ENDDATA" << endl;
		break;

	default:
		cout << "��ǰnastran�������Ϊ��" << SOL << "���������λ��void StructPart::SaveAsNastran(string fileName)" << endl;
		break;
	}
	ofs.close();
	ofs1.close();
	ofs2.close();
	ofs3.close();
	cout << "struct-2D����(*.BDF)������ɣ�(SOL = " << SOL << ")" << endl;
}

void StructPart::SaveAsNastran(string fileName, int SOL)
{
	NasPrinter np;

	// ofstream ofs, ofs1, ofs2;
	// ofs.open(fileName + "_mesh.bdf", ios::trunc);
	// ofs1.open(fileName + "_property.bdf", ios::trunc);
	// ofs2.open(fileName + "_force.bdf", ios::out);
	// string names = fileName.substr(fileName.find_last_of("/") + 1);
	// ofs << fixed << setprecision(8);//�����������Ϊ8
	// ofs1 << showpoint << setprecision(2);//ǿ�����С����
	// ofs2 << fixed << setprecision(5);//�����������Ϊ5
	int ID = 0;//�ڵ���
	int EID = 0;//��Ԫ���
	int PID = 0;//����or�ǵ�Ԫ���Ա��
	int SID = 1;//�غ�orԼ���������

	//-----------------------------��Ԫ�ڵ���Ϣ----------------------------------
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
	//λ��Լ�����
	int SID_SPC = 2;
	int SPC1num = numi * 2 - 2;
	vector<int> spc_list(0);
	spc_list.reserve(SPC1num + siteZ_2D.n_elem * (numk - 1));
	for (int n = 0; n < SPC1num; n++)//�������Ƥ��Լ��
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

	for (size_t m = 0; m < siteZ_2D.n_elem; m++)//�����ĸ���
	{
		int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//��n�����ĵ�һ����Ԫ�ĵ�Ԫid
		for (int k = 0; k < numk - 1; k++)
		{
			int n = SPC1num + m * (numk - 1) + k;//�������ѭ���е�n
			spc_list.push_back((int)elem_strcZ(firstID + k, 3) + 1);
		}
	}
	np.addSPC1(SID_SPC, 123456, spc_list);
	// ofs2 << "SPC1    " << setw(8) << SID_SPC << setw(8) << "123456";
	
	// for (int n = 0; n < SPC1num; n++)//�������Ƥ��Լ��
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

	// for (size_t m = 0; m < siteZ_2D.n_elem; m++)//�����ĸ���
	// {
	// 	int firstID = elem_strcZ.n_rows * m / siteZ_2D.n_elem;//��n�����ĵ�һ����Ԫ�ĵ�Ԫid
	// 	for (int k = 0; k < numk - 1; k++)
	// 	{
	// 		int n = SPC1num + m * (numk - 1) + k;//�������ѭ���е�n
	// 		ofs2 << setw(8) << (int)elem_strcZ(firstID + k, 3) + 1;
	// 		if (n % 8 == 5)	ofs2 << endl << "        ";
	// 	}
	// }

	//----------------element----------------
	for (size_t i = 0; i < elem_aero.n_rows; i++)//�����ǵ�Ԫ
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

	for (size_t i = 0; i < elem_strcX.n_rows; i++)//�ṹ�ǵ�Ԫ ����
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
	for (size_t i = 0; i < elem_strcZ.n_rows; i++)//�ṹ�ǵ�Ԫ ����
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
        (*it).second.MID2 = 1; // ȱ�������л������
        (*it).second.MID3 = 1; // ֮ǰ�ļ��㶼ʹ���������� Ӧ��������˶���ĸն�
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
		//ִ�п��Ʋ���
		//ofs3 << "NASTRAN parallel=8" << endl;//��Ӧ��Ϲ�����SID
		np.ssHeader << "SOL 101" << endl;
		//������������
		np.ssHeader << "CEND" << endl;
		np.ssHeader << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		np.ssHeader << "ECHO = NONE" << endl;
		np.ssHeader << "SUBCASE 1" << endl;//��һ������
		np.ssHeader << "   SUBTITLE=Default" << endl;
		np.ssHeader << "   LOAD = " << SID_F_all << endl;//��Ӧ��Ϲ�����SID
		np.ssHeader << "   SPC  = " << SID_SPC << endl;//��Ӧ��Ϲ�����SID

		np.ssHeader << "   STRESS(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//�����ԪӦ�� ,PUNCH
		//ofs3 << "   STRAIN(SORT1,PUNCH,REAL,VONMISES,BILIN)=ALL" << endl;//
		np.ssHeader << "   DISPLACEMENT(SORT1,PUNCH,REAL)=ALL" << endl;//
		 //ofs3 << "   DISPLACEMENT(SORT1,REAL)=ALL" << endl;//
		//BEGIN BULK
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "$$                                Bulk Data Cards                               $" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "BEGIN BULK" << endl;
		np.ssHeader << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		np.ssHeader << "PARAM   PRTMAXIM YES" << endl;
		//ofs3 << "PARAM,INREL,-2" << endl;//��������ͷ�

		// np.ssHeader << "include '" << (names + "_mesh.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		// np.ssHeader << "include '" << (names + "_property.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		// np.ssHeader << "include '" << (names + "_force.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		// np.ssHeader << "ENDDATA" << endl;
		break;

	case 105://��������
		// ofs3.open(fileName + "_Header.bdf", ios::trunc);
		//ִ�п��Ʋ���
		np.ssHeader << "SOL 105" << endl;
		np.ssHeader << "TIME 10000" << endl;
		//������������
		np.ssHeader << "CEND" << endl;
		np.ssHeader << "TITLE = MSC.Nastran job for Single StructPart" << endl;
		np.ssHeader << "ECHO = NONE" << endl;
		np.ssHeader << "SPC  = " << SID_SPC << endl;//��Ӧ��Ϲ�����SID
		np.ssHeader << "SUBCASE 1" << endl;//��һ������
		np.ssHeader << "   LOAD = " << SID_F_all << endl;//��Ӧ��Ϲ�����SID
		np.ssHeader << "SUBCASE 2" << endl;//�ڶ�������
		np.ssHeader << "   METHOD = 10" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "$$                                Bulk Data Cards                               $" << endl;
		np.ssHeader << "$$------------------------------------------------------------------------------$" << endl;
		np.ssHeader << "BEGIN BULK" << endl;
		np.ssHeader << "PARAM   POST     0" << endl;//����'.XBD'���Ϳ��ӻ��ļ�
		np.ssHeader << "PARAM   PRTMAXIM YES" << endl;
		// np.ssHeader << "include '" << (names + "_mesh.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		// np.ssHeader << "include '" << (names + "_property.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		// np.ssHeader << "include '" << (names + "_force.bdf") << "'" << endl;//ͷ�ļ��ڲ�ʹ�����·��
		np.ssHeader << "EIGRL   10                      5" << endl;//ģ̬������ƹؼ��� ���ǰ���ģ̬
		// np.ssHeader << "ENDDATA" << endl;
		break;

	default:
		cout << "��ǰnastran�������Ϊ��" << SOL << "���������λ��void StructPart::SaveAsNastran(string fileName)" << endl;
		break;
	}
	// ofs.close();
	// ofs1.close();
	// ofs2.close();
	// ofs3.close();
	string::size_type iPos = (fileName.find_last_of('\\') + 1) == 0 ? fileName.find_last_of('/') + 1 : fileName.find_last_of('\\') + 1;
	string datapath = fileName.substr(0, iPos);//��ȡ�ļ�·��
    string names = fileName.substr(iPos);
    iPos = (names.find_last_of('\\') + 1) == 0 ? names.find_last_of('/') + 1
                                               : names.find_last_of('\\') + 1;
    datapath += names.substr(0, iPos);
    names = names.substr(iPos);
	np.PrintBDF(datapath, names, NasPrinter::onefile);
	cout << "struct-2D����(*.BDF)������ɣ�(SOL = " << SOL << ")" << endl;
}

int StructPart::printAeroCalcMesh(string filepath)
{
	ofstream ofs;
	ofs.open(filepath, ios::trunc);
	ofs << fixed << setprecision(8);//�����������Ϊ8
	int ID = 0;//�ڵ���
	int EID = 0;//��Ԫ���
	int PID = 0;//����or�ǵ�Ԫ���Ա��

	//-----------------------------��Ԫ�ڵ���Ϣ----------------------------------
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
	for (size_t i = 0; i < elem_aero.n_rows; i++)//�����ǵ�Ԫ
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
		cout << "���� StructPart::subduction_node_f �����������б��޷���Ӧ" << endl;
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

int StructPart::ijk2ID(int i, int j, int k, int numi, int numj, int numk)//���ڼ�������������ֱ���id
{
	return k * numi * numj + j * numi + i;
}

int StructPart::updateNode(const map<int, Point&> Displacements)
{
	int state = 0;//������0˵������
	double deltaDisp = 0;
	for (size_t i = 0; i < node.size(); i++)
	{
		int ID = i + 1;
		map<int, Point&>::const_iterator it = Displacements.find(ID);
		if (it != Displacements.end())//����ڵ�λ�ƴ���
		{
			//Point<int> node_0(node(i, 0), node(i, 1), node(i, 2));
			double delta = it->second.Norm1();// / node_0.Norm1();//λ����
			if (delta > DISP_CONDITIONAL_CONVERGENCE)//1e-3
			{
				state = 1;//˵��λ������ ��Ҫ��������
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
	int state = 0;//������0˵������
	double deltaDisp = 0;
	mat disp_node_initial;//�洢��һ����������Ľڵ�����
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
		if (it != Displacements.end())//����ڵ�λ�ƴ���
		{
			disp_node(i, 0) += it->second.getX();
			disp_node(i, 1) += it->second.getY();
			disp_node(i, 2) += it->second.getZ();

			double delta = Point(disp_node.row(i) - disp_node_initial.row(i)).Norm1();
			if (delta > DISP_CONDITIONAL_CONVERGENCE)//1e-3
			{
				state = 1;//˵��λ������ ��Ҫ��������
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
		//AVL����������
		string aeropath = fileName + "_avl";
		int a = _mkdir(aeropath.c_str());
		int state = PartList[i].calcAeroForce_AVL(aeropath, exepath);
		if (state < 0)
		{
			cout << "����������ʧ�ܣ�����AVL�������" << endl;
			return -2;
		}

		PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_AvlRes.txt");
	}
	return 0;
}

int Struct2D::AeroelasticAnalysis(string fileName, string exepath, const int iteration_num)
{
	//const int maxIterationNum = iteration_num;//���������������
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
			//AVL����������
			string aeropath = fileName + "_avl";
			int a = _mkdir(aeropath.c_str());
			int state = PartList[i].calcAeroForce_AVL(aeropath, exepath);//
			if (state < 0)
			{
				cout << "����������ʧ�ܣ�����AVL�������" << endl;
				//���ش�����������ֵ
				PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_A&S_Res.txt", true);
				return -2;
			}
			//��Ԫ������������
			//PartList[i].calcAeroForce(aeropath);
			//
			//PartList[i].setProperty(p);
			PartList[i].SaveAsNastran(fileName + "_" + to_string(i), 101);
			myNastran nas;
			nas.NasPath = "E:\\Nastran2020\\bin\\nastranw.exe";
			nas.BDFPath = fileName + "_" + to_string(i) + "_Header.bdf";
			nas.CalcFilePath();
			nas.NastranCalc();
			bool isOnlyReadDisp = true;//ֻ��ȡ�ڵ�λ������
			nas.ReadResPCH(isOnlyReadDisp);


			int updatestate = PartList[i].updateDispNode(nas.GetDisplacements());
			bool if_convergent = false;
			if (iteration_num == 1)
			{// �������������Ϊ1 ��һ��ֱ��������
				if_convergent = true;
			}
			else if (updatestate == 1)
			{// �����ж�������
				cout << "����δ�������ѵ������� = " << iterationNum << endl;
			}
			else if (iterationNum > iteration_num)
			{// ���������������������������ƣ�����ٽ��
				cout << "�����ѳ����������������ѵ������� = " << iterationNum << endl;
				//���ش�����������ֵ
				PartList[i].printAnalysisRes(fileName + "_" + to_string(i) + "_A&S_Res.txt", true);
				break;
			}
			else
			{
				if_convergent = true;
			}
			// �Ƿ����������˳�����
			if (if_convergent)
			{
				////�������� ��������
				//PartList[i].SaveAsNastran(fileName + "_" + to_string(i) + "_105", 105);
				//myNastran nas105;
				//nas105.NasPath = "E:\\Nastran2020\\bin\\nastranw.exe";
				//nas105.BDFPath = fileName + "_" + to_string(i) + "_105" + "_Header.bdf";
				//nas105.CalcFilePath();
				//nas105.NastranCalc();
				////return -2;//�������һ��֮����ʱ�˳���������ѭ��
				////
				if (iteration_num == 1)
					cout << "����������ɣ� �������� = 1" << endl;
				else
					cout << "�������������������� = " << iterationNum << endl;

				//��ȡ����
				bool isOnlyReadDisp = false;//��ȡ�ڵ�λ�ơ���ԪӦ���Ӧ������
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

