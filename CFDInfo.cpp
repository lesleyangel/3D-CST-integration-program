#include "CFDInfo.h"
#include <fstream>
#include <iomanip>
#include <time.h>
#include <windows.h>
AeroInfo::AeroInfo()
{
	CalcBlockID.clear();
	Ma.clear();
	Alpha.clear();
	Beta.clear();
	Hm.clear();

	sRef = 1;
	cRef = 1;

	meshPath.clear();
	outputPath.clear();
	outputName = "unnameRes.txt";

	ifViscous = 0;//�Ƿ�����ճ�Լ���
	HeadRadius = 0.03;//ͷ���뾶
	WallTemperature = 300;//�����¶�
	WallEmissivity = 0.7;//���������
	AirType = 1;//����ģ�ͣ�0��1��
	HeatTransferCoefficient = 0;//����ϵ��
	MediumTemperature = 280;//����ϵ��
	InnerWallTemperature = 350;//�ڱ����¶�
}

int AeroCalc::CalcAero()
{
	cout << "$----------------------------------------$" << endl;
	cout << "$             ���������㿪ʼ             $" << endl;
	cout << "$-----                              -----$" << endl;
	//for (int i = 0; i < CSTall.size(); i++)
	int i = 0;
	{
		cout << "~~~~~~~~��" << i + 1 << "����������������������~~~~~~~~" << endl;
		string command = "aero_calc.exe " + CommandfileName;
		system(command.c_str());
		cout << endl;
	}
	cout << "$-----                              -----$" << endl;
	cout << "$             �������������             $" << endl;
	cout << "$----------------------------------------$" << endl;
	return 0;
}

int AeroCalc::getAeroForce(map<int, Point>& nodeForce)
{
	nodeForce.clear();
	string node_force_filepath = ai.outputPath + "/node_force.bdf";
	ifstream ifs(node_force_filepath);
	if (ifs.fail())
	{
		cout << "aerocalc����ļ������ڣ�" << endl;
		return -1;
	}
	else
	{
		cout << "�򿪳ɹ���" << endl;
	}


	string docPath = node_force_filepath.substr(0, node_force_filepath.find_last_of("/"));//��ȡ�ļ���·��
	if (docPath == node_force_filepath)//���·��
	{
		docPath = "";
	}
	string str_line;
	while (getline(ifs, str_line))
	{
		//string str_title = str_line.substr(0, 7);
		string strr = str_line.substr(0, 6);
		if (!str_line.substr(0, 6).compare("FORCE*"))
		{
			int nodeID; Field<double> force; double temp;
			stringstream  iss;
			iss << str_line.substr(24, 16); iss >> nodeID; iss.clear();
			iss << str_line.substr(56, 16); iss >> temp; iss.clear(); force.setP(temp);
			if (abs(force.getP()) > 1e-8  )
			{
				getline(ifs, str_line);
				iss << str_line.substr(8, 16); iss >> temp; iss.clear(); force.setX(temp);
				iss << str_line.substr(24, 16); iss >> temp; iss.clear(); force.setY(temp);
				iss << str_line.substr(40, 16); iss >> temp; iss.clear(); force.setZ(temp);
				pair<int, Point> nforce;
				nforce.first = nodeID;
				nforce.second = force * force.getP() * TIMES_FORCE;
				
				nodeForce.insert(nforce);
			}
		}
	}

	return 0;
}

int AeroCalc::checkData()
{
	//AeroInfo ai;
	if (ai.CalcBlockIDNum() == 0)						return -1;
	if (ai.CalcBlockArea.size() != ai.CalcBlockIDNum())	return -1;
	if (ai.CalcBlockType.size() != ai.CalcBlockIDNum())	return -1;
	if (ai.CalcBlockOrient.size() != ai.CalcBlockIDNum())	return -1;

	if (ai.MaNum() == 0)		return -2;
	if (ai.AlphaNum() == 0)	return -2;
	if (ai.BetaNum() == 0)		return -2;
	if (ai.HmNum() == 0)		return -2;

	if (ai.sRef <= 0)  return -2;
	if (ai.cRef <= 0)  return -2;

	if (ai.meshPath.size() == 0)		return -3;
	if (ai.outputPath.size() == 0)		return -3;
	if (ai.outputName.size() == 0)		return -3;
	
	return 0;
}

int AeroCalc::printCommand(string fileName)
{
	this->CommandfileName = fileName;
	int state = checkData();
	if (state < 0) return -1;
	//
	ofstream ofs;
	ofs.open(CommandfileName, ios::trunc);
	if (!ofs.is_open()) return -2;
	//
	ofs << "//BEGIN//" << endl;
	ofs << endl;
	ofs << "//////////////////����״̬(���飩///////////////////////" << endl;
	ofs << "mach	nAlpha   nBeta  nH x�᷽��(1Ϊͷ��ָ��β����-1Ϊ������)" << endl;
	ofs << left << setw(8) << ai.MaNum()
		<< left << setw(9) << ai.AlphaNum()
		<< left << setw(7) << ai.BetaNum()
		<< left << setw(3) << ai.HmNum() << "  1" << endl;

	ofs << "mach    //�����" << endl;
	for (int i = 0; i < ai.MaNum(); i++)		ofs << ai.Ma[i] << " " << endl;

	ofs << "alpha	//����(deg)" << endl;
	for (int i = 0; i < ai.AlphaNum(); i++)		ofs << ai.Alpha[i] << " " << endl;

	ofs << "beta 	//�໬��(deg)" << endl;
	for (int i = 0; i < ai.BetaNum(); i++)		ofs << ai.Beta[i] << " " << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////////���и߶�(m)////////////////////////" << endl;
	for (int i = 0; i < ai.HmNum(); i++)
		ofs << ai.Hm[i] << " " << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////////��������///////////////////////////" << endl;
	ofs << "��������� " << ai.CalcBlockIDNum() << endl;
	ofs << "//////�������ID" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockID[i];	ofs << endl;
	ofs << "//////�����������" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockArea[i];	ofs << endl;
	ofs << "//////�����������" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockType[i];	ofs << endl;
	ofs << "///////�����ⷨ�߷���" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockOrient[i];	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "//////////////////////������////////////////////////////" << endl;
	ofs << "���������� 0" << endl;
	ofs << "/////ÿ������������ռһ�У���������������ö�Ӧ���̶���ת��(0�ޣ�1X,2Y,3Z)������֮�䲻Ҫ���ո��ɶ��Ÿ�������ת����P2->P1����ʱ����תΪ������λm" << endl;
	ofs << "^^^��Ӧ����^^^^^^^��ת��(��1����)^^^^^^^��ת��(��2����)^^^^^^^�̶���ת��^^^��ƫ����^^^��ƫ������^^^^^^^^^^^" << endl;
	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////sRef�ο����//cRef�ο�����//�ο��ߴ�(m)////////" << endl;
	ofs << "1 1 1 0 0" << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////�Ƿ���ճ��(��Ϊ1����Ϊ0)/////////////" << endl;
	ofs << "����ճ�Լ���(0,1)     " << ai.ifViscous << endl;//
	ofs << "ͷ���뾶(m)           " << ai.HeadRadius << endl;//
	ofs << "�����¶�(K)           " << ai.WallTemperature << endl;//
	ofs << "���������            " << ai.WallEmissivity << endl;//
	ofs << "����ģ��(0,1)         " << ai.AirType << endl;//
	ofs << "����ϵ��              " << ai.HeatTransferCoefficient << endl;//
	ofs << "�����¶�(K)           " << ai.MediumTemperature << endl;//
	ofs << "�ڱ����¶�(K)         " << ai.InnerWallTemperature << endl;//
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "////////////////�Ƕ�������ѡ��(��Ϊ1����Ϊ0)////////////" << endl;
	ofs << "�Ƕ�������   �������:ʱ����(t)    �ڵ��ٶ�rpt�ĵ� " << endl;
	ofs << "0            0.0005                  velocity.rpt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////����ļ�//////////////////////////" << endl;
	ofs << "//�����ļ���λ(m/mm):		mm" << endl;
	ofs << "//��������                  CATIA" << endl;
	ofs << "�����ļ���                  " << ai.meshPath << endl;
	ofs << "���·����                  " << ai.outputPath << endl;
	ofs << "�������ļ�����            " << ai.outputName << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "�Ƿ����ѹ���ֲ��ļ�plot.plt(plot3D��ʽ):     YES" << endl;
	ofs << "�Ƿ������Ԫ���ֲ��ļ�face.dat:               NO" << endl;
	ofs << "�Ƿ�����ڵ����ֲ��ļ�node.dat:               YES" << endl;
	ofs << "�Ƿ�����¶ȷֲ��ļ�temperature.plt(��Ҫճ�Լ���) NO" << endl;
	ofs << "�Ƿ�������طֲ��ļ�qload.plt(��Ҫճ�Լ���)       NO" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << endl;
	ofs << endl;
	ofs << endl;
	ofs << "//END//" << endl;
	

	return 0;
}

AVLInfo::AVLInfo()
{
	name = "none";
	Ma = 0;
	Alpha = 0;
	V = 0;
	_q = 0;
	rho = 1.225;
	grav = 9.81;
	IYsym = IZsym = Zsym = 0;
	Sref = Cref = Bref = 1;
	Nchordwise = Nspanwise = 1;
	Angle = 0;
	Scale = Point(1, 1, 1);
	Orgin = Point(0, 0, 0);
	sectionInfo.clear();
}

void AVLInfo::AVLSectionInfo::printSectionFile(string path)
{
	ofstream ofs;
	string filepath = path + sectionName + ".dat";
	ofs.open(filepath, ios::trunc);
	ofs << sectionName << endl;
	for (size_t i = 0; i < pointList.size(); i++)
	{
		ofs << setw(16) << pointList[i].getX() << setw(16) << pointList[i].getY() << endl;
	}
	ofs.close();
}

int AVLServer::printCommand(string aeropath)
{
	this->workpath = aeropath;
	string inputRunPath = aeropath + "/AVL_run.run";
	ofstream ofs0;
	ofs0.open(inputRunPath, ios::trunc);
	ofs0.setf(ios::left);
	ofs0 << endl;
	ofs0 << " ---------------------------------------------" << endl;
	ofs0 << " Run case  1:  -CPPTEAST- " << endl;    
	ofs0 << " alpha        ->  alpha       =   " << ai.Alpha << endl;
	ofs0 << " beta         ->  beta        =   0.00000    " << endl;
	ofs0 << " " << endl;
	ofs0 << setw(59) << " alpha = " << ai.Alpha << endl;
	ofs0 << setw(59) << " beta = 0.00000     deg" << endl;
	ofs0 << setw(59) << " CL = 1" << endl;
	ofs0 << setw(59) << " CDo = 10.00000" << endl;
	ofs0 << setw(59) << " bank = 0.00000     deg" << endl;
	ofs0 << setw(59) << " elevation = 0.00000     deg" << endl;
	ofs0 << setw(59) << " heading = 0.00000     deg" << endl;
	ofs0 << " Mach = " << ai.Ma << endl;
	ofs0 << setw(59) << " velocity  =   300     m/s" << endl;
	ofs0 << " density   =   " << ai.rho << endl;
	ofs0 << setw(59) << " grav.acc. =   9.81000     m/s^2" << endl;
	ofs0.close();

	//input�ļ�
	string inputPath = aeropath + "/AVL_Input.avl";
	ofstream ofs1;
	ofs1.open(inputPath, ios::trunc);
	ofs1.setf(ios::left);
	ofs1 << partname << endl;
	ofs1 << setw(16) << ai.Ma << setw(16) << " Mach" << endl;
	ofs1 << setw(8) << ai.IYsym << setw(8) << ai.IZsym << setw(8) << ai.Zsym << setw(8) << " iYsym  iZsym  Zsym" << endl;
	ofs1 << setw(8) << ai.Sref << setw(8) << ai.Cref << setw(8) << ai.Bref << setw(8) << " Sref   Cref   Bref" << endl;
	ofs1 << setw(8) << 0 << setw(8) << 0 << setw(8) << 0 << setw(8) << " Xref   Yref   Zref" << endl;
	ofs1 << "0.00                   CDp " << endl;//������������ϵ�� Ĭ��Ϊ0��
	ofs1 << "SURFACE" << endl;
	ofs1 << "WING-body" << endl;
	ofs1 << setw(8) << ai.Nchordwise << setw(8) << "1.0" << setw(8) << ai.Nspanwise << setw(8) << "0" << endl;
	//ofs1 << "YDUPLICATE" << endl;//�Ƿ�����xz��Գ� �ɴ��������Ի�����жԳƴ���
	//ofs1 << "0.0" << endl;
	ofs1 << "SCALE" << endl;//�ߴ������
	ofs1 << setw(16) << ai.Scale.getX() << setw(16) << ai.Scale.getY() << setw(16) << ai.Scale.getZ() << endl;
	ofs1 << "TRANSLATE" << endl;//λ�õ�ƽ��
	ofs1 << setw(16) << ai.Orgin.getX() << setw(16) << ai.Orgin.getY() << setw(16) << ai.Orgin.getZ() << endl;
	ofs1 << "ANGLE" << endl;//���н���Ƕȵĸı�
	ofs1 << ai.Angle << endl;
	for (size_t i = 0; i < ai.sectionInfo.size(); i++)
	{
		ai.sectionInfo[i].printSectionFile(aeropath + "/" );
	ofs1 << "SECTION" << endl
			<< setw(16) << ai.sectionInfo[i].Orgin.getX() 
			<< setw(16) << ai.sectionInfo[i].Orgin.getY() 
			<< setw(16) << ai.sectionInfo[i].Orgin.getZ()
			<< setw(16) << ai.sectionInfo[i].Scale << setw(16) << ai.sectionInfo[i].Angle << endl;
		ofs1 << "AFIL" << endl;
		ofs1 << aeropath + "/" + ai.sectionInfo[i].sectionName + ".dat" << endl;
		//ofs1 << aeropath + "/" + "S00_.dat" << endl;
	}
	ofs1.close();
	//ɾ�����еĽ���ļ� ��ֹ���ɴ����ȡ��ԭ����
	ResFilePath = aeropath + "/ef.dat";
	remove(ResFilePath.c_str());
	//�������ļ�
	ofstream ofs2;
	CommandfileName = aeropath + "/AVL_cmd.cmd";
	
	ofs2.open(CommandfileName, ios::trunc);
	ofs2 << "load " << inputPath << endl;	//����������������
	ofs2 << "case " << inputRunPath << endl;//���빤����������
	ofs2 << "oper" << endl;					//����-�ؼ���
	ofs2 << "x" << endl;					//����������
	ofs2 << "fe" << endl;					//��Ԫ�ڵ���������
	ofs2 << ResFilePath << endl;			//�ڵ�������ļ�·��
	ofs2 << "o" << endl;					//����ԭʼ�ļ�
	ofs2 << endl;							//ȷ��
	ofs2 << "quit" << endl;					//�˳�����
	ofs2.close();
	return 0;
}

//ͨ�������е�������������
int AVLServer::CalcAero()
{
	cout << "$----------------------------------------$" << endl;
	cout << "$             ���������㿪ʼ             $" << endl;
	cout << "$-----                              -----$" << endl;
	//for (int i = 0; i < CSTall.size(); i++)
	int i = 0;
	{
		cout << "~~~~~~~~��" << i + 1 << "����������������������~~~~~~~~" << endl;
		string command = AVLexePath + " <" + CommandfileName;// +" exit";
		system(command.c_str());
		cout << endl;
	}
	cout << "$-----                              -----$" << endl;
	cout << "$             �������������             $" << endl;
	cout << "$----------------------------------------$" << endl;
	return 0;
}


void AVLServer::setAVLexePath(string exepath) 
{ 
	if (exepath.size() == 0)
		AVLexePath = "avl.exe"; 
	else 
		AVLexePath = exepath + "\\avl.exe";
}


int AVLServer::getAeroForce(map<int, AVLres_surface>& res)
{

	ifstream ifs(ResFilePath);
	if (ifs.fail())		return -1;	
	else 				cout << "��������ļ���" << ResFilePath << "�򿪳ɹ�" << endl;

	res.clear();
	string str_line;
	while (getline(ifs, str_line))
	{
		vector<string> strlist = split(str_line, ' ');
		if (strlist.size() != 0 && !strlist[0].compare("Surface"))
		{
			AVLres_surface res_surf;
			res_surf.SurfaceID = atoi(strlist[2].c_str());
			res_surf.SurfName = strlist[3];
			getline(ifs, str_line); //��һ������
			res_surf.Chordwise = atoi(str_line.substr(19, 6).c_str());
			res_surf.Spanwise = atoi(str_line.substr(37, 6).c_str());
			res_surf.FirstStrip = atoi(str_line.substr(57).c_str());
			getline(ifs, str_line);//�ڶ�������
			res_surf.Surface_area = atof(str_line.substr(19, 19).c_str());
			for (int i = 0; i < 6; i++)//�ն�����
			{
				getline(ifs, str_line);
			}
			getline(ifs, str_line);//�ھ�������
			res_surf.CLsurf = atof(str_line.substr(14, 15).c_str());
			res_surf.CDsurf = atof(str_line.substr(38).c_str());
			//getline(ifs, str_line);//��ʮ������
			//getline(ifs, str_line);//��ʮһ������
			int StripNum = 0;
			while (getline(ifs, str_line))
			{
				vector<string> strlist = split(str_line, ' ');
				if (strlist.size() != 0 && !strlist[0].compare("Strip"))
				{
					StripNum++;
					AVLres_surface::AVLres_strip res_strip;
					res_strip.StripID = atoi(str_line.substr(8, 8).c_str());
					res_strip.Chordwise = atoi(str_line.substr(29, 6).c_str());
					res_strip.FirstVortex = atoi(str_line.substr(49).c_str());
					getline(ifs, str_line);//��һ������
					res_strip.Orgin.setX(atof(str_line.substr(9, 14).c_str()));
					res_strip.Chord = atof(str_line.substr(37, 13).c_str());
					res_strip.Incidence = atof(str_line.substr(62, 11).c_str());
					getline(ifs, str_line);//�ڶ�������
					res_strip.Orgin.setY(atof(str_line.substr(9, 14).c_str()));
					res_strip.StripWidth = atof(str_line.substr(37, 13).c_str());
					res_strip.StripArea = atof(str_line.substr(62).c_str());
					getline(ifs, str_line);//����������
					res_strip.Orgin.setZ(atof(str_line.substr(9, 14).c_str()));
					getline(ifs, str_line);//����������
					getline(ifs, str_line);//����������
					res_strip.cl = atof(str_line.substr(9, 17).c_str());
					res_strip.cd = atof(str_line.substr(31, 16).c_str());

					while (getline(ifs, str_line))
					{
						if (str_line.size() != 0 && !str_line.compare("    I        X           Y           Z           DX        Slope        dCp"))
						{
							break;
						}
					}
					for (int n = 0; n < res_strip.Chordwise; n++)
					{
						AVLres_surface::AVLres_strip::AVLres_elem res_elem;
						res_elem.elemID = res_strip.FirstVortex + n;
						getline(ifs, str_line);//��ȡһ������
						res_elem.site.setX(atof(str_line.substr(5, 12).c_str()));
						res_elem.site.setY(atof(str_line.substr(17, 12).c_str()));
						res_elem.site.setZ(atof(str_line.substr(29, 12).c_str()));
						res_elem.Dx = atof(str_line.substr(41, 12).c_str());
						res_elem.Slope = atof(str_line.substr(53, 12).c_str());
						res_elem.dCp = atof(str_line.substr(65, 12).c_str());
						res_strip.elemRes[res_elem.elemID]= res_elem;
					}
					// res_surf.stripRes.insert(pair<int, AVLres_strip>(res_strip.StripID, res_strip));
					res_surf.stripRes[res_strip.StripID] = res_strip;
					if (StripNum == res_surf.Spanwise)//������������ָ������ʱ ����ѭ��
					{
						break;
					}
				}
			}
			res_surf.ai = ai;
			res.insert(pair<int, AVLres_surface>(res_surf.SurfaceID, res_surf));
		}
	}
	return 0;
}



int XFoilSolver::solve()
{
	cout << "XFoil ���ڼ���...";
	printShape();
	printInput();
	remove(m_info.resultPath().c_str());
	//
	const string command = ExePath + " <" + m_info.inputPath();// +" exit";
	system(command.c_str());
	//
	const clock_t start = clock();
	while (readRes() != 0)
	{
		if(clock()-start > 1000)
		{
			cout << "XFoil����ʧ�ܣ�δ�����ļ���"+m_info.resultPath() << endl;
			return -1;
		}
	}
	
	//
	cout << "����" << endl;
	return 0;
}

int XFoilSolver::printShape()
{
	ofstream ofs;
	//string path = m_info.shapePath();
	ofs.open(m_info.shapePath(), ios::trunc);
	ofs.setf(ios::left);
	//
	ofs << m_info.foilName << endl;
	for (size_t i = 0; i < m_info.foilList.size(); i++)
	{
		ofs	<< setw(20) << m_info.foilList[i].getX()
			<< setw(20) << m_info.foilList[i].getY() << endl;
	}
	ofs.close();
	return 0;
}//��ӡ�����ļ�

int XFoilSolver::printInput()
{
	ofstream ofs;
	ofs.open(m_info.inputPath(), ios::trunc);
	//
	ofs << std::left
		<< "PLOP\n"		//���ò���ʾͼ�δ���
		<< "G\n\n"		//���ò���ʾͼ�δ���
		<< "LOAD\n"
		<< m_info.shapePath() << "\n"
		<< "PPAR\n"
		<< setw(10) << "N" << m_info.outputNum << "\n\n\n"
		<< "PANE\n"
		<< "OPER\n"
		<< setw(10) << "VISC" << m_info.VISC << "\n"
		<< setw(10) << "M" << m_info.Ma << "\n"
		<< setw(10) << "ITER" << m_info.ITER << "\n"
		<< setw(10) << "ASEQ"
		<< setw(18) << get<0>(m_info.alpha)
		<< setw(18) << get<0>(m_info.alpha)
		<< setw(18) << 1.0 << "\n"
		<< "CPWR\n\n\n"
		<< "Quit";
	ofs.close();
	return 0;
}

int XFoilSolver::readRes()
{
	ifstream ifs(m_info.resultPath());
	if (ifs.fail())		return -1;
	string str_line;
	getline(ifs, str_line);
	m_res.clear();
	m_res.reserve(m_info.outputNum);//Ŀǰ����xfoil���ɵļ������Ƕ�����160�����ݵ�
	while (getline(ifs, str_line))
	{
		stringstream ss(str_line);
		double x = 0;
		double cp = 0;
		ss >> x >> cp;
		m_res.push_back(Field<double>(x, 0, 0, cp));
	}
	return 0;
}