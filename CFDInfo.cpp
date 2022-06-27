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

	ifViscous = 0;//是否启用粘性计算
	HeadRadius = 0.03;//头部半径
	WallTemperature = 300;//壁面温度
	WallEmissivity = 0.7;//壁面辐射率
	AirType = 1;//气体模型（0，1）
	HeatTransferCoefficient = 0;//换热系数
	MediumTemperature = 280;//介质系数
	InnerWallTemperature = 350;//内壁面温度
}

int AeroCalc::CalcAero()
{
	cout << "$----------------------------------------$" << endl;
	cout << "$             气动力计算开始             $" << endl;
	cout << "$-----                              -----$" << endl;
	//for (int i = 0; i < CSTall.size(); i++)
	int i = 0;
	{
		cout << "~~~~~~~~第" << i + 1 << "个气动部件的气动力计算~~~~~~~~" << endl;
		string command = "aero_calc.exe " + CommandfileName;
		system(command.c_str());
		cout << endl;
	}
	cout << "$-----                              -----$" << endl;
	cout << "$             气动力计算结束             $" << endl;
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
		cout << "aerocalc结果文件不存在！" << endl;
		return -1;
	}
	else
	{
		cout << "打开成功！" << endl;
	}


	string docPath = node_force_filepath.substr(0, node_force_filepath.find_last_of("/"));//获取文件夹路径
	if (docPath == node_force_filepath)//相对路径
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
	ofs << "//////////////////计算状态(数组）///////////////////////" << endl;
	ofs << "mach	nAlpha   nBeta  nH x轴方向(1为头部指向尾部，-1为反方向)" << endl;
	ofs << left << setw(8) << ai.MaNum()
		<< left << setw(9) << ai.AlphaNum()
		<< left << setw(7) << ai.BetaNum()
		<< left << setw(3) << ai.HmNum() << "  1" << endl;

	ofs << "mach    //马赫数" << endl;
	for (int i = 0; i < ai.MaNum(); i++)		ofs << ai.Ma[i] << " " << endl;

	ofs << "alpha	//攻角(deg)" << endl;
	for (int i = 0; i < ai.AlphaNum(); i++)		ofs << ai.Alpha[i] << " " << endl;

	ofs << "beta 	//侧滑角(deg)" << endl;
	for (int i = 0; i < ai.BetaNum(); i++)		ofs << ai.Beta[i] << " " << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////////飞行高度(m)////////////////////////" << endl;
	for (int i = 0; i < ai.HmNum(); i++)
		ofs << ai.Hm[i] << " " << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////////分区类型///////////////////////////" << endl;
	ofs << "计算分区数 " << ai.CalcBlockIDNum() << endl;
	ofs << "//////计算分区ID" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockID[i];	ofs << endl;
	ofs << "//////计算分区区域" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockArea[i];	ofs << endl;
	ofs << "//////计算分区类型" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockType[i];	ofs << endl;
	ofs << "///////分区外法线方向" << endl;
	for (int i = 0; i < ai.CalcBlockIDNum(); i++)		ofs << left << setw(8) << ai.CalcBlockOrient[i];	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "//////////////////////控制面////////////////////////////" << endl;
	ofs << "控制面总数 0" << endl;
	ofs << "/////每个控制面设置占一行，分区请与分区设置对应，固定旋转轴(0无，1X,2Y,3Z)，坐标之间不要留空格，由逗号隔开，旋转轴由P2->P1，逆时针旋转为正，单位m" << endl;
	ofs << "^^^对应分区^^^^^^^旋转轴(点1坐标)^^^^^^^旋转轴(点2坐标)^^^^^^^固定旋转轴^^^舵偏角数^^^舵偏角数组^^^^^^^^^^^" << endl;
	ofs << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////sRef参考面积//cRef参考长度//参考尺寸(m)////////" << endl;
	ofs << "1 1 1 0 0" << endl;
	ofs << "^^^^^^^^^^^" << endl;

	ofs << "/////////////////是否考虑粘性(是为1，否为0)/////////////" << endl;
	ofs << "启用粘性计算(0,1)     " << ai.ifViscous << endl;//
	ofs << "头部半径(m)           " << ai.HeadRadius << endl;//
	ofs << "壁面温度(K)           " << ai.WallTemperature << endl;//
	ofs << "壁面辐射率            " << ai.WallEmissivity << endl;//
	ofs << "气体模型(0,1)         " << ai.AirType << endl;//
	ofs << "换热系数              " << ai.HeatTransferCoefficient << endl;//
	ofs << "介质温度(K)           " << ai.MediumTemperature << endl;//
	ofs << "内壁面温度(K)         " << ai.InnerWallTemperature << endl;//
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "////////////////非定常计算选项(是为1，否为0)////////////" << endl;
	ofs << "非定常开启   计算参数:时间间隔(t)    节点速度rpt文档 " << endl;
	ofs << "0            0.0005                  velocity.rpt" << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "//////////////////////相关文件//////////////////////////" << endl;
	ofs << "//网格文件单位(m/mm):		mm" << endl;
	ofs << "//网格类型                  CATIA" << endl;
	ofs << "网格文件：                  " << ai.meshPath << endl;
	ofs << "输出路径：                  " << ai.outputPath << endl;
	ofs << "输出结果文件名：            " << ai.outputName << endl;
	ofs << "^^^^^^^^^^^" << endl;
	ofs << "是否输出压力分布文件plot.plt(plot3D格式):     YES" << endl;
	ofs << "是否输出面元力分布文件face.dat:               NO" << endl;
	ofs << "是否输出节点力分布文件node.dat:               YES" << endl;
	ofs << "是否输出温度分布文件temperature.plt(需要粘性计算) NO" << endl;
	ofs << "是否输出热载分布文件qload.plt(需要粘性计算)       NO" << endl;
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

	//input文件
	string inputPath = aeropath + "/AVL_Input.avl";
	ofstream ofs1;
	ofs1.open(inputPath, ios::trunc);
	ofs1.setf(ios::left);
	ofs1 << partname << endl;
	ofs1 << setw(16) << ai.Ma << setw(16) << " Mach" << endl;
	ofs1 << setw(8) << ai.IYsym << setw(8) << ai.IZsym << setw(8) << ai.Zsym << setw(8) << " iYsym  iZsym  Zsym" << endl;
	ofs1 << setw(8) << ai.Sref << setw(8) << ai.Cref << setw(8) << ai.Bref << setw(8) << " Sref   Cref   Bref" << endl;
	ofs1 << setw(8) << 0 << setw(8) << 0 << setw(8) << 0 << setw(8) << " Xref   Yref   Zref" << endl;
	ofs1 << "0.00                   CDp " << endl;//零升阻力修正系数 默认为0吧
	ofs1 << "SURFACE" << endl;
	ofs1 << "WING-body" << endl;
	ofs1 << setw(8) << ai.Nchordwise << setw(8) << "1.0" << setw(8) << ai.Nspanwise << setw(8) << "0" << endl;
	//ofs1 << "YDUPLICATE" << endl;//是否沿着xz面对称 由此项则代表对机翼进行对称处理
	//ofs1 << "0.0" << endl;
	ofs1 << "SCALE" << endl;//尺寸的缩放
	ofs1 << setw(16) << ai.Scale.getX() << setw(16) << ai.Scale.getY() << setw(16) << ai.Scale.getZ() << endl;
	ofs1 << "TRANSLATE" << endl;//位置的平移
	ofs1 << setw(16) << ai.Orgin.getX() << setw(16) << ai.Orgin.getY() << setw(16) << ai.Orgin.getZ() << endl;
	ofs1 << "ANGLE" << endl;//所有截面角度的改变
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
	//删除已有的结果文件 防止生成错误读取到原来的
	ResFilePath = aeropath + "/ef.dat";
	remove(ResFilePath.c_str());
	//命令行文件
	ofstream ofs2;
	CommandfileName = aeropath + "/AVL_cmd.cmd";
	
	ofs2.open(CommandfileName, ios::trunc);
	ofs2 << "load " << inputPath << endl;	//载入气动外形数据
	ofs2 << "case " << inputRunPath << endl;//载入工况参数设置
	ofs2 << "oper" << endl;					//操作-关键字
	ofs2 << "x" << endl;					//计算气动力
	ofs2 << "fe" << endl;					//单元节点力结果输出
	ofs2 << ResFilePath << endl;			//节点力输出文件路径
	ofs2 << "o" << endl;					//覆盖原始文件
	ofs2 << endl;							//确认
	ofs2 << "quit" << endl;					//退出程序
	ofs2.close();
	return 0;
}

//通过命令行调用气动求解程序
int AVLServer::CalcAero()
{
	cout << "$----------------------------------------$" << endl;
	cout << "$             气动力计算开始             $" << endl;
	cout << "$-----                              -----$" << endl;
	//for (int i = 0; i < CSTall.size(); i++)
	int i = 0;
	{
		cout << "~~~~~~~~第" << i + 1 << "个气动部件的气动力计算~~~~~~~~" << endl;
		string command = AVLexePath + " <" + CommandfileName;// +" exit";
		system(command.c_str());
		cout << endl;
	}
	cout << "$-----                              -----$" << endl;
	cout << "$             气动力计算结束             $" << endl;
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
	else 				cout << "气动结果文件：" << ResFilePath << "打开成功" << endl;

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
			getline(ifs, str_line); //第一行内容
			res_surf.Chordwise = atoi(str_line.substr(19, 6).c_str());
			res_surf.Spanwise = atoi(str_line.substr(37, 6).c_str());
			res_surf.FirstStrip = atoi(str_line.substr(57).c_str());
			getline(ifs, str_line);//第二行内容
			res_surf.Surface_area = atof(str_line.substr(19, 19).c_str());
			for (int i = 0; i < 6; i++)//空读六行
			{
				getline(ifs, str_line);
			}
			getline(ifs, str_line);//第九行内容
			res_surf.CLsurf = atof(str_line.substr(14, 15).c_str());
			res_surf.CDsurf = atof(str_line.substr(38).c_str());
			//getline(ifs, str_line);//第十行内容
			//getline(ifs, str_line);//第十一行内容
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
					getline(ifs, str_line);//第一行内容
					res_strip.Orgin.setX(atof(str_line.substr(9, 14).c_str()));
					res_strip.Chord = atof(str_line.substr(37, 13).c_str());
					res_strip.Incidence = atof(str_line.substr(62, 11).c_str());
					getline(ifs, str_line);//第二行内容
					res_strip.Orgin.setY(atof(str_line.substr(9, 14).c_str()));
					res_strip.StripWidth = atof(str_line.substr(37, 13).c_str());
					res_strip.StripArea = atof(str_line.substr(62).c_str());
					getline(ifs, str_line);//第三行内容
					res_strip.Orgin.setZ(atof(str_line.substr(9, 14).c_str()));
					getline(ifs, str_line);//第四行内容
					getline(ifs, str_line);//第五行内容
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
						getline(ifs, str_line);//读取一行内容
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
					if (StripNum == res_surf.Spanwise)//条带数量到达指定个数时 跳出循环
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
	cout << "XFoil 正在计算...";
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
			cout << "XFoil计算失败！未生成文件："+m_info.resultPath() << endl;
			return -1;
		}
	}
	
	//
	cout << "结束" << endl;
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
}//打印翼型文件

int XFoilSolver::printInput()
{
	ofstream ofs;
	ofs.open(m_info.inputPath(), ios::trunc);
	//
	ofs << std::left
		<< "PLOP\n"		//设置不显示图形窗口
		<< "G\n\n"		//设置不显示图形窗口
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
	m_res.reserve(m_info.outputNum);//目前看来xfoil生成的计算结果是定死的160个数据点
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