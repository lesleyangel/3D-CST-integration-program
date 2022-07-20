#include"Aircraft.h"
// #include<boost/foreach.hpp>
#include <iomanip>
#include <sstream>
#include <Windows.h>
#include<direct.h>
#include"Bone.h"
#include"myNastran.h"
//using namespace boost::assign;

//��ʼ��������
Aircraft::Aircraft() 
{
	Name = "testFile";
	isTecPlotFile = false;
	isAeroForce = false;
	isLiftDragRatio = false;
	isVolume = false;
	isStructInp = false;
	isBoneInp = false;
	isBoneTec = false;
	isStructBDF = false;
	isStructTec = false;
	isStruct2D_Tec = false;
	isStruct2D_Nas = false;
	isStruct2D_AeroStruct = false;
	isStruct2D_LDratio = false;
	isStruct2D_TotalMass = false;
	isStruct2D_FixedMass = false;
}

//�ṩ�����������ͷŲ���
Aircraft::~Aircraft()
{
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		if ((*it) != nullptr)
		{
			delete (*it);
			(*it) = nullptr;
		}
	}
}

//Ϊ��������Ӳ���
void Aircraft::AddShape(AbstructShape* Shape)
{
	CSTall.push_back(Shape);
}

//�������ɷ�����ģ��
void Aircraft::CalculateAircraft()
{
	int i = 1;
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		cout << "-------------------" << endl << endl;
		cout << "��" << i << "���������ڽ��м���..." << endl;
		(*it)->BuildShape();
		i++;
	}
	
	cout << "-------------------" << endl;
	cout << "Calculate Aircraft Completed!" << endl;

	
}

//���������ݵ���Ϊ�ı���ʽ
void Aircraft::SaveTecPlotAll(string name)
{
	int x = 1;
	ofstream ofs;
	ofs.open(name, ios::trunc);
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		for(CSTsurface i: (*it)->CSTsf)
		{
			toTecplot(i.GridUpp.P, i.GridUpp.E, ofs, "\"part" + to_string(x) + "'s GridUpp\"");
			toTecplot(i.GridLow.P, i.GridLow.E, ofs, "\"part" + to_string(x) + "'s GridLow\"");
		}
		x++;
	}
	ofs.close();
	cout << "-------------------" << endl;
	/*cout << "SaveTecPlotAll Completed!" << endl;*/
	cout << "����TecPlot��������ɹ���" << endl;
	cout << "-------------------" << endl;
}

inline void Aircraft::compareKeyWords(stringstream& ss, string str, bool& keywords)
{
	getline(ss, str, ',');
	if (!str.compare("on"))	keywords = true;
	else keywords = false;
}

void Aircraft::SaveToClac(string path)
{
	int num = 0;
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		num++;
		//(*it)->SaveToAero("C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh" + to_string(num) + ".dat");//����������
		string inpdatfile = path +"mesh"+ to_string(num) + ".dat";
		(*it)->SaveToAero(inpdatfile);//����������
		(*it)->MakeAeroInfo(path,num);
	}
}

void Aircraft::ReadAeroFromFile(string path)
{
	int num = 0;
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		num++;
		//��ȡ����������ļ����������������ڴ���
		ifstream ifs;
		//string path = "C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh";
		ifs.open(path + "mesh" + to_string(num) + "\\node_force.dat", ios::in);
		if (!ifs.is_open()) cout << path + "mesh" + to_string(num) + "\\node_force.dat ��ȡʧ�ܣ�" << endl;
		else cout << path + "mesh" + to_string(num) + "\\node_force.dat ��ȡ�ɹ���" << endl;

		string strLine;
		getline(ifs, strLine);//�ն�����
		getline(ifs, strLine);//�ն�����
		getline(ifs, strLine);//�ն�����
		for (size_t i = 0; i < (*it)->CSTsf.size(); i++)//�ϱ���
		{
			int pNumUpp = (*it)->CSTsf[i].GridUpp.pointNum();
			(*it)->CSTsf[i].GridUpp.nodeForce = zeros(pNumUpp, 3);
			//���ж�ȡ
			for (int j = 0; j < pNumUpp; j++)
			{
				if (getline(ifs, strLine))//�ж��Ƿ����п��Զ�ȡ
				{
					istringstream is(strLine);
					double s;
					is >> s; is >> s; is >> s; is >> s;//���� �ڵ���/CX/CY/CZ
					is >> s; (*it)->CSTsf[i].GridUpp.nodeForce(j, 0) = s;//FX
					is >> s; (*it)->CSTsf[i].GridUpp.nodeForce(j, 1) = s;//FY
					is >> s; (*it)->CSTsf[i].GridUpp.nodeForce(j, 2) = s;//FZ

				}
			}
		}
		for (size_t i = 0; i < (*it)->CSTsf.size(); i++)//�±���
		{
			int pNumLow = (*it)->CSTsf[i].GridLow.pointNum();
			(*it)->CSTsf[i].GridLow.nodeForce = zeros(pNumLow, 3);
			//���ж�ȡ
			for (int j = 0; j < pNumLow; j++)
			{
				if (getline(ifs, strLine))//�ж��Ƿ����п��Զ�ȡ
				{
					istringstream is(strLine);
					double s;
					is >> s; is >> s; is >> s; is >> s;//���� �ڵ���/CX/CY/CZ
					is >> s; (*it)->CSTsf[i].GridLow.nodeForce(j, 0) = s;//FX
					is >> s; (*it)->CSTsf[i].GridLow.nodeForce(j, 1) = s;//FY
					is >> s; (*it)->CSTsf[i].GridLow.nodeForce(j, 2) = s;//FZ
				}
			}
		}
	}
}

void Aircraft::CalcAero(string path)
{
	cout << "$----------------------------------------$" << endl;
	cout << "$             ���������㿪ʼ             $" << endl;
	cout << "$----------------------------------------$" << endl;
	for (size_t i = 0; i < CSTall.size(); i++)
	{
		cout << "~~~~~~~~��" << i + 1 << "����������������������~~~~~~~~" << endl;
		string command = "aero_calc.exe " + path + "mesh" + to_string(i+1) + "_cmd.txt";
		system(command.c_str());
		cout << endl;
	}
	cout << "$----------------------------------------$" << endl;
	cout << "$             �������������             $" << endl;
	cout << "$----------------------------------------$" << endl;
}

int Aircraft::AeroStrcStruct2D(int XSnum, int ZSnum, CSTsurface& surf, double& Stress, double& Disp, double& Twist)
{
	cout << "���ڽ��нṹ������ XSnum = " << XSnum << "  ZSnum = " << ZSnum << endl;
	Struct2D s2d;
	const int Xnum = surf.NEta;
	const int Znum = surf.NFaiL;
	ZSnum = ZSnum + 2;
	//��һ��
	int Xdelta = Xnum / (XSnum - 1);
	int Zdelta = Znum / (ZSnum - 1);
	mat Xsite = linspace(0, (double)Xdelta * ((double)XSnum - 1), XSnum);
	mat Zsite = linspace(Zdelta, (double)Zdelta * ((double)ZSnum - 2), ZSnum - 2);
	//
	StructPart SP(surf.node_2D);
	SP.setProperty(m_property);
	SP.setIsFixedMass(this->isStruct2D_FixedMass);
	SP.setSite(Xsite, Zsite);
	SP.calcStructPart(surf.Origin, surf.Rotation);
	s2d.PartList.push_back(SP);
	int state = s2d.AeroelasticAnalysis(docPath + Name, exePath);
	Stress = s2d.PartList[0].MaxElemStress.second;
	Disp = s2d.PartList[0].MaxNodeDisp.second.getY();
	Twist = s2d.PartList[0].MaxTwisting;
	return state;
}

void Aircraft::GetBone()
{
	cout << "-------------------" << endl << endl;
	int i = 1;
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{
		
		cout << "��" << i << "��������������ṹ����..." << endl;
		if ((*it)->ShapeType == 0)
		{
			bone.MainBone = (*it)->BuildBone()[0];//ͷ������
		}
		else
		{
			vector<SingleBone> SB = (*it)->BuildBone();
			for (size_t num = 0; num < SB.size(); num++)
			{
				bone.BrachBone.push_back(SB[num]);
			}
		}
		i++;
	}
	cout << "�������Ӹ������ṹ���񡣡���" << endl;
	bone.CalcBone();
	
	cout << "�ṹ�������ɳɹ�����" << endl;
	cout << "-------------------" << endl;
}

void Aircraft::calcAeroForce(string path)
{
	SaveToClac(path);						//����aero_xxxx���õ������ļ�
	CalcAero(path);							//����aero_xxxx����������
	ReadAeroFromFile(path);					//�ѽڵ�����������Ϣ�����ڴ�
}

double Aircraft::calcK(string savepath)
{
	double K = 0;//��ʼ������ֵ-�����
	double lift = 0, drag = 0;//��ʼ����������������
	for (size_t i = 0; i < CSTall.size(); i++)
	{
		for (size_t j = 0; j < CSTall[i]->CSTsf.size(); j++)
		{
			for (uword num = 0; num < CSTall[i]->CSTsf[j].GridUpp.nodeForce.n_rows; num++)
			{
				lift += CSTall[i]->CSTsf[j].GridUpp.nodeForce(num, 1);
				drag += CSTall[i]->CSTsf[j].GridUpp.nodeForce(num, 0);
			}
			for (uword num = 0; num < CSTall[i]->CSTsf[j].GridLow.nodeForce.n_rows; num++)
			{				
				lift += CSTall[i]->CSTsf[j].GridLow.nodeForce(num, 1);
				drag += CSTall[i]->CSTsf[j].GridLow.nodeForce(num, 0);
			}
		}
	}
	K = lift / drag;
	if (savepath.compare("none"))//������·�����ǿ�
	{
		ofstream ofs;
		ofs.open(savepath);
		ofs << K;
		ofs.close();
	}
	cout << "$---------  ����� = " << K << "   -----------$" << endl;
	return K;
}

double Aircraft::getVol(string savepath/* = "none"*/)
{
	double Vol = 0;
	for (size_t i = 0; i < CSTall.size(); i++)
	{
		for (size_t j = 0; j < CSTall[i]->CSTsf.size(); j++)
		{
			Vol += CSTall[i]->CSTsf[j].Info.Volume;
		}
	}
	if (savepath.compare("none"))//������·�����ǿ�
	{
		ofstream ofs;
		ofs.open(savepath);
		ofs << Vol;
		ofs.close();
	}
	return Vol;
}

void Aircraft::CalcStruct2D()
{
	Struct2D s2d;
	for (size_t n = 0; n < CSTall.size(); n++)
	{
		for (size_t i = 0; i < CSTall[n]->CSTsf.size(); i++)
		{
			CSTsurface &sf = CSTall[n]->CSTsf[i];
			sf.MakeStruct2D_byNum();
			//
			// int Xnum = sf.NEta;
			// int Znum = sf.NFaiL;
			// int XSnum = sf.Struct2dXnum;
			// int ZSnum = sf.Struct2dZnum + 2;
			// int Xdelta = Xnum / (XSnum - 1);
			// int Zdelta = Znum / (ZSnum - 1);
			// mat Xsite = linspace(0, (double)Xdelta * ((double)XSnum - 1), XSnum);
			// mat Zsite = linspace(Zdelta, (double)Zdelta * ((double)ZSnum - 2), ZSnum - 2);
			//

			
			const double eta_block = 1.0 / sf.eta_block_size;
			const int eta_struct_num = eta_block - (int)eta_block < 1e-5 ? (int)eta_block +1: (int)eta_block + 2;
			const double eta_T_ratio = (1.0 - (eta_struct_num - 2.0) * sf.eta_block_size) / (1.0 / (eta_struct_num - 1));

			const double fai_block = 1.0 / sf.fai_l_block_size;
			const int fai_struct_num = fai_block - (int)fai_block < 1e-5 ? (int)fai_block -1: (int)fai_block;
			const double fai_T_ratio = (1.0 - fai_struct_num * sf.fai_l_block_size) / (1.0 / (fai_struct_num + 1));

			if (eta_struct_num < 2)
			{
				cout << "���������������" << endl;
			}
			if (fai_struct_num < 0)
			{
				cout << "���������������" << endl;
			}
			
			vec x_site = zeros(eta_struct_num);
			for (int id = 0; id < eta_struct_num - 1; id++)
			{
				x_site(id) = id * sf.n_eta;
			}
			x_site(eta_struct_num - 1) = sf.NEta - 1;
			vec z_site = zeros(fai_struct_num);
			for (int id = 0; id < fai_struct_num; id++)
			{
				z_site(id) = (id+1) * sf.n_fai_l;
			}
			//
			StructPart SP(sf.node_2D);//����������Ϣ
			SP.setProperty(m_property);//���ò�������
			SP.setIsFixedMass(this->isStruct2D_FixedMass);//�����Ƿ�̶�����
			SP.x_T_ratio = eta_T_ratio; // ����ĩ���ṹ�������ߵĲ���ϵ��
			SP.z_T_ratio = fai_T_ratio; // ����ĩ���ṹ�������ߵĲ���ϵ��
			SP.setSite(x_site, z_site);//����eta/fai������������
			SP.calcStructPart(sf.Origin, sf.Rotation);//��������ƽ�ƺ���ת

			s2d.PartList.push_back(SP);
		}
	}
	struct2d = s2d;
	cout << "struct2D�������ɹ���" << endl;
	cout << "-------------------" << endl;
}

void Aircraft::FindOptimalStruct2D()
{
	//Struct2D s2d;
	/*for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)*/
	vector<AbstructShape*>::iterator it = CSTall.begin();
	{
		/*for (int i = 0; i < (*it)->CSTsf.size(); i++)*/
		int i = 0;
		{
			(*it)->CSTsf[i].MakeStruct2D_byNum();

			const int minXnum = 5;	//������С����
			const int maxXnum = 25;	//����������
			const int minZnum = 1;	//������С����
			const int maxZnum = 5;	//����������
			const int Xsize = maxXnum + 1;
			const int Zsize = maxZnum + 1;
			vector<vector<double>> StressXZ(Xsize);
			vector<vector<double>> DispXZ(Xsize);
			vector<vector<double>> TwistXZ(Xsize);
			for (size_t xid = 0; xid < StressXZ.size(); xid++)
			{
				StressXZ[xid].resize(Zsize);
				DispXZ[xid].resize(Zsize);
				TwistXZ[xid].resize(Zsize);
				for (size_t zid = 0; zid < StressXZ[xid].size(); zid++)
				{
					StressXZ[xid][zid] = -1;
					DispXZ[xid][zid] = -1;
					TwistXZ[xid][zid] = -1;
				}
			}
			double minStress = 3.5e8;
			const double upboundStress = 3.5e8;
			const double upboundDisp = 3.0;
			const double upboundTwist = 0.1;
			int nowXSnum = (*it)->CSTsf[i].Struct2dXnum;
			int nowZSnum = (*it)->CSTsf[i].Struct2dZnum;

			//��һ��
			int X = (*it)->CSTsf[i].Struct2dXnum;//���߳�ʼ����
			int Z = (*it)->CSTsf[i].Struct2dZnum;//������ʼ����
			AeroStrcStruct2D(X, Z, (*it)->CSTsf[i], StressXZ[X][Z], DispXZ[X][Z], TwistXZ[X][Z]);
			minStress = StressXZ[X][Z];

			vector<pair<int, int>> recordXZ;
			while (true)
			{
				recordXZ.push_back(pair<int, int>(nowXSnum, nowZSnum));
				int tempX = nowXSnum; int tempZ = nowZSnum;

				vector<pair<int, int>> XZlist(4);
				XZlist[0] = pair<int, int>(nowXSnum + 1, nowZSnum);
				XZlist[1] = pair<int, int>(nowXSnum - 1, nowZSnum);
				XZlist[2] = pair<int, int>(nowXSnum, nowZSnum + 1);
				XZlist[3] = pair<int, int>(nowXSnum, nowZSnum - 1);
				for (size_t id = 0; id < XZlist.size(); id++)
				{
					X = XZlist[id].first; Z = XZlist[id].second;
					if (minXnum <= X && X <= maxXnum && minZnum <= Z && Z <= maxZnum)
					{
						int calcState = 0;
						if (-2 < StressXZ[X][Z] && StressXZ[X][Z] < 0)//δ����״̬-1
						{
							calcState = AeroStrcStruct2D(X, Z, (*it)->CSTsf[i], StressXZ[X][Z], DispXZ[X][Z], TwistXZ[X][Z]);
							if (calcState < 0)//������� ����ֵΪ-3
								StressXZ[X][Z] = DispXZ[X][Z] = TwistXZ[X][Z] = -3;
							else
							{
								if (StressXZ[X][Z] < upboundStress && DispXZ[X][Z] < upboundDisp && TwistXZ[X][Z] < upboundTwist)//����Լ��
									if (StressXZ[X][Z] < minStress) { minStress = StressXZ[X][Z]; tempX = X; tempZ = Z; }//������СӦ���Ͷ�Ӧ�Ĳ���
							}
						}
					}
				}

				if (tempX != nowXSnum || tempZ != nowZSnum)//������ԭ�Ȳ�ͬ
				{
					nowXSnum = tempX; nowZSnum = tempZ;
				}
				else
				{
					cout << "��������������������" << endl;
					break;
				}
			}
			
		}
	}

}

//void Aircraft::SaveStruct2DAsTecplot(string savepath)
//{
//	int x = 1;
//	ofstream ofs;
//	ofs.open(savepath, ios::trunc);
//	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
//	{
//		BOOST_FOREACH(CSTsurface i, (*it)->CSTsf)
//		{
//			if (i.StructPE.P.n_rows != 0)
//			{
//				ofs << "Title = \"part" << x << "'s Struct\"" << endl;
//				ofs << "variables = \"x\",\"y\",\"z\"" << endl;
//				ofs << "zone N=" << i.StructPE.P.n_rows << ",E=" << i.StructPE.E.n_rows << ",datapacking=block" << endl;
//				ofs << "zonetype=fequadrilateral" << endl;
//				ofs << trans(i.StructPE.P) << endl;
//				ofs << i.StructPE.E + 1 << endl;
//			}
//		}
//		x++;
//	}
//	ofs.close();
//
//	cout << "struct2D�����ѳɹ�����TecPlot�ļ���" << endl;
//	cout << "-------------------" << endl;
//}

void Aircraft::CalcStruct1D()
{
	mat X = linspace(0, 1, 5);
	mat Z = linspace(0, 1, 3);
	cout << "-------------------" << endl ;
	int i = 1;
	
	for (vector<AbstructShape*>::iterator it = CSTall.begin(); it != CSTall.end(); it++)
	{

		cout << "��" << i << "����������������ǽṹ����..." << endl;
		if ((*it)->ShapeType == 0)
		{
			struct1d.MainBone = (*it)->BuildBone()[0];//ͷ������
		}
		else
		{
			vector<ShellandBeam> SB = (*it)->BuildStrc1D();
			for (size_t num = 0; num < SB.size(); num++)
			{
				struct1d.BarchStruc.push_back(SB[num]);
			}
		}
		i++; 
	}
	cout << "�������崦���������񡣡���" << endl;
	struct1d.getConnPoint();

	cout << "�����������ɳɹ�����" << endl;
	cout << "-------------------" << endl;
	
}

int Aircraft::RunFromFile(string path)
{
	cout << endl;
	cout << "�����ļ���" << path << endl;
	ifstream ifs(path);
	if (ifs.fail())
	{
		cout << "�ļ������ڣ�" << endl;
		return -1;
	}
	else
	{
		cout << "�򿪳ɹ���" << endl;
	}
	
	
	docPath = path.substr(0, path.find_last_of("/"));//��ȡ�ļ���·��
	if (docPath == path)//���·��
	{
		docPath = path.substr(0, path.find_last_of("\\"));//��ȡ�ļ���·��
		if (docPath == path)//���·��
		{
			docPath = "";
		}
		//docPath = "";
	}

	string str_line;
	while (getline(ifs, str_line))
	{
		if (str_line[0] == '$') {}//�ж��׸��ַ�,ע����
		else if (str_line[0] == '*')//�ж��׸��ַ�,������벿��������
		{
			str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//����ո���
			string partName = str_line.substr(1, str_line.length() - 1);//

			if (!partName.compare("MainBody1"))		 {	SMB1 mb1;	mb1.readFromFile(ifs);	AddShape(new ShapeMainBody1(mb1));	}
			else if (!partName.compare("MainBody2")) {	SMB2 mb2;	mb2.readFromFile(ifs);	AddShape(new ShapeMainBody2(mb2));	}
			else if (!partName.compare("MainBody3")) {	SMB3 mb3;	mb3.readFromFile(ifs);	AddShape(new ShapeMainBody3(mb3));	}
			else if (!partName.compare("ShapeWing1")){	SW1 sw1;	sw1.readFromFile(ifs);	AddShape(new ShapeWing1(sw1));		}
			else if (!partName.compare("ShapeWing2")){	SW2 sw2;	sw2.readFromFile(ifs);	AddShape(new ShapeWing2(sw2));		}
			else if (!partName.compare("ShapeTail1")){	ST1 st1;	st1.readFromFile(ifs);	AddShape(new ShapeTail1(st1));		}
			else if (!partName.compare("ShapeTail2")){	ST2 st2;	st2.readFromFile(ifs);	AddShape(new ShapeTail2(st2));		}
			else if (!partName.compare("Property"))	 {	m_property.readFromFile(ifs);		}
		}
		else if (!str_line.compare("BEGIN INPUT"))//�ж��ַ���,�������������
		{
			while (getline(ifs, str_line) && str_line.compare("END INPUT"))//����"END INPUT"ʱֹͣ�ò���
			{
				str_line.erase(remove_if(str_line.begin(), str_line.end(), isspace), str_line.end());//��������ո���
				stringstream ss(str_line);
				string str_tmp;//�洢�·��ж���ʹ�õ�һ������
				
				getline(ss, str_line, ',');//��ȡ�׸�string
				if (str_line[0] == '$') {}//����
				else if (!str_line.compare("NAME"))					getline(ss, Name, ',');
				else if (!str_line.compare("TecPlotFile"))			compareKeyWords(ss, str_tmp, isTecPlotFile);
				else if (!str_line.compare("AeroForce"))			compareKeyWords(ss, str_tmp, isAeroForce);
				else if (!str_line.compare("LiftDragRatio"))		compareKeyWords(ss, str_tmp, isLiftDragRatio);
				else if (!str_line.compare("Volume"))				compareKeyWords(ss, str_tmp, isVolume);
				else if (!str_line.compare("StructInp"))			compareKeyWords(ss, str_tmp, isStructInp);
				else if (!str_line.compare("BoneInp"))				compareKeyWords(ss, str_tmp, isBoneInp);
				else if (!str_line.compare("BoneTec"))				compareKeyWords(ss, str_tmp, isBoneTec);
				else if (!str_line.compare("StructBDF"))			compareKeyWords(ss, str_tmp, isStructBDF);
				else if (!str_line.compare("StructTec"))			compareKeyWords(ss, str_tmp, isStructTec);
				else if (!str_line.compare("Struct2D_Tec"))			compareKeyWords(ss, str_tmp, isStruct2D_Tec);
				else if (!str_line.compare("Struct2D_Nas"))			compareKeyWords(ss, str_tmp, isStruct2D_Nas);
				else if (!str_line.compare("Struct2D_AeroStruct"))	compareKeyWords(ss, str_tmp, isStruct2D_AeroStruct);
				else if (!str_line.compare("Struct2D_LDratio"))		compareKeyWords(ss, str_tmp, isStruct2D_LDratio);
				else if (!str_line.compare("Struct2D_TotalMass"))	compareKeyWords(ss, str_tmp, isStruct2D_TotalMass);
				else if (!str_line.compare("Struct2D_FixedMass"))	compareKeyWords(ss, str_tmp, isStruct2D_FixedMass);
			}
		}
		
	}
	ifs.close(); 
	//-------------------------��ʼ����--------------------------------
	if (docPath.size() != 0)
	{
		docPath += "/";
	}
	CalculateAircraft();//���ɷ���������
	if (isTecPlotFile)
	{
		SaveTecPlotAll(docPath  + Name + ".dat");//+ "/"
	}
	if (isVolume)
	{
		getVol(docPath + Name + "_Vol.txt");// + "/"
	}
	if (isAeroForce || isLiftDragRatio || isStructBDF || isStructInp || isBoneInp)
	{
		string aeropath = docPath + Name + "_aero\\";
		int a = _mkdir(aeropath.c_str());
		calcAeroForce(aeropath);//����aero_calc����������
	}
	if (isLiftDragRatio)
	{
		calcK(docPath + Name + "_K.txt");// + "/"
	}
	//------------------------StructBone-------------------------------
	if (isBoneInp || isBoneTec)
	{
		GetBone();//���ɹǼ���ģ��
	}
	if (isBoneInp)
	{
		bone.SaveAsAbaqus(docPath + Name + "_bone.inp");// + "/"
	}
	if (isBoneTec)
	{
		bone.SaveAsTcp(docPath + Name + "_bone.dat");// + "/"
	}
	//-------------------------Struct1D--------------------------------
	if (isStructBDF || isStructInp || isStructTec)
	{
		CalcStruct1D();//�������ǽ�Ͻṹ����
	}
	if (isStructTec)
	{
		struct1d.SaveAsTecplt(docPath + Name);// + "/"
	}
	if (isStructBDF)
	{
		struct1d.SaveAsNastran(docPath + Name + "_strc1d");// + "/"
	}
	if (isStructInp)
	{
		struct1d.SaveAsAbaqus(docPath + Name + "_strc1d.inp");// + "/"
	}
	//-------------------------Struct2D--------------------------------
	if (isStruct2D_Tec || isStruct2D_AeroStruct || isStruct2D_LDratio || isStruct2D_TotalMass || isStruct2D_Nas)
	{
		CalcStruct2D();
	}
	if (isStruct2D_Tec)
	{
		struct2d.SaveAsTecplot(docPath + Name);
	}
	if (isStruct2D_Nas)
	{
		struct2d.SaveAsNastran(docPath + Name, 101);
	}
	if (isStruct2D_TotalMass)
	{
		double mass = struct2d.CalcTotalMass(docPath + Name + "_Mass.txt");
	}
	if (isStruct2D_LDratio)
	{
		struct2d.AeroAnalysis(docPath + Name, exePath);
	}
	if (isStruct2D_AeroStruct)
	{
		struct2d.AeroelasticAnalysis(docPath + Name, exePath);
	}
	if (isStruct2D_FixedMass)
	{
		FindOptimalStruct2D();
	}
	return 0;
}
