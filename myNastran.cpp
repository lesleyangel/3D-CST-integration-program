#include "myNastran.h"
#include <windows.h>
#include <tlhelp32.h>
#include <comdef.h>
int find_process(string process_name)
{
	int count = 0;//���̼��� 
	PROCESSENTRY32 pe32;

	pe32.dwSize = sizeof(PROCESSENTRY32);
	HANDLE process_snapshot_handle = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);//�������̿��վ��

	if (process_snapshot_handle == INVALID_HANDLE_VALUE) return -1;//�������ʧ��

	bool is_exist = Process32First(process_snapshot_handle, &pe32);//�ҵ�һ��
	while (is_exist)
	{
		if (!_stricmp(process_name.c_str(), _bstr_t(pe32.szExeFile))) count++;//�����������ִ�Сд
		is_exist = Process32Next(process_snapshot_handle, &pe32);//����һ��
	}
	return count;
}

void myNastran::CalcFilePath()
{
	logPath = BDFPath.substr(0, BDFPath.find(".")) + ".log";
	xdbPath = BDFPath.substr(0, BDFPath.find(".")) + ".xdb";
	pchPath = BDFPath.substr(0, BDFPath.find(".")) + ".pch";
	if (ResPath.size() == 0)
	{
		ResPath = BDFPath.substr(0, BDFPath.find_last_of("/"));
	}
}

int myNastran::NastranCalc()
{
	cout << "Nastran���ڵ��á�������" << endl;
	string command = NasPath + " " + BDFPath + " out=" + ResPath;

	remove(logPath.c_str());//���֮ǰ��ͬ��log�ļ�������ɾ��ԭ�ȵ�log�ļ�
	remove(xdbPath.c_str());//���֮ǰ��ͬ��xbd�ļ�������ɾ��ԭ�ȵ�log�ļ�
	remove(pchPath.c_str());//���֮ǰ��ͬ��pch�ļ�������ɾ��ԭ�ȵ�log�ļ�
	system(command.c_str());



	//�ж�log�ļ��Ƿ����ɳɹ�
	clock_t brginTime = clock();
	while (true)
	{
		ifstream ifs(logPath);
		if (!ifs.fail())
		{
			//cout << logPath << " �ļ��ѳɹ�����" << endl;
			cout << " | Nastran���������С�������" << endl;;
			break;
		}
		else if (clock() - brginTime > 20000)
		{
			cout << " | δ����ȷ����Nastran������ "<< BDFPath <<" �ļ��Ƿ����*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//�ж�nastran.exe�Ƿ����н���
	while (true)
	{
		int count = 0;
		string process_name = "nastran.exe";

		count = find_process(process_name);
		if (count == 0)
		{
			cout << " | nastran.exe ���н���" << endl;
			break;
		}
	}
	//��ȡlog�ļ��жϷ����Ƿ�����
	bool isFinishCalc = false;
	clock_t brginTime1 = clock();
	while (true)
	{
		ifstream ifs(logPath);
		string str_line;
		while (getline(ifs, str_line))
		{
			string::size_type idx = str_line.find("Analysis complete");
			if (idx != string::npos)//���ּ�����ɵı�ʶ��----����ѭ����׼����ȡ����
			{
				isFinishCalc = true;

				ifs.close();
				break;
			}

		}
		if (isFinishCalc)
		{
			//system("C:\\Users\\yycab\\Desktop\\bdf\\nastranTESTfile\\KillNastran.bat");//dos()ɱ��nastran���� ����д �����ټӰ�
			break;
		}
		else if (clock() - brginTime1 > 8000)
		{
			cout << " | ����8s��û����log�����ɹ���Nastran��������*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//�ж�pch�ļ��Ƿ����ɳɹ�
	clock_t brginTime2 = clock();
	while (true)
	{
		ifstream ifs(pchPath);
		if (!ifs.fail())
		{
			cout << " | " << pchPath << " �ļ��ѳɹ�����" << endl;
			break;
		}
		else if (clock() - brginTime2 > 8000)
		{
			cout << " | ����8s��δ����ȷ���� '.pch'������ " << pchPath << " �ļ��Ƿ����*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//Sleep(2000);//ϵͳ��ͣ����
	cout << "Nastran�����ɹ�" << endl;
	return 0;
}

int myNastran::ReadResPCH(bool isOnlyReadDisp)
{
	//���ԭʼ����
	Displacements.clear();
	ElemStress.clear();
	MaxElemStress = pair<int, double>(0, 0);
	MaxElemStrain = pair<int, double>(0, 0);
	MaxNodeDisp = pair<int, Point>(0, Point());
	//�򿪽���ļ�
	ifstream ifs(pchPath);
	if (ifs.fail())
	{
		cout << "Nastran����ʧ�ܣ�.pch �ļ������ڣ������������" << endl;
		return -1;
	}
	else cout << "���ݽ���ļ��򿪳ɹ��� ���ڷ������ݡ���" << endl;
	int elemType = 0;
	int loadtype = 0;//1��ԪӦ�� 2��ԪӦ�� 3�ڵ�λ��

	string str_line;
	while (getline(ifs, str_line))
	{
		//string::size_type idx = str_line.find("SUBTITLE");
		if (str_line.find("SUBTITLE") != string::npos)//���ֱ�ʶ��----׼����ȡ����
		{
			//getline(ifs, str_line);
			getline(ifs, str_line);//�ն�1��
			getline(ifs, str_line);

			if (str_line.find("STRAINS") != string::npos)			loadtype = 2;//cout << "�ҵ���Ӧ������" << endl;
			else if (str_line.find("STRESSES") != string::npos)		loadtype = 1;//cout << "�ҵ���Ӧ������" << endl;
			else if (str_line.find("DISPLACEMENTS") != string::npos)loadtype = 3;//�ҵ��˽ڵ�λ������

			getline(ifs, str_line);
			getline(ifs, str_line);//�ն�����
			if (loadtype != 3) getline(ifs, str_line);

			if (str_line.find("QUAD4") != string::npos)			elemType = 4;//cout << "�ҵ����ı��ε�Ԫ" << endl;
			else if (str_line.find("TRIA3") != string::npos)	elemType = 3;//cout << "�ҵ��������ε�Ԫ" << endl;
			else if (str_line.find("BAR") != string::npos)		elemType = 2;//cout << "�ҵ���һά����Ԫ" << endl;
		}

		if (loadtype == 1)//Ӧ��
		{
			switch (elemType)
			{
			case 4://�ı��ε�Ԫ
				cout << "���ڶ�ȡ�ı��ε�ԪӦ��������";
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)//����һֱ����TITLE��ֹͣ
				{
					int id = atoi(str_line.substr(0, 17).c_str());//1 
					getline(ifs, str_line);//2
					getline(ifs, str_line);//3
					getline(ifs, str_line); double data = atof(str_line.substr(18, 18).c_str());//4

					getline(ifs, str_line);//5
					getline(ifs, str_line); data = (abs(atof(str_line.substr(54, 18).c_str())) > abs(data)) ? atof(str_line.substr(54, 18).c_str()) : data;//6
					if (abs(data) > MaxElemStress.second)
					{
						MaxElemStress = pair<int, double>(id, data);
					}
					ElemStress.insert(pair<int, double>(id, data));
					for (int i = 0; i < 24; i++)
					{
						getline(ifs, str_line);
					}
				}
				loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				cout << "���" << endl;
				break;
			case 3://�����ε�Ԫ
				break;
			case 2://����Ԫ
				break;
			}
		}
		else if (loadtype == 2)//Ӧ��
		{
			switch (elemType)
			{
			case 4://�ı��ε�Ԫ
				cout << "���ڶ�ȡ�ı��ε�ԪӦ�䡣����";
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)
				{
					int id = atoi(str_line.substr(0, 17).c_str());//1 
					getline(ifs, str_line);//2
					getline(ifs, str_line);//3
					getline(ifs, str_line); double data = atof(str_line.substr(18, 18).c_str());//4

					getline(ifs, str_line);//5
					getline(ifs, str_line); data = (abs(atof(str_line.substr(54, 18).c_str())) > abs(data)) ? atof(str_line.substr(54, 18).c_str()) : data;//6
					
					if (MaxElemStrain.second < abs(data))
					{
						MaxElemStrain = pair<int, double>(id, data);
					}
					ElemStrain.insert(pair<int, double>(id, data));
					for (int i = 0; i < 24; i++)
					{
						getline(ifs, str_line);
					}
				}
				loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				cout << "���" << endl;
				break;
			case 3://�����ε�Ԫ
				break;
			case 2://����Ԫ
				break;
			}
		}
		else if (loadtype == 3)//�ڵ�λ��
		{
			cout << "���ڶ�ȡ�ڵ�λ�ơ�����";
			while (getline(ifs, str_line))//����һֱ����TITLE��ֹͣ
			{
				if (str_line.find("TITLE") != string::npos)
				{
					break;
				}
				Point disp;
				int id = atoi(str_line.substr( 0, 17).c_str());//�ڵ���1
				disp.setX(atof(str_line.substr(23, 18).c_str()));//x����λ��
				disp.setY(atof(str_line.substr(41, 18).c_str()));//y����λ��
				disp.setZ(atof(str_line.substr(59, 18).c_str()));//z����λ��
				if (abs(MaxNodeDisp.second.getY()) < abs(disp.getY()))
				{
					MaxNodeDisp = pair<int, Point>(id, disp);
				}
				Displacements.insert(pair<int, Point>(id, disp));
				getline(ifs, str_line);
			}
			loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
			elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
			cout << "���" << endl;
			if (isOnlyReadDisp)
			{
				break;//��ǰֻŪ���ڵ�λ�ƣ���Ҫ�������������ɾ��
			}
		}
	}
	return 0;
}

int myNastran::ReadResPCH(vector<int> CQUAD4id, vector<int>CBARid)
{
	int cabinID = 0;
	int bulkhID = 0;

	ifstream ifs(pchPath);
	if (ifs.fail())
	{
		cout << "Nastran����ʧ�ܣ�.pch �ļ������ڣ������������" << endl;
		return -1;
	}
	else cout << "���ݽ���ļ��򿪳ɹ��� ���ڷ������ݡ���" << endl;

	int elemType = 0;
	int loadtype = 0;//1��ԪӦ�� 2��ԪӦ��
	//��¼���Ӧ��
	vector<double> maxstress(CQUAD4id.size());//�洢���Ӧ��ֵ
	vector<int> maxstressID(CQUAD4id.size());//�洢���Ӧ����Ԫ���
	vector<double> maxstrain(CQUAD4id.size());//�洢���Ӧ��ֵ
	vector<int> maxstrainID(CQUAD4id.size());//�洢���Ӧ�䵥Ԫ���
	//��¼�������
	vector<double> maxstress2(CBARid.size());//�洢���Ӧ��ֵ
	vector<int> maxstressID2(CBARid.size());//�洢���Ӧ����Ԫ���
	vector<double> maxstrain2(CBARid.size());//�洢���Ӧ��ֵ
	vector<int> maxstrainID2(CBARid.size());//�洢���Ӧ�䵥Ԫ���
	for (size_t i = 0; i < CQUAD4id.size(); i++)
	{
		maxstress[i] = 0;
		maxstressID[i] = 0;
		maxstrain[i] = 0;
		maxstrainID[i] = 0;
	}
	string str_line;
	while (getline(ifs, str_line))
	{
		//string::size_type idx = str_line.find("SUBTITLE");
		if (str_line.find("SUBTITLE") != string::npos)//���ֱ�ʶ��----׼����ȡ����
		{
			//getline(ifs, str_line);
			getline(ifs, str_line);//�ն�1��
			getline(ifs, str_line);

			if (str_line.find("STRAINS") != string::npos)			loadtype = 2;//cout << "�ҵ���Ӧ������" << endl;
			else if (str_line.find("STRESSES") != string::npos)		loadtype = 1;//cout << "�ҵ���Ӧ������" << endl;

			getline(ifs, str_line);
			getline(ifs, str_line);//�ն�����
			getline(ifs, str_line);

			if (str_line.find("QUAD4") != string::npos)			elemType = 4;//cout << "�ҵ����ı��ε�Ԫ" << endl;
			else if (str_line.find("TRIA3") != string::npos)	elemType = 3;//cout << "�ҵ��������ε�Ԫ" << endl;
			else if (str_line.find("BAR") != string::npos)		elemType = 2;//cout << "�ҵ���һά����Ԫ" << endl;
		}

		if (loadtype == 1)//Ӧ��
		{
			switch (elemType)
			{
			case 4://�ı��ε�Ԫ
				cabinID = 0;
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)//����һֱ����TITLE��ֹͣ
				{
					int id = atoi(str_line.substr(0, 17).c_str());//1 
					getline(ifs, str_line);//2
					getline(ifs, str_line);//3
					getline(ifs, str_line); double data = atof(str_line.substr(18, 18).c_str());//4
					//string dawdwadw = str_line.substr(18,35);
					if ((size_t)cabinID + 1 < CQUAD4id.size() && id == CQUAD4id[cabinID + 1])
					{
						cabinID++;
					}

					if (abs(maxstress[cabinID]) < abs(data))
					{
						maxstress[cabinID] = data;
						maxstressID[cabinID] = id;
					}
					for (int i = 0; i < 26; i++)
					{
						getline(ifs, str_line);
					}
				}
				loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				break;
			case 3://�����ε�Ԫ
				break;
			case 2://����Ԫ
				bulkhID = 0;
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)
				{
					double data = 0;
					int id = atoi(str_line.substr(0, 18).c_str());//1
					if ((size_t)bulkhID + 1 < CBARid.size() && id == CBARid[bulkhID + 1] + 1)
					{
						bulkhID++;
					}

					getline(ifs, str_line); double dataA = atof(str_line.substr(54, 18).c_str());//string dsadawd = str_line.substr(54,18);//2
					if (abs(maxstress2[bulkhID] < abs(dataA)))
					{
						maxstress2[bulkhID] = dataA;
						maxstressID2[bulkhID] = id;
					}
					getline(ifs, str_line);//3
					getline(ifs, str_line);//4
					getline(ifs, str_line); double dataB = atof(str_line.substr(18, 18).c_str());//5
					if (abs(maxstress2[bulkhID] < abs(dataB)))
					{
						maxstress2[bulkhID] = dataB;
						maxstressID2[bulkhID] = id;
					}
					getline(ifs, str_line);
					loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
					elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				}
				break;
			}
		}
		else if (loadtype == 2)//Ӧ��
		{
			switch (elemType)
			{
			case 4://�ı��ε�Ԫ
				cabinID = 0;
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)
				{
					int id = atoi(str_line.substr(0, 17).c_str());//1 
					getline(ifs, str_line);//2
					getline(ifs, str_line);//3
					getline(ifs, str_line); double data = atof(str_line.substr(18, 35).c_str());//4

					if ((size_t)cabinID + 1 < CQUAD4id.size() && id == CQUAD4id[cabinID + 1])
					{
						cabinID++;
					}
					if (abs(maxstrain[cabinID]) < abs(data))
					{
						maxstrain[cabinID] = data;
						maxstrainID[cabinID] = id;
					}
					for (int i = 0; i < 26; i++)
					{
						getline(ifs, str_line);
					}
				}
				loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				break;
			case 3://�����ε�Ԫ
				break;
			case 2://����Ԫ
				bulkhID = 0;
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)
				{
					double data = 0;
					int id = atoi(str_line.substr(0, 18).c_str());//1
					if ((size_t)bulkhID + 1 < CBARid.size() && id == CBARid[bulkhID + 1] + 1)
					{
						bulkhID++;
					}
					getline(ifs, str_line); double dataA = atof(str_line.substr(54, 18).c_str());//2
					if (abs(maxstrain2[bulkhID] < abs(dataA)))
					{
						maxstrain2[bulkhID] = dataA;
						maxstrainID2[bulkhID] = id;
					}
					getline(ifs, str_line);//3
					getline(ifs, str_line);//4
					getline(ifs, str_line); double dataB = atof(str_line.substr(18, 18).c_str());//5
					if (abs(maxstrain2[bulkhID] < abs(dataB)))
					{
						maxstrain2[bulkhID] = dataB;
						maxstrainID2[bulkhID] = id;
					}
					getline(ifs, str_line);
					loadtype = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
					elemType = 0;//��Ԫ��ȡ��ϣ�������ͱ���������
				}
				break;
			}
		}
	}
	this->maxstress = maxstress;//�洢���Ӧ��ֵ
	this->maxstressID = maxstressID;//�洢���Ӧ����Ԫ���
	this->maxstrain = maxstrain;//�洢���Ӧ��ֵ
	this->maxstrainID = maxstrainID;//�洢���Ӧ�䵥Ԫ���

	this->maxstress2 = maxstress2;//�洢���Ӧ��ֵ
	this->maxstressID2 = maxstressID2;//�洢���Ӧ����Ԫ���
	this->maxstrain2 = maxstrain2;//�洢���Ӧ��ֵ
	this->maxstrainID2 = maxstrainID2;//�洢���Ӧ�䵥Ԫ���
	cout << "������ɣ������ѱ��浽ָ���ļ���" << endl;
	return 0;
}
