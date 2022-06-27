#include "myNastran.h"
#include <windows.h>
#include <tlhelp32.h>
#include <comdef.h>
int find_process(string process_name)
{
	int count = 0;//进程计数 
	PROCESSENTRY32 pe32;

	pe32.dwSize = sizeof(PROCESSENTRY32);
	HANDLE process_snapshot_handle = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);//创建进程快照句柄

	if (process_snapshot_handle == INVALID_HANDLE_VALUE) return -1;//创建句柄失败

	bool is_exist = Process32First(process_snapshot_handle, &pe32);//找第一个
	while (is_exist)
	{
		if (!_stricmp(process_name.c_str(), _bstr_t(pe32.szExeFile))) count++;//进程名不区分大小写
		is_exist = Process32Next(process_snapshot_handle, &pe32);//找下一个
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
	cout << "Nastran正在调用。。。。" << endl;
	string command = NasPath + " " + BDFPath + " out=" + ResPath;

	remove(logPath.c_str());//如果之前有同名log文件，则先删除原先的log文件
	remove(xdbPath.c_str());//如果之前有同名xbd文件，则先删除原先的log文件
	remove(pchPath.c_str());//如果之前有同名pch文件，则先删除原先的log文件
	system(command.c_str());



	//判断log文件是否生成成功
	clock_t brginTime = clock();
	while (true)
	{
		ifstream ifs(logPath);
		if (!ifs.fail())
		{
			//cout << logPath << " 文件已成功生成" << endl;
			cout << " | Nastran正在运行中。。。。" << endl;;
			break;
		}
		else if (clock() - brginTime > 20000)
		{
			cout << " | 未能正确运行Nastran！请检查 "<< BDFPath <<" 文件是否存在*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//判断nastran.exe是否运行结束
	while (true)
	{
		int count = 0;
		string process_name = "nastran.exe";

		count = find_process(process_name);
		if (count == 0)
		{
			cout << " | nastran.exe 运行结束" << endl;
			break;
		}
	}
	//读取log文件判断分析是否正常
	bool isFinishCalc = false;
	clock_t brginTime1 = clock();
	while (true)
	{
		ifstream ifs(logPath);
		string str_line;
		while (getline(ifs, str_line))
		{
			string::size_type idx = str_line.find("Analysis complete");
			if (idx != string::npos)//发现计算完成的标识符----跳出循环，准备提取数据
			{
				isFinishCalc = true;

				ifs.close();
				break;
			}

		}
		if (isFinishCalc)
		{
			//system("C:\\Users\\yycab\\Desktop\\bdf\\nastranTESTfile\\KillNastran.bat");//dos()杀掉nastran进程 懒得写 后面再加吧
			break;
		}
		else if (clock() - brginTime1 > 8000)
		{
			cout << " | 超过8s还没读到log分析成功，Nastran分析有误*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//判断pch文件是否生成成功
	clock_t brginTime2 = clock();
	while (true)
	{
		ifstream ifs(pchPath);
		if (!ifs.fail())
		{
			cout << " | " << pchPath << " 文件已成功生成" << endl;
			break;
		}
		else if (clock() - brginTime2 > 8000)
		{
			cout << " | 超过8s还未能正确生成 '.pch'！请检查 " << pchPath << " 文件是否存在*******" << endl;
			//system("pause");
			return -1;
		}
	}
	//Sleep(2000);//系统暂停两秒
	cout << "Nastran分析成功" << endl;
	return 0;
}

int myNastran::ReadResPCH(bool isOnlyReadDisp)
{
	//清除原始数据
	Displacements.clear();
	ElemStress.clear();
	MaxElemStress = pair<int, double>(0, 0);
	MaxElemStrain = pair<int, double>(0, 0);
	MaxNodeDisp = pair<int, Point>(0, Point());
	//打开结果文件
	ifstream ifs(pchPath);
	if (ifs.fail())
	{
		cout << "Nastran计算失败！.pch 文件不存在！请检查输入参数" << endl;
		return -1;
	}
	else cout << "数据结果文件打开成功！ 正在分析数据。。" << endl;
	int elemType = 0;
	int loadtype = 0;//1单元应力 2单元应变 3节点位移

	string str_line;
	while (getline(ifs, str_line))
	{
		//string::size_type idx = str_line.find("SUBTITLE");
		if (str_line.find("SUBTITLE") != string::npos)//发现标识符----准备提取数据
		{
			//getline(ifs, str_line);
			getline(ifs, str_line);//空读1行
			getline(ifs, str_line);

			if (str_line.find("STRAINS") != string::npos)			loadtype = 2;//cout << "找到了应变数据" << endl;
			else if (str_line.find("STRESSES") != string::npos)		loadtype = 1;//cout << "找到了应力数据" << endl;
			else if (str_line.find("DISPLACEMENTS") != string::npos)loadtype = 3;//找到了节点位移数据

			getline(ifs, str_line);
			getline(ifs, str_line);//空读两行
			if (loadtype != 3) getline(ifs, str_line);

			if (str_line.find("QUAD4") != string::npos)			elemType = 4;//cout << "找到了四边形单元" << endl;
			else if (str_line.find("TRIA3") != string::npos)	elemType = 3;//cout << "找到了三角形单元" << endl;
			else if (str_line.find("BAR") != string::npos)		elemType = 2;//cout << "找到了一维梁单元" << endl;
		}

		if (loadtype == 1)//应力
		{
			switch (elemType)
			{
			case 4://四边形单元
				cout << "正在读取四边形单元应力。。。";
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)//读到一直出现TITLE再停止
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
				loadtype = 0;//单元读取完毕，清空类型编号以免出错
				elemType = 0;//单元读取完毕，清空类型编号以免出错
				cout << "完成" << endl;
				break;
			case 3://三角形单元
				break;
			case 2://梁单元
				break;
			}
		}
		else if (loadtype == 2)//应变
		{
			switch (elemType)
			{
			case 4://四边形单元
				cout << "正在读取四边形单元应变。。。";
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
				loadtype = 0;//单元读取完毕，清空类型编号以免出错
				elemType = 0;//单元读取完毕，清空类型编号以免出错
				cout << "完成" << endl;
				break;
			case 3://三角形单元
				break;
			case 2://梁单元
				break;
			}
		}
		else if (loadtype == 3)//节点位移
		{
			cout << "正在读取节点位移。。。";
			while (getline(ifs, str_line))//读到一直出现TITLE再停止
			{
				if (str_line.find("TITLE") != string::npos)
				{
					break;
				}
				Point disp;
				int id = atoi(str_line.substr( 0, 17).c_str());//节点编号1
				disp.setX(atof(str_line.substr(23, 18).c_str()));//x方向位移
				disp.setY(atof(str_line.substr(41, 18).c_str()));//y方向位移
				disp.setZ(atof(str_line.substr(59, 18).c_str()));//z方向位移
				if (abs(MaxNodeDisp.second.getY()) < abs(disp.getY()))
				{
					MaxNodeDisp = pair<int, Point>(id, disp);
				}
				Displacements.insert(pair<int, Point>(id, disp));
				getline(ifs, str_line);
			}
			loadtype = 0;//单元读取完毕，清空类型编号以免出错
			elemType = 0;//单元读取完毕，清空类型编号以免出错
			cout << "完成" << endl;
			if (isOnlyReadDisp)
			{
				break;//当前只弄到节点位移，需要读其他量请把他删掉
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
		cout << "Nastran计算失败！.pch 文件不存在！请检查输入参数" << endl;
		return -1;
	}
	else cout << "数据结果文件打开成功！ 正在分析数据。。" << endl;

	int elemType = 0;
	int loadtype = 0;//1单元应力 2单元应变
	//记录外壳应力
	vector<double> maxstress(CQUAD4id.size());//存储最大应力值
	vector<int> maxstressID(CQUAD4id.size());//存储最大应力单元编号
	vector<double> maxstrain(CQUAD4id.size());//存储最大应变值
	vector<int> maxstrainID(CQUAD4id.size());//存储最大应变单元编号
	//记录隔框参数
	vector<double> maxstress2(CBARid.size());//存储最大应力值
	vector<int> maxstressID2(CBARid.size());//存储最大应力单元编号
	vector<double> maxstrain2(CBARid.size());//存储最大应变值
	vector<int> maxstrainID2(CBARid.size());//存储最大应变单元编号
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
		if (str_line.find("SUBTITLE") != string::npos)//发现标识符----准备提取数据
		{
			//getline(ifs, str_line);
			getline(ifs, str_line);//空读1行
			getline(ifs, str_line);

			if (str_line.find("STRAINS") != string::npos)			loadtype = 2;//cout << "找到了应变数据" << endl;
			else if (str_line.find("STRESSES") != string::npos)		loadtype = 1;//cout << "找到了应力数据" << endl;

			getline(ifs, str_line);
			getline(ifs, str_line);//空读两行
			getline(ifs, str_line);

			if (str_line.find("QUAD4") != string::npos)			elemType = 4;//cout << "找到了四边形单元" << endl;
			else if (str_line.find("TRIA3") != string::npos)	elemType = 3;//cout << "找到了三角形单元" << endl;
			else if (str_line.find("BAR") != string::npos)		elemType = 2;//cout << "找到了一维梁单元" << endl;
		}

		if (loadtype == 1)//应力
		{
			switch (elemType)
			{
			case 4://四边形单元
				cabinID = 0;
				getline(ifs, str_line);
				while (str_line.find("TITLE") == string::npos)//读到一直出现TITLE再停止
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
				loadtype = 0;//单元读取完毕，清空类型编号以免出错
				elemType = 0;//单元读取完毕，清空类型编号以免出错
				break;
			case 3://三角形单元
				break;
			case 2://梁单元
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
					loadtype = 0;//单元读取完毕，清空类型编号以免出错
					elemType = 0;//单元读取完毕，清空类型编号以免出错
				}
				break;
			}
		}
		else if (loadtype == 2)//应变
		{
			switch (elemType)
			{
			case 4://四边形单元
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
				loadtype = 0;//单元读取完毕，清空类型编号以免出错
				elemType = 0;//单元读取完毕，清空类型编号以免出错
				break;
			case 3://三角形单元
				break;
			case 2://梁单元
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
					loadtype = 0;//单元读取完毕，清空类型编号以免出错
					elemType = 0;//单元读取完毕，清空类型编号以免出错
				}
				break;
			}
		}
	}
	this->maxstress = maxstress;//存储最大应力值
	this->maxstressID = maxstressID;//存储最大应力单元编号
	this->maxstrain = maxstrain;//存储最大应变值
	this->maxstrainID = maxstrainID;//存储最大应变单元编号

	this->maxstress2 = maxstress2;//存储最大应力值
	this->maxstressID2 = maxstressID2;//存储最大应力单元编号
	this->maxstrain2 = maxstrain2;//存储最大应变值
	this->maxstrainID2 = maxstrainID2;//存储最大应变单元编号
	cout << "分析完成！数据已保存到指定文件夹" << endl;
	return 0;
}
