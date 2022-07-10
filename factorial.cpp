#include "factorial.h"
#include <iomanip>
int factorial(int a)
{
	if (a < 0)
	{
		cout << "使用错误factorial，N必须为非负实整数" << endl;
		return 0;
	}
	else if(a == 0)
	{
		return 1;
	}
	int b = 1;
	for (int i = 0; i < a; i++)
	{
		b = b * (i + 1);
	}
	return b;
}

void toTecplot(mat node, mat elem, ofstream& ofs, string partname)
{
	ofs << "Title = \""<< partname <<"\"" << endl;
	ofs << "variables = \"x\"" << endl;
	ofs << "\"y\"" << endl;
	ofs << "\"z\"" << endl;
	ofs << "ZONE T = \"ZONE 001\"" << endl;
	ofs << "STRANDID = 0, SOLUTIONTIME = 0" << endl;
	ofs << "Nodes=" << node.n_rows << ",Elements=" << elem.n_rows << ",datapacking=block" << endl;
	ofs << "zonetype=fequadrilateral" << endl;
	ofs << "DATAPACKING=BLOCK" << endl;
	for (int m = 0; m < 3; m++)//sa
	{
		for (size_t n = 0; n < node.n_rows; n++)
		{
			if (n % 5 != 4)		ofs << setw(12) << node(n, m) << " ";
			else				ofs << setw(12) << node(n, m) << endl;
		}
		ofs << endl;
	}
	ofs << elem << endl;
}

vector<string> split(string str, char pattern)
{
	vector<string> res;//存储结果
	stringstream ss(str);//通过字符串str初始化ss
	string temp;
	//使用getline，通过pattern分割stringstream
	while (getline(ss, temp, pattern))
	{
		if (temp.size() != 0)//如果结果不是空
		{
			res.push_back(temp);//结果存入temp
		}
	}
	return res;
}

vector<Field<double>> getPlistFromfile(string filepath)
{
	ifstream ifs(filepath);
	if (ifs.fail())
	{
		return vector<Field<double>>();
	}
	vector<Field<double>> plist;
	string str_tmp;
	while (getline(ifs, str_tmp))
	{
		double x, y, z, p;
		stringstream ss(str_tmp);
		stringstream iss;
		getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> x; iss.clear();
		getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> y; iss.clear();
		getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> z; iss.clear();
		getline(ss, str_tmp, ',');	iss << str_tmp;	iss >> p; iss.clear();
		plist.push_back(Field<double>(x, y, z, p));
	}
	return plist;
}

//打印表面节点P
void PointSet::printP()
{
	cout << "thia is print P" << endl;
	cout << "X = " << endl;
	cout << X << endl;
	cout << "Y = " << endl;
	cout << Y << endl;
	cout << "Z = " << endl;
	cout << Z << endl;
	cout << "ifUse = " << endl;
	cout << ifuse << endl;
	cout << "- - - - - - - - -" << endl;
}


//template<class T>
//Point<T> cross(const Point<T>& a, const Point<T>& b)
//{
//
//}

///////////////////////////////////////

mat uniqueMat::P_new()
{
	if (selfP.size()==0)//如果self为空，则运行CleanMat计算self
	{
		CleanMat();
		return selfP;
	}
	else//如果self有值，则直接调用
	{
		return selfP;
	}
}

mat uniqueMat::getID()
{
	if (id.size() == 0)//如果self为空，则运行CleanMat计算self
	{
		CleanMat();
		return id;
	}
	else//如果self有值，则直接调用
	{
		return id;
	}
}

mat uniqueMat::E_new(int elemNum)
{
	if (E_old.size() == 0)
	{
		cout << "未传入原始单元信息！！" << endl;
		return E_old;
	}
	if (id.size() == 0)
	{
		CleanMat();//清除网格重复单元
	}
	if (selfE.size() == 0)
	{
		getNewE(elemNum);//获得新的单元集
		return selfE;
	}
	else
	{
		return selfE;
	}
}

mat uniqueMat::E_new(const mat& E, int elemNum)
{
	E_old = E;
	if (E_old.size() == 0)
	{
		cout << "未传入原始单元信息！！" << endl;
		return E_old;
	}
	if (id.size() == 0)
	{
		CleanMat();//清除网格重复单元
	}
	
	getNewE(elemNum);//获得新的单元集
	return selfE;
	
}

mat uniqueMat::renewF(const mat& force)
{
	if (force.n_rows != P_old.n_rows)
	{
		cout << "renewF函数传入参数与节点数无法对应！更新节点力错误！！" << endl;
		return mat();
	}
	mat F_new = zeros(P_new().n_rows, 6);
	for (size_t j = 0; j < F_new.n_rows; j++)
	{
		for (size_t i = 0; i < id.n_rows; i++)
		{
			if (j == id(i))
			{
				F_new.row(j) += force.row(i);
			}
		}
	}
	return F_new;
}

void uniqueMat::CleanMat()
{
	cout << "正在调用CleanMat清除重复节点......";
	id = zeros(P_old.n_rows, 1);
	//初始化id号码
	for (size_t i = 0; i < P_old.n_rows; i++)
	{
		id(i, 0) = i;
	}
	//比较和清理重复的id
	for (size_t i = 0; i < P_old.n_rows - 1; i++)
	{
		for (size_t j = i + 1; j < P_old.n_rows; j++)
		{
			if (id(i,0)<i)//如果该行id号小于行号，跳出循环
			{
				break;
			}
			if (id(j,0)==j)//如果j的id号不等于行号，不做判断
			{
				bool isEqual = true;//判断两行内容是否相等
				for (size_t n = 0; n < P_old.n_cols; n++)
				{
					if (abs(P_old(i, n) - P_old(j, n)) > 1e-3)
					{
						isEqual = false;
						break;
					}
				}
				if (isEqual)//如果i j两行相等，则i的id赋给j
				{
					id(j, 0) = id(i, 0);
				}
				//if (isEqual(P_old.row(i), P_old.row(j)))//如果i j两行相等，则i的id赋给j
				//{
				//	id(j, 0) = id(i, 0);
				//}
			}
			
		}
	}
	//id号密铺排列
	id_num = unique(id);
	for (size_t num = 0; num < id_num.n_rows; num++)
	{
		for (size_t i = 0; i < id.n_rows; i++)
		{
			if (id(i, 0) == id_num(num,0))
			{
				id(i, 0) = num;
			}
		}
	}
	//得到清除后的节点
	selfP = zeros(id_num.n_rows, P_old.n_cols);
	for (size_t i = 0; i < selfP.n_rows; i++)
	{
		selfP.row(i) = P_old.row((size_t)id_num(i));
	}
	cout << "清除成功！ 重复节点数：" << P_old.n_rows - selfP.n_rows << endl;
}

void uniqueMat::getNewE(int elemNum)
{
	if (elemNum < 2 || (size_t)elemNum > E_old.n_cols)
	{
		elemNum = E_old.n_cols;
	}
	selfE = E_old;// zeros(E_old.n_rows, E_old.n_cols);
	for (size_t i = 0; i < E_old.n_rows; i++)
	{
		for (int j = 0; j < elemNum; j++)
		{
			if (E_old(i, j) != id((size_t)E_old(i, j), 0))
			{
				selfE(i, j) = id((size_t)E_old(i, j), 0);
			}
		}
	}
}

bool uniqueMat::isEqual(const mat &A, const mat &B)//判断两个矩阵是否全等
{
	size_t a_rows = A.n_rows;
	size_t a_cols = A.n_cols;
	size_t b_rows = B.n_rows;
	size_t b_cols = B.n_cols;

	if (a_rows == b_rows && a_cols == b_cols)
	{
		double delta = 0;
		for (size_t i = 0; i < a_rows; i++)
		{
			for (size_t j = 0; j < a_cols; j++)
			{
				delta += pow(A(i, j) - B(i, j),2);
			}
		}
		if (sqrt(delta) < 1e-2)
		{
			return true;
		}
		else
		{
			return false;
		}
		//return true;
	}
	else
	{
		return false;
	}
}

bool uniqueMat::isPointEqual(const mat& p1, const mat& p2, double delta /*= 1e-2*/)//判断两个矩阵是否全等
{

	if (p1.n_elem == p2.n_elem)
	{
		double temp = 0;
		for (size_t i = 0; i < p1.n_elem; i++)
		{
			//temp += pow(p1(i) - p2(i), 2);
			temp += abs(p1(i) - p2(i));
		}
		if (temp < delta)
		{
			return true;
		}
		else
		{
			return false;
		}
		//return true;
	}
	else
	{
		return false;
	}
}

vector<int> uniqueMat::isPointEqual(const vector<mat>& p, double delta)
{
	vector<int> equallist(p.size());
	for (size_t i = 0; i < p.size(); i++)
	{
		equallist[i] = i;
	}
	for (size_t i = 0; i < p.size(); i++)
	{
		if (equallist[i] == i)//如果第i个节点不是重复节点
		{
			for (size_t j = i + 1; j < p.size(); j++)
			{
				if (equallist[j] == j)
				{
					if (isPointEqual(p[i], p[j], delta))
					{
						equallist[j] = i;
					}
				}
				
			}
		}
	}
	return equallist;
}

mat uniqueMat::isPointEqual(mat p, double delta)
{
	mat equallist(p.n_rows,1);
	for (size_t i = 0; i < p.n_rows; i++)
	{
		equallist[i] = i;
	}
	for (size_t i = 0; i < p.n_rows; i++)
	{
		if (equallist[i] == i)//如果第i个节点不是重复节点
		{
			for (size_t j = i + 1; j < p.n_rows; j++)
			{
				if (equallist[j] == j)
				{
					if (isPointEqual(p.row(i), p.row(j), delta))
					{
						equallist(j) = i;
					}
				}

			}
		}
	}
	return equallist;
}
