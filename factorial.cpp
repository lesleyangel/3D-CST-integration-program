#include "factorial.h"
#include <iomanip>
int factorial(int a)
{
	if (a < 0)
	{
		cout << "ʹ�ô���factorial��N����Ϊ�Ǹ�ʵ����" << endl;
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
	vector<string> res;//�洢���
	stringstream ss(str);//ͨ���ַ���str��ʼ��ss
	string temp;
	//ʹ��getline��ͨ��pattern�ָ�stringstream
	while (getline(ss, temp, pattern))
	{
		if (temp.size() != 0)//���������ǿ�
		{
			res.push_back(temp);//�������temp
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

//��ӡ����ڵ�P
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
	if (selfP.size()==0)//���selfΪ�գ�������CleanMat����self
	{
		CleanMat();
		return selfP;
	}
	else//���self��ֵ����ֱ�ӵ���
	{
		return selfP;
	}
}

mat uniqueMat::getID()
{
	if (id.size() == 0)//���selfΪ�գ�������CleanMat����self
	{
		CleanMat();
		return id;
	}
	else//���self��ֵ����ֱ�ӵ���
	{
		return id;
	}
}

mat uniqueMat::E_new(int elemNum)
{
	if (E_old.size() == 0)
	{
		cout << "δ����ԭʼ��Ԫ��Ϣ����" << endl;
		return E_old;
	}
	if (id.size() == 0)
	{
		CleanMat();//��������ظ���Ԫ
	}
	if (selfE.size() == 0)
	{
		getNewE(elemNum);//����µĵ�Ԫ��
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
		cout << "δ����ԭʼ��Ԫ��Ϣ����" << endl;
		return E_old;
	}
	if (id.size() == 0)
	{
		CleanMat();//��������ظ���Ԫ
	}
	
	getNewE(elemNum);//����µĵ�Ԫ��
	return selfE;
	
}

mat uniqueMat::renewF(const mat& force)
{
	if (force.n_rows != P_old.n_rows)
	{
		cout << "renewF�������������ڵ����޷���Ӧ�����½ڵ������󣡣�" << endl;
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
	cout << "���ڵ���CleanMat����ظ��ڵ�......";
	id = zeros(P_old.n_rows, 1);
	//��ʼ��id����
	for (size_t i = 0; i < P_old.n_rows; i++)
	{
		id(i, 0) = i;
	}
	//�ȽϺ������ظ���id
	for (size_t i = 0; i < P_old.n_rows - 1; i++)
	{
		for (size_t j = i + 1; j < P_old.n_rows; j++)
		{
			if (id(i,0)<i)//�������id��С���кţ�����ѭ��
			{
				break;
			}
			if (id(j,0)==j)//���j��id�Ų������кţ������ж�
			{
				bool isEqual = true;//�ж����������Ƿ����
				for (size_t n = 0; n < P_old.n_cols; n++)
				{
					if (abs(P_old(i, n) - P_old(j, n)) > 1e-3)
					{
						isEqual = false;
						break;
					}
				}
				if (isEqual)//���i j������ȣ���i��id����j
				{
					id(j, 0) = id(i, 0);
				}
				//if (isEqual(P_old.row(i), P_old.row(j)))//���i j������ȣ���i��id����j
				//{
				//	id(j, 0) = id(i, 0);
				//}
			}
			
		}
	}
	//id����������
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
	//�õ������Ľڵ�
	selfP = zeros(id_num.n_rows, P_old.n_cols);
	for (size_t i = 0; i < selfP.n_rows; i++)
	{
		selfP.row(i) = P_old.row((size_t)id_num(i));
	}
	cout << "����ɹ��� �ظ��ڵ�����" << P_old.n_rows - selfP.n_rows << endl;
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

bool uniqueMat::isEqual(const mat &A, const mat &B)//�ж����������Ƿ�ȫ��
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

bool uniqueMat::isPointEqual(const mat& p1, const mat& p2, double delta /*= 1e-2*/)//�ж����������Ƿ�ȫ��
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
		if (equallist[i] == i)//�����i���ڵ㲻���ظ��ڵ�
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
		if (equallist[i] == i)//�����i���ڵ㲻���ظ��ڵ�
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
