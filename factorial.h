#include<iostream>
#include"include/armadillo"

using namespace arma;
using namespace std;
#define DISP_CONDITIONAL_CONVERGENCE 1e-2 //λ����������
#define TIMES_FORCE 10000 //�غɷŴ���
#pragma once



//#include"factorial.h"

	//���ļ�����Ҫ�õ���һЩ��ѧʵ�֡�����������޹�

	int factorial(int a);//ʵ���������׳˵ĺ�������

	void toTecplot(mat node, mat elem, ofstream& ofs, string partname = "none");

	vector<string> split(string str, char pattern);

	//����
	class PointSet
	{
	public:
		mat X;	mat Y;	mat Z;
		mat ifuse;//Ŀǰ����2D�ǽṹģ��ʱ��ʹ�ã�Ĭ�ϲ�ʹ�øýڵ㣨����1����ʹ�ýڵ㣩
		void printP();			//��ӡ����ڵ�P****4.12���¼���
		//mat getPoint(int i, int j);
	};

	//template<class T>
	//class Type;
	
	class Point
	{
	public:
		Point() { x = 0; y = 0; z = 0;  }
		Point(double x, double y = 0, double z = 0) { SetSite(x, y, z); }
		
		Point(const mat& point);
		void SetSite(double x, double y = 0, double z = 0) { this->x = x; this->y = y; this->z = z; }
		void setX(double x) { this->x = x; }
		void setY(double y) { this->y = y; }
		void setZ(double z) { this->z = z; }
		
		double getX() const { return x; }
		double getY() const { return y; }
		double getZ() const { return z; }
		tuple<double, double, double> to_tuple() { return tuple<double, double, double>(x, y, z); }
		vec to_vec() { vec res = {x, y, z}; return res;	}

		double Norm1() const { return abs(x) + abs(y) + abs(z); }
		double Norm2() const { return sqrt(x * x + y * y + z * z); }
		double Normi(int i) const { return pow(pow(abs(x), i) + pow(abs(y), i) + pow(abs(z), i), 1.0 / (double)i); }
		
		bool operator<(const Point& a) { return (*this).x < a.x; }
		bool operator>(const Point& a) { return (*this).x > a.x; }
		bool operator==(const Point& a) { return ((*this).x == a.x && (*this).y == a.y && (*this).z == a.z); }
		Point operator+(const Point& a) { return Point((*this).x + a.x, (*this).y + a.y, (*this).z + a.z); }
		Point operator-(const Point& a) { return Point((*this).x - a.x, (*this).y - a.y, (*this).z - a.z); }
		Point operator-=(const Point& a) { (*this) = (*this) - a; return (*this); }
		Point operator+=(const Point& a) { (*this) = (*this) + a; return (*this); }
		Point operator-(const double& n) { return Point((*this).x - n, (*this).y - n, (*this).z - n); }
		Point operator/(const double& n) { return Point((*this).x / n, (*this).y / n, (*this).z / n); }
		Point operator*(const double& n) { return Point((*this).x * n, (*this).y * n, (*this).z * n); }
	protected:
		double x;
		double y;
		double z;
		
		//Type<T> a;
	};
	static bool sort_token_point(const Point& p1, const Point& p2)
	{
		return p1.getX() < p2.getX();
	}
	
	inline Point::Point(const mat& point)
	{
		if (point.n_rows == 1 && point.n_cols > 2)
			SetSite(point(0), point(1), point(2));
		else if (point.n_cols == 1 && point.n_rows > 2)
			SetSite(point(0), point(1), point(2));
		else
			SetSite(0, 0, 0);
	}
	
	inline Point cross(const Point& a, const Point& b)
	{
		double x = a.getY() * b.getZ() - a.getZ() * b.getY();
		double y = a.getZ() * b.getX() - a.getX() * b.getZ();
		double z = a.getX() * b.getY() - a.getY() * b.getX();
		return Point(x, y, z);
	}


	template<class T>
	class Field : public Point
	{
	public:
		Field() { x = 0; y = 0; z = 0; p = 0; }
		Field(double x, double y = 0, double z = 0) { SetSite(x, y, z); p = 0; }
		Field(double x, double y, double z, T p) { SetSite(x, y, z); this->p = p; }
		void setP(T p) { this->p = p; }
		T getP() const { return p; }
	private:
		T p;
	};
	
	vector<Field<double>> getPlistFromfile(string filepath);

	//����������Ϣ��
	class MeshInfo
	{
	public:
		MeshInfo();
		mat FaiU;	mat FaiL;	mat EtaU;	mat EtaL;
		mat CURyx;	mat CURyzUpp;	mat CURyzLow;
		//����������Ϣ
		double SurfaceAera;	double Volume;	double planforarea;
	};

	//������
	class Grid
	{
	public:
		mat P;			//�ڵ���ϢPoint
		mat E;			//��Ԫ�����ϢElem
		mat nodeForce;
		int pointNum() { return P.n_rows; }
		int elemtNum() { return E.n_rows; }
		Grid() { P.clear(); E.clear(); nodeForce.clear(); }
	};

	//������������ظ��㡢����Ԫ���
	class uniqueMat
	{
	public:
		uniqueMat(mat p_old)
		{
			selfP.clear(); id.clear(); P_old = p_old; E_old.clear(); selfE.clear();
		}
		uniqueMat(mat p_old, mat e_old)
		{
			selfP.clear(); id.clear(); P_old = p_old; E_old = e_old; selfE.clear();
		}
		uniqueMat(Grid gid)
		{
			selfP.clear(); id.clear(); P_old = gid.P; E_old = gid.E; selfE.clear();
		}
		mat P_new();//��ȡ�½ڵ㼯��
		mat getID();//��ȡ�½ڵ���
		mat E_new(int elemNum = 0);//��ȡ�µĵ�Ԫ���
		mat E_new(const mat& E, int elemNum = 0);
		//��elemNum������ʾ��Ԫǰ�����ǽڵ�����Ϣ Ĭ��ȫ���ǣ�
		mat renewF(const mat& force);
		static bool isPointEqual(const mat& p1, const mat& p2, double delta = 1e-2);//�ж����������Ƿ�ȫ��
		static vector<int> isPointEqual(const vector<mat>& p, double delta = 1e-2);//�ж����������Ƿ�ȫ��
		static mat isPointEqual(mat p, double delta = 1e-2);//�ж����������Ƿ�ȫ��
	private:
		bool isEqual(const mat& A, const mat& B);//�ж�����mat�Ƿ�ȫ��

		void CleanMat();//��������ظ���Ԫ
		void getNewE(int elemNum = 0);//����µĵ�Ԫ��

		mat P_old;//ԭʼ����ĵ㼯
		mat selfP;//������ظ��ڵ��ĵ㼯
		mat id; //ԭ�㼯��Ӧ�µ㼯�ı��
		mat id_num;//���̵�id��

		mat E_old;//ԭʼ����ĵ�Ԫ�ڵ�
		mat selfE;//������ظ��ڵ��ĵ�Ԫ��
	};




