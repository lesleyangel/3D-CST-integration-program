#include<iostream>
#include"include/armadillo"

using namespace arma;
using namespace std;
#define DISP_CONDITIONAL_CONVERGENCE 1e-2 //位移收敛条件
#define TIMES_FORCE 10000 //载荷放大倍数
#pragma once



//#include"factorial.h"

	//该文件是需要用到的一些数学实现、与具体内容无关

	int factorial(int a);//实现正整数阶乘的函数声明

	void toTecplot(mat node, mat elem, ofstream& ofs, string partname = "none");

	vector<string> split(string str, char pattern);

	//点类
	class PointSet
	{
	public:
		mat X;	mat Y;	mat Z;
		mat ifuse;//目前仅在2D壳结构模型时候使用，默认不使用该节点（等于1代表使用节点）
		void printP();			//打印表面节点P****4.12更新加入
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

	//网格其他信息类
	class MeshInfo
	{
	public:
		MeshInfo();
		mat FaiU;	mat FaiL;	mat EtaU;	mat EtaL;
		mat CURyx;	mat CURyzUpp;	mat CURyzLow;
		//评估网格信息
		double SurfaceAera;	double Volume;	double planforarea;
	};

	//网格类
	class Grid
	{
	public:
		mat P;			//节点信息Point
		mat E;			//单元编号信息Elem
		mat nodeForce;
		int pointNum() { return P.n_rows; }
		int elemtNum() { return E.n_rows; }
		Grid() { P.clear(); E.clear(); nodeForce.clear(); }
	};

	//该类用于清除重复点、整理单元编号
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
		mat P_new();//获取新节点集合
		mat getID();//获取新节点编号
		mat E_new(int elemNum = 0);//获取新的单元编号
		mat E_new(const mat& E, int elemNum = 0);
		//（elemNum参数表示单元前几列是节点编号信息 默认全都是）
		mat renewF(const mat& force);
		static bool isPointEqual(const mat& p1, const mat& p2, double delta = 1e-2);//判断两个坐标是否全等
		static vector<int> isPointEqual(const vector<mat>& p, double delta = 1e-2);//判断两个坐标是否全等
		static mat isPointEqual(mat p, double delta = 1e-2);//判断两个坐标是否全等
	private:
		bool isEqual(const mat& A, const mat& B);//判断两个mat是否全等

		void CleanMat();//清除网格重复单元
		void getNewE(int elemNum = 0);//获得新的单元集

		mat P_old;//原始传入的点集
		mat selfP;//清除完重复节点后的点集
		mat id; //原点集对应新点集的编号
		mat id_num;//密铺的id号

		mat E_old;//原始传入的单元节点
		mat selfE;//清除完重复节点后的单元集
	};




