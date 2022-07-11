#pragma once
#include <iostream>
#include <memory>
#include <functional>
#include "include/armadillo"
#include "factorial.h"
using namespace std;
//
enum class CurveStyle
{
    error,      //默认是错误的类型
    cst,        //CST curve
    func,       //任意曲线
    list,       //散点曲线
};

//通用一元函数类
class Curve
{
public:
    CurveStyle curve_style{CurveStyle::error};
    vector<double> vec_info;//记录参数

protected:
    double max_value{1};//计算得出的函数自身最大值
private:
    double func_size{1};//曲线最大值设置

public:
    Curve() {}
    ~Curve() {}

    void set_size(const double size) { func_size = size; } // 设置曲线最大值大小
    virtual double get(double x) = 0; // 获得坐标点
    double get_normal(double x);      // 将曲线归一化处理
    double get_diff(double x);
    vector<double> get(const vector<double> &x_list, bool is_normalization);

private:
    virtual void clac_span() = 0; //计算曲线因变量的最大值和最小值
    virtual int check() = 0;//检查输入参数是否符合规范
};
using ptr_Curve = shared_ptr<Curve>;
// Curve2D
class Curve2D
{
public:
    enum Correlations
    {
        cartesian, // 笛卡尔坐标系
        polar,     // 极坐标系
    };

    ptr_Curve pt1;
    ptr_Curve pt2;
    Correlations m_correlations{cartesian}; // 坐标系风格默认为笛卡尔坐标系

    public:
        Curve2D() = default;
        Curve2D(const ptr_Curve &curve, Correlations relation);
        void set_size(const double size) { pt2->set_size(size); } // 设置曲线最大值大小
        double get_diff(double t);
};

//----------------------------------------------------------------------
#define CURVE_FUNC                   \
public:                              \
    double get(double x);            \
private:                             \
    void clac_span();                \
    int check();


//通用函数的表达曲线类型
struct CST_info
{
    double N1;
    double N2;
    double Scale{1};

    CST_info() = default;
    CST_info(const vector<double>& vec) {  N1 = vec[0];  N2 = vec[1];  }
    CST_info(double n1,double n2, double scale = 1)  { N1 = n1; N2 = n2; Scale = scale; }
};
class Curve_cst : public Curve
{CURVE_FUNC
private:
    CST_info m_info;
public:
    Curve_cst(const CST_info &info);
    Curve_cst(CST_info &&info);
};
//
class Curve_func : public Curve
{CURVE_FUNC
private:
    std::function<double(double)> m_func;
public:
    Curve_func(const std::function<double(double)> &func);
};
//
class Curve_list : public Curve
{CURVE_FUNC
private:
    std::vector<Point> m_list;//必须是正则节点向量！
public:
    Curve_list(const vector<Point>& list);
    Curve_list(vector<Point>&& list);
};
//
#undef CURVE_FUNC
//

class GuideFunc
{
public:
    ptr_Curve pt1;   //前截面类别曲线 psi
    ptr_Curve pt2;    //后截面类别曲线 psi
    ptr_Curve pt3;   //脊线 eta
    ptr_Curve theta1; // 默认输入弧度制 
    ptr_Curve theta2; // 默认输入弧度制 
    ptr_Curve theta3; // 默认输入弧度制

    Point coordinate_position{Point(0, 0, 0)};//初始位置坐标
    Point coordinate_angle{Point(0, 0, 0)};//初始角度坐标
    Point coordinate_end_position{Point(0, 0, 0)};//结束位置坐标
    Point coordinate_end_angle{Point(0, 0, 0)};//结束角度坐标
public:
    GuideFunc();
    Point update_point(const Point &p, const double eta) const;
    Point get_theta(const double eta) const;
    void rotate(Point& p, const double eta) const;
    void translate(Point &p, const double eta) const;
};

class ClassFunc
{
public:
    Curve2D front;   //前截面类别曲线 psi
    Curve2D back;    //后截面类别曲线 psi
    Curve2D ridge;   //脊线 eta
    Curve2D lateral; //侧向类别曲线 eta
    GuideFunc guide;           //引导线 eta
public:
    Point get(const double eta, const double psi) const
    {
        Point res(0,0,0);
        if (front.pt2->curve_style == CurveStyle::cst && back.pt2->curve_style == CurveStyle::cst)
        {
            CST_info cst_info;
            cst_info.N1 = front.pt2->vec_info[0] + (back.pt2->vec_info[0] - front.pt2->vec_info[0]) * eta;
            cst_info.N2 = front.pt2->vec_info[1] + (back.pt2->vec_info[1] - front.pt2->vec_info[1]) * eta;
            Curve2D middle_eta(ptr_Curve(new Curve_cst(cst_info)), Curve2D::cartesian);

            res.setX(0);
            res.setY(0 + ridge.pt2->get_normal(eta) * middle_eta.pt2->get_normal(psi));
            res.setZ(0 + middle_eta.pt1->get(psi) * lateral.pt2->get_normal(eta));
            // res = guide.update_point(res, eta);
        }
        else
        {
            cout << "这个后期再制作！！！" << endl;
        }
        return res;
    }
private:

};
