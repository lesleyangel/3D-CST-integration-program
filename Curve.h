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
    error,      //Ĭ���Ǵ��������
    cst,        //CST curve
    func,       //��������
    list,       //ɢ������
};

//ͨ��һԪ������
class Curve
{
public:
    CurveStyle curve_style{CurveStyle::error};
    vector<double> vec_info;//��¼����

protected:
    double max_value{1};//����ó��ĺ����������ֵ
private:
    double func_size{1};//�������ֵ����

public:
    Curve() {}
    ~Curve() {}

    void set_size(const double size) { func_size = size; } // �����������ֵ��С
    virtual double get(double x) = 0; // ��������
    double get_normal(double x);      // �����߹�һ������
    double get_diff(double x);
    vector<double> get(const vector<double> &x_list, bool is_normalization);

private:
    virtual void clac_span() = 0; //������������������ֵ����Сֵ
    virtual int check() = 0;//�����������Ƿ���Ϲ淶
};
using ptr_Curve = shared_ptr<Curve>;
// Curve2D
class Curve2D
{
public:
    enum Correlations
    {
        cartesian, // �ѿ�������ϵ
        polar,     // ������ϵ
    };

    ptr_Curve pt1;
    ptr_Curve pt2;
    Correlations m_correlations{cartesian}; // ����ϵ���Ĭ��Ϊ�ѿ�������ϵ

    public:
        Curve2D() = default;
        Curve2D(const ptr_Curve &curve, Correlations relation);
        void set_size(const double size) { pt2->set_size(size); } // �����������ֵ��С
        double get_diff(double t);
};

//----------------------------------------------------------------------
#define CURVE_FUNC                   \
public:                              \
    double get(double x);            \
private:                             \
    void clac_span();                \
    int check();


//ͨ�ú����ı����������
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
    std::vector<Point> m_list;//����������ڵ�������
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
    ptr_Curve pt1;   //ǰ����������� psi
    ptr_Curve pt2;    //������������ psi
    ptr_Curve pt3;   //���� eta
    ptr_Curve theta1; // Ĭ�����뻡���� 
    ptr_Curve theta2; // Ĭ�����뻡���� 
    ptr_Curve theta3; // Ĭ�����뻡����

    Point coordinate_position{Point(0, 0, 0)};//��ʼλ������
    Point coordinate_angle{Point(0, 0, 0)};//��ʼ�Ƕ�����
    Point coordinate_end_position{Point(0, 0, 0)};//����λ������
    Point coordinate_end_angle{Point(0, 0, 0)};//�����Ƕ�����
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
    Curve2D front;   //ǰ����������� psi
    Curve2D back;    //������������ psi
    Curve2D ridge;   //���� eta
    Curve2D lateral; //����������� eta
    GuideFunc guide;           //������ eta
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
            cout << "�������������������" << endl;
        }
        return res;
    }
private:

};
