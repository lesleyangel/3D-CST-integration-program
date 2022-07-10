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
    double max_value{1};

public:
    Curve() {}
    ~Curve() {}

    virtual double get(double x) = 0; //��������
    double get_normal(double x);      //�����߹�һ������
    vector<double> get(const vector<double> &x_list, bool is_normalization);

private:
    virtual void clac_span() = 0; //������������������ֵ����Сֵ
    virtual int check() = 0;//�����������Ƿ���Ϲ淶
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

    CST_info() = default;
    CST_info(const vector<double>& vec) {  N1 = vec[0];  N2 = vec[1];  }
    CST_info(double n1,double n2)  { N1 = n1; N2 = n2; }
};
class Curve_cst : public Curve
{CURVE_FUNC

private:
    CST_info m_info;
public:
    Curve_cst(const CST_info &info);

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

};
//
#undef CURVE_FUNC
//
using ptr_Curve = shared_ptr<Curve>;
class GuideFunc
{
public:
    ptr_Curve pt1;   //ǰ����������� psi
    ptr_Curve pt2;    //������������ psi
    ptr_Curve pt3;   //���� eta
    ptr_Curve theta1; // Ĭ������Ƕ��� 
    ptr_Curve theta2; // Ĭ������Ƕ��� 
    ptr_Curve theta3; // Ĭ������Ƕ���

public:
    GuideFunc();
    Point update_point(const Point &p, const double eta) const;
    void rotate(Point& p, const double eta) const;
    void translate(Point &p, const double eta) const;
};

class ClassFunc
{
public:
    ptr_Curve front;   //ǰ����������� psi
    ptr_Curve back;    //������������ psi
    ptr_Curve ridge;   //���� eta
    ptr_Curve lateral; //����������� eta
    GuideFunc guide;           //������ eta
public:
    Point get(const double eta, const double psi) const
    {
        Point res(0,0,0);
        if (front->curve_style == CurveStyle::cst && back->curve_style == CurveStyle::cst)
        {
            CST_info cst_info;
            cst_info.N1 = front->vec_info[0] + (back->vec_info[0] - front->vec_info[0]) * eta;
            cst_info.N2 = front->vec_info[1] + (back->vec_info[1] - front->vec_info[1]) * eta;
            unique_ptr<Curve> middle_eta(new Curve_cst(cst_info));

            res.setX(0);
            res.setY(0 + ridge->get_normal(eta) * middle_eta->get_normal(psi));
            res.setZ(0 + (2*psi-1) * lateral->get_normal(eta));
            res = guide.update_point(res, eta);
        }
        else
        {
            cout << "�������������������" << endl;
        }
        return res;
    }
private:

};
