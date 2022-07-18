#include "Curve.h"
#define PI 3.1415926
vector<double> Curve::get(const vector<double> &x_list, bool is_normalization)
{
    vector<double> res(x_list.size(), 0);
    if (is_normalization)
    {
        for (size_t i = 0; i < x_list.size(); i++)
        {
            res[i] = get_normal(x_list[i]);
        }
    }
    else
    {
        for (size_t i = 0; i < x_list.size(); i++)
        {
            res[i] = get(x_list[i]);
        }
    }
    return res;
}

double Curve::get_normal(double x)
{
    double res = get(x) / abs(max_value) * func_size;
    return get(x) / abs(max_value) * func_size;
}

double Curve::get_diff(double x)
{
    double dx = 0.01;
    if (x<0 || x> 1)
    {
        cout << "输入了错误的自变量参数：x = " << x << endl;
        return -100000;
    }
    else if (x > (1-dx))
    {
        return (get(x) - get(x - dx)) / dx;
    }
    else
    {
        return (get(x + dx) - get(x)) / dx;
    }
    return 0;
}

Curve2D::Curve2D(const ptr_Curve &curve, Correlations relation)
{
    switch (relation)
    {
    case cartesian:
        // y(x) = f(x); x => 0:1
        // Convert to
        // x(t) = 2 * t - 1; t => 0:1
        // y(t) = f(t); t => 0:1
        pt1 = ptr_Curve(new Curve_list({Point(0, -1), Point(1, 1)}));
        pt2 = curve;
        break;
    case polar:
        // r(theta) = f(theta); theta => 0:pi
        // Convert to
        // x(t) = f(t*pi) * cos(t*pi); t => 0:1
        // y(t) = f(t*pi) * sin(t*pi); t => 0:1
        pt1 = ptr_Curve(new Curve_func([&](double t)
                                       { return curve->get(t * PI) * cos(t * PI); }));
        pt2 = ptr_Curve(new Curve_func([&](double t)
                                       { return curve->get(t * PI) * sin(t * PI); }));
    default:
        break;
    }
}
double Curve2D::get_diff(double t)
{
    return pt2->get_diff(t) / pt1->get_diff(t);
}
//------------------------------------
Curve_cst::Curve_cst(const CST_info &info)
{
    m_info = info;
    curve_style = CurveStyle::cst;
    clac_span();
    vec_info.resize(2);
    vec_info[0] = m_info.N1;
    vec_info[1] = m_info.N2;
}
Curve_cst::Curve_cst(CST_info &&info)
{
    m_info = info;
    curve_style = CurveStyle::cst;
    clac_span();
    vec_info.resize(2);
    vec_info[0] = m_info.N1;
    vec_info[1] = m_info.N2;
}
double Curve_cst::get(double x)
{
    double res = pow(x, m_info.N1) * pow(1 - x, m_info.N2) * ((m_info.Scale - 1) * x + 1);
    if (isnan(res))
    {
        return 0;
    }
    return res;
}
void Curve_cst::clac_span()
{
    const double n_max = m_info.N1 / (m_info.N1 + m_info.N2);
    const double fai_max = pow(n_max, m_info.N1) * pow(1 - n_max, m_info.N2);    
    max_value = isnan(fai_max) ? 1.0 : fai_max;
}
int Curve_cst::check()
{
    if (m_info.N1 < 0)
        return -1;
    if (m_info.N2 < 0)
        return -1;
    return 0;
}
//
Curve_func::Curve_func(const std::function<double(double)> &func)
{
    m_func = func;
    curve_style = CurveStyle::func;
    clac_span();
}
double Curve_func::get(double x)
{
    return m_func(x);
}
void Curve_func::clac_span()
{
    double max = abs(get(0));
    const size_t size = 10000;
    for (size_t i = 0; i < size; i++)
    {
        const double value = get((double)i / ((double)size - 1));
        if (abs(max) < abs(value))
            max = value;
    }
    max_value = max;
}
int Curve_func::check()
{
    return 0;
}
//
Curve_list::Curve_list(const std::vector<Point> &list)
{
    m_list = list;
    curve_style = CurveStyle::list;
    clac_span();
}
Curve_list::Curve_list(std::vector<Point> &&list)
{
    m_list = list;
    curve_style = CurveStyle::list;
    clac_span();
}
double Curve_list::get(double x)
{
    double x0 = m_list[0].getX();
    double x1 = m_list[0].getX();
    double y0 = m_list[0].getY();
    double y1 = m_list[0].getY();
    for (size_t i = 1; i < m_list.size(); i++)
    {
        if (x <= m_list[i].getX() + 1e-5) // +1e-5避免网格计算出现1.000002这样的误差
        {
            x0 = m_list[i - 1].getX();
            x1 = m_list[i].getX();
            y0 = m_list[i - 1].getY();
            y1 = m_list[i].getY();
            break;
        }
    }
    const double res = y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    if (isnan(res))
    {
        cout << "Curve_list::get(double x) 返回了错误的数值！nan" << endl;
    }
    return res;//??????
}
void Curve_list::clac_span()
{
    max_value = abs(m_list[0].getY());//???????
    for (size_t i = 0; i < m_list.size(); i++)
    {
        if (abs(max_value) < abs(m_list[i].getY()))
            max_value = abs(m_list[i].getY());
    }
    
}
int Curve_list::check()
{
    if (m_list.size() < 1)
    {//??????????е?
        cout << "??????????е?" << endl;
        return -1;
    }
    else if (m_list.size() == 1)
    {
        return 0;
    }
    
    for (size_t i = 1; i < m_list.size(); i++)
    {
        if (m_list[i-1].getX() > m_list[i].getX())
        {//?????i???????????????i-1??

            return -2;
        }
    }
    return 0;
}
//
//--------------------------------------------------------------
//
GuideFunc::GuideFunc()
{
    vector<Point> list_one(2);
    list_one[0].SetSite(0, 0, 0);
    list_one[1].SetSite(1, 1, 0);
    vector<Point> list_zero(2);
    list_zero[0].SetSite(0, 0, 0);
    list_zero[1].SetSite(1, 0, 0);
    pt1 = ptr_Curve(new Curve_list(list_one));
    pt2 = ptr_Curve(new Curve_list(list_zero));
    pt3 = ptr_Curve(new Curve_list(list_zero));
    theta1 = ptr_Curve(new Curve_list(list_zero));
    theta2 = ptr_Curve(new Curve_list(list_zero));
    theta3 = ptr_Curve(new Curve_list(list_zero));
}
Point GuideFunc::update_point(const Point& p, const double eta) const
{
    Point res = p;
    rotate(res, eta);
    translate(res, eta);
    return res;
}

Point GuideFunc::get_theta(const double eta)const
{
    return Point(
        theta1->get(eta)+ this->coordinate_angle.getX(),
        theta2->get(eta)+ this->coordinate_angle.getY(),
        theta3->get(eta)+ this->coordinate_angle.getZ());
}

void GuideFunc::rotate(Point& p, const double eta)const 
{
    Point theta = get_theta(eta);
    const double _theta1 = theta.getX();
    const double _theta2 = theta.getY();
    const double _theta3 = theta.getZ();
    //
    const arma::mat33 r1 = {
        {cos(_theta1), sin(_theta1), 0},
        {-sin(_theta1), cos(_theta1), 0},
        {0, 0, 1},
    };
    const arma::mat33 r2 = {
        {cos(_theta2), 0, sin(_theta2)},
        {0, 1, 0},
        {-sin(_theta2), 0, cos(_theta2)},
    };
    const arma::mat33 r3 = {
        {1, 0, 0},
        {0,cos(_theta3), sin(_theta3)},
        {0,-sin(_theta3), cos(_theta3)},
    };
    arma::vec p_mat = { {p.getX()}, {p.getY()}, {p.getZ()} };
    arma::mat res_mat = r3 * r2 * r1 * p_mat;//???????
    p.SetSite(res_mat(0), res_mat(1), res_mat(2));
}

void GuideFunc::translate(Point& p, const double eta)const
{
    const double _theta1 = coordinate_angle.getX();
    const double _theta2 = coordinate_angle.getY();
    const double _theta3 = coordinate_angle.getZ();
    //
    const arma::mat33 r1 = {
        {cos(_theta1), sin(_theta1), 0},
        {-sin(_theta1), cos(_theta1), 0},
        {0, 0, 1},
    };
    const arma::mat33 r2 = {
        {cos(_theta2), 0, sin(_theta2)},
        {0, 1, 0},
        {-sin(_theta2), 0, cos(_theta2)},
    };
    const arma::mat33 r3 = {
        {1, 0, 0},
        {0,cos(_theta3), sin(_theta3)},
        {0,-sin(_theta3), cos(_theta3)},
    };
    arma::vec pt = {
        pt1->get(eta),
        pt2->get(eta),
        pt3->get(eta),
    };
    arma::mat res_mat = r3 * r2 * r1 * pt;
    p.setX(res_mat(0) + p.getX());
    p.setY(res_mat(1) + p.getY());
    p.setZ(res_mat(2) + p.getZ());
    p += coordinate_position;
    
}