/*
插值函数V1.2.3
N维插值适用于两种标准气动文件格式的插值，支持面元法输出格式行数乱序的情况，支持有维度表头长度为1的情况，支持逆序插值，自动识别顺序逆序。
附带的一维插值和二维插值为通用格式，效率与N维插值类似,参数带vector的支持逆序插值，自动识别顺序逆序。
附带的三维和四维插值在插值点超过下限时会报错，且效率比N维插值低，不支持逆序插值。

日期：2019.12.02
*/
#pragma once
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
using namespace std;
#define DLLEXPORT __declspec(dllexport)
//  插值表数据类型


// class DLLEXPORT Interpolation
class Interpolation
{
public:
	Interpolation(void);
	~Interpolation(void);

public:
	//初始化数据，两种方式 1：读txt；2：手动设置
	bool InterpTableRead(string FileName);                  // 读取老格式插值表，准备废弃
	bool InterpTableReadFile(string FileName);			    //	读取插值表txt，新格式
	bool InterpTableReadData(string strdata);               //  读取字符串，新格式
	
	bool SetHeader(vector<vector<double>> header,vector<string> headername);//手动设置表头
	bool SetData(vector<int> cubeorder, vector<vector<double>> cubedata);   //手动设置数据
	
	//插值主函数
	double InterpolationN(const vector<double>& x, bool out = false);//	插值函数,out=1允许外插，=0不允许外插

	//获取数据
	vector<vector<double>> GetData(const vector<int>& cubeorder);//获取指定位置的数据矩阵
	int GetDim();											        //	读取维度
	const vector<int>& GetHeaderLength();							//	读取表头长度
	const vector<vector<double>>& GetHeader();						//	读取表头
	const vector<string>& GetHeaderName();							//	读取表头名称
	bool WriteTxt(string txtpath);//数据写成标准格式
	void CreateZeroInterp();//建一个为0的插值表
	int Check();  // 数据检测
public:	 //  for serializing
	bool IfOut;						//	IfOut=1允许外插，=0不允许外插
	bool Ifturn;					//	如果为ture，则输入的插值点逆序操作
	int	m_Dim;			            //	维数
	vector<int> m_HeaderLength;		//	各维列表长度
	vector<vector<double>> m_Header;//	各维列表
	vector<double> m_data;			//	数据
	vector<string> m_HeaderName;    //  表头名字
private:
	bool CheckTable();	//	检查表头是否升序排列，检查data长度与表头是否匹配
	double Interp1(double x1, double x2, double y1, double y2, double x);//	一维插值计算
	double InterpN(int start, const vector<double>& x, int Dim, const vector<vector<int>>& XNum, const vector<vector<double>>& XVal);//	n维插值计算
	void WriteCube(int n,fstream &ff); //写一个方块的数据
	
//暂时不用的，针对另一种气动数据格式
private:
	void ReadAero(const char *FileName,int n);	//读新的气动文件,n为读取用来插值数据的列数-1,针对面元法输出的气动格式
	void InterpTableRead2(const char *FileName,int n);	//支持乱序输入的ReadAreo,n为读取用来插值数据的列数-1
	int FindPosition(const vector<double>& head); //InterpTableRead2用，寻找插到m_AllData中的地方
	vector<double> Shorten(const vector<double>& a); //缩短表头用	
public:
	vector<vector<double>> m_AllData;	//InterpTableRead2读出来的数据存在这里
};

//////////////////////////////////////////////////////////////////////////
//一维差值
double Interpolation1(int n_Xs, double *Xs,double *Data, double X );
extern "C" DLLEXPORT double Interpolation1(const vector<double> &Xs,const vector<double>& Data, double X );
//二维差值
double Interpolation2(int n_Xs, double *Xs, int n_Ys, double *Ys, double *Data, double X, double Y);
extern "C" DLLEXPORT double Interpolation2(const vector<double> &Xs, const vector<double> &Ys, const vector<double> &Data, double X, double Y);
//三维差值
double Interpolation3(int n_Xs, double *Xs, int n_Ys, double *Ys, int n_Zs, double *Zs, double *Data, double X, double Y, double Z );
extern "C" DLLEXPORT double Interpolation3(const vector<double> &Xs, const vector<double> &Ys, const vector<double> &Zs, const vector<double> &Data, double X, double Y, double Z );
//四维差值                 
double Interpolation4(int n_X1, double *X1s, int n_X2, double *X2s, int n_X3, double *X3s, int n_X4s, double *X4s, double *Data, double X1, double X2, double X3, double X4);
extern "C" DLLEXPORT double Interpolation4(const vector<double> &X1s, const vector<double> &X2s, const vector<double> &X3s, const vector<double> &X4s, const vector<double> &Data, double X1, double X2, double X3, double X4);


/*
更新日志
V1.2~V1.2.1，20190325
修改了double类型的等于0判定
改进了表头长度为1时N维插值的处理方式

V1.2.1~V1.2.2，20190326
改进了InterpTableRead读取效率，减少了data长度次大小判断

V1.2.2~V1.2.3，20191202
增加了把读入的数据输出成标准气动txt的函数
*/