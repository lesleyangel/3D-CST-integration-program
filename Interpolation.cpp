//#include "StdAfx.h"
#include "interpolation.h"
#include<iostream>
#include<cmath>
#include<vector>
#include <sstream>
#include <iomanip>
#include "Interpolation.h"
#include <algorithm>
#include <functional>
using namespace std;

inline string& StrlTrimL(string &str)
{
	string::iterator p = find_if(str.begin(), str.end(), not1(ptr_fun<int, int>(isspace)));
	str.erase(str.begin(), p);
	return str;
}

inline string& StrTrimR(string &str) 
{
	string::reverse_iterator p = find_if(str.rbegin(), str.rend(), not1(ptr_fun<int , int>(isspace)));
	str.erase(p.base(), str.end());
	return str;
}

inline string& StrTrim(string &str)
{
	StrlTrimL(StrTrimR(str));
	return str;
}

inline vector<double> SplitString(string str,  bool withheader = true, string header = string() )
{
	istringstream istr(StrTrim(str));
	double dtemp = 0;
	vector<double> data;
	if (withheader)
		istr>>header;
	while(!istr.eof())
	{		istr>>dtemp;		data.push_back(dtemp);	}
	return (data);
}


Interpolation::Interpolation( void )
{
	IfOut = false;
	m_Dim = 0;
	Ifturn =false;
	//CreateZeroInterp();
}

Interpolation::~Interpolation( void )
{

}

bool Interpolation::SetHeader( vector<vector<double>> header,vector<string> headername )
{
	if (header.size() != headername.size())
	{
		return false;
	}
	if (0 == headername.size())
	{
		return false;
	}
	for(int i = 0; i < header.size(); i++)
	{
		if (0 == header[i].size())
		{
			return false;
		}
	}
	m_Dim = header.size();
	m_HeaderLength.resize(m_Dim);
	int sum = 1;
	for (int i = 0 ; i < m_Dim; i++)
	{
		m_HeaderLength[i] = header[i].size();
		sum = sum * m_HeaderLength[i];
	}
	m_Header = header;
	m_HeaderName = headername;
	m_data.resize(sum);
	
	return CheckTable();
}

bool Interpolation::SetData( vector<int> cubeorder, vector<vector<double>> cubedata )
{
	if ((cubeorder.size()+2) != m_Dim)//判断序号长度是否符合维数
	{
		return false;
	}
	if (1 == m_Dim)
	{
		if (1 != cubedata.size())//判断是否只有一行
		{
			return false;
		}
		if (cubedata[0].size() != m_HeaderLength[0])//判断列数与第一维长度是否相等
		{
			return false;
		}
	}
	else
	{
		for (size_t i = 1; i < cubedata.size(); i++)//判断数据是否为完整矩阵
		{
			if (cubedata[i-1].size() != cubedata[i].size())
			{
				return false;
			}
		}
		if ((int)cubedata.size() != m_HeaderLength[1])//判断行数与第二维长度是否相等
		{
			return false;
		}
		if ((int)cubedata[0].size() != m_HeaderLength[0])//判断列数与第一维长度是否相等
		{
			return false;
		}
	}
	//一维二维的情况
	if (0 == cubeorder.size())
	{
		for (size_t i = 0; i < cubedata.size(); i++)
		{
			for (size_t j = 0; j < cubedata[i].size(); j++)
			{
				m_data[j + i*cubedata[i].size()] = cubedata[i][j];
			}
		}
		return true;
	}
	//二维以上的情况
	int cubebegin = 0;
	int tempsum = 1;
	for (size_t i = 0; i < cubeorder.size(); i++ )
	{
		cubebegin = cubebegin + cubeorder[i] * tempsum;
		tempsum = tempsum * m_HeaderLength[i+2];
	}
	cubebegin = cubebegin * cubedata.size() * cubedata[0].size();
	for (size_t i = 0; i < cubedata.size(); i++)
	{
		for (size_t j = 0; j < cubedata[i].size(); j++)
		{
			m_data[cubebegin + j + i*cubedata[i].size()] = cubedata[i][j];
		}
	}
	return true;
}

vector<vector<double>> Interpolation::GetData(const vector<int>& cubeorder )
{
	vector<vector<double>> tempempty(0);
	if (cubeorder.size() > 0 && (cubeorder.size()+2) != m_Dim)
	{
		return tempempty;
	}
	if (cubeorder.size() == 0 && m_Dim > 2)
	{
		return tempempty;
	}
	if (0 == cubeorder.size() && 1 == m_Dim)
	{
		m_HeaderLength.resize(2);
		m_HeaderLength[1] = 1;
	}
	int cubebegin = 0;
	int tempsum = 1;
	for (int i = 0; i < cubeorder.size(); i++ )
	{
		cubebegin = cubebegin + cubeorder[i] * tempsum;
		tempsum = tempsum * m_HeaderLength[i+2];
	}
	cubebegin = cubebegin * m_HeaderLength[0] * m_HeaderLength[1];
	vector<double> temp(m_HeaderLength[0]);
	vector<vector<double>> cubedata(m_HeaderLength[1],temp);
	for (int i = 0; i < m_HeaderLength[1]; i++ )
	{
		for(int j = 0; j < m_HeaderLength[0]; j++)
		{
			cubedata[i][j] = m_data[cubebegin + j + i * m_HeaderLength[0]];
		}
	}
	return (cubedata);
}

bool Interpolation::InterpTableRead(string FileName)// 读取old插值表
{
	Ifturn = false;
	char		str1[4096];
	ifstream	fid;
	int			DataLength=1;
	m_HeaderLength.clear();		//	各维列表长度
	m_Header.clear();			//	各维列表
	m_data.clear();				//	数据
	m_HeaderName.clear();
	fid.open(FileName,ios::in);
	if ( fid)     //打开成功
	{	
		fid >> str1;  //字符串 维度：
		fid >> m_Dim; //维数值
		if (m_Dim == 0 )
		{			cout<<"interpdata error with file "<<FileName<<endl;		return false;		}
		m_HeaderLength.resize( m_Dim );
		for ( int id=0; id!=m_Dim; id++ )
		{
			fid >> str1;
			fid	>> m_HeaderLength[id];	  //读入各表头的长度
		}
		m_HeaderName.resize( m_Dim );
		m_Header.resize( m_Dim );

		for ( int id=0; id!=m_Dim; id++ )
		{
			m_Header[id].resize( m_HeaderLength[id] );
			fid >> m_HeaderName[id];
			for ( int i=0; i!=m_HeaderLength[id]; i++ )
			{
				fid >> m_Header[id][i];  //读入各维表头				
			}
			DataLength=DataLength * m_HeaderLength[id];  //记录各维离散点个数，即为气动数据的总个数
		}

		m_data.resize( DataLength );

		if (m_Dim <= 2)//一维和二维情况
		{
			fid >>str1;  fid.getline(str1,4096);//读入表头行
			for ( int i=0; i!=DataLength; i++ )
			{
				if (i%m_HeaderLength[0]==0)//每次第一列先跳过竖直表头
				{
					fid>>str1;
				}			
				fid >> m_data[i];			
			}
		} 
		else//三维及以上情况
		{
			for ( int i=0; i!=DataLength; i++ )
			{
				if ((i%(m_HeaderLength[0]*m_HeaderLength[1]))==0)//每组数据开始时
				{
					fid >> str1;  fid.getline(str1,4096);   //读入注释行
					fid >> str1;  fid.getline(str1,4096);	//读入表头行
				}
				if (i%m_HeaderLength[0]==0)//每次第一列先跳过竖直表头
				{
					fid>>str1;
				}			
				fid >> m_data[i];			
			}
		}

		fid.close();
		return CheckTable();
	}
	else      //打开文件失败
	{
		return false;
	}
}

bool Interpolation::InterpTableReadFile(string filename)// 读取插值表
{
	Ifturn = false;
	m_HeaderLength.clear();		//	各维列表长度
	m_Header.clear();			//	各维列表
	m_data.clear();				//	数据
	m_HeaderName.clear();
	int DataLength=1;
	//
	ifstream fid(filename, ios::in);
	if (!fid)
	{
		cout<<"error open interp file "<<filename<<endl;
		return false;
	}
	string templine;
	vector<double> tempdata;
	// dim
	getline(fid, templine);  //info
	getline(fid, templine);  
	tempdata = SplitString(templine);
	if (1!=tempdata.size())
	{			cout<<"error interpdata dim,  error with file "<<filename<<endl;		return false;		}
	m_Dim = tempdata[0];
	if (m_Dim == 0 )
	{			cout<<"error interpdata dim set 0,  error with file "<<filename<<endl;		return false;		}
	// header
	getline(fid, templine);  //info
	m_HeaderLength.resize(m_Dim);
	m_Header.resize(m_Dim);
	m_HeaderName.resize(m_Dim);
	for (size_t i=0; i<m_Dim; i++)
	{
		getline(fid, templine);  //info
		m_Header[i] = SplitString(templine, true, m_HeaderName[i]);
		m_HeaderLength[i] = m_Header[i].size();
		if (0==m_HeaderLength[i])
			{			cout<<"error interpdata header dim set 0,  error with file "<<filename<<endl;		return false;		}		
		DataLength = DataLength * m_HeaderLength[i];  //记录各维离散点个数，即为气动数据的总个数
	}
	// data
	getline(fid, templine);  //info
	m_data.resize(DataLength);	
	if (m_Dim <= 2)//一维和二维情况
	{
		getline(fid, templine);  //读入表头行
		for ( int i=0; i!=DataLength; i++ )
		{
			if (i%m_HeaderLength[0]==0)//每次第一列先跳过竖直表头
			{				fid>>templine;			}			
			fid >> m_data[i];			
		}
	} 
	else//三维及以上情况
	{
		int curpos = 0;
		int outsize = DataLength/m_HeaderLength[0]/m_HeaderLength[1];
		for (int iout = 0;iout<outsize;iout++)
		{
			getline(fid, templine); //读入辅助表头行
			getline(fid, templine);	//读入表头行
			for (int i = 0; i<m_HeaderLength[1];i++)
			{
				fid >> templine;
				for (int j = 0; j<m_HeaderLength[0];j++)
					fid >> m_data[curpos++];
			}
			getline(fid, templine);  //info  // 跳过换行符
		}
	}
	return true;
}

bool Interpolation::InterpTableReadData(string strdata)
{
	Ifturn = false;
	m_HeaderLength.clear();		//	各维列表长度
	m_Header.clear();			//	各维列表
	m_data.clear();				//	数据
	m_HeaderName.clear();
	int DataLength=1;
	//
	istringstream fid(strdata);
	if (!fid)
	{
		cout<<"error open interp strdata"<<endl;
		return false;
	}
	string templine;
	vector<double> tempdata;
	// dim
	getline(fid, templine);  //info
	getline(fid, templine);  
	tempdata = SplitString(templine);
	if (1!=tempdata.size())
	{			cout<<"error interpdata dim"<<endl;		return false;		}
	m_Dim = tempdata[0];
	if (m_Dim == 0 )
	{			cout<<"error interpdata dim set 0"<<endl;		return false;		}
	// header
	getline(fid, templine);  //info
	m_HeaderLength.resize(m_Dim);
	m_Header.resize(m_Dim);
	m_HeaderName.resize(m_Dim);
	for (size_t i=0; i<m_Dim; i++)
	{
		getline(fid, templine);  //info
		m_Header[i] = SplitString(templine, true, m_HeaderName[i]);
		m_HeaderLength[i] = m_Header[i].size();
		if (0==m_HeaderLength[i])
		{			cout<<"error interpdata header dim set 0"<<endl;		return false;		}
		DataLength = DataLength * m_HeaderLength[i];  //记录各维离散点个数，即为气动数据的总个数
	}
	// data
	getline(fid, templine);  //info
	m_data.resize(DataLength);	
	if (m_Dim <= 2)//一维和二维情况
	{
		getline(fid, templine);  //读入表头行
		for ( int i=0; i!=DataLength; i++ )
		{
			if (i%m_HeaderLength[0]==0)//每次第一列先跳过竖直表头
			{				fid>>templine;			}			
			fid >> m_data[i];			
		}
	} 
	else//三维及以上情况
	{
		int curpos = 0;
		int outsize = DataLength/m_HeaderLength[0]/m_HeaderLength[1];
		for (int iout = 0;iout<outsize;iout++)
		{
			getline(fid, templine); //读入辅助表头行
			getline(fid, templine);	//读入表头行
			for (int i = 0; i<m_HeaderLength[1];i++)
			{
				fid >> templine;
				for (int j = 0; j<m_HeaderLength[0];j++)
					fid >> m_data[curpos++];
			}
			getline(fid, templine);  //info  // 跳过换行符
		}
	}
	return true;
}

void Interpolation::InterpTableRead2(const char *FileName,int num)
{
	Ifturn = false;
	char str[4096];
	int info1,info2;
	m_HeaderLength.clear();		//	各维列表长度
	m_Header.clear();//	各维列表
	m_data.clear();			//	数据
	ifstream fid ( FileName ); 
	fid >> str;
	fid >> info1;
	fid >> info2;
	vector<double> a(info2,0);
	vector<vector<double>> Data(info1,a);
	for (int i = 0; i < 11; i++)
	{
		fid.getline(str,4096);
	}
	fid >> str;
	fid >> m_Dim;
	m_Header.resize(m_Dim);
	m_HeaderLength.resize(m_Dim);
	int sum=1;
	for (int i=0; i<m_Dim; i++)
	{
		fid >> str;
		fid >> m_HeaderLength[i];
		sum=sum*m_HeaderLength[i];
	}
	m_data.resize(sum);
	for (int i=0; i<m_Dim; i++)
	{
		m_Header[i].resize(m_HeaderLength[i]);
		fid >> str;
		for (int j=0; j<m_HeaderLength[i]; j++)
		{
			fid >> m_Header[i][j];
		}
	}
	fid.getline(str,4096);
	fid.getline(str,4096);//这下面开始读data
	
	for (int i=0; i<info1; i++)
	{
		for (int j=0; j<info2; j++)
		{
			fid >> Data[i][j];
		}
	}
	a.resize(sum);
	vector<vector<double>> AllData(info2,a);
	for (int i=0; i<info1; i++)
	{
		vector<double> head(m_Dim,0);
		for (int j=0; j<m_Dim; j++)
		{
			head[j] = Data[i][j];
		}
		int n = FindPosition(head);
		if (n != -1)
		{
			for (int j=0; j<info2; j++)
			{
				AllData[j][n] = Data[i][j];
			}
		}
	}
	m_data.resize(sum);
	m_data = AllData[num];

}

int Interpolation::FindPosition(const vector<double>& head)
{
	int Position = 0;
	int sum=1;
	for (int i=0; i<m_Dim; i++)
	{
		double tempmax = fabs(m_Header[i][0]);
		if (fabs(m_Header[i][m_HeaderLength[i]-1])>tempmax)//寻找表头中最大的绝对值
		{
			tempmax=fabs(m_Header[i][m_HeaderLength[i]-1]);
		}
		double Small = tempmax/1000.0;//相对小量
		if (fabs(Small) < 0.00000001)//如果表头只有一个0，则相对小量为0.001
		{
			Small = 0.001;
		}
		int id = -1;
		for (int j=0; j<m_HeaderLength[i]; j++)//寻找当前表头值对应表头序号
		{
			if (fabs(m_Header[i][j]-head[i])<Small)
			{
				id=j;
				break;
			}
		}
		if (id == -1)//如果没有找到相等的，则返回-1
		{
			return -1;
		}
		Position = Position + sum * id;
		sum=sum*m_HeaderLength[i];
	}
	return Position;
}

double Interpolation::Interp1( double x1 , double x2 , double y1 , double y2 , double x)//一维插值计算
{
	if (fabs(x1-x2)<DBL_EPSILON)
	{
		return y1;
	}
	if (IfOut) //允许外插
	{
		return (x-x1) * (y2-y1) / (x2-x1) + y1;
	} 
	else //不允许外插
	{
		if (x < x1)
		{
			return y1;
		} 
		else if (x > x2)
		{
			return y2;
		}
		else
		{
			return (x-x1) * (y2-y1) / (x2-x1) + y1;
		}
	}
}

double Interpolation::InterpN(int start, const vector<double>& x , int Dim , const vector<vector<int>>& XNum , const vector<vector<double>>& XVal)//n维插值计算//date数据，插入点，维数，插入点两侧表头序号，插入点两侧表头值
{
	int sum=1;
	for( int i=0; i<Dim-1; i++)//求一组date的总数
		sum *= m_HeaderLength[i];
	double y=0.0;
	if (Dim>1)
	{
		double y1 = InterpN(start + sum*XNum[0][Dim-1],  x, Dim-1, XNum, XVal);
		if (XNum[1][Dim-1] == 0)
		{
			y = y1;
		} 
		else
		{
			double y2 = InterpN(start + sum*XNum[1][Dim-1],  x, Dim-1, XNum, XVal);//两组date压缩成两点
			y = Interp1(XVal[0][Dim-1], XVal[1][Dim-1], y1, y2, x[Dim-1]);		
		}
	}
	else //一维的情况
	{
		y = Interp1(XVal[0][0], XVal[1][0], m_data[start + sum*XNum[0][Dim-1]], m_data[start + sum*XNum[1][Dim-1]], x[0]);
	}
	return y;
}

double Interpolation::InterpolationN(const vector<double>& xtemp, bool out )	//n维插值前处理及结果
{	
	// 为空时返回0
	if (0 == m_Dim)
		return 0;
	vector<double> x = xtemp;
	if (Ifturn)//如果是ReadAreo读取的数据，需要将x逆序
	{
		for (int i=0; i!=xtemp.size() ;i++)
			x[xtemp.size()-i-1]=xtemp[i];
	}
	IfOut=out;
	vector<int> Num(m_Dim,0);
	vector<double> Val(m_Dim,0.0);
	vector<vector<int>> XNum(2,Num);//存放插入点两侧的表头序号
	vector<vector<double>> XVal(2,Val);//存放插入点两侧的表头的值
	for( int i=0; i<m_Dim; i++)
	{
		if (m_HeaderLength[i] <= 1)
		{
			XNum[0][i] = 0;
			XNum[1][i] = 0;
			XVal[0][i] = m_Header[i][0];
			XVal[1][i] = m_Header[i][0];
		}
		for ( int j=0; j<m_HeaderLength[i]-1; j++)
			{				
				XNum[0][i] = j;
				XNum[1][i] = j+1;
				XVal[0][i] = m_Header[i][j];
				XVal[1][i] = m_Header[i][j+1];
				if (( m_Header[i][j]<= x[i] && x[i] <= m_Header[i][j+1])||( m_Header[i][j]>= x[i] && x[i] >= m_Header[i][j+1]))
					break;
				if ( x[i] <= m_Header[i][0] && m_Header[i][0] < m_Header[i][1] )
					break;
				if ( x[i] >= m_Header[i][0] && m_Header[i][0] > m_Header[i][1] )
					break;			
				else 
					continue;
			}
	}//以上完成了求插入点两侧表头序号和值
	return InterpN(0, x, m_Dim, XNum, XVal);//插值结果
}

int Interpolation::GetDim()									//读取维度
{
	return m_Dim;
}

const vector<int>& Interpolation::GetHeaderLength()					//读取表头长度
{
	return m_HeaderLength;
}

const vector<vector<double>>& Interpolation::GetHeader()			//读取表头
{
	return m_Header;
}

bool Interpolation::CheckTable()									//检查表格数据是否异常
{
	int DateLength=1;
	for ( int id=0; id<m_Dim; id++)
	{
		DateLength=DateLength*m_HeaderLength[id];
		for ( int i=0; i<m_HeaderLength[id]-1; i++)
		{
			if ( m_Header[id][i] >= m_Header[id][i+1] )
			{
				return false;
			}
		}
	}
	if ( DateLength != m_data.size() )
	{
		return false;
	}
	return true;
}



vector<double> Interpolation::Shorten(const vector<double>& a)
{
	int n = a.size();
	vector<double> b;
	b.push_back(a[0]);
	bool turn = false;
	for (int i=1; i<n; i++)
	{
		if (turn && fabs(a[i]-a[0])<DBL_EPSILON)
		{
			break;
		}
		if (fabs(a[i]-a[i-1])>DBL_EPSILON)
		{
			b.push_back(a[i]);
			turn = true;
		}
	}
	return (b);

}

void Interpolation::ReadAero(const char *FileName,int n)
{
	Ifturn = true;
	char str[4096];
	int info1,info2;
	int Dim;
	m_HeaderLength.clear();		//	各维列表长度
	m_Header.clear();//	各维列表
	m_data.clear();			//	数据
	vector<vector<double>> Header;
	ifstream fid ( FileName ); 
	fid >> str;
	fid >> info1;
	fid >> info2;
	vector<double> a(info1,0);
	vector<vector<double>> Data(info2,a);
	for (int i = 0; i < 11; i++)
	{
		fid.getline(str,4096);
	}
	fid >> str;
	fid >> Dim;
	Header.resize(Dim);

	fid.getline(str,4096);
	fid.getline(str,4096);
	
	for (int i=0; i<info1; i++)
	{
		for (int j=0; j<info2; j++)
		{
			fid >> Data[j][i];
		}
	}
	for (int i=0;i<Dim;i++)
	{
		Header[i]=Shorten(Data[i]);
	}
	m_Dim = Dim;
	m_data.resize(info1);
	m_data = Data[n];
	m_HeaderLength.resize(Dim);
	m_Header.resize(Dim);
	for (int i=0;i<Dim;i++)
	{
		m_HeaderLength[i]=Header[Dim-i-1].size();
		m_Header[i].resize(Header[Dim-i-1].size());
		m_Header[i]=Header[Dim-i-1];
	}
}

bool Interpolation::WriteTxt(string txtpath)
{
	if (0 == m_data.size())
		return false;

	fstream ff;
	ff.open(txtpath,ios::out);
	ff.precision(12);
	ff.flags(ios::left);
	//
	ff<<"%格式说明：(1)采用多个二维表表示任意高维数据；(2)请严格遵守数据顺序；(3)不允许出现空白行；(4)数据请采用空格隔开；(5)%开头为注释行，请不要删除"<<endl;
	ff<<setw(20)<<"维数维度(1-any)"<<setw(20)<<m_Dim<<endl;
	//
	ff<<"%下述表头必须和维度和各维度表头维数完全对应，增序排列，标识符支持(攻角、侧滑角、马赫数、高度、当量比、动压)"<<endl;
	for (int i = 0; i < m_Dim ;i++ )
	{
		ff<<setw(20)<<m_HeaderName[i];
		for (int j = 0; j < m_HeaderLength[i]; j++)
			ff<<setw(20)<<m_Header[i][j];
		ff<<endl;
	}
	//
	ff<<"%下述为插值数据，按照表头逆序定义成多个顺序排列二维表，一维二维情况不需要高维辅助表头"<<endl;
	if (1==m_Dim)
	{
		ff<<setw(20)<<m_HeaderName[0];
		for (int j = 0; j<m_HeaderLength[0]; j++)
			ff<<setw(20)<<m_Header[0][j];
		ff<<endl;
		for (int j = 0; j<m_HeaderLength[0]; j++)
			ff<<setw(20)<<m_data[j];
	}
	else if (2==m_Dim)
	{
		ff<<setw(20)<<m_HeaderName[1] + "\\" + m_HeaderName[0];
		for (int j = 0; j<m_HeaderLength[0]; j++)
			ff<<setw(20)<<m_Header[0][j];
		ff<<endl;
		for (int i = 0; i<m_HeaderLength[1]; i++)
		{
			ff<<setw(20)<<m_Header[1][i];
			for (int j = 0; j<m_HeaderLength[0]; j++)
				ff<<setw(20)<<m_data[j + i * m_HeaderLength[0]];
			ff<<endl;
		}
	}
	else if (3<=m_Dim)
	{
		int sum = 1;
		for (int i = 2; i<m_Dim; i++)
			sum = sum * m_HeaderLength[i];
		for (int i = 0; i<sum; i++ )
			WriteCube(i,ff);
	}
	ff.close();
	return true;
}

void Interpolation::WriteCube( int n, fstream &ff)
{
	if (!ff)
		return;
	vector<int> order(m_Dim-2);
	int tempsum = 1;
	int cubbegin = n * m_HeaderLength[0]*m_HeaderLength[1];
	for (int i = 0; i < m_Dim-2; i++ )
	{
		order[i] = (n/tempsum)%m_HeaderLength[i+2];
		tempsum = tempsum * m_HeaderLength[i+2];
	}
	ff<<"["<<" ";
	for (int i = 0; i<m_Dim - 2; i++)
		ff<<m_HeaderName[i+2]<<" = "<<m_Header[i+2][order[i]]<<" ";
	ff<<"]"<<endl;
	ff<<setw(20)<<m_HeaderName[1] + "\\" + m_HeaderName[0];
	for (int j = 0; j<m_HeaderLength[0]; j++)
		ff<<setw(20)<<m_Header[0][j];
	ff<<endl;
	for (int i = 0; i<m_HeaderLength[1]; i++)
	{
		ff<<setw(20)<<m_Header[1][i];
		for (int j = 0; j<m_HeaderLength[0]; j++)
			ff<<setw(20)<<m_data[cubbegin + j + i * m_HeaderLength[0]];
		ff<<endl;
	}
}

void Interpolation::CreateZeroInterp()
{
	m_Dim = 1;
	m_HeaderLength.resize(1,1);
	vector<double> temp(1);
	m_Header.resize(1,temp);
	m_HeaderName.resize(1,"马赫数");
	m_data.resize(1,0);
}

int Interpolation::Check()
{
	if (m_Dim<=0 || (int)m_data.size()<=0 || (int)m_Header.size()<=0 || (int)m_HeaderLength.size()<=0)
	{
		cout<<"interp data size equal to zero"<<endl;
		return -1;
	}
	if (m_Dim!=m_HeaderLength.size() || m_HeaderLength.size()!=m_Header.size())
	{
		cout<<"interp data size not equal"<<endl;
		return -2;
	}
	if (!CheckTable())
	{
		cout<<"interp data struct error"<<endl;
		return -3;
	}
	return 0;
}

const vector<string>& Interpolation::GetHeaderName()
{
	return m_HeaderName;
}

////////////////////////////////////////////////////////////////////////// 
// 更高效的插值代码，但是数据需要额外处理
//
//////////////////////////////////////////////////////////////////////////
double Interpolation1( int n_Xs, double *Xs,double *Data, double X )
{	
	// Step1：确定Xs向位置
	int id=0;
	for ( int i=0; i!=n_Xs-1; i++ )
	{
		id=i;
		if ( X>=Xs[i] && X<=Xs[i+1] )			
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double ratio=0;

	if ( X<Xs[0] )
	{
		ratio=(Xs[0]-X)/(Xs[1]-Xs[0]);
		return Data[0]-ratio*(Data[1]-Data[0]);
	}
	else if ( X>=Xs[n_Xs-1] )
	{
		ratio=(X-Xs[n_Xs-1])/(Xs[n_Xs-1]-Xs[n_Xs-2]);
		return Data[n_Xs-1]+ratio*(Data[n_Xs-1]-Data[n_Xs-2]);
	}
	else 
	{
		ratio=(X-Xs[id])/(Xs[id+1]-Xs[id]);
		return Data[id]+(Data[id+1]-Data[id])*ratio;
	}
}

double Interpolation1(const vector<double> &Xs, const vector<double>& Data, double X )
{
	int n_Xs = Xs.size();
	double tem=0;
	// Step1：确定Xs向位置
	int id =0;
	for ( int i=0; i!=n_Xs-1; i++ )
	{
		id=i;
		if ((Xs[i] <= X && X <= Xs[i+1])||(Xs[i] >= X && X >= Xs[i+1]))			
			break;
		if ( X<Xs[0] && Xs[0]<Xs[1] )
			break;
		if ( X>Xs[0] && Xs[0]>Xs[1] )
			break;
		else
			continue;
	}
	// Step2：Xs向插值
	return Data[id]+(Data[id+1]-Data[id])*(X-Xs[id])/(Xs[id+1]-Xs[id]);
}

double Interpolation2( int n_Xs, double *Xs, int n_Ys, double *Ys, double *Data, double X, double Y)
{	
	// Step1：确定Xs向位置
	int id = 0;
	for ( int i=0; i!=n_Xs-1; i++ )
	{
		id=i;
		if ( X>=Xs[i] && X<=Xs[i+1] )			
			break;
		else
			continue;
	}

	// Step2：Xs向插值
	double *Y_Data= new double[n_Ys];
	double ratio;

	if ( X<Xs[0] )
	{
		ratio=(Xs[0]-X)/(Xs[1]-Xs[0]);
		for ( int i=0; i!=n_Ys; i++ )
			Y_Data[i]=Data[i*n_Xs]-ratio*(Data[i*n_Xs+1]-Data[i*n_Xs]);
	}
	else if ( X>=Xs[n_Xs-1] )
	{
		ratio=(X-Xs[n_Xs-1])/(Xs[n_Xs-1]-Xs[n_Xs-2]);
		for ( int i=0; i!=n_Ys; i++ )
			Y_Data[i]=Data[i*n_Xs+n_Xs-1]+ratio*(Data[i*n_Xs+n_Xs-1]-Data[i*n_Xs+n_Xs-2]);
	}
	else 
	{
		ratio=(X-Xs[id])/(Xs[id+1]-Xs[id]);
		for ( int i=0; i!=n_Ys; i++ )
			Y_Data[i]=Data[i*n_Xs+id]+(Data[i*n_Xs+id+1]-Data[i*n_Xs+id])*ratio;
	}

	// Step3：确定Ys向位置
	for ( int i=0; i!=n_Ys-1; i++ )
	{
		id=i;
		if ( Y>=Ys[i] && Y<=Ys[i+1] )
			break;
		else
			continue;
	}

	// Step4：Ys向插值
	double	ResultVal = 0;
	if ( Y<Ys[0] )
	{
		ratio=(Ys[0]-Y)/(Ys[1]-Ys[0]);
		ResultVal = Y_Data[0]-ratio*(Y_Data[1]-Y_Data[0]);
	}

	else if ( Y>=Ys[n_Ys-1] )
	{
		ratio=(Y-Ys[n_Ys-1])/(Ys[n_Ys-1]-Ys[n_Ys-2]);
		ResultVal = Y_Data[n_Ys-1]+ratio*(Y_Data[n_Ys-1]-Y_Data[n_Ys-2]);
	}
	else
	{
		ratio=(Y-Ys[id])/(Ys[id+1]-Ys[id]);
		ResultVal = Y_Data[id]+(Y_Data[id+1]-Y_Data[id])*ratio;
	}
	delete [] Y_Data;
	return ResultVal;
}

double Interpolation2(const vector<double> &Xs, const vector<double> &Ys, const vector<double> &Data, double X, double Y)
{	
	int n_Xs = Xs.size();
	int n_Ys = Ys.size();
	int id=0;
	// Step1：确定Ys向位置
	for ( int i=0; i!=n_Ys-1; i++ )
	{
		id=i;
		if ((Ys[i] <= Y && Y <= Ys[i+1])||(Ys[i] >= Y && Y >= Ys[i+1]))		
			break;
		if ( Y<Ys[0] && Ys[0]<Ys[1] )
			break;
		if ( Y>Ys[0] && Ys[0]>Ys[1] )
			break;
		else
			continue;
	}
	// Step2: 由一维插值得出两个点
	vector<double> data1(n_Xs,0);
	vector<double> data2(n_Xs,0);
	for ( int i=0; i!=n_Xs; i++ )
	{
		data1[i]=Data[i+id*n_Xs];
		data2[i]=Data[i+(id+1)*n_Xs];
	}
	double Y1=Interpolation1(Xs,data1,X);
	double Y2=Interpolation1(Xs,data2,X);
	// Step3：Ys向插值
	return Y1+(Y2-Y1)*(Y-Ys[id])/(Ys[id+1]-Ys[id]);
}


//三维差值
double Interpolation3( int n_Xs, double *Xs, int n_Ys, double *Ys, int n_Zs, double *Zs, double *Data, double X, double Y, double Z )
{
	// Step1：Xs、Ys向二维插值
	double *Z_Data= new double[n_Zs];

	for ( int i=0; i!=n_Zs; i++ )
		Z_Data[i]=Interpolation2(n_Xs,Xs,n_Ys,Ys,Data+i*n_Xs*n_Ys,X,Y); 

	// Step2：确定Zs向位置
	int id = 0;
	for ( int i=0; i!=n_Zs-1; i++ )
	{
		id=i;
		if ( Z>=Zs[i] && Z<=Zs[i+1] )
			break;
		else
			continue;
	}


	// Step3：Zs向插值
	double ResultVal = 0;
	double ratio = 0;
	if ( Z<Zs[0] )
	{
		ratio=(Zs[0]-Z)/(Zs[1]-Zs[0]);
		ResultVal = Z_Data[0]-ratio*(Z_Data[1]-Z_Data[0]);
	}
	else if ( Z>=Zs[n_Zs-1] )
	{
		ratio=(Z-Zs[n_Zs-1])/(Zs[n_Zs-1]-Zs[n_Zs-2]);
		ResultVal = Z_Data[n_Zs-1]+ratio*(Z_Data[n_Zs-1]-Z_Data[n_Zs-2]);
	}

	else
	{		
		ratio=(Z-Zs[id])/(Zs[id+1]-Zs[id]);
		ResultVal = Z_Data[id]+(Z_Data[id+1]-Z_Data[id])*ratio;
	}
	delete[] Z_Data;
	return ResultVal;
}

double Interpolation3(const vector<double>& Xs, const vector<double>& Ys, const vector<double> &Zs, const vector<double> &Data, double X, double Y, double Z )

{
	int n_Xs = Xs.size();
	int n_Ys = Ys.size();
	int n_Zs = Zs.size();
	// Step1：Xs、Ys向二维插值
	vector<double> Z_Data(n_Zs,0.0);
	vector<double> New_Data;
	for ( int i=0; i!=n_Zs; i++ )
	{
		New_Data.assign(Data.begin()+i*n_Xs*n_Ys,Data.end());
		Z_Data[i]=Interpolation2(Xs,Ys,New_Data,X,Y); 
	}    

	// Step2：确定Zs向位置
	int id = 0;
	for ( int i=0; i!=n_Zs-1; i++ )
	{
		id=i;
		if (Z<Zs[0])
		{
			break;
		}
		if ( Z>=Zs[i] && Z<=Zs[i+1] )
			break;
		else
			continue;
	}


	// Step3：Zs向插值
	double ResultVal;
	double ratio;
	if ( Z<Zs[0] )
	{
		ratio=(Zs[0]-Z)/(Zs[1]-Zs[0]);
		ResultVal = Z_Data[0]-ratio*(Z_Data[1]-Z_Data[0]);
	}
	else if ( Z>=Zs[n_Zs-1] )
	{
		ratio=(Z-Zs[n_Zs-1])/(Zs[n_Zs-1]-Zs[n_Zs-2]);
		ResultVal = Z_Data[n_Zs-1]+ratio*(Z_Data[n_Zs-1]-Z_Data[n_Zs-2]);
	}
	else
	{		
		ratio=(Z-Zs[id])/(Zs[id+1]-Zs[id]);
		ResultVal = Z_Data[id]+(Z_Data[id+1]-Z_Data[id])*ratio;
	}
	return ResultVal;
}

//四维差值，采用一个三维差值组合一维差值实现
//一、二、三维维组成一组三维数据，进行差值，
//结果和第四维构成一维数据
double Interpolation4( int n_X1, double *X1s, int n_X2, double *X2s, int n_X3, double *X3s, int n_X4, double *X4s, double *Data, double X1, double X2, double X3, double X4)
{
	double w; //返回值

	double	*Temp=new double [n_X4];
	int		Location; //三维差值表内数据的总个数

	for(int i = 0;i != n_X4; i++)
	{
		Location=n_X1 * n_X2 * n_X3 * i; //各组三维差值表起始位置相对Data的偏移量

		Temp[i]=Interpolation3(n_X1,X1s,n_X2,X2s,n_X3,X3s,Data+Location,X1,X2,X3);//调用三维差值，计算给定第一、二、四维需差值数据，在第三维节点上的气动数据
	}
	w = Interpolation1(n_X4,X4s,Temp,X4);
	delete[] Temp;
	return(w);
}

double Interpolation4(const vector<double> &X1s,  const vector<double> &X2s, const vector<double> &X3s, const vector<double> &X4s, const vector<double> &Data, double X1, double X2, double X3, double X4)
{
	int n_X1 = X1s.size();
	int n_X2 = X2s.size();
	int n_X3 = X3s.size();
	int n_X4 = X4s.size();
	double w; //返回值
	vector<double> Temp(n_X4,0.0);
	int		Location; //三维差值表内数据的总个数
	vector<double> New_Data;
	for(int i = 0;i != n_X4; i++)
	{
		Location=n_X1 * n_X2 * n_X3 * i; //各组三维差值表起始位置相对Data的偏移量
		New_Data.assign(Data.begin()+Location,Data.end());
		Temp[i]=Interpolation3(X1s,X2s,X3s,New_Data,X1,X2,X3);//调用三维差值，计算给定第一、二、四维需差值数据，在第三维节点上的气动数据
	}
	w = Interpolation1(X4s,Temp,X4);
	return(w);
}