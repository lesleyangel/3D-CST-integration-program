#include"CST_Instantiation.h"
#include"Aircraft.h"
#include"Curve.h"
#include <sstream> 
using namespace arma;
using namespace std;

#define pi 3.1415926


mat copyMall(vector<vector<double>> &v)
{
	vector<vector<double>>::iterator s = v.begin();
	mat m = zeros(v.size(),(*s).size());
	int i = 0;
	//int j = 0;
	for (vector<vector<double>>::iterator it = v.begin(); it != v.end(); it++)
	{
		int j = 0;
		//(*it)---容器 vector<int>
		for (vector<double>::iterator vit = (*it).begin(); vit != (*it).end(); vit++)
		{
			//cout << *vit << " ";
			m(i, j) = *vit;
			j++;
		}
		i++;
	}
	return m;
}

//test CopyMall
void test05()
{
	vector<vector<double>>v;

	vector<double>v1;
	vector<double>v2;
	vector<double>v3;
	for (int i = 0; i < 1; i++)
	{
		v1.push_back((double)i);
		v2.push_back((double)i + 3);
		v3.push_back((double)i + 6);
	}
	v.push_back(v1);
	//v.push_back(v2);
	//v.push_back(v3);

	cout << "v1.size:" << v1.size() << endl;
	vector<vector<double>>::iterator ss = v.begin();
	cout << v.size() << " " << (*ss).size() << endl;
	
	mat AA = copyMall(v);
	cout << AA << endl;

	vector<vector<double>>BUPP;
	vector<double> va;
	//va.resize(1);
	va.push_back(1);

	//BUPP.resize(1);
	BUPP.push_back(v1);

	cout << "va.size:" << va.size() << endl;
	vector<vector<double>>::iterator s = BUPP.begin();
	cout << BUPP.size() << " " << (*s).size() << endl;
	AbstructShape::copyM(v, AA);
}

//test Aircraft
void testAircraft()
{
	SMB1 mb1;
	mb1.NFaiL = mb1.NFaiU = mb1.NEta = mb1.NHeight = 40;
	mb1.NEta = 20;
	mb1.N2 = { 0.1,0.1,0.5,0.5 };
	mb1.N1 = { 3.0,3.0,0.1,0.1 };
	// mb1.N2 = { 0.1,0.1,0.5,0.5 };
	// mb1.N1 = { 0.1,0.1,0.5,0.5 };
	// mb1.N2 = { 0.5,1.0,0.5,1.0 };
	// mb1.N1 = { 0.5,1.0,0.5,1.0 };
	//mb1.Length = {15000, 4000, 4000, 8000};
	mb1.Length = { 150, 40, 40, 80 };
	mb1.M1 = { 0.0,0,0.0,0 };
	mb1.T1 = { 0.0,0,0.0,0 };
	//mb1.M1 = { 0.0,0.1,0.0,0.1 };
	//mb1.T1 = { 0.0,0.1,0.0,0.1 };
	//mb1.M1 = { 0.5,0.0,0.5,0.0 };
	//mb1.T1 = { 0.5,0.0,0.5,0.0 };
	mb1.GridRefineType = -1;
	SMB2 mb2;
	mb2.NFaiL = mb2.NFaiU = 40;  mb2.NHeight = 4;
	mb2.N1 = { 0.2,0.2,0.1,0.1 };
	mb2.N2 = { 0.4,0.4,0.1,0.1 };
	mb2.N3 = { 0.1,0.1,0.1,0.1 };
	mb2.LHead[0] += 1000;
	mb2.LBody[0] -= 1000;
	mb2.LHead[0] *= 1.2;
	mb2.LBody[0] *= 1.2;
	mb2.LHead[3] = mb2.LBody[3] = 3500;
	mb2.LBody[1] = mb2.LHead[1] = 2300;
	
	//mb2.BUPP1 = { { 10,1,10,1,10 },
	//			{ 10,1,10,1,10 }, 
	//			{ 10,1,10,1,10 },
	//			{ 10,1,10,1,10 },
	//			{ 10,1,10,1,10 },
	//			{ 10,1,10,1,10 } };
	//mb2.BUPP2 = { { 1.5,.5,1.5,9,1,1 },
	//			{ 1.5,.5,1.5,9,1,1 } ,
	//			{ 1.5,.5,1.5,9,1,1 } ,
	//			{ 1.5,.5,1.5,9,1,1 } ,
	//			{ 1.5,.5,1.5,9,1,1 } ,
	//			{ 1.5,.5,1.5,9,1,1 } };
	//mb2.DUPP1 = { {1,2,1,2,1} };
	
	//mb2.NFaiL =  mb2.NFaiU = mb2.NEta = mb2.NHeight = 2;
	SMB3 mb3;
	//变形1
	//mb3.Length[0] *= 2.5;
	//mb3.Length[4] *= 2;
	//mb3.N1 = { 0.5,0.5,0.5,0.5 };
	//mb3.N4 = mb3.N1;
	//变形2
	//mb3.BUPP1 = { {1,2,1,2.2,2} };
	//mb3.BUPP2 = { {2,2,1,2,1,2,2} };
	//mb3.BUPP3 = { {2,2,1,2,1} };
	//变形3
	mb3.DUPP1 = mb3.DLOW1 = { {1.2,.7,1.2,.7,.7} } ;
	mb3.DUPP2 = mb3.DLOW2 = { {.7,.7,1.2,.7,1.2,.7,.7} };
	mb3.DUPP3 = mb3.DLOW3 = { {.7,.7,1.2,.7,1.2} };

	SW1 sw1; 
	sw1.Origin = { 16500*1.2,-900,1830 };
	sw1.SpanL = 5020*1.2;
	sw1.RootChordL = 6500*1.2;
	sw1.TipRootRatio = 0.268;
	sw1.Thickness = 600;
	sw1.NFaiL = sw1.NFaiU = sw1.NEta = 40;
	sw1.NHeight = 4;
	sw1.ifUseMid = false;
	
	SW2 sw2;
	sw2.NFaiL = sw2.NFaiU = 30;
	sw2.NEta_inn = sw2.NEta_out = 30;
	sw2.NHeight = 4;
	sw2.SpanL_inn = 3000;
	sw2.SweepBackAngle_inn = 65;
	sw2.ifUseMid = false;
	
	ST1 st1;
	st1.Origin = { 24000,2330,0 };
	st1.NFaiL = st1.NFaiU = st1.NEta = 40; st1.NHeight = 3;
	st1.RootChordL = 3500;
	st1.SpanL = 3500;
	st1.SweepBackAngle = 30;
	
	ST2 st2;
	st2.NFaiL = st2.NFaiU = st2.NEta = 10; st2.NHeight = 3;
	
	Aircraft A;								//实例化一个部件
	//使用AddShape函数将部件加入飞行器A
	A.AddShape(new ShapeMainBody1(mb1));
	//A.AddShape(new ShapeMainBody2(mb2));	//将头身组合体部件加入飞行器
	//A.AddShape(new ShapeMainBody2);		//将头身组合体部件加入飞行器
	//A.AddShape(new ShapeMainBody3(mb3));
	
	//A.AddShape(new ShapeWing1(sw1));
	//A.AddShape(new ShapeWing2);			//将机翼部件加入飞行器
	//A.AddShape(new ShapeWing2(sw2));
	
	//A.AddShape(new ShapeTail1);			//将单垂尾部件加入飞行器
	//A.AddShape(new ShapeTail1(st1));
	//A.AddShape(new ShapeTail2);
	//A.AddShape(new ShapeTail2(st2));		//将双垂尾部件加入飞行器
	
	
	A.CalculateAircraft();					//生成飞行器网格
	A.SaveTecPlotAll("C:/Users/yycab/Desktop/database/c++/CST/example/A.dat");		//生成可视化软件可读文本文件
	//A.calcAeroForce();						//调用aero_calc计算气动力
	//A.CalcStruct1D();							//生成梁壳结合结构网格
	//A.CalcStruct2D();							//生成壳单元结构网格
	////A.GetBone();							//生成骨架梁模型
	//
	//A.strct1d.SaveAsAbaqus("test.inp");
	//A.strct1d.SaveAsTecplt("dwadw");
	//A.bone.SaveAsTcp("bone.dat");			//存储
	//A.bone.SaveAsAbaqus("bone_abaqus.inp");

}

void ShapeX34()
{
	//X34
	SMB2 X34mb2;
	X34mb2.NFaiL = X34mb2.NFaiU = 40;  X34mb2.NHeight = 4;
	X34mb2.N1 = { 0.2,0.2,0.1,0.1 };
	X34mb2.N2 = { 0.4,0.4,0.1,0.1 };
	X34mb2.N3 = { 0.1,0.1,0.1,0.1 };
	X34mb2.LHead[0] += 1000;	X34mb2.LBody[0] -= 1000;
	X34mb2.LHead[0] *= 1.2;	X34mb2.LBody[0] *= 1.2;
	X34mb2.LHead[3] = X34mb2.LBody[3] = 3500;
	X34mb2.LBody[1] = X34mb2.LHead[1] = 2300;

	SW1 X34sw1;
	X34sw1.Origin = { 16500 * 1.2,-900,1830 };
	X34sw1.SpanL = 5020 * 1.2;
	X34sw1.RootChordL = 6500 * 1.2;
	X34sw1.TipRootRatio = 0.268;
	X34sw1.Thickness = 600;
	X34sw1.NFaiL = X34sw1.NFaiU = X34sw1.NEta = 40;
	X34sw1.NHeight = 4;
	X34sw1.ifUseMid = false;

	ST1 X34st1;
	X34st1.Origin = { 24000,2330,0 };
	X34st1.NFaiL = X34st1.NFaiU = X34st1.NEta = 40; X34st1.NHeight = 3;
	X34st1.RootChordL = 3500;
	X34st1.SpanL = 3500;
	X34st1.SweepBackAngle = 30;


	Aircraft X34;								//实例化一个部件
	X34.AddShape(new ShapeMainBody2(X34mb2));	//将头身组合体部件加入飞行器
	X34.AddShape(new ShapeWing1(X34sw1));		//将机翼部件加入飞行器
	X34.AddShape(new ShapeTail1(X34st1));		//将单垂尾部件加入飞行器

	X34.CalculateAircraft();					//生成飞行器网格
	X34.SaveTecPlotAll("AircraftData/X34.dat");	//生成可视化软件可读文本文件
	X34.calcAeroForce("C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh");//调用aero_calc计算气动力
	X34.CalcStruct1D();								//生成梁壳结合结构网格

	X34.struct1d.SaveAsAbaqus("AircraftData/X34.inp");
	X34.struct1d.SaveAsTecplt("AircraftData/X34");
	X34.bone.SaveAsTcp("AircraftData/X34-bone.dat");			//存储
	X34.bone.SaveAsAbaqus("AircraftData/X34-bone.inp");
}

void ShapeX37()
{
	string name = "X100";

	SMB3 sm3;
	sm3.N4 = { 0.4,0.4,0.1,0.1 };
	sm3.Ratio3 = 4;
	sm3.NEta = sm3.NFaiL = sm3.NFaiU = 30;
	sm3.GridRefineType = -1;
	//sm3.BUPP1 = {
	//	{1,1,1,2,2,2},
	//	{1,1,1,1,1,3},
	//	{1,2,2,2,2,2},
	//	{1,3,1,1,1,1} };
	//sm3.BUPP2 = {
	//	{1,2,1,2,1,1},
	//	{2,4,2,4,2,2},
	//	{3,6,3,6,3,3},
	//	{1,2,1,2,1,1},
	//	{2,4,2,4,2,2},
	//	{3,6,3,6,3,3} };
	//sm3.BUPP3 = {
	//	{1,2,3,1,2,3},
	//	{1,2,3,1,2,3} };
	SW1 sw1;
	sw1.Origin = { 5000 ,-900,1830 };
	sw1.SpanL = 3020;
	sw1.RootChordL = 4500;
	sw1.TipRootRatio = 0.268;
	sw1.Thickness = 600;
	sw1.NFaiL = sw1.NFaiU = sw1.NEta = 30;
	sw1.NHeight = 4;
	sw1.ifUseMid = false;

	SW2 sw2;
	sw2.Origin[0] = 5000;
	sw2.NFaiL = sw2.NFaiU = 30;
	sw2.NEta_inn = sw2.NEta_out = 30;
	sw2.NHeight = 4;
	sw2.RootChordL_inn = 7400;
	sw2.SpanL_inn = 1500;
	sw2.SweepBackAngle_inn = 72;
	sw2.SweepBackAngle_out = 30;
	sw2.SpanL_out = 3300;
	sw2.TipRootRatio_inn = 1.0 / 2.5;
	sw2.TipRootRatio_out = 0.65*0.7;
	sw2.ifUseMid = false;

	ST2 st2;
	st2.NFaiL = st2.NFaiU = st2.NEta = 30; st2.NHeight = 3;
	st2.Origin = { 16000,1200,1800 };
	st2.SpanL = 5000;
	st2.RootChordL = 3000;
	st2.SideAngle = 45;

	Aircraft A;								//实例化一个部件
	A.AddShape(new ShapeMainBody3(sm3));
	A.AddShape(new ShapeWing2(sw2));
	A.AddShape(new ShapeTail2(st2));		//将双垂尾部件加入飞行器

	A.CalculateAircraft();					//生成飞行器网格
	A.getVol("AircraftData/" + name + "_Vol.txt");
	A.SaveTecPlotAll("AircraftData/"+name+".dat");	//生成可视化软件可读文本文件
	A.calcAeroForce("C:\\Users\\yycab\\source\\repos\\test001\\test001\\mesh");	//调用aero_calc计算气动力
	A.getVol();
	A.GetBone();
	A.CalcStruct1D();								//生成梁壳结合结构网格
	A.CalcStruct2D();							//生成壳单元结构网格
	//A.strct1d.SaveAsAbaqus("AircraftData/" + name + ".inp");
	A.struct1d.SaveAsTecplt("AircraftData/" + name + "");
	//A.strct1d.SaveAsNastran("C:/Users/yycab/Desktop/bdf/X37/" + name);
	A.bone.SaveAsTcp("AircraftData/" + name + "-bone.dat");			//存储
	//A.bone.SaveAsAbaqus("AircraftData/" + name + "-bone.inp");
	
}

void testInterp()
{
	mat X = ones(1, 4);
	for (uword i = 0; i < X.n_cols; i++)
	{
		X.col(i) = X.col(i) * i;
	}
	mat Y = X;
	mat Z = { {0,2,44,6} };
	cout << X << endl;
	cout << Z << endl;
	mat X1 = { {0.2,0.3} }, Y1;
	interp1(X, trans(Z), trans(X1), Y1);
	cout << Y1 << endl;

	mat Z2d = { {1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16} };
	mat Z1;
	interp2(X, Y, Z2d, X1, Y1, Z1);
	cout << Z1 << endl;
}

void testXfoil(string path)
{
	XFoilInfo x;
	x.workPath = path;
	x.foilName = "testfoil";
	x.foilList = {
		Field<double>(0.999999999999998, -0.0157947278957458),
		Field<double>(0.994992310390415, -0.0145913714486130),
		Field<double>(0.989979224563579, -0.0133988074359859),
		Field<double>(0.984971534953996, -0.0122062434233611),
		Field<double>(0.979958449127160, -0.0110244718452371),
		Field<double>(0.974945363300324, -0.00984809648436348),
		Field<double>(0.969932277473488, -0.00867711734074263),
		Field<double>(0.964913795429405, -0.00751153441437458),
		Field<double>(0.959900709602569, -0.00635674392250969),
		Field<double>(0.954882227558480, -0.00520195343064480),
		Field<double>(0.949858349297143, -0.00405795537328307),
		Field<double>(0.944839867253055, -0.00291395731592134),
		Field<double>(0.939815988991713, -0.00178075169306277),
		Field<double>(0.934792110730377, -0.000658338504707355),
		Field<double>(0.929762836251787, 0.000464074683648059),
		Field<double>(0.924738957990446, 0.00158109165475067),
		Field<double>(0.919709683511857, 0.00268731619135012),
		Field<double>(0.914680409033267, 0.00379354072794957),
		Field<double>(0.909645738337425, 0.00488897283004586),
		Field<double>(0.904616463858836, 0.00597900871489180),
		Field<double>(0.899581793162994, 0.00706364838248494),
		Field<double>(0.894547122467151, 0.00813749561557492),
		Field<double>(0.889512451771309, 0.00921134284866244),
		Field<double>(0.884472384858214, 0.0102743976472493),
		Field<double>(0.879432317945119, 0.0113374524458361),
		Field<double>(0.874392251032024, 0.0123897148099197),
		Field<double>(0.869352184118929, 0.0134365809567506),
		Field<double>(0.864306720988586, 0.0144780508863286),
		Field<double>(0.859266654075491, 0.0155141245986563),
		Field<double>(0.854221190945144, 0.0165394058764809),
		Field<double>(0.849170331597548, 0.0175646871543054),
		Field<double>(0.844124868467205, 0.0185791759976268),
		Field<double>(0.839074009119604, 0.0195828724064426),
		Field<double>(0.834023149772009, 0.0205811725980080),
		Field<double>(0.828972290424413, 0.0215740765723231),
		Field<double>(0.823916034859560, 0.0225561881121325),
		Field<double>(0.818859779294711, 0.0235329034346916),
		Field<double>(0.813803523729863, 0.0245042225400004),
		Field<double>(0.808741871947761, 0.0254593529935507),
		Field<double>(0.803685616382913, 0.0264090872298531),
		Field<double>(0.798618568383564, 0.0273534252489003),
		Field<double>(0.793556916601462, 0.0282815746161940),
		Field<double>(0.788489868602108, 0.0292043277662373),
		Field<double>(0.783422820602754, 0.0301162884817750),
		Field<double>(0.778355772603405, 0.0310174567628095),
		Field<double>(0.773283328386798, 0.0319078326093409),
		Field<double>(0.768210884170196, 0.0327928122386219),
		Field<double>(0.763133043736341, 0.0336616032161446),
		Field<double>(0.758060599519734, 0.0345196017591665),
		Field<double>(0.752977362868627, 0.0353668078676852),
		Field<double>(0.747899522434772, 0.0362032215416984),
		Field<double>(0.742816285783664, 0.0370288427812108),
		Field<double>(0.737733049132557, 0.0378382753689673),
		Field<double>(0.732644416264196, 0.0386423117394710),
		Field<double>(0.727555783395836, 0.0394301594582211),
		Field<double>(0.722467150527480, 0.0402018185252153),
		Field<double>(0.717373121441867, 0.0409626851577039),
		Field<double>(0.712279092356259, 0.0417127593556918),
		Field<double>(0.707185063270646, 0.0424520411191765),
		Field<double>(0.702085637967784, 0.0431751342309052),
		Field<double>(0.696986212664918, 0.0438820386908780),
		Field<double>(0.691886787362057, 0.0445781507163477),
		Field<double>(0.686781965841943, 0.0452634703073141),
		Field<double>(0.681677144321830, 0.0459326012465271),
		Field<double>(0.676572322801716, 0.0465855435339841),
		Field<double>(0.671462105064349, 0.0472276933869355),
		Field<double>(0.666351887326982, 0.0478536545881334),
		Field<double>(0.661241669589615, 0.0484634271375778),
		Field<double>(0.656126055635001, 0.0490624072525165),
		Field<double>(0.651010441680381, 0.0496451987157018),
		Field<double>(0.645894827725762, 0.0502171977443814),
		Field<double>(0.640773817553894, 0.0507784043385579),
		Field<double>(0.635652807382022, 0.0513288184982337),
		Field<double>(0.630537193427407, 0.0518630440061510),
		Field<double>(0.625410787038287, 0.0523918732968180),
		Field<double>(0.620289776866415, 0.0529045139357315),
		Field<double>(0.615168766694547, 0.0534117583573922),
		Field<double>(0.610042360305422, 0.0539028141272969),
		Field<double>(0.604915953916302, 0.0543830774626985),
		Field<double>(0.599789547527182, 0.0548525483635969),
		Field<double>(0.594657744920809, 0.0553166230472450),
		Field<double>(0.589531338531688, 0.0557645090791346),
		Field<double>(0.584399535925315, 0.0562069988937764),
		Field<double>(0.579273129536195, 0.0566333000566597),
		Field<double>(0.574141326929822, 0.0570542050022927),
		Field<double>(0.569009524323449, 0.0574643175134225),
		Field<double>(0.563872325499823, 0.0578690338073020),
		Field<double>(0.558740522893450, 0.0582629576666758),
		Field<double>(0.553608720287077, 0.0586460890915465),
		Field<double>(0.548471521463451, 0.0590184280819140),
		Field<double>(0.543334322639830, 0.0593853708550312),
		Field<double>(0.538202520033457, 0.0597415211936428),
		Field<double>(0.533065321209831, 0.0600922753150040),
		Field<double>(0.527928122386205, 0.0604322370018621),
		Field<double>(0.522790923562579, 0.0607614062542170),
		Field<double>(0.517648328521706, 0.0610851792893191),
		Field<double>(0.512511129698080, 0.0614035561071684),
		Field<double>(0.507373930874454, 0.0617111404905170),
		Field<double>(0.502231335833580, 0.0620133286566127),
		Field<double>(0.497094137009954, 0.0623047243882054),
		Field<double>(0.491951541969076, 0.0625853276852924),
		Field<double>(0.486808946928202, 0.0628605347651290),
		Field<double>(0.481671748104576, 0.0631303456277153),
		Field<double>(0.476529153063702, 0.0633893640557960),
		Field<double>(0.471386558022824, 0.0636429862666264),
		Field<double>(0.466243962981950, 0.0638858160429536),
		Field<double>(0.461101367941071, 0.0641178533847776),
		Field<double>(0.455958772900198, 0.0643444945093488),
		Field<double>(0.450816177859319, 0.0645657394166697),
		Field<double>(0.445668186601192, 0.0647761918894850),
		Field<double>(0.440525591560314, 0.0649812481450499),
		Field<double>(0.435382996519440, 0.0651755119661116),
		Field<double>(0.430235005261308, 0.0653589833526702),
		Field<double>(0.425092410220435, 0.0655370585219760),
		Field<double>(0.419944418962308, 0.0657097374740314),
		Field<double>(0.414801823921429, 0.0658716239915813),
		Field<double>(0.409653832663303, 0.0660227180746304),
		Field<double>(0.404505841405171, 0.0661684159404267),
		Field<double>(0.399363246364298, 0.0663087175889702),
		Field<double>(0.394215255106171, 0.0664382268030105),
		Field<double>(0.389067263848040, 0.0665569435825477),
		Field<double>(0.383919272589913, 0.0666702641448346),
		Field<double>(0.378771281331786, 0.0667727922726182),
		Field<double>(0.373628686290908, 0.0668645279658963),
		Field<double>(0.368480695032781, 0.0669508674419240),
		Field<double>(0.363332703774655, 0.0670318107006990),
		Field<double>(0.358184712516523, 0.0671019615249707),
		Field<double>(0.353036721258397, 0.0671613199147393),
		Field<double>(0.347888730000270, 0.0672098858700047),
		Field<double>(0.342740738742139, 0.0672530556080198),
		Field<double>(0.337592747484012, 0.0672854329115293),
		Field<double>(0.332444756225886, 0.0673124139977884),
		Field<double>(0.327296764967754, 0.0673286026495419),
		Field<double>(0.322148773709627, 0.0673339988667947),
		Field<double>(0.317000782451501, 0.0673286026495419),
		Field<double>(0.311852791193369, 0.0673178102150388),
		Field<double>(0.306704799935243, 0.0672962253460325),
		Field<double>(0.301556808677116, 0.0672638480425230),
		Field<double>(0.296408817418990, 0.0672260745217607),
		Field<double>(0.291260826160858, 0.0671721123492425),
		Field<double>(0.286112834902732, 0.0671073577422235),
		Field<double>(0.280964843644605, 0.0670372069179518),
		Field<double>(0.275816852386474, 0.0669508674419240),
		Field<double>(0.270674257345600, 0.0668537355313932),
		Field<double>(0.265526266087473, 0.0667512074036095),
		Field<double>(0.260378274829342, 0.0666324906240723),
		Field<double>(0.255230283571215, 0.0665083776272848),
		Field<double>(0.250087688530337, 0.0663680759787388),
		Field<double>(0.244939697272210, 0.0662169818956921),
		Field<double>(0.239797102231336, 0.0660550953781399),
		Field<double>(0.234649110973205, 0.0658770202088341),
		Field<double>(0.229506515932331, 0.0656881526050251),
		Field<double>(0.224363920891453, 0.0654884925667106),
		Field<double>(0.219221325850579, 0.0652726438766425),
		Field<double>(0.214078730809700, 0.0650406065348185),
		Field<double>(0.208936135768826, 0.0647977767584937),
		Field<double>(0.203793540727948, 0.0645387583304106),
		Field<double>(0.198650945687074, 0.0642635512505739),
		Field<double>(0.193513746863448, 0.0639721555189813),
		Field<double>(0.188376548039822, 0.0636645711356352),
		Field<double>(0.183239349216196, 0.0633407981005306),
		Field<double>(0.178102150392575, 0.0629954401964222),
		Field<double>(0.172964951568950, 0.0626392898578106),
		Field<double>(0.167833148962576, 0.0622615546501903),
		Field<double>(0.162701346356203, 0.0618676307908164),
		Field<double>(0.157569543749830, 0.0614521220624338),
		Field<double>(0.152437741143457, 0.0610150284650473),
		Field<double>(0.147311334754337, 0.0605563499986521),
		Field<double>(0.142184928365212, 0.0600814828805009),
		Field<double>(0.137063918193345, 0.0595796346760930),
		Field<double>(0.131942908021477, 0.0590562016026763),
		Field<double>(0.126821897849605, 0.0585057874430030),
		Field<double>(0.121706283894990, 0.0579337884143209),
		Field<double>(0.116596066157623, 0.0573294120821318),
		Field<double>(0.111485848420257, 0.0566926584464283),
		Field<double>(0.106381026900143, 0.0560235275072178),
		Field<double>(0.101281601597277, 0.0553220192644953),
		Field<double>(0.0961875725116684, 0.0545773412837603),
		Field<double>(0.0911043358605607, 0.0537948897822629),
		Field<double>(0.0860264954267059, 0.0529530798909970),
		Field<double>(0.0809540512100990, 0.0520573078272128),
		Field<double>(0.0758977956452505, 0.0510913849391568),
		Field<double>(0.0708577287321555, 0.0500553112268291),
		Field<double>(0.0658338504708190, 0.0489275018212234),
		Field<double>(0.0608315570784887, 0.0477079567223371),
		Field<double>(0.0558562447724127, 0.0463804872784192),
		Field<double>(0.0509187059871014, 0.0449289048377089),
		Field<double>(0.0460189407225499, 0.0433424169657079),
		Field<double>(0.0411785338477695, 0.0415994387934049),
		Field<double>(0.0363974853627601, 0.0396783854517986),
		Field<double>(0.0317135687882767, 0.0375468796373747),
		Field<double>(0.0271483689933354, 0.0351671478293731),
		Field<double>(0.0227450557159439, 0.0324906240725262),
		Field<double>(0.0185575911286156, 0.0295011197150803),
		Field<double>(0.0146021638831186, 0.0262040309742882),
		Field<double>(0.0109057550656975, 0.0226209427191539),
		Field<double>(0.00751153441437458, 0.0187518549496798),
		Field<double>(0.00452203005692868, 0.0145643903623564),
		Field<double>(0.00213150581442154, 0.0100045867846654),
		Field<double>(0.000561206594174027, 0.00511021773736674),
		Field<double>(0, 0),
		Field<double>(0.000760866632488583, -0.00502387826133902),
		Field<double>(0.00284920270889749, -0.00965922888055453),
		Field<double>(0.00595742384588058, -0.0136955993848313),
		Field<double>(0.00974017213932699, -0.0171006124706574),
		Field<double>(0.0139114480748969, -0.0200361546555850),
		Field<double>(0.0183039689177877, -0.0226263389364043),
		Field<double>(0.0228583762782283, -0.0249089388338774),
		Field<double>(0.0275368966354589, -0.0269379165205177),
		Field<double>(0.0322909640342092, -0.0287780266033515),
		Field<double>(0.0371043898227305, -0.0304670426031338),
		Field<double>(0.0419609853492693, -0.0320265493888781),
		Field<double>(0.0468445619620624, -0.0334781318295860),
		Field<double>(0.0517605158783673, -0.0348433747942686),
		Field<double>(0.0566926584464258, -0.0361330707174266),
		Field<double>(0.0616409896662428, -0.0373472195990601),
		Field<double>(0.0666109057550660, -0.0385074063081778),
		Field<double>(0.0715862180611371, -0.0396190270620301),
		Field<double>(0.0765723228017138, -0.0406820818606144),
		Field<double>(0.0815692199767962, -0.0417073631384390),
		Field<double>(0.0865715133691264, -0.0427002671127540),
		Field<double>(0.0915792029787095, -0.0436553975663068),
		Field<double>(0.0965922888055454, -0.0445835469336005),
		Field<double>(0.101610770849634, -0.0454793189973822),
		Field<double>(0.106634649110971, -0.0463481099749073),
		Field<double>(0.111663923589560, -0.0471953160834260),
		Field<double>(0.116698594285402, -0.0480155411056857),
		Field<double>(0.121733264981244, -0.0488141812589365),
		Field<double>(0.126773331894339, -0.0495912365431835),
		Field<double>(0.131813398807434, -0.0503521031756746),
		Field<double>(0.136858861937782, -0.0510859887219065),
		Field<double>(0.141909721285378, -0.0518036856163824),
		Field<double>(0.146960580632973, -0.0524997976418521),
		Field<double>(0.152016836197822, -0.0531797210155682),
		Field<double>(0.157073091762670, -0.0538380595202755),
		Field<double>(0.162129347327524, -0.0544856055904797),
		Field<double>(0.167190999109620, -0.0551061705744273),
		Field<double>(0.172252650891721, -0.0557159431238692),
		Field<double>(0.177319698891076, -0.0563095270215576),
		Field<double>(0.182386746890430, -0.0568815260502397),
		Field<double>(0.187453794889779, -0.0574427326444162),
		Field<double>(0.192520842889133, -0.0579877505868391),
		Field<double>(0.197593287105735, -0.0585219760947565),
		Field<double>(0.202665731322342, -0.0590454091681732),
		Field<double>(0.207738175538944, -0.0595526535898338),
		Field<double>(0.212816015972804, -0.0600491055769914),
		Field<double>(0.217888460189406, -0.0605293689123929),
		Field<double>(0.222966300623260, -0.0609988398132914),
		Field<double>(0.228044141057115, -0.0614575182796866),
		Field<double>(0.233127377708223, -0.0619000080943259),
		Field<double>(0.238205218142083, -0.0623317054744620),
		Field<double>(0.243288454793190, -0.0627418179855919),
		Field<double>(0.248371691444298, -0.0631411380622185),
		Field<double>(0.253460324312653, -0.0635242694870892),
		Field<double>(0.258543560963761, -0.0638966084774567),
		Field<double>(0.263632193832121, -0.0642473625988180),
		Field<double>(0.268720826700482, -0.0645819280684232),
		Field<double>(0.273809459568842, -0.0649003048862750),
		Field<double>(0.278898092437198, -0.0652024930523708),
		Field<double>(0.283992121522811, -0.0654884925667106),
		Field<double>(0.289080754391171, -0.0657583034292969),
		Field<double>(0.294174783476780, -0.0660119256401248),
		Field<double>(0.299268812562393, -0.0662493591992016),
		Field<double>(0.304362841648001, -0.0664706041065200),
		Field<double>(0.309456870733614, -0.0666810565793377),
		Field<double>(0.314550899819223, -0.0668753204003995),
		Field<double>(0.319650325122089, -0.0670587917869581),
		Field<double>(0.324744354207697, -0.0672260745217607),
		Field<double>(0.329843779510563, -0.0673771686048074),
		Field<double>(0.334937808596171, -0.0675174702533533),
		Field<double>(0.340037233899037, -0.0676415832501409),
		Field<double>(0.345136659201899, -0.0677495075951749),
		Field<double>(0.350236084504760, -0.0678358470712026),
		Field<double>(0.355330113590373, -0.0679113941127272),
		Field<double>(0.360429538893234, -0.0679707525024958),
		Field<double>(0.365528964196095, -0.0680085260232581),
		Field<double>(0.370628389498961, -0.0680355071095147),
		Field<double>(0.375727814801822, -0.0680462995440203),
		Field<double>(0.380827240104683, -0.0680409033267675),
		Field<double>(0.385926665407550, -0.0680247146750116),
		Field<double>(0.391026090710411, -0.0679869411542517),
		Field<double>(0.396125516013272, -0.0679383751989863),
		Field<double>(0.401224941316138, -0.0678790168092153),
		Field<double>(0.406324366618999, -0.0678034697676931),
		Field<double>(0.411418395704607, -0.0677117340744126),
		Field<double>(0.416517821007473, -0.0676092059466314),
		Field<double>(0.421617246310334, -0.0674904891670942),
		Field<double>(0.426716671613196, -0.0673663761703042),
		Field<double>(0.431810700698809, -0.0672206783045079),
		Field<double>(0.436910126001670, -0.0670695842214612),
		Field<double>(0.442004155087283, -0.0668969052694058),
		Field<double>(0.447098184172891, -0.0667188301001000),
		Field<double>(0.452197609475757, -0.0665245662790383),
		Field<double>(0.457291638561366, -0.0663141138062205),
		Field<double>(0.462385667646979, -0.0660928688989021),
		Field<double>(0.467479696732587, -0.0658608315570781),
		Field<double>(0.472573725818200, -0.0656126055635006),
		Field<double>(0.477667754903809, -0.0653535871354174),
		Field<double>(0.482756387772169, -0.0650837762728335),
		Field<double>(0.487850416857782, -0.0647923805412409),
		Field<double>(0.492939049726138, -0.0644955885923955),
		Field<double>(0.498027682594498, -0.0641772117745462),
		Field<double>(0.503116315462859, -0.0638426463049410),
		Field<double>(0.508204948331219, -0.0634918921835797),
		Field<double>(0.513288184982327, -0.0631249494104625),
		Field<double>(0.518376817850687, -0.0627418179855919),
		Field<double>(0.523460054501795, -0.0623371016917148),
		Field<double>(0.528537894935650, -0.0619161967460819),
		Field<double>(0.533621131586757, -0.0614737069314401),
		Field<double>(0.538698972020612, -0.0610150284650449),
		Field<double>(0.543776812454467, -0.0605401613468961),
		Field<double>(0.548849256671074, -0.0600437093597386),
		Field<double>(0.553921700887676, -0.0595310687208275),
		Field<double>(0.558994145104283, -0.0590022394301581),
		Field<double>(0.564066589320885, -0.0584572214877376),
		Field<double>(0.569133637320239, -0.0579014111108114),
		Field<double>(0.574200685319588, -0.0573294120821293),
		Field<double>(0.579267733318942, -0.0567412244016937),
		Field<double>(0.584329385101043, -0.0561422442867550),
		Field<double>(0.589391036883145, -0.0555324717373131),
		Field<double>(0.594452688665246, -0.0549065105361127),
		Field<double>(0.599514340447342, -0.0542697569004117),
		Field<double>(0.604570596012196, -0.0536222108302075),
		Field<double>(0.609626851577044, -0.0529638723255001),
		Field<double>(0.614683107141893, -0.0522947413862872),
		Field<double>(0.619733966489488, -0.0516148180125735),
		Field<double>(0.624790222054337, -0.0509241022043542),
		Field<double>(0.629841081401938, -0.0502225939616318),
		Field<double>(0.634886544532280, -0.0495102932844062),
		Field<double>(0.639937403879881, -0.0487925963899302),
		Field<double>(0.644982867010224, -0.0480641070609511),
		Field<double>(0.650028330140572, -0.0473302215147167),
		Field<double>(0.655073793270915, -0.0465855435339817),
		Field<double>(0.660119256401262, -0.0458354693359962),
		Field<double>(0.665159323314357, -0.0450799989207556),
		Field<double>(0.670204786444700, -0.0443137360710142),
		Field<double>(0.675244853357795, -0.0435420770040200),
		Field<double>(0.680279524053637, -0.0427596255025226),
		Field<double>(0.685319590966732, -0.0419771740010253),
		Field<double>(0.690359657879823, -0.0411893262822751),
		Field<double>(0.695394328575665, -0.0403960823462746),
		Field<double>(0.700428999271507, -0.0395974421930213),
		Field<double>(0.705469066184602, -0.0387934058225176),
		Field<double>(0.710503736880444, -0.0379839732347611),
		Field<double>(0.715538407576286, -0.0371745406470071),
		Field<double>(0.720573078272128, -0.0363597118419978),
		Field<double>(0.725602352750723, -0.0355502792542414),
		Field<double>(0.730642419663813, -0.0347462428837377),
		Field<double>(0.735677090359655, -0.0339422065132340),
		Field<double>(0.740711761055497, -0.0331489625772335),
		Field<double>(0.745751827968592, -0.0323611148584834),
		Field<double>(0.750791894881687, -0.0315894557914892),
		Field<double>(0.755831961794782, -0.0308339853762509),
		Field<double>(0.760877424925125, -0.0300893073955159),
		Field<double>(0.765928284272725, -0.0293716105010375),
		Field<double>(0.770979143620321, -0.0286701022583175),
		Field<double>(0.776030002967917, -0.0279901788846013),
		Field<double>(0.781086258532765, -0.0273318403798940),
		Field<double>(0.786147910314867, -0.0267004829614433),
		Field<double>(0.791209562096968, -0.0260853141947486),
		Field<double>(0.796276610096322, -0.0254971265143130),
		Field<double>(0.801343658095671, -0.0249359199201365),
		Field<double>(0.806410706095025, -0.0243909019777135),
		Field<double>(0.811483150311632, -0.0238674689042969),
		Field<double>(0.816560990745487, -0.0233710169171418),
		Field<double>(0.821638831179342, -0.0228961497989906),
		Field<double>(0.826716671613197, -0.0224428675498481),
		Field<double>(0.831794512047051, -0.0220111701697095),
		Field<double>(0.836877748698159, -0.0216010576585822),
		Field<double>(0.841966381566519, -0.0212179262337090),
		Field<double>(0.847049618217627, -0.0208509834605943),
		Field<double>(0.852138251085987, -0.0205110217737362),
		Field<double>(0.857226883954348, -0.0201926449558845),
		Field<double>(0.862320913039956, -0.0198958530070415),
		Field<double>(0.867409545908317, -0.0196260421444577),
		Field<double>(0.872503574993930, -0.0193724199336273),
		Field<double>(0.877597604079538, -0.0191457788090536),
		Field<double>(0.882691633165151, -0.0189407225534887),
		Field<double>(0.887785662250760, -0.0187572511669326),
		Field<double>(0.892885087553626, -0.0186007608666331),
		Field<double>(0.897984512856487, -0.0184604592180872),
		Field<double>(0.903078541942095, -0.0183471386558028),
		Field<double>(0.908177967244961, -0.0182554029625223),
		Field<double>(0.913277392547822, -0.0181852521382505),
		Field<double>(0.918376817850683, -0.0181366861829851),
		Field<double>(0.923476243153550, -0.0181097050967259),
		Field<double>(0.928575668456411, -0.0181043088794756),
		Field<double>(0.933675093759272, -0.0181258937484819),
		Field<double>(0.938774519062138, -0.0181636672692442),
		Field<double>(0.943873944364999, -0.0182284218762656),
		Field<double>(0.948967973450607, -0.0183147613522908),
		Field<double>(0.954067398753473, -0.0184226856973249),
		Field<double>(0.959166824056334, -0.0185521949113677),
		Field<double>(0.964260853141948, -0.0187032889944143),
		Field<double>(0.969360278444809, -0.0188759679464698),
		Field<double>(0.974454307530422, -0.0190756279847819),
		Field<double>(0.979548336616030, -0.0192914766748500),
		Field<double>(0.984642365701643, -0.0195343064511771),
		Field<double>(0.989736394787252, -0.0197987210965106),
		Field<double>(0.994825027655612, -0.0200847206108505),
		Field<double>(0.999913660523973, -0.0203977012114494),
	};
	x.VISC = 5000000;
	x.Ma = 0.76;
	x.ITER = 150;
	x.alpha = {0.4, 0.4, 1};
	XFoilSolver xsolver(x);
	xsolver.ExePath = path + "\\xfoil.exe";
	xsolver.solve();
}

void test_curve()
{
	CST_info cst_info;
	cst_info.N1 = 0.5;
	cst_info.N2 = 0.5;
	Curve* cv = new Curve_cst(cst_info);
	vector<double> list = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
	vector<double> res = cv->get(list,false);
	vector<double> res1 = cv->get(list,true);

	auto funcc = [](double x)
	{ return x * x*4-100; };
	Curve *ccv = new Curve_func(funcc);
	auto res2 = ccv->get(list, true);

	GuideFunc gf;
	gf.pt1 = shared_ptr<Curve>(new Curve_func([](double x){ return 0; }));
	gf.pt2 = shared_ptr<Curve>(new Curve_func([](double x){ return 0; }));
	gf.pt3 = shared_ptr<Curve>(new Curve_func([](double x){ return 0; }));
	gf.theta1 = shared_ptr<Curve>(new Curve_func([](double x){ return 90 * x; }));
	gf.theta2 = shared_ptr<Curve>(new Curve_func([](double x){ return 90 * x; }));
	gf.theta3 = shared_ptr<Curve>(new Curve_func([](double x){ return 0; }));

	Point res_p = gf.update_point(Point(0, 1, 0), 0.5);
}

int main()
{
	testAircraft();

	return 0;
}

int main1(int argc, char* argv[])
{
	string exefilepath = argv[0];
	string exePath = exefilepath.substr(0, exefilepath.find_last_of("/"));//获取exe所在文件夹路径
	if (exePath == exefilepath)//相对路径
	{
		exePath = exefilepath.substr(0, exefilepath.find_last_of("\\"));//获取文件夹路径
		if (exePath == exefilepath)//相对路径
		{
			exePath = "";
		}
	//docPath = "";
	}
	//计算输入文件路径
	string datpath = "";
	switch (argc)
	{
	case 2:
		datpath = argv[1];//这里的路径必须是反斜杠才可以正确生成输出文件
		break;
	case 1:
		datpath = "E:\\Matlab-Xfoil\\cst\\inpSW.cst";//
		break;
	default:
		cout << "错误的输入参数个数" << endl;
		return -1;
	}
	

	cout << "$----------------------------------------$\n";
	cout << "$       3D-CST integration program       $\n";
	cout << "$                 v 5.0.0                $\n";
	cout << "$----------------------------------------$\n";
	cout << std::flush;

	Aircraft A;
	
	A.setEXEworkPath(exePath);
	int state = A.RunFromFile(datpath);
	//system("pause");
	return state;
}