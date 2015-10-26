#pragma once

#include "StructData.h"
#include "MEM.h"



class CEL
{
//типы элементов по геометрии
	//22 - 6 узловой плоский треугольник
	//24 - 10 узловой тетраэдр
	//25 - 20 узловой кубический

public:
	CEL(void);
public:
	~CEL(void);

public:
	MEM *MM;
	int ielobjnumber;
public:
	int eltype; //тип элемента
	int NORT; //число пространственных координат в узле
	int NORTfullcrd; //число пространственных координат в узле в полном массиве координат модели (нужно для обработки плоских граней трехмерных моделей)
	double lcrdlim[2]; //показывает пределы изменения локальных координат внутри элемента
	int NDOF; //число степеней свободы в узле
	int NN; //число узлов в элементе
	int NNE; //число степеней свободы в элементе 	NNE = NN*NDOF;
	int NNORT; //число координат узлов в элементе 	NNORT = NN*NORT;
	int NDEF; // число компонентов деформаций
	int NSTRS; // число компонентов напряжений
	
	int *ind; //массив номеров узлов элемента (указатель на соответствующую строку матрицы индексов)
	double *fullcrd; //указатель на полный массив координат узлов модели

	double **STIF; // матрица жесткости элемента
	double **MASS; //матрица масс

	double **B; // матрица градиентов
	double **D; //матрица свойств для текущего элемента

	//локальные координаты точек интегрирования по объему
	double **vlcrdint; //[номер точки][номер координаты]
	int nvlcrd; //количество локальных координат при интегрировании по объему (для тетраэдра, например, 4)
				//однако производные берутся только по независимым координатам NORT (для тетраэдра 1 производная зависима)
	int NPINT; //количество точек интегрирования по объему
	double *Wvint;//весовые коэффициенты в точках интегрирования по объему

	//локальные координаты точки интегрирования по граням
	double ***faceloccrdint; //[грань][номер точки][номер координаты]
	int Nface; //количество граней
	int NNface; //количество узлов на грани
	int Nfp; //количество точек интегрирования на грани
	int nflcrd; //количество локальных координат при интегрировании по грани
	double **Wfint;//весовые коэффициенты в точках интегрирования [грань][номер точки]
	int *faceorientcrd_num; //номер лок. координаты, которая наа грани является постоянной
	double *faceorientcrd_val; //значение локальной координаты, которая постоянна на грани

	//массивы флагов принадлежности узлов в локальной нумерации граням
	int **FaceNodes; //[грань][номер узла в лок. нумерации элемента] = 1, если данный узел присутствует в грани
	int **FaceNodesList; //[грань][номер узла в лок. нумерации грани] = локальному номеру узла в нумерации элемента 

	//функции формы и интегрирование
	double **FFv;//матрица значений функций формы в точках интегрир-ия по объему
	double ***FFf;//матрица значений функций формы в точках интегрир-ия по граням
	double ***dFFvloc;	//матрица значений производных функций формы по локальным координатам в точках интегрир-ия по объему;
					//[точка инт][номер координаты взятия производной][номер функции формы]
	double ****dFFfloc;	//матрица значений производных функций формы по локальным координатам в точках интегрир-ия по грани;
					//[грань][точка инт][номер координаты взятия производной][номер функции формы]
	
	

	int *GLOBNE; //массив глобальных номеров степеней свободы
	
	double VEL;

	INTPOINT *P; //массив параметров точек интегрирования
	MATPROP *material; //параметры материала

	//флаги

	int isVPbelowzero; //флаг показывает появление отрицательного объема точки интегрирования
	int isVELcalculated; //показывает, определялся ли объем элемента
	int isInitialized;

public:

	void Initialize(int ieltype);
	void GeneralMemInit();
	void InitGLOBNE();

	//задание основных параметров по типам элементов
	void GeneralParameters5();
	void GeneralParameters22();
	void GeneralParameters24();
	void GeneralParameters25();

	//задание локальных координат точек интегрирования по объему и весовых функций
	void Init_VolIntPointLocCrd_5();
	void Init_VolIntPointLocCrd_22();
	void Init_VolIntPointLocCrd_24();
	void Init_VolIntPointLocCrd_25();

	//задание локальных координат точек интегрирования по граням и весовых функций
	void Init_FacePointLocCrd_25();

	//поиск грани, определение площади
	void Init_FaceNodes_24(); //задания массивов флагов принадлежностей узлов граням
	void Init_FaceNodes_25(); //задания массивов флагов принадлежностей узлов граням
	int IdentifyFace(int *facenodes, int nfn);//определение соответствия грани (передаются узлы в глоб нумер)
	double FaceArea(int iface);
	void FaceNodeNums(int iface, int *NodeNums, int *NNf);

	// вычисление значений функций формы
	void FF_all(double *F, double *lcrd);
	void FF_5(double *F, double *lcrd);
	void FF_22(double *F, double *lcrd);
	void FF_24(double *F, double *lcrd);
	void FF_25(double *F, double *lcrd);

	//вычисление производных функций формы по локальным координатам
	void dFF_all(double **dF, double *lcrd);
	void dFF_5(double **dF, double *lcrd);
	void dFF_22(double **dF, double *lcrd);
	void dFF_24(double **dF, double *lcrd);
	void dFF_25(double **dF, double *lcrd);

	//определение матрицы Якоби
	void Jacoby(double **dFloc, double **J);
	void DetJacoby(double **J, double *detJ);
	void InvJacoby(double **J, double detJ);
	int GetLocCrd(double *globcrd, double *loccrd, double *FF); //определение локальных координат точки по глобальным координатам

	//расчет матрицы жесткости в точке интегрирования с накоплением в общей матрице жесткости
	void CalcBGrad_3D(double **invJ, double **dFF);
	void CalcSTIFPoint(INTPOINT *P); //вычисление матрицы жесткости для точки интегрирования

	int DINT(); //вычисление матрицы упругих свойств

	void ELSTIF(); //вычисление матрицы жесткости
	void ELMASS(); //вычисление матрицы масс
	
	void shifting( double shift); //сдвиг по частоте для типа расчета 2

	//геометрические характеристики
	void GetBoundBox(double *vmin, double *vmax);
	void GetCentr(double *vcentr);
};
