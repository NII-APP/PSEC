#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "stdlib.h"
#include "EL.h"
#include "FACE.h"
#include "L_FDPoint.h"

class FULLMODEL
{
public:
	FULLMODEL(void);
public:
	~FULLMODEL(void);

public:
	MEM *MM;

	int MAXELTYPES;
	CEL *el; //библиотека конечных элементов

	int loadtype;

	char name[256];
//	char path[256];
	char netpath[256];
	char pathmatr[256]; // путь к файлам с исходными данными по СЭ
	char pathmain[256]; // путь к файлам с исходными данными по модели и СЭ разбиению в целом

	int nmat; //количество материалов
	MATPROP *mat;

	INTPOINT *P;
	int *envP; //массив показывающий стартовый номер позиции в массиве точек интегрирования, соответствующий каждому конечному элементу
	int nIntP;

	int KORT;
	int NEL; //количество элементов
	int NN; //количество узлов
	int NNE; //количество степеней свободы
	int **IND; // матрица индексов
	double *CRD; // массив координат
	int *MTR; // массив материалов
	int *FIX; // массив закреплений
	double *UFIX; // массив заданных перемещений
	double *Mdiag; //диагональная матрица масс

	//графы
	int **NDEL;
	int **NDGR;
	int **ELGR;

	//массив флагов и список КЭ для вспомогательных целей
	int *LIST_EL;
	int nlistel;
	int *MASK_EL;


	//поверхность модели
	FACE *surfmodel;
	int nfaces; //количество поверхностных граней
	
	// векторы нагрузки и решения
	int nvect;
	double **LV;
	double **RV;
	int fl_LVRV_internal; // == 1 : вектора нагрузок и решений инициализированы внутри объекта FullModel (требуется удаление памяти)

	// набор точек приложения сил и поиска перемещений
	FDPOINT *pforce; //точки приложения сил
	int npforce;
	FDPOINT *pdisp; //точки расчета перемещений
	int npdisp;


	void ReadFullmodel();
	void ReadFromUzor();
	void ReadFromAnsys_FullGrid();
	void ReadFromUzor_NewFormat();

	void InitIntPoint();
	void AttachElement(CEL *el, int elnum);
	void AttachElement(CEL *el, FACE *face); //только для геометрических расчетов
	void MassMatrix();
	void MassAssemblingDiag(CEL *el,int elnum);
	void MassDiagMatrixFix();

	void GenSurf();
	int CheckNewFace( FACE* surfmodel_tmp, int inewface);

	void InitLIST_EL();
	void ClearLIST_EL();
	void DelLIST_EL();

	void CalcActualForcePointPosition();
	void CalcActualDispPointPosition();

	void ParaView_PrintCRD(FILE *fp);
	void ParaView_PrintIND(FILE *fp);
	void ParaView_PrintGrid(FILE *fp);
	void ParaView_PrintMaterial(FILE *fp);
	void ParaView_StartCellDataSection(FILE *fp);
	void ParaView_StartNodeDataSection(FILE *fp);
	void ParaView_SingleDispPrint(FILE *fp, char *name);
	void ParaView_SingleXYZS(FILE *fp, char *name, double *data);
	void ParaView_SingleVector(FILE *fp, char *name, double *data);
	void ParaView_SingleVector(FILE *fp, char *name, float *data);
	void ParaView_PrintSurfModel(FILE *fp);

};