#pragma once
#include "stdafx.h"
#include "StructData.h"
#include "EL.h"
#include "SE_structure.h"
#include "MEM.h"

#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <process.h>

class CSE
{
public:
	CSE(void);
public:
	~CSE(void);

public:
	double bignumber;
	MEM *MM;

public:
	int MAXELTYPES;

public:

	int NS; //количество граничных узлов, переходящих на верхний уровень
	int NI; //количество исключаемых узлов
	int NN; //количество узлов
	int NSE; //количество уравнений, соответствующих NS
	int NIE; //количество уравнений, соответствующих NI
	int NNE;
	int NEL; //количество составляющих элементов
	
	int KORT;

	SEstruct *sedat;
	//PROBLEMPROP *prob;
	MATPROP *pmat;


	INTPOINT *P;
	int *envP; //массив показывающий стартовый номер позиции в массиве точек интегрирования, соответствующий каждому конечному элементу
	int nIntP;

	//матрицы
	double *STII;
	double **STIS;
	double **STISTR;
	//double **STIS_1;
	double *STSS;
	long long int SSlength;
	long long int *STII_ENV;
	int **STIS_ENV;

	//для метода PCG
	double **STIINZ;
	double **STIIPC;
	int **STII_ENVNZ;
	int **STII_ENVPC;

	//исходные данные
	int **IND;
	double *CRD;
	int *MTR;
	int *FIX;
	double *UFIX;
	int *REN; //глобальные номера узлов
	int *REE; //глобальные номера элементов

	//матрицы векторов сил и перемещений

	double **FORSE;
	double **DISP;
	int nvect; //количество векторов сил и перемещений для одновременного решения

	//матрица жесткости в подпространстве (для поиска СЧ)
	double **KK;
	//многопоточное вычисление матрицы жесткости в подпространстве
	int *ThreadQueue;
	int inewtrnum;
	HANDLE EvNotEmptyQueue;
	int Queuesize; //должен быть хотя бы на 1 больше nthread
	int nthread;
	int blocksize;





public:

	void Initialize(SEstruct *ses);
	void InitIntPoint();
	void ReadSEinitials();
	void MainMatrCalc();
	void MainLoadCalc();
	void MainResCalc();
	void MainSubspaceCalc();
	
	void CalcSubspaceKM();
	void CalcSubspaceKMPar();

	void STIIenv ();
	void STISenv ();
	void StiffMatrix();
	void AttachElement(CEL *el, int elnum);

	void StifIIAssembling(CEL *el,int elnum);
	void StifISAssembling(CEL *el,int elnum);
	void StifSSAssembling(CEL *el,int elnum);
	void StifIIAssemblingFromSE(double *STSSse, int KU, int *NUR );
	void StifISAssemblingFromSE(double *STSSse, int KU, int *NUR );
	void StifSSAssemblingFromSE(double *STSSse, int KU, int *NUR );

	void StiffMatrixTransform();
	void MatrixFixing();
	void SScolTransform(int colstart,int colend);
	void IScolExpand(double *COLIS,int isscol);
	void SScolCalc(double *COLSS, double *COLIS, int jsscol);
	void SScolInsert(double *COLSS, int isscol);

	void LoadForse();
	void ForseTransform();
	void ForseFixing();
	void FORSEcolTransform(int colstart, int colend);
	void FORSEcolS(double *RV,int icol);

	void LoadDisp();
	void DispTransform();
	void DISPcolTransform(int colstart, int colend);

	void StifPrecondition();
	void STIINZenv (int flagfull); //0 - только II часть, 1 - полная матрица
	void StifIINZAssembling(CEL *el,int elnum,int flagfull);
	void StifIINZAssemblingSE(double *STSSse, int KU, int *NUR);
	void PrecondHolSol(double *B, double *C);
	void StifMultPC(double *C, double *B);
	void PCGSolve( double *U );

	void StiiWrite();
	void StiiRead();
	void StiiWriteNZPC();
	void StiiReadNZPC();
	void WriteNZfull();
	void ReadNZfull();
	void StisWrite();
	void StisRead();
	void StssWrite();
	void StssRead();
	void ForseRead();
	void ForseWrite();
	void ForseWriteSS();
	void DispRead();
	void DispWrite();
	void DispWriteSS();
	void IntPWrite();
	void IntPRead();
	void WriteSubspaceKM();

	void FULLDEL();
	void DeleteLevel_1(); //удаление матриц
	void DeleteLevel_2(); //удаление векторов правых частей и решений
	void DeleteLevel_3(); //удаление информации по точкам интегрирования
	void DeleteLevel_4(); //удаление базовой информации	
	
};

struct KKparDAT 
{
	int startiur;
	int endiur;

	int ithread;

	CSE *pcse;
	HANDLE start;
	HANDLE totalfin;
	HANDLE potok;
	DWORD uThreadID;

	int infcycleflag;
};