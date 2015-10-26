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

	int NS; //���������� ��������� �����, ����������� �� ������� �������
	int NI; //���������� ����������� �����
	int NN; //���������� �����
	int NSE; //���������� ���������, ��������������� NS
	int NIE; //���������� ���������, ��������������� NI
	int NNE;
	int NEL; //���������� ������������ ���������
	
	int KORT;

	SEstruct *sedat;
	//PROBLEMPROP *prob;
	MATPROP *pmat;


	INTPOINT *P;
	int *envP; //������ ������������ ��������� ����� ������� � ������� ����� ��������������, ��������������� ������� ��������� ��������
	int nIntP;

	//�������
	double *STII;
	double **STIS;
	double **STISTR;
	//double **STIS_1;
	double *STSS;
	long long int SSlength;
	long long int *STII_ENV;
	int **STIS_ENV;

	//��� ������ PCG
	double **STIINZ;
	double **STIIPC;
	int **STII_ENVNZ;
	int **STII_ENVPC;

	//�������� ������
	int **IND;
	double *CRD;
	int *MTR;
	int *FIX;
	double *UFIX;
	int *REN; //���������� ������ �����
	int *REE; //���������� ������ ���������

	//������� �������� ��� � �����������

	double **FORSE;
	double **DISP;
	int nvect; //���������� �������� ��� � ����������� ��� �������������� �������

	//������� ��������� � ��������������� (��� ������ ��)
	double **KK;
	//������������� ���������� ������� ��������� � ���������������
	int *ThreadQueue;
	int inewtrnum;
	HANDLE EvNotEmptyQueue;
	int Queuesize; //������ ���� ���� �� �� 1 ������ nthread
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
	void STIINZenv (int flagfull); //0 - ������ II �����, 1 - ������ �������
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
	void DeleteLevel_1(); //�������� ������
	void DeleteLevel_2(); //�������� �������� ������ ������ � �������
	void DeleteLevel_3(); //�������� ���������� �� ������ ��������������
	void DeleteLevel_4(); //�������� ������� ����������	
	
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