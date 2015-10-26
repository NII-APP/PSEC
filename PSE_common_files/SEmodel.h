#pragma once
#include "stdafx.h"
#include <windows.h>
#include <winbase.h>
#include "stdlib.h"
#include "MEM.h"
#include "SE_structure.h"
#include "SE.h"
#include "FullModel.h"

class SEMODEL
{
public:
	SEMODEL(void);
public:
	~SEMODEL(void);

public:
	MEM *MM;

	char name[256];
//	char path[256];
	char netpath[256];
	char pathmatr[256]; // ���� � ������ � ��������� ������� �� ��
	char pathmain[256]; // ���� � ������ � ��������� ������� �� ������ � �� ��������� � �����

//	int nmat; //���������� ����������
//	MATPROP *mat;

	FULLMODEL *pfm;

	int SEmaxLev;
	int SElevelsnum;
	int *numSEbylevels;
	SEstruct **SE;

	//��� ��������������� �������
	int **Levels; //��������� ��������������� ������� �� ��, ���������� ��� ������� ������ � ������� �� ������ �� ���������� ������ ����������� ������ �������
	FILE *fp_ase;


	

	int loadtype; //��� �������� ������
	// 1 - �� ��������� UZOR
	// 2 - �������������� ���������

	//������

	void ReadSEmodel();
	void ReadFromUzor();
	void ReadFromUzor_NewFormat();
	void ReadFromAutoSE();

	void AutoMLDivision_S();
	void AutoMLDivision_LevelOutChange(int ilev, int nse, int *RENlev, int NNlev, int nel, int **IND, int **INDnextlev, int *NNnextlev);
	void AutoMLD_L0(int *DS_);
	void AutoMLD_L0_EL(int *DS_);
	void AutoMLD_L1(int ilev, int maxboundary_size, int maxtotal_size, int **IND_OLD, int nel, int NNlev );
	//void FindNotNumNeib(int **GR, int *MASK, int iSE, int *NEIB, int *neib, int DS, int startnode);
	void FindNotNumNeib(int **GR, int *MASK, int *DEG, int iSE, int *NEIB, int *neib, int DS, int startnode);
	void FindSeBound(int ***SEBND_, int **GREL, int **NDEL, int NNlev, int nel, int ilev);
	void AutoDivision_PrintLevels_Paraview();

	void SE_struct_init(SEstruct *SE, int ilev, int ise);

	void CalcStifMatr();
	//void CalcLoadVect(double **LV, int nV);
	//void CalcResultVect(double **RV, int nV);
	void CalcLoadVect();
	void CalcResultVect();
	void CalcSubspaceMatrix(double **KK);

	void SetNewShiftSEM(double shifting);
	void SetOperationType(int type);

};