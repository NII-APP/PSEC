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
	CEL *el; //���������� �������� ���������

	int loadtype;

	char name[256];
//	char path[256];
	char netpath[256];
	char pathmatr[256]; // ���� � ������ � ��������� ������� �� ��
	char pathmain[256]; // ���� � ������ � ��������� ������� �� ������ � �� ��������� � �����

	int nmat; //���������� ����������
	MATPROP *mat;

	INTPOINT *P;
	int *envP; //������ ������������ ��������� ����� ������� � ������� ����� ��������������, ��������������� ������� ��������� ��������
	int nIntP;

	int KORT;
	int NEL; //���������� ���������
	int NN; //���������� �����
	int NNE; //���������� �������� �������
	int **IND; // ������� ��������
	double *CRD; // ������ ���������
	int *MTR; // ������ ����������
	int *FIX; // ������ �����������
	double *UFIX; // ������ �������� �����������
	double *Mdiag; //������������ ������� ����

	//�����
	int **NDEL;
	int **NDGR;
	int **ELGR;

	//������ ������ � ������ �� ��� ��������������� �����
	int *LIST_EL;
	int nlistel;
	int *MASK_EL;


	//����������� ������
	FACE *surfmodel;
	int nfaces; //���������� ������������� ������
	
	// ������� �������� � �������
	int nvect;
	double **LV;
	double **RV;
	int fl_LVRV_internal; // == 1 : ������� �������� � ������� ���������������� ������ ������� FullModel (��������� �������� ������)

	// ����� ����� ���������� ��� � ������ �����������
	FDPOINT *pforce; //����� ���������� ���
	int npforce;
	FDPOINT *pdisp; //����� ������� �����������
	int npdisp;


	void ReadFullmodel();
	void ReadFromUzor();
	void ReadFromAnsys_FullGrid();
	void ReadFromUzor_NewFormat();

	void InitIntPoint();
	void AttachElement(CEL *el, int elnum);
	void AttachElement(CEL *el, FACE *face); //������ ��� �������������� ��������
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