#pragma once

#include "StructData.h"
#include "MEM.h"



class CEL
{
//���� ��������� �� ���������
	//22 - 6 ������� ������� �����������
	//24 - 10 ������� ��������
	//25 - 20 ������� ����������

public:
	CEL(void);
public:
	~CEL(void);

public:
	MEM *MM;
	int ielobjnumber;
public:
	int eltype; //��� ��������
	int NORT; //����� ���������������� ��������� � ����
	int NORTfullcrd; //����� ���������������� ��������� � ���� � ������ ������� ��������� ������ (����� ��� ��������� ������� ������ ���������� �������)
	double lcrdlim[2]; //���������� ������� ��������� ��������� ��������� ������ ��������
	int NDOF; //����� �������� ������� � ����
	int NN; //����� ����� � ��������
	int NNE; //����� �������� ������� � �������� 	NNE = NN*NDOF;
	int NNORT; //����� ��������� ����� � �������� 	NNORT = NN*NORT;
	int NDEF; // ����� ����������� ����������
	int NSTRS; // ����� ����������� ����������
	
	int *ind; //������ ������� ����� �������� (��������� �� ��������������� ������ ������� ��������)
	double *fullcrd; //��������� �� ������ ������ ��������� ����� ������

	double **STIF; // ������� ��������� ��������
	double **MASS; //������� ����

	double **B; // ������� ����������
	double **D; //������� ������� ��� �������� ��������

	//��������� ���������� ����� �������������� �� ������
	double **vlcrdint; //[����� �����][����� ����������]
	int nvlcrd; //���������� ��������� ��������� ��� �������������� �� ������ (��� ���������, ��������, 4)
				//������ ����������� ������� ������ �� ����������� ����������� NORT (��� ��������� 1 ����������� ��������)
	int NPINT; //���������� ����� �������������� �� ������
	double *Wvint;//������� ������������ � ������ �������������� �� ������

	//��������� ���������� ����� �������������� �� ������
	double ***faceloccrdint; //[�����][����� �����][����� ����������]
	int Nface; //���������� ������
	int NNface; //���������� ����� �� �����
	int Nfp; //���������� ����� �������������� �� �����
	int nflcrd; //���������� ��������� ��������� ��� �������������� �� �����
	double **Wfint;//������� ������������ � ������ �������������� [�����][����� �����]
	int *faceorientcrd_num; //����� ���. ����������, ������� ��� ����� �������� ����������
	double *faceorientcrd_val; //�������� ��������� ����������, ������� ��������� �� �����

	//������� ������ �������������� ����� � ��������� ��������� ������
	int **FaceNodes; //[�����][����� ���� � ���. ��������� ��������] = 1, ���� ������ ���� ������������ � �����
	int **FaceNodesList; //[�����][����� ���� � ���. ��������� �����] = ���������� ������ ���� � ��������� �������� 

	//������� ����� � ��������������
	double **FFv;//������� �������� ������� ����� � ������ ��������-�� �� ������
	double ***FFf;//������� �������� ������� ����� � ������ ��������-�� �� ������
	double ***dFFvloc;	//������� �������� ����������� ������� ����� �� ��������� ����������� � ������ ��������-�� �� ������;
					//[����� ���][����� ���������� ������ �����������][����� ������� �����]
	double ****dFFfloc;	//������� �������� ����������� ������� ����� �� ��������� ����������� � ������ ��������-�� �� �����;
					//[�����][����� ���][����� ���������� ������ �����������][����� ������� �����]
	
	

	int *GLOBNE; //������ ���������� ������� �������� �������
	
	double VEL;

	INTPOINT *P; //������ ���������� ����� ��������������
	MATPROP *material; //��������� ���������

	//�����

	int isVPbelowzero; //���� ���������� ��������� �������������� ������ ����� ��������������
	int isVELcalculated; //����������, ����������� �� ����� ��������
	int isInitialized;

public:

	void Initialize(int ieltype);
	void GeneralMemInit();
	void InitGLOBNE();

	//������� �������� ���������� �� ����� ���������
	void GeneralParameters5();
	void GeneralParameters22();
	void GeneralParameters24();
	void GeneralParameters25();

	//������� ��������� ��������� ����� �������������� �� ������ � ������� �������
	void Init_VolIntPointLocCrd_5();
	void Init_VolIntPointLocCrd_22();
	void Init_VolIntPointLocCrd_24();
	void Init_VolIntPointLocCrd_25();

	//������� ��������� ��������� ����� �������������� �� ������ � ������� �������
	void Init_FacePointLocCrd_25();

	//����� �����, ����������� �������
	void Init_FaceNodes_24(); //������� �������� ������ ��������������� ����� ������
	void Init_FaceNodes_25(); //������� �������� ������ ��������������� ����� ������
	int IdentifyFace(int *facenodes, int nfn);//����������� ������������ ����� (���������� ���� � ���� �����)
	double FaceArea(int iface);
	void FaceNodeNums(int iface, int *NodeNums, int *NNf);

	// ���������� �������� ������� �����
	void FF_all(double *F, double *lcrd);
	void FF_5(double *F, double *lcrd);
	void FF_22(double *F, double *lcrd);
	void FF_24(double *F, double *lcrd);
	void FF_25(double *F, double *lcrd);

	//���������� ����������� ������� ����� �� ��������� �����������
	void dFF_all(double **dF, double *lcrd);
	void dFF_5(double **dF, double *lcrd);
	void dFF_22(double **dF, double *lcrd);
	void dFF_24(double **dF, double *lcrd);
	void dFF_25(double **dF, double *lcrd);

	//����������� ������� �����
	void Jacoby(double **dFloc, double **J);
	void DetJacoby(double **J, double *detJ);
	void InvJacoby(double **J, double detJ);
	int GetLocCrd(double *globcrd, double *loccrd, double *FF); //����������� ��������� ��������� ����� �� ���������� �����������

	//������ ������� ��������� � ����� �������������� � ����������� � ����� ������� ���������
	void CalcBGrad_3D(double **invJ, double **dFF);
	void CalcSTIFPoint(INTPOINT *P); //���������� ������� ��������� ��� ����� ��������������

	int DINT(); //���������� ������� ������� �������

	void ELSTIF(); //���������� ������� ���������
	void ELMASS(); //���������� ������� ����
	
	void shifting( double shift); //����� �� ������� ��� ���� ������� 2

	//�������������� ��������������
	void GetBoundBox(double *vmin, double *vmax);
	void GetCentr(double *vcentr);
};
