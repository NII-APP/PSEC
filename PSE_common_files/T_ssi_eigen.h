#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "FullModel.h"
#include "SEmodel.h"
#include "T_elastat.h"
#include "time.h"




class SSIEIGEN
{
public:
	SSIEIGEN(void);
public:
	~SSIEIGEN(void);

public:
	MEM *MMM;

	char name[256];
	int number; //����� ������� ������
	
	FULLMODEL *pfm;
	SEMODEL *psem;

	ELASTAT *pelast;


public:

	int endsolflag;

	int nform, nform_old, nf, nf_old, nfconv, nur, ITL, ITLlastshift, NEW_RAND;
	int NEigen;
	int *tmpi, *tmpi_old, *tmpi_old_2; //������ ��������� �� � ������� 
	//����������� �� �������, ���������� �������� � 2 �������� �����
	int *low_err;
//	int *NDL;
	double *tmp;
	double err_form,err_freq,err_freq_old,err_form_old, eps_form, eps_freq, *err_freq_mas,*err_form_mas;
	double **KK;//������� ��������� � ���������������
	double *M;//��������� ������� ���� ������ ������
	double **MM;//������� ���� � ���������������
	double *LM;//������������ ������� ����������� ����� � ���������������
	double *LMconv; //������ ��, ���������� �� ��������� �������� �� ������ ��������� ��������
	double *LMold;//������������ ������� ����������� ����� � ��������������� (���������� �������� �� ����������)
	double *LMold_2;//���� ������� ����������� �������� 2 �������� �����
	double *CR; //������ ��������� ���������� 
	double *CRold; //������ ��������� ���������� 1 �������� ����� 
	double *ALFR; //������ ������������� ���������� � ��������� ��������� ����������
	double **Q;//������� ����������� ���� � ���������������
	double *UFmaxold;//����� ����������� �������� �� ������������� �������� �� ������� ��������
	double *UFmax, **UF, **UFold, **UFconv;
	FILE *fp, *fp3, *Fp_NND;

	int nthread;
	int blocksize;

	double shifting;
	double shiftingold;
	int fl_isBlockConv; //=1, ���� ������� nf ������ � �����
	double minfr_interval;
	int LastConvNum;
	double lqp1; //������ ��� ������� ����� nform+1
	double lqp1_sum;
	int nlqp1;
	int consecutive_convfreq;
	int nfconvold;
	int *ortoseqv;

	int Init(int NEigen, double EpsEigen, int nf);
	int MainSSI();
	void UFGetMax();
	void NewForms();
	void JacobySweep();
	void JacobyOFF (double *KKoff, double *MMoff);
	void Jacoby();
	void GaussShmidt_2();
	void GaussShmidt_3();
	void GaussShmidt_2_SpecSeqv();
	void FormSort();
	void GrammShmidtCyckl();
	void NormalizeEigVect();
	void CollinearDetection();
	void GrammShmidtCyckl_SpecSeqv();
	void StandartStartingVectors();
	void PrintRez_EigenSSI_iter();
	void PrintRez_EigenSSI_full();
	void DeleteMem();
	void CopyConvForms();
	void DefineNewShift();
	void SetNewShiftSEM();	
	void ErrorCalc();
	void ConvAccelCalc();
	void GrammShmidtCheck();
	void ConvCheck();

	void statusbar(char *str);
	

};