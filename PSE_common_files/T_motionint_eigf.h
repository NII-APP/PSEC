#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "FullModel.h"

class MOTIONEIGF
{
public:
	MOTIONEIGF(void);
public:
	~MOTIONEIGF(void);

public:
	MEM *MM;

	char name[256];
	int number; //����� ������� ������

	bool is_initialized; //����, ������������, ��������������� �� ������ ������
	
	FULLMODEL *pfm;


	//����������� �����
	int nform;
	int NNE;
	float **EF;
	
	//������ ������� ���
	float *FORCE;

	//������ ������� �����������
	float *DISP;

	float *FREQ; //������ ����������� ������ � ���/c
	float *MODM; //��������� �����
	float *MODD; //��������� �������������

	// ��� ��������������
	float *VST; //������ ��������� nform*2 - � ����� ���� �� ������� (������� ������������ ��������)
	float *VSTpi; //������ ��������� - � ����� ���� �� ������� (�� ���������� ��������)
	float *VST0; //������ ��������� �� ������ ���� ��������������
	float *qf; //���� ����������� � ������ ��������� (������� �� ����, ������������ � ������������� �������� ��������������)
	float *qf_ss; // ���� ����������� � ������ ��������� - �� ������ ���� �� �������
	float *qf_es; // ���� ����������� � ������ ��������� - �� ����� ���� �� �������

	// �������� ����������
	float it_err;
	int it;
	float eps_it;
	float lvst; //������� ������ ������� ��������� �� ����������� �����������, ������� �� 2 ��������� ��������
	float dvst; //������� ������ �������� ������� ��������� �� ��������� 2� ��������� 
	
	//��������� ������ �� ��������� ����������� �������� �� ���� ���������
	float **axqf_in_1; //[x,y,z][nform]

	


	//� ���������� ����� ��������� �������� ������� � �������� ������� ����������
	void LoadEIGF();
	void InitIntegr();
	void InitConstModalDamp(float damp);
	void InitInertiaUnitModalForces();
	void AttachToFullModel(FULLMODEL *pf);

	void ResetModalForces_EndStep();
	void AddInertiaForces_EndStep(float *ax_val);
	void AddFDPointForces_EndStep();

	void MakeSingleTimeStep(float dt);
	void StartNewStep();
	void StartNewIt();
	bool IsItConverged();

	void EvaluateAllDOF();
	void EvaluateFDPoint_EndStep();
//	void OutParaviewModel();


	void TestCase_InertiaLoading();

	void FullDel();



};