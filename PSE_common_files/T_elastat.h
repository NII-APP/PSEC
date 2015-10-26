#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "FullModel.h"
#include "SEmodel.h"

class ELASTAT
{
public:
	ELASTAT(void);
public:
	~ELASTAT(void);

public:
	MEM *MM;

	char name[256];
	int number; //����� ������� ������
	
	FULLMODEL *pfm;
	SEMODEL *psem;

	//� ���������� ����� ��������� �������� ������� � �������� ������� ����������
	void ReinitLVRV(int newnvect); //��������������� ������ ��� ������� ��������/�������
	void CreateLV(); //���������� �������
	void StaticSolve();
	void StatMultSolveExt(int nv, double **R);
	void AttachToFullModel();

	void OutParaviewModel();

	void FullDel();

	//������� �������� � �������
	int nvect;
	double **LV;
	double **RV;
	int fl_LVRV_internal; // == 1 : ������� �������� � ������� ���������������� ������ ������� FullModel (��������� �������� ������) 


};
