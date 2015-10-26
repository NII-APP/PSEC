#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "FullModel.h"
#include "SEmodel.h"
#include "T_elastat.h"
#include "T_ssi_eigen.h"
#include "T_motionint_eigf.h"



class PROBLEM
{
public:
	PROBLEM(void);
public:
	~PROBLEM(void);

public:
	MEM *MM;

	char name[256];
	int number; //����� ������� ������
	int type; // ����� ���� ������ 
	// 1 - ������� �������
	// 2 - ������ ����������� ������ � ���� ���������


	int nmat; //���������� ����������
	MATPROP *mat;




	//� ���������� ����� ��������� �������� ������� � �������� ������� ����������

};
