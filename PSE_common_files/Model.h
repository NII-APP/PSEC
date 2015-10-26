#pragma once
#include "stdafx.h"
#include "MEM.h"
#include "SEmodel.h"
#include "FullModel.h"



class MODEL
{
public:
	MODEL(void);
public:
	~MODEL(void);

public:
	MEM *MM;

	char name[256];
	char path[256];
	char netpath[256];

	int nmat; //���������� ����������
	MATPROP *mat;

	SEMODEL *sem; //��������������� ������ �������
	FULLMODEL *fm; //������ �������-���������� ������ �������



};