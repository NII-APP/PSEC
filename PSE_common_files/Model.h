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

	int nmat; //количество материалов
	MATPROP *mat;

	SEMODEL *sem; //суперэлементная модель объекта
	FULLMODEL *fm; //полная конечно-элементная модель объекта



};