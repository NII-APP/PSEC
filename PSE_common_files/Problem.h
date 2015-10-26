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
	int number; //номер текущей задачи
	int type; // номер типа задачи 
	// 1 - упругая статика
	// 2 - расчет собственных частот и форм колебаний


	int nmat; //количество материалов
	MATPROP *mat;




	//в дальнейшем можно дополнить системой задания и хранения законов нагружения

};
