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
	int number; //номер текущей задачи
	
	FULLMODEL *pfm;
	SEMODEL *psem;

	//в дальнейшем можно дополнить системой задания и хранения законов нагружения
	void ReinitLVRV(int newnvect); //реинициализация памяти под векторы нагрузки/решения
	void CreateLV(); //отладочная фукнция
	void StaticSolve();
	void StatMultSolveExt(int nv, double **R);
	void AttachToFullModel();

	void OutParaviewModel();

	void FullDel();

	//векторы нагрузки и решения
	int nvect;
	double **LV;
	double **RV;
	int fl_LVRV_internal; // == 1 : вектора нагрузок и решений инициализированы внутри объекта FullModel (требуется удаление памяти) 


};
