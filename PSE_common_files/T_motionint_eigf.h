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
	int number; //номер текущей задачи

	bool is_initialized; //флаг, показывающий, инициализирован ли данный расчет
	
	FULLMODEL *pfm;


	//собственные формы
	int nform;
	int NNE;
	float **EF;
	
	//вектор узловых сил
	float *FORCE;

	//вектор узловых перемещений
	float *DISP;

	float *FREQ; //массив собственных частот в рад/c
	float *MODM; //модальные массы
	float *MODD; //модальные демпфирования

	// для интегрирования
	float *VST; //вектор состояния nform*2 - в конце шага по времени (текущее определяемое значение)
	float *VSTpi; //вектор состояния - в конце шага по времени (на предыдущей итерации)
	float *VST0; //вектор состояния на начало шага интегрирования
	float *qf; //силы приведенные к формам колебаний (средние на шаге, используется в аналитических формулах интегрирования)
	float *qf_ss; // силы приведенные к формам колебаний - на начало шага по времени
	float *qf_es; // силы приведенные к формам колебаний - на конец шага по времени

	// параметр сходимости
	float it_err;
	int it;
	float eps_it;
	float lvst; //квадрат модуля вектора состояния по компонентам перемещения, средний за 2 последние итерации
	float dvst; //квадрат модуля разности вектора состояния на последних 2х итерациях 
	
	//модальные вклады от единичной инерционной нагрзуки по осям координат
	float **axqf_in_1; //[x,y,z][nform]

	


	//в дальнейшем можно дополнить системой задания и хранения законов нагружения
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