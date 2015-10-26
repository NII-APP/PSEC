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
	int number; //номер текущей задачи
	
	FULLMODEL *pfm;
	SEMODEL *psem;

	ELASTAT *pelast;


public:

	int endsolflag;

	int nform, nform_old, nf, nf_old, nfconv, nur, ITL, ITLlastshift, NEW_RAND;
	int NEigen;
	int *tmpi, *tmpi_old, *tmpi_old_2; //массив нумерации СЧ в порядке 
	//возрастания на текущей, предыдущей итерации и 2 итерации назад
	int *low_err;
//	int *NDL;
	double *tmp;
	double err_form,err_freq,err_freq_old,err_form_old, eps_form, eps_freq, *err_freq_mas,*err_form_mas;
	double **KK;//матрица жесткости в подпространстве
	double *M;//диагональ матрицы масс полной модели
	double **MM;//матрица масс в подпространстве
	double *LM;//диагональная матрица собственных чисел в подпространстве
	double *LMconv; //массив СЧ, сошедшихся до требуемой точности на пролых итерациях шифтинга
	double *LMold;//диагональная матрица собственных чисел в подпространстве (предыдущая итерация до нормировки)
	double *LMold_2;//диаг матрица собственных значений 2 итерации назад
	double *CR; //вектор скоростей сходимости 
	double *CRold; //вектор скоростей сходимости 1 итерацию назад 
	double *ALFR; //вектор коэффициентов релаксации в алгоритме ускорения сходимости
	double **Q;//матрица собственных форм в подпространстве
	double *UFmaxold;//нормы собственных векторов по максимальному элементу на прошлой итерации
	double *UFmax, **UF, **UFold, **UFconv;
	FILE *fp, *fp3, *Fp_NND;

	int nthread;
	int blocksize;

	double shifting;
	double shiftingold;
	int fl_isBlockConv; //=1, если сошлось nf частот в блоке
	double minfr_interval;
	int LastConvNum;
	double lqp1; //оценка для частоты номер nform+1
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