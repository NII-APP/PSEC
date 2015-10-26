#include "stdafx.h"

#pragma once


class MEM 
{
public:
	MEM();
	~MEM();
public:
	double curmem;

public://функции работы с памятью
	int* MEM_NEW(int *B, int k);
	long long int* MEM_NEW(long long int *B, int k);
	float* MEM_NEW(float *B,int k);
	double* MEM_NEW(double *B,int k);
	double* MEM_NEW(double *B,long long int k);
	
	int** MEM_NEW(int **B, int str, int stb);
	float** MEM_NEW(float **B, int str, int stb);
	double** MEM_NEW(double **B, int str, int stb);
		
	int*** MEM_NEW(int ***B, int n1, int n2, int n3);
	float*** MEM_NEW(float ***B, int n1, int n2, int n3);
	double*** MEM_NEW(double ***B, int n1, int n2, int n3);
	
	float**** MEM_NEW(float ****B, int n1, int n2, int n3, int n4);
	double**** MEM_NEW(double ****B, int n1, int n2, int n3, int n4);
	
	int* MEM_DEL(int *A, int k);
	int* MEM_DEL(int *A, long long int k);
	long long int* MEM_DEL(long long int *B, int k);
	float* MEM_DEL(float *A, int k);
	double* MEM_DEL(double *A, int k);
	double* MEM_DEL(double *A, long long int k);
	
	int** MEM_DEL(int **A, int str, int stb);
	float** MEM_DEL(float **A, int str, int stb);
	double** MEM_DEL(double **A, int str, int stb);
		
	int*** MEM_DEL(int ***A, int n1, int n2,int n3);
	float*** MEM_DEL(float ***A, int n1, int n2,int n3);
	double*** MEM_DEL(double ***A, int n1, int n2,int n3);

	float**** MEM_DEL(float ****A, int n1, int n2, int n3, int n4);
	double**** MEM_DEL(double ****A, int n1, int n2, int n3, int n4);

public: //функции обнуления массивов
	int QQZERO (int *A, int k);
	int QQZERO (long long int *A, int k);
	int QQZERO (float *A, int k);
	int QQZERO (double *A, int k);
	int QQZERO (double *A, long long int k);
	
	int QQZERO (int **A, int k, int n);
	int QQZERO (float **A, int k, int n);
	int QQZERO (double **A, int k, int n);

	// работа с графом

	void GRNDEL(int ***NDEL_, int **IND, int NEL, int NN );
	void GenELGR(int ***GREL_, int NEL, int NN, int minnodes, int **NDEL, int **IND);
	void GenNDGR(int ***GRND_, int NEL, int NN, int NI, int **NDEL, int **IND);
	void sort (int *S,int s);

	//перенумерация и подсчет размеров профиля
	void GENRCM(int NN, int NI, int NEL, int NORT, double *CRD, int **IND, int *INVP, int *size_before, int *size_after, int *width_before, int *width_after, int *nbreaks );
	int RCM (int **GR,int NI, int stnode, int *PERM, int *INVP);
	void FindNOTNUMBERED (int **GR,int *MASK,int *TMP, int *temp,int *DEG,int uz);
	int FindStartingNode(double *CRD, int NI, int NORT);
	int SIZE_PROFILE (int **G, int NN, int NORT); //возвращает целое число мегабайт
	int WIDTH_PROFILE(int **G, int NN, int NORT); //максимальная ширина профиля, число уравнений

	void RenumberIND (int **IND, int NEL, int NI, int *INVP);
	void RenumberVect (int *VEK, int *INVP, int NI, int NORT);
	void RenumberVect (float *VEK, int *INVP, int NI, int NORT);
	void RenumberVect (double *VEK, int *INVP, int NI, int NORT);
};