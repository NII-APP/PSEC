#pragma once

typedef struct   /* данные о материалах */
{
	double E;
	double MU;
	double RO;
	double SPL;
	double AL;
} MATPROP;

typedef struct              /* данные о точке интегрирования */
{
	double S[6];
	double E[6];
	double VP;
} INTPOINT;

//typedef struct              /* данные о задаче */
//{
//	char name[256];
//	int number; //номер текущей задачи
//	int type; // номер типа задачи 
//	// 1 - упругая статика
//	// 2 - расчет собственных частот и форм колебаний
//
//
//	int nmat; //количество материалов
//	MATPROP *mat;
//
//
//	//в дальнейшем можно дополнить системой задания и хранения законов нагружения
//	
//} PROBLEMPROP;
//
//typedef struct
//{
//	int SElevelsnum;
//	int *numSEbylevels;
//
//} SEMODEL;


//typedef struct
//{
//	int NEL; //количество элементов
//	int NN; //количество узлов
//	int **IND; // матрица индексов
//	double *CRD; // массив координат
//	int *MAT; // массив материалов
//	int *FIX; // массив закреплений
//	int *UFIX; // массив заданных перемещений
//
//} FULLMODEL;

//
//typedef struct             
//{
//	char IPname[256];
//	char WorkFolder[256];
//
//	
//} CLIENTDATA;