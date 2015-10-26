#pragma once

typedef struct   /* ������ � ���������� */
{
	double E;
	double MU;
	double RO;
	double SPL;
	double AL;
} MATPROP;

typedef struct              /* ������ � ����� �������������� */
{
	double S[6];
	double E[6];
	double VP;
} INTPOINT;

//typedef struct              /* ������ � ������ */
//{
//	char name[256];
//	int number; //����� ������� ������
//	int type; // ����� ���� ������ 
//	// 1 - ������� �������
//	// 2 - ������ ����������� ������ � ���� ���������
//
//
//	int nmat; //���������� ����������
//	MATPROP *mat;
//
//
//	//� ���������� ����� ��������� �������� ������� � �������� ������� ����������
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
//	int NEL; //���������� ���������
//	int NN; //���������� �����
//	int **IND; // ������� ��������
//	double *CRD; // ������ ���������
//	int *MAT; // ������ ����������
//	int *FIX; // ������ �����������
//	int *UFIX; // ������ �������� �����������
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