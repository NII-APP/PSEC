#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"


void FULLMODEL::ReadFromUzor_NewFormat()
{	
	int i,j,k,k1,nn;
	int type, newtype;
	double tmpd;
	char str[256];
	char strtmp[256];
	FILE *fp;

	//массив соответствия номеров типов элементов в UZOR и PSE (аналогично PARAVIEW)
	int renumber_types[50];
	int nodenumber_newtypes[50];

	for (i=0; i<50; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}


	renumber_types[0] = 12;
	renumber_types[1] = 25;
	renumber_types[2] = 13;
	renumber_types[3] = -1;
	renumber_types[4] = 10;
	renumber_types[5] = 24;
	renumber_types[6] = -1;
	renumber_types[7] = -1;
	renumber_types[8] = 3;
	renumber_types[9] = 21;
	renumber_types[10] = 3;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;



	//считывание объема модели
	sprintf(str,"%s\\siz0000.dat",pathmatr);
	fp = fopen(str,"r");
	fscanf(fp,"%s",strtmp);
	NN = atoi(strtmp);
	fscanf(fp,"%s",strtmp);
	NEL = atoi(strtmp);
	fclose(fp);

	//считывание матрицы индексов
	
	IND = new int*[NEL];
	sprintf(str,"%s\\ind0000.dat",pathmatr);
	fp = fopen(str,"r");
	for (i=0; i<NEL; i++)
	{
		fscanf(fp,"%s",strtmp);
		type = atoi(strtmp);
		newtype = renumber_types[type];
		IND[i] = new int[nodenumber_newtypes[newtype]+2];
		IND[i][0] = newtype;
		IND[i][1] = nodenumber_newtypes[newtype];

		for (k1=0; k1<IND[i][1]; k1++)
		{
			fscanf(fp,"%s",strtmp);
			IND[i][2+k1] = atoi(strtmp) - 1;
		}
	}
	fclose(fp);

	if (IND[0][0] == 3 || IND[0][0] == 10 || IND[0][0] == 12 || IND[0][0] == 24 || IND[0][0] == 25 )
	{
		KORT = 3;
	}
	NNE = NN*KORT;

	//считывание массива координат
	CRD = new double[KORT*NN];
	sprintf(str,"%s\\crd0000.dat",pathmatr);
	fp = fopen(str,"rb");
	fread(CRD,sizeof(double),KORT*NN,fp);
	fclose(fp);

	//считывание массива материалов
	MTR = new int [NEL];
	sprintf(str,"%s\\mtr0001.sba",pathmain);
	fp = fopen(str,"rb");
	fread(MTR,sizeof(int),NEL,fp);
	fclose(fp);

	for (i=0; i<NEL; i++)
	{
		MTR[i] -= 1;
	}

	//считывание закреплений
	FIX = new int [KORT*NN];
	UFIX = new double [KORT*NN];

	sprintf(str,"%s\\fix0000.dat",pathmatr);
	fp = fopen(str,"r");
	for (i=0; i<NN; i++)
	{
		for (j=0; j<KORT; j++)
		{
			fscanf(fp,"%s",strtmp);
			tmpd = atof(strtmp);
			if ( abs(tmpd) > 1.e-21 )
			{
				FIX[i*KORT+j] = 1;
				UFIX[i*KORT+j] = tmpd;
				if ( abs(tmpd) < 2.e-20 )
				{
					UFIX[i*KORT+j] = 0.0;
				}
			}
			else
			{
				FIX[i*KORT+j] = 0;
				UFIX[i*KORT+j] = 0.0;
			}
		}
	}
	fclose(fp);

}


void FULLMODEL::ReadFromUzor()
{	
	int i,j,k,k1,nn;
	int type, newtype;
	double tmpd;
	char str[256];
	char strtmp[256];
	FILE *fp;

	//массив соответствия номеров типов элементов в UZOR и PSE (аналогично PARAVIEW)
	int renumber_types[50];
	int nodenumber_newtypes[50];

	for (i=0; i<50; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}


	renumber_types[0] = 12;
	renumber_types[1] = 25;
	renumber_types[2] = 13;
	renumber_types[3] = -1;
	renumber_types[4] = 10;
	renumber_types[5] = 24;
	renumber_types[6] = -1;
	renumber_types[7] = -1;
	renumber_types[8] = 3;
	renumber_types[9] = 21;
	renumber_types[10] = 3;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;



	//считывание объема модели
	sprintf(str,"%s\\siz0000.dat",pathmatr);
	fp = fopen(str,"r");
	fscanf(fp,"%s",strtmp);
	NN = atoi(strtmp);
	fscanf(fp,"%s",strtmp);
	NEL = atoi(strtmp);
	fclose(fp);

	//считывание матрицы индексов
	
	IND = new int*[NEL];
	sprintf(str,"%s\\ind0000.dat",pathmatr);
	fp = fopen(str,"r");
	for (i=0; i<NEL; i++)
	{
		fscanf(fp,"%s",strtmp);
		type = atoi(strtmp);
		newtype = renumber_types[type];
		IND[i] = new int[nodenumber_newtypes[newtype]+2];
		IND[i][0] = newtype;
		IND[i][1] = nodenumber_newtypes[newtype];

		for (k1=0; k1<IND[i][1]; k1++)
		{
			fscanf(fp,"%s",strtmp);
			IND[i][2+k1] = atoi(strtmp) - 1;
		}
	}
	fclose(fp);

	if (IND[0][0] == 3 || IND[0][0] == 10 || IND[0][0] == 12 || IND[0][0] == 24 || IND[0][0] == 25 )
	{
		KORT = 3;
	}
	NNE = NN*KORT;

	//считывание массива координат
	CRD = new double[KORT*NN];
	sprintf(str,"%s\\crd0000.dat",pathmatr);
	fp = fopen(str,"rb");
	fread(CRD,sizeof(double),KORT*NN,fp);
	fclose(fp);

	//считывание массива материалов
	MTR = new int [NEL];
	sprintf(str,"%s\\mtr0001.sba",pathmain);
	fp = fopen(str,"rb");
	fread(MTR,sizeof(int),NEL,fp);
	fclose(fp);

	for (i=0; i<NEL; i++)
	{
		MTR[i] -= 1;
	}

	//считывание закреплений
	FIX = new int [KORT*NN];
	UFIX = new double [KORT*NN];

	sprintf(str,"%s\\fix0000.dat",pathmatr);
	fp = fopen(str,"r");
	for (i=0; i<NN; i++)
	{
		for (j=0; j<KORT; j++)
		{
			fscanf(fp,"%s",strtmp);
			tmpd = atof(strtmp);
			if ( abs(tmpd) > 1.e-21 )
			{
				FIX[i*KORT+j] = 1;
				UFIX[i*KORT+j] = tmpd;
				if ( abs(tmpd) < 2.e-20 )
				{
					UFIX[i*KORT+j] = 0.0;
				}
			}
			else
			{
				FIX[i*KORT+j] = 0;
				UFIX[i*KORT+j] = 0.0;
			}
		}
	}
	fclose(fp);

}