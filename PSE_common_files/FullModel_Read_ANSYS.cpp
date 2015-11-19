#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <string.h>

void FULLMODEL::ReadFromAnsys_FullGrid()
{	
	int i,j,tmpi,*perenum,crdrenum[3],neltypes,*model_eltypes,eltype;
	double tmpd,cx,cy,cz,scale;
	char str[256],strtmp[256];
	FILE *fp1, *fp2;
	char str1[256],str2[256],str3[256];

	//массив соответствия номеров типов элементов в ANSYS и PSE (аналогично PARAVIEW)
	int renumber_types[500];
	int nodenumber_newtypes[500];

	for (i=0; i<500; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}

	renumber_types[187] = 24;
	renumber_types[182] = 5;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[5] = 3;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;

	printf("\nModel loading...\n");

	sprintf(str1,"%s\\ans_model.txt",pathmain);
	fp1 = fopen(str1,"r");

	fscanf(fp1,"%s",strtmp);
	KORT = atoi(strtmp);
	fscanf(fp1,"%s",strtmp);
	NN = atoi(strtmp);
	fscanf(fp1,"%s",strtmp);
	NEL = atoi(strtmp);

	NNE = NN*KORT;

	fscanf(fp1,"%s",strtmp);
	scale = atof(strtmp);

	fscanf(fp1,"%s",strtmp);
	nmat = atoi(strtmp);
	mat = new MATPROP[nmat];
	for (i=0; i<nmat; i++)
	{
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		mat[i].RO = atof(strtmp)*scale*scale*scale;
		fscanf(fp1,"%s",strtmp);
		mat[i].E = atof(strtmp)*scale*scale;
		fscanf(fp1,"%s",strtmp);
		mat[i].MU = atof(strtmp);
	}

	fscanf(fp1,"%s",strtmp);
	neltypes = atoi(strtmp);

	model_eltypes = new int[neltypes];
	for (i=0; i<neltypes; i++)
	{
		fscanf(fp1,"%s",strtmp);
		model_eltypes[i] = atoi(strtmp);
	}	

	for (i=0; i<KORT; i++)
	{
		fscanf(fp1,"%s",strtmp);
		crdrenum[i] = atoi(strtmp);
	}

	CRD = new double[NNE];
	
	for (i=0; i<NN; i++)
	{
		fscanf(fp1,"%s",strtmp);
		for (j=0; j<KORT; j++)
		{
			fscanf(fp1,"%s",strtmp);
			CRD[ i*KORT + crdrenum[j] ] = atof(strtmp)/scale;
		}
	}
	
	IND = new int*[NEL];
	MTR = new int [NEL];
	
	for (i=0; i<NEL; i++)
	{
		fscanf(fp1,"%s",strtmp);
		eltype = atoi(strtmp) - 1;
		eltype = renumber_types[ model_eltypes[ eltype ] ];
		IND[i] = new int [ nodenumber_newtypes[eltype] + 2 ];

		fscanf(fp1,"%s",strtmp);
		MTR[i] = atoi(strtmp)-1;

		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);
		fscanf(fp1,"%s",strtmp);

		IND[i][0] = eltype;
		IND[i][1] = nodenumber_newtypes[eltype]; 
		for (j=0; j<IND[i][1]; j++)
		{
			fscanf(fp1,"%s",strtmp);
			IND[i][j+2] = atoi(strtmp) - 1;
		}
	}

	//локальная перенумерация узлов
	perenum = new int[100];
	int *tmp = new int[100];

	for (i=0;i<20; i++)
	{
		perenum[i] = i;
	}

	//перенумерация локального порядка узлов для основных объемных элементов требует проверки
	perenum[0] = 0;
	perenum[1] = 1;
	perenum[2] = 2;
	perenum[3] = 3;
	perenum[4] = 4;
	perenum[5] = 5;
	perenum[6] = 6;
	perenum[7] = 7;
	perenum[8] = 8;
	perenum[9] = 9;
	perenum[10] = 10;
	perenum[11] = 11;
	perenum[12] = 16;
	perenum[13] = 17;
	perenum[14] = 18;
	perenum[15] = 19;
	perenum[16] = 12;
	perenum[17] = 13;
	perenum[18] = 14;
	perenum[19] = 15;

	for (i=0; i<NEL; i++)
	{
		for (j=0; j<IND[i][1]; j++)
		{
			tmp[j] = IND[i][j+2];
		}
		for (j=0; j<IND[i][1]; j++)
		{
			IND[i][j+2] = tmp[perenum[j]];
		}
	}

	//выделение памяти и формирование закреплений
	FIX = new int[NNE];
	UFIX = new double[NNE];
	for(i=0; i<NNE; i++) 
	{
		FIX[i] = 0;
		UFIX[i] = 0.0;
	}

	//считывание файла закреплений

	int nfix,ifix;
	fscanf(fp1,"%ld",&nfix);
	//требуется исправить код для учета заданных перемещений
	
	for(i=0; i<nfix; i++)
	{
		fscanf(fp1,"%s",strtmp);
		ifix = atoi(strtmp) - 1;
		
		for (j=0; j<KORT; j++)
		{
			FIX[ifix*KORT+j]=1;
		}
	}

	delete [] perenum;
	delete [] tmp;
	delete [] model_eltypes;

	printf("\nAnsys model was loaded...\n\n");
}


void FULLMODEL::ReadFromAnsys_inpFormat()
{	
	int i,j,tmpi,*perenum,crdrenum[3],neltypes,*model_eltypes,eltype,fl;
	double tmpd,cx,cy,cz,scale;
	char str[256],strtmp[256], *pch = NULL;
	FILE *fp1, *fp2;
	char str1[256],str2[256],str3[256];

	//массив соответствия номеров типов элементов в ANSYS и PSE (аналогично PARAVIEW)
	int renumber_types[500];
	int nodenumber_newtypes[500];

	for (i=0; i<500; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}

	renumber_types[187] = 24;
	renumber_types[182] = 5;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[5] = 3;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;

	printf("\nModel loading...\n");

	sprintf(str1,"%s\\ans_model.txt",pathmain);
	fp1 = fopen(str1,"r");

	//подсчет количества узлов
	fl = 1;
	NN = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"NBLOCK,");
		if ( pch != NULL )
		{
			fgets(str1,255,fp1);
			while (fl == 1)
			{
				fgets(str1,255,fp1);

				pch = strstr(str1,"-1,");
				if ( pch != NULL ) { fl = 0; break; }

				sscanf(str1,"%ld",&NN);
			}
		}
	}

	//подсчет количества материалов
	fl = 1;
	nmat = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"MP, EX,");
		if ( pch != NULL )
		{
			sscanf(pch,"%ld",&nmat);
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}
	}

	//подсчет количества типа элементов
	fl = 1;
	neltypes = 0;
	while (fl == 1)
	{
		//требуется проанализировать количество координат в решаемой задаче и количество степеней свободы
		fgets(str1,255,fp1);

		pch = strstr(str1,"ET,");
		if ( pch != NULL )
		{
			sscanf(pch,"%ld",&neltypes);
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}
	}

	//временное назначение числа рабочих координат
	KORT = 2;

	//подсчет количества элементов
	fl = 1;
	NEL = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);
		int t;
		pch = strstr(str1,"EBLOCK,");
		if ( pch != NULL )
		{
			fgets(str1,255,fp1);
			while (fl == 1)
			{
				fgets(str1,255,fp1);

				pch = strstr(str1,"-1");
				if ( pch != NULL ) { fl = 0; break; }

				sscanf(str1,"%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld ",
							 &t, &t, &t, &t, &t, &t, &t, &t, &t, &t, &NEL);
			}
		}
	}

	NNE = NN*KORT;

	CRD = new double[NNE];
	IND = new int*[NEL];
	MTR = new int [NEL];
	mat = new MATPROP[nmat];
	model_eltypes = new int[neltypes];
	
	fseek(fp1,0,0);

	//считывание координат узлов
	fl = 1;
	NN = 0;
	int in = 0,iort = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"NBLOCK,");
		if ( pch != NULL )
		{
			fgets(str1,255,fp1);
			while (fl == 1)
			{		
				fscanf(fp1,"%s",str1);

				pch = strstr(str1,"N");
				if ( pch != NULL ) { fl = 0; break; }

				fscanf(fp1,"%s",str1);
				fscanf(fp1,"%s",str1);
				fscanf(fp1,"%s",str1);
				CRD[in*KORT+iort] = atof(str1);
				iort++;
				fscanf(fp1,"%s",str1);
				CRD[in*KORT+iort] = atof(str1);
				iort++;
				fscanf(fp1,"%s",str1);
				
				in++;
				iort = 0;
			}
		}
	}

	mat[0].RO = 7.872000E-06;
	mat[0].E = 2.000000E+08;
	mat[0].MU = 2.500000E-01;

	model_eltypes[0] = 182;

	for (i=0; i<KORT; i++)
	{
		crdrenum[i] = i;
	}

	//подсчет количества элементов
	fl = 1;
	int iel = 0;
	in = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);
		int t;
		pch = strstr(str1,"EBLOCK,");
		if ( pch != NULL )
		{
			fgets(str1,255,fp1);
			while (fl == 1)
			{
				fscanf(fp1,"%s",str1);

				pch = strstr(str1,"-1");
				if ( pch != NULL ) { fl = 0; break; }

				fscanf(fp1,"%s",str1);
				MTR[iel] = atoi(str1) - 1;

				fscanf(fp1,"%s",strtmp);
				eltype = atoi(strtmp) - 1;
				eltype = renumber_types[ model_eltypes[ eltype ] ];
				IND[iel] = new int [ nodenumber_newtypes[eltype] + 2 ];

				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);

				IND[iel][0] = eltype;
				IND[iel][1] = nodenumber_newtypes[eltype]; 
				for (j=0; j<IND[iel][1]; j++)
				{
					fscanf(fp1,"%s",strtmp);
					IND[iel][j+2] = atoi(strtmp) - 1;
				}

				iel++;
			}
		}
	}
	

	//локальная перенумерация узлов
	perenum = new int[100];
	int *tmp = new int[100];

	for (i=0;i<20; i++)
	{
		perenum[i] = i;
	}

	//перенумерация локального порядка узлов для основных объемных элементов требует проверки
	perenum[0] = 0;
	perenum[1] = 1;
	perenum[2] = 2;
	perenum[3] = 3;
	perenum[4] = 4;
	perenum[5] = 5;
	perenum[6] = 6;
	perenum[7] = 7;
	perenum[8] = 8;
	perenum[9] = 9;
	perenum[10] = 10;
	perenum[11] = 11;
	perenum[12] = 16;
	perenum[13] = 17;
	perenum[14] = 18;
	perenum[15] = 19;
	perenum[16] = 12;
	perenum[17] = 13;
	perenum[18] = 14;
	perenum[19] = 15;

	for (i=0; i<NEL; i++)
	{
		for (j=0; j<IND[i][1]; j++)
		{
			tmp[j] = IND[i][j+2];
		}
		for (j=0; j<IND[i][1]; j++)
		{
			IND[i][j+2] = tmp[perenum[j]];
		}
	}

	//выделение памяти и формирование закреплений
	FIX = new int[NNE];
	UFIX = new double[NNE];
	for(i=0; i<NNE; i++) 
	{
		FIX[i] = 0;
		UFIX[i] = 0.0;
	}

	//считывание закреплений
	fl = 1;
	neltypes = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"COM,");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}

		pch = strstr(str1,"D,");
		if ( pch != NULL )
		{
			sscanf(pch,"%ld",&in);
			in--;
			iort = -1;
			strstr(pch,"UX,");	
			if ( pch != NULL )
			{
				iort = 0;
			}
			strstr(pch,"UY,");	
			if ( pch != NULL )
			{
				iort = 1;
			}
			/*strstr(pch,"UZ,");	
			if ( pch != NULL )
			{
				iort = 2;
			}*/
			if (iort > -1)
			{
				FIX[in*KORT + iort] = 1;
				sscanf(pch,"%e",&UFIX[in*KORT + iort]);
		}
		
	}


	delete [] perenum;
	delete [] tmp;
	delete [] model_eltypes;

	printf("\nAnsys model was loaded...\n\n");
}