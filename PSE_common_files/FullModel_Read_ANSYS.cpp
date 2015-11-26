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
	int nort_oldtypes[500];
	int similarcrdflag[3];
	double similarcrd[3];

	for (i=0; i<500; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
		nort_oldtypes[i] = -1;
	}

	renumber_types[182] = 5;
	renumber_types[187] = 24;
	
	nort_oldtypes[182] = 2;
	nort_oldtypes[187] = 3;

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
			pch+=7;
			sscanf(pch,"%ld",&nmat);
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			pch = strstr(str1,"MATERIALS");
			if (pch == NULL)
			{
				fl = 0;
				break;
			}
		}
	}

	//подсчет количества типа элементов
	fl = 1;
	neltypes = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"ET,");
		if ( pch != NULL )
		{
			pch += 3;
			sscanf(pch,"%ld",&neltypes);
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			pch = strstr(str1,"ET");
			if (pch == NULL)
			{
				fl = 0;
				break;
			}
		}
	}

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

	
	double *CRDTMP;

	CRDTMP = new double[NN*3]; // временный массив с расчетом на максимальное число координат
	IND = new int*[NEL];
	MTR = new int [NEL];
	mat = new MATPROP[nmat];
	model_eltypes = new int[neltypes];
	
	fseek(fp1,0,0);

	//считывание координат узлов
	fl = 1;
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
				CRDTMP[in*3+iort] = atof(str1);
				iort++;
				fscanf(fp1,"%s",str1);
				CRDTMP[in*3+iort] = atof(str1);
				iort++;
				fscanf(fp1,"%s",str1);
				CRDTMP[in*3+iort] = atof(str1);
				iort++;
				
				in++;
				iort = 0;
			}
		}
	}

	//считывание материалов
	fl = 1;
	int imat = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"MP, EX,");
		if ( pch != NULL )
		{
			pch+=7;
			sscanf(pch,"%ld",&imat);
			pch = strstr(pch,",");
			pch++;
			mat[imat-1].E = atof(pch);
		}
		pch = strstr(str1,"MP, NUXY,");
		if ( pch != NULL )
		{
			pch+=9;
			sscanf(pch,"%ld",&imat);
			pch = strstr(pch,",");
			pch++;
			mat[imat-1].MU = atof(pch);
		}
		pch = strstr(str1,"MP, DENS,");
		if ( pch != NULL )
		{
			pch+=9;
			sscanf(pch,"%ld",&imat);
			pch = strstr(pch,",");
			pch++;
			mat[imat-1].RO = atof(pch);
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			pch = strstr(str1,"MATERIALS");
			if (pch == NULL)
			{
				fl = 0;
				break;
			}
		}
	}


	//считывание типов элементов
	fl = 1;
	int ieltype = 0;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"ET,");
		if ( pch != NULL )
		{
			pch+=3;
			sscanf(pch,"%ld",&ieltype);
			pch = strstr(pch,"PLANE182");
			if ( pch != NULL )
			{
				model_eltypes[ieltype-1] = 182;
			}
		}
		pch = strstr(str1,"--------");
		if ( pch != NULL )
		{
			pch = strstr(str1,"ET");
			if (pch == NULL)
			{
				fl = 0;
				break;
			}
		}
	}

	// проверка числа рабочих координат, по типам элементов. 
	// В одной модели должны быть только типы элементов с одинаоквым числом рабочих координат
	KORT = nort_oldtypes[ model_eltypes[0] ]; //как минимум 1 тип элемента должен быть задан обязательно
	for (i=1; i<neltypes; i++)
	{
		if ( nort_oldtypes[ model_eltypes[i] ] != KORT )
		{
			printf("\ninappropriate mix of element types by NORT number....\n");
			break;
		}
	}

	for (i=0; i<KORT; i++)
	{
		crdrenum[i] = i;
	}

	//если задействованы менее 3х координат
	//осуществляется анализ, какие из 3х координат используются
	for (i=0; i<3; i++) 
	{ 
		similarcrdflag[i] = 0; //0 - координата не используется, 1 - координата используется  
		similarcrd[i] = CRDTMP[i]; //координаты для первого узла
	}
	double relcrdeps = 1.e-15;

	for (i=1; i<NN; i++)
	{
		for (j=0; j<3; j++)
		{
			if ( fabs( CRDTMP[i*3+j] ) > relcrdeps )
			{
				if ( fabs( (similarcrd[j] - CRDTMP[i*3+j])/CRDTMP[i*3+j] ) > relcrdeps )
				{
					similarcrdflag[j] = 1;
				}
			}
			else
			{
				if ( fabs( similarcrd[j] - CRDTMP[i*3+j] ) > relcrdeps )
				{
					similarcrdflag[j] = 1;
				}
			}
		}
	}
	KORT = 0;
	for (i=0; i<3; i++)
	{
		if ( similarcrdflag[i] == 1 ) KORT++;
	}

	NNE = NN*KORT;

	//копирование рабочих координат в окончательный используемый далее массив координат
	if ( KORT == 3 )
	{
		CRD = CRDTMP;
	}
	else
	{
		CRD = new double[NNE];
		int iort = 0;
		for (j=0; j<3; j++)
		{
			if ( similarcrdflag[j] == 1 )
			{
				for (i=0; i<NN; i++)
				{
					CRD[ i*KORT + iort ] = CRDTMP[ i*3 + j ];
				}
				iort++;
			}
		}
		delete []CRDTMP;
	}
			
		


	//считывание матрицы индексов
	fl = 1;
	int iel = 0, writtennodes = 0;
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
				writtennodes = atoi(strtmp);
				fscanf(fp1,"%s",strtmp);
				fscanf(fp1,"%s",strtmp);

				IND[iel][0] = eltype;
				IND[iel][1] = nodenumber_newtypes[eltype]; 
				for (j=0; j<writtennodes; j++)
				{
					fscanf(fp1,"%s",strtmp);
					if (j < IND[iel][1]) IND[iel][j+2] = atoi(strtmp) - 1;
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

	//!!!!провверить: при экспорте из NX нумерация узлов может требовать другой перенумерации
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

	//переход внутрь шага нагружения для последующего считывания закреплений-нагрузок
	//может быть добавлен анализ числа шагов нагружения и их видов
	fl = 1;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"Boundary Conditions");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}
	}

	//считывание закреплений
	fl = 1;
	neltypes = 0;
	char *pch2;
	while (fl == 1)
	{
		fgets(str1,255,fp1);

		pch = strstr(str1,"Loads");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}
		pch = strstr(str1,"FINISH");
		if ( pch != NULL )
		{
			fl = 0;
			break;
		}

		pch = strstr(str1,"D,");
		if ( pch != NULL )
		{
			pch2 = NULL;
			pch+=2;
			sscanf(pch,"%ld",&in);
			in--;
			iort = -1;
			pch = strstr(str1,"UX,");	
			if ( pch != NULL )
			{
				pch+=3;
				iort = 0;
				pch2 = pch;
			}
			pch = strstr(str1,"UY,");	
			if ( pch != NULL )
			{
				pch+=3;
				iort = 1;
				pch2 = pch;
			}
			pch = strstr(str1,"UZ,");	
			if ( pch != NULL )
			{
				pch+=3;
				iort = 2;
				pch2 = pch;
			}
			if ( iort > (KORT-1) ) iort = -1;
			if (iort > -1 && pch2 != NULL)
			{
				FIX[in*KORT + iort] = 1;
				UFIX[in*KORT + iort] = atof(pch2);
			}
		}
		
	}

	fclose(fp1);

	delete [] perenum;
	delete [] tmp;
	delete [] model_eltypes;

	printf("\nAnsys model was loaded...\n\n");
}