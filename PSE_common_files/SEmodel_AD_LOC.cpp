#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"
#include "SE_structure.h"

void SEMODEL::AutoMLDivision_LevelOutChange(int ilev, int nse, int *RENlev, int NNlev, int nel, int **IND, int **INDnextlev, int *NNnextlev)
{
	// ilev - номер уровня, который нужно сформировать и вывести при данном запуске
	// Levels - массив принадлежности элементов (суперэлементов уровня (ilev-1) ) новым суперэлементам формируемого уровня ilev
	// nse - количество новых суперэлементов уровня ilev
	// RENlev - массив перенумерации к глобальным номерам исходной модели от локальной нумерации узлов в пределах уровня (ilev-1) в целом
	// NNlev - количество узлов в уровне (ilev-1) в целом
	// nel - количество элементов (суперэлементов уровня (ilev-1) ) составляющих уровень ilev в целом
	// IND - матрица индексов уровня (ilev-1) в целом
	// INDnextlev - матрица индексов для СЭ формируемого уровня ilev (должна быть внешним образом выделена память под первый указатель)
	// NNnextlev - количество узлов в формируемом уровне ilev, возвращаемая величина

	int i,j,iel,in,ise,nblev,k,k1;
	char str[256];
	FILE *fp;
	
	//1 - определение граничных узлов
	int *fl_b; //флаг граничных узлов (== -2)
	fl_b = new int[NNlev];
	for (i=0; i<NNlev; i++)
	{
		fl_b[i] = -1;
	}

	for (iel = 0; iel < nel; iel++)
	{
		for (j = 0; j < IND[iel][1]; j++)
		{
			in = IND[iel][j+2];
			if ( fl_b[ in ] == -1 ) //узел встречается впервые
			{
				fl_b[ in ] = Levels[ilev][iel];
			}
			if ( fl_b[ in ] > -1 ) // узел уже встречался в одном СЭ
			{
				if ( fl_b[ in ] != Levels[ilev][iel] ) // если теперь этот узел встретился в другом суперэлементе
				{
					fl_b[ in ] = -2; //узел встретился более, чем в одном СЭ, значит является граничным
				}
			}
		}
	}
	nblev = 0; //количество граничных узлов в уровне (количество узлов в следующем уровне)
	for (i=0; i < NNlev; i++)
	{
		if ( fl_b[i] == -2 ) nblev++;
	}

	//2 - составление массива перенумерации граничных узлов (перенумерация от номеров узлов следующего уровня к номерам узлов текущего уровня)
	int *REN_B_LEV,*REN_B_LEV_INV,inb;
	REN_B_LEV = new int [nblev];
	REN_B_LEV_INV = new int [NNlev];

	inb = 0;
	for (i=0; i < NNlev; i++)
	{
		REN_B_LEV_INV[i] = -1;
		if ( fl_b[i] == -2 )
		{
			REN_B_LEV_INV[i] = inb;
			REN_B_LEV[inb] = i;
			inb++;
		}
	}

	int *REE,*REN,*REN_B;
	bool *fl_n;
	int **INDse;
	REE = new int [nel];
	REN = new int [NNlev];
	REN_B = new int [NNlev];
	int iree = 0,ireni = 0, irenb = 0;

	fl_n = new bool[NNlev];
	for (i=0; i<NNlev; i++)
	{
		fl_n[i] = false;
	}

	// создание информационной структуры о суперэлементах данного уровня
	numSEbylevels[ilev] = nse;
	SE[ilev] = new SEstruct [numSEbylevels[ilev]];


	for (ise = 0; ise < nse; ise++)
	{
		//экономичная "очистка" массивов перенумерации
		iree = 0;
		ireni = 0;
		irenb = 0;
		SE_struct_init(&SE[ilev][ise], ilev, ise);
		for (iel=0; iel<nel; iel++)
		{
			if ( Levels[ilev][iel] == ise )
			{
				//запоминаем глобальный номер составляющего элемента в рамках уровня
				REE[iree] = iel;
				iree++;

				for (j=0; j <IND[iel][1]; j++)
				{
					in = IND[iel][j+2];
					if ( fl_n[in] == false )
					{
						if ( fl_b[ in ] != -2 ) //узел является внутренним
						{
							REN[ireni] = in;
							ireni++;
						}
						else
						{
							REN_B[irenb] = in;
							irenb++;
						}
						fl_n[in] = true;
					}
				}
			}
		}
		SE[ilev][ise].NEL = iree;
		SE[ilev][ise].NS = irenb;
		SE[ilev][ise].NI = ireni;
		SE[ilev][ise].NN = SE[ilev][ise].NI + SE[ilev][ise].NS;
		SE[ilev][ise].NSE = SE[ilev][ise].NS*pfm->KORT;
		SE[ilev][ise].NIE = SE[ilev][ise].NI*pfm->KORT;
		SE[ilev][ise].NNE = SE[ilev][ise].NN*pfm->KORT;

		SE[ilev][ise].is_matrix_calculated = FALSE;
		SE[ilev][ise].is_rightvect_calculated = FALSE;
		SE[ilev][ise].is_solvect_calculated = FALSE;

		//составление строки матрицы индексов для следующего уровня
		INDnextlev[ise] = new int[SE[ilev][ise].NS + 2];
		INDnextlev[ise][0] = 100; //тип - суперэлемент
		INDnextlev[ise][1] = SE[ilev][ise].NS;
		for (j=0; j<SE[ilev][ise].NS; j++)
		{
			INDnextlev[ise][j+2] = REN_B_LEV_INV[ REN_B[j] ];
		}

		//составление общего массива перенумерации узлов от локальных номеров СЭ к глобальным номерам уровня (ilev-1)
		for (j=0; j<irenb; j++)
		{
			REN[ireni] = REN_B[j]; 
			ireni++;
		}
		//экономичная очистка fl_n
		for (j=0; j<SE[ilev][ise].NN; j++)
		{
			fl_n[ REN[j] ] = false;
		}
		//перенумерация узлов, от нумерации уровня (ilev-1) к локальной нумерации внутри СЭ
		for (j=0; j<NNlev; j++)
		{
			REN_B[j] = -1;
		}
		for (j=0; j<SE[ilev][ise].NN; j++)
		{
			REN_B[ REN[j] ] = j;
		}

		//составление матрицы индексов
		INDse = new int*[SE[ilev][ise].NEL];
		for (i=0; i<SE[ilev][ise].NEL; i++)
		{
			iel = REE[i];
			in = IND[iel][1];
			INDse[i] = new int[in+2];
			INDse[i][0] = IND[iel][0];
			INDse[i][1] = in;

			for (j=0; j<in; j++)
			{
				INDse[i][j+2] = REN_B[ IND[iel][j+2] ]; //локальная нумерация узлов внутри суперэлемента
			}
		}
		//перевод вектора перенумерации узлов к глобальным номерам узлов исходной модели
		for(i=0; i<SE[ilev][ise].NN; i++)
		{
			REN[i] = RENlev[ REN[i] ];
		}

		//создание первичного массива координат узлов для осуществления выбора стартового узла при перенумерации
		double *CRDSE;
		CRDSE = new double [SE[ilev][ise].NNE];
		for (k = 0; k<SE[ilev][ise].NN; k++)
		{
			for (k1=0; k1<pfm->KORT; k1++)
			{
				CRDSE[ k*pfm->KORT + k1 ] = pfm->CRD[ REN[k]*pfm->KORT + k1 ];
			}
		}
		
		fprintf(fp_ase,"ilev = %d\tise = %d\tNN = %d\tNI = %d\tNS = %d\t",ilev,ise, SE[ilev][ise].NN,SE[ilev][ise].NI,SE[ilev][ise].NS);
		fflush(fp_ase);
		
		// перенумерация узлов Катхила-Макки
		int *INVP;
		int sb,sa,wb,wa,nbr;
		INVP = new int [SE[ilev][ise].NI];
		MM->GENRCM(SE[ilev][ise].NN,SE[ilev][ise].NI,SE[ilev][ise].NEL,pfm->KORT,CRDSE,INDse,INVP,
			&sb,&sa,&wb,&wa,&nbr);
		//перестановка матрицы индексов и массива координат делается внутри GENRCM
		
		fprintf(fp_ase,"RCM\tsize_before = %d\tsize_after = %d\twidth_before = %d\twidth_after = %d\n",sb,sa,wb,wa);
		fflush(fp_ase);

		// перестановка массива глобальной перенумерации узлов
		MM->RenumberVect(REN,INVP,SE[ilev][ise].NI,1);

		delete []INVP;

		// вывод матрицы индексов
		sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_ind);
		fp = fopen(str,"wb");
		for (k=0; k<SE[ilev][ise].NEL; k++)
		{
			fwrite(INDse[k],sizeof(int),INDse[k][1]+2,fp);
		}
		fclose(fp);

		for (k=0; k<SE[ilev][ise].NEL; k++)
		{
			delete [] INDse[k];
		}
		delete [] INDse;


		// вывод массива перенумерации узлов
		sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_rennodes);
		fp = fopen(str,"wb");
		fwrite(REN,sizeof(int),SE[ilev][ise].NN,fp);
		fclose(fp);

		//вывод массива перенумерации элементов
		if (nse != 1) //для последнего уровня не требуется
		{
			sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_renelem);
			fp = fopen(str,"wb");
			fwrite(REE,sizeof(int),SE[ilev][ise].NEL,fp);
			fclose(fp);
		}

		//составление вектора координат узлов, закреплений,  материалов, только для СЭ нулевого уровня
		//вектор нагрузок составляется и передается в процессе решения тоже только СЭ нулевого уровня
		if (ilev == 0)
		{
			// вывод координат узлов
			sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_crd);
			fp = fopen(str,"wb");
			fwrite(CRDSE,sizeof(double),SE[ilev][ise].NNE,fp);
			fclose(fp);

			int *MTRSE;
			MTRSE = new int[SE[ilev][ise].NEL];
			for (k=0; k<SE[ilev][ise].NEL; k++)
			{
				MTRSE[k] = pfm->MTR[REE[k]];
			}
			sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_mat);
			fp = fopen(str,"wb");
			fwrite(MTRSE,sizeof(int),SE[ilev][ise].NEL,fp);
			fclose(fp);
			delete []MTRSE;

			//создание файлов закреплений СЭ
			int *FIX;
			double *UFIX;
			FIX = new int [SE[ilev][ise].NNE];
			UFIX = new double [SE[ilev][ise].NNE];
			for (k = 0; k<SE[ilev][ise].NN; k++)
			{
				for (k1=0; k1<pfm->KORT; k1++)
				{
					FIX[ k*pfm->KORT + k1 ] = pfm->FIX[ REN[k]*pfm->KORT + k1 ];
					UFIX[ k*pfm->KORT + k1 ] = pfm->UFIX[ REN[k]*pfm->KORT + k1 ];
				}
			}
			sprintf(str,"%s\\%s",pathmatr,SE[ilev][ise].name_fix);
			fp = fopen(str,"wb");
			fwrite(FIX,sizeof(int),SE[ilev][ise].NNE,fp);
			fwrite(UFIX,sizeof(double),SE[ilev][ise].NNE,fp);
			fclose(fp);
			delete []FIX;
			delete []UFIX;
		}

		delete []CRDSE;
	}

	//переформирование RENlev для перехода на следующий уровень (теперь он будет соответствовать матрице индексов INDnextlev)
	for (i=0; i<nblev; i++)
	{
		REN_B_LEV[i] = RENlev[ REN_B_LEV[i] ];
	}
	for (i=0; i<nblev; i++)
	{
		RENlev[i] = REN_B_LEV[i];
	}
	*NNnextlev = nblev;

	//очистка памяти
	delete [] fl_n;
	delete [] REE;
	delete [] REN;
	delete [] REN_B;
	delete [] REN_B_LEV_INV;
	delete [] REN_B_LEV;
	delete [] fl_b;

}