#include "StdAfx.h"
#include "SE.h"

//CLIENTDATA cldat;

CSE::CSE(void)
{
	bignumber = 1.e80;
	MM = new MEM;
	KORT = 3;

	nIntP = 0;

	sedat = NULL;
	pmat = NULL;

	P = NULL;
	envP = NULL; //массив показывающий стартовый номер позиции в массиве точек интегрирования, соответствующий каждому конечному элементу

	//мтарицы
	STII = NULL;
	STIS = NULL;
	STISTR = NULL;
	STSS = NULL;
	STII_ENV = NULL;
	STIS_ENV = NULL;

	STIINZ = NULL;
	STIIPC = NULL;
	STII_ENVNZ= NULL;
	STII_ENVPC =NULL;

	KK = NULL;

	//исходные данные
	IND = NULL;
	CRD = NULL;
	MTR = NULL;
	REN = NULL; //глобальные номера узлов
	REE = NULL; //глобальные номера элементов
	FIX = NULL;
	UFIX= NULL;

	FORSE = NULL;
	DISP = NULL;

	MAXELTYPES = 30;

	ThreadQueue = NULL;
	EvNotEmptyQueue = NULL;

	nvect = -1;

	//параметры параллельности при расчете матрицы жесткости в подпространстве собственных форм
	nthread = 4;
	Queuesize = nthread + 2; //должен быть хотя бы на 1 больше nthread
	blocksize = 200;
}

CSE::~CSE(void)
{
	FULLDEL();

	ThreadQueue = NULL;
	EvNotEmptyQueue = NULL;

	delete MM;

}


void CSE::Initialize(SEstruct *ses)
{
	sedat = ses;
	NS = ses->NS;
	NI = ses->NI;
	NN = ses->NN;
	NEL = ses->NEL;

}

void CSE::MainMatrCalc()
{
	printf("\n\nMainMatrCalc  iLevel = %d   iSE = %d....\n",sedat->SElevel,sedat->iSE);

	ReadSEinitials();
	InitIntPoint();
	STIIenv();
	if ( sedat->fl_PCG_internal == true )
	{
		STIINZenv (0);
	}
	if ( sedat->cur_operation_type == 2 && sedat->SElevel == 0)
	{//для расчета собственных частот
		STIINZenv (1);
	}
	if (NS != 0)
	{
		STISenv();
	}
	printf("Assembling....\n");
	StiffMatrix();
	IntPWrite();
	if (NS != 0)
	{
		StisWrite();
	}
	printf("Transformation....\n");
	StiffMatrixTransform();
	
	if ( sedat->cur_operation_type == 2 && sedat->SElevel == 0)
	{//для расчета собственных частот
		WriteNZfull();
	}
	
	if ( sedat->fl_PCG_internal == true )
	{
		StifPrecondition();
		StiiWriteNZPC();
	}
	else
	{
		StiiWrite();
	}
	if (NS != 0)
	{
		StssWrite();
	}
	sedat->is_matrix_calculated = true;
	printf("Finished....\n");
}

void CSE::MainLoadCalc()
{
	printf("\n\n MainLoadCalc  iLevel = %d   iSE = %d....\n",sedat->SElevel,sedat->iSE);
	ReadSEinitials();
	if ( sedat->fl_PCG_internal == true )
	{
		StiiReadNZPC();
	}
	else
	{
		StiiRead();
	}

	if ( NS != 0 )
	{
		StisRead();
	}
	printf("vector loading....\n");
	LoadForse();
	ForseTransform();
	if ( NS == 0 )
	{
		//для верхнего уровня сразу осуществляется решение СЛАУ, запись векторов результата
		DispWriteSS(); //вместе с удалением памяти
	}
	else
	{
		ForseWriteSS(); //вместе с удалением памяти
	}
	
}

void CSE::MainResCalc()
{
	printf("\n\n MainResCalc  iLevel = %d   iSE = %d....\n",sedat->SElevel,sedat->iSE);
	ReadSEinitials();
	if ( sedat->fl_PCG_internal == true )
	{
		StiiReadNZPC();
	}
	else
	{
		StiiRead();
	}
	StisRead();
	printf("vector loading....\n");
	LoadDisp();
	DispTransform();
	if ( sedat->SElevel == 0 )
	{
		//для нижнего уровня  запись векторов результата осуществляется в упрощенной форме
		DispWrite(); //вместе с удалением памяти
	}
	else
	{
		DispWriteSS(); //вместе с удалением памяти
	}
	
}


void CSE::MainSubspaceCalc()
{
	printf("\n\n MainSubspaceCalc  iLevel = %d   iSE = %d....\n",sedat->SElevel,sedat->iSE);
	ReadSEinitials();
	ReadNZfull();
	
	printf("vector loading....\n");
	LoadForse(); //через вектора сил передаются вектора собственных форм.
	ForseFixing();
	
	//CalcSubspaceKM();
	CalcSubspaceKMPar();

	WriteSubspaceKM();	
}



void CSE::InitIntPoint()
{
	int i,j,ieltype;

	if ( sedat->SElevel == 0 )
	{
		//подготовка библиотеки конечных элементов
		CEL *el;
		el = new CEL[MAXELTYPES];
		for (ieltype = 0; ieltype<MAXELTYPES; ieltype++)
		{
			el[ieltype].Initialize(ieltype);
		}

		//расчет профиля хранения точек интегрирования, подсчет их общего количества
		nIntP = 0;
		envP = MM->MEM_NEW(envP,NEL);
		envP[0] = 0;
		for (i=0; i<(NEL-1); i++)
		{
			envP[i+1] = envP[i] + el[ IND[i][0] ].NPINT;
		}
		nIntP = envP[NEL-1] + el[ IND[NEL-1][0] ].NPINT;

		delete [] el;

		//создание массива точек интегрирования

		P = new INTPOINT[nIntP];
		for (i=0; i<nIntP; i++)
		{
			P[i].VP =0.0;
			MM->QQZERO(P[i].E,6);
			MM->QQZERO(P[i].S,6);
		}
	}

}

void CSE::StiffMatrix()
{
	int i,j,k,kn,ii,ieltype,icurse,iur;
	long long int tmpi;
	char path[256];

	icurse = 0;

	//подготовка библиотеки конечных элементов
	CEL *el;
	el = new CEL[MAXELTYPES];
	for (ieltype = 0; ieltype<MAXELTYPES; ieltype++)
	{
		el[ieltype].Initialize(ieltype);
	}

	//выделение памяти под матрицы
	STII = MM->MEM_NEW(STII,STII_ENV[NIE]);
	if ( sedat->fl_PCG_internal == true )
	{
		STIINZ = new double*[NIE];
		for (k=0; k<NIE; k++)
		{
			STIINZ[k] = NULL;
			STIINZ[k] = MM->MEM_NEW(STIINZ[k],STII_ENVNZ[k][0]);
		}
	}
	if ( sedat->cur_operation_type == 2 && sedat->SElevel == 0)
	{//для расчета собственных частот
		STIINZ = new double*[NNE];
		for (k=0; k<NNE; k++)
		{
			STIINZ[k] = NULL;
			STIINZ[k] = MM->MEM_NEW(STIINZ[k],STII_ENVNZ[k][0]);
		}
	}
	if (NS != 0)
	{
		STIS = new double *[NSE];
		for (k=0; k<NSE; k++)
		{
			kn = (int)(k/KORT);
			STIS[k] = NULL;
			STIS[k] = MM->MEM_NEW(STIS[k],STIS_ENV[kn][0]*KORT);
		}
		SSlength = ((long long int)((NSE*NSE - NSE + 0.1)/2.0)) + NSE;
		STSS = MM->MEM_NEW(STSS,SSlength);
	}


	//вычисление МЖ и ассемблирование
	for (i=0; i<NEL; i++)
	{
		ieltype = IND[i][0];
		if (ieltype >= 0 && ieltype < 100) //обычный конечный элемент
		{
			if (el[ieltype].isInitialized == 0) printf("unknown element iSe = %d level = %d iel= %d\n",sedat->iSE,sedat->SElevel,i);
			AttachElement(&el[ieltype],i);

			el[ieltype].ELSTIF();

			if ( sedat->cur_operation_type == 2 && sedat->SElevel == 0)
			{//для расчета собственных часто
				el[ieltype].ELMASS();
				el[ieltype].shifting(sedat->shift);
			}
			
			StifIIAssembling(&el[ieltype],i);
			if ( sedat->fl_PCG_internal == true )
			{
				StifIINZAssembling(&el[ieltype],i,0);
			}
			if ( sedat->cur_operation_type == 2 && sedat->SElevel == 0)
			{//для расчета собственных частот
				StifIINZAssembling(&el[ieltype],i,1);
			}
			
			if (NS != 0)
			{
				StifISAssembling(&el[ieltype],i);
				StifSSAssembling(&el[ieltype],i);
			}
		}
		else //суперэлемент
		{
			//считывание матрицы SS суперэлемента
			FILE *fp;
			sprintf(path,"%s\\stss%s.bin",sedat->PathInclSEs[icurse],sedat->namenumInclSEs[icurse]);
			fp  = fopen(path,"rb");

			tmpi = -1;
			fread(&tmpi,sizeof(long long int),1,fp);  //количество элементов в полностью заполненном профиле

			double *STSSse = NULL;
			STSSse = MM->MEM_NEW(STSSse,tmpi);
			fread(STSSse,sizeof(double),tmpi,fp);

			fclose(fp);
			
			// формирование массива NUR
			int KU, *NUR = NULL;
			KU = IND[i][1]*KORT;
			NUR = MM->MEM_NEW(NUR,KU);

			for (j=0; j<IND[i][1]; j++)
			{
				for (k=0; k<KORT; k++)
				{
					iur = j*KORT + k;
					NUR[iur] = IND[i][2+j]*KORT + k;
				}
			}

			// составление матриц нового суперэлемента

			StifIIAssemblingFromSE(STSSse, KU, NUR );
			if ( sedat->fl_PCG_internal == true )
			{
				StifIINZAssemblingSE(STSSse, KU, NUR);
			}
			if (NS != 0)
			{
				StifISAssemblingFromSE(STSSse, KU, NUR );
				StifSSAssemblingFromSE(STSSse, KU, NUR );
			}

			MM->MEM_DEL(STSSse,tmpi);
			MM->MEM_DEL(NUR,KU);

			icurse++;
		}
	}

	delete [] el;

	sedat->fl_matrix_recalc = true;

//	NEW = 1;


}


void CSE::AttachElement(CEL *el, int elnum)
{
	int i,j;
//	el->current_el_number = elnum;
	el->material = &pmat[MTR[elnum]];
	el->P = &P[envP[elnum]];
	
	el->fullcrd = CRD;
	el->ind = &IND[elnum][2]; //ссылка на первый номер узла для КЭ

	el->isVELcalculated = 0;
}


void CSE::DeleteLevel_1()
{//удаление матриц СЭ (наибольшая часть объема)
	int i,j,k;

	if ( STII != NULL )
	{
		delete []STII;
		STII = NULL;
		delete []STII_ENV;
		STII_ENV = NULL;
	}
	if ( STIINZ != NULL )
	{
		if (sedat->cur_operation_type == 2)
		{
			for (i = NNE-1; i>=0; i--)
			{
				delete []STIINZ[i];
			}
		}
		else
		{
			for (i = NIE-1; i>=0; i--)
			{
				delete []STIINZ[i];
			}
		}
		delete []STIINZ;
		STIINZ = NULL;
	}
	if ( STIIPC != NULL )
	{
		for (i = NIE-1; i>=0; i--)
		{
			delete []STIIPC[i];
		}
		delete []STIIPC;
		STIIPC = NULL;
	}
	if ( STII_ENVNZ != NULL )
	{
		if (sedat->cur_operation_type == 2)
		{
			for (i = NNE-1; i>=0; i--)
			{
				delete []STII_ENVNZ[i];
			}
		}
		else
		{
			for (i = NIE-1; i>=0; i--)
			{
				delete []STII_ENVNZ[i];
			}
		}
		delete []STII_ENVNZ;
		STII_ENVNZ = NULL;
	}
	if ( STII_ENVPC != NULL )
	{
		for (i = NIE-1; i>=0; i--)
		{
			delete []STII_ENVPC[i];
		}
		delete []STII_ENVPC;
		STII_ENVPC = NULL;
	}

	if ( STIS != NULL )
	{
		for (k=0; k<NSE; k++)
		{
			delete []STIS[k];
		}
		delete []STIS;
		STIS = NULL;

		for (k=0; k<NS; k++)
		{
			delete []STIS_ENV[k];
		}
		delete []STIS_ENV;
		STIS_ENV = NULL;
	}

	if ( STSS != NULL )
	{
		delete []STSS;
		STSS = NULL;
	}
}

void CSE::DeleteLevel_2()
{
	//удаление векторов правых частей и решений, матрицы в подпространстве для поиска собственных частот

	int i,j,k;

	if (nvect > 0)
	{
		if ( FORSE != NULL )
		{			
			for (k=0; k<nvect; k++)
			{
				delete []FORSE[k];
			}
			delete []FORSE;
			FORSE = NULL;
		}
		if ( DISP != NULL )
		{
			for (k=0; k<nvect; k++)
			{
				delete []DISP[k];
			}
			delete []DISP;
			DISP = NULL;
		}
		if ( KK != NULL )
		{
			for (k=0; k<nvect; k++)
			{
				delete []KK[k];
			}
			delete []KK;
			KK = NULL;
		}
	}
}

void CSE::DeleteLevel_3()
{
	//удаление данных по точкам интегрирования
	if ( P != NULL )
	{
		delete []P;
		P = NULL;
		nIntP = 0;
	}
	if ( envP != NULL )
	{
		delete []envP;
		envP = NULL;
	}
	
}

void CSE::DeleteLevel_4()
{
	//удаление базовой информации
	int i,j,k;

	if (IND != NULL)
	{
		for (i=0; i<NEL; i++)
		{
			delete []IND[i];
		}
		delete []IND;
		IND = NULL;
	}

	if (CRD != NULL)
	{
		delete []CRD;
		CRD = NULL;
	}

	if (MTR != NULL)
	{
		delete []MTR;
		MTR = NULL;
	}

	if (REN != NULL)
	{
		delete []REN;
		REN = NULL;
	}

	if (REE != NULL)
	{
		delete []REE;
		REE = NULL;
	}

	if (FIX != NULL)
	{
		delete []FIX;
		FIX = NULL;
	}

	if (UFIX != NULL)
	{
		delete []UFIX;
		UFIX = NULL;
	}
}

void CSE::FULLDEL()
{
	DeleteLevel_1();
	DeleteLevel_2();
	DeleteLevel_3();
	DeleteLevel_4();
}
