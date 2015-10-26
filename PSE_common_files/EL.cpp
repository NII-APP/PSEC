#include "StdAfx.h"
#include "EL.h"

CEL::CEL(void)
{
	isInitialized = 0;
	isVELcalculated = 0;
	isVPbelowzero = 0;

	ind = NULL;
	fullcrd = NULL;
	STIF = NULL;
	MASS = NULL;
	B = NULL;
	D = NULL;
	vlcrdint = NULL;
	Wvint = NULL;
	faceloccrdint = NULL;
	faceorientcrd_num = NULL;
	faceorientcrd_val = NULL;
	FaceNodes = NULL;
	FaceNodesList = NULL;
	Wfint = NULL;
	FFv = NULL;
	FFf = NULL;
	dFFvloc = NULL;
	dFFfloc = NULL;
	GLOBNE = NULL;
	P = NULL;
	material = NULL;

	MM = new MEM;
}

void CEL::Initialize(int ieltype)
{
	int i,j;

	ielobjnumber = ieltype;
	isInitialized = 0;
	eltype = -1;
	if (ieltype == 5) //линейный плоский треугольник
	{
		eltype = ieltype;
		GeneralParameters5();
		GeneralMemInit();
//		Init_VolIntPointLocCrd_22();
//		Init_FaceNodes_22();
		
	}
	if (ieltype == 22) //квадратичный плоский треугольник
	{
		eltype = ieltype;
		GeneralParameters22();
		GeneralMemInit();
//		Init_VolIntPointLocCrd_22();
//		Init_FaceNodes_22();
		
	}
	if (ieltype == 24) //квадратичный тетраэдр
	{
		eltype = ieltype;
		GeneralParameters24();
		GeneralMemInit();
		Init_VolIntPointLocCrd_24();
		Init_FaceNodes_24();
				
	}
	if (ieltype == 25) //квадратичный кубик
	{
		eltype = ieltype;
		GeneralParameters25();
		GeneralMemInit();
		Init_VolIntPointLocCrd_25();
		Init_FaceNodes_25();
		Init_FacePointLocCrd_25();

		//расчет производных функций формы по локальным координатам в точках интегрирования по граням
		for (i=0; i<Nface; i++)
		{
			for (j=0; j<Nfp; j++)
			{
				dFF_25(dFFfloc[i][j],faceloccrdint[i][j]);
			}
		}
	}

	if (eltype > 0 && eltype != 22)//!!!!!!!!!!временно, нужно доопределить 22 тип
	{
		//расчет значений функций формы в точках интегрирования по объему
		for (i=0; i<NPINT; i++)
		{
			FF_all(FFv[i],vlcrdint[i]);
		}
		//расчет производных функций формы по локальным координатам в точках интегрирования по объему
		for (i=0; i<NPINT; i++)
		{
			dFF_all(dFFvloc[i],vlcrdint[i]);
		}

		isInitialized = 1;
	}
}

void CEL::FF_all(double *F, double *lcrd)
{
	switch (eltype)
	{
	case 5:
		FF_5(F,lcrd);
		break;
	case 22:
		FF_22(F,lcrd);
		break;
	case 24:
		FF_24(F,lcrd);
		break;
	case 25:
		FF_25(F,lcrd);
	default:
		break;
	}
}
void CEL::dFF_all(double **dF, double *lcrd)
{
	switch (eltype)
	{
	case 5:
		dFF_5(dF,lcrd);
		break;
	case 22:
		dFF_22(dF,lcrd);
		break;
	case 24:
		dFF_24(dF,lcrd);
		break;
	case 25:
		dFF_25(dF,lcrd);
	default:
		break;
	}
}

void CEL::GeneralMemInit()
{
	// выделение памяти для матрицы градиентов
		B = MM->MEM_NEW(B,NDEF,NNE);

		// матрица свойств 
		D = MM->MEM_NEW(D,NDEF,NDEF);
				
		// матрица жесткости элемента
		STIF = MM->MEM_NEW(STIF,NNE,NNE);

		// матрица масс элемента
		MASS = MM->MEM_NEW(MASS,NNE,NNE);

		// локальные координаты точек интегрирования по объему
		vlcrdint = MM->MEM_NEW(vlcrdint,NPINT,nvlcrd);

		// весовые коэффициенты при интегрировании по объему
		Wvint = MM->MEM_NEW(Wvint,NPINT);

		// локальные координаты в точках интегрирования по граням
		faceloccrdint = MM->MEM_NEW(faceloccrdint,Nface,Nfp,nflcrd);

		// весовые коэффициенты в точках интегрирования по граням
		Wfint = MM->MEM_NEW(Wfint,Nface,Nfp);

		// значения функций формы в точках интегрирования по объему
		FFv = MM->MEM_NEW(FFv,NPINT,NN);

		// значения функций формы в точках интегрирования по граням
		FFf = MM->MEM_NEW(FFf,Nface,Nfp,NN);

		// значения производных функций формы по локальным координатам в точках интегрирования по объему
		dFFvloc = MM->MEM_NEW(dFFvloc,NPINT,NORT,NN);

		// значения производных функций формы по локальным координатам в точках интегрирования по граням
		dFFfloc = MM->MEM_NEW(dFFfloc,Nface,NPINT,NORT,NN);

		//номер лок. координаты, которая наа грани является постоянной
		faceorientcrd_num = MM->MEM_NEW(faceorientcrd_num,Nface);

		//значение локальной координаты, которая постоянна на грани
		faceorientcrd_val = MM->MEM_NEW(faceorientcrd_val,Nface);
		
		//массивы флагов принадлежности узлов в локальной нумерации граням
		FaceNodes = MM->MEM_NEW(FaceNodes,Nface,NN);

		//локальные номера узлов граней
		FaceNodesList = MM->MEM_NEW(FaceNodesList,Nface,NNface);

		//массив номеров степеней свободы
		GLOBNE = MM->MEM_NEW(GLOBNE,NNE);
}

CEL::~CEL(void)
{
	int j;
	eltype = -1;

	if (isInitialized == 1)
	{
	// выделение памяти для матрицы градиентов
		B = MM->MEM_DEL(B,NDEF,NNE);

		// матрица свойств 
		D = MM->MEM_DEL(D,NDEF,NDEF);
				
		// матрица жесткости элемента
		STIF = MM->MEM_DEL(STIF,NNE,NNE);

		// матрица масс элемента
		MASS = MM->MEM_DEL(MASS,NNE,NNE);

		// локальные координаты точек интегрирования по объему
		vlcrdint = MM->MEM_DEL(vlcrdint,NPINT,nvlcrd);

		// весовые коэффициенты при интегрировании по объему
		Wvint = MM->MEM_DEL(Wvint,NPINT);

		// локальные координаты в точках интегрирования по граням
		faceloccrdint = MM->MEM_DEL(faceloccrdint,Nface,Nfp,nflcrd);

		// весовые коэффициенты в точках интегрирования по граням
		Wfint = MM->MEM_DEL(Wfint,Nface,Nfp);

		// значения функций формы в точках интегрирования по объему
		FFv = MM->MEM_DEL(FFv,NPINT,NN);

		// значения функций формы в точках интегрирования по граням
		FFf = MM->MEM_DEL(FFf,Nface,Nfp,NN);

		// значения производных функций формы по локальным координатам в точках интегрирования по объему
		dFFvloc = MM->MEM_DEL(dFFvloc,NPINT,NORT,NN);

		// значения производных функций формы по локальным координатам в точках интегрирования по граням
		dFFfloc = MM->MEM_DEL(dFFfloc,Nface,NPINT,NORT,NN);

		//номер лок. координаты, которая наа грани является постоянной
		faceorientcrd_num = MM->MEM_DEL(faceorientcrd_num,Nface);

		//значение локальной координаты, которая постоянна на грани
		faceorientcrd_val = MM->MEM_DEL(faceorientcrd_val,Nface);
		
		//массивы флагов принадлежности узлов в локальной нумерации граням
		FaceNodes = MM->MEM_DEL(FaceNodes,Nface,NN);

		//локальные номера узлов граней
		FaceNodesList = MM->MEM_DEL(FaceNodesList,Nface,NNface);

		//массив номеров степеней свободы
		GLOBNE = MM->MEM_DEL(GLOBNE,NNE);
	}
	delete MM;
}

void CEL::InitGLOBNE()
{
	int i,j;
	for(i=0; i<NN; i++)
	{
		for (j=0; j<NDOF; j++)
		{
			GLOBNE[i*NDOF+j] = ind[i]*NDOF + j;
		}
	}
}


void CEL::shifting( double shift)
{
	int i;

	for (i=0; i<NNE; i++)
	{
		STIF[i][i] -= shift*MASS[i][i];
	}

}


void CEL::GetBoundBox(double *vmin, double *vmax)
{
	int i,j;
	//определение минимальных и максимальных координат ограничивающего элемент параллелепипеда.
	for (j=0; j<NORTfullcrd; j++)
	{
		vmin[j] = 10000000.0;
		vmax[j] = -10000000.0;
	}
	for (i=0; i<NN; i++)
	{
		for (j=0; j<NORTfullcrd; j++)
		{
			if (fullcrd[ind[i]*NORTfullcrd+j] < vmin[j]) vmin[j] = fullcrd[ind[i]*NORTfullcrd+j];
			if (fullcrd[ind[i]*NORTfullcrd+j] > vmax[j]) vmax[j] = fullcrd[ind[i]*NORTfullcrd+j];
		}
	}
}

void CEL::GetCentr(double *vcentr)
{
	int i,j;
	//определение центра элемента, как среднего значения координат его узлов
	for (j=0; j<NORTfullcrd; j++)
	{
		vcentr[j] = 0.0;
	}
	for (i=0; i<NN; i++)
	{
		for (j=0; j<NORTfullcrd; j++)
		{
			vcentr[j] += fullcrd[ind[i]*NORTfullcrd+j];
		}
	}
	for (j=0; j<NORTfullcrd; j++)
	{
		vcentr[j] = vcentr[j]/NN;
	}
}