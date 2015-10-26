#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"

FULLMODEL::FULLMODEL()
{
	sprintf(pathmain,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	loadtype = 1;

	MAXELTYPES = 30;
	MM = NULL;
	
	mat = NULL;

	IND = NULL;
	CRD = NULL;
	MTR = NULL;
	FIX = NULL;
	UFIX = NULL;
	Mdiag = NULL;
	P = NULL;
	envP = NULL;
	LV = NULL;
	RV = NULL;
	fl_LVRV_internal = 0;
	pforce = NULL;
	pdisp = NULL;
	npforce = 0;
	npdisp = 0;

	NDEL = NULL;
	NDGR = NULL;
	ELGR = NULL;

	LIST_EL = NULL;
	MASK_EL = NULL;
	nlistel = 0;

	surfmodel = NULL;

	//подготовка библиотеки конечных элементов
	int ieltype;
	el = new CEL[MAXELTYPES];
	for (ieltype = 0; ieltype<MAXELTYPES; ieltype++)
	{
		el[ieltype].Initialize(ieltype);
	}
}

FULLMODEL::~FULLMODEL()
{
	if ( mat != NULL ) { delete []mat; mat = NULL; }

	IND = MM->MEM_DEL(IND,NEL,1); //второй номер фиктивный
	CRD = MM->MEM_DEL(CRD,NNE); 
	MTR = MM->MEM_DEL(MTR,NEL);
	FIX = MM->MEM_DEL(FIX,NNE);
	UFIX = MM->MEM_DEL(UFIX,NNE);
	Mdiag = MM->MEM_DEL(Mdiag,NNE);
	if ( P != NULL )
	{
		delete []P;
		P = NULL;
	}
	if ( envP != NULL)
	{
		envP = MM->MEM_DEL(envP,NEL);
		envP = NULL;
	}
	if (fl_LVRV_internal == 1)
	{
		LV = MM->MEM_DEL(LV,nvect,NNE);
		RV = MM->MEM_DEL(RV,nvect,NNE);
	}
	if ( pforce != NULL ) { delete []pforce; pforce = NULL; }
	if ( pdisp != NULL ) { delete []pdisp; pdisp = NULL; }
	npforce = 0;
	npdisp = 0;

	NDEL = MM->MEM_DEL(NDEL,NN,1);
	NDGR = MM->MEM_DEL(NDGR,NN,1);
	ELGR = MM->MEM_DEL(ELGR,NEL,1);

	DelLIST_EL();

	if ( surfmodel != NULL ) { delete []surfmodel; surfmodel = NULL; }
	delete [] el;

}

void FULLMODEL::InitLIST_EL()
{
	LIST_EL = MM->MEM_NEW(LIST_EL,NEL);
	MASK_EL = MM->MEM_NEW(MASK_EL,NEL);
	nlistel = 0;
}

void FULLMODEL::ClearLIST_EL()
{
	int i;
	for (i=0; i<nlistel; i++)
	{
		MASK_EL[ LIST_EL[i] ] = 0;
	}
	nlistel = 0;
}

void FULLMODEL::DelLIST_EL()
{
	LIST_EL = MM->MEM_DEL(LIST_EL,NEL);
	MASK_EL = MM->MEM_DEL(MASK_EL,NEL);
	nlistel = 0;
}

void FULLMODEL::CalcActualForcePointPosition()
{
	int i,j;

	for(i=0; i<npforce; i++)
	{
		pforce[i].TrackingPoint((void*)this,&pforce[i].crd[0]);
	}
}
void FULLMODEL::CalcActualDispPointPosition()
{
	int i,j;
	for(i=0; i<npdisp; i++)
	{
		pdisp[i].TrackingExternPoint((void*)this,&pdisp[i].crd[0]);
	}
}

void FULLMODEL::ReadFullmodel()
{	
	if (loadtype == 1)
	{
		sprintf(pathmatr,"%s\\ds",pathmain);
		ReadFromUzor();
	}
	if (loadtype == 2)
	{
		sprintf(pathmatr,"%s\\ds",pathmain);
		ReadFromUzor_NewFormat();
	}
	if (loadtype == 3)
	{
		sprintf(pathmatr,"%s",pathmain);
		ReadFromAnsys_FullGrid();
	}
}




void FULLMODEL::AttachElement(CEL *el, int elnum)
{
	int i,j;
//	el->current_el_number = elnum;
	el->material = &mat[MTR[elnum]];
	el->P = &P[envP[elnum]];
	
	el->fullcrd = CRD;
	el->ind = &IND[elnum][2]; //ссылка на первый номер узла для КЭ
	el->isVELcalculated = 0;
}

void FULLMODEL::AttachElement(CEL *el, FACE *face) 
{
	el->fullcrd = CRD;
	el->NORTfullcrd = KORT;
	el->ind = face->facenodenums; //ссылка на первый номер узла для КЭ
	el->isVELcalculated = 0;
}


void FULLMODEL::InitIntPoint()
{
	int i,j,ieltype;

	if ( envP == NULL && P == NULL ) //выполняется однократно для модели
	{
		//расчет профиля хранения точек интегрирования, подсчет их общего количества
		nIntP = 0;
		envP = MM->MEM_NEW(envP,NEL);
		envP[0] = 0;
		for (i=0; i<(NEL-1); i++)
		{
			envP[i+1] = envP[i] + el[ IND[i][0] ].NPINT;
		}
		nIntP = envP[NEL-1] + el[ IND[NEL-1][0] ].NPINT;

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