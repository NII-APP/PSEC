#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"

void FULLMODEL::MassMatrix()
{
	int i,j,k,kn,ii,ieltype,icurse,iur;
	long long int tmpi;
	char path[256];

	//выделение памяти под матрицы
	Mdiag = MM->MEM_NEW(Mdiag,NNE);

	//вычисление ММ и ассемблирование
	for (i=0; i<NEL; i++)
	{
		ieltype = IND[i][0];

		if (el[ieltype].isInitialized == 0) printf("FULLMODEL  unknown element  iel= %d\n",i);
		AttachElement(&el[ieltype],i);

		el[ieltype].ELMASS();
		
		MassAssemblingDiag(&el[ieltype],i);
	}
}

void FULLMODEL::MassDiagMatrixFix()
{
	int i,j;

	for (i=0; i<NNE; i++){
		if (FIX[i] > 0){
			Mdiag[i] = 0.0;
		}
	}
}

void FULLMODEL::MassAssemblingDiag(CEL *el,int elnum)
{
	int i,j,iur,jur,tmp;
	
	for(i=0; i<el->NNE; i++)
	{
		iur = el->GLOBNE[i];

		Mdiag[iur] += el->MASS[i][i];
	}
}
