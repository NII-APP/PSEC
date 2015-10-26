#include "StdAfx.h"
#include "EL.h"
#include "math.h"

void CEL::ELSTIF()
{
	int ipo,j;
	double **J,detJ;
	
	InitGLOBNE();

	J = NULL;
	J = MM->MEM_NEW(J,NORT,NORT);

	MM->QQZERO(STIF,NNE,NNE);
	VEL = 0.0;
	isVPbelowzero = 0;

	DINT();
	for (ipo=0; ipo<NPINT; ipo++)
	{
		Jacoby(dFFvloc[ipo],J);
		DetJacoby(J,&detJ);
		
		P[ipo].VP = detJ*Wvint[ipo];
		if ( eltype == 24) P[ipo].VP = P[ipo].VP/6;
		if ( P[ipo].VP < 0.0 ) P[ipo].VP = fabs(P[ipo].VP); 
		VEL += P[ipo].VP;

		InvJacoby(J,detJ);
		CalcBGrad_3D(J,dFFvloc[ipo]);
		CalcSTIFPoint(&P[ipo]);
	}

	isVELcalculated = 1;
	J = MM->MEM_DEL(J,NORT,NORT);
}


void CEL::ELMASS()
{
	int i,j;
	double mm;
	double **J,detJ;

	J = NULL;
	J = MM->MEM_NEW(J,NORT,NORT);

	InitGLOBNE();
	MM->QQZERO(MASS,NNE,NNE);
	
	if ( isVELcalculated == 0 )
	{
		VEL = 0;
		if ( P[0].VP > 0 )
		{
			for (i=0; i<NPINT; i++)
			{
				VEL += P[i].VP;
			}
		}
		else
		{
			for (i=0; i<NPINT; i++)
			{
				Jacoby(dFFvloc[i],J);
				DetJacoby(J,&detJ);
				P[i].VP = detJ*Wvint[i];
				if ( eltype == 24) P[i].VP = P[i].VP/6;
				if ( P[i].VP < 0.0 ) P[i].VP = fabs(P[i].VP); 
				VEL += P[i].VP;
			}
		}
		isVELcalculated = 1;
	}


	mm = VEL*material->RO;
	mm = mm/NN;

	for (i=0; i<NNE; i++)
	{
		MASS[i][i] = mm;
	}

	J = MM->MEM_DEL(J,NORT,NORT);
}

void CEL::CalcBGrad_3D(double **invJ, double **dFF)
{
	int i,j,k,z;
	double s;
	double **dFFglob;

	dFFglob = NULL;
	dFFglob = MM->MEM_NEW(dFFglob,NORT,NN);

	//вычисление производных функций формы по глобальным координатам для элемента IEL в точке интегрирования IPO

	for (i=0; i<NORT; i++)
	{
		for (j=0; j<NN; j++)
		{
			dFFglob[i][j] = 0;
			for (k=0; k<NORT; k++)
			{
				dFFglob[i][j] += invJ[i][k]*dFF[k][j];
			}
		}
	}

	//заполнение матрицы градиентов
	for(i=0; i<NN; i++)
	{
		B[0][i*NORT] = dFFglob[0][i];
		B[1][i*NORT+1] = dFFglob[1][i];
		B[2][i*NORT+2] = dFFglob[2][i];
		B[3][i*NORT] = dFFglob[1][i]; B[3][i*NORT+1] = dFFglob[0][i];
		B[4][i*NORT+1] = dFFglob[2][i]; B[4][i*NORT+2] = dFFglob[1][i];
		B[5][i*NORT] = dFFglob[2][i]; B[5][i*NORT+2] = dFFglob[0][i];
	}
	
	dFFglob = MM->MEM_DEL(dFFglob,NORT,NN);
}

void CEL::CalcSTIFPoint(INTPOINT *P)
{
	int i,j,z;
	double s;
	double **BD;
	
	BD = NULL;
	BD = MM->MEM_NEW(BD,NNE,NDEF);
	
	//B транспонированное умножаем на D
	for(i=0; i<NNE; i++)
	{
		for(j=0; j<NDEF; j++)
		{
			for(z=0; z<NDEF; z++)
			{
				BD[i][j] += B[z][i]*D[z][j];
			}
		}
	}
	//BD умножаем на B и на объем точки интегрирования. Результат помещаем в GE
	for(i=0; i<NNE; i++)
	{
		for(j=0; j<NNE; j++)
		{
			s = 0.0;
			for(z=0; z<NDEF; z++)
			{
				s+=BD[i][z]*B[z][j];
			}
			STIF[i][j] += s*P->VP;
		}
	}

	BD = MM->MEM_DEL(BD,NNE,NDEF);
}






