#include "StdAfx.h"
#include "EL.h"

void CEL::GeneralParameters5()
{
	NN = 3;
	NORT = 2;
	NORTfullcrd = NORT;
	NDOF = 2;
	NPINT = 3; // O(h^3)
	NNE = NN*NDOF;
	NNORT = NN*NORT;
	NDEF = 3;
	NSTRS = 3;
	VEL = 0.0;

	lcrdlim[0] = 0.0;
	lcrdlim[1] = 1.0;

	nvlcrd = 3; //3 локальные координаты (L-координаты) одна из них зависимая
	Nface = 0; // сам же является поверхностью
	NNface = 0; //
	Nfp = 1; //!!!!!!!!!!!!!!!!! проверить!!!
	nflcrd = 4; //!!!!!!!!!!!!!!!!! проверить!!!
}


void CEL::Init_VolIntPointLocCrd_5()
{
	int i,j;
	double a,b;

	a = 0.0;
	b = 0.5;

	for (i=0; i<NPINT; i++)
	{
		for (j=0; j<nvlcrd; j++)
		{
			vlcrdint[i][j] = b;
		}
		vlcrdint[i][i] = a;
		Wvint[i] = 1.0/3.0;
	}
}

void CEL::FF_5(double *F, double *lcrd)
{
	lcrd[2] = 1.0 - lcrd[0] - lcrd[1];

	F[0] = lcrd[0];
	F[1] = lcrd[1];
	F[2] = lcrd[2];
}

void CEL::dFF_5(double **dF, double *lcrd)
{
	lcrd[2] = 1.0 - lcrd[0] - lcrd[1];

	dF[0][0] = 1;		dF[1][0] = 0;				
	dF[0][1] = 0;		dF[1][1] = 1;				
	dF[0][2] = -1;		dF[1][2] = -1;	
}