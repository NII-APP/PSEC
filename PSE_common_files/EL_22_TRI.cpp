#include "StdAfx.h"
#include "EL.h"

void CEL::GeneralParameters22()
{
	NN = 6;
	NORT = 2;
	NORTfullcrd = NORT;
	NDOF = 2;
	NPINT = 4; //??????
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


void CEL::Init_VolIntPointLocCrd_22()
{

}

void CEL::FF_22(double *F, double *lcrd)
{
	lcrd[2] = 1.0 - lcrd[0] - lcrd[1];

	F[0] = lcrd[0]*(2*lcrd[0]-1);
	F[1] = lcrd[1]*(2*lcrd[1]-1);
	F[2] = lcrd[2]*(2*lcrd[2]-1);
	F[3] = 4*lcrd[0]*lcrd[1];
	F[4] = 4*lcrd[1]*lcrd[2];
	F[5] = 4*lcrd[2]*lcrd[0];
}

void CEL::dFF_22(double **dF, double *lcrd)
{
	lcrd[2] = 1.0 - lcrd[0] - lcrd[1];

	dF[0][0] = 4*lcrd[0]-1;		dF[1][0] = 0;				
	dF[0][1] = 0;				dF[1][1] = 4*lcrd[1]-1;				
	dF[0][2] = -4*lcrd[2]+1;	dF[1][2] = -4*lcrd[2]+1;	
	
	dF[0][3] = 4*lcrd[1];				dF[1][3] = 4*lcrd[0];				
	dF[0][4] = -4*lcrd[1];				dF[1][4] = 4*(lcrd[2]-lcrd[1]);	
	dF[0][5] = 4*(lcrd[2]-lcrd[0]);		dF[1][5] = -4*lcrd[0];				
}