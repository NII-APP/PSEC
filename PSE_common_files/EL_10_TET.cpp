#include "StdAfx.h"
#include "EL.h"

void CEL::GeneralParameters24()
{
	NN = 10;
	NORT = 3;
	NDOF = 3;
	NPINT = 4;
	NNE = NN*NDOF;
	NNORT = NN*NORT;
	NDEF = 6;
	NSTRS = 6;
	VEL = 0.0;

	nvlcrd = 4;
	Nface = 4;
	NNface = 6;
	Nfp = 1; //!!!!!!!!!!!!!!!!! проверить!!!
	nflcrd = 4; //!!!!!!!!!!!!!!!!! проверить!!!
}


void CEL::Init_VolIntPointLocCrd_24()
{
	int i,j;
	double a,b;

	a = 0.58541020;
	b = 0.13819660;

	for (i=0; i<NPINT; i++)
	{
		for (j=0; j<nvlcrd; j++)
		{
			vlcrdint[i][j] = b;
		}
		vlcrdint[i][i] = a;
		Wvint[i] = 0.25;
	}
}

void CEL::FF_24(double *F, double *lcrd)
{
	F[0] = lcrd[0]*(2*lcrd[0]-1);
	F[1] = lcrd[1]*(2*lcrd[1]-1);
	F[2] = lcrd[2]*(2*lcrd[2]-1);
	F[3] = lcrd[3]*(2*lcrd[3]-1);
	F[4] = 4*lcrd[0]*lcrd[1];
	F[5] = 4*lcrd[1]*lcrd[2];
	F[6] = 4*lcrd[2]*lcrd[0];
	F[7] = 4*lcrd[0]*lcrd[3];
	F[8] = 4*lcrd[1]*lcrd[3];
	F[9] = 4*lcrd[2]*lcrd[3];
}

void CEL::dFF_24(double **dF, double *lcrd)
{
	dF[0][0] = 4*lcrd[0]-1;		dF[1][0] = 0;				dF[2][0] = 0;
	dF[0][1] = 0;				dF[1][1] = 4*lcrd[1]-1;		dF[2][1] = 0;
	dF[0][2] = 0;				dF[1][2] = 0;				dF[2][2] = 4*lcrd[2]-1;
	dF[0][3] = -4*lcrd[3]+1;	dF[1][3] = -4*lcrd[3]+1;	dF[2][3] = -4*lcrd[3]+1;
	
	dF[0][4] = 4*lcrd[1];				dF[1][4] = 4*lcrd[0];				dF[2][4] = 0;
	dF[0][5] = 0;						dF[1][5] = 4*lcrd[2];				dF[2][5] = 4*lcrd[1];
	dF[0][6] = 4*lcrd[2];				dF[1][6] = 0;						dF[2][6] = 4*lcrd[0];
	dF[0][7] = 4*(lcrd[3]-lcrd[0]);		dF[1][7] = -4*lcrd[0];				dF[2][7] = -4*lcrd[0];
	dF[0][8] = -4*lcrd[1];				dF[1][8] = 4*(lcrd[3]-lcrd[1]);		dF[2][8] = -4*lcrd[1];
	dF[0][9] = -4*lcrd[2];				dF[1][9] = -4*lcrd[2];				dF[2][9] = 4*(lcrd[3]-lcrd[2]);
}

void CEL::Init_FaceNodes_24()
{
	int i,j,iface;
	int *p;

	p = FaceNodesList[0];
	p[0] = 0;	p[1] = 1;	p[2] = 2;	p[3] = 4;	p[4] = 5;	p[5] = 6;
	p = FaceNodesList[1];
	p[0] = 0;	p[1] = 1;	p[2] = 3;	p[3] = 4;	p[4] = 8;	p[5] = 7;
	p = FaceNodesList[2];
	p[0] = 3;	p[1] = 1;	p[2] = 2;	p[3] = 8;	p[4] = 5;	p[5] = 9;
	p = FaceNodesList[3];
	p[0] = 0;	p[1] = 3;	p[2] = 2;	p[3] = 7;	p[4] = 9;	p[5] = 6;

	for (i=0; i<Nface; i++)
	{
		for (j=0; j<NN; j++)
		{
			FaceNodes[i][j] = 0;
		}
	}
	for (i=0; i<Nface; i++)
	{
		for (j=0; j<NNface; j++)
		{
			FaceNodes[i][FaceNodesList[i][j]] = 1;
		}
	}
}