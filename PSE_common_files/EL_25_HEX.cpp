#include "StdAfx.h"
#include "EL.h"

void CEL::GeneralParameters25()
{
	NN = 20;
	NORT = 3;
	NORTfullcrd = NORT;
	NDOF = 3;
	NPINT = 8; //!!!!!!!!!!!!переделать на схему 27 точек
	NNE = NN*NDOF;
	NNORT = NN*NORT;
	NDEF = 6;
	NSTRS = 6;
	VEL = 0.0;

	lcrdlim[0] = -1.0;
	lcrdlim[1] = 1.0;

	nvlcrd = 3;
	Nface = 6;
	NNface = 8; //число узлов на грани
	Nfp = 4; //!!!!!!!!!!!!!!!!! проверить!!!
	nflcrd = 3; //!!!!!!!!!!!!!!!!! проверить!!!
}

void CEL::Init_FaceNodes_25()
{
	int i,j,iface;
	int *p;

	p = FaceNodesList[0];
	p[0] = 0;	p[1] = 1;	p[2] = 4;	p[3] = 5;	p[4] = 8;	p[5] = 12;	p[6] = 16;	p[7] = 17; //лок.крд[0] = +1  (кси)
	p = FaceNodesList[1];
	p[0] = 2;	p[1] = 3;	p[2] = 6;	p[3] = 7;	p[4] = 10;	p[5] = 14;	p[6] = 18;	p[7] = 19; //лок.крд[0] = -1
	p = FaceNodesList[2];
	p[0] = 1;	p[1] = 2;	p[2] = 5;	p[3] = 6;	p[4] = 9;	p[5] = 13;	p[6] = 17;	p[7] = 18; //лок.крд[1] = +1  (этта)
	p = FaceNodesList[3];
	p[0] = 0;	p[1] = 3;	p[2] = 4;	p[3] = 7;	p[4] = 11;	p[5] = 15;	p[6] = 16;	p[7] = 19; //лок.крд[1] = -1
	p = FaceNodesList[4];
	p[0] = 4;	p[1] = 5;	p[2] = 6;	p[3] = 6;	p[4] = 12;	p[5] = 13;	p[6] = 14;	p[7] = 15; //лок.крд[2] = +1 (дзета)
	p = FaceNodesList[5];
	p[0] = 0;	p[1] = 1;	p[2] = 2;	p[3] = 3;	p[4] = 8;	p[5] = 9;	p[6] = 10;	p[7] = 11; //лок.крд[2] = -1

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


int CEL::IdentifyFace(int *facenodes, int nfn)
{
	//идентификация осуществляется только по всем узлам грани
	//в массиве facenodes сначала должны следовать угловые узлы, затем, узлы на серединах сторон

	int i,j;
	int *nodesflag;

	nodesflag = MM->MEM_NEW(nodesflag,NN);

	//составление массива флагов узлов грани в локальной нумерации элемента
	for (i=0; i<nfn; i++)
	{
		for (j=0; j<NN; j++)
		{
			if ( facenodes[i] == ind[j] )
			{
				nodesflag[j] = 1;
			}
		}
	}

	//сопоставление с гранями элемента по произведению массива флагов
	int res;
	for (i=0; i<Nface; i++)
	{
		res = 0;
		for (j=0; j<NN; j++)
		{
			res += nodesflag[j]*FaceNodes[i][j];
		}
		if ( res == nfn ) break;
	}

	if ( i == Nface )
	{
		return(-1); //ошибка, не удалось найти грань
	}
	return(i);

	nodesflag = MM->MEM_DEL(nodesflag,NN);
}

void CEL::Init_FacePointLocCrd_25()
{
	//у элемента типа 25 - 6 граней. На каждой грани 4 точки (схема 2х2)
	//последовательность записи граней: лок.крд1 +,-;лок.крд2 +,-;лок.крд3 +,-;
	int iface;
	double a,b,c;
	int ipo;

	a = 0.577350269189626;

	iface = 0;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = 1.0; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = a;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = 1.0; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = a;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = 1.0; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = -a;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = 1.0; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = -a;

	faceorientcrd_num[iface] = 0;
	faceorientcrd_val[iface] = 1.0;

	iface = 1;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -1.0; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = a;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -1.0; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = a;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -1.0; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = -a;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -1.0; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = -a;

	faceorientcrd_num[iface] = 0;
	faceorientcrd_val[iface] = -1.0;

	iface = 2;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = 1.0; faceloccrdint[iface][ipo][2] = a;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = 1.0; faceloccrdint[iface][ipo][2] = a;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = 1.0; faceloccrdint[iface][ipo][2] = -a;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = 1.0; faceloccrdint[iface][ipo][2] = -a;

	faceorientcrd_num[iface] = 1;
	faceorientcrd_val[iface] = 1.0;

	iface = 3;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = -1.0; faceloccrdint[iface][ipo][2] = a;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = -1.0; faceloccrdint[iface][ipo][2] = a;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = -1.0; faceloccrdint[iface][ipo][2] = -a;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = -1.0; faceloccrdint[iface][ipo][2] = -a;

	faceorientcrd_num[iface] = 1;
	faceorientcrd_val[iface] = -1.0;

	iface = 4;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = 1.0;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = 1.0;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = 1.0;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = 1.0;

	faceorientcrd_num[iface] = 2;
	faceorientcrd_val[iface] = 1.0;

	iface = 5;
	ipo = 0; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = -1.0;
	ipo = 1; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = a; faceloccrdint[iface][ipo][2] = -1.0;
	ipo = 2; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = -a; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = -1.0;
	ipo = 3; Wfint[iface][ipo] = 1.0;
	faceloccrdint[iface][ipo][0] = a; faceloccrdint[iface][ipo][1] = -a; faceloccrdint[iface][ipo][2] = -1.0;

	faceorientcrd_num[iface] = 2;
	faceorientcrd_val[iface] = -1.0;

}

void CEL::Init_VolIntPointLocCrd_25()
{
	//!!!!!!!!!!!!!требуется заменить на схему 3х3х3
	double a,b,c;
	int ipo;

	a = 0.577350269189626;
	ipo = 0;
	vlcrdint[ipo][0] = a; vlcrdint[ipo][1] = -a; vlcrdint[ipo][2] = -a; Wvint[ipo] = 1.0;
	ipo = 1;
	vlcrdint[ipo][0] = a; vlcrdint[ipo][1] = a; vlcrdint[ipo][2] = -a; Wvint[ipo] = 1.0;
	ipo = 2;
	vlcrdint[ipo][0] = -a; vlcrdint[ipo][1] = a; vlcrdint[ipo][2] = -a; Wvint[ipo] = 1.0;
	ipo = 3;
	vlcrdint[ipo][0] = -a; vlcrdint[ipo][1] = -a; vlcrdint[ipo][2] = -a; Wvint[ipo] = 1.0;
	ipo = 4;
	vlcrdint[ipo][0] = a; vlcrdint[ipo][1] = -a; vlcrdint[ipo][2] = a; Wvint[ipo] = 1.0;
	ipo = 5;
	vlcrdint[ipo][0] = a; vlcrdint[ipo][1] = a; vlcrdint[ipo][2] = a; Wvint[ipo] = 1.0;
	ipo = 6;
	vlcrdint[ipo][0] = -a; vlcrdint[ipo][1] = a; vlcrdint[ipo][2] = a; Wvint[ipo] = 1.0;
	ipo = 7;
	vlcrdint[ipo][0] = -a; vlcrdint[ipo][1] = -a; vlcrdint[ipo][2] = a; Wvint[ipo] = 1.0;
}

void CEL::FF_25(double *F, double *lcrd)
{
	F[0] = (1.0 + lcrd[0])*(1.0 - lcrd[1])*(1.0 - lcrd[2])*(lcrd[0] - lcrd[1] - lcrd[2] - 2.0)/8.0;
	F[1] = (1.0 + lcrd[0])*(1.0 + lcrd[1])*(1.0 - lcrd[2])*(lcrd[0] + lcrd[1] - lcrd[2] - 2.0)/8.0;
	F[2] = (1.0 - lcrd[0])*(1.0 + lcrd[1])*(1.0 - lcrd[2])*(- lcrd[0] + lcrd[1] - lcrd[2] - 2.0)/8.0;
	F[3] = (1.0 - lcrd[0])*(1.0 - lcrd[1])*(1.0 - lcrd[2])*(- lcrd[0] - lcrd[1] - lcrd[2] - 2.0)/8.0;
	F[4] = (1.0 + lcrd[0])*(1.0 - lcrd[1])*(1.0 + lcrd[2])*(lcrd[0] - lcrd[1] + lcrd[2] - 2.0)/8.0;
	F[5] = (1.0 + lcrd[0])*(1.0 + lcrd[1])*(1.0 + lcrd[2])*(lcrd[0] + lcrd[1] + lcrd[2] - 2.0)/8.0;
	F[6] = (1.0 - lcrd[0])*(1.0 + lcrd[1])*(1.0 + lcrd[2])*(- lcrd[0] + lcrd[1] + lcrd[2] - 2.0)/8.0;
	F[7] = (1.0 - lcrd[0])*(1.0 - lcrd[1])*(1.0 + lcrd[2])*(- lcrd[0] - lcrd[1] + lcrd[2] - 2.0)/8.0;

	// узлы, для которых нулевая координата равна 0
	F[9] = (1.0 - lcrd[0]*lcrd[0])*(1.0 + lcrd[1])*(1.0 - lcrd[2])/4.0;
	F[11] = (1.0 - lcrd[0]*lcrd[0])*(1.0 - lcrd[1])*(1.0 - lcrd[2])/4.0;
	F[13] = (1.0 - lcrd[0]*lcrd[0])*(1.0 + lcrd[1])*(1.0 + lcrd[2])/4.0;
	F[15] = (1.0 - lcrd[0]*lcrd[0])*(1.0 - lcrd[1])*(1.0 + lcrd[2])/4.0;

	// узлы, для которых первая координата равна 0
	F[8] = (1.0 - lcrd[1]*lcrd[1])*(1.0 + lcrd[0])*(1.0 - lcrd[2])/4.0;
	F[10] = (1.0 - lcrd[1]*lcrd[1])*(1.0 - lcrd[0])*(1.0 - lcrd[2])/4.0;
	F[12] = (1.0 - lcrd[1]*lcrd[1])*(1.0 + lcrd[0])*(1.0 + lcrd[2])/4.0;
	F[14] = (1.0 - lcrd[1]*lcrd[1])*(1.0 - lcrd[0])*(1.0 + lcrd[2])/4.0;

	// узлы, для которых вторая координата равна 0
	F[16] = (1.0 - lcrd[2]*lcrd[2])*(1.0 + lcrd[0])*(1.0 - lcrd[1])/4.0;
	F[17] = (1.0 - lcrd[2]*lcrd[2])*(1.0 + lcrd[0])*(1.0 + lcrd[1])/4.0;
	F[18] = (1.0 - lcrd[2]*lcrd[2])*(1.0 - lcrd[0])*(1.0 + lcrd[1])/4.0;
	F[19] = (1.0 - lcrd[2]*lcrd[2])*(1.0 - lcrd[0])*(1.0 - lcrd[1])/4.0;


}

void CEL::dFF_25(double **dF, double *lcrd)
{
	//!!!!!!!!!!!! не рационаьный код, нужно избавиться от массива L и циклов
	int iform,i;
	double t,n,z;
	t = lcrd[0];
	n = lcrd[1];
	z = lcrd[2];

	double L[20][3];

	iform = 0; L[iform][0] = 1.0; L[iform][1] = -1.0; L[iform][2] = -1.0;
	iform = 1; L[iform][0] = 1.0; L[iform][1] = 1.0; L[iform][2] = -1.0;
	iform = 2; L[iform][0] = -1.0; L[iform][1] = 1.0; L[iform][2] = -1.0;
	iform = 3; L[iform][0] = -1.0; L[iform][1] = -1.0; L[iform][2] = -1.0;
	iform = 4; L[iform][0] = 1.0; L[iform][1] = -1.0; L[iform][2] = 1.0;
	iform = 5; L[iform][0] = 1.0; L[iform][1] = 1.0; L[iform][2] = 1.0;
	iform = 6; L[iform][0] = -1.0; L[iform][1] = 1.0; L[iform][2] = 1.0;
	iform = 7; L[iform][0] = -1.0; L[iform][1] = -1.0; L[iform][2] = 1.0;
	iform = 8; L[iform][0] = 1.0; L[iform][1] = 0.0; L[iform][2] = -1.0;
	iform = 9; L[iform][0] = 0.0; L[iform][1] = 1.0; L[iform][2] = -1.0;
	iform = 10; L[iform][0] = -1.0; L[iform][1] = 0.0; L[iform][2] = -1.0;
	iform = 11; L[iform][0] = 0.0; L[iform][1] = -1.0; L[iform][2] = -1.0;
	iform = 12; L[iform][0] = 1.0; L[iform][1] = 0.0; L[iform][2] = 1.0;
	iform = 13; L[iform][0] = 0.0; L[iform][1] = 1.0; L[iform][2] = 1.0;
	iform = 14; L[iform][0] = -1.0; L[iform][1] = 0.0; L[iform][2] = 1.0;
	iform = 15; L[iform][0] = 0.0; L[iform][1] = -1.0; L[iform][2] = 1.0;
	iform = 16; L[iform][0] = 1.0; L[iform][1] = -1.0; L[iform][2] = 0.0;
	iform = 17; L[iform][0] = 1.0; L[iform][1] = 1.0; L[iform][2] = 0.0;
	iform = 18; L[iform][0] = -1.0; L[iform][1] = 1.0; L[iform][2] = 0.0;
	iform = 19; L[iform][0] = -1.0; L[iform][1] = -1.0; L[iform][2] = 0.0;

	for (i = 0; i<8; i++)
	{
		dF[0][i] = L[i][0]*(1+n*L[i][1])*(1+z*L[i][2])*(2*t*L[i][0]+n*L[i][1]+z*L[i][2] - 1)/8;
		dF[1][i] = L[i][1]*(1+t*L[i][0])*(1+z*L[i][2])*(t*L[i][0]+2*n*L[i][1]+z*L[i][2] - 1)/8;
		dF[2][i] = L[i][2]*(1+t*L[i][0])*(1+n*L[i][1])*(t*L[i][0]+n*L[i][1]+2*z*L[i][2] - 1)/8;
	}

	for (i = 9; i<16; i+=2)
	{
		dF[0][i] = 0.0 - t*(1+n*L[i][1])*(1+z*L[i][2])/2;
		dF[1][i] = (1-t*t)*L[i][1]*(1+z*L[i][2])/4;
		dF[2][i] = (1-t*t)*L[i][2]*(1+n*L[i][1])/4;
	}

	for (i = 8; i<15; i+=2)
	{
		dF[0][i] = L[i][0]*(1-n*n)*(1+z*L[i][2])/4;
		dF[1][i] = 0.0 - n*(1+t*L[i][0])*(1+z*L[i][2])/2;
		dF[2][i] = (1-n*n)*L[i][2]*(1+t*L[i][0])/4;
	}

	for (i = 16; i<20; i++)
	{
		dF[0][i] = (1-z*z)*L[i][0]*(1+n*L[i][1])/4;
		dF[1][i] = (1-z*z)*L[i][1]*(1+t*L[i][0])/4;
		dF[2][i] = 0.0 - z*(1+t*L[i][0])*(1+n*L[i][1])/2;
	}

}