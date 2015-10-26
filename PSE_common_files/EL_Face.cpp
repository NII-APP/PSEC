#include "StdAfx.h"
#include "EL.h"

double CEL::FaceArea(int iface)
{
	int ipo,j;
	double **J,detJ;
	double area;
	J = NULL;
	J = MM->MEM_NEW(J,NORT,NORT);

	area = 0.0;
	for (ipo=0; ipo<Nfp; ipo++)
	{
		Jacoby(dFFfloc[iface][ipo],J);
		DetJacoby(J,&detJ);
		area += detJ*Wfint[iface][ipo];
	}

	J = MM->MEM_DEL(J,NORT,NORT);
	return(area);
}

void CEL::FaceNodeNums(int iface, int *NodeNums, int *NNf)
{
	int i;
	for (i=0; i<NNface; i++)
	{
		NodeNums[i] = ind[ FaceNodesList[iface][i] ];
	}
	*NNf = NNface;
}