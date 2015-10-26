#include <stdio.h>
#include "stdafx.h"
#include "EL.h"
#include "string.h"

int CEL::GE_MAKE(int ipo)
{
	int i,j,z;
	double s;


	for(i=0; i<KN; i++)
	{
		for (j=0; j<KORT; j++)
		{
			NUR[i*KORT+j] = ind[i]*KORT + j;
		}
	}

	for(i=0; i<KN; i++)
	{
		B[0][i*KORT] = dFrz[0][i];
		B[1][i*KORT+1] = dFrz[1][i];
		B[2][i*KORT+2] = dFrz[2][i];
		B[3][i*KORT] = dFrz[1][i]; B[3][i*KORT+1] = dFrz[0][i];
		B[4][i*KORT+1] = dFrz[2][i]; B[4][i*KORT+2] = dFrz[1][i];
		B[5][i*KORT] = dFrz[2][i]; B[5][i*KORT+2] = dFrz[0][i];
	}

	//B транспонированное умножаем на D

	MM->QQZERO(BD,KU,LE);

	for(i=0; i<KU; i++)
	{
		for(j=0; j<LE; j++)
		{
			for(z=0; z<LE; z++)
			{
				BD[i][j]+=B[z][i]*D[z][j];
			}
		}
	}
	//BD умножаем на B и на объем точки интегрирования. Результат помещаем в GE

	for(i=0; i<KU; i++)
	{
		for(j=0; j<KU; j++)
		//for(j=i; j<KU; j++)
		{
			s=0;
			for(z=0; z<LE; z++)
			{
				s+=BD[i][z]*B[z][j];
			}
			
			GE[i][j] += s*P[ipo].VP;
		}
	}

	return(0);
}

