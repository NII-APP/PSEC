#include "StdAfx.h"
#include "EL.h"
#include "math.h"

void CEL::Jacoby(double **dFloc, double **J)
{
	int i,j,k;
	//!!!!!!!! согласовать n и NN

	//вычисление матрицы Якоби
	for (i=0; i<NORT; i++)
	{
		for(j=0; j<NORT; j++)
		{
			J[i][j] = 0;
			for(k=0; k<NN; k++)
			{
				J[i][j] += dFloc[i][k]*fullcrd[ind[k]*NORTfullcrd+j];
			}
		}
	}
}

void CEL::DetJacoby(double **J, double *detJ)
{
	// вычисление определителя матрицы Якоби
	if (NORT == 3)
	{
		*detJ = J[0][0]*(J[1][1]*J[2][2] - J[2][1]*J[1][2]) - 
			J[0][1]*(J[1][0]*J[2][2] - J[2][0]*J[1][2]) + 
			J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]);
	}
	if (NORT == 2)
	{
		*detJ = J[0][0]*J[1][1] - J[1][0]*J[0][1];
	}
}

void CEL::InvJacoby(double **J, double detJ)
{
	int i,j;
	double MMM[3][3];
	// вычисление матрицы дополнительных миноров
	if ( NORT == 3 )
	{
		MMM[0][0] = J[1][1]*J[2][2] - J[2][1]*J[1][2];
		MMM[0][1] = J[1][0]*J[2][2] - J[2][0]*J[1][2];
		MMM[0][2] = J[1][0]*J[2][1] - J[2][0]*J[1][1];
		MMM[1][0] = J[0][1]*J[2][2] - J[0][2]*J[2][1];
		MMM[1][1] = J[0][0]*J[2][2] - J[2][0]*J[0][2];
		MMM[1][2] = J[0][0]*J[2][1] - J[2][0]*J[0][1];
		MMM[2][0] = J[0][1]*J[1][2] - J[1][1]*J[0][2];
		MMM[2][1] = J[0][0]*J[1][2] - J[1][0]*J[0][2];
		MMM[2][2] = J[0][0]*J[1][1] - J[1][0]*J[0][1];
	}
	if ( NORT == 2 )
	{
		MMM[0][0] = J[1][1];
		MMM[0][1] = J[1][0];
		MMM[1][0] = J[0][1];
		MMM[1][1] = J[0][0];
	}
		
	// вычисление матрицы обратной к якобиану
	for (i=0; i<NORT; i++)
	{
		for(j=0; j<NORT; j++)
		{
			J[i][j] = pow(double(-1),i+j) * MMM[j][i] / detJ;
		}
	}
}

