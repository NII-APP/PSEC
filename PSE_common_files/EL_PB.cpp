#include <stdio.h>
#include "stdafx.h"
#include "EL.h"
#include "math.h"

int CEL::PBVP (int IPO)
{
/*BXO�> IEL,IPO - HOMEPA ��EMEHTA � TO�K� �HTE�P�POBAH�� B HEM
C:PO��: - ���������� ����������� ������� ����� �� ���������� ����������� ��� ���� ��������� � ���� ����� ��������������
        - ������� ������ ����� �������������� � ������ �������� � ������� ������������� Hi,Hj
C<B�XO�< -    */
	int i,j,k;
	double J[3][3],JA[3][3],MMM[3][3],detJ;
	
//���������� ������� �����
	for(i=0; i<KORT; i++)
	{
		for(j=0; j<KORT; j++)
		{
			J[i][j]=0;
			JA[i][j]=0;
		}
	}

	for (i=0; i<KORT; i++)
	{
		for(j=0; j<KORT; j++)
		{
			J[i][j] = 0;
			for(k=0; k<KN; k++)
			{
				//J[i][j] += dFke[IPO][i][k]*KOORD[NND[IEL][k]][j];
				//J[i][j] += dFke[IPO][i][k]*NodesCRD[k*KORT+j];
				J[i][j] += dFke[IPO][i][k]*fullcrd[ind[k]*KORT+j];
			}
		}
	}
	// ���������� ������� �������������� �������
	MMM[0][0] = J[1][1]*J[2][2] - J[2][1]*J[1][2];
	MMM[0][1] = J[1][0]*J[2][2] - J[2][0]*J[1][2];
	MMM[0][2] = J[1][0]*J[2][1] - J[2][0]*J[1][1];
	MMM[1][0] = J[0][1]*J[2][2] - J[0][2]*J[2][1];
	MMM[1][1] = J[0][0]*J[2][2] - J[2][0]*J[0][2];
	MMM[1][2] = J[0][0]*J[2][1] - J[2][0]*J[0][1];
	MMM[2][0] = J[0][1]*J[1][2] - J[1][1]*J[0][2];
	MMM[2][1] = J[0][0]*J[1][2] - J[1][0]*J[0][2];
	MMM[2][2] = J[0][0]*J[1][1] - J[1][0]*J[0][1];

	detJ=J[0][0]*MMM[0][0] - J[0][1]*MMM[0][1] + J[0][2]*MMM[0][2];
		
	if(detJ <=0 )
	{
		i=i;
	}

	// ���������� ������� �������� � ��������
	for (i=0; i<KORT; i++)
	{
		for(j=0; j<KORT; j++)
		{
			JA[i][j] = pow(double(-1),i+j) * MMM[j][i] / detJ;
		}
	}

//���������� ����������� ������� ����� �� ���������� ����������� ��� �������� IEL � ����� �������������� IPO

	for (i=0; i<KORT; i++)
	{
		for (j=0; j<KN; j++)
		{
			dFrz[i][j] = 0;
			for (k=0; k<KORT; k++)
			{
				dFrz[i][j] += JA[i][k]*dFke[IPO][k][j];
			}
		}
	}

	P[IPO].VP = Hi[IPO]*detJ;
	if (KN==10) P[IPO].VP = P[IPO].VP/6;

	if (P[IPO].VP <=0 )
	{
		isVPbelowzero = 1;
		P[IPO].VP = abs(P[IPO].VP);
	}

	return(0);
}
/*********************************************************************/