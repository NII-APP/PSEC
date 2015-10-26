#include "StdAfx.h"
#include "SE.h"
#include "math.h"


void CSE::StifPrecondition()
{
	int i,j,k,isel,imax,nn,*FLAG,tmpi,adelem;
	double max,tmp,adsum,tol;

//	double pcpercent = 0.0;
	double pcpercent = 0.1;
//	double pcpercent = 2.00;
	tol = 0.1;

	adelem = 0;

	if ( STII_ENVPC != NULL )
	{
		for (i=0; i<NIE; i++)
		{
			delete [] STII_ENVPC[i];
		}
		delete [] STII_ENVPC;
		STII_ENVPC = NULL;
	}

	if ( STIIPC != NULL )
	{
		for (i=0; i<NIE; i++)
		{
			delete [] STIIPC[i];
		}
		delete [] STIIPC;
		STIIPC = NULL;
	}

	STII_ENVPC = new int*[NIE];
	STIIPC = new double*[NIE];

	FLAG = new int[NIE];
	
	for (i=0; i<NIE; i++)
	{
		FLAG[i] = 0;
	}
	isel = 0;

	for (i=0; i<NIE; i++)
	{
		if ( (i%1000) == 0 )
		{
			tmp = 100.0*(((double)i)/((double)NIE));
			tmpi = (int)tmp;
			printf( " Precond %d%% IUR %d\r",tmpi,i);
		}
		//определение высоты столбца предобусловливателя
//		nn = (int)(ENV[i]*pcpercent);
//		nn = (int)(STII_ENVNZ[i][0]*pcpercent + adelem);
		if (i > 0) 
		{
			nn = (int)( (STII_ENV[i] - STII_ENV[i-1])*pcpercent );
		}
		else
		{
			nn = 1;
		}
		if (nn < 1) nn = 1;
		if (i > 0)
		{
			if (nn > (STII_ENV[i] - STII_ENV[i-1])) nn = (STII_ENV[i] - STII_ENV[i-1]);
		}
		else
		{
			nn = 1;
		}

		if (i == 34)
		{
			i=i;
		}
		
		isel = 0;
		if ( nn != 1 )
		{
			//отбор максимальных элементов столбца
			/*for (j=0; j<ENV[i]; j++)
			{
				if ( fabs(STIF[i][j])/STIF[i][0] > tol )
				{
					FLAG[j] = 1;
					isel++;
					if (isel == nn) break;
				}
			}*/
			for (k=0; k<nn; k++)
			{
				max = -1.0;
				imax = -1;
				for(j=0; j< (STII_ENV[i] - STII_ENV[i-1]); j++)
				{
					if ( FLAG[j] == 0 )
					{
						if ( fabs(STII[ STII_ENV[i]-j ]) > max ) 
						{
							max = fabs(STII[ STII_ENV[i]-j ]);  
							imax = j;
						}
					}
				}
				if ( (max/STII[ STII_ENV[i] ]) > tol )
				{
					FLAG[imax] = 1;
					isel++;
				}
				else
				{
					break;
				}
			}
			nn = isel;
		}
		
		//выделение памяти под столбец предобусловливателя
		STII_ENVPC[i] = new int[nn+1];
		STII_ENVPC[i][0] = nn;
		STIIPC[i] = new double [nn];

		adsum = 0.0;
		if ( nn > 1 )
		{
			isel = 0;
			for (k=0; k<(STII_ENV[i] - STII_ENV[i-1]); k++)
			{
				if (FLAG[k] == 1)
				{
					STII_ENVPC[i][isel+1] = i - k;
					STIIPC[i][isel] = STII[STII_ENV[i]-k];
					isel++;
					FLAG[k] = 0;
					if (isel > nn)
					{//error
						i=i;
					}
				}
				else
				{
//					adsum += fabs(STIF[i][k]);
//					STIFPC[i-k][0] += fabs(STIF[i][k]); 
				}
			}
		}
		else
		{
			FLAG[0] = 0;
		}
		//диагональный элемент копируется в любом случае
		STIIPC[i][0] = STII[STII_ENV[i]] + adsum;
//		STIFPC[i][0] = 1.0;
		STII_ENVPC[i][1] = i;
	}

	delete [] FLAG;

	printf( "\n");
}


void CSE::PrecondHolSol(double *B, double *C)
{
	int i,j,nach,pp;
	double s;
	
	for(i=0; i<NIE; i++)
	{
		B[i] = C[i];
	}

	for(i=0; i<NIE; i++)
	{
		s=0.0;
		for(j=STII_ENVPC[i][0]; j>1; j--)
		{			
			s += STIIPC[i][j-1]*B[STII_ENVPC[i][j]];
		}
		B[i]=(B[i]-s)/STIIPC[i][0];
	}
	
	for(i=NIE-1; i>=0; i--)
	{
		B[i] = B[i]/STIIPC[i][0];
		for(j=STII_ENVPC[i][0]; j>1; j--)
		{
			B[STII_ENVPC[i][j]] -= B[i]*STIIPC[i][j-1];
		}
	}

}

void CSE::StifMultPC(double *C, double *B)
{
	int i,j,nach,pp;
	double s;
	
	for(i=0; i<NIE; i++)
	{
		C[i] = 0.0;
	}

	for(i=0; i<NIE; i++)
	{
		//недиагональные элементы участвуют в умножении дважды из-за симметрии
		for(j=2; j<=STII_ENVNZ[i][0]; j++)
		{
			//нижний треугольник
			C[i] += STIINZ[i][j-1]*B[STII_ENVNZ[i][j]];
			//верхний треугольник
			C[STII_ENVNZ[i][j]] += STIINZ[i][j-1]*B[i];
		}
		//для диагонального элемента
		C[i] += STIINZ[i][0]*B[i];
	}

}

void CSE::PCGSolve( double *U )
{
	int i,j,k,iter;
	double s,tol,eps,al,bet,s1,s2,modB;
	double *B,*R,*Rold,*Z,*Zold,*P,*tmp;

	int KUR = NIE;

	eps = 1.e-8;

	B = new double[KUR];
	R = new double[KUR];
	Rold = new double[KUR];
	Z = new double[KUR];
	Zold = new double[KUR];
	P = new double[KUR];
	tmp = new double[KUR];

	//1 - выбор начального приближения
	modB = 0.0;
	for (i=0; i<KUR; i++)
	{
		B[i] = U[i];
		U[i] = 0.0;
		modB += B[i]*B[i];
	}
	modB = sqrt(modB);

	// поиск первой невязки
	StifMultPC(R,U);
	for (i=0; i<KUR; i++)
	{
		R[i] = B[i] - R[i];
	}

	PrecondHolSol(Z,R);
	for (i=0; i<KUR; i++)
	{
		P[i] = Z[i];
	}

	//итерации
	tol = 1.0;
	iter = 0;
	while ( tol > eps )
	{
		iter++;
		if ( iter%100 == 0 )
		{
			iter = iter;
		}
		// 1
		StifMultPC(tmp,P);
		s1 = 0.0;
		s2 = 0.0;
		for (i=0; i<KUR; i++)
		{
			s1 += Z[i]*R[i];
			s2 += P[i]*tmp[i];
		}
		al = s1/s2;

		//2
		for (i=0; i<KUR; i++)
		{
			Zold[i] = Z[i];
			Rold[i] = R[i];
			U[i] += al*P[i];
			R[i] -= al*tmp[i];
		}
		
		//3
		PrecondHolSol(Z,R);

		//4
		s1 = 0.0;
		s2 = 0.0;
		for (i=0; i<KUR; i++)
		{
			s1 += Z[i]*R[i];
			s2 += Zold[i]*Rold[i];
		}
		bet = s1/s2;

		//5
		s1 = 0.0;
		for (i=0; i<KUR; i++)
		{
			P[i] = Z[i] + bet*P[i];
			s1 += R[i]*R[i];
		}
		s1 = sqrt(s1);

		tol = s1/modB;
	}

	printf("  it = %d\n",iter);

	delete [] B;
	delete [] R;
	delete [] Rold;
	delete [] Z;
	delete [] Zold;
	delete [] P;
	delete [] tmp;
	

}