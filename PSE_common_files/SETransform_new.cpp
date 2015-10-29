#include "StdAfx.h"
#include "SE.h"

extern void xolp( double *A, long long int *D, int k);
extern void xolsinglesolv( double *L, long long int *D, int k, double *B );
extern void xolmultisolv_1th( double *L, long long int *D, int k, double **B, int nv );

void CSE::StiffMatrixTransform()
{
	int i,j;
	
	//необходимо задать закрепления!!!

	if ( sedat->SElevel == 0 ) 
	{
		MatrixFixing();
	}
	xolp(STII, STII_ENV, NIE);
	if (NS != 0)
	{
		SScolTransform(0,NSE);
	}
}

void CSE::MatrixFixing()
{
	//закрепление за счет задания большого числа на диагонали (матрица IS не участвует)
	int i,j,k;

	for (i=0; i<NIE; i++)
	{
		if ( FIX[i] != 0 ) 
		{
			STII[ STII_ENV[i] ] = bignumber;
			if ( sedat->fl_PCG_internal == true )
			{
				STIINZ[i][0] = bignumber;
			}
		}
	}

	for (i=0; i<NSE; i++)
	{
		if ( FIX[NIE+i] != 0 )
		{
			if (i == 0)
			{
				STSS[0] = bignumber;
			}
			else
			{
				k = (int)(((i)*(i)-(i))/2 + 0.1);
				k += i;
				k += i;
				STSS[k] = bignumber;
			}
		}
	}
}

void CSE::SScolTransform(int colstart,int colend)
{
	int i,j,isscol,ncol;
	double **COLSS,**COLIS;
	COLSS = NULL;
	COLIS = NULL;

	ncol = colend - colstart;

	printf("SS calc....\n");

	COLSS = MM->MEM_NEW(COLSS,ncol,NSE);
	COLIS = MM->MEM_NEW(COLIS,ncol,NIE);

	if ( (isscol%10) == 0 ) printf("iss = %d // %d\r",isscol,colend);
	for (isscol = colstart; isscol<colend; isscol++)
	{
		IScolExpand(COLIS, isscol);
	}

	xolmultisolv_1th(STII, STII_ENV, NIE, COLIS, ncol);
	
	SScolCalc(COLSS,COLIS,isscol);

	//далее нужна синхронизация
	for (isscol = colstart; isscol<colend; isscol++)
	{
		SScolInsert(COLSS,isscol);
	}
	//конец синхронизации

	printf("\n");

	MM->MEM_DEL(COLSS,ncol,NSE);
	MM->MEM_DEL(COLIS,ncol,NIE);
}

void CSE::IScolExpand(double *COLIS,int isscol)
{
	int i,j,k,ip,jnode,iur;

	MM->QQZERO(COLIS,NIE);

	ip = 0;
	jnode = (int)(isscol/KORT);
	for(i=0; i<STIS_ENV[jnode][1]; i+=2)
	{
		for (j=0; j<STIS_ENV[jnode][i+2+1]; j++)
		{
			for (k=0; k<KORT; k++)
			{
				iur = STIS_ENV[jnode][i+2]*KORT + j*KORT + k;
				COLIS[iur] = STIS[isscol][ip];
				ip++;
			}
		}
	}
}

void CSE::SScolCalc(double *COLSS, double *COLIS, int jsscol)
{
	int i,j,k,ip,jnode,isscol,iur;

	MM->QQZERO(COLSS,NSE);
	for (isscol = 0; isscol <= jsscol; isscol++)
	{
		ip = 0;
		jnode = (int)(isscol/KORT);
		for(i=0; i<STIS_ENV[jnode][1]; i+=2)
		{
			for (j=0; j<STIS_ENV[jnode][i+2+1]; j++)
			{
				for (k=0; k<KORT; k++)
				{
					iur = STIS_ENV[jnode][i+2]*KORT + j*KORT + k;
					COLSS[isscol] += COLIS[iur]*STIS[isscol][ip];
					ip++;
				}
			}
		}
	}
}


void CSE::SScolInsert(double *COLSS, int isscol)
{
	int i,j,ip;

	if (isscol > 0)
	{
		ip = (int)((isscol*isscol-isscol)/2 + 0.1);
		ip += isscol;
	}
	else
	{
		ip = 0;
	}

	//необходимо вставить mutex для параллельной версии !!!

	if (ip == 0)
	{
		ip = ip;
	}
	
	for (i=0; i<=isscol; i++)
	{
		STSS[ip+i] -= COLSS[i];
	}
}