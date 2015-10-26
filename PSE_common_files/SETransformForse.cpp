#include "StdAfx.h"
#include "SE.h"

extern void xolp( double *A, long long int *D, int k);
extern void xolsinglesolv( double *L, long long int *D, int k, double *B );


void CSE::LoadForse()
{
	int i,j,k,kk;
	FILE *fp;
	char strfrom[256];

	if (sedat->SElevel == 0)
	{
		ForseRead();
	}
	else
	{
		//составление массива обратной перенумерации степеней свободы
		int *RENINV;
		RENINV = new int [sedat->fullnn];
		for (i=0; i<sedat->fullnn; i++)
		{
			RENINV[i] = -1;
		}
		for (i=0; i<NN; i++)
		{
			RENINV[ REN[i] ] = i;
		}

		//для суперэлементов старших уровней необходимо осуществить ассемблирование вектора нагрузки
		for (i=0; i<sedat->NInclSEs; i++)
		{
			int ns1,ni1,nie1,nne1,nn1,nse1,kort1,*RENse;
			double *lv;
			sprintf(strfrom,"%s\\fff%s.bin",sedat->PathInclSEs[i],sedat->namenumInclSEs[i]);
			fp = fopen(strfrom,"rb");
			fread(&nvect,sizeof(int),1,fp);
			fread(&ni1,sizeof(int),1,fp);
			fread(&ns1,sizeof(int),1,fp);
			fread(&kort1,sizeof(int),1,fp);
			nn1 = ni1+ns1;
			nne1 = nn1*kort1;
			nse1 = ns1*kort1;
			RENse = new int [nn1];
			lv = new double[nne1];
			if (i==0)
			{
				FORSE = new double*[nvect];
				for(k=0; k<nvect; k++)
				{
					FORSE[k] = new double [sedat->NNE];
					MM->QQZERO(FORSE[k],sedat->NNE);
				}
			}

			fread(RENse,sizeof(int),nn1,fp);
			for (j=0; j<nvect; j++)
			{
				fread(lv,sizeof(double),nne1,fp);
				for (k=0; k<ns1; k++)
				{
					for (kk=0; kk<kort1; kk++)
					{
						if ( (RENINV[ RENse[ni1+k] ]*KORT + kk) < 0 || (RENINV[ RENse[ni1+k] ]*KORT + kk) >= sedat->NNE  )
						{
							k=k;
						}
						FORSE[j][ RENINV[ RENse[ni1+k] ]*KORT + kk ] += lv[ (ni1+k)*kort1 + kk ];
					}
				}
			}
			fclose(fp);
			delete []lv;
			delete []RENse;
		}
		delete []RENINV;		
	}

	printf("%d vectors loaded...\n",nvect);
}


void CSE::ForseTransform()
{
	int i,j;

	if ( sedat->SElevel == 0 ) ForseFixing();

	if ( NS == 0 )
	{
		//для верхнего уровня сразу осуществляется решение СЛАУ, запись векторов результата
		DISP = new double*[nvect];
		for (i=0; i<nvect; i++)
		{
			DISP[i] = new double[sedat->NNE];
			for (j=0; j<sedat->NNE; j++)
			{
				DISP[i][j] = FORSE[i][j];
			}
		}
		printf("\n");
		for (i=0; i<nvect; i++)
		{
			
			if ( sedat->fl_PCG_internal == true )
			{
				printf("ivect = %d // %d  \r",i, nvect);
				PCGSolve(DISP[i]);
			}
			else
			{
				printf("ivect = %d // %d  \r",i, nvect);
				xolsinglesolv(STII, STII_ENV, NIE, DISP[i]);
			}
		}
	}
	else
	{
		FORSEcolTransform(0,nvect);
	}
}

void CSE::ForseFixing()
{
	int i,j,k;

	for (i=0; i<sedat->NNE; i++)
	{
		if (FIX[i] != 0)
		{
			for (j=0; j<nvect; j++)
			{
				FORSE[j][i] = UFIX[i]*bignumber;
			}
		}
	}
}

void CSE::FORSEcolTransform(int colstart, int colend)
{
	int i,j,icol;
	double *RV=NULL;

	RV = MM->MEM_NEW(RV,NIE);
	printf("\n");
	for (icol = colstart; icol<colend; icol++)
	{
		for (j=0; j<NIE; j++) RV[j] = FORSE[icol][j];
		if ( sedat->fl_PCG_internal == true )
		{
			printf("ivect = %d // %d  \r",icol, colend);
			PCGSolve(RV);
		}
		else
		{
			printf("ivect = %d // %d  \r",icol, colend);
			xolsinglesolv(STII, STII_ENV, NIE, RV);
		}
		FORSEcolS(RV,icol);
	}

	MM->MEM_DEL(RV,NIE);
}

void CSE::FORSEcolS(double *RV,int icol)
{
	int i,j,k,ip,jnode,isscol,iur;
	double *ftmp=NULL;

	ftmp = MM->MEM_NEW(ftmp,NSE);

	for (isscol = 0; isscol < NSE; isscol++)
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
					//FORSE[icol][NIE+isscol] -= RV[iur]*STIS[isscol][ip];
					ftmp[isscol] += RV[iur]*STIS[isscol][ip];
					ip++;
				}
			}
		}
	}

	//в параллельной версии следующий код нужно заключить в критическую секцию
	for (isscol = 0; isscol < NSE; isscol++)
	{
		FORSE[icol][NIE+isscol] -= ftmp[isscol];
	}
	/// конец  критической секции

	MM->MEM_DEL(ftmp,NSE);
}