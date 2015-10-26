#include "StdAfx.h"
#include "SE.h"

extern void xolp( double *A, long long int *D, int k);
extern void xolsinglesolv( double *L, long long int *D, int k, double *B );

void CSE::LoadDisp()
{
	int i,j,k,kk;
	FILE *fp;
	char strfrom[256];

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

	// необходимо осуществить считывание вектора перемещений —Ё, включающего данный —Ё
	int ns1,ni1,nie1,nne1,nn1,nse1,kort1,*RENse,inode;
	double *rv;
	sprintf(strfrom,"%s\\uuu%s.bin",sedat->PathIncludedInSE,sedat->namenumIncludedInSE);
	fp = fopen(strfrom,"rb");
	fread(&nvect,sizeof(int),1,fp);
	fread(&ni1,sizeof(int),1,fp);
	fread(&ns1,sizeof(int),1,fp);
	fread(&kort1,sizeof(int),1,fp);
	nn1 = ni1+ns1;
	nne1 = nn1*kort1;
	nse1 = ns1*kort1;
	RENse = new int [nn1];
	rv = new double[nne1];

	DISP = new double*[nvect];
	for(k=0; k<nvect; k++)
	{
		DISP[k] = new double [sedat->NNE];
		MM->QQZERO(DISP[k],sedat->NNE);
	}

	fread(RENse,sizeof(int),nn1,fp);
	for (j=0; j<nvect; j++)
	{
		fread(rv,sizeof(double),nne1,fp);
		for (k=0; k<nn1; k++)
		{
			for (kk=0; kk<kort1; kk++)
			{
				inode = RENINV[ RENse[k] ];	
				if ( inode >= 0 )
				{
					if ( inode >= sedat->NNE )
					{
						j=j;
					}
					DISP[j][ inode*KORT + kk ] = rv[ k*kort1 + kk ];
				}
			}
		}
	}
	fclose(fp);

	delete []rv;
	delete []RENse;

	//считывание вектора сил внутренних узлов данного суперэлемента
	sprintf(strfrom,"%s\\fff%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp = fopen(strfrom,"rb");
	fread(&nvect,sizeof(int),1,fp);
	fread(&ni1,sizeof(int),1,fp);
	fread(&ns1,sizeof(int),1,fp);
	fread(&kort1,sizeof(int),1,fp);
	nn1 = ni1+ns1;
	nne1 = nn1*kort1;
	nse1 = ns1*kort1;
	RENse = new int [nn1];
	rv = new double[nne1];

	fread(RENse,sizeof(int),nn1,fp);
	for (j=0; j<nvect; j++)
	{
		fread(rv,sizeof(double),nne1,fp);
		for (k=0; k<(ni1*kort1); k++)
		{
			DISP[j][k] = rv[k];
		}
	}
	fclose(fp);

	delete []rv;
	delete []RENse;
	delete []RENINV;	

	printf("%d vectors loaded...\n",nvect);

}

void CSE::DispTransform()
{
	int i,j;
	
	DISPcolTransform(0,nvect);

	//необходимо организовать сохранение!!!

}

void CSE::DISPcolTransform(int colstart, int colend)
{

	int i,j,k,ip,jnode,isscol,fcol,iur;

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
					for (fcol=colstart; fcol<colend; fcol++)
					{
						DISP[fcol][iur] -= DISP[fcol][NIE+isscol]*STIS[isscol][ip];
					}
					ip++;
				}
			}
		}
	}

	printf("\n");
	for (fcol=colstart; fcol<colend; fcol++)
	{
		if ( sedat->fl_PCG_internal == true )
		{
			printf("ivect = %d // %d  \r",fcol, colend);
			PCGSolve(DISP[fcol]);
		}
		else
		{
			printf("ivect = %d // %d  \r",fcol, colend);
			xolsinglesolv(STII, STII_ENV, NIE, DISP[fcol]);
		}
	}
}
