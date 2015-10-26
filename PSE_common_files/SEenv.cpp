#include "StdAfx.h"
#include "SE.h"


void CSE::STIIenv ()
{

	int i,j,z,min,KN;
	int TMP[100],tmp=0;

	if (STII_ENV != NULL) delete [] STII_ENV; 
	STII_ENV = new long long int [NIE+1];

	min=0;

	for(i=0; i<NI; i++)//профильная схема
	{
		min=i;
		tmp=0;
		for(j=0; j<NEL; j++)
		{
			KN = IND[j][1];
			for(z=0;z<KN;z++)
			{
				if(IND[j][z+2] == i)
				{
					TMP[tmp]=j;
					tmp++;
					break;
				}
			}
		}
		for(j=0; j<tmp; j++)
		{
			KN = IND[TMP[j]][1];
			for(z=0; z<KN; z++)
			{
				if ( IND[TMP[j]][z+2] < min )  min=IND[TMP[j]][z+2];
			}
		}
		for(j=0; j<KORT; j++)
		{
			if ( (i*KORT+j) > 0 ) STII_ENV[i*KORT+j] = STII_ENV[i*KORT+j - 1] + (i-min)*KORT + 1 + j;
			else STII_ENV[i*KORT+j] = 0;
		}
	}
	STII_ENV[NIE] = STII_ENV[NIE-1] + 1;

}


void CSE::STISenv ()
{

	int i,j,z,m,min,KN;
	int tmp=0,env=0,node,nonzeronum,counter,firstnum;
	int *TMP,*ENV,*FLAG;

	STIS_ENV = new int*[NS];


	TMP = new int[NEL];
	ENV = new int[NI*3];
	FLAG = new int[NN];

	for(i=NI; i<(NI+NS); i++)//безнулевая схема
	{
		tmp=0;
		for(j=0; j<NEL; j++)
		{
			KN = IND[j][1];
			for(z=0;z<KN;z++)
			{
				if(IND[j][z+2] == i)
				{
					TMP[tmp]=j;
					tmp++;
					break;
				}
			}
		}

		//определение всех узлов, входящих в столбец узла i
		for (j=0; j<NI; j++) FLAG[j] = 0;

		nonzeronum = 0;
		for(j=0; j<tmp; j++)
		{
			KN = IND[TMP[j]][1];
			for(z=0; z<KN; z++)
			{
				node = IND[TMP[j]][z+2];
				if ( (FLAG[node] == 0) && (node < NI) ) 
				{
					FLAG[node] = 1;
					nonzeronum++;
				}
			}
		}

		//переход к сокращенной схеме хранения профиля столбца 
		//для каждой непрерывной последовательности ненулевых элементов (по узлам, не по степеням свободы!) указывается номер первого узла и их количество
		counter = 0;
		firstnum = 0;
		env = 0;
		for (j=0; j<NI; j++)
		{
			if (FLAG[j] == 1)
			{
				if (counter == 0)
				{
					firstnum = j;
				}
				counter++;
			}
			if ( (FLAG[j] == 0 && counter > 0) || (j == (NI-1) && counter > 0) )
			{
				ENV[env] = firstnum;
				env++;
				/*if (env > 1)
				{
					ENV[env] = ENV[env-2] + counter;
					env++;
				}
				else
				{
					ENV[env] = counter;
					env++;
				}*/
				ENV[env] = counter;
				env++;
				counter = 0;
			}
		}

		//вся запись осуществляется по узлам, а не по степеням свободы
		STIS_ENV[i-NI] = new int [env + 2];
		STIS_ENV[i-NI][0] = nonzeronum; //первое число показывает количество ненулевых элементов (узлов!) в столбце IS
		STIS_ENV[i-NI][1] = env; //количество элементов в массиве ENV
		for (j=0; j<env; j++)
		{
			STIS_ENV[i-NI][2+j] = ENV[j];
		}
	}

	delete []TMP;
	delete []ENV;
	delete []FLAG;
	

}


void  CSE::STIINZenv (int flagfull)
{
	int		IEL, EL, KUL,KN;
	int		n, i, j, ik, jk, iu, ju, IL, totalnum,tmpi;
	double	tmp, totalmemstif;
	char strl[256];
	int *Flag, *NonZerNum, nzn, nzntr, numcomp,iur,jur,inode,iort;


	int NIE,NI;

	if (flagfull == 0)
	{
		NIE = this->NIE;
		NI = this->NI;
	}
	if (flagfull == 1)
	{
		NIE = this->NNE;
		NI = this->NN;
	}

	if ( STII_ENVNZ != NULL )
	{
		for (i=0; i<NIE; i++)
		{
			delete [] STII_ENVNZ[i];
		}
		delete [] STII_ENVNZ;
		STII_ENVNZ = NULL;
	}

	STII_ENVNZ = new int*[NIE];
	Flag = new int[NN];
	NonZerNum = new int[NN];

	for (i=0; i<NN; i++)
	{
		Flag[i] = 0;
		NonZerNum[i] = 0;
	}

	totalnum = 0;
	for ( inode = 0; inode < NI; inode++)
	{
		if ( (inode%1000) == 0 )
		{
			tmp = 100.0*(((double)inode)/((double)NI));
			tmpi = (int)tmp;
			printf( " ENVNZ %d%% INODE %d\r",tmpi,inode);
			//sprintf( strl," ENVNZ %d%% INODE %d",tmpi,inode);
			//statusbar( strl );
		}
		nzn = 0;
		for ( IEL = 0; IEL <NEL; IEL++ )
		{
			KN = IND[IEL][1];
			for (j=0; j<KN; j++)
			{
				if ( (IND[IEL][2+j]) == inode )
				{
					break;
				}
			}
			if (j < KN)
			{
				for (j=0; j<KN; j++)
				{
					if ( Flag[ IND[IEL][j+2] ] == 0 )
					{
						Flag[ IND[IEL][j+2] ] = 1;
						NonZerNum[nzn] = IND[IEL][j+2];
						nzn++;
					}
				}
			}
		}

		nzntr = 0;
		for (i=0; i<inode; i++)
		{
			if ( Flag[i] == 1 ) 
			{
				nzntr++;
			}
		}
		
		for (iort = 0; iort<KORT; iort++)
		{
			iur = inode*KORT + iort;
			numcomp = nzntr*KORT + iort + 1;
			STII_ENVNZ[iur] = new int [numcomp+1]; 
			STII_ENVNZ[iur][0] = numcomp;
			jur = 0;
			for (i=0; i<inode; i++)
			{
				if ( Flag[i] == 1 ) 
				{
					for (j=0; j<KORT; j++)
					{
						STII_ENVNZ[iur][numcomp-jur] = i*KORT+j;
						jur++;
						totalnum++;
					}
				}
			}
			for (j=0; j<=iort; j++) //диагональный блок
			{
				STII_ENVNZ[iur][numcomp-jur] = inode*KORT + j;
				jur++;
			}
		}

		for (i=0; i<nzn; i++)
		{
			Flag[ NonZerNum[i] ] = 0;
		}
	}

	totalmemstif = ((double)(totalnum*(4+8)))/(1024*1024);

	delete [] Flag;
	delete [] NonZerNum;

	printf("\n");

}

