#include "StdAfx.h"
#include "MEM.h"

void MEM::GENRCM(int NN, int NI, int NEL, int NORT, double *CRD, int **IND, int *INVP, int *size_before, int *size_after, int *width_before, int *width_after, int *nbreaks)
{
	int **GR = NULL,*PERM;
	int **NDEL = NULL;
	int i,j,stnode;

	PERM = new int[NI];// соответствие новых номеров узлов старым номерам

	GRNDEL(&NDEL,IND,NEL,NI);
	GenNDGR(&GR,NEL,NN,NI,NDEL,IND);

	*size_before = SIZE_PROFILE(GR,NI,NORT);
	*width_before = WIDTH_PROFILE(GR,NI,NORT);

	if (NEL > 2) //перенумерация имеет смысл только если число составляющих элементов (суперэлементов) 3 и более
	{
		stnode = FindStartingNode(CRD,NI,NORT);
		*nbreaks = RCM(GR,NI,stnode,PERM,INVP);

		GR = MEM_DEL(GR,NI,1);
		NDEL = MEM_DEL(NDEL,NI,1);

		RenumberIND(IND,NEL,NI,INVP);
		RenumberVect(CRD,INVP,NI,NORT);

		GRNDEL(&NDEL,IND,NEL,NI);
		GenNDGR(&GR,NEL,NN,NI,NDEL,IND);

		*size_after = SIZE_PROFILE(GR,NI,NORT);
		*width_after = WIDTH_PROFILE(GR,NI,NORT);
	}
	else
	{
		*size_after = *size_before;
		*width_after = *width_before;

		for (i=0; i<NI; i++)
		{
			INVP[i] = i;
		}
	}

	GR = MEM_DEL(GR,NI,1);
	NDEL = MEM_DEL(NDEL,NI,1);
	
	delete PERM;
}

int MEM::RCM (int **GR,int NI, int stnode, int *PERM, int *INVP)
{
	int i,j,k,mask,nenums,tmp;
	int *MASK,*NENUMS,*TMP,*DEG;
	int nbreaks = 1; //количество односвязных областей
	
	MASK= new int[NI];
	NENUMS= new int[NI];
	TMP= new int[NI];
	DEG= new int[NI];

	for(i=0;i<NI;i++)
	{
		MASK[i]=0;
		NENUMS[i]=-1;
		TMP[i]=-1;
		DEG[i]=-1;
	}

	PERM[0] = stnode;
	MASK[stnode] = 1;
	nenums=0;
//определение прямого упорядочения
	for (i=0; i<NI-1; i++)
	{
		FindNOTNUMBERED(GR,MASK,TMP,&tmp,DEG,PERM[i]);
		for (j=0;j<tmp;j++)
		{
			NENUMS[nenums]=TMP[j];
			MASK[TMP[j]]=1;
			nenums++;
		}

		PERM[i+1]=NENUMS[i];
		if (PERM[i+1] == -1)
		{
			nbreaks++;
			//поиск нового узла для продолжения перенумерации из числа непронумерованных
			for (j=0; j<NI; j++)
			{
				if ( MASK[j] != 1 ) break;
			}
			if ( j != NI )
			{
//				fprintf(fp,"RCM Error PERM[i+1]==-1: I_Ur_Se = %d I_Se=%d i+1=%d tmp=%d  newstarting_node = %d NI = %d\n",I_Ur_Se,I_SeL,i+1,tmp,j,NI);
				//sprintf(strl,"RCM Error PERM[i+1]==-1: I_SeL = %d I_Se=%d i+1=%d tmp=%d  newstarting_node = %d NI = %d\n",I_SeL,I_Se,i+1,tmp,j,NI);
				//msg(strl);
				NENUMS[i] = j;
				MASK[j] = 1;
				nenums++;
				PERM[i+1]=NENUMS[i];
			}
			else
			{
//				fprintf(fp,"RCM FULLError PERM[i+1]==-1: I_SeL = %d I_Se=%d i+1=%d tmp=%d  newstarting_node = %d NI = %d\n",I_Ur_Se,I_SeL,i+1,tmp,j,NI);
//				sprintf(strl,"RCM FULLError PERM[i+1]==-1: I_SeL = %d I_Se=%d i+1=%d tmp=%d  newstarting_node = %d NI = %d\n",I_Ur_Se,I_SeL,i+1,tmp,j,NI);
//				msg(strl);
				break;
			}
		}
	}

//перестановка перенумеровывается с конца для повышения эффективности
	k=NI-1;
	i=0;
	while (i<k)
	{
		j=PERM[i];
		PERM[i]=PERM[k];
		PERM[k]=j;
		i++;
		k--;
	}

//расчет обратной перестановки, которая показывает на какие новые номера нужно поставить старые узлы
	for (i=0;i<NI;i++)
	{
		INVP[PERM[i]]=i;
	}

	delete DEG;
	delete TMP;
	delete NENUMS;
	delete MASK;
	
	return(nbreaks);
}

void MEM::FindNOTNUMBERED (int **GR,int *MASK,int *TMP, int *temp,int *DEG,int uz)
{
	int i,j,k,min,i_min,ft,fd,deg;
	int tmp;

	tmp=0;
//vibiraem nenumerovannix sosedei uzla "uz"
	for (i=1; i<=GR[uz][0];i++)
	{
		if (MASK[GR[uz][i]]!=1)
		{
			TMP[tmp]=GR[uz][i];
			DEG[tmp]=GR[TMP[tmp]][0];
			tmp++;
		}
	}

//vibrannoe sortiruem po vozrastaniu stepenei
	k=0;
	i_min=-1;
	while (k<tmp)
	{
		min=100000000;
		for(i=k; i<tmp;i++)
		{
			if(DEG[i]<min)
			{
				min=DEG[i];
				i_min=i;
			}
		}
		ft=TMP[k];
		fd=DEG[k];
		TMP[k]=TMP[i_min];
		DEG[k]=DEG[i_min];
		TMP[i_min]=ft;
		DEG[i_min]=fd;
		k++;
	}
	*temp=tmp;
}

int MEM::FindStartingNode(double *CRD, int NI, int NORT)
{
	int i;
	double xmin, ymin, zmin, mindist, dist;
	xmin = 10000000000000.0;
	ymin = 10000000000000.0;
	zmin = 10000000000000.0;
	for (i=0; i<NI; i++)
	{
		if ( xmin > CRD[i*NORT ] ) xmin = CRD[i*NORT ];
		if ( ymin > CRD[i*NORT +1] ) ymin = CRD[i*NORT +1];
		if ( zmin > CRD[i*NORT +2] ) zmin = CRD[i*NORT +2];
	}
	mindist = 10000000000000000.0;
	int imindist = -1;
	for (i=0; i<NI; i++)
	{
		dist = (CRD[i*NORT ] - xmin)*(CRD[i*NORT ] - xmin) +
			(CRD[i*NORT +1] - ymin)*(CRD[i*NORT +1] - ymin)  +
			(CRD[i*NORT +2] - zmin)*(CRD[i*NORT +2] - zmin);
		if ( dist < mindist )
		{
			mindist = dist;
			imindist = i;
		}
	}
	return(imindist);
}

int MEM::SIZE_PROFILE (int **G, int NN, int NORT)
{
	int i,j,k,r2;
	long long int size;
	double tmp;

	size = 0;

	for(i=0; i<NN; i++)
	{
		if ( G[i][0] > 0 )
		{
			k=i-G[i][1];
		}
		else
		{
			k=0;
		}
		size += k*NORT*NORT;
		size += (int)( (NORT*NORT - NORT)/2 + NORT );
	}
	tmp = (double)size;
	tmp = tmp*8/(1024*1024);
	r2 = (int)tmp;
	return(r2);

}

int MEM::WIDTH_PROFILE(int **G, int NN, int NORT)
{
	int i,k,max;
	max=0;

	for(i=0; i<NN; i++)
	{
		if ( G[i][0] > 0 )
		{
			k=i-G[i][1];
			if (k>max) 
			{
				max=k;
			}
		}
	}
	return(max*NORT);
}

void MEM::RenumberIND (int **IND, int NEL, int NI, int *INVP)
{
	int i,j,k;
	
	for(i=0; i<NEL; i++)
	{
		for(j=2; j < (IND[i][1] + 2); j++)
		{
			if(IND[i][j] < NI)
			{
				IND[i][j] = INVP[ IND[i][j] ];
			}
		}
	}
}

void MEM::RenumberVect (float *VEK, int *INVP, int NI, int NORT)
{
	int i,j;
	float *A;

	A = new float[NI*NORT];

	for(i=0; i < NI; i++)
	{
		for(j=0; j < NORT; j++)
		{
			A[ NORT*INVP[i] + j ] = VEK[ NORT*i + j ];
		}
	}

	for(i=0; i<NI*NORT; i++)
	{
		VEK[i] = A[i];
	}

	delete [] A;
}

void MEM::RenumberVect (int *VEK, int *INVP, int NI, int NORT)
{
	int i,j;
	int *A;

	A = new int[NI*NORT];

	for(i=0; i < NI; i++)
	{
		for(j=0; j < NORT; j++)
		{
			A[ NORT*INVP[i] + j ] = VEK[ NORT*i + j ];
		}
	}

	for(i=0; i<NI*NORT; i++)
	{
		VEK[i] = A[i];
	}

	delete [] A;
}

void MEM::RenumberVect (double *VEK, int *INVP, int NI, int NORT)
{
	int i,j;
	double *A;

	A = new double[NI*NORT];

	for(i=0; i < NI; i++)
	{
		for(j=0; j < NORT; j++)
		{
			A[ NORT*INVP[i] + j ] = VEK[ NORT*i + j ];
		}
	}

	for(i=0; i<NI*NORT; i++)
	{
		VEK[i] = A[i];
	}

	delete [] A;
}