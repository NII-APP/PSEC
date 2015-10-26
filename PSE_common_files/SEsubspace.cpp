#include "StdAfx.h"
#include "SE.h"

	DWORD WINAPI thread_BlockReducedSTIFmatrix(void *n);
	CRITICAL_SECTION CS_BlockKK, CS_AddKK;

void CSE::CalcSubspaceKMPar()
{
	int i,j,k,iform,jform,ik,jk,tmpi,ithread, curthread;
	double *tmpL,*tmp2,kij,tmp;
	char strl[256];
	char str[128];
	DWORD dwExitCode;

	KK = MM->MEM_NEW(KK,nvect,nvect);

	struct KKparDAT *dat;
	dat = new struct KKparDAT[nthread];

	EvNotEmptyQueue = CreateEvent(NULL,TRUE,FALSE,NULL);
	ResetEvent(EvNotEmptyQueue);
	inewtrnum = 0;
	Queuesize = nthread+2;
	ThreadQueue = new int[Queuesize];

	InitializeCriticalSection(&CS_BlockKK);
	InitializeCriticalSection(&CS_AddKK);

	for (i=0; i<nthread; i++)
	{
		dat[i].ithread = i;
		dat[i].infcycleflag = 1;
		dat[i].pcse = this;
		dat[i].startiur = 0;
		dat[i].endiur = 0;

		dat[i].start = CreateEvent(NULL,TRUE,FALSE,NULL);
		ResetEvent(dat[i].start);
		dat[i].totalfin = CreateEvent(NULL,TRUE,FALSE,NULL);
		ResetEvent(dat[i].totalfin);

		dat[i].potok = CreateThread(NULL,0,thread_BlockReducedSTIFmatrix,(void*)&dat[i],0,&dat[i].uThreadID);
	}




	i = 0;
	ithread = 0;
	printf("\n");
	while ( i < NNE )
	{
		tmp = 100.0*(((double)i)/((double)NNE));
		tmpi = (int)tmp;
		printf( "SubspStif %d%%\r",tmpi);

		//проверяем не пуста ли очередь потоков
		WaitForSingleObject(EvNotEmptyQueue,INFINITE);

		//определяем, какому потоку дать задание
		EnterCriticalSection(&CS_BlockKK);

		curthread = ThreadQueue[ithread];
		ithread ++;
		if (ithread == Queuesize) ithread = 0;
		if (ithread != inewtrnum) SetEvent(EvNotEmptyQueue);
		else ResetEvent(EvNotEmptyQueue);

		LeaveCriticalSection(&CS_BlockKK);

		dat[curthread].startiur = i;
		dat[curthread].endiur = i + blocksize;
		
		if ( dat[curthread].endiur > NNE ) dat[curthread].endiur = NNE;

		i = dat[curthread].endiur;

		SetEvent(dat[curthread].start);
	}

	//завершение потоков
	for (i=0; i<nthread; i++)
	{
		//проверяем не пуста ли очередь потоков
		WaitForSingleObject(EvNotEmptyQueue,INFINITE);

		//определяем, какому потоку дать задание
		EnterCriticalSection(&CS_BlockKK);

		curthread = ThreadQueue[ithread];
		ithread ++;
		if (ithread == Queuesize) ithread = 0;
		if (ithread != inewtrnum) SetEvent(EvNotEmptyQueue);
		else ResetEvent(EvNotEmptyQueue);

		LeaveCriticalSection(&CS_BlockKK);

		//завершение потока
		dat[curthread].infcycleflag = 0;
		ResetEvent(dat[curthread].totalfin);
		SetEvent(dat[curthread].start);
		WaitForSingleObject(dat[curthread].totalfin,INFINITE);
	}

	DeleteCriticalSection(&CS_BlockKK);
	DeleteCriticalSection(&CS_AddKK);

	for (i=0; i<nthread; i++)
	{
		CloseHandle(dat[i].potok);
		CloseHandle(dat[i].start);
		CloseHandle(dat[i].totalfin);
	}

	CloseHandle(EvNotEmptyQueue);
	delete [] ThreadQueue;
	delete [] dat;

}

DWORD WINAPI thread_BlockReducedSTIFmatrix(void *n)
{

	int i,j,k,iform,jform,ik,jk,tmpi,nform;
	double *tmpL,*tmp2,**KK,kij,tmp;
	char strl[256];
	int **STII_ENVNZ;
	double  **STIINZ,**FORSE; 

	
	struct KKparDAT *p;
	p = (struct KKparDAT *) n;

	nform = p->pcse->nvect;
	STII_ENVNZ = p->pcse->STII_ENVNZ;
	STIINZ = p->pcse->STIINZ;
	FORSE = p->pcse->FORSE;

	tmpL = new double [nform];
	tmp2 = new double [nform];

	KK = new double *[nform];
	for (i=0; i<nform; i++)
	{
		KK[i] = new double [nform];
	}

	while (p->infcycleflag == 1)
	{
		EnterCriticalSection(&CS_BlockKK);

		p->pcse->ThreadQueue[p->pcse->inewtrnum] = p->ithread;
		SetEvent(p->pcse->EvNotEmptyQueue);
		p->pcse->inewtrnum++;
		if (p->pcse->inewtrnum == p->pcse->Queuesize) p->pcse->inewtrnum = 0;

		LeaveCriticalSection(&CS_BlockKK);

		WaitForSingleObject(p->start, INFINITE);
		ResetEvent(p->start);

		if (p->infcycleflag == 1)
		{
			
			for (i=0; i<nform; i++)
			{
				for (j=0; j<nform; j++)
				{
					KK[i][j] = 0.0;
				}
			}
			for (i=p->startiur; i<p->endiur; i++)
			{
				for (j=0; j<STII_ENVNZ[i][0]; j++)
				{
					kij = STIINZ[i][j];// элемент
					ik = STII_ENVNZ[i][j+1];// строка
					jk = i;// столбец
					if (ik != jk){
						for (iform=0; iform<nform; iform++){ //умножение слева на UF
							tmpL[iform] = FORSE[iform][ik]*kij; //от элемента в верхнем треугольнике
							tmp2[iform] = FORSE[iform][jk]*kij; //от симметричного элемента в нижнем треугольнике
						}

						for (iform=0; iform<nform; iform++){
							for(jform=0; jform<nform; jform++){
								KK[iform][jform] += tmpL[iform]*FORSE[jform][jk];
								KK[iform][jform] += tmp2[iform]*FORSE[jform][ik];
							}
						}
					}
					else{
						for (iform=0; iform<nform; iform++){ //умножение слева на UF
							tmpL[iform] = FORSE[iform][ik]*kij; //от элемента в верхнем треугольнике
						}
						for (iform=0; iform<nform; iform++){
							for(jform=0; jform<nform; jform++){
								KK[iform][jform] += tmpL[iform]*FORSE[jform][jk];
							}
						}
					}
				}
			}

			EnterCriticalSection(&CS_AddKK);

			for (i=0; i<nform; i++)
			{
				for (j=0; j<nform; j++)
				{
					p->pcse->KK[i][j] += KK[i][j];
				}
			}

			LeaveCriticalSection(&CS_AddKK);

		}
	}

	delete [] tmpL;
	delete [] tmp2;

	for (i=0; i<nform; i++)
	{
		delete [] KK[i];
	}
	delete [] KK;

	SetEvent(p->totalfin);

	return(0);

}

void CSE::CalcSubspaceKM()
{
	int i,j,k,nform,iform,jform,ik,jk,tmpi;
	double *tmpL,*tmp2,kij,tmp;
	char strl[256];

	tmpL = new double [nvect];
	tmp2 = new double [nvect];

	KK = MM->MEM_NEW(KK,nvect,nvect);


	for (i=0; i<NNE; i++)
	{
		/*if ( (i%100) == 0 )
		{
			tmp = 100.0*(((double)i)/((double)NNE));
			tmpi = (int)tmp;
			sprintf( strl," ITL %d KKcalc %d%% IEL %d",ITL,tmpi,i);
			statusbar( strl );
		}*/

		for (j=0; j<STII_ENVNZ[i][0]; j++)
		{
			kij = STIINZ[i][j];// элемент
			ik = STII_ENVNZ[i][j+1];// строка
			jk = i;// столбец
			if (ik != jk){
				for (iform=0; iform<nvect; iform++){ //умножение слева на UF
					tmpL[iform] = FORSE[iform][ik]*kij; //от элемента в верхнем треугольнике
					tmp2[iform] = FORSE[iform][jk]*kij; //от симметричного элемента в нижнем треугольнике
				}

				for (iform=0; iform<nvect; iform++){
					for(jform=0; jform<nvect; jform++){
						KK[iform][jform] += tmpL[iform]*FORSE[jform][jk];
						KK[iform][jform] += tmp2[iform]*FORSE[jform][ik];
					}
				}
			}
			else{
				for (iform=0; iform<nvect; iform++){ //умножение слева на UF
					tmpL[iform] = FORSE[iform][ik]*kij; //от элемента в верхнем треугольнике
				}
				for (iform=0; iform<nvect; iform++){
					for(jform=0; jform<nvect; jform++){
						KK[iform][jform] += tmpL[iform]*FORSE[jform][jk];
					}
				}
			}
		}
	}


	delete [] tmpL;
	delete [] tmp2;


}