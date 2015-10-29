#include "StdAfx.h"
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <process.h>
//#include "..\uztype.h"

extern clock_t BegClock, EndClock;
//extern double msgclock();

int msg(char *str);
//void statusbar(char *strl );
void xolp( double *A, long long int *D, int k);
void xolsinglesolv( double *L, long long int *D, int k, double *B );
//void xolmultisolv_1th( double *L, long long int *D, int k, double **B, int nv );

DWORD WINAPI element(void *n);

#define kolpotok 20
#define kol 100

struct DAT {
	double *A;
	long long int *D;
	int i;
	int j;
	int sk;	
};
HANDLE Fin[kolpotok],Start[kolpotok];
HANDLE potok_e[kolpotok];
CRITICAL_SECTION CS[kolpotok];
struct DAT *dat;
int zahod=0;

DWORD WINAPI element(void *n)
{
	long long int i,f,j,ip,nnn,curp;
	long long int zmin,z;
	double s;
	double *A,*pst1,*pst2;
	long long int *D;
	int *nn;

	f=1;
	struct DAT *p;
    nn=(int *) n;
	i=1;
	curp = *nn;

	while (i >= 0)	{
		WaitForSingleObject(Start[*nn],INFINITE);
		p=&dat[*nn];
		D=p->D;
		A=p->A;
		j = p->j;
		if(p->i >= 0){	//для реализации механизма завершения потока
			for (i= j+1 - (D[j]-D[j-1]); i<=j; i++){
				//проверка готов ли к использованию столбец i
				if ( (j-i) < kolpotok ){  //!!//есть вероятность, что разложение столбца j-1 не закончено
					ip = curp - 1;
					if (ip < 0) ip = ip + kolpotok;
					if ( ip != curp ){	// текущий поток проверять не нужно
						EnterCriticalSection(&CS[ip]);
						if ( dat[ip].j < j ){
							WaitForSingleObject(Fin[ip],INFINITE);
							SetEvent(Fin[ip]);	//требуется повторная установка события, поскольку WaitForSingleObject сбросила его после проверки
						}							
						LeaveCriticalSection(&CS[ip]);
					}
				}
				s=0;
				if ( i > 0 ) zmin=((D[j]-D[j-1]-1-j+i)<(D[i]-D[i-1]-1)?(D[j]-D[j-1]-1-j+i):(D[i]-D[i-1]-1));
				else zmin = 0;

				if (i!=j){ 
					pst1 = &A[D[i]];
					pst1--;
					pst2 = &A[D[j]-j+i-1];
					for(z=1; z<=zmin; ++z) s += *pst1--**pst2--;
				}
				else
				{
					pst1 = &A[D[i]];
					pst1--;
					pst2 = pst1;
					for (z=1; z<=zmin; ++z)
					{
						*pst1 = *pst1/A[D[i-z]];
						s += *pst1--**pst2--*A[D[i-z]];
					}
				}
				A[ D[j] - j + i ] -= s; 
			}
		}
		else{
			i = p->i;
		}
		ResetEvent(Start[*nn]);
		SetEvent(Fin[*nn]);
	}
	return(0);
}

extern void xolp( double *A, long long int *D, int k )
{
	int i,j,ii,h,h1,p,f,prov,np,pp,dlina;
	int nn[kolpotok];
	char str[128];
	char	strl[256];
	DWORD uThreadIDs[kolpotok];
	DWORD dwExitCode;
	for (i=0;i<kolpotok;i++){
		sprintf(str,"%d",i);
		Fin[i]=CreateEvent(NULL,TRUE,TRUE,NULL);
	}
			
	for (i=0;i<kolpotok;i++){
		sprintf(str,"000%d",i*2);
		Start[i]=CreateEvent(NULL,TRUE,FALSE,NULL);
		ResetEvent(Start[i]);
	}

	for (i=0;i<kolpotok;i++){
		InitializeCriticalSection(&CS[i]);
	}

	zahod++;
	dat=new struct DAT [kolpotok];	
	
	for(i=0;i<kolpotok;i++){
		dat[i].A=A;
		dat[i].D=D;
		nn[i]=i;
	}
	i=0;

	for (i=0;i<kolpotok;i++){
		potok_e[i] = CreateThread(NULL,0,element,(void*)&nn[i],0,&uThreadIDs[i]);
	}
	//первый элемент обрабатывается отдельно
	//A[0]=sqrt(A[0]);
	//основной алгоритм распределения потоков
	p=0;
	for (i=1; i<k; i++){
		if(((i+1)%1000)==0){
//			sprintf(strl," xolp %5.2f percent",float(i*100/k));
			printf(" xolp %5.2f percent A[%d] = %e\r",float(i*100/k),i-2*kolpotok,A[D[i-2*kolpotok]]);
//			statusbar( strl );
		}
		EnterCriticalSection(&CS[p]);
		WaitForSingleObject(Fin[p],INFINITE);
		
		dat[p].i = i;
		dat[p].j = i;
			
		ResetEvent(Fin[p]);
		SetEvent(Start[p]);

		LeaveCriticalSection(&CS[p]);
		p++;
		if ( p == kolpotok ) p = 0;
	}

	for(p=0; p<kolpotok; p++){
		WaitForSingleObject(Fin[p],INFINITE);
	}
	for (p=0;p<kolpotok;p++){
		dat[p].i = -1;
		dat[p].j = 0;
		dat[p].sk = 0;
		ResetEvent(Fin[p]);
		SetEvent(Start[p]);
	}
	for(p=0; p<kolpotok; p++){
		WaitForSingleObject(Fin[p],INFINITE);
	}

	delete [] dat;
	for( i = kolpotok-1; i >= 0; i-- ){
		CloseHandle(potok_e[i]);
	}
	for( i = kolpotok-1; i >= 0; i-- ){
		DeleteCriticalSection(&CS[i]);
	}
	for( i = kolpotok-1; i >= 0; i-- ){
		CloseHandle(Start[i]);
		CloseHandle(Fin[i]);
	}
	printf("\n");
}


void xolsinglesolv( double *L, long long int *D, int k, double *B )
{
	long long int  i,j,jmax;
	double bi,s,*pL,*pB;

	for(i=1; i<k; i++)
	{
		s=0;
		pL = &L[D[i-1]+1];
		pB = &B[i+1-D[i]+D[i-1]];
		jmax = D[i] - D[i-1] - 1;
		for (j=0; j<jmax; j++) s += *pL++**pB++;
		B[i] -= s;
	}
	for (i=0; i<k; i++)
	{
		B[i] = B[i]/L[D[i]];
	}

	for(i=k-1; i>=0; i--)
	{
		bi = B[i];
		pB = &B[i];
		pL = &L[D[i]];
		if (i>0)
		{
			pL--;
			pB--;
			jmax = D[i] - D[i-1] - 1;
		}
		else jmax = 0;
		for (j=0; j<jmax; j++) *pB-- -= *pL--*bi;
	}	

}


//void xolmultisolv_1th( double *L, long long int *D, int k, double **B, int nv )
//{
//	long long int  i,j,jmax;
//	double bi,s,*pL,*pB;
//
//	for(i=1; i<k; i++)
//	{
//		s=0;
//		pL = &L[D[i-1]+1];
//		pB = &B[i+1-D[i]+D[i-1]];
//		jmax = D[i] - D[i-1] - 1;
//		for (j=0; j<jmax; j++) s += *pL++**pB++;
//		B[i] -= s;
//	}
//	for (i=0; i<k; i++)
//	{
//		B[i] = B[i]/L[D[i]];
//	}
//
//	for(i=k-1; i>=0; i--)
//	{
//		bi = B[i];
//		pB = &B[i];
//		pL = &L[D[i]];
//		if (i>0)
//		{
//			pL--;
//			pB--;
//			jmax = D[i] - D[i-1] - 1;
//		}
//		else jmax = 0;
//		for (j=0; j<jmax; j++) *pB-- -= *pL--*bi;
//	}	
//
//}

