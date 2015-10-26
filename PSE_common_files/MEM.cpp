#include "stdafx.h"
#include "MEM.h"

MEM::MEM()
{
	curmem = 0.0;
}

MEM::~MEM()
{
}

int* MEM::MEM_NEW(int *B, int k)
{
	int *A;
	if (B != NULL)
	{
		A = MEM_DEL(B,k);
	}

	A=new int [k];
	QQZERO (A,k);

	curmem += 4*((double)k)/(1024*1024);

	return(A);
}

long long int* MEM::MEM_NEW(long long int *B, int k)
{
	long long int *A;
	if (B != NULL)
	{
		A = MEM_DEL(B,k);
	}

	A = new long long int [k];
	QQZERO (A,k);

	curmem += 8*((double)k)/(1024*1024);

	return(A);
}

float* MEM::MEM_NEW(float *B, int k)
{
	float *A;
	if (B != NULL)
	{
		A = MEM_DEL(B,k);
	}

	A=new float [k];
	QQZERO (A,k);

	curmem += 4*((double)k)/(1024*1024);

	return(A);
}

double* MEM::MEM_NEW(double *B,int k)
{
	double *A;
	if (B != NULL)
	{
		A = MEM_DEL(B,k);
	}
	
	A=new double [k];
	QQZERO (A,k);

	curmem += 8*((double)k)/(1024*1024);

	return(A);
}

double* MEM::MEM_NEW(double *B,long long int k)
{
	double *A;
	if (B != NULL)
	{
		A = MEM_DEL(B,k);
	}
	
	A=new double [k];
	QQZERO (A,k);

	curmem += 8*((double)k)/(1024*1024);

	return(A);
}

int** MEM::MEM_NEW(int **B, int str, int stb)
{
	int i;
	int **A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,str,stb);
	}

	A= new int *[str];

    for (i=0; i < str; i++) 
	{
		A[i]=new int [stb];
    }
	QQZERO(A,str,stb);

	curmem += 4*((double)str)*((double)stb)/(1024*1024);

	return(A);
}

float** MEM::MEM_NEW(float **B, int str, int stb)
{
	int i;
	float **A;
	if (B != NULL)
	{
		A = MEM_DEL(B,str,stb);
	}
	A= new float *[str];
    for (i=0; i < str; i++) 
	{
		A[i]=new float [stb];
    }
	QQZERO(A,str,stb);

	curmem += 4*((double)str)*((double)stb)/(1024*1024);

	return(A);
}

double** MEM::MEM_NEW(double **B, int str, int stb)
{
	int i;
	double **A;
	if (B != NULL)
	{
		A = MEM_DEL(B,str,stb);
	}
	A= new double *[str];
    for (i=0; i < str; i++) 
	{
		A[i]=new double [stb];
    }
	QQZERO(A,str,stb);

	curmem += 8*((double)str)*((double)stb)/(1024*1024);

	return(A);
}

int*** MEM::MEM_NEW(int ***B, int n1, int n2, int n3)
{
	int i,j;
	int ***A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,n1,n2,n3);
	}

	A= new int **[n1];
    for (i=0; i < n1; i++) 
	{
		A[i]=new int *[n2];
		for (j=0; j<n2; j++)
		{
			A[i][j]=new int [n3];
		}
		QQZERO(A[i],n2,n3);
	}

	curmem += 4*((double)n1)*((double)n2)*((double)n3)/(1024*1024);

	return(A);
}

float*** MEM::MEM_NEW(float ***B, int n1, int n2, int n3)
{
	int i,j;
	float ***A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,n1,n2,n3);
	}

	A= new float **[n1];
    for (i=0; i < n1; i++) 
	{
		A[i]=new float *[n2];
		for (j=0; j<n2; j++)
		{
			A[i][j]=new float [n3];
		}
		QQZERO(A[i],n2,n3);
	}

	curmem += 4*((double)n1)*((double)n2)*((double)n3)/(1024*1024);

	return(A);
}

double*** MEM::MEM_NEW(double ***B, int n1, int n2, int n3)
{
	int i,j;
	double ***A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,n1,n2,n3);
	}

	A= new double **[n1];
    for (i=0; i < n1; i++) 
	{
		A[i]=new double *[n2];
		for (j=0; j<n2; j++)
		{
			A[i][j]=new double [n3];
		}
		QQZERO(A[i],n2,n3);
	}

	curmem += 8*((double)n1)*((double)n2)*((double)n3)/(1024*1024);

	return(A);
}

float**** MEM::MEM_NEW(float ****B, int n1, int n2, int n3, int n4)
{
	int i,j,k;
	float ****A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,n1,n2,n3,n4);
	}

	A= new float ***[n1];
    for (i=0; i < n1; i++) 
	{
		A[i]=new float **[n2];
		for (j=0; j<n2; j++)
		{
			A[i][j]=new float *[n3];
			for (k=0; k<n3; k++)
			{
				A[i][j][k] = new float [n4];
				QQZERO(A[i][j][k],n4);
			}
		}
	}

	curmem += 4*((double)n1)*((double)n2)*((double)n3)*((double)n4)/(1024*1024);

	return(A);
}

double**** MEM::MEM_NEW(double ****B, int n1, int n2, int n3, int n4)
{
	int i,j,k;
	double ****A;
	
	if (B != NULL)
	{
		A = MEM_DEL(B,n1,n2,n3,n4);
	}

	A= new double ***[n1];
    for (i=0; i < n1; i++) 
	{
		A[i]=new double **[n2];
		for (j=0; j<n2; j++)
		{
			A[i][j]=new double *[n3];
			for (k=0; k<n3; k++)
			{
				A[i][j][k] = new double [n4];
				QQZERO(A[i][j][k],n4);
			}
		}
	}

	curmem += 8*((double)n1)*((double)n2)*((double)n3)*((double)n4)/(1024*1024);

	return(A);
}



int* MEM::MEM_DEL(int *A, int k)
{
	if (A != NULL) 
	{
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)k)/(1024*1024);
	return(A);
}

long long int* MEM::MEM_DEL(long long int *A, int k)
{
	if (A != NULL) 
	{
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)k)/(1024*1024);
	return(A);
}

float* MEM::MEM_DEL(float *A, int k)
{
	if (A != NULL) 
	{
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)k)/(1024*1024);
	return(A);
}

double* MEM::MEM_DEL(double *A, int k)
{
	if (A != NULL) 
	{
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)k)/(1024*1024);
	return(A);
}

double* MEM::MEM_DEL(double *A, long long int k)
{
	if (A != NULL) 
	{
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)k)/(1024*1024);
	return(A);
}

int** MEM::MEM_DEL(int **A, int str, int stb)
{
	int i;
	if (A != NULL) 
	{
		for (i=str-1; i>=0; i--) 
	    {
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)str)*((double)stb)/(1024*1024);
	return(A);
}

float** MEM::MEM_DEL(float **A, int str, int stb)
{
	int i;
	if (A != NULL) 
	{
		for (i=str-1; i>=0; i--) 
	    {
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)str)*((double)stb)/(1024*1024);
	return(A);
}

double** MEM::MEM_DEL(double **A, int str, int stb)
{
	int i;
	if (A != NULL) 
	{
		for (i=str-1; i>=0; i--) 
	    {
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)str)*((double)stb)/(1024*1024);
	return(A);
}

int*** MEM::MEM_DEL(int ***A, int n1, int n2,int n3)
{
	int i,j;
	if (A != NULL) 
	{
		for (i=n1-1; i>=0; i--) 
	    {
			for (j=n2-1; j>=0; j--)
			{
				if (A[i][j] != NULL) delete [] A[i][j];
			}
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)n1)*((double)n2)*((double)n3)/(1024*1024);
	return(A);
}

float*** MEM::MEM_DEL(float ***A, int n1, int n2,int n3)
{
	int i,j;
	if (A != NULL) 
	{
		for (i=n1-1; i>=0; i--) 
	    {
			for (j=n2-1; j>=0; j--)
			{
				if (A[i][j] != NULL) delete [] A[i][j];
			}
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)n1)*((double)n2)*((double)n3)/(1024*1024);
	return(A);
}

double*** MEM::MEM_DEL(double ***A, int n1, int n2,int n3)
{
	int i,j;
	if (A != NULL) 
	{
		for (i=n1-1; i>=0; i--) 
	    {
			for (j=n2-1; j>=0; j--)
			{
				if (A[i][j] != NULL) delete [] A[i][j];
			}
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)n1)*((double)n2)*((double)n3)/(1024*1024);
	return(A);
}

float**** MEM::MEM_DEL(float ****A, int n1, int n2, int n3, int n4)
{
	int i,j,k;
	if (A != NULL) 
	{
		for (i=n1-1; i>=0; i--) 
	    {
			for (j=n2-1; j>=0; j--)
			{
				for (k=n3-1; k>=0; k--)
				{
					if (A[i][j][k] != NULL) delete [] A[i][j][k];
				}
				if (A[i][j] != NULL) delete [] A[i][j];
			}
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 4*((double)n1)*((double)n2)*((double)n3)*((double)n4)/(1024*1024);
	return(A);
}

double**** MEM::MEM_DEL(double ****A, int n1, int n2, int n3, int n4)
{
	int i,j,k;
	if (A != NULL) 
	{
		for (i=n1-1; i>=0; i--) 
	    {
			for (j=n2-1; j>=0; j--)
			{
				for (k=n3-1; k>=0; k--)
				{
					if (A[i][j][k] != NULL) delete [] A[i][j][k];
				}
				if (A[i][j] != NULL) delete [] A[i][j];
			}
			if (A[i] != NULL) delete [] A[i];
		}
		delete [] A;
		A = NULL;
	}
	curmem -= 8*((double)n1)*((double)n2)*((double)n3)*((double)n4)/(1024*1024);
	return(A);
}


int MEM::QQZERO (int *A, int k)
{
	int i;
	for(i=0; i<k; i++) A[i]=0;
	return(0);
}
int MEM::QQZERO (long long int *A, int k)
{
	int i;
	for(i=0; i<k; i++) A[i]=0;
	return(0);
}

int MEM::QQZERO (float *A, int k)
{
	int i;
	for(i=0; i<k; i++) A[i]=0.0;
	return(0);
}

int MEM::QQZERO (double *A, int k)
{
	int i;
	for(i=0; i<k; i++) A[i]=0.0;
	return(0);
}

int MEM::QQZERO (double *A, long long int k)
{
	long long int i;
	for(i=0; i<k; i++) A[i]=0.0;
	return(0);
}

int MEM::QQZERO (int **A, int k, int n)
{
	int i,j;

	for(i=0; i<k; i++)
	{
		for(j=0; j<n; j++)
		{
			A[i][j]=0;
		}
	}
	return(0);
}

int MEM::QQZERO (float **A, int k, int n)
{
	int i,j;

	for(i=0; i<k; i++)
	{
		for(j=0; j<n; j++)
		{
			A[i][j]=0.0;
		}
	}
	return(0);
}

int MEM::QQZERO (double **A, int k, int n)
{
	int i,j;

	for(i=0; i<k; i++)
	{
		for(j=0; j<n; j++)
		{
			A[i][j]=0.0;
		}
	}
	return(0);
}

