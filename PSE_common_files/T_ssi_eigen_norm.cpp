#include "StdAfx.h"
#include "T_ssi_eigen.h"
#include "math.h"
#include "stdio.h"

void SSIEIGEN::NormalizeEigVect()
{
	int i,j,k;

	double s, ss;

	for (i=0; i<nform; i++){
		//нормировка
		s=0;
		for (k=0; k<nur; k++){
			s += UF[i][k]*M[k]*UF[i][k];
		}
		s = sqrt(s);
		for (k=0; k<nur; k++){
			UF[i][k] = (UF[i][k]/s);
		}
	}
}

void SSIEIGEN::CollinearDetection()
{
	int i,j,k;

	double s, ss;

	ss=0;
	s=0;
	fprintf(fp3,"\n\n");
	for (i=0; i<nform; i++){
		for (j=i; j<nform; j++){
			s=0;
			for (k=0; k<nur; k++){
				s += UF[i][k]*M[k]*UF[j][k];
			}

			if ( (j!=i) && (abs(s-1.0) < 0.000001) )
			{//исключаем коллинеарный вектор
				fprintf(fp,"\nCOLLINEAR DETECTED,  ..... i=%d  j=%d     s = %e\n",i,j,s);
				//srand((unsigned)time( NULL ) );
				//s = 0.0;
				//for (k=0; k<nur; k++)
				//{
				//	UF[j][k] = (((double) rand())/((double) RAND_MAX)) - 0.5;
				//	if ( M[k] < 1.e-18 ) UF[j][k] = 0.0;
				//	s += UF[j][k]*M[k]*UF[j][k];
				//}
				//s = sqrt(s);
				//for (k=0; k<nur; k++)	UF[j][k] = UF[j][k]/s; 

				////проверка
				//s=0;
				//for (k=0; k<nur; k++){
				//	s += UF[i][k]*M[k]*UF[j][k];
				//}
				
			}

			fprintf(fp3,"ITL = %d  Gramm_Shmidt_Check_main 0  ..... i=%d  j=%d     s = %e\n",ITL,i,j,s);
			fflush(fp3);

		}
	}
}

void SSIEIGEN::UFGetMax()
{
	int i, ii, j, jj, k;
	double displ=0.0;

	for (i=0; i<nform; i++)
	{
		UFmax[i] = 0.0;
		
		for (j=0; j<nur; j++)
		{
			if (UFmax[i] < abs(UF[i][j])) UFmax[i] = abs(UF[i][j]);
		}
	}
}

void SSIEIGEN::FormSort( )
{
	int i,j,k,imax;
	double m;
	int *tmpL;

	tmpL = new int [nform];
	for (i=0; i<nform; i++) tmpL[i] = 0;
	for (i=0; i<nform; i++){
		m = -1.e40;
		for (j=0; j<nform; j++)	{
			if ( m <= LM[j] && tmpL[j] == 0 ) {
				m = LM[j];
				imax = j;
			}
		}
		tmpL[imax] = 1;
		tmpi[nform - 1 - i] = imax;
	}
	delete [] tmpL;
}

void SSIEIGEN::GrammShmidtCyckl()
{
	int i,j,i_ortog,k;
	double eps_ort, err_ort,s,err_orti,sij;
	char strl[256];
	
	eps_ort = 1.e-13;
	i_ortog = 0;
	err_ort = 1;
	err_orti = 1;

	//fprintf(fp3,"\n\n\nMASS\n\n");
	//for (i=0; i<nur; i++)
	//{
	//	fprintf(fp3,"M[%d] = %e\n",i,M[i]);
	//}
	

	sprintf( strl," GrammShmidtCyckl ... strated ...");
	statusbar( strl );

	fprintf(fp3,"\n\n\nGramm_Shmidt_Cyckl .....started\n\n");
	while ( err_orti > eps_ort && i_ortog < 2 )
	{
		i_ortog++;
		GaussShmidt_2();

		s=0;
		err_orti = 0;
		for (i=0; i<nform; i++)
		{
			err_ort = 0;
			for (j=0; j<nform; j++)
			{
				s=0;
				for (k=0; k<nur; k++)
				{
					s += UF[i][k]*M[k]*UF[j][k];
				}
				fprintf(fp3,"Gramm_Shmidt_Cyckl IT %d..... i=%d  j=%d     s = %e\n",i_ortog,i,j,s);
				
				if ( (i!=j) && (abs(s) > err_ort) ) err_ort = abs(s);
				if ( i == j ) sij = abs(s);
			}
			s = err_ort/sij;
			if ( s > err_orti ) err_orti = s;
		}
		fprintf(fp3,"Gramm_Shmidt_Cyckl IT %d..... errmax = %e\n",i_ortog,err_orti);
		fflush(fp3);
		sprintf( strl," GSCyckl ...IT %d ...err %e",i_ortog,err_orti);
		statusbar( strl );
	}
	sprintf( strl," GrammShmidtCyckl ... finished ...");
	statusbar( strl );
}

void SSIEIGEN::GaussShmidt_2()
{
	int i,j,k,nprevort,nfcur;

	double *tmpL, s, ss;
	tmpL = new double[nur];

	////!!!!!!!!!!!!
	//if (ITL == 30)
	//{
	//	FILE *fpp;
	//	fpp = fopen("orto_1.bin","wb");
	//	int nfsum;
	//	nfsum = nfconv + nform;
	//	fwrite(&nfsum,sizeof(int),1,fpp);
	//	fwrite(&nur,sizeof(int),1,fpp);
	//	for (i=0; i<nfconv; i++)
	//	{
	//		fwrite(UFconv[i],sizeof(double),nur,fpp);
	//	}
	//	for (i=0; i<nform; i++)
	//	{
	//		fwrite(UF[i],sizeof(double),nur,fpp);
	//	}
	//	fclose(fpp);
	//}


	for (i=0; i<nform; i++){
		//ортогонализация
		for (k=0; k<nur; k++) tmpL[k] = UF[i][k];	
		//nprevort = i+LastConvNum+1;
		//nprevort = i+nfconv;
		nprevort = nfconv;
		for (j=0; j<nprevort; j++){
			s=0;
			ss=0;
			//if (j<(LastConvNum+1))
			if (j<(nfconv))
			{
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UFconv[j][k];
					ss += UFconv[j][k]*M[k]*UFconv[j][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UFconv[j][k];
				}
			}
			else
			{
				//nfcur = j - (LastConvNum+1);
				nfcur = j - nfconv;
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UF[nfcur][k];
					ss += UF[nfcur][k]*M[k]*UF[nfcur][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UF[nfcur][k];
				}
			}
		}
		//нормировка
		s=0;
		for (k=0; k<nur; k++){
			s += tmpL[k]*M[k]*tmpL[k];
		}
		if ( s <= 0.0 )
		{
			fprintf(fp,"\n\n ORTOGONALIZATION   i = %d s = %e\n",i,s);
			for (j=0; j<nur; j++)
			{
				fprintf(fp, "UF[%d][%d] = %e\t\t\tM = %e\n",i,j,tmpL[j],M[j]);
			}
		}
		s = sqrt(s);
		for (k=0; k<nur; k++){
//			UF[i][k] = (tmpL[k]/s)*1000.;
			UF[i][k] = (tmpL[k]/s);
		}
	}
	delete [] tmpL;
}


void SSIEIGEN::GrammShmidtCheck()
{
	int i,j,k;
	double s, ss;

	ss=0;
	s=0;
	fprintf(fp,"\n\n");
	for (i=0; i<nform; i++){
		for (j=i; j<nform; j++){
			s=0;
			for (k=0; k<nur; k++){
				s += UF[i][k]*M[k]*UF[j][k];
			}
			fprintf(fp3,"ITL = %d  Gramm_Shmidt_Check_main 1  ..... i=%d  j=%d     s = %e\n",ITL,i,j,s);
			fflush(fp3);
		}
	}
}

void SSIEIGEN::GrammShmidtCyckl_SpecSeqv()
{
	int i,j,i_ortog,k;
	double eps_ort, err_ort,s,err_orti,sij;
	char strl[256];
	
	eps_ort = 1.e-13;
	i_ortog = 0;
	err_ort = 1;
	err_orti = 1;

	sprintf( strl," GrammShmidtCyckl ... strated ...");
	statusbar( strl );

	fprintf(fp3,"\n\n\nGramm_Shmidt_Cyckl .....started\n\n");
	while ( err_orti > eps_ort && i_ortog < 2 )
	{
		i_ortog++;
		GaussShmidt_2_SpecSeqv();

		s=0;
		err_orti = 0;
		for (i=0; i<nform; i++)
		{
			err_ort = 0;
			for (j=0; j<nform; j++)
			{
				s=0;
				for (k=0; k<nur; k++)
				{
					s += UF[i][k]*M[k]*UF[j][k];
				}
				fprintf(fp3,"Gramm_Shmidt_Cyckl IT %d..... i=%d  j=%d     s = %e\n",i_ortog,i,j,s);
				
				if ( (i!=j) && (abs(s) > err_ort) ) err_ort = abs(s);
				if ( i == j ) sij = abs(s);
			}
			s = err_ort/sij;
			if ( s > err_orti ) err_orti = s;
		}
		fprintf(fp3,"Gramm_Shmidt_Cyckl IT %d..... errmax = %e\n",i_ortog,err_orti);
		fflush(fp3);
		sprintf( strl," GSCyckl ...IT %d ...err %e",i_ortog,err_orti);
		statusbar( strl );
	}
	sprintf( strl," GrammShmidtCyckl ... finished ...");
	statusbar( strl );
}

void SSIEIGEN::GaussShmidt_2_SpecSeqv()
{
	int i,j,k,nprevort,nfcur;

	double *tmpL, s, ss;
	tmpL = new double[nur];

	for (i=0; i<nform; i++){
		//ортогонализация
		for (k=0; k<nur; k++) tmpL[k] = UF[ortoseqv[i]][k];	
		//nprevort = i+LastConvNum+1;
		nprevort = i+nfconv;
		for (j=0; j<nprevort; j++){
			s=0;
			ss=0;
			//if (j<(LastConvNum+1))
			if (j<(nfconv))
			{
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UFconv[j][k];
					ss += UFconv[j][k]*M[k]*UFconv[j][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UFconv[j][k];
				}
			}
			else
			{
				//nfcur = j - (LastConvNum+1);
				nfcur = ortoseqv[j - nfconv];
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UF[nfcur][k];
					ss += UF[nfcur][k]*M[k]*UF[nfcur][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UF[nfcur][k];
				}
			}
		}
		//нормировка
		s=0;
		for (k=0; k<nur; k++){
			s += tmpL[k]*M[k]*tmpL[k];
		}
		s = sqrt(s);
		for (k=0; k<nur; k++){
//			UF[i][k] = (tmpL[k]/s)*1000.;
//			UF[i][k] = (tmpL[k]/s);
			UF[ortoseqv[i]][k] = (tmpL[k]/s);
		}
	}
	delete [] tmpL;
}

void SSIEIGEN::GaussShmidt_3()
{
	int i,j,k,nprevort,nfcur;

	double *tmpL, s, ss;
	tmpL = new double[nur];

	for (i=0; i<nform; i++){
		//ортогонализация
		for (k=0; k<nur; k++) tmpL[k] = UF[i][k];	
		//nprevort = i+LastConvNum+1;
		//nprevort = i+nfconv;
		nprevort = i;
		for (j=0; j<nprevort; j++){
			s=0;
			ss=0;
			//if (j<(LastConvNum+1))
			if (j<(nfconv))
			{
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UFconv[j][k];
					ss += UFconv[j][k]*M[k]*UFconv[j][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UFconv[j][k];
				}
			}
			else
			{
				//nfcur = j - (LastConvNum+1);
				nfcur = j - nfconv;
				for (k=0; k<nur; k++){
	//				s += UF[i][k]*UF[j][k];
					s += tmpL[k]*M[k]*UF[nfcur][k];
					ss += UF[nfcur][k]*M[k]*UF[nfcur][k];
				}
				s = s/ss;
				for (k=0; k<nur; k++){
					tmpL[k] -= s*UF[nfcur][k];
				}
			}
		}
		//нормировка
		s=0;
		for (k=0; k<nur; k++){
			s += tmpL[k]*M[k]*tmpL[k];
		}
		if ( s <= 0.0 )
		{
			fprintf(fp,"\n\n ORTOGONALIZATION   i = %d s = %e\n",i,s);
			for (j=0; j<nur; j++)
			{
				fprintf(fp, "UF[%d][%d] = %e\t\t\tM = %e\n",i,j,tmpL[j],M[j]);
			}
		}
		s = sqrt(s);
		for (k=0; k<nur; k++){
//			UF[i][k] = (tmpL[k]/s)*1000.;
			UF[i][k] = (tmpL[k]/s);
		}
	}
	delete [] tmpL;
}