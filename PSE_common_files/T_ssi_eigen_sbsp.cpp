#include "StdAfx.h"
#include "T_ssi_eigen.h"
#include "math.h"
#include "stdio.h"


void SSIEIGEN::JacobySweep()
{//реализует одно "прометание" матриц в алгоритме Якоби
	int in,jn,i,j,k,signum;
	double al,gam,kii,kjj,kk,t,x;
	double *tmpki,*tmpkj,*tmpmi,*tmpmj;
	double cfk,cfm;
	double tollerance;

	tollerance = 1.e-7;

	tmpki = new double[nform];
	tmpkj = new double[nform];
	tmpmi = new double[nform];
	tmpmj = new double[nform];

	for (in=0; in<nform-1; in++){
		for (jn=in+1; jn<nform; jn++){
			//поворот блока для in jn
			cfk = sqrt(KK[in][jn]*KK[in][jn]/(KK[in][in]*KK[jn][jn]));
			cfm = sqrt(MM[in][jn]*MM[in][jn]/(MM[in][in]*MM[jn][jn]));
			if ( abs(cfk) > tollerance || abs(cfm) > tollerance ){
				
				//вычисление параметров поворота
				kii = KK[in][in]*MM[in][jn] - MM[in][in]*KK[in][jn];
				kjj = KK[jn][jn]*MM[in][jn] - MM[jn][jn]*KK[in][jn];
				kk = KK[in][in]*MM[jn][jn] - MM[in][in]*KK[jn][jn];
				if ( (kk*kk/4 + kii*kjj) <= 0.0 )
				{
					fprintf(fp,"\n1___ in = %d  jn = %d SWEEP KKii/MMii = %e   KKjj/MMjj = %e    KKij/MMij = %e \n",in,jn,KK[in][in]/MM[in][in],KK[jn][jn]/MM[jn][jn],KK[in][jn]/MM[in][jn]);
					fprintf(fp,"\n1___ in = %d  jn = %d SWEEP KKii = %e   KKjj = %e    KKij = %e \n",in,jn,KK[in][in],KK[jn][jn],KK[in][jn]);
					fprintf(fp,"\n1___ in = %d  jn = %d SWEEP MMii = %e   MMjj = %e    MMij = %e \n",in,jn,MM[in][in],MM[jn][jn],MM[in][jn]);
					al = 0.0;
					gam = 0.0 - KK[in][jn]/KK[jn][jn];
					fprintf(fp,"\n al = %e   gam = %e cfk=%e cfm = %e\n",al,gam,cfk,cfm);
				}
				else
				{
					t = sqrt(kk*kk/4 + kii*kjj);
					signum = (kk >= 0)?(1):(-1);
					x = kk/2 + signum*t;
					gam = - kii/x;
					al = kjj/x;
				}
				//выполнение поворота в строках матриц жесткости и масс
				for (i=0; i<nform; i++){
					tmpki[i] = KK[in][i] + gam*KK[jn][i]; 
					tmpkj[i] = KK[jn][i] + al*KK[in][i];
					tmpmi[i] = MM[in][i] + gam*MM[jn][i]; 
					tmpmj[i] = MM[jn][i] + al*MM[in][i];
				}
				for (i=0; i<nform; i++){
					KK[in][i] = tmpki[i];
					KK[jn][i] = tmpkj[i];
					MM[in][i] = tmpmi[i];
					MM[jn][i] = tmpmj[i];
				}
				//выполнение поворота в столбцах матриц жесткости и масс
				for (i=0; i<nform; i++){
					tmpki[i] = KK[i][in] + gam*KK[i][jn]; 
					tmpkj[i] = KK[i][jn] + al*KK[i][in];
					tmpmi[i] = MM[i][in] + gam*MM[i][jn]; 
					tmpmj[i] = MM[i][jn] + al*MM[i][in];
				}
				for (i=0; i<nform; i++){
					KK[i][in] = tmpki[i];
					KK[i][jn] = tmpkj[i];
					MM[i][in] = tmpmi[i];
					MM[i][jn] = tmpmj[i];
				}
				//накопление собственных форм в матрице QQ - домножение столбцов
				for (i=0; i<nform; i++){
					tmpki[i] = Q[i][in] + gam*Q[i][jn]; 
					tmpkj[i] = Q[i][jn] + al*Q[i][in];
				}
				for (i=0; i<nform; i++){
					Q[i][in] = tmpki[i];
					Q[i][jn] = tmpkj[i];
				}

			}
		}
	}

	delete [] tmpmj;
	delete [] tmpmi;
	delete [] tmpkj;
	delete [] tmpki;
}

void SSIEIGEN::JacobyOFF (double *KKoff, double *MMoff)
{
	int i,j;
	double Koff, Moff, Kdiag, Mdiag;

	Koff = 0;
	Moff = 0;
	Kdiag = 0;
	Mdiag = 0;

	for (i=0; i<nform; i++){
		for (j=i; j<nform; j++){
			if (i == j){
				Kdiag += KK[i][j]*KK[i][j];
				Mdiag += MM[i][j]*MM[i][j];
			}
			else{
				Koff += KK[i][j]*KK[i][j];
				Moff += MM[i][j]*MM[i][j];
			}
//			fprintf(fp,"\ni = %d  j = %d JACOBY  init  KKoff = %e   MMoff = %e    KKdiag = %e   MMdiag = %e \n", i,j,Koff, Moff, Kdiag, Mdiag);
		}
	}
	*KKoff = Koff/Kdiag;
	*MMoff = Moff/Mdiag;
}

void SSIEIGEN::Jacoby()
{
	int in,jn,i,j,k, itersweep;
	double KKoff, MMoff,*LMoldj,errfr,epsfr,err;
	double tollerance;

	LMoldj = new double[nform];

	tollerance = 1.e-14;
	errfr = 1.0;
	epsfr = 1.e-7;

	for (i=0; i<nform; i++){
		LMoldj[i] = 0.0;
		for(j=0; j<nform; j++){
			if ( i == j ) Q[i][j]=1.0;
			else Q[i][j] = 0.0;
		}
	}

	KKoff = 0;
	MMoff = 0;
	JacobyOFF(&KKoff,&MMoff);
	fprintf(fp,"\n\n JACOBY  init  KKoff = %e   MMoff = %e  \n", KKoff, MMoff);
	itersweep = 0;
	while ( ( KKoff > tollerance || MMoff > tollerance || errfr > epsfr ) && itersweep < 25 ){
		JacobySweep();
		JacobyOFF( &KKoff, &MMoff );

		//определение собственных чисел
		for (i=0; i<nform; i++){
			LMoldj[i] = LM[i];
			LM[i] = KK[i][i]/MM[i][i];
	//		LM[i] = KK[i][i]/MM[i][i] + shifting;
		}
				
		if (itersweep > 2)
		{
			errfr = 0.0;
			for (i=0; i<nform; i++)
			{
				err = (abs(LM[i] - LMoldj[i]))/LM[i];
				if ( err > errfr ) errfr = err;
			}
		}

		//вывод KK
		fprintf(fp,"\n\nITsweep=%d   QQ JACOBY   KKoff = %e   MMoff = %e  errfr = %e\n",itersweep, KKoff, MMoff,errfr);
		
		fflush(fp);
		itersweep++;
	}

	delete [] LMoldj;
}

void SSIEIGEN::NewForms()
{//с релаксацией
	int i,j,k;
	double *tmpL;

	tmpL = new double[nform];

	for (i=0; i<nur; i++){
		for(j=0; j<nform; j++){
			tmpL[j] = 0.0;
			for (k=0; k<nform; k++){
				tmpL[j] += UF[k][i]*Q[k][j];
			}
		}
		for (j=0; j<nform; j++){
//			UF[j][i] = UFold[j][i] + (tmpL[j] - UFold[j][i])*ALFR[j];
			if ( ITL > (ITLlastshift+2) )
			{
				UF[tmpi[j]][i] = UFold[tmpi_old[j]][i] + (tmpL[tmpi[j]] - UFold[tmpi_old[j]][i])*ALFR[j];
			}
			else
			{
				UF[j][i] = tmpL[j];
			}
			
		}
	}

	delete [] tmpL;
}