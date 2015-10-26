#include "StdAfx.h"
#include "T_ssi_eigen.h"
#include "math.h"
#include "stdio.h"


void SSIEIGEN::ErrorCalc()
{
	int i,j;
	double s_freq,s_form;

	char tmpbuf[256];

	if (ITL > ITLlastshift)
	{
		err_freq_old = err_freq;
		err_form_old = err_form;

		err_freq = 0;
		err_form = 0;
		for (i=0; i<nform; i++)
		{
			err_freq_mas[i] = abs((LM[tmpi[i]] - LMold[tmpi_old[i]])/LM[tmpi[i]]);
			err_form_mas[i] = abs((UFmax[tmpi[i]] - UFmaxold[tmpi_old[i]])/UFmax[tmpi[i]]);
		}
		for (i=0; i<nf; i++)
		{
			s_freq = abs((LM[tmpi[i]] - LMold[tmpi_old[i]])/LM[tmpi[i]]);
			s_form = abs((UFmax[tmpi[i]] - UFmaxold[tmpi_old[i]])/UFmax[tmpi[i]]);
			if ( low_err[tmpi[i]] == 0 && err_freq < s_freq) err_freq = s_freq;
			fprintf(fp,"low_err %d = %d\n",i,low_err[tmpi[i]]);
			fprintf(fp,"freq %d = %e   err_freq = %e  \n",i,LM[tmpi[i]],s_freq);
			
			if ( low_err[tmpi[i]] == 0 && err_form < s_form) err_form = s_form;
			fprintf(fp,"form %d = %e   err_form = %e  \n",i,UFmax[tmpi[i]],s_form);
			fprintf(fp,"\n");
			fflush(fp);
		}
	}
	_strtime_s( tmpbuf, 128 );
	fprintf(fp," TOLLERANCE  ITL = %d    freq = %f  form = %f      OS time:\t%s\n\n\n",ITL,err_freq,err_form,tmpbuf);
	fflush(fp);

//определение количества сошедшихся до требуемой точности частот, следующих подряд
	consecutive_convfreq = 0;
	for (i=0; i<nform; i++)
	{
		if ( err_freq_mas[i] < eps_freq )
		{
			consecutive_convfreq++;
		}
		else
		{
			break;
		}
	}

}

void SSIEIGEN::ConvAccelCalc()
{
	int i,j;
	double rconv;


	fprintf(fp,"\n\n CONVERGENCE RATE AND ACCEL   ITL = %d  ITLlastshift  = %d\n",ITL,ITLlastshift);

	if (ITL > (ITLlastshift+1))
	{
		//for (i=0; i<nform; i++)
		for (i=0; i<nf; i++)
		{
			CR[i] = abs((LM[tmpi[i]] - LMold[tmpi_old[i]])/(LMold[tmpi_old[i]] - LMold_2[tmpi_old_2[i]]));
			fprintf(fp,"CR %d = %e\n",i,CR[i]);
		}
	}

	fprintf(fp,"\n");

	nlqp1 = 0;
	lqp1_sum = 0.0;
//		lqp1 = 0.0;
	if (ITL > (ITLlastshift+2))
	{
		//for (i=0; i<nform; i++)
		for (i=0; i<nf; i++)
		{
			rconv = abs(CR[i] - CRold[i])/abs(CR[i]);
			fprintf(fp,"rconv %d = %e         err_freq_mas = %e \n",i,rconv,err_freq_mas[i]);
			if (  (rconv < 0.25) && (err_freq_mas[i] < 0.001) )
			{
				lqp1_sum += (LM[tmpi[i]]+shifting)/sqrt(CR[i]);
				nlqp1++;
			}
		}

		if (nlqp1 > 1) lqp1 = lqp1_sum/((double)nlqp1);

//			nlqp1 = 1;
		//if (nfconv == 0) lqp1 = LM[tmpi[0]]+LM[tmpi[nform-1]];
		//else lqp1 = LMconv[0]+LM[tmpi[nform-1]];
		//if ( (nf+5) < nform) lqp1 = LM[tmpi[0]]+LM[tmpi[nf+5]];
		//else lqp1 = LM[tmpi[0]]+LM[tmpi[nform-1]];
		

	}

	fprintf(fp,"lqp1 = %e    nlqp1 = %d\n",lqp1,nlqp1);

	if (nlqp1 > 1)
	{
		for (i=0; i<nform; i++) 
		{
			ALFR[i] = 1.0/(1.0 - (LM[tmpi[i]]+shifting)/lqp1);
			if (ALFR[i] > 10.0) ALFR[i] = 10.0;
			if (ALFR[i] < 0.0) ALFR[i] = 1.0;
			
			/*rconv = abs(CR[i] - CRold[i])/abs(CR[i]);
			if (  (rconv < 0.25) && (err_freq_mas[i] < 0.001) )
			{
				ALFR[i] = 1.0/(1.0 - (LM[tmpi[i]]+shifting)/lqp1);
				if (ALFR[i] > 10.0) ALFR[i] = 10.0;
				if (ALFR[i] < 0.0) ALFR[i] = 1.0;
			}
			else
			{
				ALFR[i] = 1.0;
			}*/
		}
	}
	else
	{
		for (i=0; i<nform; i++) ALFR[i] = 1.0;
	}

	for (i=0; i<nform; i++) fprintf(fp,"ALFR %d = %e \n",i,ALFR[i]);

	fflush(fp);

	for (i=0; i<nform; i++){
		LMold_2[i] = LMold[i];
		tmpi_old_2[i] = tmpi_old[i];
		CRold[i] = CR[i];
		UFmaxold[i] = UFmax[i];
		LMold[i] = LM[i];
		tmpi_old[i] = tmpi[i];
	}
}

void SSIEIGEN::ConvCheck()
{
	int i,j;

	if( (err_freq <= eps_freq) || 
		( (abs(ITLlastshift - ITL) > 15) && (consecutive_convfreq > (int)(nf*0.3)) ) ||
		(abs(ITLlastshift - ITL) > 30) )
	{
		CopyConvForms(); //копирование форм для сошедшихся частот в массив для постоянного хранения
		PrintRez_EigenSSI_full();
		if ( (nfconv < NEigen)  ){
		//if ( (SSIL->nfconv < NEigen) && ( abs(SSIL->ITLlastshift - SSIL->ITL) <= 20 ) ){
//					SSIL->NewRandStartVect(); //замена 
			DefineNewShift();//определение нового смещения по частоте
			if ( shifting > shiftingold )
			{
				psem->SetNewShiftSEM(shifting);
			}
		}
	}

}

void SSIEIGEN::CopyConvForms()
{
	int i,j,k,nrepeate,nfcur;
	double s,tol;
	int *fl_nocopy;

	tol = 1.e-3;

	fl_nocopy = new int[nf];
	for (i=0; i<nf; i++) fl_nocopy[i] = 0;

	nrepeate = nfconv - LastConvNum;

	//поиск вновом блоке СФ, повторяющих старое решение
	for (i = LastConvNum+1; i<nfconv; i++)
	{
		for (j=0; j<nf; j++)
		{
			if (fl_nocopy[j] == 0)
			{
				s = 0.0;
				for (k=0; k<nur; k++)
				{
					s += UF[tmpi[j]][k]*M[k]*UFconv[i][k];
				}
				if (abs(s)>tol) fl_nocopy[j] = 1; 
				fprintf(fp,"\niconvnum = %d     infnum = %d     s = %e   nocopy = %d\n",i,j,s,fl_nocopy[j]);
				fflush(fp);
			}
		}
	}

	nfconvold = nfconv;

	for (j=0; j<nform; j++) ortoseqv[j] = -1;
	
	srand((unsigned)time( NULL ) );
	for (i=0; i<nf; i++)
	{
		if ( err_freq_mas[i] > eps_freq ) break;
		if (fl_nocopy[i] == 0)
		{
			LMconv[nfconv] = LM[tmpi[i]]+shifting;
//			LMconv[nfconv] = LM[tmpi[i]];
			for (j=0; j<nur; j++)
			{
				UFconv[nfconv][j] = UF[tmpi[i]][j];
			}
			nfconv++;

			s=0.0;
			nfcur = tmpi[i];
			for (j=0; j<nur; j++){
				UF[nfcur][j] = (((double) rand())/((double) RAND_MAX)) - 0.5;
				s += UF[nfcur][j]*UF[nfcur][j];
			}
			s = sqrt(s);
			for (j=0; j<nur; j++)	UF[nfcur][j] = UF[nfcur][j]/s;
			
		}
		ortoseqv[nform-1-i] = i;
		if (nfconv == NEigen) break;
	}
	for (j=i; j<nform; j++)
	{
		ortoseqv[j-i] = j;
	}

	for (j=0; j<nform; j++) fprintf(fp,"ortoseqv %d = %d\n",j,ortoseqv[j]);

	fprintf(fp,"\n\nNEWNFCONV = %d\n\n",nfconv);
	fflush(fp);

	//ортогонализация до требуемой степени ортогональности
	//GrammShmidtCyckl();
	GrammShmidtCyckl_SpecSeqv();

	delete []fl_nocopy;
}


void SSIEIGEN::DefineNewShift()
{
	int i,j;
	double frrelint,alowsh;
	double shifting_treshold;

	shifting_treshold = 5.0;


	shiftingold = shifting;

//	alowsh = LMconv[nfconvold] + 0.3*(LMconv[nfconv-1] - LMconv[nfconvold]);
//	fprintf(fp,"Allowable shift = %e\n\n",alowsh);

	for ( i = nfconv-1; i > 0; i-- )
	//i = nfconv;
	//while ( shifting > alowsh )
	{
//		i--;
		frrelint = (LMconv[i] - LMconv[i-1])/ LMconv[i];
		if (frrelint > minfr_interval)
		{
			shifting = (LMconv[i] + LMconv[i-1])/2;
			break;
		}
	}

	
	if ( shifting > 0.0000001 )
	{
		if ( sqrt( (shifting - shiftingold)/shifting ) < minfr_interval ) shifting = shiftingold;
	}
	if ( (ITL - ITLlastshift) < 7 ) shifting = shiftingold; //скорость сходимости пока достаточная, не обязательно делать сдвиг

//	shifting = 0.0;


	LastConvNum = i-1;

	fprintf(fp,"LastConvNum = %d\nShift = %e\n\n",i-1,shifting);
	fflush(fp);

	ITLlastshift = ITL;
	err_freq = 1.0;
	err_form = 1.0;

	for (i=0; i<nform; i++){
		ALFR[i] = 1.0;
		CR[i] = 0.0;
		CRold[i] = 0.0;
	}

	
	lqp1 = 0.0;
	lqp1_sum = 0.0;
	nlqp1 = 0;

}