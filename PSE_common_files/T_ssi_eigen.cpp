#include "StdAfx.h"
#include "T_ssi_eigen.h"
#include "math.h"
#include "stdio.h"

SSIEIGEN::SSIEIGEN()
{
	endsolflag = 0;
	MMM = NULL;
	pfm = NULL;
	psem = NULL;
	pelast = NULL;
	
	tmpi = NULL;
	tmpi_old = NULL;
	tmpi_old_2 = NULL;
	low_err = NULL;
//	NDL = NULL;
	tmp = NULL;
	err_freq_mas = NULL;
	err_form_mas = NULL;
	KK = NULL;
	M = NULL;
	MM = NULL;
	LM = NULL;
	LMconv = NULL;
	LMold = NULL;
	LMold_2 = NULL;
	CR = NULL;
	CRold = NULL;
	ALFR = NULL;
	Q  =NULL;
	UFmaxold = NULL;
	UFmax = NULL;
	UF = NULL;
	UFold = NULL;
	UFconv = NULL;
	ortoseqv = NULL;
	
	fp = NULL;
}

SSIEIGEN::~SSIEIGEN()
{
	tmpi = MMM->MEM_DEL(tmpi,nur);
	tmpi_old = MMM->MEM_DEL(tmpi_old,nform);
	tmpi_old_2 = MMM->MEM_DEL(tmpi_old_2,nform);
	low_err = MMM->MEM_DEL(low_err,nform);
	tmp = MMM->MEM_DEL(tmp,nur);
	err_freq_mas = MMM->MEM_DEL(err_freq_mas,nform);
	err_form_mas = MMM->MEM_DEL(err_form_mas,nform);
	KK = MMM->MEM_DEL(KK,nform,nform);

	MM = MMM->MEM_DEL(MM,nform,nform);
	LM = MMM->MEM_DEL(LM,nform);
	LMconv = MMM->MEM_DEL(LMconv,NEigen);
	LMold = MMM->MEM_DEL(LMold,nform);
	LMold_2 = MMM->MEM_DEL(LMold_2,nform);
	CR = MMM->MEM_DEL(CR,nform);
	CRold = MMM->MEM_DEL(CRold,nform);
	ALFR = MMM->MEM_DEL(ALFR,nform);
	Q  = MMM->MEM_DEL(Q,nform,nform);
	UFmaxold = MMM->MEM_DEL(UFmaxold,nform);
	UFmax = MMM->MEM_DEL(UFmax,nform);
	UF = MMM->MEM_DEL(UF,nform,nur);
	UFold = MMM->MEM_DEL(UFold,nform,nur);
	UFconv = MMM->MEM_DEL(UFconv,NEigen,nur);
	ortoseqv = MMM->MEM_DEL(ortoseqv,nform);
	
	delete [] pelast;

	if (fp != NULL) { fclose(fp); fp = NULL;}
}

void SSIEIGEN::statusbar(char *str)
{
	printf ("\n%s\n",str);
}


int SSIEIGEN::MainSSI()
{
	int		i, j, k, *tmpiL;
	double	s, ss, s_freq, s_form,rconv;
	char	strl[256], nameX[256], tmpbuf[128];
	FILE	*fp2, *fpnf;

	endsolflag = 0;

	while (endsolflag != 1)
	{
		// вызвать упругое решение для всех векторов UF, необходимо умножить этивектора на матрицу масс
		for (i=0; i<nform; i++)
		{
			for (j=0; j<nur; j++)
			{
				UF[i][j] = UF[i][j]*M[j];
			}
		}

		pelast->StatMultSolveExt(nform,UF);

		NormalizeEigVect();

		sprintf( strl," ITL.......GrammShmidt_Check 0......",ITL);
		statusbar( strl );

		CollinearDetection();
		
		// определение строк матриц жесткости и масс в подпространстве
		sprintf( strl," ITL..%d.......make subspace matrix......",ITL);
		statusbar( strl );
		_strtime_s( tmpbuf, 128 );
		fprintf(fp,"%sOS time:\t\t\t\t%s\n",strl,tmpbuf);
		fflush(fp);
		
		//KKmaincalc();
		//KKmaincalcByEl();
		psem->CalcSubspaceMatrix(KK);

		//вывод KK
		fprintf(fp,"\n\nIT=%d   KK _after KKmaincalcByEl !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",ITL);
		for (i=0; i<nform; i++)
		{
			for (j=0; j<nform; j++)
			{
				fprintf(fp,"%e\t",KK[i][j]);
			}
			fprintf(fp,"\n");
		}
		
		//расчет MM
		for (i=0; i<nform; i++){
			for (j=0; j<nur; j++){
				tmp[j] = UF[i][j]*M[j];
			}
			
			for (j=0; j<nform; j++){
				s=0;
				for (k=0; k<nur; k++){
					s += tmp[k]*UF[j][k];
				}
				MM[i][j] = s;
			}
		}


		sprintf( strl," ITL.......Jacoby started......",ITL);
		statusbar( strl );
		
		
		//определение собственных векторов и собственных частот в подпространстве
		Jacoby();


		//вывод KK
		fprintf(fp,"\n\nIT=%d   KK _after JACOBY\n",ITL);
		for (i=0; i<nform; i++)
		{
			for (j=0; j<nform; j++)
			{
				fprintf(fp,"%e\t",KK[i][j]);
			}
			fprintf(fp,"\n");
		}
		fflush(fp);
		//вывод MM
		fprintf(fp,"\n\nIT=%d   MM _after JACOBY\n",ITL);
		for (i=0; i<nform; i++)
		{
			for (j=0; j<nform; j++)
			{
				fprintf(fp,"%e\t",MM[i][j]);
			}
			fprintf(fp,"\n");
		}
		fflush(fp);

		

//проверка сходимости и рассчет скоростей сходимости

		//определение норм (максимальных компонент) собственных векторов
		UFGetMax();

		//сортировка по возрастанию частот
		FormSort( );

		fprintf(fp,"\n-------------LM------------------\n");
		for (i=0; i<nform; i++){
			double frhertz = sqrt((LM[tmpi[i]] + shifting))/(2*3.1415926);
			fprintf(fp,"LM[%d] = %e  + shift = %e  fr_hertz = %e\n",i,LM[tmpi[i]],(LM[tmpi[i]] + shifting),frhertz);
		}
		fprintf(fp,"\n-------------LM------------------\n");
		fflush(fp);

		//определение параметров сходимости итерационного процесса
		ErrorCalc();

		//определение скоростей сходимости и коэффициентов ускорения

		ConvAccelCalc();

		//определение нового приближения к собственным формам
		NewForms();

		sprintf( strl," ITL.......GrammShmidt 1......",ITL);
		printf( strl );
		GrammShmidtCyckl();


		sprintf( strl," ITL.......GrammShmidt_Check 1......",ITL);
		printf( strl );

		GrammShmidtCheck();		
		
		for (i=0; i<nform; i++){
			for (j=0; j<nur; j++) UFold[i][j] = UF[i][j];
		}
		

		ITL++;
		//запись для продолжения решения
		for( i = 0; i < 256; i++ ) strl[i] = 0;
		sprintf( strl,"%s\\UFtmp.dat\0", pfm->pathmain );

		fp2 = fopen ( strl ,"wb");
		fwrite(&ITL,sizeof(int),1,fp2);
		fwrite(&err_freq,sizeof(double),1,fp2);
		fwrite(&err_form,sizeof(double),1,fp2);
		fwrite(UFmaxold,sizeof(double),nform,fp2);
		fwrite(LMold,sizeof(double),nform,fp2);
		fwrite(tmpi_old,sizeof(int),nform,fp2);
		fwrite(low_err,sizeof(int),nform,fp2);
		for(i=0; i<nform; i++){
			fwrite(UF[i],sizeof(double),nur,fp2);
		}		 
		fclose(fp2);
		for( i = 0; i < 256; i++ ) strl[i] = 0;
		sprintf( strl,"%s\\nf_tmp.txt\0", pfm->pathmain );
		
		fpnf = fopen ( strl ,"w");
		fprintf(fpnf,"%d %d",nf,nform);
		fclose(fpnf);



		PrintRez_EigenSSI_iter();

// Проверка необходимости осуществления шифтинга
// Проверка условия сходимости по частотам и формам в блоке для nf частот

		ConvCheck();

		// Проверка окончания решения

		if( nfconv == NEigen)
		{
			endsolflag = 1;
		}

		sprintf( strl," ITER  SSI end ....ITL %d.............",ITL);
		statusbar( strl );
	}

//	PrintRez_EigenSSI_full();

	return 0;
}


void SSIEIGEN::PrintRez_EigenSSI_iter()
{
	//вывод информации
	char strl[256], nameX[256], tmpbuf[128];
	int m3, m4, m5, i1, i2, i3, dirlen, nvalue, in, n, k, kk, jj, j, i;
	int *tmpiL;
	float *tus8;
	double s;
	double pi2,lamda,omega,frc,*tus8d;
	double form_mashtab = 1;
	FILE *fpL;

	int KND = pfm->NN; 

	
//	tus8 = new float[nur];
//	tmpiL = new int[nur];
	pi2 = 3.1415926 * 2.0;

	sprintf( strl," ITER  SSI print ....IT %d.............",ITL);
	statusbar( strl );
	FormSort( );	// заполняется массив tmpi

//	for( k = 0; k < nf; k++ ){
//		kk = k + 1;
//		n = 0;
//		s=0;
//		for (jj=0; jj<nur; jj++)  s += UF[tmpi[k]][jj]*M[jj]*UF[tmpi[k]][jj];
//		s = sqrt(s); 
//		for ( j = 0; j < KORT; j++){
//			for ( i = 0; i < KND; i++){
//				in = i * KORT + j;
//				tus8[n] = (float) UF[tmpi[k]][in]/s; 
//				tus8[n] = tus8[n]*form_mashtab;
//				n++;
//			}
//		}
//		for( i = 0; i < 256; i++ ) nameX[i] = 0;
//		sprintf( nameX,"it_forma0000.sba\0");
//		for( i = 0; i < 256; i++ ) strl[i] = 0;
//		sprintf( strl,"%s\0", FName[1] );
//		strcat( strl, nameX );
//		for( i = 0; i < strlen(strl); i++ ) nameX[i] = strl[i];
//		dirlen = strlen( FName[1] );
//		m3 = dirlen + 9;
//		m4 = dirlen + 10;
//		m5 = dirlen + 11;
//		i1 = (kk) % 10;
//		i2 = ((kk) / 10) % 10;
//		i3 = (kk) / 100;
//		nameX[m3] = Cff[i3]; nameX[m4] = Cff[i2]; nameX[m5] = Cff[i1];
//		if(( fpL = fopen( nameX,"wb" ) ) == NULL ) goto e1;
//		fwrite( tus8, sizeof(float), nur, fpL );
//		fclose( fpL );
//	}
//e1:	{};
//	delete [] tus8;
//
//	tus8d = new double [nur];
//
//	for( i = 0; i < 256; i++ ) strl[i] = 0;
//	sprintf( strl,"%s\0", FName[1] );
//	strcat( strl, "it_evalues.wrm" );
//	if(( fpL = fopen( strl,"wb" ) ) == NULL ){
//		fprintf( FLST,"PrintRez_EigenSSI: Wrong file name=%s (wb)\n", strl );
//		goto e2;
//	}
//	fwrite( &nf, sizeof(int), 1, fpL );
//	for( j = 0; j < nf; j++ ){
//		lamda = LM[tmpi[j]];
//		omega = sqrt(LM[tmpi[j]]);
//		frc = omega/pi2;
//		fwrite( &lamda, sizeof(double), 1, fpL );
//		fwrite( &omega, sizeof(double), 1, fpL );
//		fwrite( &frc, sizeof(double), 1, fpL );
//		s=0;
//		for (k=0; k<nur; k++) s += UF[tmpi[j]][k]*M[k]*UF[tmpi[j]][k];
//		s = sqrt(s);
//		for (k=0; k<nur; k++) tus8d[k] = (UF[tmpi[j]][k]/s)*sqrt(M[k]);
//
//		fwrite( tus8d, sizeof(double), nur, fpL );
//	}
//	for (i=0; i<nur; i++) tus8d[i] = sqrt(M_BN[i]);
//	fwrite( tus8d, sizeof(double), nur, fpL ); 
//	fclose( fpL );
//e2:	{};
		for( i = 0; i < 256; i++ ) strl[i] = 0;
		sprintf( strl,"%s\\it_freq.lst\0", pfm->pathmain );
		
	fpL = fopen ( strl ,"w");
	fprintf(fpL,"\nkk=%d Eigen Values were calculated in ITL = %d\n\n", nf, ITL );
	for( j = 0; j < nf; j++ ){
		lamda = LM[tmpi[j]];
		omega = sqrt(LM[tmpi[j]]);
		frc = omega/pi2;
		fprintf( fpL,"j=%d  lamda=%e omega=%e frc=%e\n",j+1, lamda, omega, frc );
	}
	fclose(fpL);
//	delete [] tmpiL;
//	delete [] tus8d;

	
}

void SSIEIGEN::PrintRez_EigenSSI_full()
{
	//вывод информации
	char strl[256];
	int  j,i,k;
	double maxd,s;
	double modm;
	double pi2,lamda,omega,frc;
	FILE *fpL;
	int KND;

	KND = pfm->NN;
	pi2 = 3.1415926 * 2.0;

	sprintf( strl," ITER  SSI print ....IT %d.............",ITL);
	statusbar( strl );

	sprintf( strl,"%s\\EIGF.bin\0", pfm->pathmain );
	fpL = fopen(strl,"wb");
	fwrite(&pfm->NNE,sizeof(int),1,fpL);
	fwrite(&nfconv,sizeof(int),1,fpL);
	for (j=0; j < nfconv; j++)
	{
		omega = sqrt(LMconv[j]);
		fwrite(&omega,sizeof(double),1,fpL);
		modm = 0.0;
		for (k=0; k < pfm->NNE; k++)
		{
			modm += UFconv[j][k]*M[k]*UFconv[j][k];
		}
		fwrite(&modm,sizeof(double),1,fpL);
		fwrite(UFconv[j],sizeof(double),pfm->NNE,fpL);
	}
	fclose(fpL);

	for( i = 0; i < 256; i++ ) strl[i] = 0;
	sprintf( strl,"%s\\freq.lst\0", pfm->pathmain );
		
	fpL = fopen ( strl ,"w");
	fprintf(fpL,"\nkk=%d Eigen Values were calculated in ITL = %d\n\n", nfconv, ITL );
	for( j = 0; j < nfconv; j++ ){
		lamda = LMconv[j];
		omega = sqrt(LMconv[j]);
		frc = omega/pi2;
		fprintf( fpL,"j=%d  lamda=%e omega=%e frc=%e\n",j+1, lamda, omega, frc );
	}
	fclose(fpL);

	//печать для Paraview
	for (i=0; i<nfconv; i++)
	{
		sprintf( strl,"%s\\EIGF_view_%d.vtk\0", pfm->pathmain, i+1 );
		fpL = fopen(strl,"w");
		pfm->ParaView_PrintGrid(fpL);
		pfm->ParaView_StartCellDataSection(fpL);
		pfm->ParaView_PrintMaterial(fpL);
		pfm->ParaView_StartNodeDataSection(fpL);
		
		//для вывода собственной формы осуществляется прямая и обратная нормировка к единичному максимальному отклонению узла
		maxd = 0.0;
		for (j=0; j<pfm->NN; j++)
		{
			s = 0.0;
			for (k=0; k < pfm->KORT; k++)
			{
				s += UFconv[i][j*pfm->KORT + k]*UFconv[i][j*pfm->KORT + k];
			}
			if ( s > maxd ) maxd = s;
		}
		maxd = sqrt(maxd);
		for  (j=0; j < pfm->NNE; j++)
		{
			UFconv[i][j] = UFconv[i][j]/maxd;
		}

		sprintf(strl,"eigen_form");
//		pfm->ParaView_SingleXYZS(fpL,strl,UFconv[i]);
		pfm->ParaView_SingleVector(fpL,strl,UFconv[i]);
		fclose(fpL);

		for  (j=0; j < pfm->NNE; j++)
		{
			UFconv[i][j] = UFconv[i][j]*maxd;
		}
	}

}