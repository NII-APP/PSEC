#include "StdAfx.h"
#include "T_ssi_eigen.h"
#include "math.h"
#include "stdio.h"

int SSIEIGEN::Init(int NEigen, double EpsEigen, int nf)
{
	int		i, j, k, n, mmm, nu, IP, *tmpiL, ING, IOBR;
	long long int	ii;
	float	aa, bb, cc;
	char	strl[256], nameX[256], tmpbuf[128];
	FILE	*fp2, *fpnf, *Fp;

	for( i = 0; i < 256; i++ ) strl[i] = 0;
	sprintf( strl,"%s\\Eigen_SSI_SE_listing.lst", pfm->pathmain );
	fp = fopen ( strl ,"w");
	for( i = 0; i < 256; i++ ) strl[i] = 0;
	sprintf( strl,"%s\\ortogonalization.lst", pfm->pathmain );
	fp3 = fopen ( strl ,"w");

	_tzset();

	_strtime_s( tmpbuf, 128 );
	fprintf(fp,"SSI  initializetion........................OS time:\t\t\t\t%s\n",tmpbuf);
	fflush(fp);

	sprintf( strl," ITER  SSI start ...");
	statusbar( strl );

	nur = pfm->NNE;

	fprintf(fp,"nur = %d\n",nur);
	fflush(fp);

	////количество потоков

//	nthread = NTHN+2;
	nthread = 4;
	blocksize = 100;

//	nthread = 1;
//	blocksize = 400000;


	//задается в MainFrm
//	nf = 6; //количество частот, поиск которых осуществляется в одном блоке перед шифтингом
	//nf = NEigen;

	

	eps_freq = EpsEigen;
	eps_form = eps_freq*10;
	if (nf>NEigen) nf = NEigen;
	this->nf = nf;
	this->NEigen = NEigen;

	nfconv = 0;
	nfconvold = 0;
	nf_old = nf;
	nform = (2*nf < nf+8)?(2*nf):(nf+8);
	shifting = 0.0;
	shiftingold = 0.0;

	LastConvNum = -1;
	ITL=0;
	ITLlastshift = 0;
	err_freq = 1;
	err_form = 1;

	nlqp1 = 0;
	lqp1 = 0.0;
	lqp1_sum = 0.0;

	fl_isBlockConv = 0;
	minfr_interval = 0.01;

	low_err = new int[nform];
	err_freq_mas = new double[nform];
	err_form_mas = new double[nform];
	tmp = new double[nur];
	tmpi = new int[nur];
	tmpi_old = new int[nform];
	tmpi_old_2 = new int[nform];
	ortoseqv = new int[nform];
	UF = new double*[nform];
	for (i=0; i<nform; i++) UF[i] = new double[nur];
	UFold = new double*[nform];
	for (i=0; i<nform; i++) UFold[i] = new double[nur];
	UFconv = new double*[NEigen];
	for (i=0; i<NEigen; i++) UFconv[i] = new double[nur];

	KK = new double*[nform];
	for (i=0; i<nform; i++)	KK[i] = new double[nform];
	MM = new double*[nform];
	for (i=0; i<nform; i++)	MM[i] = new double[nform];
	LM = new double[nform];
	LMconv = new double[NEigen];
	LMold = new double[nform];
	LMold_2 = new double[nform];
	CR = new double[nform];
	CRold = new double[nform];
	ALFR = new double[nform];


	UFmaxold = new double[nform];
	UFmax = new double[nform];
	Q = new double*[nform];
	for (i=0; i<nform; i++) Q[i] = new double[nform];
	for (i=0; i<nform; i++){
		LMold[i] = 0.0;
		UFmaxold[i] = 0.0;
		ALFR[i] = 1.0;
		LM[i] = 0.0;
		LMold_2[i] = 0.0;
		CR[i] = 0.0;
		CRold[i] = 0.0;
		ortoseqv[i] = i;
		low_err[i] = 0;
	}

	////найти диагональную матрицу масс в полной модели
	//
	pfm->InitIntPoint();
	pfm->MassMatrix();
	pfm->MassDiagMatrixFix();
	M = pfm->Mdiag;

	//1 - подготовка начального приближения собственных векторов конструкции

	StandartStartingVectors();

	psem->SetNewShiftSEM(shifting);
	psem->SetOperationType(2);//при расчете собственных частот требуется сохранять полнлную матрицу СЭ в безнулевой форме

	pelast = new ELASTAT[1];
	pelast->pfm = pfm;
	pelast->psem = psem;

	_strtime_s( tmpbuf, 128 );
	fprintf(fp,"ITerations started........................OS time:\t\t\t\t%s\n",tmpbuf);
	fflush(fp);

	return 0;
}

void SSIEIGEN::StandartStartingVectors()
{
	char	strl[256], strL[256], tmpbuf[128];
	int		i,j,jj,k,ff,inum,*tmpiL;
	int		DOFstep, DOFm,nrandforces,icur;
	double	s, maxold,tmp;
	FILE*	fp2;

	tmpiL = new int [nur];

	_strtime_s( tmpbuf, 128 );
	fprintf(fp,"StandartStartingVector........................OS time:\t\t\t\t%s\n\n\n",tmpbuf);
	fflush(fp);

	srand((unsigned)time( NULL ) );

	//nrandforces = int( nur/(5*nform) ); //!!!!!!!!!!!!!!!!!
	//if (nrandforces < 1) nrandforces = 1;
	nrandforces = 10;
	if ( nrandforces*nform*2 > nur ) nrandforces = 1; 

	for (i=0; i<nur; i++)
	{
		tmpiL[i] = 0;
	}

	for ( i=0; i<nform; i++)
	{
		for (j=0; j<nur; j++)
		{
			UF[i][j] = 0.0;
		}
	}

	int irf = 0;
	for (i=0; i<(nform-2); i++)
	{
		irf = 0;
		while (irf < nrandforces)
		{
			tmp = (((double) rand())/((double) RAND_MAX))*nur;
			icur = (int)tmp;
			if ( (M[icur] > 1.e-18) && tmpiL[icur] == 0 )
			{
				UF[i][icur] = 1.0;
				tmpiL[icur] = 1;
				irf++;
			}
			else
			{
				irf = irf;
			}
		}
	}

	//оставшиеся собственные векторы заполняем случайными числами
	srand((unsigned)time( NULL ) );
	for (i=nform-2; i<nform; i++){
	//for (i=0; i<nform; i++){
		s=0;
		for (j=0; j<nur; j++){
			UF[i][j] = (((double) rand())/((double) RAND_MAX)) - 0.5;
			s += UF[i][j]*UF[i][j];
		}
		s = sqrt(s);
		for (j=0; j<nur; j++)	UF[i][j] = UF[i][j]/s; 
	}

	//ортогонализация до требуемой степени ортогональности
	GaussShmidt_3();

	delete [] tmpiL;
}

