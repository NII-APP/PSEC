#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"


//тест интегрирования движения
int usermain_motionint_eigf()
{
	int i,j,nmat;
	char strpath[256];

	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE_3levels");
	//sprintf(strpath,"C:\\Temp\\Test_Shevchenko");

//	sprintf(strpath,"C:\\Temp\\UGM_LAST\\STATOR");
	//sprintf(strpath,"C:\\Temp\\Reshetka");
sprintf(strpath,"C:\\Temp\\simple1_out_model");

	MEM *pmem = new MEM[1];

	FULLMODEL *pfm;
	SEMODEL *psem;
	MATPROP *pm;
	SSIEIGEN *peig;

	//nmat = 1;
	//pm = new MATPROP[nmat];
	//pm->E = 2.1e11;
	//pm->MU = 0.3;
	//pm->RO = 7800;

	pfm = new FULLMODEL[1];
	pfm->MM = pmem;
	sprintf(pfm->pathmain,"%s",strpath);
//	pfm->nmat = nmat;
//	pfm->mat = pm;
//	pfm->loadtype = 2; //UZOR новый формат
	pfm->loadtype = 3; //Ansys Full Grid

	pfm->ReadFullmodel();
	pfm->GenSurf();


	psem = new SEMODEL[1];
	psem->MM = pmem;
	psem->pfm = pfm;
	//psem->loadtype = 1; //из UZOR
	psem->loadtype = 2; //автоматическое деление	
	psem->ReadSEmodel();


	peig = new SSIEIGEN[1];
	peig->MMM = pmem;
	peig->pfm = pfm;
	peig->psem = psem;
	peig->Init(5,0.00001,15);
	peig->MainSSI();
	
	MOTIONEIGF *pmef;
	pmef = new MOTIONEIGF[1];
	pmef->MM = pmem;
	pmef->AttachToFullModel(pfm);
	pmef->LoadEIGF();
	pmef->InitConstModalDamp(0.1);
	pmef->InitIntegr();
	pmef->InitInertiaUnitModalForces();
	pmef->TestCase_InertiaLoading();

	i=0;

	delete []peig;
	delete []psem;
	delete []pfm;
//	delete []pm;
	delete []pmem;
	delete []pmef;
	
	return 0;
}
