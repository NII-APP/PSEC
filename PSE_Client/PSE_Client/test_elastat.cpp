#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"


int usermain_elastat()
{
	int i,j,nmat;
	char strpath[256];

	sprintf(strpath,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	//sprintf(strpath,"C:\\Temp\\UGM_LAST\\STATOR");
	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE_3levels");
	//sprintf(strpath,"C:\\Temp\\Test_Shevchenko");

	MEM *pmem = new MEM[1];

	FULLMODEL *pfm;
	SEMODEL *psem;
	MATPROP *pm;
	ELASTAT *pest;

	nmat = 1;
	pm = new MATPROP[nmat];
	pm->E = 2.1e11;
	pm->MU = 0.3;
	pm->RO = 7800;

	pfm = new FULLMODEL[1];
	pfm->MM = pmem;
	sprintf(pfm->pathmain,"%s",strpath);
	pfm->nmat = nmat;
	pfm->mat = pm;
	pfm->loadtype = 2;
	pfm->ReadFullmodel();


	psem = new SEMODEL[1];
	psem->MM = pmem;
	psem->pfm = pfm;
	psem->loadtype = 2;
	
	psem->ReadSEmodel();
	
	//psem->AutoMLDivision_S();

	pest = new ELASTAT[1];
	pest->MM = pmem;
	pest->pfm = pfm;
	pest->psem = psem;
	pest->CreateLV();
	pest->StaticSolve();


	i=0;

	delete []pest;
	delete []psem;
	delete []pfm;
	delete []pm;
	delete []pmem;
	
	return 0;
}
