#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"

extern char glob_str_taskpath[256];

//тест (расчет СЭ модели и собственных частот)
int usermain_eigf()
{
	int i,j,nmat;
	char strpath[256];

//sprintf(strpath,"C:\\Temp\\simple1_out_model");
	sprintf(strpath,"%s",glob_str_taskpath);

	MEM *pmem = new MEM[1];

	FULLMODEL *pfm;
	SEMODEL *psem;
	MATPROP *pm;
	SSIEIGEN *peig;

	pfm = new FULLMODEL[1];
	pfm->MM = pmem;
	sprintf(pfm->pathmain,"%s",strpath);
	pfm->loadtype = 4; //NX .inp file for ANSYS

	pfm->ReadFullmodel();

	psem = new SEMODEL[1];
	psem->MM = pmem;
	psem->pfm = pfm;
	psem->loadtype = 2; //автоматическое деление	
	psem->ReadSEmodel();

	peig = new SSIEIGEN[1];
	peig->MMM = pmem;
	peig->pfm = pfm;
	peig->psem = psem;
	peig->Init(2,0.00001,15);
	peig->MainSSI();
	

	i=0;

	delete []peig;
	delete []psem;
	delete []pfm;
	delete []pmem;
	
	return 0;
}