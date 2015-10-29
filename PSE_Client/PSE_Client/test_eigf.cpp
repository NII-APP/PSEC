#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"


//���� (������ �� ������ � ����������� ������)
int usermain_eigf()
{
	int i,j,nmat;
	char strpath[256];

sprintf(strpath,"C:\\Temp\\simple1_out_model");

	MEM *pmem = new MEM[1];

	FULLMODEL *pfm;
	SEMODEL *psem;
	MATPROP *pm;
	SSIEIGEN *peig;

	pfm = new FULLMODEL[1];
	pfm->MM = pmem;
	sprintf(pfm->pathmain,"%s",strpath);
	pfm->loadtype = 3; //Ansys Full Grid

	pfm->ReadFullmodel();

	psem = new SEMODEL[1];
	psem->MM = pmem;
	psem->pfm = pfm;
	psem->loadtype = 2; //�������������� �������	
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