#include "StdAfx.h"
#include "EL.h"

CEL::CEL(void)
{
	isInitialized = 0;
	isVELcalculated = 0;
	isVPbelowzero = 0;

	ind = NULL;
	fullcrd = NULL;
	STIF = NULL;
	MASS = NULL;
	B = NULL;
	D = NULL;
	vlcrdint = NULL;
	Wvint = NULL;
	faceloccrdint = NULL;
	faceorientcrd_num = NULL;
	faceorientcrd_val = NULL;
	FaceNodes = NULL;
	FaceNodesList = NULL;
	Wfint = NULL;
	FFv = NULL;
	FFf = NULL;
	dFFvloc = NULL;
	dFFfloc = NULL;
	GLOBNE = NULL;
	P = NULL;
	material = NULL;

	MM = new MEM;
}

void CEL::Initialize(int ieltype)
{
	int i,j;

	ielobjnumber = ieltype;
	isInitialized = 0;
	eltype = -1;
	if (ieltype == 5) //�������� ������� �����������
	{
		eltype = ieltype;
		GeneralParameters5();
		GeneralMemInit();
//		Init_VolIntPointLocCrd_22();
//		Init_FaceNodes_22();
		
	}
	if (ieltype == 22) //������������ ������� �����������
	{
		eltype = ieltype;
		GeneralParameters22();
		GeneralMemInit();
//		Init_VolIntPointLocCrd_22();
//		Init_FaceNodes_22();
		
	}
	if (ieltype == 24) //������������ ��������
	{
		eltype = ieltype;
		GeneralParameters24();
		GeneralMemInit();
		Init_VolIntPointLocCrd_24();
		Init_FaceNodes_24();
				
	}
	if (ieltype == 25) //������������ �����
	{
		eltype = ieltype;
		GeneralParameters25();
		GeneralMemInit();
		Init_VolIntPointLocCrd_25();
		Init_FaceNodes_25();
		Init_FacePointLocCrd_25();

		//������ ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		for (i=0; i<Nface; i++)
		{
			for (j=0; j<Nfp; j++)
			{
				dFF_25(dFFfloc[i][j],faceloccrdint[i][j]);
			}
		}
	}

	if (eltype > 0 && eltype != 22)//!!!!!!!!!!��������, ����� ������������ 22 ���
	{
		//������ �������� ������� ����� � ������ �������������� �� ������
		for (i=0; i<NPINT; i++)
		{
			FF_all(FFv[i],vlcrdint[i]);
		}
		//������ ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		for (i=0; i<NPINT; i++)
		{
			dFF_all(dFFvloc[i],vlcrdint[i]);
		}

		isInitialized = 1;
	}
}

void CEL::FF_all(double *F, double *lcrd)
{
	switch (eltype)
	{
	case 5:
		FF_5(F,lcrd);
		break;
	case 22:
		FF_22(F,lcrd);
		break;
	case 24:
		FF_24(F,lcrd);
		break;
	case 25:
		FF_25(F,lcrd);
	default:
		break;
	}
}
void CEL::dFF_all(double **dF, double *lcrd)
{
	switch (eltype)
	{
	case 5:
		dFF_5(dF,lcrd);
		break;
	case 22:
		dFF_22(dF,lcrd);
		break;
	case 24:
		dFF_24(dF,lcrd);
		break;
	case 25:
		dFF_25(dF,lcrd);
	default:
		break;
	}
}

void CEL::GeneralMemInit()
{
	// ��������� ������ ��� ������� ����������
		B = MM->MEM_NEW(B,NDEF,NNE);

		// ������� ������� 
		D = MM->MEM_NEW(D,NDEF,NDEF);
				
		// ������� ��������� ��������
		STIF = MM->MEM_NEW(STIF,NNE,NNE);

		// ������� ���� ��������
		MASS = MM->MEM_NEW(MASS,NNE,NNE);

		// ��������� ���������� ����� �������������� �� ������
		vlcrdint = MM->MEM_NEW(vlcrdint,NPINT,nvlcrd);

		// ������� ������������ ��� �������������� �� ������
		Wvint = MM->MEM_NEW(Wvint,NPINT);

		// ��������� ���������� � ������ �������������� �� ������
		faceloccrdint = MM->MEM_NEW(faceloccrdint,Nface,Nfp,nflcrd);

		// ������� ������������ � ������ �������������� �� ������
		Wfint = MM->MEM_NEW(Wfint,Nface,Nfp);

		// �������� ������� ����� � ������ �������������� �� ������
		FFv = MM->MEM_NEW(FFv,NPINT,NN);

		// �������� ������� ����� � ������ �������������� �� ������
		FFf = MM->MEM_NEW(FFf,Nface,Nfp,NN);

		// �������� ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		dFFvloc = MM->MEM_NEW(dFFvloc,NPINT,NORT,NN);

		// �������� ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		dFFfloc = MM->MEM_NEW(dFFfloc,Nface,NPINT,NORT,NN);

		//����� ���. ����������, ������� ��� ����� �������� ����������
		faceorientcrd_num = MM->MEM_NEW(faceorientcrd_num,Nface);

		//�������� ��������� ����������, ������� ��������� �� �����
		faceorientcrd_val = MM->MEM_NEW(faceorientcrd_val,Nface);
		
		//������� ������ �������������� ����� � ��������� ��������� ������
		FaceNodes = MM->MEM_NEW(FaceNodes,Nface,NN);

		//��������� ������ ����� ������
		FaceNodesList = MM->MEM_NEW(FaceNodesList,Nface,NNface);

		//������ ������� �������� �������
		GLOBNE = MM->MEM_NEW(GLOBNE,NNE);
}

CEL::~CEL(void)
{
	int j;
	eltype = -1;

	if (isInitialized == 1)
	{
	// ��������� ������ ��� ������� ����������
		B = MM->MEM_DEL(B,NDEF,NNE);

		// ������� ������� 
		D = MM->MEM_DEL(D,NDEF,NDEF);
				
		// ������� ��������� ��������
		STIF = MM->MEM_DEL(STIF,NNE,NNE);

		// ������� ���� ��������
		MASS = MM->MEM_DEL(MASS,NNE,NNE);

		// ��������� ���������� ����� �������������� �� ������
		vlcrdint = MM->MEM_DEL(vlcrdint,NPINT,nvlcrd);

		// ������� ������������ ��� �������������� �� ������
		Wvint = MM->MEM_DEL(Wvint,NPINT);

		// ��������� ���������� � ������ �������������� �� ������
		faceloccrdint = MM->MEM_DEL(faceloccrdint,Nface,Nfp,nflcrd);

		// ������� ������������ � ������ �������������� �� ������
		Wfint = MM->MEM_DEL(Wfint,Nface,Nfp);

		// �������� ������� ����� � ������ �������������� �� ������
		FFv = MM->MEM_DEL(FFv,NPINT,NN);

		// �������� ������� ����� � ������ �������������� �� ������
		FFf = MM->MEM_DEL(FFf,Nface,Nfp,NN);

		// �������� ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		dFFvloc = MM->MEM_DEL(dFFvloc,NPINT,NORT,NN);

		// �������� ����������� ������� ����� �� ��������� ����������� � ������ �������������� �� ������
		dFFfloc = MM->MEM_DEL(dFFfloc,Nface,NPINT,NORT,NN);

		//����� ���. ����������, ������� ��� ����� �������� ����������
		faceorientcrd_num = MM->MEM_DEL(faceorientcrd_num,Nface);

		//�������� ��������� ����������, ������� ��������� �� �����
		faceorientcrd_val = MM->MEM_DEL(faceorientcrd_val,Nface);
		
		//������� ������ �������������� ����� � ��������� ��������� ������
		FaceNodes = MM->MEM_DEL(FaceNodes,Nface,NN);

		//��������� ������ ����� ������
		FaceNodesList = MM->MEM_DEL(FaceNodesList,Nface,NNface);

		//������ ������� �������� �������
		GLOBNE = MM->MEM_DEL(GLOBNE,NNE);
	}
	delete MM;
}

void CEL::InitGLOBNE()
{
	int i,j;
	for(i=0; i<NN; i++)
	{
		for (j=0; j<NDOF; j++)
		{
			GLOBNE[i*NDOF+j] = ind[i]*NDOF + j;
		}
	}
}


void CEL::shifting( double shift)
{
	int i;

	for (i=0; i<NNE; i++)
	{
		STIF[i][i] -= shift*MASS[i][i];
	}

}


void CEL::GetBoundBox(double *vmin, double *vmax)
{
	int i,j;
	//����������� ����������� � ������������ ��������� ��������������� ������� ���������������.
	for (j=0; j<NORTfullcrd; j++)
	{
		vmin[j] = 10000000.0;
		vmax[j] = -10000000.0;
	}
	for (i=0; i<NN; i++)
	{
		for (j=0; j<NORTfullcrd; j++)
		{
			if (fullcrd[ind[i]*NORTfullcrd+j] < vmin[j]) vmin[j] = fullcrd[ind[i]*NORTfullcrd+j];
			if (fullcrd[ind[i]*NORTfullcrd+j] > vmax[j]) vmax[j] = fullcrd[ind[i]*NORTfullcrd+j];
		}
	}
}

void CEL::GetCentr(double *vcentr)
{
	int i,j;
	//����������� ������ ��������, ��� �������� �������� ��������� ��� �����
	for (j=0; j<NORTfullcrd; j++)
	{
		vcentr[j] = 0.0;
	}
	for (i=0; i<NN; i++)
	{
		for (j=0; j<NORTfullcrd; j++)
		{
			vcentr[j] += fullcrd[ind[i]*NORTfullcrd+j];
		}
	}
	for (j=0; j<NORTfullcrd; j++)
	{
		vcentr[j] = vcentr[j]/NN;
	}
}