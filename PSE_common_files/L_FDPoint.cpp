#include "StdAfx.h"
#include "L_FDPoint.h"
#include "FullModel.h"
#include "math.h"

FDPOINT::FDPOINT(void)
{
	ConnectedNodes = new int[20]; //������ �������� � �������, �� ����������� ����� ����� �� ����������� ����� ���������
	FF = new double[20];
	ncn = 0;
	connected_el = -1;
	isactive = false;
}

FDPOINT::~FDPOINT(void)
{
	delete [] ConnectedNodes;
	delete [] FF;
}

void FDPOINT::TrackingPoint( void *vpfm, double *gcrd)
{
	if (isactive == true) AttachToFullModel(vpfm,gcrd,connected_el,&lcrd[0]);
}

void FDPOINT::TrackingPoint( void *vpfm, float *gcrdf)
{
	int k;
	double gcrd[3];
	for (k=0; k<3; k++) gcrd[k] = (double)gcrdf[k];
	if (isactive == true) AttachToFullModel(vpfm,&gcrd[0],connected_el,&lcrd[0]);
}

int FDPOINT::AttachToFullModel( void *vpfm, double *gcrd, int posel, double *poslcrd)
{
	//����������
	// 1 - ���� ����� ������� � �������� posel;
	// 2 - ���� ����� ������� � ������ ��������;
	// 3 - ���� ����� �� ������ �� � ���� �������, �� ���������������� �� ��������� ����� � ��������� � ��� ���������
	// -1 - ���� ����� ������ �� �������

	int i;
	int iel,ieltype;
	int glcretflag,fl;
	double FFtmp[20];
	double vmin[3], vmax[3], vcentr[3];
	double mindist,dist;
	int ielmindist;

	FULLMODEL *pfm;
	pfm = (FULLMODEL *)vpfm;

	for (i=0; i<pfm->KORT; i++)
	{
		crd[i] = gcrd[i];
	}

	//����� ���������� � �������� "posel", ������� �� ������� ������� ����� ��������� �������� ��������� ������ ��������� ����� ���������� ���� (����������� �����������)
	//�� ���� ����� ���������� �����
	iel = 0;
	if ( posel >= 0 )
	{
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = poslcrd[i];
		}
		iel = posel;
	}
	else
	{
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = 0.0;
		}
	}

	pfm->ClearLIST_EL();

	pfm->LIST_EL[pfm->nlistel] = iel;
	pfm->MASK_EL[iel] = 1;
	pfm->nlistel++;

	int ilistel = 0; //������� ������������ ����� �������� � ������

	ielmindist = -1;
	mindist = 100000.0;
	while ( pfm->nlistel < pfm->NEL || ilistel < pfm->nlistel )
	{
		iel = pfm->LIST_EL[ilistel];
		
		ieltype = pfm->IND[iel][0];
		pfm->AttachElement(&pfm->el[ieltype],iel);
		pfm->el[ieltype].GetCentr(&vcentr[0]);
		dist = 0.0;
		for (i=0; i<pfm->KORT; i++)
		{
			dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
		}
		dist = sqrt(dist);
		if ( dist < mindist ) 
		{
			mindist = dist;
			ielmindist = iel;
		}
		pfm->el[ieltype].GetBoundBox(&vmin[0],&vmax[0]);

		//�������� �������� �� ������� ����� � ��������������, �������������� ������ �������
		//���� �������� - ����� ���� ����� ��������� ��������� ������
		fl = 0;
		for (i=0; i<pfm->KORT; i++)
		{
			if ( crd[i] < vmin[i] || crd[i] > vmax[i] )
			{
				fl = 1;
			}
		}
		if ( fl == 0 )
		{
			if ( ilistel > 0 ) //���� ��� ������� ������� �������� (�������� ����������) ����� ������� �� ����
			{//�� ��� ��������� ��������� ��������� ����������� ��������� ��������� ����������
				for (i=0; i<pfm->KORT; i++)
				{
					lcrd[i] = 0.0;
				}
			}
			glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
			
			if ( glcretflag > 0 )
			{
				//����� ����������������
				connected_el = iel;
				ncn = pfm->el[ieltype].NN;
				//ConnectedNodes = new int[ncn];
				//FF = new double[ncn];
				for (i=0; i<ncn; i++)
				{
					ConnectedNodes[i] = pfm->IND[iel][i+2];
					FF[i] = FFtmp[i];
				}
				if ( ilistel == 0 ) return(1);
				else return(2);
			}
		}

		//���������� � ������ ��� ������� ������� �������� �������������� ��������
		for (i=0; i<pfm->ELGR[iel][0]; i++)
		{
			if ( pfm->MASK_EL[ pfm->ELGR[iel][i+1] ] == 0 )
			{
				pfm->LIST_EL[ pfm->nlistel ] = pfm->ELGR[iel][i+1];
				pfm->nlistel++;
				pfm->MASK_EL[ pfm->ELGR[iel][i+1] ] = 1;
			}
		}
		ilistel++; //������� � ���������� �������� �� ������
	}
	// ���� ��������� �� ����������� �� ����� ����� - ������ ������� �� ��������� �� � ����� �� ��
	// �� �������������� ����� ������������ ��������� ���������� ����� ��� ��������, � �������� ��� ����� ����� �����������

	iel = ielmindist;
	ieltype = pfm->IND[iel][0];
	pfm->AttachElement(&pfm->el[ieltype],iel);
	for (i=0; i<pfm->KORT; i++)
	{
		lcrd[i] = 0.0;
	}
	glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
	connected_el = iel;
	ncn = pfm->el[ieltype].NN;
	//ConnectedNodes = new int[ncn];
	//FF = new double[ncn];
	for (i=0; i<ncn; i++)
	{
		ConnectedNodes[i] = pfm->IND[iel][i+2];
		FF[i] = FFtmp[i];
	}
	return(3);
}

void FDPOINT::TrackingExternPoint( void *vpfm, double *gcrd)
{
	int i,j;
	int iel,ieltype;
	double FFtmp[20];
	double vcentr[3];
	double mindist,dist;
	int ielmindist, ielmindist_old;
	int glcretflag;

	FULLMODEL *pfm;
	pfm = (FULLMODEL *)vpfm;

	for (i=0; i<pfm->KORT; i++)
	{
		crd[i] = gcrd[i];
	}

	//��������� ����������� ���������� �� � ��������� ���������
	if (connected_el == -1) AttachToFullModel(vpfm,gcrd,connected_el,&lcrd[0]);
	else
	{
		ielmindist_old = -2;
		ielmindist = connected_el;
		
		while (ielmindist_old != ielmindist) //�������� ����������� ��������� ��� �� ��� ������, ���� ������� ����� ������������� �������, ��� �� ���� ������� �������
		{//����������� ���������� ��������� ��������
			mindist = 100000.0;
			ielmindist_old = ielmindist;
			//���������� �� ��, ����������� ��������� � ������� ���
			iel = ielmindist_old;
			ieltype = pfm->IND[iel][0];
			pfm->AttachElement(&pfm->el[ieltype],iel);
			pfm->el[ieltype].GetCentr(&vcentr[0]);
			dist = 0.0;
			for (i=0; i<pfm->KORT; i++)
			{
				dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
			}
			dist = sqrt(dist);
			if ( dist < mindist ) 
			{
				mindist = dist;
				ielmindist = iel;
			}
			//����� ������������ ���������� �� �������� ��������� ����� ������� ����������� ���������� ��
			for (j=1; j< (pfm->ELGR[ielmindist_old][0] + 1); j++)
			{
				iel = pfm->ELGR[ielmindist_old][j];
				ieltype = pfm->IND[iel][0];
				pfm->AttachElement(&pfm->el[ieltype],iel);
				pfm->el[ieltype].GetCentr(&vcentr[0]);
				dist = 0.0;
				for (i=0; i<pfm->KORT; i++)
				{
					dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
				}
				dist = sqrt(dist);
				if ( dist < mindist ) 
				{
					mindist = dist;
					ielmindist = iel;
				}
			}
		}
		//��� ��������, ����������� ���������, ���������� ��������� ���������� ������� �����
		iel = ielmindist;
		ieltype = pfm->IND[iel][0];
		pfm->AttachElement(&pfm->el[ieltype],iel);
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = 0.0;
		}
		glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
		connected_el = iel;
		ncn = pfm->el[ieltype].NN;
		//ConnectedNodes = new int[ncn];
		//FF = new double[ncn];
		for (i=0; i<ncn; i++)
		{
			ConnectedNodes[i] = pfm->IND[iel][i+2];
			FF[i] = FFtmp[i];
		}
	}
}