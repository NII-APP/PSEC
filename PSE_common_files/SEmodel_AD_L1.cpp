#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"

void SEMODEL::AutoMLD_L1(int ilev, int maxboundary_size, int maxtotal_size, int **IND_OLD, int nel, int NNlev )
{
	int i,j,k,kk;
	int nbn,ntotaln,nin;
	int istartse = 0;
	int icurnewse = 0;
	int nincludedse = 0;
	int maxbn, imaxbn, isemaxbn,curbn,ise;
	int **GREL = NULL;
	int **NDEL = NULL;
	int **SEBND = NULL;
	int *MASK = NULL;
	int *NEIB = NULL;
	int neib;

	// ������ ����� �� ��������������
	
	MM->GRNDEL(&NDEL,IND_OLD,nel,NNlev);
	MM->GenELGR(&GREL,nel,NNlev,1,NDEL,IND_OLD);

	//����������� ����� ��������� ����� ��� ������� ������������� � ������ ���� ������������ �������������� �������� ������
	FindSeBound(&SEBND,GREL,NDEL,NNlev,nel,ilev-1);

	//������ ������������ ��:
	// ���������� ��������� ��. � �������� ������������ ����� �� � ���� �������������� �� �� ������ ��� ������� 1��� � ����� ������� �������
	// ��� ������������� ���������� � ������ ������� �� �� �� ������ �������, ������� ����� ������ ��������� ����� � ����� ����������� ��
	// ����������� ����������� �� ���������� ��������� ����� (�, ��������, �� ������ ���������� �����)
	// ����������� �� ���������� ��������� ����� ������ �� ����� ��� ����� (����� � 2.5 - 3 ����) ��������� ������������ ����� ��������� ����� ����� �� �������� ������  

	// � �������� ������� �� ����� ����� ������� �� � ���������� ����������� ��������� �����

	MASK = new int[nel];
	NEIB = new int[nel];
	for (i=0; i<nel; i++)
	{
		MASK[i] = -1; //����� ����� ���������� � ����� ����� ����� �� ��� ������� ������ ������������
	}
	neib = 0; //� NEIB ����� ��������� ������� ������ ��, ���������� ��� ������������ ������ ��  
	nincludedse = 0;
	while ( nincludedse < nel )
	{
		//������ ������������ ���������� ������ ��
		neib = 0;
		NEIB[neib] = istartse;
		neib++;
		nincludedse++;
		MASK[istartse] = icurnewse;

		nbn = 0; // ������� ���������� ��������� �����
		ntotaln = 0; //������� ������ ����� �����
		nin = 0; //������� ���������� ���������� �����
		imaxbn = 1;
		while ( nbn < maxboundary_size && ntotaln < maxtotal_size && imaxbn != -1 )
		{
			// ������ ��������� ����� � �������� ���� �������
			maxbn = 0;
			imaxbn = -1;
			
			for ( i=0; i < neib; i++)
			{
				for ( j=1; j < (GREL[ NEIB[i] ][0] + 1); j++)
				{
					ise = GREL[ NEIB[i] ][j];
					if ( MASK[ ise ] == -1 ) //���� ���� ������������ ��� ������ �� �������
					{
						//���������� ����� ��� �����, ���������� � ����������� ��
						curbn = 0;
						for (k=0; k<neib; k++)
						{
							for (kk=1; kk<(GREL[ise][0]+1); kk++)
							{
								if ( GREL[ise][kk] == NEIB[k] )
								{
									curbn += SEBND[ise][kk];
								}
							}
						}
						if ( curbn > maxbn )
						{
							maxbn = curbn;
							imaxbn = ise;
						}
					}
				}
			}
			if ( imaxbn != -1 ) //�� �������� ����������������� �� � ������� ��������� � ����� ����������� ��
			{
				// ��������� �� �������� � NEIB � ���������������� ����� ��������� ����� ������ ��
				NEIB[neib] = imaxbn;
				neib++;
				nincludedse++;
				MASK[imaxbn] = icurnewse;

				nbn = 0;
				ntotaln = 0;
				nin = 0;
				//������ ���������� ��������� ����� �������������� ���:
				// ������� ����������� ����� �������� ���������� ��������� ����� ����� ������ ����� �������� ��
				for (i=0; i<neib; i++)
				{
					for (j=i+1; j<neib; j++)
					{
						for (k=1; k<(GREL[ NEIB[i] ][0] + 1); k++)
						{
							if ( GREL[ NEIB[i] ][k] == NEIB[j] )
							{
								nin += SEBND[ NEIB[i] ][k];
							}
						}
					}
					ntotaln += SEBND[ NEIB[i] ][0];
				}
				ntotaln -= nin;
				nbn = ntotaln - nin;
			}
		}

		// ����� ����, ��� ����� ������������ ������, ���� ����� ��������� ��
		istartse = -1;
		
		while ( istartse == -1 && nincludedse < nel )
		{
			//������� ��������� ������ ��, �������� � ������� ��, �������� � ������ ���������� ��������������� �� ������ ������
			for (i=0; i<neib; i++)
			{
				for (k=1; k<(GREL[ NEIB[i] ][0] + 1); k++)
				{
					if ( MASK[ GREL[ NEIB[i] ][k] ] == -1 && istartse == -1 )
					{
						istartse = GREL[ NEIB[i] ][k];
					}
				}
			}
			//���� ����� ������ ������� �� ������� �� ������ ������������������ �������������
			// ��������� ������ ����� ������ ����������������� �������
			if (istartse == -1)
			{
				for (i=0; i<nel; i++)
				{
					if ( MASK[i] == -1 ) break;
				}
				if ( i < nel ) //����������������� ��� ��������
				{
					istartse = i;
				}
			}

			//��������� ���������, �� �������� �� ����� ��������� �� �������������
			if (istartse != -1)
			{
				for (j=1; j < (GREL[istartse][0] + 1); j++)
				{
					if ( MASK[ GREL[istartse][j] ] == -1 ) break;
				}
				if ( j == (GREL[istartse][0] + 1) )
				{
					//���� ������� ����������, ��������� ��� ������������ � ������ �� �������� � ��� ����� ��������������
					// ����� ����������������� - ������������ ���, ����� ���������� ����������� ����� ��������� �����
					// � ������ ������ ������������ � ���� �� � ������� ������ ����� ��������� �����.
					curbn = 0;
					imaxbn = -1;
					for ( k=1; k < (GREL[istartse][0]+1); k++)
					{
						if ( SEBND[istartse][k] > curbn )
						{
							curbn = SEBND[istartse][k];
							imaxbn = k;
						}
					}
					//������������ i ������������ � ���� �� �������������, � �������� ����������� GREL[istartse][k]
					MASK[ istartse ] = MASK[ GREL[istartse][imaxbn] ];
					nincludedse++;
					istartse = -1; //������ ��������� ����� ������ ���������� ��
				}
			}
		}

		// � ������� ������ �������������� ��� ��
		neib = 0;

		icurnewse++;
	}

	//� icurnewse ���������� ����� �������������� �������������� ������ ������
	numSEbylevels[ilev] = icurnewse;
	//����� ��������� � ������� � ���������� ������
	Levels[ilev] = new int[nel];
	for (i=0; i<nel; i++)
	{
		Levels[ilev][i] = MASK[i];
	}
	
	delete []NEIB;
	delete []MASK;
	SEBND = MM->MEM_DEL(SEBND,nel,1);
	GREL = MM->MEM_DEL(GREL,nel,1);
	NDEL = MM->MEM_DEL(NDEL,NNlev,1);

}

		//if ( neib == 1 ) //��������� �� �������� ������������ ����������������� � ������� �������. ��� ��������� ������������ � ������ ������ ��������������� ��������� ��
		//{
		//	for (k=1; k<(GREL[ NEIB[0] ][0] + 1); k++)
		//	{
		//		if ( MASK[ GREL[ NEIB[0] ][k] ] != -1 ) break;
		//	}
		//	if ( k != (GREL[ NEIB[0] ][0] + 1) )
		//	{
		//		MASK[ NEIB[0] ] = MASK[ GREL[ NEIB[0] ][k] ];
		//	}
		//	else
		//	{
		//		k=k; //������, ��� ��������������� �������� ��
		//	}
		//}

void SEMODEL::FindSeBound(int ***SEBND_, int **GREL, int **NDEL, int NNlev, int nel, int ilev)
{
	int i,j,k,in,ise1,ise2;
	int **SEBND=NULL;

	//������ ���������� � ��������� ����� ������������� �� ��������� � ������ �������������� ��������� � ������ �� ��
	SEBND = new int *[nel];
	for (i=0; i<nel; i++)
	{
		SEBND[i] = new int [GREL[i][0]+1]; // ������ ����� � ������ - ����� ����� ��������� ����� �������������
		for (j=0; j< (GREL[i][0]+1); j++)
		{
			SEBND[i][j] = 0;
		}
		SEBND[i][0] = SE[ilev][i].NS;
	}

	for (in = 0; in < NNlev; in++)
	{
		//��� ������� ���� �������������� ������ ������ ������� �� �� NDEL
		//��� ������ ���� �� �� ������ ��������� ����������� 1 � ��������������� ������ SEBND

		for ( i = 2; i < (NDEL[in][0]+2); i++)
		{
			for (j = i+1; j < (NDEL[in][0]+2); j++)
			{
				ise1 = NDEL[in][i];
				ise2 = NDEL[in][j];
				
				for ( k = 1; k < (GREL[ise1][0]+1); k++)
				{
					if ( GREL[ise1][k] == ise2 ) break;
				}
				if ( k == (GREL[ise1][0]+1) )
				{
					k=k; //������: ������������ �� ������ � �����, ���� ����� ����� ���� � ������� ��
				}
				else
				{
					SEBND[ise1][k]++;
				}
				//���������� ������� ��������
				for ( k = 1; k < (GREL[ise2][0]+1); k++)
				{
					if ( GREL[ise2][k] == ise1 ) break;
				}
				if ( k == (GREL[ise2][0]+1) )
				{
					k=k; //������: ������������ �� ������ � �����, ���� ����� ����� ���� � ������� ��
				}
				else
				{
					SEBND[ise2][k]++;
				}
			}
		}
	}
	*SEBND_ = SEBND;
}
