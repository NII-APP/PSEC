#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"


void SEMODEL::AutoMLD_L0_EL(int *DS_)
{
	//������� �� ������������� �������� ������
	//������� �� ��������� ������� � ���, ��� ����� �� �������������� ������� � ������ �� ����� �� ������ ������� ����� ������

	int i,j,k, NN, NEL, numNodes, startnode, neib,node;
	int **GR = NULL,**NDEL = NULL, **IND = NULL;
	int *MASK;
	int *NEIB;
	int *DEG;
	int DS = *DS_;

	NN = pfm->NN;
	IND = pfm->IND;
	NEL = pfm->NEL;

	// 1 -  ����������� ������� ����� �� ������
	MM->GRNDEL(&NDEL,IND,NEL,NN);
	MM->GenELGR(&GR,NEL,NN,4,NDEL,IND);
	
	//�������� �� ������� ���������, �� ������� ����� ������ � ������� ����������
	int fl = 0;
	for (i=0; i<NEL; i++)
	{
		if ( GR[i][0] == 0 ) fl++;
	}
	if (fl > 0)
	{
		printf("\nthere is %d elements without neibourghs, exiting......\n",fl);
		exit(-4);
	}

	// 2 - ����������� ����� (��) ������� ������ 
	
	MASK = new int [NEL];
	DEG = new int [NEL];
	NEIB = new int [NEL];
	for (i=0; i<NEL; i++)
	{
		MASK[i] = -1;
		DEG[i] = 0;
	}

	if ( DS > (1.05*NEL/2) ) DS = (int)(1.05*NEL/2);

	// 2.1 - ����� ���������� ��������� ��������
	startnode = MM->FindStartingNode(pfm->CRD,NN,pfm->KORT);
	startnode = NDEL[startnode][2];

	// ����������� �����
	int iSE = 0;
	numNodes = 0;
	while ( numNodes < NEL )
	{
		printf("\niSE = %d starting node determination ....",iSE);
		// ����� ���������� ���� ��� ���� ��, ����� �������
		if (iSE > 0)
		{
			//�������� ������������ � ������� ��� ���� �� ����� ������ 
			for (i = neib-1; i>=0; i--)
			{
				node = NEIB[i];
				// �������� �� ���� ��� ������� � �������� ��� �������� ����, ������� ��� �� ��������� (�� ����������� �������� ��)
				for (j = 1; j< (GR[node][0]+1); j++)
				{
					if ( MASK[ GR[node][j] ] == -1 )
					{
						i = -99;  //!!!!!
						break;
					}
				}
			}

			if (i == -100) //���� ������ !!!!!!!!!!!!!!!
			{
				startnode = GR[node][j];
			}
			else
			{
				// ��� ����, ������������� � �������� �� ��� �������� � ������ ��
				// �������� ������ ���������� ����������������� ����
				for (i=0; i<NEL; i++)
				{
					if ( MASK[i] == -1 )
					{
						break;
					}
				}
				startnode = i;
			}
		}

		printf("\nstartnode = %d ",startnode);

		//��������� ������������ �������� �� �������������� ����� � ������������ �� �������
		neib = 0;
		//FindNotNumNeib(GR, MASK, iSE, NEIB, &neib, DS, startnode);
		FindNotNumNeib(GR, MASK, DEG, iSE, NEIB, &neib, DS, startnode);
		printf("\nSE prepared (neib = %d) ",neib);

		//������ ����������� ��������. ���� ������ ������ ������, ��� ��������� ������������ � ������ �� ���������� � ��� ��
		// � ��������� ������ ��������� ����� �� � ��� ����� ��� ��� ���������� ����� ��������� � MASK

		if ( neib > (0.3*DS) )
		{
			printf("\nSize normal, simple copy.... ");
			for (i=0; i<neib; i++)
			{
				MASK[ NEIB[i] ] = iSE;
				numNodes++;
			}
			iSE++;
		}
		else
		{
			printf("\nSize small, find neighbour SE.... ");
			//��������� ���������������� � ������ �� �������� ������ ��, � ������������ ��� � ����, � ������� ������ �������� �����
			int *nNodesNeibToSe,isupel;
			int nbse = 0, maxneib,imax;
			nNodesNeibToSe = new int[iSE+1]; //������ - ����� ��, ����� �� ������� - ���������� ����� ����� ��, ������������� �� ����� ����������� ��
			//���� �������� 0 - ��������� ���

			for (i=0; i<(iSE+1); i++)
			{
				nNodesNeibToSe[i] = 0;
			}

			for (i=0; i< neib; i++)
			{
				for (j=1; j< (GR [ NEIB[i] ][0] + 1); j++)
				{
					isupel = MASK[ GR[ NEIB[i] ][j] ];
					if ( isupel >= 0 )
					{
						nNodesNeibToSe[isupel]++;
					}
				}
			}
			maxneib = 0;
			imax = -1;
			for (i=0; i<(iSE+1); i++)
			{
				if ( maxneib < nNodesNeibToSe[i] )
				{
					imax = i;
					maxneib = nNodesNeibToSe[i];
				}
			}
			// � ���������� imax - ����� ���� ��, � �������� ����� ������������ ��������� ������ �����
			// �������� �� ��� ��, ��� � ���������� ������, �� � ���� ��������� ������� "������������" ��
			// ���������� iSE � ������ ������ �� �������������, ��� ��� ����� �� �� ������

			if (imax >= 0)
			{
				for (i=0; i<neib; i++)
				{
					MASK[ NEIB[i] ] = imax;
					numNodes++;
				}
			}
			else
			{
				//������!!! �� ������� �������� ��, � �� ��������� ���������� �������. �������� �������� ������, �� �������� ����������� ������?
				i=i;
				printf("\n Can not connect small SE, there is no neibourgh SEs: FEM model is not consistent.... \n");
			}

			delete [] nNodesNeibToSe;

			printf("\nneibourgh SE = %d    neibnodes = %d.... \n\n", imax, maxneib);
		}
	}
	
	Levels[0] = new int[NEL];
	for (i=0; i<NEL; i++)
	{
		Levels[0][i] = MASK[i];
	}

	numSEbylevels[0] = iSE;

	MASK = MM->MEM_DEL(MASK,NEL);
	NEIB = MM->MEM_DEL(NEIB,NEL);
	DEG = MM->MEM_DEL(DEG,NEL);
	GR = MM->MEM_DEL(GR,NEL,1);
	NDEL = MM->MEM_DEL(NDEL,NN,1);

	*DS_ = DS;
}


void SEMODEL::AutoMLD_L0(int *DS_)
{
	//������� �� ������������� �������� ������
	//������� �� ��������� ������� � ���, ��� ����� �� �������������� ������� � ������ �� ����� �� ������ ������� ����� ������

	int i,j,k, NN, NEL, numNodes, startnode, neib,node;
	int **GR = NULL,**NDEL = NULL, **IND = NULL;
	int *MASK;
	int *NEIB;
	int *DEG;
	int DS = *DS_;

	NN = pfm->NN;
	IND = pfm->IND;
	NEL = pfm->NEL;

	// 1 -  ����������� ������� ����� �� ������
	MM->GRNDEL(&NDEL,IND,NEL,NN);
	MM->GenNDGR(&GR,NEL,NN,NN,NDEL,IND);
	NDEL = MM->MEM_DEL(NDEL,NN,1);

	// 2 - ����������� ����� (��) ������� ������ 
	
	MASK = new int [NN];
	DEG = new int [NN];
	NEIB = new int [NN];
	for (i=0; i<NN; i++)
	{
		MASK[i] = -1;
		DEG[i] = 0;
	}

	if ( DS > (1.05*NN/2) ) DS = (int)(1.05*NN/2);

	// 2.1 - ����� ������� ����, ��� ���� � ���������� ����������� �� ���� ��������������� ������ ��������������� � ������������ ������������ �� ���� 3� ����

	startnode = MM->FindStartingNode(pfm->CRD,NN,pfm->KORT);

	// ����������� �����
	int iSE = 0;
	numNodes = 0;
	while ( numNodes < NN )
	{
		printf("\niSE = %d starting node determination ....",iSE);
		// ����� ���������� ���� ��� ���� ��, ����� �������
		if (iSE > 0)
		{
			//�������� ������������ � ������� ��� ���� �� ����� ������ 
			for (i = neib-1; i>=0; i--)
			{
				node = NEIB[i];
				// �������� �� ���� ��� ������� � �������� ��� �������� ����, ������� ��� �� ��������� (�� ����������� �������� ��)
				for (j = 1; j< (GR[node][0]+1); j++)
				{
					if ( MASK[ GR[node][j] ] == -1 )
					{
						i = -99;  //!!!!!
						break;
					}
				}
			}

			if (i == -100) //���� ������ !!!!!!!!!!!!!!!
			{
				startnode = GR[node][j];
			}
			else
			{
				// ��� ����, ������������� � �������� �� ��� �������� � ������ ��
				// �������� ������ ���������� ����������������� ����
				for (i=0; i<NN; i++)
				{
					if ( MASK[i] == -1 )
					{
						break;
					}
				}
				startnode = i;
			}
		}

		printf("\nstartnode = %d ",startnode);

		//��������� ������������ �������� �� �������������� ����� � ������������ �� �������
		neib = 0;
		//FindNotNumNeib(GR, MASK, iSE, NEIB, &neib, DS, startnode);
		FindNotNumNeib(GR, MASK, DEG, iSE, NEIB, &neib, DS, startnode);
		printf("\nSE prepared (neib = %d) ",neib);

		//������ ����������� ��������. ���� ������ ������ ������, ��� ��������� ������������ � ������ �� ���������� � ��� ��
		// � ��������� ������ ��������� ����� �� � ��� ����� ��� ��� ���������� ����� ��������� � MASK

		if ( neib > (0.3*DS) )
		{
			printf("\nSize normal, simple copy.... ");
			for (i=0; i<neib; i++)
			{
				MASK[ NEIB[i] ] = iSE;
				numNodes++;
			}
			iSE++;
		}
		else
		{
			printf("\nSize small, find neighbour SE.... ");
			//��������� ���������������� � ������ �� �������� ������ ��, � ������������ ��� � ����, � ������� ������ �������� �����
			int *nNodesNeibToSe,isupel;
			int nbse = 0, maxneib,imax;
			nNodesNeibToSe = new int[iSE+1]; //������ - ����� ��, ����� �� ������� - ���������� ����� ����� ��, ������������� �� ����� ����������� ��
			//���� �������� 0 - ��������� ���

			for (i=0; i<(iSE+1); i++)
			{
				nNodesNeibToSe[i] = 0;
			}

			for (i=0; i< neib; i++)
			{
				for (j=1; j< (GR [ NEIB[i] ][0] + 1); j++)
				{
					isupel = MASK[ GR[ NEIB[i] ][j] ];
					if ( isupel >= 0 )
					{
						nNodesNeibToSe[isupel]++;
					}
				}
			}
			maxneib = 0;
			imax = -1;
			for (i=0; i<(iSE+1); i++)
			{
				if ( maxneib < nNodesNeibToSe[i] )
				{
					imax = i;
					maxneib = nNodesNeibToSe[i];
				}
			}
			// � ���������� imax - ����� ���� ��, � �������� ����� ������������ ��������� ������ �����
			// �������� �� ��� ��, ��� � ���������� ������, �� � ���� ��������� ������� "������������" ��
			// ���������� iSE � ������ ������ �� �������������, ��� ��� ����� �� �� ������

			if (imax >= 0)
			{
				for (i=0; i<neib; i++)
				{
					MASK[ NEIB[i] ] = imax;
					numNodes++;
				}
			}
			else
			{
				//������!!! �� ������� �������� ��, � �� ��������� ���������� �������. �������� �������� ������, �� �������� ����������� ������?
				i=i;
			}
			delete [] nNodesNeibToSe;

			printf("\nneibourgh SE = %d    neibnodes = %d.... \n\n", imax, maxneib);
		}
	}

	NEIB = MM->MEM_DEL(NEIB,NN);
	GR = MM->MEM_DEL(GR,NN,1);

	FILE *fp;
	char strname[256];
	int KORT = pfm->KORT,KN;
	sprintf(strname,"%s\\SELEVELS_0_NODE.vtk",pfm->pathmain);

	fp = fopen(strname,"w");

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"%s\n",strname);
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d float\n",pfm->NN);
	for (i=0; i<pfm->NN; i++)
	{
		fprintf(fp,"%f %f %f\n", pfm->CRD[KORT*i], pfm->CRD[KORT*i+1], pfm->CRD[KORT*i+2]);
	}

	// � ������ �������� �������� ������ ����� �������, ��������� �� ������������� ���� ��
	KN = pfm->IND[0][1];

	if ( KN == 10 )
	{
		fprintf(fp,"CELLS %d %d\n",NEL,11*NEL);
		for (i=0; i<NEL; i++)
		{
			fprintf(fp,"10 \n");
			for (j=0; j<KN; j++)
			{
				fprintf(fp,"%d ",pfm->IND[i][j+2]);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"CELL_TYPES %d\n",NEL);
		for (i=0; i<NEL; i++)
		{
			fprintf(fp,"24\n");
		}
	}
	if( KN == 20 )
	{
		//NEL_OLD = NEL;

		//NEL = 1;
		//fprintf(fp,"CELLS %d %d\n",NEL,21*NEL);
		//for (i=0; i<NEL; i++)
		//{
		//	fprintf(fp,"20 \n");
		//	for (j=0; j<KN; j++)
		//	{
		//		fprintf(fp,"%d ",pfm->IND[i][j]);
		//	}
		//	fprintf(fp,"\n");
		//}
		//fprintf(fp,"CELL_TYPES %d\n",NEL);
		//for (i=0; i<NEL; i++)
		//{
		//	fprintf(fp,"25\n");
		//}

		//NEL = NEL_OLD;

		fprintf(fp,"CELLS %d %d\n",NEL,21*NEL);
		for (i=0; i<NEL; i++)
		{
			fprintf(fp,"20 \n");
			/*for (j=0; j<KN; j++)
			{
				fprintf(fp,"%d ",pfm->IND[i][j+2]);
			}*/
			fprintf(fp,"%d ",pfm->IND[i][2+0]);
			fprintf(fp,"%d ",pfm->IND[i][2+1]);
			fprintf(fp,"%d ",pfm->IND[i][2+2]);
			fprintf(fp,"%d ",pfm->IND[i][2+3]);
			fprintf(fp,"%d ",pfm->IND[i][2+4]);
			fprintf(fp,"%d ",pfm->IND[i][2+5]);
			fprintf(fp,"%d ",pfm->IND[i][2+6]);
			fprintf(fp,"%d ",pfm->IND[i][2+7]);
			fprintf(fp,"%d ",pfm->IND[i][2+8]);
			fprintf(fp,"%d ",pfm->IND[i][2+9]);
			fprintf(fp,"%d ",pfm->IND[i][2+10]);
			fprintf(fp,"%d ",pfm->IND[i][2+11]);
			fprintf(fp,"%d ",pfm->IND[i][2+16]);
			fprintf(fp,"%d ",pfm->IND[i][2+17]);
			fprintf(fp,"%d ",pfm->IND[i][2+18]);
			fprintf(fp,"%d ",pfm->IND[i][2+19]);
			fprintf(fp,"%d ",pfm->IND[i][2+12]);
			fprintf(fp,"%d ",pfm->IND[i][2+13]);
			fprintf(fp,"%d ",pfm->IND[i][2+14]);
			fprintf(fp,"%d ",pfm->IND[i][2+15]);

			fprintf(fp,"\n");
		}
		fprintf(fp,"CELL_TYPES %d\n",NEL);
		for (i=0; i<NEL; i++)
		{
			fprintf(fp,"25\n");
		}
	}

	fprintf(fp,"POINT_DATA %d\n",pfm->NN);
	fprintf(fp,"SCALARS Level0node int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<pfm->NN; i++)
	{
		fprintf(fp,"%d\n",MASK[i]);
	}

	fclose(fp);



	//�� ����������� �� �� �����, � �� ��. ������� �������������� ������ ����� ������������� � ������ �� ���������. 
	// ���� ������� �������� ����, ����������� � ������ ��, �� ������� ������� � ���� ��, ��� � ���� ������ �����
	int iLev,NSE,isupel,imax;
	int *relSE, *numrelSE;
	int irse,tmpi;
	double tmp;

	
	Levels[0] = new int[NEL];
	for (i=0; i<NEL; i++)
	{
		Levels[0][i] = -1;
	}

	numSEbylevels[0] = iSE;

	iSE = 0;
	iLev = 0;
	NSE = numSEbylevels[0];

	relSE = new int[NSE];
	irse = 0;
	numrelSE = new int[NSE];
	for(i=0; i<NSE; i++)
	{
		numrelSE[i] = 0;
	}

	printf("\nreforming SE model from nodes to elements formulation..... NEL = %d \n", NEL);

	for (i=0; i<NEL; i++)
	{
		if ( (i%1000) == 0 )
		{
			tmp = 100.0*(((double)i)/((double)NEL));
			tmpi = (int)tmp;
			
			printf(" SE model reform %d%% IEL %d\r",tmpi,i);
		}
		for (j=2; j< (IND[i][1]+2); j++)
		{
			isupel = MASK[ IND[i][j] ];
			if (numrelSE[isupel] == 0)
			{
				relSE[irse] = isupel;
				irse++;
			}
			numrelSE[isupel]++;
		}
		if (irse == 0)
		{
			//!!!!������, ������� �� ��������� �� � ������ ��
			i=i;
		}
		imax = 0;
		isupel = -1;
		for(j=0; j<irse; j++)
		{
			if ( imax < numrelSE[ relSE[j] ] )
			{
				imax = numrelSE[ relSE[j] ];
				isupel = relSE[j];
			}
		}

		Levels[0][i] = isupel;

		//����������� ������� ��������������� ��������
		for(j=0; j<irse; j++)
		{
			numrelSE[ relSE[j] ] = 0;
		}
		irse = 0;
	}

	relSE = MM->MEM_DEL(relSE,NSE);
	numrelSE = MM->MEM_DEL(numrelSE,NSE);
	MASK = MM->MEM_DEL(MASK,NN);


	*DS_ = DS;


	////������: ����� ������� � ���������  ������� � ��������� ����, � ����� �������� ����� ��� Paraview

	//printf("\n\n print  Level001.sba  .....");

	//FILE *fp;
	//char str[256];
	//sprintf(str,"%s\\Level001.sba",pfm->pathmain);
	//fp = fopen(str,"wb");
	//for (i=0; i<NEL; i++)
	//{
	//	Levels[0][i] ++;
	//}
	//fwrite(Levels[0],sizeof(int),NEL,fp);
	//for (i=0; i<NEL; i++)
	//{
	//	Levels[0][i] --;
	//}
	//fclose(fp);

}

//void SEMODEL::FindNotNumNeib(int **GR, int *MASK, int iSE, int *NEIB, int *neib, int DS, int startnode)
//{
//	int i,j,k,nb,nban,ib;
//
//	NEIB[0] = startnode;
//	MASK[ NEIB[0] ] = -2;
//	nb = 1;
//
//	for (ib = 0; ib < nb; ib++)
//	{
//		for (j=1; j< (GR[ NEIB[ib] ][0] + 1); j++)
//		{
//			if ( MASK [ GR[ NEIB[ib] ][j] ] == -1 ) //���� ���� ������ �� ���������
//			{
//				NEIB[nb] = GR[ NEIB[ib] ][j];
//				MASK [ NEIB[nb] ] = -2; //���� ���������� ���������� � ����� ������������
//				nb++;
//			}
//		}
//		if ( nb > DS )
//		{
//			break;
//		}
//	}
//	
//	*neib = nb;
//
//}


//void SEMODEL::FindNotNumNeib(int **GR, int *MASK, int *DEG, int iSE, int *NEIB, int *neib, int DS, int startnode)
//{
//	int i,j,k,nb,ib,node;
//	int *PSN,ipsn_start,npsn,ipsn;
//	int mindeg, imindeg, loc_mindeg, iloc_mindeg;
//
//	//�������� ������������ � ������������� � ������ ������� ����, ������� ����������� ���������� �������, �� ����������� � ������� ��
//	//��� ������� ���� � ��� ������� ������� ������ ����������� ���������� ������� DEG, ������������ ���������� ������� ������� ����, �� ������������� ������� ��
//
//	PSN = new int[pfm->NN];
//	npsn = 0;
//	ipsn_start = 0;
//	for (ipsn = 0; ipsn< pfm->NN; ipsn++)
//	{
//		PSN[ipsn] = 0;
//	}
//
//	NEIB[0] = startnode;
//	MASK[ NEIB[0] ] = -2;
//	nb = 1;
//	DEG[ NEIB[0] ] = GR[ NEIB[0] ][0]; //��� ������ ������� ���� ���� ����������� ������ �� ��� �� �������������
//	
//	mindeg = 100000000;
//	imindeg = -1;
//
//	if (iSE == 10)
//	{
//		node = 0;
//		for (i=0; i< pfm->NN; i++)
//		{
//			if (MASK[i] == -1) node++;
//		}
//		node = 0;
//		for (i=0; i< pfm->NN; i++)
//		{
//			if (MASK[i] == -2) node++;
//		}
//		node = 0;
//		for (i=0; i< pfm->NN; i++)
//		{
//			if (MASK[i] == -3) node++;
//		}
//		node = 0;
//		for (i=0; i< pfm->NN; i++)
//		{
//			if (MASK[i] >= 0) node++;
//		}
//	}
//
//
//
//	for (ib = 1; ib < DS; ib++)
//	{
//		if (ib == 37 && iSE == 2)
//		{
//			ib = ib;
//		}
//		imindeg = -1;
//		for (i=1; i<(GR[ NEIB[ib-1] ][0] + 1); i++)
//		{
//			node = GR[ NEIB[ib-1] ][i]; //��������� ����� ���������� ��������������� ����
//			//��������� �������� ��� �����-������� ���������� ��������������� ����
//			DEG[node] = GR[ node ][0];
//			for (j=1; j<(GR[ node ][0] + 1); j++) //���� ��� ������ ����������� ������ ������� ��, ������� �����������
//			{
//				if ( MASK[ GR[node][j] ] == -2 ) DEG[node]--; 
// 			}
//			//����� ����� ���������� �� �������������
//			if ( MASK[ node ] == -1 )
//			{
//				PSN[npsn] = node;
//				npsn++;
//				MASK[ node ] = -3; //���� ���������� �� �������������
//			}
//			//���������� ������ ���� � ����������� �������� (������ ����)
//			if ( DEG[node] > 0 && MASK[node] == -3 ) 
//			{
//				if ( DEG[node] <= mindeg )
//				{
//					mindeg = DEG[node];
//					imindeg = node;
//				}
//			}
//		}
//
//		//���� ���� � ����������� �������� �� ��� ������ ����� ������� ���������� ��������������� ����
//		//��������� ������������� ������ PSN ��� ����������� ���� � ����������� �������� � ������ �������� �����
//		if ( imindeg == -1 && ipsn_start != npsn ) 
//		{
//			loc_mindeg = 100000000;
//			iloc_mindeg = -1;
//			for (ipsn = ipsn_start; ipsn<npsn; ipsn++)
//			{
//				if ( MASK[ PSN[ipsn] ] == -3 )
//				{
//					if ( DEG[ PSN[ipsn] ] < loc_mindeg )
//					{
//						loc_mindeg = DEG[ PSN[ipsn] ];
//						iloc_mindeg = PSN[ipsn];
//					}
//					if ( DEG[ PSN[ipsn] ] <= mindeg ) 
//					{
//						mindeg = DEG[ PSN[ipsn] ];
//						imindeg = PSN[ipsn];
//						break; //��������� ������� � ���� ����� ��������� �� ������ �����������
//					}
//				}
//			}
//			if ( ipsn == npsn )
//			{
//				mindeg = loc_mindeg;
//				imindeg = iloc_mindeg;
//			}
//		}
//
//		if (imindeg == -1)
//		{
//			//���� �����������, ��� �������������
//			break;
//		}
//		//������������� ���� � ����������� ��������
//		if (MASK[imindeg] == -1 || MASK[imindeg] == -3)
//		{
//			NEIB[nb] = imindeg;
//			nb++;
//			MASK[ imindeg ] = -2;
//		}
//		else
//		{
//			nb=nb;
//		}
//
//		//�������� ������ ������� PSN �� ���� �������������� ����������
//		if ( ipsn_start < npsn )
//		{
//			while ( MASK[ PSN[ ipsn_start ] ] == -2  ||  MASK[ PSN[ ipsn_start ] ] >=0 )
//			{
//				ipsn_start++;
//				if (ipsn_start == npsn) break;
//			}
//		}
//
//	}
//
//	//����������� ���� -1 ���� ����� - ����������, �� ���������� � ������� ��
//	for (ipsn = ipsn_start; ipsn<npsn; ipsn++)
//	{
//		if ( MASK[ PSN[ipsn] ] == -3 ) MASK[ PSN[ipsn] ] = -1;
//	}
//	
//	*neib = nb;
//
//	delete []PSN;
//}



void SEMODEL::FindNotNumNeib(int **GR, int *MASK, int *DEG, int iSE, int *NEIB, int *neib, int DS, int startnode)
{
	int i,j,k,nb,ib,node;
	int *PSN,ipsn_start,npsn,ipsn;
	int mindeg, imindeg, loc_mindeg, iloc_mindeg;

	//�������� ������������ � ������������� � ������ ������� ����, ������� ����������� ���������� �������, �� ����������� � ������� ��
	//��� ������� ���� � ��� ������� ������� ������ ����������� ���������� ������� DEG, ������������ ���������� ������� ������� ����, �� ������������� ������� ��

	PSN = new int[10*DS];
	npsn = 0;
	ipsn_start = 0;
	for (ipsn = 0; ipsn< 10*DS; ipsn++)
	{
		PSN[ipsn] = 0;
	}

	NEIB[0] = startnode;
	MASK[ NEIB[0] ] = -2;
	nb = 1;
	DEG[ NEIB[0] ] = GR[ NEIB[0] ][0]; //��� ������ ������� ���� ���� ����������� ������ �� ��� �� �������������
	
	mindeg = -100000000;
	imindeg = -1;

	printf("\n");
	for (ib = 1; ib < DS; ib++)
	{
		if ( ib%100 == 0)
		{
			printf("%d%%\r",(int)(ib*100.0/DS));
		}
		imindeg = -1;
		for (i=1; i<(GR[ NEIB[ib-1] ][0] + 1); i++)
		{
			node = GR[ NEIB[ib-1] ][i]; //��������� ����� ���������� ��������������� ����
			//��������� �������� ��� �����-������� ���������� ��������������� ����
			DEG[node] = GR[ node ][0];
			for (j=1; j<(GR[ node ][0] + 1); j++) //���� ��� ������ �� ����������� ������ ������� ��, ������� �����������
			{
				if ( MASK[ GR[node][j] ] != -2 ) DEG[node]--; 
 			}
			//����� ����� ���������� �� �������������
			if ( MASK[ node ] == -1 )
			{
				PSN[npsn] = node;
				npsn++;
				MASK[ node ] = -3; //���� ���������� �� �������������
			}
			//���������� ������ ���� � ������������ ��������
			if ( MASK[node] == -3 ) 
			{
				if ( DEG[node] >= mindeg )
				{
					mindeg = DEG[node];
					imindeg = node;
				}
			}
		}

		//���� ���� � ����������� �������� �� ��� ������ ����� ������� ���������� ��������������� ����
		//��������� ������������� ������ PSN ��� ����������� ���� � ����������� �������� � ������ �������� �����
		if ( imindeg == -1 && ipsn_start != npsn ) 
		{
			loc_mindeg = -100000000;
			iloc_mindeg = -1;
			for (ipsn = ipsn_start; ipsn<npsn; ipsn++)
			{
				if ( MASK[ PSN[ipsn] ] == -3 )
				{
					if ( DEG[ PSN[ipsn] ] > loc_mindeg )
					{
						loc_mindeg = DEG[ PSN[ipsn] ];
						iloc_mindeg = PSN[ipsn];
					}
					if ( DEG[ PSN[ipsn] ] >= mindeg ) 
					{
						mindeg = DEG[ PSN[ipsn] ];
						imindeg = PSN[ipsn];
						break; //��������� ������� � ���� ����� ��������� �� ������ �����������
					}
				}
			}
			if ( ipsn == npsn )
			{
				mindeg = loc_mindeg;
				imindeg = iloc_mindeg;
			}
		}

		if (imindeg == -1)
		{
			//���� �����������, ��� �������������
			break;
		}
		//������������� ���� � ������������ ��������
		if (MASK[imindeg] == -1 || MASK[imindeg] == -3)
		{
			NEIB[nb] = imindeg;
			nb++;
			MASK[ imindeg ] = -2;
		}
		else
		{
			nb=nb;
		}

		//�������� ������ ������� PSN �� ���� �������������� ����������
		if ( ipsn_start < npsn )
		{
			while ( MASK[ PSN[ ipsn_start ] ] == -2  ||  MASK[ PSN[ ipsn_start ] ] >=0 )
			{
				ipsn_start++;
				if (ipsn_start == npsn) break;
			}
		}

	}

	//����������� ���� -1 ���� ����� - ����������, �� ���������� � ������� ��
	for (ipsn = ipsn_start; ipsn<npsn; ipsn++)
	{
		if ( MASK[ PSN[ipsn] ] == -3 ) MASK[ PSN[ipsn] ] = -1;
	}
	
	*neib = nb;

	delete []PSN;
}
