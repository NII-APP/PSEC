#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"

void SEMODEL::AutoDivision_PrintLevels_Paraview()
{
	int i,ii,j,k,KN;
	int **ELSELEV = NULL; //дл€ каждого  Ё, показывает в какой —Ё каждого уровн€ он входит
	int NEL = pfm->NEL;
	int KORT = pfm->KORT;

	char strtmp[256];
	char strname[256];
	FILE *fp;

	ELSELEV = MM->MEM_NEW(ELSELEV,NEL,SElevelsnum);

	for (i=0; i<NEL; i++)
	{
		for (j=0; j<SElevelsnum; j++)
		{
			ii=i;
			for (k=0; k<j; k++) // можно оптимизировать и исключить этот цикл?
			{
				ii = Levels[k][ii];
			}
			ELSELEV[i][j] = Levels[j][ii]%15; //фикитивные номера с приведением к 15 цветам
		}
	}


	sprintf(strname,"%s\\SELEVELS.vtk",pfm->pathmain);

	fp = fopen(strname,"w");

	pfm->ParaView_PrintGrid(fp);

	//fprintf(fp,"# vtk DataFile Version 3.0\n");
	//fprintf(fp,"%s\n",strname);
	//fprintf(fp,"ASCII\n");
	//fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	//fprintf(fp,"POINTS %d float\n",pfm->NN);
	//for (i=0; i<pfm->NN; i++)
	//{
	//	fprintf(fp,"%f %f %f\n", pfm->CRD[KORT*i], pfm->CRD[KORT*i+1], pfm->CRD[KORT*i+2]);
	//}

	//// в данном варианте возможен только вывод моделей, состо€щих из единственного типа  Ё
	//KN = pfm->IND[0][1];

	//if ( KN == 10 )
	//{
	//	fprintf(fp,"CELLS %d %d\n",NEL,11*NEL);
	//	for (i=0; i<NEL; i++)
	//	{
	//		fprintf(fp,"10 \n");
	//		for (j=0; j<KN; j++)
	//		{
	//			fprintf(fp,"%d ",pfm->IND[i][j+2]);
	//		}
	//		fprintf(fp,"\n");
	//	}
	//	fprintf(fp,"CELL_TYPES %d\n",NEL);
	//	for (i=0; i<NEL; i++)
	//	{
	//		fprintf(fp,"24\n");
	//	}
	//}
	//if( KN == 20 )
	//{
	//	//NEL_OLD = NEL;

	//	//NEL = 1;
	//	//fprintf(fp,"CELLS %d %d\n",NEL,21*NEL);
	//	//for (i=0; i<NEL; i++)
	//	//{
	//	//	fprintf(fp,"20 \n");
	//	//	for (j=0; j<KN; j++)
	//	//	{
	//	//		fprintf(fp,"%d ",pfm->IND[i][j]);
	//	//	}
	//	//	fprintf(fp,"\n");
	//	//}
	//	//fprintf(fp,"CELL_TYPES %d\n",NEL);
	//	//for (i=0; i<NEL; i++)
	//	//{
	//	//	fprintf(fp,"25\n");
	//	//}

	//	//NEL = NEL_OLD;

	//	fprintf(fp,"CELLS %d %d\n",NEL,21*NEL);
	//	for (i=0; i<NEL; i++)
	//	{
	//		fprintf(fp,"20 \n");
	//		/*for (j=0; j<KN; j++)
	//		{
	//			fprintf(fp,"%d ",pfm->IND[i][j+2]);
	//		}*/
	//		fprintf(fp,"%d ",pfm->IND[i][2+0]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+1]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+2]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+3]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+4]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+5]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+6]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+7]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+8]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+9]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+10]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+11]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+16]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+17]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+18]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+19]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+12]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+13]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+14]);
	//		fprintf(fp,"%d ",pfm->IND[i][2+15]);

	//		fprintf(fp,"\n");
	//	}
	//	fprintf(fp,"CELL_TYPES %d\n",NEL);
	//	for (i=0; i<NEL; i++)
	//	{
	//		fprintf(fp,"25\n");
	//	}
	//}

	fprintf(fp,"CELL_DATA %d\n",NEL);
	for (j=0; j<SElevelsnum; j++)
	{
		fprintf(fp,"SCALARS Level_%d int 1\n",j);
		fprintf(fp,"LOOKUP_TABLE default\n");
		for (i=0; i<NEL; i++)
		{
			fprintf(fp,"%d\n",ELSELEV[i][j]);
		}
	}

	fclose(fp);

	ELSELEV = MM->MEM_DEL(ELSELEV,NEL,SElevelsnum);
}