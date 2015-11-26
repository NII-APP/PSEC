#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"

void FULLMODEL::ParaView_PrintCRD(FILE *fp)
{
	int i;

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"%s\n",pathmain);
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d float\n",NN);
	if (KORT == 3)
	{
		for (i=0; i<NN; i++)
		{
			fprintf(fp,"%f %f %f\n", CRD[KORT*i], CRD[KORT*i+1], CRD[KORT*i+2]);
		}
	}
	if (KORT == 2)
	{
		for (i=0; i<NN; i++)
		{
			fprintf(fp,"%f %f 0.000000\n", CRD[KORT*i], CRD[KORT*i+1]);
		}
	}
}

void FULLMODEL::ParaView_PrintIND(FILE *fp)
{
	int i,j,k,cellsize;
	// расчет количества чисел в поле CELLS
	// (для каждого элемента суммируется его число узлов + 1
	cellsize = 0;
	for (i=0; i < NEL; i++)
	{
		cellsize += (IND[i][1] + 1);
	}
	//вывод матрицы индексов
	fprintf(fp,"CELLS %d %d\n",NEL,cellsize);
	for (i=0; i<NEL; i++)
	{
		fprintf(fp,"%d ",IND[i][1]);
		if ( IND[i][0] == 24 || IND[i][0] == 5 )
		{
			for (j=0; j<IND[i][1]; j++)
			{
				fprintf(fp,"%d ",IND[i][j+2]);
			}
		}
		if ( IND[i][0] == 25 )
		{
			fprintf(fp,"%d ",IND[i][2+0]);
			fprintf(fp,"%d ",IND[i][2+1]);
			fprintf(fp,"%d ",IND[i][2+2]);
			fprintf(fp,"%d ",IND[i][2+3]);
			fprintf(fp,"%d ",IND[i][2+4]);
			fprintf(fp,"%d ",IND[i][2+5]);
			fprintf(fp,"%d ",IND[i][2+6]);
			fprintf(fp,"%d ",IND[i][2+7]);
			fprintf(fp,"%d ",IND[i][2+8]);
			fprintf(fp,"%d ",IND[i][2+9]);
			fprintf(fp,"%d ",IND[i][2+10]);
			fprintf(fp,"%d ",IND[i][2+11]);
			fprintf(fp,"%d ",IND[i][2+16]);
			fprintf(fp,"%d ",IND[i][2+17]);
			fprintf(fp,"%d ",IND[i][2+18]);
			fprintf(fp,"%d ",IND[i][2+19]);
			fprintf(fp,"%d ",IND[i][2+12]);
			fprintf(fp,"%d ",IND[i][2+13]);
			fprintf(fp,"%d ",IND[i][2+14]);
			fprintf(fp,"%d ",IND[i][2+15]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"CELL_TYPES %d\n",NEL);
	for (i=0; i<NEL; i++)
	{
		fprintf(fp,"%d\n",IND[i][0]);
	}
}

void FULLMODEL::ParaView_PrintGrid(FILE *fp)
{	
	ParaView_PrintCRD(fp);
	ParaView_PrintIND(fp);
}

void FULLMODEL::ParaView_PrintMaterial(FILE *fp)
{
	int i;
	fprintf(fp,"SCALARS material_number int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NEL; i++)
	{
		fprintf(fp,"%d\n",MTR[i]);
	}
}

void FULLMODEL::ParaView_StartCellDataSection(FILE *fp)
{
	fprintf(fp,"CELL_DATA %d\n",NEL);
}

void FULLMODEL::ParaView_StartNodeDataSection(FILE *fp)
{
	fprintf(fp,"POINT_DATA %d\n",NN);
}

void FULLMODEL::ParaView_SingleXYZS(FILE *fp, char *name, double *data)
{
	int i,j;
	fprintf(fp,"SCALARS %s_X float 1\n",name);
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e\n",(float)data[i*KORT]);
	}

	fprintf(fp,"SCALARS %s_Y float 1\n",name);
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e\n",(float)data[i*KORT+1]);
	}

	fprintf(fp,"SCALARS %s_Z float 1\n",name);
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e\n",(float)data[i*KORT+2]);
	}

	fprintf(fp,"SCALARS %s_SUM float 1\n",name);
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e\n",(float)(sqrt( data[i*KORT]*data[i*KORT] + data[i*KORT+1]*data[i*KORT+1] + data[i*KORT+2]*data[i*KORT+2] )));
	}
}

void FULLMODEL::ParaView_SingleVector(FILE *fp, char *name, double *data)
{
	int i,j;
	fprintf(fp,"VECTORS %s float\n",name);
//	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e %e %e\n",(float)data[i*KORT],(float)data[i*KORT+1],(float)data[i*KORT+2]);
	}
}

void FULLMODEL::ParaView_SingleVector(FILE *fp, char *name, float *data)
{
	int i,j;
	fprintf(fp,"VECTORS %s float\n",name);
//	fprintf(fp,"LOOKUP_TABLE default\n");
	for (i=0; i<NN; i++)
	{
		fprintf(fp,"%e %e %e\n",data[i*KORT],data[i*KORT+1],data[i*KORT+2]);
	}
}

void FULLMODEL::ParaView_PrintSurfModel(FILE *fp)
{
	int i,j,k,cellsize,surfcelltypes[9];

	for (i=0; i<9; i++) surfcelltypes[i] = 0;
	surfcelltypes[3] = 5; //устанавливает соответствие между количеством узлов ячейки поверхности и ее номером типа в VTK
	surfcelltypes[4] = 9;
	surfcelltypes[6] = 22;
	surfcelltypes[8] = 23;

	ParaView_PrintCRD(fp);

	// расчет количества чисел в поле CELLS
	// (для каждого элемента суммируется его число узлов + 1
	cellsize = 0;
	for (i=0; i < nfaces; i++)
	{
		cellsize += (surfmodel[i].NNface + 1);
	}
	//вывод матрицы индексов
	fprintf(fp,"CELLS %d %d\n",nfaces,cellsize);
	for (i=0; i<nfaces; i++)
	{
		fprintf(fp,"%d ",surfmodel[i].NNface);
		for (j=0; j<surfmodel[i].NNface; j++)
		{
			fprintf(fp,"%d ",surfmodel[i].facenodenums[j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"CELL_TYPES %d\n",nfaces);
	for (i=0; i<nfaces; i++)
	{
		fprintf(fp,"%d\n",surfcelltypes[ surfmodel[i].NNface ]);
	}

}