#include "StdAfx.h"
#include "T_elastat.h"
#include "math.h"
#include "stdio.h"

ELASTAT::ELASTAT()
{
	MM = NULL;
	pfm = NULL;
	psem = NULL;

	nvect = 0;
	LV = NULL;
	RV = NULL;
	fl_LVRV_internal = 0;
}

ELASTAT::~ELASTAT()
{
	FullDel();
}

void ELASTAT::FullDel()
{
	if (fl_LVRV_internal == 1)
	{
		MM->MEM_DEL(LV,nvect,pfm->NNE);
		MM->MEM_DEL(RV,nvect,pfm->NNE);
	}
	nvect = 0;
}

void ELASTAT::AttachToFullModel()
{
	int i,j;

	pfm->nvect = nvect;
	pfm->LV = LV;
	pfm->RV = RV;

}

void ELASTAT::ReinitLVRV(int newnvect)
{
	//выделение пам€ти под векторы нагрузки
	if ( newnvect != nvect )
	{
		if (fl_LVRV_internal == 1)
		{
			LV = MM->MEM_DEL(LV,nvect,pfm->NNE);
			RV = MM->MEM_DEL(RV,nvect,pfm->NNE);
		}
		nvect = newnvect;
		LV = MM->MEM_NEW(LV,nvect,pfm->NNE);
		RV = MM->MEM_NEW(RV,nvect,pfm->NNE);
		fl_LVRV_internal = 1;
	}	
}

void ELASTAT::CreateLV()
{
	int i,j;

	nvect = 1;

	LV = MM->MEM_NEW(LV,nvect,pfm->NNE);
	RV = MM->MEM_NEW(RV,nvect,pfm->NNE);

	for (i=0; i<pfm->NN; i++)
	{
		LV[0][i*pfm->KORT+2] = 1.0; //дл€ теста
	}

	//LV[0][5000*pfm->KORT+0] = 1.0;

}

void ELASTAT::StaticSolve()
{
	AttachToFullModel();
	
	psem->CalcStifMatr();
	//psem->CalcLoadVect(LV, nvect);
	//psem->CalcResultVect(RV, nvect);
	psem->CalcLoadVect();
	psem->CalcResultVect();

//	OutParaviewModel();

}

void ELASTAT::StatMultSolveExt(int nv, double **R)
{
	nvect = nv;
	LV = R;
	RV = LV;
	AttachToFullModel();
	
	psem->CalcStifMatr();
	//psem->CalcLoadVect(LV, nvect);
	//psem->CalcResultVect(RV, nvect);
	psem->CalcLoadVect();
	psem->CalcResultVect();

//	OutParaviewModel();

}

void ELASTAT::OutParaviewModel()
{
	//int i,j,k,KEL_OLD,KORT,KEL,**NND,KND,KN;
	//double *adddisp;
	//char strtmp[256];
	//char strname[256];
	//FILE *fp;

	//KORT = pfm->KORT;
	//KEL = pfm->NEL;
	//NND = pfm->IND;
	//KND = pfm->NN;

	//double maxrezmod = -0.01;
	//for (i=0; i<KND*KORT; i++)
	//{
	//	if ( fabs(RV[0][i]) > maxrezmod ) maxrezmod = fabs(RV[0][i]);
	//}
	//double maxcoord = -100000000.0, mincoord = 100000000.0;
	//for (i=0; i<KND*KORT; i++)
	//{
	//	if ( pfm->CRD[i] > maxcoord ) maxcoord = pfm->CRD[i];
	//	if ( pfm->CRD[i] < mincoord ) mincoord = pfm->CRD[i];
	//}
	//double koefdisp;
	//koefdisp = 0.15*(maxcoord-mincoord)/maxrezmod;

	//sprintf(strname,"%s\\outrez.vtk",pfm->pathmain);

	//fp = fopen(strname,"w");

	//fprintf(fp,"# vtk DataFile Version 3.0\n");
	//fprintf(fp,"%s\n",strname);
	//fprintf(fp,"ASCII\n");
	//fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	//fprintf(fp,"POINTS %d float\n",pfm->NN);
	//for (i=0; i<pfm->NN; i++)
	//{
	//	fprintf(fp,"%f %f %f\n", pfm->CRD[KORT*i] + koefdisp*RV[0][KORT*i], pfm->CRD[KORT*i+1] + koefdisp*RV[0][KORT*i+1], pfm->CRD[KORT*i+2] + koefdisp*RV[0][KORT*i+2]);
	//}

	//KN = 20;

	//if ( KN == 10 )
	//{
	//	fprintf(fp,"CELLS %d %d\n",KEL,11*KEL);
	//	for (i=0; i<KEL; i++)
	//	{
	//		fprintf(fp,"10 \n");
	//		for (j=0; j<KN; j++)
	//		{
	//			fprintf(fp,"%d ",NND[i][j+2]);
	//		}
	//		fprintf(fp,"\n");
	//	}
	//	fprintf(fp,"CELL_TYPES %d\n",KEL);
	//	for (i=0; i<KEL; i++)
	//	{
	//		fprintf(fp,"24\n");
	//	}
	//}
	//if( KN == 20 )
	//{
	//	//KEL_OLD = KEL;

	//	//KEL = 1;
	//	//fprintf(fp,"CELLS %d %d\n",KEL,21*KEL);
	//	//for (i=0; i<KEL; i++)
	//	//{
	//	//	fprintf(fp,"20 \n");
	//	//	for (j=0; j<KN; j++)
	//	//	{
	//	//		fprintf(fp,"%d ",NND[i][j]);
	//	//	}
	//	//	fprintf(fp,"\n");
	//	//}
	//	//fprintf(fp,"CELL_TYPES %d\n",KEL);
	//	//for (i=0; i<KEL; i++)
	//	//{
	//	//	fprintf(fp,"25\n");
	//	//}

	//	//KEL = KEL_OLD;

	//	fprintf(fp,"CELLS %d %d\n",KEL,21*KEL);
	//	for (i=0; i<KEL; i++)
	//	{
	//		fprintf(fp,"20 \n");
	//		/*for (j=0; j<KN; j++)
	//		{
	//			fprintf(fp,"%d ",NND[i][j+2]);
	//		}*/
	//		fprintf(fp,"%d ",NND[i][2+0]);
	//		fprintf(fp,"%d ",NND[i][2+1]);
	//		fprintf(fp,"%d ",NND[i][2+2]);
	//		fprintf(fp,"%d ",NND[i][2+3]);
	//		fprintf(fp,"%d ",NND[i][2+4]);
	//		fprintf(fp,"%d ",NND[i][2+5]);
	//		fprintf(fp,"%d ",NND[i][2+6]);
	//		fprintf(fp,"%d ",NND[i][2+7]);
	//		fprintf(fp,"%d ",NND[i][2+8]);
	//		fprintf(fp,"%d ",NND[i][2+9]);
	//		fprintf(fp,"%d ",NND[i][2+10]);
	//		fprintf(fp,"%d ",NND[i][2+11]);
	//		fprintf(fp,"%d ",NND[i][2+16]);
	//		fprintf(fp,"%d ",NND[i][2+17]);
	//		fprintf(fp,"%d ",NND[i][2+18]);
	//		fprintf(fp,"%d ",NND[i][2+19]);
	//		fprintf(fp,"%d ",NND[i][2+12]);
	//		fprintf(fp,"%d ",NND[i][2+13]);
	//		fprintf(fp,"%d ",NND[i][2+14]);
	//		fprintf(fp,"%d ",NND[i][2+15]);

	//		fprintf(fp,"\n");
	//	}
	//	fprintf(fp,"CELL_TYPES %d\n",KEL);
	//	for (i=0; i<KEL; i++)
	//	{
	//		fprintf(fp,"25\n");
	//	}
	//}

	//fprintf(fp,"CELL_DATA %d\n",KEL);
	//fprintf(fp,"SCALARS material_number int 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KEL; i++)
	//{
	//	fprintf(fp,"%d\n",pfm->MTR[i]);
	//}
	//
	//
	///*fprintf(fp,"CELL_DATA %d\n",KEL);
	//fprintf(fp,"SCALARS SE_number int 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KEL; i++)
	//{
	//	fprintf(fp,"%d\n",psem->Levels[0][i]);
	//}*/

	////fprintf(fp,"POINT_DATA %d\n",KND);
	////fprintf(fp,"VECTORS displacement_vector double\n");
	////for (i=0; i<KND; i++)
	////{
	////	fprintf(fp,"%e %e %e\n",DU[i*KORT],DU[i*KORT+1],DU[i*KORT+2]);
	////}
	//


	//fprintf(fp,"POINT_DATA %d\n",KND);
	//fprintf(fp,"SCALARS displacement_X float 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KND; i++)
	//{
	//	fprintf(fp,"%e\n",(float)(RV[0][i*KORT]*1.e12));
	//}

	////fprintf(fp,"POINT_DATA %d\n",KND);
	//fprintf(fp,"SCALARS displacement_Y float 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KND; i++)
	//{
	//	fprintf(fp,"%e\n",(float)(RV[0][i*KORT+1]*1.e12));
	//}

	////fprintf(fp,"POINT_DATA %d\n",KND);
	//fprintf(fp,"SCALARS displacement_Z float 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KND; i++)
	//{
	//	fprintf(fp,"%e\n",(float)(RV[0][i*KORT+2]*1.e12));
	//}

	////fprintf(fp,"POINT_DATA %d\n",KND);
	//fprintf(fp,"SCALARS displacement_SUM float 1\n");
	//fprintf(fp,"LOOKUP_TABLE default\n");
	//for (i=0; i<KND; i++)
	//{
	//	fprintf(fp,"%e\n",(float)(sqrt(RV[0][i*KORT]*RV[0][i*KORT] + RV[0][i*KORT+1]*RV[0][i*KORT+1] + RV[0][i*KORT+2]*RV[0][i*KORT+2])*1.e12));
	//}


	//fclose(fp);
}


