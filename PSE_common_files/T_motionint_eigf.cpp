#include "StdAfx.h"
#include "T_motionint_eigf.h"
#include "math.h"
#include "stdio.h"

MOTIONEIGF::MOTIONEIGF()
{
	MM = NULL;
	pfm = NULL;

	NNE = 0;
	nform = 0;
	
	EF = NULL;
	FORCE = NULL;
	DISP = NULL;
	FREQ = NULL;
	MODM = NULL;
	MODD = NULL;
	VST = NULL;
	VSTpi = NULL;
	VST0 = NULL;
	qf = NULL;
	qf_ss = NULL;
	qf_es = NULL;
	axqf_in_1 = NULL;

	is_initialized = false;
}

MOTIONEIGF::~MOTIONEIGF()
{
	FullDel();
}

void MOTIONEIGF::AttachToFullModel(FULLMODEL *pf)
{
	pfm = pf;
}

void MOTIONEIGF::FullDel()
{
	EF = MM->MEM_DEL(EF,nform,NNE);
	FORCE = MM->MEM_DEL(FORCE,NNE);
	DISP = MM->MEM_DEL(DISP,NNE);
	FREQ = MM->MEM_DEL(FREQ,nform);
	MODM = MM->MEM_DEL(MODM,nform);
	MODD = MM->MEM_DEL(MODD,nform);
	VST = MM->MEM_DEL(VST,2*nform);
	VSTpi = MM->MEM_NEW(VSTpi,2*nform);
	VST0 = MM->MEM_DEL(VST0,2*nform);
	qf = MM->MEM_DEL(qf,nform);
	qf_ss = MM->MEM_DEL(qf_ss,nform);
	qf_es = MM->MEM_DEL(qf_es,nform);
	if (pfm != NULL) axqf_in_1 = MM->MEM_DEL(axqf_in_1,pfm->KORT,nform);
	nform = 0;
	NNE = 0;
}

void MOTIONEIGF::LoadEIGF()
{
	int i,j,nformloc;
	double tmp,*tmpa;
	FILE *fp;
	char strl[256];

	EF = MM->MEM_DEL(EF,nform,NNE);
	FORCE = MM->MEM_DEL(FORCE,NNE);
	DISP = MM->MEM_DEL(DISP,NNE);
	FREQ = MM->MEM_DEL(FREQ,nform);
	MODM = MM->MEM_DEL(MODM,nform); //удаление прошлых загруженных данных

	sprintf( strl,"%s\\EIGF.bin\0", pfm->pathmain );
	fp = fopen(strl,"rb");
	fread(&NNE,sizeof(int),1,fp);
	if ( NNE != pfm->NNE )
	{
		printf("\nERROR: LoadEIGF  NNE is not equal pfm->NNE; check data-model pair\n");
	}
	fread(&nformloc,sizeof(int),1,fp);
	if ( nform > 0 && nformloc != nform )
	{
		printf("\nERROR: LoadEIGF  new nform is not equal old nform (uncorrect eigf reloading)\n");
	}
	else
	{
		nform = nformloc;
	}

	EF = MM->MEM_NEW(EF,nform,NNE);
	FORCE = MM->MEM_NEW(FORCE,NNE);
	DISP = MM->MEM_NEW(DISP,NNE);
	FREQ = MM->MEM_NEW(FREQ,nform);
	MODM = MM->MEM_NEW(MODM,nform);

	tmpa = new double[NNE];

	for (j=0; j < nform; j++)
	{
		fread(&tmp,sizeof(double),1,fp);
		FREQ[j] = (float)tmp;
		fread(&tmp,sizeof(double),1,fp);
		MODM[j] = (float)tmp;
		fread(tmpa,sizeof(double),NNE,fp);
		for (i=0; i<NNE; i++)
		{
			EF[j][i] = (float)tmpa[i];
		}
	}
	delete []tmpa;

	fclose(fp);
}

void MOTIONEIGF::InitIntegr()
{
	int i;

	VST = MM->MEM_DEL(VST,2*nform);
	VSTpi = MM->MEM_DEL(VSTpi,2*nform);
	VST0 = MM->MEM_DEL(VST0,2*nform);
	qf = MM->MEM_DEL(qf,nform);
	qf_ss = MM->MEM_DEL(qf_ss,nform);
	qf_es = MM->MEM_DEL(qf_es,nform);

	if (nform > 0)
	{
		VST = MM->MEM_NEW(VST,2*nform);
		VSTpi = MM->MEM_NEW(VSTpi,2*nform);
		VST0 = MM->MEM_NEW(VST0,2*nform);
		qf = MM->MEM_NEW(qf,nform);
		qf_ss = MM->MEM_NEW(qf_ss,nform);
		qf_es = MM->MEM_NEW(qf_es,nform);

		for (i=0; i < nform; i++)
		{
			MODD[i] = MODD[i]*FREQ[i];
			FREQ[i] = sqrt( FREQ[i]*FREQ[i] - MODD[i]*MODD[i] );
		}
	}
}

void MOTIONEIGF::InitConstModalDamp(float damp)
{
	int i;

	MODD = MM->MEM_DEL(MODD,nform);
	MODD = MM->MEM_NEW(MODD,nform);
	for (i=0; i < nform; i++)
	{
		MODD[i] = damp;
	}
}

void MOTIONEIGF::InitInertiaUnitModalForces()
{
	int i,j,k;
	axqf_in_1 = MM->MEM_DEL(axqf_in_1,pfm->KORT,nform);
	axqf_in_1 = MM->MEM_NEW(axqf_in_1,pfm->KORT,nform);

	for (i = 0; i < pfm->KORT; i++)
	{
		for (k = 0; k < nform; k++)
		{
			axqf_in_1[i][k] = 0.0;
			for (j = 0; j < pfm->NN; j++)
			{
				axqf_in_1[i][k] -= (float)pfm->Mdiag[j*pfm->KORT + i]*EF[k][j*pfm->KORT + i];
			}
		}
	}
}

void MOTIONEIGF::ResetModalForces_EndStep()
{
	int i;
	for (i=0; i < nform; i++) qf_es[i] = 0.0;
}

void MOTIONEIGF::StartNewStep()
{
	int i;
	it = 0;
	it_err = 1.0;
	for (i=0; i < nform; i++) qf_ss[i] = qf_es[i];
	ResetModalForces_EndStep();
	for (i=0; i < 2*nform; i++) VST0[i] = VST[i];
}

void MOTIONEIGF::StartNewIt()
{
	int i;
	it++;
	for (i=0; i < 2*nform; i++) VSTpi[i] = VST[i];
	ResetModalForces_EndStep();
}

bool MOTIONEIGF::IsItConverged()
{
	int i;
	double lvst, dvst;
	if (it > 1)
	{
		lvst = 0.0;
		dvst = 0.0;
		for (i=0; i<nform; i++)
		{
			lvst += (VST[2*i]*VST[2*i] + VSTpi[2*i]*VSTpi[2*i])/4;
			dvst += (VST[2*i] - VSTpi[2*i])*(VST[2*i] - VSTpi[2*i]);
		}
		if ( lvst > 0.0 )
		{
			it_err = (float)sqrt(dvst/lvst);
		}
		else
		{
			it_err = 0.0;
		}
		if ( it_err < eps_it )
		{
			return(true);
		}
		else
		{
			return(false);
		}
	}
	else
	{
		return(false);
	}
}

void MOTIONEIGF::AddInertiaForces_EndStep(float *ax_val)
{
	int i,k;
	
	for (k=0; k<pfm->KORT; k++)
	{
		for (i=0; i < nform; i++)
		{
			qf_es[i] += axqf_in_1[k][i]*ax_val[k];
		}
	}
}

void MOTIONEIGF::AddFDPointForces_EndStep()
{
	int i,j,k,iform,node;
	double ff;

	for (i=0; i<pfm->npforce; i++)
	{
		if ( pfm->pforce[i].isactive == true )
		{
			for (j=0; j<pfm->pforce[i].ncn; j++)
			{
				node = pfm->pforce[i].ConnectedNodes[j];
				ff = pfm->pforce[i].FF[j];
				for (iform = 0; iform < nform; iform++)
				{
					for (k=0; k< pfm->KORT; k++)
					{
						qf_es[iform] += pfm->pforce[i].Force[k]*ff*EF[iform][node*pfm->KORT + k];
					}
				}
			}
		}
	}

}

void MOTIONEIGF::MakeSingleTimeStep(float dt)
{
	int i;
	for (i = 0; i < nform; i++)
	{ 
		qf[i] = (qf_ss[i] + qf_es[i])/2.0;
	}
	for (i=0; i<nform; i++)
	{
		VST[2*i] = exp(-MODD[i]*dt)*( VST0[2*i]*cos(FREQ[i]*dt) + (VST0[2*i+1] + MODD[i]*VST0[2*i])*sin(FREQ[i]*dt)/FREQ[i] )
			+ (qf[i]/(MODM[i]*FREQ[i]*(MODD[i]*MODD[i] + FREQ[i]*FREQ[i])))*(FREQ[i] - exp(-MODD[i]*dt)*(FREQ[i]*cos(FREQ[i]*dt) + MODD[i]*sin(FREQ[i]*dt)) );

		VST[2*i+1] = exp(-MODD[i]*dt)*( VST0[2*i+1]*( cos(FREQ[i]*dt) - MODD[i]*sin(FREQ[i]*dt)/FREQ[i] ) - VST0[2*i]*( MODD[i]*MODD[i]*sin(FREQ[i]*dt)/FREQ[i] + FREQ[i]*sin(FREQ[i]*dt) ) ) 
			+ qf[i]*sin(FREQ[i]*dt)*exp(-MODD[i]*dt)/(MODM[i]*FREQ[i]);
	}
}

void MOTIONEIGF::EvaluateAllDOF()
{
	int i,j;
	for (i=0; i<NNE; i++)
	{
		DISP[i] = 0.0;
	}
	for (j=0; j<nform; j++)
	{
		for (i=0; i<NNE; i++)
		{
			DISP[i] += VST[2*j]*EF[j][i];
		}
	}
}

void MOTIONEIGF::EvaluateFDPoint_EndStep()
{
	int i,j,k,iform,node;
	double ff;

	for (i=0; i<pfm->npdisp; i++)
	{
		if ( pfm->pdisp[i].isactive == true )
		{
			for (k=0; k<pfm->KORT; k++)
			{
				pfm->pdisp[i].Disp[k] = 0.0;
			}
			for (j=0; j<pfm->pdisp[i].ncn; j++)
			{
				node = pfm->pdisp[i].ConnectedNodes[j];
				ff = pfm->pdisp[i].FF[j];
				for (iform = 0; iform < nform; iform++)
				{
					for (k=0; k< pfm->KORT; k++)
					{
						pfm->pdisp[i].Disp[k] += ff*VST[iform*2]*EF[iform][node*pfm->KORT + k];
					}
				}
			}
		}
	}

}

void MOTIONEIGF::TestCase_InertiaLoading()
{
	//int i,j,k,nstep,printstep;
	//float dt;
	//float t,tend,taxstart;
	//float ax[3],axold[3],axstep[3];
	//char strl[256];
	//FILE *fp;

	////программа реализует внезапное приложение ускорения

	//dt = 0.001; //шаг по времени
	//tend = 0.2; //полное время моделирования
	//taxstart = 0.1*tend;

	//nstep = (int) (tend/dt);
	//printstep = 1;

	//for (k=0;k<3; k++) ax[k] = 0.0;

	//printf("\nMotionInt   TestCase_InertiaLoading()\n");
	//for (i=1; i<=nstep; i++)
	//{
	//	t = i*dt; //время на конец шага
	//	for (k=0;k<3; k++) axold[k] = ax[k];
	//	if ( t > taxstart )
	//	{
	//		ax[0] = 1.0;
	//		ax[1] = 1.0;
	//		ax[2] = 1.0;
	//	}
	//	else
	//	{
	//		ax[0] = 0.0;
	//		ax[1] = 0.0;
	//		ax[2] = 0.0;
	//	}

	//	for (k=0;k<3; k++) axstep[k] = (ax[k] + axold[k])/2;

	//	ResetModalForces();
	//	AddInertiaForces(axstep);

	//	for (k=0; k<2*nform; k++)
	//	{
	//		VST0[k] = VST[k];
	//	}

	//	MakeSingleTimeStep(dt);

	//	if ( i%printstep == 0 )
	//	{
	//		printf("\ri = %d \\ %d",i,nstep);
	//		EvaluateAllDOF();
	//		sprintf( strl,"%s\\MotionInt_DISP_%d.vtk\0", pfm->pathmain, i );
	//		fp = fopen(strl,"w");
	//		pfm->ParaView_PrintGrid(fp);
	//		pfm->ParaView_StartNodeDataSection(fp);
	//		sprintf(strl,"deflection_by_time");
	//		pfm->ParaView_SingleVector(fp,strl,DISP);
	//		fclose(fp);
	//	}
	//}

}
	
	