#include "StdAfx.h"
#include "SE.h"

//extern CLIENTDATA cldat;

void CSE::StiiWrite()
{
	int i,j,itmp;
	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stii);
	fp  = fopen(path,"wb");
	itmp = NIE+1;
	fwrite(&itmp,sizeof(int),1,fp);
	//запись профил€
	fwrite(STII_ENV,sizeof(long long int),NIE+1,fp);
	//запись матрицы
	fwrite(STII,sizeof(double),STII_ENV[NIE],fp);
	fclose(fp);

	//!!!!!!!!!!!!!!!!требуетс€ дополнительный анализ ситуации, когда количество элементов в STII превышает число int. ѕринимает ли fread int64?
}

void CSE::StiiRead()
{
	int i,j,itmp;
	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stii);
	fp  = fopen(path,"rb");
	fread(&itmp,sizeof(int),1,fp);
	STII_ENV = MM->MEM_NEW(STII_ENV,NIE+1);
	//чтение профил€
	fread(STII_ENV,sizeof(long long int),NIE+1,fp);
	STII = MM->MEM_NEW(STII,STII_ENV[NIE]);
	//чтение матрицы
	fread(STII,sizeof(double),STII_ENV[NIE],fp);
	fclose(fp);

	//!!!!!!!!!!!!!!!!требуетс€ дополнительный анализ ситуации, когда количество элементов в STII превышает число int. ѕринимает ли fread int64?
}

void CSE::StiiWriteNZPC()
{
	//чтение II матрицы в разреженном виде с предобусловливателем
	int i,j,itmp;
	char path[512];

	FILE *fp;

	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stii);
	fp  = fopen(path,"wb");
	for  (i=0; i<NIE; i++)
	{
		fwrite(STII_ENVNZ[i],sizeof(int),STII_ENVNZ[i][0]+1,fp);
	}
	for  (i=0; i<NIE; i++)
	{
		fwrite(STII_ENVPC[i],sizeof(int),STII_ENVPC[i][0]+1,fp);
	}
	for  (i=0; i<NIE; i++)
	{
		fwrite(STIINZ[i],sizeof(double),STII_ENVNZ[i][0],fp);
	}
	for  (i=0; i<NIE; i++)
	{
		fwrite(STIIPC[i],sizeof(double),STII_ENVPC[i][0],fp);
	}

	fclose(fp);

}

void CSE::StiiReadNZPC()
{
	//запись II матрицы в разреженном виде с предобусловливателем
	int i,j,itmp;
	char path[512];

	FILE *fp;

	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stii);
	fp  = fopen(path,"rb");
	
	STII_ENVNZ = new int*[NIE];
	STII_ENVPC = new int*[NIE];
	STIINZ = new double*[NIE];
	STIIPC = new double*[NIE];

	for  (i=0; i<NIE; i++)
	{
		fread(&itmp,sizeof(int),1,fp);
		STII_ENVNZ[i] = MM->MEM_NEW(STII_ENVNZ[i],itmp+1);
		STII_ENVNZ[i][0] = itmp;
		fread(&STII_ENVNZ[i][1],sizeof(int),itmp,fp);
	}
	for  (i=0; i<NIE; i++)
	{
		fread(&itmp,sizeof(int),1,fp);
		STII_ENVPC[i] = MM->MEM_NEW(STII_ENVPC[i],itmp+1);
		STII_ENVPC[i][0] = itmp;
		fread(&STII_ENVPC[i][1],sizeof(int),itmp,fp);
	}
	for  (i=0; i<NIE; i++)
	{
		STIINZ[i] =  MM->MEM_NEW(STIINZ[i],STII_ENVNZ[i][0]);
		fread(STIINZ[i],sizeof(double),STII_ENVNZ[i][0],fp);
	}
	for  (i=0; i<NIE; i++)
	{
		STIIPC[i] =  MM->MEM_NEW(STIIPC[i],STII_ENVPC[i][0]);
		fread(STIIPC[i],sizeof(double),STII_ENVPC[i][0],fp);
	}

	fclose(fp);

}

void CSE::WriteNZfull()
{
	//запись полной матрицы —Ё в разреженном виде
	int i,j,itmp;
	char path[512];

	FILE *fp;

	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stnz);
	fp  = fopen(path,"wb");
	for  (i=0; i<NNE; i++)
	{
		fwrite(STII_ENVNZ[i],sizeof(int),STII_ENVNZ[i][0]+1,fp);
	}
	for  (i=0; i<NNE; i++)
	{
		fwrite(STIINZ[i],sizeof(double),STII_ENVNZ[i][0],fp);
	}

	fclose(fp);

}

void CSE::ReadNZfull()
{
	//чтение полной матрицы —Ё в разреженном виде
	int i,j,itmp;
	char path[512];

	FILE *fp;

	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stnz);
	fp  = fopen(path,"rb");
	
	STII_ENVNZ = new int*[NNE];
	STIINZ = new double*[NNE];

	for  (i=0; i<NNE; i++)
	{
		fread(&itmp,sizeof(int),1,fp);
		STII_ENVNZ[i] = NULL;
		STII_ENVNZ[i] = MM->MEM_NEW(STII_ENVNZ[i],itmp+1);
		STII_ENVNZ[i][0] = itmp;
		fread(&STII_ENVNZ[i][1],sizeof(int),itmp,fp);
	}
	for  (i=0; i<NNE; i++)
	{
		STIINZ[i] = NULL;
		STIINZ[i] =  MM->MEM_NEW(STIINZ[i],STII_ENVNZ[i][0]);
		fread(STIINZ[i],sizeof(double),STII_ENVNZ[i][0],fp);
	}

	fclose(fp);

}


void CSE::StisWrite()
{
	int i,j,tmpi,inode;
	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stis);
	fp  = fopen(path,"wb");
	fwrite(&NS,sizeof(int),1,fp);
	//запись профил€
	for (i=0; i<NS; i++)
	{
		tmpi = STIS_ENV[i][1] + 2;
		fwrite(&tmpi,sizeof(int),1,fp);
		fwrite(STIS_ENV[i],sizeof(int),tmpi,fp);
	}
	//запись матрицы
	for (i=0; i<NSE; i++)
	{
		inode = (int)(i/KORT);
		fwrite(STIS[i],sizeof(double),STIS_ENV[inode][0]*KORT,fp);
	}
	fclose(fp);
}

void CSE::StisRead()
{
	int i,j,tmpi,inode;
	char path[512];
	double locmem = 0.0;

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stis);
	fp  = fopen(path,"rb");
	fread(&NS,sizeof(int),1,fp);
	STIS_ENV = new int*[NS];
	//чтение профил€
	for (i=0; i<NS; i++)
	{
		fread(&tmpi,sizeof(int),1,fp);
		STIS_ENV[i] = new int[tmpi];
		fread(STIS_ENV[i],sizeof(int),tmpi,fp);

		locmem += 4*((double)tmpi)/(1024*1024);
	}
	//чтение матрицы
	STIS = new double*[NSE];
	for (i=0; i<NSE; i++)
	{
		inode = (int)(i/KORT);
		tmpi = STIS_ENV[inode][0]*KORT;
		STIS[i] = new double[tmpi];
		fread(STIS[i],sizeof(double),tmpi,fp);
		
		locmem += 8*((double)tmpi)/(1024*1024);
	}
	fclose(fp);

	//!!!!!!! требуетс€ создать критическую секцию дл€ суммировани€ объема пам€ти
	MM->curmem += locmem;
}


void CSE::StssWrite()
{
	int i,j,tmpi,inode;
	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stss);
	fp  = fopen(path,"wb");

//	tmpi = ((int)((NSE*NSE - NSE)/2.0)) + NSE; //количество элементов в полностью заполненном профиле
	fwrite(&SSlength,sizeof(long long int),1,fp);

	//запись матрицы
	fwrite(STSS,sizeof(double),SSlength,fp);

	fclose(fp);
}


void CSE::StssRead()
{
	int i,j,tmpi,inode;
	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\%s",sedat->fullnetfolder,sedat->name_stss);
	fp  = fopen(path,"rb");

	fread(&tmpi,sizeof(long long int),1,fp);  //количество элементов в полностью заполненном профиле

	STSS = MM->MEM_NEW(STSS,tmpi);
	//запись матрицы
	fread(STSS,sizeof(double),tmpi,fp);

	fclose(fp);
}

void CSE::ForseWrite()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\fff%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"wb");
	fwrite(&nvect,sizeof(int),1,fp);
	for(i=0; i<nvect; i++)
	{
		fwrite(FORSE[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
	for(i=0; i<nvect; i++)
	{
		delete []FORSE[i];
	}
	delete []FORSE;
	FORSE = 0;
}

void CSE::ForseWriteSS()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\fff%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"wb");
	fwrite(&nvect,sizeof(int),1,fp);
	fwrite(&NI,sizeof(int),1,fp);
	fwrite(&NS,sizeof(int),1,fp);
	fwrite(&KORT,sizeof(int),1,fp);
	fwrite(REN,sizeof(int),NN,fp);	
	for(i=0; i<nvect; i++)
	{
		fwrite(FORSE[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
	for(i=0; i<nvect; i++)
	{
		delete []FORSE[i];
	}
	delete []FORSE;
	FORSE = 0;
}

void CSE::ForseRead()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\fff%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"rb");
	fread(&nvect,sizeof(int),1,fp);
	FORSE = new double*[nvect];
	for(i=0; i<nvect; i++)
	{
		FORSE[i] = new double [sedat->NNE];
		fread(FORSE[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
}

void CSE::WriteSubspaceKM()
{

	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\sbsp%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"wb");
	fwrite(&nvect,sizeof(int),1,fp);
	
	for(i=0; i<nvect; i++)
	{
		fwrite(KK[i],sizeof(double),nvect,fp);
	}
	fclose(fp);

}

void CSE::DispWrite()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\uuu%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"wb");
	fwrite(&nvect,sizeof(int),1,fp);
	for(i=0; i<nvect; i++)
	{
		fwrite(DISP[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
	for(i=0; i<nvect; i++)
	{
		delete []DISP[i];
	}
	delete []DISP;
	DISP = 0;
}

void CSE::DispWriteSS()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\uuu%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"wb");
	fwrite(&nvect,sizeof(int),1,fp);
	fwrite(&NI,sizeof(int),1,fp);
	fwrite(&NS,sizeof(int),1,fp);
	fwrite(&KORT,sizeof(int),1,fp);
	fwrite(REN,sizeof(int),NN,fp);	
	for(i=0; i<nvect; i++)
	{
		fwrite(DISP[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
	for(i=0; i<nvect; i++)
	{
		delete []DISP[i];
	}
	delete []DISP;
	DISP = 0;
}

void CSE::DispRead()
{
	int i,j;

	char path[512];

	FILE *fp;

	//sprintf(path,"%s\\%s\\%s\\%s",cldat.IPname,cldat.WorkFolder,prob->name,sedat->name_stii);
	sprintf(path,"%s\\uuu%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"rb");
	fread(&nvect,sizeof(int),1,fp);
	DISP = new double*[nvect];
	for(i=0; i<nvect; i++)
	{
		DISP[i] = new double [sedat->NNE];
		fread(DISP[i],sizeof(double),sedat->NNE,fp);
	}
	fclose(fp);
}

void CSE::IntPWrite()
{
	int i,j;
	char path[512];
	FILE *fp;

	if ( sedat->SElevel == 0 )
	{
		sprintf(path,"%s\\intp%s.bin",sedat->fullnetfolder,sedat->namenum);
		fp  = fopen(path,"wb");
		fwrite(&nIntP,sizeof(int),1,fp);
		fwrite(envP,sizeof(int),NEL,fp);
		fwrite(P,sizeof(INTPOINT),nIntP,fp);
		fclose(fp);
	}
}

void CSE::IntPRead()
{
	int i,j;
	char path[512];
	FILE *fp;

	if ( (envP != NULL) || (P != 0) ) DeleteLevel_3();

	sprintf(path,"%s\\intp%s.bin",sedat->fullnetfolder,sedat->namenum);
	fp  = fopen(path,"rb");
	fread(&nIntP,sizeof(int),1,fp);
	envP = MM->MEM_NEW(envP,NEL);
	fread(envP,sizeof(int),NEL,fp);
	P = new INTPOINT[nIntP];
	fread(P,sizeof(INTPOINT),nIntP,fp);
	fclose(fp);
}


void CSE::ReadSEinitials()
{
	int i,j,k;
	int type,nn;
	char str[256],strtmp[256];
	FILE *fp;

	//считывание матрицы индексов
	sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_ind);
	fp = fopen(str,"rb");

	IND = new int*[NEL];
	for (i=0; i<NEL; i++)
	{
		fread(&type,sizeof(int),1,fp);
		fread(&nn,sizeof(int),1,fp);
		IND[i] = new int[nn+2];
		IND[i][0] = type;
		IND[i][1] = nn;
		fread(&IND[i][2],sizeof(int),nn,fp);
	}
	fclose(fp);
	
	if (IND[0][0] == 3 || IND[0][0] == 10 || IND[0][0] == 12 || IND[0][0] == 24 || IND[0][0] == 25 )
	{
		KORT = 3;
	}

	NIE = NI*KORT;
	NSE = NS*KORT;
	NNE = NN*KORT;
	sedat->NNE = NNE;

	if (sedat->SElevel == 0)
	{
		//считывание матрицы координат
		sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_crd);
		fp = fopen(str,"rb");
		CRD = new double [NN*KORT];
		fread(CRD,sizeof(double),NN*KORT,fp);
		fclose(fp);
		//считывание матрицы материалов
		sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_mat);
		fp = fopen(str,"rb");
		MTR = new int [NEL];
		fread(MTR,sizeof(int),NEL,fp);
		fclose(fp);

		//считывание вектора закреплений
		sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_fix);
		fp = fopen(str,"rb");
		FIX = new int [NN*KORT];
		UFIX = new double [NN*KORT];
		fread(FIX,sizeof(int),NN*KORT,fp);
		fread(UFIX,sizeof(double),NN*KORT,fp);
		fclose(fp);
	}

	//считывание массивов перенумерации
	sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_rennodes);
	fp = fopen(str,"rb");
	REN = new int [NN];
	fread(REN,sizeof(int),NN,fp);
	fclose(fp);
	if ( NS != 0) //дл€ верхнего уровн€ нет перенумерации элементов
	{
		sprintf(str,"%s\\%s",sedat->fullnetfolder,sedat->name_renelem);
		fp = fopen(str,"rb");
		REE = new int [NEL];
		fread(REE,sizeof(int),NEL,fp);
		fclose(fp);
	}


}