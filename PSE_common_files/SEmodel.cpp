#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"

char glob_str_path[256];


void startCalcMatr(SEstruct *SE);
void startLoadVect(SEstruct *SE);
void startResVect(SEstruct *SE);
void startSubspaceMatr(SEstruct *SE);

//временный указатель для передачи данных
MATPROP *pmglob;

SEMODEL::SEMODEL(void)
{
	MM = NULL;
	pfm = NULL;
	numSEbylevels = NULL;
	SE = NULL;
	Levels = NULL;
	fp_ase = NULL;

	SEmaxLev = 10;
	loadtype = 1;
	SElevelsnum = 0;
//	sprintf(glob_str_path,"C:\\PSEtemp\\testeigen\\stifmatr");
}

SEMODEL::~SEMODEL(void)
{
	int i,j,k;
	if (SE != NULL)
	{
		for (i=0; i<SEmaxLev; i++)
		{
			for (j=0; j<numSEbylevels[i]; j++)
			{
				for (k=0; k<SE[i][j].NInclSEs; k++)
				{
					if (SE[i][j].PathInclSEs[k] != NULL) { delete [] SE[i][j].PathInclSEs[k]; SE[i][j].PathInclSEs[k] = NULL;}
					if (SE[i][j].namenumInclSEs[k] != NULL) { delete [] SE[i][j].namenumInclSEs[k]; SE[i][j].namenumInclSEs[k] = NULL;}
				}
				if (SE[i][j].PathIncludedInSE != NULL) { delete [] SE[i][j].PathIncludedInSE; SE[i][j].PathIncludedInSE = NULL;}
				if (SE[i][j].namenumIncludedInSE != NULL) { delete [] SE[i][j].namenumIncludedInSE; SE[i][j].namenumIncludedInSE = NULL;}
			}
			if (SE[i] != NULL) {delete []SE[i]; SE[i] = NULL;}
		}
		delete []SE;
	}
	
	Levels = MM->MEM_DEL(Levels,SEmaxLev,1);
	numSEbylevels = MM->MEM_DEL(numSEbylevels,SEmaxLev);
}


void SEMODEL::ReadSEmodel()
{
	sprintf(pathmain,"%s",pfm->pathmain);
	
	
	
	if (loadtype == 1)
	{
		sprintf(pathmatr,"%s",pfm->pathmatr);
		pmglob = pfm->mat;
		ReadFromUzor();
	}
	if (loadtype == 2)
	{
		sprintf(pathmatr,"%s",pfm->pathmatr);
		pmglob = pfm->mat;
		AutoMLDivision_S();
		AutoDivision_PrintLevels_Paraview();
	}
}

void SEMODEL::SE_struct_init(SEstruct *SE, int ilev, int ise)
{
	int i,j,num,nn;
	char namenum[10];

	SE->iSE = ise;
	SE->SElevel = ilev;
	SE->fullnn = pfm->NN;

	namenum[0] = '\0';
	num = ise+1;
	//формирование числовой части имени
	sprintf(namenum,"%s%d",namenum,ilev); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
	nn = (int)(num/1000);
	num -= nn*1000;
	sprintf(namenum,"%s%d",namenum,nn);
	nn = (int)(num/100);
	num -= nn*100;
	sprintf(namenum,"%s%d",namenum,nn);
	nn = (int)(num/10);
	num -= nn*10;
	sprintf(namenum,"%s%d",namenum,nn);
	sprintf(namenum,"%s%d",namenum,num);
	sprintf(SE->namenum,"%s",namenum);

	//запись имен
	sprintf(SE->name_crd,"crd%s.bin",namenum);
	sprintf(SE->name_ind,"ind%s.bin",namenum);
	sprintf(SE->name_renelem,"ree%s.bin",namenum);  
	sprintf(SE->name_rennodes,"ren%s.bin",namenum);

	sprintf(SE->name_mat,"mat%s.bin",namenum);
	sprintf(SE->name_fix,"fix%s.bin",namenum);

	sprintf(SE->name_state,"state%s.dat",namenum);
	sprintf(SE->name_stii,"stii%s.bin",namenum);
	sprintf(SE->name_stnz,"stnz%s.bin",namenum);
	sprintf(SE->name_stis,"stis%s.bin",namenum);
	sprintf(SE->name_stise,"stise%s.bin",namenum);
	sprintf(SE->name_stss,"stss%s.bin",namenum);


	SE->NInclSEs = 0;
	SE->PathIncludedInSE = NULL;
	SE->namenumIncludedInSE = NULL;
	SE->PathInclSEs = NULL;
	SE->namenumInclSEs = NULL;
}

void SEMODEL::ReadFromAutoSE()
{
	int i,j,k,nn,num,type,newtype,k1;
	char str[256];
	char namenum[10],namenum2[10],namenum_level[10];
	char strtmp[256];
	FILE *fp;

	//массив соответствия номеров типов элементов в UZOR и PSE (аналогично PARAVIEW)
	int renumber_types[50];
	int nodenumber_newtypes[50];

	for (i=0; i<50; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}


	renumber_types[0] = 12;
	renumber_types[1] = 25;
	renumber_types[2] = 13;
	renumber_types[3] = -1;
	renumber_types[4] = 10;
	renumber_types[5] = 24;
	renumber_types[6] = -1;
	renumber_types[7] = -1;
	renumber_types[8] = 3;
	renumber_types[9] = 21;
	renumber_types[10] = 3;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;


	/*sprintf(pathmain,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	sprintf(pathmatr,"%s\\ds",pathmain);*/

	sprintf(str,"%s\\semodel.dat",pfm->pathmain);
	fp = fopen(str,"r");
	fscanf(fp,"%s",strtmp);
	SElevelsnum = atoi(strtmp);

	numSEbylevels = new int [SElevelsnum];
	SE = new SEstruct *[SElevelsnum];


	for (i=0; i<(SElevelsnum); i++)
	{
		fscanf(fp,"%s",strtmp);
		numSEbylevels[i] = atoi(strtmp);
		SE[i] = new SEstruct [numSEbylevels[i]];
		//for (j=0; j<numSEbylevels[i]; j++)
		//{
		//	fscanf(fp,"%s",strtmp); //холостое считывание данных о повторяемости (при построении модели она должна быть отключена)
		//}
	}
	fclose(fp);
	////для верхнего уровня файл setpln не содержит информации
	//i = SElevelsnum-1;
	//numSEbylevels[i] = 1;
	//SE[i] = new SEstruct [numSEbylevels[i]];

	//int **RENlev;
	//int *renlevsize;
	//RENlev = new int*[SElevelsnum-1];
	//renlevsize = new int[SElevelsnum-1];
	//for (i=0; i<(SElevelsnum-1); i++)
	//{
	//	RENlev[i] = NULL;
	//	renlevsize[i] = 0;
	//}



	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			namenum[0] = '\0';
			num = j;
			//формирование числовой части имени
			sprintf(namenum,"%s%d",namenum,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			nn = (int)(num/1000);
			num -= nn*1000;
			sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum,"%s%d",namenum,nn);
			sprintf(namenum,"%s%d",namenum,num);
			sprintf(SE[i][j].namenum,"%s",namenum);

			////имена UZOR
			//namenum2[0] = '\0';
			//if (i == (SElevelsnum-1) )
			//{
			//	num = j;
			//}
			//else
			//{
			//	num = j+1;
			//}
			////формирование числовой части имени
			//sprintf(namenum2,"%s%d",namenum2,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			////nn = (int)(num/1000);
			////num -= nn*1000;
			////sprintf(namenum,"%s%d",namenum,nn);
			//nn = (int)(num/100);
			//num -= nn*100;
			//sprintf(namenum2,"%s%d",namenum2,nn);
			//nn = (int)(num/10);
			//num -= nn*10;
			//sprintf(namenum2,"%s%d",namenum2,nn);
			//sprintf(namenum2,"%s%d",namenum2,num);

			//num = 0;
			//namenum_level[0] = '\0';
			////формирование числовой части имени уровня
			//sprintf(namenum_level,"%s%d",namenum_level,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			//
			//nn = (int)(num/100);
			//num -= nn*100;
			//sprintf(namenum_level,"%s%d",namenum_level,nn);
			//nn = (int)(num/10);
			//num -= nn*10;
			//sprintf(namenum_level,"%s%d",namenum_level,nn);
			//sprintf(namenum_level,"%s%d",namenum_level,num);

			
			//считывание характеристик СЭ
			sprintf(str,"%s\\siz%s.dat",pathmatr,namenum);
			fp = fopen(str,"r");
			fscanf(fp,"%s",strtmp);
			SE[i][j].NN = atoi(strtmp);
			fscanf(fp,"%s",strtmp);
			SE[i][j].NEL = atoi(strtmp);
			fscanf(fp,"%s",strtmp);
			SE[i][j].NS = atoi(strtmp);

			SE[i][j].fullnn = pfm->NN;
			
			fclose(fp);

			SE[i][j].NI = SE[i][j].NN - SE[i][j].NS; 
			SE[i][j].SElevel = i;
			SE[i][j].iSE = j;

			//запись имен
			sprintf(SE[i][j].name_crd,"crd%s.bin",namenum);
			sprintf(SE[i][j].name_ind,"ind%s.bin",namenum);
			sprintf(SE[i][j].name_renelem,"ree%s.bin",namenum);  
			sprintf(SE[i][j].name_rennodes,"ren%s.bin",namenum);

			sprintf(SE[i][j].name_mat,"mat%s.bin",namenum);
			sprintf(SE[i][j].name_fix,"fix%s.bin",namenum);

			sprintf(SE[i][j].name_state,"state%s.dat",namenum);
			sprintf(SE[i][j].name_stii,"stii%s.bin",namenum);
			sprintf(SE[i][j].name_stnz,"stnz%s.bin",namenum);
			sprintf(SE[i][j].name_stis,"stis%s.bin",namenum);
			sprintf(SE[i][j].name_stise,"stise%s.bin",namenum);
			sprintf(SE[i][j].name_stss,"stss%s.bin",namenum);

			//преобразование исходных файлов
			if (i == 0)
			{
				sprintf(strtmp,"%s\\crd%s.bin",pathmatr,namenum);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_crd);
				CopyFileA(strtmp,str,false);
			}
			
			sprintf(strtmp,"%s\\ren%s.bin",pathmatr,namenum);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
			CopyFileA(strtmp,str,false);
			fp = fopen(str, "rb");
			int *REN;
			REN = new int[SE[i][j].NN];
			fread(REN,sizeof(int),SE[i][j].NN,fp);
			fclose(fp);


			if ( i != (SElevelsnum - 1) )
			{
				sprintf(strtmp,"%s\\ree%s.bin",pathmatr,namenum);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_renelem);
				CopyFileA(strtmp,str,false);		
			}

			//if ( i == 0 )
			//{
			//	if (j == 0)
			//	{

			//	}
			//}

			
			sprintf(strtmp,"%s\\ind%s.bin",pathmatr,namenum);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_ind);
			CopyFileA(strtmp,str,false);

			
			//создание файлов материалов для СЭ
			if (i == 0)
			{
				sprintf(strtmp,"%s\\mtr%s.bin",pathmatr,namenum);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_mat);
				CopyFileA(strtmp,str,false);
		
				//создание файлов закреплений СЭ
				int *FIX;
				double *UFIX;
				FIX = new int [pfm->KORT*SE[i][j].NN];
				UFIX = new double [pfm->KORT*SE[i][j].NN];
				for (k = 0; k<SE[i][j].NN; k++)
				{
					for (k1=0; k1<pfm->KORT; k1++)
					{
						FIX[ k*pfm->KORT + k1 ] = pfm->FIX[ REN[k]*pfm->KORT + k1 ];
						UFIX[ k*pfm->KORT + k1 ] = pfm->UFIX[ REN[k]*pfm->KORT + k1 ];
					}
				}
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_fix);
				fp = fopen(str,"wb");
				fwrite(FIX,sizeof(int),pfm->KORT*SE[i][j].NN,fp);
				fwrite(UFIX,sizeof(double),pfm->KORT*SE[i][j].NN,fp);
				fclose(fp);
				delete []FIX;
				delete []UFIX;
			}

			SE[i][j].is_matrix_calculated = false;

			delete [] REN;
			//!!!!!!!!!!!!!!!

		}
	}

	//for (i=0; i<(SElevelsnum-1); i++)
	//{
	//	delete [] RENlev[i];
	//}
	//delete []renlevsize;

}

void SEMODEL::ReadFromUzor()
{
	int i,j,k,nn,num,type,newtype,k1;
	char str[256];
	char namenum[10],namenum2[10],namenum_level[10];
	char strtmp[256];
	FILE *fp;

	//массив соответствия номеров типов элементов в UZOR и PSE (аналогично PARAVIEW)
	int renumber_types[50];
	int nodenumber_newtypes[50];

	for (i=0; i<50; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}


	renumber_types[0] = 12;
	renumber_types[1] = 25;
	renumber_types[2] = 13;
	renumber_types[3] = -1;
	renumber_types[4] = 10;
	renumber_types[5] = 24;
	renumber_types[6] = -1;
	renumber_types[7] = -1;
	renumber_types[8] = 3;
	renumber_types[9] = 21;
	renumber_types[10] = 3;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;


	/*sprintf(pathmain,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	sprintf(pathmatr,"%s\\ds",pathmain);*/

	sprintf(str,"%s\\setpln",pathmain);
	fp = fopen(str,"r");
	fscanf(fp,"%s",strtmp);
	SElevelsnum = atoi(strtmp)+1;

	numSEbylevels = new int [SElevelsnum];
	SE = new SEstruct *[SElevelsnum];


	for (i=0; i<(SElevelsnum-1); i++)
	{
		fscanf(fp,"%s",strtmp);
		numSEbylevels[i] = atoi(strtmp);
		SE[i] = new SEstruct [numSEbylevels[i]];
		for (j=0; j<numSEbylevels[i]; j++)
		{
			fscanf(fp,"%s",strtmp); //холостое считывание данных о повторяемости (при построении модели она должна быть отключена)
		}
	}
	fclose(fp);
	//для верхнего уровня файл setpln не содержит информации
	i = SElevelsnum-1;
	numSEbylevels[i] = 1;
	SE[i] = new SEstruct [numSEbylevels[i]];

	int **RENlev;
	int *renlevsize;
	RENlev = new int*[SElevelsnum-1];
	renlevsize = new int[SElevelsnum-1];
	for (i=0; i<(SElevelsnum-1); i++)
	{
		RENlev[i] = NULL;
		renlevsize[i] = 0;
	}



	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			namenum[0] = '\0';
			num = j+1;
			//формирование числовой части имени
			sprintf(namenum,"%s%d",namenum,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			nn = (int)(num/1000);
			num -= nn*1000;
			sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum,"%s%d",namenum,nn);
			sprintf(namenum,"%s%d",namenum,num);
			sprintf(SE[i][j].namenum,"%s",namenum);

			//имена UZOR
			namenum2[0] = '\0';
			if (i == (SElevelsnum-1) )
			{
				num = j;
			}
			else
			{
				num = j+1;
			}
			//формирование числовой части имени
			sprintf(namenum2,"%s%d",namenum2,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			//nn = (int)(num/1000);
			//num -= nn*1000;
			//sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum2,"%s%d",namenum2,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum2,"%s%d",namenum2,nn);
			sprintf(namenum2,"%s%d",namenum2,num);

			num = 0;
			namenum_level[0] = '\0';
			//формирование числовой части имени уровня
			sprintf(namenum_level,"%s%d",namenum_level,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum_level,"%s%d",namenum_level,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum_level,"%s%d",namenum_level,nn);
			sprintf(namenum_level,"%s%d",namenum_level,num);

			
			//считывание характеристик СЭ
			sprintf(str,"%s\\siz%s.dat",pathmatr,namenum2);
			fp = fopen(str,"r");
			fscanf(fp,"%s",strtmp);
			SE[i][j].NN = atoi(strtmp);
			fscanf(fp,"%s",strtmp);
			SE[i][j].NEL = atoi(strtmp);

			SE[i][j].fullnn = pfm->NN;
			
			if ( i != (SElevelsnum - 1) )
			{
				fscanf(fp,"%s",strtmp);
				SE[i][j].NS = atoi(strtmp);
			}
			else
			{
				SE[i][j].NS = 0;
			}
			fclose(fp);

			SE[i][j].NI = SE[i][j].NN - SE[i][j].NS; 
			SE[i][j].SElevel = i;
			SE[i][j].iSE = j;

			//запись имен
			sprintf(SE[i][j].name_crd,"crd%s.bin",namenum);
			sprintf(SE[i][j].name_ind,"ind%s.bin",namenum);
			sprintf(SE[i][j].name_renelem,"ree%s.bin",namenum);  
			sprintf(SE[i][j].name_rennodes,"ren%s.bin",namenum);

			sprintf(SE[i][j].name_mat,"mat%s.bin",namenum);
			sprintf(SE[i][j].name_fix,"fix%s.bin",namenum);

			sprintf(SE[i][j].name_state,"state%s.dat",namenum);
			sprintf(SE[i][j].name_stii,"stii%s.bin",namenum);
			sprintf(SE[i][j].name_stnz,"stnz%s.bin",namenum);
			sprintf(SE[i][j].name_stis,"stis%s.bin",namenum);
			sprintf(SE[i][j].name_stise,"stise%s.bin",namenum);
			sprintf(SE[i][j].name_stss,"stss%s.bin",namenum);

			//преобразование исходных файлов
			if (i == 0)
			{
				sprintf(strtmp,"%s\\crd%s.dat",pathmatr,namenum2);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_crd);
				CopyFileA(strtmp,str,false);
			}
			

			
			if ( (i > 0) ) //для нулевого уровня массив перенумерации уровня в целом не требуется
			{
				//считывание матрицы перенумерации уровня
				sprintf(strtmp,"%s\\siz%s.dat",pathmatr,namenum_level);
				fp = fopen(strtmp,"r");
				fscanf(fp,"%s",strtmp);
				renlevsize[i-1] = atoi(strtmp);
				fclose(fp);
				sprintf(strtmp,"%s\\msv%s.dat",pathmatr,namenum_level);
				RENlev[i-1] = new int[renlevsize[i-1]];
				fp = fopen(strtmp,"rb");
				fread(RENlev[i-1],sizeof(int),renlevsize[i-1],fp);
				if ( i > 1 )
				{
					for (k=0; k<renlevsize[i-1]; k++) 
					{
						RENlev[i-1][k] = RENlev[i-2][ RENlev[i-1][k] - 1 ];
					}
				}
				else
				{
					for (k=0; k<renlevsize[i-1]; k++) 
					{
						RENlev[i-1][k] = RENlev[i-1][k] - 1;
					}
				}
				fclose(fp);
			}

			//разобраться с перенумерацией узлов на старших уровнях

			sprintf(strtmp,"%s\\msv%s.dat",pathmatr,namenum2);
			int *REN;
			REN = new int[SE[i][j].NN];
			fp = fopen(strtmp,"rb");
			fread(REN,sizeof(int),SE[i][j].NN,fp);
			if ( i == (SElevelsnum - 1) )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] = RENlev[i-1][ k ];
			}
			if ( (i > 0) && (i != (SElevelsnum - 1)) )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] = RENlev[i-1][ REN[k] - 1 ];
			}
			if ( i == 0 )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] =  REN[k] - 1;
			}
			fclose(fp);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
			fp = fopen(str,"wb");
			fwrite(REN,sizeof(int),SE[i][j].NN,fp);
			fclose(fp);
			
			int *REE;
			if ( i != (SElevelsnum - 1) )
			{
				sprintf(str,"%s\\mne%s.dat",pathmatr,namenum2);
				
				REE = new int[SE[i][j].NEL];
				fp = fopen(str,"r");
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					REE[k] = atoi(strtmp) - 1;
				}
				fclose (fp);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_renelem);
				fp = fopen(str,"wb");
				fwrite(REE,sizeof(int),SE[i][j].NEL,fp);
				fclose(fp);
			}
			

			sprintf(str,"%s\\ind%s.dat",pathmatr,namenum2);
			int **IND;
			IND = new int *[SE[i][j].NEL];
			
			fp = fopen(str,"r");
			if (i == 0)
			{
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					type = atoi(strtmp);
					newtype = renumber_types[type];
					IND[k] = new int[nodenumber_newtypes[newtype]+2];
					IND[k][0] = newtype;
					IND[k][1] = nodenumber_newtypes[newtype];

					for (k1=0; k1<IND[k][1]; k1++)
					{
						fscanf(fp,"%s",strtmp);
						IND[k][2+k1] = atoi(strtmp) - 1;
					}
				}
			}
			else
			{
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					nn = atoi(strtmp);
					IND[k] = new int[nn+2];
					IND[k][0] = 100; //тип указывающий на СЭ
					IND[k][1] = nn;

					for (k1=0; k1< IND[k][1]; k1++)
					{
						fscanf(fp,"%s",strtmp);
						IND[k][2+k1] = atoi(strtmp) - 1;
					}
				}
			}
			fclose (fp);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_ind);
			fp = fopen(str,"wb");
			for (k=0; k<SE[i][j].NEL; k++)
			{
				fwrite(IND[k],sizeof(int),IND[k][1]+2,fp);
			}
			fclose(fp);

			for (k=0; k<SE[i][j].NEL; k++)
			{
				delete [] IND[k];
			}
			delete [] IND;

			
			//создание файлов материалов для СЭ
			if (i == 0)
			{
				int *MTRSE;
				MTRSE = new int[SE[i][j].NEL];
				for (k=0; k<SE[i][j].NEL; k++)
				{
					MTRSE[k] = pfm->MTR[REE[k]];
				}
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_mat);
				fp = fopen(str,"wb");
				fwrite(MTRSE,sizeof(int),SE[i][j].NEL,fp);
				fclose(fp);
				delete []MTRSE;
			

				//создание файлов закреплений СЭ
				int *FIX;
				double *UFIX;
				FIX = new int [pfm->KORT*SE[i][j].NN];
				UFIX = new double [pfm->KORT*SE[i][j].NN];
				for (k = 0; k<SE[i][j].NN; k++)
				{
					for (k1=0; k1<pfm->KORT; k1++)
					{
						FIX[ k*pfm->KORT + k1 ] = pfm->FIX[ REN[k]*pfm->KORT + k1 ];
						UFIX[ k*pfm->KORT + k1 ] = pfm->UFIX[ REN[k]*pfm->KORT + k1 ];
					}
				}
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_fix);
				fp = fopen(str,"wb");
				fwrite(FIX,sizeof(int),pfm->KORT*SE[i][j].NN,fp);
				fwrite(UFIX,sizeof(double),pfm->KORT*SE[i][j].NN,fp);
				fclose(fp);
				delete []FIX;
				delete []UFIX;
			}

			SE[i][j].is_matrix_calculated = false;

			//!!!!!!!!!!!!!!!

			if ( i != (SElevelsnum - 1) ) delete [] REE;
			delete [] REN;
		}
	}

	for (i=0; i<(SElevelsnum-1); i++)
	{
		delete [] RENlev[i];
	}
	delete []renlevsize;

}

void SEMODEL::ReadFromUzor_NewFormat()
{
	int i,j,k,nn,num,type,newtype,k1;
	char str[256];
	char namenum[10],namenum2[10],namenum_level[10];
	char strtmp[256];
	FILE *fp;

	//массив соответствия номеров типов элементов в UZOR и PSE (аналогично PARAVIEW)
	int renumber_types[50];
	int nodenumber_newtypes[50];

	for (i=0; i<50; i++)
	{
		renumber_types[i] = -1;
		nodenumber_newtypes[i] = -1;
	}


	renumber_types[0] = 12;
	renumber_types[1] = 25;
	renumber_types[2] = 13;
	renumber_types[3] = -1;
	renumber_types[4] = 10;
	renumber_types[5] = 24;
	renumber_types[6] = -1;
	renumber_types[7] = -1;
	renumber_types[8] = 3;
	renumber_types[9] = 21;
	renumber_types[10] = 3;

	nodenumber_newtypes[3] = 2;
	nodenumber_newtypes[10] = 4;
	nodenumber_newtypes[12] = 8;
	nodenumber_newtypes[13] = 6;
	nodenumber_newtypes[21] = 3;
	nodenumber_newtypes[24] = 10;
	nodenumber_newtypes[25] = 20;


	/*sprintf(pathmain,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
	sprintf(pathmatr,"%s\\ds",pathmain);*/

	sprintf(str,"%s\\setpln",pathmain);
	fp = fopen(str,"r");
	fscanf(fp,"%s",strtmp);
	SElevelsnum = atoi(strtmp)+1;

	numSEbylevels = new int [SElevelsnum];
	SE = new SEstruct *[SElevelsnum];


	for (i=0; i<(SElevelsnum-1); i++)
	{
		fscanf(fp,"%s",strtmp);
		numSEbylevels[i] = atoi(strtmp);
		SE[i] = new SEstruct [numSEbylevels[i]];
		for (j=0; j<numSEbylevels[i]; j++)
		{
			fscanf(fp,"%s",strtmp); //холостое считывание данных о повторяемости (при построении модели она должна быть отключена)
		}
	}
	fclose(fp);
	//для верхнего уровня файл setpln не содержит информации
	i = SElevelsnum-1;
	numSEbylevels[i] = 1;
	SE[i] = new SEstruct [numSEbylevels[i]];

	int **RENlev;
	int *renlevsize;
	RENlev = new int*[SElevelsnum-1];
	renlevsize = new int[SElevelsnum-1];
	for (i=0; i<(SElevelsnum-1); i++)
	{
		RENlev[i] = NULL;
		renlevsize[i] = 0;
	}



	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			namenum[0] = '\0';
			SE_struct_init(&SE[i][j], i, j);
			sprintf(namenum,"%s",SE[i][j].namenum);

			//имена UZOR
			namenum2[0] = '\0';
			if (i == (SElevelsnum-1) )
			{
				num = j;
			}
			else
			{
				num = j+1;
			}
			//формирование числовой части имени
			sprintf(namenum2,"%s%d",namenum2,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			//nn = (int)(num/1000);
			//num -= nn*1000;
			//sprintf(namenum,"%s%d",namenum,nn);
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum2,"%s%d",namenum2,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum2,"%s%d",namenum2,nn);
			sprintf(namenum2,"%s%d",namenum2,num);

			num = 0;
			namenum_level[0] = '\0';
			//формирование числовой части имени уровня
			sprintf(namenum_level,"%s%d",namenum_level,i); //номер уровня (предполагается, что уровней не больше 10 (от нуля)
			
			nn = (int)(num/100);
			num -= nn*100;
			sprintf(namenum_level,"%s%d",namenum_level,nn);
			nn = (int)(num/10);
			num -= nn*10;
			sprintf(namenum_level,"%s%d",namenum_level,nn);
			sprintf(namenum_level,"%s%d",namenum_level,num);

			
			//считывание характеристик СЭ
			sprintf(str,"%s\\siz%s.dat",pathmatr,namenum2);
			fp = fopen(str,"r");
			fscanf(fp,"%s",strtmp);
			SE[i][j].NN = atoi(strtmp);
			fscanf(fp,"%s",strtmp);
			SE[i][j].NEL = atoi(strtmp);

			SE[i][j].fullnn = pfm->NN;
			
			if ( i != (SElevelsnum - 1) )
			{
				fscanf(fp,"%s",strtmp);
				SE[i][j].NS = atoi(strtmp);
			}
			else
			{
				SE[i][j].NS = 0;
			}
			fclose(fp);

			SE[i][j].NI = SE[i][j].NN - SE[i][j].NS; 
			SE[i][j].SElevel = i;
			SE[i][j].iSE = j;

			//преобразование исходных файлов
			if (i == 0)
			{
				sprintf(strtmp,"%s\\crd%s.dat",pathmatr,namenum2);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_crd);
				CopyFileA(strtmp,str,false);
			}
			
			if ( (i > 0) ) //для нулевого уровня массив перенумерации уровня в целом не требуется
			{
				//считывание матрицы перенумерации уровня
				sprintf(strtmp,"%s\\siz%s.dat",pathmatr,namenum_level);
				fp = fopen(strtmp,"r");
				fscanf(fp,"%s",strtmp);
				renlevsize[i-1] = atoi(strtmp);
				fclose(fp);
				sprintf(strtmp,"%s\\msv%s.dat",pathmatr,namenum_level);
				RENlev[i-1] = new int[renlevsize[i-1]];
				fp = fopen(strtmp,"rb");
				fread(RENlev[i-1],sizeof(int),renlevsize[i-1],fp);
				if ( i > 1 )
				{
					for (k=0; k<renlevsize[i-1]; k++) 
					{
						RENlev[i-1][k] = RENlev[i-2][ RENlev[i-1][k] - 1 ];
					}
				}
				else
				{
					for (k=0; k<renlevsize[i-1]; k++) 
					{
						RENlev[i-1][k] = RENlev[i-1][k] - 1;
					}
				}
				fclose(fp);
			}

			//разобраться с перенумерацией узлов на старших уровнях

			sprintf(strtmp,"%s\\msv%s.dat",pathmatr,namenum2);
			int *REN;
			REN = new int[SE[i][j].NN];
			fp = fopen(strtmp,"rb");
			fread(REN,sizeof(int),SE[i][j].NN,fp);
			if ( i == (SElevelsnum - 1) )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] = RENlev[i-1][ k ];
			}
			if ( (i > 0) && (i != (SElevelsnum - 1)) )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] = RENlev[i-1][ REN[k] - 1 ];
			}
			if ( i == 0 )
			{
				for (k=0; k<SE[i][j].NN; k++) REN[k] =  REN[k] - 1;
			}
			fclose(fp);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
			fp = fopen(str,"wb");
			fwrite(REN,sizeof(int),SE[i][j].NN,fp);
			fclose(fp);
			
			int *REE;
			if ( i != (SElevelsnum - 1) )
			{
				sprintf(str,"%s\\mne%s.dat",pathmatr,namenum2);
				
				REE = new int[SE[i][j].NEL];
				fp = fopen(str,"r");
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					REE[k] = atoi(strtmp) - 1;
				}
				fclose (fp);
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_renelem);
				fp = fopen(str,"wb");
				fwrite(REE,sizeof(int),SE[i][j].NEL,fp);
				fclose(fp);
			}
			

			sprintf(str,"%s\\ind%s.dat",pathmatr,namenum2);
			int **IND;
			IND = new int *[SE[i][j].NEL];
			
			fp = fopen(str,"r");
			if (i == 0)
			{
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					type = atoi(strtmp);
					newtype = renumber_types[type];
					IND[k] = new int[nodenumber_newtypes[newtype]+2];
					IND[k][0] = newtype;
					IND[k][1] = nodenumber_newtypes[newtype];

					for (k1=0; k1<IND[k][1]; k1++)
					{
						fscanf(fp,"%s",strtmp);
						IND[k][2+k1] = atoi(strtmp) - 1;
					}
				}
			}
			else
			{
				for (k=0; k<SE[i][j].NEL; k++)
				{
					fscanf(fp,"%s",strtmp);
					nn = atoi(strtmp);
					IND[k] = new int[nn+2];
					IND[k][0] = 100; //тип указывающий на СЭ
					IND[k][1] = nn;

					for (k1=0; k1< IND[k][1]; k1++)
					{
						fscanf(fp,"%s",strtmp);
						IND[k][2+k1] = atoi(strtmp) - 1;
					}
				}
			}
			fclose (fp);
			sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_ind);
			fp = fopen(str,"wb");
			for (k=0; k<SE[i][j].NEL; k++)
			{
				fwrite(IND[k],sizeof(int),IND[k][1]+2,fp);
			}
			fclose(fp);

			for (k=0; k<SE[i][j].NEL; k++)
			{
				delete [] IND[k];
			}
			delete [] IND;

			
			//создание файлов материалов для СЭ
			if (i == 0)
			{
				int *MTRSE;
				MTRSE = new int[SE[i][j].NEL];
				for (k=0; k<SE[i][j].NEL; k++)
				{
					MTRSE[k] = pfm->MTR[REE[k]];
				}
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_mat);
				fp = fopen(str,"wb");
				fwrite(MTRSE,sizeof(int),SE[i][j].NEL,fp);
				fclose(fp);
				delete []MTRSE;
			

				//создание файлов закреплений СЭ
				int *FIX;
				double *UFIX;
				FIX = new int [pfm->KORT*SE[i][j].NN];
				UFIX = new double [pfm->KORT*SE[i][j].NN];
				for (k = 0; k<SE[i][j].NN; k++)
				{
					for (k1=0; k1<pfm->KORT; k1++)
					{
						FIX[ k*pfm->KORT + k1 ] = pfm->FIX[ REN[k]*pfm->KORT + k1 ];
						UFIX[ k*pfm->KORT + k1 ] = pfm->UFIX[ REN[k]*pfm->KORT + k1 ];
					}
				}
				sprintf(str,"%s\\%s",pathmatr,SE[i][j].name_fix);
				fp = fopen(str,"wb");
				fwrite(FIX,sizeof(int),pfm->KORT*SE[i][j].NN,fp);
				fwrite(UFIX,sizeof(double),pfm->KORT*SE[i][j].NN,fp);
				fclose(fp);
				delete []FIX;
				delete []UFIX;
			}

			SE[i][j].is_matrix_calculated = false;

			//!!!!!!!!!!!!!!!

			if ( i != (SElevelsnum - 1) ) delete [] REE;
			delete [] REN;
		}
	}

	for (i=0; i<(SElevelsnum-1); i++)
	{
		delete [] RENlev[i];
	}
	delete []renlevsize;

}

void SEMODEL::CalcStifMatr()
{
	int i,j,k,nse;
	char strcp[256],strfrom[256],strto[256],strtmp[256];
	FILE *fp;

	//sprintf(strcp,"C:\\PSEtemp\\testeigen\\stifmatr");
	sprintf(strcp,"%s",glob_str_path);


	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			if (SE[i][j].is_matrix_calculated == false)
			{

				//заполнение места расположения СЭ
				sprintf(SE[i][j].fullnetfolder,"%s",strcp);
				//1 - копирование исходных файлов СЭ

				
				sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_ind);
				sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_ind);
				CopyFileA(strfrom,strto,false);
				if ( i != (SElevelsnum - 1) )
				{
					sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_renelem);
					sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_renelem);
					CopyFileA(strfrom,strto,false);
				}
				sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
				sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_rennodes);
				CopyFileA(strfrom,strto,false);
				
				
				if (i == 0)
				{
					sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_crd);
					sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_crd);
					CopyFileA(strfrom,strto,false);
					sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_fix);
					sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_fix);
					CopyFileA(strfrom,strto,false);
					sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_mat);
					sprintf(strto,"%s\\%s",SE[i][j].fullnetfolder,SE[i][j].name_mat);
					CopyFileA(strfrom,strto,false);
				}

				//заполнение расположения составляющих СЭ (для СЭ > 0 уровня)
				if (i > 0)
				{
					SE[i][j].NInclSEs = SE[i][j].NEL;
					SE[i][j].PathInclSEs = new char *[SE[i][j].NInclSEs];
					SE[i][j].namenumInclSEs = new char *[SE[i][j].NInclSEs];
					int *REE;
					REE = new int[SE[i][j].NEL];

					if ( i != (SElevelsnum-1) )
					{
						sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_renelem);
						fp = fopen(strfrom,"rb");
						fread(REE,sizeof(int),SE[i][j].NEL,fp);
						fclose(fp);
					}
					else
					{
						for (k=0; k<SE[i][j].NEL; k++)
						{
							REE[k] = k;
						}
					}
				
					for (k=0; k<SE[i][j].NEL; k++)
					{
						SE[i][j].PathInclSEs[k] = new char [256];
						SE[i][j].namenumInclSEs[k] = new char [10];
						sprintf(SE[i][j].PathInclSEs[k],"%s",SE[i-1][REE[k]].fullnetfolder);
						sprintf(SE[i][j].namenumInclSEs[k],"%s",SE[i-1][REE[k]].namenum);

						SE[i-1][REE[k]].PathIncludedInSE = new char[256];
						SE[i-1][REE[k]].namenumIncludedInSE = new char[10];
						sprintf(SE[i-1][REE[k]].PathIncludedInSE,"%s",SE[i][j].fullnetfolder);
						sprintf(SE[i-1][REE[k]].namenumIncludedInSE,"%s",SE[i][j].namenum);					
					}

					delete []REE;
				}

				//SE[i][j].fl_PCG_internal = true;
				SE[i][j].fl_PCG_internal = false;

				//передача информации о материалах модели (однократно для каждого клиента)

				//!!!!!!!!!!!!!!!

				//передача структуры СЭ


				// запуск расчета
				startCalcMatr(&SE[i][j]);
			}

		}
	}

}

void startCalcMatr(SEstruct *SE)
{
	int i,j;

	CSE *pSE;

	pSE = new CSE[1];

	pSE->pmat = pmglob;

	pSE->Initialize(SE);
	pSE->MainMatrCalc();
	pSE->FULLDEL();

	delete []pSE;

}

void SEMODEL::CalcLoadVect()
{
	int i,j,k,nse,ii,jj,iur;
	char strcp[256],strfrom[256],strto[256],strtmp[256];
	FILE *fp;

	//sprintf(strcp,"C:\\PSEtemp\\testeigen\\stifmatr");
	sprintf(strcp,"%s",glob_str_path);

	int *FLAG; //для учета приложения силы в граничном узле
	FLAG = new int[pfm->NNE];
	for (i=0; i<pfm->NNE; i++)
	{
		FLAG[i] = 0;
	}


	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			if (i == 0)
			{
				//1 - загрузка вектора перенумерации узлов
				sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
				fp = fopen(strfrom,"rb");
				int *REN;
				REN = new int [SE[i][j].NN];
				fread(REN,sizeof(int),SE[i][j].NN,fp);
				fclose(fp);

				//2 - создание векторов нагрузки
				int nvect;
				double **LVSE;
				nvect = pfm->nvect;
				LVSE = new double*[nvect];
				for (k=0; k<nvect; k++)
				{
					LVSE[k] = new double [SE[i][j].NN*pfm->KORT];
				}

				for (ii=0; ii<SE[i][j].NN; ii++)
				{
					for (jj=0; jj<pfm->KORT; jj++)
					{
						iur = REN[ii]*pfm->KORT + jj;
						if (FLAG[iur] == 0)
						{
							for (k=0; k<nvect; k++)
							{
								LVSE[k][ii*pfm->KORT+jj] = pfm->LV[k][ iur ];				
							}
							FLAG[iur] = 1;
						}
						else
						{
							for (k=0; k<nvect; k++)
							{
								LVSE[k][ii*pfm->KORT+jj] = 0.0;				
							}
						}
					}
				}

				sprintf(strto,"%s\\fff%s.bin",SE[i][j].fullnetfolder,SE[i][j].namenum);
				fp = fopen(strto,"wb");
				fwrite(&nvect,sizeof(int),1,fp);
				for (k=0; k<nvect; k++)
				{
					fwrite(LVSE[k],sizeof(double),SE[i][j].NN*pfm->KORT,fp);
				}
				fclose(fp);

				delete []REN;
				for (k=0; k<nvect; k++)
				{
					delete []LVSE[k];
				}
				delete []LVSE;

			}

			

			//передача структуры СЭ


			// запуск расчета
			startLoadVect(&SE[i][j]);

		}
	}

	delete []FLAG;

}

void startLoadVect(SEstruct *SE)
{
	int i,j;

	CSE *pSE;

	pSE = new CSE[1];

	pSE->pmat = pmglob;

	pSE->Initialize(SE);
	pSE->MainLoadCalc();
	pSE->FULLDEL();

	delete []pSE;

}


void SEMODEL::CalcResultVect()
{
	int i,j,k,nse,tmpi,ii,jj;
	char strcp[256],strfrom[256],strto[256],strtmp[256];
	FILE *fp;

	//sprintf(strcp,"C:\\PSEtemp\\testeigen\\stifmatr");
	sprintf(strcp,"%s",glob_str_path);


	for (i=SElevelsnum-2; i>=0; i--)
	{//цикл по уровням сверху-вниз
		//для верхнего уровня обработка (расчет вектора решения) была произведена при формировании векторов нагрузок
		for (j=0; j<numSEbylevels[i]; j++)
		{
			//передача структуры СЭ

			// запуск расчета
			// клиенты самостоятельно считывают вектора перемещений внутренних узлов соответствующих СЭ верхнего уровня
			
			startResVect(&SE[i][j]);

			if (i == 0)
			{ //для СЭ 0 уровня необходимо "собрать" результат для всех узлов полной модели
				

				//1 - загрузка вектора перенумерации узлов
				sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
				fp = fopen(strfrom,"rb");
				int *REN;
				REN = new int [SE[i][j].NN];
				fread(REN,sizeof(int),SE[i][j].NN,fp);
				fclose(fp);

				
				//2 - создание векторов решения
				int nvect;
				double **RVSE;
				nvect = pfm->nvect;
				RVSE = new double*[nvect];
				for (k=0; k<nvect; k++)
				{
					RVSE[k] = new double [SE[i][j].NN*pfm->KORT];
				}
				//загрузка векторов решения
				sprintf(strfrom,"%s\\uuu%s.bin",SE[i][j].fullnetfolder,SE[i][j].namenum);
				fp = fopen(strfrom,"rb");
				fread(&tmpi,sizeof(int),1,fp);
				for (k=0; k<nvect; k++)
				{
					fread(RVSE[k],sizeof(double),SE[i][j].NN*pfm->KORT,fp);
				}
				fclose(fp);


				for (ii=0; ii<SE[i][j].NN; ii++)
				{
					for (jj=0; jj<pfm->KORT; jj++)
					{
						for (k=0; k<nvect; k++)
						{
							pfm->RV[k][ REN[ii]*pfm->KORT + jj ] = RVSE[k][ii*pfm->KORT+jj]; 
						}
					}
				}

				delete []REN;
				for (k=0; k<nvect; k++)
				{
					delete []RVSE[k];
				}
				delete []RVSE;
			}
	
		}
	}

}

void startResVect(SEstruct *SE)
{
	int i,j;

	CSE *pSE;

	pSE = new CSE[1];

	pSE->pmat = pmglob;

	pSE->Initialize(SE);
	pSE->MainResCalc();
	pSE->FULLDEL();

	delete []pSE;

}


void SEMODEL::CalcSubspaceMatrix(double **KK)
{
	int i,j,k,nse,ii,jj,iur,itmp;
	char strcp[256],strfrom[256],strto[256],strtmp[256];
	double **KK_1;
	FILE *fp;

	//sprintf(strcp,"C:\\PSEtemp\\testeigen\\stifmatr");
	sprintf(strcp,"%s",glob_str_path);

	KK_1 = new double *[pfm->nvect];
	for (i=0; i<pfm->nvect; i++)
	{
		KK_1[i] = new double [pfm->nvect];
	}
	for (i=0; i<pfm->nvect; i++)
	{
		for (j=0; j<pfm->nvect; j++)
		{
			KK[i][j] = 0.0;
		}
	}

	i = 0; //обработка ведется только для нулевого уровня СЭ
	for (j=0; j<numSEbylevels[i]; j++)
	{

		//1 - загрузка вектора перенумерации узлов
		sprintf(strfrom,"%s\\%s",pathmatr,SE[i][j].name_rennodes);
		fp = fopen(strfrom,"rb");
		int *REN;
		REN = new int [SE[i][j].NN];
		fread(REN,sizeof(int),SE[i][j].NN,fp);
		fclose(fp);

		//2 - создание векторов нагрузки
		int nvect;
		double **LVSE;
		nvect = pfm->nvect;
		LVSE = new double*[nvect];
		for (k=0; k<nvect; k++)
		{
			LVSE[k] = new double [SE[i][j].NN*pfm->KORT];
		}

		for (ii=0; ii<SE[i][j].NN; ii++)
		{
			for (jj=0; jj<pfm->KORT; jj++)
			{
				iur = REN[ii]*pfm->KORT + jj;

				for (k=0; k<nvect; k++)
				{
					LVSE[k][ii*pfm->KORT+jj] = pfm->LV[k][ iur ];				
				}

			}
		}

		sprintf(strto,"%s\\fff%s.bin",SE[i][j].fullnetfolder,SE[i][j].namenum);
		fp = fopen(strto,"wb");
		fwrite(&nvect,sizeof(int),1,fp);
		for (k=0; k<nvect; k++)
		{
			fwrite(LVSE[k],sizeof(double),SE[i][j].NN*pfm->KORT,fp);
		}
		fclose(fp);

		delete []REN;
		for (k=0; k<nvect; k++)
		{
			delete []LVSE[k];
		}
		delete []LVSE;

		//передача структуры СЭ


		// запуск расчета
		startSubspaceMatr(&SE[i][j]);

		//считывание найденной матрицы жесткости в подпространстве
		sprintf(strto,"%s\\sbsp%s.bin",SE[i][j].fullnetfolder,SE[i][j].namenum);
		fp  = fopen(strto,"rb");
		fread(&nvect,sizeof(int),1,fp);
		for(ii=0; ii<nvect; ii++)
		{
			fread(KK_1[ii],sizeof(double),nvect,fp);
		}
		fclose(fp);

		for (ii=0; ii<nvect; ii++)
		{
			for (jj=0; jj<nvect; jj++)
			{
				KK[ii][jj] += KK_1[ii][jj];
			}
		}

	}

	for (i=(pfm->nvect-1); i>=0; i--)
	{
		delete []KK_1[i];
	}
	delete []KK_1;

}

void startSubspaceMatr(SEstruct *SE)
{
	int i,j;

	CSE *pSE;

	pSE = new CSE[1];

	pSE->pmat = pmglob;

	pSE->Initialize(SE);
	pSE->MainSubspaceCalc();
	pSE->FULLDEL();

	delete []pSE;

}

void SEMODEL::SetNewShiftSEM(double shifting)
{
	int i,j;

	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			SE[i][j].shift = shifting;
			SE[i][j].is_matrix_calculated = false;
		}
	}
}

void SEMODEL::SetOperationType(int type)
{
	int i,j;

	for (i=0; i<SElevelsnum; i++)
	{
		for (j=0; j<numSEbylevels[i]; j++)
		{
			SE[i][j].cur_operation_type = type;
		}
	}
}