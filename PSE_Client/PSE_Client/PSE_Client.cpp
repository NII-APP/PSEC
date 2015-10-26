// PSE_Client.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"
//
//int _tmain(int argc, _TCHAR* argv[])
//{
//	int i,j,nmat;
//	char strpath[256];
//
//	sprintf(strpath,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
//	//sprintf(strpath,"C:\\Temp\\UGM_LAST\\STATOR");
//	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE_3levels");
//	//sprintf(strpath,"C:\\Temp\\Test_Shevchenko");
//
//	MEM *pmem = new MEM[1];
//
//	FULLMODEL *pfm;
//	SEMODEL *psem;
//	MATPROP *pm;
//	ELASTAT *pest;
//
//	nmat = 1;
//	pm = new MATPROP[nmat];
//	pm->E = 2.1e11;
//	pm->MU = 0.3;
//	pm->RO = 7800;
//
//	pfm = new FULLMODEL[1];
//	pfm->MM = pmem;
//	sprintf(pfm->pathmain,"%s",strpath);
//	pfm->nmat = nmat;
//	pfm->mat = pm;
//	pfm->loadtype = 2;
//	pfm->ReadFullmodel();
//
//
//	psem = new SEMODEL[1];
//	psem->MM = pmem;
//	psem->pfm = pfm;
//	psem->loadtype = 2;
//	
//	psem->ReadSEmodel();
//	
//	//psem->AutoMLDivision_S();
//
//	pest = new ELASTAT[1];
//	pest->MM = pmem;
//	pest->pfm = pfm;
//	pest->psem = psem;
//	pest->CreateLV();
//	pest->StaticSolve();
//
//
//	i=0;
//
//	delete []pest;
//	delete []psem;
//	delete []pfm;
//	delete []pm;
//	delete []pmem;
//	
//	return 0;
//}
//

////тест (расчет СЭ модели и собственных частот)
//int _tmain(int argc, _TCHAR* argv[])
//{
//	int i,j,nmat;
//	char strpath[256];
//
//	sprintf(strpath,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
//	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE_3levels");
//	//sprintf(strpath,"C:\\Temp\\Test_Shevchenko");
//
////	sprintf(strpath,"C:\\Temp\\UGM_LAST\\STATOR");
//	//sprintf(strpath,"C:\\Temp\\Reshetka");
//sprintf(strpath,"C:\\Temp\\simple1_out_model");
////sprintf(strpath,"C:\\Temp\\3DCUT_TEMP_1_dyntestMC\\FEMCUT");
//
//	MEM *pmem = new MEM[1];
//
//	FULLMODEL *pfm;
//	SEMODEL *psem;
//	MATPROP *pm;
//	SSIEIGEN *peig;
//
//	//nmat = 1;
//	//pm = new MATPROP[nmat];
//	//pm->E = 2.1e11;
//	//pm->MU = 0.3;
//	//pm->RO = 7800;
//
//	pfm = new FULLMODEL[1];
//	pfm->MM = pmem;
//	sprintf(pfm->pathmain,"%s",strpath);
////	pfm->nmat = nmat;
////	pfm->mat = pm;
////	pfm->loadtype = 2; //UZOR новый формат
//	pfm->loadtype = 3; //Ansys Full Grid
//
//	pfm->ReadFullmodel();
//
//
//	psem = new SEMODEL[1];
//	psem->MM = pmem;
//	psem->pfm = pfm;
//	//psem->loadtype = 1; //из UZOR
//	psem->loadtype = 2; //автоматическое деление	
//	psem->ReadSEmodel();
//
//
//	peig = new SSIEIGEN[1];
//	peig->MMM = pmem;
//	peig->pfm = pfm;
//	peig->psem = psem;
//	peig->Init(2,0.00001,15);
//	peig->MainSSI();
//	
//
//
//	i=0;
//
//	delete []peig;
//	delete []psem;
//	delete []pfm;
////	delete []pm;
//	delete []pmem;
//	
//	return 0;
//}

////тест интегрирования движения
//int _tmain(int argc, _TCHAR* argv[])
//{
//	int i,j,nmat;
//	char strpath[256];
//
//	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE\\Test_eigen_SE");
//	//sprintf(strpath,"C:\\Temp\\Test_eigen_SE_3levels");
//	//sprintf(strpath,"C:\\Temp\\Test_Shevchenko");
//
////	sprintf(strpath,"C:\\Temp\\UGM_LAST\\STATOR");
//	//sprintf(strpath,"C:\\Temp\\Reshetka");
//sprintf(strpath,"C:\\Temp\\simple1_out_model");
//
//	MEM *pmem = new MEM[1];
//
//	FULLMODEL *pfm;
//	SEMODEL *psem;
//	MATPROP *pm;
//	SSIEIGEN *peig;
//
//	//nmat = 1;
//	//pm = new MATPROP[nmat];
//	//pm->E = 2.1e11;
//	//pm->MU = 0.3;
//	//pm->RO = 7800;
//
//	pfm = new FULLMODEL[1];
//	pfm->MM = pmem;
//	sprintf(pfm->pathmain,"%s",strpath);
////	pfm->nmat = nmat;
////	pfm->mat = pm;
////	pfm->loadtype = 2; //UZOR новый формат
//	pfm->loadtype = 3; //Ansys Full Grid
//
//	pfm->ReadFullmodel();
//	pfm->GenSurf();
//
//
//	psem = new SEMODEL[1];
//	psem->MM = pmem;
//	psem->pfm = pfm;
//	//psem->loadtype = 1; //из UZOR
//	psem->loadtype = 2; //автоматическое деление	
//	psem->ReadSEmodel();
//
//
//	peig = new SSIEIGEN[1];
//	peig->MMM = pmem;
//	peig->pfm = pfm;
//	peig->psem = psem;
//	peig->Init(5,0.00001,15);
//	peig->MainSSI();
//	
//	MOTIONEIGF *pmef;
//	pmef = new MOTIONEIGF[1];
//	pmef->MM = pmem;
//	pmef->AttachToFullModel(pfm);
//	pmef->LoadEIGF();
//	pmef->InitConstModalDamp(0.1);
//	pmef->InitIntegr();
//	pmef->InitInertiaUnitModalForces();
//	pmef->TestCase_InertiaLoading();
//
//	i=0;
//
//	delete []peig;
//	delete []psem;
//	delete []pfm;
////	delete []pm;
//	delete []pmem;
//	delete []pmef;
//	
//	return 0;
//}

//управление программой через Pipe
int _tmain(int argc, _TCHAR* argv[])
{
	int i,j,nmat,fwhile;
	char strpath[256];
	DWORD cb;

	MEM *pmem = new MEM[1];

	FULLMODEL *pfm;
	SEMODEL *psem;
	SSIEIGEN *peig;
	ELASTAT *pest;

	pest = NULL;
	psem = NULL;
	pfm = NULL;
	peig = NULL;

	char strpipecmd[256];
	char strpipedata[256];
	for (i=0; i<256; i++) 
	{
		strpipecmd[i] = '\0';
		strpipedata[i] = '\0';
	}

	sprintf(strpipecmd,"\\\\.\\pipe\\PSEcmdpipe");

	sprintf(strpipedata,"\\\\.\\pipe\\PSEdatapipe");

	//PSE создает канал первой
	HANDLE hPSEcmdpipe;
	HANDLE hPSEdatapipe;
	
	hPSEcmdpipe = CreateNamedPipeA(LPCSTR(strpipecmd),
		PIPE_ACCESS_DUPLEX, 
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
		PIPE_UNLIMITED_INSTANCES,
		512,512,5000,NULL);
	hPSEdatapipe = CreateNamedPipeA(LPCSTR(strpipedata),
		PIPE_ACCESS_DUPLEX, 
		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
		PIPE_UNLIMITED_INSTANCES,
		512,512,5000,NULL);

	 if(hPSEcmdpipe == INVALID_HANDLE_VALUE)
  {
    printf("CreateNamedPipe: Error %ld\n", 
      GetLastError());
    getch();
    return 0;
  }

	 	int message,message100 = 100;
	//Ожидание соединения со стороны внешней программы, определяющей задания для расчета
	BOOL fl_is_connected_cmd = FALSE, fl_is_connected_data = FALSE;
	fl_is_connected_cmd = ConnectNamedPipe(hPSEcmdpipe, NULL);
	//сигнал об окончании процесса соединения по данному каналу
	WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);

	fl_is_connected_data = ConnectNamedPipe(hPSEdatapipe, NULL);
	//сигнал об окончании процесса соединения по данному каналу
	WriteFile(hPSEdatapipe,&message100,sizeof(int),&cb,NULL);

	printf("\n1");
	if ( fl_is_connected_cmd == TRUE && fl_is_connected_data == TRUE )
	{
		fwhile = 1;
		while (fwhile == 1)
		{
			//основной цикл обработки сообщений
			printf("\nWaiting new command.....\n");
			ReadFile(hPSEcmdpipe,&message,sizeof(int),&cb,NULL);
			printf("\nmessage = %d\n",message);
			switch (message)
			{
			case 1: //создание и считывание модели
				int nchar;
				ReadFile(hPSEdatapipe,&nchar,sizeof(int),&cb,NULL);//длина строки
				ReadFile(hPSEdatapipe,&strpath[0],nchar*sizeof(char),&cb,NULL);// путь к папке, где нужно осуществить расчет
				printf("Model receiving: %s\n",strpath);
				pfm = new FULLMODEL[1];
				pfm->MM = pmem;
				sprintf(pfm->pathmain,"%s",strpath);
				
				ReadFile(hPSEdatapipe,&pfm->loadtype,sizeof(int),&cb,NULL);//тип загрузки данных
				
				pfm->ReadFullmodel();
				
				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 2: //создание суперэлементной модели автоматически
				psem = new SEMODEL[1];
				psem->MM = pmem;
				psem->pfm = pfm;

				printf("SE model generation\n");
				//psem->loadtype = 1; //из UZOR
				psem->loadtype = 2; //автоматическое деление	
				psem->ReadSEmodel();

				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 3: //инициализация статического расчета
				pest = new ELASTAT[1];
				pest->MM = pmem;
				pest->pfm = pfm;
				pest->psem = psem;
				printf("Static Calculation intit....\n");
				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 20: //задание вектора заданных перемещений
				printf("Fixed displacement receiving.....\n");
				//передача вектора флагов закреплений
				ReadFile(hPSEdatapipe,pfm->FIX,pfm->NNE*sizeof(int),&cb,NULL); 

				//передача вектора значений заданных перемещений
				ReadFile(hPSEdatapipe,pfm->UFIX,pfm->NNE*sizeof(double),&cb,NULL); 

				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 21: //задание векторов узловых сил
				int newnvect;
				//передача количества векторов сил
				ReadFile(hPSEdatapipe,&newnvect,sizeof(int),&cb,NULL); 
				printf("Forces receiving.....nvect = %d\n",newnvect);
				pest->ReinitLVRV(newnvect);

				//передача векторов сил
				for (i=0; i<pest->nvect; i++)
				{
					ReadFile(hPSEdatapipe,pest->LV[i],pfm->NNE*sizeof(double),&cb,NULL); 
				}

				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 22: //обратная передача векторов узловых перемещений после расчета
				printf("Nodal displacements give back.....\n");
				for (i=0; i<pfm->nvect; i++)
				{
					WriteFile(hPSEdatapipe,pest->RV[i],pfm->NNE*sizeof(double),&cb,NULL); 
				}
				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 40: //запуск статического расчета
				printf("Starting static solution.....\n");
				pest->StaticSolve();
				
				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 41: //выполнение расчета на собственные частоты
				int nfr,nblock;
				double eps;

				//передача количества частот, которые требуется найти
				ReadFile(hPSEdatapipe,&nfr,sizeof(int),&cb,NULL); 
				//передача точности
				ReadFile(hPSEdatapipe,&eps,sizeof(double),&cb,NULL); 
				//передача количества частот в блоке
				ReadFile(hPSEdatapipe,&nblock,sizeof(int),&cb,NULL); 
				
				printf("Starting Eigen solution  nfr = %d eps = %e nblock = %d\n",nfr,eps,nblock);

				peig = new SSIEIGEN[1];
				peig->MMM = pmem;
				peig->pfm = pfm;
				peig->psem = psem;
				peig->Init(nfr,eps,nblock);
				peig->MainSSI();
				delete []peig;
				peig = NULL;
				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 80: //очистка памяти

				printf("Memory full cleaning.....\n");
				if (pest != NULL) delete []pest;
				if (psem != NULL) delete []psem;
				if (pfm != NULL) delete []pfm;

				pest = NULL;
				psem = NULL;
				pfm = NULL;

				//сигнал об окончании команды
				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
				break;
			case 1000: //выход из программы
				fwhile = 0;
				printf("Exiting.....\n");
				break;
			default:
				break;
			}
		}
	}


	if (pest != NULL) delete []pest;
	if (psem != NULL) delete []psem;
	if (pfm != NULL) delete []pfm;
	if (peig != NULL) delete []peig;

	pest = NULL;
	psem = NULL;
	pfm = NULL;
	peig = NULL;

	//сигнал об окончании команды
	WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);

	DisconnectNamedPipe(hPSEcmdpipe);
	DisconnectNamedPipe(hPSEdatapipe);

	CloseHandle(hPSEcmdpipe);
	CloseHandle(hPSEdatapipe);

	if (pmem != NULL) delete []pmem;
	pmem = NULL;

	return 0;
}
