// PSE_Client.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"

extern void usermain();


int _tmain(int argc, _TCHAR* argv[])
{
	usermain();	
	return 0;
}

////���������� ���������� ����� Pipe
//int _tmain(int argc, _TCHAR* argv[])
//{
//	int i,j,nmat,fwhile;
//	char strpath[256];
//	DWORD cb;
//
//	MEM *pmem = new MEM[1];
//
//	FULLMODEL *pfm;
//	SEMODEL *psem;
//	SSIEIGEN *peig;
//	ELASTAT *pest;
//
//	pest = NULL;
//	psem = NULL;
//	pfm = NULL;
//	peig = NULL;
//
//	char strpipecmd[256];
//	char strpipedata[256];
//	for (i=0; i<256; i++) 
//	{
//		strpipecmd[i] = '\0';
//		strpipedata[i] = '\0';
//	}
//
//	sprintf(strpipecmd,"\\\\.\\pipe\\PSEcmdpipe");
//
//	sprintf(strpipedata,"\\\\.\\pipe\\PSEdatapipe");
//
//	//PSE ������� ����� ������
//	HANDLE hPSEcmdpipe;
//	HANDLE hPSEdatapipe;
//	
//	hPSEcmdpipe = CreateNamedPipeA(LPCSTR(strpipecmd),
//		PIPE_ACCESS_DUPLEX, 
//		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
//		PIPE_UNLIMITED_INSTANCES,
//		512,512,5000,NULL);
//	hPSEdatapipe = CreateNamedPipeA(LPCSTR(strpipedata),
//		PIPE_ACCESS_DUPLEX, 
//		PIPE_TYPE_BYTE | PIPE_READMODE_BYTE | PIPE_WAIT,
//		PIPE_UNLIMITED_INSTANCES,
//		512,512,5000,NULL);
//
//	 if(hPSEcmdpipe == INVALID_HANDLE_VALUE)
//  {
//    printf("CreateNamedPipe: Error %ld\n", 
//      GetLastError());
//    getch();
//    return 0;
//  }
//
//	 	int message,message100 = 100;
//	//�������� ���������� �� ������� ������� ���������, ������������ ������� ��� �������
//	BOOL fl_is_connected_cmd = FALSE, fl_is_connected_data = FALSE;
//	fl_is_connected_cmd = ConnectNamedPipe(hPSEcmdpipe, NULL);
//	//������ �� ��������� �������� ���������� �� ������� ������
//	WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//
//	fl_is_connected_data = ConnectNamedPipe(hPSEdatapipe, NULL);
//	//������ �� ��������� �������� ���������� �� ������� ������
//	WriteFile(hPSEdatapipe,&message100,sizeof(int),&cb,NULL);
//
//	printf("\n1");
//	if ( fl_is_connected_cmd == TRUE && fl_is_connected_data == TRUE )
//	{
//		fwhile = 1;
//		while (fwhile == 1)
//		{
//			//�������� ���� ��������� ���������
//			printf("\nWaiting new command.....\n");
//			ReadFile(hPSEcmdpipe,&message,sizeof(int),&cb,NULL);
//			printf("\nmessage = %d\n",message);
//			switch (message)
//			{
//			case 1: //�������� � ���������� ������
//				int nchar;
//				ReadFile(hPSEdatapipe,&nchar,sizeof(int),&cb,NULL);//����� ������
//				ReadFile(hPSEdatapipe,&strpath[0],nchar*sizeof(char),&cb,NULL);// ���� � �����, ��� ����� ����������� ������
//				printf("Model receiving: %s\n",strpath);
//				pfm = new FULLMODEL[1];
//				pfm->MM = pmem;
//				sprintf(pfm->pathmain,"%s",strpath);
//				
//				ReadFile(hPSEdatapipe,&pfm->loadtype,sizeof(int),&cb,NULL);//��� �������� ������
//				
//				pfm->ReadFullmodel();
//				
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 2: //�������� ��������������� ������ �������������
//				psem = new SEMODEL[1];
//				psem->MM = pmem;
//				psem->pfm = pfm;
//
//				printf("SE model generation\n");
//				//psem->loadtype = 1; //�� UZOR
//				psem->loadtype = 2; //�������������� �������	
//				psem->ReadSEmodel();
//
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 3: //������������� ������������ �������
//				pest = new ELASTAT[1];
//				pest->MM = pmem;
//				pest->pfm = pfm;
//				pest->psem = psem;
//				printf("Static Calculation intit....\n");
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 20: //������� ������� �������� �����������
//				printf("Fixed displacement receiving.....\n");
//				//�������� ������� ������ �����������
//				ReadFile(hPSEdatapipe,pfm->FIX,pfm->NNE*sizeof(int),&cb,NULL); 
//
//				//�������� ������� �������� �������� �����������
//				ReadFile(hPSEdatapipe,pfm->UFIX,pfm->NNE*sizeof(double),&cb,NULL); 
//
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 21: //������� �������� ������� ���
//				int newnvect;
//				//�������� ���������� �������� ���
//				ReadFile(hPSEdatapipe,&newnvect,sizeof(int),&cb,NULL); 
//				printf("Forces receiving.....nvect = %d\n",newnvect);
//				pest->ReinitLVRV(newnvect);
//
//				//�������� �������� ���
//				for (i=0; i<pest->nvect; i++)
//				{
//					ReadFile(hPSEdatapipe,pest->LV[i],pfm->NNE*sizeof(double),&cb,NULL); 
//				}
//
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 22: //�������� �������� �������� ������� ����������� ����� �������
//				printf("Nodal displacements give back.....\n");
//				for (i=0; i<pfm->nvect; i++)
//				{
//					WriteFile(hPSEdatapipe,pest->RV[i],pfm->NNE*sizeof(double),&cb,NULL); 
//				}
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 40: //������ ������������ �������
//				printf("Starting static solution.....\n");
//				pest->StaticSolve();
//				
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 41: //���������� ������� �� ����������� �������
//				int nfr,nblock;
//				double eps;
//
//				//�������� ���������� ������, ������� ��������� �����
//				ReadFile(hPSEdatapipe,&nfr,sizeof(int),&cb,NULL); 
//				//�������� ��������
//				ReadFile(hPSEdatapipe,&eps,sizeof(double),&cb,NULL); 
//				//�������� ���������� ������ � �����
//				ReadFile(hPSEdatapipe,&nblock,sizeof(int),&cb,NULL); 
//				
//				printf("Starting Eigen solution  nfr = %d eps = %e nblock = %d\n",nfr,eps,nblock);
//
//				peig = new SSIEIGEN[1];
//				peig->MMM = pmem;
//				peig->pfm = pfm;
//				peig->psem = psem;
//				peig->Init(nfr,eps,nblock);
//				peig->MainSSI();
//				delete []peig;
//				peig = NULL;
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 80: //������� ������
//
//				printf("Memory full cleaning.....\n");
//				if (pest != NULL) delete []pest;
//				if (psem != NULL) delete []psem;
//				if (pfm != NULL) delete []pfm;
//
//				pest = NULL;
//				psem = NULL;
//				pfm = NULL;
//
//				//������ �� ��������� �������
//				WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//				break;
//			case 1000: //����� �� ���������
//				fwhile = 0;
//				printf("Exiting.....\n");
//				break;
//			default:
//				break;
//			}
//		}
//	}
//
//
//	if (pest != NULL) delete []pest;
//	if (psem != NULL) delete []psem;
//	if (pfm != NULL) delete []pfm;
//	if (peig != NULL) delete []peig;
//
//	pest = NULL;
//	psem = NULL;
//	pfm = NULL;
//	peig = NULL;
//
//	//������ �� ��������� �������
//	WriteFile(hPSEcmdpipe,&message100,sizeof(int),&cb,NULL);
//
//	DisconnectNamedPipe(hPSEcmdpipe);
//	DisconnectNamedPipe(hPSEdatapipe);
//
//	CloseHandle(hPSEcmdpipe);
//	CloseHandle(hPSEdatapipe);
//
//	if (pmem != NULL) delete []pmem;
//	pmem = NULL;
//
//	return 0;
//}
