#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"


void SEMODEL::AutoMLDivision_S()
{
	int i,j,k, NN, numNodes, startnode, neib,node,numSElastlev;
	int **GR = NULL,**NDEL = NULL;
	int *MASK;
	int *NEIB;
	int NNlev, nel, **IND_OLD;
	int *RENlev = NULL;
	int **INDnextlev = NULL; //матрица индексов следующего уровня 
	int NNnextlev; //число узлов в следующем уровне
	int ilev;
	char strname[256];

	sprintf(strname,"%s\\SELEVELS_listing.dat",pfm->pathmain);
	fp_ase = fopen(strname,"w");

	NN = pfm->NN;

	SEmaxLev = 10;
	int DS = 2000; //характерное число элементов в СЭ нулевого уровня
	int DSНB = 1000; //характерное число граничных узлов СЭ старших уровней
	int DSHT = 5000; //характерное полное число узлов СЭ старших уровней
	

	//создание информационной структуры о суперэлементах по максимально допустимому числу уровней, не все из них впоследствии будут заполнены
	numSEbylevels = new int [SEmaxLev];
	for (i=0; i<SEmaxLev; i++)
	{
		numSEbylevels[i] = 0;
	}
	SE = new SEstruct *[SEmaxLev];
	for (i=0; i<SEmaxLev; i++) SE[i] = NULL;
	Levels = new int*[SEmaxLev];
	for (i=0; i<SEmaxLev; i++) Levels[i] = NULL;
	RENlev = new int[NN]; //для всех уровней, максимальный размер.
	
	ilev = 0;
	numSElastlev = 2;
	while ( ilev < 10 && numSElastlev > 1 )
	{

		if (ilev == 0)
		{
			NNlev = pfm->NN;
			nel = pfm->NEL;
			IND_OLD = pfm->IND;
			for (i=0; i<NN; i++) //для нулевого уровня номера узлов уровня совпадают с номерами узлов исходной конструкции
			{
				RENlev[i] = i;
			}

			//AutoMLD_L0(&DS); //деление нулевого уровня (заполнение соответствующей строки в Levels)
			AutoMLD_L0_EL(&DS); //деление нулевого уровня (заполнение соответствующей строки в Levels)
		}
		else
		{
			IND_OLD = MM->MEM_DEL(IND_OLD,nel,1);
			nel = numSEbylevels[ilev-1];
			NNlev = NNnextlev;
			NNnextlev = 0;
			IND_OLD = INDnextlev;
			INDnextlev = NULL;

			AutoMLD_L1(ilev,DSНB,DSHT,IND_OLD,nel,NNlev); //деление старших уровней (заполнение соответствующей строки в Levels)
		}

		//вывод исходных данных по уровню и переход к следующему уровню
		NNnextlev = 0;
		INDnextlev = new int *[ numSEbylevels[ilev] ];
		for (i=0; i<numSEbylevels[ilev]; i++) INDnextlev[i] = NULL;

		//вывод сформированного уровня суперэлементов, на основе информации из Levels
		//расчет матрицы индексов конструкции, составленной из этих суперэлементов, и вспомогательных массивов перенумерации
		AutoMLDivision_LevelOutChange(ilev, numSEbylevels[ilev], RENlev, NNlev, nel, IND_OLD, INDnextlev, &NNnextlev);

		if (ilev == 0) IND_OLD = NULL; //память по этому указателю выделена в pfm->IND для нулевого уровня и ее не надо удалять

		numSElastlev = numSEbylevels[ilev];
		//приращение номера уровня
		ilev++;
	}

	SElevelsnum = ilev; //количество сформированных уровней

	IND_OLD = MM->MEM_DEL(IND_OLD,nel,1);
	INDnextlev = MM->MEM_DEL(INDnextlev,numSElastlev,1);
	delete [] RENlev;	


	fclose(fp_ase);
}






//void SEMODEL::FindNotNumNeib(int **GR, int *MASK, int iSE, int *NEIB, int *neib, int DS, int startnode)
//{
//	int i,j,k,nb,nban,ib;
//	int *NEIBCAND,*DEGREE;
//	int nc;
//
//	NEIBCAND = new int[pfm->NN];
//	DEGREE = new int[pfm->NN];
//	nc = 0;
//	for (i=0; i<pfm->NN; i++)
//	{
//		DEGREE[i] = 0;
//	}
//
//	NEIBCAND[0] = startnode;
//	NEIB[0] = startnode;
//	MASK[ NEIB[0] ] = -2;
//	nb = 1;
//
//	cur_neib_num = 0;
//
//	for (ib = 0; ib < nb; ib++)
//	{
//		for (j=1; j< (GR[ NEIBCAND[cur_neib_num] ][0] + 1); j++)
//		{
//			if ( MASK [ GR[ NEIBCAND[cur_neib_num] ][j] ] == -1 ) //узел пока никуда не относится
//			{
//				NEIBCAND[nc] = GR[ NEIBCAND[cur_neib_num] ][j];
//				MASK [ NEIBCAND[nc] ] = -3; //узел становится кандидатом в новый суперэлемент
//				nc++;
//			}
//		}
//		//из узлов-кандидатов нужно выбрать для включения в СЭ тот, который обладает наибольшим числом соседей из числа узлов, включенных в СЭ ранее
//		if ( nb > DS )
//		{
//			break;
//		}
//	}
//	
//	*neib = nb;
//
//}