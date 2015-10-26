#include "StdAfx.h"
#include "SEmodel.h"
#include "math.h"
#include "stdio.h"

void SEMODEL::AutoMLD_L1(int ilev, int maxboundary_size, int maxtotal_size, int **IND_OLD, int nel, int NNlev )
{
	int i,j,k,kk;
	int nbn,ntotaln,nin;
	int istartse = 0;
	int icurnewse = 0;
	int nincludedse = 0;
	int maxbn, imaxbn, isemaxbn,curbn,ise;
	int **GREL = NULL;
	int **NDEL = NULL;
	int **SEBND = NULL;
	int *MASK = NULL;
	int *NEIB = NULL;
	int neib;

	// расчет графа по суперэлементам
	
	MM->GRNDEL(&NDEL,IND_OLD,nel,NNlev);
	MM->GenELGR(&GREL,nel,NNlev,1,NDEL,IND_OLD);

	//определение числа граничных узлов для каждого суперэлемента и каждой пары составляющих суперэлементов прошлого уровня
	FindSeBound(&SEBND,GREL,NDEL,NNlev,nel,ilev-1);

	//логика формирования СЭ:
	// выбирается стартовый СЭ. В процессе формирования новых СЭ к нему присоединяются СЭ из списка его соседей 1ого и более высоких уровней
	// для присоединения выбираются в первую очередь те СЭ из списка соседей, которые имеют больше граничных узлов с ранее отобранными СЭ
	// проверяется ограничение на количество граничных узлов (и, возможно, на полное количество узлов)
	// ограничение на количество граничных узлов должно не менее чем вдвое (лучше в 2.5 - 3 раза) превышать максимальное число граничных узлов среди СЭ прошлого уровня  

	// в качестве первого СЭ лучше всего выбрать СЭ с наибольшим количеством граничных узлов

	MASK = new int[nel];
	NEIB = new int[nel];
	for (i=0; i<nel; i++)
	{
		MASK[i] = -1; //целое число показывает в какой новый номер СЭ уже включен данный суперэлемент
	}
	neib = 0; //в NEIB будут храниться текущие номера СЭ, отобранных для формирования нового СЭ  
	nincludedse = 0;
	while ( nincludedse < nel )
	{
		//начало формирования очередного нового СЭ
		neib = 0;
		NEIB[neib] = istartse;
		neib++;
		nincludedse++;
		MASK[istartse] = icurnewse;

		nbn = 0; // текущее количество граничных узлов
		ntotaln = 0; //текущее полное число узлов
		nin = 0; //текущее количество внутренних узлов
		imaxbn = 1;
		while ( nbn < maxboundary_size && ntotaln < maxtotal_size && imaxbn != -1 )
		{
			// анализ граничных узлов с соседями всех уровней
			maxbn = 0;
			imaxbn = -1;
			
			for ( i=0; i < neib; i++)
			{
				for ( j=1; j < (GREL[ NEIB[i] ][0] + 1); j++)
				{
					ise = GREL[ NEIB[i] ][j];
					if ( MASK[ ise ] == -1 ) //если этот суперэлемент еще никуда не включен
					{
						//определяем число его узлов, граничащих с отобранными СЭ
						curbn = 0;
						for (k=0; k<neib; k++)
						{
							for (kk=1; kk<(GREL[ise][0]+1); kk++)
							{
								if ( GREL[ise][kk] == NEIB[k] )
								{
									curbn += SEBND[ise][kk];
								}
							}
						}
						if ( curbn > maxbn )
						{
							maxbn = curbn;
							imaxbn = ise;
						}
					}
				}
			}
			if ( imaxbn != -1 ) //не осталось непронумерованных СЭ в области связанной с новым формируемым СЭ
			{
				// выбранный СЭ включаем в NEIB и перерассчитываем число граничных узлов нового СЭ
				NEIB[neib] = imaxbn;
				neib++;
				nincludedse++;
				MASK[imaxbn] = icurnewse;

				nbn = 0;
				ntotaln = 0;
				nin = 0;
				//расчет количества граничных узлов осуществляется так:
				// сначала вычисляется сумма значений количества граничных узлов между каждой парой входящих СЭ
				for (i=0; i<neib; i++)
				{
					for (j=i+1; j<neib; j++)
					{
						for (k=1; k<(GREL[ NEIB[i] ][0] + 1); k++)
						{
							if ( GREL[ NEIB[i] ][k] == NEIB[j] )
							{
								nin += SEBND[ NEIB[i] ][k];
							}
						}
					}
					ntotaln += SEBND[ NEIB[i] ][0];
				}
				ntotaln -= nin;
				nbn = ntotaln - nin;
			}
		}

		// после того, как новый суперэлемент собран, ищем новый стартовый СЭ
		istartse = -1;
		
		while ( istartse == -1 && nincludedse < nel )
		{
			//сначала проверяем список СЭ, соседних с набором СЭ, вошедших в состав последнего сформированного СЭ нового уровня
			for (i=0; i<neib; i++)
			{
				for (k=1; k<(GREL[ NEIB[i] ][0] + 1); k++)
				{
					if ( MASK[ GREL[ NEIB[i] ][k] ] == -1 && istartse == -1 )
					{
						istartse = GREL[ NEIB[i] ][k];
					}
				}
			}
			//если среди списка соседей не нашлось ни одного непронумерованного суперэлемента
			// требуется искать любой другой непронумерованный элемент
			if (istartse == -1)
			{
				for (i=0; i<nel; i++)
				{
					if ( MASK[i] == -1 ) break;
				}
				if ( i < nel ) //непронумерованные еще остались
				{
					istartse = i;
				}
			}

			//требуется проверить, не является ли новый стартовый СЭ изолированным
			if (istartse != -1)
			{
				for (j=1; j < (GREL[istartse][0] + 1); j++)
				{
					if ( MASK[ GREL[istartse][j] ] == -1 ) break;
				}
				if ( j == (GREL[istartse][0] + 1) )
				{
					//если элемент изолирован, требуется его присоединить к одному из соседних с ним новых суперэлементов
					// можно усовершенствовать - присоединять так, чтобы получилось минимальное число граничных узлов
					// в данной версии присоеднияем к тому СЭ с которым больше всего граничных узлов.
					curbn = 0;
					imaxbn = -1;
					for ( k=1; k < (GREL[istartse][0]+1); k++)
					{
						if ( SEBND[istartse][k] > curbn )
						{
							curbn = SEBND[istartse][k];
							imaxbn = k;
						}
					}
					//суперэлемент i присоединяем к тому же суперэлементу, к которому присоединен GREL[istartse][k]
					MASK[ istartse ] = MASK[ GREL[istartse][imaxbn] ];
					nincludedse++;
					istartse = -1; //теперь требуется поиск нового стартового СЭ
				}
			}
		}

		// и очищаем список сформировавших его СЭ
		neib = 0;

		icurnewse++;
	}

	//в icurnewse содержится число сформированных суперэлементов нового уровня
	numSEbylevels[ilev] = icurnewse;
	//вывод разбиения и переход к следующему уровню
	Levels[ilev] = new int[nel];
	for (i=0; i<nel; i++)
	{
		Levels[ilev][i] = MASK[i];
	}
	
	delete []NEIB;
	delete []MASK;
	SEBND = MM->MEM_DEL(SEBND,nel,1);
	GREL = MM->MEM_DEL(GREL,nel,1);
	NDEL = MM->MEM_DEL(NDEL,NNlev,1);

}

		//if ( neib == 1 ) //стартовый СЭ оказался единственным непронумерованным в текущей области. Его треубется присоединить к любому новому сформированному соседнему СЭ
		//{
		//	for (k=1; k<(GREL[ NEIB[0] ][0] + 1); k++)
		//	{
		//		if ( MASK[ GREL[ NEIB[0] ][k] ] != -1 ) break;
		//	}
		//	if ( k != (GREL[ NEIB[0] ][0] + 1) )
		//	{
		//		MASK[ NEIB[0] ] = MASK[ GREL[ NEIB[0] ][k] ];
		//	}
		//	else
		//	{
		//		k=k; //ошибка, нет пронумерованных соседних СЭ
		//	}
		//}

void SEMODEL::FindSeBound(int ***SEBND_, int **GREL, int **NDEL, int NNlev, int nel, int ilev)
{
	int i,j,k,in,ise1,ise2;
	int **SEBND=NULL;

	//массив информации о граничных узлах соответствует по структуре и должен использоваться совместно с графом по СЭ
	SEBND = new int *[nel];
	for (i=0; i<nel; i++)
	{
		SEBND[i] = new int [GREL[i][0]+1]; // первое число в строке - общее число граничных узлов суперэлемента
		for (j=0; j< (GREL[i][0]+1); j++)
		{
			SEBND[i][j] = 0;
		}
		SEBND[i][0] = SE[ilev][i].NS;
	}

	for (in = 0; in < NNlev; in++)
	{
		//для каждого узла осуществляется анализ списка смежных СЭ из NDEL
		//для каждой пары СЭ из списка смежности добавляется 1 в соответствующие ячейки SEBND

		for ( i = 2; i < (NDEL[in][0]+2); i++)
		{
			for (j = i+1; j < (NDEL[in][0]+2); j++)
			{
				ise1 = NDEL[in][i];
				ise2 = NDEL[in][j];
				
				for ( k = 1; k < (GREL[ise1][0]+1); k++)
				{
					if ( GREL[ise1][k] == ise2 ) break;
				}
				if ( k == (GREL[ise1][0]+1) )
				{
					k=k; //ошибка: суперэлемент не найден в графе, хотя имеет общий узел с текущим СЭ
				}
				else
				{
					SEBND[ise1][k]++;
				}
				//вычисление парного элемента
				for ( k = 1; k < (GREL[ise2][0]+1); k++)
				{
					if ( GREL[ise2][k] == ise1 ) break;
				}
				if ( k == (GREL[ise2][0]+1) )
				{
					k=k; //ошибка: суперэлемент не найден в графе, хотя имеет общий узел с текущим СЭ
				}
				else
				{
					SEBND[ise2][k]++;
				}
			}
		}
	}
	*SEBND_ = SEBND;
}
