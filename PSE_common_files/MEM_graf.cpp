#include "StdAfx.h"
#include "MEM.h"

void MEM::GRNDEL(int ***NDEL_, int **IND, int NEL, int NN )
{
	//создание массива, показывающего в какие конечные элементы входит каждый узел
	int i,j,k;
	int *ttt,ittt;
	int **NDEL;
	int tmpi;

	NDEL = new int*[NN];

	for(j=0;j<NN;j++)
	{
		NDEL[j] = new int[12];
		NDEL[j][0] = 0; //текущее количество зарегистрированных элементов
		NDEL[j][1] = 10; //текущее максимальное количество элементов
	}

	for(j=0;j<NEL;j++)
	{
		for(k=2; k<(IND[j][1]+2); k++)
		{
			if ( IND[j][k] < NN )
			{
				NDEL[IND[j][k]][ NDEL[IND[j][k]][0] + 2 ] = j;
				NDEL[IND[j][k]][0] ++;
				if ( NDEL[IND[j][k]][0] == NDEL[IND[j][k]][1] ) 
				{	//необходимо нарастить количество выделенной памяти, т.к. список элементов оказался длинным
					//без двойного копирования
					NDEL[IND[j][k]][1] += 10; //длина увеличивается на 10 
					ttt = new int [ NDEL[IND[j][k]][1] + 2 ]; 
					ittt = 0;
					for (i=0; i< NDEL[IND[j][k]][0]+2; i++)
					{
						ttt[ittt] = NDEL[IND[j][k]][i];
						ittt++;
					}
					delete [] NDEL[IND[j][k]];
					NDEL[IND[j][k]] = ttt; 
				}
			}
		}
	}

	*NDEL_ = NDEL;
}

void MEM::GenELGR(int ***GREL_, int NEL, int NN, int minnodes, int **NDEL, int **IND)
{	//создание несортированного графа по конечным элементам
	int i,j,k,iel,node,el;
	int **GREL;
	//создания массива флагов выбора и соответствующего массива выбранных элементов
	int *MASK,*MASKN;
	int *MASKLIST;
	MASK = new int [NEL];
	MASKLIST = new int [NEL];
	int imsl = 0;
	
	for(i=0; i<NEL; i++)
	{
		MASK[i] = 0;
	}

	//выделение памяти под граф по элементам
	GREL = new int*[NEL];

	//проход по конечным элементам
	for (iel = 0; iel<NEL; iel++)
	{
		//для каждого элемента составляем список других элементов, имеющих с ним хотя бы 1 общий узел
		//используя массив NDEL
		MASK[iel] = 1; //в список смежности не должен попасть сам номер текущего элемента
		for (i = 2; i < IND[iel][1]+2; i++)
		{
			node = IND[iel][i];
			for (j=2; j < NDEL[node][0]+2; j++)
			{ //выбираем элементы, в которые входит данный узел текущего КЭ
				el = NDEL[node][j];
				if ( MASK[el] == 0 )
				{
					MASK[el] = 1;
					MASKLIST[imsl] = el;
					imsl++;
				}
			}
		}
		
		GREL[iel] = new int[imsl+1];
		GREL[iel][0] = imsl;
		for (i=0; i<imsl; i++)
		{
			GREL[iel][i+1] = MASKLIST[i];
			MASK[MASKLIST[i]] = 0; //экономичная очистка флагов
		}
		MASK[iel] = 0; //очистка флага для текущего элемента
		imsl = 0; //сброс списка соседних элементов
	}

	if (minnodes > 1)
	{
		MASKN = new int[NN];
		for(i=0; i<NN; i++)
		{
			MASKN[i] = 0;
		}

		//оставляем в списках соседей только те элементы, которые имеют не менее minnodes общих узлов
		int curminnodes;
		for (iel = 0; iel<NEL; iel++)
		{
			for (j=2; j< (IND[iel][1] + 2); j++)
			{
				MASKN[ IND[iel][j] ] = 1;
			}
			//проверка на необходимое число общих узлов всех соседей из списка смежности
			imsl = 0;
			for (j=1; j< (GREL[iel][0] + 1); j++)
			{
				el = GREL[iel][j];
				curminnodes = 0;
				for (k=2; k< (IND[el][1] + 2); k++)
				{
					if ( MASKN[ IND[el][k] ] == 1 ) 
					{
						curminnodes++;
						if (curminnodes == minnodes) break;
					}
				}
				if (curminnodes >= minnodes)
				{
					MASKLIST[imsl] = el;
					imsl++;
				}
			}
			//обратное копирование оставшихся соседей
			delete []GREL[iel];
			GREL[iel] = new int[imsl+1];
			GREL[iel][0] = imsl;
			for (j=0; j<imsl; j++)
			{
				GREL[iel][j+1] = MASKLIST[j];
			}
			//очистка MASK
			for (j=2; j< (IND[iel][1] + 2); j++)
			{
				MASKN[ IND[iel][j] ] = 0;
			}
		}
		delete [] MASKN;
	}


	delete [] MASK;
	delete [] MASKLIST;

	*GREL_ = GREL;
}

void MEM::GenNDGR(int ***GRND_, int NEL, int NN, int NI, int **NDEL, int **IND)
{	
	//создание сортированного графа по NI первых узлов из NN узлов в сетке, определяемой матрицей индексов IND 
	int i,j,k,iel,node,el;
	int **GRND;
	//создания массива флагов выбора и соответствующего массива выбранных элементов
	int *MASK;
	int *MASKLIST;
	MASK = new int [NN];
	MASKLIST = new int [NN];
	int imsl = 0;
	
	for(i=0; i<NN; i++)
	{
		MASK[i] = 0;
	}

	//выделение памяти под граф 
	GRND = new int*[NI];

	for(i=0;i<NI;i++)
	{
		// для каждого узла проходим по всем элементам, в которые он входит (массив NDEL) и составляем список смежных узлов
		MASK[i] = 1; //текущий узел не должен попасть в список смежности
		imsl = 0;
		for (j = 2; j < (NDEL[i][0] + 2); j++)
		{
			el = NDEL[i][j];
			for (k = 2; k < (IND[el][1] + 2); k++)
			{
				node = IND[el][k];
				if ( MASK[node] == 0 && node < NI)
				{
					MASK[node] = 1;
					MASKLIST[imsl] = node;
					imsl++;
				}
			}
		}
		MASK[i] = 0; //текущий узел не должен попасть в список смежности

		//выделение памяти под строку графа
		GRND[i] = new int[imsl+1];
		GRND[i][0] = imsl;

		if ( ((float)imsl)/((float)NI) < 0.01 )
		{
			//если смежных узлов мало по сравнению с общим числом узлов модели, то выполняется сортировка смежных узлов и их запись в строку графа
			//обычно выполняется для моделей составленных из конечных элементов (не из суперэлементов)
			sort(MASKLIST,imsl);
			for (j = 0; j<imsl; j++)
			{
				GRND[i][j+1] = MASKLIST[j];
				MASK[ MASKLIST[j] ] = 0; //экономичная очистка массива флагов
			}
		}
		else
		{
			//если смежных узлов много - сортировка не выгодна
			//в этом случае порядок по возрастанию устанавливается однократным сквозным проходом по массиву MASK
			// эта ситуация типична в случае, когда модель составлена из суперэлементов
			imsl = 0;
			for ( j=0; j < NI; j++)
			{
				if ( MASK[j] == 1 )
				{
					GRND[i][imsl+1] = j;
					imsl++;
					MASK[j] = 0;
				}
			}
		}
	}


	delete [] MASK;
	delete [] MASKLIST;

	*GRND_ = GRND;
}

void MEM::sort (int *S,int s)
{
	int i,j,f,ns,min,i_min;
	
	ns=0;
	f=0;
	i_min=0;
	
	while (ns<s)
	{
		min=100000000;
		for(j=ns;j<s;j++)
		{
			if(S[j]<min) 
			{
				min=S[j];
				i_min=j;
			}
		}
		f=S[ns];
		S[ns]=min;
		S[i_min]=f;
		ns++;
	}
}
