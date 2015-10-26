#include "StdAfx.h"
#include "L_FDPoint.h"
#include "FullModel.h"
#include "math.h"

FDPOINT::FDPOINT(void)
{
	ConnectedNodes = new int[20]; //память выделена с запасом, по наибольшему числу узлов из применяемых типов элементов
	FF = new double[20];
	ncn = 0;
	connected_el = -1;
	isactive = false;
}

FDPOINT::~FDPOINT(void)
{
	delete [] ConnectedNodes;
	delete [] FF;
}

void FDPOINT::TrackingPoint( void *vpfm, double *gcrd)
{
	if (isactive == true) AttachToFullModel(vpfm,gcrd,connected_el,&lcrd[0]);
}

void FDPOINT::TrackingPoint( void *vpfm, float *gcrdf)
{
	int k;
	double gcrd[3];
	for (k=0; k<3; k++) gcrd[k] = (double)gcrdf[k];
	if (isactive == true) AttachToFullModel(vpfm,&gcrd[0],connected_el,&lcrd[0]);
}

int FDPOINT::AttachToFullModel( void *vpfm, double *gcrd, int posel, double *poslcrd)
{
	//возвращает
	// 1 - если точка найдена в элементе posel;
	// 2 - если точка найдена в другом элементе;
	// 3 - если точка не попала ни в один элемент, но идентифицировано ее положение рядом с ближайшим к ней элементом
	// -1 - если точка совсем не найдена

	int i;
	int iel,ieltype;
	int glcretflag,fl;
	double FFtmp[20];
	double vmin[3], vmax[3], vcentr[3];
	double mindist,dist;
	int ielmindist;

	FULLMODEL *pfm;
	pfm = (FULLMODEL *)vpfm;

	for (i=0; i<pfm->KORT; i++)
	{
		crd[i] = gcrd[i];
	}

	//поиск начинается с элемента "posel", который из истории расчета может считаться наиболее вероятным местом положения точки приложения силы (определения перемещений)
	//за счет этого ускоряется поиск
	iel = 0;
	if ( posel >= 0 )
	{
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = poslcrd[i];
		}
		iel = posel;
	}
	else
	{
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = 0.0;
		}
	}

	pfm->ClearLIST_EL();

	pfm->LIST_EL[pfm->nlistel] = iel;
	pfm->MASK_EL[iel] = 1;
	pfm->nlistel++;

	int ilistel = 0; //текущий обработанный номер элемента в списке

	ielmindist = -1;
	mindist = 100000.0;
	while ( pfm->nlistel < pfm->NEL || ilistel < pfm->nlistel )
	{
		iel = pfm->LIST_EL[ilistel];
		
		ieltype = pfm->IND[iel][0];
		pfm->AttachElement(&pfm->el[ieltype],iel);
		pfm->el[ieltype].GetCentr(&vcentr[0]);
		dist = 0.0;
		for (i=0; i<pfm->KORT; i++)
		{
			dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
		}
		dist = sqrt(dist);
		if ( dist < mindist ) 
		{
			mindist = dist;
			ielmindist = iel;
		}
		pfm->el[ieltype].GetBoundBox(&vmin[0],&vmax[0]);

		//проверка попадает ли искомая точка в параллелепипед, ограничивающий данный элемент
		//если попадает - тогда есть смысл проводить подробный анализ
		fl = 0;
		for (i=0; i<pfm->KORT; i++)
		{
			if ( crd[i] < vmin[i] || crd[i] > vmax[i] )
			{
				fl = 1;
			}
		}
		if ( fl == 0 )
		{
			if ( ilistel > 0 ) //если при анализе первого элемента (наиболее вероятного) точка найдена не была
			{//то для следующих элементов начальное приближение локальных координат обнуляется
				for (i=0; i<pfm->KORT; i++)
				{
					lcrd[i] = 0.0;
				}
			}
			glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
			
			if ( glcretflag > 0 )
			{
				//точка идентифицирована
				connected_el = iel;
				ncn = pfm->el[ieltype].NN;
				//ConnectedNodes = new int[ncn];
				//FF = new double[ncn];
				for (i=0; i<ncn; i++)
				{
					ConnectedNodes[i] = pfm->IND[iel][i+2];
					FF[i] = FFtmp[i];
				}
				if ( ilistel == 0 ) return(1);
				else return(2);
			}
		}

		//добавление в список для анализа соседей текущего анализируемого элемента
		for (i=0; i<pfm->ELGR[iel][0]; i++)
		{
			if ( pfm->MASK_EL[ pfm->ELGR[iel][i+1] ] == 0 )
			{
				pfm->LIST_EL[ pfm->nlistel ] = pfm->ELGR[iel][i+1];
				pfm->nlistel++;
				pfm->MASK_EL[ pfm->ELGR[iel][i+1] ] = 1;
			}
		}
		ilistel++; //переход к следующему элементу из списка
	}
	// если программа не завершилась до этого места - значит элемент не находится ни в одном из КЭ
	// на заключительном этапе определяются локальные координаты точки для элемента, к которому она ближе всего расположена

	iel = ielmindist;
	ieltype = pfm->IND[iel][0];
	pfm->AttachElement(&pfm->el[ieltype],iel);
	for (i=0; i<pfm->KORT; i++)
	{
		lcrd[i] = 0.0;
	}
	glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
	connected_el = iel;
	ncn = pfm->el[ieltype].NN;
	//ConnectedNodes = new int[ncn];
	//FF = new double[ncn];
	for (i=0; i<ncn; i++)
	{
		ConnectedNodes[i] = pfm->IND[iel][i+2];
		FF[i] = FFtmp[i];
	}
	return(3);
}

void FDPOINT::TrackingExternPoint( void *vpfm, double *gcrd)
{
	int i,j;
	int iel,ieltype;
	double FFtmp[20];
	double vcentr[3];
	double mindist,dist;
	int ielmindist, ielmindist_old;
	int glcretflag;

	FULLMODEL *pfm;
	pfm = (FULLMODEL *)vpfm;

	for (i=0; i<pfm->KORT; i++)
	{
		crd[i] = gcrd[i];
	}

	//первичное обнаружение ближайшего КЭ и локальных координат
	if (connected_el == -1) AttachToFullModel(vpfm,gcrd,connected_el,&lcrd[0]);
	else
	{
		ielmindist_old = -2;
		ielmindist = connected_el;
		
		while (ielmindist_old != ielmindist) //процедра запускается несколько раз на тот случай, если искомая точка переместилась сильнее, чем на один уровень соседей
		{//предыдущего ближайщего конечного элемента
			mindist = 100000.0;
			ielmindist_old = ielmindist;
			//расстояние до КЭ, являвшегося ближайшим в прошлый раз
			iel = ielmindist_old;
			ieltype = pfm->IND[iel][0];
			pfm->AttachElement(&pfm->el[ieltype],iel);
			pfm->el[ieltype].GetCentr(&vcentr[0]);
			dist = 0.0;
			for (i=0; i<pfm->KORT; i++)
			{
				dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
			}
			dist = sqrt(dist);
			if ( dist < mindist ) 
			{
				mindist = dist;
				ielmindist = iel;
			}
			//поиск минимального расстояния до конечных элементов среди соседей предыдущего ближайшего КЭ
			for (j=1; j< (pfm->ELGR[ielmindist_old][0] + 1); j++)
			{
				iel = pfm->ELGR[ielmindist_old][j];
				ieltype = pfm->IND[iel][0];
				pfm->AttachElement(&pfm->el[ieltype],iel);
				pfm->el[ieltype].GetCentr(&vcentr[0]);
				dist = 0.0;
				for (i=0; i<pfm->KORT; i++)
				{
					dist += (vcentr[i] - crd[i])*(vcentr[i] - crd[i]);
				}
				dist = sqrt(dist);
				if ( dist < mindist ) 
				{
					mindist = dist;
					ielmindist = iel;
				}
			}
		}
		//для элемента, являющегося ближайшим, определяем локальные координаты внешней точки
		iel = ielmindist;
		ieltype = pfm->IND[iel][0];
		pfm->AttachElement(&pfm->el[ieltype],iel);
		for (i=0; i<pfm->KORT; i++)
		{
			lcrd[i] = 0.0;
		}
		glcretflag = pfm->el[ieltype].GetLocCrd(&crd[0],&lcrd[0],&FFtmp[0]);
		connected_el = iel;
		ncn = pfm->el[ieltype].NN;
		//ConnectedNodes = new int[ncn];
		//FF = new double[ncn];
		for (i=0; i<ncn; i++)
		{
			ConnectedNodes[i] = pfm->IND[iel][i+2];
			FF[i] = FFtmp[i];
		}
	}
}