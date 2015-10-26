#include "StdAfx.h"
#include "FullModel.h"
#include "math.h"
#include "stdio.h"

void FULLMODEL::GenSurf()
{
	//отбор граней основан на проверке первых 3х узлов, входящих в перечень грани
	int i,j,k;
	int nmaxfaces,inewface,eltype;

	InitIntPoint();

	nmaxfaces = 0;

	for (i=0; i<NEL; i++)
	{
		eltype = IND[i][0];
		nmaxfaces += el[eltype].Nface;
	}

	FACE *surfmodel_tmp;

	surfmodel_tmp = new FACE[nmaxfaces];

	inewface = 0;
	for (i=0; i<NEL; i++)
	{
		eltype = IND[i][0];
		AttachElement(&el[eltype],i);
		for (j=0; j<el[eltype].Nface; j++)
		{
			el[eltype].FaceNodeNums(j,&surfmodel_tmp[inewface].facenodenums[0],&surfmodel_tmp[inewface].NNface);

			if ( CheckNewFace(surfmodel_tmp,inewface) == 1 ) //такая ячейка не встречалась раньше
			{
				surfmodel_tmp[inewface].iel = i;
				surfmodel_tmp[inewface].isouter = true; //сначала по умолчанию ячейка считается внешней
				inewface++;
			}
		}
	}

	nfaces = 0;
	for (i=0; i<inewface; i++)
	{
		if ( surfmodel_tmp[i].isouter == true )
		{
			nfaces++;
		}
	}

	surfmodel = new FACE[nfaces];
	int iface = 0;

	for (i=0; i<inewface; i++)
	{
		if ( surfmodel_tmp[i].isouter == true )
		{ 
			surfmodel[iface].NNface = surfmodel_tmp[iface].NNface;
			for (j=0; j<surfmodel[iface].NNface; j++) surfmodel[iface].facenodenums[j] = surfmodel_tmp[i].facenodenums[j];
			surfmodel[iface].iel = surfmodel_tmp[i].iel;
			surfmodel[iface].isouter = surfmodel_tmp[i].isouter;			
			iface++;
		}
	}
	delete []surfmodel_tmp;

	char strl[256];
	FILE *fp;
	sprintf( strl,"%s\\SurfModel.vtk\0", pathmain );
	fp = fopen(strl,"w");
	ParaView_PrintSurfModel(fp);
	fclose(fp);
}

int FULLMODEL::CheckNewFace( FACE* surfmodel_tmp, int inewface)
{
	//возвращает 1 если новая грань не встречается в списке
	//0 если грань встретилась в списке
	int i,j,k,inode,nchn;
	i=0;

	if ( surfmodel_tmp[inewface].NNface > 4 ) nchn = (int) (surfmodel_tmp[inewface].NNface/2 + 0.1);
	else nchn = surfmodel_tmp[inewface].NNface; //количество сличаемых узлов (угловые узлы грани). Годится для линейных и квадратичных кубиков и тетраэдров.

	//программа может некорректно работать, если используются разные типы элементов одновременно

	for (i=0; i<inewface; i++)
	{
		if (surfmodel_tmp[i].isouter == true)
		{
			for (k=0; k<nchn; k++)
			{
				inode = surfmodel_tmp[inewface].facenodenums[k];
				for (j=0; j<nchn; j++)
				{
					if ( inode == surfmodel_tmp[i].facenodenums[j] ) break;
				}
				if (j == nchn) break; //узел не найден в списке узлов грани, значит новая грань точно не совпадает с гранью "i"
			}
			if ( k == nchn )
			{
				// грани совпали.
				// в этом случае ни новая грань, ни грань "i" не являются поверхностными
				surfmodel_tmp[i].isouter = false;
				return(0);
			}
		}
	}

	if (i == inewface) //новая грань не была найдена в списке
	{
		return(1);
	}

	//в случае неожиданной ошибки
	return(-1);
}
