#include "StdAfx.h"
#include "SE.h"

//CLIENTDATA cldat;


void CSE::StifIIAssembling(CEL *el,int elnum)
{
	int i,j,iur,jur,tmp;
	
	for(i=0; i<el->NNE; i++)
	{
		for(j=i; j<el->NNE; j++)
		{
			iur = el->GLOBNE[i];
			jur = el->GLOBNE[j];

			//if (iur <= jur) STIF[jur][jur-iur] += el->GE[i][j];
			//else STIF[iur][iur-jur] += el->GE[i][j];
			if ( (iur < NIE) && (jur < NIE) )
			{
				if (iur <= jur) STII[ STII_ENV[jur] - (jur-iur) ] += el->STIF[i][j];
				else STII[ STII_ENV[iur] - (iur-jur) ] += el->STIF[i][j];	
			}
		}
	}
}


void CSE::StifIINZAssembling(CEL *el,int elnum, int flagfull)
{
	int ii,jj,i,j,in,jn,ik,jk,k,iur,jur,tmpi;

	int NIE,NI;

	if (flagfull == 0)
	{
		NIE = this->NIE;
		NI = this->NI;
	}
	if (flagfull == 1)
	{
		NIE = this->NNE;
		NI = this->NN;
	}

	for(i=0; i<el->NNE; i++)
	{
		for(j=i; j<el->NNE; j++)
		{
			iur = el->GLOBNE[i];// строка
			jur = el->GLOBNE[j];// столбец

			if (jur < iur)
			{
				tmpi = jur;
				jur = iur;
				iur = tmpi;
			}

			if ( (iur < NIE) && (jur < NIE) )
			{
				//поиск позиции в соответствующем столбце матрицы жесткости
				for (k=0; k<STII_ENVNZ[jur][0]; k++)
				{
					if ( STII_ENVNZ[jur][k+1] == iur ) break;
				}

				if (k == STII_ENVNZ[jur][0])
				{
					// Error
					k = k;
				}
				else
				{
					STIINZ[jur][k] += el->STIF[i][j];
				}
			}
		}
	}

}


void CSE::StifISAssembling(CEL *el,int elnum)
{
	int i,j,k,env,iur,jur,tmp,inode,jnode,nnode,iurenv,nnodesum;
	
	for(i=0; i<el->NNE; i++)
	{
		for(j=i; j<el->NNE; j++)
		{
			iur = el->GLOBNE[i];
			jur = el->GLOBNE[j];

			if ( iur > jur ) //заполняем верхний треугольник
			{
				tmp = iur;
				iur = jur;
				jur = tmp;
			}

			if ( (jur >= NIE) && (iur < NIE) ) //проверка попадания в IS часть
			{
				inode = (int)(iur/KORT);
				jnode = (int)((jur-NIE)/KORT);

				nnodesum = 0;

				env = STIS_ENV[jnode][1];
				for (k=0; k<env; k+=2)
				{
					/*if ( k > 0 )
					{
						nnode = STIS_ENV[jnode][k+2+1] - STIS_ENV[jnode][k+1];
					}
					else
					{
						nnode = STIS_ENV[jnode][k+2+1];
					}*/
					nnode = STIS_ENV[jnode][k+2+1];
					if ( (inode >= STIS_ENV[jnode][k+2]) && (inode < (STIS_ENV[jnode][k+2] + nnode)) ) break; //найден нужный интервал ненулей в столбце
					nnodesum += nnode;
				}

				/*if ( k > 0 )
				{
					iurenv = STIS_ENV[jnode][k+1]*KORT;
				}
				else
				{
					iurenv = 0;
				}
				iurenv =  iurenv + (iur - STIS_ENV[jnode][k+2]*KORT);*/
				iurenv =  nnodesum*KORT + (iur - STIS_ENV[jnode][k+2]*KORT);

				if (iurenv >= STIS_ENV[jnode][0]*KORT)
				{
					iurenv = iurenv;
				}

				STIS[jur-NIE][iurenv] += el->STIF[i][j];
			}
		}
	}
}

void CSE::StifSSAssembling(CEL *el,int elnum)
{
	//заполняется полный треугольник
	int i,j,iur,jur,tmp,ip,iurl, jurl;
	
	for(i=0; i<el->NNE; i++)
	{
		for(j=i; j<el->NNE; j++)
		{
			iur = el->GLOBNE[i];
			jur = el->GLOBNE[j];

			if ( iur > jur ) //заполняем верхний треугольник
			{
				tmp = iur;
				iur = jur;
				jur = tmp;
			}
			iurl = iur - NIE;
			jurl = jur - NIE;

			if ( (iurl >= 0) && (jurl >= 0) )
			{
				if (jurl > 0)
				{
					ip = (int)((jurl*jurl-jurl)/2 + 0.01);
					ip += jurl;
					ip += iurl;
				}
				else
				{
					ip = 0;
				}

				STSS[ip] += el->STIF[i][j];
			}
		}
	}
}


void CSE::StifIIAssemblingFromSE(double *STSSse, int KU, int *NUR )
{
	int i,j,iur,jur,tmp,ip;
	
	for(i=0; i<KU; i++)
	{
		for(j=i; j<KU; j++)
		{
			iur = NUR[i];
			jur = NUR[j];


			if (j > 0)
			{
				ip = (int)((j*j-j)/2 + 0.01);
				ip += j;
				ip += i;
			}
			else
			{
				ip = 0;
			}

			//if (iur <= jur) STIF[jur][jur-iur] += el->GE[i][j];
			//else STIF[iur][iur-jur] += el->GE[i][j];
			if ( (iur < NIE) && (jur < NIE) )
			{
				if (iur <= jur) 
				{
					STII[ STII_ENV[jur] - (jur-iur) ] += STSSse[ip];
				}
				else 
				{
					STII[ STII_ENV[iur] - (iur-jur) ] += STSSse[ip];	
				}
			}
		}
	}
}


void CSE::StifIINZAssemblingSE(double *STSSse, int KU, int *NUR)
{
	int ii,jj,i,j,in,jn,ik,jk,k,iur,jur,tmpi,ip;

	for(i=0; i<KU; i++)
	{
		for(j=i; j<KU; j++)
		{
			iur = NUR[i];// строка
			jur = NUR[j];// столбец

			if (jur < iur)
			{
				tmpi = jur;
				jur = iur;
				iur = tmpi;
			}

			if (j > 0)
			{
				ip = (int)((j*j-j)/2 + 0.01);
				ip += j;
				ip += i;
			}
			else
			{
				ip = 0;
			}

			if ( (iur < NIE) && (jur < NIE) )
			{
				//поиск позиции в соответствующем столбце матрицы жесткости
				for (k=0; k<STII_ENVNZ[jur][0]; k++)
				{
					if ( STII_ENVNZ[jur][k+1] == iur ) break;
				}

				if (k == STII_ENVNZ[jur][0])
				{
					// Error
					k = k;
				}
				else
				{
					STIINZ[jur][k] += STSSse[ip];
				}
			}
		}
	}

}


void CSE::StifISAssemblingFromSE(double *STSSse, int KU, int *NUR )
{
	int i,j,env,k,iur,jur,tmp,inode,jnode,nnode,nnodesum,iurenv,ip;
	
	for(i=0; i<KU; i++)
	{
		for(j=i; j<KU; j++)
		{
			iur = NUR[i];
			jur = NUR[j];


			if (j > 0)
			{
				ip = (int)((j*j-j)/2 + 0.01);
				ip += j;
				ip += i;
			}
			else
			{
				ip = 0;
			}


			if ( iur > jur ) //заполняем верхний треугольник
			{
				tmp = iur;
				iur = jur;
				jur = tmp;
			}

			if ( (jur >= NIE) && (iur < NIE) ) //проверка попадания в IS часть
			{
				inode = (int)(iur/KORT);
				jnode = (int)((jur-NIE)/KORT);

				nnodesum = 0;
				env = STIS_ENV[jnode][1];
				for (k=0; k<env; k+=2)
				{
					/*if ( k > 0 )
					{
						nnode = STIS_ENV[jnode][k+2+1] - STIS_ENV[jnode][k+1];
					}
					else
					{
						nnode = STIS_ENV[jnode][k+2+1];
					}*/
					nnode = STIS_ENV[jnode][k+2+1];
					if ( (inode >= STIS_ENV[jnode][k+2]) && (inode < (STIS_ENV[jnode][k+2] + nnode)) ) break; //найден нужный интервал ненулей в столбце
					nnodesum += nnode;
				}

				/*if ( k > 0 )
				{
					iurenv = STIS_ENV[jnode][k+1]*KORT;
				}
				else
				{
					iurenv = 0;
				}
				iurenv =  iurenv + (iur - STIS_ENV[jnode][k+2]*KORT);*/
				iurenv =  nnodesum*KORT + (iur - STIS_ENV[jnode][k+2]*KORT);

				STIS[jur-NIE][iurenv] += STSSse[ip];
			}
		}
	}
}

void CSE::StifSSAssemblingFromSE(double *STSSse, int KU, int *NUR )
{
	//заполняется полный треугольник
	int i,j,iur,jur,tmp,ip,ip2,iurl, jurl;
	
	for(i=0; i<KU; i++)
	{
		for(j=i; j<KU; j++)
		{
			iur = NUR[i];
			jur = NUR[j];

			if (j > 0) //вычисление номера элемента в профиле STSSse
			{
				ip2 = (int)((j*j-j)/2 + 0.01);
				ip2 += j;
				ip2 += i;
			}
			else
			{
				ip2 = 0;
			}

			if ( iur > jur ) //заполняем верхний треугольник
			{
				tmp = iur;
				iur = jur;
				jur = tmp;
			}
			iurl = iur - NIE;
			jurl = jur - NIE;

			if ( (iurl >= 0) && (jurl >= 0) )
			{
				if (jurl > 0)
				{
					ip = (int)((jurl*jurl-jurl)/2 + 0.01);
					ip += jurl;
					ip += iurl;
				}
				else
				{
					ip = 0;
				}
				STSS[ip] += STSSse[ip2]; 
			}
			
		}
	}
}