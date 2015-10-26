#include "StdAfx.h"
#include "EL.h"
#include "math.h"

int CEL::GetLocCrd(double *globcrd, double *extloccrd, double *FF)
{
	// ����������: -2  -  ����� ��� �������� (�� ���. ���.)
	// -1  -  ��������� ���������� ����� ��������
	// 2...14 - ����� �������� � �������� �����������

	int i,j,it,fl;
	double **J,detJ;
	double **dFloc; //����������� ������ ����� �� ��������� �����������
	double *loccrd,*loccrd_old, *globcrd_it;
	double err,eps;
	
	loccrd = new double[nvlcrd];
	loccrd_old = new double[nvlcrd];
	globcrd_it = new double[NORT];
	J = NULL;
	J = MM->MEM_NEW(J,NORT,NORT);
	dFloc = NULL;
	dFloc = MM->MEM_NEW(dFloc,NORT,NN);

	for (i=0; i<NORT; i++) loccrd[i] = extloccrd[i];

	err = 1.0;
	eps = 1.e-6;
	it = 0;
	while ( err > eps && it < 15 )
	{
		for (i=0; i<NORT; i++) loccrd_old[i] = loccrd[i];
		//������ ������� ����� � �����-����������� � �������� ���������� ������������
		FF_all(FF,loccrd);
		//������ ����������� ������� ����� �� ��������� �����������
		dFF_all(dFloc,loccrd);
		//������ ���������� ��������� �����-�����������
		for (i=0; i<NORT; i++) globcrd_it[i] = 0.0;
		for (i=0; i<NORT; i++)
		{
			for (j=0; j<NN; j++) globcrd_it[i] += fullcrd[ind[j]*NORTfullcrd+i]*FF[j];
		}
		//������ ������� ���������� ��������� �����-�����������
		for (i=0; i<NORT; i++) globcrd_it[i] = globcrd[i] - globcrd_it[i];
		//������ ������� ����� �������������� ���������
		Jacoby(dFloc, J);
		DetJacoby(J,&detJ);
		if (detJ == 0.0) //��������� ������������ ������ ����������� ����������� Z, ������� ����� ���������, ����������� ��� 
		{
			break;
		}
		InvJacoby(J,detJ);
		//�������������� ������� ���������� ��������� � ���������� ��������� ���������
		//������ ����� ��������� ���������
		for (i=0; i<NORT; i++)
		{
			for (j=0; j<NORT; j++) loccrd[i] += J[j][i]*globcrd_it[j];
		}
		//������ ������� � ��������� ����������� �� ��������
		//���������� �� ��������� �.�. ���. ��������� � ��� �����������
		err = 0.0;
		for (i=0; i<NORT; i++)
		{
			if ( fabs( loccrd[i] - loccrd_old[i] ) > err ) err = fabs( loccrd[i] - loccrd_old[i] );
		}
		it++;
	}
	if ( it > 0 )
	{
		for (i=0; i<NORT; i++) extloccrd[i] = loccrd[i];
		FF_all(FF,loccrd);
		//�������� ��������� ����� ������ ��������
		fl = 0;
		for (i=0; i<nvlcrd; i++)
		{
			if ( loccrd[i] < lcrdlim[0] || loccrd[i] > lcrdlim[1] )
			{
				fl++;
			}
		}
		if ( fl > 0 ) it = -2; 
		if ( it == 15 ) it = -1;
	}
	dFloc = MM->MEM_DEL(dFloc,NORT,NN);
	J = MM->MEM_DEL(J,NORT,NORT);
	delete [] loccrd;
	delete [] loccrd_old;
	delete [] globcrd_it;

	return( it );
}