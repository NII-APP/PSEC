#include <stdio.h>
#include "stdafx.h"
#include "EL.h"

int CEL::DINT ()

{
int i;
double C,c1,c2,E,mu;	


if ( (eltype == 24) || (eltype == 25) )
{
	E = material->E;
	mu = material->MU;
	C = E*(1-mu)/((1+mu)*(1-2*mu));
	c1 = mu/(1-mu);
	c2 = (1-2*mu)/(2*(1-mu));
	D[0][0]=C;
	D[0][1]=C*c1;
	D[0][2]=C*c1;
	D[0][3]=0;
	D[0][4]=0;
	D[0][5]=0;
	D[1][0]=C*c1;
	D[1][1]=C;
	D[1][2]=C*c1;
	D[1][3]=0;
	D[1][4]=0;
	D[1][5]=0;
	D[2][0]=C*c1;
	D[2][1]=C*c1;
	D[2][2]=C;
	D[2][3]=0;
	D[2][4]=0;
	D[2][5]=0;
	D[3][0]=0;
	D[3][1]=0;
	D[3][2]=0;
	D[3][3]=C*c2;
	D[3][4]=0;
	D[3][5]=0;
	D[4][0]=0;
	D[4][1]=0;
	D[4][2]=0;
	D[4][3]=0;
	D[4][4]=C*c2;
	D[4][5]=0;
	D[5][0]=0;
	D[5][1]=0;
	D[5][2]=0;
	D[5][3]=0;
	D[5][4]=0;
	D[5][5]=C*c2;
}

if (eltype == 5) 
{
	E = material->E;
	mu = material->MU;
	C = E/(1-mu*mu);
	D[0][0]=C*1;
	D[0][1]=C*mu;
	D[0][2]=0;
	D[1][0]=C*mu;
	D[1][1]=C*1;
	D[1][2]=0;
	D[2][0]=0;
	D[2][1]=0;
	D[2][2]=C*(1-mu)/2;
	
}

	 return(0);

}
