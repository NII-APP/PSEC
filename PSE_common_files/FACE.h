#pragma once

#include "StructData.h"
#include "MEM.h"



class FACE
{
public:
	FACE(void);
public:
	~FACE(void);
	
public:
	int facenodenums[8]; //������������ ������� (��������)
	int NNface; //����� ����� �����
	double normal[3];
	int iel;
	bool isouter;
};