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
	int facenodenums[8]; //квадратичный элемент (максимум)
	int NNface; //число узлов грани
	double normal[3];
	int iel;
	bool isouter;
};