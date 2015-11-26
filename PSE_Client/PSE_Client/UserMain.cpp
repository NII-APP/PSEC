#include "stdafx.h"
#include "StructData.h"
#include "SE_structure.h"
#include "Problem.h"
#include "stdio.h"
#include "conio.h"
#include "Windows.h"

extern int usermain_motionint_eigf();
extern int usermain_eigf();
extern int usermain_elastat();

extern char glob_str_path[256];
char glob_str_taskpath[256];

void usermain()
{
	sprintf(glob_str_path,"C:\\PSEtemp\\testeigen\\stifmatr1");

	sprintf(glob_str_taskpath,"C:\\Temp\\TestNX\\AnsFormat1");

	//usermain_motionint_eigf();
	usermain_eigf();
	//usermain_elastat();

}