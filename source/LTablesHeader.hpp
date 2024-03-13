#ifndef HEADERFILE_LTABLESHEADER
#define HEADERFILE_LTABLESHEADER

#include "Grids.hpp"

//defining global variables to be used in all source files:

extern string LT_pName;		//particle name
extern double LT_xB;        //xB value
extern GridPoints LT_Grids; //grids
extern int LdndxMaxPoints;  //maximal number of points for Ldndx integration
extern int LCollMaxPoints;  //maximal number of points for collisional integration

#endif