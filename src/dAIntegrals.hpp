#ifndef HEADERFILE_DAHEADER
#define HEADERFILE_DAHEADER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FdA integrals definitions:
void FdAHaltonSeqInit(int FdAMaxPts);										 //generates Halton sequences for FdA integrals
double dAp410(double ph, interpFun &normint);								 //dAp410 integral definition
double FdA411(double ph, double dp, interpFun &normint, interpFun &dndxint); //FdA411 integral definition
double FdA412(double ph, double dp, interpFun &normint, interpFun &dndxint); //FdA412 integral definition
double FdA413(double ph, double dp, interpFun &normint, interpFun &dndxint); //FdA413 integral definition
double FdA414(double ph, double dp, interpFun &normint, interpFun &dndxint); //FdA414 integral definition
double FdA415(double ph, double dp, interpFun &normint, interpFun &dndxint); //FdA415 integral definition

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//dA integrals definitions:
void dAHaltonSeqInit(int dAMaxPts); 							 //generates Halton sequences for dA integrals
double dA410(double ph, interpFun &normint); 					 //dA410 integral definition
double dA411(double ph, interpFun &normint, interpFun &dndxint); //dA411 integral definition
double dA412(double ph, interpFun &normint, interpFun &dndxint); //dA412 integral definition
double dA413(double ph, interpFun &normint, interpFun &dndxint); //dA413 integral definition
double dA414(double ph, interpFun &normint, interpFun &dndxint); //dA414 integral definition
double dA415(double ph, interpFun &normint, interpFun &dndxint); //dA415 integral definition
double dA416(double ph, interpFun &normint, interpFun &dndxint); //dA416 integral definition
double dA417(double ph, interpFun &normint, interpFun &dndxint); //dA417 integral definition

#endif