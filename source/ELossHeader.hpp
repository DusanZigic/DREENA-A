#ifndef HEADERFILE_ELOSSHEADER
#define HEADERFILE_ELOSSHEADER

#include "LinearInterpolation.hpp"
#include "Grids.hpp"

//defining global variables to be used in all source files:

extern string pName; 					    		  //particle name
extern string centrality;				    		  //centrality class
extern double xB;						    		  //xB value
extern GridPoints Grids;				    		  //grid points
extern interpFun dsdpti2; 				    		  //initial pT distribution
extern interpFun LNorm, Ldndx, LColl;	    		  //interpolated L tables
extern interpFun TProfile, BinCollDensity;			  //temperature evolution and binary collision density
extern double tau0;									  //thermalization time
extern double mgC, MC;					   			  //constant particle and gluon masses used for dA integrals
extern double TCollConst;				    		  //constant temperature used for Gauss filter integration
extern int xGridN, yGridN, phiGridN;	    		  //initial position grid points and angle number
extern vector<double> xGridPts, yGridPts, phiGridPts; //defining vectors that store initial position points and angles
extern double TIMESTEP, TCRIT; 					   	  //time step and critical temperature
extern string pTinit_path, temp_path, bcd_path;	   	  //initial pT distribution, temperature evolution and binary collision density paths

#endif