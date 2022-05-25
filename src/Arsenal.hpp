#ifndef HEADERFILE_ARSENALHEADER
#define HEADERFILE_ARSENALHEADER

#include "LinearInterpolation.hpp"

int GetHeaderLength(ifstream &myfile); //function that determines header length of file

double ProductLog(double x); 		 //ProductLog function
double UnitStep(double x); 			 //UnitStep for double
long double UnitStep(long double x); //UnitStep for long double

double HaltonSequence(int index, int base); //generates numbers of Halton sequence

double LinearIntegrate(vector<double> xdata, vector<double> fdata);							    //analytical integration of 1rd order polynomial in it's domain
double LinearIntegrate(vector<double> xdata, vector<double> fdata, double lim_l, double lim_h); //analytical integration of 1rd order polynomial in given range
double  CubicIntegrate(vector<double> xdata, vector<double> fdata);							    //analytical integration of 3rd order polynomial in it's domain
double  CubicIntegrate(vector<double> xdata, vector<double> fdata, double lim_l, double lim_h); //analytical integration of 3rd order polynomial in given range

int GenerateInitPosPoints(); //generates initial positions points and angles

void GaussFilterIntegrate(vector<double> radiativeRAA1, vector<vector<double>> radiativeRAA2, vector<double> collisionalEL, vector<double> &singRAA1, vector<vector<double>> &singRAA2);						  //function that performs Gauss filter integration - modefied pT integration algorithm
void GaussFilterIntegrate(interpFun dsdpti2lquark, vector<double> radiativeRAA1, vector<vector<double>> radiativeRAA2, vector<double> collisionalEL, vector<double> &singRAA1, vector<vector<double>> &singRAA2); //function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
void GaussFilterIntegrate(vector<double> radiativeRAA, vector<double> collisionalEL, vector<double> &singRAA); 																	         						  //function that performs Gauss filter integration - default algorithm

void CalcObservables(vector<vector<double>> RAApTphi, vector<double> &RAApT, vector<double> &v2pT); 	   //function that calculates observables RAA(pT) and v2(pT)
void CalcAvgPLT(double *pl_dist_a, double *temp_dist_a, vector<double> &avg_pl, vector<double> &avg_temp); //function that calculates path-lengths and temperatures

#endif