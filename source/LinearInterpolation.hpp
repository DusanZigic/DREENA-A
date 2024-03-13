#ifndef HEADERFILE_LINEARINTERPOLATION
#define HEADERFILE_LINEARINTERPOLATION

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <vector>
using namespace std;

class interpFun {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//public functions:
public:

	//CONSTRUCTORS:

	interpFun();

	//input is 2 1D arrays:
	interpFun(double *XData, double *FData, int NofElements);
	void SetData(double *XData, double *FData, int NofElements);

	//input is 2 1D vectors:
	interpFun(vector<double> XData, vector<double> FData);
	void SetData(vector<double> XData, vector<double> FData);

	//input is 3 1D arrays:
	interpFun(double *X1Data, double *X2Data, double *FData, int NofElements);
	void SetData(double *X1Data, double *X2Data, double *FData, int NofElements);

	//input is 3 1D vectors:
	interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> FData);
	void SetData(vector<double> X1Data, vector<double> X2Data, vector<double> FData);

	//input is 2 1D vectors (grids) and 1 2d vector (function values):
	interpFun(vector<double> X1Data, vector<double> X2Data, vector<vector<double>> FData);
	void SetData(vector<double> X1Data0, vector<double> X2Data0, vector<vector<double>> FData0);

	//input is 4 1D arrays:
	interpFun(double *X1Data, double *X2Data, double *X3Data, double *FData, int NofElements);
	void SetData(double *X1Data, double *X2Data, double *X3Data, double *FData, int NofElements);

	//input is 4 1D vectors:
	interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> FData);
	void SetData(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> FData);

	//input is 5 1D arrays:
	interpFun(double *X1Data, double *X2Data, double *X3Data, double *X4Data, double *FData, int NofElements);
	void SetData(double *X1Data, double *X2Data, double *X3Data, double *X4Data, double *FData, int NofElements);

	//input is 5 1D vectors:
	interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> X4Data, vector<double> FData);
	void SetData(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> X4Data, vector<double> FData);

	//DESTRUCTOR:
	~interpFun();

	//INTERPOLATION FUNCTIONS:
	//1D interpolation
	double interp(double pointValue);

	//2D interpolation
	double interp(double pointValue1, double pointValue2);

	//3D interpolation
	double interp(double pointValue1, double pointValue2, double pointValue3);

	//4D interpolation
	double interp(double pointValue1, double pointValue2, double pointValue3, double pointValue4);

	//miscellaneous FUNCTIONS:

	//function that returns domains:
	vector<vector<double>> domain();

	//function that returns codomain:
	vector<double> codomain();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private variables and functions:
private:

	int dataLength; //total length of the data
	vector<vector<double> > data; //2D vector that stores all the data (dimension: (No variables+1) x dataLength)
	int dim; //No of variables
	int gridLengths[4]; //1D array that stores lengths of specific axes
	int relPosition[4]; //1D array that stores relative positions of specific axes
	double maxValues[4], minValues[4]; //1D arrays that store maximum and minimum values od variables data

	//function tha determines axes lengths
	void gridLengthsF();

	//function that determinens maximum and minumum values of data
	void MaxMinDataF();

	//function that locates points
	void locatePointF(double *points, int *position);

	//1D linear interpolation
	double lin1DInterpolation(double *x, double *y, double xx);

	//1D interpolation (full function)
	double interpolation1D(double pt);

	//2D interpolation
	double interpolation2D(double pt1, double pt2);

	//3D interpolation
	double interpolation3D(double pt1, double pt2, double pt3);

	//4D interpolation
	double interpolation4D(double pt1, double pt2, double pt3, double pt4);
};

#endif