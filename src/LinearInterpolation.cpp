#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
using namespace std;

#include "LinearInterpolation.hpp"

//CONSTRUCTORS:
interpFun::interpFun() {

}

//input is 2 1D arrays:
interpFun::interpFun(double *XData, double *FData, int NofElements)
{
	SetData(XData, FData, NofElements);
}

void interpFun::SetData(double *XData, double *FData, int NofElements)
{
	dim = 1;
	dataLength = NofElements;

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = vector<double> (XData, XData + dataLength);
    data[1] = vector<double> (FData, FData + dataLength);

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 2 1D vectors:
interpFun::interpFun(vector<double> XData, vector<double> FData)
{
	SetData(XData, FData);
}

void interpFun::SetData(vector<double> XData, vector<double> FData)
{
	dim = 1;
	dataLength = FData.size();

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = XData;
    data[1] = FData;

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 3 1D arrays:
interpFun::interpFun(double *X1Data, double *X2Data, double *FData, int NofElements)
{
	SetData(X1Data, X2Data, FData, NofElements);
}

void interpFun::SetData(double *X1Data, double *X2Data, double *FData, int NofElements)
{
	dim = 2;
	dataLength = NofElements;

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = vector<double> (X1Data, X1Data + dataLength);
	data[1] = vector<double> (X2Data, X2Data + dataLength);
    data[2] = vector<double> (FData, FData + dataLength);

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 3 1D vectors:
interpFun::interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> FData)
{
	SetData(X1Data, X2Data, FData);
}

void interpFun::SetData(vector<double> X1Data, vector<double> X2Data, vector<double> FData)
{
	dim = 2;
	dataLength = FData.size();

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = X1Data;
	data[1] = X2Data;
    data[2] = FData;

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 2 1D vectors (grids) and 1 2d vector (function values):
interpFun::interpFun(vector<double> X1Data, vector<double> X2Data, vector<vector<double>> FData)
{
	SetData(X1Data, X2Data, FData);
}

void interpFun::SetData(vector<double> X1Data0, vector<double> X2Data0, vector<vector<double>> FData0)
{
	vector<double> X1Data, X2Data, FData;

	for (int i=0; i<X1Data0.size(); i++)
	{
		for (int j=0; j<X2Data0.size(); j++)
		{
			X1Data.push_back(X1Data0[i]);
			X2Data.push_back(X2Data0[j]);
			 FData.push_back(FData0[i][j]);
		}
	}

	dim = 2;
	dataLength = FData.size();

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = X1Data;
	data[1] = X2Data;
    data[2] = FData;

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 4 1D arrays:
interpFun::interpFun(double *X1Data, double *X2Data, double *X3Data, double *FData, int NofElements)
{
	SetData(X1Data, X2Data, X3Data, FData, NofElements);
}

void interpFun::SetData(double *X1Data, double *X2Data, double *X3Data, double *FData, int NofElements)
{
	dim = 3;
	dataLength = NofElements;

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = vector<double> (X1Data, X1Data + dataLength);
	data[1] = vector<double> (X2Data, X2Data + dataLength);
	data[2] = vector<double> (X3Data, X3Data + dataLength);
    data[3] = vector<double> (FData, FData + dataLength);

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 4 1D vectors:
interpFun::interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> FData)
{
	SetData(X1Data, X2Data, X3Data, FData);
}

void interpFun::SetData(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> FData)
{
	dim = 3;
	dataLength = FData.size();

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = X1Data;
	data[1] = X2Data;
	data[2] = X3Data;
    data[3] = FData;

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 5 1D arrays:
interpFun::interpFun(double *X1Data, double *X2Data, double *X3Data, double *X4Data, double *FData, int NofElements)
{
	SetData(X1Data, X2Data, X3Data, X4Data, FData, NofElements);
}

void interpFun::SetData(double *X1Data, double *X2Data, double *X3Data, double *X4Data, double *FData, int NofElements)
{
	dim = 4;
	dataLength = NofElements;

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = vector<double> (X1Data, X1Data + dataLength);
	data[1] = vector<double> (X2Data, X2Data + dataLength);
	data[2] = vector<double> (X3Data, X3Data + dataLength);
	data[3] = vector<double> (X4Data, X4Data + dataLength);
    data[4] = vector<double> (FData, FData + dataLength);

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//input is 5 1D vectors:
interpFun::interpFun(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> X4Data, vector<double> FData)
{
	SetData(X1Data, X2Data, X3Data, X4Data, FData);
}

void interpFun::SetData(vector<double> X1Data, vector<double> X2Data, vector<double> X3Data, vector<double> X4Data, vector<double> FData)
{
	dim = 4;
	dataLength = FData.size();

	//resizing vector for dim:
	data.resize(dim+1);

	//setting values of coordinates and function values:
	data[0] = X1Data;
	data[1] = X2Data;
	data[2] = X3Data;
	data[3] = X4Data;
    data[4] = FData;

	//determining lengts of grids and minimum and maximum values of data
	gridLengthsF();
	MaxMinDataF();

	//error message (if there is not enough data for interpolation):
	bool err = false;
	for (int i = 0; i<dim; i++) { if (gridLengths[i] < 2) err = true; }
	if (err) printf("Error: Interpolation not possible (not enough data).");
}

//DESTRUCTORS:
interpFun::~interpFun() {

}

//INTERPOLATION FUNCTIONS:
//1D interpolation
double interpFun::interp(double pointValue)
{
	if (dim >1) {
		printf("Error: Not enough points for interpolation.\n");
		return 0;
	}
	else {
		return interpolation1D(pointValue);
	}
}

//2D interpolation
double interpFun::interp(double pointValue1, double pointValue2)
{
	if (dim < 2) {
		printf("Error: Too much points for interpolation.\n");
		return 0;
	}
	else if (dim >2) {
		printf("Error: Not enough points for interpolation.\n");
		return 0;
	}
	else {
		return interpolation2D(pointValue1, pointValue2);
	}
}

//3D interpolation
double interpFun::interp(double pointValue1, double pointValue2, double pointValue3)
{
	if (dim < 3) {
		printf("Error: Too much points for interpolation.\n");
		return 0;
	}
	else if (dim >3) {
		printf("Error: Not enough points for interpolation.\n");
	}
	else {
		return interpolation3D(pointValue1, pointValue2, pointValue3);
	}
	return 0.0;
}

//4D interpolation
double interpFun::interp(double pointValue1, double pointValue2, double pointValue3, double pointValue4)
{
	if (dim < 4) {
		printf("Error: Too much points for interpolation.\n");
		return 0;
	}
	else if (dim >4) {
		printf("Error: Not enough points for interpolation.\n");
	}
	else {
		return interpolation4D(pointValue1, pointValue2, pointValue3, pointValue4);
	}
	return 0.0;
}

//DOMAIN FUNCTIONS:
vector<vector<double>> interpFun::domain()
{
	vector<vector<double>> domain_vector;

	domain_vector.resize(dim);
	for (int i=0; i<dim; i++)
	{
		domain_vector[i].push_back(minValues[i]); domain_vector[i].push_back(maxValues[i]);
	}

	return domain_vector;
}

vector<double> interpFun::codomain()
{
	vector<double> codomain_vector;
	codomain_vector.push_back(*min_element(begin(data[dim]), end(data[dim])));
	codomain_vector.push_back(*max_element(begin(data[dim]), end(data[dim])));
	return codomain_vector;
}

//function tha determines axes lengths
void interpFun::gridLengthsF()
{
	//determing grid lengths
	int cnt;
	bool pomocnaBool;
	double pomData;
	for (int i = 0; i<dim; i++) { gridLengths[i] = 1; }
	int gLM = 1;

	for (int i = dim - 1; i>0; i--)
	{
		pomocnaBool = true;
		cnt = 0;
		pomData = data[i][0];

		while (pomocnaBool)
		{
			cnt++;
			if (pomData == data[i][gLM*cnt]) {
				pomocnaBool = false;
			}
		}
		gridLengths[i] = cnt;
		gLM *= gridLengths[i];
	}
	gridLengths[0] = dataLength;
	for (int i = dim - 1; i>0; i--) { gridLengths[0] /= gridLengths[i]; }

	//determining relative postions (first different value of a point on the grid)
	for (int i = 0; i<dim; i++)
	{
		relPosition[i] = 1;
		for (int j = dim - 1; j>i; j--)
		{
			relPosition[i] *= gridLengths[j];
		}
	}
}

//function that determinens maximum and minumum values of data
void interpFun::MaxMinDataF()
{
	double mmVal;
	//determining maximum values of data
	for (int i = 0; i<dim; i++)
	{
		mmVal = data[i][0];
		for (int j = 1; j<gridLengths[i]; j++)
		{
			if (mmVal < data[i][j*relPosition[i]]) {
				mmVal = data[i][j*relPosition[i]];
			}
		}
		maxValues[i] = mmVal;
	}
	//determining minimum values of data
	for (int i = 0; i<dim; i++)
	{
		mmVal = data[i][0];
		for (int j = 1; j<gridLengths[i]; j++)
		{
			if (mmVal > data[i][j*relPosition[i]]) {
				mmVal = data[i][j*relPosition[i]];
			}
		}
		minValues[i] = mmVal;
	}
}

//function that locates points
void interpFun::locatePointF(double *points, int *position)
{
	//checking if values lie within the range
	for (int i = 0; i<dim; i++)
	{
		if ((points[i] > maxValues[i]) || points[i] < minValues[i]) {
			printf("ERROR: Value lies outside the range.\n");
			for (int i = 0; i<dim; i++) position[i] = -1;
			return;
		}
	}

	//finding the position of points
	int ju, jm, jl, mm, n;
	int mmmin; //helping variable
	bool ascnd;
	for (int i = 0; i<dim; i++)
	{
		jl = 0;
		ju = gridLengths[i] - 1;
		mm = 2;
		ascnd = (data[i][dataLength - 1] >= data[i][0]);

		while ((ju - jl) >1) {
			jm = (ju + jl) >> 1;
			if ((points[i] >= data[i][jm*relPosition[i]]) == ascnd) {
				jl = jm;
			}
			else {
				ju = jm;
			}
		}
		n = gridLengths[i];
		mmmin = n - mm < jl - ((mm - 2) >> 1) ? n - mm : jl - ((mm - 2) >> 1);
		position[i] = 0 > mmmin ? 0 : mmmin;
	}
}

//1D linear interpolation
double interpFun::lin1DInterpolation(double *x, double *y, double xx)
{
	//x - points values
	//y - function values
	//xx - point in which the value will be determined

	double res = y[0] + (xx - x[0])*(y[1] - y[0]) / (x[1] - x[0]);

	return res;
}

//1D interpolation (full function)
double interpFun::interpolation1D(double pt)
{
	//searching for position
	double *pts = new double[dim];
	pts[0] = pt;
	int *pos = new int[dim];
	locatePointF(pts, pos);

	//setting x and Q values
	double xPts[2];
	double Q[2];
	for (int i = 0; i<2; i++)
	{
		xPts[i] = data[0][pos[0] + i];
		Q[i] = data[1][pos[0] + i];
		//cout<<xPts[i]<<" "<<Q[i]<<"\n";
	}

	//interpolation
	double res = lin1DInterpolation(xPts, Q, pts[0]);

	//deleting allocated memory for pts, pos
	delete[] pts; delete[] pos;

	//returning interpolated value
	return res;
}

//2D interpolation
double interpFun::interpolation2D(double pt1, double pt2)
{
	//searching for position
	double pts[2];
	pts[0] = pt1;
	pts[1] = pt2;
	int pos[2];
	locatePointF(pts, pos);

	//setting x1 values
	double x1Pts[2];
	for (int i = 0; i<2; i++) { x1Pts[i] = data[0][(pos[0] + i)*relPosition[0]]; }

	//setting x2 values
	double x2Pts[2];
	for (int i = 0; i<2; i++) { x2Pts[i] = data[1][(pos[1] + i)*relPosition[1]]; }

	double Q1[2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			Q1[i][j] = data[2][(pos[0] + i)*relPosition[0] + (pos[1] + j)*relPosition[1]];
		}
	}

	double Q2[2];
	for (int i = 0; i<2; i++)
	{
		Q2[i] = lin1DInterpolation(x2Pts, Q1[i], pts[1]);
	}

	//obtaining final result by interpolating in x1 axis
	double res = lin1DInterpolation(x1Pts, Q2, pts[0]);

	//returning interpolated value
	return res;
}
//3D interpolation
double interpFun::interpolation3D(double pt1, double pt2, double pt3)
{
	//searching for position
	double pts[3];
	pts[0] = pt1;
	pts[1] = pt2;
	pts[2] = pt3;
	int pos[3];
	locatePointF(pts, pos);

	//setting x1 values
	double x1Pts[2];
	for (int i = 0; i<2; i++) { x1Pts[i] = data[0][(pos[0] + i)*relPosition[0]]; }

	//setting x2 values
	double x2Pts[2];
	for (int i = 0; i<2; i++) { x2Pts[i] = data[1][(pos[1] + i)*relPosition[1]]; }

	//setting x3 values
	double x3Pts[2];
	for (int i = 0; i<2; i++) { x3Pts[i] = data[2][(pos[2] + i)*relPosition[2]]; }

	double Q1[2][2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			for (int k = 0; k<2; k++)
			{
				Q1[i][j][k] = data[3][(pos[0] + i)*relPosition[0] + (pos[1] + j)*relPosition[1] + (pos[2] + k)*relPosition[2]];
			}
		}
	}

	double Q2[2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			Q2[i][j] = lin1DInterpolation(x3Pts, Q1[i][j], pts[2]);
		}
	}

	double Q3[2];
	for (int i = 0; i<2; i++)
	{
		Q3[i] = lin1DInterpolation(x2Pts, Q2[i], pts[1]);
	}

	//obtaining final result by interpolating in x1 axis
	double res = lin1DInterpolation(x1Pts, Q3, pts[0]);

	//returning interpolated value
	return res;
}

//4D interpolation
double interpFun::interpolation4D(double pt1, double pt2, double pt3, double pt4)
{
	//searching for position
	double pts[4];
	pts[0] = pt1;
	pts[1] = pt2;
	pts[2] = pt3;
	pts[3] = pt4;
	int pos[4];
	locatePointF(pts, pos);

	//setting x1 values
	double x1Pts[2];
	for (int i = 0; i<2; i++) { x1Pts[i] = data[0][(pos[0] + i)*relPosition[0]]; }

	//setting x2 values
	double x2Pts[2];
	for (int i = 0; i<2; i++) { x2Pts[i] = data[1][(pos[1] + i)*relPosition[1]]; }

	//setting x3 values
	double x3Pts[2];
	for (int i = 0; i<2; i++) { x3Pts[i] = data[2][(pos[2] + i)*relPosition[2]]; }

	//setting x4 values
	double x4Pts[2];
	for (int i = 0; i<2; i++) { x4Pts[i] = data[3][(pos[3] + i)*relPosition[3]]; }

	double Q1[2][2][2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			for (int k = 0; k<2; k++)
			{
				for (int l = 0; l<2; l++)
				{
					Q1[i][j][k][l] = data[4][(pos[0] + i)*relPosition[0] + (pos[1] + j)*relPosition[1] + (pos[2] + k)*relPosition[2] + (pos[3] + l)*relPosition[3]];
				}
			}
		}
	}

	double Q2[2][2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			for (int k = 0; k<2; k++)
			{
				Q2[i][j][k] = lin1DInterpolation(x4Pts, Q1[i][j][k], pts[3]);
			}
		}
	}
	
	double Q3[2][2];
	for (int i = 0; i<2; i++)
	{
		for (int j = 0; j<2; j++)
		{
			Q3[i][j] = lin1DInterpolation(x3Pts, Q2[i][j], pts[2]);
		}
	}

	double Q4[2];
	for (int i = 0; i<2; i++)
	{
		Q4[i] = lin1DInterpolation(x2Pts, Q3[i], pts[1]);
	}

	//obtaining final result by interpolating in x1 axis
	double res = lin1DInterpolation(x1Pts, Q4, pts[0]);

	//returning interpolated value
	return res;
}