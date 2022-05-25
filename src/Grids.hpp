#ifndef HEADERFILE_GRIDPOINTS
#define HEADERFILE_GRIDPOINTS

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

class GridPoints {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//public functions:
public:

	//CONSTRUCTORS:
	GridPoints();
	GridPoints(string particle_name);
	void SetGridPoints(string particle_name);

	//DESTRUCTOR:
	~GridPoints();

	//GRID FUNCTIONS:
	vector<double> tauPts();
	double tauPts(int i);
	int tauPtsLength();

	vector<double> pPts();
	double pPts(int i);
	int pPtsLength();

	vector<double> xPts();
	double xPts(int i);
	int xPtsLength();

	vector<double> TPts();
	double TPts(int i);
	int TPtsLength();

	vector<double> FdpPts();
	double FdpPts(int i);
	int FdpPtsLength();

	vector<double> RadPts();
	double RadPts(int i);
	int RadPtsLength();

	vector<double> pCollPts();
	double pCollPts(int i);
	int pCollPtsLength();

	vector<double> TCollPts();
	double TCollPts(int i);
	int TCollPtsLength();

	vector<double> finPts();
	double finPts(int i);
	int finPtsLength();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private variables and functions:
private:

	vector<double> grid_tauPts, grid_pPts, grid_TPts, grid_xPts, grid_RadPts, grid_FdpPts, grid_pCollPts, grid_TCollPts, grid_finPts; //vectors that store grids

	double LinearIntegrate(vector<double> data_x, vector<double> data_f, double xH); //linear integration function
	vector<double> GenerateGrids(vector<vector<double>> density, int numpts);		 //function that generates grid points
};

#endif