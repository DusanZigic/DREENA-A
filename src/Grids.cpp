#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

#include "Grids.hpp"
#include "Arsenal.hpp"
#include "LinearInterpolation.hpp"

//CONSTRUCTORS:

GridPoints::GridPoints() {

}

GridPoints::GridPoints(string particle_name)
{
	SetGridPoints(particle_name);
}

void GridPoints::SetGridPoints(string particle_name)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (particle_name == "Bottom") {

		//tauPts:
		int taugridn = 21;
		vector<vector<double>> tauden{{0.0, 10}, {20, 10}};
		grid_tauPts = GenerateGrids(tauden, taugridn);

		//pPts:
		int pgridn = 25;
		vector<vector<double>> pden{{1, 8}, {20, 7}, {30, 3}, {60, 5}, {200, 1}};
		grid_pPts = GenerateGrids(pden, pgridn);

		//TPts:
		int Tgridn = 40;
		vector<vector<double>> Tden{{0.01, 10}, {2.000, 10}};
		grid_TPts = GenerateGrids(Tden, Tgridn);

		//xPts:
		double T = grid_TPts[0], nf = 3.0, lambda = 0.2;
		double mu = 0.197*sqrt((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda / ProductLog((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda));
		double mg = mu/sqrt(2.0);
		double M = 4.75;
		double MAXP = 200.0;
		int xgridn = 30;
		double xmin = mg/(MAXP + sqrt(MAXP*MAXP + M*M));
		for (int i=0; i<xgridn; i++)
			grid_xPts.push_back(exp(log(xmin) - log(xmin)/(xgridn-1.0)*i));

		//RadPts:
		int Radgridn = 20;
		vector<vector<double>> Radden{{2, 10}, {21.8, 10.0}, {44.5, 1.050001}, {170, 1}};
		grid_RadPts = GenerateGrids(Radden, Radgridn);

		//FdpPts:
		double mgC = 0.5140824927;
		int Fdpgridn = 16;
		vector<vector<double>> Fdpden = {{5.0*mgC/2.0, 10}, {12.0, 5}, {30, 0.}};
		grid_FdpPts = GenerateGrids(Fdpden, Fdpgridn-4);
		grid_FdpPts.insert(grid_FdpPts.begin(), 4.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 3.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 2.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 1.0*mgC/2.0);

		//pCollPts:
		int pCollgridn = 20;
		vector<vector<double>> pCollden{{1, 10}, {4, 10}, {9, 2.5}, {30, 0.6}, {60, 0.5}, {170, 0.3}};
		grid_pCollPts = GenerateGrids(pCollden, pCollgridn);

		//TCollPts:
		int TCollgridn = 40;
		vector<vector<double>> TCollden{{0.01, 10}, {2.000, 10}};
		grid_TCollPts = GenerateGrids(TCollden, TCollgridn);

		//finpts:
		int fingridn = 30;
		vector<vector<double>> finden{{5, 10}, {50, 10}, {70, 5}, {150, 3}};
		grid_finPts = GenerateGrids(finden, fingridn);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (particle_name == "Charm") {

		//tauPts:
		int taugridn = 21;
		vector<vector<double>> tauden{{0.0, 10}, {5.0, 10}, {10.0, 10}, {15.0, 10}, {20, 10}};
		grid_tauPts = GenerateGrids(tauden, taugridn);
		
		//pPts:
		int pgridn = 25;
		vector<vector<double>> pden{{1, 8}, {20, 7}, {30, 3}, {60, 5}, {200, 1}};
		grid_pPts = GenerateGrids(pden, pgridn);

		//TPts:
		int Tgridn = 40;
		vector<vector<double>> Tden{{0.01, 10}, {2.000, 10}};
		grid_TPts = GenerateGrids(Tden, Tgridn);

		//xPts:
		double T = grid_TPts[0], nf = 3.0, lambda = 0.2;
		double mu = 0.197*sqrt((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda / ProductLog((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda));
		double mg = mu/sqrt(2.0);
		double M = 1.2;
		double MAXP = 200.0;
		int xgridn = 30;
		double xmin = mg/(MAXP + sqrt(MAXP*MAXP + M*M));
		for (int i=0; i<xgridn; i++)
			grid_xPts.push_back(exp(log(xmin) - log(xmin)/(xgridn-1.0)*i));

		//RadPts:
		int Radgridn = 20;
		vector<vector<double>> Radden{{2, 10}, {21.8, 10.0}, {44.5, 1.050001}, {170, 1}};
		grid_RadPts = GenerateGrids(Radden, Radgridn);

		//FdpPts:
		double mgC = 0.5140824927;
		int Fdpgridn = 16;
		vector<vector<double>> Fdpden = {{5.0*mgC/2.0, 10}, {12.0, 5}, {30, 0.}};
		grid_FdpPts = GenerateGrids(Fdpden, Fdpgridn-4);
		grid_FdpPts.insert(grid_FdpPts.begin(), 4.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 3.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 2.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 1.0*mgC/2.0);

		//pCollPts:
		int pCollgridn = 20;
		vector<vector<double>> pCollden{{1, 10}, {4, 10}, {9, 2.5}, {30, 0.6}, {60, 0.5}, {170, 0.3}};
		grid_pCollPts = GenerateGrids(pCollden, pCollgridn);

		//TCollPts:
		int TCollgridn = 40;
		vector<vector<double>> TCollden{{0.01, 10}, {2.000, 10}};
		grid_TCollPts = GenerateGrids(TCollden, TCollgridn);

		//finpts:
		int fingridn = 30;
		vector<vector<double>> finden{{5, 10}, {50, 10}, {70, 5}, {150, 3}};
		grid_finPts = GenerateGrids(finden, fingridn);
		
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else if (particle_name == "Gluon") {

		//tauPts:
		int taugridn = 21;
		vector<vector<double>> tauden{{0.0, 10}, {20, 10}};
		grid_tauPts = GenerateGrids(tauden, taugridn);

		//pPts:
		int pgridn = 50;
		vector<vector<double>> pden{{1, 8}, {20, 7}, {40, 3}, {100, 5}, {450, 1}};
		grid_pPts = GenerateGrids(pden, pgridn);

		//TPts:
		int Tgridn = 40;
		vector<vector<double>> Tden{{0.01, 10}, {2.000, 10}};
		grid_TPts = GenerateGrids(Tden, Tgridn);

		//xPts:
		double T = grid_TPts[0], nf = 3.0, lambda = 0.2;
		double mu = 0.197*sqrt((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda / ProductLog((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda));
		double mg = mu/sqrt(2.0);
		double M = mu/sqrt(2.0);
		double MAXP = 450.0;
		int xgridn = 50;
		double xmin = mg/(MAXP + sqrt(MAXP*MAXP + M*M));
		for (int i=0; i<xgridn; i++)
			grid_xPts.push_back(exp(log(xmin) - log(xmin)/(xgridn-1.0)*i));

		//RadPts:
		int Radgridn = 40;
		vector<vector<double>> Radden{{2, 10}, {50, 10}, {70, 1}, {420, 1}};
		grid_RadPts = GenerateGrids(Radden, Radgridn);

		//FdpPts:
		double mgC = 0.5140824927;
		int Fdpgridn = 22;
		vector<vector<double>> Fdpden = {{5.0*mgC/2.0, 10}, {12.0, 5}, {30, 0.}};
		grid_FdpPts = GenerateGrids(Fdpden, Fdpgridn-4);
		grid_FdpPts.insert(grid_FdpPts.begin(), 4.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 3.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 2.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 1.0*mgC/2.0);

		//pCollPts:
		int pCollgridn = 40;
		vector<vector<double>> pCollden{{1, 10}, {4, 10}, {9, 2.5}, {30, 0.6}, {60, 0.5}, {420, 0.3}};
		grid_pCollPts = GenerateGrids(pCollden, pCollgridn);

		//TCollPts:
		int TCollgridn = 40;
		vector<vector<double>> TCollden{{0.01, 10}, {2.000, 10}};
		grid_TCollPts = GenerateGrids(TCollden, TCollgridn);

		//finpts:
		int fingridn = 50;
		vector<vector<double>> finden{{5, 10}, {50, 10}, {70, 5}, {400, 3}};
		grid_finPts = GenerateGrids(finden, fingridn);
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	else {

		//tauPts:
		int taugridn = 21;
		vector<vector<double>> tauden{{0.0, 10}, {20, 10}};
		grid_tauPts = GenerateGrids(tauden, taugridn);

		//pPts:
		int pgridn = 50;
		vector<vector<double>> pden{{1, 8}, {20, 7}, {40, 3}, {100, 5}, {450, 1}};
		grid_pPts = GenerateGrids(pden, pgridn);

		//TPts:
		int Tgridn = 40;
		vector<vector<double>> Tden{{0.01, 10}, {2.000, 10}};
		grid_TPts = GenerateGrids(Tden, Tgridn);

		//xPts:
		double T = grid_TPts[0], nf = 3.0, lambda = 0.2;
		double mu = 0.197*sqrt((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda / ProductLog((-8.0*(6 + nf)*M_PI*M_PI*T*T) / (2.0*nf - 33.0) / lambda / lambda));
		double mg = mu/sqrt(2.0);
		double M = mu/sqrt(6.0);
		double MAXP = 450.0;
		int xgridn = 50;
		double xmin = mg/(MAXP + sqrt(MAXP*MAXP + M*M));
		for (int i=0; i<xgridn; i++)
			grid_xPts.push_back(exp(log(xmin) - log(xmin)/(xgridn-1.0)*i));

		//RadPts:
		int Radgridn = 40;
		vector<vector<double>> Radden{{2, 10}, {50, 10}, {70, 1}, {420, 1}};
		grid_RadPts = GenerateGrids(Radden, Radgridn);

		//FdpPts:
		double mgC = 0.5140824927;
		int Fdpgridn = 22;
		vector<vector<double>> Fdpden = {{5.0*mgC/2.0, 10}, {12.0, 5}, {30, 0.}};
		grid_FdpPts = GenerateGrids(Fdpden, Fdpgridn-4);
		grid_FdpPts.insert(grid_FdpPts.begin(), 4.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 3.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 2.0*mgC/2.0);
		grid_FdpPts.insert(grid_FdpPts.begin(), 1.0*mgC/2.0);

		//pCollPts:
		int pCollgridn = 40;
		vector<vector<double>> pCollden{{1, 10}, {4, 10}, {9, 2.5}, {30, 0.6}, {60, 0.5}, {420, 0.3}};
		grid_pCollPts = GenerateGrids(pCollden, pCollgridn);

		//TCollPts:
		int TCollgridn = 40;
		vector<vector<double>> TCollden{{0.01, 10}, {2.000, 10}};
		grid_TCollPts = GenerateGrids(TCollden, TCollgridn);

		//finpts:
		int fingridn = 50;
		vector<vector<double>> finden{{5, 10}, {50, 10}, {70, 5}, {400, 3}};
		grid_finPts = GenerateGrids(finden, fingridn);
	}

	//rounding grids to 10 decimal points:
	for (int i=0; i<grid_tauPts.size(); i++) grid_tauPts[i] = round(grid_tauPts[i]*1e10)/1e10;
	for (int i=0; i<grid_pPts.size(); i++)   grid_pPts[i]   = round(grid_pPts[i]*1e10)/1e10;
	for (int i=0; i<grid_TPts.size(); i++)   grid_TPts[i]   = round(grid_TPts[i]*1e10)/1e10;
	for (int i=0; i<grid_xPts.size(); i++)   grid_xPts[i]   = round(grid_xPts[i]*1e10)/1e10;
	for (int i=0; i<grid_RadPts.size(); i++) grid_RadPts[i] = round(grid_RadPts[i]*1e10)/1e10;
	for (int i=0; i<grid_FdpPts.size(); i++) grid_FdpPts[i] = round(grid_FdpPts[i]*1e10)/1e10;

	for (int i=0; i<grid_pCollPts.size(); i++) grid_pCollPts[i] = round(grid_pCollPts[i]*1e10)/1e10;
	for (int i=0; i<grid_TCollPts.size(); i++) grid_TCollPts[i] = round(grid_TCollPts[i]*1e10)/1e10;

	for (int i=0; i<grid_finPts.size(); i++) grid_finPts[i] = round(grid_finPts[i]*1e10)/1e10;
}

//DESTRUCTORS:
GridPoints::~GridPoints() {

}

//GRID FUNCTIONS:

vector<double> GridPoints::tauPts() {return grid_tauPts;}
double GridPoints::tauPts(int i) {if (i < 0) return grid_tauPts[grid_tauPts.size()-abs(i)]; return grid_tauPts[i];}
int GridPoints::tauPtsLength() {return grid_tauPts.size();}

vector<double> GridPoints::pPts() {return grid_pPts;}
double GridPoints::pPts(int i) {if (i < 0) return grid_pPts[grid_pPts.size()-abs(i)]; return grid_pPts[i];}
int GridPoints::pPtsLength() {return grid_pPts.size();}

vector<double> GridPoints::TPts() {return grid_TPts;}
double GridPoints::TPts(int i) {if (i < 0) return grid_TPts[grid_TPts.size()-abs(i)]; return grid_TPts[i];}
int GridPoints::TPtsLength() {return grid_TPts.size();}

vector<double> GridPoints::xPts() {return grid_xPts;}
double GridPoints::xPts(int i) {if (i < 0) return grid_xPts[grid_xPts.size()-abs(i)]; return grid_xPts[i];}
int GridPoints::xPtsLength() {return grid_xPts.size();}

vector<double> GridPoints::RadPts() {return grid_RadPts;}
double GridPoints::RadPts(int i) {if (i < 0) return grid_RadPts[grid_RadPts.size()-abs(i)]; return grid_RadPts[i];}
int GridPoints::RadPtsLength() {return grid_RadPts.size();}

vector<double> GridPoints::FdpPts() {return grid_FdpPts;}
double GridPoints::FdpPts(int i) {if (i < 0) return grid_FdpPts[grid_FdpPts.size()-abs(i)]; return grid_FdpPts[i];}
int GridPoints::FdpPtsLength() {return grid_FdpPts.size();}

vector<double> GridPoints::pCollPts() {return grid_pCollPts;}
double GridPoints::pCollPts(int i) {if (i < 0) return grid_pCollPts[grid_pCollPts.size()-abs(i)]; return grid_pCollPts[i];}
int GridPoints::pCollPtsLength() {return grid_pCollPts.size();}

vector<double> GridPoints::TCollPts() {return grid_TCollPts;}
double GridPoints::TCollPts(int i) {if (i < 0) return grid_TCollPts[grid_TCollPts.size()-abs(i)]; return grid_TCollPts[i];}
int GridPoints::TCollPtsLength() {return grid_TCollPts.size();}

vector<double> GridPoints::finPts() {return grid_finPts;}
double GridPoints::finPts(int i) {if (i < 0) return grid_finPts[grid_finPts.size()-abs(i)]; return grid_finPts[i];}
int GridPoints::finPtsLength() {return grid_finPts.size();}

//PRIVATE FUNCTIONS:

double GridPoints::LinearIntegrate(vector<double> data_x, vector<double> data_f, double xH)
{
	vector<double> k, c;
	for (int i=0; i<(data_x.size()-1); i++)
	{
		k.push_back((data_f[i+1]-data_f[i])/(data_x[i+1]-data_x[i]));
		c.push_back(data_f[i]-k.back()*data_x[i]);
	}

	int xHi = 0; while (xH > data_x[xHi]) xHi++; xHi--;

	double sum = 0.0;

	for (int i=0; i<xHi; i++)
	{
		sum += 0.5*k[i]*(data_x[i+1]*data_x[i+1] - data_x[i]*data_x[i]) + c[i]*(data_x[i+1] - data_x[i]);
	}

	sum += 0.5*k[xHi]*(xH*xH - data_x[xHi]*data_x[xHi]) + c[xHi]*(xH - data_x[xHi]);

	return sum;
}

vector<double> GridPoints::GenerateGrids(vector<vector<double>> density, int numpts)
{
	vector<double> density_x, density_f;
	for (int i=0; i<density.size(); i++) {density_x.push_back(density[i][0]); density_f.push_back(density[i][1]);}

	vector<double> inttab_x, inttab_f;

	double xxx = density_x.front();
	inttab_x.push_back(0.0);
	inttab_f.push_back(xxx);

	for (int i=1; i<20; i++)
	{
		xxx = density_x.front() + (density_x.back()-density_x.front())*i/19.0;
		inttab_x.push_back(LinearIntegrate(density_x, density_f, xxx));
		inttab_f.push_back(xxx);
	}

	interpFun inttab_int(inttab_x, inttab_f);

	vector<double> gridpoints;

	gridpoints.push_back(density_x.front());

	for (int i=1; i<numpts; i++)
	{
		double a = 0.0 + inttab_x.back()*i/(numpts-1.0);
		gridpoints.push_back(inttab_int.interp(a));
	}

	return gridpoints;
}