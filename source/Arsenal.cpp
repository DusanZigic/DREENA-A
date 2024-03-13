#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cfloat>

using namespace std;

#include "ELossHeader.hpp"
#include "Arsenal.hpp"
#include "ImportExport.hpp"
#include "LinearInterpolation.hpp"

//function that determines header length of file:
int GetHeaderLength(ifstream &myfile)
//myfile - file to determine header length of <- input
//return value: header length
{
	int header_len = 0;
	
	string line;

	while (getline(myfile, line))
	{
		int cnt = 0;
    	while (isspace(line[cnt])) cnt++;

    	if (isdigit(line[cnt])) break;

    	header_len++;
	}

	//returning input file to begining:
	myfile.clear(); myfile.seekg(0);

	return header_len;
}

//ProductLog function (needed for mu)
double ProductLog(double x)
{
	if (x == 0) {
		return 0;
	}

	double w0, w1;
	if (x > 0) {
		w0 = log(1.2 * x / log(2.4 * x / log1p(2.4 * x)));
	}
	else {
		double v = 1.4142135623730950488 * sqrt(1 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (true) {
		double e = exp(w0);
		double f = w0 * e - x;
		w1 = w0 - f / ((e * (w0 + 1) - (w0 + 2) * f / (w0 + w0 + 2)));
		if (fabs(w0 / w1 - 1) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}
	return w1;
}

//UnitStep function (for double and overloaded for long double)
double UnitStep(double x)
{
	if (x < 0.0) {
		return 0.0;
	}
	else {
		return 1.0;
	}
}
long double UnitStep(long double x)
{
	if (x < 0.0L) {
		return 0.0;
	}
	else {
		return 1.0;
	}
}

//function that generates Halton sequence of quasi random numbers
double HaltonSequence(int index, int base)
{
	double f = 1;
	double res = 0;

	while (index > 0) {
		f = f / base;
		res += f * (index % base);
		index = index / base; // integer division
	}

	return res;
}

//function that locates point in vector:
int LocatePoint(vector<double> data, double x, int int_ord)
//data    - vector in which to search for <- input
//x       - point value                   <- input
//int_ord - interpolation order 		  <- input
//return value: point position
//function taken from recipes in C++
{
	int ju, jm, jl;
	int mm = int_ord + 1;
	int n = data.size();
	int mmmin;
	bool ascnd = (data.back() >= data.front());
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1)
	{
		jm = (ju + jl) >> 1;
		if ((x >= data[jm]) == ascnd) {
			jl = jm;
		}
		else {
			ju = jm;
		}
	}
	mmmin = n - mm < jl - ((mm - 2) >> 1) ? n - mm : jl - ((mm - 2) >> 1);
	return (0 > mmmin ? 0 : mmmin);
}

//function that findes coefficients of the polynomial (order of the polynomila is data length -1)
void PolynomialCoeff(vector<double> data_x, vector<double> data_f, vector<double> &coeff)
//data_x - points values, data_f - function values <-input
//coeff - polynomial coefficients <- output
//function taken from recipes in C++
{
	int n = data_x.size();

	coeff.resize(n); fill(coeff.begin(), coeff.end(), 0.0);

	int k, j, i;
	double phi, ff, b;
	
	vector<double> s(n); fill(s.begin(), s.end(), 0.0);
	
	s[n - 1] = -data_x[0];

	for (i = 1; i<n; i++)
	{
		for (j = n - 1 - i; j<n - 1; j++)
		{
			s[j] -= data_x[i] * s[j + 1];
		}
		s[n - 1] -= data_x[i];
	}

	for (j = 0; j<n; j++)
	{
		phi = n;
		
		for (k = n - 1; k>0; k--)
		{
			phi = k * s[k] + data_x[j] * phi;
		}
		
		ff = data_f[j] / phi;
		
		b = 1.0;
		
		for (k = n - 1; k >= 0; k--)
		{
			coeff[k] += b * ff;
			b = s[k] + data_x[j] * b;
		}
	}
}

//function that integrates linear interpolated polynomial:
double LinearIntegrate(vector<double> xdata, vector<double> fdata)
//xdata - points, fdata - function values <- input
//return value: integrated value
{
	if (xdata.size() < 2) return 0.0;

	vector<double> k, c;
	for (int i=0; i<(xdata.size()-1); i++)
	{
		k.push_back((fdata[i+1]-fdata[i])/(xdata[i+1]-xdata[i]));
		c.push_back(fdata[i]-k.back()*xdata[i]);
	}

	double res = 0.0;

	for (int i=0; i<(xdata.size()-1); i++)
		res += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	return res;
}

//function that integrates linear interpolated polynomial within given range:
double LinearIntegrate(vector<double> xdata, vector<double> fdata, double lim_l, double lim_h)
//xdata - points values 		  <- input
//fdata - function values 		  <- input
//lim_l - lower integration limit <- input
//lim_h - higer integration limit <- input
//return value: value of the integral
{
	if (xdata.size() < 2) return 0.0;

	vector<double> k, c;
	for (int i=0; i<(xdata.size()-1); i++)
	{
		k.push_back((fdata[i+1]-fdata[i])/(xdata[i+1]-xdata[i]));
		c.push_back(fdata[i]-k.back()*xdata[i]);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0;

	for (int i=0; i<(xdata.size()-1); i++)
		sum += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	int lim_l_pos = LocatePoint(xdata, lim_l, 1);

	double sum_l = 0.0;

	for (int i=0; i<lim_l_pos; i++)
		sum_l += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	sum_l += 0.5*k[lim_l_pos]*(lim_l*lim_l - xdata[lim_l_pos]*xdata[lim_l_pos]) + c[lim_l_pos]*(lim_l - xdata[lim_l_pos]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	int lim_h_pos = LocatePoint(xdata, lim_h, 1);

	double sum_h = 0.0;

	sum_h += 0.5*k[lim_h_pos]*(xdata[lim_h_pos+1]*xdata[lim_h_pos+1] - lim_h*lim_h) + c[lim_h_pos]*(xdata[lim_h_pos+1] - lim_h);

	for (int i=lim_h_pos+1; i<xdata.size()-1; i++)
		sum_h += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - sum_h - sum_l);
}

//function that integrates 3rd order polynomial in it's domain:
double CubicIntegrate(vector<double> xdata, vector<double> fdata)
//xdata - points values 		  <- input
//fdata - function values 		  <- input
//return value: value of the integral
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	vector<vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (int i=0; i<coefficents.size(); i++)
	{
		int pt_loc = LocatePoint(xdata, xdata[i], 3);
		vector<double> xdatatemp(xdata.begin()+pt_loc, xdata.begin()+pt_loc+4);
		vector<double> fdatatemp(fdata.begin()+pt_loc, fdata.begin()+pt_loc+4);
		PolynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//calculating value of integral:
	double sum = 0.0;
	for (int i=0; i<xdata.size()-1; i++)
		for (int j=0; j<coefficents[i].size(); j++)
			sum += 1.0/(j+1)*coefficents[i][j]*(pow(xdata[i+1], j+1) - pow(xdata[i], j+1));

	return sum;
}

//function that integrates 3rd order polynomial within limits:
double CubicIntegrate(vector<double> xdata, vector<double> fdata, double lim_l, double lim_h)
//xdata - points values 		  <- input
//fdata - function values 		  <- input
//lim_l - lower integration limit <- input
//lim_h - higer integration limit <- input
//return value: value of the integral
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	vector<vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (int i=0; i<coefficents.size(); i++)
	{
		int pt_loc = LocatePoint(xdata, xdata[i], 3);
		vector<double> xdatatemp(xdata.begin()+pt_loc, xdata.begin()+pt_loc+4);
		vector<double> fdatatemp(fdata.begin()+pt_loc, fdata.begin()+pt_loc+4);
		PolynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0;
	for (int i=0; i<xdata.size()-1; i++)
		for (int j=0; j<coefficents[i].size(); j++)
			sum += 1.0/(j+1)*coefficents[i][j]*(pow(xdata[i+1], j+1) - pow(xdata[i], j+1));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	int lim_l_pos = LocatePoint(xdata, lim_l, 1);

	double sum_l = 0.0;

	for (int i=0; i<lim_l_pos; i++)
		for (int j=0; j<coefficents[i].size(); j++)
			sum_l += 1.0/(j+1)*coefficents[i][j]*(pow(xdata[i+1], j+1) - pow(xdata[i], j+1));

	for (int j=0; j<coefficents[lim_l_pos].size(); j++)
		sum_l += 1.0/(j+1)*coefficents[lim_l_pos][j]*(pow(lim_l, j+1) - pow(xdata[lim_l_pos], j+1));	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	int lim_h_pos = LocatePoint(xdata, lim_h, 1);

	double sum_h = 0.0;

	for (int j=0; j<coefficents[lim_h_pos].size(); j++)
		sum_h += 1.0/(j+1)*coefficents[lim_h_pos][j]*(pow(xdata[lim_h_pos+1], j+1) - pow(lim_h, j+1));

	for (int i=lim_h_pos+1; i<xdata.size()-1; i++)
		for (int j=0; j<coefficents[i].size(); j++)
			sum_h += 1.0/(j+1)*coefficents[i][j]*(pow(xdata[i+1], j+1) - pow(xdata[i], j+1));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - sum_h - sum_l);
}

//function that generates initial position points and angles:
int GenerateInitPosPoints()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	if (yGridN == -2) {
		//if yGridN is set to -2, MonteCarlo method is used to generate initial position points and angles
		//number of x-y initial position points is equal to xGridN

		interpFun BCDensInt; if (LoadBinCollDensity(BCDensInt) != 1) return -1; //loading binary collisio density
		
		vector<vector<double>> BCDensDomain = BCDensInt.domain(); //getting binary collision density domain
		vector<double> BCDensCoDomain = BCDensInt.codomain();	  //getting binary collision density codomain

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating x and y points:

		random_device rd1; mt19937 mt1(rd1()); uniform_real_distribution<double> distX(BCDensDomain[0][0], nextafter(BCDensDomain[0][1], DBL_MAX)); //defining random generator for x points
		random_device rd2; mt19937 mt2(rd2()); uniform_real_distribution<double> distY(BCDensDomain[1][0], nextafter(BCDensDomain[1][1], DBL_MAX)); //defining random generator for y points
		random_device rd3; mt19937 mt3(rd3()); uniform_real_distribution<double> distZ(BCDensCoDomain[0],  nextafter(BCDensCoDomain[1],  DBL_MAX)); //defining random generator for z points
	
		double x, y, z; //defining buffers

		for (int i=0; i<xGridN; i++) //generating initial x and y points
		{
			do {
				x = distX(mt1);
				y = distY(mt2);
				z = distZ(mt3);
			}
			while (z > BCDensInt.interp(x, y));

			xGridPts.push_back(x); yGridPts.push_back(y); //setting x and y points values
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating phi points:

		random_device rd4; mt19937 mt4(rd4()); uniform_real_distribution<double> distPhi(0.0, nextafter(2.0*M_PI, DBL_MAX)); //defining random generator for phi points

		double phi; //defining phi buffer

		phi = 0.0; 		phiGridPts.push_back(phi); sort(phiGridPts.begin(), phiGridPts.end()); //artificially adding 0 and 2Pi to the list
		phi = 2.0*M_PI; phiGridPts.push_back(phi); sort(phiGridPts.begin(), phiGridPts.end());

		//generating other points:
		for (int i=2; i<phiGridN; i++)
		{
			do {phi = distPhi(mt4);} while (binary_search(phiGridPts.begin(), phiGridPts.end(), phi)); //checking if value is already in the list
			phiGridPts.push_back(phi); sort(phiGridPts.begin(), phiGridPts.end());					   //adding value and sorting vector
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating binary collision density with function values 1:

		vector<double> bcd_x, bcd_y, bcd_f;

		for (int i_x=0; i_x<=10; i_x++)
		{
			for (int i_y=0; i_y<=10; i_y++)
			{
				bcd_x.push_back(BCDensDomain[0][0] + (BCDensDomain[0][1]-BCDensDomain[0][0])*i_x/10.0); //setting x points
				bcd_y.push_back(BCDensDomain[1][0] + (BCDensDomain[1][1]-BCDensDomain[1][0])*i_y/10.0); //setting y points
				bcd_f.push_back(1.0);																	//setting bc dens values
			}
		}

		BinCollDensity.SetData(bcd_x, bcd_y, bcd_f); //setting binary collision density interpolated function
	}
	else if (yGridN == -1) {
		//if yGridN is set to -1, MonteCarlo method is used to generate initial position points, while angles are on equidistant grid
		//number of x-y initial position points is equal to xGridN

		interpFun BCDensInt; if (LoadBinCollDensity(BCDensInt) != 1) return -2; //loading binary collisio density
		
		vector<vector<double>> BCDensDomain = BCDensInt.domain(); //getting binary collision density domain
		vector<double> BCDensCoDomain = BCDensInt.codomain();	  //getting binary collision density codomain

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating x and y points:

		random_device rd1; mt19937 mt1(rd1()); uniform_real_distribution<double> distX(BCDensDomain[0][0], nextafter(BCDensDomain[0][1], DBL_MAX)); //defining random generator for x points
		random_device rd2; mt19937 mt2(rd2()); uniform_real_distribution<double> distY(BCDensDomain[1][0], nextafter(BCDensDomain[1][1], DBL_MAX)); //defining random generator for y points
		random_device rd3; mt19937 mt3(rd3()); uniform_real_distribution<double> distZ(BCDensCoDomain[0],  nextafter(BCDensCoDomain[1],  DBL_MAX)); //defining random generator for z points
	
		double x, y, z; //defining buffers

		for (int i=0; i<xGridN; i++) //generating initial x and y points
		{
			do {
				x = distX(mt1);
				y = distY(mt2);
				z = distZ(mt3);
			}
			while (z > BCDensInt.interp(x, y));

			xGridPts.push_back(x); yGridPts.push_back(y); //setting x and y points values
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating phi points:

		for (int i_phi=0; i_phi<phiGridN; i_phi++) phiGridPts.push_back(2.0*M_PI*i_phi/(phiGridN-1.0));

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating binary collision density with function values 1:

		vector<double> bcd_x, bcd_y, bcd_f;

		for (int i_x=0; i_x<=10; i_x++)
		{
			for (int i_y=0; i_y<=10; i_y++)
			{
				bcd_x.push_back(BCDensDomain[0][0] + (BCDensDomain[0][1]-BCDensDomain[0][0])*i_x/10.0); //setting x points
				bcd_y.push_back(BCDensDomain[1][0] + (BCDensDomain[1][1]-BCDensDomain[1][0])*i_y/10.0); //setting y points
				bcd_f.push_back(1.0);																	//setting bc dens values
			}
		}

		BinCollDensity.SetData(bcd_x, bcd_y, bcd_f); //setting binary collision density interpolated function
	}
	else if (yGridN == 0) {
		//if yGridN is set to 0, initial position points are randomly selected from a list, while angles are on equidistant grid
		//number of x-y initial position points is equal to xGridN

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//selecting x and y points:

		vector<vector<double>> bcpts; if (LoadBinCollPoints(bcpts) != 1) return -3; //importing binary collision points

		if ((xGridN == bcpts.size()) || (xGridN == 0)) { //take all points if xGridN is equal to total number of points or 0
			for (int i=0; i<xGridN; i++)
			{
				xGridPts.push_back(bcpts[i][0]); yGridPts.push_back(bcpts[i][1]);
			}
		}

		random_device rd; auto rng = default_random_engine { rd() }; //defining random device
		shuffle(begin(bcpts), end(bcpts), rng);				 		 //shuffling vector

		for (int i=0; i<xGridN; i++) //setting points values
		{
			xGridPts.push_back(bcpts[i][0]); yGridPts.push_back(bcpts[i][1]);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating phi points:

		for (int i_phi=0; i_phi<phiGridN; i_phi++) phiGridPts.push_back(2.0*M_PI*i_phi/(phiGridN-1.0));

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating binary collision density with function values 1:

		vector<double> bcd_x, bcd_y, bcd_f;

		for (int i_x=0; i_x<=10; i_x++)
		{
			for (int i_y=0; i_y<=10; i_y++)
			{
				bcd_x.push_back(-20.0 + 40.0*i_x/10.0); //setting x points
				bcd_y.push_back(-20.0 + 40.0*i_y/10.0); //setting y points
				bcd_f.push_back(1.0);					//setting bc dens values
			}
		}

		BinCollDensity.SetData(bcd_x, bcd_y, bcd_f); //setting binary collision density interpolated function

	}
	else {
		//if yGridN is larger than 0, initial position points and angles are generated on equidistant grid
		//number of x-y initial position points is equal to (xGridN+1)*(yGridN+1)

		if (LoadBinCollDensity() != 1) return -4; //loading binary collision density

		vector<vector<double>> BCDensDomain = BinCollDensity.domain(); //getting binary collision density domain

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating x and y points:

		//creating range for x-y grids:
		double init_grid_range;
		if (abs(BCDensDomain[0][0]) > BCDensDomain[0][1]) init_grid_range = abs(BCDensDomain[0][0]) - 0.5;
		else init_grid_range = abs(BCDensDomain[0][1]) - 0.5;

		double x, y; //setting x and y buffers

		for (int i_x=0; i_x<=xGridN; i_x++)
			for (int i_y=0; i_y<=yGridN; i_y++)
				{ xGridPts.push_back(-1.0*init_grid_range + 2.0*i_x*init_grid_range/xGridN); yGridPts.push_back(-1.0*init_grid_range + 2.0*i_y*init_grid_range/yGridN);}
	

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//generating phi points:

		for (int i_phi=0; i_phi<phiGridN; i_phi++) phiGridPts.push_back(2.0*M_PI*i_phi/(phiGridN-1.0));
	}

	return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gauss filter integration functions:

//function that generates sampling points for Gaussian integration:
void GenerateGaussTab(vector<double> &qGTab, vector<double> &fGTab)
//qGTab, fGTab - vectors that store sampling point <- output
{	
	double sigmaNum = 3.5; //setting sigma
	double sigmaStep = 0.25; //setting step
	int GTabLen = 2 * (int)(sigmaNum / sigmaStep) + 1; //setting length of sampling points
	
	double GaussTabSum = 0.0; //setting normalization sum to zero
	
	for (int g_i=0; g_i<GTabLen; g_i++) //calculating sampling points
	{
		qGTab.push_back(-1.0*sigmaNum + g_i*sigmaStep);		  //setting qGaussTab values
		fGTab.push_back(exp(-qGTab.back()*qGTab.back()/2.0)); //setting fGaussTab values
		GaussTabSum += fGTab.back();						  //adding to normalization sum
	}
	
	for (int g_i=0; g_i<GTabLen; g_i++)  //normalizing
	{
		fGTab[g_i] /= GaussTabSum; //dividing fGaussTab values with total sum
	}
}

//function that performs Gauss filter integration - modefied pT integration algorithm
void GaussFilterIntegrate(vector<double> radiativeRAA1, vector<vector<double>> radiativeRAA2, vector<double> collisionalEL, vector<double> &singRAA1, vector<vector<double>> &singRAA2)
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
	interpFun muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	vector<double> qGaussTabOG, fGaussTabOG; 	//defining vectors that will store original Gauss filter sampling points
	GenerateGaussTab(qGaussTabOG, fGaussTabOG); //generating sampling points and settin number of sampling poins

	vector<double> qGaussTab, fGaussTab; 		//defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
		interpFun RadRelInt(Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (auto pT : Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interp(pT);

			sigmaColl = sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						   //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						   //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (int g_i=0; g_i<qGaussTab.size(); g_i++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[g_i];			
				GFSum += (dsdpti2.interp(pT + dppT)*RadRelInt.interp(pT + dppT)*(pT + dppT) / pT * fGaussTab[g_i]);
			}

			singRAA1.push_back(1.0 / dsdpti2.interp(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpFun RadRelInt(Grids.RadPts(), Grids.FdpPts(), radiativeRAA2); //creating radiative RAA2 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (auto pT : Grids.finPts())
		{
			singRAA2.push_back(vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interp(pT);

			sigmaColl = sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						   //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						   //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (auto dpT : Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (int g_i=0; g_i<qGaussTab.size(); g_i++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[g_i];
					GFSum += (dsdpti2.interp(pT + dpT + dppT)*RadRelInt.interp(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[g_i]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2.interp(pT) * GFSum);
			}
		}
	}
}

//function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm:
void GaussFilterIntegrate(interpFun dsdpti2lquark, vector<double> radiativeRAA1, vector<vector<double>> radiativeRAA2, vector<double> collisionalEL, vector<double> &singRAA1, vector<vector<double>> &singRAA2)
//dsdpti2lquark - light quark initial pT distribution      						  <- input
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
	interpFun muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	vector<double> qGaussTabOG, fGaussTabOG; 	//defining vectors that will store original Gauss filter sampling points
	GenerateGaussTab(qGaussTabOG, fGaussTabOG); //generating sampling points and settin number of sampling poins

	vector<double> qGaussTab, fGaussTab; 		//defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
		interpFun RadRelInt(Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (auto pT : Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interp(pT);

			sigmaColl = sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						   //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						   //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (int g_i=0; g_i<qGaussTab.size(); g_i++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[g_i];			
				GFSum += (dsdpti2lquark.interp(pT + dppT)*RadRelInt.interp(pT + dppT)*(pT + dppT) / pT * fGaussTab[g_i]);
			}

			singRAA1.push_back(1.0 / dsdpti2lquark.interp(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpFun RadRelInt(Grids.RadPts(), Grids.FdpPts(), radiativeRAA2); //creating radiative RAA2 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (auto pT : Grids.finPts())
		{
			singRAA2.push_back(vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interp(pT);

			sigmaColl = sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						   //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						   //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	   //setting rescaling factor
				for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (auto dpT : Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (int g_i=0; g_i<qGaussTab.size(); g_i++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[g_i];
					GFSum += (dsdpti2lquark.interp(pT + dpT + dppT)*RadRelInt.interp(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[g_i]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2lquark.interp(pT) * GFSum);
			}
		}
	}
}

//function that performs Gauss filter integration - default algorithm:
void GaussFilterIntegrate(vector<double> radiativeRAA, vector<double> collisionalEL, vector<double> &singRAA)
//radiativeRAA  - raditive RAA 							   <- input
//collisionalEL - collisional energy loss				   <- input
//singRAA 		- RAA array after Gauss filter integration <- output
{
	interpFun RadRelInt(Grids.RadPts(),   radiativeRAA);  //creating radiative RAA interpolated function
	interpFun muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	vector<double> qGaussTabOG, fGaussTabOG; 	//defining vectors that will store original Gauss filter sampling points
	GenerateGaussTab(qGaussTabOG, fGaussTabOG); //generating sampling points and settin number of sampling poins

	vector<double> qGaussTab, fGaussTab; 		//defining vectors that will store Gauss filter sampling points

	long double GFSum; //defining sum variable for Gauss filter

	long double pT, dpT; //defining pT and dpT variables

	long double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value

	long double sigmaColl; //defining variable for collisional sigma
	
	//Gauss filter
	//for (auto pT : Grids.finPts())
	for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
	{
		pT = Grids.finPts(f_i);

		GFSum = 0.0L;

		muCollCurrVal = muCollInt.interp(pT);

		sigmaColl = sqrt(2.0*TCollConst*muCollCurrVal);

		qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

		if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						   //checking if Gauss is out of bound on lower bound
			double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	   //setting rescaling factor
			for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}		
		
		if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						   //checking if Gauss is out of bound on upper bound
			double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	   //setting rescaling factor
			for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}
		
		//calculating Gauss filter
		for (int g_i=0; g_i<qGaussTab.size(); g_i++)
		{
			dpT = muCollCurrVal + sigmaColl * qGaussTab[g_i];			
			GFSum += (dsdpti2.interp(pT + dpT)*RadRelInt.interp(pT + dpT)*(pT + dpT) / pT * fGaussTab[g_i]);
		}

		singRAA.push_back(1.0 / dsdpti2.interp(pT) * GFSum);
	}
}

//function that calculates observables RAA(pT) and v2(pT):
void CalcObservables(vector<vector<double>> RAApTphi, vector<double> &RAApT, vector<double> &v2pT)
//RAApTphi 	  - RAA(pT,phi) distribution <- input
//RAApT, v2pT - RAA(pT), v2(pT)          <- output
{
	for (int p_i=0; p_i<Grids.finPtsLength(); p_i++)
	{
		vector<double> RAAphi(RAApTphi[p_i]);
		
		double norm = CubicIntegrate(phiGridPts, RAAphi);

		RAApT.push_back(norm/2.0/M_PI);

		interpFun RAAphiInt(phiGridPts, RAAphi);

		int phiIntN = 100;

		vector<double> phi_integrate(phiIntN), psi2s_integrate(phiIntN), psi2c_integrate(phiIntN), v2_integrate(phiIntN);

		/*
		phi_integrate.push_back(phiGridPts.front());
		v2_integrate.push_back(RAAphiInt.interp(phi_integrate.back())*cos(2.0*phi_integrate.back()));

		for (int phi_i=1; phi_i<(phiIntN-1); phi_i++)
		{
			phi_integrate.push_back(phiGridPts.front() + (phiGridPts.back()-phiGridPts.front())*phi_i/(phiIntN-1.0));
			 v2_integrate.push_back(RAAphiInt.interp(phi_integrate.back())*cos(2.0*phi_integrate.back()));
		}
		*/

		for (int phi_i=0; phi_i<phiIntN; phi_i++)
		{
			  phi_integrate[phi_i] = phiGridPts.front() + 1.0*phi_i/(phiIntN-1.0)*(phiGridPts.back()-phiGridPts.front());
			psi2s_integrate[phi_i] = RAAphiInt.interp(phi_integrate[phi_i])*sin(2.0*phi_integrate[phi_i]);
			psi2c_integrate[phi_i] = RAAphiInt.interp(phi_integrate[phi_i])*cos(2.0*phi_integrate[phi_i]);
		}

		double psi2 = 1/2.0*atan(CubicIntegrate(phi_integrate, psi2s_integrate)/CubicIntegrate(phi_integrate, psi2c_integrate));
		
		for (int phi_i=0; phi_i<phiIntN; phi_i++)
		{
			phi_integrate[phi_i] = phiGridPts.front() + 1.0*phi_i/(phiIntN-1.0)*(phiGridPts.back()-phiGridPts.front());
			 v2_integrate[phi_i] = RAAphiInt.interp(phi_integrate[phi_i])*cos(2.0*phi_integrate[phi_i] - 2.0*psi2);
		}
		
		v2pT.push_back(CubicIntegrate(phi_integrate, v2_integrate)/norm);
	}
	
}

//function that calculates path-lengths and temperatures:
void CalcAvgPLT(double *pl_dist_a, double *temp_dist_a, vector<double> &avg_pl, vector<double> &avg_temp)
//pl_dist_a   - path-length angle distribution 						  <- input
//temp_dist_a - temperature angle distribution 						  <- input
//avg_pl      - angular-averaged, in-plane and out-plane path-lengths <- output
//avg_temp    - angular-averaged, in-plane and out-plane temperatures <- output
{
	vector<double> pl_dist(pl_dist_a, pl_dist_a+phiGridPts.size());											//setting path-length distribution vector
	interpFun pl_dist_int(phiGridPts, pl_dist);																//setting path-length distribution interpolated function
	avg_pl.push_back(CubicIntegrate(phiGridPts, pl_dist)/2.0/M_PI);											//calculating angular-averaged path-length
	avg_pl.push_back((pl_dist_int.interp(phiGridPts.front()) + pl_dist_int.interp(phiGridPts.back()))/2.0); //calculating in-plane path-length
	avg_pl.push_back((pl_dist_int.interp(M_PI/2.0)           + pl_dist_int.interp(3.0*M_PI/2.0))/2.0); 		//calculating out-of-plane path-length

	vector<double> temp_dist(temp_dist_a, temp_dist_a+phiGridPts.size());										  //setting temperature distribution vector
	interpFun temp_dist_int(phiGridPts, temp_dist);																  //setting temperature distribution interpolated function
	avg_temp.push_back(CubicIntegrate(phiGridPts, temp_dist)/2.0/M_PI);											  //calculating angular-averaged temperature
	avg_temp.push_back((temp_dist_int.interp(phiGridPts.front()) + temp_dist_int.interp(phiGridPts.back()))/2.0); //calculating in-plane temperature
	avg_temp.push_back((temp_dist_int.interp(M_PI/2.0)           + temp_dist_int.interp(3.0*M_PI/2.0))/2.0); 	  //calculating out-of-plane temperature
}