#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cctype>

using namespace std;

#include "ELossHeader.hpp"
#include "LTablesHeader.hpp"
#include "Arsenal.hpp"
#include "LinearInterpolation.hpp"
#include "Grids.hpp"

//function that loads dsdpti2 table and interpolates it:
int Loaddsdpti2()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	//setting file path:
	string path_in;
	if (pTinit_path.length() == 0) path_in = "pTinitDists/pTinitDist_" + pName + ".dat";
	else  						   path_in = pTinit_path;

	ifstream file_in(path_in); 												   //opening file
	if (!file_in.is_open()) {												   //checking if file is open
		cerr << "Error: unable to open initial pT distribution file." << endl; //if not, priting error and
		return -1;															   //and return -1
	}

	vector<double> pTdist_x, pTdist_f; //defining vectors that store dsdpti2 values

	string line; //string type buffer
	double d;	 //double type buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	while (getline(file_in, line))
	{
		stringstream ss(line);
		ss >> d; pTdist_x.push_back(d);
		ss >> d; pTdist_f.push_back(d); //importing values into vectors
	}

	dsdpti2.SetData(pTdist_x, pTdist_f); //setting dsdpti2 interpolated function

	file_in.close(); //closing file

	return 1;
}

//function that loads dsdpti2 table and interpolates it for given particle:
int Loaddsdpti2(string pname, interpFun &dsdpti2int)
//this function is used for algorithm that calculates energy loss for all light quarks
//pname - particle name 						  <- input
//dsdpti2in - interpolated intial pT distribution <- output
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string path_in = "pTinitDists/pTinitDist_" + pname + ".dat"; //setting file path

	ifstream file_in(path_in); 												   //opening file
	if (!file_in.is_open()) {												   //checking if file is open
		cerr << "Error: unable to open initial pT distribution file." << endl; //if not, priting error and
		return -1;															   //and return -1
	}

	vector<double> pTdist_x, pTdist_f; //defining vectors that store dsdpti2 values

	string line; //string type buffer
	double d;	 //double type buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	while (getline(file_in, line))
	{
		stringstream ss(line);
		ss >> d; pTdist_x.push_back(d);
		ss >> d; pTdist_f.push_back(d); //importing values into vectors
	}

	dsdpti2int.SetData(pTdist_x, pTdist_f); //setting dsdpti2 interpolated function

	file_in.close(); //closing file

	return 1;
}

//function that loads Ldndx table and interpolates it:
int LoadLdndx()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string part_name;
	if (pName == "Bottom") part_name = "Bottom";
	else if (pName == "Charm") part_name = "Charm";
	else if (pName == "Gluon") part_name = "Gluon";
	else part_name = "LQuarks";

	stringstream xb; xb << fixed << setprecision(1) << xB;

	string path_in = "LTables/LdndxTbl_" + part_name + "_xB=" + xb.str() + ".dat"; //setting file path

	ifstream file_in(path_in); 									   //opening file
	if (!file_in.is_open()) {									   //checking if file is open
		cerr << "Error: unable to open Ldndx table file." << endl; //if not, priting error and
		return -1;												   //and return -1
	}

	vector<double> Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f; //defining vectors that store Ldndx table values

	string line; //string type buffer
	double d;    //double type buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	while (getline(file_in, line))
	{
		stringstream ss(line);
		ss >> d; Ldndx_tau.push_back(d);
		ss >> d; Ldndx_p.push_back(d);
		ss >> d; Ldndx_T.push_back(d);
		ss >> d; Ldndx_x.push_back(d);
		ss >> d; Ldndx_f.push_back(d); //importing values into vectors
	}

	Ldndx.SetData(Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f); //setting Ldndx interpolated function

	file_in.close(); //closing file

	return 1;
}

//function that loads LNorm table and interpolates it:
int LoadLNorm()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string part_name;
	if (pName == "Bottom") part_name = "Bottom";
	else if (pName == "Charm") part_name = "Charm";
	else if (pName == "Gluon") part_name = "Gluon";
	else part_name = "LQuarks";

	stringstream xb; xb << fixed << setprecision(1) << xB;

	string path_in = "LTables/LNormTbl_" + part_name + "_xB=" + xb.str() + ".dat"; //setting file path

	ifstream file_in(path_in); 									   //opening file
	if (!file_in.is_open()) {									   //checking if file is open
		cerr << "Error: unable to open LNorm table file." << endl; //if not, priting error and
		return -1;												   //and return -1
	}

	vector<double> LNorm_tau, LNorm_p, LNorm_T, LNorm_f; //defining vectors that store LNorm table values

	string line; //string type buffer
	double d;    //double type buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	while (getline(file_in, line))
	{
		stringstream ss(line);
		ss >> d; LNorm_tau.push_back(d);
		ss >> d; LNorm_p.push_back(d);
		ss >> d; LNorm_T.push_back(d);
		ss >> d; LNorm_f.push_back(d); //importing values into vectors
	}

	LNorm.SetData(LNorm_tau, LNorm_p, LNorm_T, LNorm_f); //setting LNorm interpolated function

	file_in.close(); //closing file

	return 1;
}

//function that loads LColl table and interpolates it:
int LoadLColl()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string part_name;
	if (pName == "Bottom") part_name = "Bottom";
	else if (pName == "Charm") part_name = "Charm";
	else if (pName == "Gluon") part_name = "Gluon";
	else part_name = "LQuarks";

	string path_in = "LTables/LCollTbl_" + part_name + ".dat"; //setting file path

	ifstream file_in(path_in); 									   //opening file
	if (!file_in.is_open()) {									   //checking if file is open
		cerr << "Error: unable to open LColl table file." << endl; //if not, priting error and
		return -1;												   //and return -1
	}

	vector<double> LColl_p, LColl_T, LColl_f; //defining vectors that store LColl table values

	string line; //string type buffer
	double d;    //double type buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	while (getline(file_in, line))
	{
		stringstream ss(line);
		ss >> d; LColl_p.push_back(d);
		ss >> d; LColl_T.push_back(d);
		ss >> d; LColl_f.push_back(d); //importing values into vectors
	}

	LColl.SetData(LColl_p, LColl_T, LColl_f); //setting LColl interpolated function

	file_in.close(); //closing file

	return 1;
}

//function that loads TProfile and interpolates it:
int LoadTProfile()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	//setting file path
	string path_in;
	if (temp_path.length() == 0) path_in = "TProfiles/TProfile_cent=" + centrality + ".dat";
	else                         path_in = temp_path;

	ifstream file_in(path_in); 									//opening file
	if (!file_in.is_open()) {									//checking if file is open
		cerr << "Error: unable to open TProfile file." << endl; //if not, priting error and
		return -1;												//and return -1
	}

	string line; //string buffer
	double d;    //double buffer

	//checking header length:
	int header_len = GetHeaderLength(file_in);

	//skiping header:
	for (int i=0; i<header_len; i++) getline(file_in, line);
	
	//checking how many columns evolution file has:
	int column_cnt = 0;
	getline(file_in, line);
	stringstream first_line(line); while (first_line >> d) column_cnt++;

	file_in.clear(); file_in.seekg(0);
	for (int i=0; i<header_len; i++) getline(file_in, line); //returning to first line after header

	vector<double> temp_tau, temp_x, temp_y, temp_f1, temp_f2, temp_T; //vectors that store values from file

	if (column_cnt == 4) { //evolution file has 4 columns (just temperature)

		while (getline(file_in, line)) //importing values from file
		{
			stringstream ss(line);
			ss >> d; temp_tau.push_back(d);
			ss >> d; temp_x.push_back(d);
			ss >> d; temp_y.push_back(d);
			ss >> d; temp_T.push_back(d); //importing values into vectors
		}
	}
	else if (column_cnt == 5) { //evolution file has 5 columns (energy density and temperature)

		while (getline(file_in, line)) //importing values from file
		{
			stringstream ss(line);
			ss >> d; temp_tau.push_back(d);
			ss >> d; temp_x.push_back(d);
			ss >> d; temp_y.push_back(d);
			ss >> d; temp_f1.push_back(d);
			ss >> d; temp_f2.push_back(d); //importing values into vectors
		}

		double temp_f1_max = *max_element(begin(temp_f1), end(temp_f1)); //finding maximal values of columns 4 and 5 to see which one is 
		double temp_f2_max = *max_element(begin(temp_f2), end(temp_f2)); //energy density and which one is temperature since energy density has larger values

		if (temp_f1_max > temp_f2_max) temp_T = temp_f2; //setting temp_T to 5th column
		else temp_T = temp_f1; 							 //setting temp_T to 4th column
	}
	else { //evolution file is not suitable for interpolation

		cerr << "Error: not enough columns for temperature evolution interpolation." << endl;
		return -2;
	}


	if (((*min_element(temp_x.begin(), temp_x.end())) >= 0.0) && ((*min_element(temp_y.begin(), temp_y.end())) >= 0.0)) {
		//if temperature evolution is defined only in the first quadrant:

		//creating tau grid:
		vector<double> temp_tau_grid = temp_tau;
		sort(temp_tau_grid.begin(), temp_tau_grid.end());
		temp_tau_grid.erase(unique(temp_tau_grid.begin(), temp_tau_grid.end()), temp_tau_grid.end());

		//creating full x grid:
		vector<double> temp_x_grid = temp_x;
		double factor = -1.0; for_each(temp_x_grid.begin(), temp_x_grid.end(), [factor](double &c){ c *= factor; });
		temp_x_grid.insert(temp_x_grid.end(), temp_x.begin(), temp_x.end());
		sort(temp_x_grid.begin(), temp_x_grid.end());
		temp_x_grid.erase(unique(temp_x_grid.begin(), temp_x_grid.end()), temp_x_grid.end());

		//creating full y grid:
		vector<double> temp_y_grid = temp_y;
		factor = -1.0; for_each(temp_y_grid.begin(), temp_y_grid.end(), [factor](double &c){ c *= factor; });
		temp_y_grid.insert(temp_y_grid.end(), temp_y.begin(), temp_y.end());
		sort(temp_y_grid.begin(), temp_y_grid.end());
		temp_y_grid.erase(unique(temp_y_grid.begin(), temp_y_grid.end()), temp_y_grid.end());

		interpFun TProfileFQ(temp_tau, temp_x, temp_y, temp_T); //creating interpolated temperature evolution defined in first quadrant

		//creating full temperature evolution table:
		vector<double> temp_tau_full, temp_x_full, temp_y_full, temp_T_full;
		for (auto tau : temp_tau_grid)
			for (auto x : temp_x_grid)
				for (auto y : temp_y_grid)
				{
					temp_tau_full.push_back(tau);
					  temp_x_full.push_back(x);
					  temp_y_full.push_back(y);
					  temp_T_full.push_back(TProfileFQ.interp(tau, abs(x), abs(y)));
				}

		TProfile.SetData(temp_tau_full, temp_x_full, temp_y_full, temp_T_full); //setting temperature evolution interpolated function
	}
	else {
		//if not, creating interpolated function with values from file:

		TProfile.SetData(temp_tau, temp_x, temp_y, temp_T); //setting temperature evolution interpolated function
	}

	tau0 = TProfile.domain()[0][0]; //setting value of tau0

	file_in.close(); //closing file

	return 1;
}

//function that loads binary collision density and interpolates it:
int LoadBinCollDensity()
//return value: 1 if execution is succesfull, negative integer otherwise
{
	//setting file path:
	string path_in;
	if (bcd_path.length() == 0) path_in = "BinaryCollDensities/BinaryCollDensity_cent=" + centrality + ".dat";
	else                        path_in = bcd_path;

	ifstream file_in(path_in); 													//opening file
	if (!file_in.is_open()) {													//checking if file is open
		cerr << "Error: unable to open binary collision density file." << endl; //if not, priting error and
		return -1;																//and return -1
	}

	string line; //string buffer
	double d;    //double buffer

	vector<double> bcd_x, bcd_y, bcd_f; //vectors that store values from file

	while (getline(file_in, line)) //importing values from file
	{
		stringstream ss(line);
		ss >> d; bcd_x.push_back(d);
		ss >> d; bcd_y.push_back(d);
		ss >> d; bcd_f.push_back(d); //importing values into vectors
	}

	file_in.close(); //closing file

	if (((*min_element(bcd_x.begin(), bcd_x.end())) >= 0.0) && ((*min_element(bcd_y.begin(), bcd_y.end())) >= 0.0)) {
		//if binary collision density is defined only in the first quadrant:

		//creating full x grid:
		vector<double> bcd_x_grid = bcd_x;
		double factor = -1.0; for_each(bcd_x_grid.begin(), bcd_x_grid.end(), [factor](double &c){ c *= factor; });
		bcd_x_grid.insert(bcd_x_grid.end(), bcd_x.begin(), bcd_x.end());
		sort(bcd_x_grid.begin(), bcd_x_grid.end());
		bcd_x_grid.erase(unique(bcd_x_grid.begin(), bcd_x_grid.end()), bcd_x_grid.end());

		//creating full x grid:
		vector<double> bcd_y_grid = bcd_y;
		factor = -1.0; for_each(bcd_y_grid.begin(), bcd_y_grid.end(), [factor](double &c){ c *= factor; });
		bcd_y_grid.insert(bcd_y_grid.end(), bcd_y.begin(), bcd_y.end());
		sort(bcd_y_grid.begin(), bcd_y_grid.end());
		bcd_y_grid.erase(unique(bcd_y_grid.begin(), bcd_y_grid.end()), bcd_y_grid.end());

		interpFun BinCollDensFQ(bcd_x, bcd_y, bcd_f); //creating interpolated binary collision density defined in first quadrant

		//creating full binary collision density table:
		vector<double> bcd_x_full, bcd_y_full, bcd_f_full;
		for (auto x : bcd_x_grid)
			for (auto y : bcd_y_grid)
			{
				bcd_x_full.push_back(x);
				bcd_y_full.push_back(y);
				bcd_f_full.push_back(BinCollDensFQ.interp(abs(x), abs(y)));
			}

		BinCollDensity.SetData(bcd_x_full, bcd_y_full, bcd_f_full); //setting binary collision density interpolated function
	}
	else {
		//if not, creating interpolated function with values from file:

		BinCollDensity.SetData(bcd_x, bcd_y, bcd_f); //setting binary collision density interpolated function
	}

	return 1;
}

//function that loads binary collision density and interpolates it for given interpFun object:
int LoadBinCollDensity(interpFun &BinCollDensInt)
//BinCollDensInt - binary collision density interpolated function <- output
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string path_in = "BinaryCollDensities/BinaryCollDensity_cent=" + centrality + ".dat";

	ifstream file_in(path_in); 													//opening file
	if (!file_in.is_open()) {													//checking if file is open
		cerr << "Error: unable to open binary collision density file." << endl; //if not, priting error and
		return -1;																//and return -1
	}

	string line; //string buffer
	double d;    //double buffer

	vector<double> bcd_x, bcd_y, bcd_f; //vectors that store values from file

	while (getline(file_in, line)) //importing values from file
	{
		stringstream ss(line);
		ss >> d; bcd_x.push_back(d);
		ss >> d; bcd_y.push_back(d);
		ss >> d; bcd_f.push_back(d); //importing values into vectors
	}

	file_in.close(); //closing file

	if (((*min_element(bcd_x.begin(), bcd_x.end())) >= 0.0) && ((*min_element(bcd_y.begin(), bcd_y.end())) >= 0.0)) {
		//if binary collision density is defined only in the first quadrant:

		//creating full x grid:
		vector<double> bcd_x_grid = bcd_x;
		double factor = -1.0; for_each(bcd_x_grid.begin(), bcd_x_grid.end(), [factor](double &c){ c *= factor; });
		bcd_x_grid.insert(bcd_x_grid.end(), bcd_x.begin(), bcd_x.end());
		sort(bcd_x_grid.begin(), bcd_x_grid.end());
		bcd_x_grid.erase(unique(bcd_x_grid.begin(), bcd_x_grid.end()), bcd_x_grid.end());

		//creating full x grid:
		vector<double> bcd_y_grid = bcd_y;
		factor = -1.0; for_each(bcd_y_grid.begin(), bcd_y_grid.end(), [factor](double &c){ c *= factor; });
		bcd_y_grid.insert(bcd_y_grid.end(), bcd_y.begin(), bcd_y.end());
		sort(bcd_y_grid.begin(), bcd_y_grid.end());
		bcd_y_grid.erase(unique(bcd_y_grid.begin(), bcd_y_grid.end()), bcd_y_grid.end());

		interpFun BinCollDensFQ(bcd_x, bcd_y, bcd_f); //creating interpolated binary collision density defined in first quadrant

		//creating full binary collision density table:
		vector<double> bcd_x_full, bcd_y_full, bcd_f_full;
		for (auto x : bcd_x_grid)
			for (auto y : bcd_y_grid)
			{
				bcd_x_full.push_back(x);
				bcd_y_full.push_back(y);
				bcd_f_full.push_back(BinCollDensFQ.interp(abs(x), abs(y)));
			}

		BinCollDensInt.SetData(bcd_x_full, bcd_y_full, bcd_f_full); //setting binary collision density interpolated function
	}
	else {
		//if not, creating interpolated function with values from file:

		BinCollDensInt.SetData(bcd_x, bcd_y, bcd_f); //setting binary collision density interpolated function
	}

	return 1;
}

//function that loads binary collision points:
int LoadBinCollPoints(vector<vector<double>> &bcpoints)
//bcpoints - vector that store binary collision points <- output
//return value: 1 if execution is succesfull, negative integer otherwise
{
	string path_in = "BinaryCollPoints/BinaryCollPoints_cent=" + centrality + ".dat"; //setting file path

	ifstream file_in(path_in);												   //opening file
	if (!file_in.is_open()) {												   //checking if file is open
		cerr << "Error: unable to open binary collision points file." << endl; //if not, priting error and
		return -1;															   //and return -1
	}

	string line; //string buffer
	double d;    //double buffer

	vector<double> xpoints, ypoints; //vectors that sore points

	while (getline(file_in, line)) //importing values from file
	{
		stringstream ss(line);
		ss >> d; xpoints.push_back(d);
		ss >> d; ypoints.push_back(d); //importing values into vectors
	}

	bcpoints.resize(xpoints.size()); //resizing output vector

	for (int i=0; i<xpoints.size(); i++)
	{
		bcpoints[i].push_back(xpoints[i]); bcpoints[i].push_back(ypoints[i]);
	}

	file_in.close(); //closing file

	return 1;
}

//function that exports RAA(pT,phi) distribution and RAA(pT), v2(pT) observables:
int ExportResults(vector<vector<double>> RAApTphi, vector<double> RAApT, vector<double> v2pT, vector<double> avg_pl, vector<double> avg_temp)
//RAApTphi - RAA(pT,phi) distribution           <- input
//RAApT, v2pT - RAA(pT) and v2(pT) observables  <- input
//avg_pl, avg_t - path-lengths and temperatures <- input
//return value: 1 if execution is succesfull, negative integer otherwise
{
	//creating header:
	vector<string> header;
	header.push_back("collision_system: Pb+Pb");
	header.push_back("collision_energy: 5.02TeV");
	header.push_back("particle_type: " + pName);
	header.push_back("centrality: " + centrality);

	stringstream xb_str; xb_str << fixed << setprecision(1) << xB;
	header.push_back("xB = " + xb_str.str());

	stringstream avg_pl_str[3]; for (int i=0; i<3; i++) avg_pl_str[i] << fixed << setprecision(6) << avg_pl[i];
	header.push_back("average_path-lengths: " + avg_pl_str[0].str() + ", " + avg_pl_str[1].str() + ", " + avg_pl_str[2].str());

	stringstream avg_temp_str[3]; for (int i=0; i<3; i++) avg_temp_str[i] << fixed << setprecision(6) << avg_temp[i];
	header.push_back("average_temperatures: " + avg_temp_str[0].str() + ", " + avg_temp_str[1].str() + ", " + avg_temp_str[2].str());
	
	if (yGridN <= 0) {
		header.push_back("number_of_angles:    " + to_string(phiGridN));
		header.push_back("number_of_xy_points: " + to_string(xGridN));
	}
	else {
		header.push_back("number_of_angles:      " + to_string(phiGridN));
		header.push_back("number_of_grid_points: " + to_string(xGridN) + ", " + to_string(yGridN));
	}

	{//exporting distribution to a file:
		
		string path_out = "CResults/CResults_" + pName + "/" + pName + "_5TeV_cent=" + centrality + "_xB=" + xb_str.str() + "_dist.dat";

		ofstream file_out(path_out);

		if (!file_out.is_open()) {
			cerr << "Error: unable to open RAA(pT,phi) distribution file." << endl;
			return -1;
		}

		for (auto head : header) file_out << head << endl;

		file_out << "--------------------------------------------------------" << endl;
		file_out << "    pT [GeV]       phi          R_AA   " << endl;

		for (int f_i= 0; f_i<Grids.finPtsLength(); f_i++)
			for (int phi_i=0; phi_i<phiGridN; phi_i++)
				file_out << fixed << setw(14) << setprecision(10) <<    Grids.finPts(f_i) << " "
						 << fixed << setw(12) << setprecision(10) <<    phiGridPts[phi_i] << " "
						 << fixed << setw(12) << setprecision(10) << RAApTphi[f_i][phi_i] << endl;

		file_out.close();
	}

	{//exporting observables to a file:
		
		string path_out = "CResults/CResults_" + pName + "/" + pName + "_5TeV_cent=" + centrality + "_xB=" + xb_str.str() + "_obs.dat";

		ofstream file_out(path_out);

		if (!file_out.is_open()) {
			cerr << "Error: unable to open RAA(pT) and v2(pT) observables file." << endl;
			return -2;
		}

		for (auto head : header) file_out << head << endl;
		
		file_out << "--------------------------------------------------------" << endl;
		file_out << "    pT [GeV]       R_AA         v_2" << endl;

		for (int f_i= 0; f_i<Grids.finPtsLength(); f_i++)
			file_out << fixed << setw(14) << setprecision(10) << Grids.finPts(f_i) << " "
					 << fixed << setw(12) << setprecision(10) <<        RAApT[f_i] << " "
					 << fixed << setw(12) << setprecision(10) <<         v2pT[f_i] << endl;

		file_out.close();
	}

	return 1;
}

//function that exports RAA(pT,phi) distribution and RAA(pT), v2(pT) observables with given particle name - used in lquarks algorithm:
int ExportResults(string part_name, vector<vector<double>> RAApTphi, vector<double> RAApT, vector<double> v2pT, vector<double> avg_pl, vector<double> avg_temp)
//part_name     - particle name 				 <- input
//RAApTphi      - RAA(pT,phi) distribution       <- input
//RAApT, v2pT   - RAA(pT) and v2(pT) observables <- input
//avg_pl, avg_t - path-lengths and temperatures  <- input
//return value: 1 if execution is succesfull, negative integer otherwise
{
	//creating header:
	vector<string> header;
	header.push_back("collision_system: Pb+Pb");
	header.push_back("collision_energy: 5.02TeV");
	header.push_back("particle_type: " + part_name);
	header.push_back("centrality: " + centrality);

	stringstream xb_str; xb_str << fixed << setprecision(1) << xB;
	header.push_back("xB = " + xb_str.str());

	stringstream avg_pl_str[3]; for (int i=0; i<3; i++) avg_pl_str[i] << fixed << setprecision(6) << avg_pl[i];
	header.push_back("average_path-lengths: " + avg_pl_str[0].str() + ", " + avg_pl_str[1].str() + ", " + avg_pl_str[2].str());

	stringstream avg_temp_str[3]; for (int i=0; i<3; i++) avg_temp_str[i] << fixed << setprecision(6) << avg_temp[i];
	header.push_back("average_temperatures: " + avg_temp_str[0].str() + ", " + avg_temp_str[1].str() + ", " + avg_temp_str[2].str());
	
	if (yGridN <= 0) {
		header.push_back("number_of_angles:    " + to_string(phiGridN));
		header.push_back("number_of_xy_points: " + to_string(xGridN));
	}
	else {
		header.push_back("number_of_angles:      " + to_string(phiGridN));
		header.push_back("number_of_grid_points: " + to_string(xGridN) + ", " + to_string(yGridN));
	}

	{//exporting distribution to a file:
		
		string path_out = "CResults/CResults_" + part_name + "/" + part_name + "_5TeV_cent=" + centrality + "_xB=" + xb_str.str() + "_dist.dat";

		ofstream file_out(path_out);

		if (!file_out.is_open()) {
			cerr << "Error: unable to open RAA(pT,phi) distribution file." << endl;
			return -1;
		}

		for (auto head : header) file_out << head << endl;

		file_out << "--------------------------------------------------------" << endl;
		file_out << "    pT [GeV]       phi          R_AA   " << endl;

		for (int f_i= 0; f_i<Grids.finPtsLength(); f_i++)
			for (int phi_i=0; phi_i<phiGridN; phi_i++)
				file_out << fixed << setw(14) << setprecision(10) <<    Grids.finPts(f_i) << " "
						 << fixed << setw(12) << setprecision(10) <<    phiGridPts[phi_i] << " "
						 << fixed << setw(12) << setprecision(10) << RAApTphi[f_i][phi_i] << endl;

		file_out.close();
	}

	{//exporting observables to a file:
		
		string path_out = "CResults/CResults_" + part_name + "/" + part_name + "_5TeV_cent=" + centrality + "_xB=" + xb_str.str() + "_obs.dat";

		ofstream file_out(path_out);

		if (!file_out.is_open()) {
			cerr << "Error: unable to open RAA(pT) and v2(pT) observables file." << endl;
			return -2;
		}

		for (auto head : header) file_out << head << endl;
		
		file_out << "--------------------------------------------------------" << endl;
		file_out << "    pT [GeV]       R_AA         v_2" << endl;

		for (int f_i= 0; f_i<Grids.finPtsLength(); f_i++)
			file_out << fixed << setw(14) << setprecision(10) << Grids.finPts(f_i) << " "
					 << fixed << setw(12) << setprecision(10) <<        RAApT[f_i] << " "
					 << fixed << setw(12) << setprecision(10) <<         v2pT[f_i] << endl;

		file_out.close();
	}

	return 1;
}

//function that export LTables:
int ExportLTables(vector<double> &ldndxtbl, vector<double> &lnormtbl, vector<double> &lcolltbl)
//ldndxtbl, lnormtbl, lcolltbl - LTables to be exported to files    <- input
//vectors are passed by reference to avoid copying because of their size
//return value: 1 if execution is succesfull, negative integer otherwise
{
	stringstream xb_str; xb_str << fixed << setprecision(1) << LT_xB; //setting xB string

	//opening Ldndx export file:
	string path_out_ldndx = "LTables/LdndxTbl_" + LT_pName + "_xB=" + xb_str.str() + ".dat";
	ofstream file_out_ldndx(path_out_ldndx);
	if (!file_out_ldndx.is_open()) {
		cerr << "Error: unable to open Ldndx table export file." << endl;
		return -1;
	}

	//opening LNorm export file:
	string path_out_lnorm = "LTables/LNormTbl_" + LT_pName + "_xB=" + xb_str.str() + ".dat";
	ofstream file_out_lnorm(path_out_lnorm);
	if (!file_out_lnorm.is_open()) {
		cerr << "Error: unable to open LNorm table export file." << endl;
		return -1;
	}

	for (int tau_i=0; tau_i<LT_Grids.tauPtsLength(); tau_i++)
	{
		for (int p_i=0; p_i<LT_Grids.pPtsLength(); p_i++)
		{
			for (int T_i=0; T_i<LT_Grids.TPtsLength(); T_i++)
			{
				for (int x_i=0; x_i<LT_Grids.xPtsLength(); x_i++)
				{
					//calculating ldndx index:
					int ldndx_index = tau_i*LT_Grids.xPtsLength()*LT_Grids.TPtsLength()*LT_Grids.pPtsLength()
									+ p_i*LT_Grids.xPtsLength()*LT_Grids.TPtsLength()
									+ T_i*LT_Grids.xPtsLength()
									+ x_i;

					//printing ldndx to file:
					file_out_ldndx << fixed << setprecision(10) << LT_Grids.tauPts(tau_i) << " "
								   << fixed << setprecision(10) << LT_Grids.pPts(p_i) << " "
								   << fixed << setprecision(10) << LT_Grids.TPts(T_i) << " "
								   << fixed << setprecision(10) << LT_Grids.xPts(x_i) << " "
								   << fixed << setprecision(10) << ldndxtbl[ldndx_index] << endl;
				}

				//calculating lnorm index:
				int lnorm_index = tau_i*LT_Grids.TPtsLength()*LT_Grids.pPtsLength()
								+ p_i*LT_Grids.TPtsLength()
								+ T_i;

				//printing lnorm to file:
				file_out_lnorm << fixed << setprecision(10) << LT_Grids.tauPts(tau_i) << " "
							   << fixed << setprecision(10) << LT_Grids.pPts(p_i) << " "
							   << fixed << setprecision(10) << LT_Grids.TPts(T_i) << " "
							   << fixed << setprecision(10) << lnormtbl[lnorm_index] << endl;
			}
		}
	}

	file_out_ldndx.close(); file_out_lnorm.close(); //closing Ldndx and LNorm files

	//opening LColl export file:
	string path_out_lcoll = "LTables/LCollTbl_" + LT_pName + ".dat";
	ofstream file_out_lcoll(path_out_lcoll);
	if (!file_out_lcoll.is_open()) {
		cerr << "Error: unable to open LColl table export file." << endl;
		return -1;
	}

	//printing lcoll to file:
	for (int p_i=0; p_i<LT_Grids.pCollPtsLength(); p_i++)
	{
		for (int T_i=0; T_i<LT_Grids.TCollPtsLength(); T_i++)
		{
			file_out_lcoll << fixed << setprecision(10) << LT_Grids.pCollPts(p_i) << " "
						   << fixed << setprecision(10) << LT_Grids.TCollPts(T_i) << " "
						   << fixed << setprecision(10) << lcolltbl[p_i*LT_Grids.TCollPtsLength() + T_i] << endl;
		}
	}

	file_out_lcoll.close(); //closing LColl file
}