#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

#include "InputParser.hpp"

int GetInputs(vector<string> inputs, string &pname, string &cent, double &xb, int &xptsn, int &yptsn, int &phiptsn, string &pTinitpath, string &temppath, string &bcdpath, double &timestep, double &tcrit)
{
	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		cout << "default parameter values: -pName=Charm -centrality=30-40% -xB=0.4 -xGridN=40 -yGridN=40 -phiGridN=50 -TIMESTEP=0.1 -TCRIT=0.155 -pTinit_path='' -temp_path='' -bcd_path=''" << endl;
		return 0;
	}

	vector<string> key, value;
	for (auto in : inputs)
	{
		key.push_back(in.substr(0, in.find("=")));
		value.push_back(in.substr(in.find("=")+1, in.length()));
	}

	pname = "Charm";
	for (int i=0; i<key.size(); i++) if (key[i] == "-pName") {pname=value[i]; break;}

	cent = "30-40%";
	for (int i=0; i<key.size(); i++) if (key[i] == "-centrality") {cent=value[i]; break;}

	xb = 0.4;
	for (int i=0; i<key.size(); i++) if (key[i] == "-xB") {stringstream ss(value[i]); ss >> xb; break;}

	xptsn = 40;
	for (int i=0; i<key.size(); i++) if (key[i] == "-xGridN") {stringstream ss(value[i]); ss >> xptsn; break;}

	yptsn = 40;
	for (int i=0; i<key.size(); i++) if (key[i] == "-yGridN") {stringstream ss(value[i]); ss >> yptsn; break;}

	phiptsn = 50;
	for (int i=0; i<key.size(); i++) if (key[i] == "-phiGridN") {stringstream ss(value[i]); ss >> phiptsn; break;}

	pTinitpath = "";
	for (int i=0; i<key.size(); i++) if (key[i] == "-pTinit_path") {pTinitpath=value[i]; break;}

	temppath = "";
	for (int i=0; i<key.size(); i++) if (key[i] == "-temp_path") {temppath=value[i]; break;}

	bcdpath = "";
	for (int i=0; i<key.size(); i++) if (key[i] == "-bcd_path") {bcdpath=value[i]; break;}

	timestep = 0.1;
	for (int i=0; i<key.size(); i++) if (key[i] == "-TIMESTEP") {stringstream ss(value[i]); ss >> timestep; break;}

	tcrit = 0.155;
	for (int i=0; i<key.size(); i++) if (key[i] == "-TCRIT") {stringstream ss(value[i]); ss >> tcrit; break;}

	return 1;
}

int GetInputs(vector<string> inputs, string &pname, double &xb, int &LdndxMaxPts, int &LCollMaxPts)
{
	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		cout << "default values: -pName=Charm -xB=0.4 -LdndxMaxPoints=100000 -LCollMaxPoints=10000" << endl;
		return 0;
	}

	vector<string> key, value;
	for (auto in : inputs)
	{
		key.push_back(in.substr(0, in.find("=")));
		value.push_back(in.substr(in.find("=")+1, in.length()));
	}

	pname = "Charm";
	for (int i=0; i<key.size(); i++) if (key[i] == "-pName") {pname=value[i]; break;}

	xb = 0.4;
	for (int i=0; i<key.size(); i++) if (key[i] == "-xB") {stringstream ss(value[i]); ss >> xb; break;}

	LdndxMaxPts = 500000;
	for (int i=0; i<key.size(); i++) if (key[i] == "-LdndxMaxPoints") {stringstream ss(value[i]); ss >> LdndxMaxPts; break;}

	LCollMaxPts = 10000;
	for (int i=0; i<key.size(); i++) if (key[i] == "-LCollMaxPoints") {stringstream ss(value[i]); ss >> LCollMaxPts; break;}

	return 1;
}