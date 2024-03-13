#include "energyloss.hpp"
#include "grids.hpp"
#include "linearinterpolation.hpp"
#include "polyintegration.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <tuple>
#include <fstream>
#include <cmath>
#include <limits>
#include <iomanip>

energyLoss::energyLoss(int argc, const char *argv[])
{
	m_error = false;

	std::vector<std::string> inputs; for (int i=2; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: --collsys=PbPb --sNN=5020GeV --pName=Charm --centrality=30-40% --xB=0.6 --xGridN=50 --yGridN=50 --phiGridN=25 --TIMESTEP=0.1 --TCRIT=0.155" << std::endl;
		m_error = true;
	}

	std::map<std::string, std::string> inputparams;
	for (const auto &in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

	//checking if configuration file is provided:
	std::map<std::string, std::string> inputparams_f;
	if (inputparams.count("c") > 0) {
		std::ifstream file_in(inputparams["c"]);
		if (!file_in.is_open()) {
			std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
			m_error = true;
		}
		std::string line, key, sep, val;
		while (std::getline(file_in, line))
		{
			std::stringstream ss(line);
			ss >> key; ss >> sep; ss >> val;
			inputparams_f[key] = val;
		}
		file_in.close();
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
	m_collsys = "PbPb"; if (inputparams_f.count("collsys") > 0) m_collsys = inputparams_f["collsys"];
						if (  inputparams.count("collsys") > 0) m_collsys =   inputparams["collsys"];
	
	m_sNN = "5020GeV"; if (inputparams_f.count("sNN") > 0) m_sNN = inputparams_f["sNN"];
					   if (  inputparams.count("sNN") > 0) m_sNN =   inputparams["sNN"];

	m_pName = "Charm"; if (inputparams_f.count("pName") > 0) m_pName = inputparams_f["pName"];
					   if (  inputparams.count("pName") > 0) m_pName =   inputparams["pName"];

	m_centrality = "30-40%"; if (inputparams_f.count("centrality") > 0) m_centrality = inputparams_f["centrality"];
						     if (  inputparams.count("centrality") > 0) m_centrality =   inputparams["centrality"];

	m_xB = 0.6; if (inputparams_f.count("xB") > 0) m_xB = stod(inputparams_f["xB"]);
				if (  inputparams.count("xB") > 0) m_xB = stod(  inputparams["xB"]);

    m_xGridN = 25; if (inputparams_f.count("xGridN") > 0) m_xGridN = stoi(inputparams_f["xGridN"]);
				   if (  inputparams.count("xGridN") > 0) m_xGridN = stoi(  inputparams["xGridN"]);
    
    m_yGridN = 25; if (inputparams_f.count("yGridN") > 0) m_yGridN = stoi(inputparams_f["yGridN"]);
				   if (  inputparams.count("yGridN") > 0) m_yGridN = stoi(  inputparams["yGridN"]);

	m_phiGridN = 25; if (inputparams_f.count("phiGridN") > 0) m_phiGridN = stoi(inputparams_f["phiGridN"]);
					 if (  inputparams.count("phiGridN") > 0) m_phiGridN = stoi(  inputparams["phiGridN"]);

	m_TIMESTEP = 0.1; if (inputparams_f.count("TIMESTEP") > 0) m_TIMESTEP = stod(inputparams_f["TIMESTEP"]);
					  if (  inputparams.count("TIMESTEP") > 0) m_TIMESTEP = stod(  inputparams["TIMESTEP"]);

	m_TCRIT = 0.155; if (inputparams_f.count("TCRIT") > 0) m_TCRIT = stod(inputparams_f["TCRIT"]);
					 if (  inputparams.count("TCRIT") > 0) m_TCRIT = stod(  inputparams["TCRIT"]);

	//checking if provided value of sNN is an option:
	if ((m_sNN != "5440GeV") && (m_sNN != "5020GeV") && (m_sNN != "2760GeV") && (m_sNN != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		m_error = true;
	}

	m_nf = m_sNN == "200GeV" ? 2.5 : 3.0;
	double T = 3.0 / 2.0*m_TCRIT;
	double mu = 0.197*std::sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
	m_mgC = mu / std::sqrt(2.0);
	if (m_pName == "Bottom") m_MC = 4.75;
	else if (m_pName == "Charm") m_MC = 1.2;
	else if (m_pName == "Gluon") m_MC = mu/std::sqrt(2.0);
	else m_MC = mu/sqrt(6.0);
	m_TCollConst = T;
}

energyLoss::~energyLoss() {}

void energyLoss::runEnergyLoss()
{
    if (m_error) return;

	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

	if (loadLdndx() != 1) return;
	if (loadLNorm() != 1) return;
	if (loadLColl() != 1) return;

    if (loadTempEvol()          != 1) return;
    if (generateInitPosPoints() != 1) return;
}

double energyLoss::productLog(double x) const
{
	if (x == 0.0) {
		return 0.0;
	}

	double w0, w1;
	if (x > 0.0) {
		w0 = std::log(1.2 * x / std::log(2.4 * x / std::log1p(2.4 * x)));
	}
	else {
		double v = 1.4142135623730950488 * std::sqrt(1.0 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (true) {
		double e = std::exp(w0);
		double f = w0 * e - x;
		w1 = w0 - f / ((e * (w0 + 1.0) - (w0 + 2.0) * f / (w0 + w0 + 2.0)));
		if (std::abs(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}
	return w1;
}

int energyLoss::loaddsdpti2()
{
	const std::string path_in = "./ptDists/ptDist" + m_sNN + "/ptDist_" + m_sNN + "_" + m_pName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file. Aborting..." << std::endl;
		return -1;
	}

	std::vector<double> pTdistX, pTdistF; //defining vectors that store dsdpti2 values

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; pTdistX.push_back(buffer);
		ss >> buffer; pTdistF.push_back(buffer);
	}

	m_dsdpti2.setData(pTdistX, pTdistF);

	file_in.close();

	return 1;
}

int energyLoss::loaddsdpti2(const std::string &pname, interpolationF<double> &dsdpti2int)const 
{
	const std::string path_in = "./ptDists/ptDist" + m_sNN + "/ptDist_" + m_sNN + "_" + pname + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file. Aborting..." << std::endl;
		return -1;
	}

	std::vector<double> pTdistX, pTdistF;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; pTdistX.push_back(buffer);
		ss >> buffer; pTdistF.push_back(buffer);
	}

	dsdpti2int.setData(pTdistX, pTdistF);

	file_in.close();

	return 1;
}

int energyLoss::loadLdndx()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/ldndx_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open Ldndx table file. Aborting..." << std::endl;
		return -1;
	}

	std::vector<double> Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; Ldndx_tau.push_back(buffer);
		ss >> buffer; Ldndx_p.push_back(buffer);
		ss >> buffer; Ldndx_T.push_back(buffer);
		ss >> buffer; Ldndx_x.push_back(buffer);
		ss >> buffer; Ldndx_f.push_back(buffer);
	}

	file_in.close();

	m_Ldndx.setData(Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f);

	std::vector<std::vector<double>> domain = m_Ldndx.domain();
	if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.xPts(0)    < domain[3][0]) {std::cerr << "Error:   x grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.xPts(-1)   > domain[3][1]) {std::cerr << "Error:   x grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::loadLNorm()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lnorm_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LNorm table file. Aborting..." << std::endl;
		return -1;
	}

	std::vector<double> LNorm_tau, LNorm_p, LNorm_T, LNorm_f; //defining vectors that store LNorm table values

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; LNorm_tau.push_back(buffer);
		ss >> buffer; LNorm_p.push_back(buffer);
		ss >> buffer; LNorm_T.push_back(buffer);
		ss >> buffer; LNorm_f.push_back(buffer);
	}

	file_in.close();

	m_LNorm.setData(LNorm_tau, LNorm_p, LNorm_T, LNorm_f);

	std::vector<std::vector<double>> domain = m_LNorm.domain();
	if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::loadLColl()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lcoll_nf=" + nfss.str() + "_" + partName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LColl table file. Aborting..." << std::endl;
		return -1;
	}

	std::vector<double> LColl_p, LColl_T, LColl_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;
            
		std::stringstream ss(line);
		ss >> buffer; LColl_p.push_back(buffer);
		ss >> buffer; LColl_T.push_back(buffer);
		ss >> buffer; LColl_f.push_back(buffer);
	}

	file_in.close();

	m_LColl.setData(LColl_p, LColl_T, LColl_f);

	std::vector<std::vector<double>> domain = m_LColl.domain();
	if (m_Grids.pCollPts(0)  < domain[0][0]) {std::cerr << "Error: p grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pCollPts(-1) > domain[0][1]) {std::cerr << "Error: p grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TCollPts(0)  < domain[1][0]) {std::cerr << "Error: T grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TCollPts(-1) > domain[1][1]) {std::cerr << "Error: T grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::loadTempEvol()
{
	std::string path_in = "./evols/tempevol_cent=" + m_centrality + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open temperature evolution file. Aborting..." << std::endl;
		return -1;
	}

	std::string line; double buffer;

    while (std::getline(file_in, line)) { // skiping header lines that start with '#'
        if (line.at(0) == '#')
            continue;
        break;
    }

	//checking how many columns evolution file has:
	size_t columnCnt = 0;
	std::stringstream lineSStr(line); while (lineSStr >> buffer) columnCnt++;

	file_in.clear(); file_in.seekg(0); //return to the begining of file

	std::vector<double> tempTau, tempX, tempY, tempDataA, tempDataB, tempT;

	if (columnCnt == 4) { //evolution file has 4 columns (just temperature)

		while (std::getline(file_in, line)) {
			std::stringstream ss(line);
			ss >> buffer; tempTau.push_back(buffer);
			ss >> buffer; tempX.push_back(buffer);
			ss >> buffer; tempY.push_back(buffer);
			ss >> buffer; tempT.push_back(buffer);
		}
	}
	else if (columnCnt == 5) { //evolution file has 5 columns (energy density and temperature)

		while (std::getline(file_in, line)) {
			std::stringstream ss(line);
			ss >> buffer; tempTau.push_back(buffer);
			ss >> buffer; tempX.push_back(buffer);
			ss >> buffer; tempY.push_back(buffer);
			ss >> buffer; tempDataA.push_back(buffer);
			ss >> buffer; tempDataB.push_back(buffer);
		}

		double tempDataAMax = *std::max_element(tempDataA.begin(), tempDataA.end());
		double tempDataBMax = *std::max_element(tempDataB.begin(), tempDataB.end());

        if (tempDataAMax < tempDataBMax)
            tempT.assign(tempDataA.begin(), tempDataA.end()); // 4th column is temperature
        else
            tempT.assign(tempDataB.begin(), tempDataB.end()); // 5th column is temperature

	}
	else { //evolution file is not suitable for interpolation

		std::cerr << "Error: number of columns is not appropriate for temperature evolution interpolation. Aborting..." << std::endl;
		return -2;
	}

    file_in.close();

    double tempXMin = *std::min_element(tempX.begin(), tempX.end());
    double tempYMin = *std::min_element(tempY.begin(), tempY.end());

	if ((tempXMin >= 0.0) && (tempYMin >= 0.0)) {// if temperature evolution is defined only in the first quadrant

		//creating tau grid:
		std::vector<double> tempTauGrid(tempTau.begin(), tempTau.end());
		std::sort(tempTauGrid.begin(), tempTauGrid.end());
		tempTauGrid.erase(unique(tempTauGrid.begin(), tempTauGrid.end()), tempTauGrid.end());

		//creating full x grid:
		std::vector<double> tempXGrid(tempX.begin(), tempX.end());
        size_t sizeX = tempXGrid.size();
        tempXGrid.reserve(sizeX * 2);
        for (size_t i=0; i<sizeX; ++i)
            tempXGrid.push_back(-1.0*tempXGrid[i]);
		sort(tempXGrid.begin(), tempXGrid.end());
		tempXGrid.erase(unique(tempXGrid.begin(), tempXGrid.end()), tempXGrid.end());

		//creating full y grid:
		std::vector<double> tempYGrid(tempY.begin(), tempY.end());
        size_t sizeY = tempYGrid.size();
        tempYGrid.reserve(sizeY * 2);
        for (size_t i=0; i<sizeY; ++i)
            tempYGrid.push_back(-1.0*tempYGrid[i]);
		sort(tempYGrid.begin(), tempYGrid.end());
		tempYGrid.erase(unique(tempYGrid.begin(), tempYGrid.end()), tempYGrid.end());

        // temperature volution interpolated function in first quadtrant:
		interpolationF<double> tempEvolFirstQuadrant(tempTau, tempX, tempY, tempT);

		//creating full temperature evolution table:
		std::vector<double> tempTauFull, tempXFull, tempYFull, tempTFull;
		for (const auto &tau : tempTauGrid) {
			for (const auto &x : tempXGrid) {
				for (const auto &y : tempYGrid) {
					tempTauFull.push_back(tau);
					  tempXFull.push_back(x);
					  tempYFull.push_back(y);
					  tempTFull.push_back(tempEvolFirstQuadrant.interpolation(tau, std::abs(x), std::abs(y)));
				}
            }
        }

		m_tempEvol.setData(tempTauFull, tempXFull, tempYFull, tempTFull);
	}
	else {// if not, creating interpolated function with values from file:

		m_tempEvol.setData(tempTau, tempX, tempY, tempT);
	}

	m_tau0 = m_tempEvol.domain()[0][0];

	return 1;
}

int energyLoss::loadBinCollDensity(interpolationF<double> &binCollDensity)
{
	std::string path_in = "binarycolldensities/binarybolldensity_cent=" + m_centrality + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open binary collision density file." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	std::vector<double> bcdX, bcdY, bcdData;

	while (std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> buffer;    bcdX.push_back(buffer);
		ss >> buffer;    bcdY.push_back(buffer);
		ss >> buffer; bcdData.push_back(buffer);
	}

	file_in.close();

    double bcdXMin = *std::min_element(bcdX.begin(), bcdX.end());
    double bcdYMin = *std::min_element(bcdY.begin(), bcdY.end());

	if ((bcdXMin >= 0.0) && (bcdYMin >= 0.0)) {// if binary collision density is defined only in the first quadrant:

		//creating full x grid:
		std::vector<double> bcdXGrid(bcdX.begin(), bcdX.end());
        size_t sizeX = bcdXGrid.size();
        bcdXGrid.reserve(sizeX * 2);
        for (size_t i=0; i<sizeX; ++i)
            bcdXGrid.push_back(-1.0*bcdXGrid[i]);
		sort(bcdXGrid.begin(), bcdXGrid.end());
		bcdXGrid.erase(unique(bcdXGrid.begin(), bcdXGrid.end()), bcdXGrid.end());

		//creating full y grid:
		std::vector<double> bcdYGrid(bcdY.begin(), bcdY.end());
        size_t sizeY = bcdYGrid.size();
        bcdYGrid.reserve(sizeY * 2);
        for (size_t i=0; i<sizeY; ++i)
            bcdYGrid.push_back(-1.0*bcdYGrid[i]);
		sort(bcdYGrid.begin(), bcdYGrid.end());
		bcdYGrid.erase(unique(bcdYGrid.begin(), bcdYGrid.end()), bcdYGrid.end());

        // creating interpolated binary collision density defined in first quadrant:
		interpolationF<double> binCollDensityFirstQuadrant(bcdX, bcdY, bcdData);

		//creating full binary collision density table:
		std::vector<double> bcdXFull, bcdYFull, bcdDataFull;
		for (const auto &x : bcdXGrid) {
			for (const auto &y : bcdYGrid) {
				bcdXFull.push_back(x);
				bcdYFull.push_back(y);
				bcdDataFull.push_back(binCollDensityFirstQuadrant.interpolation(std::abs(x), std::abs(y)));
			}
        }

		binCollDensity.setData(bcdXFull, bcdYFull, bcdDataFull);
	}
	else { // if not, creating interpolated function with values from file:

		binCollDensity.setData(bcdX, bcdY, bcdData);
	}

	return 1;
}

int energyLoss::loadPhiPoints(std::vector<double> &phiPoints)
{
	std::string path_in = "./phiGaussPts/phiptsgauss" + std::to_string(m_phiGridN) + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open phi points file. Aborting..." << std::endl;
		return -1;
	}

	phiPoints.resize(0);

	std::string line; double buffer;

	while (std::getline(file_in, line)) {
		if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; phiPoints.push_back(buffer);
	}

	file_in.close();

	return 1;
}

int energyLoss::loadBinCollPoints(std::vector<std::vector<double>> &bcPoints)
{
	std::string path_in = "binarycollpoints/binarycollpoints_cent=" + m_centrality + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open binary collision points file." << std::endl;
		return -1;
	}

	bcPoints.resize(0);

	std::string line; double bufferX, bufferY;

	while (std::getline(file_in, line)) {
		std::stringstream ss(line);
		ss >> bufferX;
		ss >> bufferY;
        bcPoints.push_back({bufferX, bufferY});
	}

	file_in.close();

	return 1;
}

int energyLoss::generateInitPosPoints()
{
	if (m_yGridN == -2) {
		//if yGridN is set to -2, MonteCarlo method is used to generate initial position points and angles
		//number of x-y initial position points is equal to m_xGridN

		interpolationF<double> binCollDensity; if (loadBinCollDensity(binCollDensity) != 1) return -1;
		
		std::vector<std::vector<double>> bcDensDomain = binCollDensity.domain();
		std::vector<double> bcDensCoDomain = binCollDensity.codomain();

		//generating x and y points:
		std::random_device rdX; std::mt19937 mtX(rdX());
        std::uniform_real_distribution<double> distX(bcDensDomain[0][0], std::nextafter(bcDensDomain[0][1], std::numeric_limits<double>::max()));
        std::random_device rdY; std::mt19937 mtY(rdY());
        std::uniform_real_distribution<double> distY(bcDensDomain[1][0], std::nextafter(bcDensDomain[1][1], std::numeric_limits<double>::max()));
        std::random_device rdZ; std::mt19937 mtZ(rdZ());
        std::uniform_real_distribution<double> distZ(bcDensCoDomain[0],  std::nextafter(bcDensCoDomain[1],  std::numeric_limits<double>::max()));
	
		double x, y, z;

		for (size_t iXY=0; iXY<m_xGridN; iXY++) {
			do {
				x = distX(mtX);
				y = distY(mtY);
				z = distZ(mtZ);
			}
			while (z > binCollDensity.interpolation(x, y));

			m_xGridPts.push_back(x); m_yGridPts.push_back(y);
		}

		std::random_device rdPhi; std::mt19937 mtPhi(rdPhi());
        std::uniform_real_distribution<double> distPhi(0.0, std::nextafter(2.0*M_PI, std::numeric_limits<double>::max()));

		double phi;

        // artificially adding 0 and 2Pi to the list:
		phi = 0.0; 		m_phiGridPts.push_back(phi); std::sort(m_phiGridPts.begin(), m_phiGridPts.end());
		phi = 2.0*M_PI; m_phiGridPts.push_back(phi); std::sort(m_phiGridPts.begin(), m_phiGridPts.end());

		//generating other points:
		for (size_t iPhi=2; iPhi<m_phiGridN; iPhi++) {
			do {
                phi = distPhi(mtPhi);
            }
            while (std::binary_search(m_phiGridPts.begin(), m_phiGridPts.end(), phi));
			m_phiGridPts.push_back(phi); std::sort(m_phiGridPts.begin(), m_phiGridPts.end());
		}

		//generating binary collision density with function values 1:
		std::vector<double> bcdX, bcdY, bcdData;
		for (size_t iX=0; iX<=10; iX++){
			for (size_t iY=0; iY<=10; iY++){
				   bcdX.push_back(bcDensDomain[0][0] + (bcDensDomain[0][1]-bcDensDomain[0][0])*static_cast<double>(iX)/10.0);
				   bcdY.push_back(bcDensDomain[1][0] + (bcDensDomain[1][1]-bcDensDomain[1][0])*static_cast<double>(iY)/10.0);
				bcdData.push_back(1.0);
			}
		}
		m_binCollDensity.setData(bcdX, bcdY, bcdData);

	}
	else if (m_yGridN == -1) {
		//if yGridN is set to -1, MonteCarlo method is used to generate initial position points, while angles are on equidistant grid
		//number of x-y initial position points is equal to xGridN

		interpolationF<double> binCollDensity; if (loadBinCollDensity(binCollDensity) != 1) return -1;
		
		std::vector<std::vector<double>> bcDensDomain = binCollDensity.domain();
		std::vector<double> bcDensCoDomain = binCollDensity.codomain();

		//generating x and y points:
		std::random_device rdX; std::mt19937 mtX(rdX());
        std::uniform_real_distribution<double> distX(bcDensDomain[0][0], std::nextafter(bcDensDomain[0][1], std::numeric_limits<double>::max()));
        std::random_device rdY; std::mt19937 mtY(rdY());
        std::uniform_real_distribution<double> distY(bcDensDomain[1][0], std::nextafter(bcDensDomain[1][1], std::numeric_limits<double>::max()));
        std::random_device rdZ; std::mt19937 mtZ(rdZ());
        std::uniform_real_distribution<double> distZ(bcDensCoDomain[0],  std::nextafter(bcDensCoDomain[1],  std::numeric_limits<double>::max()));
	
		double x, y, z;

		for (size_t iXY=0; iXY<m_xGridN; iXY++) {
			do {
				x = distX(mtX);
				y = distY(mtY);
				z = distZ(mtZ);
			}
			while (z > binCollDensity.interpolation(x, y));

			m_xGridPts.push_back(x); m_yGridPts.push_back(y);
		}

		if (loadPhiPoints(m_phiGridPts) != 1) return -3;

		//generating binary collision density with function values 1:
		std::vector<double> bcdX, bcdY, bcdData;
		for (size_t iX=0; iX<=10; iX++){
			for (size_t iY=0; iY<=10; iY++){
				   bcdX.push_back(bcDensDomain[0][0] + (bcDensDomain[0][1]-bcDensDomain[0][0])*static_cast<double>(iX)/10.0);
				   bcdY.push_back(bcDensDomain[1][0] + (bcDensDomain[1][1]-bcDensDomain[1][0])*static_cast<double>(iY)/10.0);
				bcdData.push_back(1.0);
			}
		}
		m_binCollDensity.setData(bcdX, bcdY, bcdData);

	}
	else if (m_yGridN == 0) {
		//if yGridN is set to 0, initial position points are randomly selected from a list, while angles are on equidistant grid
		//number of x-y initial position points is equal to xGridN

		std::vector<std::vector<double>> bcPoints; if (loadBinCollPoints(bcPoints) != 1) return -4;

		if ((m_xGridN == bcPoints.size()) || (m_xGridN == 0)) {// take all points if xGridN is equal to total number of points or 0
			for (size_t iXY=0; iXY<bcPoints.size(); iXY++) {
				m_xGridPts.push_back(bcPoints[iXY][0]);
                m_yGridPts.push_back(bcPoints[iXY][1]);
			}
		}
        else {// randomly select from imported points
            std::random_device rd; auto rng = std::default_random_engine{rd()};
	        std::shuffle(bcPoints.begin(), bcPoints.end(), rng);
            for (size_t iXY=0; iXY<m_xGridN; iXY++) {
				m_xGridPts.push_back(bcPoints[iXY][0]);
                m_yGridPts.push_back(bcPoints[iXY][1]);
			}
        }

		if (loadPhiPoints(m_phiGridPts) != 1) return -5;

        //generating binary collision density with function values 1 in domain [-20, 20]:
		std::vector<double> bcdX, bcdY, bcdData;
		for (size_t iX=0; iX<=10; iX++) {
			for (size_t iY=0; iY<=10; iY++) {
				   bcdX.push_back(-20.0 + 40.0*static_cast<double>(iX)/10.0);
				   bcdY.push_back(-20.0 + 40.0*static_cast<double>(iY)/10.0);
				bcdData.push_back(1.0);
			}
		}

		m_binCollDensity.setData(bcdX, bcdY, bcdData);

	}
	else {
		//if yGridN is larger than 0, initial position points and angles are generated on equidistant grid
		//number of x-y initial position points is equal to (xGridN+1)*(yGridN+1)

		if (loadBinCollDensity(m_binCollDensity) != 1) return -6; //loading binary collision density

		std::vector<std::vector<double>> bcDensDomain = m_binCollDensity.domain();

		double initGridRange = 0.0;
		if (std::abs(bcDensDomain[0][0]) > bcDensDomain[0][1])
            initGridRange = std::abs(bcDensDomain[0][0]) - 0.5;
		else
            initGridRange = std::abs(bcDensDomain[0][1]) - 0.5;

		for (size_t iX=0; iX<=m_xGridN; iX++) {
			for (size_t iY=0; iY<=m_yGridN; iY++) {
                m_xGridPts.push_back(-1.0*initGridRange + 2.0*iX*initGridRange/static_cast<double>(m_xGridN));
                m_yGridPts.push_back(-1.0*initGridRange + 2.0*iY*initGridRange/static_cast<double>(m_yGridN));
            }
        }

		if (loadPhiPoints(m_phiGridPts) != 1) return -7;
	}

	return 1;
}