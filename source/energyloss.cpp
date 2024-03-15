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
	std::map<std::string, std::string> inputParamsFile;
	if (inputparams.count("c") > 0) {
		if (loadInputsFromFile(inputparams.at("c"), inputParamsFile) != 1) {
			m_error = true;
		}
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
	m_collsys = "PbPb"; if (inputParamsFile.count("collsys") > 0) m_collsys = inputParamsFile["collsys"];
						if (  inputparams.count("collsys") > 0) m_collsys =   inputparams["collsys"];
	
	m_sNN = "5020GeV"; if (inputParamsFile.count("sNN") > 0) m_sNN = inputParamsFile["sNN"];
					   if (  inputparams.count("sNN") > 0) m_sNN =   inputparams["sNN"];

	m_pName = "Charm"; if (inputParamsFile.count("pName") > 0) m_pName = inputParamsFile["pName"];
					   if (  inputparams.count("pName") > 0) m_pName =   inputparams["pName"];

	m_centrality = "30-40%"; if (inputParamsFile.count("centrality") > 0) m_centrality = inputParamsFile["centrality"];
						     if (  inputparams.count("centrality") > 0) m_centrality =   inputparams["centrality"];

	m_xB = 0.6; if (inputParamsFile.count("xB") > 0) m_xB = stod(inputParamsFile["xB"]);
				if (  inputparams.count("xB") > 0) m_xB = stod(  inputparams["xB"]);

    m_xGridN = 25; if (inputParamsFile.count("xGridN") > 0) m_xGridN = stoi(inputParamsFile["xGridN"]);
				   if (  inputparams.count("xGridN") > 0) m_xGridN = stoi(  inputparams["xGridN"]);
    
    m_yGridN = 25; if (inputParamsFile.count("yGridN") > 0) m_yGridN = stoi(inputParamsFile["yGridN"]);
				   if (  inputparams.count("yGridN") > 0) m_yGridN = stoi(  inputparams["yGridN"]);

	m_phiGridN = 25; if (inputParamsFile.count("phiGridN") > 0) m_phiGridN = stoi(inputParamsFile["phiGridN"]);
					 if (  inputparams.count("phiGridN") > 0) m_phiGridN = stoi(  inputparams["phiGridN"]);

	m_TIMESTEP = 0.1; if (inputParamsFile.count("TIMESTEP") > 0) m_TIMESTEP = stod(inputParamsFile["TIMESTEP"]);
					  if (  inputparams.count("TIMESTEP") > 0) m_TIMESTEP = stod(  inputparams["TIMESTEP"]);

	m_TCRIT = 0.155; if (inputParamsFile.count("TCRIT") > 0) m_TCRIT = stod(inputParamsFile["TCRIT"]);
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

int energyLoss::loadInputsFromFile(const std::string &filePath, std::map<std::string, std::string> &inputParamsFile)
{
	std::ifstream file_in(filePath);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
		return -1;
	}
	std::string line, key, sep, val;
	while (std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> key; ss >> sep; ss >> val;
		inputParamsFile[key] = val;
	}
	file_in.close();
	return 1;
}

energyLoss::~energyLoss() {}

void energyLoss::runEnergyLoss()
{
    if (m_error) return;

	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

	if (loadLdndx() != 1) return;
	if (loadLNorm() != 1) return;
	if (loadLColl() != 1) return;

    if (generateInitPosPoints() != 1) return;
    if (loadTempEvol()          != 1) return;

	if ((m_pName == "Bottom") || (m_pName == "Charm")) {
		runELossHeavyFlavour();
	}
	else if (m_pName == "LQuarks") {
		runELossLightQuarks();
	}
	else {
		runELossLightFlavour();
	}
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
			if (line.at(0) == '#')
				continue;
			
			std::stringstream ss(line);
			ss >> buffer; tempTau.push_back(buffer);
			ss >> buffer; tempX.push_back(buffer);
			ss >> buffer; tempY.push_back(buffer);
			ss >> buffer; tempT.push_back(buffer);
		}
	}
	else if (columnCnt == 5) { //evolution file has 5 columns (energy density and temperature)

		while (std::getline(file_in, line)) {
			if (line.at(0) == '#')
				continue;
			
			std::stringstream ss(line);
			ss >> buffer; tempTau.push_back(buffer);
			ss >> buffer; tempX.push_back(buffer);
			ss >> buffer; tempY.push_back(buffer);
			ss >> buffer; tempDataA.push_back(buffer);
			ss >> buffer; tempDataB.push_back(buffer);
		}

		double tempDataAMax = *std::max_element(tempDataA.begin(), tempDataA.end());
		double tempDataBMax = *std::max_element(tempDataB.begin(), tempDataB.end());

        if (tempDataAMax < tempDataBMax) {
            tempT.assign(tempDataA.begin(), tempDataA.end()); // 4th column is temperature
		} else {
            tempT.assign(tempDataB.begin(), tempDataB.end()); // 5th column is temperature
		}

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
		if (line.at(0) == '#')
			continue;
		
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
	m_xGridPts.resize(0); m_yGridPts.resize(0); m_phiGridPts.resize(0);

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
			for (long iY=0; iY<=m_yGridN; iY++) {
                m_xGridPts.push_back(-1.0*initGridRange + 2.0*iX*initGridRange/static_cast<double>(m_xGridN));
                m_yGridPts.push_back(-1.0*initGridRange + 2.0*iY*initGridRange/static_cast<double>(m_yGridN));
            }
        }

		if (loadPhiPoints(m_phiGridPts) != 1) return -7;
	}

	return 1;
}


void energyLoss::radCollEnergyLoss(double x, double y, double phi, std::vector<double> &radRAA1, std::vector<std::vector<double>> &radRAA2, std::vector<double> &collEL, double &pathLength, double &temperature) const
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (modefied pT integration algorithm)
//x, y, phi     - inital position and angle 					  		     <- input
//radiativeRAA1 - radiative RAA for single trajectory (dA410)	  		     <- output
//radiativeRAA2 - radiative RAA for single trajectory (rest of dA integrals) <- output
//collisionalEL - collisional energy loss for single trajectory   		     <- output
//pathLength    - path-length for single trajectory							 <- output
//temperature   - temperature for single trajectory							 <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = m_tau0, currTemp; // defining current path-length (time) and temperature

	while ((currTemp = m_tempEvol.interpolation(t, x + t*std::cos(phi), y + t*std::sin(phi))) > m_TCRIT) {// calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += m_TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) {// calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); // LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											  // table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); // Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	  // table for currDndxInterp

		for (const auto &p : m_Grids.pPts()) //loop over ppts
		{
			for (size_t l=0; l<currLTTabL.size(); l++) {// loop over current path-length and temperature table
				currNormTabTau[l] = currLTTabL[l];											//setting path-lengths
				currNormTabVal[l] = m_LNorm.interpolation(currLTTabL[l], p, currLTTabT[l]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);													  //setting p of current norm table
			NormSparseV.push_back(poly::linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : m_Grids.xPts()) {// loop over xpts
				for (size_t l=0; l<currLTTabL.size(); l++) {// loop over current path-length and temperature table
					currDndxTabTau[l] = currLTTabL[l];											   //setting path-lengths
					currDndxTabVal[l] = m_Ldndx.interpolation(currLTTabL[l], p, currLTTabT[l], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p);													  //setting p of current dndx table
				dndxSparseX.push_back(x);													  //setting x of current dndx table
				dndxSparseV.push_back(poly::linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF<double> currNorm(NormSparseP, NormSparseV);				//constructing interpolated current norm
		interpolationF<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		

		for (const auto &ph : m_Grids.RadPts()) {// loop over Radpts
			radRAA1.push_back(dAp410(ph, currNorm));

			radRAA2.push_back(std::vector<double>());
			for (const auto &Fdp : m_Grids.FdpPts())
				radRAA2.back().push_back(FdA(ph, Fdp, currNorm, currDndx));

		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto &p : m_Grids.pCollPts()) {// loop over pCollPts
			for (size_t l=0; l<currLTTabL.size(); l++) {// loop over current path-length and temperature table
				currCollTabTau[l] = currLTTabL[l];							 //setting path-lengths
				currCollTabVal[l] = m_LColl.interpolation(p, currLTTabT[l]); //setting LColl values
			}

			collEL.push_back(poly::linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLength = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temperature = 0.0;
		for (size_t l=0; l<currLTTabL.size(); l++) temperature += currLTTabT[l];
		temperature /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLength  = 0.0; //setting path-length and temperature
		temperature = 0.0;
	}
}

void energyLoss::radCollEnergyLoss(double x, double y, double phi, std::vector<double> &radRAA, std::vector<double> &collEL, double &pathLenght, double &temperature) const
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (standard algorithm)
//x, y, phi   - inital position and angle 					  <- input
//radRAA      - radiative RAA for single trajectory 		  <- output
//collEL      - collisional energy loss for single trajectory <- output
//pathLenght  - path-length for single trajectory			  <- output
//temperature - temperature for single trajectory 			  <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = m_tau0, currTemp; //defining current path-length (time) and temperature

	while ((currTemp = m_tempEvol.interpolation(t, x + t*std::cos(phi), y + t*std::sin(phi))) > m_TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += m_TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											  //table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	  //table for currDndxInterp

		for (const auto &p : m_Grids.pPts()) //loop over ppts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currNormTabTau[iL] = currLTTabL[iL];										   //setting path-lengths
				currNormTabVal[iL] = m_LNorm.interpolation(currLTTabL[iL], p, currLTTabT[iL]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);													  //setting p of current norm table
			NormSparseV.push_back(poly::linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : m_Grids.xPts()) //loop over xpts
			{
				for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
				{
					currDndxTabTau[iL] = currLTTabL[iL]; 									          //setting path-lengths
					currDndxTabVal[iL] = m_Ldndx.interpolation(currLTTabL[iL], p, currLTTabT[iL], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 													  //setting p of current dndx table
				dndxSparseX.push_back(x);													  //setting x of current dndx table
				dndxSparseV.push_back(poly::linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF<double> currNorm(NormSparseP, NormSparseV); 			    //constructing interpolated current norm
		interpolationF<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (const auto &p : m_Grids.RadPts())
			radRAA.push_back(dA41(p, currNorm, currDndx)/m_dsdpti2.interpolation(p)); //calculating radiative RAA
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto &p : m_Grids.pCollPts()) //loop over pCollPts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currCollTabTau[iL] = currLTTabL[iL]; 				            //setting path-lengths
				currCollTabVal[iL] = m_LColl.interpolation(p, currLTTabT[iL]); //setting LColl values
			}

			collEL.push_back(poly::linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLenght = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temperature = 0.0;
		for (size_t iL=0; iL<currLTTabL.size(); iL++) temperature += currLTTabT[iL];
		temperature /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLenght   = 0.0; //setting path-length and temperature
		temperature  = 0.0;
	}
}


void energyLoss::generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab) const
//function that generates sampling points for Gaussian integration
//qGTab, fGTab - vectors that store sampling point <- output
{	
	double sigmaNum = 3.5; //setting sigma
	double sigmaStep = 0.25; //setting step
	size_t GTabLen = 2 * static_cast<size_t>(sigmaNum / sigmaStep) + 1; //setting length of sampling points
	
	double GaussTabSum = 0.0; //setting normalization sum to zero
	
	for (size_t iG=0; iG<GTabLen; iG++) //calculating sampling points
	{
		qGTab.push_back(-1.0*sigmaNum + static_cast<double>(iG)*sigmaStep); //setting qGaussTab values
		fGTab.push_back(std::exp(-qGTab.back()*qGTab.back()/2.0));          //setting fGaussTab values
		GaussTabSum += fGTab.back();                                        //adding to normalization sum
	}
	
	for (size_t iG=0; iG<GTabLen; iG++)  //normalizing
	{
		fGTab[iG] /= GaussTabSum; //dividing fGaussTab values with total sum
	}
}

void energyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (m_dsdpti2.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();            //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (m_dsdpti2.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
			}
		}
	}
}

void energyLoss::gaussFilterIntegrate(const interpolationF<double> &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
//dsdpti2lquark - light quark initial pT distribution      						  <- input
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (dsdpti2lquark.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (dsdpti2lquark.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
			}
		}
	}
}

void energyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA) const
//function that performs Gauss filter integration - default algorithm
//radiativeRAA  - raditive RAA 							   <- input
//collisionalEL - collisional energy loss				   <- input
//singRAA 		- RAA array after Gauss filter integration <- output
{
    interpolationF<double> RadRelInt(m_Grids.RadPts(),   radiativeRAA);  //creating radiative RAA interpolated function
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	double GFSum; //defining sum variable for Gauss filter

	double dpT; //defining pT and dpT variables

	double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value

	double sigmaColl; //defining variable for collisional sigma
	
	//Gauss filter
	for (const auto &pT : m_Grids.finPts())
	{
		GFSum = 0.0L;

		muCollCurrVal = muCollInt.interpolation(pT);

		sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

		qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

		if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
			double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}		
		
		if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
			double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}
		
		//calculating Gauss filter
		for (size_t iG=0; iG<qGaussTab.size(); iG++)
		{
			dpT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
			GFSum += (m_dsdpti2.interpolation(pT + dpT)*RadRelInt.interpolation(pT + dpT)*(pT + dpT) / pT * fGaussTab[iG]);
		}

		singRAA.push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
	}
}


void energyLoss::calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp) const
{
	interpolationF<double> pathLenghDistInt(m_phiGridPts, pathLenghDist);
	avgPathLength[0] = poly::cubicIntegrate(m_phiGridPts, pathLenghDist)/2.0/M_PI;
	avgPathLength[1] = (pathLenghDistInt.interpolation(m_phiGridPts.front()) + pathLenghDistInt.interpolation(m_phiGridPts.back()))/2.0;
	avgPathLength[2] = (pathLenghDistInt.interpolation(M_PI/2.0)             + pathLenghDistInt.interpolation(3.0*M_PI/2.0))       /2.0;

	interpolationF<double> temperatureDistInt(m_phiGridPts, temperatureDist);
	avgTemp[0] = poly::cubicIntegrate(m_phiGridPts, temperatureDist)/2.0/M_PI;
	avgTemp[1] = (temperatureDistInt.interpolation(m_phiGridPts.front()) + temperatureDistInt.interpolation(m_phiGridPts.back()))/2.0;
	avgTemp[2] = (temperatureDistInt.interpolation(M_PI/2.0)             + temperatureDistInt.interpolation(3.0*M_PI/2.0))       /2.0;
}

int energyLoss::exportResults(const std::string &pName, const std::vector<std::vector<double>> &RAADist, const std::vector<double> avgPathLength, const std::vector<double> avgTemp)
{
	std::vector<std::string> header;
	header.push_back("#collision_system: " + m_collsys);
	header.push_back("#collision_energy: " + m_sNN);
	header.push_back("#particle_type: " + pName);
	header.push_back("#centrality: " + m_centrality);

	std::stringstream xbSStr; xbSStr << std::fixed << std::setprecision(1) << m_xB;
	header.push_back("#xB = " + xbSStr.str());

	std::stringstream avgPathLengthSStr[3];
    for (size_t i=0; i<3; i++) avgPathLengthSStr[i] << std::fixed << std::setprecision(6) << avgPathLength[i];
	header.push_back("#average_path-lengths: " + avgPathLengthSStr[0].str() + ", " + avgPathLengthSStr[1].str() + ", " + avgPathLengthSStr[2].str());

	std::stringstream avgTempSStr[3];
    for (size_t i=0; i<3; i++) avgTempSStr[i] << std::fixed << std::setprecision(6) << avgTemp[i];
	header.push_back("#average_temperatures: " + avgTempSStr[0].str() + ", " + avgTempSStr[1].str() + ", " + avgTempSStr[2].str());
	
	if (m_yGridN <= 0) {
		header.push_back("#number_of_angles:    " + std::to_string(m_phiGridN));
		header.push_back("#number_of_xy_points: " + std::to_string(m_xGridN));
	}
	else {
		header.push_back("#number_of_angles:      " + std::to_string(m_phiGridN));
		header.push_back("#number_of_grid_points: " + std::to_string(m_xGridN) + ", " + std::to_string(m_yGridN));
	}

	header.push_back("#-------------------------------------------------------");
	header.push_back("#   pT [GeV]       phi          R_AA   ");

	const std::string path_out = "./results/results" + pName + "/" + pName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xbSStr.str() + "_dist.dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable to open RAA(pT,phi) distribution file. Aborting..." << std::endl;
		return -1;
	}

	for (const auto &h : header) file_out << h << "\n";

	for (size_t ipT= 0; ipT<m_Grids.finPtsLength(); ipT++)
		for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
			file_out << std::fixed << std::setw(14) << std::setprecision(10) <<   m_Grids.finPts(ipT) << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) <<    m_phiGridPts[iPhi] << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) << RAADist[ipT][iPhi] << "\n";
		}

	file_out.close();

	return 1;
}


void energyLoss::runELossHeavyFlavour()
{
	if (loaddsdpti2() != 1) return;

	FdAHaltonSeqInit(150);

	std::vector<std::vector<double>> RAADist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

	std::vector<double> pathLengthDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

	#pragma omp declare reduction(vectorDoublePlus : std::vector<double> : \
                              	  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    			  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	
	#pragma omp parallel for reduction(vectorDoublePlus : pathLengthDist, temperatureDist) schedule(dynamic)
	for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
		double phi = m_phiGridPts[iPhi];

		std::vector<double> sumRAA1(m_Grids.finPtsLength(), 0.0);

		std::vector<std::vector<double>> sumRAA2(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0));

		double weightsumEL = 0.0, weightsumPLT = 0.0; //energy and path-length and temperature loss weightsum

		for (size_t iXY=0; iXY<m_xGridPts.size(); iXY++) {
			double x = m_xGridPts[iXY], y = m_yGridPts[iXY];
			double binCollDensity = m_binCollDensity.interpolation(x, y);

			if (binCollDensity > 0) {
				weightsumEL += binCollDensity;

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				radCollEnergyLoss(x, y, phi, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					for (auto &cEL : collEL) cEL += 1e-12; //modifying collEL to prevent division by 0

							 weightsumPLT += binCollDensity;
					 pathLengthDist[iPhi] += (pathLength*binCollDensity);
					temperatureDist[iPhi] += (temperature*binCollDensity);

					std::vector<double> singleRAA1; std::vector<std::vector<double>> singleRAA2;
					gaussFilterIntegrate(radRAA1, radRAA2, collEL, singleRAA1, singleRAA2);

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
						sumRAA1[iFinPts] += singleRAA1[iFinPts]*binCollDensity;
						for (size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++) {
							sumRAA2[iFinPts][iFdp] += singleRAA2[iFinPts][iFdp]*binCollDensity;
						}
					}
				}
				else {// if path length is smaller than tau0:						
					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {//multiplying RAA1 (which is 1) with binary collision function as weigth and adding to RAA sum; RAA2 is 0 in this case
						sumRAA1[iFinPts] += binCollDensity;
					}
				}
			}
		}

		std::for_each(sumRAA1.begin(), sumRAA1.end(), [weightsumEL](double &c){c/=weightsumEL;});
		for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
			std::for_each(sumRAA2[iFinPts].begin(), sumRAA2[iFinPts].end(), [weightsumEL](double &c){c/=weightsumEL;});
		}
		for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
			RAADist[iFinPts][iPhi] = sumRAA1[iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iFinPts])/m_Grids.finPts(iFinPts);
		}
		pathLengthDist[iPhi] /= weightsumPLT; temperatureDist[iPhi] /= weightsumPLT;
	}

	std::vector<double> avgPathLength(3, 0.0), avgTemp(3, 0.0);
	calculateAvgPathlenTemps(pathLengthDist, temperatureDist, avgPathLength, avgTemp);

	if (exportResults(m_pName, RAADist, avgPathLength, avgTemp) != 1) return;
}

void energyLoss::runELossLightQuarks()
{
	const std::vector<std::string> lightQuarksList{"Down", "DownBar", "Strange", "Up", "UpBar"};

	std::vector<interpolationF<double>> dsdpti2LightQuarks(lightQuarksList.size());
	for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
		if (loaddsdpti2(lightQuarksList[iLQ], dsdpti2LightQuarks[iLQ]) != 1) return;

	FdAHaltonSeqInit(100);

	std::vector<std::vector<std::vector<double>>> RAADist(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN)));

	std::vector<double> pathLengthDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

	#pragma omp declare reduction(vectorDoublePlus : std::vector<double> : \
                              	  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    			  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	
	#pragma omp parallel for reduction(vectorDoublePlus : pathLengthDist, temperatureDist) schedule(dynamic)
	for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
		double phi = m_phiGridPts[iPhi];

		std::vector<std::vector<double>> sumRAA1(lightQuarksList.size(), std::vector<double>(m_Grids.finPtsLength(), 0.0));		
		std::vector<std::vector<std::vector<double>>> sumRAA2(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0)));
	
		double weightsumEL = 0.0, weightsumPLT = 0.0; //energy and path-length and temperature loss weightsum

		for (size_t iXY=0; iXY<m_xGridPts.size(); iXY++) {
			double x = m_xGridPts[iXY], y = m_yGridPts[iXY];
			double binCollDensity = m_binCollDensity.interpolation(x, y);

			if (binCollDensity > 0) {
				weightsumEL += binCollDensity;

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				radCollEnergyLoss(x, y, phi, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					for (auto &cEL : collEL) cEL += 1e-12; //modifying collEL to prevent division by 0

							 weightsumPLT += binCollDensity;
					 pathLengthDist[iPhi] += (pathLength*binCollDensity);
					temperatureDist[iPhi] += (temperature*binCollDensity);

					std::vector<std::vector<double>> singleRAA1(lightQuarksList.size());
					std::vector<std::vector<std::vector<double>>> singleRAA2(lightQuarksList.size());
					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						gaussFilterIntegrate(dsdpti2LightQuarks[iLQ], radRAA1, radRAA2, collEL, singleRAA1[iLQ], singleRAA2[iLQ]);
					
					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
						for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
							sumRAA1[iLQ][iFinPts] += singleRAA1[iLQ][iFinPts]*binCollDensity;
							for (size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++)
								sumRAA2[iLQ][iFinPts][iFdp] += singleRAA2[iLQ][iFinPts][iFdp]*binCollDensity;
						}
					}
				}
				else {
					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
							sumRAA1[iLQ][iFinPts] += binCollDensity;
				}
			}			
		}

		for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
			std::for_each(sumRAA1[iLQ].begin(), sumRAA1[iLQ].end(), [weightsumEL](double &c){ c/=weightsumEL; });
			for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
				std::for_each(sumRAA2[iLQ][iFinPts].begin(), sumRAA2[iLQ][iFinPts].end(), [weightsumEL](double &c){ c/=weightsumEL; });
			}
		}

		//setting RAA(pT,phi) value by integrating over p:
		for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
			for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
				RAADist[iLQ][iFinPts][iPhi] = sumRAA1[iLQ][iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iLQ][iFinPts])/m_Grids.finPts(iFinPts);
			}
		}

		pathLengthDist[iPhi] /= weightsumPLT; temperatureDist[iPhi] /= weightsumPLT;
	}

	std::vector<double> avgPathLength(3, 0.0), avgTemp(3, 0.0);
	calculateAvgPathlenTemps(pathLengthDist, temperatureDist, avgPathLength, avgTemp);

	for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
		if (exportResults(lightQuarksList[iLQ], RAADist[iLQ], avgPathLength, avgTemp) != 1) return;
	}
}

void energyLoss::runELossLightFlavour()
{
	if (loaddsdpti2() != 1) return;

	dAHaltonSeqInit(1000);

	std::vector<std::vector<double>> RAADist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

	std::vector<double> pathLengthDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

	#pragma omp declare reduction(vectorDoublePlus : std::vector<double> : \
                              	  std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    			  initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	
	#pragma omp parallel for reduction(vectorDoublePlus : pathLengthDist, temperatureDist) schedule(dynamic)
	for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
		double phi = m_phiGridPts[iPhi];

		std::vector<double> sumRAA(m_Grids.finPtsLength(), 0.0);

		double weightsumEL = 0.0, weightsumPLT = 0.0; //energy and path-length and temperature loss weightsum

		for (size_t iXY=0; iXY<m_xGridPts.size(); iXY++) {
			double x = m_xGridPts[iXY], y = m_yGridPts[iXY];
			double binCollDensity = m_binCollDensity.interpolation(x, y);

			if (binCollDensity > 0) {
				weightsumEL += binCollDensity;

				std::vector<double> radRAA; std::vector<double> collEL;
				double pathLength, temperature;
				radCollEnergyLoss(x, y, phi, radRAA, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					for (auto &cEL : collEL) cEL += 1e-12; //modifying collEL to prevent division by 0

							 weightsumPLT += binCollDensity;
					 pathLengthDist[iPhi] += (pathLength*binCollDensity);
					temperatureDist[iPhi] += (temperature*binCollDensity);

					std::vector<double> singleRAA;
					gaussFilterIntegrate(radRAA, collEL, singleRAA);

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
						sumRAA[iFinPts] += singleRAA[iFinPts]*binCollDensity;
					}
				}
				else {// if path length is smaller than tau0:						
					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {//multiplying RAA1 (which is 1) with binary collision function as weigth and adding to RAA sum
						sumRAA[iFinPts] += binCollDensity;
					}
				}
			}
		}

		for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
			RAADist[iFinPts][iPhi] = sumRAA[iFinPts]/weightsumEL;
		}

		pathLengthDist[iPhi] /= weightsumPLT; temperatureDist[iPhi] /= weightsumPLT;
	}

	std::vector<double> avgPathLength(3, 0.0), avgTemp(3, 0.0);
	calculateAvgPathlenTemps(pathLengthDist, temperatureDist, avgPathLength, avgTemp);

	if (exportResults(m_pName, RAADist, avgPathLength, avgTemp) != 1) return;	
}