#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
using namespace std;

#include <omp.h>

#include "MainHeader.hpp"
#include "ELossHeader.hpp"
#include "Arsenal.hpp"
#include "LinearInterpolation.hpp"
#include "Grids.hpp"
#include "ImportExport.hpp"
#include "dAIntegrals.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//defining global variables to be used in all source files:

string pName; 					    		   //particle name
string centrality;				    		   //centrality class
double xB;						    		   //xB value
GridPoints Grids;				    		   //grid points
interpFun dsdpti2; 				    		   //initial pT distribution
interpFun LNorm, Ldndx, LColl;	    		   //interpolated L tables
interpFun TProfile, BinCollDensity; 		   //temperature evolution and binary collision density
double tau0;								   //thermalization time
double mgC, MC;					   			   //constant particle and gluon masses used for dA integrals
double TCollConst;				    		   //constant temperature used for Gauss filter integration
int xGridN, yGridN, phiGridN;	    		   //initial position grid points and angle number
vector<double> xGridPts, yGridPts, phiGridPts; //defining vectors that store initial position points and angles
double TIMESTEP, TCRIT; 					   //time step and critical temperature
string pTinit_path, temp_path, bcd_path;	   //initial pT distribution, temperature evolution and binary collision density paths

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//defining global variables for this source file:

static double nf = 3.0;		  //effective number of flavours
static double lambda = 0.2;   //QCD scale

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FdA and dA INTEGRALS FUNCTIONS:

//FdA function (number of FdA integrals depends on particle type) - used for heavy flavour
double FdA(double ph, double dp, interpFun& currNormInt, interpFun& currDndxInt)
{
	return (FdA411(ph, dp, currNormInt, currDndxInt) + FdA412(ph, dp, currNormInt, currDndxInt) + FdA413(ph, dp, currNormInt, currDndxInt) +
				FdA414(ph, dp, currNormInt, currDndxInt) + FdA415(ph, dp, currNormInt, currDndxInt));
}

//dA41 function (number of dA integrals depends on particle type) - used for light flavour
double dA41(double ph, interpFun &currNormInt, interpFun &currDndxInt)
{
	if (pName == "Gluon") { //gluon needs 7 dA integrals
		
		return (dA410(ph, currNormInt) + dA411(ph, currNormInt, currDndxInt) + dA412(ph, currNormInt, currDndxInt) +dA413(ph, currNormInt, currDndxInt) +
					dA414(ph, currNormInt, currDndxInt) + dA415(ph, currNormInt, currDndxInt) + dA416(ph, currNormInt, currDndxInt) +
					dA417(ph, currNormInt, currDndxInt));
	}
	else { //light quarks need 5 dA integrals
		
		return (dA410(ph, currNormInt) + dA411(ph, currNormInt, currDndxInt) + dA412(ph, currNormInt, currDndxInt) + dA413(ph, currNormInt, currDndxInt) +
					dA414(ph, currNormInt, currDndxInt) + dA415(ph, currNormInt, currDndxInt));
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ENERGY LOSS FUNCTION:

//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (modefied pT integration algorithm)
void RadCollEL(double X0, double Y0, double phi0, vector<double> &radiativeRAA1, vector<vector<double>> &radiativeRAA2, vector<double> &collisionalEL, double &pathL, double &temp)
//X0, Y0, phi0  - inital position and angle 					  		     <- input
//radiativeRAA1 - radiative RAA for single trajectory (dA410)	  		     <- output
//radiativeRAA2 - radiative RAA for single trajectory (rest of dA integrals) <- output
//collisionalEL - collisional energy loss for single trajectory   		     <- output
//pathL, temp - path-length and temperature for single trajectory 		     <- output
{
	vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = tau0, curr_temp; //defining current path-length (time) and temperature

	while ((curr_temp = TProfile.interp(t, X0 + t*cos(phi0), Y0 + t*sin(phi0))) > TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(curr_temp);
		t += TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		vector<double> NormSparseP, NormSparseV;											 //table for currNormInterp
		
		vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	 //table for currDndxInterp

		for (auto p : Grids.pPts()) //loop over ppts
		{
			for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
			{
				currNormTabTau[l] = currLTTabL[l]; 								   //setting path-lengths
				currNormTabVal[l] = LNorm.interp(currLTTabL[l], p, currLTTabT[l]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(LinearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (auto x : Grids.xPts()) //loop over xpts
			{
				for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
				{
					currDndxTabTau[l] = currLTTabL[l]; 									  //setting path-lengths
					currDndxTabVal[l] = Ldndx.interp(currLTTabL[l], p, currLTTabT[l], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(LinearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpFun currNorm(NormSparseP, NormSparseV); 			   //constructing interpolated current norm
		interpFun currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		//radiativeRAA2.resize(Grids.RadPtsLength()); //resizing radiative RAA 2d vector

		for (auto ph : Grids.RadPts()) //loop over Radpts
		{
			radiativeRAA1.push_back(dAp410(ph, currNorm)); //calculating radiative energy loss for dA410

			radiativeRAA2.push_back(vector<double>()); //resizing radiativeRAA2 2d vector

			for (auto Fdp : Grids.FdpPts())
				radiativeRAA2.back().push_back(FdA(ph, Fdp, currNorm, currDndx)); //calculating radiative energy loss for rest of the dA integrals

		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (auto p : Grids.pCollPts()) //loop over pCollPts
		{
			for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
			{
				currCollTabTau[l] = currLTTabL[l]; 				    //setting path-lengths
				currCollTabVal[l] = LColl.interp(p, currLTTabT[l]); //setting LColl values
			}

			collisionalEL.push_back(LinearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathL = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (int l=0; l<currLTTabL.size(); l++) temp += currLTTabT[l];
		temp /= currLTTabL.size();
	}
	else { //if path-length is smaller than thermalization time:

		pathL = 0.0; //setting path-length and temperature
		temp  = 0.0;
	}
}

//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (standard algorithm)
void RadCollEL(double X0, double Y0, double phi0, vector<double> &radiativeRAA, vector<double> &collisionalEL, double &pathL, double &temp)
//X0, Y0, phi0  - inital position and angle 					  <- input
//radiativeRAA  - radiative RAA for single trajectory 			  <- output
//collisionalEL - collisional energy loss for single trajectory   <- output
//pathL, temp - path-length and temperature for single trajectory <- output
{
	vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = tau0, curr_temp; //defining current path-length (time) and temperature

	while ((curr_temp = TProfile.interp(t, X0 + t*cos(phi0), Y0 + t*sin(phi0))) > TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(curr_temp);
		t += TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		vector<double> NormSparseP, NormSparseV;											 //table for currNormInterp
		
		vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	 //table for currDndxInterp

		for (auto p : Grids.pPts()) //loop over ppts
		{
			for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
			{
				currNormTabTau[l] = currLTTabL[l]; 								   //setting path-lengths
				currNormTabVal[l] = LNorm.interp(currLTTabL[l], p, currLTTabT[l]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(LinearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (auto x : Grids.xPts()) //loop over xpts
			{
				for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
				{
					currDndxTabTau[l] = currLTTabL[l]; 									  //setting path-lengths
					currDndxTabVal[l] = Ldndx.interp(currLTTabL[l], p, currLTTabT[l], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(LinearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpFun currNorm(NormSparseP, NormSparseV); 			   //constructing interpolated current norm
		interpFun currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (auto p : Grids.RadPts())
			radiativeRAA.push_back(dA41(p, currNorm, currDndx)/dsdpti2.interp(p)); //calculating radiative RAA
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (auto p : Grids.pCollPts()) //loop over pCollPts
		{
			for (int l=0; l<currLTTabL.size(); l++) //loop over current path-length and temperature table
			{
				currCollTabTau[l] = currLTTabL[l]; 				    //setting path-lengths
				currCollTabVal[l] = LColl.interp(p, currLTTabT[l]); //setting LColl values
			}

			collisionalEL.push_back(LinearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathL = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (int l=0; l<currLTTabL.size(); l++) temp += currLTTabT[l];
		temp /= currLTTabL.size();
	}
	else { //if path-length is smaller than thermalization time:

		pathL = 0.0; //setting path-length and temperature
		temp  = 0.0;
	}
}

//function that sets energy loss parameters:
void SetELParameters()
{
	//setting energy loss parameters:
	double T = 3.0 / 2.0*TCRIT;
	double mu = 0.197*sqrt((-8.0*(6+nf)*M_PI*M_PI*T*T)/(2.0*nf-33.0)/lambda/lambda/ProductLog((-8.0*(6+nf)*M_PI*M_PI*T*T)/(2.0*nf-33.0)/lambda/lambda));
	mgC = mu / sqrt(2.0);
	if (pName == "Bottom") MC = 4.75;
	else if (pName == "Charm") MC = 1.2;
	else if (pName == "Gluon") MC = mu/sqrt(2.0);
	else MC = mu/sqrt(6.0);
	TCollConst = T;
}

//function that calculates averaged energy loss:
void AverageEL()
{
	SetELParameters(); //setting energy loss parameters

	Grids.SetGridPoints(pName); //generating grids

	if (LoadLdndx()    != 1) return; //loading Ldndx table
	if (LoadLNorm()    != 1) return; //loading LNorm table
	if (LoadLColl()    != 1) return; //loading LColl table
	
	if (LoadTProfile() != 1) return; //loading TProfile

	if (GenerateInitPosPoints() != 1) return; //generating initial position points
	
	
	double PLDist[1000] = {0.0}, TDist[1000] = {0.0}; //deffining path-length and temperature distribution arrays


	if ((pName == "Bottom") || (pName == "Charm")) {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//HEAVY FLAVOR CALCULATION:

		if (Loaddsdpti2()  != 1) return; //loading dsdpti2

		vector<vector<double>> RAAdist(Grids.finPtsLength(), vector<double>(phiGridN)); //defining RAA(pT,phi) distribution 2d vector

		FdAHaltonSeqInit(150); //initialozing Halton sequences

		#pragma omp parallel for reduction(+:PLDist[:1000],TDist[:1000]) schedule(dynamic)
		for (int phi_i=0; phi_i<phiGridN; phi_i++) //loop over angles
		{
			double phi = phiGridPts[phi_i]; //setting value of phi

			//vectors that stores RAA(pT) for single value of phi:
			vector<double> sumRAA1(Grids.finPtsLength()); fill(sumRAA1.begin(), sumRAA1.end(), 0.0);

			vector<vector<double>> sumRAA2(Grids.finPtsLength(), vector<double>(Grids.FdpPtsLength()));
			for (int i=0; i<Grids.finPtsLength(); i++) fill(sumRAA2[i].begin(), sumRAA2[i].end(), 0.0);

			double  weightsumEL = 0.0; //defining energy loss weightsum
			double weightsumPLT = 0.0; //defining path-length and temperature weightsum

			for (int xy_i=0; xy_i<xGridPts.size(); xy_i++) //loop over x and y initial position points
			{
				double x = xGridPts[xy_i], y = yGridPts[xy_i]; //setting values of x and y

				if (BinCollDensity.interp(x,y) > 0.0) { //checking if binary collision density is larger than 0

					weightsumEL += BinCollDensity.interp(x,y); //adding to weightsum

					vector<double> radRAA1;	vector<vector<double>> radRAA2; vector<double> collEL; //defining vectors for raditive RAA and collisional enerfy loss for single trajectory
					double path_length, temperature;											   //defining path-length and temperature for single trajectory
					RadCollEL(x, y, phi, radRAA1, radRAA2, collEL, path_length, temperature);	   //calculating energy loss for single trajectory

					if (path_length > tau0) { //checking if path-length is larger than thermalization time

						for (int pc_i=0; pc_i<Grids.pCollPtsLength(); pc_i++) collEL[pc_i] += 1e-12; //modifying collEL to prevent division by 0

						weightsumPLT  += BinCollDensity.interp(x,y);			   //calculating path-length and temperature distributions
						PLDist[phi_i] += (path_length*BinCollDensity.interp(x,y));
						 TDist[phi_i] += (temperature*BinCollDensity.interp(x,y));

						vector<double> singleRAA1; vector<vector<double>> singleRAA2;
						GaussFilterIntegrate(radRAA1, radRAA2, collEL, singleRAA1, singleRAA2); //performing Gauss filter integration

						//multiplying RAA with bin coll density as weigth and adding to RAA sum
						for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
						{
							sumRAA1[f_i] += singleRAA1[f_i]*BinCollDensity.interp(x,y);

							for (int dp_i=0; dp_i<Grids.FdpPtsLength(); dp_i++)
								sumRAA2[f_i][dp_i] += singleRAA2[f_i][dp_i]*BinCollDensity.interp(x,y);
						}
					}
					else {
						//multiplying RAA1 (which is 1) with binary collision function as weigth and adding to RAA sum; RAA2 is 0 in this case:
						for (int f_i=0; f_i<Grids.finPtsLength(); f_i++) sumRAA1[f_i] += BinCollDensity.interp(x,y);
					}

				}
			}

			for_each(sumRAA1.begin(), sumRAA1.end(), [weightsumEL](double &c){ c/=weightsumEL; }); //dividing sumRAA1 with weight sum

			for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
				for_each(sumRAA2[f_i].begin(), sumRAA2[f_i].end(), [weightsumEL](double &c){ c/=weightsumEL; }); //dividing sumRAA2 with weight sum

			//setting RAA(pT,phi) value by integrating over p:
			for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
				RAAdist[f_i][phi_i] = sumRAA1[f_i] + CubicIntegrate(Grids.FdpPts(), sumRAA2[f_i])/Grids.finPts(f_i);

			PLDist[phi_i] /= weightsumPLT; TDist[phi_i] /= weightsumPLT; //norming path-length and temperature evolutions
		}

		vector<double> RAA, v2;
		CalcObservables(RAAdist, RAA, v2); //calculating observables: RAA(pT) and v2(pT)

		vector<double> avgPLength, avgTemp;
		CalcAvgPLT(PLDist, TDist, avgPLength, avgTemp); //calculating path-lengths and temperatures

		if (ExportResults(RAAdist, RAA, v2, avgPLength, avgTemp) != 1) return; //exporting results to file
	
	}
	else if (pName == "LQuarks") {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ALL LIGHT QUARKS CALCULATION:

		vector<string> LQuarksList{"Down", "DownBar", "Strange", "Up", "UpBar"}; //defining light quarks particle name list

		vector<vector<vector<double>>> RAAdist(LQuarksList.size(), vector<vector<double>>(Grids.finPtsLength(), vector<double>(phiGridN))); //defining RAA(pT,phi) distribution 3d vector

		FdAHaltonSeqInit(100); //initializing Halton sequences

		vector<interpFun> dsdpti2lq(5); //defining initial pT distributions for all light quarks

		for (int i=0; i<LQuarksList.size(); i++) if (Loaddsdpti2(LQuarksList[i], dsdpti2lq[i]) != 1) return; //loading light quarks initial pT distributions

		#pragma omp parallel for reduction(+:PLDist[:1000],TDist[:1000]) schedule(dynamic)
		for (int phi_i=0; phi_i<phiGridN; phi_i++) //loop over angles
		{
			double phi = phiGridPts[phi_i]; //setting value of phi

			//vectors that stores RAA(pT) for single value of phi:
			vector<vector<double>> sumRAA1(LQuarksList.size(), vector<double>(Grids.finPtsLength()));
			for (int i=0; i<sumRAA1.size(); i++) fill(sumRAA1[i].begin(), sumRAA1[i].end(), 0.0);

			vector<vector<vector<double>>> sumRAA2(LQuarksList.size(), vector<vector<double>>(Grids.finPtsLength(), vector<double>(Grids.FdpPtsLength())));
			for (int i=0; i<sumRAA2.size(); i++)
				for (int j=0; j<sumRAA2[i].size(); j++)
					fill(sumRAA2[i][j].begin(), sumRAA2[i][j].end(), 0.0);

			double  weightsumEL = 0.0; //defining energy loss weightsum
			double weightsumPLT = 0.0; //defining path-length and temperature weightsum

			for (int xy_i=0; xy_i<xGridPts.size(); xy_i++) //loop over x and y initial position points
			{
				double x = xGridPts[xy_i], y = yGridPts[xy_i]; //setting values of x and y

				if (BinCollDensity.interp(x,y) > 0.0) { //checking if binary collision density is larger than 0

					weightsumEL += BinCollDensity.interp(x,y); //adding to weightsum

					vector<double> radRAA1;	vector<vector<double>> radRAA2; vector<double> collEL; //defining vectors for raditive RAA and collisional energy loss for single trajectory
					double path_length, temperature;											   //defining path-length and temperature for single trajectory
					RadCollEL(x, y, phi, radRAA1, radRAA2, collEL, path_length, temperature);	   //calculating energy loss for single trajectory

					if (path_length > tau0) { //checking if path-length is larger than thermalization time

						for (int pc_i=0; pc_i<Grids.pCollPtsLength(); pc_i++) collEL[pc_i] += 1e-12; //modifying collEL to prevent division by 0

						weightsumPLT  += BinCollDensity.interp(x,y);			   //calculating path-length and temperature distributions
						PLDist[phi_i] += (path_length*BinCollDensity.interp(x,y));
						 TDist[phi_i] += (temperature*BinCollDensity.interp(x,y));

						vector<vector<double>> singleRAA1(LQuarksList.size());
						vector<vector<vector<double>>> singleRAA2(LQuarksList.size());

						for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
							GaussFilterIntegrate(dsdpti2lq[lq_i], radRAA1, radRAA2, collEL, singleRAA1[lq_i], singleRAA2[lq_i]); //performing Gauss filter integration

						//multiplying RAA with bin coll density as weigth and adding to RAA sum:
						for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
						{
							for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
							{
								sumRAA1[lq_i][f_i] += singleRAA1[lq_i][f_i]*BinCollDensity.interp(x,y);

								for (int dp_i=0; dp_i<Grids.FdpPtsLength(); dp_i++)
									sumRAA2[lq_i][f_i][dp_i] += singleRAA2[lq_i][f_i][dp_i]*BinCollDensity.interp(x,y);
							}
						}
					}
					else {
						//multiplying RAA1 (which is 1) with binary collision function as weigth and adding to RAA sum; RAA2 is 0 in this case:
						for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
							for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
								sumRAA1[lq_i][f_i] += BinCollDensity.interp(x,y);
					}

				}	
			}

			for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
			{
				for_each(sumRAA1[lq_i].begin(), sumRAA1[lq_i].end(), [weightsumEL](double &c){ c/=weightsumEL; }); //dividing sumRAA1 with weight sum

				for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
					for_each(sumRAA2[lq_i][f_i].begin(), sumRAA2[lq_i][f_i].end(), [weightsumEL](double &c){ c/=weightsumEL; }); //dividing sumRAA2 with weight sum
			}

			//setting RAA(pT,phi) value by integrating over p:
			for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
				for (int f_i=0; f_i<Grids.finPtsLength(); f_i++)
					RAAdist[lq_i][f_i][phi_i] = sumRAA1[lq_i][f_i] + CubicIntegrate(Grids.FdpPts(), sumRAA2[lq_i][f_i])/Grids.finPts(f_i);

			PLDist[phi_i] /= weightsumPLT; TDist[phi_i] /= weightsumPLT; //norming path-length and temperature evolutions
		}

		vector<vector<double>> RAA(LQuarksList.size()), v2(LQuarksList.size());
		for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++) CalcObservables(RAAdist[lq_i], RAA[lq_i], v2[lq_i]); //calculating observables: RAA(pT) and v2(pT)

		vector<double> avgPLength, avgTemp;
		CalcAvgPLT(PLDist, TDist, avgPLength, avgTemp); //calculating path-lengths and temperatures

		for (int lq_i=0; lq_i<LQuarksList.size(); lq_i++)
			if (ExportResults(LQuarksList[lq_i], RAAdist[lq_i], RAA[lq_i], v2[lq_i], avgPLength, avgTemp) != 1) return; //exporting results to file

	}
	else {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//LIGHT FLAVOR CALCULATION (INDIVIDUAL PARTICLES):

		if (Loaddsdpti2()  != 1) return; //loading dsdpti2

		vector<vector<double>> RAAdist(Grids.finPtsLength(), vector<double>(phiGridN)); //defining RAA(pT,phi) distribution 2d vector

		dAHaltonSeqInit(1000); //initializing Halton sequences for dA integration

		#pragma omp parallel for reduction(+:PLDist[:1000],TDist[:1000]) schedule(dynamic)
		for (int phi_i=0; phi_i<phiGridN; phi_i++) //loop over angles
		{
			double phi = phiGridPts[phi_i]; //setting value of phi

			//vector that stores RAA(pT) for single value of phi:
			vector<double> sumRAA(Grids.finPtsLength()); fill(sumRAA.begin(), sumRAA.end(), 0.0);

			double  weightsumEL = 0.0; //defining energy loss weightsum
			double weightsumPLT = 0.0; //defining path-length and temperature weightsum

			for (int xy_i=0; xy_i<xGridPts.size(); xy_i++) //loop over x and y initial position points
			{
				double x = xGridPts[xy_i], y = yGridPts[xy_i]; //setting values of x and y

				if (BinCollDensity.interp(x,y) > 0.0) { //checking if binary collision density is larger than 0

					weightsumEL += BinCollDensity.interp(x,y); //adding to weightsum

					vector<double> radRAA, collEL; 									//radiative RAA and collisional energy loss vectors
					double path_length, temperature;								//defining path-length and temperature for single trajectory
					RadCollEL(x, y, phi, radRAA, collEL, path_length, temperature); //calculating energy loss for single trajectory

					if (path_length > tau0) { //checking if path-length is larger than thermalization time

						for (int pc_i=0; pc_i<Grids.pCollPtsLength(); pc_i++) collEL[pc_i] += 1e-12; //modifying collEL to prevent division by 0

						weightsumPLT  += BinCollDensity.interp(x,y);			   //calculating path-length and temperature distributions
						PLDist[phi_i] += (path_length*BinCollDensity.interp(x,y));
						 TDist[phi_i] += (temperature*BinCollDensity.interp(x,y));
						
						vector<double> singleRAA;
						GaussFilterIntegrate(radRAA, collEL, singleRAA); //performing Gauss filter integration
						
						//multiplying RAA with bin coll density as weigth and adding to RAA sum
						for (int f_i=0; f_i<Grids.finPtsLength(); f_i++) sumRAA[f_i] += singleRAA[f_i]*BinCollDensity.interp(x,y);
					}
					else {
						//multiplying RAA (which is 1) with binary collision function as weigth and adding to RAA sum
						for (int f_i=0; f_i<Grids.finPtsLength(); f_i++) sumRAA[f_i] += BinCollDensity.interp(x,y);
					}
				}
			}

			for (int f_i= 0; f_i<Grids.finPtsLength(); f_i++) RAAdist[f_i][phi_i] = sumRAA[f_i]/weightsumEL; //setting value of RAA(pT,phi)

			PLDist[phi_i] /= weightsumPLT; TDist[phi_i] /= weightsumPLT; //norming path-length and temperature evolutions
		}

		vector<double> RAA, v2;
		CalcObservables(RAAdist, RAA, v2); //calculating observables: RAA(pT) and v2(pT)

		vector<double> avgPLength, avgTemp;
		CalcAvgPLT(PLDist, TDist, avgPLength, avgTemp); //calculating path-lengths and temperatures

		if (ExportResults(RAAdist, RAA, v2, avgPLength, avgTemp) != 1) return; //exporting results to file
	}
}