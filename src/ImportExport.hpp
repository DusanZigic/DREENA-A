#ifndef HEADERFILE_IMPORTEXPORTHEADER
#define HEADERFILE_IMPORTEXPORTHEADER

int Loaddsdpti2();								   		 //function that loads dsdpti2
int Loaddsdpti2(string pname, interpFun &dsdpti2int);	 //function that loads dsdpti2 for given particle
int LoadLdndx();								   		 //function that loads Ldndx
int LoadLNorm();								   		 //function that loads LNorm
int LoadLColl();								   		 //function that loads LColl
int LoadTProfile(); 	  						   		 //function that loads TProfile
int LoadBinCollDensity();						   		 //function that loads binary collision density
int LoadBinCollDensity(interpFun &BinCollDensInt); 		 //function that loads binary collision density for given interpFun object
int LoadBinCollPoints(vector<vector<double>> &bcpoints); //function that loads binary collision points

int ExportResults(vector<vector<double>> RAApTphi, vector<double> RAApT, vector<double> v2pT, vector<double> avg_pl, vector<double> avg_temp);					 //function that exports results
int ExportResults(string part_name, vector<vector<double>> RAApTphi, vector<double> RAApT, vector<double> v2pT, vector<double> avg_pl, vector<double> avg_temp); //function that exports results - lquarks algorithm
int ExportLTables(vector<double> &ldndxtbl, vector<double> &lnormtbl, vector<double> &lcolltbl); 																 //function that export LTables

#endif