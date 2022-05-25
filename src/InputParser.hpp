#ifndef HEADERFILE_INPUTHEADER
#define HEADERFILE_INPUTHEADER

int GetInputs(vector<string> inputs, string &pname, string &cent, double &xb, int &xptsn, int &yptsn, int &phiptsn, string &pTinitpath, string &temppath, string &bcdpath, double &timestep, double &tcrit); //function that gets inputs for energy loss calculations
int GetInputs(vector<string> inputs, string &pname, double &xb, int &LdndxMaxPts, int &LCollMaxPts);				 //function that gets inputs for LTables calculations

#endif