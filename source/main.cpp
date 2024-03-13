/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

																		    DREENA-A

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "MainHeader.hpp"
#include "InputParser.hpp"
#include "ELossHeader.hpp"
#include "LTablesHeader.hpp"

int main(int argc, char const *argv[])
{
	if (argc == 1) {
		cout << "type -h for help" << endl;
		return -1;
	}

	if (argv[1] ==  string("-h")) {

		cout << "posible calculations: AverageEL, LTables" << endl;
		cout << "type: 'calculation' -h for extra help" << endl;
		return -2;
	}
	else if (argv[1] ==  string("AverageEL")) {
		
		vector<string> argv_inputs;
		for (int i=2; i<argc; i++) argv_inputs.push_back(argv[i]);

		if (GetInputs(argv_inputs, pName, centrality, xB, xGridN, yGridN, phiGridN, pTinit_path, temp_path, bcd_path, TIMESTEP, TCRIT) == 0) return -3;

		
		AverageEL();
	}
	else if (argv[1] == string("LTables")) {

		vector<string> argv_inputs;
		for (int i=2; i<argc; i++) argv_inputs.push_back(argv[i]);

		if (GetInputs(argv_inputs, LT_pName, LT_xB, LdndxMaxPoints, LCollMaxPoints) == 0) return -4;

		
		GenerateLTables();
	}

	return 1;
}