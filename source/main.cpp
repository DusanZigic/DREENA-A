/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

																		DREENA-A

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

#include "ltables.hpp"
#include "energyloss.hpp"

#include <iostream>
#include <string>

int main(int argc, char const *argv[])
{
	if (argc == 1) {
		std::cout << "use -h for help" << std::endl;
		return 0;
	}

    if (argv[1] ==  std::string("-h")) {

		std::cout << "posible calculations: AverageEL, LTables" << std::endl;
		std::cout << "use: 'calculation' -h for extra help" << std::endl;
		return 0;
	}
    else if (argv[1] ==  std::string("AverageEL")) {
        energyLoss EnergyLoss(argc, argv);
		EnergyLoss.runEnergyLoss();
    }
    else if (argv[1] == std::string("LTables")) {
		lTables LTables(argc, argv);
		LTables.runLTables();
	}
	else {
		std::cerr << "posible calculations: AverageEL, LTables" << std::endl;
		return -3;
	}

	return 1;
}