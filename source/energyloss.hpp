#ifndef HEADERFILE_ELOSSHEADER
#define HEADERFILE_ELOSSHEADER

#include "grids.hpp"
#include "linearinterpolation.hpp"

#include <string>
#include <vector>

class energyLoss {

public:
    energyLoss(int argc, const char *argv[]);
    ~energyLoss();
    void runEnergyLoss();

private:
    bool m_error; //flag that checks if previous calculation is done properly

    std::string m_collsys;                 // collision system
    std::string m_sNN;                     // collision energy
    std::string m_pName;                   // particle name
    std::string m_centrality;              // centrality class
    double m_xB;                           // xB value
    size_t m_xGridN, m_yGridN, m_phiGridN; // initial position grid points and angle number
    double m_TIMESTEP, m_TCRIT;	           // time step and critical temperature

    double m_nf;			     // effective number of flavours
    const double m_lambda = 0.2; // QCD scale
    double m_mgC, m_MC;          // constant particle and gluon masses used for dA integrals
    double m_TCollConst;         // constant temperature used for Gauss filter integration
    double m_tau0;               // thermalization time

    gridPoints m_Grids; // grid points
    

    double productLog(double x) const;

    interpolationF<double> m_dsdpti2; // initial pT distribution interpolated function
    int loaddsdpti2();
    int loaddsdpti2(const std::string &pname, interpolationF<double> &dsdpti2int) const;
    
    interpolationF<double> m_LNorm, m_Ldndx, m_LColl; // interpolated L tables
    int loadLdndx();
    int loadLNorm();
    int loadLColl();

    interpolationF<double> m_tempEvol; // temperature evolution interpolated function
    int loadTempEvol();
    
    interpolationF<double> m_binCollDensity;                  // binary collision density interpolated function
    int loadBinCollDensity(interpolationF<double> &binCollDensity);
    int loadPhiPoints(std::vector<double> &phipoints);
    int loadBinCollPoints(std::vector<std::vector<double>> &bcpoints);
    std::vector<double> m_xGridPts, m_yGridPts, m_phiGridPts; // vectors that store initial position points and angles
    int generateInitPosPoints();
};

#endif