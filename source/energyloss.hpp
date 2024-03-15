#ifndef HEADERFILE_ELOSSHEADER
#define HEADERFILE_ELOSSHEADER

#include "grids.hpp"
#include "linearinterpolation.hpp"

#include <string>
#include <vector>
#include <map>

class energyLoss {

public:
    energyLoss(int argc, const char *argv[]);
    ~energyLoss();
    void runEnergyLoss();

private:
    bool m_error; //flag that checks if previous calculation is done properly

    std::string m_collsys;      // collision system
    std::string m_sNN;          // collision energy
    std::string m_pName;        // particle name
    std::string m_centrality;   // centrality class
    double m_xB;                // xB value
    size_t m_xGridN;            // initial position grid points and angle number
    long m_yGridN;              // initial position grid points and angle number
    size_t m_phiGridN;          // initial position grid points and angle number
    double m_TIMESTEP, m_TCRIT;	// time step and critical temperature

    double m_nf;			     // effective number of flavours
    const double m_lambda = 0.2; // QCD scale
    double m_mgC, m_MC;          // constant particle and gluon masses used for dA integrals
    double m_TCollConst;         // constant temperature used for Gauss filter integration
    double m_tau0;               // thermalization time

    gridPoints m_Grids; // grid points
    

    int loadInputsFromFile(const std::string &filePath, std::map<std::string, std::string> &inputParamsFile);

    double productLog(double x) const;

    interpolationF<double> m_dsdpti2; // initial pT distribution interpolated function
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


    double haltonSequence(int index, int base) const;

    size_t m_FdAMaxPoints2, m_FdAMaxPoints3, m_FdAMaxPoints4, m_FdAMaxPoints5; //number of points for FdA integration
    std::vector<double> m_FdAHS2, m_FdAHS3, m_FdAHS4, m_FdAHS5;                //vectors that store Halton sequences for FdA integrals
    void FdAHaltonSeqInit(size_t FdAMaxPts);
    double dAp410(double ph, const interpolationF<double> &norm) const;
    double FdA411(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double FdA412(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double FdA413(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double FdA414(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double FdA415(double ph, double dp, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double FdA(double ph, double dp, const interpolationF<double> &currnorm, const interpolationF<double> &currdndx) const;
    
    size_t m_dAMaxPoints1, m_dAMaxPoints2, m_dAMaxPoints3, m_dAMaxPoints4, m_dAMaxPoints5, m_dAMaxPoints6, m_dAMaxPoints7; //number of points for dA integration
    std::vector<double> m_dAHS1, m_dAHS2, m_dAHS3, m_dAHS4, m_dAHS5, m_dAHS6, m_dAHS7; 								 	   //vectors that store Halton sequences for dA integrals
    void dAHaltonSeqInit(size_t dAMaxPts);
    double dA410(double ph, const interpolationF<double> &norm) const;
    double dA411(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA412(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA413(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA414(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA415(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA416(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA417(double ph, const interpolationF<double> &norm, const interpolationF<double> &dndx) const;
    double dA41(double ph, interpolationF<double> &currnorm, interpolationF<double> &currdndx) const;

    void radCollEnergyLoss(double x, double y, double phi, std::vector<double> &radRAA1, std::vector<std::vector<double>> &radRAA2, std::vector<double> &collEL, double &pathLength, double &temperature) const;
    void radCollEnergyLoss(double x, double y, double phi, std::vector<double> &radRAA, std::vector<double> &collEL, double &pathLenght, double &temperature) const;

    void generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab) const;
    void gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const;
    void gaussFilterIntegrate(const interpolationF<double> &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const;
    void gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA) const;

    void calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp) const;
    
    int exportResults(const std::string &pName, const std::vector<std::vector<double>> &RAADist, const std::vector<double> avgPathLength, const std::vector<double> avgTemp);

    void runELossHeavyFlavour();
    void runELossLightQuarks();
    void runELossLightFlavour();
};

#endif