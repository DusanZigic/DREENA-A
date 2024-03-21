<h1><img src="logo/qgptomography.png" alt="logo" width='215' align="right"/> DREENA-A</h1>

DREENA-A is computational framework for generating high-pT predictions based on a dynamical energy loss formalism. The framework can include 
any, in principle arbitrary, temperature evolution within the dynamical energy loss formalism. This version is generalized to account for both LHC and RHIC energies and collision systems.


## < 1 > compilation

Compilation of the source code, performed using gcc compiler:

g++ source/*.cpp -fopenmp -O3 -o ebeDREENA

a) -fopenmp is necessary to enable parallelization using OpenMP;  
b) -O3 optimization is not necessary, but recommended;  
c) -o DREENAA is optional and if omitted, the output of the compilation will be placed in a.out;  


## < 2 > prerequisite files

All prerequisite files need to be textual tables. They can have a different number of columns depending on the file. They can contain headers that start with '#'.

#### a) initial pT distributions

initial pT distributions file should have 2 columns in format:

pT | dsigma/d(pT^2) |
--- | --- |
... | ... |

for heavy flavor, initial pT distributions can be obtained from this [web site](http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html);  
> [!CAUTION]
> note that the default output of this web interface is dsigma/dpT, while ebeDREENA initial pT distribution input needs to be dsigma/d(pT^2), so these distributions need to be modified;

to avoid modifying different parameters that are hard-coded, initial pT distribution for heavy flavour should be in the range of 1GeV to at least 200GeV and for light flavour from 1GeV to at least 450GeV; for most distributions pT step of 1GeV seems to be sufficient;  

initial pT distribution file path relative to the executable should be: `./pTDists/ptDists[sNN]/ptDist_[sNN]_[particleName].dat`, where *sNN* is collision energy that can be 200GeV, 2760GeV, 5020GeV and *particleName* is the name of the particle that can be Bottom, Charm, Down, DownBar, Gluon, Strange, Up, UpBar;

#### b) binary collision density

binary collision density file contains jet creation probability in transverse plane as a function of x and y; x and y need to form ordered grid with run order y > x, while probability is the 3rd column in the table;  
the binary collision density can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;  
binary collision density needs to be consistent with initial conditions used to generate temperature evolution;  
binary collision density provided in this demo (./binarycolldensities/binarycolldensity_cent=30-40%.dat) is generated using optical Glauber model (see [arxiv:2110.01544](https://inspirehep.net/literature/2606181) for more details) and it corresponds to provided temperature evolution;

#### c) temperature evolution

temperature evolution file contains temperature as a function of proper time, tau, and x and y spatial coordinates in the transverse plane in that order; tau, x and y need to form an ordered grid with run order y > x > tau;  
temperature evolution file can contain an additional column with energy density evolution; if the file contains temperature and energy density evolution, the order of these two columns can be arbitrary (energy density can be 4th column and temperature 5th, or vice versa);  
the evolution can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;  
the evolution provided in this demo (./evols/tempevol_cent=30-40%.dat) is generated evolvong optical Glauber initial conditions using 3D hydro model (see [arxiv:2110.01544](https://inspirehep.net/literature/2606181) for more details);

#### d) LTables

LTables files contain pre-generated radiated gluon rates and collisional energy loss; for radiative energy loss, there are 2 tables lnorm table with column format:  
tau | p | temp | LNorm
--- | --- | --- | --- |
... | ... | ... | ... |

and ldndx table with column format:  
tau | p | temp | x | Ldndx
--- | --- | --- | --- | --- |
... | ... | ... | ... | ... |

these tables are also a function of, effective number of flavours, (nf=2.5 for sNN=200GeV and nf=3.0 for higher collision energies), particle type and magnetic to electric mass ratio, xB, which figure in file names;  
for collisional energy loss, there is one table lcoll with column format:
p | temp | LColl
--- | --- | --- |
... | ... | ... |

this table is also a function of particle mass (particle name is in file name) and effective number of flavours, nf;

LTables files path relative to the executable should be:

+ `./ltables/lnorm_nf=[nf]_[particleName]_xB=[xB].dat`,
+ `./ltables/ldndx_nf=[nf]_[particleName]_xB=[xB].dat` and
+ `./ltables/lcoll_nf=[nf]_[particleName].dat`,

where *nf* is the effective number of flavours that can be 2.5 for 200GeV or 3.0 for LHC energy collisions, *particleName* is the name of the particle that can be Bottom, Charm, LQuarks or Gluon (all light quarks are taken to have same mass, so their LTables are the same), and *xB* is the chromo-magnetic to electric mass ratio;

unlike previous files, DREENA-A calculates LTables; however, these tables need to be calculated only once and can be reused while calculating high-pT energy loss with different temperature evolution backgrounds;

within this repository there is an example of LTables files for Charm quark, nf=3.0 and xB=0.6;

#### e) phiGausPts

in ./phiGaussPts/ directory are textual tables containing jet's direction angles and weights that correspond to Gaussian quadrature integration method in range [0, 2Pi];  
jet's direction angles are sampled in these points, so that afterwards, when integrating final pT,phi distribution over phi, which is nedeed to obtain R_AA and v_2 there is no need for angle resampling;

## < 3 > run DREENA-A

There are two possible calculation options within DREENA-A framework: LTables calculation and energy loss calculation. Since DREENA-A is parallelized using OpenMP, set OMP_NUM_THREADS environmental variable to desired value before running calculations.

#### a) LTables calculation

to see all parameters and their default values run:  
```
./DREENAA LTables -h
```
parameters for LTables calculations are:    
+ **sNN** parameter: case sensitive string with possible options: 200GeV, 2760GeV, 5020GeV, 5440GeV;  
*default value: PbPb*

+ **pName** parameter: case sensitive string with possible options: Bottom, Charm, LQuarks, Gluon, where LQuarks stand for light quarks, since all light quarks (down, down-bar, strange, strange-bar, up and up-bar) use the same LTables;  
*default value: Charm*

+ **xB** parameter: float that represents magnetic to electric mass ratio; based on lattice calculation: xB=0.6;  
*default value: 0.6*

additional parameters are number of points used for QuasiMonteCarlo integration LdndxMaxPoints and LCollMaxPoints, with default values 500000 and 10000, respectively;  

to generate LTables provided in this demo, that are for 5020GeV collision energy, charm quark and for xB value of 0.6, use:

```
./DREENAA LTables --sNN=5020GeV --pName=Charm --xB=0.6

```

or just:

```
./DREENAA LTables

```

since these are all default parameter values;

#### b) energy loss calculation

to see all parameters and their default values run:  
```
./DREENAA AverageEL -h
```
parameters for energy loss calculations are:  
+ **collsys** parameter: case sensitive string with possible options: AuAu, PbPb, XeXe,...  
*default value: PbPb*

+ **sNN** parameter: case sensitive string with possible options: 200GeV, 2760GeV, 5020GeV, 5440GeV;  
*default value: 5020GeV*

+ **pName** parameter: case sensitive string with possible options: Bottom, Charm, Gluon, Down, DownBar, Strange, Up, UpBar, LQuarks;  
the calculation can be done for each parton individualy, however there is an option to calculate all light quarks at the same time using modified algorithm since the only thing differentiating light quarks is initial pT distribution; this leads to 4x speed-up of the calculation time compared to calculating each light quark individually;  
*default value: Charm*

+ **centrality** parameter: string in format 'xx-xx%' (ie 0-5%, 10-20%,...);  
*default value: 30-40%*

+ **xB** parameter: float representing magnetic to electric mass ratio; based on latest lattice calculation: xB=0.6;  
*default value: 0.6*

+ **xGriN** and **yGridN** parameters: positive integers representing number of x and y points on the equidistant grid in transverse plane that will be used as initial position points for jets;  
*default values: 50, 50*

+ **phiGridN** parameter: positive integer representing number of sampled angles between jet's trajectory and the x-axis in the transverse plane;  
*default value: 25*

+ **TIMESTEP** parameter: positive float representing timestep along jet's trajectory in fm;  
*default value: 0.1*

+ **TCRIT** parameter: positive float representing critical temperature in GeV, i.e. temperature value for which jet's energy loss stops;  
*default value: 0.155*

to generate results provided in this demo (*./results/resultsCharm/Charm_PbPb_sNN=5020GeV_cent=30-40%_xB=0.6_dist.dat*) provided in this demo, that are for Pb+Pb 5020GeV collisions, charm quarks, 30-40% centrality class and for xB value of 0.6, use:

```
./DREENAA AverageEL collsys=PbPb --sNN=5020GeV --pName=Charm --centrality=30-40% --xB=0.6 --xGridN=50 --yGridN=50 --phiGridN=25

```

or just:

```
./DREENAA AverageEL

```

since these are all default parameter values;


## < 4 > output

#### a) LTables output

Output of LTables run is in ./ltables/ directory; for single run, there are three output files - two for radiative and one for collisional part of energy loss.  

These tables are determined by collision energy throught effective number of flavours, nf, particle type and value of xB.  

These tables will later be used for energy loss calculations on hydro background and they can be reused for different hydro backgrounds.

As an example, tables for 5020GeV collisions, charm quark and xB value of 0.4 are given.

#### b) energy loss output

Output of energy loss run is in *./results/* directory;  

These are R_AA(pT,phi) distributions that are later used to calculate R_AA and v_2; for more details on how R_AA(pT) and v_2(pT) are calculated see equations (16) and (17) in [arxiv:2110.01544](https://inspirehep.net/literature/2606181).  

These files also have headers, that contain various informations about the energy loss calculation such as number of angles and points in x-y plane used as jet's initial positions, average jet path-length and temperature jets experience along the trajectory,...  

For file name pattern see file provided in this repository.

As an example file contaning RAA(pT,phi) for charm quark at 30-40% centrality class for Pb+Pb 5020GeV collisions and xB value of 0.6 is provided. This RAA(pT,phi) distribution is obtained calculating energy loss on hydro background provided in this demo (./evols/tempevol_cent=30-40%.dat), while the jet's intial positions are generated based on probability also provided in this demo (./binarycolldensities/binarycolldensity_cent=30-40%.dat)

## < 5 > disclaimer

This version of the DREENA-A framework is meant for high-pT energy loss calculation on smooth hydro background, i.e., it cannot be used for event-by-event calculations.

DREENA-A framework provides results for bare quarks and gluons. To obtain results comparable to experimental data, fragmentation functions such as DSS, BCFY and KLP have to be used.

Calculation time for LTables calculation is large (up to about one hour on 112 cores). However, they need to be calculated only once. Energy loss calculation time is mainly dependent on the number of sampled jet's trajectory and for default values of parameters (xGridN=50, yGridN=50, phiGridN=25) whose values have shown to be sufficient for multiple different temperature evolutions and binary collision densities, and it is about several minutes on 4 cores, which means it can easily be run locally.

-------------------------------------------------------

for questions contact Zigic Dusan at zigic@ipb.ac.rs
