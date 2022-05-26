===========================================
DREENA-A
===========================================

DREENA-A is computational framework for generating high-pT predictions based on a dynamical energy loss formalism. The framework can include 
any, in principle arbitrary, temperature profile within the dynamical energy loss formalism.


-----------------------------
< 1 > compilation
-----------------------------

Compilation of the source code is performed using gcc v7.5.0 on ubuntu 18.04 inside terminal opened in the main directory with the command:

g++ src/*.cpp -fopenmp -O3 -o DREENAA

a) -fopenmp is necessary to enable parallelization using OpenMP;
b) -O3 optimization is not necessary, but recommended;
c) -o DREENAA is optional and if omitted, the output of the compilation will be placed in a.out;


-----------------------------
< 2 > prerequisite files
-----------------------------

All prerequisite files need to be textual tables. They can have a different number of columns depending on the file. They can contain headers in which
lines do not start with a number or a minus sign, "-". Indentation for numbers is accounted for.

a) initial pT distributions
   
   initial pT distributions file should have 2 columns in format: pT | dsigma/d(pT^2);  

   for heavy flavor, initial pT distributions can be obtained from http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html
   *IMPORTANT* note that the default output of this web interface is dsigma/dpT, while DREENA-A initial pT distribution input needs to
               be dsigma/d(pT^2), so these distributions need to be modified;
               to avoid modifying different parameters that are hard-coded, initial pT distribution for heavy flavour should be in the range of 1GeV
               to at least 200GeV and for light flavour from 1GeV to at least 450GeV;
               for most distributions pT step of 1GeV seems to be sufficient;

b) temperature evolution

   temperature evolution file contains temperature as a function of proper time, tau, and x and y spatial coordinates in the transverse plane
   in that order; tau, x and y need to form an ordered grid with run order y > x > tau;
   temperature evolution file can contain an additional column with energy density evolution; if the file contains temperature and energy
   density evolution, the order of these two columns can be arbitrary (energy density can be 4th column and temperature 5th, or vice versa);
   the evolution can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;

   *the evolution provided in this demo (./TProfiles/TProfile_cent=30-40%.dat) corresponds to 'Glauber' evolution outlined in the manuscript;   

c) binary collision density

   binary collision density file contains jet creation probability in transverse plane as a function of x and y; x and y need to form ordered
   grid with run order y > x, while probability is the 3rd column in the table;
   the binary collision density can be given only in the 1st quadrant of the transverse plane (x>=0, y>=0) if initial conditions have this symmetry;
   binary collision density needs to be consistent with initial conditions used to generate temperature evolution;

   *binary collision density provided in this demo (./BinaryCollDensities/BinaryCollDensity_cent=30-40%.dat) corresponds to provided temperature evolution;

d) LTables

   LTables files contain pre-generated radiated gluon rates and collisional energy loss; for radiative energy loss, there are 2 tables LNorm table with
   column format: tau | p | temp | LNorm and Ldndx table with column format: tau | p | temp | x | Ldndx; these tables are also a function of particle
   type and magnetic to electric mass ratio, xB, which figure in file names; for collisional energy loss, there is one table LColl with column format:
   p | temp | LColl; this table is also a function of particle type, which figures in the filename; 
   unlike previous files, DREENA-A calculates LTables; however, these tables need to be calculated only once and can be reused while calculating
   high-pT energy loss with different temperature evolution backgrounds;

   *LTables provided in this demo (in ./LTables/ directory) are for charm quark and for xB=0.4;


-----------------------------
< 3 > run DREENA-A
-----------------------------

a) LTables calculation

   to perform basic LTables calculation use command: ./DREENAA LTables -pName=[particle_name] -xB=[xB]

   particle_name parameter: case sensitive string with possible options: Bottom, Charm, LQuarks, Gluon, where LQuarks stand for light quarks, since
                            all light quarks (down, down-bar, strange, strange-bar, up and up-bar) use the same LTables;
   xB parameter:            float that represents magnetic to electric mass ratio; based on lattice calculation: 0.4<=xB<=0.6;

   all parameters are optional; to see all parameters that can be changed and their default values (values that will be used if another value is not
   provided) use: ./DREENAA LTables -h

   since LTables calculation is time consuming, pregenerated LTables for charm quark with xB=0.4 is provided in ./LTables/ directory;

b) energy loss calculation

   to perform basic energy loss calculation use command: ./DREENAA AverageEL -pName=[particle_name] -centrality=[centrality] -xB=[xB]

   particle_name parameter: case sensitive string with possible options: Bottom, Charm, Down, DownBar, Gluon, Strange, Up, UpBar, or LQuarks, where
                            LQuarks stands for light quarks; the calculation for all light quarks is performed at the same time (preferable
                            if all light flavor particles are needed for fragmentation functions);
                            pregenerated LTables need to match particle name parameter to run energy loss calculation;
   centrality:              string in format 'xx-xx%' (ie 0-5%, 10-20%,...) that represents centrality class parameter
   xB parameter:            float that represents magnetic to electric mass ratio; based on lattice calculation: 0.4<=xB<=0.6;
                            pregenerated LTables need to match xB parameter to run energy loss calculation;
   
   additional parameters are:
   -xGriN, yGridN: number of x and y points on the equidistant grid in transverse plane that will be used as initial position points for jets (default values:
                   -xGridN=40 -yGridN=40);
   -phiGridN:      number of sampled angles between jet's trajectory and the x-axis in the transverse plane (default value: -phiGridN=50)
   -pTinit_path:   absolute path or path relative to executable of initial pT distribution; if omitted, the default path is used: ./pTinitDists/pTinitDist_[particle_name].dat
   -temp_path:     absolute path or path relative to executable of temperature evolution file; if omitted, the default path is used: ./TProfiles/TProfile_cent=[centrality].dat
   -bcd_path:      absolute path or path relative to executable of binary collision density file; if omitted, the default path is used: ./BinaryCollDensities/BinaryCollDensity_cent=[centrality].dat
   -TIMESTEP:      timestep along jet's trajectory (default value: 0.1; unit: fm)
   -TCRIT:         critical temperature (default value: 0.155; unit: GeV)

   all parameters are optional; to see all parameters that can be changed and their default values (values that will be used if another value is not
   provided) use: ./DREENAA AverageEL -h

   code is parallelized using OpenMP, which means the number of cores used by the program can be controlled with the OMP_NUM_THREADS environmental variable; for example,
   to run energy loss on 4 cores for charm and for 30-40% centrality, command would be: export OMP_NUM_THREADS=4; ./DREENAA AverageEL -pName=Charm -centrality=30-40%;


-----------------------------
< 4 > output
-----------------------------

a) LTables output

   output of LTables run is ./LTables directory; for single run, there are three output files - two for radiative and one for collisional; these tables are
   determined by particle type and value of xB;
   these tables will later be used for energy loss calculations on hydro background and they can be reused for different hydro bacgrounds;
   as an example, tables for charm quark are given for xB value of 0.4;

b) energy loss output

   energy loss run outputs two files: one containing RAA(pT, phi) and the other one containing RAA(pT) and v_2(pT), where RAA(pT) and v_2(pT) are calculated
   from RAA(pT, phi);
   files are contained in ./CResults/CResults_[particle_name] directory and file name format is: [particle_name]_5TeV_cent=[centrality]_xB=[xB]_dist.dat
   for file that contains RAA(pT, phi) and [particle_name]_5TeV_cent=[centrality]_xB=[xB]_obs.dat for file that contains RAA(pT) and v_2(pT);
   both files have headers that contain additional data such as: collision system, collision energy, particle type, centrality, xB value, average path lengths
   and temperatures both contaning three values: angular-averaged, in-plane and out-of-plane, as well as number of sampled angles and initial (x,y) points and
   column description of the data that follows;
   both files have three columns of data after header; file that containes RAA(pT, phi) (*_dist.dat file) has: pT | phi | RAA(pT,phi), while the file that
   contains RAA(pT) and v_2(pT) (*_obs.dat file) has: pT | RAA(pT) | v_2(pT);
   as an example dist (RAA(pT,phi)), and obs (RAA(pT) and v2(pT)) files for charm quark are given; energy loss is calculated using initial momentum distributions
   obtained using web interface mentioned above and on hydro background that is provided in this demo;


-----------------------------
< 5 > disclaimer
-----------------------------

This version of the DREENA-A framework is meant for high-pT energy loss calculation on smooth hydro background, i.e., it cannot be used for event-by-event
calculations.

DREENA-A framework provides results for bare quarks and gluons. To obtain results comparable to experimental data, fragmentation functions such as DSS,
BCFY and KLP have to be used.

This version of DREENA-A framework is tuned to LHC energies.

Calculation time for LTables calculation is large (up to about one hour on 112 cores). However, they need to be calculated only once. Energy loss calculation
time is mainly dependent on the number of sampled jet's trajectory and for default values of parameters (xGridN=40, yGridN=40, phiGridN=50) whose values
have shown to be sufficient for multiple different temperature evolutions and binary collision densities, and it is about several minutes on 4 cores,
which means it can easily be run locally.

-------------------------------------------------------

for questions contact Zigic Dusan at zigic@ipb.ac.rs
