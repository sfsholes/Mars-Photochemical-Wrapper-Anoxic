# Mars-Photochemical-Wrapper-Anoxic
Wrapper for the Mars photochemical model used in Sholes et al. (2017)
doi:10.1029/2018JE005837

This is an old Python wrapper for running the Mars photochemical FORTRAN code for various different volcanic outgassing ratios. 
The photochemical model is based off the Smith et al. (2014) code from the Catling group at the University of Washington. The history of the code and other versions can be found at: https://github.com/VirtualPlanetaryLaboratory/atmos

This version of the wrapper uses Python 2.7 and is fairly clunky. An improved model was built for the updated Sholes et al. (2019, Astrobiology) paper which also introduced an optimization module, and an even smoother version is currently being overhauled (3/2020) which focuses more on editing the actual FORTRAN code to produce better output files that can be more easily read in and plotted with matplotlib. These improvements will rely less on having to deal with the unreliable and finicky nature of the main output and input files for the FORTRAN code. 

To run, copy the files into the same directory as the TOTCdev model and run the sfs_mars1.py file from the terminal. 
