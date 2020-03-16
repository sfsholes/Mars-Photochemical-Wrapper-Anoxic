#  sfs_species.py
#  Written by: Steven F. Sholes (sfsholes@uw.edu)
#  Last Updated: October 1, 2014
#  
#  This file is used by sfs_mars.py to call the global variables for each of
#  the species. Make sure that the numbering is equivalent to the order the
#  species appear in the out.out output file.
#
#  This also includes the atomic weights of variable species.
#  
#  If you change the number of species, make sure you change NSP

#SPECIES IN ORDER THEY APPEAR IN SPECIES.DAT
O = 0
O2 = 1
H2O = 2
H = 3
OH = 4
HO2 = 5
H2O2 = 6
H2 = 7
CO = 8
HCO = 9
H2CO = 10
CH4 = 11
CH3 = 12
C2H6 = 13
NO = 14
NO2 = 15
HNO = 16
H2S = 17
HS = 18
S = 19
SO = 20
SO2 = 21
H2SO4 = 22
HSO = 23
S2 = 24
S4 = 25
S8 = 26
SO3 = 27
OCS = 28
S3 = 29
O3 = 30
HNO3 = 31
N = 32
HNO4 = 33
NO3 = 34
SO4AER = 35
S8AER = 36

#ATOMIC WEIGHTS OF SPECIES
wtH = 1.0
wtC = 12.0
wtN = 14.0
wtO = 16.0
wtS = 32.0

wtH2O = 2*(wtH) + 1*(wtO)
wtH2 = 2*(wtH)
wtH2S = 2*(wtH) + 1*(wtS)
wtSO2 = 1*(wtS) + 2*(wtO)
wtS2 = 2*(wtS)
wtCO = 1*(wtC) + 1*(wtO)
wtCO2 = 1*(wtC) + 2*(wtO)

NSP = S8AER     #Change this if you change the number of species. Just put the highest number here
NSP1 = NSP + 1
