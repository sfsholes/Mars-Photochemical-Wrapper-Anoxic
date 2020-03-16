#  sfs_mars1.py
#  Written by: Steven F. Sholes (sfsholes@uw.edu)
#  Last Updated: June 12, 2017
#
#  This is used to run and plot the TOTCdev model for the photochemistry of Mars
#  for incremental steps in increasing the amount of volcanic gasses.
#  NOTE: THIS VERSION WORKS FOR PYTHON 2.7 (Has not been updated for Python 3)
#
#  The function run() is used to run the model, it will ask for a run number
#    (used to create a unique directory), notes (used to store information in a
#    readme file in the directory), and information about which Gaillard et al.
#    2013 magma melt values you wish to model (% water content, fugacity, and
#    pressure). 

#  run() calls Gbuffer() which looks up what the corresponding mixtures are for
#    the different species (SO2, H2S, H2, CO, S2) and then calls mulltirun().

#  multirun() then calculates what the flux is for the given parameters and runs
#    modelrun() which modifies the INPUTFILES/species.dat file and runs the
#    photochemical code (./TOTCdev). It loops over a reasonable volcanic flux
#    range, edits the upwards fluxes for the volcanic gases, makes the file,
#    saves the model outputs into the created directory, removes the temporary
#    outputs, and copies the model-produced steady-state mixing ratios as the 
#    new starting point for the next ramped-up flux model run (in.dist).

#  Finally, modelplot() is run to actually plot out what happens in the model
#    runs. It goes through and finds the ground-level mixing ratios of all
#    species as well as the corresponding fluxes (it is only designed for up to
#    89 species in the NSP list, more than that will require additional edits to
#    ensure all the mixing ratios/fluxes are found). It will then plot the data
#    with the following plots: 1) volcanic flux vs. mixing ratio/flux (with -
#    legend70.png - and without - results100.png - a legend; 2) Plots just for
#    testing that show the mixingratios - mixingratios.png - and fluxes -
#    fluxes.png - of only the 'important' species; 3) Redox state of the atmos.
#    as a function of volcanic flux (pOx = 2pO2 - pCO - pH2 - 3pOCS) -
#    redoxstate.png; 4) Plots that are the same as number one, but have been
#    normalized (since the photochemical code does not ensure that total
#    mixing ratios are = 1) - NormalizedMixingRatios.png; 5) Plots for testing
#    that show the mixing ratios of all species grouped into oxygen species -
#    Ospecies.png - carbon and nitrogen species - Cspecies.png - and sulfur
#    species - Sspecies.png; 6) Plots that show the redox state along with the
#    mixing ratios of the dominant redox species - RedoxFigPaper.png; 7) A test
#    plot (currently set to plotting the aerosol fluxes), that can be changed
#    according to what you need debugging.
#
#  Some useful variables to change:
#    - xmaximum:   Highest volcanic flux tested (in km3/yr)
#    - xminimum:   Lowest volcanic flux tested (in km3/yr)
#    - increments: Multiplier step for ramping up volcanism (1.5 is usually good)
#    - volden:     Density of melt (g/km3)
#    - ch4ratio:   Upwards flux of methane (molecules/cm2/s), currently set to
#                    1.0e7 in adherence to hypervelocity impact CO removal
#    - marstemp:   Global mean Mars temperature (209.51 K is good)
#    - NSTEPS:     Same as NSTEPS in TOTCtester (maximum # of time steps)
#
#  *Also, before running the code double check the following files:
#    - sfs_species.py   Assigns numbers to each of the different species. Only
#                         need to change this if you add any species or change
#                         their order in species.dat.
#    - species.dat      Located in INPUTFILES/, Make sure that the deposition
#                         velocities are what you want (0.02 cm/s is good), that
#                         all the species you want are there, and that the
#                         LBOUND, DISTH, MBOUND, SMFLUX are correct. *Important:*
#                         If you change what is coming out of your volcanos make
#                         sure that you change the LBOUND and SGFLUX both here
#                         and below in the modelrun() part of the code.
#
#  Important: make sure sfs_species.py is in the same folder as this file
#  The species' position in the output file are referenced there as globals
#
#  In order to run, while in the folder with TOTCdev, type into the command line:
#  python sfs_mars1.py


# IMPORT SCIENCE/MATH PACKAGES
import subprocess
import numpy
import pylab
import matplotlib.pyplot as plt
from time import strftime, time
from sfs_species import *

subprocess.call("source ~/.profile", shell=True)     #This is to make sure the bashrc is correct
#subprocess.call("echo $PATH", shell=True)           #Used to test that the bashrc is correct, it should contain the correct folders and paths

def checkconv(volflux, runnum):
    """This function just checks through each of the model-run writeouts and
    will print out the maximum step number, to spot-check whether or not the
    model run converged. 
    
    volflux: the current volcanic flux (km3/yr), a float
    runnum: a unique number for the output file name (e.g. 1)"""
    
    j=0                 #Counter for start looking for convergence number
    bottomline = 0      #Counter to find the last line of run printout from the model output
    outputfile = "Output" + `runnum` + "/" + str('%.2E' % volflux) + "/ModelOutput.txt"     #Goes into each of subdirectories to import the output file
    inputsfile = file(outputfile,'r')           #Reads the output files as a string
    lines = inputsfile.readlines()              #Breaks it into lines
    while j < len(lines):                       #Search the output for the the O2 flux line 
        if lines[j].find("Sprod:  O2") > 0 or lines[j].find("Sdep:   O2") > 0:     #This is the 2nd line after the model has ran. It may either say Sdep or Sprod, so check for both. 
            bottomline = j - O2 - 1             #Set the bottomline as above all the productions up to O2
        j += 1
    
    lastrun = lines[bottomline][0:8].strip()    #Strip away everything but the convergence number
    
    print "Converged at line: " + `lastrun`
    return int(lastrun)

def changetemp(diddle):
    """Used to change the temperature profile in the model.
    diddle: ground temperature (K) you want to use for global mean Mars temp."""
    
    indist = file('in.dist','r')   #import in.dist (temp profile stored here)
    lines = indist.readlines()
    indist.close()
    temps = lines
    
    i = 0
    while i < 50:
        T = diddle - 1.4*i     #This is the profile used in Zahnle et al. 2008
        repeat = T
        #print temps[440+i][4:18]
        others = temps[440+i][17:]
        write = "   " + str('%.8E' % T) + others
        temps[440+i] = write
        #print temps[440+i]
        i += 1
    
    while i < 110:
        others = temps[440+i][17:]
        #print temps[440+i][4:18]
        write = "   " + str('%.8E' % repeat) + others
        temps[440+i] = write
        #print temps[440+i]
        i += 1
    
    newtemps = "".join(temps)       #Join all the lines back together
    #print "###############"
    #print newtemps
    #print "Temperature: " + `diddle`
    outdist = file('in.dist','w')     #Defines output file, note that this apparently deletes all the original file so careful if you move it
    writetemps = newtemps      #Set the original lines file to the new file
    outdist.write(writetemps)     #Write the output file with the new file data
    outdist.close()      #Close the file

def modelrun(Fmagma, so2ratio, h2sratio, coratio, h2ratio, s2ratio, ch4ratio, incremental, iterations, runnum):
    """Runs the TOTCdev a number of times with incremental steps in increasing
    the volcanic flux (and the ratio between the SO2/H2S/CO/H2/CH4 species).
    
    Fmagma:      the starting point for the volcanic flux [scientific notation]
    so2ratio:    SO2 value based on Eqn. 4 in Sholes et al. 2017 (w/o V)
    h2sratio:    H2S value    ""     ""
    coratio:     CO value    ""     ""
    h2ratio:     H2 value    ""     ""
    s2ratio:     S2 value    ""     ""
    ch4ratio:    fixed flux for CH4 based on hypervelocity CO removal
    incremental: multiplier to ramp up flux (e.g. 1.5)
    iterations:  how many times you want to run the model (e.g. 30)
    runnum:      a unique integer for the output file names (e.g. 1)"""
    
    marstemp = 209.51   #Temp on Mars
    NSTEPS = 5000       #Max number of time steps
    volflux = Fmagma       #sets the base volcanic flux
    convint = 0                 #Counter for finding first non-converged run
    i = 0                       #used as the base increment for the loop, keep as zero
    
    changetemp(marstemp)   #change the temperature
    
    while i <= iterations:                                        #Start the loop for increasing the flux
        input_species = file('INPUTFILES/species.dat','r') #import the species.dat file for reading
        
        #THIS IS EDITING THE SPECIES.dat FILE
        lines = input_species.readlines()       #break species.dat into lines
        species = lines                         #did this to keep the data and delete the old version, i.e. it kept adding on the additional species 
        
        #Reset the fluxes
        h2sflux = volflux * h2sratio    #Initial H2S flux
        coflux = volflux * coratio      #Initial CO flux
        h2flux = volflux * h2ratio      #Initial H2 flux
        ch4flux = ch4ratio    #Initial CH4 flux
        s2flux = volflux * s2ratio      #Initial S2 flux
        so2flux = volflux * so2ratio    #Initial SO2 Flux
        
        #print species[17]       #Used for testing; This prints all the column labels
        
        ###To FIX: Eventually want to rework it so that the you can change the mbound conditions in species.dat, these are hardcoded in
        ch4print = str('%.3E' % ch4flux)
        species[18 + CH4] = species[18 + CH4][:53] + ch4print + ' ' + '20.     0      0.      0.        !Open up sfs_mars.py to change MBOUND stuff\n'
        print "CH4 Flux: ", ch4print
        
        h2sprint = str('%.3E' % h2sflux)                 #Converts the flux into a string in scientific notation with one decimal point
        species[18 + H2S] = species[18 + H2S][:53] + h2sprint + ' ' + '20.     0      0.      0.        !Open up sfs_mars.py to change MBOUND stuff \n'     #Rewrites the line [32] so that it has the right number of spaces
        print "H2S Flux: ", h2sprint                                      #Used for testing 
        
        coprint = str('%.3E' % coflux)                 #Converts the flux into a string in scientific notation with one decimal point
        species[18 + CO] = species[18 + CO][:53] + coprint + ' ' + '20.     2      -2.0E7  0.         !Open up sfs_mars.py to change MBOUND stuff \n'     #Rewrites the line so that it has the right number of spaces
        print "CO Flux: ", coprint                                        #Used for testing
        
        h2print = str('%.3E' % h2flux)                 #Converts the flux into a string in scientific notation with one decimal point
        species[18 + H2] = species[18 + H2][:53] + h2print + ' ' + '20.     0      0.      0.         !Open up sfs_mars.py to change MBOUND stuff\n'     #Rewrites the line so that it has the right number of spaces
        print "H2 Flux: ", h2print                                        #Used for testing 
        
        so2print = str('%.3E' % so2flux)                 #Converts the flux into a string in scientific notation with one decimal point
        species[18 + SO2] = species[18 + SO2][:53] + so2print + ' ' + '20.     0      0.      0.        !Open up sfs_mars.py to change MBOUND stuff\n'     #Rewrites the line so that it has the right number of spaces
        print "SO2 Flux: ", so2print                                        #Used for testing
        
        s2print = str('%.3E' % s2flux)
        species[18 + S2] = species[18 + S2][:53] + s2print + " " + '20.     0      0.      0.         !Open up sfs_mars.py to change MBOUND stuff\n'
        print "S2 Flux", s2print
        print "------------"        #Used for testing
        
        newspecies = "".join(species)       #Join all the lines back together
        output_species = file('INPUTFILES/species.dat','w')     #Defines output file, note that this apparently deletes all the original file so careful if you move it
        lines = newspecies      #Set the original lines file to the new file
        output_species.write(lines)     #Write the output file with the new file data
        output_species.close()      #Close the file
        
        #THIS IS FOR RUNNING THE MODEL
        filename = "mkdir Output" + `runnum` + "/" + str('%.2E' % volflux)      #Make a subdirectory for this particular run
        access = "Output" + `runnum` + "/" + str('%.2E' % volflux)      #Used this for ease below
        subprocess.call(filename, shell=True)
        
        #subprocess.call("make clean", shell=True)      #Makes the model clean, will take longer to run but use this if there are problems
        subprocess.call("make", shell=True)     #Make sure the TOTCdev is up to date
        print "#################################"
        print "Starting model-run " + `i + 1` + " / " + `iterations + 1` + " at:       " + strftime("%H:%M")
        print "Vol Flux:  " + str('%.2E' % volflux)
        subprocess.call("./TOTCdev > " + access + "/ModelOutput.txt", shell=True)        #Run the model and save the terminal output
        
        #THIS IS FOR COPYING THE FILES TO A NEW SUBDIRECTORY
        #Don't need most, so cleanup() removes them, but kept in case anyone
        #wants to save them in the future
        #subprocess.call("cp INPUTFILES/species.dat " + access, shell=True)
        #subprocess.call("mv out.trs " + access, shell=True)
        #subprocess.call("mv out.time " + access, shell=True)
        #subprocess.call("mv out.tim " + access, shell=True)
        #subprocess.call("mv out.terse " + access, shell=True)
        #subprocess.call("mv out.tau " + access, shell=True)
        #subprocess.call("mv out.so2 " + access, shell=True)
        #subprocess.call("mv out.redox " + access, shell=True)
        #subprocess.call("mv out.rates " + access, shell=True)
        #subprocess.call("mv out.raingc " + access, shell=True)
        #subprocess.call("mv out.rad " + access, shell=True)
        #subprocess.call("mv out.prod " + access, shell=True)
        #subprocess.call("mv out.params " + access, shell=True)
        subprocess.call("mv out.out " + access, shell=True)   #most important one
        #subprocess.call("mv out.NOprates " + access, shell=True)
        #subprocess.call("mv out.gridz " + access, shell=True)
        #subprocess.call("mv out.gridw " + access, shell=True)
        #subprocess.call("mv out.flux " + access, shell=True)
        #subprocess.call("mv out.flow " + access, shell=True)
        #subprocess.call("mv out.finalden " + access, shell=True)
        #subprocess.call("mv out.error " + access, shell=True)
        
        #UPDATE in.dist TO BETTER-CONVERGED VERSION
        subprocess.call("cp in.dist " + access + "/in.dist", shell=True)
        subprocess.call("cp out.dist in.dist", shell=True)
        subprocess.call("mv out.dist " + access, shell=True)
        
        #subprocess.call("mv out.densities " + access, shell=True)
        #subprocess.call("mv ISOinert.dist " + access, shell=True)
        #subprocess.call("mv ISOin.dist " + access, shell=True)
        #subprocess.call("mv out.cl " + access, shell=True)
        #subprocess.call("mv out.converge " + access, shell=True)
        #subprocess.call("mv out.xsec " + access, shell=True)
        #subprocess.call("mv out.oxygencolumnden " + access, shell=True)
        #subprocess.call("mv out.perchloratemr " + access, shell=True)
        #subprocess.call("mv out.pna " + access, shell=True)
        #subprocess.call("mv out.temprun " + access, shell=True)
        #subprocess.call("mv out.intrates " + access, shell=True)
        #subprocess.call("mv fort.13 " + access, shell=True)
        #subprocess.call("mv fort.22 " + access, shell=True)
        #subprocess.call("mv fort.20 " + access, shell=True)
        #subprocess.call("mv inertmixingratios.out " + access, shell=True)
        
        print "Done running model-run " + `i + 1` + " / " + `iterations + 1` + "  at:  " + strftime("%H:%M")
        lastrun = checkconv(volflux, runnum)
        print "#################################"
        
        if (lastrun >= NSTEPS) and (i > 1):              #Find the first non-converging line that isn't the first or second run
            if convint < 1:
                convint += 1
                global finalconv
                global finalvol
                global finalvolfloat
                finalconv = i - 1               #Don't want i+1 because i is from 0-iterations, and we need 0 (and the previous)
                finalvol = str('%.2E' % (volflux / incremental))
                finalvolfloat = volflux / incremental
        else:
            global finalso2float
            finalvolfloat = volflux
            finalconv = i - 1
            finalvol = "All converged"
        
        #INCREMENT THE STEP
        volflux *= incremental                         #Increment the H2S flux
        i += 1                                       #Increment the while loop
    
    pass

def modelplot(Fmagma, incremental, iterations, runnum, xminimum, xmaximum):
    """This process will plot the results from a model run using the same
    input as the modelrun process. It will plot and save the following graphs in
    the runnum output folder (full list above):
    
    legend70.png  : The main results with a legend
    results100.png: The main results with no lenged
    redoxstate.png: The redox state of the atm.
    normalizedmixingratios.png: Main results normalized to 1
    
    Fmagma:      the starting point for the volcanic flux [scientific notation]
    incremental: multiplier to ramp up flux (e.g. 1.5)
    iterations:  how many times you want to run the model (e.g. 30)
    runnum:      a unique integer for the output file names (e.g. 1)
    xminimum:    minimum volcanic flux for plotting
    xmaximum:    maximum volcanic flux for plotting"""
    
    #SET UP THE 'EMPTY' LISTS
    volflux = Fmagma                                       #This just sets the base SO2 flux to be changed
    global mixinglist                                           #Set globals so we can call them in other functions
    global fluxlist
    length = iterations + 1      
    mixinglist = numpy.ones(shape=[NSP1,length])*-9999      #Creates an array to be filled in with the mixing ratios for each species. -9999 is used to test for errors
    vollist = numpy.ones(shape=[length,1])*-9999            #Creates an array to be filled in with the used SO2 fluxes
    fluxlist = numpy.ones(shape=[NSP1,length])*-9999        #Creates array for ground level fluxes of species
    deplist = numpy.ones(shape=[NSP1,length])*-9999         #Creates array for ground deposition of species
    losslist = numpy.ones(shape=[NSP1,length])*-9999        #Creates array for top of atmosphere loss of species
    #NSP1 is the number of species, referenced in sfs_species.py
    
    redox = []  #(out - in) / in
    redox2 = [] # out - in
    
    #CHECK FOR HOW MANY BATCHES OF OUTPUTS THERE ARE
    #This is to check how many blocks of mixing ratios etc. there are
    if NSP < 72:
        #Only 4 batches
        mixing5 = ""
        depline5 = ""
        fluxline5 = ""
        lossline5 = ""
    if NSP < 54:
        #Only 3 batches
        mixing4 = ""
        depline4 = ""
        fluxline4 = ""
        lossline4 = ""
    if NSP < 36:
        #Only 2 batches
        mixing3 = ""
        depline3 = ""
        fluxline3 = ""
        lossline3 = ""
    
    #RUN THE LOOP TO EXTRACT DATA FROM OUT.OUT
    i = 0                   #Used for the while loop for finding folders
    resultsection = 0       #The line where we find the true results once the model reaches a steady state, set to zero before we find it
    while i <= iterations:
        outputfile = "Output" + `runnum` + "/" + str('%.2E' % volflux) + "/out.out"     #Goes into each of subdirectories to import the output file
        inputspecies = file(outputfile,'r')                                             #Reads the output files as a string
        text = inputspecies.readlines()                                                 #Breaks the output files into lines
        
        #*************************************
        #FIND THE MIXING RATIOS
        aa = 0           #Used for the while loop
        while aa < len(text):        #This loop finds all TIMEY = ***, the last of which shows the actual ressults
            if text[aa].find("NPHOT") > 0:      #I'm using this to judge where the final output lines are. Should be good, but may need to hardcode a better tracer
                resultsection = aa
            aa += 1
        
        bb = resultsection     #Going to now find the first mixing ratios
        while bb < len(text):        #Need to find the Z after NPHOT which shows the mixingratios
            if text[bb].find("Z") > 0:       #Goes through each line in the text starting at the results section to find a Z
                mixing1 = text[bb + 1]
                mixing1 = mixing1[11:].strip("\n")
                cc = bb + 1                     #Go to the next line to start the next search
                bb = len(text) + 1       #If it finds one set the counter to larger than the text length to shut off the loop
            bb +=1
        
        while cc < len(text):   #Find the next batch of mixing ratios
            if text[cc].find("Z") > 0:
                mixing2 = text[cc + 1]
                mixing2 = mixing2[11:].strip("\n")
                dd = cc + 1
                cc = len(text) + 1
            cc += 1
        
        if NSP > 35:        #These are set to batches of 18 species per line of output. If there are more species, then find them. 
            while dd < len(text):   #Find the next batch of mixing ratios
                if text[dd].find("Z") > 0:
                    mixing3 = text[dd + 1]
                    mixing3 = mixing3[11:].strip("\n")
                    ee = dd + 1
                    dd = len(text) + 1
                dd += 1
        
        if NSP > 53:
            while ee < len(text):   #Find the next batch of mixing ratios
                if text[ee].find("Z") > 0:
                    mixing4 = text[ee + 1]
                    mixing4 = mixing4[11:].strip("\n")
                    ff = ee + 1
                    ee = len(text) + 1
                ee += 1
        
        if NSP > 71:
            while ff < len(text):   #Find the next batch of mixing ratios
                if text[ff].find("Z") > 0:
                    mixing5 = text[ff + 1]
                    mixing5 = mixing5[11:].strip("\n")
                    gg = ff + 1
                    ff = len(text) + 1
                ff += 1
        
        mixingratios = mixing1 + " " + mixing2 + " " + mixing3 + " " + mixing4 + " " + mixing5      #Adds all the mixing ratio species into one line
        
        #SEPARATE THE VALUES IN THE MIXING OUTPUT
        #Need to do this b/c values don't always have a separator
        rrr = ''
        sss = 0
        ttt = 8     #The length of the values in scientific notation
        while sss <= len(mixingratios) + 1:
            rrr = rrr + ' ' + mixingratios[sss:ttt]
            sss = ttt
            ttt += 9
        zempty = []
        zlist = rrr.split()
        
        zfull = []
        for zitem in zlist:         #This loop is to set values really small if they are less than 9.99E-99
            zitem = zitem.strip()   #This is because it will write anything smaller without the 'E'
            if zitem.find('E') < 1: #e.g. 1.00-100
                zitem = '1.00E-99'
            zfull.append(zitem)
            
        mixingratios = " ".join(str(zitem) for zitem in zfull)
        array = numpy.fromstring(mixingratios, dtype = float, sep=" ")  #Create a numpy array from the mixing ratios (in float because it can't read scientific notation otherwise)
        mixinglist[:,i]=array     #Add the mixing ratios for this species at this SO2 flux to he array
        vollist[i,0] = volflux      #Add the SO2 flux value into its array
        
        #*************************************
        #FIND THE FLUXES
        
        if NSP < 72:            #Do this in case the variable gg is not made because the if loop was not initiated
            gg = resultsection
        
        while gg < len(text):
            if text[gg].find("FLUXES OF LONG-LIVED SPECIES") > 0:
                hh = gg + 1
                gg = len(text) + 1
            gg += 1
        
        while hh < len(text):               #Find the first batch of fluxes, after the mixing ratios
            if text[hh].find("Z") > 0:
                fluxline1 = text[hh + 1]
                fluxline1 = fluxline1[11:].strip("\n")
                jj = hh + 1
                hh = len(text) + 1
            hh += 1
        
        while jj < len(text):               #Find the second batch of fluxes
            if text[jj].find("Z") > 0:
                fluxline2 = text[jj + 1]
                fluxline2 = fluxline2[11:].strip("\n")
                kk = jj + 1
                jj = len(text) + 1
            jj += 1
        
        if NSP > 35:
            while kk < len(text):               #Find the third batch of fluxes
                if text[kk].find("Z") > 0:
                    fluxline3 = text[kk + 1]
                    fluxline3 = fluxline3[11:].strip("\n")
                    ll = kk + 1
                    kk = len(text) + 1
                kk += 1
        
        if NSP > 53:
            while ll < len(text):               #Find the fourth batch of fluxes
                if text[ll].find("Z") > 0:
                    fluxline4 = text[ll + 1]
                    fluxline4 = fluxline4[11:].strip("\n")
                    mm = ll + 1
                    ll = len(text) + 1
                ll += 1
        
        if NSP > 71:
            while mm < len(text):               #Find the fifth batch of fluxes
                if text[mm].find("Z") > 0:
                    fluxline5 = text[mm + 1]
                    fluxline5 = fluxline5[11:].strip("\n")
                    nn = mm + 1
                    mm = len(text) + 1
                mm += 1
        
        fluxes = fluxline1 + " " + fluxline2 + " " + fluxline3 + " " + fluxline4 + " " + fluxline5
        
        #TAKES THE ABSOLUTE VALUES FOR THE FLUXES
        #Can't plot (-) values in log space so take the abs()
        #This basically takes the string, converts to list, evaluates, then back to string
        fempty = []
        fother = []
        flist = fluxes.split()
        for fitem in flist:
            if fitem.find('E') < 1:
                fitem = '1.00E-99'      #If the value is too small (<1E-100) then just call it 1e-99
            #Ideally I would want to add in something that checks to see if the flux is too great
            #i.e. 3.1+140 and print out an error statement but adding fitem.find('-')>1
            #Doesn't seem to work when there is multiple negatives. Needs fixing eventually.
            fempty.append(abs(float(fitem)))

        ffull = " ".join(str(gitem) for gitem in fempty)
        
        #print str('%.1E' % so2flux) + "   " + str('%.1E' % fempty[S8])     #Used for testing, prints out the so2flux and corresponding flux for S8
        
        array2 = numpy.fromstring(ffull, dtype=float, sep=" ")
        fluxlist[:,i] = array2
        
        #*************************************
        #FIND THE DEPOSITIONS AND LOSSES
        #Should be able to consolidate the flux section with this section (same info is included)
        
        if NSP < 72:            #Do this in case nn is not defined because of if loop not initiated
            nn = resultsection
        
        while nn < len(text):                   #Find where the deposition and loss values are
            if text[nn].find("PHIDEP") > 0:
                oo = nn + 1
                nn = len(text) + 1
            nn += 1
        
        while oo < len(text):
            if text[oo].find("Z") > 0:  #Find the first batch
                #Note: [hh + 0] is the labels
                #       [hh + 1] are the rainout rates
                #       [hh + 2] are the depositions
                #       [hh + 3] are the losses at the top
                #       [hh + 4] are the lower boundary conditions (0-3)
                #       [hh + 5] are the deposition velocities
                #       [hh + 7] are the total production (column integrated)
                #       [hh + 8] are the total loss (column integrated)
                #       [hh + 9] are the flux at the upper boundary
                #       [hh + 10] are the flux at the lower boundary
                #       [hh + 11] are the convergence (should be relatively small)
                depline1 = text[oo + 2]                 #These are the depositions
                depline1 = depline1[11:].strip("\n")
                lossline1 = text[oo + 9]                #These are the losses
                lossline1 = lossline1[11:].strip("\n")
                pp = oo + 1
                oo = len(text) + 1
            oo += 1
        
        while pp < len(text):       #Find the second batch
            if text[pp].find("Z") > 0:
                depline2 = text[pp + 2]
                depline2 = depline2[11:].strip("\n")
                lossline2 = text[pp + 9]
                lossline2 = lossline2[11:].strip("\n")
                qq = pp + 1
                pp = len(text) + 1
            pp += 1
        
        if NSP > 35:
            while qq < len(text):       #Find the third batch
                if text[qq].find("Z") > 0:
                    depline3 = text[qq + 2]
                    depline3 = depline3[11:].strip("\n")
                    lossline3 = text[qq + 9]
                    lossline3 = lossline3[11:].strip("\n")
                    rr = qq + 1
                    qq = len(text) + 1
                qq += 1
        
        if NSP > 53:
            while rr < len(text):       #Find the fourth batch
                if text[rr].find("Z") > 0:
                    depline4 = text[rr + 2]
                    depline4 = depline4[11:].strip("\n")
                    lossline4 = text[rr + 9]
                    lossline4 = lossline4[11:].strip("\n")
                    ss = rr + 1
                    rr = len(text) + 1
                rr += 1
        
        if NSP > 71:
            while ss < len(text):       #Find the fifth batch
                if text[ss].find("Z") > 0:
                    depline5 = text[ss + 2]
                    depline5 = depline5[11:].strip("\n")
                    lossline5 = text[ss + 9]
                    lossline5 = lossline5[11:].strip("\n")
                    tt = ss + 1
                    ss = len(text) + 1
                ss += 1
        
        deps = depline1 + " " + depline2 + " " + depline3 + " " + depline4 + " " + depline5
        
        #SEPARATE THE VALUES IN THE DEPS OUTPUT
        #Need to do this b/c values don't always have a separator
        xxx = ''
        yyy = 0
        zzz = 8
        while yyy <= len(deps) + 1:
            xxx = xxx + ' ' + deps[yyy:zzz]
            yyy = zzz
            zzz += 9
        yempty = []
        ylist = xxx.split()
        for yitem in ylist:
            yitem = yitem.strip()
            if yitem.find('E') < 1:
                yitem = '1.00E-99'      #If it's too small set it to this
            yempty.append(yitem)
        deps = " ".join(str(yitem) for yitem in yempty)
        
        array3 = numpy.fromstring(deps, dtype=float, sep=" ")
        deplist[:,i] = array3
        
        loss = lossline1 + " " + lossline2 + " " + lossline3 + " " + lossline4 + " " + lossline5
        
        #SEPARATE THE VALUES IN THE LOSS OUTPUT
        #Need to do this b/c values have to common separator
        ooo = ''
        ppp = 0
        qqq = 8
        while ppp <= len(loss) + 1:
            ooo = ooo + ' ' + loss[ppp:qqq]
            ppp = qqq
            qqq += 9
        sempty = []
        slist = ooo.split()
        for sitem in slist:
            sitem = sitem.strip()
            if sitem.find('E') < 1:
                sitem = '1.00E-99'      #If it's too small set it to this
            sempty.append(sitem)
        loss = " ".join(str(sitem) for sitem in sempty)
        
        array4 = numpy.fromstring(loss, dtype=float, sep=" ")
        losslist[:,i] = array4
        
        #**************************************
        #  FIND THE REDOX CONSERVATION VALUES
        tt = resultsection
        while tt < len(text):                   #Find where the deposition and loss values are
            if text[tt].find("redox budget:") > 0:
                vv = tt + 1
                tt = len(text) + 1
            tt += 1
        
        redoxbudget = text[vv]
        o_in = float(redoxbudget[14:27])
        o_out = float(redoxbudget[39:52])
        r_in = float(redoxbudget[62:75])
        r_out = float(redoxbudget[86:99])
        r_rain = float(redoxbudget[111:121])
        o_rain = float(redoxbudget[133:143])
        red_in = r_in - o_in
        red_out = r_out - o_out + r_rain - o_rain
        redox_conv = (red_out - red_in) / red_in
        redox.append(redox_conv)
        redox2.append(red_out - red_in)
        
        #INCREMENT THE LOOPS
        volflux *= incremental      #Increase the volanic flux the same way the directories were set up
        i += 1
    
    print redox
    
    #Testing section
    #print "\nMixing Ratios at " + str('%.1E' % so2baseflux) + ":"
    #print "CO: %.2f" % (mixinglist[CO,0] * 100) + " %"
    #print "O2: %.2e" % mixinglist[O2,0]
    #print "Total:   " + `sum(mixinglist[:,0]) + 0.95 + 0.016`
    #print "Mixing Ratios at " + finalso2 + ":"
    #print "CO: %.2f" % (mixinglist[CO,finalconv] * 100) + " %"
    #print "O2: %.2e" % mixinglist[O2,finalconv]
    #print "Total:   " + `sum(mixinglist[:,finalconv]) + 0.95 + 0.016` + "\n"
    #
    #*************************************
    #     NOW FOR PLOTTING THE DATA
    #*************************************
    
    #For now I'm plotting mixing ratios with . and fluxes/losses with -- and deps with -
    #These colors exactly match up with the plotting colors of Kevin's unpublished paper
    mixcolors = ["#000000", "#008011", "#DD0806", "#02ABEA", "#BC9A05", "#999999", "#FF69B4"]   #Has pink as backup
    fluxcolors = ["#9816D4", "#FF8B1B"]
    depcolors = ["#0000D4", "#27FC20", "#02ABEA", "#FF69B4", "#FCF305"]    #Has yellow as backup
    losscolors = ["#008011"]
    mixmarkers = ["o", "s", "D", "^", "v", "d", "p"]  #o is circle, s is square, D is diamond, ^ is up triangle, v is down triangle, has + as backup
    
    #The species you want to include 
    mixspecies = ["$\mathregular{O_2}$", "$\mathregular{H_2}$", "CO", "$\mathregular{SO_2}$", "OCS", "$\mathregular{CH_4}$"]
    mixnums = [O2, H2, CO, SO2, OCS, CH4]
    fluxspecies = ["$\mathregular{{H_2SO_4}_{aer}}$ flux", "$\mathregular{S_{{8}_{aer}}}$ flux"]
    fluxnums = [SO4AER, S8AER]
    depspecies = ["$\mathregular{H_2O_2}$ dep", "$\mathregular{O_3}$ dep", "$\mathregular{SO_2}$ dep", "SO dep"]
    depnums = [H2O2, O3, SO2, SO]   #sfs removed SO
    lossspecies = ["$\mathregular{H_2}$ loss"]
    lossnums = losslist[H2]
    
    #CHANGE YOUR PLOT BOUNDARIES HERE
    xmin = 0.001
    xmax = 1.0
    y1min = 1e-7
    y1max = 1e0
    y2min = 5e6
    y2max = 2e10
    
    counter = 0
    scounter = 0
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    ax1.set_ylabel("Mixing Ratio")
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(y1min, y1max)
    ax1.set_xlim(xmin, xmax)
    
    #Just plotting a line to show modern Earth values
    #earthso2flux = 2.0e9        #This is to create a vertical line showing SO2 flux on Earth for comparison
    #ax1.plot([earthso2flux, earthso2flux], [y1min, y1max], color = '#888888', linestyle='-', linewidth=2)
    #ax1.text(earthso2flux - earthso2flux*0.2, y1min - y1min*0.45, "Earth", rotation=0)
    
    #Plot mixing ratioson the y1 axis
    for zz in mixnums:
        ax1.plot(vollist[:], mixinglist[zz,:], color=mixcolors[scounter], marker=mixmarkers[scounter], linestyle='none', markersize=5.0, label=mixspecies[scounter])
        scounter += 1
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("Fluxes [molecules $\mathregular{cm^{-2} s^{-1}}$]")
    ax2.set_yscale('log')
    ax2.set_ylim(y2min, y2max)
    ax2.set_xlim(xmin, xmax)
    
    #Plot fluxes on the y2 axis
    scounter = 0
    for yy in fluxnums:
        ax2.plot(vollist[:], fluxlist[yy,:], color=fluxcolors[scounter], lw=2.5, label=fluxspecies[scounter])
        scounter += 1
    
    #Plot deposition on the y2 axis
    scounter = 0
    for xx in depnums:
        ax2.plot(vollist[:], deplist[xx,:], color=depcolors[scounter], lw=2.5, label=depspecies[scounter])
        scounter += 1
    
    #Plot loss on the y2 axis
    ax2.plot(vollist[:], lossnums, color=losscolors[0], lw=2.5, label=lossspecies[0])
    
    plt.savefig("Output" + `runnum` + "/results100.png")       #Save the figure, was designed for without the legend
    
    ax1.set_position([0.1,0.1,0.6,0.8])     #Used to set the left, bottom, width, height of the actual plot. Used to make room for the legend.
    ax2.set_position([0.1,0.1,0.6,0.8])     #Same as above but for the other axes
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    ax2.legend(loc='lower left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    plt.savefig("Output" + `runnum` + "/legend70.png")       #Save the figure with legends
    #plt.show()      #Show the plot
    
    plt.close()
    
    #********************************
    #PLOT JUST IMPORTANT MIXING RATIOS
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(y1min, y1max)
    ax1.set_ylabel("Mixing Ratio")
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    counter = 0
    scounter = 0
    for zz in mixnums:
        ax1.plot(vollist[:], mixinglist[zz,:], color=mixcolors[scounter], linestyle='-', linewidth=3.5, label=mixspecies[scounter])
        scounter += 1
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.01,1.01), fancybox=True, shadow=True)
    plt.savefig("Output" + `runnum` + "/mixingratios.png")
    plt.close()
    
    #********************************
    #PLOT JUST IMPORTANT FLUXES
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(y2min, y2max)
    ax1.set_ylabel("Fluxes [molecules $\mathregular{cm^{-2} s^{-1}}$]")
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    counter = 0
    scounter = 0
    for yy in fluxnums:
        ax1.plot(vollist[:], fluxlist[yy,:], color=fluxcolors[scounter], lw=3.5, label=fluxspecies[scounter])
        scounter += 1
    
    #Plot deposition on the y2 axis
    scounter = 0
    for xx in depnums:
        ax1.plot(vollist[:], deplist[xx,:], color=depcolors[scounter], lw=3.5, label=depspecies[scounter])
        scounter += 1
    
    #Plot loss on the y2 axis
    ax1.plot(vollist[:], lossnums, color=losscolors[0], lw=3.5, label=lossspecies[0])
    
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.01,1.01), fancybox=True, shadow=True)
    plt.savefig("Output" + `runnum` + "/fluxes.png")
    plt.close()
    
    #********************************
    #PLOT THE REDOX STATE OF THE ATM.
    #pOx = 2pO2 - pCO - pH2
    #Look back to pg.4 in Kevin's unpublished paper
    #Modern Mars has pOx approx. 10 ubar
    plt.figure()
    
    plt.xscale('log')
    plt.xlim(xmin, xmax)
    plt.ylim(-750, 100)
    plt.ylabel("pOx [$\mathregular{\mu}$bar]")
    plt.xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    
    #Just plotting a line to show modern Earth values
    #plt.axvline(x=earthso2flux, ymin=0, ymax=1, color = 'g', linestyle='-', linewidth=2)
    
    plt.plot([xmin, xmax], [0, 0], color = '#888888', linestyle='-', linewidth=2)  #Plot the redox neutral boundary line
    pOx = (2*mixinglist[O2,:]-mixinglist[CO,:]-mixinglist[H2,:]-3*mixinglist[OCS,:])*6500     #Multiply by 6500 to get partial pressures in ubar
    plt.plot(vollist[:], pOx, color='k', linewidth=2)
    
    plt.savefig("Output" + `runnum` + "/redoxstate.png")
    #plt.show()
    plt.close()
    
    testytest = redoxswitch(mixinglist, vollist)
    
    #*******************************
    # Redox and mixing ratios
    
    normalizedmr = numpy.copy(mixinglist)
    co2list = numpy.ones(shape=[iterations, 1])*0.95
    
    font = 14
    
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(y1min, y1max)
    ax1.set_ylabel("Mixing Ratio", fontsize=font)
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]", fontsize=font)
    plt.tick_params(labelsize=font)
    
    mrsums = numpy.sum(normalizedmr, axis = 0)
    for item in range(0,len(mrsums)):
        mrsums[item] += (0.95 + 0.016)
    for item in range(0,len(co2list)):
        co2list[item] /= mrsums[item]
    for item in range(0, len(normalizedmr)):
        for element in range(0, len(normalizedmr[item])):
            normalizedmr[item][element] /= mrsums[element]
    
    mixcolors2 = ['k','r','g']
    mixmarkers2 = ['o','d','s']
    mixspecies2 = ["O2","CO","H2"]
    scounter2 = 0
    for uuu in [O2,CO,H2]:
        ax1.plot(vollist[:], normalizedmr[uuu,:], color=mixcolors2[scounter2], marker=mixmarkers2[scounter2], linestyle='--', markersize=5.0)
        scounter2 += 1
    
    plt.plot([testytest, testytest], [y1min, y1max], color = '#888888', linestyle='-', linewidth=2)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Redox State [$\si{\micro \meter}]", fontsize=font)
    ax2.set_ylim(-300, 100)
    ax2.set_xlim(xmin, xmax)
    ax2.tick_params(labelsize=font)
    
    plt.plot(vollist[:], pOx, color='k', linewidth=2)
    plt.plot([xmin, xmax], [0, 0], color = '#888888', linestyle='-', linewidth=2)  #Plot the redox neutral boundary line
    
    ax1.set_position([0.1,0.1,0.6,0.8])     #Used to set the left, bottom, width, height of the actual plot. Used to make room for the legend.
    ax2.set_position([0.1,0.1,0.6,0.8])     #Same as above but for the other axes
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    ax2.legend(loc='lower left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    
    plt.savefig("Output" + `runnum` + "/RedoxFigPaper.png")
    plt.close()
    
    #********************************
    #PLOT ALL MIXING RATIOS
    olist = [O, O2, H2O, H, OH, HO2 ,H2O2, H2]
    clist = [CO, HCO, H2CO, CH4, CH3, C2H6, NO, NO2, HNO, HNO3, N, HNO4]
    slist = [H2S, HS, S, S2, S3, S4, S8AER, SO, HSO, SO4AER, SO3, OCS, SO2]
    oleg = ["O", "O2", "H2O", "H", "OH", "HO2", "H2O2", "H2"]
    cleg = ["CO", "HCO", "H2CO", "CH4", "CH3", "C2H6", "NO", "NO2", "HNO", "HNO3", "N", "HNO4"]
    sleg = ["H2S", "HS", "S", "S2", "S3", "S4", "S8AER", "SO", "HSO", "SO4AER", "SO3", "OCS", "SO2"]
    plt3colors = ["b", "g", "r", "c", "m", "k", "b", "g", "r", "c", "m", "k", "b"]
    plt3styles = ["-", "-", "-", "-", "-", "-", ":", ":", ":", ":", ":", ":", "--"]
    
    #Plot O species
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
#    ax1.set_ylim(1e-15, 1e-2)
    ax1.set_ylabel("Mixing Ratio")
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    counter = 0
    for item in olist:
        ax1.plot(vollist[:], mixinglist[item,:], color=plt3colors[counter], linestyle=plt3styles[counter], label=oleg[counter])
        counter += 1
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.01,1.01), fancybox=True, shadow=True)
    plt.savefig("Output" + `runnum` + "/Ospecies.png")
    plt.close()
    
    #Plot C and N species
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
#    ax1.set_ylim(1e-26, 1e-2)
    ax1.set_ylabel("Mixing Ratio")
    ax1.set_xlabel("Volcanic Magma Flux [molecules $\mathregular{cm^{-2} s^{-1}}$]")
    counter = 0
    for item in clist:
        ax1.plot(vollist[:], mixinglist[item,:], color=plt3colors[counter], linestyle=plt3styles[counter], label=cleg[counter])
        counter += 1
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.01,1.01), fancybox=True, shadow=True)
    plt.savefig("Output" + `runnum` + "/Cspecies.png")
    plt.close()
    
    #Plot S species
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
#    ax1.set_ylim(1e-40, 1e-7)
    ax1.set_ylabel("Mixing Ratio")
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    counter = 0
    for item in slist:
        ax1.plot(vollist[:], mixinglist[item,:], color=plt3colors[counter], linestyle=plt3styles[counter], label=sleg[counter])
        counter += 1
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.01,1.01), fancybox=True, shadow=True)
    plt.savefig("Output" + `runnum` + "/Sspecies.png")
    plt.close()
    
    #********************************
    #        TEST PLOT AREA
    #Use this area to test any certain species you want to test
    #Currently Testing:     S8aer Fluxes
    plt.figure()
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(xmin, xmax)
    plt.ylim(1e1, 1e9)
    plt.ylabel("Fluxes")
    plt.xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]")
    
    #Just plotting a line to show modern Earth values
    #plt.axvline(x=earthso2flux, ymin=0, ymax=1, color = 'g', linestyle='-', linewidth=2)
    plt.plot(vollist[:], fluxlist[SO4AER], color='k', linewidth=2, label="H2SO4aer")
    plt.plot(vollist[:], fluxlist[S8AER], color='b', linewidth=2, label="S8aer")
    
    plt.savefig("Output" + `runnum` + "/AER_fluxes.png")
    plt.close()
    
    #*******************************
    #    NORMALIZED MIXING RATIOS
    #This is just a test to normalize the mixing ratios so they add up to 1
    #Not sure how realistic this is, will think about it
    
    normalizedmr = numpy.copy(mixinglist)
    co2list = numpy.ones(shape=[iterations, 1])*0.95
    
    font = 14
    
    fig, ax1 = plt.subplots()
    ax1.set_position([0.1,0.1,0.7,0.8])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(y1min, y1max)
    ax1.set_ylabel("Mixing Ratio", fontsize=font)
    ax1.set_xlabel("Volcanic Magma Flux [$\mathregular{km^{3} yr^{-1}}$]", fontsize=font)
    plt.tick_params(labelsize=font)
    
    mrsums = numpy.sum(normalizedmr, axis = 0)
    for item in range(0,len(mrsums)):
        mrsums[item] += (0.95 + 0.016)
    for item in range(0,len(co2list)):
        co2list[item] /= mrsums[item]
    for item in range(0, len(normalizedmr)):
        for element in range(0, len(normalizedmr[item])):
            normalizedmr[item][element] /= mrsums[element]
    
    scounter = 0
    for uu in mixnums:
        ax1.plot(vollist[:], normalizedmr[uu,:], color=mixcolors[scounter], marker=mixmarkers[scounter], linestyle='--', markersize=5.0, label=mixspecies[scounter])
        scounter += 1
    
    #ax1.plot(vollist[:], co2list[:], color=mixcolors[scounter], marker=mixmarkers[scounter], linestyle='none', markersize=5.0, label="CO2")
    
    plt.plot([testytest, testytest], [y1min, y1max], color = '#888888', linestyle='-', linewidth=2)
    #plt.text(earthso2flux - earthso2flux*0.2, y1min - y1min*0.45, "Earth", rotation=0)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Fluxes [molecules $\mathregular{cm^{-2} s^{-1}}$]", fontsize=font)
    ax2.set_yscale('log')
    ax2.set_ylim(y2min, y2max)
    ax2.set_xlim(xmin, xmax)
    ax2.tick_params(labelsize=font)
    
    scounter = 0
    for yy in fluxnums:
        ax2.plot(vollist[:], fluxlist[yy,:], color=fluxcolors[scounter], lw=2.5, label=fluxspecies[scounter])
        scounter += 1
    scounter = 0
    for xx in depnums:
        ax2.plot(vollist[:], deplist[xx,:], color=depcolors[scounter], lw=2.5, label=depspecies[scounter])
        scounter += 1
    ax2.plot(vollist[:], lossnums, color=losscolors[0], lw=2.5, label=lossspecies[0])
    ax1.set_position([0.1,0.1,0.6,0.8])     #Used to set the left, bottom, width, height of the actual plot. Used to make room for the legend.
    ax2.set_position([0.1,0.1,0.6,0.8])     #Same as above but for the other axes
    ax1.legend(loc='upper left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    ax2.legend(loc='lower left', ncol=1, bbox_to_anchor = (1.17,0.5), fancybox=True, shadow=True, prop={'size':10})
    
    #plt.gcf().subplots_adjust(bottom=0.15)
    
    plt.savefig("Output" + `runnum` + "/NormalizedMixingRatios")
    plt.close()
    
    #testytest is the magma flux at which the atmosphere switches
    return testytest

def convwrite(Fmagma, incremental, iterations, runnum):
    """THIS VERSION IS FOR WRITING THE README FILE
    
    This function just checks through each of the model-run writeouts and
    will print out the maximum step number, to check whether or not the
    model run converged.
    
    Fmagma:      the starting point for the volcanic flux [scientific notation]
    incremental: multiplier to ramp up flux (e.g. 1.5)
    iterations:  how many times you want to run the model (e.g. 30)
    runnum:      a unique integer for the output file names (e.g. 1)"""
    
    i = 0
    j = 0
    volflux = Fmagma
    bottomline = 0
    while i <= iterations:
        outputfile = "Output" + `runnum` + "/" + str('%.2E' % volflux) + "/ModelOutput.txt"     #Goes into each of subdirectories to import the output file
        inputsfile = file(outputfile,'r')           #Reads the output files as a string
        lines = inputsfile.readlines()              #Breaks it into lines
        while j < len(lines):
            if lines[j].find("Sprod:  O2") > 0 or lines[j].find("Sdep:   O2") > 0:     #This is the 2nd line after the model has ran. It may either say Sdep or Sprod, so check for both. 
                bottomline = j - O2 - 1
            j += 1
        
        lastrun = lines[bottomline][0:8].strip()
        openreadme = "Output" + `runnum` + "/" + "00_readme.txt"
        openreadmefile = file(openreadme, 'a')
        openreadmefile.write(str('%.2E' % volflux) + "     " + lastrun + "\n")
        openreadmefile.close()
        
        j = 0       #You need this here or else it will present a list index error
        volflux *= incremental
        i += 1

def redoxswitch(mixinglist, vollist):
    """This is used to find where the redox state of the atmospheres switches
    from being oxic to anoxic. It is based on the equation pOx = 2pO2 - pCO - pH2.
    
    To do this is calculates the redox state and finds when it goes from being
    positive to being negative. It then calculates a simple line between the two
    points and finds the volcanic flux level where that line is zero.
    
    mixinglist: this is the list of mixing ratios that is created in modelplot()
    vollist   : This is the list of volcanic flux values being plotted"""
    
    pOx = (2*mixinglist[O2,:]-mixinglist[CO,:]-mixinglist[H2,:])*6500       #Multiply by 6500 to get values in ubar
    print "Starting Redox Switch Calculations"
    i = 0
    j = 0
    while i < len(pOx):
        print "pOx: ", pOx[i]
        if pOx[i] < 0:
            j = i
            i = len(pOx) + 1
        i += 1
    
    #Simple line building
    m = (pOx[j] - pOx[j-1]) / (vollist[j] - vollist[j-1])
    b = pOx[j-1] - m*vollist[j-1]
    
    print "Slope: ", m[0]
    print "b    : ", b[0]
    print "pOx Before: ", pOx[j-1]
    print "pOx After : ", pOx[j]
    print "SO2 After : ", vollist[j]
    print " y = " + `m[0]` + " x + " + `b[0]`
    
    eqn = [m[0], b[0]]
    root = numpy.roots(eqn)
    answer = root[0]
    
    print "SWITCH VALUE:  ", str('%.3e' % answer)
    return answer

def multirun(so2_i, h2_i, h2s_i, s2_i, co_i, runnum, notes):
    """Calculates Eqn. 4 in Sholes et al. 2017 (without V) and runs it through
    modelrun(). 
    
    so2_i:     SO2 input value for (from Gaillard et al. 2013)
    h2_i:      H2 input    ""    ""
    h2s_i:     H2S input    ""    ""
    s2_i:      S2 input    ""    ""
    CO_i:      CO input    ""    ""
    runnum:    a unique integer for the output file names (e.g. 1)
    notes:     user-inputted notes for each model run"""
    
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!!!!!!!!!!STARTING!!!!!!!!!!!"
    print "Running TOTCdev for varying levels of volcanic gas fluxes"
    
    start_time = time()   #for timing how long it takes
    xmaximum = 2.0        #maximum volcanic flux for testing (km3/yr)
    xminimum = 0.001      #minimum volcanic flux for testing (km3/yr)
    increments = 1.5      #multiplier to ramp up flux (e.g. 1.5)
    
    #Calculates how many iterations are required to ramp up volcanic flux by increments
    iterations = int(numpy.log10((xmaximum)/(xminimum))/numpy.log10(increments))+1
    
    volden =2.9e15         #g/km3
    volnorm = xminimum      #km3/yr
    avogadro = 6.022e23     #molecules/mol
    speryr = 3.2e-7         #yr/s
    samars = 1.4e18         #cm2
    ppm = 1e-6
    
    #Get the values in ppmv (Gaillard has them in ppm wt)
    so2_g = (so2_i / wtSO2) * ppm
    h2_g =  (h2_i / wtH2) * ppm
    h2s_g = (h2s_i / wtH2S) * ppm
    s2_g =  (s2_i / wtS2) * ppm
    co_g =  (co_i / wtCO) * ppm
    
    #Flux of A = ppmvA * 1e-6 * Fmagma
    #They will be multiplied by Fmagma in modelrun()
    h2ratio = h2_g * volden * avogadro * speryr / samars
    h2sratio = h2s_g * volden * avogadro * speryr / samars
    so2ratio = so2_g * volden * avogadro * speryr / samars
    s2ratio = s2_g * volden * avogadro * speryr / samars
    coratio = co_g * volden * avogadro * speryr / samars
    ch4ratio = 1.0e7   #Fixed (see Appendix B in Sholes et al. 2017)
    
    subprocess.call("mkdir " + "Output" + `runnum`, shell=True)     #Make a new directory to store the model runs
    
    #incopy.dist is a well-converged list of mixing ratios, whenever you start modeling
    #make sure to replace the latest copy of in.dist with a copy of incopy.dist renamed
    #this is b/c in.dist is changes every time to make it more well-convereged
    #but this doesn't help when you go from 2 km/3 to 0.0001 km3/yr
    subprocess.call("rm in.dist", shell=True)
    subprocess.call("cp incopy.dist in.dist", shell=True)   #incopy.dist is a v
    modelrun(volnorm, so2ratio, h2sratio, coratio, h2ratio, s2ratio, ch4ratio, increments, iterations, runnum)
    elapsed_time = time() - start_time
    print ("*** Elapsed Time:  %.2f min ***" % ((time() - start_time)/60.))
    elapsed_time_str = "%.2f min" % ((time() - start_time)/60.)
    
    volcswitch = modelplot(volnorm, increments, iterations, runnum, xminimum, xmaximum)
    
    #PRINT OUT A README FILE WITH THE PARAMETERS TO REFER BACK TO
    text_readme = open("Output" + `runnum` + "/00_readme.txt", "w")
    printtext = "Parameters Used for Output Run " + `runnum` + ':\n' + \
                notes + \
                "\n\nGaillard et al. 2009 Values: " + \
                "\nH2:  " + `h2_g` + "\nH2S: " + `h2s_g` + "\nSO2: " + `so2_g` + \
                "\nS2:  " + `s2_g` + "\nCO:  " + `co_g` + \
                "\n\nMagma Base Flux: " + str('%.2E' % volnorm) + \
                "\nH2S Ratio: " + str('%.3f' % h2sratio) + "\nH2 Ratio:  " + \
                str('%.3f' % h2ratio) + "\nCO Ratio:  " + str('%.3f' % coratio) + \
                "\nCH4 Ratio: " + str('%.3f' % ch4ratio) + \
                "\nS2 Ratio:  " + str('%.3f' % s2ratio) + \
                "\n**Flux is in molecules/cm2/s" + \
                "\nNumber of runs:      " + `iterations + 1` + "\nMagma Flux Multiplier: " + \
                `increments` + "\n\nMixing Ratios at " + str('%.1E' % volnorm) + ":" + \
                "\nCO: " + str('%.2f' % (mixinglist[CO,0] * 100)) + " %" + \
                "\nO2: " + str('%.2f' % (mixinglist[O2,0] * 100)) + " %" + \
                "\nTotal: " + `sum(mixinglist[:,0]) + 0.95 + 0.016` + \
                "\n\nMixing Ratios at " + `finalvol` + ":" + \
                "\nCO: " + str('%.2f' % (mixinglist[CO,finalconv] * 100)) + " %" + \
                "\nO2: " + str('%.2f' % (mixinglist[O2,finalconv] * 100)) + " %" + \
                "\nTotal: " + `sum(mixinglist[:,finalconv]) + 0.95 + 0.016` + \
                "\n\nTotal Run Time:   " + elapsed_time_str + \
                "\n\nSO2 Flux    Converged Line\n"
    text_readme.write(printtext)
    text_readme.close()
    
    convwrite(volnorm, increments, iterations, runnum)
    
    return volcswitch  #returns the redox-switch value

def Gbuffer(FMQ, water, pressure):
    """Use this to define which set of Gaillard et al. 2013 values you want to use depending
    on the type of buffer."""
    
    print FMQ
    print water
    print pressure
    
    if FMQ == 0.5 and water == 0.4 and pressure == 0.01:
        h2_i =  87
        h2s_i = 85
        so2_i = 2667
        s2_i =  1404
        co_i =  318
        
    elif FMQ == 0.5 and water == 0.4 and pressure == 0.1:
        h2_i =  61
        h2s_i = 171
        so2_i = 1981
        s2_i =  1137
        co_i =  243
    
    elif FMQ == 0.5 and water == 0.4 and pressure == 1:
        h2_i =  36
        h2s_i = 264
        so2_i = 1357
        s2_i =  704
        co_i =  170
        
    elif FMQ == 0.5 and water == 0.2 and pressure == 0.01:
        h2_i =  46
        h2s_i = 60
        so2_i = 1766
        s2_i =  1542
        co_i =  338
        
    elif FMQ == 0.5 and water == 0.2 and pressure == 0.1:
        h2_i =  30
        h2s_i = 108
        so2_i = 1387
        s2_i =  1095
        co_i =  251
        
    elif FMQ == 0.5 and water == 0.2 and pressure == 1:
        h2_i =  15
        h2s_i = 140
        so2_i = 975
        s2_i =  599
        co_i =  168
        
    elif FMQ == 0.5 and water == 0.1 and pressure == 0.01:
        h2_i =  23
        h2s_i = 36
        so2_i = 1319
        s2_i =  1521
        co_i =  345
        
    elif FMQ == 0.5 and water == 0.1 and pressure == 0.1:
        h2_i =  14
        h2s_i = 60
        so2_i = 1095
        s2_i =  1004
        co_i =  250
        
    elif FMQ == 0.5 and water == 0.1 and pressure == 1:
        h2_i =  6
        h2s_i = 63
        so2_i = 783
        s2_i =  503
        co_i =  160
        
    elif FMQ == 1.4 and water == 0.4 and pressure == 0.01:
        h2_i =  96
        h2s_i = 102
        so2_i = 1910
        s2_i =  1562
        co_i =  320
        
    elif FMQ == 1.4 and water == 0.4 and pressure == 0.1:
        h2_i =  69
        h2s_i = 200
        so2_i = 1324
        s2_i =  1140
        co_i =  251
        
    elif FMQ == 1.4 and water == 0.4 and pressure == 1:
        h2_i =  42
        h2s_i = 298
        so2_i = 781
        s2_i =  601
        co_i =  184
        
    elif FMQ == 1.4 and water == 0.2 and pressure == 0.01:
        h2_i =  52
        h2s_i = 72
        so2_i = 1053
        s2_i =  1544
        co_i =  309
        
    elif FMQ == 1.4 and water == 0.2 and pressure == 0.1:
        h2_i =  35
        h2s_i = 126
        so2_i = 771
        s2_i =  966
        co_i =  237
        
    elif FMQ == 1.4 and water == 0.2 and pressure == 1:
        h2_i =  18
        h2s_i = 154
        so2_i = 431
        s2_i =  417
        co_i =  169
        
    elif FMQ == 1.4 and water == 0.1 and pressure == 0.01:
        h2_i =  26
        h2s_i = 44
        so2_i = 625
        s2_i =  1338
        co_i =  252
        
    elif FMQ == 1.4 and water == 0.1 and pressure == 0.1:
        h2_i =  16
        h2s_i = 69
        so2_i = 480
        s2_i =  744
        co_i =  188
        
    elif FMQ == 1.4 and water == 0.1 and pressure == 1:
        h2_i =  7
        h2s_i = 63
        so2_i = 233
        s2_i =  251
        co_i =  132
        
    elif FMQ == 1.4 and water == 0.01 and pressure == 0.01:
        h2_i =  2
        h2s_i = 4
        so2_i = 404
        s2_i =  1078
        co_i =  365
        
    elif FMQ == 1.4 and water == 0.01 and pressure == 0.1:
        h2_i =  1
        h2s_i = 5
        so2_i = 332
        s2_i =  548
        co_i =  262
        
    elif FMQ == 1.4 and water == 0.01 and pressure == 1:
        h2_i =  0
        h2s_i = 2
        so2_i = 166
        s2_i =  176
        co_i =  181
        
    elif FMQ == 3.5 and water == 0.4 and pressure == 0.01:
        h2_i =  115
        h2s_i = 131
        so2_i = 1027
        s2_i =  1595
        co_i =  140
        
    elif FMQ == 3.5 and water == 0.4 and pressure == 0.1:
        h2_i =  86
        h2s_i = 244
        so2_i = 616
        s2_i =  957
        co_i =  114
        
    elif FMQ == 3.5 and water == 0.4 and pressure == 1:
        h2_i =  58
        h2s_i = 334
        so2_i = 238
        s2_i =  339
        co_i =  93
        
    elif FMQ == 3.5 and water == 0.1 and pressure == 0.01:
        h2_i =  33
        h2s_i = 50
        so2_i = 156
        s2_i =  770
        co_i =  164
        
    elif FMQ == 3.5 and water == 0.1 and pressure == 0.1:
        h2_i =  23
        h2s_i = 70
        so2_i = 74
        s2_i =  262
        co_i =  138
        
    elif FMQ == 3.5 and water == 0.1 and pressure == 1:
        h2_i =  14
        h2s_i = 49
        so2_i = 6
        s2_i =  23
        co_i =  134
        
    elif FMQ == 3.5 and water == 0.01 and pressure == 0.01:
        h2_i =  3
        h2s_i = 4
        so2_i = 20
        s2_i =  165
        co_i =  177
        
    elif FMQ == 3.5 and water == 0.01 and pressure == 0.1:
        h2_i =  1
        h2s_i = 3
        so2_i = 3
        s2_i =  21
        co_i =  167
        
    elif FMQ == 3.5 and water == 0.01 and pressure == 1:
        h2_i =  0
        h2s_i = 1
        so2_i = 0
        s2_i =  2
        co_i =  168
    
    returnset = []
    returnset.append(so2_i)
    returnset.append(h2s_i)
    returnset.append(h2_i)
    returnset.append(s2_i)
    returnset.append(co_i)
    
    return returnset

def run():
    """Runs the model, asks for buffer parameters and run number"""
    
    othernotes = raw_input("Notes:  ")
    runnum = int(raw_input("Run Number:  "))
    FMQ = float(raw_input("FMQ-0.5 (Enter '0.5')\nFMQ-1.4 (Enter '1.4')\nIW (Enter '3.5')\n"))
    water = float(raw_input("Water Content Percent (e.g. '0.4' for 0.4%)"))
    pressure = float(raw_input("Pressure in bar (0.01, 0.1, or 1):"))
    
    notes =  "FMQ-" + `FMQ` + ", " + `water` + "% H2O, " + `pressure` + " bar" + othernotes
    
    so2_i = Gbuffer(FMQ, water, pressure)[0]
    h2s_i = Gbuffer(FMQ, water, pressure)[1]
    h2_i = Gbuffer(FMQ, water, pressure)[2]
    s2_i = Gbuffer(FMQ, water, pressure)[3]
    co_i = Gbuffer(FMQ, water, pressure)[4]
    multirun(so2_i, h2_i, h2s_i, s2_i, co_i, runnum, notes)

def cleanup():
    """Removes all the temporary output files to cleanup the main directory"""
    subprocess.call("rm out.trs ", shell=True)
    subprocess.call("rm out.time ", shell=True)
    subprocess.call("rm out.tim ", shell=True)
    subprocess.call("rm out.terse ", shell=True)
    subprocess.call("rm out.tau ", shell=True)
    subprocess.call("rm out.so2 ", shell=True)
    subprocess.call("rm out.redox ", shell=True)
    subprocess.call("rm out.rates ", shell=True)
    subprocess.call("rm out.raingc ", shell=True)
    subprocess.call("rm out.rad ", shell=True)
    subprocess.call("rm out.prod ", shell=True)
    subprocess.call("rm out.params ", shell=True)
    subprocess.call("rm out.NOprates ", shell=True)
    subprocess.call("rm out.gridz ", shell=True)
    subprocess.call("rm out.gridw ", shell=True)
    subprocess.call("rm out.flux ", shell=True)
    subprocess.call("rm out.flow ", shell=True)
    subprocess.call("rm out.finalden ", shell=True)
    subprocess.call("rm out.error ", shell=True)
    subprocess.call("rm out.densities ", shell=True)
    subprocess.call("rm ISOinert.dist ", shell=True)
    subprocess.call("rm ISOin.dist ", shell=True)
    subprocess.call("rm out.cl ", shell=True)
    subprocess.call("rm out.converge ", shell=True)
    subprocess.call("rm out.xsec ", shell=True)
    subprocess.call("rm out.oxygencolumnden ", shell=True)
    subprocess.call("rm out.perchloratemr ", shell=True)
    subprocess.call("rm out.pna ", shell=True)
    subprocess.call("rm out.temprun ", shell=True)
    subprocess.call("rm out.intrates ", shell=True)
    subprocess.call("rm fort.13 ", shell=True)
    subprocess.call("rm fort.22 ", shell=True)
    subprocess.call("rm fort.20 ", shell=True)
    subprocess.call("rm inertmixingratios.out ", shell=True)

#To run properly want run() and cleanup()
run()
cleanup()
#multirun(0,0,0,0,0,64,"Testing with no volcanic fluxes")
