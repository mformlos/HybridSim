import sys
import glob 
import os 
import random
import shutil
import subprocess 
from collections import namedtuple 

################################
### THIS FILE IS DEPRECATED ####
################################

Parameter = namedtuple("Parameters", "Box, Shear, Temperature, MPCRho, MPCStep, MDStep, SimTime, EquilTime, SCNPNumber")

Box = namedtuple("Box", ["Lx", "Ly", "Lz"])

ParameterSets = []

ParameterSets.append(Parameter(Box(50, 50, 50), [0.01], 1.0, 10, 0.1, 0.001, 200000.0, 0.0,[18, 27]))
#ParameterSets.append(Parameter(Box(50, 50, 50), [0.0, 0.00001, 0.00005, 0.0001], 1.0, 10, 0.1, 0.01, 8010000.0, 0.0,[i for i in range(50)]))
#ParameterSets.append(Parameter(Box(50, 50, 50), [0.0], 1.0, 5, 0.1, 0.01, 1000.0, 0.0,[i for i in range(5)]))

submit_files =open("runs_to_submit.dat", "w")
if not os.path.exists("/scratch-new/formanek/HYBRIDSIM/"): 
    os.makedirs("/scratch-new/formanek/HYBRIDSIM/")
if not os.path.exists("/scratch-new/formanek/HYBRIDSIM/runs/"):
    os.makedirs("/scratch-new/formanek/HYBRIDSIM/runs/")

for paramSet in ParameterSets: 
    Lx = paramSet.Box.Lx
    Ly = paramSet.Box.Ly
    Lz = paramSet.Box.Lz
    Temperature = paramSet.Temperature
    Rho = paramSet.MPCRho
    MPCStep = paramSet.MPCStep
    MDStep = paramSet.MDStep
    SimTime = paramSet.SimTime
    EquilTime = paramSet.EquilTime
    for Shear in paramSet.Shear: 
        for SCNP in paramSet.SCNPNumber: 
            run_name = "SCNP-"+str(SCNP)+"-Shear-"+str(Shear)
            directory = "/scratch-new/formanek/HYBRIDSIM/runs/"+run_name
            
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.makedirs(directory) 
            os.makedirs(directory+"/configs")
            shutil.copyfile("SCNPs/SCNP-"+str(SCNP)+"-config", directory+"/config")
            shutil.copyfile("SCNPs/SCNP-"+str(SCNP)+"-bonds", directory+"/bonds")    
            shutil.copyfile("SCNPs/SCNP-"+str(SCNP)+"-chain", directory+"/chain")
            shutil.copyfile("filetimeoutput", directory+"/steps")
            submit_files.write(directory+"\n") 
            with open("parameter_template.dat", "r") as inF:
                lines = inF.readlines()
            with open(directory+"/parameters.dat", "w") as outF:
                outF.write("BoxX = "+str(Lx)+"\n")
                outF.write("BoxY = "+str(Ly)+"\n")
                outF.write("BoxZ = "+str(Lz)+"\n\n")
                outF.write("Shear = "+str(Shear)+"\n")    
                outF.write("Temperature = "+str(Temperature)+"\n\n")
                outF.write("MPCRho = "+str(Rho)+"\n\n")
                outF.write("MPCStep = "+str(MPCStep)+"\n")  
                outF.write("MDStep = "+str(MDStep)+"\n\n") 
                outF.write("StartTime = "+str(0.0)+"\n")
                outF.write("SimTime = "+str(SimTime)+"\n")
                outF.write("EquilTime = "+str(EquilTime)+"\n\n")
                lineInd = 0
                for currLine in lines: 
                    if lineInd >= 16 : 
                        outF.write(currLine)
                    lineInd += 1
        
                
                
                
            

