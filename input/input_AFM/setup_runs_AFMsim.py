import sys
import glob 
import os 
import random
import shutil
import subprocess 
from collections import namedtuple 


random.seed()
RunsPerSCNP = 5
Parameter = namedtuple("Parameters", "Box, Temperature, Gamma, MDStep, SimTime, SurfaceEnergyStart, SurfaceEnergyEquil, AnchorTime, SCNPNumber")

Box = namedtuple("Box", ["Lx", "Ly", "Lz"])

ParameterSets = []

ParameterSets.append(Parameter(Box(50, 50, 200), 1.0, 0.05, 0.01, 1000000.0, 2.0, 1.5, 500., [1, 17]))

submit_files =open("runs_to_submit.dat", "w")
if not os.path.exists("/scratch-new/formanek/AFMSIM/"): 
    os.makedirs("/scratch-new/formanek/AFMSIM/")
if not os.path.exists("/scratch-new/formanek/AFMSIM/runs/"):
    os.makedirs("/scratch-new/formanek/AFMSIM/runs/")
    
for paramSet in ParameterSets: 
    Lx = paramSet.Box.Lx
    Ly = paramSet.Box.Ly
    Lz = paramSet.Box.Lz
    Temperature = paramSet.Temperature
    Gamma = paramSet.Gamma
    MDStep = paramSet.MDStep
    SimTime = paramSet.SimTime
    AnchorTime = paramSet.AnchorTime
    SurfaceEnergyStart = paramSet.SurfaceEnergyStart
    SurfaceEnergyEquil = paramSet.SurfaceEnergyEquil    
    for SCNP in paramSet.SCNPNumber: 
        for i in range(RunsPerSCNP): 
            Seed = random.randint(0,999999)
            run_name = "AFM-SCNP-"+str(SCNP)+"-REPL-"+str(i)
            directory = "/scratch-new/formanek/AFMSIM/runs/"+run_name
            
            if os.path.exists(directory):
                shutil.rmtree(directory)
            os.makedirs(directory) 
            os.makedirs(directory+"/configs")
            shutil.copyfile("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-"+str(SCNP)+"-config", directory+"/config")
            shutil.copyfile("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-"+str(SCNP)+"-bonds", directory+"/bonds")    
            shutil.copyfile("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-"+str(SCNP)+"-chain", directory+"/chain")
            shutil.copyfile("filetimeoutput", directory+"/steps")
            shutil.copyfile("forceupdatefile", directory+"/force")
            submit_files.write(directory+"\n") 
            with open("parameter_template_afm.dat", "r") as inF:
                lines = inF.readlines()
            with open(directory+"/parameters.dat", "w") as outF:
                outF.write("BoxX = "+str(Lx)+"\n")
                outF.write("BoxY = "+str(Ly)+"\n")
                outF.write("BoxZ = "+str(Lz)+"\n\n")   
                outF.write("Temperature = "+str(Temperature)+"\n")
                outF.write("Gamma = "+str(Gamma)+"\n\n") 
                outF.write("MDStep = "+str(MDStep)+"\n\n") 
                outF.write("StartTime = "+str(0.0)+"\n")
                outF.write("SimTime = "+str(SimTime)+"\n\n")
                outF.write("SurfaceEnergyStart = "+str(SurfaceEnergyStart)+"\n")
                outF.write("SurfaceEnergyEquil = "+str(SurfaceEnergyEquil)+"\n\n")
                outF.write("AnchorTime = "+str(AnchorTime)+"\n\n")
                outF.write("Seed = "+str(AnchorTime)+"\n\n")
                
                lineInd = 0
                for currLine in lines: 
                    if lineInd >= 19 : 
                        outF.write(currLine)
                    lineInd += 1

