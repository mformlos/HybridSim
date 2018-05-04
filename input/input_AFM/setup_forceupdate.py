import sys
import glob 
import os 
import random
import shutil

ForceUpdateInterval = 10000000
ForceIncrement = 0.01
ForceMax = 5.0
Force = 0.0

StepInterval = 0
Stepcurrent = 0 

Stepin = open("filetimetemplate", "r")
Stepout = open("filetimeoutput", "w")
Forceout = open("forceupdatefile", "w")

while (Force <= ForceMax):
    Forceout.write(str(StepInterval)+" 0.0 0.0 "+str(Force)+" \n") 
    while (Stepcurrent < ForceUpdateInterval): 
        Stepcurrent = int(Stepin.readline())
        Stepout.write(str(StepInterval+Stepcurrent)+ " \n")
    Stepin.seek(0)    
    Stepcurrent = 0
    StepInterval += ForceUpdateInterval
    Force += ForceIncrement
        
        
    

