import sys
import numpy as np 

if len(sys.argv) < 2:
    print "usage: python force_extension_average.py dirname begintime maxtime " 
    exit()
    
    
dirname = str(sys.argv[1])
filename = dirname+"extension"

begintime = float(sys.argv[2])
maxtime = float(sys.argv[3])

outputf = open(dirname+"extension_average", "w")
inputf = open(filename, "r")

extension = float(inputf.readline().split()[1])
print extension
 

intervaltime = 0.0
force = 0.0

n = 0

for i, line in enumerate(inputf):
    line = line.split()
    extension_now = float(line[1])
    time_now = float(line[0])
    if (extension_now != extension):
        force = force/n 
        outputf.write(str(extension)+" "+str(force)+" "+str(n)+"\n")
        force = 0.0
        n = 0
        extension = extension_now
        intervaltime = time_now    
    if (time_now > begintime+intervaltime):
        force += float(line[2]) 
        n += 1
