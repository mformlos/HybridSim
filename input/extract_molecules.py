
molecule_file = open("/home/formanek/HYBRIDSIM/input/filebondunlinked","r") 
bond_file = open("/home/formanek/HYBRIDSIM/input/bondlinked0","r")
config_file = open("/home/formanek/HYBRIDSIM/input/config0", "r")

current_molecule = 0; 
mono_count = 0

molecule = []
bond1 = []
bond2 = []
molecule_sizes = []

##CHAIN TOPOLOGY

for i, line in enumerate(molecule_file):
    if i == 0: 
        continue
    data = line.split()
    molecule.append(int(data[0]))
    bond1.append(int(data[1]))
    bond2.append(int(data[2]))    
    
last_molecule = molecule[-1]

for i in range(1,last_molecule+1): 
    molecule_sizes.append(molecule.count(i))

    
for i in range(len(molecule_sizes)): 
    filename = "SCNPs/SCNP-"+str(i)+"-chain" 
    with open(filename, "w") as outF: 
        outF.write(" %10d \n" %molecule_sizes[i])
        for j in range(mono_count,mono_count+molecule_sizes[i]): 
            outF.write(" %10d %10d %10d \n" %(1, bond1[j], bond2[j]))
    mono_count += molecule_sizes[i]

#
## SCNP TOPOLOGY
#


molecule =[]
bond1 = []
bond2 = []
molecule_bonds = []
mono_count = 0

for i, line in enumerate(bond_file): 
    if i == 0: 
        continue
    data = line.split()
    molecule.append(int(data[0]))
    bond1.append(int(data[1]))
    bond2.append(int(data[2]))

last_molecule = molecule[-1]

for i in range(1,last_molecule+1): 
    molecule_bonds.append(molecule.count(i))
    
for i in range(len(molecule_bonds)): 
    filename = "SCNPs/SCNP-"+str(i)+"-bonds"
    with open(filename, "w") as outF: 
        outF.write(" %10d \n" %molecule_bonds[i])
        for j in range(mono_count,mono_count+molecule_bonds[i]): 
            outF.write(" %10d %10d %10d \n" %(1, bond1[j], bond2[j])) 
    mono_count += molecule_bonds[i]    
    
    
## CONFIG 

current_mol = 0
count = 0
outF = open("SCNPs/SCNP-"+str(current_mol)+"-config", "w")
for line in config_file: 
    print current_mol
    if (count > molecule_sizes[current_mol]):
        outF.close()
        current_mol += 1
        outF = open("SCNPs/SCNP-"+str(current_mol)+"-config", "w")
        count = 0
    count += 1
    outF.write(line) 
    
outF.close()
    
        


