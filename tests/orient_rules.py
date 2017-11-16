import sys
from openbabel import *
from ElectrostaticTools import *

filename = str(sys.argv[1])

print ("Using file : ", filename)

mol = GeneralMultipoledMolecule()
mol.readMe(filename)

mol.printMe()

rules = mol.getOrientMatches()

k = 0

for i in rules:
  print(str(k), " - ", i)
  k+=1 
