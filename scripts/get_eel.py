#!/usr/bin/python3

import pybel
from openbabel import *
from ElectrostaticTools import *
import argparse
import itertools

parser = argparse.ArgumentParser(description='Get intra- and intermolecular electrostatic energies'\
                                             ' for a set of molecules in mmol format.')
parser.add_argument('-I', '--input', type=argparse.FileType('r'), nargs='+', required=True, help="Molecule file(s)")
parser.add_argument('--weight-1-2', type=float, default=0.0)
parser.add_argument('--weight-1-3', type=float, default=0.0)
parser.add_argument('--weight-1-4', type=float, default=1/1.2)
parser.add_argument('--weight-intra', type=float, default=1.0)
parser.add_argument('--weight-inter', type=float, default=1.0)


args = parser.parse_args()

mols = []
energies = []
labels = []

for f in args.input:
  mol = GeneralMultipoledMolecule()
  mol.readMe(f.name)
  mols.append(mol)

for f in args.input:
  labels.append("intra_{}".format(f.name))
for f in args.input:
  labels.append("inter_{}".format(f.name))

# intramolecular
 
weights = [args.weight_intra]*100
weights[0:3] = [0.0, args.weight_1_2*args.weight_intra, args.weight_1_3*args.weight_intra, args.weight_1_4*args.weight_intra]

for mol in mols:
  energy = 0.0
  for i in range(1, mol.NumAtoms()+1):
    for atom, depth in OBMolAtomBFSIter(mol, i):
      e  = 0.0
      dist = depth - 1
      if dist != 0:
        e = mol.GetMultipoles(i).energy( mol.GetMultipoles(atom.GetIdx()) ) * weights[dist] * 0.5
      energy = energy + e
  energies.append(energy*627.503)

# intermolecular  

i = 0
for molA in mols:
  selector = [1]*len(mols)
  selector[i] = 0
  energy = 0.0
  for molB in itertools.compress(mols, selector):
    energy = energy + molA.E_El(molB)*627.503
  i = i + 1
  energies.append(energy)

print (" ".join(labels))
print (" ".join(format(x, "10.7f") for x in energies))

