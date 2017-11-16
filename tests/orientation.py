from ElectrostaticTools import *

mol = GeneralMultipoledMolecule()
mol.readMe("h2co.mmol")

mol_bak = mol

print ("Initial multipoles:")
print (mol.GetMultipoles(2).toString())

x = Vector3(1,0,0)
y = Vector3(0,1,0)
z = Vector3(0,0,1)

mat = Matrix3(x,y,z)
mol = mol_bak
mol.Rotate(mat)

print ("Rotation matrix: ")
print (mat.toString())

print ("Rotated multipoles: ")
print (mol.GetMultipoles(2).toString())

mat = Matrix3(y,x,z)
mol = mol_bak
mol.Rotate(mat)

print ("Rotation matrix: ")
print (mat.toString())

print ("Rotated multipoles: ")
print (mol.GetMultipoles(2).toString())

mat = Matrix3(z,y,x)
mol = mol_bak
mol.Rotate(mat)

print ("Rotation matrix: ")
print (mat.toString())

print ("Rotated multipoles: ")
print (mol.GetMultipoles(2).toString())
