from openbabel import *
from ElectrostaticTools import *

print ("This is a simple test to check whether multipole energies and potentials are calculated correctly")
print ("If you see lots of 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 6.0, -6.0 and the values does not change in a single row - then its ok.")

plus = Vector3(bohr/2,0,0)
minus = Vector3(-bohr/2,0,0)

mult1 = Multipoles(plus)
mult2 = Multipoles(minus)

mult1.setMonopole(Monopole(1.0))
mult2.setMonopole(Monopole(1.0))
print ("m-m,m-m 1,1 : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult2.setMonopole(Monopole(-1.0))
print ("m-m,m-m 1,-1 : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult2.deleteMonopole()
mult2.setDipole(Dipole(1,0,0))
print ("m-d,d-m 1,(1,0,0) : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult2.setDipole(Dipole(-1,0,0))
print ("m-d,d-m 1,(-1,0,0) : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult2.deleteDipole()
mult2.setQuadrupole(Quadrupole(1, -0.5, -0.5))
print ("m-q,q-m 1,(1.0,-0.5,-0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult2.setQuadrupole(Quadrupole(-1, 0.5, 0.5))
print ("m-q,q-m 1,(-1.0,0.5,0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1), " ", mult2.potential(plus))

mult1.deleteMonopole()
mult1.setDipole(Dipole(1.0,0.0,0.0))
mult2.deleteQuadrupole()
mult2.setDipole(Dipole(1.0,0.0,0.0))
print ("d-d,d-d (1,0,0),(1,0,0) : ", mult1.energy(mult2), " ", mult2.energy(mult1))

mult2.setDipole(Dipole(-1.0,0.0,0.0))
print ("d-d,d-d (1,0,0),(-1,0,0) : ", mult1.energy(mult2), " ", mult2.energy(mult1))

mult2.deleteDipole()
mult2.setQuadrupole(Quadrupole(1, -0.5, -0.5))
print ("d-q,q-d (1,0,0),(1.0,-0.5,-0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1))

mult2.setQuadrupole(Quadrupole(-1, 0.5, 0.5))
print ("d-q,q-d (1,0,0),(-1.0,0.5,0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1))

mult1.deleteDipole()
mult1.setQuadrupole(Quadrupole(1, -0.5, -0.5))
print ("q-q,q-q (1,-0.5,-0.5),(-1.0,0.5,0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1))

mult2.setQuadrupole(Quadrupole(1, -0.5, -0.5))
print ("q-q,q-q (1,-0.5,-0.5),(1.0,-0.5,-0.5) : ", mult1.energy(mult2), " ", mult2.energy(mult1))
