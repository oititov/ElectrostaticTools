#!/usr/bin/python3

from openbabel import *
from ElectrostaticTools import *
import argparse

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Get multipole parameters from molecule')
parser.add_argument('-I', '--input', type=argparse.FileType('r'), nargs='+', \
                      required=True, help="Input molecule file")
parser.add_argument('-M', '--mask', required=True, help="Atom SMARTS mask")
parser.add_argument("-m", '--monopole', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get charge information.")
parser.add_argument("-d", '--dipole', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get dipole information.")
parser.add_argument('--dx', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get Dx information.")
parser.add_argument('--dy', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get Dy information.")
parser.add_argument('--dz', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get Dz information.")
parser.add_argument("-q", '--quadrupole', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get quadrupole information.")
parser.add_argument('--qxx', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get quadrupole information.")
parser.add_argument('--qyy', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get quadrupole information.")
parser.add_argument('--qzz', type=str2bool, nargs='?', const=True, default=False, \
                      help="Get quadrupole information.")
parser.add_argument('-b', '--fill-blanks', nargs='?', const='0', default=None, \
                      help="Use a specified spacer, if no multipole found.")
parser.add_argument('--no-header', type=str2bool, nargs='?', const=True, default=False, \
                      help="Hide a header.")

args = parser.parse_args()


header="filename"
if not args.no_header :
  if args.monopole :
    header = header + " m"
  if args.dipole :
    header = header + " ( dx, dy, dz )"
  if args.dx :
    header = header + " dx"
  if args.dy :
    header = header + " dy"
  if args.dz :
    header = header + " dz"
  if args.quadrupole :
    header = header + " ( Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz )"
  if args.qxx :
    header = header + " Qxx"
  if args.qyy :
    header = header + " Qyy"
  if args.qzz :
    header = header + " Qzz"
  print(header)  

for inmol in args.input:
  mol = GeneralMultipoledMolecule()
  mol.readMe(inmol.name)
  result=inmol.name + " "

  for atom in OBMolAtomIter(mol):
    if atom.MatchesSMARTS(args.mask) :
      mult = mol.GetRawMultipoles(atom.GetIdx())
      if args.monopole :
        result = result + "{:.5f} ".format(mult.monopole().value())
      dipole = Vector3()
      if mult.hasDipole() :
        dipole = mult.dipole().value()
        if args.dipole :
          result = result + "{} ".format(dipole.toString())
        if args.dx :
          result = result + "{:.5f} ".format(dipole.x())
        if args.dy :
           result = result + "{:.5f} ".format(dipole.y())
        if args.dz :
          result = result + "{:.5f} ".format(dipole.z())
      elif args.fill_blanks is not None:
        if args.dipole :
          result = result + "( {}, {}, {} ) ".format(args.fill_blanks, args.fill_blanks, args.fill_blanks)
        if args.dx :
          result = result + "{} ".format(args.fill_blanks)
        if args.dy :
          result = result + "{} ".format(args.fill_blanks)
        if args.dz :
          result = result + "{} ".format(args.fill_blanks)
      if mult.hasQuadrupole() :
        quad = mult.quadrupole().value()
        if args.quadrupole :
          result = result + "{} ".format(quad.toString())
        if args.qxx :
          result = result + "{:.5f} ".format(quad(0,0))
        if args.qyy :
          result = result + "{:.5f} ".format(quad(1,1))
        if args.qzz :
          result = result + "{:.5f} ".format(quad(2,2))
      elif args.fill_blanks is not None:
        if args.quadrupole :
          result = result + "( {}, {}, {}, {}, {}, {}, {}, {}, {} ) ".format( \
                                             args.fill_blanks, args.fill_blanks, args.fill_blanks, \
                                             args.fill_blanks, args.fill_blanks, args.fill_blanks, \
                                             args.fill_blanks, args.fill_blanks, args.fill_blanks \
                                                                          )
        if args.qxx :
          result = result + "{} ".format(args.fill_blanks)
        if args.qyy :
          result = result + "{} ".format(args.fill_blanks)
        if args.qzz :
          result = result + "{} ".format(args.fill_blanks)
      print(result)
      
