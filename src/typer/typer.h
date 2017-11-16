#ifndef TYPER_HEADER
#define TYPER_HEADER

#include <string>

#include <openbabel/mol.h>

class AtomTyper
{
public:
  virtual std::string type(OpenBabel::OBAtom* a) const = 0;
  virtual ~AtomTyper() {}
};

class BondTyper
{
public:
  virtual std::string type(OpenBabel::OBBond* b) const = 0;
  virtual ~BondTyper() {}
};


#endif // TYPERHEADER

