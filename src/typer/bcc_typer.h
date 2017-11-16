#ifndef BCC_TYPER_HEADER
#define BCC_TYPER_HEADER

#include <string>
#include <vector>

#include <openbabel/mol.h>

#include "typer.h"

class BCCAtomTyper : public AtomTyper
{
public:
  BCCAtomTyper(const std::vector<std::string>& types) : t(types) {}
  std::string type(OpenBabel::OBAtom* a) const { return t[a->GetIdx() - 1]; } // a stub 
private:
  const std::vector<std::string>& t;
};

class BCCBondTyper : public BondTyper
{
public:
  BCCBondTyper(const BCCAtomTyper& t) : atyper(t) {}
  std::string type(OpenBabel::OBBond* bond) const
  {
    OpenBabel::OBAtom* pa = bond->GetBeginAtom();
    OpenBabel::OBAtom* pb = bond->GetEndAtom();
    std::string bo, a = atyper.type(pa), b = atyper.type(pb);
    if (bond->IsSingle())
    {
      if (bond->IsAromatic()) bo = "07"; else
      {
        int ach = pa->GetFormalCharge(), bch = pb->GetFormalCharge();
        if ( (ach > 0 && bch < 0) || ( ach < 0 && bch > 0) ) bo = "06"; else
        if ( (ach == 0 && bch < 0) || ( ach < 0 && bch == 0) ) bo = "09"; else
        bo = "01";
      }
    } else if (bond->IsDouble())
    {
      if (bond->IsAromatic()) bo = "08"; else
      bo = "09";
    } else if (bond->IsTriple())
    {
      bo = "03";
    } else if (bond->IsAromatic())
    {
      bo = "09";
    } else throw std::domain_error("Unknown bond detected");
    return ( a > b ) ? b + bo + a : a + bo + b;
  }
private:
  const BCCAtomTyper& atyper;
};

#endif // TYPERHEADER

