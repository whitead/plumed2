#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "core/ActionAtomistic.h"
#include "tools/Pbc.h"
#include "tools/LinkCells.h"
#include "tools/AtomNumber.h"

//for MPI communicator
#include "core/PlumedMain.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace Virial{

//+PLUMEDOC MCOLVAR DISTANCES
/*
Calculate ?

\par Examples


\verbatim
?
\endverbatim
(See also \ref ?).

\par TODO

1. histogram RDF and validate
3. Get rij outter product loop
2. Our function is derivative of log of RDF, which I will add to second grid at each step (averaging it).
4. Get rj via unwrapping
5. Compute sum


*/
//+ENDPLUMEDOC


class Virial : public Colvar {
private:
  //if we are computing self-virial
  bool b_self_virial;
  //Neighbor list linkscells to accelerate
  LinkCells linkcells;
  //atoms
  vector<AtomNumber> group_1;
  vector<AtomNumber> group_2;
  //used to store neighboring atoms
  vector<unsigned int> neighs;

  void setupLinkCells();  

public:
  static void registerKeywords( Keywords& keys );
  explicit Virial(const ActionOptions&);
// active methods:
  virtual void calculate();
  virtual void prepare();
};

PLUMED_REGISTER_ACTION(Virial,"VIRIAL")

void Virial::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);

  keys.add("compulsory", "CUTOFF", "The maximum distance to consider for pairwise calculations");

  //  keys.reset_style("ATOMS","atoms");
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
                              "the atoms in GROUPB. This must be used in conjuction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
                              "in GROUPB. This must be used in conjuction with GROUPA.");
}

Virial::Virial(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
b_self_virial(false),
neighs(0),
linkcells(comm)
{
  //used for link list
  double cutoff;
  // Read in the atoms
  vector< vector<AtomNumber> > all_atoms(3);

  parse("CUTOFF", cutoff);
  parseAtomList("GROUP", all_atoms[0]);
  parseAtomList("GROUPA", all_atoms[1]);
  parseAtomList("GROUPB", all_atoms[2]);
  
  if(all_atoms[0].size() < 2) {
    if(all_atoms[1].size() < 2 && all_atoms[2].size() < 2) {
      error("Must specify both GROUPA and GROUPB or GROUP");
    }
    b_self_virial = false;    
    group_1 = all_atoms[1];
    group_2 = all_atoms[2];
  } else {
    b_self_virial = true;    
    group_1 = all_atoms[0];
    group_2 = all_atoms[0];
  }

  //make neighbors big enough for all neighbors case
  neighs.resize(group_1.size() + group_2.size());
			      
  // And check everything has been read in correctly
  checkRead();


  linkcells.setCutoff( cutoff );

  addValueWithDerivatives();
  setNotPeriodic();

}



  void Virial::prepare() {
    requestAtoms(group_1);
    requestAtoms(group_2);
  }
  
void Virial::setupLinkCells(){
  std::vector<Vector> ltmp_pos( group_1.size() + group_2.size() * b_self_virial ); 
  std::vector<unsigned> ltmp_ind( group_1.size() + group_2.size() * b_self_virial ); 


  for(unsigned int i = 0; i < group_1.size(); i++) {
    ltmp_pos[i] = getPosition(group_1[i]);
    ltmp_ind[i] = i;
  }  
  unsigned int j;
  for(unsigned int i = group_1.size(); i < ltmp_pos.size(); i++) {
    j = i - group_1.size();
    ltmp_pos[i] = getPosition(group_2[j]);
    ltmp_ind[i] = j;
  }

  // Build the lists for the link cells
  linkcells.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );
}



void Virial::calculate()  {
  setupLinkCells();

  //number of neighbors
  unsigned int nn;
  Vector rij;
  double r;
  for(unsigned int i = 0; i < group_1.size(); ++i) {

    nn = 1;
    neighs[0] = i;
    linkcells.retrieveNeighboringAtoms( getPosition(group_1[i]), nn, neighs);
    for(unsigned int j = 0; j < nn; ++j) {
      if(!b_self_virial && neighs[j] < group_1.size())
	continue;
      rij = delta(getPosition(group_1[i]), getPosition(group_2[neighs[j]]));
      r = rij.modulo();
      //do thing tiwht r
      log.printf("pairwise = %f\n", r);
    }
  }

  //get pairs 
  setValue(0);

}


}
}

