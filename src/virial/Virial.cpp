#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "core/ActionAtomistic.h"
#include "tools/Pbc.h"
#include "tools/LinkCells.h"
#include "tools/AtomNumber.h"
#include "tools/Grid.h"

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
  double rdf_bw_;
  double cutoff_;

  vector<unsigned int> rdf_kernel_support_;
  //atoms
  vector<AtomNumber> group_1;
  vector<AtomNumber> group_2;
  //used to store neighboring atoms
  vector<unsigned int> neighs;
  //Neighbor list linkscells to accelerate
  LinkCells linkcells;
  OFile rdffile_;
  OFile mf_file_;


  Grid* gr_grid_;
  Grid* dgr_grid_;


  void setup_link_cells_();
  double eval_rdf_(const double x, const double r, const double scale, vector<double>& der, vector<double>& der2);

public:
  static void registerKeywords( Keywords& keys );
  explicit Virial(const ActionOptions&);
  ~Virial();
// active methods:
  virtual void calculate();
  virtual void prepare();
};

PLUMED_REGISTER_ACTION(Virial,"VIRIAL")

void Virial::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);

  keys.add("compulsory", "CUTOFF", "The maximum distance to consider for pairwise calculations");
  keys.add("compulsory", "RDF_MIN", "The maximum distance to consider for pairwise calculations");
  keys.add("compulsory", "RDF_MAX", "The maximum distance to consider for pairwise calculations");
  keys.add("compulsory", "RDF_NBINS", "The maximum distance to consider for pairwise calculations");
  keys.add("compulsory", "RDF_BW", "The maximum distance to consider for pairwise calculations");
  keys.addFlag("RDF_SPLINE", true, "The maximum distance to consider for pairwise calculations");

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
linkcells(comm),
gr_grid_(NULL),
dgr_grid_(NULL)
{
  // Read in the atoms
  vector< vector<AtomNumber> > all_atoms(3);

  parse("CUTOFF", cutoff_);

  string funcl = getLabel() + ".gr";
  string dfuncl = getLabel() + ".dgr";
  vector<string> gmin, gmax;
  vector<unsigned int> nbin;
  bool b_spline;
  parseVector("RDF_MIN", gmin);
  parseVector("RDF_MAX", gmax);
  parseVector("RDF_NBINS", nbin);
  parse("RDF_BW", rdf_bw_);
  parseFlag("RDF_SPLINE", b_spline);

  vector<string> gnames;
  vector<bool> gpbc;
  gnames.push_back("virial_rdf");
  gpbc.push_back(false);

  //grid wants to have some kind of 
  //string to look up the periodic domain
  //if it is periodic. Since it's not, I'll just 
  //pass some string vectors I have laying around
  //to satisfy the constructor.
  gr_grid_ = new Grid(funcl, gnames, gmin, gmax, nbin, b_spline, true, true, gpbc, gnames, gnames);
  dgr_grid_ = new Grid(dfuncl, gnames, gmin, gmax, nbin, b_spline, true, true, gpbc, gnames, gnames);

  rdf_kernel_support_.push_back(ceil(rdf_bw_ / gr_grid_->getDx()[0]));
  

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


  linkcells.setCutoff( cutoff_ );

  addValueWithDerivatives();
  setNotPeriodic();


  rdffile_.link(*this);
  mf_file_.link(*this);
  rdffile_.open("rdf.dat");
  mf_file_.open("drdf.dat");

}

  Virial::~Virial(){

    if(gr_grid_) {
      free(gr_grid_);
      free(dgr_grid_);
    }

       
    
  }

  void Virial::prepare() {
    requestAtoms(group_1);
    requestAtoms(group_2);
  }
  
void Virial::setup_link_cells_(){
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

  double Virial::eval_rdf_(const double x,
			   const double r,
			   const double scale,
			   vector<double>& der,
    			   vector<double>& der2) {
    const double A = 4 * pi * 35. / 32. * scale;
    //code for these was autogenerated from a symbolic math engine
    der[0] = 2*A*pow(pow(rdf_bw_, 2) - pow(r - x, 2), 2)*(-pow(rdf_bw_, 2) - 3*r*(r - x) + pow(r - x, 2))/(pow(rdf_bw_, 6)*pow(r, 3));
    der2[0] = 6*A*(pow(rdf_bw_, 2) - pow(r - x, 2))*(pow(r, 2)*(-pow(rdf_bw_, 2) + 5*pow(r - x, 2)) + 4*r*(pow(rdf_bw_, 2) - pow(r - x, 2))*(r - x) + pow(pow(rdf_bw_, 2) - pow(r - x, 2), 2))/(pow(rdf_bw_, 6)*pow(r, 4));
    return A*pow(1 - pow(-r + x, 2)/pow(rdf_bw_, 2), 3)/pow(r, 2);
  }

void Virial::calculate()  {
  setup_link_cells_();

  //number of neighbors
  unsigned int nn;
  Vector rij;
  double r, dx, gr;
  Grid::index_t ineigh;
  vector<double> rvec(1);
  vector<Grid::index_t> grid_neighs;
  vector<double> grid_der(1);
  vector<double> grid_der2(1);
  vector<double> x(1);


  double count = 0;
  dx = gr_grid_->getDx()[0];
  
  for(unsigned int i = 0; i < group_1.size(); ++i) {

    nn = 1;
    neighs[0] = i;
    linkcells.retrieveNeighboringAtoms( getPosition(group_1[i]), nn, neighs);
    for(unsigned int j = 0; j < nn; ++j) {
      if(!b_self_virial && neighs[j] < group_1.size())
	continue;
      rij = delta(getPosition(group_1[i]), getPosition(group_2[neighs[j]]));
      r = rij.modulo();
      rvec[0] = r;
      
      grid_neighs = gr_grid_->getNeighbors(rvec, rdf_kernel_support_);
      for(unsigned k = 0; k < grid_neighs.size();++k){
        ineigh = grid_neighs[k];
        gr_grid_->getPoint(ineigh,x);
        gr = eval_rdf_(x[0], r, dx, grid_der, grid_der2);
        gr_grid_->addValueAndDerivatives(ineigh, gr, grid_der);
	dgr_grid_->addValueAndDerivatives(ineigh, grid_der[0], grid_der2);
      }
      count++;
      //do thing tiwht r
      log.printf("pairwise = %f\n", r);
    }
  }

  setValue(0);
  
  rdffile_.rewind();
  mf_file_.rewind();
  gr_grid_->scaleAllValuesAndDerivatives(1 / count);
  dgr_grid_->scaleAllValuesAndDerivatives(4. / 3 * pi * pow(cutoff_,3) / count);
  gr_grid_->writeToFile(rdffile_);
  dgr_grid_->writeToFile(mf_file_);
}


}
}


