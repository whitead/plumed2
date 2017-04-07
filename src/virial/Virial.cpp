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
  double smoothing_;
  double pairwise_density_;

  double kbt_;
  unsigned int pair_count_;
  unsigned int rdf_samples_;

  vector<unsigned int> rdf_kernel_support_;
  //atoms
  vector<AtomNumber> group_1;
  vector<AtomNumber> group_2;
  //used to store neighboring atoms
  vector<unsigned int> neighs;
  //Neighbor list linkscells to accelerate
  LinkCells linkcells;
  OFile rdf_file_;
  OFile drdf_file_;
  OFile virial_file_;


  Grid* rdf_grid_;
  Grid* drdf_grid_;
  Grid* virial_grid_;


  void setup_link_cells_();
  double eval_virial_(const double r,
			     const double rdf,
			     const double drdf,
			     const double d2rdf,
			      double* force) const;
  double eval_rdf_(const double x, 
		   const double r, 
		   const double scale, 
		   vector<double>& der, 
		   vector<double>& der2) const;
  void compute_virial_grid_() const;

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
  keys.add("optional", "TEMPERATURE", "The maximum distance to consider for pairwise calculations");
  keys.add("optional", "AVG_TIME_CONST", "The maximum distance to consider for pairwise calculations");
  keys.add("optional", "VOLUME", "The maximum distance to consider for pairwise calculations");
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
kbt_(1),
rdf_samples_(0),
neighs(0),
linkcells(comm),
rdf_grid_(NULL),
drdf_grid_(NULL),
virial_grid_(NULL)
{
  // Read in the atoms
  vector< vector<AtomNumber> > all_atoms(3);
  double T = 0;
  double time_const = 1;
  double vol = 0;

  parse("CUTOFF", cutoff_);

  string funcl = getLabel() + ".gr";
  string dfuncl = getLabel() + ".dgr";
  string ffuncl = getLabel() + ".force";
  vector<string> gmin, gmax;
  vector<unsigned int> nbin;
  bool b_spline;
  parseVector("RDF_MIN", gmin);
  parseVector("RDF_MAX", gmax);
  parseVector("RDF_NBINS", nbin);
  parse("RDF_BW", rdf_bw_);
  parse("TEMPERATURE", T);
  parse("AVG_TIME_CONST", time_const);
  parse("VOLUME", vol);
  parseFlag("RDF_SPLINE", b_spline);

  if(T>=0.0) kbt_=plumed.getAtoms().getKBoltzmann()*T;
  else kbt_ = plumed.getAtoms().getKbT();

  log.printf("  kBT taken to be %4.2f", kbt_);

  if(time_const == 0) {
    time_const = 1;
  }
  smoothing_ = 1.0 - exp(-1.0 / time_const);
  log.printf("  Using exponential moving average with time constant %4.2f (smoothing = %4.2f)\n", time_const, smoothing_);

  vector<string> gnames;
  vector<bool> gpbc;
  gnames.push_back("virial_rdf");
  gpbc.push_back(false);

  //grid wants to have some kind of 
  //string to look up the periodic domain
  //if it is periodic. Since it's not, I'll just 
  //pass some string vectors I have laying around
  //to satisfy the constructor.
  rdf_grid_ = new Grid(funcl, gnames, gmin, gmax, nbin, b_spline, true, true, gpbc, gnames, gnames);
  drdf_grid_ = new Grid(dfuncl, gnames, gmin, gmax, nbin, b_spline, true, true, gpbc, gnames, gnames);
  virial_grid_ = new Grid(ffuncl, gnames, gmin, gmax, nbin, b_spline, true, true, gpbc, gnames, gnames);

  //rdf_kernel_support_.push_back(ceil(rdf_bw_ / rdf_grid_->getDx()[0]));
  rdf_kernel_support_.push_back(4 * ceil(rdf_bw_ / rdf_grid_->getDx()[0]));
  

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
    log.printf("  Computing virial between two groups\n");
  } else {
    b_self_virial = true;    
    group_1 = all_atoms[0];
    group_2 = all_atoms[0];
    log.printf("  Computing self-virial\n");
  }

  //make neighbors big enough for all neighbors case
  neighs.resize(group_1.size() + group_2.size());
			      
  // And check everything has been read in correctly
  checkRead();

  if(vol == 0) {
    vol = getBox().determinant();
    if(vol == 0) {
      log.printf("  Warning: Unable to determine box, so using cutoff for volume estimate.\n");
      log.printf("    This will lead to an overestimation of density.\n");
      log.printf("    Conisder setting volume (in plumed units) manually\n");
      vol = 4. / 3 * pi * pow(cutoff_, 3);
      log.printf("  Estimating volume as %4.2f\n", vol);    
    }
  }
  //calculate pairwise density
  //divide by 2 since we don't double count
  pairwise_density_ = (group_1.size() * (group_2.size() - 1)) / 2.0f / vol;
  log.printf("  Setting pairwise density to %4.2f\n", pairwise_density_);



  linkcells.setCutoff( cutoff_ );

  addValueWithDerivatives();
  setNotPeriodic();


  rdf_file_.link(*this);
  drdf_file_.link(*this);
  virial_file_.link(*this);
  rdf_file_.open("rdf.dat");
  drdf_file_.open("drdf.dat");
  virial_file_.open("virial.dat");

}

  Virial::~Virial(){

    if(rdf_grid_) {
      free(rdf_grid_);
      free(drdf_grid_);
    }
    if(virial_grid_)
      free(virial_grid_);
       
    
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


  /*double Virial::eval_rdf_(const double x,
			   const double r,
			   const double scale,
			   vector<double>& der,
    			   vector<double>& der2) const {
    const double A = 4 * pi * 35. / 32. * scale;
    //code for these was autogenerated from a symbolic math engine
    der[0] = 2*A*pow(pow(rdf_bw_, 2) - pow(r - x, 2), 2)*(-pow(rdf_bw_, 2) - 3*r*(r - x) + pow(r - x, 2))/(pow(rdf_bw_, 6)*pow(r, 3));
    der2[0] = 6*A*(pow(rdf_bw_, 2) - pow(r - x, 2))*(pow(r, 2)*(-pow(rdf_bw_, 2) + 5*pow(r - x, 2)) + 4*r*(pow(rdf_bw_, 2) - pow(r - x, 2))*(r - x) + pow(pow(rdf_bw_, 2) - pow(r - x, 2), 2))/(pow(rdf_bw_, 6)*pow(r, 4));
    return A*pow(1 - pow(-r + x, 2)/pow(rdf_bw_, 2), 3)/pow(r, 2);
    }
  */

    double Virial::eval_rdf_(const double x,
			   const double r,
			   const double scale,
			   vector<double>& der,
    			   vector<double>& der2) const {
      const double A = 1. / (4 * pi * scale) / sqrt(2 * pi) / rdf_bw_ / pairwise_density_;
      const double u = (x - r) * (x - r) / rdf_bw_ / rdf_bw_;
      const double g = A * exp(-u / 2) / r/ r;
      der[0] = g * (r - x) / rdf_bw_ / rdf_bw_;
      der2[0] = g  * (u - 1) / rdf_bw_ / rdf_bw_;
      return g;
    }
  
  double Virial::eval_virial_(const double r,
			     const double rdf,
			     const double drdf,
			     const double d2rdf,
			     double* force) const {

    *force = (kbt_ / r * (d2rdf / rdf - drdf*drdf / rdf / rdf - drdf / r / rdf));
    return kbt_ / rdf * drdf;
  }
    
			     

  void Virial::compute_virial_grid_() const {
    //compute and write out pairwise force
    vector<double> x(1);
    vector<double> der(1);
    vector<double> der2(1);
    double v, f;
    for(Grid::index_t i = 0; i < virial_grid_->getSize(); ++i) {
      rdf_grid_->getPoint(i, x);
      v = rdf_grid_->getValueAndDerivatives(i, der);
      drdf_grid_->getValueAndDerivatives(i, der2);
      //get virial and reuse  v
      v = eval_virial_(x[0], v, der[0], der2[0], &f);
      //reuse der
      der[0] = f;
      virial_grid_->setValueAndDerivatives(i, v, der);
    }

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
  vector<double> temp(1);
  vector<double> x(1);


  dx = rdf_grid_->getDx()[0];
  
  for(unsigned int i = 0; i < group_1.size(); ++i) {

    nn = 1;
    neighs[0] = i;
    linkcells.retrieveNeighboringAtoms( getPosition(group_1[i]), nn, neighs);
    for(unsigned int j = 0; j < nn; ++j) {
      //don't double count
      if(!b_self_virial && neighs[j] < group_1.size())
	continue;
      rij = delta(getPosition(group_1[i]), getPosition(group_2[neighs[j]]));
      r = rij.modulo();
      rvec[0] = r;
      
      grid_neighs = rdf_grid_->getNeighbors(rvec, rdf_kernel_support_);
      for(unsigned k = 0; k < grid_neighs.size();++k){
	//get position of neighbor
        ineigh = grid_neighs[k];
        rdf_grid_->getPoint(ineigh,x);
	
	//eval the rdf and derivatives
        gr = eval_rdf_(x[0], r, dx, grid_der, grid_der2);

	//add them with smoothing factor (exp moving avg)
	gr = gr * smoothing_ + (1 - smoothing_) * rdf_grid_->getValueAndDerivatives(ineigh, temp);
	grid_der[0] = grid_der[0] * smoothing_ + (1 - smoothing_) * temp[0];
        rdf_grid_->addValueAndDerivatives(ineigh, gr, grid_der);

	drdf_grid_->getValueAndDerivatives(ineigh, temp);
	grid_der2[0] = grid_der2[0] * smoothing_ + (1 - smoothing_) * temp[0];
	drdf_grid_->addValueAndDerivatives(ineigh, grid_der[0], grid_der2);
      }
      rdf_samples_++;
    }
  }

  setValue(0);
  
  rdf_file_.rewind();
  drdf_file_.rewind();
  virial_file_.rewind();
  compute_virial_grid_();
  rdf_grid_->writeToFile(rdf_file_);
  drdf_grid_->writeToFile(drdf_file_);
  virial_grid_->writeToFile(virial_file_);

}


}
}


