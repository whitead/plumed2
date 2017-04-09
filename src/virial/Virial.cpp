#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "tools/Pbc.h"
#include "tools/LinkCells.h"
#include "tools/AtomNumber.h"
#include "tools/Grid.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Units.h"

//for MPI communicator
#include "core/PlumedMain.h"

#include <string>
#include <cmath>
#include <iostream>

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

  class Virial : public ActionAtomistic, public ActionWithValue {
  private:
    //if we are computing self-virial
    bool b_self_virial_;
    bool b_rdf_first_update_;
    double rdf_bw_;
    double rdf_cutoff_;
    double virial_cutoff_;
    double smoothing_;
    double pairwise_density_;    
    double temperature_;
    double kbt_;
    double virial_scaling_;

    unsigned int rdf_stride_;
    unsigned int grid_w_stride_;
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
    
    void write_grid_();
    double do_calc_();
    
  public:
    static void registerKeywords( Keywords& keys );
    explicit Virial(const ActionOptions&);
    ~Virial();
    // active methods:
    void calculate();
    void prepare();
    void apply();
    virtual unsigned int getNumberOfDerivatives();
  };
  
PLUMED_REGISTER_ACTION(Virial,"VIRIAL")
  
  void Virial::registerKeywords( Keywords& keys ){
    Action::registerKeywords( keys );
    ActionWithValue::registerKeywords( keys );
    ActionAtomistic::registerKeywords( keys );
    
    keys.add("compulsory", "CUTOFF", "The maximum distance to consider for pairwise calculations");
    keys.add("compulsory", "RDF_MIN", "The maximum grid value for the running estimate of RDF");
    keys.add("compulsory", "RDF_MAX", "The minimum grid value for the running estimate of RDF");
    keys.add("compulsory", "RDF_NBINS", "The number of bins to use for the running estimate of the RDF");
    keys.add("optional", "RDF_BW", "The kernel bandwidth for running estimate of RDF");
    keys.add("optional", "TEMPERATURE", "The temperature used for calculating Boltzmann's constant");
    keys.add("optional", "AVG_TIME_CONST", "The exponential moving average RDF time constant. "
	     "Corresponds to how long old RDF datapoints influence it.");
    keys.add("optional", "VOLUME", "Used to manually specify the volume for computing number density. "
	     "If unspecified, the simulation box volume will be used if available and failing that, "
	     "a sphere with radius equal to CUTOFF is used." );
    keys.add("optional", "RDF_STRIDE", "How often add to the RDF.");
    keys.add("optional", "GRID_WRITE_STRIDE", 
	     "How often to write the RDF, its derivative, and pairwise virial grids."
	     "Currently their names are fixed. ");
    keys.addFlag("RDF_SPLINE", true, "Whether or not to use spline interpolation for computing RDF estimate");
    keys.addFlag("MEAN_FIELD", false, 
		 "Instead of computing pairwise distances at each time step, "
		 "replace them with the mean-field approximate radial density");
    
    
    keys.add("atoms-1","GROUP","Calculate the virial between each distinct pair of atoms in the group");
    keys.add("atoms-2","GROUPA","Calculate the virial between all the atoms in GROUPA and all "
	     "the atoms in GROUPB. This must be used in conjuction with GROUPB.");
    keys.add("atoms-2","GROUPB","Calculate the virial between all the atoms in GROUPA and all the atoms "
	     "in GROUPB. This must be used in conjuction with GROUPA.");
    
    keys.addOutputComponent("virial", "default", "The value of the virial collective variable (can be biased)");
    keys.addOutputComponent("pressure", "default", "The pressure computed from the virial (only available for within atom group virial)");
  }
  
  Virial::Virial(const ActionOptions&ao):
    Action(ao),
    ActionAtomistic(ao),
    ActionWithValue(ao),
    b_self_virial_(false),
    b_rdf_first_update_(true),
    rdf_bw_(0),
    temperature_(0),
    kbt_(1),
    rdf_stride_(250),
    grid_w_stride_(1000),
    rdf_samples_(0),
    neighs(0),
    linkcells(comm),
    rdf_grid_(NULL),
    drdf_grid_(NULL),
    virial_grid_(NULL)
  {
    // Read in the atoms
    vector< vector<AtomNumber> > all_atoms(3);
    double time_const = 1;
    double vol = 0;
    
    parse("CUTOFF", virial_cutoff_);
    
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
    parse("RDF_STRIDE", rdf_stride_);
    parse("GRID_WRITE_STRIDE", grid_w_stride_);
    parse("TEMPERATURE", temperature_);
    parse("AVG_TIME_CONST", time_const);
    parse("VOLUME", vol);
    parseFlag("RDF_SPLINE", b_spline);
    
    log.printf("  will compute RDF every %d steps\n", rdf_stride_);
    if(grid_w_stride_ > 0)
      log.printf("  will write grids every %d steps\n", grid_w_stride_);
    
    if(temperature_ >= 0) kbt_ = plumed.getAtoms().getKBoltzmann() * temperature_;
    else kbt_ = plumed.getAtoms().getKbT();
    
    log.printf("  kBT taken to be %4.2f\n", kbt_);
    
    if(rdf_bw_ == 0) {
      log.printf("  Warning: Setting bandwidth to zero means no forces/virial will be calculated.\n"
		 "  Only valid RDF will be output\n");
    }
    if(time_const == 0) {
      time_const = 5;
    }
    smoothing_ = 1.0 - exp(-1.0 / time_const);
    log.printf("  using exponential moving average with time constant %4.2f (smoothing = %4.2f)\n", time_const, smoothing_);
    
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
    //3 -> gaussian stds
    rdf_kernel_support_.push_back(3 * ceil(rdf_bw_ / rdf_grid_->getDx()[0]));
    //we need to go a little bit farther to prevent edge effects in virial
    rdf_cutoff_ = virial_cutoff_ + 3 * rdf_bw_;
    
    
    
    parseAtomList("GROUP", all_atoms[0]);
    parseAtomList("GROUPA", all_atoms[1]);
    parseAtomList("GROUPB", all_atoms[2]);
    
    if(all_atoms[0].size() < 2) {
      if(all_atoms[1].size() < 2 && all_atoms[2].size() < 2) {
	error("Must specify both GROUPA and GROUPB or GROUP");
      }
      b_self_virial_ = false;    
      group_1 = all_atoms[1];
      group_2 = all_atoms[2];
      log.printf("  computing virial between two groups\n");
    } else {
      b_self_virial_ = true;    
      group_1 = all_atoms[0];
      group_2 = all_atoms[0];
      log.printf("  computing self-virial\n");
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
	vol = 4. / 3 * pi * pow(rdf_cutoff_, 3);
	log.printf("  estimating volume...");    
      }
    }
    log.printf("  volume is %4.2f\n", vol);    
    //calculate pairwise density
    //divide by 2 since we don't double count
    pairwise_density_ = (group_1.size() * (group_2.size() - b_self_virial_)) / 2.0f / vol;
    log.printf("  setting pairwise density to %4.2f\n", pairwise_density_);
    
    //compute virial scaling which is kT  /(3 V)
    virial_scaling_ = kbt_ / (3 * vol);
    
    
    
    linkcells.setCutoff( rdf_cutoff_ );
    
    addComponentWithDerivatives("virial");
    if(b_self_virial_){
      addComponentWithDerivatives("pressure");
      componentIsNotPeriodic("pressure");
    }

    componentIsNotPeriodic("virial");
    
    
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
    rdf_file_.close();
    drdf_file_.close();
    virial_file_.close();
    
  }
  
  void Virial::prepare() {
    requestAtoms(group_1);
    requestAtoms(group_2);    
  }
  
  void Virial::setup_link_cells_(){
    std::vector<Vector> ltmp_pos( group_1.size() + group_2.size() * !b_self_virial_ ); 
    std::vector<unsigned> ltmp_ind( group_1.size() + group_2.size() * !b_self_virial_ ); 
    
    
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
    //A contains (1) the guassian denominator (2) local volume
    //pairwise density (3) per-particle-origin averaging prefactor (4) density
    // (5) scale factor for disc grid (which cancles with local vol dx term)
    double A = 1. / (4 * pi * r * r * sqrt(2 * pi)* rdf_bw_ * pairwise_density_);
    double u = (x - r) * (x - r);
    if(rdf_bw_ == 0)
      A = 1 / (4 * pi * r * r * scale * pairwise_density_);
    else
      u /= rdf_bw_ * rdf_bw_;
    const double g = A * exp(-u / 2);
    der[0] = g * (r - x) / rdf_bw_ / rdf_bw_;
    der2[0] = g  * (u - 1) / rdf_bw_ / rdf_bw_;
    return g;
  }
  
  double Virial::eval_virial_(const double r,
			      const double rdf,
			      const double drdf,
			      const double d2rdf,
			      double* force) const {
    
    *force = (virial_scaling_ / r * (d2rdf / rdf - drdf*drdf / rdf / rdf - drdf / r / rdf));
    //positive because we take it for rij vector
    return virial_scaling_ / rdf * drdf;
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
  
  double Virial::do_calc_() {
    
    //number of neighbors
    unsigned int nn;
    Vector ri, rj, rij, force;
    double r, dx, sr, gr, tmp, virial = 0, smooth = smoothing_;
    bool rdf_update = getStep() % rdf_stride_ == 0;
    
    Grid::index_t ineigh;
    vector<Grid::index_t> grid_neighs;
    
    vector<double> rvec(1);
    vector<double> grid_der(1);
    vector<double> grid_der2(1);
    vector<double> temp(1);
    vector<double> x(1);

    //we need this in case our CV is biased
    //to propogate chain rule
    double colvar_force = getPntrToComponent(0)->getForce() + getPntrToComponent(1)->getForce();

    //check if we are on the first update
    if(b_rdf_first_update_) {
      smooth = 1.0;
      b_rdf_first_update_ = false;
    }

    vector<Vector>& forces(modifyForces());
    Tensor virial_tensor(modifyVirial());
    
    if(rdf_update) {
      dx = rdf_grid_->getDx()[0];    
      rdf_grid_->scaleAllValuesAndDerivatives(1 - smoothing_);
      drdf_grid_->scaleAllValuesAndDerivatives(1 - smoothing_);
      rdf_samples_ = 0;
    }
    
    
    for(unsigned int i = 0; i < group_1.size(); ++i) {
      
      nn = 1;
      neighs[0] = i;
      ri = getPosition(group_1[i]);
      linkcells.retrieveNeighboringAtoms( ri, nn, neighs);
      //start at 1 since the first is self
      for(unsigned int j = 1; j < nn; ++j) {
	//don't double count
	if((!b_self_virial_ && neighs[j] < group_1.size()) ||
	   (b_self_virial_ && neighs[j] < i))
	  continue;
	rj = getPosition(group_2[neighs[j]]);
	rij = pbcDistance(ri, rj);
	r = rij.modulo();
	if(r >= rdf_cutoff_)
	  continue;

	rvec[0] = r;
	
	if(r < virial_cutoff_) {
	  //only apply if particle is away from edge (rdf vs virial cutoff)

	  //get mean force radial component
	  sr = 1 / r * virial_grid_->getValueAndDerivatives(rvec, grid_der);

	  //compute forces and virial term
	  //and virial sum terms to my own virial
	  for(unsigned int vi = 0; vi < 3; ++vi) {
	    force[vi] = 0;
	    for(unsigned int vj = 0; vj < 3; ++vj) {
	      //get virial derivative and apply chain rule
	      tmp = colvar_force * rij[vi] * (sr * (vi == vj) + rij[vi] * rij[vj] / r * grid_der[0]);
	      force[vi] += tmp;
	      virial_tensor(vi,vj) += tmp;
	      if(vi == vj) {
		virial += rij[vi] * sr;
	      }
	    }
	  }

	  forces[group_1[i].index()] -= force;
	  forces[group_2[neighs[j]].index()] += force;

	}

	if(rdf_update) {	  
	  grid_neighs = rdf_grid_->getNeighbors(rvec, rdf_kernel_support_);
	  for(unsigned k = 0; k < grid_neighs.size();++k){
	    //get position of neighbor
	    ineigh = grid_neighs[k];
	    rdf_grid_->getPoint(ineigh,x);
	    
	    //eval the rdf and derivatives
	    gr = eval_rdf_(x[0], r, dx, grid_der, grid_der2);
	    
	    //add them with smoothing factor (exp moving avg)
	    //exponential moving average
	    gr *= smooth;
	    grid_der[0] *= smooth;
	    grid_der2[0] *= smooth;
	    
	    rdf_grid_->addValueAndDerivatives(ineigh, gr, grid_der);
	    drdf_grid_->addValueAndDerivatives(ineigh, grid_der[0], grid_der2);
	  }
	  rdf_samples_++;
	}
      }
    }
    
    return virial;
  }
  
  void Virial::write_grid_() {

    rdf_file_ << "#! TIME " << getTime() << "\n";
    drdf_file_ << "#! TIME " << getTime() << "\n";
    virial_file_ << "#! TIME " << getTime() << "\n";
    
    rdf_grid_->writeToFile(rdf_file_);
    drdf_grid_->writeToFile(drdf_file_);
    virial_grid_->writeToFile(virial_file_);
    
  }
  
  void Virial::calculate()  {
    
    setup_link_cells_();
    
    double v = do_calc_();
    getPntrToComponent("virial")->set(v);
    //add ideal gas contribution
    if(b_self_virial_)
      getPntrToComponent("pressure")->set(v + 3 * virial_scaling_ * group_1.size());
    
    if(getStep() % rdf_stride_ == 0) {
      //only changes when rdf is updated
      compute_virial_grid_();	
    }
    if(getStep() % grid_w_stride_ == 0)
      write_grid_();
    
  }
  
  void Virial::apply() {
    //nothing. Handled in actionatomistic
  }

  unsigned int Virial::getNumberOfDerivatives(){
    //plus 9 for virial (?)
    return 3*getNumberOfAtoms() + 9;
  }
  
}
}
