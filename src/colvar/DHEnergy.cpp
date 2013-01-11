/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "CoordinationBase.h"
#include "tools/SwitchingFunction.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR DHENERGY
/*
Calculate Debye-Huckel interaction energy among GROUPA and GROUPB. The two groups should
be disjointed.

\par Examples
\verbatim
# this is printing the electrostatic interaction between two groups of atoms
dh: DHEN GROUPA=1-10 GROUPB=11-20 EPS=80.0 I=0.1 TEMP=300.0
PRINT ARG=dh
\endverbatim
(see also \ref PRINT)

*/
//+ENDPLUMEDOC
   
class DHEnergy : public CoordinationBase{
  double k; // Inverse Debye screening length
  double constant;
  double epsilon;

public:
  DHEnergy(const ActionOptions&);
// active methods:
  static void registerKeywords( Keywords& keys );
  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const;
};

PLUMED_REGISTER_ACTION(DHEnergy,"DHENERGY")

void DHEnergy::registerKeywords( Keywords& keys ){
  CoordinationBase::registerKeywords(keys);
  keys.add("compulsory","I","1.0","Ionic strength (M)");
  keys.add("compulsory","TEMP","300.0","Simulation temperature (K)");
  keys.add("compulsory","EPSILON","80.0","Dielectric constant of solvent");
}

  /*
  Global constants in SI unit used in this calculation:
        N_A = 6.0221412927 * 10^(23) mol^(-1) : Avogadro number
        q = 1.60217656535 * 10^(-19) C : proton charge
        e_0 = 8.854187817620 * 10^(-12) C^2/(N*m^2) : vacuum's dielectric constant         
        k_B = 1.380648813 * 10^(-23) N*m/K : Boltzmann constant
  In SI unit, Debye Huckel CV is defined as:
        DHen = \sum_i\sum_j (q_i*q_j*q^2*N_A)/(4*pi*eps*e_0) * exp(-k*|f_ij|)/(|f_ij|) 
               + \sum_i\sum_j (q_i*q_j*q^2*N_A)/(4*pi*epp*e_0) * (1/|r_ij| - 1/|f_ij|)
             = (q^2*N_A)/(4*pi*e_0) * \sum_i\sum_j q_i*q_j * (exp(-k*|f_ij|)/(eps*|f_ij|) + 1/epp*(1/|r_ij| - 1/|f_ij|))
  (in which |f_ij| = \sqrt(|r_ij|^2+\sigma_i*\sigma_j*exp(-|r_ij|^2/4*\sigma_i*\sigma_j)), 
   \sigma_i and \sigma_j are the effective Born radius.)   
  For an efficient calculation, we group constants and variables into groups:
        constant = (q^2*N_A)/(4*pi*e_0)
        tmp = 1/eps*exp(-k*|f_ij|)/(|f_ij|) + 1/epp*(1/|r_ij| - 1/|f_ij|)

  To speed up the loop calculation, constant can be modified as followed:
        constant= (q^2*N_A)/(4*pi*e_0*10^(-9))*10^(-3) (kJ/mol)
                = ((1.60217656535*10^(-19))^2*6.0221412927*10^(23)*10^(-3))/(4*3.14159265*8.854187817620*10^(-12)*10^(-9))
                = 138.935458111 (kJ/mol)
  
  */

DHEnergy::DHEnergy(const ActionOptions&ao):
Action(ao),
CoordinationBase(ao),
k(0.0),
constant(0.0)
{
  double I,T;
  parse("I",I);
  parse("TEMP",T);
  parse("EPSILON",epsilon);
  checkRead();
  plumed_assert(!plumed.getAtoms().usingNaturalUnits());
  atoms.getUnits().getLength();
  constant=138.935458111/atoms.getUnits().getEnergy();
  k=sqrt(I/(epsilon*T))*502.903741125;
  checkRead();
  log<<"  with solvent dielectric constant "<<epsilon<<"\n";
  log<<"  at temperature "<<T<<" K\n";
  log<<"  at ionic strength "<<I<< "M\n";
  log<<"  Bibliography "<<plumed.cite("Trang, Carloni, Varani and Bussi, submitted (2013)");
}

double DHEnergy::pairing(double distance,double&dfunc,unsigned i,unsigned j)const{
  if(getAbsoluteIndex(i)==getAbsoluteIndex(j)){
    dfunc=0.0;
    return 0.0;
  }
  double invdistance=1.0/distance;
  double tmp=exp(-k*distance)*invdistance*constant*getCharge(i)*getCharge(j)/epsilon;
  double dtmp=-(k+invdistance)*tmp;
  dfunc=dtmp*invdistance;
  return tmp;
}

}

}
