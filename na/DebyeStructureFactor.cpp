/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"

#include <cmath>
#include <sstream> //std::ostringstream
#include <string.h> //strcasecmp
//#include "Compton_consts.h" //source: http://lammps.sandia.gov/doc/compute_xrd.html

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DEBYE_STRUCTURE_FACTOR
/*
The Debye Structure Factor as defined on [wikipedia](https://en.wikipedia.org/wiki/Structure_factor)
\f[
  DS(q) = 1+ \frac{2}{\sum^N_i f_i^2(q)} \sum^N_{i<j} f_i(q)f_j(q) \frac{sin(q r_{ij})}{q r_{ij}}w(r_{ij})
\f]
\f$r_{ij}\f$ is the distance between atom \f$i\f$ and \f$j\f$, and \f$f_i(q)\f$ is the Compton scattering factor, to be considered in case of X-ray diffraction.
To calculate this factor we referenced to the [LAMMPS module USER-DIFFRACTION](http://lammps.sandia.gov/doc/compute_xrd.html).
The supported atom types can be found there.

Two different kind of Debye Structure Factors are available.
By default periodic boundary conditions are expected and it is compulsory to provide a `CUTOFF` radius \f$R_c\f$ which should be smaller than half the box edge.
Only \r_{ij}<R_c\f$ will contribute and a window function is used to correct for this truncation:
\f[
  w(r_{ij})=\frac{sin(\pi r_{ij}/R_c)}{\pi r_{ij}/R_c}
\f]
If instead the `NOPBC` keyword is specified, no cutoff radius is used, distances \f$r_{ij}\f$ are calculated as if the box was not periodic and the window function is considered constant \f$w(r_{ij})=1\f$.
In this way the pbc symmetry is broken, but a more accurate Debye Structure Factor is calculated, which is useful for certain applications.

The Debye Structure Factor can be obtained also as a function of the scattering angle \f$2\theta\f$, via Bragg law for elastic scattering:
\f[
  \frac{2\pi}{\lambda}sin(\theta)=\frac{q}{2}
\f]
When a value for the diffraction `LAMBDA` is given, the output will be given as a function of the \f$2\theta\f$ angle instead of \f$q\f$.

This CV can work in two modes:
- manual pick active \f$q\f$ or \f$2\theta\f$ (useful when biasing)
- grid of certain \f$q\f$ or \f$2\theta\f$ values (mainly for building the full structure factor)
The grid spacing is automatically chosen using the given `BOX_EDGE` which is taken as reference, doesn't have to be exact.
This spacing can be fine tuned through the `RESOLUTION` keyword.

\par Examples
some usage examples
\plumedfile
manual: DEBYE_STRUCTURE_FACTOR LAMBDA=1.5406 ACTIVE_2THETA=13.5,25 CUTOFF=8
grid_q: DEBYE_STRUCTURE_FACTOR NOBPC BOX_EDGE=16 MAX=9
\endplumedfile

*/
//+ENDPLUMEDOC

struct indexes_pair //auxiliary class for double loop parallelization
{
  unsigned i,j;
  indexes_pair(unsigned _i,unsigned _j) : i(_i),j(_j) {}
};

class DebyeStructureFactor : public Colvar {

private:
  bool no_deriv_;
  unsigned NumParallel_; //number of parallel tasks

  unsigned NumAtom_;
  std::vector< std::vector<indexes_pair> > AtomPair_;

  bool pbc_;
  double cutoff_;

  bool grid_mode_;
  unsigned n_min_;
  double q_const_;

  std::vector<double> active_q_;
  std::vector< std::vector<double> > norm_form_factor_;

  std::vector<Value*> valueDS;

public:
  DebyeStructureFactor(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& );
};

PLUMED_REGISTER_ACTION(DebyeStructureFactor,"DEBYE_STRUCTURE_FACTOR")

void DebyeStructureFactor::registerKeywords(Keywords& keys)
{
  Colvar::registerKeywords(keys);
  keys.add("optional","CUTOFF","cutoff distance, should not be bigger than half of the box size");
  keys.add("numbered","ATOMS","groups of atoms to be used. ATOMS1, ATOMS2,... and then specify the different ATOM_TYPES. Default is all atoms in group 1");
  keys.reset_style("ATOMS","atoms");
  keys.add("optional","ATOM_TYPE","add effect of Compton scattering for given atom type, use only for x-ray experiments");

//set considered q by manually choosing them. useful for biasing
  keys.add("optional","ACTIVE_Q","manually set which q frequencies will be considered");
  keys.add("optional","ACTIVE_2THETA","manually set which frequencies will be considered, by setting the angle in degrees");
  keys.add("optional","LAMBDA","wavelength of incident radiation. Compulsory when using angles instead of frequencies");
//set considered q through a uniform grid. used mainly when running the driver
  keys.add("optional","BOX_EDGE","set a reference value for the edge L of the simulation box. Compulsory for 'grid mode'");
  keys.add("optional","MAX","maximum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
  keys.add("optional","MIN","minimum grid value calculated. If LAMBDA is set, 2theta value is expected, otherwise q value. Default is reasonable when few primitive cells are simulated");
  keys.add("optional","RESOLUTION","change q grid resolution, default is 3. Actual resolution depends on the size of the simulation box");

//some flags to toggle
  keys.addFlag("SERIAL",false,"perform the calculation in serial even if multiple tasks are available");
  keys.addFlag("NO_DERIV",false,"do not calculate derivatives");
  //the NOPBC flag is present by default

//output components
  keys.add("optional","NAME_PRECISION","set the number of digits used for components name");
  keys.addOutputComponent("ds","default","the instantaneous Debye Structure Factor at a given frequency q (or angle 2theta)");
  ActionWithValue::useCustomisableComponents(keys); //needed to have an unknown number of components
}

DebyeStructureFactor::DebyeStructureFactor(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
//parse and initialize:
//- get pbc and cutoff radius
  bool no_pbc=false;
  parseFlag("NOPBC",no_pbc);
  pbc_=!no_pbc;
  cutoff_=-1;
  parse("CUTOFF",cutoff_);
  if (pbc_)
  {
    plumed_massert(cutoff_>0,"when using PBC a CUTOFF must be set");
    log.printf(" -- PBC: distances will be calculated using the minimal image convention\n");
    log.printf("  distance CUTOFF radius: %g\n",cutoff_);
  }
  else
  {
    plumed_massert(cutoff_=-1,"when using NOPBC the CUTOFF cannot be set");
    log.printf(" -- NOPBC: distances will be calculated ignoring pbc\n");
  }

//- get active q frequencies
  double lambda=-1;
  parse("LAMBDA",lambda);
  const auto from_2theta_to_q=[&lambda](const double _2theta)
    { plumed_massert(_2theta>0 && _2theta<=180,"2theta must be between 0 and 180 degrees");
      return 4*PLMD::pi/lambda*sin(_2theta*PLMD::pi/360); };
  std::vector<double> active_2theta;
  parseVector("ACTIVE_2THETA",active_2theta);
  parseVector("ACTIVE_Q",active_q_);
  if (active_q_.size()>0)
  {
    plumed_massert(lambda==-1,"when manually setting ACTIVE_Q, no LAMBDA is needed");
    plumed_massert(active_2theta.size()==0,"either set ACTIVE_Q or ACTIVE_2THETA or none");
  }
  else if (active_2theta.size()>0)
  {
    plumed_massert(lambda!=-1,"a LAMBDA is needed in order to set ACTIVE_2THETA");
    active_q_.resize(active_2theta.size());
    for (unsigned q=0; q<active_q_.size(); q++)
      active_q_[q]=from_2theta_to_q(active_2theta[q]);
    log.printf("  converting 2theta values to q=4*pi/lambda*sin(2theta*pi/360), using LAMBDA = %g\n",lambda);
  }

  if (active_q_.size()>0)
  {
    //print info
    log.printf("  using only %d manually selected q:",active_q_.size());
    for (unsigned q=0; q<active_q_.size(); q++)
      log.printf("  %g",active_q_[q]);
    log.printf("\n");
    grid_mode_=false;
  }
  else
    grid_mode_=true;

//- get q frequency grid
  double box_edge=-1;
  double max=-1;
  double min=-1;
  unsigned resolution=0;
  parse("BOX_EDGE",box_edge);
  parse("MIN",min);
  parse("MAX",max);
  parse("RESOLUTION",resolution);
  if (grid_mode_)
  {
    plumed_massert(box_edge!=-1,"use either specific ACTIVE_Q or a grid (thus at least BOX_EDGE is needed)");
    if (resolution==0)
      resolution=3; //default resolution
    //set the grid
    q_const_=PLMD::pi/box_edge/resolution;
    if (min==-1)
      min=4*PLMD::pi/box_edge; //below this is not physical
    else if (lambda!=-1)
      min=from_2theta_to_q(min);
    n_min_=std::ceil(min/q_const_);
    unsigned n_max;
    if (max==-1)
      n_max=n_min_+149; //good guess for small simulations
    else
    {
      if (lambda!=-1)
        max=from_2theta_to_q(max);
      n_max=std::floor(max/q_const_);
    }
    active_q_.resize(n_max-n_min_+1,q_const_);
    for (unsigned q=0; q<active_q_.size(); q++)
      active_q_[q]*=(q+n_min_);
    //print grid info
    log.printf("  using a grid on q space:\n");
    log.printf("    reference BOX_EDGE L = %g\n",box_edge);
    log.printf("    MIN q = %g [should be greater than 2pi/L=%g]\n",q_const_*n_min_,2*PLMD::pi/box_edge);
    log.printf("    MAX q = %g\n",q_const_*n_max);
    log.printf("    q RESOLUTION = %d --> %d grid points\n",resolution,active_q_.size());
  }
  else
    plumed_massert(box_edge==-1 && min==-1 && max==-1 && resolution==0,"if specific ACTIVE_Q are given, no grid parameter (BOX_EDGE,MIN,MAX,RESOLUTION) should be set");

//add colvar components
  valueDS.resize(active_q_.size());
  std::string val="q";
  if (lambda!=-1 && active_2theta.size()==0)
  {
    active_2theta.resize(active_q_.size());
    for (unsigned q=0; q<active_q_.size(); q++)
      active_2theta[q]=360/PLMD::pi*asin(active_q_[q]*lambda/(4*PLMD::pi));
    val="2theta";
  }
  std::ostringstream oss;
  unsigned name_precision=7;
  parse("NAME_PRECISION",name_precision);
  oss.precision(name_precision);
  log.printf("  components name are %s values, with NAME_PRECISION = %d\n",val.c_str(),name_precision);
  for (unsigned q=0; q<active_q_.size(); q++)
  {
    oss.str("");
    if (lambda!=-1)
      oss<<"ds-"<<active_2theta[q];
    else
      oss<<"ds-"<<active_q_[q];
    addComponentWithDerivatives(oss.str());
    componentIsNotPeriodic(oss.str());
    valueDS[q]=getPntrToComponent(oss.str());
  }

//get the atoms groups
  std::vector<AtomNumber> atoms;
  std::vector<unsigned> group_sizes;
  for (unsigned a=1;; a++)
  {
    std::vector<AtomNumber> atoms_group;
    parseAtomList("ATOMS",a,atoms_group);
    if (!atoms_group.empty())
    {
      atoms.insert(atoms.end(),atoms_group.begin(),atoms_group.end());
      group_sizes.push_back(atoms_group.size());
    }
    else
      break;
  }
  NumAtom_=atoms.size();
  if (NumAtom_==0)
  {
    NumAtom_=plumed.getAtoms().getNatoms();
    atoms.resize(NumAtom_);
    group_sizes.push_back(NumAtom_);
    for(unsigned j=0; j<NumAtom_; j++)
      atoms[j].setIndex(j);
  }
  requestAtoms(atoms);//this must stay after the addComponentWithDerivatives otherwise segmentation violation
  const unsigned NumTypes=group_sizes.size();
  const auto group_index=[&NumTypes](unsigned m,unsigned n)
    { if (m>n) std::swap(m,n); //just for safety
      return n+m*(NumTypes-1)-m*(m-1)/2; };

//get parallelization stuff
  NumParallel_=comm.Get_size();
  unsigned rank=comm.Get_rank();
  bool serial=false;
  parseFlag("SERIAL",serial);
  if (serial)
  {
    log.printf(" -- SERIAL: running without loop parallelization\n");
    NumParallel_=1;
    rank=0;
  }
//initialize the array of atoms pairs
  AtomPair_.resize(NumTypes*(NumTypes+1)/2);
  std::vector<unsigned> tot_pairs(NumTypes*(NumTypes+1)/2);
  std::vector<unsigned> begin_group(NumTypes,0);
  for (unsigned m=1; m<NumTypes; m++)
    begin_group[m]=begin_group[m-1]+group_sizes[m-1];
  unsigned pivot=0;
  for (unsigned m=0; m<NumTypes; m++)
  {
    tot_pairs[group_index(m,m)]=group_sizes[m]*(group_sizes[m]-1)/2;
    AtomPair_[group_index(m,m)].reserve(tot_pairs[group_index(m,m)]/NumParallel_+1);
    for (unsigned i=0; i<group_sizes[m]; i++)
      for (unsigned j=i+1; j<group_sizes[m]; j++) //same group
      {
        pivot++;
        if ((pivot+rank)%NumParallel_==0)
          AtomPair_[group_index(m,m)].emplace_back(begin_group[m]+i,begin_group[m]+j);
      }
    for (unsigned n=m+1; n<NumTypes; n++)
    {
      tot_pairs[group_index(m,n)]=group_sizes[m]*group_sizes[n];
      AtomPair_[group_index(m,n)].reserve(tot_pairs[group_index(m,n)]/NumParallel_+1);
      for (unsigned i=0; i<group_sizes[m]; i++)
        for (unsigned j=0; j<group_sizes[n]; j++)
        {
          pivot++;
          if ((pivot+rank)%NumParallel_==0)
            AtomPair_[group_index(m,n)].emplace_back(begin_group[m]+i,begin_group[n]+j);
        }
    }
  }
  log.printf("  over a total of N_tot=%d, considering a number of atoms N=%d\n",plumed.getAtoms().getNatoms(),NumAtom_);
  log.printf("    for a total number of ordered pairs equal to %d\n",pivot);
  if(NumParallel_>1)
    log.printf("    redistributed over %d processors, each dealing with maximum %d pairs\n",NumParallel_,pivot/NumParallel_+1);
  if (NumTypes>1)
  {
    log.printf("    separated in %d different groups, based on %d atom types\n",AtomPair_.size(),NumTypes);
    log.printf("    in detail:  type_i  type_j  tot_pairs_ij  max_per_processor\n");
    for (unsigned m=0; m<NumTypes; m++) for (unsigned n=m; n<NumTypes; n++)
      log.printf("               %6d  %6d  %11d  %10d\n",m,n,tot_pairs[group_index(m,n)],tot_pairs[group_index(m,n)]/NumParallel_+1);
  }

//get scale corrections
  log.printf("  Form factor:\n");
  std::vector< std::vector<double> > form_factor(AtomPair_.size()); //f_i*f_j
  for (unsigned g=0; g<AtomPair_.size(); g++)
    form_factor[g].resize(active_q_.size(),1);
  //Compton scattering (to be used only for x-rays)
  std::vector<std::string> atom_type;
  parseVector("ATOM_TYPE",atom_type);
  if(atom_type.size()==0)
  {
    log.printf("    no ATOM_TYPE specified, no Compton scattering factor: neutron diffraction is assumed\n");
    plumed_massert(NumTypes==1,"if different ATOMS groups are given, their ATOM_TYPE should be specified");
  }
  else
  {
    error("'Compton_consts.h' is not included");
/*
    plumed_massert(atom_type.size()==NumTypes,"for each atom type specify both ATOM_TYPE and ATOMSn");
    log.printf("    considering Compton scattering factor (x-ray only):\n");
    std::vector< std::vector<double> > compton_scat(NumTypes);
    for (unsigned m=0; m<NumTypes; m++)
    {
      compton_scat[m].resize(active_q_.size(),0);
      unsigned a_type=XRDmaxType;//check if provided atom_type is known
      for (unsigned x=0; x<XRDmaxType; x++)
      {
        if (strcasecmp(atom_type[m].c_str(),XRDtypeList[x])==0)
        {
          a_type=x;
          break;
        }
      }
      plumed_massert(a_type<XRDmaxType,"no match found for ATOM_TYPE="+atom_type[m]);
      for (unsigned q=0; q<active_q_.size(); q++)
      {
        const double argument=pow(active_q_[q]/(4*PLMD::pi),2);
        for (unsigned C=0; C<8; C+=2)
          compton_scat[m][q]+=ASFXRD[a_type][C]*exp(-1*ASFXRD[a_type][C+1]*argument);
        compton_scat[m][q]+=ASFXRD[a_type][8];
      }
      log.printf("      for ATOM_TYPE=%s\n",XRDtypeList[a_type]);
    }
    for (unsigned m=0; m<NumTypes; m++)
      for (unsigned n=m; n<NumTypes; n++)
        for (unsigned q=0; q<active_q_.size(); q++)
          form_factor[group_index(m,n)][q]*=compton_scat[m][q]*compton_scat[n][q];
*/
  }

  //for the structure factor the shift is 1, here must be properly scaled
  std::vector<double> normalization(active_q_.size(),0);
  for (unsigned q=0; q<active_q_.size(); q++)
    for (unsigned m=0; m<NumTypes; m++)
      normalization[q]+=form_factor[group_index(m,m)][q]*group_sizes[m]/2; //ordered pairs count twice

  //apply normalization
  norm_form_factor_.resize(AtomPair_.size());
  for (unsigned g=0; g<AtomPair_.size(); g++)
  {
    norm_form_factor_[g].resize(active_q_.size());
    for (unsigned q=0; q<active_q_.size(); q++)
      norm_form_factor_[g][q]=form_factor[g][q]/normalization[q];
  }

  //print everything, it might be useful
  log.printf("   ~|# q  normalization  form_factor(s)\n");
  for (unsigned q=0; q<active_q_.size(); q++)
  {
    log.printf("   ~|%9.5f  %10f",active_q_[q],normalization[q]*2);
    for (unsigned g=0; g<AtomPair_.size(); g++)
      log.printf("  %10f",form_factor[g][q]);
    log.printf("\n");
  }

//get last flag
  no_deriv_=false;
  parseFlag("NO_DERIV",no_deriv_);
  if (no_deriv_)
    log.printf("  --- NO_DERIV: caution, derivatives are turned off! no bias expected\n");

//parsing finished
  checkRead();
}

void DebyeStructureFactor::calculate()
{
//Calculate the Debye Structure Factor components
  std::vector<double> DebyeS(active_q_.size(),0.);
  std::vector<double> d_DebyeS;
  std::vector<double> virial;
  if (!doNotCalculateDerivatives() && !no_deriv_)
  {
    d_DebyeS.resize(3*NumAtom_*active_q_.size(),0.);
    virial.resize(3*3*active_q_.size(),0.);
  }
  for (unsigned g=0; g<AtomPair_.size(); g++) //loop over possible combinations of atoms group (each has a different norm_form_factor_)
  {
    std::vector<double> DebyeS_group(active_q_.size(),0.);
    for (unsigned h=0; h<AtomPair_[g].size(); h++)
    {
      Vector vR_ij;
      if(pbc_)
        vR_ij=pbcDistance(getPosition(AtomPair_[g][h].j),getPosition(AtomPair_[g][h].i));
      else
        vR_ij=getPosition(AtomPair_[g][h].i)-getPosition(AtomPair_[g][h].j);
      const double R_ij=vR_ij.modulo();
      double window=1;
      double d_window=0;
      if (pbc_)
      {
        if (R_ij>=cutoff_) //TODO implement neighbor list
          continue;
        //window correction function
        const double window_arg=PLMD::pi*R_ij/cutoff_;
        window=sin(window_arg)/window_arg;
        d_window=(window_arg*cos(window_arg)-sin(window_arg))/window_arg/R_ij;
      }
      //trigonometric tricks
      double base_cos=0;
      double base_sin=0;
      double prev_cos=0;
      double prev_sin=0;
      if (grid_mode_)
      {
        base_cos=cos(q_const_*R_ij);
        base_sin=sin(q_const_*R_ij);
        prev_cos=cos(q_const_*R_ij*(n_min_-1));//room for improvements?
        prev_sin=sin(q_const_*R_ij*(n_min_-1));
      }
      for (unsigned q=0; q<active_q_.size(); q++)
      {
        const double QR=active_q_[q]*R_ij;
        double cosQR;
        double sinQR;
        if (grid_mode_)
        {
          cosQR=base_cos*prev_cos-base_sin*prev_sin;
          sinQR=base_cos*prev_sin+base_sin*prev_cos;
          prev_cos=cosQR;
          prev_sin=sinQR;
        }
        else
        {
          sinQR=sin(QR);
          cosQR=cos(QR);
        }
        DebyeS_group[q]+=window*sinQR/QR;
        //get derivatives
        if (!doNotCalculateDerivatives() && !no_deriv_)
        {
          const Vector vDeriv_ij=vR_ij*(norm_form_factor_[g][q]*(window*(QR*cosQR-sinQR)/R_ij+d_window*sinQR)/QR/R_ij);
          for (unsigned l=0; l<3; l++)
            d_DebyeS[3*(q*NumAtom_+AtomPair_[g][h].i)+l]+=vDeriv_ij[l];
          for (unsigned l=0; l<3; l++)
            d_DebyeS[3*(q*NumAtom_+AtomPair_[g][h].j)+l]-=vDeriv_ij[l];
          if (pbc_)
          {
            for(unsigned l=0; l<3; l++)
              for(unsigned m=0; m<3; m++)
                virial[q*9+l*3+m]-=vR_ij[l]*vDeriv_ij[m];
          }
        }
      }
    }
    for (unsigned q=0; q<active_q_.size(); q++)
      DebyeS[q]+=norm_form_factor_[g][q]*DebyeS_group[q];
  }
  if (NumParallel_>1)
  {
    comm.Sum(DebyeS);
    if (!doNotCalculateDerivatives() && !no_deriv_)
    {
      comm.Sum(d_DebyeS);
      if(pbc_)
        comm.Sum(virial);
    }
  }

  //set the components values
  for (unsigned q=0; q<active_q_.size(); q++)
  {
    valueDS[q]->set(1+DebyeS[q]);
    if (!doNotCalculateDerivatives() && !no_deriv_)
    {
      for(unsigned ii=0; ii<3*NumAtom_; ii++)
        valueDS[q]->setDerivative(ii,d_DebyeS[3*NumAtom_*q+ii]);
      if (pbc_)
      {
        for (unsigned ll=0; ll<3*3; ll++)
          valueDS[q]->setDerivative(3*NumAtom_+ll,virial[q*9+ll]);
      }
      else
        setBoxDerivativesNoPbc(valueDS[q]);
    }
  }
}

} //colvar namespace
} //PLMD namespace
