/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "bias/Bias.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "core/Atoms.h"
#include "tools/Communicator.h"
#include "tools/File.h"

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS OPES
/*
\par Examples

OPES ...
  LABEL=opes
  ARG=cv
  MODE=CONVERGE
  PACE=500
  SIGMA=0.2
  BARRIER=15
... OPES


*/
//+ENDPLUMEDOC

class OPES : public Bias {

private:
  bool isFirstStep_;
  bool afterCalculate_;
  unsigned NumParallel_;
  unsigned rank_;
  unsigned NumWalkers_;
  unsigned walker_rank_;
  unsigned ncv_;
  unsigned long counter_;

  double kbt_;
  double biasfactor_;
  double bias_prefactor_;
  unsigned stride_;
  std::vector<double> sigma0_;
  std::vector<double> sigma_min_;
  bool fixed_sigma_;
  double epsilon_;
  double prob_norm_;
  double av_weights_;
  double av_weights2_;
  double max_max_core_p_;
  double max_core_p_;
  double current_bias_;
  double tau_;
  bool exploration_;

  double threshold2_;
  bool recursive_merge_;
//kernels for now are diagonal truncated Gaussians
  struct kernel
  {
    double height;
    std::vector<double> center;
    std::vector<double> sigma;

    inline void merge_me_with(const kernel & );
    kernel(double h, const std::vector<double> & c,const std::vector<double> & s):
      height(h),center(c),sigma(s) {}
  };
  double cutoff2_;
  double val_at_cutoff_;
  inline double evaluateKernel(const kernel&,const std::vector<double>&) const;
  inline double evaluateKernel(const kernel&,const std::vector<double>&,std::vector<double>&);
  std::vector<kernel> kernels_;
  OFile kernelsOfile_;

  double work_;
  double old_prob_norm_;
  std::vector<kernel> delta_kernels_;

  OFile probOfile_;
  int wProbStride_;
  bool storeOldProb_;

public:
  OPES(const ActionOptions&);
  void calculate() override;
  void update() override;
  double getProbAndDerivatives(const std::vector<double>&,std::vector<double>&);
  void addKernel(const kernel&,const bool);
  unsigned getMergeableKernel(const std::vector<double>&,const unsigned);
  void dumpProbToFile();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(OPES,"OPES")

void OPES::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","MODE","choose between CONVERGE or EXPLORE mode");
  keys.add("compulsory","TEMP","-1","temperature. If not specified tries to get it from MD engine");
  keys.add("compulsory","PACE","the frequency for kernel addition");
  keys.add("compulsory","SIGMA","the initial widths of the kernels");
  keys.add("compulsory","BARRIER","0","the free energy barrier to be overcome. It is used to set BIASFACTOR, EPSILON, and KERNEL_CUTOFF to reasonable values");
  keys.add("compulsory","COMPRESSION_THRESHOLD","1","merge kernels if closer than this threshold (ignored if a grid is used)");
//extra options
  keys.add("optional","BIASFACTOR","the \\f$\\gamma\\f$ bias factor used for well-tempered target \\f$p(\\mathbf{s})\\f$."
           " Set to 0 for non-tempered flat target");
  keys.add("optional","EPSILON","the value of the regularization constant for the probability");
  keys.add("optional","KERNEL_CUTOFF","truncate kernels at this distance (in units of sigma)");
  keys.add("optional","FES_MAX","the max value of free energy you want to sample");
  keys.add("optional","SIGMA_MIN","never go below this value of sigma");
  keys.addFlag("FIXED_SIGMA",false,"do not decrease sigma as simulation goes on");
  keys.addFlag("RECURSIVE_MERGE_OFF",false,"do not recursevely attempt kernel merging when a new one is added. This add an overhead, but keeps the total number of kernels lower");
  keys.add("optional","TAU","converge variant: inital exploration time; exploration variant: exponential decaying average time");
//kernels file
  keys.add("compulsory","FILE","KERNELS","a file in which the list of added kernels is stored");
  keys.add("optional","FMT","specify format for KERNELS file");
//save probability estimate (compressed kernels)
  keys.add("optional","PROB_WFILE","the file on which to write the estimated probability");
  keys.add("optional","PROB_WSTRIDE","write the estimated probability to a file every N steps");
  keys.addFlag("STORE_PROB",false,"store all the estimated probability files the calculation generates. They will be deleted if this keyword is not present");
//miscellaneous
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers");
  keys.addFlag("SERIAL",false,"perform calculations in serial. It is faster for small number of kernels e.g. 1D systems");
  keys.use("RESTART");

//output components
  componentsAreNotOptional(keys);
  keys.addOutputComponent("work","default","work done by the last kernel added"); //calculating this maybe is only a useless overhead...
  keys.addOutputComponent("ci","default","convergence indicator: \\f$-\\log \\lange e^{\\beta V} \\rangle\\f$");
  keys.addOutputComponent("maxprob","default","max value of the probability estimate");
  keys.addOutputComponent("neff","default","effective number of samples");
  keys.addOutputComponent("nker","default","total number of compressed kernels employed");
}

OPES::OPES(const ActionOptions&ao)
  : PLUMED_BIAS_INIT(ao)
  , isFirstStep_(true)
  , afterCalculate_(false)
  , counter_(1)
  , stride_(500)
  , prob_norm_(1)
  , av_weights_(1)
  , av_weights2_(1)
  , max_core_p_(0)
  , work_(0)
  , old_prob_norm_(1)
{
  ncv_=getNumberOfArguments();
//set mode
  std::string mode;
  parse("MODE",mode);
  if(mode=="CONVERGE")
    exploration_=false;
  else if(mode=="EXPLORE")
    exploration_=true;
  else
    error("must set either MODE=CONVERGE or MODE=EXPLORE, setting MODE="+mode+" is not an option");

//set kbt_
  const double Kb=plumed.getAtoms().getKBoltzmann();
  double temp=-1;
  parse("TEMP",temp);
  kbt_=Kb*temp;
  if(kbt_<0)
  {
    kbt_=plumed.getAtoms().getKbT();
    plumed_massert(kbt_>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
  }

//other compulsory input
  parse("PACE",stride_);
  parseVector("SIGMA",sigma0_);
  plumed_massert(sigma0_.size()==ncv_,"number of SIGMA parameters does not match number of arguments");
  biasfactor_=0;
  epsilon_=0;
  double cutoff=0;
  double barrier=0;
  parse("BARRIER",barrier);
  if(barrier!=0)
  {
    biasfactor_=barrier/kbt_;
    cutoff=sqrt(2.*barrier/kbt_);
    epsilon_=std::exp(-barrier/kbt_);
  }
  bias_prefactor_=kbt_;
  parse("BIASFACTOR",biasfactor_);
  if(exploration_)
  {
    plumed_massert(biasfactor_>1,"BIASFACTOR must be greater than one, when using EXPLORATION");
    bias_prefactor_*=(biasfactor_-1.);
    for(unsigned i=0; i<ncv_; i++)
      sigma0_[i]*=std::sqrt(biasfactor_); //the sigma of the target is broader F_t(s)=1/gamma*F(s)
    epsilon_=std::exp(-1);
  }
  else
  {
    plumed_massert(biasfactor_==0 || biasfactor_>1,"BIASFACTOR must be zero (for uniform target) or greater than one");
    if(biasfactor_!=0)
      bias_prefactor_*=(1-1./biasfactor_);
  }
  parse("EPSILON",epsilon_);
  plumed_massert(epsilon_>0,"you must choose a value for EPSILON if you do not use the BARRIER keyword");
  parse("KERNEL_CUTOFF",cutoff);
  plumed_massert(cutoff>0,"you must choose a value for KERNEL_CUTOFF if you do not use the BARRIER keyword");
  cutoff2_=cutoff*cutoff;
  val_at_cutoff_=std::exp(-0.5*cutoff2_);
  if(val_at_cutoff_>epsilon_)
    log.printf("+++ WARNING +++ kernels might be trunchated too much for the given epsilon");

//other options
  threshold2_=1;
  parse("COMPRESSION_THRESHOLD",threshold2_);
  threshold2_*=threshold2_;
  if(threshold2_!=0)
    plumed_massert(threshold2_>0 && threshold2_<cutoff2_,"COMPRESSION_THRESHOLD cannot be smaller than KERNEL_CUTOFF");
  max_max_core_p_=0;
  double max_fes=0;
  parse("FES_MAX",max_fes);
  if(max_fes!=0)
  {
    max_max_core_p_=std::exp(max_fes/kbt_)-1;
    if(exploration_)
      max_max_core_p_=std::exp(max_fes/kbt_/biasfactor_)-1;
  }
  parseVector("SIGMA_MIN",sigma_min_);
  plumed_massert(sigma_min_.size()==0 || sigma_min_.size()==ncv_,"number of SIGMA_MIN does not match number of arguments");
  fixed_sigma_=false;
  parseFlag("FIXED_SIGMA",fixed_sigma_);
  bool recursive_merge_off=false;
  parseFlag("RECURSIVE_MERGE_OFF",recursive_merge_off);
  recursive_merge_=!recursive_merge_off;
  tau_=0;
  parse("TAU",tau_);
  if(exploration_)
    tau_=tau_/(getTimeStep()*stride_);

//kernels file
  std::string kernelsFileName("KERNELS");
  parse("FILE",kernelsFileName);
  std::string fmt;
  parse("FMT",fmt);

//save used bias
  std::string probFileName;
  parse("PROB_WFILE",probFileName);
  wProbStride_=0;
  parse("PROB_WSTRIDE",wProbStride_);
  storeOldProb_=false;
  parseFlag("STORE_PROB",storeOldProb_);
  if(wProbStride_!=0 || storeOldProb_)
    plumed_massert(probFileName.length()>0,"filename for estimated probability not specified, use PROB_WFILE");
  if(probFileName.length()>0 && wProbStride_==0)
    wProbStride_=-1;//will print only on CPT events

//multiple walkers //TODO implement also external mw for cp2k
  bool walkers_mpi=false;
  parseFlag("WALKERS_MPI",walkers_mpi);
  if(walkers_mpi)
  {
    if(comm.Get_rank()==0)//multi_sim_comm works on first rank only
    {
      NumWalkers_=multi_sim_comm.Get_size();
      walker_rank_=multi_sim_comm.Get_rank();
    }
    if(comm.Get_size()>1) //if each walker has more than one processor update them all
    {
      comm.Bcast(NumWalkers_,0);
      comm.Bcast(walker_rank_,0);
    }
  }
  else
  {
    NumWalkers_=1;
    walker_rank_=0;
  }

//parallelization stuff
  NumParallel_=comm.Get_size();
  rank_=comm.Get_rank();
  bool serial=false;
  parseFlag("SERIAL",serial);
  if(serial)
  {
    log.printf(" -- SERIAL: running without loop parallelization\n");
    NumParallel_=1;
    rank_=0;
  }

  checkRead();

//TODO add option to restart from dumped PROB file
//restart if needed
  if(getRestart())
  {
    plumed_massert(max_max_core_p_==0 && tau_==0 && sigma_min_.size()==0,"RESTART is not supported when this options are used");
    IFile ifile;
    ifile.link(*this);
    if(NumWalkers_>1)
      ifile.enforceSuffix("");
    if(ifile.FileExist(kernelsFileName))
    {
      ifile.open(kernelsFileName);
      log.printf("  RESTART - make sure all used options are compatible\n");
      log.printf("    Restarting from: %s\n",kernelsFileName.c_str());
      std::string old_mode;
      ifile.scanField("mode",old_mode);
      plumed_massert(old_mode==mode,"Cannot restart directly form a different MODE. Kernels file should be properly prepared");
      double old_biasfactor;
      ifile.scanField("biasfactor",old_biasfactor);
      if(old_biasfactor!=biasfactor_)
        log.printf(" +++ WARNING +++ previous bias factor was %g while now it is %g. diff = %g\n",old_biasfactor,biasfactor_,biasfactor_-old_biasfactor);
      if(exploration_)
        plumed_massert(old_biasfactor==biasfactor_,"in MODE=EXPLORE restarting form different bias factor is not supported");
      double old_epsilon;
      ifile.scanField("epsilon",old_epsilon);
      if(old_epsilon!=epsilon_)
        log.printf(" +++ WARNING +++ previous epsilon was %g while now it is %g. diff = %g\n",old_epsilon,epsilon_,epsilon_-old_epsilon);
      double old_cutoff;
      ifile.scanField("kernel_cutoff",old_cutoff);
      if(old_cutoff!=cutoff)
        log.printf(" +++ WARNING +++ previous kernel_cutoff was %g while now it is %g. diff = %g\n",old_cutoff,cutoff,cutoff-old_cutoff);
      double old_threshold;
      const double threshold=sqrt(threshold2_);
      ifile.scanField("compression_threshold",old_threshold);
      if(old_threshold!=threshold)
        log.printf(" +++ WARNING +++ previous compression_threshold was %g while now it is %g. diff = %g\n",old_threshold,threshold,threshold-old_threshold);
      for(unsigned i=0; i<ncv_; i++)
      {
        if(getPntrToArgument(i)->isPeriodic())
        {
          std::string arg_min,arg_max;
          getPntrToArgument(i)->getDomain(arg_min,arg_max);
          std::string file_min,file_max;
          ifile.scanField("min_"+getPntrToArgument(i)->getName(),file_min);
          ifile.scanField("max_"+getPntrToArgument(i)->getName(),file_max);
          plumed_massert(file_min==arg_min,"mismatch between restart and ARG periodicity");
          plumed_massert(file_max==arg_max,"mismatch between restart and ARG periodicity");
        }
      }
      ifile.allowIgnoredFields(); //this allows for multiple restart, but without checking for consistency between them!
      double time;
      while(ifile.scanField("time",time))
      {
        std::vector<double> center(ncv_);
        std::vector<double> sigma(ncv_);
        double height;
        double logweight;
        for(unsigned i=0; i<ncv_; i++)
          ifile.scanField(getPntrToArgument(i)->getName(),center[i]);
        for(unsigned i=0; i<ncv_; i++)
          ifile.scanField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
        ifile.scanField("height",height);
        ifile.scanField("logweight",logweight);
        ifile.scanField();
        kernel new_kernel(height,center,sigma);
        addKernel(new_kernel,false);
        counter_++;
        const double weight=std::exp(logweight);
        av_weights_+=(weight-av_weights_)/counter_;
        av_weights2_+=(weight*weight-av_weights2_)/counter_;
        if(!exploration_)
          prob_norm_+=weight;
      }
      if(exploration_)
        prob_norm_=counter_;
      log.printf("    A total of %d kernels where read, and compressed to %d\n",counter_-1,kernels_.size());
      ifile.reset(false);
      ifile.close();
    }
    else
      log.printf("+++ WARNING +++ restart requested, but file '%s' was not found!\n",kernelsFileName.c_str());
    //sync all walkers and treads. Not sure is mandatory but is no harm
    comm.Barrier();
    if(comm.Get_rank()==0)
      multi_sim_comm.Barrier();
  }

//setup output kernels file
  kernelsOfile_.link(*this);
  if(NumWalkers_>1)
  {
    if(walker_rank_>0)
      kernelsFileName="/dev/null"; //only first walker writes on file
    kernelsOfile_.enforceSuffix("");
  }
  kernelsOfile_.open(kernelsFileName);
  if(fmt.length()>0)
    kernelsOfile_.fmtField(" "+fmt);
  kernelsOfile_.setHeavyFlush(); //do I need it?
  //define and set const fields
  kernelsOfile_.addConstantField("mode");
  kernelsOfile_.addConstantField("biasfactor");
  kernelsOfile_.addConstantField("epsilon");
  kernelsOfile_.addConstantField("kernel_cutoff");
  kernelsOfile_.addConstantField("compression_threshold");
  for(unsigned i=0; i<ncv_; i++)
    kernelsOfile_.setupPrintValue(getPntrToArgument(i));
  kernelsOfile_.printField("mode",mode);
  kernelsOfile_.printField("biasfactor",biasfactor_);
  kernelsOfile_.printField("epsilon",epsilon_);
  kernelsOfile_.printField("kernel_cutoff",sqrt(cutoff2_));
  kernelsOfile_.printField("compression_threshold",sqrt(threshold2_));

//open file for storing estimated probability
  if(wProbStride_!=0)
  {
    probOfile_.link(*this);
    if(NumWalkers_>1)
    {
      if(walker_rank_>0)
        probFileName="/dev/null"; //only first walker writes on file
      probOfile_.enforceSuffix("");
    }
    probOfile_.open(probFileName);
    if(fmt.length()>0)
      probOfile_.fmtField(" "+fmt);
  }

//add and set output components
  addComponent("work"); componentIsNotPeriodic("work");
  addComponent("ci"); componentIsNotPeriodic("ci");
  getPntrToComponent("ci")->set(-std::log(av_weights_));
  addComponent("maxprob"); componentIsNotPeriodic("maxprob");
  getPntrToComponent("maxprob")->set((max_core_p_+epsilon_)/prob_norm_);
  addComponent("neff"); componentIsNotPeriodic("neff");
  getPntrToComponent("neff")->set(av_weights_*av_weights_/av_weights2_);
  addComponent("nker"); componentIsNotPeriodic("nker");

//printing some info
  if(exploration_)
    log.printf("  MODE=EXPLORE - the estimated probability will be the biased one, not the reweighted one\n");
  else
    log.printf("  MODE=CONVERGE - the estimated probability will be the unbiased one, via reweighting\n");
  log.printf("  Adding new kernels with PACE = %d\n",stride_);
  log.printf("  Kernels have initial sigma = ");
  for(unsigned i=0; i<ncv_; i++)
    log.printf(" %g",sigma0_[i]);
  log.printf("\n");
  if(exploration_)
    log.printf("    the sigma is SIGMA*sqrt(BIASFACTOR)\n");
  if(barrier!=0)
    log.printf("  Expected barrier is %g\n",barrier);
  log.printf("  Kernels are truncated with cutoff = %g\n",cutoff);
  if(cutoff<3.5)
    log.printf("+++ WARNING +++ probably kernels are truncated too much!\n");
  log.printf("      value at cutoff is = %g\n",val_at_cutoff_);
  log.printf("  Regularization epsilon = %g\n",epsilon_);
  if(max_fes!=0)
    log.printf("  The fes will be explored up to FES_MAX = %g  (max_max_core_p_=%g)\n",max_fes,max_max_core_p_);
  if(fixed_sigma_)
    log.printf(" -- FIXED_SIGMA: sigma will not decrease\n");
  if(threshold2_!=0)
    log.printf("  Kernel compression will be used, with threshold = %g\n",sqrt(threshold2_));
  if(!recursive_merge_)
    log.printf(" -- RECURSIVE_MERGE_OFF: only one merge for each new kernel will be attempted. This is faster only if total number of kernels does not grow too much\n");
  if(tau_!=0)
    log.printf("  Increasing exploration by using TAU = %g\n",exploration_?tau_*getTimeStep()*stride_:tau_); //FIXME make it two different keywords?
  if(biasfactor_!=0)
    log.printf("  Using target distribution with bias factor gamma = %g\n",biasfactor_);
  else
    log.printf("  Using flat target distribution, no well-tempering\n");
  log.printf("  Temperature T = %g\n",kbt_/Kb);
  log.printf("  Beta = %g\n",1./kbt_);
  if(wProbStride_!=0 && walker_rank_==0)
    log.printf("  Probability estimate is written on file %s with stride %d\n",probFileName.c_str(),wProbStride_);
  if(walkers_mpi)
    log.printf("  -- WALKERS_MPI: if present, multiple simulations will communicate\n");
  if(NumWalkers_>1)
  {
    log.printf("  Using multiple walkers\n");
    log.printf("    number of walkers: %d\n",NumWalkers_);
    log.printf("    walker rank: %d\n",walker_rank_);
  }
  if(NumParallel_>1)
    log.printf("  Using multiple threads per simulation: %d\n",NumParallel_);
  log.printf(" Bibliography ");
  log<<plumed.cite("Invernizzi and Parrinello, arXiv:1909.07250 (2019)");
}

void OPES::calculate()
{
  std::vector<double> cv(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    cv[i]=getArgument(i);

  std::vector<double> der_prob(ncv_,0);
  const double prob=getProbAndDerivatives(cv,der_prob);
  current_bias_=bias_prefactor_*std::log(prob);
  if(exploration_)
    current_bias_-=kbt_*std::log(av_weights_);

  setBias(current_bias_);
  for(unsigned i=0; i<ncv_; i++)
    setOutputForce(i,-bias_prefactor_/prob*der_prob[i]);

//calculate work
  double tot_delta=0;
  for(unsigned d=0; d<delta_kernels_.size(); d++)
    tot_delta+=evaluateKernel(delta_kernels_[d],cv);
  work_-=bias_prefactor_*std::log(prob_norm_/old_prob_norm_*(1.-tot_delta/(prob_norm_*prob)));

  if( (wProbStride_>0 && getStep()%wProbStride_==0) || (wProbStride_==-1 && getCPT()) )
    dumpProbToFile();
  afterCalculate_=true;
}

void OPES::update() //FIXME should check what happens when one wants to do a rerun for post processing, e.g. to get V_t(s)
{
  if(getStep()%stride_!=0)
    return;
  if(isFirstStep_)//same in MetaD, but is it needed?
  {
    isFirstStep_=false;
    return;
  }
  plumed_massert(afterCalculate_,"OPES::update() must be called after OPES::calculate() to work properly");
  afterCalculate_=false;

//work done by the bias in one iteration, uses as zero reference a point at inf -> V(s)>=0
  const double min_shift=bias_prefactor_*std::log(old_prob_norm_/prob_norm_);
  getPntrToComponent("work")->set(work_-stride_*min_shift);
  work_=0;
  delta_kernels_.clear();
  old_prob_norm_=prob_norm_;

//get new kernel height
  double height=std::exp(current_bias_/kbt_); //this assumes that calculate() always runs before update()
  if(tau_!=0 && !exploration_)
    height*=(1-std::exp(-getTime()/tau_));

//update av_weights_ and neff
  double sum_heights=height;
  double sum_heights2=height*height;
  if(NumWalkers_>1)
  {
    if(comm.Get_rank()==0)
    {
      multi_sim_comm.Sum(sum_heights);
      multi_sim_comm.Sum(sum_heights2);
    }
    if(comm.Get_size()>1)
    {
      comm.Bcast(sum_heights,0);
      comm.Bcast(sum_heights2,0);
    }
  }
  counter_+=NumWalkers_;
  av_weights_+=(sum_heights-av_weights_*NumWalkers_)/counter_;
  av_weights2_+=(sum_heights2-av_weights2_*NumWalkers_)/counter_;
  double neff=av_weights_*av_weights_/av_weights2_*counter_;
  getPntrToComponent("ci")->set(-std::log(av_weights_));
  getPntrToComponent("neff")->set(neff);
  if(exploration_)
  {
    height=1.;
    prob_norm_=counter_;
    if(tau_!=0 && prob_norm_>tau_)
      prob_norm_=tau_; //exponentially decaying average
    neff=prob_norm_;
  }
  else
    prob_norm_+=sum_heights;
  getPntrToComponent("maxprob")->set((max_core_p_+epsilon_)/prob_norm_);

//if needed, rescale sigma and height
  double s_rescaling=1;
  std::vector<double> sigma(ncv_);
  if(!fixed_sigma_)
    s_rescaling=std::pow(neff*(ncv_+2.)/4.,-1./(4+ncv_));
  if(sigma_min_.size()==0)
  {
    for(unsigned i=0; i<ncv_; i++)
      sigma[i]=sigma0_[i]*s_rescaling;
    height/=std::pow(s_rescaling,ncv_);
  }
  else
  {
    for(unsigned i=0; i<ncv_; i++)
    {
      sigma[i]=std::max(sigma0_[i]*s_rescaling,sigma_min_[i]);
      height*=sigma0_[i]/sigma[i];
    }
  }

//get new kernel center
  std::vector<double> center(ncv_);
  for(unsigned i=0; i<ncv_; i++)
    center[i]=getArgument(i);

//add new kernel(s)
  if(NumWalkers_>1)
  {
    std::vector<double> all_height(NumWalkers_,0.0);
    std::vector<double> all_center(NumWalkers_*ncv_,0.0);
    std::vector<double> all_sigma(NumWalkers_*ncv_,0.0);
    if(rank_==0)
    {
      multi_sim_comm.Allgather(height,all_height); //room for improvements: heights are communicated twice!
      multi_sim_comm.Allgather(center,all_center);
      multi_sim_comm.Allgather(sigma,all_sigma);
    }
    if(comm.Get_size()>1)
    {
      comm.Bcast(all_height,0);
      comm.Bcast(all_center,0);
      comm.Bcast(all_sigma,0);
    }
    for(unsigned w=0; w<NumWalkers_; w++)
    {
      std::vector<double> center_w(all_center.begin()+ncv_*w,all_center.begin()+ncv_*(w+1));
      std::vector<double> sigma_w(all_sigma.begin()+ncv_*w,all_sigma.begin()+ncv_*(w+1));
      kernel new_kernel(all_height[w],center_w,sigma_w);
      addKernel(new_kernel,true);
    }
  }
  else
  {
    kernel new_kernel(height,center,sigma);
    addKernel(new_kernel,true);
  }
  getPntrToComponent("nker")->set(kernels_.size());
}

double OPES::getProbAndDerivatives(const std::vector<double> &cv,std::vector<double> &der_prob)
{
  double prob=0.0;
  for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_) //TODO add neighbor list
    prob+=evaluateKernel(kernels_[k],cv,der_prob);
  if(NumParallel_>1)
  {
    comm.Sum(prob);
    comm.Sum(der_prob);
  }

  //set maxprob
  if(prob>max_core_p_)
  {
    max_core_p_=prob;
    if(max_max_core_p_>0 && max_core_p_>max_max_core_p_) //FIXME should say somthing to the user when this happens?
      epsilon_*=max_core_p_/max_max_core_p_; //this assumes min_core_p=epsilon
  }

  //regularize the probability
  prob=(prob+epsilon_)/prob_norm_;
  for(unsigned i=0; i<ncv_; i++)
    der_prob[i]/=prob_norm_;
  return prob;
}

void OPES::addKernel(const kernel &new_kernel,const bool write_to_file)
{
  bool no_match=true;
  if(threshold2_!=0)
  {
    unsigned taker_k=getMergeableKernel(new_kernel.center,kernels_.size());
    if(taker_k<kernels_.size())
    {
      no_match=false;
      delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
      kernels_[taker_k].merge_me_with(new_kernel);
      delta_kernels_.push_back(kernels_[taker_k]);
      if(recursive_merge_) //is the overhead worth it?
      {
      //TODO: this second check could run only through the kernels closer than, say, 2*threshold
      //      the function getMergeableKernel could return a list of such neighbors
        unsigned giver_k=taker_k;
        taker_k=getMergeableKernel(kernels_[giver_k].center,giver_k);
        while(taker_k<kernels_.size())
        {
          delta_kernels_.pop_back();
          delta_kernels_.emplace_back(-1*kernels_[taker_k].height,kernels_[taker_k].center,kernels_[taker_k].sigma);
          if(taker_k>giver_k) //saves time when erasing
            std::swap(taker_k,giver_k);
          kernels_[taker_k].merge_me_with(kernels_[giver_k]);
          delta_kernels_.push_back(kernels_[taker_k]);
          kernels_.erase(kernels_.begin()+giver_k);
          giver_k=taker_k;
          taker_k=getMergeableKernel(kernels_[giver_k].center,giver_k);
        }
      }
    }
  }
  if(no_match)
  {
    kernels_.push_back(new_kernel);
    delta_kernels_.push_back(new_kernel);
  }

//write to file
  if(write_to_file)
  {
    kernelsOfile_.printField("time",getTime());
    for(unsigned i=0; i<ncv_; i++)
      kernelsOfile_.printField(getPntrToArgument(i),new_kernel.center[i]);
    for(unsigned i=0; i<ncv_; i++)
      kernelsOfile_.printField("sigma_"+getPntrToArgument(i)->getName(),new_kernel.sigma[i]);
    kernelsOfile_.printField("height",new_kernel.height);
    kernelsOfile_.printField("logweight",current_bias_/kbt_);
    kernelsOfile_.printField();
  }
}

unsigned OPES::getMergeableKernel(const std::vector<double> &giver_center,const unsigned giver_k)
{ //returns kernels_.size() if no match is found
  unsigned min_k=kernels_.size();
  double min_dist2=threshold2_;
  for(unsigned k=rank_; k<kernels_.size(); k+=NumParallel_) //TODO add neighbor list
  {
    if(k==giver_k) //a kernel should not be merged with itself
      continue;
    double dist2=0;
    for(unsigned i=0; i<ncv_; i++)
    { //TODO implement merging on the border for periodic CVs
      const double d=(kernels_[k].center[i]-giver_center[i])/kernels_[k].sigma[i];
      dist2+=d*d;
      if(dist2>=min_dist2)
        break;
    }
    if(dist2<min_dist2)
    {
      min_dist2=dist2;
      min_k=k;
    }
  }
  if(NumParallel_>1)
  {
    std::vector<double> all_min_dist2(NumParallel_);
    std::vector<unsigned> all_min_k(NumParallel_);
    comm.Allgather(min_dist2,all_min_dist2);
    comm.Allgather(min_k,all_min_k);
    const unsigned best=std::distance(std::begin(all_min_dist2),std::min_element(std::begin(all_min_dist2),std::end(all_min_dist2)));
    if(all_min_dist2[best]<threshold2_)
      min_k=all_min_k[best];
  }
  return min_k;
}

void OPES::dumpProbToFile() //FIXME decide what to print in the header XXX
{
  if(storeOldProb_)
    probOfile_.clearFields();
  else if(walker_rank_==0)
    probOfile_.rewind();

  probOfile_.addConstantField("mode");
  probOfile_.addConstantField("biasfactor");
  probOfile_.addConstantField("epsilon");
  probOfile_.addConstantField("kernel_cutoff");
  probOfile_.addConstantField("compression_threshold");
  probOfile_.addConstantField("ci");
  for(unsigned i=0; i<ncv_; i++) //print periodicity of CVs
    probOfile_.setupPrintValue(getPntrToArgument(i));
  std::string mode("CONVERGE");
  if(exploration_)
    mode="EXPLORE";
  probOfile_.printField("mode",mode);
  probOfile_.printField("biasfactor",biasfactor_);
  probOfile_.printField("epsilon",epsilon_);
  probOfile_.printField("kernel_cutoff",sqrt(cutoff2_));
  probOfile_.printField("compression_threshold",sqrt(threshold2_));
  probOfile_.printField("ci",-std::log(av_weights_));
  for(unsigned k=0; k<kernels_.size(); k++)
  {
    probOfile_.printField("time",getTime());
    for(unsigned i=0; i<ncv_; i++)
      probOfile_.printField(getPntrToArgument(i),kernels_[k].center[i]);
    for(unsigned i=0; i<ncv_; i++)
      probOfile_.printField("sigma_"+getPntrToArgument(i)->getName(),kernels_[k].sigma[i]);
    probOfile_.printField("height",kernels_[k].height);
    probOfile_.printField();
  }
  if(!storeOldProb_)
    probOfile_.flush();
}

inline double OPES::evaluateKernel(const kernel& G,const std::vector<double>& x) const
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  for(unsigned i=0; i<ncv_; i++)
  {
    const double diff_i=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=diff_i*diff_i;
    if(norm2>=cutoff2_)
      return 0;
  }
  return G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
}

inline double OPES::evaluateKernel(const kernel& G,const std::vector<double>& x, std::vector<double> & acc_der)
{ //NB: cannot be a method of kernel class, because uses external variables (for cutoff)
  double norm2=0;
  std::vector<double> diff(ncv_);
  for(unsigned i=0; i<ncv_; i++)
  {
    diff[i]=difference(i,G.center[i],x[i])/G.sigma[i];
    norm2+=diff[i]*diff[i];
    if(norm2>=cutoff2_)
      return 0;
  }
  const double val=G.height*(std::exp(-0.5*norm2)-val_at_cutoff_);
  for(unsigned i=0; i<ncv_; i++)
    acc_der[i]-=diff[i]/G.sigma[i]*val; //NB: we accumulate the derivative into der
  return val;
}

inline void OPES::kernel::merge_me_with(const kernel & other)
{
  const double h=height+other.height;
  for(unsigned i=0; i<center.size(); i++)
  {
    const double c_i=(height*center[i]+other.height*other.center[i])/h;
    const double s_my_part=height*(sigma[i]*sigma[i]+center[i]*center[i]);
    const double s_other_part=other.height*(other.sigma[i]*other.sigma[i]+other.center[i]*other.center[i]);
    const double s2_i=(s_my_part+s_other_part)/h-c_i*c_i;
    center[i]=c_i;
    sigma[i]=sqrt(s2_i);
  }
  height=h;
}

}
}
