
// mclimit_csm -- use a chisquared minimized over nuisance parameters
// to compute cls.

// Class to run the TMinuit minimization of T. Devlin's chisquared
// defined in CDF 3126, minimized over the nuisance parameters.

#ifndef CSM_H
#define CSM_H

#include "TH1.h"
#include <vector>

using std::vector;

struct svstruct_s
{
  Int_t itemplate;   // which template histo this syst. variation applies to
  char *sysname;     // name of nuisance parameter this corresponds to
  Double_t sysfracl;
  Double_t sysfrach;
  // if there is no shape uncertainty associated with this dependence of a particular
  // template on a nuisance parameter, then these pointers should be set to zero.
  TH1 *lowshape;     // for shape uncertainty -- low histogram shape histo id
  TH1 *highshape;    // for shape uncertainty -- high histogram shape
  Double_t xsiglow;  // how many sigma low lowshape corresponds to (should be a negative number)
  Double_t xsighigh; // how many sigma high highshape corresponds to. (should be a positive number)
};

typedef struct svstruct_s svstruct;

// Constraint equations between nuisance parameters.  Sometimes there are fewer degrees of freedom
// than there are nuisance parameters.  Example:  MET vs. Iso evaluation of the non-W
// background in a signal region, using the 4-sector method.  A*C/B=D relates four
// nuisance parameters down to three: D is computed from A, B, and C.  Note:  multiple constraints between nuisance
// parameters will not be solved for, they will just be evaluated in order.  In pseudexperiments,
// nuisance parameters that are constrained are not randomly chosen, but rather are calculated
// based on the other nuisance parameters.  The arbitrary function here means that no checking
// can be done to make sure that the computed nuisance parameters are physical.

struct npcstruct_s
{
  Int_t ninput;       // number of nuisance paramters input to the calculation of a constrained one
  char **pnameinput;  // the names of the input nuisance parameters
  char *pnameoutput;  // the name of the output parameter
  Double_t (*f)(Double_t*);  // a function taking an array of all the input parameters and computing the output parameter
};

typedef struct npcstruct_s npcstruct;

// one-sided or two-sided 3-sigma or 5-sigma

#define MCLIMIT_CSM_TWOSIDED

// what to do about 1-sided or 2-sided 2-sigmas?

#ifdef MCLIMIT_CSM_TWOSIDED
#define MCLIMIT_CSM_2S 0.02275
#define MCLIMIT_CSM_3S 1.349898E-3
#define MCLIMIT_CSM_5S 2.866516E-7
#else
#define MCLIMIT_CSM_2S 0.0455
#define MCLIMIT_CSM_3S 2.6998E-3
#define MCLIMIT_CSM_5S 5.7330E-7
#endif

// cumulative probabilities for defining bands on test statistic
// and CL plots

#define MCLIMIT_CSM_MCLM2S 0.02275
#define MCLIMIT_CSM_MCLM1S 0.16
#define MCLIMIT_CSM_MCLMED 0.5
#define MCLIMIT_CSM_MCLP1S 0.84
#define MCLIMIT_CSM_MCLP2S 0.97725

// some messages to pass around inside for the s95 calculator

#define MCLIMIT_CSM_CLS 1
#define MCLIMIT_CSM_CLSM2 2
#define MCLIMIT_CSM_CLSM1 3
#define MCLIMIT_CSM_CLSMED 4
#define MCLIMIT_CSM_CLSP1 5
#define MCLIMIT_CSM_CLSP2 6

#define MCLIMIT_CSM_LUMI95 1
#define MCLIMIT_CSM_LUMI3S 2
#define MCLIMIT_CSM_LUMI5S 3

// horizontal is the interpolation style which calls csm_pvmorph (or csm_pvmorph_2d)
// vertical interpolation just does a linear interpolation bin-by-bin, but not
// letting any bin go below zero

typedef enum {
  CSM_INTERP_HORIZONTAL,
  CSM_INTERP_VERTICAL
} INTERPSTYLE;


/* a full model for a particular channel.  csm_model is just a collection
   of channel models and associated names */

class csm_channel_model
{
 public:
   csm_channel_model();
   ~csm_channel_model();
   void add_template( TH1 *,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      TH1 *[],    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation
                      TH1 *[],    // array of high histogram shapes, one for each nuisance param
		      Double_t *, // number of sigma high for each shape variation
                      Int_t,      // Poisson flag -- 1 if Poisson, 0 of not.  There is a split between the
                                  // interpretation of this flag when the model is an ensemble model and when
                                  // it is a test model.  In an ensemble model, the bin contents are treated
                                  // as Poisson means and the values are fluctuated according to Poisson statistics.
                                  // These random components are then used with their scale factors as components
                                  // of another Poisson mean for generating pseudodata. 
                                  // In the model used to fit to the data, these values are treated as integer
                                  // measurements from subsidiary experiments.
                      Int_t);     // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.


   csm_channel_model* Clone();   // make an exact copy of this channel model -- all the internal
                                 // cloned histograms are cloned again.  Better than just
                                 // assigning a new model to this one because the
                                 // destructors won't delete the same memory twice.
   void nuisance_response(Int_t, char *[], Double_t []); // update the internal copies of the varied histograms
                                                         // and scale factors according to the nuisance parameters supplied
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values

   void print();                 // print out some details of the model
   void plotwithdata(TH1*);      // compare data with a model.
   double kstest(TH1*);          // get ROOT's raw KS prob
   double kstest_px(TH1*);       // get ROOT's raw KS px prob
   // adding two models together, and multiplying a model by a scalar

   // note that all of these operations on models create a new model (which must be cleaned up later)

   csm_channel_model* add(csm_channel_model&); // adds the components of two models together
   csm_channel_model* scale(Double_t);        // scales all parts of the model
   csm_channel_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"
   csm_channel_model* scale_err(Double_t);  // scales rates up with scale factor, and scales
                                            // systematic errors down with scale factor.  n.b. --
                                            // MC statistical errors on model histos cannot be scaled
   Double_t chisquared1(TH1 *);  // Inputs a pointer to data histogram -- computes the chisquared of this model
                                 // compared to the supplied data histogram, without minimizing over nuisance
                                 // parameters, a la Tom Devlin's note

   Int_t checkneg();             // check for negative bins

   void set_interpolation_style(INTERPSTYLE);  // either CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL
                                               // horizontal (csm_pvmorph) is the default if this is not called.

   friend class csm_model;       // so we can access systematic error names and limits
   friend class mclimit_csm;

 private:
  vector<TH1*> histotemplate;
  vector<TH1*> histotemplate_varied;
  vector<Double_t> sft;
  vector<Double_t> sft_varied;
  vector<Int_t> poissflag;
  vector<Int_t> scaleflag;
  vector<svstruct> syserr;
  INTERPSTYLE chan_istyle;
};

class csm_model
{
 public:
   csm_model();
   ~csm_model();
   void add_template( TH1 *,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      TH1 *[],    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation (should be negative!)
                      TH1 *[],    // array of high histogram shapes, one for each nuisance param
		      Double_t *, // number of sigma high for each shape variation (should be positive!)
                      Int_t,      // Poisson flag -- 1 if Poisson, 0 of not.  There is a split between the
                                  // interpretation of this flag when the model is an ensemble model and when
                                  // it is a test model.  In an ensemble model, the bin contents are treated
                                  // as Poisson means and the values are fluctuated according to Poisson statistics.
                                  // These random components are then used with their scale factors as components
                                  // of another Poisson mean for generating pseudodata. 
                                  // In the model used to fit to the data, these values are treated as integer
                                  // measurements from subsidiary experiments.
                      Int_t,      // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.
                      char *);    // Channel name

   void set_interpolation_style(char *,INTERPSTYLE);  //  sets the interpolation style 
            //  for a particlar channel -- first arg: channel name
            // of channel -- second arg:  interpolation style:  CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL

   void add_chanmodel(csm_channel_model*,
                      char*); // instead of adding a template at a time, let's add a whole channel's model

   void add_npcons(Int_t, char**, char*, Double_t (*f)(Double_t*)); // for a nuisance parameter which
	    // can be computed given the values of other nuisance parameters, this allows
	    // the user to specify such a funtion.   Arguments;  number of n.p.'s the one to calculate
            // depends on, their names, the name of the n.p. to calculate, and a pointer to the
            // function which does the job.

   void plotwithdata(char*,TH1*); // plot a named channel's model with the supplied data histogram
   double kstest(char*,TH1*); // get raw ROOT ks prob
   double kstest_px(char*,TH1*); // get ROOT's ks test prob with px simulation (stat only)

   csm_model* Clone();   // make an exact copy of this model -- all the internal
                                 // cloned histograms are cloned again.  Better than just
                                 // assigning a new model to this one because the
                                 // destructors won't delete the same memory twice.

   void print();                 // print out some details of the model
   // adding two models together, and multiplying a model by a scalar

   // note that all of these operations on models create a new model (which must be cleaned up later)

   csm_model* add(csm_model&); // adds the components of two models together
   csm_model* scale(Double_t);        // scales all parts of the model
   csm_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"
   csm_model* scale_err(Double_t);  // scales rates up with scale factor, and scales
                                            // systematic errors down with scale factor.  n.b. --
                                            // MC statistical errors on model histos cannot be scaled
   void nuisance_response(Int_t, char *[], Double_t []); // updates the fluctuated version of the histogram templates
                                 // and scale factors inside the channel model
                                 // according to the nuisance parameters provided.  Inputs: number of nuisance
                                 // parameters and their names
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values
   void varysyst();   // randomly choose nuisance parameter values and call nuisance_response

   void single_pseudoexperiment(TH1 *[]); // generate pseudodata for all the channels in this model.
   void list_nparams(vector<char *> *, vector<Double_t> *, vector<Double_t> *); // get a list of
      // unique nuisance parameter names for all the channels in this model and their most
      // restrictive lower and upper bounds.
   Double_t chisquared1(TH1**); // calls the chisquared routine for each channel model.

   friend class mclimit_csm;
   friend class csm;

 private:
  vector<char*> channame;
  vector<csm_channel_model*> chanmodel;
  Int_t lookup_add_channame(char *);  // look up the channel name in the channame vector.  If it's
                                       // not there, add it.
  vector<npcstruct> npcm;  // constraint equations between nuisance parameters
};


class csm
{
  public:
    	csm();
	~csm();
        void set_htofit(TH1*,char*); //histogram and channel name
        void set_modeltofit(csm_model*); // a set of template histograms to fit to the data histograms
	Double_t chisquared();           // calculates chisquared
        Int_t ndof();                    // calculates an (approximate!) NDOF
	Int_t    getnparams();
        Double_t getparam(Int_t);
        Double_t getperror(Int_t);
        char* getpname(Int_t);
        csm_model* getbestmodel(); // a model with the template histograms all normalized to the best fit values.
                                   // suitable for plotting.  This varies the input model (given in set_modeltofit)
                                   //  with the best fit nuisance parameters and returns a pointer to the same model
                                   // given in set_modeltofit.
        void plotcompare(char *); // make a stacked plot of data in the named channel against the central
                                              // value model provided.
 private:
        vector<Double_t> fitparam;   // parameters of the fit
	vector<Double_t> fiterror;   // errors on fit parameters
        vector<char*> fitparamname; //  names of fit parameters
};

void csm_interpolate_histogram(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

void csm_pvmorph(Int_t *nb1, Double_t *xmin1, Double_t *xmax1, Double_t *dist1,
                 Int_t *nb2, Double_t *xmin2, Double_t *xmax2, Double_t *dist2,
                 Int_t *nbn, Double_t *xminn, Double_t *xmaxn, Double_t *distn,
                 Double_t *par1, Double_t *par2, Double_t *parn);

//maximum number of iterations in the quadratic system solver
#define MAXITER 100
//if all rates change fractionally by less than this, then we declare
//the system to be solved.
#define PREC1 1.0e-8

#define PVMORPH_MAXBINS 5000
#define CSM_DEBUGPRINT 1

/* routines needed for 2D histogram interpolation from Alex Read */

void csm_acnvec(Double_t* vec, Int_t* n);

void csm_pvmorph_2d(Int_t* nx1, Double_t* xmin1, Double_t* xmax1,
                    Int_t* ny1, Double_t* ymin1, Double_t* ymax1, 
                    Double_t* xydist1, 
                    Int_t* nx2, Double_t* xmin2, Double_t* xmax2, 
                    Int_t* ny2, Double_t* ymin2, Double_t* ymax2,
                    Double_t* xydist2, 
                    Int_t* nx3, Double_t* xmin3, Double_t* xmax3,
                    Int_t* ny3, Double_t* ymin3, Double_t* ymax3,
                    Double_t* xydist3, 
                    Double_t* par1, Double_t* par2, Double_t* par3);

void csm_ypvscat(Double_t* ydist, Double_t* xydist, Int_t* nx, Int_t* ny);

void csm_getycont(Double_t* ydist1, Int_t* ny1, 
                  Double_t* ydist2, Int_t* ny2, 
                  Double_t* ydist3, Int_t* ny3, 
	          Double_t* alpha);

class mclimit_csm
{
 public:
   mclimit_csm();
   ~mclimit_csm();
   void set_datahist(TH1 *,char *); /* data histogram and channel name */
   // the hypothesis set routines do not make clones of their inputs, they just
   // store pointers to the models. Best practice -- fully define a model (that
   // is, add all templates, before calling these routines and do not update
   // the models before the pseudoexperiments are run.

   void set_null_hypothesis(csm_model *);
   void set_test_hypothesis(csm_model *);
   void set_null_hypothesis_pe(csm_model *);
   void set_test_hypothesis_pe(csm_model *);
   void set_npe(Int_t);      // sets the number of pseudoexperiments to do.  
                             // The default is set in the constructor to 10000
   Int_t get_npe();          // returns the value set in set_npe

   void set_chisquarehistos(TH1 *,TH1 *,TH1 *,TH1 *);
   // set pointers to histograms to accumulate chisquare distributions for
   // null hyp chisquare distrib in null hyp pseudoexperiments
   // test hyp chisquare distrib in null hyp pseudoexperiments
   // null hyp chisquare distrib in test hyp pseudoexperiments
   // test hyp chisquare distrib in test hyp pseudoexperiments

   void run_pseudoexperiments(); // see set_npe() above to determine how many pseudoexperiments to run

   // output retrieval:

   Double_t cls();
   Double_t clsw();  // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clsb();
   Double_t clsbw();  // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clb();
   Double_t omclb(); // "1-clb" computed as a p-value, including the probability of the exact outcome
                     // observed in the data.  Computed with null hypothesis pseudoexperiments
   Double_t omclbw(); // Same as above, but using test hypothesis pseudoexperiments, reweighted
   // with the likelihood ratio to approximate the null hypothesis distribution

   Double_t ts();  // test statistic -- a delta chisquared -- computed for the observed data

   Double_t tsbm2(); // distributions of test statistic in null hyp pseudoexperiments 2 sigma low edge
   Double_t tsbm1(); // 1 sigma low edge
   Double_t tsbmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tsbp1();  // 1 sigma upper edge
   Double_t tsbp2();  // 2 sigma upper edge

   Double_t tssm2(); // distributions of test statistic in test hyp pseudoexperiments 2 sigma low edge
   Double_t tssm1(); // 1 sigma low edge
   Double_t tssmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tssp1();  // 1 sigma upper edge
   Double_t tssp2();  // 2 sigma upper edge
 
   Double_t clsexpbm2(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmed(); // Expected cls in null hypothesis -- median
   Double_t clsexpbp1(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsbexpbm2(); // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbm1(); // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbmed(); // Expected clsb in null hypothesis -- median
   Double_t clsbexpbp1(); // Expected clsb in null hypothesis -- 1 sigma upper edge
   Double_t clsbexpbp2(); // Expected clsb in null hypothesis -- 2 sigma upper edge


   // computed with null hypothesis px's reweighted -- good for small expected CLs's

   Double_t clsexpbm2w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmedw(); // Expected cls in null hypothesis -- median
   Double_t clsexpbp1w(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2w(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsexpsm2(); // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsm1(); // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsmed(); // Expected cls in test hypothesis -- median
   Double_t clsexpsp1(); // Expected cls in test hypothesis -- 1 sigma upper edge
   Double_t clsexpsp2(); // Expected cls in test hypothesis -- 2 sigma upper edge
 
   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)  Not to be used
   // for discovery significance!
   Double_t clbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsmed(); // Expected clb in test hypothesis -- median
   Double_t clbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t clbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome.  These are computed
   // using null hypothesis px's to compute 1-CLb's.
   Double_t omclbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmed(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpsm2w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmedw(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1w(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2w(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)
   Double_t clbexpbm2(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbm1(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbmed(); // Expected clb in null hypothesis -- median
   Double_t clbexpbp1(); // Expected clb in null hypothesis -- 1 sigma upper edge
   Double_t clbexpbp2(); // Expected clb in null hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome
   Double_t omclbexpbm2(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmed(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpbm2w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmedw(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1w(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2w(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   //these three below are computed with test hyp. px's reweighted to get the small
   //tails of the null hyp distribution better modeled with low px stats
   Double_t p2sigmat(); // Probability of a 2-sigma evidence assuming test hyp. is true
   Double_t p3sigmat(); // Probability of a 3-sigma evidence assuming test hyp. is true
   Double_t p5sigmat(); // Probability of a 5-sigma discovery assuming test hyp. is true

   //these are computed with the null hyp px's -- can be unreliable for low px stats.
   //hopefully the 1-CLb p-value is uniformly distributed between 0 and 1, but
   // these routines, as well as omclbexpb* are designed to quantify just how
   // uniform that is.
   Double_t p2sigman(); // Probability of a 2-sigma evidence assuming null hyp. is true
   Double_t p3sigman(); // Probability of a 3-sigma evidence assuming null hyp. is true
   Double_t p5sigman(); // Probability of a 5-sigma discovery assuming null hyp. is true

   Double_t calc_chi2(csm_model *, TH1 *[]);  // interface to chisquare calculator
                                                      // using our models

   Double_t weightratio(csm_model *testhyp, csm_model *nullhyp, TH1 *[]); // weight factor
   // for reweighting MC px's using the varied templates.

   // rate limite calculators
   Double_t s95();      // scale factor on signal which is excluded at exactly 95% CL
   Double_t s95m2();   // variation around the median expected s95 in the bg hypothesis, -2 sigma
   Double_t s95m1();   // variation around the median expected s95 in the bg hypothesis, -1 sigma
   Double_t s95med();  // median expected s95 in the background hypothesis
   Double_t s95p1();   // variation around the median expected s95 in the bg hypothesis, +1 sigma
   Double_t s95p2();   // variation around the median expected s95 in the bg hypothesis, +2 sigma

   Double_t lumi95(); // calculates the lumi needed for a median experiment to exclude at 95% CL
                      // what's returned is a multiplicative factor on whatever luminosity was used
                      // to construct the test and null hypotheses.
   Double_t lumi3s(); // calculates the lumi needed for a median experiment to discover at 3 sigma
   Double_t lumi5s(); // calculates the lumi needed for a median experiment to discover at 5 sigma

   void tshists(TH1*,TH1*);  // fills histograms with test statisic values in the pseudoexperiments
                             // (you define the binning).  First histo: test hypothesis, second histo:
                             // null hypothesis

   // Call Joel Heinrich's (CDF 7587) genlimit
   //Bayesian limit calculator.  First arg:  credibility level:  e.g., 0.95.  Second arg, 
   //scale factor on signal to produce the limit.  Third arg, uncertainty on the limit scale factor.
   //as with the s95 routines above, this assumes that the signal adds incoherently to the
   //background.  Requires set_test_hypothesis_pe to be called first in order to make the
   //"Prior ensemble".  The size of the prior ensemble is set with set_npe()
   //also call the relevant set_datahist methods too.

   void bayes_heinrich(Double_t beta, Double_t* sflimit, Double_t* unc);

   // Routine to call Joel Heinrich's (CDF 7587) genlimit, a Bayesian limit calculator,
   // but to repeat the calculation for pseudoexpeirments drawn from the null hypothesis
   // (be sure to call both set_test_hypothesis_pe and set_null_hypothesis_pe before using
   // this).  This computes the observed and expected limits.
   // Arguments:  1:  beta (credibility level, e.g. 0.95)
   // Argument 2:  observed limit
   // Agrument 3:  error on observed limit
   // Argument 4:  npx pseudoexperiments to run to compute expected limits
   // Arguments 5-9: Expected limits.  Median +-1, 2 sigma expectations

   void bayes_heinrich_withexpect(Double_t beta, Double_t* sflimit, Double_t* unc,
				  Int_t npx,
                                  Double_t* sm2, Double_t* sm1, Double_t* smed, Double_t* sp1,
                                  Double_t* sp2);  // Call Joel Heinrich's (CDF 7587) genlimit

   void bayes_heinrich_coverage_check(Double_t beta,
                                      Double_t* sflimit,
                                      Double_t* unc,
				      Int_t npx,
                                      Double_t* falsex);  // checks the false exclusion rate (coverage)
 private:
   csm_model *null_hypothesis;
   csm_model *test_hypothesis;
   csm_model *null_hypothesis_pe;
   csm_model *test_hypothesis_pe;

   Int_t nmc;  // number of pseudoexperiments which have been done
               // this is set to zero by the constructor and by
               // anything that modifies the hypotheses, which is the
               // indication that the pseudoexperiments need to be redone.

   Int_t nmc_req;  // number of pseudoexperiments to do

   Int_t recalctsflag; // 1 if we need to redo the data test statistic

   Double_t *tsb;  // test statistic in null hypothesis -- one per pseudoexperiment
   Double_t *tss;  // test statistic in test hypothesis -- one per pseudoexperiment
   Double_t *wtsb;  // weight for converting the corresponding tsb into a s+b px probability
   Double_t *wtss;  // weight for converting the corresponding tss into a b px probability
   Int_t *itsb; // sort order for tsb
   Int_t *itss; // sort order for tss
   Double_t tsd;                     // test statistic for data
   vector<TH1*> datahist;
   vector<char*> dhname;            // channel names for each data histogram -- must match
                                    // the channel names used in the models.

   Double_t s95aux(Int_t); // s95 calculator using a function you pass in.
   Double_t lumipaux(Int_t); // luminosity threshold calculator
   Double_t clsaux(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsauxw(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsbaux(Double_t);
   Double_t clsbauxw(Double_t);
   Double_t clbaux(Double_t);
   Double_t omclbaux(Double_t);
   Double_t omclbauxw(Double_t);

   TH1 *nullnullchisquare;
   TH1 *nulltestchisquare;
   TH1 *testnullchisquare;
   TH1 *testtestchisquare;
};

#endif

#ifndef GENLIMIT

typedef struct {
  float e,b;
}EB;

typedef enum {
  flat=10,
  corr=20
} PRIOR;

double cspdf(double s,double norm,
             int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior);

double csint(double s0,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty);

void csint2(double s1,double s2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22);

void csint2cut(double s1,double s2,double shi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22);

double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
               int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty);

double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty);

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens);

void gausslaguerre(double x[],double lw[],int n,double alpha);

#define GENLIMIT 1

#endif
