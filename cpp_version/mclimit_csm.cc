/* This is a set of routines that works like mclimit.f but uses csm.c to compute
   the chisquared of the data histogram compared to a sum of models, each of which
   may (or may not) have sensitivities to nuisance parameters in their shapes,
   normalizations, and even statistical uncertainty in each bin */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include<assert.h>
#include <stddef.h>
#include <algorithm>
#include "TRandom.h"
#include "TMinuit.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "mclimit_csm.h"
#include "TString.h"

using namespace std;

#define max(max_a,max_b) (((max_a)>(max_b)) ? (max_a) : (max_b))
#define min(min_a,min_b) (((min_a)>(min_b)) ? (min_b) : (min_a))

// Minuit ugliness -- data communication with the function to fit is either
// via member functions of a new class inherited from TObject (not done here),
// or in global data (at least global to this source file) storage, which is
// the method chosen here because it's easier.

static vector<TH1*> datatofit;
static vector<char*> datatofitname;
static csm_model *modeltofit;
static vector<Int_t> constrainedfitparam;
static vector<char*> npfitname;

void csm_minuit_fcn(Int_t &npar, double *gin, double &f,
                        double *par, Int_t iflag);

double csint0(double xlo,double logscale,
		     int nchan,const int nobs[],const EB chan[],
		     int ngl,const double xgl[],const double lwgl[],
		     PRIOR prior);

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
	          int nchan,const int nobs[],const EB chan[],
		  int ngl,const double xgl[],const double lwgl[],PRIOR prior,
	          double* int1,double* int2);

void gameansigma(double *mean,double *sigma,
			int nchan,int nens,const int nobs[],const EB* ens);
double arcfreq(double y);
#define freq(x) (0.5*erfc(-0.707106781186547524*(x)))

void csint02(double xlo1,double xlo2,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2);

/*----------------------------------------------------------------------------*/

// constructor

mclimit_csm::mclimit_csm()
{
  nmc = 0;
  nmc_req = 10000;
  recalctsflag = 1;
  // set null pointers to our cumulative histograms -- if we
  // don't get any from the user, don't bother filling them.
  nullnullchisquare = 0;
  nulltestchisquare = 0;
  testnullchisquare = 0;
  testtestchisquare = 0;
  // null pointers to the test statistic arrays -- need to allocate memory when
  // we know how many to do
  tss = 0;
  tsb = 0;
  wtss = 0;
  wtsb = 0;
  itss = 0;
  itsb = 0;
}

/*----------------------------------------------------------------------------*/

// destructor

mclimit_csm::~mclimit_csm()
{
  Int_t i;

  // deallocate cloned input data histograms and their names.

  for (i=0; i<(Int_t) datahist.size(); i++)
    {
      delete datahist[i];
      delete[] dhname[i]; 
    }
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel data histograms and channel names that is sorted by channel */
/* name during the building process.  Use the same sorting procedure as used in */
/* csm_model::lookup_add_channame so that our data histograms are stored in */
/* the same order as our models */

void mclimit_csm::set_datahist(TH1 *h, char *cname)
{
  Int_t i,ifound,j,jfound;
  char *s;
  vector<char*>::iterator nhi;
  vector<TH1*>::iterator dhi;

  recalctsflag = 1;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) dhname.size(); i++)
    {
      j = (Int_t) strcmp(cname,dhname[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings.  If the name is on the 
     list, replace the existing data histogram with a clone of the one supplied. */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      strcpy(s,cname);
      if (jfound == -1)
	{
          dhname.push_back(s);
          datahist.push_back((TH1*) h->Clone());
	}
      else
	{
	  nhi = dhname.begin() + jfound;
	  dhname.insert(nhi,s);
	  dhi = datahist.begin() + jfound;
	  datahist.insert(dhi,(TH1*) h->Clone());
	}
    }
  else
    {
      delete datahist[ifound];
      datahist[ifound] = (TH1*) h->Clone();
    }
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_npe(Int_t nperequest)
{
  if (nperequest < 0)
    {
      cout << "mclimit_csm::set_npe: Invalid pseudoexperiment request: " << nperequest << endl;
      exit(0);
    }
  nmc_req = nperequest;
}

/*----------------------------------------------------------------------------*/

Int_t mclimit_csm::get_npe()
{
  return(nmc_req);
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis(csm_model *model)
{
  null_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis(csm_model *model)
{
  test_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis_pe(csm_model *model)
{
  null_hypothesis_pe = model;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis_pe(csm_model *model)
{
  test_hypothesis_pe = model;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

// test statistic of observed data

Double_t mclimit_csm::ts()
{
  Int_t i;

  if (recalctsflag) 
    {
      // copy the data histogram pointers into a flat array
      // for the chisquare calculator.
      TH1** darray = new TH1*[datahist.size()];
      for (i=0;i<(Int_t)datahist.size();i++)
	{
	  darray[i] = datahist[i];
	}
      tsd = calc_chi2(test_hypothesis,darray) -
	calc_chi2(null_hypothesis,darray);
      recalctsflag = 0;
      delete[] darray;
    }
  return(tsd);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/
// confidence levels for an arbitrary test statistic.  Used for computing
// actual and expected confidence levels.
// may need to work on this a bit because of finite MC statistics -- it's hard to
// compute the confidence levels on the tails unless the histograms are properly
// filled out there.  Particularly values of clb very close to 1 need to be computed
// with extreme care.
// This is addressed with the reweighting procedure suggested by Alex Read --
// the likelihood ratio can be used to reweight test hypothesis pseudoexperiments
// to model the null hypotheis background distribution, and vice versa.

Double_t mclimit_csm::clsbaux(Double_t tsaux)
{
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clsbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tss[itss[i]] < tsaux)
	{ 
	  clsbloc = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(1-clsbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clsbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] >= tsaux)
	{ 
          clsbloc += wtsb[itsb[i]];
	}
    }
  clsbloc /= ((Double_t) nmc);
  return(clsbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbaux(Double_t tsaux)
{
  Int_t i;
  Double_t clbloc;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] < tsaux)
	{ 
	  clbloc = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(1-clbloc);
}

/*----------------------------------------------------------------------------*/

// compute 1-CLb (as a p-value, including the outcomes with exactly the -2lnQ
// observed) using null hypothesis px's.

Double_t mclimit_csm::omclbaux(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 1;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] <= tsaux)
	{ 
	  omclbloc = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/

/*  Compute 1-CLb using reweighted test hypothesis pseudoexperiments.  Reweight
    using the inverese of the likelihood ratio, 1/(p(data|test)/p(data|null)) */

Double_t mclimit_csm::omclbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc = 0;

  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tss[itss[i]] <= tsaux)
	{ 
          omclbloc += wtss[itss[i]];
	}
    }
  omclbloc /= ((Double_t) nmc);
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsaux(Double_t tsaux)
{
  Double_t clbloc,clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0)
    {
      clsloc = clsbaux(tsaux)/clbloc;
    }
  else
    {
      clsloc = 1;
    }
  return(clsloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsauxw(Double_t tsaux)
{
  Double_t clbloc,clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0)
    {
      clsloc = clsbauxw(tsaux)/clbloc;
    }
  else
    {
      clsloc = 1;
    }
  return(clsloc);
}

/*----------------------------------------------------------------------------*/
// confidence levels using the data test statistic.  Recompute the data test
// statistic if need be.

Double_t mclimit_csm::cls()
{
  return(clsaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsb()
{
  return(clsbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsw()
{
  return(clsauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbw()
{
  return(clsbauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clb()
{
  return(clbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclb()
{
  return(omclbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbw()
{
  return(omclbauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2()
{
  return(clsaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1()
{
  return(clsaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmed()
{
  return(clsaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1()
{
  return(clsaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2()
{
  return(clsaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2w()
{
  return(clsauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1w()
{
  return(clsauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmedw()
{
  return(clsauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1w()
{
  return(clsauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2w()
{
  return(clsauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm2()
{
  return(clsbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm1()
{
  return(clsbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbmed()
{
  return(clsbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp1()
{
  return(clsbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp2()
{
  return(clsbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm2()
{
  return(clsaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm1()
{
  return(clsaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsmed()
{
  return(clsaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp1()
{
  return(clsaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp2()
{
  return(clsaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm2()
{
  return(clbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm1()
{
  return(clbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbmed()
{
  return(clbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp1()
{
  return(clbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp2()
{
  return(clbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm2()
{
  return(clbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm1()
{
  return(clbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsmed()
{
  return(clbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp1()
{
  return(clbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp2()
{
  return(clbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2()
{
  return(omclbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1()
{
  return(omclbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmed()
{
  return(omclbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1()
{
  return(omclbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2()
{
  return(omclbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2()
{
  return(omclbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1()
{
  return(omclbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmed()
{
  return(omclbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1()
{
  return(omclbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2()
{
  return(omclbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2w()
{
  return(omclbauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1w()
{
  return(omclbauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmedw()
{
  return(omclbauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1w()
{
  return(omclbauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2w()
{
  return(omclbauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2w()
{
  return(omclbauxw(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1w()
{
  return(omclbauxw(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmedw()
{
  return(omclbauxw(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1w()
{
  return(omclbauxw(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2w()
{
  return(omclbauxw(tssp2()));
}



/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigmat()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigmat()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigmat()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigman()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigman()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigman()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95()
{
  return(s95aux(MCLIMIT_CSM_CLS));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m2()
{
  return(s95aux(MCLIMIT_CSM_CLSM2));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m1()
{
  return(s95aux(MCLIMIT_CSM_CLSM1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95med()
{
  return(s95aux(MCLIMIT_CSM_CLSMED));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p1()
{
  return(s95aux(MCLIMIT_CSM_CLSP1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p2()
{
  return(s95aux(MCLIMIT_CSM_CLSP2));
}

Double_t mclimit_csm::s95aux(Int_t itype)
{
  Double_t sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t j,foundit;
  csm_model *testhypsave,*testhyppesave;

  sfl = 0;
  sfh = 0;
  cltest = 0;

  // this hypothesis never gets destroyed, but we lose the pointer to it

  testhypsave = test_hypothesis;
  testhyppesave = test_hypothesis_pe;

  sf = 1.0;
  cll = 0.0;
  cla = 0.0;
  clh = 0.0; 
  foundit = 0;

  for (j=0;j<32;j++)
  {
    //cout << "in s95aux seek " << j << " scale: " << sf << endl;
    if (foundit == 0)
      {
	csm_model* scaledsignal = testhypsave->scalesignal(sf);
        csm_model* scaledsignalpe = testhyppesave->scalesignal(sf);
        test_hypothesis = scaledsignal;
        test_hypothesis_pe = scaledsignalpe;
	run_pseudoexperiments();
        recalctsflag = 1;
        if (itype == MCLIMIT_CSM_CLS)
	  {
          cltest = cls(); 
	  }
        else if (itype == MCLIMIT_CSM_CLSM2)
	  {
	    cltest = clsexpbm2();
	  }
        else if (itype == MCLIMIT_CSM_CLSM1)
	  {
	    cltest = clsexpbm1();
	  }
        else if (itype == MCLIMIT_CSM_CLSMED)
	  {
	    cltest = clsexpbmed();
	  }
        else if (itype == MCLIMIT_CSM_CLSP1)
	  {
	    cltest = clsexpbp1();
	  }
        else if (itype == MCLIMIT_CSM_CLSP2)
	  {
	    cltest = clsexpbp2();
	  }
        delete scaledsignal;
        delete scaledsignalpe;
	if (j==0)
	  {
	    cla = cltest;
	  }
	if (cltest<0.05)
	  {
	    if (cla>0.05)
	      {
		sfh = sf;
		clh = cltest;
		sfl = sf/2.0;
		cll = cla;
		foundit = 1;
	      }
	    sf /= 2.0;
	  }
	else if (cltest>0.05)
	  {
	    if (cla<0.05)
	      {
		sfl = sf;
		cll = cltest;
		sfh = sf*2.0;
		clh = cla;
		foundit = 1;
	      }
	    sf *= 2.0;
	  }
	else
	  {
	    test_hypothesis = testhypsave;
	    test_hypothesis_pe = testhyppesave;
	    return(sf);
	  }
	cla = cltest;
      }
  } // end of loop over 32 powers of 2 in search of a signal scale factor which brackets
  // 95% CL exclusion

  //cout << "done with seek loop " << sf << endl;
  sf = sfh;
  if (foundit == 0)
    { 
      cout << "mclimit_csm::s95** could not find s95 within 2**32 of original guess" << endl;
      sf = 0;
    }
  else
    {
      // From Tom Wright -- speed up by doing a deterministic five more
      // calcs of CL and a linear fit of log(CL) vs. sf.

      // find error on 0.05 CL for number of PEs
      double dcl=sqrt(0.05*0.95/nmc);

      // put in some protection against logarithms of negative numbers
      // makes sure -5*dcl + 0.05 is not negative.

      dcl = min(dcl,0.0099);

      // try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
      // increment stuff used for linear fit of ln(CL) vs sf
      double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
      for( int j=-5; j<6; j+=2 )
	{
	  sf = sfl + (log(0.05+j*dcl) - log(cll))*
	    (sfl-sfh)/(log(cll)-log(clh));

	  double clsbtest=0;
	  double clbtest=0;

	  // calculate CL for this sf
	  csm_model* scaledsignal = testhypsave->scalesignal(sf);
          csm_model* scaledsignalpe = testhyppesave->scalesignal(sf);
          test_hypothesis = scaledsignal;
          test_hypothesis_pe = scaledsignalpe;
          run_pseudoexperiments();
          recalctsflag = 1;
          if (itype == MCLIMIT_CSM_CLS)
	    {
	      cltest = cls();
	      clbtest = clb();
	      clsbtest = clsb();
  	    }
          else if (itype == MCLIMIT_CSM_CLSM2)
	    {
	      cltest = clsexpbm2();
	      clbtest = clbexpbm2();
	      clsbtest = clsbexpbm2();
	    }
          else if (itype == MCLIMIT_CSM_CLSM1)
	    {
	      cltest = clsexpbm1();
	      clbtest = clbexpbm1();
	      clsbtest = clsbexpbm1();
	    }
          else if (itype == MCLIMIT_CSM_CLSMED)
	    {
	      cltest = clsexpbmed();
	      clbtest = clbexpbmed();
	      clsbtest = clsbexpbmed();
	    }
          else if (itype == MCLIMIT_CSM_CLSP1)
	    {
	      cltest = clsexpbp1();
	      clbtest = clbexpbp1();
	      clsbtest = clsbexpbp1();
	    }
          else if (itype == MCLIMIT_CSM_CLSP2)
	    {
	      cltest = clsexpbp2();
	      clbtest = clbexpbp2();
	      clsbtest = clsbexpbp2();
	    }
	  delete scaledsignal;
	  delete scaledsignalpe;

	  // double dcltest=sqrt(cltest*(1-cltest)/nmc);
// 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
// 				       (1-clsbtest)/clsbtest/nmc);
	  double dcltest=sqrt(clsbtest*(1-clsbtest)/nmc)/clbtest;

          //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

	  double lcl = log(cltest);
	  double dlcl = dcltest/cltest;

	  lf_a += sf/dlcl/dlcl;
	  lf_b += 1/dlcl/dlcl;
	  lf_c += lcl/dlcl/dlcl;
	  lf_d += sf*sf/dlcl/dlcl;
	  lf_e += sf*lcl/dlcl/dlcl;
	  lf_f += lcl*lcl/dlcl/dlcl;
	}

      // Find fit parameters for log(CL)=p1+p2*sf
      double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
      double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

      //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
      //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
      //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;

      //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

      //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
      //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
      //printf("chisuare/dof: %f\n",lf_x2/4);

      // invert to get sf at 0.05 and its error
      // assuming 100% anticorrelation
      double lf_sf = (log(0.05)-lf_p1)/lf_p2;
      //printf("CL variation at %f: %f %f %f\n",lf_sf,
      //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
      //     exp(lf_p1+lf_p2*lf_sf),
      //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));

      //double lf_dsf1 = (log(0.05)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
      //double lf_dsf2 = (log(0.05)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);

      //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
      sf = lf_sf;

    }

  test_hypothesis = testhypsave;
  test_hypothesis_pe = testhyppesave;
  recalctsflag = 1;
  return(sf);
}


/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi95()
{
  return(lumipaux(MCLIMIT_CSM_LUMI95));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi3s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI3S));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi5s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI5S));
}

/*----------------------------------------------------------------------------*/
// compute median amounts of luminosity needed for 95% CL exclusion, 3 sigma
// evidence, or 5 sigma discovery -- scale the systematic errors with 1/sqrt(lumi/lumi_0)

Double_t mclimit_csm::lumipaux(Int_t itype)
{
  Double_t sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t foundit,j;
  csm_model *testhypsave,*nullhypsave;
  csm_model *testhyppesave,*nullhyppesave;
  Double_t resdes;

  resdes = 0.5; // do median luminosity thresholds

  sfl = 0;
  sfh = 0;
  cltest = 0;

  sf = 1.0;
  cll = 0.0;
  cla = 0.0;
  clh = 0.0; 
  foundit = 0;

  testhypsave = test_hypothesis;
  nullhypsave = null_hypothesis;
  testhyppesave = test_hypothesis_pe;
  nullhyppesave = null_hypothesis_pe;
  
  for (j=0;j<32;j++)
    {
      if (foundit == 0)
	{
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }
	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  if (j==0)
	    {
	      cla = cltest;
	    }
	  if (cltest < resdes)
	    {
	      if (cla > resdes)
	        {
		  sfl = sf;
	  	  cll = cltest;
		  sfh = sf*2.0;
		  clh = cla;
		  foundit = 1;
	        }
	      sf *= 2.0;
	    }
  	  else if (cltest > resdes)
	    {
	      if (cla < resdes)
	        {
 		  sfh = sf;
		  clh = cltest;
		  sfl = sf/2.0;
		  cll = cla;
		  foundit = 1;
	        }
	      sf /= 2.0;
	    }
	  else
	    {
	      test_hypothesis = testhypsave;
	      null_hypothesis = nullhypsave;
	      test_hypothesis_pe = testhyppesave;
	      null_hypothesis_pe = nullhyppesave;
	      return(sf);
	    }
	  cla = cltest;
        }
    } // end of loop over 32 powers of 2 in search of a luminosity scale factor which
      // brackets the desired sensitvity

  sf = sfl;
  if (foundit == 0)
    { 
      cout << "mclimit_csm::lumipaux** could not find s95 within 2**32 of original guess" << endl;
      sf = 0;
    }
  else
    {

      // From Tom Wright -- speed up by doing a deterministic five more
      // calcs of CL and a linear fit of log(CL) vs. sf.

      // find error on resdes CL for number of PEs
      double dcl=sqrt(resdes*(1.0-resdes)/nmc);

      // put in some protection against logarithms of negative numbers
      // makes sure -5*dcl + resdes is not negative.

      dcl = min(dcl,resdes/5 - 0.0001);

      // try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
      // increment stuff used for linear fit of ln(CL) vs sf
      double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
      for( int j=-5; j<6; j+=2 )
	{
	  sf = sfl + (log(resdes+j*dcl) - log(cll))*
	    (sfl-sfh)/(log(cll)-log(clh));

	  //double clsbtest, clbtest;

	  // calculate CL for this sf
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }

	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  // double dcltest=sqrt(cltest*(1-cltest)/nmc);
// 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
// 				       (1-clsbtest)/clsbtest/nmc);
	  double dcltest=sqrt(cltest*(1-cltest)/nmc)/cltest;

          //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

	  double lcl = log(cltest);
	  double dlcl = dcltest/cltest;

	  lf_a += sf/dlcl/dlcl;
	  lf_b += 1/dlcl/dlcl;
	  lf_c += lcl/dlcl/dlcl;
	  lf_d += sf*sf/dlcl/dlcl;
	  lf_e += sf*lcl/dlcl/dlcl;
	  lf_f += lcl*lcl/dlcl/dlcl;
	}

      // Find fit parameters for log(CL)=p1+p2*sf
      double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
      double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

      //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
      //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
      //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;

      //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

      //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
      //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
      //printf("chisuare/dof: %f\n",lf_x2/4);

      // invert to get sf at resdes and its error
      // assuming 100% correlation
      double lf_sf = (log(resdes)-lf_p1)/lf_p2;
      //printf("CL variation at %f: %f %f %f\n",lf_sf,
      //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
      //     exp(lf_p1+lf_p2*lf_sf),
      //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));

      //double lf_dsf1 = (log(resdes)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
      //double lf_dsf2 = (log(resdes)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);

      //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
      sf = lf_sf;

    }
   
  test_hypothesis = testhypsave;
  null_hypothesis = nullhypsave;
  test_hypothesis_pe = testhyppesave;
  null_hypothesis_pe = nullhyppesave;
  return(sf);
}

void mclimit_csm::tshists(TH1* testhypts, TH1* nullhypts)
{
  int i;
  testhypts->Reset();
  nullhypts->Reset();
  for (i=0;i<nmc;i++)
    { 
      testhypts->Fill(tss[i]);
      nullhypts->Fill(tsb[i]);

    }
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_chisquarehistos(TH1 *nn, TH1 *nt, TH1 *tn, TH1 *tt)
{
  nullnullchisquare = nn;
  nulltestchisquare = nt;
  testnullchisquare = tn;
  testtestchisquare = tt;

}

/*----------------------------------------------------------------------------*/
// run pseudoexperiments, allocating and filling the tss and tsb arrays
// also fill in wtss and wtsb, for use in reweighting pseudoexperiments, using
// the varied nuisance parameters.  Sort tss and tsb, and keep the sort order
// arrays itss and itsb around so the corresponding wtss and wtsb arrays can
// be used with them

void mclimit_csm::run_pseudoexperiments()
{
  Int_t i;
  char *pdname;
  // Double_t tmp;
  Double_t csnull,cstest;

  // make some histograms to store the pseudodata.

  TH1** pdarray = new TH1*[null_hypothesis_pe->channame.size()];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      pdname = new char[strlen(null_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) null_hypothesis_pe->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete pdname;
    }

  // allocate memory for test statistic and weight and sort order storage
  if (tss != 0) delete[] tss;
  if (tsb != 0) delete[] tsb;
  if (wtss != 0) delete[] wtss;
  if (wtsb != 0) delete[] wtsb;
  if (itss != 0) delete[] itss;
  if (itsb != 0) delete[] itsb;
  tss = new Double_t[nmc_req];
  tsb = new Double_t[nmc_req];
  wtss = new Double_t[nmc_req];
  wtsb = new Double_t[nmc_req];
  itss = new Int_t[nmc_req];
  itsb = new Int_t[nmc_req];

  for(i=0;i<nmc_req;i++)
    {
      //      cout << "At pseudoexperiment " << i << " out of " << nmc_req << endl;
      null_hypothesis_pe->single_pseudoexperiment(pdarray);
      wtsb[i] = weightratio(test_hypothesis_pe,null_hypothesis_pe,pdarray);
      csnull = calc_chi2(null_hypothesis,pdarray);
      cstest = calc_chi2(test_hypothesis,pdarray);
      if (nullnullchisquare != 0)
	{ nullnullchisquare->Fill(csnull); }
      if (nulltestchisquare != 0)
	{ nulltestchisquare->Fill(cstest); }
      tsb[i] = cstest - csnull;
      // cout << "null hyp chisquared: " << csnull << " test hyp chisquared: " << cstest << endl;
      // cout << cstest-csnull << endl;

      test_hypothesis_pe->single_pseudoexperiment(pdarray);
      wtss[i] = weightratio(null_hypothesis_pe,test_hypothesis_pe,pdarray);
      csnull = calc_chi2(null_hypothesis,pdarray);
      cstest = calc_chi2(test_hypothesis,pdarray);
      if (testnullchisquare != 0)
	{ testnullchisquare->Fill(csnull); }
      if (testtestchisquare != 0)
	{ testtestchisquare->Fill(cstest); }
      tss[i] = cstest - csnull;
      // cout << "null hyp chisquared: " << csnull << " test hyp chisquared: " << cstest << endl;
      // cout << " " << cstest-csnull << endl;
    }

  TMath::Sort(nmc_req,tss,itss,0);
  TMath::Sort(nmc_req,tsb,itsb,0);

  for (i=0;i<(Int_t)null_hypothesis_pe->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  nmc = nmc_req;
}

/*----------------------------------------------------------------------------*/
// update the internal representation of the nuisance response of the model.
// Use the method for doing this with each
// channel separately.

void csm_model::nuisance_response(Int_t nparams,
                                  char *paramname[],
                                  Double_t paramvalue[])
{
  Int_t i,nchans;
  Int_t ipar,icons,j,k,ifound,jfound;
  Double_t *cinput;

  /*
  cout << "in model nuisance response: " << endl;
  for (i=0;i < (Int_t) nparams; i++)
  {
    cout << "param: " << i << " name: " << paramname[i] << endl;
  }
  */

  // compute the nuisance parameters which are functions of the others
  // loop over constraints and replace the parameter values with the constrained ones
  // do not overwrite the input parameters but make our own copy

  Double_t *parloc = new Double_t[nparams];
  for (i=0;i<nparams;i++)
    {
      parloc[i] = paramvalue[i];
    }
  for (icons=0;icons<(Int_t)npcm.size();icons++)
    {
      jfound = 0;
      for (ipar=0;ipar<nparams;ipar++)
	{
          if (strcmp(npcm[icons].pnameoutput,paramname[ipar])==0)
	    {
	      jfound = 1;
	      cinput = new Double_t[npcm[icons].ninput];
	      for (j=0;j<(Int_t)npcm[icons].ninput;j++)
	        {
	          ifound = 0;
	          for (k=0;k<nparams;k++)
		    {
		      if (strcmp(npcm[icons].pnameinput[j],paramname[k])==0)
		        {
		          cinput[j] = parloc[k];
		          ifound = 1;
		        }
		    }
	          if (ifound == 0)
		    {
		      cout << "Didn't find parameter name: " << 
                          npcm[icons].pnameinput[j] << 
                          " in the list of nuisance parameters" << endl;
		      exit(0);
		    }
	        }
	      parloc[ipar] = npcm[icons].f(cinput);
	      delete[] cinput;
	    }
	}
      if (jfound == 0) 
	{
	  cout << "Constraint equation found for nuisance parameter not on our list: " 
               << npcm[icons].pnameoutput << endl;
          exit(0);
	}
    }

  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->nuisance_response(nparams,paramname,parloc);
    }
  delete[] parloc;
}

/*----------------------------------------------------------------------------*/

void csm_model::undo_nuisance_response()
{
  Int_t i,nchans;
  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->undo_nuisance_response();
    }
}

/*----------------------------------------------------------------------------*/

// Create a fluctuated channel model -- input a list of nuisance parameter names
// and values, and return a pointer to a new channel model which has 
// responded to those nuisance parameters.  Be sure to delete it when done.
// todo -- make sure that the shape errors accumulate.  Suggestion of John
// Zhou: average all shape error interpolations.

void csm_channel_model::nuisance_response(Int_t nparams,
                                          char *paramname[],
                                          Double_t paramvalue[])
{
  Int_t i,j,itpl,nsys,ntemplates;

  /*
  cout << "in channel nuisance response: " << endl;
  for (i=0;i < (Int_t) syserr.size(); i++)
  {
    cout << "error source: " << i << " name: " << syserr[i].sysname << endl;
  }
  for (i=0;i < (Int_t) nparams; i++)
  {
    cout << "param: " << i << " name: " << paramname[i] << endl;
  }
  */

  undo_nuisance_response();
  ntemplates = (Int_t) histotemplate.size();
  Int_t *avgcount = new Int_t[ntemplates];
  for (itpl=0;itpl<ntemplates;itpl++)
    {
      avgcount[itpl] = 0;
    }

  TH1* hcl = 0;

  nsys = (Int_t) syserr.size();
  for (i=0;i < nsys; i++)
    {
      for (j=0; j<nparams; j++)
	{
	  if (strcmp(syserr[i].sysname,paramname[j]) == 0)
	    {
	      itpl = syserr[i].itemplate;
	      sft_varied[itpl] *= max(1E-8,( 
                  (syserr[i].sysfrach+syserr[i].sysfracl)*paramvalue[j]*paramvalue[j]/2.0 +
                  (syserr[i].sysfrach-syserr[i].sysfracl)*paramvalue[j]/2.0 + 1.0));
	      if (paramvalue[j]>0)
		{
		  if (syserr[i].highshape != 0)
		    {
		      if (hcl == 0)
			{
                           hcl = (TH1*) histotemplate[itpl]->Clone();
			}
		      csm_interpolate_histogram(histotemplate[itpl],0.0,
                                                syserr[i].highshape,syserr[i].xsighigh,
                                                hcl,paramvalue[j],chan_istyle);
                      histotemplate_varied[itpl]->Add(histotemplate_varied[itpl],hcl,
						      ((Double_t) avgcount[itpl])/((Double_t) (avgcount[itpl]+1)),
						      1/((Double_t) (avgcount[itpl]+1)));
		      avgcount[itpl]++;
		      //cout << "did a +interpolation " << i << " " << j << " param: " << paramvalue[j] <<  " " << avgcount[itpl] << endl;
		    }
		}
	      else
		{
		  if (syserr[i].lowshape != 0)
		    {
		      if (hcl == 0)
			{
                           hcl = (TH1*) histotemplate[itpl]->Clone();
			}
		      csm_interpolate_histogram(histotemplate[itpl],0.0,
                                                syserr[i].lowshape, -syserr[i].xsiglow,
                                                hcl, -paramvalue[j],chan_istyle);
                      histotemplate_varied[itpl]->Add(histotemplate_varied[itpl],hcl,
						      ((Double_t) avgcount[itpl])/((Double_t) (avgcount[itpl]+1)),
						      1/((Double_t) (avgcount[itpl]+1)));
		      avgcount[itpl]++;
		      //cout << "did a -interpolation " << i << " " << j << " param: " << paramvalue[j] <<  " " << avgcount[itpl] << endl;
		    }
		}
	    }
	}
    }
  if (hcl != 0)
    {
      delete hcl;
    }
  delete[] avgcount;
}

/*----------------------------------------------------------------------------*/
// resets all the varied histograms and scales to their unvaried states.

void csm_channel_model::undo_nuisance_response()
{
  Int_t i,ntemplates,nbinsx,nbinsy,ix,iy;

  ntemplates = (Int_t) histotemplate.size();
  for (i=0;i<ntemplates;i++)
    {
      sft_varied[i] = sft[i];
      nbinsx = histotemplate[i]->GetNbinsX();
      nbinsy = histotemplate[i]->GetNbinsY();
      if (nbinsy==1)
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      histotemplate_varied[i]->SetBinContent(ix,histotemplate[i]->GetBinContent(ix));
	    }
	}
      else
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      for (iy=1;iy<=nbinsy;iy++)
		{
	          histotemplate_varied[i]->SetBinContent(ix,iy,histotemplate[i]->GetBinContent(ix,iy));
		}
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
// check to see if any bin has a total negative prediction in this channel

Int_t csm_channel_model::checkneg()
{
  cout << "csm_channel_model::checkneg() to be written" << endl;
  return(0);
}

/*----------------------------------------------------------------------------*/
// collect all nuisance parameter names and upper and lower bounds for this model
// Where the bounds come from -- do not allow extrapolation on histogram shapes.  
// (all histogram extrapolation should be done and verified by the user)
// Also do not allow any contribution to signal or background to go negative.

void csm_model::list_nparams(vector<char *> *npn, vector<Double_t> *nplb, vector<Double_t> *nphb)
{
  Int_t i,j,k,ifound;
  csm_channel_model* cm;
  Double_t nplb_tmp,nphb_tmp,a,b,c,disc,xp,xm,xht,xlt;
  npn->clear();
  nplb->clear();
  nphb->clear();

  for (i=0;i<(Int_t) channame.size();i++)
    {
      cm = chanmodel[i];
      for (j=0;j<(Int_t) cm->syserr.size();j++)
	{
	  //cout << "sys error item channel: " << i << 
          //" error index: " << j << " " << cm->syserr[j].sysname << " " <<
	  //(cm->syserr[j].highshape != 0) << " " <<  
          //(cm->syserr[j].lowshape != 0) << " " <<
          //cm->syserr[j].xsiglow << " " << cm->syserr[j].xsighigh << endl;  

          // question -- do we need to consider nuisance parameter variations beyond 20 sigma?
          // probably not if we only need 5-sigma discovery significance.

          nplb_tmp = -20;
	  nphb_tmp = 20;

          // Require the user to supply shape variations out to the number of sigma
          // we will investigate here.  This program won't do shape extrapolations internally,
          // but the csm_pvmorph subroutine supplied will in fact extrapolate.  Users should
          // look at what they get when extrapolating histograms, though -- check and validate.

          if (cm->syserr[j].lowshape != 0)
	    { nplb_tmp = max(nplb_tmp,cm->syserr[j].xsiglow); }
          if (cm->syserr[j].highshape != 0)
	    { nphb_tmp = min(nphb_tmp,cm->syserr[j].xsighigh); }

          // limit the nuisance paramters also so that individual scale factors do not go negative.
          // There's protection in the fit function, but we need the pseudoexperiments also to
          // be sensible -- this is the equivalent (using the asymmetric errors supplied) of the
          // truncated Gaussian

	  a = (cm->syserr[j].sysfrach + cm->syserr[j].sysfracl)/2.0;
	  b = (cm->syserr[j].sysfrach - cm->syserr[j].sysfracl)/2.0;
          c = 1;
	  if (a == 0)
	    {
	      if (b > 0)
		{ nplb_tmp = max(nplb_tmp,-1.0/b); }
	      if (b < 0)
		{ nphb_tmp = min(nphb_tmp,-1.0/b); }
	    }
	  else
	    {
	      disc = b*b - 4.0*a*c;
	      if (disc > 0)
		{ 
		  xp = (-b + sqrt(disc))/(2.0*a);
		  xm = (-b - sqrt(disc))/(2.0*a);
		  xht = max(xp,xm);
		  xlt = min(xp,xm); 
		  nphb_tmp = min(nphb_tmp,xht);
		  nplb_tmp = max(nplb_tmp,xlt);
		}
	    }

	  ifound = -1;
	  for (k=0;k<(Int_t) npn->size();k++)
	    {
	      if (strcmp(cm->syserr[j].sysname,(*npn)[k]) == 0) { ifound = k; }
	    }
	  if (ifound == -1)
	    {
	      npn->push_back(cm->syserr[j].sysname);
	      nplb->push_back(nplb_tmp);
	      nphb->push_back(nphb_tmp);
	      //cout << "sysname: " << cm->syserr[j].sysname << " assigned ranges: " << nplb_tmp << " " << nphb_tmp << endl;
	    }
	  else
	    {
	      (*nplb)[ifound] = max((*nplb)[ifound],nplb_tmp);
	      (*nphb)[ifound] = min((*nphb)[ifound],nphb_tmp); 
	      //cout << "sysname: " << cm->syserr[j].sysname << " reassigned ranges: " << nplb_tmp << " " << nphb_tmp <<  " " <<
	      // (*nplb)[ifound] << " " << (*nphb)[ifound] << endl;
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
/* a splitoff from single_pseudoexperiment -- just vary the templates but do  */
/* not generate pseudodata.   Useful for interfacing with Joel's program      */
/*----------------------------------------------------------------------------*/

void csm_model::varysyst()
{
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Int_t i;

 // systematically fluctuate our model

  list_nparams(&npn, &nplb, &nphb);
  //cout << " in pe: npn.size " << npn.size() << endl;
  Double_t *npval = new Double_t[npn.size()];
  char **npnloc = new char *[npn.size()];

  for (i=0;i<(Int_t)npn.size();i++)
    {
      npnloc[i] = npn[i];
      do 
	{
          npval[i] = gRandom->Gaus(0,1);
          //cout << i << " " << nplb[i] << " " << npval[i] << " " << nphb[i] << endl;
	}
      while (npval[i] < nplb[i] || npval[i] > nphb[i]);
    }
  nuisance_response(npn.size(),npnloc,npval);
  delete[] npval;
  delete[] npnloc;
}

/*----------------------------------------------------------------------------*/

/* Generate a single pseudoexperiment from a model -- fluctuate all nuisance parameters
   with their uncertainties -- the pseudodata histograms are in the same order as 
   the channels in the model description, with the same binning assumed.  The psuedodata
   histograms should be allocated in the calling routine.  That way the histograms don't
   have to be continually created and destroyed for each pseudoexperiment but can be
   re-used.*/

void csm_model::single_pseudoexperiment(TH1 *pseudodata[])
{
  Int_t ichan,itpl,ibinx,ibiny,nbinsx,nbinsy,nchans,ntemplates;
  csm_channel_model* cm;
  Double_t bintot;
  TH1* ht;
  Double_t r,gbc;

  // call nuisance_response with random nuisance parameters

  varysyst();
 
  // generate random pseudodata.  Randomly fluctuate the Poisson subsidiary
  // experiments (a "systematic effect") to figure out what the proper mean
  // is for the main experiment.

  nchans = (Int_t) channame.size();
  for (ichan=0;ichan<nchans;ichan++)
    {
      cm = chanmodel[ichan];
      nbinsx = cm->histotemplate[0]->GetNbinsX();
      nbinsy = cm->histotemplate[0]->GetNbinsY(); 
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      bintot = 0;
	      ntemplates = (Int_t) cm->histotemplate.size();
              for (itpl=0;itpl<ntemplates;itpl++)
	        {
	          ht = cm->histotemplate_varied[itpl];
	          if (cm->poissflag[itpl])
		    {
		      if (nbinsy == 1) 
			{ gbc = ht->GetBinContent(ibinx+1); }
		      else
			{ gbc = ht->GetBinContent(ibinx+1,ibiny+1); }
		      r = (Double_t) gRandom->Poisson(gbc);
		    }
	          else
	            {
		      if (nbinsy == 1) 
			{ r = ht->GetBinContent(ibinx+1); }
		      else
			{ r = ht->GetBinContent(ibinx+1,ibiny+1); }
		    }
	          r *= cm->sft_varied[itpl];
	          bintot += r;
		}
              r = gRandom->Poisson(bintot);
	      if (nbinsy == 1)
		{ pseudodata[ichan]->SetBinContent(ibinx+1,r); }
	      else
		{ pseudodata[ichan]->SetBinContent(ibinx+1,ibiny+1,r); }
	    } // end loop over binsy
	} // end loop over binsx
    } // end loop over channels
}

/*----------------------------------------------------------------------------*/

// calculate the ratio p(data|nmodel)/p(data|dmodel)

Double_t mclimit_csm::weightratio(csm_model *nmodel, csm_model *dmodel, TH1 *hist[])
{
  int nchans = nmodel->channame.size();
  int ichan;
  csm_channel_model *ncm;
  csm_channel_model *dcm;
  int nbinsx,nbinsy,ibin,jbin,ic;
  double wr,pn,pd;
  int dtb,ncc,dcc;

  wr = 0;
  for (ichan=0;ichan<nchans;ichan++)
    {
      ncm = nmodel->chanmodel[ichan];
      ncc = ncm->histotemplate.size();
      dcm = dmodel->chanmodel[ichan];
      dcc = dcm->histotemplate.size();

      nbinsx = hist[ichan]->GetNbinsX();
      nbinsy = hist[ichan]->GetNbinsY();
      for (ibin=0;ibin<nbinsx;ibin++)
	{
	  for (jbin=0;jbin<nbinsy;jbin++)
	    {
	      dtb = 0;
	      if (nbinsy == 1)
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1)));
		}
	      else
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1,jbin+1)));
		}

	      pn = 0; // prediction for numerator model
	      for (ic=0;ic<ncc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
                              ncm->sft_varied[ic];
		    }
		  else
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
                              ncm->sft_varied[ic];
		    }
		}
	      pd = 0; // prediction for denominator model
	      for (ic=0;ic<dcc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
                              dcm->sft_varied[ic];
		    }
		  else
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
                              dcm->sft_varied[ic];
		    }
		}
	      if (pd > 0)
		{
	          wr += dtb*log(pn/pd) - pn + pd;
		}
	    }
	}
    }
  // about the limit of exponentials
  if (wr>680.) 
    { wr = 680.; }
  wr = exp(wr);
  return(wr);
}

/*----------------------------------------------------------------------------*/

// minimize the chisquared function over the nuisance parameters

Double_t mclimit_csm::calc_chi2(csm_model *model, TH1 *hist[])
{
  Int_t i;
  Double_t chisquared;
  csm* mycsm = new csm;
  mycsm->set_modeltofit(model);

  // assume (check?) that the model channel names match up with the data channel names

  for (i=0;i<(Int_t)model->channame.size();i++)
    {
      mycsm->set_htofit(hist[i],model->channame[i]);
    }
  chisquared = mycsm->chisquared();
  delete mycsm;
  return(chisquared);
}

Double_t csm::chisquared()
{
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Double_t arglist[20];
  Int_t ierflag = 0;
  Int_t i,j,icons;
  Double_t cresult;
  Double_t param,paramerror;

  modeltofit->list_nparams(&npn, &nplb, &nphb);
  if (npn.size() > 0)
    {
      
      TMinuit *mnp = new TMinuit(npn.size()+1);

      mnp->SetFCN(csm_minuit_fcn);
      arglist[0] = -1;
      mnp->mnexcm("SET PRINT", arglist, 1, ierflag);
      arglist[0] = 2; 
      mnp->mnexcm("SET STRATEGY", arglist, 1, ierflag);	
      mnp->mnexcm("SET NOW",arglist,1,ierflag);
      mnp->mnexcm("SET NOG",arglist,1,ierflag);
      char npname[10];

      for (i=0;i < (Int_t) npn.size();i++)
        {
          sprintf(npname,"np%d",i);
	  //cout << "setting minuit parameter: " << npname << " " << nplb[i] << " " << nphb[i] << endl;
	  TString npname2 = npname;
          mnp->mnparm(i,npname2,0.0,0.5,nplb[i],nphb[i],ierflag);
          icons = 1;
	  for (j=0;j<(Int_t) modeltofit->npcm.size();j++)
	    {
	      if (strcmp(modeltofit->npcm[j].pnameoutput,npn[i])==0)
		{
		  ierflag = mnp->FixParameter(i);
		  icons = 0;
		}
	    }
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          npfitname.push_back(s);  // this copy is in static global storage so the minuit function knows about it
          if (strstr(npn[i],"UNCONSTRAINED") != 0)
  	    {
	      icons = 0;
	    }
          constrainedfitparam.push_back(icons);
        }

      arglist[0] = 1;
      mnp->mnexcm("SET ERR", arglist ,1,ierflag); 
      ierflag = 0;
      arglist[0] = 500;
      arglist[1] = 1.;

      mnp->mnexcm("SIMPLEX", arglist ,2,ierflag);
      mnp->mnexcm("MIGRAD", arglist ,2,ierflag);

      //cout << "Number of function calls in Minuit: " << mnp->fNfcn << endl;

      // copy best fit parameters for outside use

      cresult = mnp->fAmin;
      //      cout << mnp->fAmin << " fAmin" << endl;
      fitparam.clear();
      fiterror.clear();
      for (i=0;i<(Int_t) fitparamname.size();i++)
        {
          delete[] fitparamname[i];
        }
      fitparamname.clear();
      for (i=0;i < (Int_t) npn.size();i++)
        {
          mnp->GetParameter(i,param,paramerror);
          fitparam.push_back(param);
          fiterror.push_back(paramerror);
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          fitparamname.push_back(s); // this copy's part of the class private members
        }

      delete mnp;
    }
  else
    {
      i = 0;
      csm_minuit_fcn(i,0,cresult,0,0);
    }
  if (cresult < 0)
    { 
      //cout << "chisquared less than zero: " << cresult << " setting it to zero" << endl;
      cresult = 0;
    }

  for (i=0;i<(Int_t) npfitname.size();i++)
    {
      delete[] npfitname[i];
    }
  npfitname.clear();
  constrainedfitparam.clear();

  return(cresult);
}

/*----------------------------------------------------------------------------*/

// A model is a collection of channel models and names

csm_model::csm_model()
{
}

/*----------------------------------------------------------------------------*/

csm_model::~csm_model()
{
  Int_t i,j;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      delete[] channame[i];
      delete chanmodel[i];
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      for (j=0;j<npcm[i].ninput;j++)
	{
	  delete[] npcm[i].pnameinput[j];
	}
      delete[] npcm[i].pnameinput;
      delete[] npcm[i].pnameoutput;
    }

  /* The vectors themselves are deleted when the class instance is deleted */
}

/*----------------------------------------------------------------------------*/

void csm_model::add_template(TH1 *template_hist, //Poisson or non-Poisson histogram
                                Double_t sf,        //scale factor to multiply template by to compare w/ data 
                                                    //(e.g., (data_lum/MC_lum) for a MC Poisson histogram
                                Int_t nnp,          // number of nuisance parameters (Gaussian of unit width)
                                char* npname[],     // nuisance parameter names 
                                Double_t *nps_low,  // fractional uncertainty on sf due to each nuisance parameter -- low side
                                Double_t *nps_high, // fractional uncertainty on sf due to each nuisance parameter -- high side
		                                    // typically nps_low and nps_high are input with opposite signs -- if opposite
                                                    // variations of the nuisance parameter create opposite changes in sf.  The sign
                                                    // is retained in the calculation in case both variations of a nuisance parameter
                                                    // shift the normalization in the same way (either both + or both -)
                                TH1 *lowshape[],    // array of low hisogram shapes, one for each nuisance param (null if no shape error)
                                Double_t *lowsigma, // number of sigma low for each nuisance parameter shape variation
                                TH1 *highshape[],   // array of high histogram shapes, one for each nuisance param (null if no shape error)
	                        Double_t *highsigma, // number of sigma high for each shape variation
                                Int_t pflag,         // Poisson flag -- 1 if Poisson, 0 of not.
                                Int_t sflag,         // scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
                                char *cname)
{
  Int_t i;
  i = lookup_add_channame(cname);
  chanmodel[i]->add_template(template_hist,sf,nnp,npname,nps_low,
                             nps_high,lowshape,lowsigma,highshape,
                             highsigma,pflag,sflag);
}

/*----------------------------------------------------------------------------*/
// add a whole channel's model to the total set of models.

void csm_model::add_chanmodel(csm_channel_model *cm, char *cname)
{
  Int_t ichan;

  ichan = lookup_add_channame(cname);
  chanmodel[ichan] = cm->Clone();
}

/*----------------------------------------------------------------------------*/
// add a constraint function between nuisance parameters.  Make our own copies of all
// the names and the function pointer.

void csm_model::add_npcons(Int_t nparin,char **parin, char *parout, Double_t (*f)(Double_t*))
{
  Int_t i,j;
  npcstruct npc;
  char *s;

  for (i=0;i<(Int_t) npcm.size();i++)
    {
      if (strcmp(parout,npcm[i].pnameoutput)==0)
	{
	  cout << "Warning: Two constraint functions for the same nuisance parameter: " << parout << " defined" << endl;
          exit(0); // bad enough to crash
	}
      for (j=0;j<npcm[i].ninput;j++)
	{
	  if (strcmp(parout,npcm[i].pnameinput[j]) == 0)
	    {
	      cout << "Warning: nuisance parameter: " << npcm[i].pnameinput[j] << " depends on " << parout << endl;
              cout << "but " << parout << "is computed itsef by a constraint after " << npcm[i].pnameinput[j] << "is computed." << endl;
              exit(0); // bad enough to crash
	    }
	}
    }

  npc.ninput = nparin;
  npc.pnameinput = new char*[nparin];
  for (i=0;i<nparin;i++)
    {
      s = new char[strlen(parin[i])+1];
      strcpy(s,parin[i]);
      if (strcmp(s,parout)==0)
	{
	  cout << "Constraint function for nuisance parameter: " << s 
               << " depends on nuisance parameter " << s << endl;
	  exit(0);
	}
      npc.pnameinput[i] = s;
    }
  s = new char[strlen(parout)+1];
  npc.pnameoutput = s;
  npc.f = f;
  npcm.push_back(npc);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::Clone()
{
  Int_t i;
  csm_model* mclone = new csm_model;

  for (i=0;i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    } 
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/
// includes all the new contributing histograms, and also collects together all constraint
// relationships between nuisance parameters.

csm_model* csm_model::add(csm_model &a)
{
  Int_t i;
  csm_model* mclone = a.Clone();
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scale(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scalesignal(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scalesignal(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/

csm_model* csm_model::scale_err(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale_err(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel models and channel names that is sorted by channel */
/* name during the building process.  Return the vector index to use to refer */
/* to this particular channel. */

Int_t csm_model::lookup_add_channame(char *cname)
{
  Int_t i,ifound,j,jfound;
  char *s;
  csm_channel_model *cm; 
  vector<char*>::iterator cni;
  vector<csm_channel_model*>::iterator cmi;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      j = (Int_t) strcmp(cname,channame[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      cm = new csm_channel_model;
      strcpy(s,cname);
      if (jfound == -1)
	{
          ifound = channame.size();
          channame.push_back(s);
          chanmodel.push_back(cm);
	}
      else
	{
	  ifound = jfound;
	  cni = channame.begin() + jfound;
	  channame.insert(cni,s);
	  cmi = chanmodel.begin() + jfound;
	  chanmodel.insert(cmi,cm);
	}
    }

  return(ifound);
}

/*----------------------------------------------------------------------------*/

/* A channel model is a sum of template histograms along with systematic errors */

// constructor

csm_channel_model::csm_channel_model()
{
  chan_istyle = CSM_INTERP_HORIZONTAL;  //  defaults to csm_pvmorph interpolation
}

/*----------------------------------------------------------------------------*/

// destructor

csm_channel_model::~csm_channel_model()
{
  Int_t i;

  //cout << "Called csm_channel_model destructor" << endl;

  // deallocate memory used to save sytematic error names

  for (i=0;i < (Int_t) syserr.size();i++)
    {
      delete[] syserr[i].sysname;
    }
  // deallocate cloned input histograms
  for (i=0;i < (Int_t) histotemplate.size();i++)
    {
      delete histotemplate[i];
      delete histotemplate_varied[i];
    }
  for (i=0;i < (Int_t) syserr.size();i++)
    {
      if (syserr[i].lowshape !=0)
	{
	  delete syserr[i].lowshape;
	}
      if (syserr[i].highshape !=0)
	{
	  delete syserr[i].highshape;
	}
    }
}

void csm_model::print()
{
  Int_t i,j;

  cout << "csm_model::print  -- printing out model information" << endl;
  for (i=0;i<(Int_t) channame.size();i++)
    {
      cout << "Channel: " << i << " Name: " << channame[i] << endl;
      chanmodel[i]->print();
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      cout << "-------------------" << endl;
      cout << "Constraint equation:  " << npcm[i].pnameoutput << " is computed from " << endl;
      for (j=0;j<npcm[i].ninput;j++)
	{
	  cout << npcm[i].pnameinput[j] << endl;
	}
    }
  cout << "-------------------" << endl;
}

void csm_channel_model::print()
{
  Int_t i,j;
  Double_t ssum = 0;
  Double_t bsum = 0;

  cout << "Begin-----------------csm_channel_model::print()------------" << endl;
  for(i=0;i < (Int_t) histotemplate.size();i++)
    {
      cout << "Template " << i << endl;
      cout << "  sft: " << sft[i] << endl;
      cout << "  poissflag: " << poissflag[i] << endl;
      cout << "  signalflag: " << scaleflag[i] << endl;
      cout << "  Integral: " << histotemplate[i]->Integral() << endl;
      cout << "  Scaled Integral: " << histotemplate[i]->Integral()*sft[i] << endl;
      if (scaleflag[i]) 
	{
	  ssum += histotemplate[i]->Integral()*sft[i];
	}
      else
	{
	  bsum += histotemplate[i]->Integral()*sft[i];
	}

      //histotemplate[i]->Print("all");
      Double_t errtotup = 0;
      Double_t errtotdown = 0;
      for (j=0;j< (Int_t) syserr.size();j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      cout << "Syst: " << syserr[j].sysname << endl;
	      cout << "  Up rate error: " << syserr[j].sysfrach << endl;
	      errtotup += syserr[j].sysfrach*syserr[j].sysfrach;
	      errtotdown += syserr[j].sysfracl*syserr[j].sysfracl;
	      cout << "  Down rate error: " << syserr[j].sysfracl << endl;
	      if (syserr[j].highshape != 0)
		{
		  cout << "  Up shape error provided " << syserr[j].xsighigh << endl;
		  cout << "  Down shape error provided " << syserr[j].xsiglow << endl;
		}
	    }
	}
      errtotup = sqrt(errtotup);
      errtotdown = sqrt(errtotdown);
      cout << "Total relative error on this template (up): " << errtotup << endl;
      cout << "Total relative error on this template (down): " << errtotdown << endl;
    }
  cout << "Total signal this channel in this model: " << ssum << endl;
  cout << "Total background this channel in this model: " << bsum << endl;
  cout << "End-------------------csm_channel_model::print()------------" << endl;
}
/*----------------------------------------------------------------------------*/

void csm_channel_model::add_template
                               (TH1 *template_hist, //Poisson or non-Poisson histogram
                                Double_t sf,        //scale factor to multiply template by to compare w/ data 
                                                    //(e.g., (data_lum/MC_lum) for a MC Poisson histogram
                                Int_t nnp,          // number of nuisance parameters (Gaussian of unit width)
                                char* npname[],     // nuisance parameter names 
                                Double_t *nps_low,  // fractional uncertainty on sf due to each nuisance parameter -- low side
                                Double_t *nps_high, // fractional uncertainty on sf due to each nuisance parameter -- high side
		                                    // typically nps_low and nps_high are input with opposite signs -- if opposite
                                                    // variations of the nuisance parameter create opposite changes in sf.  The sign
                                                    // is retained in the calculation in case both variations of a nuisance parameter
                                                    // shift the normalization in the same way (either both + or both -)
                                TH1 *lowshape[],   // array of low hisogram shapes, one for each nuisance param (null if no shape error)
                                Double_t *lowsigma,  // number of sigma low for each nuisance parameter shape variation
                                TH1 *highshape[],   // array of high histogram shapes, one for each nuisance param (null if no shape error)
	                        Double_t *highsigma, // number of sigma high for each shape variation
                                Int_t pflag,         // Poisson flag -- 1 if Poisson, 0 of not.
                                Int_t sflag)        // scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
{
  int i;
  svstruct ses;
  char *s;

  sft.push_back(sf);
  sft_varied.push_back(sf);
  poissflag.push_back(pflag);
  scaleflag.push_back(sflag);
  for (i=0;i<nnp;i++)
    {
       s = new char[strlen(npname[i])+1];
       strcpy(s,npname[i]);
       ses.sysname = s;
       ses.itemplate = histotemplate.size();
       ses.sysfracl = nps_low[i];
       ses.sysfrach = nps_high[i];
       if (lowshape[i] !=0)
	 {
           ses.lowshape = (TH1*) lowshape[i]->Clone();
 	 }
       else
	 {
	   ses.lowshape = 0;
	 }
       if (highshape[i] !=0)
	 {
           ses.highshape = (TH1*) highshape[i]->Clone();
	 }
       else
	 {
	   ses.highshape = 0;
	 }
       ses.xsiglow = lowsigma[i];
       ses.xsighigh = highsigma[i];
       syserr.push_back(ses);

       if (ses.highshape != 0)
	 {
          if (ses.highshape->GetNbinsX() != template_hist->GetNbinsX())
            {
              cout << "Chisquared minmization:  histo template and high shape have different bin counts." << endl;
              cout << template_hist->GetNbinsX() << " != " << ses.highshape->GetNbinsX() << endl;
              exit(0);
	    }
	  if (pflag>0)
	    {
	      if (ses.highshape->GetEntries() != template_hist->GetEntries())
		{
		  // cout << "csm_channel_model::add_template:" << endl;
		  // cout << "warning -- putting in a shape uncertainty histogram" << endl;
		  // cout << "With a different number of entries than the original template histogram" << endl;
		  // cout << "invalidates the Poisson statistics assumption when they are interpolated" <<endl;
		  // cout << "Poisson Template histo and high syst. shape have different numbers of entries." << endl;
                  // cout << "Template: " << template_hist->GetEntries() << " High syst. shape: " 
                  //      << ses.highshape->GetEntries() << endl;
		  // cout << "Continuing anyhow." << endl;
		}
	    }
        }
      if (ses.lowshape != 0)
	{
          if (ses.lowshape->GetNbinsX() != template_hist->GetNbinsX())
            {
              cout << "Chisquared minmization:  histo template and low shape have different bin counts." << endl;
              cout <<  template_hist->GetNbinsX() << " != " << ses.lowshape->GetNbinsX() << endl;
              exit(0);
	    }
	  if (pflag>0)
	    {
	      if (ses.lowshape->GetEntries() != template_hist->GetEntries())
		{
		  // cout << "csm_channel_model::add_template:" << endl;
		  // cout << "putting in a shape uncertainty histogram" << endl;
		  // cout << "With a different number of entries than the original template histogram" << endl;
		  // cout << "invalidates the Poisson statistics assumption when they are interpolated" <<endl;
		  // cout << "Poisson Template histo and low syst. shape have different numbers of entries." << endl;
                  // cout << "Template: " << template_hist->GetEntries() << " Low syst. shape: " 
                  //     << ses.lowshape->GetEntries() << endl;
		  //cout << "Continuing anyhow." << endl;
	        }
	    }
	}
    }
  histotemplate.push_back((TH1*) template_hist->Clone());
  histotemplate_varied.push_back((TH1*) template_hist->Clone());
  //cout << "model::add_template: " << histotemplate.size() << endl;
  //gDirectory->ls();
}

/*----------------------------------------------------------------------------*/

// make a copy of this model by adding the templates over again.
// that way the clone can be deleted by itself, and the destructor
// won't try to delete allocated memory twice

csm_channel_model* csm_channel_model::Clone()
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1 *[syserr.size()];
  TH1 **highshape = new TH1 *[syserr.size()];
  char **ename = new char *[syserr.size()];

  csm_channel_model* mclone = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      //cout << "Model clone adding template " << i << endl;
      mclone->add_template(histotemplate[i],sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }

  mclone->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mclone);
}

/*----------------------------------------------------------------------------*/

// addition of two models makes a new model with the sum of the templates

csm_channel_model* csm_channel_model::add(csm_channel_model &a)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* mclone = a.Clone();

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      mclone->add_template(histotemplate[i],sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  mclone->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mclone);
}

/*----------------------------------------------------------------------------*/

// multiplication of a model and a scalar

csm_channel_model* csm_channel_model::scale(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/

// multiplication of a model and a scalar -- scale the systematic
// errors down with 1/sqrt(coefficient)

csm_channel_model* csm_channel_model::scale_err(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t escale;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  char **ename = new char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  escale = 1.0/sqrt(coefficient);

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i< ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j< nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl*escale;
	      nps_high[nnp] = syserr[j].sysfrach*escale;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow/escale;
	      highsigma[nnp] = syserr[j].xsighigh/escale;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/

// multiplication of only parts of a model and a scalar

csm_channel_model* csm_channel_model::scalesignal(Double_t coefficient)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  char **ename = new char*[syserr.size()];
  Double_t sc1;

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      sc1 = sft[i];
      if (scaleflag[i] != 0)
	{
	  sc1 = coefficient*sft[i];
	}
      smodel->add_template(histotemplate[i],sc1,nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

/*----------------------------------------------------------------------------*/


// Use TMinuit to minimize the chisquared in T. Devlin's note CDF 3126 wrt the
// nuisance parameters.

// updated here -- do a joint minimization over shared nuisance parameters
// for several histograms (channels).

// global (in this file) declarations are at the top of the source file

// constructor
csm::csm()
{
  //  cout << "Chisquared Minimizer Constructor called\n";
  datatofit.clear();
  datatofitname.clear();
  constrainedfitparam.clear();
  npfitname.clear();
}

//destuctor

csm::~csm()
{
  Int_t i;
  //cout << "Chisquared Minimizer Destructor called\n";

  // clear out static global variables

  for (i=0;i<(Int_t) datatofit.size(); i++)
    {
      delete datatofit[i];
      delete[] datatofitname[i];
    }
  datatofit.clear();
  datatofitname.clear();

  // clear out allocated memory pointed to by our private members

  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      delete[] fitparamname[i];
    }
}

// put in the data histograms in the same order we built up the model histograms

void csm::set_htofit(TH1 *h, char *cname)
{
  Int_t i,ifound,j,jfound;
  vector<char*>::iterator dni;
  vector<TH1*>::iterator dfi;
  char *s;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) datatofitname.size(); i++)
    {
      j = (Int_t) strcmp(cname,datatofitname[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings.  If the name is on the 
     list, replace the existing data histogram with a clone of the one supplied. */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      strcpy(s,cname);
      if (jfound == -1)
	{
          datatofitname.push_back(s);
          datatofit.push_back((TH1*) h->Clone());
	}
      else
	{
          dni = datatofitname.begin() + jfound;
	  datatofitname.insert(dni,s);
	  dfi = datatofit.begin() + jfound;
	  datatofit.insert(dfi,(TH1*) h->Clone());
	}
    }
  else
    {
      delete datatofit[ifound];
      datatofit[ifound] = (TH1*) h->Clone();
    }
}

void csm::set_modeltofit(csm_model* mtf)
{
  modeltofit = mtf;
}

csm_model* csm::getbestmodel()
{
  Int_t i;

  // make a local array of pointers to nuisance parameter names
  
  char **fpnameloc = new char *[fitparamname.size()];
  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      fpnameloc[i] = fitparamname[i];
      //cout << "in getbestmodel, paramname: " << fpnameloc[i] << endl;
    }
  Double_t *parloc = new Double_t[fitparam.size()];
  for (i=0;i<(Int_t) fitparam.size();i++)
    {
      parloc[i] = fitparam[i];
      //cout << "in getbestmodel, param: " << parloc[i] << endl;
    }

  modeltofit->nuisance_response(fitparam.size(),fpnameloc,parloc);
  delete[] fpnameloc;
  delete[] parloc;
  return(modeltofit);
}

void csm::plotcompare(char *cname)
{
  Int_t i;
  for (i=0;i<(Int_t)datatofitname.size();i++)
    {
      if (strcmp(datatofitname[i],cname)==0)
	{
	  modeltofit->plotwithdata(cname,datatofit[i]);
	}
    }
}

// Number of degrees of freedom -- this is approximately true for
// large statistics (in fact, the whole chisquared idea is only approximately
// true in cases of large statistics where distributions are approximately Gaussian)
// Degrees of freedom "freeze out" as the expected number of events gets small
// (<5 or so).  A bin with no data and no expectation shouldn't contribute either
// to the chisquared or the number of degrees of freedom, and neither really should
// a bin with 1E-6 expected and no observed events.  This routine won't draw the
// line (and even interpolated histograms can have variable numbers of bins with
// zero expectation).  This routine's very naive and just counts bins, filled or not.

Int_t csm::ndof()
{
  Int_t ndofl,i;
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;

  ndofl = 0;
  for (i=0;i<(Int_t) datatofit.size();i++)
    {
      ndofl += datatofit[i]->GetNbinsX()*datatofit[i]->GetNbinsY();
    }
  modeltofit->list_nparams(&npn, &nplb, &nphb);
  ndofl -= npn.size();
  cout << "nDOF isn't very clearly defined here... todo" << endl;
  return(ndofl);
}


Int_t csm::getnparams()
{
  return(fitparam.size());
}

Double_t csm::getparam(Int_t iparam)
{
  return(fitparam[iparam]);
}

Double_t csm::getperror(Int_t iparam)
{
  return(fiterror[iparam]);
}

char* csm::getpname(Int_t iparam)
{
  return(fitparamname[iparam]);
}

// Call the individual channel chisquared calculators inside here.
// the parameters par are labeled by their names npfitname in the static global vector.

void csm_minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Int_t i;

  //adjust the model according to the nuisance paramters
  
  char **fpnameloc = new char *[npfitname.size()];
  for (i=0;i<(Int_t) npfitname.size();i++)
    {
      fpnameloc[i] = npfitname[i];
      //cout << "in minuit fit fcn: " << i << " " << npfitname[i] << endl;
    }
  modeltofit->nuisance_response(npfitname.size(),fpnameloc,par);

  TH1** dfloc = new TH1*[datatofit.size()];
  for (i=0;i<(Int_t) datatofit.size();i++)
    {
      dfloc[i] = datatofit[i];
    }

  //cout << "In minuit function: printing out the model" << endl;
  //mfluct->print();

  f = modeltofit->chisquared1(dfloc);

  // Gaussian constraints for variables which are constrained.

  for (i=0;i<npar;i++) 
    {
      if (constrainedfitparam[i] != 0)
        {
	  //cout << "In fcn: " << i << " " << par[i] << endl;
          f += par[i]*par[i];
        }
    }

  //cout << "end of computation of f in minuit_fit_fcn: " << f << endl;

  delete[] fpnameloc;
  delete[] dfloc;
}

/*------------------------------------------------------------------------*/
// make a plot of the results, along with some data

void csm_channel_model::plotwithdata(TH1* dh)
{
  Int_t i,ntemplates,nbinsy;
  Double_t stackmax,datamax,plotmax;
  THStack *hs = new THStack("hs",dh->GetTitle());
  ntemplates = (Int_t) histotemplate.size();
  TLegend *slegend = (TLegend*) new TLegend(0.7,0.6,0.89,0.89);

  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i]->Clone();
      htl->Scale(sft_varied[i]);
      htl->SetFillColor(i+40);
      htl->SetFillStyle(1001);
      hs->Add(htl);
    }

  TList *hlist = hs->GetHists();
  TObjLink *lnk = hlist->LastLink();          
  while (lnk)
    {  slegend->AddEntry(lnk->GetObject(),lnk->GetObject()->GetName(),"F");
       lnk = lnk->Prev();                       
    }     
  // make sure the plot is big enough to fit the data, the model stack,
  // and the data error bars with a little room to spare
  stackmax = hs->GetMaximum();
  datamax = dh->GetMaximum();
  nbinsy = dh->GetNbinsY();
  datamax += sqrt(datamax);
  plotmax = max(datamax,stackmax);
  hs->SetMaximum(plotmax);
  if (nbinsy==1)
    {
      hs->Draw("HIST");
      dh->SetMarkerStyle(20);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("E0SAME");
    }
  else
    {
      hs->Draw();
      dh->SetMarkerStyle(20);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("LEGO,E0,SAME");
    }
  slegend->AddEntry(dh,dh->GetName(),"P");
  slegend->SetHeader(dh->GetTitle());
  slegend->Draw();
}


double csm_channel_model::kstest(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  TH1* hsum = (TH1*) histotemplate_varied[0]->Clone();
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh);
  delete hsum;
  return(tout);
}

double csm_channel_model::kstest_px(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  TH1* hsum = (TH1*) histotemplate_varied[0]->Clone();
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh,"X");
  delete hsum;
  return(tout);
}

/*------------------------------------------------------------------------*/

// and a method to allow an object of type csm_model to plot up one of its
// channels with some data compared.

void csm_model::plotwithdata(char* cname, TH1* dh)
{
  Int_t i;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  chanmodel[i]->plotwithdata(dh);
	}
    }
}
 
/*------------------------------------------------------------------------*/

double csm_model::kstest(char* cname, TH1* dh)
{
  Int_t i;
  double ksresult=0; 
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest(dh);
	}
    }
  return(ksresult);
}
 
/*------------------------------------------------------------------------*/

double csm_model::kstest_px(char* cname, TH1* dh)
{
  Int_t i;
  double ksresult=0;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest_px(dh);
	}
    }
  return(ksresult);
}
 
/*------------------------------------------------------------------------*/

/* chisquared1 evaluates a chisquared function in the style of T. Devlin's CDF 3126, eq's 9 and 10
   The signal is a sum of signal contributions and the background is a sum of
   background contributions.  This chisquared function is meant to be minimized 
   with respect to the free nuisance parameters 

   This version does one 1D or 2D histogram at a time.

   This version allows for multiple sources of signal and multiple sources of background,
   some of each of which are estimated using finite MC or data statistics in each bin.
   This routine does not distinguish between a signal source and a background source --
   finding the chisquared of a data distribution to a sum of models does not need a distinction
   at this level.  Instead, one may compare the chisquared of the same data against collections
   of models that include signals and those that do not include signals, calling this routine
   twice (or more times).

   CDF 3126 describes how to minimize the chisquared function over each bin's uncertain
   Poisson-constrained rates.  When multiple sources are allowed to be estimated from Poisson
   distributed subsidiary measurements, the quadratic polynomial to be solved for turns
   into a system of quadratic equations which is solved here iteratively.

   This function is meant to be part of a MINUIT minimization over the nuisance
   parameters.


            input: TH1 *dh -- data histogram to compare the channel's model against
            
            output:  chi squared, the function value.

   Update 5 July 2006 -- Reading Barlow and Beeston about bins with zero MC prediction in one
   or more source.  Take the one with the strongest contribution (here taken from the normalization
   scale factors), and set the others to zero when solving the n coupled quadratic equations.
*/

Double_t csm_channel_model::chisquared1(TH1 *dh)
{
  Double_t chi2;
  Int_t ip1,ic,ibinx,ibiny,iter,iprec;
  Double_t A,B,C,D;
  Double_t csum,cpsum,gbc;
  Int_t nbinsx,nbinsy;
  Int_t nsubs;
  Int_t dtb;  // data observed in a single bin
  Int_t nc;

  // number of template histograms
  nc = (Int_t) histotemplate.size();

  nbinsx = dh->GetNbinsX();
  nbinsy = dh->GetNbinsY();

  // allocate rho1 and rho2 for all templates, even though we're only going to need
  // them for the Poisson-distributed ones

  Double_t* rho1 = new Double_t[nc];
  Double_t* rho2 = new Double_t[nc];
  int *zlist = new int[nc];

  chi2 = 0;

  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
        {
	  if (nbinsy==1)
	    { dtb = (Int_t) dh->GetBinContent(ibinx+1); }
	  else
	    { dtb = (Int_t) dh->GetBinContent(ibinx+1,ibiny+1); }
          //cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;

	  /* the sum of non-Poisson contributions in this bin, varied by the nuisance parameters */
          csum = 0;
          for (ic=0;ic<nc;ic++)
 	    {
	      if (poissflag[ic] == 0)
	        {
		  if (nbinsy == 1)
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
	          csum += gbc*sft_varied[ic];
	        }
	    }
	  if (csum < 0) 
	    { chi2 += 1E10; }

	  /* solve for the rho's in each bin for each source of Poisson-estimated model rate
	     rho1 is the current estimate used to compute the rho2's.  On each iteration,
	     copy the previous iteration's rho2 into the rho1 array and re-solve for rho2.
	     start with nominal central values from the subsidiary measurements */

	  int haszero = 0;
          int im1=0;
          double xm1=0;
          for (ic=0;ic<nc;ic++) 
	    { 
	      if (poissflag[ic] != 0)
		{
		  if (nbinsy == 1)
		    { gbc = max(0,nearbyint(histotemplate_varied[ic]->GetBinContent(ibinx+1))); }
		  else
		    { gbc = max(0,nearbyint(histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1))); }
		  rho2[ic] = gbc*sft_varied[ic];
		  rho1[ic] = rho2[ic];
		  if (gbc == 0)
		    { 
		      haszero = 1;
		      zlist[ic] = 1;
		      if (sft_varied[ic]>xm1)
			{ 
			  xm1 = sft_varied[ic]; 
			  im1 = ic;
			}
		    }
		  else
		    {
		      zlist[ic] = 0;
		    }
		}
	      else
		{ 
		  rho2[ic] = 0;
		  rho1[ic] = 0;
		}
	    }
	  if (haszero != 0)
	    {
	      zlist[im1] = 0;
	    }

          for (iter=0;iter<MAXITER;iter++)
	    {
              for (ic=0;ic<nc;ic++) 
                { 
                  rho1[ic] = rho2[ic];
		  if (zlist[ic] == 1)
		    { 
		      rho1[ic] = 0;
		    }
                }

	      for (ic=0;ic<nc;ic++)
	        {
	          if (poissflag[ic] != 0)
	            {
		      if (zlist[ic] == 0)
			{
	                  D = csum;
                          for (ip1=0;ip1<nc;ip1++)
		            { if (poissflag[ip1] && ip1 != ic) D += rho1[ip1]; }
	                  A = 1.0 + 1.0/sft_varied[ic];
		          if (nbinsy == 1)
			    { gbc = max(0,nearbyint(histotemplate_varied[ic]->GetBinContent(ibinx+1))); }
		          else
			    { gbc = max(0,nearbyint(histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1))); }
		          B = A*D - dtb - gbc;
		          C = -gbc*D;
	                  rho2[ic] = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
			}
		      else // a la Barlow and Beeston, set only one prediction to nonzero
			{                    // if we have zero MC -- the "strongest" one among all the contributions
			  rho2[ic] = 0;     // with zero MC prediction
			}
		    }
	        }
	      iprec = 0;

	      for (ic=0;ic<nc;ic++)
	        {
	          if (poissflag[ic] != 0)
	            {
	              if (fabs(rho1[ic]) < PREC1)
		        { 
		          if (fabs(rho2[ic]-rho1[ic]) > PREC1)
		            { 
                              iprec = 1;
		              break;
		            }
		        }
	              else
		        {
		          if ( fabs((rho2[ic]-rho1[ic])/rho1[ic])>PREC1 )
		            {
		              iprec = 1;
		              break;
		            }
		        }
	            }
	        }
              if (iprec == 0) break;
  	    }  /* end loop over iterations to compute the rho's.  rho2 is the computed array */

	  /*
	  if (CSM_DEBUGPRINT >0 && iprec ==1)
	    {
	      // cout << "csm_chisquared1: iterations failed to converge " << endl;
	      cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;
	      for (ic=0;ic<nc;ic++)
		{
		  if (nbinsy == 1)
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
		  if (poissflag[ic] != 0)
		    {
		      cout << "Poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
		    }
		  else
		    {
		      cout << "Non-poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
		    }
		}
	    }
	  */
	  // When the iterations fail to converge, it is usually an oscillatory
	  // solution.  Pick the rho1 or the rho2 array which minimizes chisquare
	  // first compute the chisquare using the rho2 array, and if we need to,
          // redo it with the rho1 array, and pick the smaller of the two.

	  Double_t chi2a = chi2;

	  cpsum = csum;
	  for (ic=0;ic<nc;ic++)
	    {
	      if (poissflag[ic] != 0)
		{
		  cpsum += rho2[ic];
		}
	    }
	  if (cpsum < 0)
	    { chi2a += 1E10; }
	  else
	    {
	      chi2a += cpsum;
	      chi2a -= dtb;
	      if (dtb>0) {chi2a -= dtb*log(cpsum/((Double_t) dtb));}
	    }

	  for (ic=0;ic<nc;ic++)
	    {
	      if (poissflag[ic] != 0)
		{
		  if (nbinsy == 1)
		    { nsubs = (Int_t) histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { nsubs = (Int_t) histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }

		  /*
		    cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
		    rho2[ic] << " " << sft_varied[ic] << " " << nsubs << endl;
		  */
		  chi2a += (rho2[ic]/sft_varied[ic] - nsubs);
		  if (nsubs > 0 && rho2[ic] > 0)
		    chi2a -= ((Double_t) nsubs)*log(rho2[ic]/(sft_varied[ic]*((Double_t) nsubs)));
		}
	    }

	  Double_t chi2b = chi2a;

	  if (iprec == 1)
	    {
	      cpsum = csum;
	      for (ic=0;ic<nc;ic++)
		{
		  if (poissflag[ic] != 0)
		    {
		      cpsum += rho1[ic];
		    }
		}
	      if (cpsum < 0)
		{ chi2b += 1E10; }
	      else
		{
		  chi2b += cpsum;
		  chi2b -= dtb;
		  if (dtb>0) {chi2b -= dtb*log(cpsum/((Double_t) dtb));}
		}

	      for (ic=0;ic<nc;ic++)
		{
		  if (poissflag[ic] != 0)
		    {
		      if (nbinsy == 1)
			{ nsubs = (Int_t) histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		      else
			{ nsubs = (Int_t) histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }

		      /*
			cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
			rho1[ic] << " " << sft_varied[ic] << " " << nsubs << endl;
		      */
		      chi2b += (rho1[ic]/sft_varied[ic] - nsubs);
		      if (nsubs > 0 && rho1[ic] > 0)
			chi2b -= ((Double_t) nsubs)*log(rho1[ic]/(sft_varied[ic]*((Double_t) nsubs)));
		    }
		}
	    }
	  chi2 = min(chi2a,chi2b);
        } /* end loop over binsy */
    } /* end loop over binsx */

  chi2 *= 2.0;

  // cout << "chisquared calc: " << chi2 << endl;

  delete[] rho1;
  delete[] rho2;
  delete[] zlist;
  return(chi2);

}

Double_t csm_model::chisquared1(TH1 **dh)
{
  Int_t i;
  Double_t cs;
  cs = 0;
  for (i=0;i<(Int_t)chanmodel.size();i++)
    {
      cs += chanmodel[i]->chisquared1(dh[i]);
    }
  return(cs);
}

/*------------------------------------------------------------------------*/
//  Set the interpolation style for a particular channel.  Two methods    
//  one for channel models, and one if you just have a pointer to a csm_model
/*------------------------------------------------------------------------*/


void csm_model::set_interpolation_style(char *cname, INTERPSTYLE istyle)
{
  Int_t i;
  for (i=0;i<(Int_t) channame.size(); i++)
    {
      if (strcmp(channame[i],cname)==0)
	{
	  chanmodel[i]->set_interpolation_style(istyle);
	}
    }
}

void csm_channel_model::set_interpolation_style(INTERPSTYLE istyle)
{
  chan_istyle = istyle;
}

/*------------------------------------------------------------------------*/

// interpolate 1D histograms and 2D histograms
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output

void csm_interpolate_histogram(TH1* a, Double_t xa, 
                               TH1* b, Double_t xb,
                               TH1* c, Double_t xc,
                               INTERPSTYLE istyle)
{
  Double_t xmina,xminb,xmaxa,xmaxb,xminc,xmaxc,hnorma,hnormb,hnormc,hnormci;
  Int_t nbinsa,nbinsb,nbinsc,i,j;
  Int_t nbinsya,nbinsyb,nbinsyc;
  Double_t ymina,yminb,ymaxa,ymaxb,yminc,ymaxc;
  Double_t gbc;

  nbinsa = a->GetNbinsX();
  nbinsb = b->GetNbinsX();
  nbinsc = c->GetNbinsX();
  nbinsya = a->GetNbinsY();
  nbinsyb = a->GetNbinsY();
  nbinsyc = a->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      cout << "nbins mismatch1 in csm_interpolate_histograms: " << nbinsa << " " << nbinsb << endl;
    }
  if (nbinsb != nbinsc)
    {
      cout << "nbins mismatch2 in csm_interpolate_histograms: " << nbinsb << " " << nbinsc << endl;
    }
  if (nbinsya != nbinsyb)
    {
      cout << "nbinsy mismatch1 in csm_interpolate_histograms: " << nbinsya << " " << nbinsyb << endl;
    }
  if (nbinsyb != nbinsyc)
    {
      cout << "nbinsy mismatch2 in csm_interpolate_histograms: " << nbinsyb << " " << nbinsyc << endl;
    }

  if (xb == xa)
    {
      cout << "xb == xa in csm_interpolate_histogram " << xa << endl;
      cout << "fatal error -- exiting." << endl;
      exit(0);
    }

  if (a->Integral()<=0 || b->Integral()<=0)
    { 
      c->Reset();
      return;
    }
    
  if (nbinsya == 1)
    {
      xmina = a->GetXaxis()->GetXmin();
      xminb = b->GetXaxis()->GetXmin();
      xminc = c->GetXaxis()->GetXmin();
      xmaxa = a->GetXaxis()->GetXmax();
      xmaxb = b->GetXaxis()->GetXmax();
      xmaxc = c->GetXaxis()->GetXmax();

      Double_t *dista = new Double_t[nbinsa];
      Double_t *distb = new Double_t[nbinsb];
      Double_t *distc = new Double_t[nbinsc];

      hnorma = 0;
      hnormb = 0;
      for (i=0;i<nbinsa;i++)
        { dista[i] = a->GetBinContent(i+1); hnorma += dista[i]; }
      for (i=0;i<nbinsb;i++)
        { distb[i] = b->GetBinContent(i+1); hnormb += distb[i]; }

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);

      if (istyle == CSM_INTERP_HORIZONTAL)
	{
           csm_pvmorph(&nbinsa, &xmina, &xmaxa, dista,
                       &nbinsb, &xminb, &xmaxb, distb,
                       &nbinsc, &xminc, &xmaxc, distc,
                       &xa, &xb, &xc);

           hnormci = 0;
           for (i=0;i<nbinsc;i++) { hnormci += distc[i]; }

           for (i=0;i<nbinsc;i++)
           {
             c->SetBinContent(i+1,distc[i]*hnormc/hnormci);
           }
	}
      else if (istyle == CSM_INTERP_VERTICAL)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = dista[i] + ((xc-xa)/(xb-xa))*(distb[i]-dista[i]);
	      if (gbc < 0) 
		{
		  gbc = 0;
		}
	      c->SetBinContent(i+1,gbc);
	    }
	}
      else
	{
	  cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << endl;
	  exit(0);
	}

      //cout << xa << " " << xb << " " << xc << endl;

      delete[] dista;
      delete[] distb;
      delete[] distc;
    }
  else         // 2d case
    {
      xmina = a->GetXaxis()->GetXmin();
      xminb = b->GetXaxis()->GetXmin();
      xminc = c->GetXaxis()->GetXmin();
      xmaxa = a->GetXaxis()->GetXmax();
      xmaxb = b->GetXaxis()->GetXmax();
      xmaxc = c->GetXaxis()->GetXmax();

      ymina = a->GetYaxis()->GetXmin();
      yminb = b->GetYaxis()->GetXmin();
      yminc = c->GetYaxis()->GetXmin();
      ymaxa = a->GetYaxis()->GetXmax();
      ymaxb = b->GetYaxis()->GetXmax();
      ymaxc = c->GetYaxis()->GetXmax();

      Double_t *distxya = new Double_t[nbinsa*nbinsya];
      Double_t *distxyb = new Double_t[nbinsb*nbinsyb];
      Double_t *distxyc = new Double_t[nbinsc*nbinsyc];

      hnorma = 0;
      hnormb = 0;
      for (j=0;j<nbinsya;j++)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = a->GetBinContent(i+1,j+1);
	      distxya[i+nbinsa*j] = gbc;
	      hnorma += gbc;
	    }
	}
      hnormb = 0;
      for (j=0;j<nbinsyb;j++)
	{
	  for (i=0;i<nbinsb;i++)
	    {
	      gbc = b->GetBinContent(i+1,j+1);
	      distxyb[i+nbinsb*j] = gbc;
	      hnormb += gbc;
	    }
	}

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);

      if (istyle == CSM_INTERP_HORIZONTAL)
	{
          csm_pvmorph_2d(&nbinsa,  &xmina, &xmaxa,
                         &nbinsya, &ymina, &ymaxa,
                         distxya,
                         &nbinsb,  &xminb, &xmaxb,
                         &nbinsyb, &yminb, &ymaxb,
                         distxyb,
                         &nbinsc,  &xminc, &xmaxc,
                         &nbinsyc, &yminc, &ymaxc,
                         distxyc,
                         &xa, &xb, &xc);

          hnormci = 0;
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          hnormci += distxyc[i+nbinsc*j];
	        }
	    }
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          c->SetBinContent(i+1,j+1,distxyc[i+nbinsc*j]*hnormc/hnormci);
	        }
	    }
	}
      else if (istyle == CSM_INTERP_VERTICAL)
	{
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
		  gbc = distxya[i+nbinsc*j] + ((xc-xa)/(xb-xa))*(distxyb[i+nbinsc*j]-distxya[i+nbinsc*j]);
		  if (gbc < 0)
		    {
		      gbc = 0;
		    }
	          c->SetBinContent(i+1,j+1,gbc);
	        }
	    }
	}
      else
	{
	  cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << endl;
	  exit(0);
	}

      delete[] distxya;
      delete[] distxyb;
      delete[] distxyc;
    }

}

/*------------------------------------------------------------------------*/


/* csm_pvmorph.c -- original C'd from d_pvmorph.f using f2c 
   -- trj 5 July 2005 -- remove f2c's ugly i/o and replace
   with printf's */


/* .================================================================== */
void csm_pvmorph(Int_t *nb1, Double_t *xmin1, Double_t *xmax1, Double_t *dist1,
                 Int_t *nb2, Double_t *xmin2, Double_t *xmax2, Double_t *dist2,
                 Int_t *nbn, Double_t *xminn, Double_t *xmaxn, Double_t *distn,
                 Double_t *par1, Double_t *par2, Double_t *parn)
{
    /* f2c generated locals */
    Int_t i__1, i__2, i__3;
    Double_t d__1, d__2, d__3;

    /* Local variables */
    static Int_t allzero1, allzero2;
    static Double_t xdis[(PVMORPH_MAXBINS + 1)];
    static Int_t i__;
    static Double_t x, y, xmind, xmaxd, total, yprev, x1, x2, xmin1d, 
	    xmin2d, xmax1d, xmax2d;
    static Int_t i12;
    static Double_t x10, x20, x21, dx, y20, y21;
    static Int_t idebug;
    static Double_t x11;
    static Int_t ix;
    static Double_t y10, y11, sigdis[20004]	/* was [(PVMORPH_MAXBINS + 1)][4] */, dx1, dx2;
    static Int_t ix1, ix2, ix3, nx3, ixf;
    static Double_t wta, wtb;
    static Int_t ixl, ix1l, ix2l;

   
/* -------------------------------------------------------------------------- */
/* Author           : Alex Read */
/* Version 2.2 */
/* Date             : 28 November 2001 */
/* Last modification: 03 January 2005 (previous 05 April 2004) */

/*      Perform a linear interpolation between two probability distribution */
/*      histogram as a function of the characteristic parameter of the */
/*      distribution. */

/*      The algorithm is described in Read, A. L., "Linear Interpolation */
/*      of Histograms",NIM A 425 (1999) 357-360. */

/*      This is a completely new implementation which does not have the */
/*      limitation of ignoring structures below a bin probability of 10**-5. */
/*      In addition to changing the implementation, the use of double */
/*      precision allows pdf's to be accurately interpolated down to */
/*      something like 10**-15. */

/*      In this particular test-version it is required that the two */
/*      input histograms have identical binning. I have identified a */
/*      strategy which will allow precise interpolations down to */
/*      ~10**-36 but have not implemented it yet. */

/*      05.04.04 Protection against interpolating distributions with more */
/*               than nxmax bins (currently PVMORPH_MAXBINS). Insure endpoint of final */
/*               cdf is set correctly even if GE fails due to numerical */
/*               precision. */

/*      03.01.05 Maarten.Boonekamp@cern.ch and Vanina.Ruhlmann@in2p3.fr */
/*               found problems with sparse distributions. It occurred that */
/*               some bins end up with negative contents. After a careful */
/*               reading */
/*               of the code I found that some code treating empty bins was */
/*               irrelevant and wrong. Tests of distributions with a small */
/*               number of bins with contents revealed other problems related */
/*               to numerical precision and if-tests. */

/* nb1,nb2,nbn are INTEGER (int), all other arguments are REAL*4 (float). */

/* par1,par2 : The values of the linear parameter that characterises the */
/*             histograms (e.g. a particle mass). */
/* parn      : The value of the linear parameter we wish to interpolate to. */
/* nb1,xmin1,xmax1 : The number of bins and the range of the first input pdf. */
/*                 : Maximum number of bins is PVMORPH_MAXBINS. */
/* nb2,xmind,xmaxd : The number of bins and the range of the second input pdf. */
/*                 : Maximum number of bins is PVMORPH_MAXBINS. */
/* nbn,xminn,xmaxn : The number of bins and the range of the output */
/*                   interpolated pdf. On input nbn specifies the maximum */
/*                   size of the output array distn and returns the actual */
/*                   size of distn. */
/* dist1,dist2 : Arrays containing the binned inputs pdf's. */
/* distn : An array containing the output interpolated pdf. */

/* ------------------------------------------------------------------------ */

/* ......Arguments */


/* ......Local variables */

/* nxmax bins have */
/* nxmax+1 edges */

/* ......External functions */


/* ......Expert debugging on=1. */

    /* Parameter adjustments */
    --dist1;
    --dist2;
    --distn;

    /* Function Body */
    idebug = 0;

    /*
    cout << "nbins: " << *nb1 << " " << *nb2 << " " << *nbn << endl;
    cout << "xmin: " << *xmin1 << " " << *xmin2 << " " << *xminn << endl;
    cout << "xmax: " << *xmax1 << " " << *xmax2 << " " << *xmaxn << endl;
    cout << "par: " << *par1 << " " << *par2 << " " << *parn << endl;
    for (Int_t ipr1=0;ipr1<*nb1;ipr1++)
      {
	cout << "bin: " << ipr1 << " " << dist1[ipr1] << " " << dist2[ipr1] << endl;
      }
    */

/* ......Enforce limit on number of input bins. */

    if (max(*nb1,*nb2) > PVMORPH_MAXBINS) {
      printf("-E-phmorph More than PVMORPH_MAXBINS bins in input histograms\n");
	goto L999;
    }

/* From Tom Junk (LEP Higgs Working Group): */
/* special case -- zero histogram contents -- sometimes comes up in 2D */
/* interpolations */

    allzero1 = 1;
    allzero2 = 1;
    i__1 = *nb1;
    for (ix = 1; ix <= i__1; ++ix) {
	if (dist1[ix] > 0) {
	    allzero1 = 0;
	}
    }
    i__1 = *nb2;
    for (ix = 1; ix <= i__1; ++ix) {
	if (dist2[ix] > 0) {
	    allzero2 = 0;
	}
    }

    if ( (allzero1 !=0) || (allzero2 !=0)) {
	*nbn = max(*nb1,*nb2);
	*xminn = fmin(*xmin1,*xmin2);
	*xmaxn = fmax(*xmax1,*xmax2);
	i__1 = *nbn;
	for (ix = 1; ix <= i__1; ++ix) {
	    distn[ix] = 0;
	}
	goto L999;
    }

/* ......The weights (wta,wtb) are the "distances" between the values of the */
/*      parameters at the histograms and the desired interpolation point. */
/*      Check that they make sense. If par1=par2 then we can choose any */
/*      valid set of wta,wtb so why not take the average? */

    if (*par2 != *par1) {
	wta = 1. - (*parn - *par1) / (*par2 - *par1);
	wtb = (*parn - *par2) / (*par2 - *par1) + 1.;
    } else {
	wta = .5;
	wtb = .5;
    }
    d__1 = 1. - (wta + wtb);
    if (wta < 0.
        || wta > 1.
        || wtb < 0.
        || wtb > 1.
        || fabs(d__1) > 1e-4) {
      printf("-W-csm_pvmorph This is an extrapolation,\n"); 
      printf(" weights are: %lf %lf  sum=%lf\n",wta,wtb,wta+wtb);
    }

/* ......Perform interpolation of histogram bin parameters. Use */
/*      assignments instead of computation when input binnings */
/*      are identical to assure best possible floating point precision. */

    if (min(wta,wtb) >= 0.) {
	if (*xmin1 == *xmin2) {
	    *xminn = *xmin1;
	} else {
	    *xminn = wta * *xmin1 + wtb * *xmin2;
	}
	if (*xmax1 == *xmax2) {
	    *xmaxn = *xmax1;
	} else {
	    *xmaxn = wta * *xmax1 + wtb * *xmax2;
	}
	if (*nb1 == *nb2) {
	    *nbn = *nb1;
	} else {
	    *nbn = (int) (wta * *nb1 + wtb * *nb2);
	}
    } else {

/* ......If both weights are zero, then do something simple but */
/*      but reasonable with the histogram bin parameters. */

	*xminn = fmin(*xmin1,*xmin2);
	*xmaxn = fmax(*xmax1,*xmax2);
	*nbn = max(*nb1,*nb2);
    }

/* ......Make double precision copy of histogram bin parameters to */
/*      use in interpolation computation. */

    xmin1d = *xmin1;
    xmin2d = *xmin2;
    xmind = *xminn;
    xmax1d = *xmax1;
    xmax2d = *xmax2;
    xmaxd = *xmaxn;
    dx1 = (xmax1d - xmin1d) / (*nb1);
    dx2 = (xmax2d - xmin2d) / (*nb2);
    dx = (xmaxd - xmind) / (*nbn);
    if (idebug >= 1) {
      printf("PVMORPH: %lf %lf %lf %d\n",xmin1d,xmax1d,dx1,*nb1);
      printf("PVMORPH: %lf %lf %lf %d\n",xmin2d,xmax2d,dx2,*nb2);
      printf("PVMORPH: %lf %lf %lf %d\n",xmind,xmaxd,dx,*nbn);
    }

/* ......Extract the single precision pdf's into double precision arrays */
/*      for interpolation computation. The offset is because sigdis(i) */
/*      describes edge i (there are nbins+1 of them) while dist1/2 */
/*      describe bin i. */

    i__1 = *nb1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigdis[i__] = dist1[i__];
    }
    i__1 = *nb2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigdis[i__ + (PVMORPH_MAXBINS + 1)] = dist2[i__];
    }

/* ......Normalize pdf's to 1 and integrate to find cdf's. First */
/*      bin (edges 1 and 2) require special treatment. */

    total = 0.;
    i__1 = *nb1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	total += sigdis[i__];
    }
    sigdis[0] = 0.;
    sigdis[1] /= total;
    i__1 = *nb1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	sigdis[i__] = sigdis[i__] / total + sigdis[i__ - 1];
    }

    total = 0.;
    i__1 = *nb2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	total += sigdis[i__ + (PVMORPH_MAXBINS + 1)];
    }
    sigdis[(PVMORPH_MAXBINS + 1)] = 0.;
    sigdis[(PVMORPH_MAXBINS + 2)] /= total;
    i__1 = *nb2;
    for (i__ = 2; i__ <= i__1; ++i__) {
	sigdis[i__ + (PVMORPH_MAXBINS + 1)] = sigdis[i__ + (PVMORPH_MAXBINS + 1)] / total + sigdis[i__ + PVMORPH_MAXBINS];
    }

/* ......We are going to step through all the edges of both input */
/*      cdf's ordered by increasing value of y. We start at the */
/*      lower edge, but first we should identify the ends of the */
/*      curves. */

/* ......Find the ends of the curves (the first point from the right */
/*      that has the same integral as the last edge. */

    ix1l = *nb1 + 1;
    ix2l = *nb2 + 1;
    while(sigdis[ix1l - 2] >= sigdis[ix1l - 1]) {
	--ix1l;
    }
    while(sigdis[ix2l + (PVMORPH_MAXBINS - 1)] >= sigdis[ix2l + PVMORPH_MAXBINS]) {
	--ix2l;
    }

/* ......Step up to the beginnings of the curves. These are the */
/*      last non-zero points from the left. */

    ix1 = 1;
    ix2 = 1;
    while(sigdis[ix1] <= sigdis[0]) {
	++ix1;
    }
    while(sigdis[ix2 + (PVMORPH_MAXBINS + 1)] <= sigdis[(PVMORPH_MAXBINS + 1)]) {
	++ix2;
    }

/* ......The first interpolated point is computed now. */

    nx3 = 1;
    x1 = xmin1d + (double) (ix1 - 1) * dx1;
    x2 = xmin2d + (double) (ix2 - 1) * dx2;
    x = wta * x1 + wtb * x2;
    xdis[nx3 - 1] = x;
    sigdis[nx3 + (2*PVMORPH_MAXBINS + 1)] = 0.;
    if (idebug >= 1) {
      printf("First,last edges: (%d , %d) (%d , %d)\n",ix1,ix1l,ix2,ix2l);
      printf("First interpolated point: (%lf , %lf)\n",xdis[0],sigdis[(2*PVMORPH_MAXBINS + 2)]);
/* Computing MAX */
	i__2 = ix1l + 1, i__3 = ix2l + 1;
	i__1 = max(i__2,i__3);
	for (i__ = min(ix1,ix2); i__ <= i__1; ++i__) {
	  printf(" %d %lf %lf\n",i__,sigdis[i__ - 1],sigdis[i__ + PVMORPH_MAXBINS]);
	}
    }

/* ......Loop over the remaining point in both curves. Getting the last */
/*      points is a bit tricky due to limited floating point precision. */

    if (idebug >= 1) {
      printf("----------Double loop------------\n");
    }
    yprev = -1.;
    while(ix1 < ix1l || ix2 < ix2l) {

/* ......Increment to the next lowest point. Step up to the next */
/*      kink in case there are several empty (flat in the integral) */
/*      bins. */

	if ((sigdis[ix1] < sigdis[ix2 + (PVMORPH_MAXBINS + 1)] || ix2 == ix2l) && ix1 < ix1l) {
/* * 03.01.05         IF((sigdis(ix1+1,1).LE.sigdis(ix2+1,2).OR.ix2.EQ.ix2l) */
	    ++ix1;
/* * 03.01.05           DO WHILE(sigdis(ix1+1,1).LE.sigdis(ix1,1) */
	    while(sigdis[ix1] < sigdis[ix1 - 1] && ix1 < ix1l) {
		++ix1;
	    }
	    i12 = 1;
	} else if (ix2 < ix2l) {
	    ++ix2;
/* * 03.01.05           DO WHILE(sigdis(ix2+1,2).LE.sigdis(ix2,2) */
	    while(sigdis[ix2 + (PVMORPH_MAXBINS + 1)] < sigdis[ix2 + PVMORPH_MAXBINS] && ix2 < ix2l) {
		++ix2;
	    }
	    i12 = 2;
	}
	if (idebug >= 1) {
	  printf("Pair %d %d %d\n",ix1,ix2,i12);
	}
	if (i12 == 1) {
	    if (idebug >= 1) {
	      printf("   %lf %lf %lf\n",sigdis[ix2 + PVMORPH_MAXBINS],sigdis[ix1 - 1],sigdis[ix2 + (PVMORPH_MAXBINS + 1)]);
	    }
	    x1 = xmin1d + (double) (ix1 - 1) * dx1;
	    y = sigdis[ix1 - 1];
	    x20 = (double) (ix2 - 1) * dx2 + xmin2d;
	    x21 = x20 + dx2;
	    y20 = sigdis[ix2 + PVMORPH_MAXBINS];
	    y21 = sigdis[ix2 + (PVMORPH_MAXBINS + 1)];
	    if (y21 > y20) {
		x2 = x20 + (x21 - x20) * (y - y20) / (y21 - y20);
	    } else {
		x2 = x20;
	    }
	} else {
	    if (idebug >= 1) {
	      printf("   %lf %lf %lf\n",sigdis[ix1 - 1],sigdis[ix2 + PVMORPH_MAXBINS],sigdis[ix1]);
	    }
	    x2 = xmin2d + (double) (ix2 - 1) * dx2;
	    y = sigdis[ix2 + PVMORPH_MAXBINS];
	    x10 = (double) (ix1 - 1) * dx1 + xmin1d;
	    x11 = x10 + dx1;
	    y10 = sigdis[ix1 - 1];
	    y11 = sigdis[ix1];
	    if (y11 > y10) {
		x1 = x10 + (x11 - x10) * (y - y10) / (y11 - y10);
	    } else {
		x1 = x10;
	    }
	}
	x = wta * x1 + wtb * x2;
/* * 03.01.05        IF(y.GT.yprev) THEN */
	if (y >= yprev) {
	    if (nx3 < PVMORPH_MAXBINS) {
		++nx3;
	    } else {
	      printf("CSM_PVMORPH - error - Interpolation\n");
	      printf("\aborted, too many points (max 5k).n");
		*nbn = 0;
		*xminn = 0.;
		*xmaxn = 0.;
		goto L999;
	    }
	    if (idebug >= 1) {
              d__1 = 1. - y;
	      printf("                    %d %d %lf %lf\n",i12,nx3,x,d__1);
	    }
	    yprev = y;
	    xdis[nx3 - 1] = x;
	    sigdis[nx3 + (2*PVMORPH_MAXBINS + 1)] = y;
	    if (idebug >= 1) {
                d__1 = 1. - sigdis[ix1 - 1];
                d__2 = 1. - sigdis[ix2 + PVMORPH_MAXBINS];
                d__3 = 1. - sigdis[nx3 + (2*PVMORPH_MAXBINS + 1)];
	        printf("aaa %d %d %d %lf %lf %d %lf %lf\n",ix1,ix2,i12,d__1,d__2,nx3,x,d__3);
	    }
	}
    }

/* ......Now we loop over the edges of the bins of the interpolated */
/*      histogram and find out where the interpolated cdf 3 */
/*      crosses them. This defines the result. */

/* ......We set all the bins following the final edge to the value */
/*      of the final edge. */

    ix = *nbn + 1;
    x = xmind + (double) (ix - 1) * dx;
    sigdis[ix + (3*PVMORPH_MAXBINS + 2)] = sigdis[nx3 + (2*PVMORPH_MAXBINS + 1)];
/* Always endpoint, even if GE fails */
    while(x >= xdis[nx3 - 1]) {
	sigdis[ix + (3*PVMORPH_MAXBINS + 2)] = sigdis[nx3 + (2*PVMORPH_MAXBINS + 1)];
	if (idebug >= 1) {
	  printf("Setting additional final bins %d %lf %lf\n",ix,x,sigdis[ix + (3*PVMORPH_MAXBINS + 2)]);
	}
	--ix;
	x = xmind + (double) (ix - 1) * dx;
    }
    ixl = ix;

/* ......The beginning may be empty, so we have to step up to the first */
/*      edge where the result is nonzero. */

    ix = 1;
    x = xmind + (double) (ix - 1) * dx;
    sigdis[ix + (3*PVMORPH_MAXBINS + 2)] = sigdis[(2*PVMORPH_MAXBINS + 2)];
/* Always endpoint, even if LE fails */
    while(x <= xdis[0]) {
	sigdis[ix + (3*PVMORPH_MAXBINS + 2)] = sigdis[(2*PVMORPH_MAXBINS + 2)];
	if (idebug >= 1) {
	  printf("Setting additional initial bins %d %lf %lf\n",ix,x,sigdis[ix + (3*PVMORPH_MAXBINS + 2)]);
	}
	++ix;
	x = xmind + (double) (ix - 1) * dx;
    }
    ixf = ix;
    if (idebug >= 1) {
      printf("Bins left to loop over: %d - %d\n",ixf,ixl);
    }

/* ......Also the end (y~1.0) often comes */
/*      before the last bin so we have to set the following to 1.0 as well. */

    if (idebug >= 1) {
      printf("-------------------------------\n");
      printf("ix3,x3,sigdis3\n");
	i__1 = nx3;
	for (ix3 = 1; ix3 <= i__1; ++ix3) {
	  printf("%d %lf %lf %lf %lf\n",ix3,xdis[ix3 - 1],sigdis[ix3 - 1],sigdis[ix3 + PVMORPH_MAXBINS],sigdis[ix3 + (2*PVMORPH_MAXBINS + 1)]);
	}
      printf("-------------------------------\n");
    }
    ix3 = 1;
    i__1 = ixl;
    for (ix = ixf; ix <= i__1; ++ix) {
	x = xmind + (double) (ix - 1) * dx;
	if (x < xdis[0]) {
	    y = 0.;
	} else if (x > xdis[nx3 - 1]) {
	    y = 1.;
	} else {
	    while(xdis[ix3] <= x && ix3 < *nbn << 1) {
		++ix3;
	    }
/* ......ALR 03-JAN-2005 As far as I can tell this empty bin treatment */
/*                      is not relevant. If the interpolated cdf has been */
/*                      formed properly it shouldn't matter if 2 points */
/*                      on cdf-3 straddle more than a single bin edge on the */
/*                      x-axis. If this happens it probably means a hole in */
/*                      the distribution does not match perfectly a bin and the */
/*                      empty bin is "smeared" across 2 bins in the final pdf. */

/* $$$            IF(xdis(ix3+1)-x.GT.1.1*dx2) THEN ! Empty bin treatment */
/* $$$               y = sigdis(ix3+1,3) */
/* $$$            ELSE IF(xdis(ix3+1).GT.xdis(ix3)) THEN ! Normal bins */
	    if (xdis[ix3] > xdis[ix3 - 1]) {
/* Normal bins */
		y = sigdis[ix3 + (2*PVMORPH_MAXBINS + 1)] + (sigdis[ix3 + (2*PVMORPH_MAXBINS + 2)] - sigdis[ix3 + 
			(2*PVMORPH_MAXBINS + 1)]) * (x - xdis[ix3 - 1]) / (xdis[ix3] - xdis[ix3 
			- 1]);
	    } else {
/* This should never happen */
		y = 0.;
		// printf("csm_pvmorph - ERROR: infinite slope in cdf-3\n");
	    }
	}
	sigdis[ix + (3*PVMORPH_MAXBINS + 2)] = y;
	if (idebug >= 1) {
          d__1 = 1. - sigdis[ix3 + (2*PVMORPH_MAXBINS + 1)];
          d__2 = 1. - y;
	  d__3 = 1. - sigdis[ix3 + (2*PVMORPH_MAXBINS + 2)];

	  printf("%d %d %lf %lf %lf %lf %lf %lf\n",ix,ix3,xdis[ix3 - 1],x,xdis[ix3],d__1,d__2,d__3);
	}
/* L6600: */
    }

/* ......Differentiate to get pdf from cdf and place in single precision */
/*      output array. */
    for (ix = *nbn + 1; ix >= 2; --ix) {
	y = sigdis[ix + (3*PVMORPH_MAXBINS + 2)] - sigdis[ix + (3*PVMORPH_MAXBINS + 1)];
	distn[ix - 1] = y;
    }

L999:
    ;
} /* csm_pvmorph */

/*--------------------------------------------------------------------------------------------*/

/* csm_pvmoprh_2d:  original C'd from Alex Read's d_pvmorph_2d in d_higgs.f using f2c
  -- trj 8 Nov 2005 -- removed f2c's ugly i/o and replace with
  cout's.  There are still a lot of hardwired buffer sizes here.*/


/* .============================================================== */
void csm_acnvec(Double_t* vec, Int_t* n)
{
    /* System generated locals */
    Int_t i__1;

    /* Local variables */
    static Int_t i__;

/* --------------------------------------------------------------- */
/* Integrate the probability distribution in array vec and normalize */
/* the integral to 1.0 */
/* --------------------------------------------------------------- */


/* ......Integrate the probability distribution. */

    /* Parameter adjustments */
    --vec;

    /* Function Body */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	vec[i__] += vec[i__ - 1];
    }

    if (vec[*n] == (float)0.) {
/*         WRITE (6,*) 'csm_acnvec: total is zero' */
	return;
    }

/* ......Normalize the cummulative probability distribution to 1. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vec[i__] /= vec[*n];
    }

} /* csm_acnvec */

/* .================================================================== */
void csm_pvmorph_2d(Int_t* nx1, Double_t* xmin1, Double_t* xmax1,
                    Int_t* ny1, Double_t* ymin1, Double_t* ymax1, 
                    Double_t* xydist1, 
                    Int_t* nx2, Double_t* xmin2, Double_t* xmax2, 
                    Int_t* ny2, Double_t* ymin2, Double_t* ymax2,
                    Double_t* xydist2, 
                    Int_t* nx3, Double_t* xmin3, Double_t* xmax3,
                    Int_t* ny3, Double_t* ymin3, Double_t* ymax3,
                    Double_t* xydist3, 
                    Double_t* par1, Double_t* par2, Double_t* par3)
{
    /* System generated locals */
    Int_t i__1, i__2, i__3;

    /* Local variables */
    static Double_t alpha[45000]	/* was [150][150][2] */, xtemp[450]	/* 
	    was [150][3] */;
    static Double_t ydist1[150], ydist2[150], ydist3[150];
    static Int_t ix, iy;
    static Int_t iyc;


/* -------------------------------------------------------------------------- */
/* Author           : Alex Read 06-NOV-1998 */
/* CVS identification: $Id: d_pvsmorph.f,v 1.3 1998/11/06 11:22:59 read Exp $ */
/*                     $Source: /mn/axuoep1/u1/read/CVSROOT/methodA/d_pvsmorph.f,v $ */

/* ......Do a linear interpolation between two two-dimensional */
/*      probability distributions (scatterplots) as a function */
/*      of the characteristic parameter of the distribution. */

/*      This is a generalization of d_pvmorph. The 2d distribution */
/*      can be move around and be stretched or squeezed in two */
/*      dimenions but finite rotations (changes in the correlation) */
/*      are poorly approximated. */

/* Takes two input scatterplots (1 and 2) in the form of 1d arrays and */
/* binning constants and compute the linearly interpolated or morphed */
/* distribution (3) and its binning constants. */

/* nx1,nx2,nx3        : Number of x-bins in the input and output distributions. */
/* ny1,ny2,ny3        : Number of y-bins in the input and output distributions. */
/* xmin1,xmin2,xmin3  : Lower x-edges of histograms */
/* xmax1,xmax2,xmax3  : Upper x-edges of histograms */
/* ymin1,ymin2,ymin3  : Lower y-edges of histograms */
/* ymax1,ymax2,ymax3  : Upper y-edges of histograms */
/* xydist1,xydist2,xydist3  : Bin contents of scatterplots. The arrays should be */
/*                            packed with the index running fastest over the x */
/*                            dimension. */
/* par1,par2,par3     : Values of the linear parameters that characterise the */
/*                      histograms (e.g. the Higgs mass). */

/* Inputs: par1,par2,par3,nx1,nx2,ny1,ny2,xmin1,xmin2,ymin1,ymin2 */
/*         xydist1,xydist2 */
/* Output: nx3,ny3,xmin3,xmax3,ymin3,ymax3,xydist3 */

/* The scatterplots are limited in size to 110x110 bins, this is the limit */
/* of temporary arrays used to manipulate copies of the input scatterplots. */
/* ------------------------------------------------------------------------ */


/* ----------------------------------------------------------------------- */
/*  18.06.98 - Generalize the interpolation up to 100x100 bins */
/*  27.09.01 - Increase number of bins to 110x110 */
/*  01.07.02 - Increase number of bins to 150x150 */
/* ----------------------------------------------------------------------- */

/* ......Arguments */


/* ......Local variables */


/* ......Check that the interpolation can be performed (01.07.02) */

    /* Parameter adjustments */
    --xydist3;
    --xydist2;
    --xydist1;

    /* Function Body */
    if (max(*ny1,*ny2) > 150) {
	i__1 = max(*ny1,*ny2);
        cout << "*** Fatal error in csm_pvmorph_2d - " << endl;
	cout << "exceeded limit of " << 150 << " bins on y-axis: " << i__1 << endl;
	exit(0);
    }

    if (max(*nx1,*nx2) > 150) {
	i__1 = max(*nx1,*nx2);
        cout << "*** Fatal error in csm_pvmorph_2d - " << endl;
	cout << "exceeded limit of " << 150 << " bins on x-axis: " << i__1 << endl;
	exit(0);
    }

/* ......Extract the probability distributions projected onto */
/*      the y-axis. */

    csm_ypvscat(ydist1, &xydist1[1], nx1, ny1);

    csm_ypvscat(ydist2, &xydist2[1], nx2, ny2);

/* ......Interpolate the y-projections of the 2d probability */
/*      distribution. */

    *ny3 = 150;
/* input array size, output result size */
    csm_pvmorph(ny1, ymin1, ymax1, ydist1, ny2, ymin2, ymax2, ydist2, ny3, 
	    ymin3, ymax3, ydist3, par1, par2, par3);

/* ......Find out which y bins of the inputs (1,2) contribute to the */
/*      interpolated y-projection (3). */

    csm_getycont(ydist1, ny1, ydist2, ny2, ydist3, ny3, alpha);

/* ......Extract the x-distributions in the y-slice determined above */
/*      and interpolate them. */

    i__1 = *ny3;
    for (iy = 1; iy <= i__1; ++iy) {
/* Loop over resulting bins */
      for (int ivz=0; ivz<450; ivz++) { xtemp[ivz] = 0; }
	i__2 = *ny1;
	for (iyc = 1; iyc <= i__2; ++iyc) {
/* Loop over contributing bins */
	    i__3 = *nx1;
	    for (ix = 1; ix <= i__3; ++ix) {
/* Extract contributions to x-distributions */
		xtemp[ix - 1] += alpha[iyc + (iy + 150) * 150 - 22651] * 
			xydist1[ix + (iyc - 1) * *nx1];
		xtemp[ix + 149] += alpha[iyc + (iy + 300) * 150 - 22651] * 
			xydist2[ix + (iyc - 1) * *nx2];
	    }
	}

/* ......Interpolate the x-distributions. */

	*nx3 = 150;
	csm_pvmorph(nx1, xmin1, xmax1, xtemp, nx2, xmin2, xmax2, &xtemp[150], 
		nx3, xmin3, xmax3, &xtemp[300], par1, par2, par3);

/* ......Insert the interpolated x-distribution into the final scatterplot */
/*      array. */

	i__2 = *nx3;
	for (ix = 1; ix <= i__2; ++ix) {
	    xydist3[ix + (iy - 1) * *nx3] = xtemp[ix + 299] * ydist3[iy - 1];
	}
    }

} /* csm_pvmorph_2d */

/* .================================================================== */
void csm_ypvscat(Double_t* ydist, Double_t* xydist, Int_t* nx, Int_t* ny)
{
    /* System generated locals */
    Int_t xydist_dim1, xydist_offset, i__1, i__2;

    /* Local variables */
    static Double_t total;
    static Int_t ix, iy;

/* :------------------------------------------------------------------ */
/* Author: Alex Read    06-NOV-1998 */
/* CVS identification: */
/* $Id: d_pvsmorph.f,v 1.3 1998/11/06 11:22:59 read Exp $ */
/* $Source: /mn/axuoep1/u1/read/CVSROOT/methodA/d_pvsmorph.f,v $ */

/* Project a scatterplot onto the y-axis (called from d_pvsmorph).The */
/* projection is normalized so that the sum of the bin contents is 1.0. */

/* nx,ny    : Number of bins in the scatterplot for the x and y coordindates. */
/*            The projection is done onto <ny> bins. */
/* xydist   : The 2-dimensional array of the probabilities */
/* ydist    : The 1-dimensional array of the 2d probabilities projected onto */
/*            the y-axis. */

/* Inputs : nx,ny,xydist */
/* Outputs: ny,ydist (ny is unchanged from input to output) */
/* .------------------------------------------------------------------- */

/* ......Arguments */


/* ......Local variables */


/* ......Project onto y-axis and keep track of total */

    /* Parameter adjustments */
    xydist_dim1 = *nx;
    xydist_offset = xydist_dim1 + 1;
    xydist -= xydist_offset;
    --ydist;

    /* Function Body */
    for (int ivz=1; ivz<=*ny; ivz++) { ydist[ivz] = 0; }
    total = (float)0.;
    i__1 = *ny;
    for (iy = 1; iy <= i__1; ++iy) {
	i__2 = *nx;
	for (ix = 1; ix <= i__2; ++ix) {
	    ydist[iy] += xydist[ix + iy * xydist_dim1];
	}
	total += ydist[iy];
    }

/* ......Normalize y-projection to 1.0 */

    i__1 = *ny;
    for (iy = 1; iy <= i__1; ++iy) {
	ydist[iy] /= total;
    }

/* L999: */
} /* csm_ypvscat */


/* .============================================================== */
void csm_getycont(Double_t* ydist1, Int_t* ny1, 
                  Double_t* ydist2, Int_t* ny2, 
                  Double_t* ydist3, Int_t* ny3, 
	          Double_t* alpha)
{
    /* System generated locals */
    Int_t i__1, i__2, i__3, i__4, i__5;
    Double_t r__1, r__2;

    /* Local variables */
    static Double_t xmin, xmax;
    static Int_t ibin0, ibin1, i__, iedge;
    static Double_t x[151], y[453]	/* was [151][3] */;//, xedge;
    static Int_t ihist;
    static Double_t x0, x1, x2, dx;
    static Int_t ix;
    static Double_t xl, xr, sigdis[450]	/* was [150][3] */, xcontr[302]	/* 
	    was [151][2] */;
    static Int_t nch;


/* :------------------------------------------------------------------ */
/* Author: Alex Read    06-NOV-1998 */
/* CVS identification: */
/* $Id: d_higgs.f,v 1.4 1998/12/22 10:58:21 read Exp $ */
/* $Source: */
/* /afs/cern.ch/delphi/tasks/lep200/searches/limits/alrmc/CVSROOT */
/*  ..../methodA/delphi/d_higgs.f,v $ */

/* Bugfix 11-MAY-1999 (may have affected distributions with holes) */

/* This routine is called by <d_psvmorph>. */

/* <ydist1> and <ydist2> are the projections on the y-axis of two */
/* scatterplots which are going to be interpolated. <ydist3> is */
/* the interpolated 1d distribution which represent the projection */
/* of the interpolated scatterplot on the y-axis. This routine determines */
/* which bins of <ydist1> and <ydist2> contribute and by what amount to */
/* each bin of <ydist3>. This information is used in d_psvmorph to */
/* determine the input distributions in the x-direction of each */
/* y-bin: these are then interpolated and accumulated in the interpolated */
/* 2d distribution. */

/* Inputs : ny1,ny2,ny3,ydist1,ydist2,ydist3 */
/* Outputs: alpha */
/* .---------------------------------------------------------------- */

/* ......Arguments */


/* ......Local variables */


/* ......Verify that input arguments are acceptable. */

    /* Parameter adjustments */
    --ydist1;
    --ydist2;
    --ydist3;
    alpha -= 22651;

    /* Function Body */
    if (*ny1 != *ny2 || *ny1 != *ny3 || *ny1 > 150) {
      cout << "d_getycont: Input arrays must all have same length" << endl;
      cout << "            and size LE " << 150 << ": " << 
            *ny1 << " " << *ny2 << " " << *ny3 << endl;
      cout << "            Rest of processing skipped!" << endl;
    }

/* ......Copy input arrays to temporary storage and normalize them to 1.0 */

    for (int iuc=1;iuc<=*ny1;iuc++) { sigdis[iuc-1] = ydist1[iuc]; }
    for (int iuc=1;iuc<=*ny2;iuc++) { sigdis[149+iuc] = ydist2[iuc]; }
    for (int iuc=1;iuc<=*ny3;iuc++) { sigdis[299+iuc] = ydist3[iuc]; }
    /* ucopy_(&ydist1[1], sigdis, ny1);
       ucopy_(&ydist2[1], &sigdis[150], ny2);
       ucopy_(&ydist3[1], &sigdis[300], ny3); */

    csm_acnvec(sigdis, ny1);
    csm_acnvec(&sigdis[150], ny2);
    csm_acnvec(&sigdis[300], ny2);
/* $$$      DO i=1,ny1 */
/* $$$         print *,'getycont',i,(sigdis(i,j),j=1,3) */
/* $$$      ENDDO */

/* ......Make arrays to describe straight-line approximation */
/*      of the cummulative distribution. The x-axis is fictive. */

    xmin = (float)0.;
    xmax = (float)1.;
    dx = (float)1. / (Double_t) (*ny1);
    x[0] = xmin;
    y[0] = (float)0.;
    y[151] = (float)0.;
    y[302] = (float)0.;
    i__1 = *ny1;
    for (ix = 1; ix <= i__1; ++ix) {
	x[ix] = (Double_t) ix * dx;
	y[ix] = sigdis[ix - 1];
	y[ix + 151] = sigdis[ix + 149];
	y[ix + 302] = sigdis[ix + 299];
/* $$$         print *,'xy',ix,x(ix),x(ix+1),y(ix,2),y(ix+1,2) */
/* $$$     &        ,y(ix+1,2)-y(ix,2) */
    }

/* ......Find out where the edges of the interpolated distribution */
/*     cross the input distributions. Start with a brute-force */
/*     approach! */


/* ......Loop over internal edges to find out where the contributions */
/*     to the morphed histogram came from. */

/* ......Special treatment for empty bins! */

    xcontr[0] = xmin;
    xcontr[151] = xmin;
    xcontr[*ny1] = xmax;
    xcontr[*ny2 + 151] = xmax;
    i__1 = *ny1;
    for (iedge = 2; iedge <= i__1; ++iedge) {

/* ......Bug fixed 11-MAY-1999 by A. Read. The IF statement is meant to */
/*      detect when the cummulative distribution is constant, i.e., the */
/*      differential distribution has a bin with zero entries. In this */
/*      case the contributions from the two input distributions to the */
/*      interpolated distribution is zero. */

/*      More protection. Distributions varying over many orders of magnitude */
/*      give numerical problems when an infinitessimally narrow region makes */
/*      an inifinitessimally small contribution. The g77 compiler doesn't */
/*      signal a FP exception for the inevitable divide-by-zero but the */
/*      DIGITAL UNIX FORTRAN compiler does. Add protection against this */
/*      on 24-OCT-2000 by A. Read. */

/* $$$         IF(y(iedge,3).GT.y(iedge-2,3)) THEN ! bugfix 11-may-1999 */
	if (y[iedge + 301] > y[iedge + 300]) {
/* by A. Read */
	  //xedge = (Double_t) (iedge - 1) * dx;
	    ix = 1;
	    while(y[ix - 1] < y[iedge + 301] && ix < *ny1) {
		++ix;
	    }
/* Computing MAX */
	    r__1 = (float)1e-30, r__2 = y[ix - 1] - y[ix - 2];
	    x1 = x[ix - 2] + (y[iedge + 301] - y[ix - 2]) * dx / max(r__1,
		    r__2);
/* protection 24.10.00 A */
	    xcontr[iedge - 1] = x1;

	    ix = 1;
	    while(y[ix + 150] < y[iedge + 301] && ix < *ny1) {
		++ix;
	    }
/* Computing MAX */
	    r__1 = (float)1e-30, r__2 = y[ix + 150] - y[ix + 149];
	    x2 = x[ix - 2] + (y[iedge + 301] - y[ix + 149]) * dx / max(r__1,
		    r__2);
/* protection 24.10.00 A */
	    xcontr[iedge + 150] = x2;
	} else {
	    xcontr[iedge - 1] = xcontr[iedge - 2];
/* 0. ! bugfix 11-MAY-199 */
	    xcontr[iedge + 150] = xcontr[iedge + 149];
/* 0. ! bugfix 11-MAY-199 */
	}
    }

/* ......Loop over the bins and construct an array which tells which */
/*      fraction of each bin contributes. */

    for (int ivz=22651; ivz < 22651+45000; ivz++)
      { alpha[ivz] = 0; }
    /* f2c made this after adjusting the origin of alpha's index : vzero_(&alpha[22651], &45000); */
    for (ihist = 1; ihist <= 2; ++ihist) {
	i__1 = *ny1;
	for (iedge = 1; iedge <= i__1; ++iedge) {
	    nch = *ny1;
/* temp */
	    x0 = xcontr[iedge + ihist * 151 - 152];
	    x1 = xcontr[iedge + 1 + ihist * 151 - 152];
	    dx = (xmax - xmin) / (Double_t) nch;
	    if (dx > (float)0.) {
/* Computing MAX */
/* Computing MIN */
		i__4 = nch, i__5 = (Int_t) ((x0 - xmin) / dx) + 1;
		i__2 = 1, i__3 = min(i__4,i__5);
		ibin0 = max(i__2,i__3);
/* Computing MAX */
/* Computing MIN */
		i__4 = nch, i__5 = (Int_t) ((x1 - xmin) / dx) + 1;
		i__2 = 1, i__3 = min(i__4,i__5);
		ibin1 = max(i__2,i__3);

/* ......Mark the fraction of the single bin that contributes */

		if (ibin0 == ibin1) {
/* both edges inside same bin */
		    alpha[ibin0 + (iedge + ihist * 150) * 150] = (x1 - x0) / 
			    dx;
		} else {

/* ......Mark all the complete bins with contribution 1.0 */

		    i__ = ibin0 + 1;
		    while(i__ <= nch && i__ < ibin1) {
			alpha[i__ + (iedge + ihist * 150) * 150] = (float)1.;
			++i__;
		    }

/* ......Mark the fraction of the first partial bin. */

		    xr = xmin + dx * (Double_t) ibin0;
		    alpha[ibin0 + (iedge + ihist * 150) * 150] = (xr - x0) / 
			    dx;

/* ......Mark the fraction of the last partial bin. */

		    xl = xmin + dx * (Double_t) (ibin1 - 1);
		    alpha[ibin1 + (iedge + ihist * 150) * 150] = (x1 - xl) / 
			    dx;
		}
/* L6500: */
	    }
        }
    }

/* ......At this point we have in hand all the elements we need to */
/*      begin interpolating the two scatterplots. Using the */
/*      contribution coefficients above we need to slice out */
/*      the 1d probability distributions from the scatterplots, */
/*      interpolate them, and accumulate the resulting distribution */
/*      into the morphed scatterplot. */

} /* csm_getycont */


/*-------------------------------------------------------------------------*/
/*    Interface to Joel's genlimit.c program for a Bayesian calcualtion    */
/*    of an upper limit.  Also run pseudoexperiments if need be to compute */
/*    expected limits.                                                     */
/*    Bayesian limit calculation uses test_hypothesis_pe to compute the    */
/*    "Bayesian ensemble" because it should have signal and background     */
/*    components marked, and because it should have all systematic         */
/*    errors included.  It uses null_hypothesis_pe in order to generate    */
/*    pseudoexperiments to compute expected limits however.                */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
               sflimit:  the observed limit
               unc:      MC statistical unc. on observed limit
               npx:      Number of pseudoexperiments to run to compute expected limits
               sm2:      -2 sigma expected limit      *put in null pointers for all
               sm1:      -1 sigma expected limit      *five of these to skip the
               smed:     median expected limit        *calculation and speed it up.
               sp1:      +1 sigma expected limit
               sp2:      +2 sigma expected limit
*/

void mclimit_csm::bayes_heinrich_withexpect(Double_t beta,
                                            Double_t* sflimit,
                                            Double_t* unc,
					    Int_t npx,
                                            Double_t* sm2,
                                            Double_t* sm1,
                                            Double_t* smed,
                                            Double_t* sp1,
                                            Double_t* sp2)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<500) {nglmax = 500;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;
  *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);

  // compute expected limits

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      pdname = new char[strlen(test_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) null_hypothesis_pe->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      null_hypothesis_pe->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
                   if (nbinsy==1)
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
	           else
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		   nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
	           ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
       nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
       delete[] xgl;
       delete[] lwgl;
       xgl = new double[nglmax];
       lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++)
    {
      if (ipx>0)
	{ if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]])
	  { 
	    ngl = 0;
	  }
	}
      cslist.push_back((Double_t) cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,unc));
    }
  std::sort(cslist.begin(),cslist.end());
  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm2 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *smed = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp2 = cslist[i];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
}

/*-------------------------------------------------------------------------*/
//  Same thing as above, but only compute observed limit (much quicker)
/*-------------------------------------------------------------------------*/

void mclimit_csm::bayes_heinrich(Double_t beta,
                                 Double_t* sflimit,
                                 Double_t* unc)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.

  nglmax = nobstot;
  if (nglmax<500) {nglmax = 500;}
  double* xgl = new double[nglmax];
  double* lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;
  *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);

  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
}

/*-------------------------------------------------------------------------*/
/* Coverage checker for Joel's Bayesian limit calc.  Based on              */
/* bayes_heinrich_withexpect, but now uses test_hypothesis_pe scaled       */
/* so the signal is at the observed 95% exclusion rate.  The px's are done */
/* assuming the signal+background is present, and the false exclusion rate */
/* is computed.                                                            */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
               sflimit:  the observed limit
               unc:      MC statistical unc. on observed limit
               npx:      Number of pseudoexperiments to run to compute fales exclusion rate
	       falsex:   false exclusion rate:  Should be no more than 1-beta.

*/

void mclimit_csm::bayes_heinrich_coverage_check(Double_t beta,
                                                Double_t* sflimit,
                                                Double_t* unc,
					        Int_t npx,
                                                Double_t* falsex)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=flat;
  //  const PRIOR prior=exp;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<500) {nglmax = 500;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;
  *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);

  // Run signal+background pseudoexperiments at the observed limit, and see what
  // the distribution of limits we get out is.  Limits are computed using the
  // same Bayesian ensemble with the unscaled test_hypothesis_pe and so the
  // limits that are more restrictive than *sflimit are false exclusions.

  csm_model* testhyppescale = test_hypothesis_pe->scalesignal(*sflimit);

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      pdname = new char[strlen(testhyppescale->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,testhyppescale->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) testhyppescale->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      testhyppescale->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
                   if (nbinsy==1)
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
	           else
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		   nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
	           ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
       nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
       delete[] xgl;
       delete[] lwgl;
       xgl = new double[nglmax];
       lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++)
    {
      if (ipx>0)
	{ if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]])
	  { 
	    ngl = 0;
	  }
	}
      cslist.push_back((Double_t) cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,unc));
    }

  Int_t nfalse = 0;
  for (i = 0; i<npx; i++)
    {
      if (cslist[i] <= *sflimit)
	{
	  nfalse++;
	}
    }

  *falsex = ((Double_t) nfalse)/((Double_t) npx);

  // clean up

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
  delete testhyppescale;
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*    genlimit Bayesian code from Joel                                     */
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

/*  Joel Heinrich  8 April 2005

Returns cross section posterior p.d.f. evaluated at s.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/


double cspdf(double s,double norm,
	     int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior) {
  int i,k;
  double sum = 0, lgp = 0;
  const EB* p = ens;

  assert(nens>0);
  assert(prior==flat || prior==corr );// || prior==exp);
  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*p->e + p->b;
      const int n = nobs[k];
      esum += p->e;
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    sum += (prior==flat) ? exp(t) : esum*exp(t);
  }
  //cout << "pdf: " << s << " " << sum/nens << endl;
  return sum/(norm*nens);
}

/*  Joel Heinrich  8 April 2005

Function returns integral from xlo to infinity.

*uncertainty (if not null pointer) returned with uncertainty of
integral due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/


double csint(double xlo,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty) {
  int i,ntot=0;
  double sum=0, sum2=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    const double t = csint0(xlo,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior);
    sum += t;
    sum2 += t*t;
    p+=nchan;
  }
  sum /= nens;
  sum2 /= nens;
  if(uncertainty)
    *uncertainty = (nens>1) ? sqrt((sum2-sum*sum)/(nens-1)) : 1.0 ;
  return sum;
}

double csint0(double xlo,double logscale,
		     int nchan,const int nobs[],const EB chan[],
		     int ngl,const double xgl[],const double lwgl[],
		     PRIOR prior) {
  int i,k;
  double sum=0, esum=0, bsum=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum += chan[i].b + xlo*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo;
    double t = logscale-bsum, v;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum += v = exp(lwgl[k]+t);
    if(v<DBL_EPSILON*sum) break;
  }
  if(prior==flat)
    sum *= resum;
  return sum;
}


/*      Joel Heinrich  8 April 2005

    returns:
      *int1 = integral from xlo1 to xhi
      *int2 = integral from xlo2 to xhi
      *v11  = variance of *int1
      *v12  = covariance between *int1 and *int2
      *v22  = variance of *int2

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/



void csint2cut(double xlo1,double xlo2,double xhi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02cut(xlo1,xlo2,xhi,logscale,nchan,nobs,p,*ngl,xgl,lwgl,
	       prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, sum3=0, esum=0, bsum1=0, bsum2=0, bsum3=0, resum;

  for(i=0;i<nchan;++i) {
    const double ee=chan[i].e, bb=chan[i].b;
    esum += ee;
    bsum1 += bb + xlo1*ee;
    bsum2 += bb + xlo2*ee;
    bsum3 += bb + xhi*ee;
  }

  if(esum==0 && prior==corr) {
    *int1 = *int2 = 0;
    return;
  }

  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2);
    if(v2<DBL_EPSILON*sum2) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xhi;
    double t3 = logscale-bsum3, v3;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t3 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum3 += v3 = exp(lwgl[k]+t3);
    if(v3<DBL_EPSILON*sum3) break;
  }

  if(prior==flat) {
    *int1 = (sum1-sum3)*resum;
    *int2 = (sum2-sum3)*resum;
  } else {
    *int1 = sum1-sum3;
    *int2 = sum2-sum3;
  }
  return;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/


double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty) {
  const double eps=1.0e-6;

  /*
  cout << "nchan: " << nchan << endl;
  cout << " nens: " << nens << endl;
  cout << " prior: " << prior << endl;
  int i;
  for (i=0;i<nchan;i++)
    {
      cout << "nobs(" << i << ") = " << nobs[i] << endl;
    }
  */

  double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  double limit = galim(beta,nchan,nens,nobs,ens);
  double dl=limit, rpdf=0;
  double lo=0, hi=1e200;

  while(fabs(dl)>1.0e-10*limit) {
    const double pbeta =
      1-csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)/norm;
    rpdf = 1/cspdf(limit,norm,nchan,nens,nobs,ens,prior);
    if(pbeta>beta) {
      hi=limit*(1+eps);
    } else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }
  }

  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;
    csint2(0,limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	   &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/norm   ;
  }


  return limit;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/




double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty) {
  const double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double tail = csint(smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double eps=1.0e-6;
  double limit = galim(beta,nchan,nens,nobs,ens), rpdf=0;
  double dl=limit, lo=0, hi=smax;

  if(beta<=0) return 0;
  if(beta>=1) return smax;

  if(limit>smax || limit<0) dl = limit = 0.5*smax;

  while(fabs(dl)>1.0e-10*limit) {
    double pbeta =
      1-(csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)-tail)/
      (norm-tail);
    rpdf = 1/cspdf(limit,norm-tail,nchan,nens,nobs,ens,prior);
    if(pbeta>beta) {
      hi=limit*(1+eps);
    } else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }

  }


  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;

    csint2cut(0,limit,smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	      &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/(norm-tail)   ;
  }


  return limit;
}


/*         Joel Heinrich  24 March 2005

   returns crude Gaussian approximation to upper limit.
   For use as a starting point.
*/

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens) {
  double mean=0,sigma=0;
  gameansigma(&mean,&sigma,nchan,nens,nobs,ens);
  return mean-sigma*arcfreq( (1-freq(-mean/sigma))*(1-beta) );
}



void gameansigma(double *mean,double *sigma,
		 int nchan,int nens,const int nobs[],const EB* ens) {

  double sum=0,sum2=0,vsum=0;
  const EB* p = ens;
  int i,j;
  for(i=0;i<nens;++i) {
    double s=0,s2=0;
    for(j=0;j<nchan;++j) {
      const int n = nobs[j];
      const double eps = p->e;
      s += (n-p->b)*eps/(n+1);
      s2 += eps*eps/(n+1);
      ++p;
    }
    s /= s2;
    vsum += 1/s2;
    sum += s;
    sum2 += s*s;
  }
  
  *mean = sum/nens;
  *sigma = sqrt(vsum/nens + sum2/nens - (*mean)*(*mean));
  return;
}

#define rdfreq(x) (exp(0.5*(x)*(x))*2.50662827463100050242)

#define C0 2.515517
#define C1 0.802853
#define C2 0.010328
#define D0 1.0
#define D1 1.432788
#define D2 0.189269
#define D3 0.001308

double arcfreq(double y) {
  const double yy = (y>0.5) ? 1-y : y, t = sqrt(-2*log(yy));
  double x = (C0+t*(C1+t*C2))/(D0+t*(D1+t*(D2+t*D3))) - t;
  x -= (freq(x) - yy)*rdfreq(x);
  x -= (freq(x) - yy)*rdfreq(x);
  return (y>0.5) ? -x : x;
}


/*

   Joel Heinrich
   February 10 2005

Returns Gauss-Laguerre quadrature abscissas and log(weights) which can
be used to approximate

      integral u=0 to infinity pow(u,alpha)*exp(-u)*f(u) du
as
      sum k=0 to n-1  exp(lw[k])*f(x[k])

or equivalently

      sum k=0 to n-1  exp(lw[k]+log(f(x[k])))

The quadrature is exact for polynomial f of degree 2n-1 or less.

*/


void gausslaguerre(double x[],double lw[],int n,double alpha){
  const int nshift = 20;
  const double shift = 1<<nshift, rshift=1/shift;
  int i;
  double z=0;
  
  for(i=0;i<n;++i) {
    int j=0, k=2, nscale=0;
    double dz=0.0, p1=0, p2=0;
    if(i==0) {
      z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
    } else if(i==1) {
      z += (15.0+6.25*alpha)/(1.0+2.5*n+0.9*alpha);
    } else if(i==2) {
      const double ai=i-1;
      z += ( (1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0+3.5*ai) )*
	(z-x[i-2])/(1.0+0.3*alpha);
    } else if(i==3) {
      z = 3.0*(x[2]-x[1])+x[0];
    } else if(i==4) {
      z = 4.0*x[3] - 6.0*x[2] + 4.0*x[1] - x[0];
    } else if(i==5) {
      z = 5.0*x[4] - 10.0*x[3] + 10.0*x[2] - 5.0*x[1] + x[0];
    } else {
      z = 6.0*x[i-1] - 15.0*x[i-2] + 20.0*x[i-3] -
	15.0*x[i-4] + 6.0*x[i-5] - x[i-6];
    }
    while(k>0) {
      p1=1;
      p2=0;
      nscale=0;
      z -= dz;
      for(j=1;j<=n;++j){
	const double p3=p2;
	p2=p1;
	p1=((2*j-1+alpha-z)*p2 - (j-1+alpha)*p3)/j;
	if(fabs(p2)>shift) {
	  ++nscale;
	  p1 *= rshift;
	  p2 *= rshift;
	}
      }
      dz = p1*z/(n*p1-(n+alpha)*p2);
      if(fabs(dz)<1.0e-10*z)--k;
    }
    x[i]=z;
    lw[i] = log(z/(p2*p2)) - 2*nshift*nscale*M_LN2 ;
  }
  
  {
    double t = 0.0;
    for(i=n-1;i>=0;--i)
      t += exp(lw[i]);
    t = lgamma(alpha+1)-log(t);
    for(i=0;i<n;++i)
      lw[i] += t;
  }

  return;
}

/*    Joel Heinrich  8 April 2005

   returns:
      *int1 = integral from xlo1 to infinity
      *int2 = integral from xlo2 to infinity
      *v11  = variance of *int1
      *v12  = covariance between *int1 and *int2
      *v22  = variance of *int2


See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/



void csint2(double xlo1,double xlo2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02(xlo1,xlo2,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02(double xlo1,double xlo2,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, esum=0, bsum1=0, bsum2=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum1 += chan[i].b + xlo1*chan[i].e;
    bsum2 += chan[i].b + xlo2*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2);
    if(v2<DBL_EPSILON*sum2) break;
  }

  if(prior==flat) {
    sum1 *=resum;
    sum2 *=resum;
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}
