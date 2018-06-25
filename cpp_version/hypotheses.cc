
#include "TFile.h"
#include <TApplication.h>
#include <TGClient.h>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric> //for vector.accumulate()
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "mclimit_csm.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include <stdlib.h>
#include "math.h"
#include "TMath.h"
#include "TArrow.h"
#include <cassert>

using namespace std;

vector<double> read_in(char *infile){

  vector<double> yvals;
  int line;
  std::ifstream fin(infile, ifstream::in);
  while(!fin.eof()){
    line++;
    double xlo,xi,y,dy;
    fin >> xlo >> xi >> y ; //>> dy;
    dy=0;
    yvals.push_back(y);
  }
  return yvals;
}

double factorial(double n){
  return tgamma(n+1);
}

double Poisson(double s, double b, double d){
  return exp(-(s+b))*(pow(s+b,d))/factorial(d);
  
}

int main(int argc, char* argv[]){

  if (argc != 6){ cout << "Enter sigfile bkgfile background_syst signal_syst lumi [fb^-1] as input" << endl; return 1;}

  char* sigfile = (char*)(argv[1]); 
  char* bkgfile = (char*)(argv[2]); 
  double bkgSys=atof(argv[3]);
  double sigSys=atof(argv[4]);
  double lumi= 1e3*atof(argv[5]); //in pb^-1

  //Templates and input Data
  TH1F *hs;
  TH1F *hb;
  TH1F *hd;
  double scale_s=1; 
  double scale_b=1; 

  vector<double> sig_vals, bkg_vals;
  bkg_vals = read_in(bkgfile);
  sig_vals = read_in(sigfile);

  assert(sig_vals.size() == bkg_vals.size());

  bool rate = false; //flag for shape vs cut-and-count
  if (rate) cout << "Only using total-rate information -i.e. counting analysis" << endl;
  
  //use ALL non-empty bins in fit
  if(!rate){
    hs = new TH1F("hs","",sig_vals.size(),0,sig_vals.size()); //S
    hb = new TH1F("hb","",bkg_vals.size(),0,bkg_vals.size()); //B
    hd = new TH1F("hd","",bkg_vals.size(),0,bkg_vals.size());
    for (int i=1; i<=bkg_vals.size(); i++){
      //skip empty bins
      if (bkg_vals[i-1] == 0 || sig_vals[i-1] == 0) continue;

      hb->SetBinContent(i,  bkg_vals[i-1]*scale_b );
      //      hb->SetBinError(i,  sqrt(bkg_vals[i-1]*scale_b) );

      hs->SetBinContent(i,  sig_vals[i-1]*scale_s );
      //      hs->SetBinError(i,  sqrt(sig_vals[i-1]*scale_s) );

      hd->SetBinContent(i, scale_b*bkg_vals[i-1] );  //pseudodata = B for exclusion
      //      hd->SetBinError(i, sqrt(scale_b*bkg_vals[i-1]) ); 
    }

    for (int i=1; i<=bkg_vals.size(); i++){
      //    hd->SetBinContent(i,  scale_s*sig_vals[i-1]+scale_b*bkg_vals[i-1] ); //pseudodata = S+B for discovery
      //    hd->SetBinError(i,  sqrt(scale_s*sig_vals[i-1]+scale_b*bkg_vals[i-1]));
    }
  }
  else {
    //use only integral of distributions: i.e. cut-and-count
    hs = new TH1F("hs","",1,0,1); //S
    hb = new TH1F("hb","",1,0,1); //B
    hd = new TH1F("hd","",1,0,1); //B
    
    double tot_S = std::accumulate(sig_vals.begin(), sig_vals.end(), 0.);
    double tot_B = std::accumulate(bkg_vals.begin(), bkg_vals.end(), 0.);
    hs->SetBinContent(1,  scale_s*tot_S );
    hs->SetBinError(1,  sqrt(scale_s*tot_S) );
    hb->SetBinContent(1,  scale_b*tot_B );
    hb->SetBinError(1,  sqrt(scale_b*tot_B) );
    hd->SetBinContent(1,  scale_b*tot_B );
    hd->SetBinError(1,  sqrt(scale_b*tot_B) );
  }
   
  //data
  sig_vals.clear();
  bkg_vals.clear();


  for (int i=1; i<hb->GetSize()-2; i++){
    double N_h0 = hb->GetBinContent(i);
    double N_h1 = hb->GetBinContent(i) + hs->GetBinContent(i);
    double x = hb->GetBinContent(i);
    double dx = hb->GetBinError(i);
    double chisq_h1 = pow(N_h1-N_h0,2)/(dx*dx);
    cout << "Bin " << i << " chi-squared " << chisq_h1 << "  " << hb->GetBinContent(i) << endl;
  }
  
  cout << "Before scaling: " << endl;
  cout << hs->GetBinContent(5) << "\t" << hs->GetBinError(5) << endl;
  cout << hb->GetBinContent(5) << "\t" << hb->GetBinError(5) << endl;
  hs->Scale(lumi);
  hb->Scale(lumi);
  hd->Scale(lumi);
  cout << "After scaling: " << endl;
  cout << hs->GetBinContent(5) << "\t" << hs->GetBinError(5) << endl;
  cout << hb->GetBinContent(5) << "\t" << hb->GetBinError(5) << endl;

  hs->Sumw2();
  hb->Sumw2();
  hd->Sumw2();
  
  

  
  //declare nusiance parameters and shape function
  char *ename[15];
  double nps_low[15];
  double nps_high[15];
  double lowsigma[15];
  double highsigma[15];
  TH1 *lowshape[15];
  TH1 *highshape[15];

  for(int i=0;i<15;i++)
    {
      nps_low[i]  = 0;
      nps_high[i] = 0;
      lowsigma[i] = 0;
      highsigma[i]= 0;
      lowshape[i] = 0;
      highshape[i]= 0;
    }

  //========================================
  //Construct test/null hypothesis for fitting
  csm_model* nullhyp = new csm_model();
  csm_model* testhyp = new csm_model();
  csm_model* nullhyp_pe = new csm_model();
  csm_model* testhyp_pe = new csm_model();

  //add background templates
  int nps_count=0;
  ename   [nps_count] = (char*)"BkgSys";
  nps_low [nps_count] = -bkgSys;
  nps_high[nps_count] =  bkgSys;
  nps_count++;

  //null hypothesis = B
  nullhyp->add_template(hb,scale_b,nps_count,ename,nps_low,nps_high,lowshape,lowsigma,highshape,highsigma,1,0,(char*)"Xtt");
  testhyp->add_template(hb,scale_b,nps_count,ename,nps_low,nps_high,lowshape,lowsigma,highshape,highsigma,1,0,(char*)"Xtt");

  //add signal templates
  nps_count=0;
  ename   [nps_count] = (char*)"SigSys";
  nps_low [nps_count] = -sigSys;
  nps_high[nps_count] =  sigSys;
  nps_count++;

  //test hypothesis = S+B
  testhyp->add_template(hs,scale_s,nps_count,ename,nps_low,nps_high, lowshape,lowsigma,highshape,highsigma,1,0,(char*)"Xtt");

  //========================================
  //Construct test/null hypothesis for pseudo-experiments.
  //They can be different from test/null hypothesis for fitting in order to evaluate the wrong fitting functions
  //For now, set fitting and pseudo-experiments templates the same
  nullhyp_pe = nullhyp->Clone();
  testhyp_pe = testhyp->Clone();

  //Debug purpose
  testhyp_pe->print();


  //========================================
  // Have a visualization of how fitting looks like
  // It has nothing to do with limit calculation
  TCanvas * mycanvas = (TCanvas *) new TCanvas("Canvas1","Canvas1",0,0,600,800);
  mycanvas->Divide(1,2);

  mycanvas->cd(1);
  csm* mycsm = new csm();
  mycsm->set_htofit(hd,(char*)"Xtt");
  mycsm->set_modeltofit(testhyp);
  double chisq = mycsm->chisquared();

  csm_model* bestnullfit = mycsm->getbestmodel();
  bestnullfit->plotwithdata((char*)"Xtt",hd);
  cout << "chisq from fitter " << chisq << endl;
  delete mycsm;


  //======================================
  //Construct limit calculator - mclimit_csm
  //  Sensitivity, Significance calculation

  mclimit_csm* mymclimit = (mclimit_csm*) new mclimit_csm();

  //print out pseudo-experiments details
  //mymclimit->setpxprintflag(1);

  cout << "setting hypotheses" << endl;
  mymclimit->set_null_hypothesis(nullhyp);
  mymclimit->set_test_hypothesis(testhyp);
  mymclimit->set_null_hypothesis_pe(nullhyp_pe);
  mymclimit->set_test_hypothesis_pe(testhyp_pe);
  mymclimit->set_datahist(hd,(char*)"Xtt");

  mymclimit->set_npe(10000);
  mymclimit->run_pseudoexperiments();

  //Model
  cout <<" Bkg "<< hb->Integral() <<" Sig "<< hs->Integral() <<" Data "<< hd->Integral()<<endl;

  // Sensitivity of test hypothesis
  cout << "NULL HYPOTHESIS PROBABILITIES" << endl;
  cout << ">2 sigma probability: " << mymclimit->p2sigman() << endl; // Probability of a 3-sigma evidence assuming test hyp. is true
  cout << ">3 sigma probability: " << mymclimit->p3sigman() << endl; // Probability of a 3-sigma evidence assuming test hyp. is true
  cout << ">5 sigma probability: " << mymclimit->p5sigman() << endl; // Probability of a 5-sigma evidence assuming test hyp. is true
  cout << "----------------------" << endl;
  cout << "TEST HYPOTHESIS PROBABILITIES" << endl;
  cout << ">2 sigma probability: " << mymclimit->p2sigmat() << endl; // Probability of a 3-sigma evidence assuming test hyp. is true
  cout << ">3 sigma probability: " << mymclimit->p3sigmat() << endl; // Probability of a 3-sigma evidence assuming test hyp. is true
  cout << ">5 sigma probability: " << mymclimit->p5sigmat() << endl; // Probability of a 5-sigma evidence assuming test hyp. is true

  cout << "----------------------" << endl;
  cout << "CLS values : " << endl;
  cout << "NULL HYPOTHESIS:" << endl;
  cout << "clb_exp_bmed NULL: " << mymclimit->clbexpbmed() << endl;
  cout << "clsb_exp_bmed NULL: " << mymclimit->clsbexpbmed() << endl;
  cout << "cls_exp_bmed NULL: " << mymclimit->clsexpbmed() << endl;
  cout << "----------------------" << endl;
  cout << "TEST HYPOTHESIS: " << endl;
  cout << " clb_exp_bmed TEST   : " << mymclimit->clbexpsmed() << endl;
  //  cout << " clsb_exp_smed TEST: " << mymclimit->clsbexpsmed() << endl;
  cout << " cls_exp_smed TEST: " << mymclimit-> clsexpsmed() << endl;
  cout << "----------------------" << endl;
  cout << "Bottom line " << endl;
  cout << "CLs: " << mymclimit->cls() << endl;
  cout << "CLs+b: " << mymclimit->clsb() << endl;
  cout << "CLb: " << mymclimit->clb() << endl;
  cout << "----------------------" << endl;
  
  //Null hypothesis
  cout << "Null Hypothesis:"
       << "Median of test statistics: " << mymclimit->tsbmed()
       << " p-value=" <<mymclimit->omclbexpbmed()
       << " significance="<< TMath::ErfcInverse(mymclimit->omclbexpbmed()*2)*sqrt(2)
       <<endl;

  //Expected significance
  cout << "Test Hypothesis:"
       << "Median of test statistics: " << mymclimit->tssmed()
       << " p-value=" <<mymclimit->omclbexpsmed()
       << " significance=" << TMath::ErfcInverse(mymclimit->omclbexpsmed()*2)*sqrt(2)
       <<endl;

  // Significance of input data
  double ts_data = mymclimit->ts();
  double pval=mymclimit->omclb() ;
  double significance=TMath::ErfcInverse(pval*2)*sqrt(2); //convert pval to one-side gaussian significance

  cout << " Input Data: test statistics = " << ts_data
       << " p-val = " << pval
       << " significance= " << significance << endl;

  delete mymclimit;

  return 0;
}
