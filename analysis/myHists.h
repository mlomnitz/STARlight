#ifndef __myHists__
#define __myHists__
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TGraph.h"
//#include <string>
//#include <map>
//#include <iostream>
//
namespace ana_const{
  // -- - - -- 
  int const n_detectors = 3;
  double const detector_rap[n_detectors] = {1., 2., 3.};
  int const plot_color[n_detectors] = { 600, 418, 633};
  int const marker[n_detectors] = {20, 22, 33};
  int const plot_limits[2] = { -1, 5};
  int const plot_bins = 300;
  //enum particle{ phi, rho, j_psi, upsilon};
  string const particle_label[4] = { "#phi", "#rho", "J/#psi", "#Upsilon" };
  double const particle_mass[4] = { 1.02, 0.77, 3.097, 9.46};
  // -- - - -- 
  int const qbins  = 4;
  int w_color[qbins] = { 600, 418, 633, 618};
  int w_marker[qbins] = {20, 22, 33, 24};
  double Q2Bins[qbins][2] = { {3.2,3.4}, {3.,4.},
			      {4.,5.}, {5.,10.} };
  int const wbins = 6;
  double WBins[wbins +1] = {40,45,50,55,60,70,80};
  //
  int const ntBins = 7;
  double tBins[ntBins+1] = {0,0.1, 0.2, 0.3, 0.4, 0.6,
			    0.8,1.};
  //
};
//
class ana_hists
{
 public:
  ana_hists(std::string const accel_name, std::string const particle_name, double const xsec);
  ~ana_hists();
  void fill_detector_hists(double d1_y, double d2_y, double p_y, double p_pt);
  void make_detector_plot();
  void fill_W_hist(double W, double Egamma, double Q2, double t);
  void make_W_plot();
  void make_t_plot();
 private:
  std::string accel;
  std::string particle;
  double sample_xsec;
  // 
  THn* hWEgammaQ2t;
  TH2* hParentYPt[ana_const::n_detectors+1];
  TGraph* gDataXsec_vs_Q2;
  TGraph* gDataXsec_vs_W;
  //
  void init_histos();
  void delete_histos();
  void load_data();
  double photonFlux(double Egamma, double Q2);
  double photonFlux_HERA(double Egamma, double Q2);
  double get_scale(TGraph* data, TGraph* simu);
};

#endif //__myHists__
