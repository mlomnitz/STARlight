#ifndef __myHists__
#define __myHists__
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include <string>
//#include <map>
//#include <iostream>
//
namespace ana_const{
  // -- - - --
  int const max_n_plots = 5;
  int const marker_type[max_n_plots] = {20, 22, 33, 24, 21};
  int const marker_color[max_n_plots] = {600, 418, 633, 618, 797}; //kBlue, kGreen + 2, kRed+1, kMagenta+2, kOrange -3
  // -- - - --
  int const nYBins = 1000;
  int const yPlotLim[2] = {-10,10};
  // -- - - -- 
  int const n_detectors = 3;
  double const detector_rap[n_detectors] = {1., 2., 3.};
  // -- - - --
  int const n_EgammaBins = 5;
  double const egammaBins[n_EgammaBins +1 ] = {0,20,60,300,1800,100000};
};
//
class ana_hists
{
 public:
  ana_hists(std::string const accel_name, std::string const particle_name, double const xsec);
  ~ana_hists();
  void fill_detector_hists(double parentY, double d1_y, double d2_y, double eGamma);
  void make_plots();
 private:
  std::string accel;
  std::string particle;
  double sample_xsec;
  //
  TH1* hParentY_all;
  TH1* hParentY_acceptance[ana_const::n_detectors];
  TH1* hParentY_Egamma[ana_const::n_EgammaBins];
  TH2* hEtheta_vs_gammaQ2;
  //
  void init_histos();
  void delete_histos();
  void make_detector_plot();
  void make_gammaE_plot();
  // - - - Other utilities
  double get_scale(TGraph* data, TGraph* simu);
};

#endif //__myHists__
