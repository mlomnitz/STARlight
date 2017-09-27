#include <iostream>
#include "TF1.h"
#include "myHists.h"
#include "StyleUtilities.h"
using namespace std;
using namespace ana_const;
ana_hists::ana_hists(std::string const accel_name, std::string const particle_name, double const xsec)
{
  accel = accel_name;
  particle = particle_name;
  sample_xsec = xsec;
  init_histos();
} 
ana_hists::~ana_hists()
{
  cout<<"Closing out"<<endl;
  //delete_histos();
}
void ana_hists::init_histos()
{
  setGraphicsStyle();
  hParentW_PhotonEQ2 = new TH3D("ParentW_vs_PhotonEQ","Hadronic invariant mass vs photon E,q",200,0,200,
				750,0,1.5,600,0,60);
  for( int ii = 0 ; ii<n_detectors+1; ++ii)
    hParentYPt[ii] = new TH2D(Form("detector_%i",ii),Form("detector_%i",ii),
			      plot_bins, plot_limits[0], plot_limits[1], 500,0,5);
  load_data();
  
}
void ana_hists::delete_histos()
{
  if(hParentW_PhotonEQ2) delete hParentW_PhotonEQ2;
  for( int ii = 0 ; ii < n_detectors+1; ++ii)
    if(hParentYPt[ii]) delete hParentYPt[ii];
  //
  if(gDataXsec_vs_Q2) delete gDataXsec_vs_Q2;
  if(gDataXsec_vs_W) delete gDataXsec_vs_W;
}
void ana_hists::load_data()
{
  gDataXsec_vs_Q2 = new TGraph("HERA_comparison/H1_phi.csv","%lg , %lg");
  int const n_points = 5;
  double h1_x[n_points] = {41, 54, 67, 80, 93};
  double h1_y[n_points] = {41.2, 55.1, 49.2, 57.5, 69.6};
  gDataXsec_vs_W = new TGraph(n_points,h1_x,h1_y);
  setStyle(gDataXsec_vs_W,21, kBlack);
}
void ana_hists::fill_detector_hists(double d1_y, double d2_y, double p_y, double p_pt)
{
  for( int ii = 0 ; ii < n_detectors; ++ii){
    if( std::fabs(d1_y) < detector_rap[ii] && std::fabs(d2_y) < detector_rap[ii] )
      hParentYPt[ii] -> Fill( p_y, p_pt);
  }
  hParentYPt[n_detectors] -> Fill( p_y, p_pt ); //all parents
}
void ana_hists::make_detector_plot()
{
  float const all_parents = hParentYPt[n_detectors] ->GetEntries();
  float const rapididty_bin_width = hParentYPt[n_detectors] ->GetXaxis()->GetBinWidth(1);
  float const ratio = sample_xsec/all_parents/rapididty_bin_width;
  TH1* parent_y[n_detectors+1];
  parent_y[n_detectors] = hParentYPt[n_detectors]->ProjectionX( Form("this_hist%i",n_detectors) );
  parent_y[n_detectors]->Scale(ratio);
  // -- 
  TCanvas *cv = canvas(false);
  cv->cd();
  TLegend *leg = legend(Form("eSTARlight %s @ %s ",particle.c_str(),accel.c_str()),
			0.75,0.93,0.15,0.3);
  //
  setStyle(parent_y[n_detectors], 29, kBlack);
  parent_y[n_detectors]->Draw("l");
  leg->AddEntry(parent_y[n_detectors], "Generated parents","l");
  // -- 
  for( int ii = 0 ; ii < n_detectors ; ++ii){
    parent_y[ii] = hParentYPt[ii]->ProjectionX( Form("this_hist%i",ii) );
    float eff = parent_y[ii]->GetEntries()/all_parents;      
    parent_y[ii]->Scale(ratio);
    setStyle(parent_y[ii], marker[ii],plot_color[ii]);
    parent_y[ii]->Draw("same:l");
    leg->AddEntry(parent_y[ii],Form("|#eta_{i}|<%1.0lf, eff.=%.2lf ",detector_rap[ii],eff),"l");
  }
  leg->Draw();
  cout<<"Lomitz here"<<endl;
  cv->SaveAs( Form("scaled_%s_y.eps",particle.c_str()) );
}
void ana_hists::fill_W_hist(double W, double Egamma, double Q2)
{
  hParentW_PhotonEQ2->Fill(W,Egamma, Q2);
}
void ana_hists::make_W_plot()
{
  double mv = particle_mass[0];
  TGraph* w_plots[qbins];
  TCanvas* w_cv = canvas(false);
  w_cv->cd();
  w_cv->SetLogy();
  gDataXsec_vs_W->Draw("a:p");
  gDataXsec_vs_W->GetYaxis()->SetRangeUser(1e-3,1e4);
  gDataXsec_vs_W->GetXaxis()->SetLimits(35,250);
  gDataXsec_vs_W->GetYaxis()->SetTitle(" #sigma_{#gamma p #rightarrow #phi p} [nb]");
  for( int i = 0 ; i< qbins; ++i){
    int minQ2 = hParentW_PhotonEQ2->GetZaxis()->FindBin(Q2Bins[i][0]+1e-6);
    int maxQ2 = hParentW_PhotonEQ2->GetZaxis()->FindBin(Q2Bins[i][1]-1e-6);
    hParentW_PhotonEQ2->GetZaxis()->SetRange(minQ2,maxQ2);
    double this_plot_x[wbins] = {0};
    double this_plot_y[wbins] = {0};
    double q2 = hParentW_PhotonEQ2->GetMean(3);
    double Egamma = hParentW_PhotonEQ2->GetMean(2);
    double scale = gDataXsec_vs_Q2->Eval(q2+1.019*1.019)/(Q2Bins[i][1]-Q2Bins[i][0]); // no weighting
    cout<<Q2Bins[i][0]<< " < Q2 < "<<Q2Bins[i][1]<<" : <Q2> = "<<hParentW_PhotonEQ2->GetMean(3)
	<<" <Egamma> = "<<Egamma<<endl;
    double bin_entries = 0;
    for( int j = 0 ; j < wbins; ++j){
      int min = hParentW_PhotonEQ2->GetXaxis()->FindBin(WBins[j]+1e-6);
      int max  = hParentW_PhotonEQ2->GetXaxis()->FindBin(WBins[j+1]-1e-6);
      //
      hParentW_PhotonEQ2->GetZaxis()->SetRange(minQ2,maxQ2); 
      hParentW_PhotonEQ2->GetXaxis()->SetRange(min,max); 
      double Egamma = hParentW_PhotonEQ2->GetMean(2);
      TH1D* this_bin  = (TH1D*)(hParentW_PhotonEQ2->Project3D("x")->Clone(Form("this_bin_%i",j)));
      //
      this_plot_x[j] = hParentW_PhotonEQ2->GetMean(1);
      double integral =  0 ;
      for( int a = 0 ; a<(max-min); ++a){
	integral+=this_bin->GetBinContent(a+1);
      }
      double d_int = integral/(photonFlux(Egamma,q2)*Egamma)/(WBins[j+1]-WBins[j]); //kinda works      
      //double d_int = integral/(photonFlux(Egamma,q2))/(WBins[j+1]-WBins[j]); //kinda works      
      this_plot_y[j] = d_int;
      bin_entries+= d_int;
      //      cout<<" W bin" <<j <<" : "<<bin_entries<<endl;
    }
    scale = scale/bin_entries;
    w_plots[i] = new TGraph(wbins,this_plot_x,this_plot_y);
    setStyle(w_plots[i],w_marker[i],w_color[i]);    
    double temp_scale = get_scale(gDataXsec_vs_W,w_plots[i]);
    cout<<" TEsting "<<scale<<" vs. "<<temp_scale<<endl;
    for( int jj = 0 ; jj<w_plots[i]->GetN(); ++jj) {
      //w_plots[i]->GetY()[jj] *= scale;
      w_plots[i]->GetY()[jj] *= temp_scale;
    }
    if( i == 0 )
      w_plots[i]->Draw("same:p");
      //w_plots[i]->Draw("same:p");
    //if( i == 0 ){
    //  w_plots[i]->Draw("a:p");
    //  w_plots[i]->GetYaxis()->SetRangeUser(1,1e3);
    //}
    //else
    //      w_plots[i]->Draw("same:p");
  }
  w_cv->SaveAs("w_plot.eps");
}
double ana_hists::photonFlux(double Egamma, double Q2)
{
  //Some constants
  double mel = 0.000510998928; //electron mass
  double rap1 = acosh(53921); //27.5 GeV electrons HERA
  double rap2 = -acosh(981);  //920 GeV protons HERA
  double alpha    = 1/137.035999074;
  double pi       = 3.141592654;
  double targetLorentzGamma = cosh(rap1-rap2);
  double electronEnergy         = targetLorentzGamma * mel;
  double const ratio = Egamma/electronEnergy;
  double const minQ2 = std::pow( mel*Egamma,2.0) / (electronEnergy*(electronEnergy - Egamma));
  double to_ret = alpha/(pi) *( 1- ratio + ratio*ratio/2. - (1-ratio)*( fabs(minQ2/Q2)) );
  return to_ret/( Egamma*fabs(Q2) );
}
double ana_hists::get_scale(TGraph* data, TGraph* simu)
{
  TF1* data_f = new TF1("data_f","[0]*x**[1]",40,100);
  data->Fit(data_f,"0");
  //
  TF1* simu_f = new TF1("simu_f","[0]*x**[1]",40,100);
  simu_f->FixParameter(1,data_f->GetParameter(1));
  simu->Fit(simu_f,"0");
  return data_f->GetParameter(0)/simu_f->GetParameter(0);
  
}
