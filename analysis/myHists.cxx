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
  // - All parents Y
  hParentY_all = new TH1D("hParentY_all", "hParentY_all", nYBins,yPlotLim[0],yPlotLim[1]);
  setStyle(hParentY_all, 22, kBlack);
  hParentY_all->SetLineColor(1);
  hParentY_all->GetXaxis()->SetTitle("VM rapidity");
  hParentY_all->GetYaxis()->SetTitle("# of events [arb. units]");
  // - Detector acceptance plot
  for( int iDet = 0; iDet < ana_const::n_detectors; ++iDet){
    hParentY_acceptance[iDet] = new TH1D(Form("acceptance_plot_%i",iDet), Form("acceptance_plot_%i",iDet),
					 nYBins,yPlotLim[0],yPlotLim[1]);
    setStyle(hParentY_acceptance[iDet], marker_type[iDet], marker_color[iDet]);
    hParentY_acceptance[iDet]->SetLineColor(marker_color[iDet]);
  }
  // - Photon energy band plot
  for( int iGamma = 0 ; iGamma < ana_const::n_EgammaBins; ++iGamma){
    hParentY_Egamma[iGamma] = new TH1D(Form("gammaBand_%i",iGamma), Form("gammaBand_%i",iGamma), 
				       nYBins,yPlotLim[0],yPlotLim[1]);
    setStyle(hParentY_Egamma[iGamma], marker_type[iGamma],marker_color[iGamma]);
    hParentY_Egamma[iGamma]->SetLineColor(marker_color[iGamma]);
    hParentY_Egamma[iGamma]->SetFillColor(marker_color[iGamma]);
    hParentY_Egamma[iGamma]->SetFillStyle(3001);
  }
}
void ana_hists::delete_histos()
{
  //
}
void ana_hists::fill_detector_hists(double parentY, double d1_y, double d2_y, double eGamma)
{
  //All parents
  hParentY_all->Fill(parentY);
  // - Detector plot
  for( int ii = 0 ; ii < ana_const::n_detectors; ++ii){
    if( std::fabs(d1_y) < detector_rap[ii] && std::fabs(d2_y) < detector_rap[ii] )
      hParentY_acceptance[ii]->Fill(parentY);
  }
  // - Egamma bands
  for( int ii = 0 ; ii<ana_const::n_EgammaBins; ++ii){
    if( eGamma >= egammaBins[ii] && eGamma < egammaBins[ii+1] )
      hParentY_Egamma[ii]->Fill(parentY);
  }
}
void ana_hists::make_plots()
{
  make_detector_plot();
  make_gammaE_plot();
}
void ana_hists::make_detector_plot()
{
  TLine *zero = new TLine(-6,20,10,20);
  zero->SetLineWidth(3);
  TCanvas* cv_1 = canvas(false);
  cv_1->cd();
  cv_1->SetLogy();
  TLegend* leg1 = legend("Detector acceptance",0.6,0.9,0.65,0.9);
  hParentY_all->Draw("l");
  hParentY_all->GetYaxis()->SetRangeUser(20,8E3);
  hParentY_all->GetXaxis()->SetRangeUser(-6,10);
  leg1->AddEntry(hParentY_all,"All VM","l");
  for( int ii = 0 ; ii < ana_const::n_detectors; ++ii){
    hParentY_acceptance[ii]->Draw("same:l");
    leg1 ->AddEntry(hParentY_acceptance[ii], Form("|#eta| < %i = %2.2lf", ii+1,
					      hParentY_acceptance[ii]->GetEntries()/hParentY_all->GetEntries()), "l");
  }
  zero->Draw("same");
  leg1->Draw();
  cv_1->SaveAs("detector_acceptance.eps");  
}
void ana_hists::make_gammaE_plot()
{
  TLine *zero = new TLine(-6,20,10,20);
  zero->SetLineWidth(3);
  TCanvas* cv_1 = canvas(false);
  cv_1->cd();
  cv_1->SetLogy();
  TLegend* leg1 = legend("Photon energy [GeV]",0.6,0.9,0.65,0.9);
  hParentY_all->Draw("l");
  for( int ii = 0 ; ii < ana_const::n_EgammaBins; ++ii){
    hParentY_Egamma[ii]->Draw("same:c:CF");
    if( ii == n_EgammaBins-1 )
      leg1 ->AddEntry(hParentY_Egamma[ii],Form("%0.0lf < k ",egammaBins[ii]) , "l");
    else
      leg1 ->AddEntry(hParentY_Egamma[ii],Form("%0.0lf < k < %0.0lf",egammaBins[ii],egammaBins[ii+1]) , "l");
  }
  zero->Draw("same");
  leg1->Draw();
  cv_1->SaveAs("photon_energy_bands.eps");  
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
