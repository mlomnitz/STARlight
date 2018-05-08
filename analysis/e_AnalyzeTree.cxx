// This macro reads the starlight.root file produced by convertStarlightAsciiToTree.C, 
// which contains TLorentzVectors for the parents and a TClonesArray of TLorentzVectors 
// for the daughters. 
//
// It creates histograms of the p_T and rapidity of the daughters, as well as the p_T, 
// rapidity and mass of the parent.  While the parents may have been created as the 
// vector sum of any number of daughter particles, this macro currently produces 
// histograms for only the first two daughter particles.  The daughter histograms are
// called D1Pt, D2Pt, D1Rapidity, and D1Rapidity.  Parent histograms are 
// named ParentPt, ParentRapidity, and ParentMass.  The histograms are stored in 
// starlight_histos.root.  
//
// To use this macro, you must first run convertStarlightAsciiToTree.C to produce the 
// starlight.root file.  If needed, modify the file AnalyzeTree.h to call your input file
// (as downloaded, it calls starlight.root).  Then open root and type .x anaTree.C .
#ifndef e_AnalyzeTree_cxx
#define e_AnalyzeTree_cxx
#include "e_AnalyzeTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include "myHists.h"
//
ana_hists my_hist("HERA", "phi",17.);
//
void e_AnalyzeTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MyClass.C
//      Root > MyClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   double max = 0 ;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if( jentry%(nentries/10) == 0 ) 
       cout<<"Working on "<<jentry<<" / "<<nentries<<endl;
     daughters->Clear();
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     TLorentzVector *D1 = (TLorentzVector*)daughters->At(0);
     TLorentzVector *D2 = (TLorentzVector*)daughters->At(1);
     // if desired, acceptance or analysis cuts can be applied here, before filling histograms
     // -- Fill detector acceptance plots     
     my_hist.fill_detector_hists( parent->Rapidity(), D1->PseudoRapidity(), D2->PseudoRapidity(), Egamma);
     
   }// jentry
   my_hist.make_plots();
}
#endif
