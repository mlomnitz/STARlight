// this macro compiles and runs AnalyzeTree.cxx, which takes as input the 
// starlight.root file produced by convertStarlightAsciiToTree.cxx
// output histograms are stored in starlight_histos.root 
//
#include <TSystem.h>
//#include "e_AnalyzeTree.h"

//class e_AnalyzeTree;

void e_AnaTree(string filename = "starlight.root"){
  //gROOT->SetMacroPath("/Users/michaellomnitz/Documents/");
  //gROOT->ProcessLine(".L ./e_AnalyzeTree.cxx++");
  //  gSystem->Load("e_AnalyzeTree");
  std::cout<<"Lomnitz::: "<<filename<<std::endl;
  e_AnalyzeTree* l = new e_AnalyzeTree(filename);
  l->Loop();
}
