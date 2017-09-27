//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 17 14:39:53 2014 by ROOT version 5.34/14
// from TTree starlightTree/starlightTree
// found on file: starlight.root
//////////////////////////////////////////////////////////

#ifndef e_AnalyzeTree_h
#define e_AnalyzeTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class e_AnalyzeTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   float           q2;
   float           Egamma;
   TLorentzVector  *source;
   TLorentzVector  *parent;
   TClonesArray    *daughters;

   // List of branches
   //TBranch        *b_q2;
   //TBranch        *b_Egamma;
   TBranch        *b_source;
   TBranch        *b_parent;   //!
   TBranch        *b_daughters;   //!

   e_AnalyzeTree(string filename = "starlight.root", TTree *tree=0);
   virtual ~e_AnalyzeTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef e_AnalyzeTree_cxx
e_AnalyzeTree::e_AnalyzeTree(string filename, TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("starlight.root");
      if (!f || !f->IsOpen()) {
	std::cout<<"Lomnitz"<<filename<<std::endl;
	f = new TFile(filename.c_str());
      }
      f->GetObject("starlightTree",tree);

   }
   Init(tree);
}

e_AnalyzeTree::~e_AnalyzeTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t e_AnalyzeTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t e_AnalyzeTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void e_AnalyzeTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   source = 0 ;
   parent = 0;
   daughters = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Egamma", &Egamma);
   fChain->SetBranchAddress("Q2", &q2);
   fChain->SetBranchAddress("source", &source, &b_source);
   fChain->SetBranchAddress("parent", &parent, &b_parent);
   fChain->SetBranchAddress("daughters", &daughters, &b_daughters);
   Notify();
}

Bool_t e_AnalyzeTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void e_AnalyzeTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t e_AnalyzeTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef e_AnalyzeTree_cxx
