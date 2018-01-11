// This macro reads a starlight output file (default name slight.out) and creates a root file 
// with  TLorentzVectors for the parent and a TClonesArray of TLorentzVectors for the daughter 
// particles.  The output is stored in a root file (default name starlight.root) with one branch 
// labeled "parent" and the other labeled "daughters". Any number of daughter tracks can be 
// accomodated.  Daughter species currently accomodated are:  electrons, muons, charged or neutral 
// pions, charged or neutral kaons, and protons.  
//
// To use this macro, open root and then 
// type .x convertStarlightAsciiToTree.C("inputfilename", "outputfilename")
// 
// The root file produced can be examined in a root TBrowser.
//
// A macro to read this root file and make some standard plots is also provided.  This macro is 
// called AnalyzeTree.cxx; it can be compiled and run with the anaTree.C macro by opening root 
// and typing .x anaTree.C()

#include <iostream>
#include <fstream>
#include <sstream>

#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"


using namespace std;
double IDtoMass(int particleCode);

struct source_electron
{
  TLorentzVector* p;
  float gamma_mass;
  void set(float px, float py, float pz, float E, float Q2){
    p = new TLorentzVector(px,py,pz,E);
    gamma_mass = Q2;
  };
};

void ConvertStarlightAsciiToTree(const char* inFileName  = "slight.out",
                        const char* outFileName = "starlight.root")
{

	// create output tree
	TFile* outFile = new TFile(outFileName, "RECREATE");
	if (!outFile) {
		cerr << "    error: could not create output file '" << outFileName << "'" << endl;
		return;
	}
	TTree*          outTree           = new TTree("starlightTree", "starlightTree");
	TLorentzVector* parentParticle    = new TLorentzVector();
  	TClonesArray*   daughterParticles = new TClonesArray("TLorentzVector");
	TLorentzVector* source            = new TLorentzVector();
	TLorentzVector* target            = new TLorentzVector();
	Float_t        q2, Egamma;

	outTree->Branch("Egamma",         &Egamma, "Egamma/F");
	outTree->Branch("Q2",         &q2, "q2/F");
	outTree->Branch("Target","TLorentzVector", &target,            32000, -1);
	outTree->Branch("source",    "TLorentzVector", &source,            32000, -1);
	outTree->Branch("parent",    "TLorentzVector", &parentParticle,    32000, -1);
	outTree->Branch("daughters", "TClonesArray",   &daughterParticles, 32000, -1);

	ifstream inFile;
	inFile.open(inFileName);
	unsigned int countLines = 0;
	int tot_events=0;
	while (inFile.good()) {
		string       line;
		stringstream lineStream;
		
		// read EVENT
		string label;
		int    eventNmb, nmbTracks;
		if (!getline(inFile, line))
			break;
		++countLines;
		++tot_events;
		lineStream.str(line);
		assert(lineStream >> label >> eventNmb >> nmbTracks);
		if (!(label == "EVENT:"))
			continue;
		
		lineStream.clear();

		// read vertex
		if (!getline(inFile, line))
			break;
		++countLines;
		lineStream.str(line);
		assert(lineStream >> label);
		assert(label == "VERTEX:");
		
		*parentParticle = TLorentzVector(0, 0, 0, 0);

		lineStream.clear();		
		
		//read gamma
		if( !getline(inFile,line))
		  break;
		lineStream.str(line);
		//cout<<line.c_str()<<endl;
		assert(lineStream >> label >> Egamma >> q2 );
		assert(label == "GAMMA:");
		lineStream.clear();
 
		// read target
		if(!getline(inFile, line))
		  break;
		++countLines;
		lineStream.str(line);
		//cout<<line.c_str()<<endl;
		float tpx=0., tpy=0., tpz=0., tE=0.;
		assert(lineStream >> label >> tpx >> tpy >> tpz >> tE) ;
		assert(label == "TARGET:");
		*target = TLorentzVector(tpx, tpy, tpz, tE);

		lineStream.clear();
		// read source
		if(!getline(inFile, line))
		  break;
		++countLines;
		lineStream.str(line);
		//cout<<line.c_str()<<endl;
		float px=0., py=0., pz=0., E=0.;
		assert(lineStream >> label >> px >> py >> pz >> E) ;
		assert(label == "SOURCE:");
		*source = TLorentzVector(px, py, pz, E);

		lineStream.clear();
		//
		for (int i = 0; i < nmbTracks; ++i) {
			// read tracks
			int    particleCode;
			double momentum[3];
			if (!getline(inFile, line))
			  break;
			++countLines;
			lineStream.str(line);
			assert(lineStream >> label >> particleCode >> momentum[0] >> momentum[1] >> momentum[2]);
			assert(label == "TRACK:");
			Double_t daughterMass = IDtoMass(particleCode);
			if (daughterMass < 0) {break;}
			const double E = sqrt(  momentum[0] * momentum[0] + momentum[1] * momentum[1]
			                      + momentum[2] * momentum[2] + daughterMass * daughterMass);
			new ( (*daughterParticles)[i] ) TLorentzVector(momentum[0], momentum[1], momentum[2], E);
			*parentParticle += *(static_cast<TLorentzVector*>(daughterParticles->At(i)));

			lineStream.clear();
		}
		daughterParticles->Compress();
		outTree->Fill();
	}
	outTree->Write();
	if (outFile) {
		outFile->Close();
		delete outFile;
	}
	cout<<"Processed "<<tot_events<<" events"<<endl;
}

	double IDtoMass(int particleCode){
		double mass;
		if (particleCode == 2 || particleCode==3) {mass = 0.00051099907;} // electron
		else if (particleCode == 5 || particleCode==6) {mass = 0.105658389;} // muon
		else if (particleCode == 8 || particleCode==9)  {mass = 0.13956995;} // charged pion
		else if (particleCode == 7) {mass = 0.1345766;} // neutral pion
		else if (particleCode == 11|| particleCode==12) {mass = 0.493677;} // charged kaon
		else if (particleCode == 10 || particleCode == 16)  {mass = 0.497614;} // neutral kaon
		else if (particleCode == 14)	{mass = 0.93827231;} // proton
		else {
			cout << "unknown daughter particle (ID = " << particleCode << "), please modify code to accomodate" << endl;
			mass = -1.0;
//			exit(0); 
		     } 

		return mass;
	}
