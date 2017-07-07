///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 264                         $: revision of last commit
// $Author:: jseger                   $: author of last commit
// $Date:: 2016-06-06 21:05:12 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////
#define _makeGammaPQ2_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "e_narrowResonanceCrossSection.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
e_narrowResonanceCrossSection::e_narrowResonanceCrossSection(const inputParameters& inputParametersInstance, const beamBeamSystem&  bbsystem)
	:photonNucleusCrossSection(inputParametersInstance, bbsystem)
{
	_Ep         = inputParametersInstance.protonEnergy();	
	//this is in target frame
	_electronEnergy = inputParametersInstance.electronEnergy();
	//_target_beamLorentz = inputParametersInstance.beam2LorentzGamma();
	_target_beamLorentz = inputParametersInstance.beamLorentzGamma();
	_boost = std::acosh(inputParametersInstance.beam1LorentzGamma())
	  -std::acosh(inputParametersInstance.beam2LorentzGamma());
	_boost = _boost/2;
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();
	// Now saving the photon energy limits
	_cmsMaxPhotonEnergy=inputParametersInstance.cmsMaxPhotonEnergy();
	_cmsMinPhotonEnergy=inputParametersInstance.cmsMinPhotonEnergy();
	//
	_VMnumEgamma = inputParametersInstance.nmbEnergyBins();
}


//______________________________________________________________________________
e_narrowResonanceCrossSection::~e_narrowResonanceCrossSection()
{ }


//______________________________________________________________________________
void
e_narrowResonanceCrossSection::crossSectionCalculation(const double)  // _bwnormsave (unused)
{
        // This subroutine calculates the vector meson cross section assuming
        // a narrow resonance.  For reference, see STAR Note 386.
  
        double W, dEgamma, minEgamma;
	double ega[3] = {0};
	double int_r,dR;
	double int_r2, dR2;
	int    iEgamma, nEgamma,beam;
	
	//Integration is done with exponential steps, in target frame
	//nEgamma = _VMnumEgamma;
	nEgamma = 1000;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/nEgamma;
	minEgamma = std::log(_targetMinPhotonEnergy);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	//
        printf(" gamma+nucleon threshold (CMS): %e GeV \n", _cmsMinPhotonEnergy);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();

	if( A_2 == 0 && A_1 >= 1 ){
	  // pA, first beam is the nucleus and is in this case the target  
	  beam = 1; 
	} else if( A_1 ==0 && A_2 >= 1){       
	  // photon energy in Target frame
	  beam = 2; 
	} else {
	  // photon energy in Target frame
	  beam = 2; 
	}
  
 	int_r=0.;
	int_r2= 0;
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	// These might be useful in a future iteration
	//double target_cm = -acosh(_beamLorentzGamma);
	//double target_cm = acosh(_target_beamLorentz);
	//double coshy = cosh(target_cm);
	//double sinhy = sinh(target_cm);
	int nQ2 = 1000;
        printf(" gamma+nucleon threshold (Target): %e GeV \n", _targetMinPhotonEnergy);
	for(iEgamma = 0 ; iEgamma < nEgamma; ++iEgamma){    // Integral over photon energy
	  // Target frame energies
	  ega[0] = exp(minEgamma + iEgamma*dEgamma );
	  ega[1] = exp(minEgamma + (iEgamma+1)*dEgamma );
	  ega[2] = 0.5*(ega[0]+ega[1]);

	  // Integral over Q2		
	  double dndE[3] = {0}; // Effective photon flux
	  double full_int[3] = {0}; // Full e+X --> e+X+V.M. cross section
	  //
	  for( int iEgaInt = 0 ; iEgaInt < 3; ++iEgaInt){    // Loop over the energies for the three points to integrate over Q2
	    //
	    double Q2_min = std::pow(starlightConstants::mel*ega[iEgaInt],2.0)/(_electronEnergy*(_electronEnergy-ega[iEgaInt]));
	    double Q2_max = 4.*_electronEnergy*(_electronEnergy-ega[iEgaInt]);
	    double lnQ2ratio = std::log(Q2_max/Q2_min)/nQ2;
	    double lnQ2_min = std::log(Q2_min);
	    //
	    int q2end = -99;
	    for( int iQ2 = 0 ; iQ2 < nQ2; ++iQ2){     // Integral over photon virtuality
	      //
	      double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
	      //double cms_ega1 = ega[iEgaInt]*(coshy - sinhy) - q2_1*sinhy/(2.*_electronEnergy);
	      //
	      double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
	      //double cms_ega2 = ega[iEgaInt]*(coshy - sinhy) - q2_2*sinhy/(2.*_electronEnergy);
	      //
	      double q2_12 = (q2_2+q2_1)/2.;
	      //double cms_ega12 = ega[iEgaInt]*(coshy - sinhy) - q2_12*sinhy/(2.*_electronEnergy);
	      //if( cms_ega1 < 0 || cms_ega2 < 0 || cms_ega12 < 0 )
	      //	continue;
	      //if( cms_ega1 > _cmsMaxPhotonEnergy || cms_ega2 > _cmsMaxPhotonEnergy || cms_ega12 > _cmsMaxPhotonEnergy)
	      //continue;
	      // Integrating the effectiv photon flux
	      dndE[iEgaInt] +=(q2_2-q2_1)*( getcsgA_Q2_dep(q2_1)*photonFlux(ega[iEgaInt],q2_1)
					    +getcsgA_Q2_dep(q2_2)*photonFlux(ega[iEgaInt],q2_2)
					    +4.*getcsgA_Q2_dep(q2_12)*photonFlux(ega[iEgaInt],q2_12) );
	      //
	      full_int[iEgaInt] += (q2_2-q2_1)*( g(ega[iEgaInt],q2_1)*e_getcsgA(ega[iEgaInt],q2_1,W,beam) 
						 + g(ega[iEgaInt],q2_2)*e_getcsgA(ega[iEgaInt],q2_2,W,beam)
						 + 4.*g(ega[iEgaInt],q2_12)*e_getcsgA(ega[iEgaInt],q2_12,W,beam) );
	      //
	    }
	    //cout<<" Done with energy "<<ega[iEgaInt]<<" finished at "<<q2end<<" / "<<nQ2<<endl;
	    q2end = -99 ;
	    // Finish the Q2 integration for the three end-points (Siumpsons rule)
	    dndE[iEgaInt] = dndE[iEgaInt]/6.; 
	    full_int[iEgaInt] = full_int[iEgaInt]/6.;
	  }	    
	  // Finishing cross-section integral 
	  dR = full_int[0];
	  dR += full_int[1];
	  dR += 4.*full_int[2];
	  dR = dR*(ega[1]-ega[0])/6.;
	  
	  // Finishing integral over the effective photon flux
	  dR2 = dndE[0];
	  dR2 += dndE[1];
	  dR2 += 4.*dndE[2];
	  dR2 = dR2*(ega[1]-ega[0])/6.;
	  //
	  int_r = int_r+dR;
	  int_r2 = int_r2 + dR2;
	}
	//
	double csga_int = 0 ; 
	double ratio  = std::log(_cmsMaxPhotonEnergy/_cmsMinPhotonEnergy)/nEgamma;
	double ega_min  = std::log(_cmsMinPhotonEnergy);
	csga_int = 0;
	cout<<_cmsMinPhotonEnergy<<" - "<<_cmsMaxPhotonEnergy<<endl;
	for( int ii =0 ; ii<nEgamma; ++ii){
	  double ega1 = exp(ega_min + iEgamma*ratio );
	  double ega2 = exp(ega_min + (1+iEgamma)*ratio );
	  double ega12 = 0.5*(ega2+ega1);
	  csga_int += (ega2-ega1)*(ega1*getcsgA(ega1,W,beam) 
				   + ega2*getcsgA(ega2,W,beam) 
				   +4.*ega12*getcsgA(ega12,W,beam))/6.;
	}
	//
	//
	//int_r = int_r/int_r2;
	//
	cout<<int_r*1E4<<" / "<<int_r2*1E4<<endl;
	cout<<endl;      
	printCrossSection(" Total cross section: ",int_r);
	printCrossSection(" Total effective photon flux ", int_r2);
	printCrossSection(" gamma+X --> VM+X ", int_r/int_r2);
       	printCrossSection(" Just csga ", csga_int);
	//
	cout<<endl;
	setPhotonNucleusSigma(0.01*int_r);
	#ifdef _makeGammaPQ2_
	makeGammaPQ2dependence();
	#endif
}


//______________________________________________________________________________
void
e_narrowResonanceCrossSection::makeGammaPQ2dependence()
{
	// This subroutine calculates the Q2 dependence of 
        // gamma+X -> VM + X cross section for a narrow resonance
  
        int const nQ2bins = 19;
	double const q2Edge[nQ2bins+1] = { 0.001,1.,2.,3., 4.,5.,
					   6.,7.,8.,9.,10.,
					   11.,12.,13.,14.,15.,
					   20.,30.,40.,50.};
	//
	double full_x_section[nQ2bins] = {0};
	double effective_flux[nQ2bins] = {0};
	double gamma_x_section[nQ2bins] = {0};
	//
	ofstream  w_file, y_file, q2_file;
	//
	w_file.open("estarlight_gammap_vs_w.csv");
	y_file.open("estarlight_gammap_vs_y.csv");
	q2_file.open("estarlight_gammap_vs_q2.csv");
	//
        double W, dEgamma, minEgamma;
	double ega[3] = {0};
	double int_r,dR;
	double int_r2, dR2;
	int    iEgamma, nEgamma,beam;
	int_r = 0;
	int_r2 = 0;
	//Integration is done with exponential steps, in target frame
	//nEgamma = _VMnumEgamma;
	nEgamma = 1000;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy)/nEgamma;
	minEgamma = std::log(_targetMinPhotonEnergy);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	//
        printf(" gamma+nucleon threshold (CMS): %e GeV \n", _cmsMinPhotonEnergy);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();

	if( A_2 == 0 && A_1 >= 1 ){
	  // pA, first beam is the nucleus and is in this case the target  
	  beam = 1; 
	} else if( A_1 ==0 && A_2 >= 1){       
	  // photon energy in Target frame
	  beam = 2; 
	} else {
	  // photon energy in Target frame
	  beam = 2; 
	}
        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	//double target_cm = -acosh(_beamLorentzGamma);
	double target_cm = acosh(_target_beamLorentz);
	double coshy = cosh(target_cm);
	double sinhy = sinh(target_cm);
	int nQ2 = 500;
        printf(" gamma+nucleon threshold (Target): %e GeV \n", _targetMinPhotonEnergy);
	for( int iQ2Bin  = 0 ; iQ2Bin < nQ2bins; ++iQ2Bin){
	  for(iEgamma = 0 ; iEgamma < nEgamma; ++iEgamma){    // Integral over photon energy
	    // Target frame energies
	    ega[0] = exp(minEgamma + iEgamma*dEgamma );
	    ega[1] = exp(minEgamma + (iEgamma+1)*dEgamma );
	    ega[2] = 0.5*(ega[0]+ega[1]);
	    //
	    //
	    // Integral over Q2		
	    double dndE[3] = {0}; // Effective photon flux
	    double full_int[3] = {0}; // Full e+X --> e+X+V.M. cross section
	    double dsigmadE[3] = {0};
	    //
	    for( int iEgaInt = 0 ; iEgaInt < 3; ++iEgaInt){    // Loop over the energies for the three points to integrate over Q2	   
	      //
	      double Q2_min = q2Edge[iQ2Bin];
	      double Q2_max = q2Edge[iQ2Bin+1];
	      double lnQ2ratio = std::log(Q2_max/Q2_min)/nQ2;
	      double lnQ2_min = std::log(Q2_min);
	      for( int iQ2 = 0 ; iQ2 < nQ2; ++iQ2){     // Integral over photon virtuality
		//
		double q2_1 = exp( lnQ2_min + iQ2*lnQ2ratio);	    
		//		double x_1 =  lnQ2_min + iQ2*lnQ2ratio;
		double cms_ega1 = ega[iEgaInt]*(coshy - sinhy) - q2_1*sinhy/(2.*_electronEnergy);
		//
		double q2_2 = exp( lnQ2_min + (iQ2+1)*lnQ2ratio);
		//		double x_2 =  lnQ2_min + (1+iQ2)*lnQ2ratio;
		double cms_ega2 = ega[iEgaInt]*(coshy - sinhy) - q2_2*sinhy/(2.*_electronEnergy);
		//	      
		double q2_12 = (q2_2+q2_1)/2.;
		double cms_ega12 = ega[iEgaInt]*(coshy - sinhy) - q2_12*sinhy/(2.*_electronEnergy);
		//		double x_12 = 0.5*(x_2-x_1);
		//q2_12 = exp(x_12);
		//
		//double limit = starlightConstants::protonMass*(2*ega[iEgaInt]+protonMass);
		//if( q2_1 > limit || q2_2 > limit)
		//  continue;
		if( cms_ega1 < 0 || cms_ega2 < 0 || cms_ega12 < 0 ) //lazy solution to the integrals
		  //break;
		  continue;
		// Integrating the effectiv photon flux
		dndE[iEgaInt] +=(q2_2-q2_1)*( photonFlux(ega[iEgaInt],q2_1)
					      +photonFlux(ega[iEgaInt],q2_2)
					      +4.*photonFlux(ega[iEgaInt],q2_12) );
	      }
	      // Finish the Q2 integration for the three end-points (Siumpsons rule)
	      dndE[iEgaInt] = dndE[iEgaInt]/6.; 
	      full_int[iEgaInt] = full_int[iEgaInt]/6.;
	      dsigmadE[iEgaInt] = dsigmadE[iEgaInt]/6.;
	    }	    
	    // Finishing cross-section integral 
	    dR = full_int[0];
	    dR += full_int[1];
	    dR += 4.*full_int[2];
	    dR = dR*(ega[1]-ega[0])/6.;
	    
	    // Finishing integral over the effective photon flux
	    dR2 = dndE[0];
	    dR2 += dndE[1];
	    dR2 += 4.*dndE[2];
	    dR2 = dR2*(ega[1]-ega[0])/6.;
	    //
	    //	    int_r = int_r+dR;
	    //	    int_r2 = int_r2 + dR2;
	    full_x_section[iQ2Bin] += dR;
	    effective_flux[iQ2Bin] += dR2;	      
	    gamma_x_section[iQ2Bin] += (ega[1]-ega[0])*( dsigmadE[0] + dsigmadE[1]+4.*dsigmadE[2])/6.;
	  }
	  cout<<q2Edge[iQ2Bin]<<" - "<<q2Edge[iQ2Bin+1]<<"\t : "<<full_x_section[iQ2Bin]*10000000
	      <<"\t "<<effective_flux[iQ2Bin] <<"\t "<<full_x_section[iQ2Bin]/effective_flux[iQ2Bin]*10000000<<endl;
	  cout<<gamma_x_section[iQ2Bin]/effective_flux[iQ2Bin]*1E7<<endl;
	  q2_file<<(q2Edge[iQ2Bin+1]+q2Edge[iQ2Bin])/2.+W*W<<","<<gamma_x_section[iQ2Bin]/effective_flux[iQ2Bin]*1E7<<endl;
	}
	//
	//
	y_file.close();
	w_file.close();
	q2_file.close();
}

void e_narrowResonanceCrossSection::printCrossSection(const string name, const double x_section)
{
  if (0.01*x_section > 1.){
    cout<< name.c_str() <<0.01*x_section<<" barn."<<endl;
  } else if (10.*x_section > 1.){
    cout<< name.c_str() <<10.*x_section<<" mb."<<endl;
  } else if (10000.*x_section > 1.){
    cout<< name.c_str() <<10000.*x_section<<" microb."<<endl;
  } else if (10000000.*x_section > 1.){
    cout<< name.c_str() <<10000000.*x_section<<" nanob."<<endl;
  } else if (1.E10*x_section > 1.){
    cout<< name.c_str() <<1.E10*x_section<<" picob."<<endl;
  } else {
    cout<< name.c_str() <<1.E13*x_section<<" femtob."<<endl;
  }
  cout<<endl;
}  
