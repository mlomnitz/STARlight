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
	_narrowYmax = inputParametersInstance.maxRapidity();
	_narrowYmin = -1.0*_narrowYmax;
	_narrowNumY = inputParametersInstance.nmbRapidityBins();
	_Ep         = inputParametersInstance.protonEnergy();	
	_electronEnergy = inputParametersInstance.electronEnergy();
	//_target_beamLorentz = inputParametersInstance.beam2LorentzGamma();
	_target_beamLorentz = inputParametersInstance.beamLorentzGamma();
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
  
	double W,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double g_Eg1, g_Eg2, g_Eg12;
	double csgA1,csgA2,csgA12,int_r,dR;
	double Eth;
	int    J,NY,beam;
  
	NY   =  _narrowNumY;
	dY   = (_narrowYmax-_narrowYmin)/double(NY);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	//Lomnitz old used for XX
	Eth=0.5*(((W+protonMass)*(W+protonMass)-
	          protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
	// Adapted for eX
	//Eth=0.5*(((W+starlightConstants::mel)*(W +starlightConstants::mel)-
	//	  starlightConstants::mel*starlightConstants::mel)/(_electronEnergy+sqrt(_electronEnergy*_electronEnergy-starlightConstants::mel*starlightConstants::mel))); 
	// cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
        printf(" gamma+nucleon threshold (CMS): %e GeV \n", Eth);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();
  
 	int_r=0.;

        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	//double target_cm = -acosh(_beamLorentzGamma);
	double target_cm = -acosh(_target_beamLorentz);
	// another - sign from subraction in addition rule
	double exp_target_cm = exp(-target_cm);
	///
        printf(" gamma+nucleon threshold (Target): %e GeV \n", Eth*exp_target_cm);
	for(J=0;J<=(NY-1);J++){
    
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
		double target_ega1, target_ega2, target_ega12;
                if( A_2 == 0 && A_1 >= 1 ){
                  // pA, first beam is the nucleus and is in this case the target  
                  ega1  = 0.5*W*exp(-y1);
                  ega2  = 0.5*W*exp(-y2);
                  ega12 = 0.5*W*exp(-y12);
                  beam = 1; 
                } else if( A_1 ==0 && A_2 >= 1){
                  // pA, second beam is the nucleus and is in this case the target 
                  ega1  = 0.5*W*exp(y1);
                  ega2  = 0.5*W*exp(y2);
                  ega12 = 0.5*W*exp(y12);
		  // photon energy in Target frame
                  beam = 2; 
                } else {
                  ega1  = 0.5*W*exp(y1);
                  ega2  = 0.5*W*exp(y2);
                  ega12 = 0.5*W*exp(y12);
		  // photon energy in Target frame
                  beam = 2; 
                }
    		  // photon energy in Target frame 
		target_ega1 = ega1*exp_target_cm;
		target_ega2 = ega2*exp_target_cm;
		target_ega12 = ega12*exp_target_cm;
		//
		if(ega1 < Eth || ega2 < Eth)   
			continue;
		if(ega2 > maxPhotonEnergy() || ega1 > maxPhotonEnergy() ) 
			continue;
		//			
		g_Eg1 = integrated_Q2_dep(target_ega1);
		csgA1=getcsgA(ega1,W,beam);
		//         >> Middle Point                      =====>>>
		g_Eg12 = integrated_Q2_dep(target_ega12);
		csgA12=getcsgA(ega12,W,beam);         
		//         >> Second Point                      =====>>>
		g_Eg2 = integrated_Q2_dep(target_ega2);
		csgA2=getcsgA(ega2,W,beam);
		//testing - should be correct
		dR  = g_Eg1*csgA1*target_ega1;
		dR  = dR + 4.*g_Eg12*csgA12*target_ega12;
		dR  = dR + g_Eg2*csgA2*target_ega2;
		dR  = dR*(dY/6.);
		//>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
		//dR  = ega1*g_Eg1*csgA1;
		//dR  = dR + 4.*ega12*g_Eg12*csgA12;
		//dR  = dR + ega2*g_Eg2*csgA2;
		//dR  = dR*(dY/6.);
		
		// cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12,beam)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;

		int_r = int_r+dR;
	}


	cout<<endl;
	if (0.01*int_r > 1.){
	  cout<< " Total cross section: "<<0.01*int_r<<" barn."<<endl;
	} else if (10.*int_r > 1.){
	  cout<< " Total cross section: " <<10.*int_r<<" mb."<<endl;
        } else if (10000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000.*int_r<<" microb."<<endl;
        } else if (10000000.*int_r > 1.){
	  cout<< " Total cross section: " <<10000000.*int_r<<" nanob."<<endl;
        } else if (1.E10*int_r > 1.){
	  cout<< " Total cross section: "<<1.E10*int_r<<" picob."<<endl;
        } else {
	  cout<< " Total cross section: " <<1.E13*int_r<<" femtob."<<endl;
        }
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
	double const q2Edge[nQ2bins+1] = { 0.,1.,2.,3., 4.,5.,
					   6.,7.,8.,9.,10.,
					   11.,12.,13.,14.,15.,
					   20.,30.,40.,50.};

	double W,dY;
	double y1,y2,y12,ega1,ega2,ega12;
	double q2_cor1, q2_cor12, q2_cor2;
	double csgA1,csgA2,csgA12,int_r,dR;
	double dR2;
	double Eth;
	int    J,NY,beam;
  
	NY   =  _narrowNumY;
	dY   = (_narrowYmax-_narrowYmin)/double(NY);
  
	cout<<" Using Narrow Resonance ..."<<endl;
  
	W = getChannelMass();
	//Lomnitz old used for XX
	Eth=0.5*(((W+protonMass)*(W+protonMass)-
	          protonMass*protonMass)/(_Ep+sqrt(_Ep*_Ep-protonMass*protonMass)));
	// Adapted for eX
	//Eth=0.5*(((W+starlightConstants::mel)*(W +starlightConstants::mel)-
	//	  starlightConstants::mel*starlightConstants::mel)/(_electronEnergy+sqrt(_electronEnergy*_electronEnergy-starlightConstants::mel*starlightConstants::mel))); 
	// cout<<" gamma+nucleon  Threshold: "<<Eth<<endl;
        printf(" gamma+nucleon threshold: %e GeV \n", Eth);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();
  

        // Do this first for the case when the first beam is the photon emitter 
        // The variable beam (=1,2) defines which nucleus is the target 
	// Target beam ==2 so repidity is negative. Can generalize later
	///
	cout<<" Lomnitz debug :: sigma_gamma_p --> VM_p "<<endl;
	cout<<" Q2+MV2 \t \t"<<" sigma_gamma_p --> VM_p (nanob)"<<endl;
	double target_cm = acosh(_target_beamLorentz);
	// another - sign from subraction in addition rule
	double exp_target_cm = exp(-target_cm);
	double int_r2;
	for( int iQ2 = 0 ; iQ2 < nQ2bins; ++iQ2){
	  int_r=0.;
	  int_r2=0.;
	  for(J=0;J<=(NY-1);J++){
	    
	    y1  = _narrowYmin + double(J)*dY;
	    y2  = _narrowYmin + double(J+1)*dY;
	    y12 = 0.5*(y1+y2);
	    double target_ega1, target_ega2, target_ega12;
	    if( A_2 == 0 && A_1 >= 1 ){
	      // pA, first beam is the nucleus and is in this case the target  
	      ega1  = 0.5*W*exp(-y1);
	      ega2  = 0.5*W*exp(-y2);
	      ega12 = 0.5*W*exp(-y12);
	      beam = 1; 
	    } else if( A_1 ==0 && A_2 >= 1){
	      // pA, second beam is the nucleus and is in this case the target 
	      ega1  = 0.5*W*exp(y1);
	      ega2  = 0.5*W*exp(y2);
	      ega12 = 0.5*W*exp(y12);
	      // photon energy in Target frame
	      beam = 2; 
	    } else {
	      ega1  = 0.5*W*exp(y1);
	      ega2  = 0.5*W*exp(y2);
	      ega12 = 0.5*W*exp(y12);
	      // photon energy in Target frame
	      beam = 2; 
	    }
	    //
	    if(ega1 < Eth || ega2 < Eth)   
	      continue;
	    if(ega2 > maxPhotonEnergy() || ega1 > maxPhotonEnergy() ) 
	      continue;
	    target_ega1 = ega1*exp_target_cm;
	    target_ega12 = ega12*exp_target_cm;
	    target_ega2 = ega2*exp_target_cm;
	    //cout<<"Nortmalizations "<<integrated_x_section(ega1,0,50)<<endl;
	    //		
	    csgA1=getcsgA(ega1,W,beam);
	    double full_range_1 = integrated_x_section(target_ega1);
	    //         >> Middle Point                      =====>>>
	    csgA12=getcsgA(ega12,W,beam);         
	    double full_range_12 = integrated_x_section(target_ega12);
	    //         >> Second Point                      =====>>>
	    csgA2=getcsgA(ega2,W,beam);
	    double full_range_2 = integrated_x_section(target_ega2);
	    //
		
	    
	    //q2_cor1 = integrated_x_section(ega1,q2Edge[iQ2],q2Edge[iQ2+1]);
	    //q2_cor12 = integrated_x_section(ega12,q2Edge[iQ2],q2Edge[iQ2+1]);
	    //q2_cor2 = integrated_x_section(ega2,q2Edge[iQ2],q2Edge[iQ2+1]);
	    q2_cor1 = getcsgA_Q2_dep( (q2Edge[iQ2+1] + q2Edge[iQ2])/2. );
	    q2_cor12 = q2_cor1;
	    q2_cor2 = q2_cor1;
	    //testing - should be correct
	    dR  = csgA1*q2_cor1;
	    dR  = dR + 4.*csgA12*q2_cor12;
	    dR  = dR + csgA2*q2_cor2;
	    dR  = dR*(dY/6.);
	    //
	    dR2  = csgA1*full_range_1;
	    dR2  = dR2 + 4.*csgA12*full_range_12;
	    dR2  = dR2 + csgA2*full_range_2;
	    dR2  = dR2*(dY/6.);
	    //>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
	    //dR  = ega1*g_Eg1*csgA1;
	    //dR  = dR + 4.*ega12*g_Eg12*csgA12;
	    //dR  = dR + ega2*g_Eg2*csgA2;
	    //dR  = dR*(dY/6.);
	    
	    // cout<<" y: "<<y12<<" egamma: "<<ega12<<" flux: "<<photonFlux(ega12,beam)<<" sigma_gA: "<<10000000.*csgA12<<" dsig/dy (microb): "<<10000.*dR/dY<<endl;
	    
	    int_r = int_r+dR;
	    int_r2 = int_r2 +dR2; 
	  }
	  //cout<<(q2Edge[iQ2+1]+q2Edge[iQ2])/2.+W*W<<" ,  "<<10000000.*int_r/(q2Edge[iQ2+1]-q2Edge[iQ2])<<endl;
	  if( iQ2 ==0 )
	    cout<<"Full range "<<int_r2*10000000<<endl;
	  cout<<(q2Edge[iQ2+1]+q2Edge[iQ2])/2.+W*W<<" ,  "<<10000000.*int_r<<endl;
	}
}

