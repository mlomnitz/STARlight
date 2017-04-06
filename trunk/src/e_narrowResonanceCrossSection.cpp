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
        printf(" gamma+nucleon threshold: %e GeV \n", Eth);

        int A_1 = getbbs().beam1().A(); 
        int A_2 = getbbs().beam2().A();
  
 	int_r=0.;

        // Do this first for the case when the first beam is the photon emitter 
        // Treat pA separately with defined beams 
        // The variable beam (=1,2) defines which nucleus is the target 
	for(J=0;J<=(NY-1);J++){
    
		y1  = _narrowYmin + double(J)*dY;
		y2  = _narrowYmin + double(J+1)*dY;
		y12 = 0.5*(y1+y2);
    
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
                  beam = 2; 
                } else {
                  ega1  = 0.5*W*exp(y1);
                  ega2  = 0.5*W*exp(y2);
                  ega12 = 0.5*W*exp(y12);
                  beam = 2; 
                }
    
		if(ega1 < Eth || ega2 < Eth)   
			continue;
		if(ega2 > maxPhotonEnergy() || ega1 > maxPhotonEnergy() ) 
			continue;

		//			
		g_Eg1 = integrated_Q2_dep(ega1);
		csgA1=getcsgA(ega1,W,beam);
		//         >> Middle Point                      =====>>>
		g_Eg12 = integrated_Q2_dep(ega12);
		csgA12=getcsgA(ega12,W,beam);         
		//         >> Second Point                      =====>>>
		g_Eg2 = integrated_Q2_dep(ega2);
		csgA2=getcsgA(ega2,W,beam);
		//>> Sum the contribution for this W,Y. The 2 accounts for the 2 beams
		dR  = ega1*g_Eg1*csgA1;
		dR  = dR + 4.*ega12*g_Eg12*csgA12;
		dR  = dR + ega2*g_Eg2*csgA2;
		dR  = dR*(dY/6.);

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
}
