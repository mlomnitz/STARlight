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
#include <fstream>
#include <cmath>

#include "inputParameters.h"
#include "beambeamsystem.h"
#include "beam.h"
#include "starlightconstants.h"
#include "nucleus.h"
#include "bessel.h"
#include "gammaeluminosity.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
photonElectronLuminosity::photonElectronLuminosity(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem)
  : photonNucleusCrossSection(inputParametersInstance, bbsystem)
  ,_protonEnergy(inputParametersInstance.protonEnergy())
  ,_electronEnergy(inputParametersInstance.electronEnergy())
  ,_beamLorentzGamma(inputParametersInstance.beamLorentzGamma())
  ,_baseFileName(inputParametersInstance.baseFileName())
  ,_maxW(inputParametersInstance.maxW())
  ,_minW(inputParametersInstance.minW())
  ,_nmbWBins(inputParametersInstance.nmbWBins())
  ,_maxRapidity(inputParametersInstance.maxRapidity())
  ,_nmbRapidityBins(inputParametersInstance.nmbRapidityBins())
  ,_minGammaQ2(inputParametersInstance.minGammaQ2())
  ,_maxGammaQ2(inputParametersInstance.maxGammaQ2())
  ,_nmbGammaQ2Bins(inputParametersInstance.nmbGammaQ2Bins()) 
  ,_productionMode(inputParametersInstance.productionMode())
  ,_beamBreakupMode(inputParametersInstance.beamBreakupMode())
{
  cout <<"Creating Luminosity Tables."<<endl;
  photonNucleusDifferentialLuminosity();
  cout <<"Luminosity Tables created."<<endl;
}


//______________________________________________________________________________
photonElectronLuminosity::~photonElectronLuminosity()
{ }


//______________________________________________________________________________
void photonElectronLuminosity::photonNucleusDifferentialLuminosity()
{
  double W,dW,dY;
  double Egamma,Y;
  double testint;
  // 
  double f_WY, g_E;
  double csgA;
  double C;  
  int beam; 

  std::string wyFileName;
  wyFileName = _baseFileName +".txt";

  ofstream wylumfile;
  wylumfile.precision(15);
  
  double  bwnorm,Eth;

  dW = (_wMax-_wMin)/_nWbins;
  dY  = (_yMax-(-1.0)*_yMax)/_nYbins;
    
  // Write the values of W used in the calculation to slight.txt.  
  wylumfile.open(wyFileName.c_str());
  wylumfile << getbbs().beam1().Z() <<endl;
  wylumfile << getbbs().beam1().A() <<endl;
  wylumfile << getbbs().beam2().Z() <<endl;
  wylumfile << getbbs().beam2().A() <<endl;
  wylumfile << _beamLorentzGamma <<endl;
  wylumfile << _maxW <<endl;
  wylumfile << _minW <<endl;
  wylumfile << _nmbWBins <<endl;
  wylumfile << _maxRapidity <<endl;
  wylumfile << _nmbRapidityBins <<endl;
  wylumfile << _productionMode <<endl;
  wylumfile << _beamBreakupMode <<endl;
  wylumfile << starlightConstants::deuteronSlopePar <<endl;
  
  //     Normalize the Breit-Wigner Distribution and write values of W to slight.txt
  testint=0.0;
  //Grabbing default value for C in the breit-wigner calculation
  C=getDefaultC();
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {
    W = _wMin + double(i)*dW + 0.5*dW;
    testint = testint + breitWigner(W,C)*dW;
    wylumfile << W << endl;
  }
  bwnorm = 1./testint;
  
  //     Write the values of Y used in the calculation to slight.txt.
  for(unsigned int i = 0; i <= _nYbins - 1; ++i) {
    Y = -1.0*_yMax + double(i)*dY + 0.5*dY;
    wylumfile << Y << endl;
  }
    
  //Eth=0.5*(((_wMin+starlightConstants::mel)*(_wMin +starlightConstants::mel)-starlightConstants::mel*starlightConstants::mel)/(_electronEnergy+sqrt(_electronEnergy*_electronEnergy-starlightConstants::mel*starlightConstants::mel))); 
  Eth=0.5*(((W+protonMass)*(W+protonMass)-
	    protonMass*protonMass)/(_protonEnergy+sqrt(_protonEnergy*_protonEnergy-protonMass*protonMass)));

  int A_1 = getbbs().beam1().A(); 
  int A_2 = getbbs().beam2().A();

  // Do this first for the case when the first beam is the photon emitter 
  // Treat pA separately with defined beams 
  // The variable beam (=1,2) defines which nucleus is the target 
  for(unsigned int i = 0; i <= _nWbins - 1; ++i) {

    W = _wMin + double(i)*dW + 0.5*dW;
    
    for(unsigned int j = 0; j <= _nYbins - 1; ++j) { 

      Y = -1.0*_yMax + double(j)*dY + 0.5*dY;

      if( A_2 == 0 && A_1 != 0 ){
        // eA, first beam is the nucleus and is in this case the target  
        Egamma = 0.5*W*exp(-Y); 
        beam = 1; 
      } else if( A_1 ==0 && A_2 != 0){
        // pA, second beam is the nucleus and is in this case the target 
        Egamma = 0.5*W*exp(Y); 
        beam = 2; 
      } else {
        Egamma = 0.5*W*exp(Y);        
        beam = 2; 
      }

      f_WY = 0.; 
      g_E = 0;

      if( Egamma > Eth && Egamma < maxPhotonEnergy() ){

	csgA=getcsgA(Egamma,W,beam);
        f_WY = Egamma*csgA*breitWigner(W,bwnorm);
	g_E = integrated_Q2_dep(Egamma);
      }

      wylumfile << f_WY << endl;
      wylumfile << g_E << endl;

    }
  }

  wylumfile << bwnorm << endl;
  wylumfile.close();
   
}



