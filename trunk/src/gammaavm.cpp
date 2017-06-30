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
// $Rev:: 277                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-09-14 10:55:55 +0100 #$: date of last commit
//
// Description:
//    Added incoherent t2-> pt2 selection.  Following pp selection scheme
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "gammaavm.h"
//#include "photonNucleusCrossSection.h"
#include "wideResonanceCrossSection.h"
#include "narrowResonanceCrossSection.h"
#include "incoherentVMCrossSection.h"
//
#include "e_narrowResonanceCrossSection.h"
#include "e_wideResonanceCrossSection.h"

using namespace std;


//______________________________________________________________________________
Gammaavectormeson::Gammaavectormeson(const inputParameters& inputParametersInstance, beamBeamSystem& bbsystem):eventChannel(inputParametersInstance, bbsystem), _phaseSpaceGen(0)
{
	_VMNPT=inputParametersInstance.nmbPtBinsInterference();
	_VMWmax=inputParametersInstance.maxW();
	_VMWmin=inputParametersInstance.minW();
	_VMYmax=inputParametersInstance.maxRapidity();
	_VMYmin=-1.*_VMYmax;
	_VMnumw=inputParametersInstance.nmbWBins();
	_VMnumy=inputParametersInstance.nmbRapidityBins();
	_VMnumega=inputParametersInstance.nmbEnergyBins();
	_VMgamma_em=inputParametersInstance.beamLorentzGamma();
	_VMinterferencemode=inputParametersInstance.interferenceEnabled();
	_VMbslope=0.;//Will define in wide/narrow constructor
        _bslopeDef=inputParametersInstance.bslopeDefinition();
	_bslopeVal=inputParametersInstance.bslopeValue();
	_pEnergy= inputParametersInstance.protonEnergy();
	// electron energy in CMS frame
	_eEnergy= inputParametersInstance.electronEnergy();
	_VMpidtest=inputParametersInstance.prodParticleType();
	_VMptmax=inputParametersInstance.maxPtInterference();
	_VMdpt=inputParametersInstance.ptBinWidthInterference();
        _ProductionMode=inputParametersInstance.productionMode();
	_targetMaxPhotonEnergy=inputParametersInstance.targetMaxPhotonEnergy();
	_targetMinPhotonEnergy=inputParametersInstance.targetMinPhotonEnergy();
	// Now saving the photon energy limits
	_cmsMaxPhotonEnergy=inputParametersInstance.cmsMaxPhotonEnergy();
	_cmsMinPhotonEnergy=inputParametersInstance.cmsMinPhotonEnergy();
	_beamLorentzGamma = inputParametersInstance.beamLorentzGamma();
        N0 = 0; N1 = 0; N2 = 0; 
	if (_VMpidtest == starlightConstants::FOURPRONG){
	  // create n-body phase-spage generator
	  _phaseSpaceGen = new nBodyPhaseSpaceGen(_randy);
	}
	if(_ProductionMode == 11 || _ProductionMode == 12) //Need to define later, for full eSTARlight
	  _dummy_pncs = new photonNucleusCrossSection(inputParametersInstance, bbsystem);
	//
	const double r_04_00 = 0.674;
	const double cos_delta = 0.925;
	for( int ii = 0; ii < 100; ++ii){ //epsilon 0-1
	  double epsilon = 0.01*ii;
	  const double R = (1./epsilon)*r_04_00/(1.-r_04_00);
	  for(int jj = 0; jj < 200; ++jj){ //psi 0 - 2pi
	    double psi = jj*starlightConstants::pi/100.;
	    double max_bin;
	    for( int kk = 0; kk < 200; ++kk){ //temp
	      double theta = kk*starlightConstants::pi/100.;
	      // Fin max
	      double this_test = std::pow(std::sin(theta),2.)*(1+epsilon*cos(2.*psi)) + 2.*epsilon*R*std::pow(std::cos(theta),2.)
		+std::sqrt(2.*epsilon*(1+epsilon))*cos_delta*std::sin(2.*theta)*std::cos(psi);
	      if(this_test >  max_bin)
		max_bin = this_test;
	    }
	    _angular_max[ii][jj] = max_bin;
	  }
	}
}


//______________________________________________________________________________
Gammaavectormeson::~Gammaavectormeson()
{
	if (_phaseSpaceGen)
		delete _phaseSpaceGen;
	if (_dummy_pncs)
	  delete _dummy_pncs;
}


//______________________________________________________________________________
void Gammaavectormeson::pickwy(double &W, double &Y)
{
        double dW, dY, xw,xy,xtest,btest;
	int  IW,IY;
  
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	dY = (_VMYmax-_VMYmin)/double(_VMnumy);
  
 L201pwy:

	xw = _randy.Rndom();
	W = _VMWmin + xw*(_VMWmax-_VMWmin);

	if (W < 2 * starlightConstants::pionChargedMass)
		goto L201pwy;
  
	IW = int((W-_VMWmin)/dW);
	xy = _randy.Rndom();
	Y = _VMYmin + xy*(_VMYmax-_VMYmin);
	IY = int((Y-_VMYmin)/dY); 
	xtest = _randy.Rndom();

	if( xtest > _Farray[IW][IY] )
		goto L201pwy;

        N0++; 
	// Determine the target nucleus 
	// For pA this is given, for all other cases use the relative probabilities in _Farray1 and _Farray2 
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
           if( _ProductionMode == 2 || _ProductionMode ==3){
	     _TargetBeam = 2;
	   } else {
             _TargetBeam = 1;
           }
        } else if(  _bbs.beam1().A() != 1 && _bbs.beam2().A()==1 ){
           if( _ProductionMode == 2 || _ProductionMode ==3){
	     _TargetBeam = 1;
	   } else {
             _TargetBeam = 2;
           }
        } else {
          btest = _randy.Rndom();
	  if ( btest < _Farray1[IW][IY]/_Farray[IW][IY] ){
            _TargetBeam = 2;
            N2++;
          }  else {
            _TargetBeam = 1;
            N1++; 
          }
        }
}         


//______________________________________________________________________________                                               
void Gammaavectormeson::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                     double  W,
                                     double  px0, double  py0, double  pz0,
                                     double& px1, double& py1, double& pz1,
                                     double& px2, double& py2, double& pz2,
                                     int&    iFbadevent)
{
	// This routine decays a particle into two particles of mass mdec,
	// taking spin into account

	double pmag;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double mdec=0.0;
	double E1=0.0,E2=0.0;

	//    set the mass of the daughter particles
	mdec=getDaughterMass(ipid);

	//  calculate the magnitude of the momenta
	if(W < 2*mdec){
		cout<<" ERROR: W="<<W<<endl;
		iFbadevent = 1;
		return;
	}
	pmag = sqrt(W*W/4. - mdec*mdec);
  
	//  pick an orientation, based on the spin
	//  phi has a flat distribution in 2*pi
	phi = _randy.Rndom()*2.*starlightConstants::pi;
                                                                                                                
	//  find theta, the angle between one of the outgoing particles and
	//  the beamline, in the frame of the two photons

	theta=getTheta(ipid);
 
	//  compute unboosted momenta
	px1 = sin(theta)*cos(phi)*pmag;
	py1 = sin(theta)*sin(phi)*pmag;
	pz1 = cos(theta)*pmag;
	px2 = -px1;
	py2 = -py1;
	pz2 = -pz1;

	Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
	E1 = sqrt(mdec*mdec+px1*px1+py1*py1+pz1*pz1);
	E2 = sqrt(mdec*mdec+px2*px2+py2*py2+pz2*pz2);

	betax = -(px0/Ecm);
	betay = -(py0/Ecm);
	betaz = -(pz0/Ecm);

	transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
	transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);

	if(iFbadevent == 1)
	   return;

}


//______________________________________________________________________________                                               
void Gammaavectormeson::twoBodyDecay(starlightConstants::particleTypeEnum &ipid,
                                     double  W,
                                     double  px0, double  py0, double  pz0,
				     double  e_phi, double polarization,
                                     double& px1, double& py1, double& pz1,
                                     double& px2, double& py2, double& pz2,
                                     int&    iFbadevent)
{
	// This routine decays a particle into two particles of mass mdec,
	// taking spin into account

	double pmag;
	double phi,theta,Ecm;
	double betax,betay,betaz;
	double mdec=0.0;
	double E1=0.0,E2=0.0;

	//    set the mass of the daughter particles
	mdec=getDaughterMass(ipid);

	//  calculate the magnitude of the momenta
	if(W < 2*mdec){
		cout<<" ERROR: W="<<W<<endl;
		iFbadevent = 1;
		return;
	}
	pmag = sqrt(W*W/4. - mdec*mdec);
  
	//  pick an orientation, based on the spin
	//  phi has a flat distribution in 2*pi
	pair<double,double>* angles = getThetaPsi(e_phi,polarization);
	theta = angles->first;
	phi = angles->second;

	//phi = _randy.Rndom()*2.*starlightConstants::pi;
                                                                                                                
	//  find theta, the angle between one of the outgoing particles and
	//  the beamline, in the frame of the two photons

	//theta=getTheta(ipid);
 
	//  compute unboosted momenta
	px1 = sin(theta)*cos(phi)*pmag;
	py1 = sin(theta)*sin(phi)*pmag;
	pz1 = cos(theta)*pmag;
	px2 = -px1;
	py2 = -py1;
	pz2 = -pz1;

	Ecm = sqrt(W*W+px0*px0+py0*py0+pz0*pz0);
	E1 = sqrt(mdec*mdec+px1*px1+py1*py1+pz1*pz1);
	E2 = sqrt(mdec*mdec+px2*px2+py2*py2+pz2*pz2);

	betax = -(px0/Ecm);
	betay = -(py0/Ecm);
	betaz = -(pz0/Ecm);

	transform (betax,betay,betaz,E1,px1,py1,pz1,iFbadevent);
	transform (betax,betay,betaz,E2,px2,py2,pz2,iFbadevent);

	if(iFbadevent == 1)
	   return;

}


//______________________________________________________________________________                                               
// decays a particle into four particles with isotropic angular distribution
bool Gammaavectormeson::fourBodyDecay
(starlightConstants::particleTypeEnum& ipid,
 const double                  ,           // E (unused)
 const double                  W,          // mass of produced particle
 const double*                 p,          // momentum of produced particle; expected to have size 3
 lorentzVector*                decayVecs,  // array of Lorentz vectors of daughter particles; expected to have size 4
 int&                          iFbadevent)
{
	const double parentMass = W;

	// set the mass of the daughter particles
	const double daughterMass = getDaughterMass(ipid);
	if (parentMass < 4 * daughterMass){
		cout << " ERROR: W=" << parentMass << " GeV too small" << endl;
		iFbadevent = 1;
		return false;
	}

	// construct parent four-vector
	const double        parentEnergy = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
	                                        + parentMass * parentMass);
	const lorentzVector parentVec(p[0], p[1], p[2], parentEnergy);

	// setup n-body phase-space generator
	assert(_phaseSpaceGen);
	static bool firstCall = true;
	if (firstCall) {
		const double m[4] = {daughterMass, daughterMass, daughterMass, daughterMass};
		_phaseSpaceGen->setDecay(4, m);
		// estimate maximum phase-space weight
		_phaseSpaceGen->setMaxWeight(1.01 * _phaseSpaceGen->estimateMaxWeight(_VMWmax));
		firstCall = false;
	}

	// generate phase-space event
	if (!_phaseSpaceGen->generateDecayAccepted(parentVec))
		return false;

	// set Lorentzvectors of decay daughters
	for (unsigned int i = 0; i < 4; ++i)
		decayVecs[i] = _phaseSpaceGen->daughter(i);
	return true;
}


//______________________________________________________________________________
double Gammaavectormeson::getDaughterMass(starlightConstants::particleTypeEnum &ipid)
{
	//This will return the daughter particles mass, and the final particles outputed id...
	double mdec=0.;
  
	switch(_VMpidtest){
	case starlightConstants::RHO:
	case starlightConstants::RHOZEUS:
	case starlightConstants::FOURPRONG:
	case starlightConstants::OMEGA:
		mdec = starlightConstants::pionChargedMass;
		ipid = starlightConstants::PION;
		break;
	case starlightConstants::PHI:
		mdec = starlightConstants::kaonChargedMass;
		ipid = starlightConstants::KAONCHARGE;
		break;
	case starlightConstants::JPSI:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 
	case starlightConstants::JPSI2S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break; 
	case starlightConstants::JPSI2S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break; 

	case starlightConstants::JPSI2S:
	case starlightConstants::UPSILON:
	case starlightConstants::UPSILON2S:
	case starlightConstants::UPSILON3S:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;
		break;
	case starlightConstants::UPSILON_ee:
	case starlightConstants::UPSILON2S_ee:
	case starlightConstants::UPSILON3S_ee:
		mdec = starlightConstants::mel;
		ipid = starlightConstants::ELECTRON;
		break;
	case starlightConstants::UPSILON_mumu:
	case starlightConstants::UPSILON2S_mumu:
	case starlightConstants::UPSILON3S_mumu:
		mdec = starlightConstants::muonMass;
		ipid = starlightConstants::MUON;   
		break;
	default: cout<<"No daughtermass defined, gammaavectormeson::getdaughtermass"<<endl;
	}
  
	return mdec;
}


//______________________________________________________________________________
double Gammaavectormeson::getTheta(starlightConstants::particleTypeEnum ipid)
{
	//This depends on the decay angular distribution
	//Valid for rho, phi, omega.
	double theta=0.;
	double xtest=0.;
	double dndtheta=0.;

 L200td:
    
	theta = starlightConstants::pi*_randy.Rndom();
	xtest = _randy.Rndom();
	//  Follow distribution for helicity +/-1
	//  Eq. 19 of J. Breitweg et al., Eur. Phys. J. C2, 247 (1998)
	//  SRK 11/14/2000
  
	switch(ipid){
	  
	case starlightConstants::MUON:
	case starlightConstants::ELECTRON:
		//primarily for upsilon/j/psi.  VM->ee/mumu
		dndtheta = sin(theta)*(1.+((cos(theta))*(cos(theta))));
		break;
    
	case starlightConstants::PION:
	case starlightConstants::KAONCHARGE:
		//rhos etc
		dndtheta= sin(theta)*(1.-((cos(theta))*(cos(theta))));
		break;
    
	default: cout<<"No proper theta dependence defined, check gammaavectormeson::gettheta"<<endl;
	}//end of switch
  
	if(xtest > dndtheta)
		goto L200td;
  
	return theta;
  
}


//______________________________________________________________________________
pair<double,double>* Gammaavectormeson::getThetaPsi(double const e_plane_angle, double const epsilon)
{
  int i_plane[2] = { floor(epsilon*100.) , ceil(epsilon*100.) }; 
  // constants
  //using r_{1-1}^1 = 0.122, r_{00}^{04}=0.674 and cos(delta) = 0.925 from HERA arXiv:hep-ex/9902019
  const double r_1_1_1 = 0.122;
  const double r_04_00 = 0.674;
  const double R = (1./epsilon)*r_04_00/(1.-r_04_00);
  const double cos_delta = 0.925;
  //sample intgrated function for psi
  double theta, psi;
  double xtest = 999.;
  double this_test = 0;
  //
  while( xtest > this_test){
    psi = 2.*starlightConstants::pi*_randy.Rndom();
    this_test = (1.+2.*epsilon*r_1_1_1*std::cos(2.*psi))/(2.*starlightConstants::pi);
    xtest = _randy.Rndom();
    if( xtest > this_test )
      continue;
    int i_psi[2] = { floor( psi*100./starlightConstants::pi ), ceil( psi*100./starlightConstants::pi ) }; 
    // Got good psi, now sample for costheta
    theta = 2.*starlightConstants::pi*_randy.Rndom();
    //quick test
    this_test = std::pow(std::sin(theta),2.)*(1+epsilon*cos(2.*psi)) + 2.*epsilon*R*std::pow(std::cos(theta),2.)
      +std::sqrt(2.*epsilon*(1+epsilon))*cos_delta*std::sin(2.*theta)*std::cos(psi);
    this_test = this_test /_angular_max[ i_plane[0] ][ i_psi[0] ];

    xtest = _randy.Rndom();
    /*if( xtest < this_test)
      break;*/
  }
  std::pair<double,double>* to_ret = new std::pair<double,double>(theta, psi-e_plane_angle);
  return to_ret;
}



//______________________________________________________________________________
double Gammaavectormeson::getWidth()
{
	return _width;
}


//______________________________________________________________________________
double Gammaavectormeson::getMass()
{
	return _mass;
}


//______________________________________________________________________________
double Gammaavectormeson::getSpin()
{
	return 1.0; //VM spins are the same
}


//______________________________________________________________________________
void Gammaavectormeson::momenta(double W,double Y,double &E,double &px,double &py,double &pz,int &tcheck)
{
	//     This subroutine calculates momentum and energy of vector meson
	//     given W and Y,   without interference.  Subroutine vmpt handles
	//     production with interference
 
	double Egam,Epom,tmin,pt1,pt2,phi1,phi2;
	double px1,py1,px2,py2;
	double pt,xt,xtest,ytest;
	double t2;

  
	//Find Egam,Epom in CM frame
        if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
          // This is pA
          if( _ProductionMode == 2 || _ProductionMode ==3 ){
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp(-Y);
          }else{
    	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }  
        } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
          // This is Ap
          if( _ProductionMode == 2 || _ProductionMode == 3 ){
  	    Egam = 0.5*W*exp(-Y);
  	    Epom = 0.5*W*exp(Y);
          }else{
    	    Egam = 0.5*W*exp(Y);
  	    Epom = 0.5*W*exp(-Y);
          }
	} else {
          // This is pp or AA 
          if( _TargetBeam == 1 ){
            Egam = 0.5*W*exp(-Y);
	    Epom = 0.5*W*exp(Y);
	  }
          else {
            Egam = 0.5*W*exp(Y);
	    Epom = 0.5*W*exp(-Y);
	  }
	}

	//        } else if( _ProductionMode == 2 || _ProductionMode==3){
	//	  Egam = 0.5*W*exp(-Y);
	//	  Epom = 0.5*W*exp(Y);
	//        } else { 
	//          Egam = 0.5*W*exp(Y);
	//	  Epom = 0.5*W*exp(-Y);
	//	 }

        pt1 = pTgamma(Egam);  
	phi1 = 2.*starlightConstants::pi*_randy.Rndom();

	if( (_bbs.beam1().A()==1 && _bbs.beam2().A()==1) || 
            (_ProductionMode == 4) ) {
	    if( (_VMpidtest == starlightConstants::RHO) || (_VMpidtest == starlightConstants::RHOZEUS) || (_VMpidtest == starlightConstants::OMEGA)){
	      // Use dipole form factor for light VM
	      L555vm:
	      xtest = 2.0*_randy.Rndom();
              double ttest = xtest*xtest; 
              ytest = _randy.Rndom();
              double t0 = 1./2.23; 
              double yprob = xtest*_bbs.beam1().dipoleFormFactor(ttest,t0)*_bbs.beam1().dipoleFormFactor(ttest,t0); 
              if( ytest > yprob ) goto L555vm; 
              t2 = ttest; 
              pt2 = xtest;              
	    }else{
		//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
                double bslope_tdist = _VMbslope; 
		double Wgammap = 0.0; 
                switch(_bslopeDef){
		  case 0:
		    //This is the default, as before
		    bslope_tdist = _VMbslope;
		    break;
		  case 1:
		    //User defined value of bslope. BSLOPE_VALUE default is 4.0 if not set. 
                    bslope_tdist = _bslopeVal;
		    if( N0 <= 1 )cout<<" ATTENTION: Using user defined value of bslope = "<<_bslopeVal<<endl;
                    break; 
		  case 2:
                    //This is Wgammap dependence of b from H1 (Eur. Phys. J. C 46 (2006) 585)
		    Wgammap = sqrt(4.*Egam*_pEnergy); 
		    bslope_tdist = 4.63 + 4.*0.164*log(Wgammap/90.0);
		    if( N0 <= 1 )cout<<" ATTENTION: Using energy dependent value of bslope!"<<endl; 
		    break;
		  default:
		    cout<<" Undefined setting for BSLOPE_DEFINITION "<<endl;
		}

	        xtest = _randy.Rndom(); 
		// t2 = (-1./_VMbslope)*log(xtest);
		t2 = (-1./bslope_tdist)*log(xtest);
		pt2 = sqrt(1.*t2);
	    }
	} else {
	    // >> Check tmin
	    tmin = ((Epom/_VMgamma_em)*(Epom/_VMgamma_em));
	
	    if(tmin > 0.5){
		cout<<" WARNING: tmin= "<<tmin<<endl;
                cout<< " Y = "<<Y<<" W = "<<W<<" Epom = "<<Epom<<" gamma = "<<_VMgamma_em<<endl; 
		cout<<" Will pick a new W,Y "<<endl;
		tcheck = 1;
		return;
	    }
 L203vm:
	    xt = _randy.Rndom(); 
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
              }else{
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
              }   
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
              }else{
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
              }  
            } else if (_TargetBeam==1) {
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
            } else {
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();  
            }

	    xtest = _randy.Rndom();
	    t2 = tmin + pt2*pt2;

	    double comp=0.0; 
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
              }else{
                comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
              }   
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
              }else{
                comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
              }  
            } else if (_TargetBeam==1) {
              comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
            } else {
              comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2; 
            }	      
            if( xtest > comp ) goto L203vm;
       		
	}//else end from pp

	phi2 = 2.*starlightConstants::pi*_randy.Rndom();

	px1 = pt1*cos(phi1);
	py1 = pt1*sin(phi1);
	px2 = pt2*cos(phi2);
	py2 = pt2*sin(phi2);
        
	// Compute vector sum Pt = Pt1 + Pt2 to find pt for the vector meson
	px = px1 + px2;
	py = py1 + py2;
	pt = sqrt( px*px + py*py );
       
	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);

}


//______________________________________________________________________________
void Gammaavectormeson::momenta(double W,double Egam,double Q2, double gamma_pz, double gamma_pt,
				double &Y,double &E,double &px,double &py,double &pz,
				double &e_phi,int &tcheck)
{
	//     This subroutine calculates momentum and energy of vector meson
	//     given W and Y,   without interference.  Subroutine vmpt handles
	//     production with interference
 
	double Epom,tmin,pt2,phi1,phi2;
	double px1,py1,px2,py2;
	double pt,xt,xtest,ytest;
	double t2;
      
	// Calculations done in CMS frame (need to recalculate electron energy)
	phi1 = 2.*starlightConstants::pi*_randy.Rndom();
	e_phi = starlightConstants::pi+phi1;
	// Pomeron energy now included photon virtuality and finite transverse momenta
	Epom = 0.5*(W*W+Q2)/(Egam + std::sqrt(Egam*Egam+Q2));
	//Epom = 0.5*(W*W+Q2)/(Egam + gamma_pz);
	while( e_phi > 2.*starlightConstants::pi ) e_phi-= 2.*starlightConstants::pi;
	//
	if( (_bbs.beam1().A()==1 || _bbs.beam2().A()==1) || 
            (_ProductionMode == 4) ) {
	    if( (_VMpidtest == starlightConstants::RHO) || (_VMpidtest == starlightConstants::RHOZEUS) || (_VMpidtest == starlightConstants::OMEGA)){
	      // Use dipole form factor for light VM
	    L613vm:
	      xtest = 2.0*_randy.Rndom();
              double ttest = xtest*xtest; 
              ytest = _randy.Rndom();
              double t0 = 1./2.23; 
              double yprob = xtest*_bbs.beam1().dipoleFormFactor(ttest,t0)*_bbs.beam1().dipoleFormFactor(ttest,t0); 
              if( ytest > yprob ) goto L613vm; 
              t2 = ttest; 
              pt2 = xtest;              
	    }else{
		//Use dsig/dt= exp(-_VMbslope*t) for heavy VM
                double bslope_tdist = _VMbslope; 
		double Wgammap = 0.0; 
                switch(_bslopeDef){
		  case 0:
		    //This is the default, as before
		    bslope_tdist = _VMbslope;
		    break;
		  case 1:
		    //User defined value of bslope. BSLOPE_VALUE default is 4.0 if not set. 
                    bslope_tdist = _bslopeVal;
		    if( N0 <= 1 )cout<<" ATTENTION: Using user defined value of bslope = "<<_bslopeVal<<endl;
                    break; 
		  case 2:
                    //This is Wgammap dependence of b from H1 (Eur. Phys. J. C 46 (2006) 585)
		    Wgammap = sqrt(4.*Egam*_pEnergy); 
		    bslope_tdist = 4.63 + 4.*0.164*log(Wgammap/90.0);
		    if( N0 <= 1 )cout<<" ATTENTION: Using energy dependent value of bslope!"<<endl; 
		    break;
		  default:
		    cout<<" Undefined setting for BSLOPE_DEFINITION "<<endl;
		}

	        xtest = _randy.Rndom(); 
		// t2 = (-1./_VMbslope)*log(xtest);
		t2 = (-1./bslope_tdist)*log(xtest);
		pt2 = sqrt(1.*t2);
	    }
	} else {
	    // >> Check tmin
	    tmin = ((Epom/_VMgamma_em)*(Epom/_VMgamma_em));

	    if(tmin > 0.5){
		cout<<" WARNING: tmin= "<<tmin<<endl;
                cout<< " Y = "<<Y<<" W = "<<W<<" Epom = "<<Epom<<" gamma = "<<_VMgamma_em<<endl; 
		cout<<" Will pick a new W,Y "<<endl;
		tcheck = 1;
		return;
	    }
 L663vm:
	    xt = _randy.Rndom(); 
            if( _bbs.beam1().A()==0 && _bbs.beam2().A() != 1){ 
	      pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam2().nuclearRadius();
	    } 
	    else if( _bbs.beam2().A()==0 && _bbs.beam1().A() != 1 ){
                pt2 = 8.*xt*starlightConstants::hbarc/_bbs.beam1().nuclearRadius();  
              }
	    else{
	      std::cout<<"Can't find the electron for eX"<<std::endl;
	    }  

	    xtest = _randy.Rndom();
	    t2 = tmin + pt2*pt2;

	    double comp=0.0; 
            if( _bbs.beam1().A()==0 && _bbs.beam2().A() != 1){ 
	      comp = _bbs.beam2().formFactor(t2)*_bbs.beam2().formFactor(t2)*pt2;
	    }
	    else if( _bbs.beam2().A()==0 && _bbs.beam1().A() != 1 ){
	      comp = _bbs.beam1().formFactor(t2)*_bbs.beam1().formFactor(t2)*pt2;
            }       
	    else 
	      std::cout<<"Can't find the electron for eX"<<std::endl;
            if( xtest > comp ) goto L663vm;
       		
	}//else end from pp
	phi2 = 2.*starlightConstants::pi*_randy.Rndom();

	px1 = gamma_pt*cos(phi1);
	py1 = gamma_pt*sin(phi1);
	px2 = pt2*cos(phi2);
	py2 = pt2*sin(phi2);
        
	// Compute vector sum Pt = Pt1 + Pt2 to find pt for the vector meson
	px = px1 + px2;
	py = py1 + py2;
	pt = sqrt( px*px + py*py );
       

	//cout<<" \t 0.5*std::log( (Egam+gamma_pz) / Epom ) = "<<Y<<endl;

	E = Egam + Epom;
	pz = gamma_pz - Epom;
	// Testing different methods to extract Y
	//Y = 0.5*std::log( (Egam+gamma_pz) / Epom ); 
	Y = 0.5*std::log( (E+fabs(pz))/(E-fabs(pz)) );
	  
}


//______________________________________________________________________________
double Gammaavectormeson::pTgamma(double E)
{
    // returns on random draw from pp(E) distribution
    double ereds =0.,Cm=0.,Coef=0.,x=0.,pp=0.,test=0.,u=0.;
    double singleformfactorCm=0.,singleformfactorpp1=0.,singleformfactorpp2=0.;
    int satisfy =0;
        
    ereds = (E/_VMgamma_em)*(E/_VMgamma_em);
    //sqrt(3)*E/gamma_em is p_t where the distribution is a maximum
    Cm = sqrt(3.)*E/_VMgamma_em;
    // If E is very small, the drawing of a pT below is extremely slow. 
    // ==> Set pT = sqrt(3.)*E/_VMgamma_em for very small E. 
    // Should have no observable consequences (JN, SRK 11-Sep-2014)
    if( E < 0.0005 )return Cm; 
 
    //the amplitude of the p_t spectrum at the maximum

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3 ){
         singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
         singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
      }else{
         singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
      }  
    } else if (_TargetBeam == 1) {
      singleformfactorCm=_bbs.beam2().formFactor(Cm*Cm+ereds);
    } else {
      singleformfactorCm=_bbs.beam1().formFactor(Cm*Cm+ereds);
    }

    Coef = 3.0*(singleformfactorCm*singleformfactorCm*Cm*Cm*Cm)/((2.*(starlightConstants::pi)*(ereds+Cm*Cm))*(2.*(starlightConstants::pi)*(ereds+Cm*Cm)));
        
    //pick a test value pp, and find the amplitude there
    x = _randy.Rndom();

    if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
      if( _ProductionMode == 2 || _ProductionMode ==3){
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
      }  
    } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
      if( _ProductionMode == 2 || _ProductionMode ==3){
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam2().formFactor(pp*pp+ereds);
      }else{
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
      }  
    } else if (_TargetBeam == 1) {
        pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
    } else {
        pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
        singleformfactorpp1=_bbs.beam1().formFactor(pp*pp+ereds);
    }

    test = (singleformfactorpp1*singleformfactorpp1)*pp*pp*pp/((2.*starlightConstants::pi*(ereds+pp*pp))*(2.*starlightConstants::pi*(ereds+pp*pp)));

    while(satisfy==0){
	u = _randy.Rndom();
	if(u*Coef <= test)
	{
	    satisfy =1;
	}
	else{
	    x =_randy.Rndom();
            if( _bbs.beam1().A()==1 && _bbs.beam2().A() != 1){ 
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);
              }  
            } else if( _bbs.beam2().A()==1 && _bbs.beam1().A() != 1 ){
              if( _ProductionMode == 2 || _ProductionMode ==3){
                pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam2().formFactor(pp*pp+ereds);
              }else{
                pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
                singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
              }  
            } else if (_TargetBeam == 1) {
              pp = x*5.*starlightConstants::hbarc/_bbs.beam2().nuclearRadius(); 
              singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
            } else {
              pp = x*5.*starlightConstants::hbarc/_bbs.beam1().nuclearRadius(); 
              singleformfactorpp2=_bbs.beam1().formFactor(pp*pp+ereds);
            }
	    test = (singleformfactorpp2*singleformfactorpp2)*pp*pp*pp/(2.*starlightConstants::pi*(ereds+pp*pp)*2.*starlightConstants::pi*(ereds+pp*pp));
	}
    }

    return pp;
}


//______________________________________________________________________________
void Gammaavectormeson::vmpt(double W,double Y,double &E,double &px,double &py, double &pz,
                             int&) // tcheck (unused)
{
	//    This function calculates momentum and energy of vector meson
	//    given W and Y, including interference.
	//    It gets the pt distribution from a lookup table.
	double dY=0.,yleft=0.,yfract=0.,xpt=0.,pt1=0.,ptfract=0.,pt=0.,pt2=0.,theta=0.;
	int IY=0,j=0;
  
	dY  = (_VMYmax-_VMYmin)/double(_VMnumy);
  
	//  Y is already fixed; choose a pt
	//  Follow the approach in pickwy
	//  in  _fptarray(IY,pt) IY=1 corresponds to Y=0, IY=numy/2 corresponds to +y
 	// Changed,  now works -y to +y.
	IY=int((Y-_VMYmin)/dY);
	if (IY > (_VMnumy)-1){
        	IY=(_VMnumy)-1;
	}

	yleft=(Y-_VMYmin)-(IY)*dY;

	yfract=yleft*dY;
  
	xpt=_randy.Rndom();
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY][j]) goto L60;
	}
 L60:
  
	//  now do linear interpolation - start with extremes
  	if (j == 0){
		pt1=xpt/_fptarray[IY][j]*_VMdpt/2.;
		goto L80;
	}
	if (j == _VMNPT){
		pt1=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY][j])/(1.-_fptarray[IY][j]);
		goto L80;
	}
  
	//  we're in the middle
  	ptfract=(xpt-_fptarray[IY][j])/(_fptarray[IY][j+1]-_fptarray[IY][j]);
	pt1=(j+1)*_VMdpt+ptfract*_VMdpt;
  
	//  at an extreme in y?
	if (IY == (_VMnumy/2)-1){
		pt=pt1;
		goto L120;
	}
 L80:

	//  interpolate in y repeat for next fractional y bin      
	for(j=0;j<_VMNPT+1;j++){
		if (xpt < _fptarray[IY+1][j]) goto L90;
	}
 L90:
  
	//  now do linear interpolation - start with extremes
	if (j == 0){
		pt2=xpt/_fptarray[IY+1][j]*_VMdpt/2.;
		goto L100;
	}
	if (j == _VMNPT){
		pt2=(_VMptmax-_VMdpt/2.) + _VMdpt/2.*(xpt-_fptarray[IY+1][j])/(1.-_fptarray[IY+1][j]);
		goto L100;
	}
  
	//  we're in the middle
	ptfract=(xpt-_fptarray[IY+1][j])/(_fptarray[IY+1][j+1]-_fptarray[IY+1][j]);
	pt2=(j+1)*_VMdpt+ptfract*_VMdpt;
 L100:

	//  now interpolate in y  
	pt=yfract*pt2+(1-yfract)*pt1;
 L120:

	//  we have a pt 
	theta=2.*starlightConstants::pi*_randy.Rndom();
	px=pt*cos(theta);
	py=pt*sin(theta);

	E  = sqrt(W*W+pt*pt)*cosh(Y);
	pz = sqrt(W*W+pt*pt)*sinh(Y);
	//      randomly choose to make pz negative 50% of the time
	if(_randy.Rndom()>=0.5) pz = -pz;
}


//______________________________________________________________________________
starlightConstants::event Gammaavectormeson::produceEvent(int&)
{
	//Note used; return default event
	return starlightConstants::event();
}


//______________________________________________________________________________
upcEvent Gammaavectormeson::produceEvent()
{
	// The new event type
	upcEvent event;

	int iFbadevent=0;
	int tcheck=0;
	starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
        starlightConstants::particleTypeEnum vmpid = starlightConstants::UNKNOWN; 

	if (_VMpidtest == starlightConstants::FOURPRONG) {
		double        comenergy = 0;
		double        mom[3]    = {0, 0, 0};
		double        E         = 0;
		lorentzVector decayVecs[4];
		do {
			double rapidity = 0;
			pickwy(comenergy, rapidity);
			if (_VMinterferencemode == 0)
				momenta(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
			else if (_VMinterferencemode==1)
				vmpt(comenergy, rapidity, E, mom[0], mom[1], mom[2], tcheck);
		} while (!fourBodyDecay(ipid, E, comenergy, mom, decayVecs, iFbadevent));
		if ((iFbadevent == 0) and (tcheck == 0))
			for (unsigned int i = 0; i < 4; ++i) {
				starlightParticle daughter(decayVecs[i].GetPx(),
				                           decayVecs[i].GetPy(),
				                           decayVecs[i].GetPz(),
				                           starlightConstants::UNKNOWN,  // energy 
				                           starlightConstants::UNKNOWN,  // _mass
				                           ipid,
				                           (i < 2) ? -1 : +1);
				event.addParticle(daughter);
			}
	} else {
		double comenergy = 0.;
		double rapidity = 0.;
		double E = 0.;
		double momx=0.,momy=0.,momz=0.;

		double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
		bool accepted = false;
		do{
			pickwy(comenergy,rapidity);

			if (_VMinterferencemode==0){
				momenta(comenergy,rapidity,E,momx,momy,momz,tcheck);
			
			} else if (_VMinterferencemode==1){
				vmpt(comenergy,rapidity,E,momx,momy,momz,tcheck);
			}
	   
			_nmbAttempts++;

                        vmpid = ipid; 
			twoBodyDecay(ipid,comenergy,momx,momy,momz,px1,py1,pz1,px2,py2,pz2,iFbadevent);
			double pt1chk = sqrt(px1*px1+py1*py1);
			double pt2chk = sqrt(px2*px2+py2*py2);
			double eta1 = pseudoRapidity(px1, py1, pz1);
			double eta2 = pseudoRapidity(px2, py2, pz2);
                        

			if(_ptCutEnabled && !_etaCutEnabled){
				if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(!_ptCutEnabled && _etaCutEnabled){
				if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
					accepted = true;
					_nmbAccepted++;
				}
			}
			else if(_ptCutEnabled && _etaCutEnabled){
				if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
					if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
						accepted = true;
						_nmbAccepted++;
					}
				}
			}
			else if(!_ptCutEnabled && !_etaCutEnabled)
				_nmbAccepted++;
		}while((_ptCutEnabled || _etaCutEnabled) && !accepted);
		if (iFbadevent==0&&tcheck==0) {
			int q1=0,q2=0;
                        int ipid1,ipid2=0;

			double xtest = _randy.Rndom(); 
			if (xtest<0.5)
				{
					q1=1;
					q2=-1;
				}
			else {
				q1=-1;
				q2=1;
			}

                        if ( ipid == 11 || ipid == 13 ){
                          ipid1 = -q1*ipid;
                          ipid2 = -q2*ipid;
                        } else {
                          ipid1 = q1*ipid;
                          ipid2 = q2*ipid;
                        }

			double md = getDaughterMass(vmpid); 
                        double Ed1 = sqrt(md*md+px1*px1+py1*py1+pz1*pz1); 
			starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, ipid1, q1);
			event.addParticle(particle1);

                        double Ed2 = sqrt(md*md+px2*px2+py2*py2+pz2*pz2); 
			starlightParticle particle2(px2, py2, pz2, Ed2, starlightConstants::UNKNOWN, ipid2, q2);
			event.addParticle(particle2);


		}
	}

	return event;

}
double Gammaavectormeson::pseudoRapidity(double px, double py, double pz)
{
	double pT = sqrt(px*px + py*py);
	double p = sqrt(pz*pz + pT*pT);
	double eta = -99.9; if((p-pz) != 0){eta = 0.5*log((p+pz)/(p-pz));}
	return eta;
}

//______________________________________________________________________________
Gammaanarrowvm::Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaanarrowvm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaanarrowvm()"<<endl;
	narrowResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaanarrowvm::~Gammaanarrowvm()
{ }


//______________________________________________________________________________
Gammaaincoherentvm::Gammaaincoherentvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
        cout<<"Reading in luminosity tables. Gammaainkoherentvm()"<<endl;
        read();
        cout<<"Creating and calculating crosssection. Gammaaincoherentvm()"<<endl;
        incoherentVMCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
        _VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
Gammaaincoherentvm::~Gammaaincoherentvm()
{ }


//______________________________________________________________________________
Gammaawidevm::Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. Gammaawidevm()"<<endl;
	read();
	cout<<"Creating and calculating crosssection. Gammaawidevm()"<<endl;
	wideResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
Gammaawidevm::~Gammaawidevm()
{ }


//______________________________________________________________________________
e_Gammaanarrowvm::e_Gammaanarrowvm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. e_Gammaanarrowvm()"<<endl;
	e_read();
	cout<<"Creating and calculating crosssection. e_Gammaanarrowvm()"<<endl;
	e_narrowResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter(); 
}


//______________________________________________________________________________
e_Gammaanarrowvm::~e_Gammaanarrowvm()
{ }


//______________________________________________________________________________
e_Gammaawidevm::e_Gammaawidevm(const inputParameters& input, beamBeamSystem& bbsystem):Gammaavectormeson(input, bbsystem)
{
	cout<<"Reading in luminosity tables. e_Gammaawidevm()"<<endl;
	e_read();
	cout<<"Creating and calculating crosssection. e_Gammaawidevm()"<<endl;
	e_wideResonanceCrossSection sigma(input, bbsystem);
	sigma.crossSectionCalculation(_bwnormsave);
	setTotalChannelCrossSection(sigma.getPhotonNucleusSigma());
	_VMbslope=sigma.slopeParameter();
}


//______________________________________________________________________________
e_Gammaawidevm::~e_Gammaawidevm()
{ }


//______________________________________________________________________________
void Gammaavectormeson::pickwEgamq2(double &W, double &cmsEgamma, double &targetEgamma, 
				 double &Q2, double &gamma_pz, double &gamma_pt,//photon in target frame
				 double &E_prime, double &theta_e //electron
				 )
{
        double dW, dEgamma;
	double xw,xEgamma, xQ2, xtest, q2test, btest;
	int  IW,IGamma, IQ2;
	// ---------
	//	int egamma_draws = 0, cms_egamma_draws =0, q2_draws =0 ;
	// ---------
	dW = (_VMWmax-_VMWmin)/double(_VMnumw);
	//
	//std::chrono::steady_clock::time_point begin_evt = std::chrono::steady_clock::now();
	bool pick_state = false;
	dEgamma = std::log(_targetMaxPhotonEnergy/_targetMinPhotonEnergy);
	double dEgamma_cms = std::log(_cmsMaxPhotonEnergy/_cmsMinPhotonEnergy);
	while( pick_state == false ){
	  //cout<<"Lomnitz start pick"<<endl;
	  xw = _randy.Rndom();
	  W = _VMWmin + xw*(_VMWmax-_VMWmin);
	  if (W < 2 * starlightConstants::pionChargedMass)
	    continue;
	  IW = int((W-_VMWmin)/dW);
	  xEgamma = _randy.Rndom();
	  //
	  targetEgamma = std::exp(std::log(_targetMinPhotonEnergy) + xEgamma*(dEgamma));
	  IGamma = int(_VMnumega*xEgamma);
	  // Holds Q2 and integrated Q2 dependence. Array is saved in target frame
	  std::pair< double, std::vector<double> > this_energy = _g_EQ2array->operator[](IGamma);
	  double intgrated_q2 = this_energy.first;

	  xtest = _randy.Rndom();
	  if( xtest > intgrated_q2 ){
	    //egamma_draws+=1;
	    continue;
	  }
	  N0++; 
	  btest = _randy.Rndom();
	  
	  std::vector<double> photon_flux = this_energy.second;
	  //std::vector<double> photon_flux = _g_EQ2array->operator[](Egamma_tag);
	  double VMQ2min = photon_flux[0];
	  double VMQ2max = photon_flux[1];
	  //
	  double ratio = std::log(VMQ2max/VMQ2min);
 	  double ln_min = std::log(VMQ2min);
	  
	  xQ2 = _randy.Rndom();
	  Q2 = std::exp(ln_min+xQ2*ratio);
	  IQ2 = int(100*xQ2);	
	  // Load from look-up table. Use linear interpolation to evaluate at Q2
	  double y_1 = photon_flux[IQ2+2];
	  double y_2 = photon_flux[IQ2+3];
	  double x_1 = std::exp(ln_min+IQ2*ratio/100);
	  double x_2 = std::exp(ln_min+(1+IQ2)*ratio/100);
	  double m = (y_2 - y_1)/(x_2 - x_1);
	  double c = y_1-m*x_1;
	  double y = m*Q2+c;
	  q2test = _randy.Rndom();
	  if( y < q2test ){
	    //q2_draws++;
	    continue;
	  }
	  // -- Generate electron and photon in Target frame
	  // -- Then boost to CMS frame and check full csg
	  E_prime = _eEnergy - targetEgamma;
	  double cos_theta_e = 1. - Q2/(2.*_eEnergy*E_prime);
	  theta_e = acos(cos_theta_e);
	  double beam_y = acosh(_beamLorentzGamma);	
	  gamma_pt = E_prime*sin(theta_e);
	  double temp_pz = sqrt( targetEgamma*targetEgamma + Q2 - gamma_pt*gamma_pt);
	  // Now boost to Target frame
	  gamma_pz = temp_pz*cosh(beam_y) - targetEgamma*sinh(beam_y); 
	  cmsEgamma = targetEgamma*cosh(beam_y) - temp_pz*sinh(beam_y);
	  // Simple checkl, should not be needed but used for safety
	  if( cmsEgamma < _cmsMinPhotonEnergy ){
	      continue;
	  }
	  // -- Now in CMS frame, check csga (inherited from STARlight)
	  int iGamma_cms = int( log(cmsEgamma/_cmsMinPhotonEnergy)*_VMnumega/dEgamma_cms );
	  if( iGamma_cms < 0 || iGamma_cms > _VMnumega ){
	    //cms_egamma_draws+=1;
	    continue;
	  }
	  xtest = _randy.Rndom();
	  if( _f_WYarray[IW][iGamma_cms] < xtest ){
	    //cms_egamma_draws+=1;
	    continue;
	  }
	  pick_state = true;
	}
	return;
}


//______________________________________________________________________________
eXEvent Gammaavectormeson::e_produceEvent()
{
	// The new event type
	eXEvent event;
	
	int iFbadevent=0;
	int tcheck=0;
	starlightConstants::particleTypeEnum ipid = starlightConstants::UNKNOWN;
        starlightConstants::particleTypeEnum vmpid = starlightConstants::UNKNOWN; 
	// at present 4 prong decay is not implemented
	double comenergy = 0.;
	double rapidity = 0.;
	double Q2 = 0;
	double E = 0.;
	double momx=0.,momy=0.,momz=0.;
	double targetEgamma = 0, cmsEgamma = 0 ;
	double gamma_pz = 0 , gamma_pt = 0, e_theta = 0;
	double px2=0.,px1=0.,py2=0.,py1=0.,pz2=0.,pz1=0.;
	double e_E=0., e_phi=0;
	bool accepted = false;
	do{
	  pickwEgamq2(comenergy,cmsEgamma, targetEgamma, 
		   Q2, gamma_pz, gamma_pt, //photon infor in CMS frame
		   e_E, e_theta);	 //electron info in target frame  
	  //
	  momenta(comenergy,cmsEgamma, Q2, gamma_pz, gamma_pt, //input
		  rapidity, E, momx, momy, momz, //VM
		  e_phi,tcheck); //
	  //
	  // inelasticity: used for angular distributions
	  double col_y = 1. - (e_E/_eEnergy)*std::pow(std::cos(e_theta/2.),2.);
	  double col_polarization = (1 - col_y)/(1-col_y+col_y*col_y/2.);
	  _nmbAttempts++;
	  
	  vmpid = ipid; 
	  // Two body dedcay in eSTARlight includes the angular corrections due to finite virtuality
	  twoBodyDecay(ipid,comenergy,momx,momy,momz,e_phi,col_polarization,
		       px1,py1,pz1,px2,py2,pz2,iFbadevent);
	  double pt1chk = sqrt(px1*px1+py1*py1);
	  double pt2chk = sqrt(px2*px2+py2*py2);
	  double eta1 = pseudoRapidity(px1, py1, pz1);
	  double eta2 = pseudoRapidity(px2, py2, pz2);
                        

	  if(_ptCutEnabled && !_etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(!_ptCutEnabled && _etaCutEnabled){
	    if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
	      accepted = true;
	      _nmbAccepted++;
	    }
	  }
	  else if(_ptCutEnabled && _etaCutEnabled){
	    if(pt1chk > _ptCutMin && pt1chk < _ptCutMax &&  pt2chk > _ptCutMin && pt2chk < _ptCutMax){
	      if(eta1 > _etaCutMin && eta1 < _etaCutMax && eta2 > _etaCutMin && eta2 < _etaCutMax){
		accepted = true;
		_nmbAccepted++;
	      }
	    }
	  }
	  else if(!_ptCutEnabled && !_etaCutEnabled)
	    _nmbAccepted++;
	}while((_ptCutEnabled || _etaCutEnabled) && !accepted);
	if (iFbadevent==0&&tcheck==0) {
	  int q1=0,q2=0;
	  int ipid1,ipid2=0;
	  
	  double xtest = _randy.Rndom(); 
	  if (xtest<0.5)
	    {
	      q1=1;
	      q2=-1;
	    }
	  else {
	    q1=-1;
	    q2=1;
	  }
	  
	  if ( ipid == 11 || ipid == 13 ){
	    ipid1 = -q1*ipid;
	    ipid2 = -q2*ipid;
	  } else {
	    ipid1 = q1*ipid;
	    ipid2 = q2*ipid;
	  }

	  // - Outgoing electron - target frame - update later
	  double e_px = e_E*sin(e_theta)*cos(e_phi);
	  double e_py = e_E*sin(e_theta)*sin(e_phi);
	  double e_pz = e_E*cos(e_theta);
	  lorentzVector electron(e_px, e_py, e_pz, e_E);
	  event.addSourceElectron(electron);
	  // - Generated photon - CMS frame
	  double gamma_x = gamma_pt*cos(e_phi+starlightConstants::pi);
	  double gamma_y = gamma_pt*sin(e_phi+starlightConstants::pi);
	  lorentzVector gamma(gamma_x,gamma_y,gamma_pz,cmsEgamma);
	  event.addGamma(gamma, cmsEgamma, Q2);
	  //
	  double md = getDaughterMass(vmpid); 
	  double Ed1 = sqrt(md*md+px1*px1+py1*py1+pz1*pz1); 
	  starlightParticle particle1(px1, py1, pz1, Ed1, starlightConstants::UNKNOWN, ipid1, q1);
	  event.addParticle(particle1);
	  
	  double Ed2 = sqrt(md*md+px2*px2+py2*py2+pz2*pz2); 
	  starlightParticle particle2(px2, py2, pz2, Ed2, starlightConstants::UNKNOWN, ipid2, q2);
	  event.addParticle(particle2);
	  
	  
	}
	return event;

}
string Gammaavectormeson::gammaTableParse(int ii, int jj)
{
  ostringstream tag1, tag2;
  tag1<<ii;
  tag2<<jj;
  string to_ret = tag1.str()+","+tag2.str();
  return to_ret;
}
