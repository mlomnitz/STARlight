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
// $Rev:: 265                         $: revision of last commit
// $Author:: mlomnitz                   $: author of last commit
// $Date:: 2017-03-12 21:37:26 +0100 #$: date of last commit
//
// Description:
//    Added 18->19 for reading in the luminosity table
//    Incoherent factor added to table --Joey
//
//
//
///////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "readin_e_luminosity.h"
#include "starlightconstants.h"
#include "inputParameters.h"

using namespace std;


//______________________________________________________________________________
e_readLuminosity::e_readLuminosity(const inputParameters& inputParametersInstance)
  : _Warray(0), _Yarray(0), _f_WYarray(0), g_Earray,
  , _ReadInputnumy(inputParametersInstance.nmbRapidityBins())
  , _ReadInputnumw(inputParametersInstance.nmbWBins())
  , _baseFileName(inputParametersInstance.baseFileName())
{

}


//______________________________________________________________________________
e_readLuminosity::~e_readLuminosity()
{ 
  if(_Warray) delete [] _Warray;
  if(_Yarray) delete [] _Yarray;
  if(_Farray) delete [] _f_WYarray;
  if(_Farray1) delete [] _g_Earray;
}


//______________________________________________________________________________
void e_readLuminosity::e_read()
{
  
  if(!_Warray) _Warray = new double[_ReadInputnumw];
  if(!_Yarray) _Yarray = new double[_ReadInputnumy];
  if(!_f_WYmax) 
  {
    _f_WYarray = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _f_WYarra[i] = new double[_ReadInputnumy];
    }
  }
  if(!_g_Earray) 
  {
    _g_Earray = new double*[_ReadInputnumw];
    for(int i = 0; i < _ReadInputnumw; i++)
    {
      _g_Earray[i] = new double[_ReadInputnumy];
    }
  }

  double dummy[13]; //number of lines used to read in input parameters saved to lookup table[slight.txt].


  std::string wyFileName;
  wyFileName = _baseFileName +".txt";
  
//  cout << "wyFileName being read in" << wyFileName << endl;

  double fpart =0.;
  double fptsum=0.;
  ifstream wylumfile;

  _f_WYmax=0.0;
  _g_Emax=0.0;

  wylumfile.open(wyFileName.c_str());

  for(int i=0;i < 13;i++){ 
    wylumfile >> dummy[i];
  }
  int A_1 = dummy[1];
  int A_2 = dummy[3];

  for(int i=0;i<_ReadInputnumw;i++){
    wylumfile >> _Warray[i];
  }
  for(int i=0;i<_ReadInputnumy;i++){
    wylumfile >> _Yarray[i];
  }

  if( (A_2 == 0 && A_1 >= 1) || (A_1 ==0 && A_2 >= 1) ){ 
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        wylumfile >> _f_WYarray[i][j];
        if( _f_WYarray[i][j] > _f_WYmax ) _f_WYmax=_f_WYarray[i][j];
	//
	wylumfile >> _g_Earray[i][j];
	if( _g_Earray[i][j] > _g_Emax ) _g_Emax = _g_Earray[i][j];
      }
    }
    //Normalize f_WY array, g does not need to be normalized, it is used for normalization
    for(int i=0;i<_ReadInputnumw;i++){
      for(int j=0;j<_ReadInputnumy;j++){
        _f_WYarray[i][j] = _f_WYarray[i][j]/_f_WYmax;
	//_g_Earray[i][j] = _g_Earray[i][j]/_g_Emax;
      }
    }
  }
  wylumfile >> _bwnormsave;

  wylumfile.close();
  return;
}
