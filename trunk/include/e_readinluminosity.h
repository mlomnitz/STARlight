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
// $Rev:: 213                         $: revision of last commit
// $Author:: mlomnitz                   $: author of last commit
// $Date:: 2017-03-12 22:08:02 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef READIN_E_LUMINOSITY_H
#define READIN_E_LUMINOSITY_H


#include "inputParameters.h"
#include "starlightlimits.h"


class e_readLuminosity
{
 public:
  e_readLuminosity(const inputParameters& input);
  ~e_readLuminosity();
  
  void e_read();
  double *_Warray;
  double *_Yarray;
  double **_f_WYarray;
  double **_g_Earray;
  
  double _f_WYmax;
  double _g_Emax;

  double _bwnormsave;

 protected:
  const int _ReadInputnumy;
  const int _ReadInputnumw;
  const std::string _baseFileName;
};


#endif  // READINLUMINOSITY_H
