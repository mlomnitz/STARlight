///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016
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
// $Rev:: 0                         $: revision of last commit
// $Author:: mlomnitz                $: author of last commit
// $Date:: 2016-11-23 10:29:00 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef ELECTRON_H
#define ELECTRON_H


#include <cmath>

//This class holds the information for a target electron
class electron
{

 public:
  electron();
  electron(const int productionMode) { _productionMode = productionMode}
  ~electron();
  
  void init();
  int productionMode() const { return _productionMode;} //need to review if this will be necesarry 

  //need to review if this is necessary, will make beam class easier for photon density calculation 
  int Z() const { return 1 } 
  int A() const { return 0 } 
  double nuclearRadius() const { return 0.; }
  double thickness() const { return 0.; }

 private:
  
  int _productionMode;
}

#endif //ELECTRON_H
