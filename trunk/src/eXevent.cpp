///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2017
//
//    This file is part of estarlight.
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
// $Rev:: 263                         $: revision of last commit
// $Author:: mlomnitz                  $: author of last commit
// $Date:: 02/28/2017 #$: date of last commit
//
// Description:
//
// Container for eX eventrs
//
//
///////////////////////////////////////////////////////////////////////////


#include "eXevent.h"


eXEvent::eXEvent() :
        _particles(0)
        ,_vertices(0)
{ }

eXEvent::eXEvent(starlightConstants::event &ev) :
        _particles(0)
        ,_vertices(0)
{
  for(int i = 0; i < ev._numberOfTracks; i++)
    {
      starlightParticle p(
			  ev.px[i], 
			  ev.py[i], 
			  ev.pz[i], 
			  starlightConstants::UNKNOWN, 
			  starlightConstants::UNKNOWN, 
			  ev._fsParticle[i],
			  ev._charge[i]
			  );
      addParticle(p);
    }
}

eXEvent::~eXEvent()
{ }


eXEvent& eXEvent::operator=(const eXEvent& rhs)
{

  if(this != &rhs)
  {
    this->_particles = rhs._particles;
    this->_vertices = rhs._vertices;
    this->_gammaEnergies = rhs._gammaEnergies;
  }
  return *this;
}

eXEvent& eXEvent::operator+(const eXEvent& ev)
{
  for(unsigned int n = 0; n < ev._particles.size(); n++)
  {
    this->_particles.push_back(ev._particles.at(n));
  }
  for(unsigned int n = 0; n < ev._vertices.size(); n++)
  {
    this->_vertices.push_back(ev._vertices.at(n));
  }
 for(unsigned int n = 0; n < ev._gammaEnergies.size(); n++)
  {
    this->_gammaEnergies.push_back(ev._gammaEnergies.at(n));
  }
  return *this;
}

void eXEvent::boost(double rapidity)
{
    vector3 boostVector(0, 0, tanh(rapidity));
    std::vector<starlightParticle>::iterator part = _particles.begin();
      
    for (part = _particles.begin(); part != _particles.end(); part++)
    {
      (*part).Boost(boostVector);
    }
}
