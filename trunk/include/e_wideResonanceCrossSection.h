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
// $Rev:: 211                         $: revision of last commit
// $Author:: mlomnitz                   $: author of last commit
// $Date:: 2017-03-14 03:05:09 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef e_WIDERESONANCECROSSSECTION_H
#define e_WIDERESONANCECROSSSECTION_H


#include "photonNucleusCrossSection.h"
#include "inputParameters.h"

class e_wideResonanceCrossSection : public photonNucleusCrossSection {

public:

	e_wideResonanceCrossSection(const inputParameters& input, const beamBeamSystem& bbsystem);
	~e_wideResonanceCrossSection();

	void crossSectionCalculation(const double bwnormsave);

private:

	double _Ep;  // Proton Energy
	double _wideWmax;
	double _wideWmin;
	double _wideYmax;
	double _wideYmin;		

};


#endif  // WIDERESONANCECROSSSECTION_H
