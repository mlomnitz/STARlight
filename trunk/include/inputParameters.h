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
// $Rev:: 276                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-09-13 19:54:42 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef INPUTPARAMETERS_H
#define INPUTPARAMETERS_H


#include "starlightconstants.h"
#include "inputParser.h"
#include <string>
#include <ostream>
#include <vector>
#include <sstream>

class parameterbase;


class parameterlist
{
public:

    parameterlist() : _parameters(0) {}

    void add(parameterbase* p) {
        _parameters.push_back(p);
    }

    // Returns a string with a key of the current state of the parameter list
    // only
    inline std::string validationKey();
    

private:

    std::vector<parameterbase*> _parameters;

};

// Base class for parameters, needed to keep a list of parameters
class parameterbase
{
public:

    // Add this to parameter list
    parameterbase()
    {
        _parameters.add(this);
    }
    virtual std::string validationkey() = 0;

    template<typename T>
    std::string toString(T v)
    {
        std::stringstream s;
        s << v;
        return s.str();
    }
    inline friend std::ostream& operator<<(std::ostream& os, const parameterbase& par);
 
    // List of all parameters
    static parameterlist _parameters;


   
};
// Need to init the static variable
// parameterlist parameterbase::_parameters;


// The actual parameter class
// validate parameter specifies if the parameter should be a part of the validity check of the current parameters
template<typename T, bool validate>
class parameter : public parameterbase
{
public:

    // Constructor
    parameter(const std::string &name, T value, bool required = true) :parameterbase(),_name(name), _value(value), _validate(validate), _required(required) {}


    parameter &operator=(T v) { _value = v; return *this;}
    T* ptr() const {
        return const_cast<T*>(&_value);
    }
    
    T value() const { return _value; }
    
    std::string name() const { return _name;}
    
    bool required() const { return _required; }
    
    void setValue(T v) { _value = v; }
    
    void setName(std::string name) { _name = name; }
    
    void setRequired(bool r) { _required = r; }
    
    // Validation key for this parameter
    std::string validationkey()
    {
        return (_validate ? _name + ":" + toString(_value) + "-" : std::string(""));
    }

    template<typename S, bool v>
    inline friend std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par);



private:
    std::string _name;

    T _value; // Value
    bool _validate; // true if a change in the parameter invalidates x-sec tables
    bool _required; // true if this is required option.

    parameter();
};

template<typename S, bool v>
inline std::ostream& operator<<(std::ostream& os, const parameter<S,v>& par)
{
    os << par._value;
    return os;
}

std::ostream& operator<<(std::ostream& os, const parameterbase& par)
{
    os << par._parameters.validationKey(); 
    return os;
}
std::string parameterlist::validationKey()
{
    std::stringstream s;
    for(unsigned int i = 0; i < _parameters.size(); ++i)
    {
        s << _parameters[i]->validationkey(); // Will print names and values of validation parameters
    }
    return s.str();
}

class inputParameters {

public:
	inputParameters();
	~inputParameters();

	bool init();
	bool configureFromFile(const std::string &configFileName = "./config/slight.in");

        std::string  baseFileName          () const { return _baseFileName.value();           }
	unsigned int beam1Z                () const { return _beam1Z.value();                 }  ///< returns atomic number of beam particle 1
	unsigned int beam1A                () const { return _beam1A.value();                 }  ///< returns atomic mass number of beam particle 1
	unsigned int beam2Z                () const { return _beam2Z.value();                 }  ///< returns atomic number of beam particle 2
	unsigned int beam2A                () const { return _beam2A.value();                 }  ///< returns atomic mass number of beam particle 2
	double       targetLorentzGamma    () const { return _targetLorentzGamma;             }  ///< returns Lorentz gamma factor of source in target frame
	double       beamLorentzGamma      () const { return _beamLorentzGamma;       	      }  ///< returns Lorentz gamma factor of both beams in beam CMS frame
	double       beam1LorentzGamma     () const { return _beam1LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 1 in collider frame
	double       beam2LorentzGamma     () const { return _beam2LorentzGamma.value();      }  ///< returns Lorentz gamma factor of beam 2 in collider frame
	double       maxW                  () const { return _maxW.value();                   }  ///< returns maximum mass W of produced hadronic system [GeV/c^2]
	double       minW                  () const { return _minW.value();                   }  ///< returns minimum mass W of produced hadronic system [GeV/c^2]
	unsigned int nmbWBins              () const { return _nmbWBins.value();               }  ///< returns number of W bins in lookup table
	double       maxRapidity           () const { return _maxRapidity.value();            }  ///< returns maximum absolute value of rapidity
	unsigned int nmbRapidityBins       () const { return _nmbRapidityBins.value();        }  ///< returns number of rapidity bins in lookup table
	bool         ptCutEnabled          () const { return _ptCutEnabled.value();           }  ///< returns cut in pt
	double       ptCutMin              () const { return _ptCutMin.value();               }  ///< returns minimum pt
	double       ptCutMax              () const { return _ptCutMax.value();               }  ///< returns maximum pt
	bool         etaCutEnabled         () const { return _etaCutEnabled.value();          }  ///< returns cut in eta
	double       etaCutMin             () const { return _etaCutMin.value();              }  ///< returns minimum eta
	double       etaCutMax             () const { return _etaCutMax.value();              }  ///< returns maximum eta
	int          productionMode        () const { return _productionMode.value();         }  ///< returns production mode
	unsigned int nmbEvents             () const { return _nmbEventsTot.value();           }  ///< returns total number of events to generate
	int          prodParticleId        () const { return _prodParticleId.value();         }  ///< returns PDG particle ID of produced particle
	int          randomSeed            () const { return _randomSeed.value();             }  ///< returns seed for random number generator
	int          beamBreakupMode       () const { return _beamBreakupMode.value();        }  ///< returns breakup mode for beam particles
	bool         interferenceEnabled   () const { return _interferenceEnabled.value();    }  ///< returns whether interference is taken into account
	double       interferenceStrength  () const { return _interferenceStrength.value();   }  ///< returns percentage of interference
	double       maxPtInterference     () const { return _maxPtInterference.value();      }  ///< returns maximum p_T for interference calculation [GeV/c]
	int          nmbPtBinsInterference () const { return _nmbPtBinsInterference.value();  }  ///< returns number of p_T bins for interference calculation
	double       ptBinWidthInterference() const { return _ptBinWidthInterference.value(); }  ///< returns width of p_T bins for interference calculation [GeV/c]
	double 	     minGammaEnergy        () const { return _minGammaEnergy.value();         }  ///< returns minimum gamma energy in case of photo nuclear processes [GeV]
	double       maxGammaEnergy        () const { return _maxGammaEnergy.value();         }  ///< returns maximum gamma energy in case of photo nuclear processes [GeV]
	double       minGammaQ2            () const { return _minGammaQ2.value();             }  ///< return minimum gamma virtuality 
	double       maxGammaQ2            () const { return _maxGammaQ2.value();             }  ///< return maximum gamma virtuality 
	unsigned int nmbGammaQ2Bins        () const { return _nmbGammaQ2Bins.value();         }  ///< return number of gamma q2 bins
	std::string  pythiaParams          () const { return _pythiaParams.value();           }  ///< returns parameters to be passed to pythia
	bool         pythiaFullEventRecord () const { return _pythiaFullEventRecord.value();  }  ///< returns if the full pythia event record should be printed
	int	     xsecCalcMethod        () const { return _xsecCalcMethod.value();         }  ///< returns the method used for the x-sec calculation
        double       axionMass             () const { return _axionMass.value();              }  ///< returns axion mass //AXION HACK
	int          bslopeDefinition      () const { return _bslopeDefinition.value();       }  ///< returns the definition of b-slope
	double       bslopeValue           () const { return _bslopeValue.value();            }  ///< returns the value of b-slope
	starlightConstants::particleTypeEnum    prodParticleType     () const { return _particleType;    }  ///< returns type of produced particle
	starlightConstants::decayTypeEnum       prodParticleDecayType() const { return _decayType;       }  ///< returns decay type of produced particle
	starlightConstants::interactionTypeEnum interactionType      () const { return _interactionType; }  ///< returns interaction type
	double protonEnergy                () const { return _protonEnergy.value(); }
	double electronEnergy              () const { return _electronEnergy.value(); }
        double inputBranchingRatio         () const { return _inputBranchingRatio; }
	
        void setBaseFileName          (std::string v )  {  _baseFileName = v;     }
	void setBeam1Z                (unsigned int v)  {  _beam1Z = v;           }  ///< sets atomic number of beam particle 1
	void setBeam1A                (unsigned int v)  {  _beam1A = v;           }  ///< sets atomic mass number of beam particle 1
	void setBeam2Z                (unsigned int v)  {  _beam2Z = v;           }  ///< sets atomic number of beam particle 2
	void setBeam2A                (unsigned int v)  {  _beam2A = v;           }  ///< sets atomic mass number of beam particle 2
	void setTargetLorentzGamma    (double v)  {  _targetLorentzGamma = v;     }  ///< sets Lorentz gamma factor of both beams in beam CMS frame
	void setBeamLorentzGamma      (double v)  {  _beamLorentzGamma = v;       }  ///< sets Lorentz gamma factor of both beams in beam CMS frame
	void setBeam1LorentzGamma     (double v)  {  _beam1LorentzGamma = v;      }  ///< sets Lorentz gamma factor of beam 1 in collider frame
	void setBeam2LorentzGamma     (double v)  {  _beam2LorentzGamma = v;      }  ///< sets Lorentz gamma factor of beam 2 in collider frame
	void setMaxW                  (double v)  {  _maxW = v;                   }  ///< sets maximum mass W of produced hadronic system [GeV/c^2]
	void setMinW                  (double v)  {  _minW = v;                   }  ///< sets minimum mass W of produced hadronic system [GeV/c^2]
	void setNmbWBins              (unsigned int v)  {  _nmbWBins = v;         }  ///< sets number of W bins in lookup table
	void setMaxRapidity           (double v)  {  _maxRapidity = v;            }  ///< sets maximum absolute value of rapidity
	void setNmbRapidityBins       (unsigned int v)  {  _nmbRapidityBins = v;  }  ///< sets number of rapidity bins in lookup table
	void setPtCutEnabled          (bool v)  {  _ptCutEnabled = v;             }  ///< sets cut in pt
	void setPtCutMin              (double v)  {  _ptCutMin = v;               }  ///< sets minimum pt
	void setPtCutMax              (double v)  {  _ptCutMax = v;               }  ///< sets maximum pt
	void setEtaCutEnabled         (bool v)  {  _etaCutEnabled = v;            }  ///< sets cut in eta
	void setEtaCutMin             (double v)  {  _etaCutMin = v;              }  ///< sets minimum eta
	void setEtaCutMax             (double v)  {  _etaCutMax = v;              }  ///< sets maximum eta
	void setProductionMode        (int v)  {  _productionMode = v;            }  ///< sets production mode
	void setNmbEvents             (unsigned int v)  {  _nmbEventsTot = v;     }  ///< sets total number of events to generate
	void setProdParticleId        (int v)  {  _prodParticleId = v;            }  ///< sets PDG particle ID of produced particle
	void setRandomSeed            (int v)  {  _randomSeed = v;                }  ///< sets seed for random number generator
	void setBeamBreakupMode       (int v)  {  _beamBreakupMode = v;           }  ///< sets breakup mode for beam particles
	void setInterferenceEnabled   (bool v)  {  _interferenceEnabled = v;      }  ///< sets whether interference is taken into account
	void setInterferenceStrength  (double v)  {  _interferenceStrength = v;   }  ///< sets percentage of interference
	void setMaxPtInterference     (double v)  {  _maxPtInterference = v;      }  ///< sets maximum p_T for voiderference calculation [GeV/c]
	void setNmbPtBinsInterference (int v)  {  _nmbPtBinsInterference = v;     }  ///< sets number of p_T bins for interference calculation
	void setPtBinWidthInterference(double v)  {  _ptBinWidthInterference = v; }  ///< sets width of p_T bins for voiderference calculation [GeV/c]
	void setMinGammaEnergy        (double v)  {  _minGammaEnergy = v;         }  ///< sets minimum gamma energy in case of photo nuclear processes [GeV]
	void setMaxGammaEnergy        (double v)  {  _maxGammaEnergy = v;         }  ///< sets maximum gamma energy in case of photo nuclear processes [GeV]
	void setPythiaParams          (std::string v)  {  _pythiaParams = v;      }  ///< sets parameters to be passed to pythia
	void setPythiaFullEventRecord (bool v)  {  _pythiaFullEventRecord = v;    }  ///< sets if the full pythia event record should be prvoided
	void setXsecCalcMethod        (int v)  {  _xsecCalcMethod = v;            }  ///< sets the method used for the x-sec calculation
	void setAxionMass        (double v)  {  _axionMass = v;                   }  ///< sets axion mass    //AXION HACK
	void setbslopeDefinition      (int v)  {  _bslopeDefinition = v;          }  ///< sets the definition of b slope
        void setbslopeValue           (double v)  {  _bslopeValue = v;            }  ///< sets the value of b slope

	void setProdParticleType      (starlightConstants::particleTypeEnum v)    { _particleType = v;    }  ///< sets type of produced particle
	void setProdParticleDecayType (starlightConstants::decayTypeEnum v)       { _decayType = v;       }  ///< sets decay type of produced particle
	void setInteractionType       (starlightConstants::interactionTypeEnum v) { _interactionType = v; }  ///< sets interaction type
	 
	void setProtonEnergy        (double v)    { _protonEnergy = v;            }  ///< sets proton energy 
	void setElectronEnergy      (double v)    { _electronEnergy = v;          }  ///< sets electron energy
		
	template<typename T>
	inline bool setParameter(std::string expression);
	
	std::ostream& print(std::ostream& out) const;  ///< prints parameter summary
	std::ostream& write(std::ostream& out) const;  ///< writes parameters back to an ostream
	
	std::string parameterValueKey() const; ///< Generates key for the current parameters

  
private:

    
// To indicate if the crossection table should be re-calculated if parameter changes
#define VALIDITY_CHECK true
#define NO_VALIDITY_CHECK false
	
	std::string _configFileName;  ///< path to configuration file (default = ./config/slight.in)

	// config file parameters
        parameter<std::string,NO_VALIDITY_CHECK>   _baseFileName;
	parameter<unsigned int,VALIDITY_CHECK>     _beam1Z;                  ///< atomic number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam1A;                  ///< atomic mass number of beam particle 1
	parameter<unsigned int,VALIDITY_CHECK>     _beam2Z;                  ///< atomic number of beam particle 2
	parameter<unsigned int,VALIDITY_CHECK>     _beam2A;                  ///< atomic mass number of beam particle 2
	parameter<double, VALIDITY_CHECK>          _beam1LorentzGamma;       ///< Lorentz gamma factor of beam 1 in collider frame
	parameter<double, VALIDITY_CHECK>          _beam2LorentzGamma;       ///< Lorentz gamma factor of beam 2 in collider frame
	parameter<double, VALIDITY_CHECK>          _maxW;                    ///< maximum mass W of produced hadronic system [GeV/c^2]
	parameter<double, VALIDITY_CHECK>          _minW;                    ///< minimum mass W of produced hadronic system; if set to -1 default value is taken [GeV/c^2]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbWBins;                ///< number of W bins in lookup table
	parameter<double, VALIDITY_CHECK>          _maxRapidity;             ///< maximum absolute value of rapidity
	parameter<unsigned int, VALIDITY_CHECK>    _nmbRapidityBins;         ///< number of rapidity bins in lookup table
	parameter<bool, VALIDITY_CHECK>            _ptCutEnabled;            ///< en/disables cut in pt
	parameter<double, VALIDITY_CHECK>          _ptCutMin;                ///< minimum pt, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _ptCutMax;                ///< maximum pt, if cut is enabled
	parameter<bool, VALIDITY_CHECK>            _etaCutEnabled;           ///< en/disables cut in eta
	parameter<double, VALIDITY_CHECK>          _etaCutMin;               ///< minimum eta, if cut is enabled
	parameter<double, VALIDITY_CHECK>          _etaCutMax;               ///< maximum eta, if cut is enabled
	parameter<unsigned int, VALIDITY_CHECK>    _productionMode;          ///< \brief production mode
	                                                                     ///<
	                                                                     ///< 1 = photon-photon fusion,
	                                                                     ///< 2 = narrow vector meson resonance in photon-Pomeron fusion,
	                                                                     ///< 3 = Breit-Wigner vector meson resonance in photon-Pomeron fusion
	parameter<unsigned int, VALIDITY_CHECK>    _nmbEventsTot;            ///< total number of events to generate
	parameter<unsigned int, VALIDITY_CHECK>    _prodParticleId;          ///< PDG particle ID of produced particle
	parameter<unsigned int, VALIDITY_CHECK>    _randomSeed;              ///< seed for random number generator
	                                                                     ///<
	                                                                     ///< 1 = ASCII
	                                                                     ///< 2 = GSTARtext,
	                                                                     ///< 3 = PAW ntuple (not working)
	parameter<unsigned int, VALIDITY_CHECK>    _beamBreakupMode;         ///< \brief breakup mode for beam particles
	                                                                     ///<
	                                                                     ///< 1 = hard sphere nuclei (b > 2R),
	                                                                     ///< 2 = both nuclei break up (XnXn),
	                                                                     ///< 3 = a single neutron from each nucleus (1n1n),
	                                                                     ///< 4 = neither nucleon breaks up (with b > 2R),
	                                                                     ///< 5 = no hadronic break up (similar to option 1, but with the actual hadronic interaction)
	parameter<bool, VALIDITY_CHECK>            _interferenceEnabled;     ///< if VALIDITY_CHECK, interference is taken into account
	parameter<double, VALIDITY_CHECK>          _interferenceStrength;    ///< percentage of interference: from 0 = none to 1 = full
	parameter<double, VALIDITY_CHECK>          _maxPtInterference;       ///< maximum p_T for interference calculation [GeV/c]
	parameter<unsigned int, VALIDITY_CHECK>    _nmbPtBinsInterference;   ///< number of p_T bins for interference calculation
	parameter<double, VALIDITY_CHECK>          _ptBinWidthInterference;  ///< width of p_T bins for interference calculation [GeV/c]
	parameter<double, VALIDITY_CHECK>          _protonEnergy;
	parameter<double, VALIDITY_CHECK>          _electronEnergy;
	parameter<double, VALIDITY_CHECK>          _minGammaEnergy;          ///< minimum gamma energy in case of photo nuclear processes [GeV]
	parameter<double, VALIDITY_CHECK>          _maxGammaEnergy;          ///< maximum gamma energy in case of photo nuclear processes [GeV]
	parameter<double, VALIDITY_CHECK>          _minGammaQ2;              ///< minimum gamma Q2 in case of photo nuclear processes
	parameter<double, VALIDITY_CHECK>          _maxGammaQ2;              ///< maximum gamma Q2 in case of photo nuclear processes	
	parameter<unsigned int, VALIDITY_CHECK>     _nmbGammaQ2Bins;          ///< number of gamma q2 bins
	parameter<std::string,NO_VALIDITY_CHECK>   _pythiaParams;            ///< semi-colon separated parameters to pass to pythia, e.g. "mstj(1)=0;paru(13)=0.1" 
	parameter<bool, NO_VALIDITY_CHECK>         _pythiaFullEventRecord;   ///< if the full pythia event record should be in the output
	parameter<unsigned int, VALIDITY_CHECK>    _xsecCalcMethod;	     ///< Select x-sec calc method. (0 is standard starlight method, 1 must be used for assym. collisions (e.g. p-A), but is slow)	
        parameter<double, VALIDITY_CHECK>          _axionMass;               ///Axion mass//AXION HACK
        parameter<unsigned int, VALIDITY_CHECK>    _bslopeDefinition;        ///< Optional parameter to set different values of slope parameter
        parameter<double, VALIDITY_CHECK>          _bslopeValue;             ///< Value of slope parameter when _bslopeDefinition is set to 1

	starlightConstants::particleTypeEnum       _particleType;
	starlightConstants::decayTypeEnum          _decayType;
	starlightConstants::interactionTypeEnum    _interactionType;

	double                         _targetLorentzGamma;       ///< Lorentz gamma factor of the source in target frame, not an input parameter
	double                         _beamLorentzGamma;         ///< Lorentz gamma factor of the beams in CMS frame, not an input parameter
	double                         _inputBranchingRatio;      ///< Branching ratio defined for each channel
	
	inputParser _ip;
	
};


template<typename T>
inline 
bool inputParameters::setParameter(std::string expression)
{
   
    return _ip.parseString(expression);
   
   
}

inline
std::ostream&
operator <<(std::ostream&          out,
            const inputParameters& par)
{
	return par.print(out);
}
 

#endif  // INPUTPARAMETERS_H
