/*//////////////////////////////////////////////////////////////////                                                   
////     The SKIRT project -- advanced radiative transfer       ////
////       © Junkai Zhang, University of Bath                  ////
///////////////////////////////////////////////////////////////// */

#include "clumpSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> clumpSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::age(), SnapshotParameter::metallicity(), SnapshotParameter::custom("DTS", string(), string()), SnapshotParameter::initialMass()}; 	    
}

////////////////////////////////////////////////////////////////////

void clumpSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, filename(), "lambda(m),Age(Myr),metal(1),DTS(1),Mstar(Msun)", "Llambda(W/m)", false, false);
}

////////////////////////////////////////////////////////////////////

Range clumpSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double clumpSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double age = parameters[0] / Constants::year() /1000000;
    double metal = parameters[1];
    double DTS = parameters[2];
    double Mstar = parameters[3] / Constants::Msun(); 

    return Mstar * _table(wavelength, age, metal, DTS, Mstar);
}

////////////////////////////////////////////////////////////////////

double clumpSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                    const Array& parameters) const
{

    double age = parameters[0] / Constants::year() /1000000;
    double metal = parameters[1]; 
    double DTS = parameters[2];             
    double Mstar = parameters[3] / Constants::Msun(); 

    return Mstar * _table.cdf(lambdav, pv, Pv, wavelengthRange, age, metal, DTS, Mstar);
}

//////////////////////////////////////////////////////////////////// 



