//! \file  keAcquisition.h
//! \brief mean kinetic energy budget object declaration

#ifndef _KEACQUISITION_H
#define _KEACQUISITION_H

//! \brief Struct defining the ke fields
struct keFields
{
    // fields
    Vec           Error;           //!< error on the equation
    Vec           lEm;             //!< mechanical energy
    Vec           D;               //!< mechanical energy divergence
    Vec           F;               //!< net energy flux due to turbulent transport
    Vec           Pinf;            //!< mean background pressure
    Vec           Pf;              //!< wind turbine power extraction
    Vec           Ptheta;          //!< kinetic to potential energy conversion
    Vec           Eps;             //!< turbulent dissipation

    // working fields
    Vec           lavgUm;           //!< <u> (vector)
    Vec           lavgPm;           //!< <p> (scalar)
    Vec           lavgUpUp;         //!< <u'u'> (symmetric tensor)
    Vec           lavgUpUpUp;       //!< <u_i' u_i' u'> (vector)
    Vec           lavgUmTauSGS;     //!< <u_i TauSGS_ij> (vector)
    Vec           lavgUpPp;         //!< <u'p'> (vector)

};

#endif
