//! \file  averageAcquisition.h
//! \brief average fields object declaration

#ifndef _AVERAGEACQUISITION_H
#define _AVERAGEACQUISITION_H

//! \brief Struct defining the average fields
struct avgFields
{
    // source terms
    Vec           Coriolis;                   //!< coriolis source
    Vec           Driving;                    //!< driving pressure gradient source
    Vec           xDamping;                   //!< x damping layer (x fringe region source)
    Vec           SideForce;                  //!< side force for fringe region testing

    Vec           windFarmForce;              //!< wind farm body force

    // turbulence criterions
    Vec           Q;                          //!< Q criterion field
    Vec           L2;                         //!< L2 criterion field (not yet implemented)

    // averaging
    Vec           avgU;                       //!< sum of u, v, w over time
    Vec           avgP;                       //!< sum of p over time
    Vec           avgUU;                      //!< sum of resolved R stresses over time
    Vec           avgP2;                      //!< sum of p^2 over time
    Vec           avgNut;                     //!< sun of nut over time
    Vec           avgCs;                      //!< sum of cs over time
    Vec           avgOmega;                   //!< sum of vorticity over time
    Vec           avgOmegaOmega;              //!< sum of resolved vorticity stresses over time
    Vec           avgUdotGradP;               //!< sum of scalar u*dpdx + v*dpdy + w*dpdz over time
    Vec           avgMagGradU;                //!< sum of velocity grad magnitude over time
    Vec           avgMagUU;                   //!< sum of vector (u^2+v^2+w^2)*ui over time

    // phase averaging
    Vec           pAvgU;                      //!< phase sum of u, v, w over time
    Vec           pAvgP;                      //!< phase sum of p over time
    Vec           pAvgUU;                     //!< phase sum of resolved R stresses over time
    Vec           pAvgP2;                     //!< phase sum of p^2 over time
    Vec           pAvgNut;                    //!< phase sun of nut over time
    Vec           pAvgCs;                     //!< phase sum of cs over time
    Vec           pAvgOmega;                  //!< phase sum of vorticity over time
    Vec           pAvgOmegaOmega;             //!< phase sum of resolved vorticity stresses over time
    Vec           pAvgUdotGradP;              //!< phase sum of scalar u*dpdx + v*dpdy + w*dpdz over time
    Vec           pAvgMagGradU;               //!< phase sum of velocity grad magnitude over time
    Vec           pAvgMagUU;                  //!< phase sum of vector (u^2+v^2+w^2)*ui over time
};

#endif
