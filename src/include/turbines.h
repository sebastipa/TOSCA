//! \file  ADM.h
//! \brief Wind Farm Objects,
//! for the arrays of std::string and arrays of structs we use dynamic allocation
//! through array of pointers to each std::string or struct respectively.

#ifndef FARM_H
#define FARM_H

//! \brief Airfoil info
typedef struct
{
    word                    name;   //!< name of the airfoil
    PetscInt                size;   //!< number of points in the look up tables
    PetscReal               *aoa;   //!< angles of attack
    PetscReal                *cl;   //!< lift coefficients
    PetscReal                *cd;   //!< drag coefficients
}foilInfo;

//! \brief Blade aerodynamic properties (used to build the AD/AL models)
typedef struct
{
    PetscInt                size;   //!< number of points
    PetscReal            *radius;   //!< radius stations
    PetscReal             *chord;   //!< chord distribution
    PetscReal             *twist;   //!< twist distribution
    PetscInt            *foilIds;   //!< airfoil ids as provided in the dict, start from 0
}bladeAeroInfo;

//! \brief Structure containing a coarser AD mesh located 2.5 D upstrem each turbine for velocity sampling
typedef struct
{
    // mesh level preoperties
    PetscInt             nPoints;   //!< total number of points
    Cmpnts               *points;   //!< array containing the upstream point coordinates
    PetscReal                *dA;   //!< array containing the element area at each upstream sample point (sums up to rotor area)

    Cmpnts                center;   //!< center point of the sampling mesh

    PetscInt   thisRigControlled;   //!< true if this processor has controller cells, zero otherwise
    MPI_Comm            UPW_COMM;   //!< communicator for this sampling (TRB_COMM is a subset of UPW_COMM: all data contained in UPW_COMM are accessible from TRB_COMM)

    cellIds     *controlledCells;   //!< labels of the background mesh cells influenced by this sample points in this processor
    PetscInt         nControlled;   //!< size of controlledCells

    // time-varying variables
    PetscInt   *thisPtControlled;   //!< flags telling if a point is controlled by this processor
    cellIds        *closestCells;   //!< indices of the closest cells to this sample points

    PetscReal               Uref;   //!< wind velocity averaged on the sample points

} upSampling;

//! \brief Actuator Farm Model
typedef struct
{
    Cmpnts                 point;   //!< point coordinates
    PetscReal               Uref;   //!< reference velocity to compute CtInf (only used for data writing)
    PetscReal                 Ct;   //!< imposed thrust coefficient

    word              sampleType;   //!< velocity sampling type ("rotorDisk" or "givenVelocity")
    PetscReal     rtrUFilterFreq;   //!< frequency of the single-pole low pass filter for the rotor wind velocity (if sampling type is "rotorDisk")

    // time-varying variables
    cellIds          closestCell;   //!< indices of the closest cells to this turbine AL points
    Cmpnts                     U;   //!< flow velocity at the AF point
    Cmpnts                     B;   //!< body force at the AF point
    PetscReal             axialF;   //!< rotor axial force at the AF point

    PetscReal          rtrThrust;   //!< total rotor thrust
    PetscReal            aeroPwr;   //!< total rotor aero power

    PetscInt    thisPtControlled;   //!< flag telling if this processor controls the AF point
    PetscInt          searchDone;   //!< flag telling if the closest cell search has been done

    // debug switch
    PetscInt                 dbg;   //!< prints a lot of information

} AFM;

//! \brief Actuator Line Model
typedef struct
{
    // mesh level preoperties
    PetscInt             nPoints;   //!< total number of points
    Cmpnts               *points;   //!< array containing the AL point coordinates
    PetscReal                *dr;   //!< array containing the radial mesh size

    PetscInt             nRadial;   //!< number of AL points in the radial direction
    PetscInt            nAzimuth;   //!< number of AL points in the azimuthal direction

    PetscReal               Uref;   //!< reference velocity to compute CtInf (only used for data writing)

    PetscReal             *chord;   //!< array containing the chord at each AL point
    PetscReal             *twist;   //!< array containing the twist at each AL point
    PetscReal          *solidity;   //!< array containing the solidity at each AL point
    PetscInt           **foilIds;   //!< labels of the 2 airfoils closest to each AL point
    PetscReal               **iw;   //!< interp. weights for the the 2 airfoils closest to each AL point

    // time-varying variables
    PetscInt   *thisPtControlled;   //!< flags telling if a point is controlled by this processor
    cellIds        *closestCells;   //!< indices of the closest cells to this turbine AL points
    PetscReal                *Cd;   //!< drag coefficient at each point of the AL
    PetscReal                *Cl;   //!< lift coefficient at each point of the AL
    PetscReal             *alpha;   //!< angle of attack at each point of the AL
    Cmpnts                    *U;   //!< flow velocity at each point of the AL (relative to the blade)
    Cmpnts                    *B;   //!< body force at each point of the AL mesh
    PetscReal            *axialF;   //!< rotor axial force at each point of the AL mesh
    PetscReal            *tangtF;   //!< rotor tangential force at each point of the AL mesh

    PetscReal         rtrAvgMagU;   //!< average velocity on the rotor (includes induction from CFD)
    PetscReal          rtrTorque;   //!< total rotor torque
    PetscReal          rtrThrust;   //!< total rotor thrust
    PetscReal            aeroPwr;   //!< total rotor aero power

    PetscReal            azimuth;   //!< azimuthal angle in deg

    // debug switch
    PetscInt                 dbg;   //!< prints a lot of information

} ALM;

//! \brief Actuator Disk Model
typedef struct
{
    // mesh level preoperties
    PetscInt             nPoints;   //!< total number of points
    Cmpnts               *points;   //!< array containing the AD point coordinates
    PetscReal                *dr;   //!< array containing the radial mesh size

    PetscInt             nRadial;   //!< number of AD points in the radial direction
    PetscInt            nAzimuth;   //!< number of AD points in the azimuthal direction

    PetscReal               Uref;   //!< reference velocity to compute CtInf (only used for data writing)

    PetscReal             *chord;   //!< array containing the chord at each AD point
    PetscReal             *twist;   //!< array containing the twist at each AD point
    PetscReal          *solidity;   //!< array containing the solidity at each AD point
    PetscInt           **foilIds;   //!< labels of the 2 airfoils closest to each AD point
    PetscReal               **iw;   //!< interp. weights for the the 2 airfoils closest to each AD point

    // time-varying variables
    PetscInt   *thisPtControlled;   //!< flags telling if a point is controlled by this processor
    cellIds        *closestCells;   //!< indices of the closest cells to this turbine AD points
    PetscReal                *Cd;   //!< drag coefficient at each point of the AD
    PetscReal                *Cl;   //!< lift coefficient at each point of the AD
    PetscReal             *alpha;   //!< angle of attack at each point of the AD
    Cmpnts                    *U;   //!< flow velocity at each point of the AD (relative to the blade)
    Cmpnts                    *B;   //!< body force at each point of the AD mesh
    PetscReal            *axialF;   //!< rotor axial force at each point of the AD mesh
    PetscReal            *tangtF;   //!< rotor tangential force at each point of the AD mesh

    PetscReal         rtrAvgMagU;   //!< average velocity on the rotor (includes induction from CFD)
    PetscReal          rtrTorque;   //!< total rotor torque
    PetscReal          rtrThrust;   //!< total rotor thrust
    PetscReal            aeroPwr;   //!< total rotor aero power

    // debug switch
    PetscInt                 dbg;   //!< prints a lot of information

} ADM;

//! \brief Uniform Actuator Disk model
typedef struct
{
    // mesh level preoperties
    PetscInt             nPoints;   //!< total number of points
    Cmpnts               *points;   //!< array containing the AD point coordinates
    PetscReal                *dA;   //!< array containing the element area at each AD point (sums up to rotor area)

    PetscInt             nRadial;   //!< number of AD points in the radial direction
    PetscInt            nAzimuth;   //!< number of AD points in the azimuthal direction

    word              sampleType;   //!< velocity sampling type ("rotorUpstream" or "givenVelocity")
    PetscReal               Uref;   //!< reference velocity to make the provided Ct dimensional
    PetscReal                 Ct;   //!< imposed thrust coefficient
    PetscReal             axiInd;   //!< axial induction factor

    // time-varying variables
    PetscInt   *thisPtControlled;   //!< flags telling if a point is controlled by this processor
    cellIds        *closestCells;   //!< indices of the closest cells to this turbine AD points
    Cmpnts                    *U;   //!< flow velocity at each point of the AD
    Cmpnts                    *B;   //!< body force at each point of the AD mesh
    PetscReal            *axialF;   //!< rotor axial force at each point of the AD mesh

    PetscReal         rtrAvgMagU;   //!< average velocity on the rotor (includes induction from CFD)
    PetscReal          rtrThrust;   //!< total rotor thrust
    PetscReal            aeroPwr;   //!< total rotor aero power

    // debug switch
    PetscInt                 dbg;   //!< prints a lot of information
} UADM;

//! \brief Wind turbine tower actuator line model
typedef struct
{
    PetscReal              rBase;   //!< radius at the base of the tower
    PetscReal               rTop;   //!< radius at the top of the tower
    PetscReal                 Cd;   //!< tower drag coefficient
    PetscReal                eps;   //!< spreading width of the gaussian projection function (good is 0.035 * hTower)
    PetscReal          prjNSigma;   //!< confidence interval as number of std deviations for the projection function (hardcoded to 2.7)

    PetscInt             nPoints;   //!< number of tower points in the linear direction

    Cmpnts               *points;   //!< total number of points
    PetscReal                *dA;   //!< array containing the frontal area at each tower point

    cellIds     *controlledCells;   //!< labels of the background mesh cells influenced by this tower in this processor
    PetscInt         nControlled;   //!< size of controlledCells

    PetscInt   *thisPtControlled;   //!< flags telling if a point is controlled by this processor
    cellIds        *closestCells;   //!< indices of the closest cells to this turbine tower points

    Cmpnts                    *U;   //!< flow velocity at each point of the tower
    Cmpnts                    *B;   //!< body force at each point of the tower mesh
    PetscReal             *tangF;   //!< tower tangential force at each point of the tower mesh

    PetscReal          twrThrust;   //!< total tower thrust

    MPI_Comm            TWR_COMM;   //!< communicator for this tower
    PetscMPIInt        nProcsTwr;   //!< size of the TWR_COMM communicator

} towerModel;

//! \brief Wind turbine nacelle mounted anemometer
typedef struct
{
    Cmpnts            samplePoint;  //!< position of the anemometer
    cellIds           closestCell;  //!< cell from which to sample the velocity
    PetscInt anemometerControlled;  //!< flag telling if this processor controls the anemometer
    Cmpnts                      U;  //!< instantaneous velocity
    Cmpnts         anemFilterFreq;  //!< frequency of the single-pole low pass filter for U
    Cmpnts        acquisitionFreq;  //!< anemometer acquisition frequency
} anemometer;

//! \brief Wind turbine parameters
typedef struct
{
    // global parameters (all turbine models)
    word                      id;   //!< id of the wind turbine
    word                    type;   //!< type of the wind turbine
    PetscInt             nBlades;   //!< number of turbine blades
    PetscReal               rTip;   //!< tip radius from CoR
    PetscReal               rHub;   //!< hub radius from CoR
    PetscReal               hTwr;   //!< tower height
    PetscReal            ovrHang;   //!< nacelle overhang in the rotor direction (facing the wind)
    PetscReal            precone;   //!< blade precone (equal for all blades)
    Cmpnts                twrDir;   //!< unit vector pointing from base to tower top
    PetscReal             upTilt;   //!< nacell up-tilt (positive when blades get far from tower)
    PetscReal             genEff;   //!< electrical generator eficiency

    // detailed parameters (ADM and ALM models)
    word                  rotDir;   //!< rotation dir as seen lookingthrough the WT from the front
    Cmpnts             rotCenter;   //!< rotor rotation center
    PetscInt              nFoils;   //!< number of airfoils in the database
    word             **foilNames;   //!< array of pointers of size(n-foils in turbine) to airfoil names
    foilInfo             **foils;   //!< array of pointers of size(n-foils in turbine) to each airfoil's info
    bladeAeroInfo          blade;   //!< blade properties

    anemometer              WDAS;   //!< wind data acquisition system (nacelle mounted anemometer)

    // time-varying variables
    Cmpnts                rtrDir;   //!< (all) unit vector pointing to the rotor orientation in non-tilted position (facing the wind)
    Cmpnts               rtrAxis;   //!< (all) unit vector pointing to the rotor orientation in tilted position (facing the wind)
    Cmpnts             omega_hat;   //!< (AD/AL) turbine angular velocity unit vector (directed as rtrAxis, pointed according to rotDir)
    PetscReal           rtrOmega;   //!< (AD/AL) turbine angular velocity in rad/sec

    // wind turbine models
    ADM                      adm;   //!< actuator disk model
    UADM                    uadm;   //!< unform actuator disk model
    ALM                      alm;   //!< actuator line model
    AFM                      afm;   //!< actuator farm model

    // tower model
    towerModel               twr;   //!< actuator line tower model
    PetscInt          includeTwr;   //!< flag telling if tower is included

    // sampling points 2.5 RD upstream
    upSampling         *upPoints;   //!< struct containing the upstream sampling points information

    // useful conversion factors
    PetscReal            deg2rad;   //!< degrees to radiants conversion factor
    PetscReal            rad2deg;   //!< radiants to degrees conversion factor
    PetscReal         rpm2RadSec;   //!< RPM to rad/s conversion factor

    // flags
    PetscInt   turbineControlled;   //!< flag which tells if this proc controls this wind turbine

    // communication color
    MPI_Comm            TRB_COMM;   //!< communicator for this turbine
    PetscMPIInt       writerRank;   //!< label of master rank of the TRB_COMM communicator in the MPI_COMM_WORLD rank list
    PetscMPIInt        nProcsTrb;   //!< size of the TRB_COMM communicator

    // torque controller (ADM and ALM models)
    word       genControllerType;   //!< name of torque controller (if none preserves intial omega, else reads from control/genControllerType)
    PetscReal           genOmega;   //!< turbine angular velocity
    PetscReal       rtrOmegaFilt;   //!< turbine filtered angular velocity
    PetscReal   rtrSpdFilterFreq;   //!< frequency of the single-pole low pass filter for the rotor angular frequency
    PetscReal        cutInGenSpd;   //!< cut in generator speed
    PetscReal         cutInGenTq;   //!< cut in generator torque
    PetscReal  regTwoStartGenSpd;   //!< generator speed at the start of control region 2
    PetscReal    regTwoEndGenSpd;   //!< generator speed at the end of control region 2
    PetscReal         ratedGenTq;   //!< generator torque at the rated wind speed
    PetscReal            omegaKP;   //!< proportional gain of the generator torque controller

    PetscReal          genTorque;   //!< total generator torque
    PetscReal             genPwr;   //!< total rotor gen power

    PetscInt       tqRateLimiter;   //!< activate torque rate limiter (1 yes, 0 no)
    PetscInt       rtrSpdLimiter;   //!< activate generator speed limiter (1 yes, 0 no)
    PetscReal          tqMaxRate;   //!< maximum torque variation rate allowed for the generator torque controller
    PetscReal      ratedRotorSpd;   //!< rotor speed at rated wind speed

    // rotor dynamics (if torque controller is active)
    PetscReal  driveTrainInertia;    //!< sum of all the inertias attached to the shaft
    PetscReal         genInertia;    //!< generator inertia
    PetscReal         hubInertia;    //!< hub intertia
    PetscReal         bldInertia;    //!< blade intertia
    PetscReal        gbxRatioG2R;    //!< gearbox generator-to-rotor ratio
    PetscReal             gbxEff;    //!< gearbox mechanical efficiency

    // pitch controller (ADM and ALM models)
    word     pitchControllerType;   //!< name of torque controller (if none preserves intial omega, else reads from control/pitchControllerType)
    PetscReal             *pitch;   //!< blade pitch (array of size n-blades), could vary blade by blade
    PetscReal          collPitch;   //!< collective pitch (same for all blades)
    PetscReal            pitchKP;   //!< pitch controller proportional gain
    PetscReal            pitchKI;   //!< pitch controller integral gain
    PetscReal            pitchKD;   //!< pitch controller derivative gain
    PetscReal           pitchS2R;   //!< pitch at which the sensit. of power to pitch variations has doubled w.r.t. rated position
    PetscReal             errPID;   //!< error of the PID controller
    PetscReal          intErrPID;   //!< integrated error of the PID controller

    PetscInt    pitchRateLimiter;   //!< activate pitch rate limiter (1 yes, 0 no)
    PetscInt   pitchAngleLimiter;   //!< activate pitch angle limiter (1 yes, 0 no)
    PetscReal       pitchMaxRate;   //!< max rate of blade pitch motion
    PetscReal           pitchMin;   //!< minimum pitch angle
    PetscReal           pitchMax;   //!< maximum pitch angle

    // yaw controller (ADM and ALM models)
    word       yawControllerType;   //!< name of torque controller (if none preserves intial omega, else reads from control/yawControllerType)
    word         yawSamplingType;   //!< mode of velocity sampling (hubUpDist, anemometer)
    PetscReal   yawAverageWindow;   //!< time window used for misalignment averaging
    PetscReal    yawAllowedError;   //!< allows +-yawAllowedError flow misalignment in degrees
    PetscReal             yawMin;   //!< minimum yaw angle
    PetscReal             yawMax;   //!< maximum yaw angle
    PetscReal           yawAngle;   //!< actual yaw angle wrt xyz background ref frame
    PetscReal          flowAngle;   //!< actual averaged flow angle wrt xyz background ref frame
    PetscReal           yawError;   //!< yaw misalignment error
    PetscReal           yawSpeed;   //!< yaw speed in degs/s
    cellIds         yawSampleIds;   //!< ids of the point where the velocity for misalignment computation must be sampled
    PetscInt          yawChanged;   //!< flag telling if must do the search on the turbine points due to yaw change

    // numerical parameters
    cellIds     *controlledCells;   //!< labels of the background mesh cells influenced by this turbine in this processor
    PetscInt         nControlled;   //!< size of controlledCells
    PetscReal                eps;   //!< spreading width of the gaussian projection function (good is 0.035 * dBlade)
    PetscReal              eps_x;   //!< x spreading width of the gaussian projection function (for AFM)
    PetscReal              eps_y;   //!< y spreading width of the gaussian projection function (for AFM)
    PetscReal              eps_z;   //!< z spreading width of the gaussian projection function (for AFM)
    PetscReal          prjNSigma;   //!< confidence interval as number of std deviations for the projection function (hardcoded to 2.7)

    // debug switch
    PetscInt                 dbg;   //!< this turbines info at the begining of the simulation and at each iteration

} windTurbine;

//! \brief Wind farm parameters
struct farm_
{
    // main parametes
    word                    name;   //!< wind farm name
    PetscInt                size;   //!< number of turbines in the farm
    Cmpnts                 *base;   //!< base coordinates of each turbine

    word          **turbineTypes;   //!< array of pointers to type of the wind turbines in the farm
    word            **turbineIds;   //!< array of pointers to ID of the wind turbines in the farm
    word         **turbineModels;   //!< array of pointers to name of the wind turbine models in the farm

    windTurbine             **wt;   //!< array of pointers to wind turbines

    // output parameters
    PetscReal          timeStart;  //!< start time of acquisition system
    word            intervalType;  //!< timeStep: sample at every (timeInterval) iter, adjustableTime sample at every (timeInterval) seconds
    PetscReal       timeInterval;  //!< acquisition time interval (overrides simulation time step if smaller and adjustableTime active)
    PetscInt         writeNumber;  //!< number of mesh files written up to now

    // debug switch
    PetscInt                 dbg;  //!< global debug switch, checks for point discrimination algorithm (slower)

    Vec           lsourceFarmCat,  //!< cartesian wind farm body force
                  sourceFarmCont;  //!< contravariant wind farm body force

    // access database
    access_              *access;

    // CFL control for ALM: blade tip cannot move more than one cell
    PetscInt       checkCFL;       //!< at least one actuator line model is present in the wind farm
    PetscReal      maxTipSpeed;    //!< maximum tip speed among all rotors
};

#endif

// Top level functions
// -----------------------------------------------------------------------------

//! \brief Initialize the wind farm
PetscErrorCode InitializeWindFarm(farm_ *farm);

//! \brief Update wind turbines
PetscErrorCode UpdateWindTurbines(farm_ *farm);

// Actuator disk model routines.
// -----------------------------------------------------------------------------

//! \brief Solve rotor dynamics and compute filtered rot speed
PetscErrorCode computeRotSpeed(farm_ *farm);

//! \brief Rotate blades (only for ALM)
PetscErrorCode rotateBlades(windTurbine *wt, PetscReal angle);

//! \brief Compute generator torque with 5-regions control system model
PetscErrorCode controlGenSpeed(farm_ *farm);

//! \brief Compute blade pitch using the PID control
PetscErrorCode controlBldPitch(farm_ *farm);

//! \brief Compute nacelle yaw
PetscErrorCode controlNacYaw(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each rotor are controlled by this processor
PetscErrorCode findControlledPointsRotor(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each tower are controlled by this processor
PetscErrorCode findControlledPointsTower(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each sample rig are controlled by this processor
PetscErrorCode findControlledPointsSample(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on the rotor
PetscErrorCode checkPointDiscriminationRotor(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on the tower
PetscErrorCode checkPointDiscriminationTower(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on sample points
PetscErrorCode checkPointDiscriminationSample(farm_ *farm);

//! \brief Compute wind velocity at the rotor mesh points
PetscErrorCode computeWindVectorsRotor(farm_ *farm);

//! \brief Compute wind velocity at the tower mesh points
PetscErrorCode computeWindVectorsTower(farm_ *farm);

//! \brief Compute wind velocity at the sample mesh points
PetscErrorCode computeWindVectorsSample(farm_ *farm);

//! \brief Compute aerodynamic forces at the turbine mesh points
PetscErrorCode computeBladeForce(farm_ *farm);

//! \brief Project the wind turbine forces on the background mesh
PetscErrorCode projectBladeForce(farm_ *farm);

//! \brief Compute and project the tower forces on the background mesh
PetscErrorCode projectTowerForce(farm_ *farm);

//! \brief Transform the cartesian body force to contravariant
PetscErrorCode bodyForceCartesian2Contravariant(farm_ *farm);

//! \brief Write output if applicable
PetscErrorCode windTurbinesWrite(farm_ *farm);

//! \brief Write checkpoint file
PetscErrorCode windTurbinesWriteCheckpoint(farm_ *farm);

//! \brief Read checkpoint file and prepare wind turbines
PetscErrorCode windTurbinesReadCheckpoint(farm_ *farm);

//! \brief Print wind farm information
PetscErrorCode printFarmProperties(farm_ *farm);

//! \brief Write the wind farm AD mesh to ucd file
PetscErrorCode writeFarmADMesh(farm_ *wf);

//! \brief Write the wind farm AL mesh to ucd file
PetscErrorCode writeFarmALMesh(farm_ *farm);

//! \brief Write the wind farm tower mesh to ucd file
PetscErrorCode writeFarmTwrMesh(farm_ *wf);

//! \brief Initialize sphere cells and see what proc controls which turbine
PetscErrorCode initControlledCells(farm_ *farm);

//! \brief Initialize sphere cells and see what proc controls which sample point
PetscErrorCode initSampleControlledCells(farm_ *farm);

// \brief Initializes the ADM reading from files and allocating memory
PetscErrorCode initADM(windTurbine *wt, Cmpnts &base, const word meshName);

//! \brief Initialize the Uniform Actuator Disk Model
PetscErrorCode initUADM(windTurbine *wt, Cmpnts &base, const word meshName);

//! \brief Initialize the upstream sample points data structure
PetscErrorCode initSamplePoints(windTurbine *wt, Cmpnts &base, const word meshName);

//! \brief Initializes the ALM reading from files and allocating memory
PetscErrorCode initALM(windTurbine *wt, Cmpnts &base, const word meshName);

//! \brief Initializes the AFM reading from files and allocating memory
PetscErrorCode initAFM(windTurbine *wt, Cmpnts &base, const word meshName);

//! \brief Initializes the tower model
PetscErrorCode initTwrModel(windTurbine *wt, Cmpnts &base);

//! \brief Compute max tip speed and activate CFL control flag
PetscErrorCode computeMaxTipSpeed(farm_ *farm);

// Reading input files functions
// -----------------------------------------------------------------------------

//! \brief Read the farm_Properties file and fill the structs
PetscErrorCode readFarmProperties(farm_ *farm);

//! \brief Read the wind turbines inside the farm_Properties file
PetscErrorCode readTurbineArray(farm_ *wf);

//! \brief Fill the windTurbine struct reading in the file named as the wind turbine type
PetscErrorCode readTurbineProperties(windTurbine *wt, const char *dictName, const word meshName, const word modelName);

//! \brief Reads the names of the airfoil used in the turbine (given in the airfoils subdict inside the file named as the wind turbine type)
PetscErrorCode readAirfoilProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the blades aero properties used in the turbine (given in the bladeData subdict inside the file named as the wind turbine type)
PetscErrorCode readBladeProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the tower properties used in the turbine (given in the towerData subdict inside the file named as the wind turbine)
PetscErrorCode readTowerProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the 2D airfoil tables given in the file named as the airfoil
PetscErrorCode readAirfoilTable(foilInfo *af, const char *tableName);

//! \brief Reads generator torque controller parameters
PetscErrorCode readGenControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Reads blade pitch controller parameters
PetscErrorCode readPitchControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Reads nacelle yaw controller parameters
PetscErrorCode readYawControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

// Mathematical operations
// -----------------------------------------------------------------------------

// \brief Returns interpolation weights of the two closest values to pval and their indices in pvec
void findInterpolationWeigths(PetscReal *weights, PetscInt *labels, PetscReal *pvec, PetscInt npts, PetscReal pval);
