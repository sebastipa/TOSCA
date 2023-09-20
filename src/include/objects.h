//! \file objects.h
//! \brief Contains forward object declarations

//! \brief Access database context
struct access_;

//! \brief Solution flags context
struct flags_;

//! \brief Defined in domain.h: domain context
struct domain_;

//! \brief Defined in mesh.h: mesh context
struct mesh_;

//! \brief Defined in io.h: write settings context
struct io_;

//! \brief Defined in ueqn.h: momentum equation context
struct ueqn_;

//! \brief Defined in peqn.h: pressure equation context
struct peqn_;

//! \biref Defined in teqn.h: potential temperature transport equation context
struct teqn_;

//! \brief Defined in lesmdl.h: LES turbulence model context
struct les_;

//! \brief Defined in ADM.h: wind farm/turbine models context
struct farm_;

//! \brief Defined in abl.h: atmospheric boundary layer context
struct abl_;

//! \brief Defined in acquisition.h: simulation data gathering context
struct acquisition_;

//! \brief Defined in ibm.h: immersed boundary contect
struct ibm_;

//! \brief Defined in overset.h: overset mesh contect
struct overset_;

//! \brief Defined in precursor.h: precursor domain context
struct precursor_;
