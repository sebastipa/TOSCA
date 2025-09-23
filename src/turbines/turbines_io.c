#include "turbines_io.h"
#include "turbines_phys.h"
#include "turbines_openfast.h"

//***************************************************************************************************************//

PetscErrorCode windTurbinesWrite(farm_ *farm)   
{
    mesh_ *mesh = farm->access->mesh;

    PetscMPIInt rank, t;

    clock_ *clock = farm->access->clock;
    io_    *io    = farm->access->io;

    word   turbineFolderName     = "./postProcessing/" + mesh->meshName + "/turbines";
    word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // create/initialize turbines directory (at simulation start only)
    if
    (
        clock->it == clock->itStart
    )
    {
        if(!rank)
        {
            errno = 0;
            PetscInt dirRes;

            dirRes = mkdir(turbineFolderName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create %s directory\n", turbineFolderName.c_str());
                fatalErrorInFunction("windTurbinesWrite",  error);
            }

            dirRes = mkdir(turbineFolderTimeName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create %s directory\n", turbineFolderTimeName.c_str());
                fatalErrorInFunction("windTurbinesWrite",  error);
            }

            // if directory already exist remove everything inside
            if(errno == EEXIST)
            {
                remove_subdirs(farm->access->mesh->MESH_COMM, turbineFolderTimeName.c_str());
            }
        }

        // ensure folder is there for every processor
        MPI_Barrier(mesh->MESH_COMM);

        // create turbines files
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(rank == wt->writerRank)
            {
                FILE *f;
                char fileName[80];
                sprintf(fileName, "%s/%s", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                f = fopen(fileName, "w");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName);
                    fatalErrorInFunction("windTurbinesWrite",  error);
                }
                else
                {

                    int width = -20;

                    word w0 = "time [s]";

                    word w1  = "rtrAvgMagU [m/s]";
                    word w2  = "rtrAvgUpMagU [m/s]";
                    word w3  = "rtrThrust [kN]";
                    word w4  = "aeroPwr [MW]";

                    word w5  = "CtInf [-]";
                    word w6  = "CtLoc [-]";
                    word w7  = "CtUp [-]";

                    word w8  = "rtrTorque [kNm]";
                    word w9  = "rtrOmega [rpm]";

                    word w10 = "genTorque [kNm]";
                    word w11 = "genPwr [MW]";
                    word w12 = "genOmega [rpm]";

                    word w13 = "collPitch [deg]";

                    word w14 = "flowAngle [deg]";
                    word w15 = "yawAngle [deg]";

                    word w16 = "azimuth [deg]";

                    // actuator disk model
                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str(), width, w8.c_str(), width, w9.c_str());

                        if(wt->genControllerType != "none")
                        {
                            fprintf(f, "%*s %*s %*s ", width, w10.c_str(), width, w11.c_str(), width, w12.c_str());
                        }

                        if(wt->pitchControllerType != "none")
                        {
                            fprintf(f, "%*s ", width, w13.c_str());
                        }
                    }
                    // uniform actuator disk model
                    else if((*farm->turbineModels[t]) == "uniformADM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str());
                    }
                    // actuator farm model
                    else if((*farm->turbineModels[t]) == "AFM")
                    {
                        fprintf(f, "%*s %*s %*s ", width, w0.c_str(), width, w3.c_str(), width, w4.c_str());
                    }
                    else if((*farm->turbineModels[t]) == "ALM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str(), width, w8.c_str(), width, w9.c_str(), width, w16.c_str());

                        if(wt->genControllerType != "none")
                        {
                            fprintf(f, "%*s %*s %*s ", width, w10.c_str(), width, w11.c_str(), width, w12.c_str());
                        }

                        if(wt->pitchControllerType != "none")
                        {
                            fprintf(f, "%*s ", width, w13.c_str());
                        }
                    }

                    // yaw controller is not model specific
                    if(wt->yawControllerType != "none")
                    {
                        fprintf(f, "%*s %*s ", width, w14.c_str(), width, w15.c_str());
                    }

                    fprintf(f, "\n");

                    fclose(f);
                }

                if((*farm->turbineModels[t]) == "ALM")
                {
                    //create file to write the turbine radial point angle of attack values
                    sprintf(fileName, "%s/%s_aoa", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                        char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point relative velocity magnitude
                    sprintf(fileName, "%s/%s_relVelMag", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point rel velocity angle phi
                    sprintf(fileName, "%s/%s_phi", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in axial direction
                    sprintf(fileName, "%s/%s_uAxl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in radial direction
                    sprintf(fileName, "%s/%s_uRad", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in tangential direction
                    sprintf(fileName, "%s/%s_uTan", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point Cl
                    sprintf(fileName, "%s/%s_Cl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point Cd
                    sprintf(fileName, "%s/%s_Cd", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point axial force
                    sprintf(fileName, "%s/%s_axialF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point tangential force
                    sprintf(fileName, "%s/%s_tangF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }
                }
            }
        }
    }

    // see if must write to file
    PetscInt  turbinesWrite = 0;
    word      intervalType  = farm->intervalType;
    PetscReal timeInterval  = farm->timeInterval;
    PetscReal timeStart     = farm->timeStart;

    // write every "timeInterval" seconds
    if
    (
        (intervalType == "adjustableTime") &&
        (
            (clock->time - timeStart) / timeInterval -
            std::floor
            (
                (clock->time - timeStart) / timeInterval
            ) < 1e-10
        )
    )
    {
        turbinesWrite = 1;
    }
    // write every "timeInterval" iterations
    else if
    (
        (clock->it > 0) &&
        (intervalType == "timeStep") &&
        (
            clock->it / timeInterval -
            std::floor
            (
                clock->it / timeInterval
            ) < 1e-10
        )
    )
    {
        turbinesWrite = 1;
    }

    if(turbinesWrite)
    {
        constants_ *constants = farm->access->constants;

        // loop over wind turbines
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->turbineControlled)
            {
                PetscReal Ar          = M_PI * pow(wt->rTip,2.0);

                // UADM/ADM
                PetscReal rtrAvgMagU  = 0.0;
                PetscReal rtrAvgUpMagU= 0.0;
                PetscReal rtrThrust   = 0.0;
                PetscReal aeroPwr     = 0.0;
                PetscReal ctInf       = 0.0;
                PetscReal ctLoc       = 0.0;
                PetscReal ctUp        = 0.0;

                // ADM
                PetscReal rtrTorque   = 0.0;
                PetscReal rtrOmega    = 0.0;

                // ALM
                PetscReal azimuth     = 0.0;

                // controls
                PetscReal genTorque   = 0.0;
                PetscReal genPwr      = 0.0;
                PetscReal genOmega    = 0.0;
                PetscReal collPitch   = 0.0;

                // do turbine level scatter first
                if((*farm->turbineModels[t]) == "ADM")
                {
                    MPI_Reduce(&(wt->adm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->adm.rtrTorque),  &rtrTorque,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->rtrOmega),       &rtrOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genOmega),    &genOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        genTorque    = genTorque  / wt->nProcsTrb / 1.0e3;
                        genPwr       = genPwr     / wt->nProcsTrb / 1.0e6;
                        genOmega     = genOmega   / wt->nProcsTrb / wt->rpm2RadSec;
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        collPitch    = collPitch  / wt->nProcsTrb;
                    }

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    rtrTorque    = rtrTorque  / wt->nProcsTrb / 1.0e3;
                    rtrOmega     = rtrOmega   / wt->nProcsTrb / wt->rpm2RadSec;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->adm.Uref,   2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    MPI_Reduce(&(wt->uadm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->uadm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->uadm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->uadm.Uref,  2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }
                else if((*farm->turbineModels[t]) == "AFM")
                {
                    MPI_Reduce(&(wt->afm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->afm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;
                }
                else if((*farm->turbineModels[t]) == "ALM")
                {
                    MPI_Reduce(&(wt->alm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->alm.rtrTorque),  &rtrTorque,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->rtrOmega),       &rtrOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->alm.azimuth),    &azimuth,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genOmega),    &genOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        genTorque    = genTorque  / wt->nProcsTrb / 1.0e3;
                        genPwr       = genPwr     / wt->nProcsTrb / 1.0e6;
                        genOmega     = genOmega   / wt->nProcsTrb / wt->rpm2RadSec;
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        collPitch    = collPitch  / wt->nProcsTrb;
                    }

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    rtrTorque    = rtrTorque  / wt->nProcsTrb / 1.0e3;
                    rtrOmega     = rtrOmega   / wt->nProcsTrb / wt->rpm2RadSec;

                    azimuth      = azimuth    / wt->nProcsTrb;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->alm.Uref,   2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }

                // write to file
                if(rank == wt->writerRank)
                {
                    FILE *f;
                    char fileName[80];
                    sprintf(fileName, "%s/%s", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "a");

                    if(!f)
                    {
                        char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = -20;

                        // actuator disk model
                        if((*farm->turbineModels[t]) == "ADM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp, width, rtrTorque, width, rtrOmega);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "%*.4f %*.4f %*.4f ", width, genTorque, width, genPwr, width, genOmega);
                            }

                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "%*.4f ", width, collPitch*wt->rad2deg);
                            }
                        }
                        // uniform actuator disk model
                        else if((*farm->turbineModels[t]) == "uniformADM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp);
                        }
                        // actuator farm model
                        else if((*farm->turbineModels[t]) == "AFM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f ", width, clock->time, width, rtrThrust, width, aeroPwr);
                        }
                        else if((*farm->turbineModels[t]) == "ALM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp, width, rtrTorque, width, rtrOmega, width, azimuth);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "%*.4f %*.4f %*.4f ", width, genTorque, width, genPwr, width, genOmega);
                            }

                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "%*.4f ", width, collPitch*wt->rad2deg);
                            }
                        }

                        // yaw controller is not model specific
                        if(wt->yawControllerType != "none")
                        {
                            fprintf(f, "%*.4f %*.4f ", width, wt->flowAngle, width, wt->yawAngle);
                        }

                        fprintf(f, "\n");

                        fclose(f);
                    }

                    if((*farm->turbineModels[t]) == "ALM")
                    {
                        sprintf(fileName, "%s/%s_aoa", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sumAlpha = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sumAlpha += wt->alm.alpha[nAz*ri + ai];
                                }

                                sumAlpha/=nAz;

                                fprintf(f, "%*.4f", width, sumAlpha);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_relVelMag", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += nMag(wt->alm.U[nAz*ri + ai]);
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_phi", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai, p;
                            Cmpnts point_p, r_p;
                            Cmpnts xb_hat, yb_hat, zb_hat;
                            PetscReal ub_x, ub_y, phi;

                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    p =  nAz*ri + ai;

                                    // save this point locally for speed
                                    point_p = wt->alm.points[p];

                                    // this point position from COR
                                    r_p  = nSub(point_p, wt->rotCenter);

                                    if(wt->rotDir == "cw")
                                    {
                                        zb_hat = nUnit(r_p);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }
                                    else if(wt->rotDir == "ccw")
                                    {
                                        zb_hat = nUnit(r_p);
                                                 mScale(-1.0, zb_hat);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }

                                    ub_x = nDot(wt->alm.U[p], xb_hat);
                                    ub_y = nDot(wt->alm.U[p], yb_hat);

                                    phi = wt->rad2deg * std::atan2(ub_x, ub_y);
                                    sum += phi;
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        FILE *f1, *f2, *f3;
                        char fileName1[80], fileName2[80], fileName3[80];
                        sprintf(fileName1, "%s/%s_uAxl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f1 = fopen(fileName1, "a");
                        sprintf(fileName2, "%s/%s_uRad", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f2 = fopen(fileName2, "a");
                        sprintf(fileName3, "%s/%s_uTan", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f3 = fopen(fileName3, "a");
                        if(!f1 || !f2 || !f3)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName1);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f1, "%.4f", clock->time);
                            fprintf(f2, "%.4f", clock->time);
                            fprintf(f3, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai, p;
                            Cmpnts point_p, r_p;
                            Cmpnts xb_hat, yb_hat, zb_hat;
                            PetscReal ub_x, ub_y, ub_z, phi;

                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sumAxl = 0.0, sumTan = 0.0, sumRad = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    p =  nAz*ri + ai;

                                    // save this point locally for speed
                                    point_p = wt->alm.points[p];

                                    // this point position from COR
                                    r_p  = nSub(point_p, wt->rotCenter);

                                    if(wt->rotDir == "cw")
                                    {
                                        zb_hat = nUnit(r_p);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }
                                    else if(wt->rotDir == "ccw")
                                    {
                                        zb_hat = nUnit(r_p);
                                                 mScale(-1.0, zb_hat);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }

                                    ub_x = nDot(wt->alm.gWind[p], xb_hat);
                                    ub_y = nDot(wt->alm.gWind[p], yb_hat);
                                    ub_z = nDot(wt->alm.gWind[p], zb_hat);

                                    sumAxl += ub_x;
                                    sumTan += ub_y;
                                    sumRad += ub_z;
                                }

                                sumAxl/=nAz;
                                sumTan/=nAz;
                                sumRad/=nAz;

                                fprintf(f1, "%*.4f", width, sumAxl);
                                fprintf(f2, "%*.4f", width, sumRad);
                                fprintf(f3, "%*.4f", width, sumTan);
                            }

                            fprintf(f1, "\n");
                            fprintf(f2, "\n");
                            fprintf(f3, "\n");

                            fclose(f1);
                            fclose(f2);
                            fclose(f3);
                        }

                        sprintf(fileName, "%s/%s_Cl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.Cl[nAz*ri + ai];
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_Cd", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.Cd[nAz*ri + ai];
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_axialF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;
                                PetscReal sumDr = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.axialF[nAz*ri + ai];
                                    sumDr += wt->alm.dr[nAz*ri + ai];
                                }

                                sum/=nAz;
                                sumDr/=nAz;

                                //to compute the force per unit length divide by element radial length
                                sum/=sumDr;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_tangF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;
                                PetscReal sumDr = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.tangtF[nAz*ri + ai];
                                    sumDr += wt->alm.dr[nAz*ri + ai];
                                }

                                sum/=nAz;
                                sumDr/=nAz;

                                //to compute the force per unit length divide by element radial length
                                sum/=sumDr;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }
                    }

                }
            }
        }
    }

    // write wind farm mesh according to runTimeWrite option
    if(io->runTimeWrite)
    {
        if(farm->writeNumber == 0)
        {
            // write the wind farm mesh
            writeFarmTwrMesh(farm);
        }

        // write AD mesh
        writeFarmADMesh(farm);

        // write AL mesh
        writeFarmALMesh(farm);

        // increase write number counter
        farm->writeNumber++;
    }

    // write checkpoint file
    windTurbinesWriteCheckpoint(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windTurbinesWriteCheckpoint(farm_ *farm)
{
    // writes in fields directory which has been already created at sim start

    clock_      *clock = farm->access->clock;
    io_         *io    = farm->access->io;
    mesh_       *mesh  = farm->access->mesh;

    word        timeName;
    word        path, turbineFolderName;

    PetscInt    t;

    PetscMPIInt rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set time folder name
    timeName          = getTimeName(clock);
    path              = "./fields/" + mesh->meshName + "/turbines/" + timeName;
    turbineFolderName = "./fields/" + mesh->meshName + "/turbines";

    // create/initialize fields/turbines directory (at simulation start only)
    if
    (
        clock->it == clock->itStart && !rank
    )
    {
        errno = 0;
        PetscInt dirRes = mkdir(turbineFolderName.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create fields/turbines directory");
            fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
        }

        // if directory already exist remove everything inside except the start time (safe)
        if(errno == EEXIST)
        {
            word startTimeName = getStartTimeName(clock);
            remove_subdirs_except_keep_n(mesh->MESH_COMM, turbineFolderName.c_str(), startTimeName.c_str(), io->purgeWrite-1);
        }
    }

    // test if must write at this time step
    if(io->runTimeWrite)
    {
        // creates time folder
        if(!rank)
        {
            PetscInt dirRes = mkdir(path.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory\n", path.c_str());
                fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
            }
        }

        // ensure folder is there for every processor
        MPI_Barrier(mesh->MESH_COMM);

        // write the checkpoint for each wind turbine
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->turbineControlled)
            {
                // first scatter the necessary information on the master node of this TRB_COMM
                PetscReal rtrOmega     = 0.0;
                PetscReal rtrOmegaFilt = 0.0;
                PetscReal genTorque    = 0.0;
                PetscReal genPwr       = 0.0;
                PetscReal collPitch    = 0.0;
                PetscReal errPID       = 0.0;
                PetscReal intErrPID    = 0.0;
                PetscReal azimuth      = 0.0;

                if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                {
                    MPI_Reduce(&(wt->rtrOmega), &rtrOmega,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->rtrOmegaFilt),&rtrOmegaFilt,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->errPID),      &errPID,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->intErrPID),   &intErrPID,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }

                    if((*farm->turbineModels[t]) == "ALM")
                    {
                        MPI_Reduce(&(wt->alm.azimuth), &azimuth,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }
                }

                // now only the master of TRB_COMM writes
                if(rank == wt->writerRank)
                {
                    FILE *f;
                    char fileName[80];
                    sprintf(fileName, "%s/%s", path.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
                    }
                    else
                    {
                        // turbine level properties
                        fprintf(f, "turbineLevelProperties\n{\n");
                        fprintf(f, "    rtrDir             (%lf %lf %lf)\n", wt->rtrDir.x, wt->rtrDir.y, wt->rtrDir.z);
                        fprintf(f, "    rtrAxis            (%lf %lf %lf)\n", wt->rtrAxis.x, wt->rtrAxis.y, wt->rtrAxis.z);

                        if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                        {
                            fprintf(f, "    omega_hat          (%lf %lf %lf)\n", wt->omega_hat.x, wt->omega_hat.y, wt->omega_hat.z);
                            fprintf(f, "    rtrOmega           %lf\n",         rtrOmega / wt->nProcsTrb);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "    rtrOmegaFilt       %lf\n",         rtrOmegaFilt / wt->nProcsTrb);
                                fprintf(f, "    genTorque          %lf\n",         genTorque    / wt->nProcsTrb);
                                fprintf(f, "    genPwr             %lf\n",         genPwr       / wt->nProcsTrb);
                            }
                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "    collPitch          %lf\n",         collPitch    / wt->nProcsTrb);
                                fprintf(f, "    errPID             %lf\n",         errPID       / wt->nProcsTrb);
                                fprintf(f, "    intErrPID          %lf\n",         intErrPID    / wt->nProcsTrb);
                            }

                            if((*farm->turbineModels[t]) == "ALM")
                            {
                                fprintf(f, "    azimuth          %lf\n",           azimuth    / wt->nProcsTrb);
                            }

                        }
                        if(wt->yawControllerType != "none")
                        {
                            fprintf(f, "    yawAngle           %lf\n",         wt->yawAngle);
                            fprintf(f, "    flowAngle          %lf\n",         wt->flowAngle);
                            fprintf(f, "    yawError           %lf\n",         wt->yawError);
                        }

                        fprintf(f, "}\n\n");
                    }

                    fclose(f);
                }
            }
        }

        // ensure all files are closed
        MPI_Barrier(mesh->MESH_COMM);

        // remove old checkpoint files except last one after all files are written (safe)
        if(!rank && mesh->access->io->purgeWrite) remove_subdirs_except_keep_n(mesh->MESH_COMM, turbineFolderName.c_str(), timeName, mesh->access->io->purgeWrite-1);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windTurbinesReadCheckpoint(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    clock_      *clock = farm->access->clock;

    word        timeName;
    PetscInt         t, p;

    PetscReal      initialYaw, initialYawError;

    timeName = "./fields/" + mesh->meshName + "/turbines/" + getTimeName(clock);

    // check if this folder exists
    if(dir_exist(timeName.c_str()))
    {
        PetscPrintf(mesh->MESH_COMM, "   found checkpoint data %s, restoring configuration...\n",timeName.c_str());

        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // each processor reads data and stores
            char dictName[80];
            sprintf(dictName, "%s/%s", timeName.c_str(), (*farm->turbineIds[t]).c_str());

            // save initial yaw before overwriting it
            initialYaw = std::atan2(-1.0 * wt->rtrDir.y, -1.0 * wt->rtrDir.x) * wt->rad2deg;

            readSubDictVector(dictName, "turbineLevelProperties", "rtrDir", &(wt->rtrDir));
            readSubDictVector(dictName, "turbineLevelProperties", "rtrAxis", &(wt->rtrAxis));

            if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
            {
                readSubDictVector(dictName, "turbineLevelProperties", "omega_hat", &(wt->omega_hat));
                readSubDictDouble(dictName, "turbineLevelProperties", "rtrOmega", &(wt->rtrOmega));

                if(wt->genControllerType != "none")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "rtrOmegaFilt", &(wt->rtrOmegaFilt));
                    readSubDictDouble(dictName, "turbineLevelProperties", "genTorque", &(wt->genTorque));
                    readSubDictDouble(dictName, "turbineLevelProperties", "genPwr", &(wt->genPwr));
                }
                if(wt->pitchControllerType != "none")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "collPitch", &(wt->collPitch));
                    readSubDictDouble(dictName, "turbineLevelProperties", "errPID", &(wt->errPID));
                    readSubDictDouble(dictName, "turbineLevelProperties", "intErrPID", &(wt->intErrPID));
                }
                if((*farm->turbineModels[t]) == "ALM")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "azimuth", &(wt->alm.azimuth));
                }
            }
            if(wt->yawControllerType != "none")
            {
                readSubDictDouble(dictName, "turbineLevelProperties", "yawAngle", &(wt->yawAngle));
                readSubDictDouble(dictName, "turbineLevelProperties", "flowAngle", &(wt->flowAngle));
                readSubDictDouble(dictName, "turbineLevelProperties", "yawError", &(wt->yawError));
            }
        }

        // yaw turbines
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // test if this turbine has yaw control
            if(wt->yawControllerType != "none")
            {
                PetscReal yawError = wt->yawAngle - initialYaw;
                PetscReal rotationSign;

                if(yawError > 0.0)       rotationSign =  1.0;
                else if(yawError < -0.0) rotationSign = -1.0;
                else                     rotationSign =  0.0;

                // compute rotation angle in radiants
                PetscReal angle = rotationSign * yawError * wt->deg2rad;

                if(angle != 0.0)
                {
                    // actuator disk model
                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->adm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->adm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->adm.points[p], angle);
                            mSum(wt->adm.points[p], wt->rotCenter);
                        }
                    }
                    // uniform actuator disk model
                    else if((*farm->turbineModels[t]) == "uniformADM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->uadm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->uadm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->uadm.points[p], angle);
                            mSum(wt->uadm.points[p], wt->rotCenter);
                        }
                    }
                    // actuator line model
                    else if((*farm->turbineModels[t]) == "ALM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->alm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->alm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->alm.points[p], angle);
                            mSum(wt->alm.points[p], wt->rotCenter);
                        }
                    }
                    else if((*farm->turbineModels[t]) == "AFM")
                    {
                        // nothing to rotate
                    }

                    // rotate up-sampling points
                    {
                        upSampling *upPoints = farm->wt[t]->upPoints;

                        // rotate center
                        mSub(upPoints->center, wt->rotCenter);
                        mRot(wt->twrDir, upPoints->center, angle);
                        mSum(upPoints->center, wt->rotCenter);

                        // number of points in the sample mesh
                        PetscInt npts_t = upPoints->nPoints;

                        // loop over the sample mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(upPoints->points[p], wt->rotCenter);
                            mRot(wt->twrDir, upPoints->points[p], angle);
                            mSum(upPoints->points[p], wt->rotCenter);
                        }
                    }
                }
            }

            // rotate blades
            if((*farm->turbineModels[t]) == "ALM")
            {
                PetscInt updateAzimuth = 0;
                rotateBlades(wt, wt->alm.azimuth * wt->deg2rad, updateAzimuth);
            }
        }
    }
    else
    {
        PetscPrintf(mesh->MESH_COMM, "   no checkpoint data found, setting reference configuration...\n");
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readActuatorModelParameters(const word actuatorModel, windTurbine *wt, const word meshName)
{
    if(actuatorModel == "ADM")
    {
        readADM(wt, meshName);
    }
    else if(actuatorModel == "uniformADM")
    {
        readUADM(wt, meshName);
    }
    else if(actuatorModel == "ALM")
    {
        readALM(wt, meshName);  
    }
    else if(actuatorModel == "AFM")
    {
        readAFM(wt, meshName);    
    }
    else
    {
        char error[512];
        sprintf(error, "unknown wind turbine model for turbine %s\n", wt->id.c_str());
        fatalErrorInFunction("readActuatorModelParameters",  error);
    }
    
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readADM(windTurbine *wt, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(ADM), &(wt->adm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AD parameters
    readDictInt(descrFile.c_str(), "nRadPts", &(wt->adm.nRadial));
    readDictInt(descrFile.c_str(), "nAziPts", &(wt->adm.nAzimuth));
    readDictDouble(descrFile.c_str(), "Uref",    &(wt->adm.Uref));

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->adm.dbg));

    // set total numer of points in the mesh
    wt->adm.nPoints = wt->adm.nRadial * wt->adm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->adm.rtrThrust = 0.0;
    wt->adm.rtrTorque = 0.0;
    wt->adm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->adm.rtrAvgMagU = 0.0;

    // allocate memory for the ADM parameters
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.points));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.dr));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.chord));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.twist));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.solidity));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscInt*), &(wt->adm.foilIds));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal*), &(wt->adm.iw));

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->adm.nPoints*sizeof(PetscInt), &(wt->adm.thisPtControlled));
    PetscMalloc(wt->adm.nPoints*sizeof(cellIds),&(wt->adm.closestCells));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.Cd));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.Cl));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.alpha));
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.U));
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.B));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.axialF));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.tangtF));

    return(0);
}

PetscErrorCode meshADM(windTurbine *wt)
{
    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);

    // set rtrAxis and omega_hat turbine param. here
    {
        // set the rotor rotation unit vector
        if(wt->rotDir == "cw")
        {
            wt->omega_hat = nScale(-1.0, wt->rtrAxis);
        }
        else if(wt->rotDir == "ccw")
        {
            wt->omega_hat = nSet(wt->rtrAxis);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown rotationDir, avilable options are 'cw' or 'ccw'\n");
            fatalErrorInFunction("meshADM",  error);
        }
    }

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->adm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->adm.nAzimuth;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr;

    for(PetscInt ri=0; ri<wt->adm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->adm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

		// compute radius
		PetscReal  rMag = nMag(rvec);

        // add the initial hub radius to rvec (not aligned if precone != 0)
        mSum(rvec, rHub);

		// add hub radius (as if all was aligned)
		rMag += nMag(rHub);

        // interpolate blade propertes (only depend on r)
        PetscReal  weights[2];
        PetscInt   labels[2];

        findInterpolationWeights(weights, labels, wt->blade.radius, wt->blade.size, rMag);

        for(PetscInt ai=0; ai<wt->adm.nAzimuth; ai++)
        {
            // set interpolation variables
            PetscReal w1 = weights[0]; PetscInt l1 = labels[0];
            PetscReal w2 = weights[1]; PetscInt l2 = labels[1];

            // allocate memory for the 2 closest foil ids
            PetscMalloc(2*sizeof(PetscInt), &(wt->adm.foilIds[pi]));

            // allocate memory for the 2 interpolation weights
            PetscMalloc(2*sizeof(PetscReal), &(wt->adm.iw[pi]));

            // set airfoil chord
            wt->adm.chord[pi] = w1*wt->blade.chord[l1] + w2*wt->blade.chord[l2];

            // set airfoil twist
            wt->adm.twist[pi] = w1*wt->blade.twist[l1] + w2*wt->blade.twist[l2];

            // set airfoil ids
            wt->adm.foilIds[pi][0] = wt->blade.foilIds[l1];
            wt->adm.foilIds[pi][1] = wt->blade.foilIds[l2];

            // set interpolation weights
            wt->adm.iw[pi][0] = w1;
            wt->adm.iw[pi][1] = w2;

            // set the rotor solidity (only depends on radius)
            wt->adm.solidity[pi] = (PetscReal)wt->nBlades / (PetscReal)wt->adm.nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(wt->rtrAxis, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, wt->rotCenter);

            // set the point value
            mSet(wt->adm.points[pi], point);

            // set the dr value (uniform for now)
            wt->adm.dr[pi] = dr;

            pi++;
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readUADM(windTurbine *wt, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(UADM), &(wt->uadm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AD parameters
    readDictInt(descrFile.c_str(),    "nRadPts",    &(wt->uadm.nRadial));
    readDictInt(descrFile.c_str(),    "nAziPts",    &(wt->uadm.nAzimuth));
    readDictDouble(descrFile.c_str(), "Ct",         &(wt->uadm.Ct));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->uadm.sampleType));
    readDictDouble(descrFile.c_str(), "Uref",       &(wt->uadm.Uref));

    // check sample type
    if
    (
        wt->uadm.sampleType != "rotorUpstream" &&  // samples velocity on 2.5D upstream disk (needs Ct)
        wt->uadm.sampleType != "givenVelocity" &&  // uses input velocity (needs Ct, good for isolated turbine)
        wt->uadm.sampleType != "rotorDisk"         // samples velocity on rotor disk (needs CtPrime)
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are givenVelocity, rotorUpstream or rotorDisk");
        fatalErrorInFunction("readUADM",  error);
    }

    // compute axial induction factor
    if(wt->uadm.sampleType != "rotorDisk")
    {
        // check that Ct does not make induction complex or negative
        if(wt->uadm.Ct <= 0.0 || wt->uadm.Ct >= 1.0)
        {
            char error[512];
            sprintf(error, "provided thrust coefficient ouside of bounds ([0 1] excluded). Change or switch to sampleType = rotorDisk");
            fatalErrorInFunction("readUADM",  error);
        }

        // Ct
        wt->uadm.axiInd = (1.0 - sqrt(1 - wt->uadm.Ct)) / 2.0;
    }
    else
    {
        // CtPrime
        wt->uadm.axiInd = 0.0;
    }

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->uadm.dbg));

    // set total numer of points in the mesh
    wt->uadm.nPoints = wt->uadm.nRadial * wt->uadm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->uadm.rtrThrust = 0.0;
    wt->uadm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->uadm.rtrAvgMagU = 0.0;

    // build the AD mesh

    // allocate memory for the ADM parameters
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.points));
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscReal), &(wt->uadm.dA));

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscInt), &(wt->uadm.thisPtControlled));
    PetscMalloc(wt->uadm.nPoints*sizeof(cellIds),&(wt->uadm.closestCells));
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.U));
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.B));
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscReal), &(wt->uadm.axialF));

    return(0);
}

PetscErrorCode meshUADM(windTurbine *wt)
{
    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->uadm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->uadm.nAzimuth;

    // area of the circular crown from hub to tip radius
    PetscReal bladeSweptAreaNoHub = M_PI * (wt->rTip*wt->rTip - wt->rHub*wt->rHub);

    // hub area
    PetscReal hubFrontArea = M_PI * wt->rHub * wt->rHub;

    // hub area distribution coeff
    PetscReal hCoeff = hubFrontArea / bladeSweptAreaNoHub;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr, crownArea, r = wt->rHub;

    for(PetscInt ri=0; ri<wt->uadm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->uadm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // compute r+dr crown area
        crownArea = M_PI * ((r+dr)*(r+dr) - r*r);

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

        // add the initial hub radius
        mSum(rvec, rHub);

        for(PetscInt ai=0; ai<wt->uadm.nAzimuth; ai++)
        {
            // set cell area
            wt->uadm.dA[pi] = (1.0 + hCoeff) * crownArea / wt->uadm.nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(wt->rtrAxis, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, wt->rotCenter);

            // set the point value
            mSet(wt->uadm.points[pi], point);

            pi++;
        }

        r += dr;
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode readALM(windTurbine *wt, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(ALM), &(wt->alm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AL parameters
    readDictWord(descrFile.c_str(),   "projection", &(wt->alm.projectionType));
    readDictInt(descrFile.c_str(),    "nRadPts", &(wt->alm.nRadial));
    readDictDouble(descrFile.c_str(), "Uref",    &(wt->alm.Uref));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->alm.sampleType));

    // check sample type
    if
    (
        wt->alm.sampleType != "rotorDisk" &&
        wt->alm.sampleType != "integral"
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are rotorDisk or integral");
        fatalErrorInFunction("readALM",  error);
    }

    // check projection type
    if(wt->alm.projectionType!="isotropic" && wt->alm.projectionType!="anisotropic")
    {
        char error[512];
        sprintf(error, "unknown ALM projection, available possibilities are:\n    1. isotropic\n    2. anisotropic\n");
        fatalErrorInFunction("readALM",  error);
    }

    //read projection radius
    if(wt->alm.projectionType=="anisotropic")
    {
        if(wt->useOpenFAST)
        {
            char error[512];
            sprintf(error, "anisotropic projection not compatible with OpenFAST coupling at turbine %s", wt->id.c_str());
            fatalErrorInFunction("readALM",  error);
        }
        readDictDouble(descrFile.c_str(), "epsilonFactor_x",     &(wt->eps_x));
        readDictDouble(descrFile.c_str(), "epsilonFactor_y",     &(wt->eps_y));
        readDictDouble(descrFile.c_str(), "epsilonFactor_z",     &(wt->eps_z));
    }
    else
    {
        readDictDouble(descrFile.c_str(), "epsilon",     &(wt->eps));
    }

    wt->alm.nAzimuth = wt->nBlades;

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->alm.dbg));

    // set number of blade force points in OpenFAST 
    if(wt->useOpenFAST)
    {
        #if USE_OPENFAST
        // here we force the number of force points to be defined by TOSCA
        wt->nBladeForcePtsOF = wt->alm.nRadial;
        #endif
    }

    // set total numer of points in the mesh
    wt->alm.nPoints = wt->alm.nRadial * wt->alm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->alm.rtrThrust = 0.0;
    wt->alm.rtrTorque = 0.0;
    wt->alm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->alm.rtrAvgMagU = 0.0;

    // set initial azimuth
    wt->alm.azimuth = 0.0;

    // allocate memory for the ALM parameters
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.points));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.dr));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.chord));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.twist));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.solidity));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscInt*), &(wt->alm.foilIds));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal*), &(wt->alm.iw));

    if(wt->alm.projectionType=="anisotropic")
    {
        PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.thick));
    }

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->alm.nPoints*sizeof(PetscInt), &(wt->alm.thisPtControlled));
    PetscMalloc(wt->alm.nPoints*sizeof(cellIds),&(wt->alm.closestCells));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.Cd));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.Cl));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.alpha));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.U));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.gWind));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.B));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.axialF));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.tangtF));

    return(0);
}

PetscErrorCode meshALM(windTurbine *wt)
{
    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);

    // set omega_hat turbine param. here
    {
        // set the rotor rotation unit vector
        if(wt->rotDir == "cw")
        {
            wt->omega_hat = nScale(-1.0, wt->rtrAxis);
        }
        else if(wt->rotDir == "ccw")
        {
            wt->omega_hat = nSet(wt->rtrAxis);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown rotationDir, avilable options are 'cw' or 'ccw'\n");
            fatalErrorInFunction("meshALM",  error);
        }
    }

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->alm.nAzimuth;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr;

    for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->alm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

        // compute radius
        PetscReal  rMag = nMag(rvec);

        // add the initial hub radius to rvec (not aligned if precone != 0)
        mSum(rvec, rHub);

        // add hub radius (as if all was aligned)
        rMag += nMag(rHub);

        // interpolate blade propertes (only depend on r)
        PetscReal  weights[2];
        PetscInt   labels[2];

        findInterpolationWeights(weights, labels, wt->blade.radius, wt->blade.size, rMag);

        for(PetscInt ai=0; ai<wt->alm.nAzimuth; ai++)
        {
            // set interpolation variables
            PetscReal w1 = weights[0]; PetscInt l1 = labels[0];
            PetscReal w2 = weights[1]; PetscInt l2 = labels[1];

            // allocate memory for the 2 closest foil ids
            PetscMalloc(2*sizeof(PetscInt), &(wt->alm.foilIds[pi]));

            // allocate memory for the 2 interpolation weights
            PetscMalloc(2*sizeof(PetscReal), &(wt->alm.iw[pi]));

            // set airfoil chord
            wt->alm.chord[pi] = w1*wt->blade.chord[l1] + w2*wt->blade.chord[l2];

            // set airfoil twist
            wt->alm.twist[pi] = w1*wt->blade.twist[l1] + w2*wt->blade.twist[l2];

            // set airfoil thickness
            if(wt->alm.projectionType=="anisotropic")
            {
                wt->alm.thick[pi] = w1*wt->blade.thick[l1] + w2*wt->blade.thick[l2];
            }

            // set airfoil ids
            wt->alm.foilIds[pi][0] = wt->blade.foilIds[l1];
            wt->alm.foilIds[pi][1] = wt->blade.foilIds[l2];

            // set interpolation weights
            wt->alm.iw[pi][0] = w1;
            wt->alm.iw[pi][1] = w2;

            // set the rotor solidity (it is one for the ALM)
            wt->alm.solidity[pi] = 1.0;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(wt->rtrAxis, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, wt->rotCenter);

            // set the point value
            mSet(wt->alm.points[pi], point);

            // set the dr value (uniform for now)
            wt->alm.dr[pi] = dr;

            pi++;
        }
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode readAFM(windTurbine *wt, const word meshName)
{
    // allocate memory for the AFM
    PetscMalloc(sizeof(AFM), &(wt->afm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    readDictDouble(descrFile.c_str(), "Ct",         &(wt->afm.Ct));
    readDictDouble(descrFile.c_str(), "Uref",       &(wt->afm.Uref));
    readDictWord(descrFile.c_str(),   "projection", &(wt->afm.projectionType));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->afm.sampleType));

    // check projection type
    if(wt->afm.projectionType!="gaussexp" && wt->afm.projectionType!="anisotropic")
    {
        char error[512];
        sprintf(error, "unknown AFM projection, available possibilities are:\n    1. gaussexp\n    2. anisotropic\n");
        fatalErrorInFunction("readAFM",  error);
    }

    //read projection radius
    if(wt->afm.projectionType=="anisotropic")
    {
        readDictDouble(descrFile.c_str(), "epsilon_x",     &(wt->eps_x));
        readDictDouble(descrFile.c_str(), "epsilon_y",     &(wt->eps_y));
        readDictDouble(descrFile.c_str(), "epsilon_z",     &(wt->eps_z));
    }
    else if(wt->afm.projectionType=="gaussexp")
    {
        readDictDouble(descrFile.c_str(), "gaussexp_x",     &(wt->eps_x));
        readDictDouble(descrFile.c_str(), "gaussexp_r",     &(wt->r12));
        readDictDouble(descrFile.c_str(), "gaussexp_f",     &(wt->flat));

        // compute normalization factor
        wt->I = - 2.0*wt->eps_x*pow(M_PI,3.0/2.0)*wt->flat*wt->flat*polyLog2(-std::exp(wt->r12/wt->flat));
    }

    // check sample type
    if
    (
        wt->afm.sampleType != "rotorDisk" &&
        wt->afm.sampleType != "momentumTheory" &&
        wt->afm.sampleType != "integral"
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are momentumTheory, rotorDisk or integral");
        fatalErrorInFunction("readAFM",  error);
    }

    // check input
    if(wt->afm.sampleType == "momentumTheory" && wt->afm.Ct > 0.999)
    {
        char error[512];
        sprintf(error, "Ct coefficient must be less than one if sampling is of momentumTheory type");
        fatalErrorInFunction("readAFM",  error);
    }

    // read sampling frequency
    PetscReal windSpeedFilterPrd;
    readDictDouble(descrFile.c_str(), "windSpeedFilterPrd", &windSpeedFilterPrd);
    wt->afm.rtrUFilterFreq = 1.0 / windSpeedFilterPrd;

    wt->afm.searchDone = 0;

    return(0);
}

PetscErrorCode meshAFM(windTurbine *wt)
{
    wt->afm.point  = wt->rotCenter;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initSamplePoints(windTurbine *wt, const word meshName)
{
    // allocate memory for this turbine sample points struct
    wt->upPoints = (upSampling*)malloc(sizeof(upSampling));
    upSampling* upPoints = wt->upPoints;

    // discretization is not very fine since only area-weighted
    // average must be performed
    PetscInt nAzimuth = 12;
    PetscInt nRadial  = 12;

    upPoints->nPoints = nAzimuth * nRadial;

    // allocate objects memory
    PetscMalloc(upPoints->nPoints*sizeof(Cmpnts),    &(upPoints->points));
    PetscMalloc(upPoints->nPoints*sizeof(cellIds),   &(upPoints->closestCells));
    PetscMalloc(upPoints->nPoints*sizeof(PetscInt),  &(upPoints->thisPtControlled));
    PetscMalloc(upPoints->nPoints*sizeof(PetscReal), &(upPoints->dA));

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    // (do not rotate with uptilt)
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);

    // set the upstream sample distance
    Cmpnts upDist = nScale(2.5*2.0*wt->rTip, xr_hat);

    // set the sampling rotor center
    Cmpnts center = wt->rotCenter;
    mSum(center, upDist);
    upPoints->center = center;

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / nAzimuth;

    // area of the circular crown from hub to tip radius
    PetscReal bladeSweptAreaNoHub = M_PI * (wt->rTip*wt->rTip - wt->rHub*wt->rHub);

    // hub area
    PetscReal hubFrontArea = M_PI * wt->rHub * wt->rHub;

    // hub area distribution coeff
    PetscReal hCoeff = hubFrontArea / bladeSweptAreaNoHub;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr, crownArea, r = wt->rHub;

    for(PetscInt ri=0; ri<nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==nRadial-1) dr = drval / 2;
        else dr = drval;

        // compute r+dr crown area
        crownArea = M_PI * ((r+dr)*(r+dr) - r*r);

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, yr_hat);

        // add the initial hub radius
        mSum(rvec, rHub);

        for(PetscInt ai=0; ai<nAzimuth; ai++)
        {
            // set cell area
            upPoints->dA[pi] = (1.0 + hCoeff) * crownArea / nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(xr_hat, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, center);

            // set the point value
            mSet(upPoints->points[pi], point);

            pi++;
        }

        r += dr;
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode initTwrModel(windTurbine *wt, Cmpnts &base)
{
    // set tower thrust to zero
    wt->twr.twrThrust = 0.0;

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->twr.prjNSigma = sqrt(std::log(1.0/0.001));

    // allocate memory
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.points));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscReal), &(wt->twr.dA));
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.U));
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.B));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscReal), &(wt->twr.tangF));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscInt), &(wt->twr.thisPtControlled));
    PetscMalloc(wt->twr.nPoints*sizeof(cellIds), &(wt->twr.closestCells));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);

    // linear element distance
    PetscReal dhval  = nMag(tower) / (wt->twr.nPoints - 1);

    // rastremation coefficient
    PetscReal rastCoeff = (wt->twr.rTop - wt->twr.rBase) / wt->hTwr;

    // varying variables
    PetscReal dh, h = 0;

    for(PetscInt pi=0; pi<wt->twr.nPoints; pi++)
    {
        if(pi==0 || pi==(wt->twr.nPoints-1)) dh = dhval / 2;
        else dh = dhval;

        // bottom and top radiuses of this tower segment
        PetscReal r_lo = rastCoeff*h + wt->twr.rBase;
        PetscReal r_hi = rastCoeff*(h+dh) + wt->twr.rBase;

        // compute area using trapezoidal formula for rastremation
        wt->twr.dA[pi] = 0.5 * (r_lo + r_hi) * dh;

        // set current tower point
        Cmpnts twr_point = nScale(h, wt->twrDir);
                           mSum(twr_point, base);

        // save tower point
        wt->twr.points[pi] = nSet(twr_point);

        h += dh;
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initNacModel(windTurbine *wt, Cmpnts &base)
{
    // set tower thrust to zero
    wt->nac.nacThrust = 0.0;

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->nac.prjNSigma = sqrt(std::log(1.0/0.001));

    // set nacelle point
    wt->nac.point     = nSum(nScale(wt->hTwr, wt->twrDir), base);

    // set nacelle area
    wt->nac.A         = M_PI*wt->rHub*wt->rHub;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeFarmADMesh(farm_ *farm)
{
    mesh_     *mesh = farm->access->mesh;
    clock_    *clock = farm->access->clock;

    PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    if(!rank)
    {
        // compute number of points and cells
        PetscInt npts = 0;
        PetscInt ncll = 0;
        for(PetscInt t=0; t<farm->size; t++)
        {
            PetscInt npts_t = 0, nrc_t = 0, nac_t = 0;

            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // number of points in this turbine's AD mesh
                npts_t = farm->wt[t]->adm.nPoints;

                // number of radial cells in this turbine's AD mesh
                nrc_t  = farm->wt[t]->adm.nRadial - 1;

                // number of azimuthal cells in this turbine's AD mesh
                nac_t  = farm->wt[t]->adm.nAzimuth;
            }
            // uniform actuator disk model
            if((*farm->turbineModels[t]) == "uniformADM")
            {
                // number of points in this turbine's AD mesh
                npts_t = farm->wt[t]->uadm.nPoints;

                // number of radial cells in this turbine's AD mesh
                nrc_t  = farm->wt[t]->uadm.nRadial - 1;

                // number of azimuthal cells in this turbine's AD mesh
                nac_t  = farm->wt[t]->uadm.nAzimuth;
            }

            npts += npts_t;
            ncll += nrc_t*nac_t;
        }

        if(npts>0)
        {
            word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

            char fileName[256];
            sprintf(fileName, "%s/ADMesh_%ld.inp", turbineFolderTimeName.c_str(), farm->writeNumber);

            PetscInt width = -20;

            FILE *f = fopen(fileName, "w");

            // header
            PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");

            // number of points, number of cells
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

            // write coordinates
            npts = 1;

            for(PetscInt t=0; t<farm->size; t++)
            {
                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    // number of points in this turbine's AD mesh
                    PetscInt npts_t = farm->wt[t]->adm.nPoints;

                    for(PetscInt pi=0; pi<npts_t; pi++)
                    {
                        Cmpnts point_i;

                        point_i.x = farm->wt[t]->adm.points[pi].x;
                        point_i.y = farm->wt[t]->adm.points[pi].y;
                        point_i.z = farm->wt[t]->adm.points[pi].z;

                        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, npts, width, point_i.x, width, point_i.y, width, point_i.z);

                        npts++;
                    }
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    // number of points in this turbine's AD mesh
                    PetscInt npts_t = farm->wt[t]->uadm.nPoints;

                    for(PetscInt pi=0; pi<npts_t; pi++)
                    {
                        Cmpnts point_i;

                        point_i.x = farm->wt[t]->uadm.points[pi].x;
                        point_i.y = farm->wt[t]->uadm.points[pi].y;
                        point_i.z = farm->wt[t]->uadm.points[pi].z;

                        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, npts, width, point_i.x, width, point_i.y, width, point_i.z);

                        npts++;
                    }
                }
            }

            // write connectivity

            // storage for the labels
            PetscInt cellPtLabels[ncll][4];

            npts = 0;
            ncll = 0;

            // build connectivity
            for(PetscInt t=0; t<farm->size; t++)
            {
                PetscInt nRadCells_t, nAziCells_t;
                PetscInt nRadPts_t, nAziPts_t;

                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    // radial and azimuthal cells
                    nRadCells_t = farm->wt[t]->adm.nRadial - 1;
                    nAziCells_t = farm->wt[t]->adm.nAzimuth;

                    // radial and azimuthal points
                    nRadPts_t = farm->wt[t]->adm.nRadial;
                    nAziPts_t = farm->wt[t]->adm.nAzimuth;
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    // radial and azimuthal cells
                    nRadCells_t = farm->wt[t]->uadm.nRadial - 1;
                    nAziCells_t = farm->wt[t]->uadm.nAzimuth;

                    // radial and azimuthal points
                    nRadPts_t = farm->wt[t]->uadm.nRadial;
                    nAziPts_t = farm->wt[t]->uadm.nAzimuth;
                }

                // number of points per turbine
                PetscInt npt = nRadPts_t*nAziPts_t;

                for(PetscInt ri=0; ri<nRadPts_t; ri++)
                {
                    for(PetscInt ai=0; ai<nAziPts_t; ai++)
                    {
                        if(ri < nRadCells_t)
                        {
                            // all azi cells except last
                            if(ai < (nAziCells_t-1))
                            {
                                cellPtLabels[ncll][0] = (t*npt) + (ri * nAziPts_t + ai + 1);
                                cellPtLabels[ncll][1] = (t*npt) + (ri * nAziPts_t + ai + 2);
                                cellPtLabels[ncll][2] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 2);
                                cellPtLabels[ncll][3] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 1);

                                ncll++;
                            }
                            // last cell: pts 0,3 must be connected to pts 0,3 of first azi cell
                            //            so pt 1,2 become 0,3 of first azi cell
                            else if (ai == (nAziCells_t-1))
                            {
                                cellPtLabels[ncll][0] = (t*npt) + (ri * nAziPts_t + ai + 1);
                                cellPtLabels[ncll][1] = (t*npt) + (ri * nAziPts_t + 1);
                                cellPtLabels[ncll][2] = (t*npt) + ((ri + 1) * nAziPts_t + 1);
                                cellPtLabels[ncll][3] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 1);

                                ncll++;
                            }
                        }
                        npts++;
                    }
                }
            }

            width = -5;

            // write
            for(PetscInt ci=0; ci<ncll; ci++)
            {
                PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d %*d %*d\n", width, ci+1, width, 0, width, "quad", width, cellPtLabels[ci][0], width, cellPtLabels[ci][1], width, cellPtLabels[ci][2], width, cellPtLabels[ci][3]);
            }

            fclose(f);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeFarmTwrMesh(farm_ *farm)
{
    mesh_  *mesh  = farm->access->mesh;
    clock_ *clock = farm->access->clock;

    PetscMPIInt           rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    if(!rank)
    {
        // compute number of points and edges
        // each tower is only one line
        PetscInt npts = 2*farm->size;
        PetscInt ncll = farm->size;

        word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

        char fileName[256];
        sprintf(fileName, "%s/twrMesh.inp", turbineFolderTimeName.c_str());

        PetscInt width = -20;

        FILE *f = fopen(fileName, "w");

        // header
        PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
        PetscFPrintf(mesh->MESH_COMM, f, "#\n");
        PetscFPrintf(mesh->MESH_COMM, f, "#\n");

        // number of points, number of cells
        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

        for(PetscInt t=0; t<farm->size; t++)
        {
            // set tower top point
            Cmpnts towerPt  = nScale(farm->wt[t]->hTwr, farm->wt[t]->twrDir);
                              mSum(towerPt, farm->base[t]);

            // set tower base point
            Cmpnts basePt   = nSet(farm->base[t]);

            // set point labels
            PetscInt id_base = 2*(t + 1) - 1;
            PetscInt id_twr  = 2*(t + 1);

            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, id_base, width, basePt.x, width, basePt.y, width, basePt.z);
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, id_twr, width, towerPt.x, width, towerPt.y, width, towerPt.z);
        }

        // write connectivity
        width = -5;

        // build connectivity
        for(PetscInt t=0; t<farm->size; t++)
        {
            // set point labels
            PetscInt id_base = 2*(t + 1) - 1;
            PetscInt id_twr  = 2*(t + 1);

            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d\n", width, t+1, width, 0, width, "line", width, id_base, width, id_twr);
        }

        fclose(f);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeFarmALMesh(farm_ *farm)
{
    mesh_      *mesh = farm->access->mesh;
    clock_    *clock = farm->access->clock;

    PetscMPIInt rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    PetscInt    point_id, cell_id;

    if(!rank)
    {
        // compute number of points and cells
        PetscInt npts = 0;
        PetscInt ncll = 0;
        for(PetscInt t=0; t<farm->size; t++)
        {
            if((*farm->turbineModels[t]) == "ALM")
            {
                npts += 6;
                ncll += 3;
            }
        }

        if(npts>0)
        {
            word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

            char fileName[256];
            sprintf(fileName, "%s/ALMesh_%ld.inp", turbineFolderTimeName.c_str(), farm->writeNumber);

            PetscInt width = -20;

            FILE *f = fopen(fileName, "w");

            // header
            PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");

            // number of points, number of cells
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

            // write coordinates
            point_id = 0;

            for(PetscInt t=0; t<farm->size; t++)
            {
                // define tip/root blade points in reference (vertical) configuration
                Cmpnts bldTip  = nScale(farm->wt[t]->rTip, farm->wt[t]->twrDir);
                Cmpnts bldRoot = nScale(farm->wt[t]->rHub, farm->wt[t]->twrDir);

                Cmpnts *points;
                PetscMalloc(6*sizeof(Cmpnts), &points);

                PetscReal angle = farm->wt[t]->alm.azimuth;

                for(PetscInt b=0; b<3; b++)
                {
                    Cmpnts bldTip_b  = nRot(farm->wt[t]->omega_hat, bldTip,  angle*farm->wt[t]->deg2rad);
                    Cmpnts bldRoot_b = nRot(farm->wt[t]->omega_hat, bldRoot, angle*farm->wt[t]->deg2rad);

                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, point_id, width, bldTip_b.x, width, bldTip_b.y, width, bldTip_b.z);
                    point_id++;
                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, point_id, width, bldRoot_b.x, width, bldRoot_b.y, width, bldRoot_b.z);
                    point_id++;

                    angle = angle + 120;
                }
            }

            // write connectivity
            width = -5;

            point_id = 0;
            cell_id  = 1;

            // build connectivity
            for(PetscInt t=0; t<farm->size; t++)
            {
                for(PetscInt b=0; b<3; b++)
                {
                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d\n", width, cell_id, width, 0, width, "line", width, point_id, width, point_id+1);
                    cell_id++;
                    point_id = point_id + 2;
                }
            }

            fclose(f);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readFarmProperties(farm_ *farm)
{
    // The dictionary to read is windFarmProperties,
    // located inside ./turbines/ folder. It has a
    // first subdictionary for the write settings:
    //
    // writeSettings
    // {
    //     timeStart           10.0     // time at which the output starts to be written
    //     intervalType        timeStep // adjustableTime
    //     timeInterval        1        // output every second if adjustableTime or iteration if timeStep
    // }
    //
    // Then it contans the wind farm name and the wind farm
    // specification (onebyone or grid), if onebyone each turbine is
    // specified as follows
    //
    // WT1
    // {
    //     turbineType         NREL5MW
    //     turbineModel        ADM
    //     base                (5.0191 0.0 -90.0)
    //     windFarmController  1
    // }
    //
    // If grid we still have to implement it.

    // get mesh pointer
    mesh_ *mesh = farm->access->mesh;

    word   windFarmPropertiesFile = "./turbines/" + mesh->meshName + "/windFarmProperties";

    // allocate memory for wind farm body force
    VecDuplicate(mesh->lCent, &(farm->lsourceFarmCat));
    VecDuplicate(mesh->Cent,  &(farm->sourceFarmCont));

    // read the write settings
    readSubDictDouble(windFarmPropertiesFile.c_str(),"writeSettings","timeStart", &(farm->timeStart));
    readSubDictWord  (windFarmPropertiesFile.c_str(),"writeSettings","intervalType", &(farm->intervalType));
    readSubDictDouble(windFarmPropertiesFile.c_str(),"writeSettings","timeInterval", &(farm->timeInterval));

    // initialize write number
    farm->writeNumber = 0;

    // read debug flag
    readDictInt(windFarmPropertiesFile.c_str(), "debug", &(farm->dbg));

    // read wind farm name
    readDictWord(windFarmPropertiesFile.c_str(), "windFarmName", &(farm->name));

    // read array specification
    word arraySpec;
    readDictWord(windFarmPropertiesFile.c_str(), "arraySpecification", &arraySpec);

#if USE_OPENFAST

    // use OpenFAST turbine-specific flag
    readDictInt(windFarmPropertiesFile.c_str(),    "nFastSubSteps",&(farm->nFastSubSteps));
    
#endif

    // read the wind turbines one by one until end of file
    if(arraySpec=="onebyone")
    {
        readTurbineArray(farm);
    }
    else if(arraySpec=="grid")
    {
        char error[512];
        sprintf(error, "array specification type %s not yet implemented\n", arraySpec.c_str());
        fatalErrorInFunction("readFarmProperties",  error);
    }
    else
    {
       char error[512];
        sprintf(error, "unknown array specification type %s\n", arraySpec.c_str());
        fatalErrorInFunction("readFarmProperties",  error);
    }

    // Allocate memory for the wind turbines, this is an array of pointers.
    // Each pointer will point to a turbine. Advantages: the size is the
    // one of a pointer for the allocation and the objects can be not
    // linear in memory.

    farm->wt = new windTurbine* [farm->size];

    // go through each wind turbine and read its description file
    for(PetscInt t=0; t<farm->size; t++)
    {
        // dynamic memory allocation (make the t-th pointer point to the struct)
        farm->wt[t] = new windTurbine;

        // set this turbine type and ID
        farm->wt[t]->type = *(farm->turbineTypes[t]);
        farm->wt[t]->id   = *(farm->turbineIds[t]);

        // allocate memory for the rotor axis (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->rtrAxis));

        // allocate memory for the rotor angular velocity unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->omega_hat));

        // allocate memory for the tower direction unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->twrDir));

        // allocate memory for the rotor direction unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->rtrDir));

        // this turbine description file (can be shared if type is the same)
        word descrFile = "./turbines/" + mesh->meshName + "/" + (*farm->turbineTypes[t]);

        // read turbine properties
        readTurbineProperties(farm->wt[t], descrFile.c_str(), mesh->meshName, (*farm->turbineModels[t]), farm->base[t]);

        // read wind farm control table for this wind turbine
        if(farm->farmControlActive[t])
        {
            readWindFarmControlTable(farm->wt[t]);
        }

        // set turbine index 
        farm->wt[t]->index = t;
    }

    // set CFL checking flag to zero
    farm->checkCFL = 0;

    printFarmProperties(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTurbineArray(farm_ *farm)
{
    // get mesh pointer
    mesh_ *mesh = farm->access->mesh;

    // define the local variables
    std::vector<std::string>   turbineTypes;
    std::vector<std::string>   turbineIds;
    std::vector<std::string>   turbineModels;
    std::vector<Cmpnts>        base;
    std::vector<PetscInt>      windFarmController;

    std::string turbineTypes_i,
                turbineIds_i,
                turbineModels_i;
    Cmpnts      base_i;
    PetscInt    controller_i;

    PetscInt    nturbines = 0;

    // pointer for strtod and strtol
    char        *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];
    std::string token;

    // dictionary name
    std::string dictName = "./turbines/" + mesh->meshName + "/windFarmProperties" ;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName.c_str());
        fatalErrorInFunction("readTurbineArray",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "turbineArray",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read all turbines until end of file (should find "}" before eof)
                    while(!indata.eof())
                    {
                        // from here start to read the wind turbines

                        // read WT ID
                        indata >> word;
                        turbineIds_i = word;

                        // check if have hit the end of the list
                        if(trim(turbineIds_i)=="}")
                        {
                            if(nturbines>0)
                            {
                                // close the file
                                indata.close();

                                // allocate memory and store the variables
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineTypes));
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineIds));
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineModels));
                                PetscMalloc(nturbines*sizeof(Cmpnts),  &(farm->base));
                                PetscMalloc(nturbines*sizeof(PetscInt),&(farm->farmControlActive));

                                // base locations
                                for(PetscInt p=0; p<nturbines; p++)
                                {
                                    PetscMalloc(sizeof(std::string), &(farm->turbineTypes[p]));
                                    PetscMalloc(sizeof(std::string), &(farm->turbineIds[p]));
                                    PetscMalloc(sizeof(std::string), &(farm->turbineModels[p]));

                                    // assign the pointers to the singly created variables in memory
                                    farm->turbineTypes[p]  = new std::string(turbineTypes[p]);
                                    farm->turbineIds[p]    = new std::string(turbineIds[p]);
                                    farm->turbineModels[p] = new std::string(turbineModels[p]);

                                    //PetscPrintf(mesh->MESH_COMM,"%s\n",(*farm->turbineTypes[p]).c_str());

                                    farm->base[p].x = base[p].x;
                                    farm->base[p].y = base[p].y;
                                    farm->base[p].z = base[p].z;

                                    farm->farmControlActive[p] = windFarmController[p];
                                }

                                // wind farm size
                                farm->size = nturbines;

                                // clear the local variables
                                std::vector<std::string>   ().swap(turbineTypes);
                                std::vector<std::string>   ().swap(turbineIds);
                                std::vector<std::string>   ().swap(turbineModels);
                                std::vector<Cmpnts>        ().swap(base);

                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected at least one turbine in subdictionary turbineArray of %s dictionary\n", dictName.c_str());
                                fatalErrorInFunction("readTurbineArray",  error);
                            }
                        }

                        // check if have hit another list (this was not closed with "}")
                        if(trim(turbineIds_i)=="}")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token in subdictionary turbineArray of %s dictionary\n", dictName.c_str());
                            fatalErrorInFunction("readTurbineArray",  error);
                        }

                        // store the turbine ID
                        turbineIds.push_back(turbineIds_i);

                        // parameters are enclosed by '()', so read the first
                        indata >> word;
                        token = word;
                        if(trim(token)=="(")
                        {
                            // read turbineType keyword
                            indata >> word;
                            token = word;
                            if(trim(token)=="turbineType")
                            {
                                // read the turbine type
                                indata >> word;
                                turbineTypes_i = word;
                                turbineTypes.push_back(turbineTypes_i);

                                // read turbineModel keyword
                                indata >> word;
                                token = word;
                                if(trim(token)=="turbineModel")
                                {
                                    // read the turbine model
                                    indata >> word;
                                    turbineModels_i = word;
                                    turbineModels.push_back(turbineModels_i);

                                    // read baseLocation keyword
                                    indata >> word;
                                    token = word;
                                    if(trim(token)=="baseLocation")
                                    {
                                        // read the turbine base vector

                                        // start reading vector ------------------------------------------------------------------------------------------------------------------------------------------------

                                        // the vector is in (x y z) format, so
                                        // 1. read the first parenthesis
                                        // 2. read the 3 doubles
                                        // 3. look for the closing parethesis

                                        // read the first component (contains "(" character)
                                        indata >> word;

                                        std::string first(word);
                                        if (first.find ("(") != std::string::npos)
                                        {
                                           // remove "("" character from the first component
                                           PetscInt l1 = first.size();
                                           for(PetscInt i=0;i<l1;i++)
                                           {
                                               // save the first component
                                               word[i] = word[i+1];
                                           }

                                           base_i.x = std::strtod(word, &eptr);

                                           // check if the first component is a PetscReal, throw error otherwise
                                           std::string cmp1(word);

                                           if(isNumber(cmp1))
                                           {
                                               if (cmp1.find ('.') == std::string::npos)
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                              char error[512];
                                               sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }

                                           // read the second component
                                           indata >> word;
                                           base_i.y = std::strtod(word, &eptr);

                                           // check if the second component is a PetscReal, throw error otherwise
                                           std::string cmp2(word);

                                           if(isNumber(cmp2))
                                           {
                                               if (cmp2.find ('.') == std::string::npos)
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                              char error[512];
                                               sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }

                                           // read the third component (contains ")" character)
                                           indata >> word;

                                           std::string last(word);
                                           if (last.find (")") != std::string::npos)
                                           {
                                               // remove ") character from the last component and store
                                               base_i.z = std::strtod(word, &eptr);

                                               // check if the first component is a PetscReal, throw error otherwise
                                               std::string cmp3(word);

                                               if(isNumber(cmp3))
                                               {
                                                   if (cmp3.find ('.') == std::string::npos)
                                                   {
                                                      char error[512];
                                                       sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                       fatalErrorInFunction("readTurbineArray",  error);
                                                   }
                                               }
                                               else
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }

                                               // save the base location
                                               base.push_back(base_i);

                                               // read windFarmController keyword
                                               indata >> word;
                                               token = word;
                                               if(trim(token)=="windFarmController")
                                               {
                                                   // read wether controller is active or not
                                                   indata >> word;
                                                   controller_i = (PetscInt)std::strtol(word, &eptr, 10);
                                                   windFarmController.push_back(controller_i);

                                                   // increase turbine counter
                                                   nturbines++;

                                                   // read the closing ')'
                                                   indata >> word;
                                                   token = word;
                                                   if(trim(token)!=")")
                                                   {
                                                       char error[512];
                                                       sprintf(error, "expected <)>  at end of subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                       fatalErrorInFunction("readTurbineArray",  error);
                                                   }
                                               }
                                               else
                                               {
                                                   char error[512];
                                                   sprintf(error, "expected windFarmController keyword after %s in subdictionary %s of %s dictionary, found '%s'\n", last.c_str(), turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                               char error[512];
                                               sprintf(error, "expected <(>  after vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }
                                       }
                                       else
                                       {
                                           char error[512];
                                           sprintf(error, "expected <(>  after keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                           fatalErrorInFunction("readTurbineArray",  error);
                                       }

                                       // End of reading vector ------------------------------------------------------------------------------------------------------------------------------------------------
                                    }
                                    else
                                    {
                                        char error[512];
                                        sprintf(error, "expected turbineModel keyword after %s in subdictionary %s of %s dictionary, found '%s'\n", turbineModels_i.c_str(), turbineIds_i.c_str(), dictName.c_str(), word);
                                        fatalErrorInFunction("readTurbineArray",  error);
                                    }
                                }
                                else
                                {
                                   char error[512];
                                    sprintf(error, "expected turbineModel keyword after %s in sub-dictionary %s, found '%s'\n", turbineTypes_i.c_str(), turbineIds_i.c_str(), word);
                                    fatalErrorInFunction("readTurbineArray",  error);
                                }
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected turbineType keyword after ( in sub-dictionary %s, found '%s'\n", turbineIds_i.c_str(), word);
                                fatalErrorInFunction("readTurbineArray",  error);
                            }
                        }
                        else
                        {
                           char error[512];
                            sprintf(error, "expected <(> token after keyword %s in dictionary %s, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                            fatalErrorInFunction("readTurbineArray",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of turbineArray subdictionary in %s dictionary\n", dictName.c_str());
                    fatalErrorInFunction("readTurbineArray",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword turbineArray in dictionary %s, found '%s'\n", dictName.c_str(), word);
                    fatalErrorInFunction("readTurbineArray",  error);
                }

            }
        }

       char error[512];
        sprintf(error, "could not find subdictionary turbineArray in dictionary %s\n", dictName.c_str());
        fatalErrorInFunction("readTurbineArray",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTurbineProperties(windTurbine *wt, const char *dictName, const word meshName, const word modelName, Cmpnts &base)
{
    // local vectors (to be normalised)
    Cmpnts twrDirVec, rtrDirVec;

    // read main parameters
    readDictDouble(dictName, "rTip",        &(wt->rTip));
    readDictDouble(dictName, "rHub",        &(wt->rHub));
    readDictDouble(dictName, "hTower",      &(wt->hTwr));
    readDictDouble(dictName, "overHang",    &(wt->ovrHang));
    readDictDouble(dictName, "precone",     &(wt->precone));
    readDictVector(dictName, "towerDir",    &(twrDirVec));
    readDictVector(dictName, "rotorDir",    &(rtrDirVec));
    readDictDouble(dictName, "upTilt",      &(wt->upTilt));

    if(modelName == "ADM" || modelName == "uniformADM")
    {
        readDictDouble(dictName, "epsilon",     &(wt->eps));
    }
    else if(modelName == "AFM")
    {
        // do nothing - read later in initAFM based on projection type
    }
    else if(modelName == "ALM")
    {
        // do nothing - read later in initALM based on projection type
    }

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->prjNSigma = sqrt(std::log(1.0/0.001));

    // tower model flag
    readDictInt(dictName,    "includeTower",&(wt->includeTwr));

    // nacelle model flag
    readDictInt(dictName,    "includeNacelle",&(wt->includeNacelle));

    // set OpenFAST usage flag to zero by default
    wt->useOpenFAST = 0;

#if USE_OPENFAST

    // use OpenFAST turbine-specific flag
    readDictInt(dictName,    "useOpenFAST",&(wt->useOpenFAST));

    // check that model is compatible with OpenFAST
    if(wt->useOpenFAST)
    {
        if(modelName != "ALM")
        {
            char error[512];
            sprintf(error, "turbine model %s not compatible with useOpenFAST option in turbine type %s\n", modelName.c_str(), wt->type.c_str());
            fatalErrorInFunction("readTurbineProperties",  error);
        }
    }
    
    // initialize number of blade and tower points in OpenFAST (they are set when reading actuator models) 
    wt->nBladeForcePtsOF  = 0;
    wt->nTwrForcePtsOF    = 0;
    wt->nBladeVelPtsOF    = 0;
    wt->nTwrVelPtsOF      = 0;
    
#endif

    // debug switch
    readDictInt(dictName,    "debug",     &(wt->dbg));

    // store normalized vectors
    wt->twrDir = nUnit(twrDirVec);
    wt->rtrDir = nUnit(rtrDirVec);

    // conversion factors
    wt->deg2rad = M_PI / 180.0;
    wt->rad2deg = 180.0 / M_PI;
    wt->rpm2RadSec = 2 * M_PI / 60.0;

    // initialize rotor omega to zero (for CFL check access)
    wt->rtrOmega = 0.0;

    // set rotor center point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);
    wt->rotCenter = center;

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

    // set rotor axis (up-tilted, from nacell back to cone)
    wt->rtrAxis = xr_hat;

    // read parameters for AD/AL models
    if(modelName == "ADM" || modelName == "ALM")
    {
        PetscReal initOmega;

        readDictInt(dictName,    "nBlades",     &(wt->nBlades));
        readDictWord(dictName,   "rotationDir", &(wt->rotDir));
        readDictDouble(dictName, "initialOmega",&(initOmega));
        wt->rtrOmega = initOmega * wt->rpm2RadSec;

        // read torque controller
        readDictWord(dictName,   "genControllerType", &(wt->genControllerType));

        // read pitch controller type
        readDictWord(dictName,   "pitchControllerType", &(wt->pitchControllerType));

        // read controllers input parameters
        if(wt->genControllerType   != "none") readGenControllerParameters(wt,   wt->genControllerType.c_str(), meshName.c_str());
        if(wt->pitchControllerType != "none") readPitchControllerParameters(wt, wt->pitchControllerType.c_str(), meshName.c_str());

        // read airfoil types used in this turbine
        readAirfoilProperties(wt, dictName);

        // Allocate memory for the foil info, this is an array of pointers.
        // Each pointer will point to an airfoil. Advantages: the size is the
        // one of a pointer for the allocation and the objects can be not
        // linear in memory.
        wt->foils = (foilInfo**)malloc(wt->nFoils*sizeof(foilInfo*));

        // read the property file for each airfoil
        for(PetscInt f=0; f<wt->nFoils; f++)
        {
            // dynamic memory allocation (make the f-th pointer point to the struct)
            wt->foils[f] = (foilInfo*)malloc(sizeof(foilInfo));

            // set the name of the airfoil
            PetscMalloc(sizeof(word), &(wt->foils[f]->name));
            wt->foils[f]->name = *(wt->foilNames[f]);

            // set the path to the airfoil data file
            word name2af = "./turbines/" + meshName + "/airfoils/" + *(wt->foilNames[f]);

            // set the reset of the variables by reading the table
            readAirfoilTable(wt->foils[f], name2af.c_str());
        }

        PetscInt readThickness = 0;

        if(modelName == "ALM")
        {
            word projectionType;
            readDictWord(dictName, "projection", &projectionType);

            if(projectionType=="anisotropic")
            {
                readThickness = 1;
            }
        }

        // read blade properties
        readBladeProperties(wt, dictName, readThickness);

        // check that the max airfoil label in the blade properties actually
        // matches the number if provided airfoils
        PetscInt max_id = 0;

        // find max id
        for(PetscInt i=0; i<wt->blade.size; i++)
        {
            if(wt->blade.foilIds[i] > max_id)
            {
                max_id = wt->blade.foilIds[i];
            }
        }

        if(max_id+1>wt->nFoils)
        {
            char error[512];
            sprintf(error, "requested more airfoils than the number provided in turbine type %s (airfoils and bladeData mismatch)\n",wt->type.c_str());
            fatalErrorInFunction("readTurbineProperties",  error);
        }
    }

    // read yaw controller parameters (all models: ADM/ALM/UADM/AFM)
    readDictWord(dictName,   "yawControllerType", &(wt->yawControllerType));

    if(wt->yawControllerType   != "none") readYawControllerParameters(wt,   wt->yawControllerType.c_str(), meshName.c_str());

    if(wt->includeTwr)
    {
        readTowerProperties(wt, dictName);
    }

    if(wt->includeNacelle)
    {
        readNacelleProperties(wt, dictName);
    }

    // set yaw changed parameter to 1 at initialization
    // this makes sure all models do the first cell to point search
    // ALM always does the search as blades are rotating
    wt->yawChanged = 1;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readAirfoilProperties(windTurbine *wt, const char *dictName)
{
    // First a list of names is read from the dictionary passed as argument.
    // This must be defined as follows:
    //
    // airfoils
    // {
    //     name1
    //     ...
    //     nameN
    // }
    //
    // Then the each airfoil look up tables are read using the readAirfoilTable
    // function.

    // define the local variables
    std::vector<std::string>   foilNames;

    std::string foilNames_i;
    PetscInt         nfoils = 0;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readAirfoilProperties",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "airfoils",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the airfoil name
                        indata >> word;

                        // check that did not hit end of subdictionary "}"
                        std::string token(word);
                        if(trim(token)=="}")
                        {
                            if(nfoils>0)
                            {
                                // close file
                                indata.close();

                                // allocate memory
                                PetscMalloc(nfoils*sizeof(std::string*), &(wt->foilNames));

                                // store the data
                                for(PetscInt f=0; f<nfoils; f++)
                                {
                                    PetscMalloc(sizeof(std::string), &(wt->foilNames[f]));

                                    wt->foilNames[f] = new std::string(foilNames[f]);
                                }

                                wt->nFoils = nfoils;

                                // clear memory
                                std::vector<std::string> ().swap(foilNames);

                                // exit
                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 1 entry as reading airfoils subdictionary in %s dictionary\n", dictName);
                                fatalErrorInFunction("readAirfoilProperties",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of airfoils subdictionary in %s dictionary\n", dictName);
                            fatalErrorInFunction("readAirfoilProperties",  error);
                        }

                        // store the airfoil name
                        foilNames_i = word;
                        foilNames.push_back(foilNames_i);

                        nfoils++;

                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of airfoils subdictionary in %s dictionary\n", dictName);
                    fatalErrorInFunction("readAirfoilProperties",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword airfoils in dictionary %s, found '%s'\n", dictName, word);
                    fatalErrorInFunction("readAirfoilProperties",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword airfoils in dictionary %s\n", dictName);
        fatalErrorInFunction("readAirfoilProperties",  error);

    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readAirfoilTable(foilInfo *af, const char *tableName)
{
    // Airfoil tables must be provided as follows:
    //
    // table
    // {
    //     (alpha_deg1  cl1  cd1)
    //               ...
    //     (alpha_degN  clN  cdN)
    // }
    //
    // Suggestion: in order for the AD an AL models to work
    // properly, provide the curves between -180 to +180 degs.
    // The number of points is arbitrary, minimum is 3.

    // define the local variables
    std::vector<PetscReal>   aoa;
    std::vector<PetscReal>   cl;
    std::vector<PetscReal>   cd;

    PetscReal aoa_i, cl_i, cd_i;
    PetscInt    nPoints = 0;

    // pointer for the strtod function
    char   *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(tableName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s file\n", tableName);
        fatalErrorInFunction("readAirfoilTable",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "table",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the first component, contains "(" character
                        indata >> word;

                        // put the word in a string for performing checks
                        std::string token(word);

                        std::string first(word);
                        if (first.find ("(") != std::string::npos)
                        {
                            // remove "( character from the first component
                            PetscInt l1 = first.size();
                            for(PetscInt i=0;i<l1;i++)
                            {
                                word[i] = word[i+1];
                            }

                            // store the first variable (aoa)
                            aoa_i = std::strtod(word, &eptr);
                            aoa.push_back(aoa_i);

                            // read the second component (cl)
                            indata >> word;

                            // store second component
                            cl_i = std::strtod(word, &eptr);
                            cl.push_back(cl_i);

                            // read the third component (cd), contains ")" character
                            indata >> word;

                            std::string last(word);
                            if (last.find (")") != std::string::npos)
                            {
                                // remove ") character from the last component and store
                                cd_i = std::strtod(word, &eptr);
                                cd.push_back(cd_i);

                                // increament point counter
                                nPoints++;
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected <)>  at end of line as reading table in %s dictionary\n", tableName);
                                fatalErrorInFunction("readAirfoilTable",  error);
                            }
                        }
                        else if(trim(token)=="}")
                        {
                            if(nPoints >= 3)
                            {
                                // close file
                                indata.close();

                                // allocate memory
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->aoa));
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->cl));
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->cd));

                                // store the table
                                for(PetscInt i=0; i<nPoints; i++)
                                {
                                    af->aoa[i] = aoa[i];
                                    af->cl[i]  = cl[i];
                                    af->cd[i]  = cd[i];
                                }

                                // store number of points in the table
                                af->size = nPoints;

                                // clear memory
                                std::vector<PetscReal>   ().swap(aoa);
                                std::vector<PetscReal>   ().swap(cl);
                                std::vector<PetscReal>   ().swap(cd);

                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 3 data points as reading table in %s dictionary\n", tableName);
                                fatalErrorInFunction("readAirfoilTable",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of table in %s dictionary\n", tableName);
                            fatalErrorInFunction("readAirfoilTable",  error);
                        }
                        // we are at a new line and neither '(' nor '}' were found: throws error
                        else
                        {
                             char error[512];
                              sprintf(error, "expected either <(>  or <}> at new line as reading table in %s dictionary\n", tableName);
                              fatalErrorInFunction("readAirfoilTable",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of table in %s dictionary\n", tableName);
                    fatalErrorInFunction("readAirfoilTable",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword table in dictionary %s, found '%s'\n", tableName, word);
                    fatalErrorInFunction("readAirfoilTable",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword 'table' in file %s\n", tableName);
        fatalErrorInFunction("readAirfoilTable",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readBladeProperties(windTurbine *wt, const char *dictName, const PetscInt readThickness)
{
    // The dictionary must be located inside a
    // file named as the turbine type and located
    // inside ./turbines/, e.g.: ./turbines/NREL5MW.dat.
    // It must be defined as follows:
    //
    // bladeData
    // {
    //     (radius1 chord1 twist1 foilIds1)
    //     (radius2 chord2 twist2 foilIds2)
    //                   ...
    //     (radiusN chordN twistN foilIdsN)
    // }
    //
    // two is the minimum number of data points.
    // no checks for doubles or integers are performed,
    // only the format is checked.
    // For anisotropic gaussian ALM projection also blade thickness should be provided in the last column

    // define the local variables
    std::vector<PetscReal> radius;
    std::vector<PetscReal> chord;
    std::vector<PetscReal> twist;
    std::vector<PetscReal> thick;
    std::vector<PetscInt> foilIds;

    PetscReal   radius_i, chord_i, twist_i, thick_i;
    PetscInt    foilIds_i;
    PetscInt    nlines = 0;

    // pointer for strtod and strtol
    char   *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readBladeProperties",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "bladeData",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);
                std::string token2;

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the first component, contains "(" character
                        indata >> word;

                        // put the word in a string for performing checks
                        std::string token3(word);

                        std::string first(word);
                        if (first.find ("(") != std::string::npos)
                        {
                            // remove "( character from the first component
                            PetscInt l1 = first.size();
                            for(PetscInt i=0;i<l1;i++)
                            {
                                word[i] = word[i+1];
                            }

                            // store the first variable (radius)
                            radius_i = std::strtod(word, &eptr);
                            radius.push_back(radius_i);

                            // read the second component (chord)
                            indata >> word;

                            // store second component
                            chord_i = std::strtod(word, &eptr);
                            chord.push_back(chord_i);

                            // read the third component (twist)
                            indata >> word;

                            // store third component
                            twist_i = std::strtod(word, &eptr);
                            twist.push_back(twist_i);

                            if(readThickness)
                            {
                                // read the thickness
                                indata >> word;

                                // store thickness component
                                thick_i = std::strtod(word, &eptr);
                                thick.push_back(thick_i);
                            }

                            // read the fourth component (foilIds), contains ")" character
                            indata >> word;

                            std::string last(word);
                            if (last.find (")") != std::string::npos)
                            {
                                // remove ") character from the last component and store
                                foilIds_i = std::strtol(word, &eptr, 10);
                                foilIds.push_back(foilIds_i);

                                // increament line counter
                                nlines++;
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected <)>  at end of line as reading readData table in %s dictionary\n", dictName);
                                fatalErrorInFunction("readBladeProperties",  error);
                            }

                        }
                        // look for the terminating "}" if found: close the file, store the data and exit
                        // if not found: may be another line
                        else if(trim(token3)=="}")
                        {
                            if(nlines >= 2)
                            {
                                // close file
                                indata.close();

                                // allocate memory for the blade properties in the current wind turbine
                                PetscMalloc(sizeof(bladeAeroInfo), &(wt->blade));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.radius));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.chord));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.twist));
                                PetscMalloc(nlines*sizeof(PetscInt), &(wt->blade.foilIds));

                                if(readThickness) PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.thick));

                                // store the blade properties
                                wt->blade.size = nlines;

                                for(PetscInt p=0; p<wt->blade.size; p++)
                                {
                                    wt->blade.radius[p]  = radius[p];
                                    wt->blade.chord[p]   = chord[p];
                                    wt->blade.twist[p]   = twist[p];
                                    wt->blade.foilIds[p] = foilIds[p];

                                    if(readThickness) wt->blade.thick[p]   = thick[p];
                                }

                                // clean the local variables
                                std::vector<PetscReal> ().swap(radius);
                                std::vector<PetscReal> ().swap(chord);
                                std::vector<PetscReal> ().swap(twist);
                                std::vector<PetscReal> ().swap(thick);
                                std::vector<PetscInt>    ().swap(foilIds);

                                // exit
                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 2 data points as reading readData table in %s dictionary\n", dictName);
                                fatalErrorInFunction("readBladeProperties",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token3)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of readData table in %s dictionary\n", dictName);
                            fatalErrorInFunction("readBladeProperties",  error);
                        }
                        // we are at a new line and neither '(' nor '}' were found: throws error
                        else
                        {
                             char error[512];
                              sprintf(error, "expected either <(>  or <}> at new line as reading readData table in %s dictionary\n", dictName);
                              fatalErrorInFunction("readBladeProperties",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of readData table in %s dictionary\n", dictName);
                    fatalErrorInFunction("readBladeProperties",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword bladeData in dictionary %s, found '%s'\n", dictName, word);
                    fatalErrorInFunction("readBladeProperties",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword bladeData in dictionary %s\n", dictName);
        fatalErrorInFunction("readBladeProperties",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTowerProperties(windTurbine *wt, const char *dictName)
{
    // allocate memory for the tower model
    PetscMalloc(sizeof(towerModel), &(wt->twr));

    readSubDictDouble(dictName, "towerData", "Cd",      &(wt->twr.Cd));
    readSubDictDouble(dictName, "towerData", "epsilon", &(wt->twr.eps));
    readSubDictInt(dictName,    "towerData", "nLinPts", &(wt->twr.nPoints));
    readSubDictDouble(dictName, "towerData", "rBase",   &(wt->twr.rBase));
    readSubDictDouble(dictName, "towerData", "rTop",    &(wt->twr.rTop));

    if(wt->useOpenFAST)
    {
        #if USE_OPENFAST
        // here we force the number of force points to be defined by TOSCA
        wt->nTwrForcePtsOF = wt->twr.nPoints;
        #endif 
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readNacelleProperties(windTurbine *wt, const char *dictName)
{
    // allocate memory for the tower model
    PetscMalloc(sizeof(nacelleModel), &(wt->nac));

    readSubDictDouble(dictName, "nacelleData", "Cd",      &(wt->nac.Cd));
    readSubDictDouble(dictName, "nacelleData", "epsilon", &(wt->nac.eps));

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readGenControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // rotor dynamics parameters
    readDictDouble(path2dict, "genInertia",         &(wt->genInertia));
    readDictDouble(path2dict, "hubInertia",         &(wt->hubInertia));
    readDictDouble(path2dict, "bldIntertia",        &(wt->bldInertia));
    readDictDouble(path2dict, "gbxRatioG2R",        &(wt->gbxRatioG2R));
    readDictDouble(path2dict, "gbxEfficiency",      &(wt->gbxEff));

    // initial conditions
    readDictDouble(path2dict, "gbxRatioG2R",        &(wt->gbxRatioG2R));
    readDictDouble(path2dict, "initialGenTorque",   &(wt->genTorque));

    // generator electrical efficiency
    readDictDouble(path2dict,  "genEff",            &(wt->genEff));

    // controller parameters
    readSubDictDouble(path2dict, "genTqControllerParameters","genSpeedFilterFreq", &(wt->rtrSpdFilterFreq));
    readSubDictDouble(path2dict, "genTqControllerParameters","cutInGenSpeed",      &(wt->cutInGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","cutInGenTorque",     &(wt->cutInGenTq));
    readSubDictDouble(path2dict, "genTqControllerParameters","regTwoStartGenSpeed",&(wt->regTwoStartGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","regTwoEndGenSpeed",  &(wt->regTwoEndGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","ratedGenTorque",     &(wt->ratedGenTq));
    readSubDictDouble(path2dict, "genTqControllerParameters","controllerPGain",    &(wt->omegaKP));

    // limits
    readSubDictInt(path2dict,    "genTqControllerParameters","torqueRateLimiter",  &(wt->tqRateLimiter));
    readSubDictInt(path2dict,    "genTqControllerParameters","rtrSpeedLimiter",    &(wt->rtrSpdLimiter));
    readSubDictDouble(path2dict, "genTqControllerParameters","torqueMaxRate",      &(wt->tqMaxRate));
    readSubDictDouble(path2dict, "genTqControllerParameters","ratedRotorSpeed",    &(wt->ratedRotorSpd));

    // convert rpm to rad/s
    wt->ratedRotorSpd     *= wt->rpm2RadSec;
    wt->regTwoEndGenSpd   *= wt->rpm2RadSec;
    wt->regTwoStartGenSpd *= wt->rpm2RadSec;
    wt->cutInGenSpd       *= wt->rpm2RadSec;
    wt->omegaKP           /= (wt->rpm2RadSec*wt->rpm2RadSec);
    wt->rtrSpdFilterFreq  /= 60;

    // total inertia on HS LS shafts
    wt->driveTrainInertia
    =
    wt->nBlades*wt->bldInertia +
    wt->hubInertia +
    wt->gbxRatioG2R*wt->gbxRatioG2R*wt->genInertia;

    // generator speed in rpm
    wt->genOmega  = wt->rtrOmega * wt->gbxRatioG2R;

    // fitered rotor omega
    wt->rtrOmegaFilt = wt->rtrOmega;

    // generator power
    wt->genPwr    = 0.0;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readPitchControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // allocate memory for collective pitch
    PetscMalloc(wt->nBlades*sizeof(PetscReal), &(wt->pitch));

    // PID controller paramters
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerPGain",   &(wt->pitchKP));
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerIGain",   &(wt->pitchKI));
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerDGain",   &(wt->pitchKD));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchS2R",          &(wt->pitchS2R));

    // limits
    readSubDictInt(path2dict,    "pitchControllerParameters","pitchRateLimiter",  &(wt->pitchRateLimiter));
    readSubDictInt(path2dict,    "pitchControllerParameters","pitchAngleLimiter", &(wt->pitchAngleLimiter));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMaxRate",      &(wt->pitchMaxRate));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMin",          &(wt->pitchMin));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMax",          &(wt->pitchMax));

    // transform all in radiants
    wt->pitchMaxRate *= wt->deg2rad;
    wt->pitchMax     *= wt->deg2rad;
    wt->pitchMin     *= wt->deg2rad;
    wt->pitchKP      *= wt->deg2rad;
    wt->pitchKI      *= wt->deg2rad;
    wt->pitchKD      *= wt->deg2rad;
    wt->pitchS2R     *= wt->deg2rad;

    // check and validate
    if(fabs(wt->pitchKI) < 1e-5)
    {
       char error[512];
        sprintf(error, "integral gain 'controlledIGain' cannot be zero\n");
        fatalErrorInFunction("readPitchControllerParameters",  error);
    }
    if(fabs(wt->pitchS2R) < 1e-5)
    {
       char error[512];
        sprintf(error, "point at which the power sensitivity to pitch is twice the one at rated conditions 'pitchS2R' cannot be zero\n");
        fatalErrorInFunction("readPitchControllerParameters",  error);
    }

    // set collective and individual blade pitch
    wt->collPitch = 0.0;
    for(PetscInt bld_i=0; bld_i<wt->nBlades; bld_i++)
    {
        wt->pitch[bld_i] = 0.0;
    }

    // set initial speed error to zero
    wt->errPID = 0.0;

    // set integrated error such that if the pitch is correct at the startTime
    // it will not change at the next time step
    PetscReal G   = 1.0 / (1.0 + wt->collPitch / wt->pitchS2R);
    wt->intErrPID = wt->collPitch / (G * wt->pitchKI);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readYawControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // controller paramters
    readSubDictWord(path2dict,   "yawControllerParameters","sampleType",  &(wt->yawSamplingType));
    readSubDictDouble(path2dict, "yawControllerParameters","avgWindow",   &(wt->yawAverageWindow));
    readSubDictDouble(path2dict, "yawControllerParameters","yawMin",      &(wt->yawMin));
    readSubDictDouble(path2dict, "yawControllerParameters","yawMax",      &(wt->yawMax));
    readSubDictDouble(path2dict, "yawControllerParameters","yawSpeed",    &(wt->yawSpeed));
    readSubDictDouble(path2dict, "yawControllerParameters","allowedError",&(wt->yawAllowedError));

    // set initial flow angle
    readSubDictDouble(path2dict, "yawControllerParameters","initialFlowAngle", &(wt->flowAngle));

    // make sure yawSpeed is positive
    wt->yawSpeed = fabs(wt->yawSpeed);

    // check that sampling type is known
    if
    (
        wt->yawSamplingType != "hubUpDistance" &&
        wt->yawSamplingType != "anemometer"
    )
    {
       char error[512];
        sprintf(error, "unknown yaw sampleType %s. Known types are hubUpDistance and anemometer", wt->yawSamplingType.c_str());
        fatalErrorInFunction("readYawControllerParameters",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readWindFarmControlTable(windTurbine *wt)
{
    word tableName = "turbines/control/" + wt->id;

    std::vector<std::vector<PetscReal>> table;

    // open file stream
    std::ifstream indata;
    indata.open(tableName);

    // word by word read
    char word[256];

    // buffer for read data
    PetscReal buffer;

    // time counter
    PetscInt ntimes;

    if(!indata)
    {
        char error[512];
        sprintf(error, "cannot open file %s\n", tableName.c_str());
        fatalErrorInFunction("readWindFarmControlTable",  error);
    }
    else
    {
        std::string tmpStr;

        // read lines and get number of saved times
        ntimes = 0;
        for (int t = 0; std::getline(indata, tmpStr); t++)
        {
            if (!tmpStr.empty())
            {
                ntimes++;
            }
        }

        // first line is header
        ntimes--;

        // save the number of times
        wt->wfControlNData  = ntimes;
        wt->currentCloseIdx = 0;

        // go back on top of file
        indata.close();
        indata.open(tableName);

        // skip header line
        std::getline(indata, tmpStr);

        // resize the source table
        table.resize(ntimes);

        for(PetscInt t=0; t<ntimes; t++)
        {
            // read along the line: time | value
            for(PetscInt i=0; i<2; i++)
            {
                table[t].resize(2);

                indata >> word;
                std::sscanf(word, "%lf", &buffer);

                table[t][i] = buffer;
            }

        }

        indata.close();
    }

    // now store the source  and free the temporary variable
    PetscMalloc(sizeof(PetscReal) * ntimes, &(wt->wfControlTimes));
    PetscMalloc(sizeof(PetscReal) * ntimes, &(wt->wfControlValues));

	// hard-coded difference between the sim. start and control action start
	PetscReal initialShift = 100000;

    for(PetscInt t=0; t<ntimes; t++)
    {
        wt->wfControlTimes[t]  = table[t][0] + initialShift;
        wt->wfControlValues[t] = table[t][1];

		//printf("time %.1f, ct: %.5f\n", wt->wfControlTimes[t], wt->wfControlValues[t]);
    }

    // clean the temporary variables
    for(PetscInt t=0; t<ntimes; t++)
    {
        std::vector<PetscReal> ().swap(table[t]);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode printFarmProperties(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    // check if at least one turbine has its debug flag activated
    PetscInt atLeastOneDebug = 0;

    for(PetscInt i=0;i<farm->size;i++)
    {
        if(farm->wt[i]->dbg)
        {
            atLeastOneDebug++;
        }
    }

    if(atLeastOneDebug)
    {
        PetscPrintf(mesh->MESH_COMM,"Wind Farm Properties for debugging\n\n");
    }

    for(PetscInt i=0;i<farm->size;i++)
    {
        if(farm->wt[i]->dbg)
        {
            // main turbine properties
            PetscPrintf(mesh->MESH_COMM,"Turbine %ld\n\n", i);
            PetscPrintf(mesh->MESH_COMM," ID      : %s\n", (*farm->turbineIds[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Type    : %s\n", (*farm->turbineTypes[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Model   : %s\n", (*farm->turbineModels[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Base    : (%lf, %lf, %lf)\n", farm->base[i].x, farm->base[i].y, farm->base[i].z);

            if((*farm->turbineModels[i]) != "uniformADM" && (*farm->turbineModels[i]) != "AFM")
            {
                PetscPrintf(mesh->MESH_COMM," Airfoils: ");

                for(PetscInt j=0;j<farm->wt[i]->nFoils; j++)
                {
                    PetscPrintf(mesh->MESH_COMM,"%s, ", (*farm->wt[i]->foilNames[j]).c_str());
                }

                PetscPrintf(mesh->MESH_COMM,"\n\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                PetscPrintf(mesh->MESH_COMM," Blade properties for this turbine:\n");

                // blade properties
                {
                    PetscInt n_bldpts = farm->wt[i]->blade.size;

                    PetscInt nch = -15;

                    PetscPrintf(mesh->MESH_COMM,"    | %*s | %*s | %*s| %*s|\n", nch, "radius [m]", nch, "chord [m]", nch, "twist [deg]", nch, "foil ID [-]");
                    for(PetscInt pti=0; pti<n_bldpts; pti++)
                    {
                        PetscReal radius   = farm->wt[i]->blade.radius[pti];
                        PetscReal chord    = farm->wt[i]->blade.chord[pti];
                        PetscReal twist    = farm->wt[i]->blade.twist[pti];
                        PetscInt    foilIds  = farm->wt[i]->blade.foilIds[pti];
                        word   foilName = (*farm->wt[i]->foilNames[foilIds]);

                        PetscPrintf(mesh->MESH_COMM,"    | %*.4f | %*.4f | %*.4f| %*s|\n", nch, radius, nch, chord, nch, twist, nch, foilName.c_str());
                    }
                }

                PetscPrintf(mesh->MESH_COMM,"\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                PetscPrintf(mesh->MESH_COMM," Breakdown of this turbine airfoil tables:\n");

                // airfoil properties
                for(PetscInt j=0;j<farm->wt[i]->nFoils; j++)
                {
                    // this airfoil table size
                    PetscInt n_afpts = farm->wt[i]->foils[j]->size;

                    PetscInt nch = -15;

                    PetscPrintf(mesh->MESH_COMM,"  \n%s:\n", (*farm->wt[i]->foilNames[j]).c_str());
                    PetscPrintf(mesh->MESH_COMM,"    | %*s | %*s | %*s|\n", nch, "alpha [deg]", nch, "Cl [-]", nch, "Cd [-]");
                    for(PetscInt pti=0; pti<n_afpts; pti++)
                    {
                        PetscReal alpha = farm->wt[i]->foils[j]->aoa[pti];
                        PetscReal cl    = farm->wt[i]->foils[j]->cl[pti];
                        PetscReal cd    = farm->wt[i]->foils[j]->cd[pti];
                        PetscPrintf(mesh->MESH_COMM,"    | %*.4f | %*.4f | %*.4f|\n", nch, alpha, nch, cl, nch, cd);
                    }
                }

                PetscPrintf(mesh->MESH_COMM,"\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");

                // controller properties
                if(farm->wt[i]->genControllerType != "none")
                {
                    PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                    PetscPrintf(mesh->MESH_COMM," Generator torque controller:\n");
                    PetscPrintf(mesh->MESH_COMM, "  genControllerType : %s\n", farm->wt[i]->genControllerType.c_str());
                    PetscPrintf(mesh->MESH_COMM, "  rtrOmega          : %lf\n", farm->wt[i]->rtrOmega);
                    PetscPrintf(mesh->MESH_COMM, "  genOmega          : %lf\n", farm->wt[i]->genOmega);
                    PetscPrintf(mesh->MESH_COMM, "  rtrOmegaFilt      : %lf\n", farm->wt[i]->rtrOmegaFilt);
                    PetscPrintf(mesh->MESH_COMM, "  rtrSpdFilterFreq  : %lf\n", farm->wt[i]->rtrSpdFilterFreq);
                    PetscPrintf(mesh->MESH_COMM, "  cutInGenSpd       : %lf\n", farm->wt[i]->cutInGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  cutInGenTq        : %lf\n", farm->wt[i]->cutInGenTq);
                    PetscPrintf(mesh->MESH_COMM, "  regTwoStartGenSpd : %lf\n", farm->wt[i]->regTwoStartGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  regTwoEndGenSpd   : %lf\n", farm->wt[i]->regTwoEndGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  ratedGenTq        : %lf\n", farm->wt[i]->ratedGenTq);
                    PetscPrintf(mesh->MESH_COMM, "  omegaKP           : %lf\n", farm->wt[i]->omegaKP);
                    PetscPrintf(mesh->MESH_COMM, "  genTorque         : %lf\n", farm->wt[i]->genTorque);
                    PetscPrintf(mesh->MESH_COMM, "  genPwr            : %lf\n", farm->wt[i]->genPwr);
                    PetscPrintf(mesh->MESH_COMM, "  tqRateLimiter     : %ld\n", farm->wt[i]->tqRateLimiter);
                    PetscPrintf(mesh->MESH_COMM, "  rtrSpdLimiter     : %ld\n", farm->wt[i]->rtrSpdLimiter);
                    PetscPrintf(mesh->MESH_COMM, "  tqMaxRate         : %lf\n", farm->wt[i]->tqMaxRate);
                    PetscPrintf(mesh->MESH_COMM, "  ratedRotorSpd     : %lf\n", farm->wt[i]->ratedRotorSpd);
                    PetscPrintf(mesh->MESH_COMM, "  driveTrainInertia : %lf\n", farm->wt[i]->driveTrainInertia);
                    PetscPrintf(mesh->MESH_COMM, "  genInertia        : %lf\n", farm->wt[i]->genInertia);
                    PetscPrintf(mesh->MESH_COMM, "  hubInertia        : %lf\n", farm->wt[i]->hubInertia);
                    PetscPrintf(mesh->MESH_COMM, "  bldInertia        : %lf\n", farm->wt[i]->bldInertia);
                    PetscPrintf(mesh->MESH_COMM, "  gbxRatioG2R       : %lf\n", farm->wt[i]->gbxRatioG2R);
                    PetscPrintf(mesh->MESH_COMM, "  gbxEff            : %lf\n", farm->wt[i]->gbxEff);

                }
            }

            PetscPrintf(mesh->MESH_COMM,"\n");
        }
    }

    MPI_Barrier(mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkTurbineMesh(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, p, c;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    Cmpnts           ***cent;
    PetscReal        ***aj;
    PetscReal        lMaxCell, gMaxCell, lMinCell, gMinCell;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        Cmpnts rtrAxis  = farm->wt[t]->rtrAxis;

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            lMaxCell = 0.0, gMaxCell = 0.0;
            lMinCell = 1e20, gMinCell = 1e20;

            // loop in sphere points
            for(c=0; c<wt->nControlled; c++)
            {
                // cell indices
                PetscInt i = wt->controlledCells[c].i,
                    j = wt->controlledCells[c].j,
                    k = wt->controlledCells[c].k;

                if(lMaxCell < 1/aj[k][j][i])
                {
                    lMaxCell = 1/aj[k][j][i];
                }

                if(lMinCell > 1/aj[k][j][i])
                {
                    lMinCell = 1/aj[k][j][i];
                }

            }

            MPI_Allreduce(&lMaxCell, &gMaxCell, 1, MPIU_REAL, MPIU_MAX, wt->TRB_COMM);
            MPI_Allreduce(&lMinCell, &gMinCell, 1, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

            gMaxCell = std::pow(gMaxCell, 1.0/3.0);
            gMinCell = std::pow(gMinCell, 1.0/3.0);


            if
            (
                (*farm->turbineModels[t]) != "AFM" &&
                2.0 * wt->rTip < 8.0 * gMaxCell
            )
            {
                char warning[512];
                sprintf(warning, "turbine diameter (%lf) < 8 mesh cells (%lf), not resolved properly. Revise mesh to improve resolution.\n", 2.0 * wt->rTip, 10.0 * gMaxCell);
                warningInFunction("checkTurbineMesh",  warning);
            }

            //additional checks for Advanced actuator line method
            if(((*farm->turbineModels[t]) == "ALM"))
            {
                if(wt->alm.projectionType == "anisotropic")
                {
                    //check that there are enough radial elements
                    PetscReal cellAvg;

                    cellAvg = 0.5 * (gMaxCell + gMinCell);

                    PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                    PetscReal drvalOpt = 2.0 *  cellAvg;

                    PetscInt nRadial = PetscInt((wt->rTip - wt->rHub)/drvalOpt) + 1;

                    if(drval > drvalOpt)
                    {
                        char error[512];
                        sprintf(error, "Not enough radial elements in the Actuator line mesh. Increase the radial resolution to have %ld points\n", nRadial);
                        fatalErrorInFunction("checkTurbineMesh",  error);
                    }

                    if(drval < gMinCell)
                    {
                        char error[512];
                        sprintf(error, "Too many radial elements. Radial element size smaller than mesh size. Decrease the radial resolution to have %ld points or make the mesh finer\n", nRadial);
                        fatalErrorInFunction("checkTurbineMesh",  error);
                    }

                    //check the mesh size with respect to the projection radius
                                    // number of points in the AL mesh
                    PetscInt npts_t = wt->alm.nPoints;

                    // create temporary vectors
                    std::vector<PetscReal> lminDist(npts_t);
                    std::vector<PetscReal> gminDist(npts_t);
                    std::vector<Cmpnts> perturb(npts_t);

                    // loop over the AD mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        // initialize min dists to a big value
                        lminDist[p] = 1e20;
                        gminDist[p] = 1e20;

                        // set point perturbation
                        perturb[p].x =  procContrib;
                        perturb[p].y =  procContrib;
                        perturb[p].z =  procContrib;

                        // save this point locally for speed
                        Cmpnts point_p = wt->alm.points[p];

                        // perturb the point position
                        mSum(point_p, perturb[p]);

                        // find the closest cell center
                        PetscReal  r_c_minMag = 1e20;
                        cellIds closestCell;

                        // loop over the sphere cells
                        for(c=0; c<wt->nControlled; c++)
                        {
                            // cell indices
                            PetscInt i = wt->controlledCells[c].i,
                                    j = wt->controlledCells[c].j,
                                    k = wt->controlledCells[c].k;

                            // compute distance from mesh cell to AD point
                            Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                            // compute magnitude
                            PetscReal r_c_mag = nMag(r_c);

                            if(r_c_mag < r_c_minMag)
                            {
                                r_c_minMag = r_c_mag;
                                closestCell.i = i;
                                closestCell.j = j;
                                closestCell.k = k;
                            }
                        }

                        // save closest cell indices
                        wt->alm.closestCells[p].i = closestCell.i;
                        wt->alm.closestCells[p].j = closestCell.j;
                        wt->alm.closestCells[p].k = closestCell.k;

                        // save min dist
                        lminDist[p] = r_c_minMag;
                    }

                    MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->alm.nPoints, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

                    for(p=0; p<npts_t; p++)
                    {
                        // point is controlled
                        if(lminDist[p] == gminDist[p])
                        {
                            wt->alm.thisPtControlled[p] = 1;
                        }
                        // point is not controlled
                        else
                        {
                            wt->alm.thisPtControlled[p] = 0;
                        }
                    }

                    // clean memory
                    std::vector<PetscReal> ().swap(lminDist);
                    std::vector<PetscReal> ().swap(gminDist);
                    std::vector<Cmpnts> ().swap(perturb);

                    // loop over the AL mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->alm.thisPtControlled[p])
                        {
                            // get the closest cell center
                            PetscInt i = wt->alm.closestCells[p].i,
                                     j = wt->alm.closestCells[p].j,
                                     k = wt->alm.closestCells[p].k;

                            PetscReal cellsize = std::pow(1/aj[k][j][i], 1.0/3.0);

                            PetscInt nRadial = PetscInt(0.8 * wt->alm.nRadial);

                            PetscInt radPt = PetscInt (p /wt->alm.nAzimuth);

                            PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                            PetscReal eps_x =     wt->alm.chord[p]*wt->eps_x,
                                      eps_y =     wt->alm.thick[p]*wt->eps_y,
                                      eps_z =     drval * wt->eps_z;

                            // PetscPrintf(PETSC_COMM_SELF, "radial point = %ld, nRadial = %ld, eps_x = %lf \n", radPt, nRadial, eps_x);
                            // excluding the tip points for this check
                            if(radPt <= nRadial)
                            {
                                if(eps_x < 1.2 * cellsize)
                                {
                                    char error[512];
                                    sprintf(error, "Fluid Mesh size not optimal for AALM simulation. At radial distance %lf m (%0.2lf %%) from blade root, %lf * chordLength (%lf) < 1.5 * cell size (%lf). Refine mesh or increase epsilonFactor_x (Note: epsilonFactor_x optimal <= 1)\n", drval * radPt, drval * radPt * 100/(wt->rTip - wt->rHub), wt->eps_x, wt->alm.chord[p], cellsize);
                                    fatalErrorInFunction("checkTurbineMesh",  error);
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    return(0);
}