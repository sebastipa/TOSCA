// include base headers and objects
#include "include/base.h"
#include "include/domain.h"
#include "include/initialization.h"
#include "include/inline.h"

#if USE_CATALYST
#include "include/catalystAdaptor.h"
#endif

static char head[] = "TOSCA - Developed at UBC Okanagan CFDLab, Kelowna";

int main(int argc, char **argv)
{
    // initialize PETSc
    PetscInitialize(&argc, &argv, (char *)0, head);

    // domains array
    domain_ *domain;

    // simulation clock
    clock_   clock;
    ReadTimeControls(&clock);

    // simulation flags (initialize)
    flags_ flags;
    SetSimulationFlags(&flags);

    // simulation info (initialize)
    simInfo_ info;
    SetSimulationInfo(&info);

    // simulation time
    PetscReal solutionTimeStart, solutionTimeEnd;
    PetscReal iterationTimeStart, iterationTimeEnd;
    PetscTime(&solutionTimeStart);

    // initialize simulation
    simulationInitialize(&domain, &clock, &info, &flags);

    if(flags.isPvCatalystActive)
    {
        #if USE_CATALYST
        catalystInitialize(domain);
        #endif
    }

    while(clock.endTime - clock.time > 1e-10)
    {
        PetscTime(&iterationTimeStart);

        // adjust time step
        adjustTimeStep(domain);

        // re-read input file
        RereadIO(domain);

        setRunTimeWrite(domain);

        for(PetscInt d=0; d<info.nDomains; d++)
        {
            // reset flags based on domain preferences
            flags = domain[d].flags;

            //update the IBM position and interpolate based on it before the start of solution
            if(flags.isIBMActive)
            {
                UpdateIBM(domain[d].ibm);
            }

            // copy old fields - contains all the bc for this time step - after all the boundary updates
            VecCopy(domain[d].ueqn->Ucont, domain[d].ueqn->Ucont_o);

            if(flags.isTeqnActive)
            {
                VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);

                if(domain[d].teqn->pTildeFormulation)
                {
                    ghGradRhoK(domain[d].teqn);
                }

                Buoyancy(domain[d].ueqn, 1.0);
            }

            if(domain[d].ueqn->centralUpwindDiv || flags.isTeqnActive)
            {
                UpdateFluxLimiter(domain[d].ueqn);
            }

            if(flags.isAblActive)
            {
                if(domain[d].abl->controllerActive)
                {
                    CorrectSourceTerms(domain[d].ueqn, 1);
                }
                if(domain[d].abl->controllerActiveT && flags.isTeqnActive)
                {
                    CorrectSourceTermsT(domain[d].teqn, 1);
                }
            }

            if(flags.isLesActive)
            {
                UpdateCs (domain[d].les);
                UpdateNut(domain[d].les);
                UpdateWallModelsU(domain[d].ueqn);

                if(flags.isIBMActive)
                {
                    if(domain[d].ibm->wallShearOn)
                    {
                        findIBMWallShearChester(domain[d].ibm);
                    }
                }
            }

            if(flags.isXDampingActive || flags.isZDampingActive)
            {
                correctDampingSources(domain[d].ueqn);
            }

            // update wind turbines
            if(flags.isWindFarmActive)
            {
                UpdateWindTurbines(domain[d].farm);
            }

            // compute pressure gradient term
            GradP(domain[d].peqn);

            // finish the y-damping layer processor mapping 
            if(flags.isYDampingActive)
            {
                mapYDamping(domain[d].ueqn);
            }

            // Predictor Step
            SolveUEqn(domain[d].ueqn);

            // Pressure Correction
            SolvePEqn(domain[d].peqn);

            // transform contravariant to cartesian
            contravariantToCartesian(domain[d].ueqn);

            // temperature step
            if(flags.isTeqnActive)
            {
                UpdateWallModelsT(domain[d].teqn);

                SolveTEqn(domain[d].teqn);
            }

            MPI_Barrier(domain[d].mesh->MESH_COMM);

            // print time step continuity errors (slower)
            // ContinuityErrors(domain[d].peqn);

            // print time step continuity errors (optimized)
            ContinuityErrorsOptimized(domain[d].peqn);

            // save momentum right hand side
            if(domain[d].ueqn->ddtScheme=="backwardEuler")
            {
                VecSet(domain[d].ueqn->Rhs_o, 0.0);

                //interpolate IBM cells before computing the forces and moments on the IBM
                if(flags.isIBMActive)
                {

                    if (domain[d].ibm->IBInterpolationModel == "CURVIB")
                    {

                        if(domain[d].ibm->wallShearOn)
                        {
                            CurvibInterpolationInternalCell(domain[d].ibm);
                        }
                        else 
                        {
                            if(domain[d].ibm->curvibType == "CurvibTrilinear")
                            {
                                if(domain[d].ibm->curvibOrder == "linear")
                                {
                                    CurvibInterpolation(domain[d].ibm);
                                }
                                else if(domain[d].ibm->curvibOrder == "quadratic")
                                {
                                    CurvibInterpolationQuadratic(domain[d].ibm);
                                }
                                else
                                {
                                    char error[512];
                                    sprintf(error, "wrong interpolation order chosen. Available options are linear and quadratic\n");
                                    fatalErrorInFunction("main",  error);
                                }
                            }
                            else if(domain[d].ibm->curvibType == "CurvibTriangular")
                            {
                                CurvibInterpolationTriangular(domain[d].ibm);
                            }
                            else
                            {
                                char error[512];
                                sprintf(error, "wrong curvib interpolation type\n");
                                fatalErrorInFunction("main", error);
                            }
                        }
                    }

                    if(domain[d].ibm->wallShearOn)
                    {
                        findIBMWallShearChester(domain[d].ibm);
                    }

                    UpdateImmersedBCs(domain[d].ibm);
                }

                FormU (domain[d].ueqn, domain[d].ueqn->Rhs_o, 1.0);
            }

            if(flags.isTeqnActive)
            {
                UpdateTemperatureBCs(domain[d].teqn);
            }

            // update cartesian BC
            UpdateCartesianBCs(domain[d].ueqn);

            // update contravariant BC
            UpdateContravariantBCs(domain[d].ueqn);

            if(flags.isIBMActive)
            {
                if(domain[d].ibm->computeForce)
                {
                    if(domain[d].ibm->dynamic)
                    {
                        IBMProjectionProcessorTransfer(domain[d].ibm);
                    }

                    ComputeForceMoment(domain[d].ibm);
                }
            }

            MPI_Barrier(domain[d].mesh->MESH_COMM);
        }

        WriteAcquisition(domain);

        if(flags.isPvCatalystActive)
        {
            #if USE_CATALYST
            catalystExecute(domain);
            #endif
        }

        if(flags.isOversetActive)
        {
            UpdateOversetInterpolation(domain);
        }

        clock.it ++;

        PetscTime(&iterationTimeEnd);

        PetscPrintf(PETSC_COMM_WORLD, "Total iteration time = %lf s\n", iterationTimeEnd - iterationTimeStart);

        MPI_Barrier(PETSC_COMM_WORLD);

    }

    #if USE_CATALYST
    catalystFinalize();
    #endif

    PetscTime(&solutionTimeEnd);

    PetscPrintf(PETSC_COMM_WORLD, "\n\nIterations = %ld, Cpu Time = %lf s, Finalizing parallel run\n", clock.it-1, solutionTimeEnd - solutionTimeStart);

    PetscPrintf(PETSC_COMM_WORLD, "\n\nEnd\n\n");

    PetscFinalize();

    return(0);
}
