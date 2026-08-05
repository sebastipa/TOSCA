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

    // uncomment to enable PETSc peak-memory tracking
    // PetscMemorySetGetMaximumUsage();

    // declare simulation data
    domain_ *domain;

    // simulation clock
    clock_   clock;
    ReadTimeControls(&clock);

    // initialize simulation flags
    flags_ flags;
    SetSimulationFlags(&flags);

    // initialize simulation info
    simInfo_ info;
    SetSimulationInfo(&info);

    // initialize simulation time
    PetscReal solutionTimeStart, solutionTimeEnd;
    PetscReal iterationTimeStart, iterationTimeEnd;
    PetscTime(&solutionTimeStart);

    // initialize simulation data
    simulationInitialize(&domain, &clock, &info, &flags);

    if(flags.isPvCatalystActive)
    {
        #if USE_CATALYST
        catalystInitialize(domain);
        #endif
    }

    PetscPrintf(PETSC_COMM_WORLD, "\nStarting time loop:\n");
    PetscPrintf(PETSC_COMM_WORLD, "******************************************************************\n\n");

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
            if(flags.isOversetActive)
                PetscPrintf(PETSC_COMM_WORLD, "\nDomain: %ld\n\n", *(domain[d].access.domainID));

            // locally copythis domain flags
            flags = domain[d].flags;

            // update the IBM position and interpolate based on new positions
            if(flags.isIBMActive)
            {
                UpdateIBM(domain[d].ibm);
            }

            // save velocity at old time step
            VecCopy(domain[d].ueqn->Ucont, domain[d].ueqn->Ucont_o);

            if(flags.isTeqnActive)
            {
                // save temperature at old time step
                if(domain[d].teqn->ddtScheme == "BDF2")
                {
                    VecCopy(domain[d].teqn->Tmprt_o, domain[d].teqn->Tmprt_oo); 
                    VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
                }
                else
                {
                    VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
                }
            }

            // update flux limiter
            if(domain[d].ueqn->centralUpwindDiv || domain[d].ueqn->centralUpwindWDiv)
            {
                UpdateFluxLimiter(domain[d].ueqn);
            }

            // solve the alpha water equation
            if(flags.isAeqnActive)
            {
                // alpha water sub-cycling: take multiple smaller timesteps
                PetscInt nSubCycles  = domain[d].aeqn->nAlphaSubCycles;
                PetscReal dtOriginal = domain[d].clock->dt;
                PetscReal dtAlpha    = dtOriginal / nSubCycles;
                
                for (PetscInt subcycle = 0; subcycle < nSubCycles; subcycle++)
                {
                    // temporarily reduce timestep for alpha equation
                    domain[d].clock->dt = dtAlpha;

                    // save alpha at old time step
                    if(domain[d].aeqn->ddtScheme == "BDF2")
                    {
                        VecCopy(domain[d].aeqn->Alpha_o, domain[d].aeqn->Alpha_oo);
                        VecCopy(domain[d].aeqn->Alpha, domain[d].aeqn->Alpha_o);
                    }
                    else
                    {
                        VecCopy(domain[d].aeqn->Alpha, domain[d].aeqn->Alpha_o);  
                    }

                    // advance alpha by one sub-step
                    SolveAEqn(domain[d].aeqn);
                }

                // restore original timestep
                domain[d].clock->dt = dtOriginal;

                // bound alpha water between 0 and 1 & print min max of alpha water for logging
                PetscReal alphaMinPre; VecMin(domain[d].aeqn->Alpha, NULL, &alphaMinPre); 
                PetscReal alphaMaxPre; VecMax(domain[d].aeqn->Alpha, NULL, &alphaMaxPre);
                 
                boundAlpha(domain[d].aeqn);

                PetscReal alphaMinPost; VecMin(domain[d].aeqn->Alpha, NULL, &alphaMinPost);
                PetscReal alphaMaxPost; VecMax(domain[d].aeqn->Alpha, NULL, &alphaMaxPost);
                
                PetscPrintf(domain[d].mesh->MESH_COMM, "Alpha-Water pre/post correction min = %.5f/%.5f, max = %.5f/%.5f\n", alphaMinPre, alphaMinPost, alphaMaxPre, alphaMaxPost);

                // save old time step density before updating rho
                VecCopy(domain[d].aeqn->lRhoFace, domain[d].aeqn->lRhoFace_o);

                // update rho 
                UpdateRho(domain[d].aeqn);

                // compute density gradient
                GradRho(domain[d].aeqn);
            }

            // update SGS models
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

            // correct damping layers/fringe region
            if(flags.isXDampingActive || flags.isZDampingActive)
            {
                correctDampingSources(domain[d].ueqn);
            }

            // update flow controllers source terms
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

            // update wind turbines
            if(flags.isWindFarmActive)
            {
                UpdateWindTurbines(domain[d].farm);
            }

            // compute pressure gradient term
            if(domain[d].ueqn->central4Div)
            {
                GradP4thOrder(domain[d].peqn);
            }
            else 
            {
                GradP(domain[d].peqn);
            }

            // update y-damping layer processor mapping 
            if(flags.isYDampingActive)
            {
                mapYDamping(domain[d].ueqn);
            }

            // buoyancy term
            if(flags.isTeqnActive)
            {
                // save old buoyancy term for AB2 formulation
                if(clock.it > clock.itStart)
                {
                    if(domain[d].teqn->pTildeFormulation)
                        VecCopy(domain[d].teqn->ghGradRhok, domain[d].teqn->ghGradRhok_o);
                    else
                        VecCopy(domain[d].ueqn->bTheta, domain[d].ueqn->bTheta_o);
                }

                // compute new buoyancy term
                if(domain[d].teqn->pTildeFormulation)
                    ghGradRhoK(domain[d].teqn);
                else
                    Buoyancy(domain[d].ueqn, 1.0);
            }

            // solve the momentum equation
            SolveUEqn(domain[d].ueqn);

            // solve the pressure equation
            SolvePEqn(domain[d].peqn);

            // update cartesian velocity
            contravariantToCartesian(domain[d].ueqn);

            // potential temperature step
            if(flags.isTeqnActive)
            {
                // update SGS fields 
                UpdateCsk(domain[d].les);
                UpdatekT(domain[d].les);

                // update wall models
                UpdateWallModelsT(domain[d].teqn);
                
                // advance temperature to n+1
                SolveTEqn(domain[d].teqn);
            }

            MPI_Barrier(domain[d].mesh->MESH_COMM);

            // print time step continuity errors (slower)
            // ContinuityErrors(domain[d].peqn);

            // print time step continuity errors (optimized)
            ContinuityErrorsOptimized(domain[d].peqn);

            // save momentum right hand side
            if(domain[d].ueqn->ddtScheme=="CN")
            {
                // interpolate IBM cells before computing the forces and moments on the IBM
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

                // save full right hand side for next iteration
                VecSet(domain[d].ueqn->Rhs_o, 0.0);
                FormU (domain[d].ueqn, domain[d].ueqn->Rhs_o, 1.0);
            }

            // update temperature BC
            if(flags.isTeqnActive)
            {
                UpdateTemperatureBCs(domain[d].teqn);
            }

            // update alpha water BC
            if(flags.isAeqnActive)
            {
                UpdateAlphaWaterBCs(domain[d].aeqn);
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

        // perform paraview catalyst actions
        if(flags.isPvCatalystActive)
        {
            #if USE_CATALYST
            catalystExecute(domain);
            #endif
        }

        // perform overset interpolation
        if(flags.isOversetActive)
        {
            UpdateOversetInterpolation(domain);
        }

        // remove gauge pressure and sync
        SyncPressureAcrossDomains(domain);

        // write output files
        WriteAcquisition(domain);

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
