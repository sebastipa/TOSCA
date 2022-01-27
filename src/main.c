// include base headers and objects
#include "include/base.h"

// object declarations
#include "include/domain.h"

#include "include/initialization.h"
#include "include/inline.h"

static char head[] = "TOSCA - Developed at UBC Okanagan CFDLab, Kelowna";

int main(int argc, char **argv)
{
    // initialize PETSc
    PetscInitialize(&argc, &argv, (char *)0, head);

    // domains array
    domain_ *domain;

    // simulation clock
    clock_  clock;
    ReadTimeControls(&clock);

    // simulation flags
    flags_ flags;
    SetSimulationFlags(&flags);

    // simulation info
    simInfo_ info;
    SetSimulationInfo(&info);

    // simulation time
    PetscReal solutionTimeStart, solutionTimeEnd;
    PetscReal iterationTimeStart, iterationTimeEnd;
    PetscTime(&solutionTimeStart);

    simulationInitialize(&domain, &clock, &info, &flags);

    while(clock.time < clock.endTime)
    {
        PetscTime(&iterationTimeStart);

        // adjust time step
        adjustTimeStep(domain);

        setRunTimeWrite(domain);

        for(PetscInt d=0; d<info.nDomains; d++)
        {
            // create old fields
            VecCopy(domain[d].ueqn->Ucont, domain[d].ueqn->Ucont_o);

            if(flags.isTeqnActive)
            {
                VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
            }

            if(domain[d].ueqn->centralUpwindDiv || flags.isTeqnActive)
            {
                UpdateFluxLimiter(domain[d].ueqn);
            }

            if(flags.isLesActive)
            {
                UpdateCs (domain[d].les);
                UpdateNut(domain[d].les);
                UpdateWallModels(domain[d].ueqn);
            }

            if(flags.isAblActive)
            {
                CorrectSourceTerms(domain[d].ueqn, 1);
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

            // Predictor Step
            SolveUEqn(domain[d].ueqn);

            // Pressure Correction
            SolvePEqn(domain[d].peqn);

            // transform contravariant to cartesian
            contravariantToCartesian(domain[d].ueqn);

            // temperature step
            if(flags.isTeqnActive)
            {
                SolveTEqn(domain[d].teqn);

                // save temperature equation right hand side
                VecSet(domain[d].teqn->Rhs_o, 0.0);
                FormT (domain[d].teqn, domain[d].teqn->Rhs_o, 1.0);
            }

            MPI_Barrier(domain[d].mesh->MESH_COMM);

            // print time step continuity errors
            ContinuityErrors(domain[d].peqn);

            // save momentum right hand side
            VecSet(domain[d].ueqn->Rhs_o, 0.0);
            FormU (domain[d].ueqn, domain[d].ueqn->Rhs_o, 1.0);

            if(flags.isIBMActive)
            {
                if(atleastOneIBM(domain[d].ibm))
                {
                    PetscReal ibmTimeStart, ibmTimeEnd;
                    PetscTime(&ibmTimeStart);

                    ibmUpdate(domain[d].ibm);

                    PetscTime(&ibmTimeEnd);
                    PetscPrintf(domain[d].mesh->MESH_COMM, "ibm update time = %lf s\n", ibmTimeEnd - ibmTimeStart);
                }
            }

            if(flags.isTeqnActive)
            {
                UpdateTemperatureBCs(domain[d].teqn);
            }

            // update cartesian BC
            UpdateCartesianBCs(domain[d].ueqn);

            // update contravariant BC
            UpdateContravariantBCs(domain[d].ueqn);

        }

        WriteAcquisition(domain);

        if(flags.isOversetActive)
        {
            UpdateOversetInterpolation(domain);
        }

        clock.it ++;

        PetscTime(&iterationTimeEnd);

        PetscPrintf(PETSC_COMM_WORLD, "Total iteration time = %lf s\n", iterationTimeEnd - iterationTimeStart);

    }

    MPI_Barrier(PETSC_COMM_WORLD);

    PetscTime(&solutionTimeEnd);

    PetscPrintf(PETSC_COMM_WORLD, "\n\nIterations = %ld, Cpu Time = %lf s, Finalizing parallel run\n", clock.it-1, solutionTimeEnd - solutionTimeStart);

    PetscPrintf(PETSC_COMM_WORLD, "\n\nEnd\n\n");

    PetscFinalize();

    return(0);
}
