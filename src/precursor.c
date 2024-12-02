//! \file  precursor.c
//! \brief Contains top to bottom level routines for the concurrent precursor method.

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"
#include "include/initialization.h"
#include "include/initialField.h"

//***************************************************************************************************************//

PetscErrorCode SetSolutionFlagsPrecursor(domain_ *domain)
{
    flags_ *flags = &(domain->flags);


    // read from control file
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    flags->isOversetActive               = 0;
    flags->isLesActive                   = 1;
    flags->isTeqnActive                  = 0;
    flags->isWindFarmActive              = 0;
    flags->isAquisitionActive            = 0;
    flags->isAblActive                   = 0;
    flags->isIBMActive                   = 0;
    flags->isZDampingActive              = 0;
    flags->isXDampingActive              = 0;
    flags->isYDampingActive              = 0;
    flags->isKLeftRayleighDampingActive  = 0;
    flags->isKRightRayleighDampingActive = 0;
    flags->isAdvectionDampingActive      = 0;
    flags->isCanopyActive                = 0;
    flags->isConcurrentPrecursorActive   = 0;
    flags->isPrecursorSpinUp             = 0;
    flags->isPvCatalystActive            = 0;
    flags->isGravityWaveModelingActive   = 0;
    flags->isNonInertialFrameActive      = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-les",            &(flags->isLesActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-potentialT",     &(flags->isTeqnActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-abl",            &(flags->isAblActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-adjustTimeStep", &(flags->isAdjustableTime), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-precursorSpinUp",&(flags->isPrecursorSpinUp), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-precursorIbm",             &(flags->isIBMActive), PETSC_NULL);

    // set acquisition flags
    PetscInt isProbesActive         = 0;
    PetscInt isSectionsActive       = 0;
    PetscInt isAverageABLActive     = 1;
    PetscInt isAverage3LMActive     = 0;
    PetscInt isPhaseAveragingActive = 0;
    PetscInt isAveragingActive      = 0;
    PetscInt isKEBudgetsActive      = 0;
    PetscInt isQCritActive          = 0;
    PetscInt isL2CritActive         = 0;
    PetscInt isPerturbABLActive     = 0;

    flags->isAquisitionActive
    =
    PetscMin((PetscInt)
    (
        isProbesActive + isSectionsActive + isAverageABLActive + isAverage3LMActive + isKEBudgetsActive +
        isAveragingActive + isPhaseAveragingActive + isQCritActive + isL2CritActive + isPerturbABLActive),
        1
    );

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode concurrentPrecursorInitialize(abl_ *abl)
{
    PetscPrintf(abl->access->mesh->MESH_COMM, "Initializing concurrent precursor:\n");

    // allocate memory
    abl->precursor = new precursor_;
    precursor_ *precursor = abl->precursor;

    precursor->domain = new domain_;
    domain_    *domain    = precursor->domain;

    // copy simulation physical constants
    domain->constants.Pr    = abl->access->constants->Pr;
    domain->constants.nu    = abl->access->constants->nu;
    domain->constants.rho   = abl->access->constants->rho;

    // set solution flags
    SetSolutionFlagsPrecursor(domain);

    // set simulation info
    SetSimulationInfo(&(domain->info));

    // allocate domain memory
    SetDomainMemory(domain);

    // set precursor clock parameters
    SetStartTimePrecursor(domain, abl);

    // set access database pointers
    SetAccessPointers(domain);

    // initialize mesh and create mapping
    InitializeMeshPrecursor(abl);

    // only fringe processors proceed
    if(precursor->thisProcessorInFringe)
    {
		// activate precursor flag
	    abl->access->flags->isConcurrentPrecursorActive = 1;

        // set I/O
        InitializeIO(domain->io);

        // override any possible IO
        domain->io->averaging      = 0;
        domain->io->phaseAveraging = 0;
        domain->io->keBudgets      = 0;
        domain->io->qCrit          = 0;
        domain->io->l2Crit         = 0;
        domain->io->windFarmForce  = 0;
        domain->io->sources        = 0;

        // set wall models
        SetWallModels(domain->ueqn);

        // overwrite mesh name after initialization
        domain->mesh->meshName = "precursor";

        // initialize ibm
        InitializeIBM(domain->ibm);

        // copy ABL information into precursor and inflowFunction
        ABLInitializePrecursor(domain);

        // initialize equations
        InitializeUEqn(domain->ueqn);
        InitializePEqn(domain->peqn);

        if(domain->flags.isTeqnActive)
        {
            InitializeTEqn(domain->teqn);
        }
        if(domain->flags.isLesActive)
        {
            InitializeLES(domain->les);
        }

        // initialize acquisition
        InitializeAcquisitionPrecursor(domain);
    }

    MPI_Barrier(abl->access->mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetStartTimePrecursor(domain_ *domain, abl_ *abl)
{
    // set precursor clock pointer
    domain->clock = abl->access->clock;

    // check that the startTime is present in the available times
    if(domain->flags.isPrecursorSpinUp == 0)
    {
        word location = "./fields/precursor/";
        std::vector<PetscReal>        timeSeries;
        PetscInt                      ntimes;
        getTimeList(location.c_str(), timeSeries, ntimes);

        if(timeSeries[ntimes-1] != domain->clock->startTime)
        {
            char error[512];
            sprintf(error, "startTime value %f not available in ./fields/precursor/. Found %f", domain->clock->startTime, timeSeries[ntimes-1]);
            fatalErrorInFunction("SetStartTimePrecursor",  error);
        }

        std::vector<PetscReal> ().swap(timeSeries);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeMeshPrecursor(abl_ *abl)
{
    precursor_ *precursor = abl->precursor;
    domain_    *domain    = precursor->domain;
    flags_     *flags     = precursor->domain->access.flags;

    // get precursor and successor meshes
    mesh_         *mesh_s = abl->access->mesh;
    mesh_         *mesh_p = domain->mesh;

    // set mesh name and type equal to successor (name will be changed after initialization)
    mesh_p->meshName     = mesh_s->meshName;
    mesh_p->meshFileType = "cartesian";

    // set boundary conditions
    SetBoundaryConditionsPrecursor(mesh_p);

    // get domain context information
    DM            da = mesh_s->da, fda = mesh_s->fda;
    DMDALocalInfo info = mesh_s->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscMPIInt   rank, rankP = -1, nProcs, nProcsP, pI;
    PetscInt      i, j, k, commColor = 2;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***cent, ***coor;

    Vec           Coor;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    PetscPrintf(mesh_s->MESH_COMM, "    - creating successor/precursor 1 to 1 processor mapping\n");

    MPI_Comm_rank(mesh_s->MESH_COMM, &rank);
    MPI_Comm_size(mesh_s->MESH_COMM, &nProcs);

    std::vector<PetscMPIInt> lfringeControlled(nProcs);
    std::vector<PetscMPIInt> gfringeControlled(nProcs);

    for(pI=0; pI<nProcs; pI++)
    {
        lfringeControlled[pI] = 0;
        gfringeControlled[pI] = 0;
    }

    // find baricentric i and j ids on this processor
    PetscInt iP        = std::floor((lxs + lxe) / 2.0);
    PetscInt jP        = std::floor((lys + lye) / 2.0);

    // x damping layer parameters
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;

    DMDAVecGetArray(fda, mesh_s->lCent, &cent);

    // see if this processor is in the fringe
    PetscReal xSP    = cent[lzs  ][jP][iP].x;
    PetscReal xEP    = cent[lze-1][jP][iP].x;

    // min k and max k indices
    PetscInt lmink   = 0, gmink;
    PetscInt lmaxk   = 0, gmaxk;

    // number of ranks in the new communicator
    PetscMPIInt lranksp = 0, granksp;

    if
    (
        //(xSP >= xS && xSP <= xE) // the left-most cell center is in/on the fringe
        //||
        //(xEP >= xS && xEP <= xE) // the right-most cell center is in/on the fringe
        xSP <= xE && xEP >= xS
    )
    {
        lfringeControlled[rank]    = 1;
        commColor                  = 1;

        lmink                      = lzs;
        lmaxk                      = lze;

        lranksp                    = 1;
    }

    // scatter control information
    MPI_Allreduce(&lfringeControlled[0], &gfringeControlled[0], nProcs, MPI_INT, MPI_SUM, mesh_s->MESH_COMM);
    precursor->thisProcessorInFringe = lfringeControlled[rank];

    // scatter min and max ids
    MPI_Allreduce(&lmink, &gmink, 1, MPIU_INT, MPI_MIN, mesh_s->MESH_COMM);
    MPI_Allreduce(&lmaxk, &gmaxk, 1, MPIU_INT, MPI_MAX, mesh_s->MESH_COMM);

    // scatter number of ranks
    MPI_Allreduce(&lranksp, &granksp, 1, MPI_INT, MPI_SUM, mesh_s->MESH_COMM);

    // create new communicator
    PetscSubcomm  psubcomm;
    PetscSubcommCreate(mesh_s->MESH_COMM, &(psubcomm));
    PetscSubcommSetTypeGeneral(psubcomm, commColor, rank);
    mesh_p->MESH_COMM  = PetscSubcommChild(psubcomm);

    PetscPrintf(mesh_p->MESH_COMM,"    - creating precursor communicator\n");

    // build the new DMDA info
    PetscInt M, N, P = 0;
    PetscInt m, n, p;

    DMBoundaryType bx, by,       bz = DM_BOUNDARY_PERIODIC;
    if(flags->isPrecursorSpinUp) bz = DM_BOUNDARY_GHOSTED;
    DMDAGetInfo(da, PETSC_NULL, &M, &N, PETSC_NULL, &m, &n, PETSC_NULL, PETSC_NULL, PETSC_NULL, &bx, &by, PETSC_NULL, PETSC_NULL);

    // set number of processors in the k direction
    p = granksp / (m * n);

    // get each processor's ownership range in the successor domain
    const PetscInt *lx, *ly, *lzsucc, *lzprec;
    DMDAGetOwnershipRanges(da, &lx, &ly, &lzsucc);

    // initialize ownership range in k direction
    PetscInt* llz; PetscMalloc(sizeof(PetscInt)*p, &(llz));
    PetscInt* glz; PetscMalloc(sizeof(PetscInt)*p, &(glz));

    for(pI=0; pI<p; pI++)
    {
        llz[pI] = 0;
    }

    // build the ownership range by fixing j and i
    jP = std::floor(my/2.0);
    iP = std::floor(mx/2.0);

    // find min and max fringe k-indices
    PetscReal lminDistS = 1e30, gminDistS = 1e30;
    PetscReal lminDistE = 1e30, gminDistE = 1e30;
    PetscInt lfringeMinK = 0;
    PetscInt lfringeMaxK = 0;
    // access the line of processors in k that have i and j indices equal to mx/2 and my/2
    if
    (
        (iP < xe && iP >= xs)
        &&
        (jP < ye && jP >= ys)
    )
    {
        for(k=zs; k<lze; k++)
        {
            PetscReal distS = fabs(coor[k][jP][iP].x - xS);
            PetscReal distE = fabs(coor[k][jP][iP].x - xE);

            if(distS < lminDistS)
            {
                lminDistS   = distS;
                lfringeMinK = k;
            }

            if(distE < lminDistE)
            {
                lminDistE   = distE;
                lfringeMaxK = k;
            }
        }
    }

    MPI_Allreduce(&lminDistS, &gminDistS, 1, MPIU_REAL, MPI_MIN, mesh_s->MESH_COMM);
    MPI_Allreduce(&lminDistE, &gminDistE, 1, MPIU_REAL, MPI_MIN, mesh_s->MESH_COMM);

    // compare to see which processor got the right k
    if(gminDistS != lminDistS) lfringeMinK = 0;
    if(gminDistE != lminDistE) lfringeMaxK = 0;

    MPI_Allreduce(&lfringeMinK, &(precursor->map.kStart), 1, MPIU_INT, MPI_MAX, mesh_s->MESH_COMM);
    MPI_Allreduce(&lfringeMaxK, &(precursor->map.kEnd), 1, MPIU_INT, MPI_MAX, mesh_s->MESH_COMM);

    // access the line of processors in k that have i and j indices equal to mx/2 and my/2
    if
    (
        (iP < xe && iP >= xs)
        &&
        (jP < ye && jP >= ys)
    )
    {
        // test if this processor is entirely or partially contained in the fringe
        if
        (
            xSP <= xE && xEP >= xS
        )
        {
            MPI_Comm_rank(mesh_p->MESH_COMM, &rankP);

            // recast this rank number into an id into the k line
            PetscInt processorIdx = std::floor(rankP / (m * n));

            PetscInt ghosts = 0;

            // test if processor entirely inside the fringe inside the fringe
            if
            (
                coor[zs][jP][iP].x    < xE && coor[zs][jP][iP].x    > xS &&
                coor[lze-1][jP][iP].x < xE && coor[lze-1][jP][iP].x > xS
            )
            {
                // ownership equal to successor (no physical ghosts)
                llz[processorIdx] = lzsucc[processorIdx];
            }
            // test if left intersection
            else if(coor[zs][jP][iP].x <= xS && coor[lze-1][jP][iP].x < xE && coor[lze-1][jP][iP].x >= xS)
            {
                // ownership range equal to intersection between fringe and processor range + 1 physical ghost node
                llz[processorIdx] = ze - precursor->map.kStart + 1;
            }
            // test if right intersection
            else if(coor[lze-1][jP][iP].x >= xE && coor[zs][jP][iP].x <= xE && coor[zs][jP][iP].x > xS)
            {
                // ownership rage equal to intersection between fringe and processor range + 1 physical ghost node
                llz[processorIdx] = precursor->map.kEnd - zs + 1;
            }
            // it covers two cases: fringe entirely contained in processor, even if coincident
            else
            {
                // size of the fringe region + two ghosts
                llz[processorIdx] = precursor->map.kEnd - precursor->map.kStart + 2;
            }

        }
    }

    // scatter z ownership range (must use MPI_INT64_T since we are scattering PetscInt)
    MPI_Allreduce(&llz[0], &glz[0], p, MPIU_INT, MPI_SUM, mesh_s->MESH_COMM);

    // free memory and swap pointers
    PetscFree(llz);
    lzprec = glz; glz = NULL;

    // compute number of elements in k
    for(pI=0; pI<p; pI++)
    {
        P = P + lzprec[pI];
    }

    // set n points in GCC directions
    mesh_p->IM = M;
    mesh_p->JM = N;
    mesh_p->KM = P;

    // print some information:
    PetscPrintf(mesh_s->MESH_COMM, "      -> breakdown of processors and their concurrent solution flag\n");
    for(pI=0; pI<nProcs; pI++)
    {
        PetscPrintf(mesh_s->MESH_COMM, "         processor %ld: %ld, kStart = %ld\n", pI, gfringeControlled[pI], precursor->map.kStart);
    }

    PetscPrintf(mesh_s->MESH_COMM, "      -> ownership range with (m, n, p) = (%ld, %ld, %ld) precursor proc decomposition in (M, N, P) = (%ld, %ld, %ld) points:\n", m, n, p, M, N, P);
    for(pI=0; pI<m; pI++)
    {
        PetscPrintf(mesh_s->MESH_COMM, "         m[%ld] = %ld\n", pI, lx[pI]);
    }
    for(pI=0; pI<n; pI++)
    {
        PetscPrintf(mesh_s->MESH_COMM, "         n[%ld] = %ld\n", pI, ly[pI]);
    }
    for(pI=0; pI<p; pI++)
    {
        PetscPrintf(mesh_s->MESH_COMM, "         p[%ld] = %ld\n", pI, lzprec[pI]);
    }

    DMDAVecRestoreArray(fda, mesh_s->lCent, &cent);
    DMDAVecRestoreArray(fda, Coor, &coor);

    MPI_Barrier(mesh_s->MESH_COMM);

    // test if this processor is in the fringe
    if(precursor->thisProcessorInFringe)
    {
        DMDACreate3d
        (
            mesh_p->MESH_COMM,
            bx, by, bz,               // boundary type
            DMDA_STENCIL_BOX,         // stencil type
            M, N, P,                  // global points in x,y,z
            m, n, p,                  // processes in x,y,z
            1,                        // dofs
            3,                        // ghost stencil width
            lx, ly, lzprec,           // lx,ly,lz
            &(mesh_p->da)
        );

        // set up da
        DMSetUp(mesh_p->da);

        // set coordinates on da (xmin,xmax,ymin,ymax,zmin,zmax)
        DMDASetUniformCoordinates(mesh_p->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

        // create fda: distributed array data structure for vectors
        DMGetCoordinateDM(mesh_p->da, &(mesh_p->fda));

        // get information about da and this processor location in it
        DMDAGetLocalInfo(mesh_p->da, &(mesh_p->info));

        // create da: distributed array data structure for symmetric tensors
        DMDACreate3d
        (
            mesh_p->MESH_COMM,
            bx, by, bz,               // boundary type
            DMDA_STENCIL_BOX,         // stencil type
            M, N, P,                  // global points in x,y,z
            m, n, p,                  // processes in x,y,z
            6,                        // dofs
            3,                        // ghost stencil width
            lx, ly, lzprec,           // lx,ly,lz
            &(mesh_p->sda)
        );

        // set up sda
        DMSetUp(mesh_p->sda);

        // initialize coordinates and allocate distributed arrays
        {
            PetscPrintf(mesh_p->MESH_COMM, "      -> initializing mesh...");

            DM       da = mesh_p->da, fda = mesh_p->fda;
            DMDALocalInfo info = mesh_p->info;
            PetscInt xs = info.xs, xe = info.xs + info.xm;
            PetscInt ys = info.ys, ye = info.ys + info.ym;
            PetscInt zs = info.zs, ze = info.zs + info.zm;
            PetscInt mx = info.mx, my = info.my, mz = info.mz;

            Vec      lCoorP, lCoorS, gCoorP;
            Cmpnts   ***lcoorp, ***lcoors, ***gcoorp;

            PetscInt i, j, k;
            PetscInt lxs, lxe, lys, lye, lzs, lze;

            lxe = xe; if (xe == mx) lxe = xe - 1;
            lye = ye; if (ye == my) lye = ye - 1;
            lze = ze; if (ze == mz) lze = ze - 1;

            // global working vectors for writing mesh file
            std::vector<PetscReal> lxpoints(mesh_p->KM-1), lypoints(mesh_p->IM-1), lzpoints(mesh_p->JM-1),
                                   gxpoints(mesh_p->KM-1), gypoints(mesh_p->IM-1), gzpoints(mesh_p->JM-1);

            // get successor data
            DMGetCoordinatesLocal(mesh_s->da, &lCoorS);
            DMDAVecGetArray(mesh_s->fda, lCoorS, &lcoors);

            // get precursor data
            DMGetCoordinates(da, &gCoorP);
            DMGetCoordinatesLocal(da, &lCoorP);

            // set global coordinates on precursor mesh
            VecSet(gCoorP, 0.);
            DMDAVecGetArray(fda, gCoorP, &gcoorp);

            // mapping info
            PetscInt kStart = precursor->map.kStart;
            PetscInt kEnd   = precursor->map.kEnd;

            // read x coords - loop over k,i,j with global indexing
            for (k=0; k<mesh_p->KM-1; k++)
            {
                for (j=0; j<mesh_p->JM-1; j++)
                {
                    for (i=0; i<mesh_p->IM-1; i++)
                    {
                        if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                        {
                            gcoorp[k][j][i].x = lcoors[k+kStart][j][i].x;
                            gcoorp[k][j][i].y = lcoors[k+kStart][j][i].y;
                            gcoorp[k][j][i].z = lcoors[k+kStart][j][i].z;

                            // save coordinates to write the mesh
                            if(i==0 && j==0)
                            {
                                lxpoints[k] = gcoorp[k][j][i].x;
                            }
                            if(j==0 && k==0)
                            {
                                lypoints[i] = gcoorp[k][j][i].y;
                            }
                            if(k==0 && i==0)
                            {
                                lzpoints[j] = gcoorp[k][j][i].z;
                            }
                        }
                    }
                }
            }

            DMDAVecRestoreArray(fda, gCoorP, &gcoorp);

            // fill fda data structure with values from all processes (scatter/gather)
            DMGlobalToLocalBegin(fda, gCoorP, INSERT_VALUES, lCoorP);
            DMGlobalToLocalEnd  (fda, gCoorP, INSERT_VALUES, lCoorP);

            // set local coordinates into da
            DMSetCoordinatesLocal(da, lCoorP);

            // restore successor arrays
            DMDAVecRestoreArray(mesh_s->fda, lCoorS, &lcoors);

            // write mesh
            MPI_Reduce(&lxpoints[0], &gxpoints[0], mesh_p->KM-1, MPIU_REAL, MPIU_SUM, 0, mesh_p->MESH_COMM);
            MPI_Reduce(&lypoints[0], &gypoints[0], mesh_p->IM-1, MPIU_REAL, MPIU_SUM, 0, mesh_p->MESH_COMM);
            MPI_Reduce(&lzpoints[0], &gzpoints[0], mesh_p->JM-1, MPIU_REAL, MPIU_SUM, 0, mesh_p->MESH_COMM);

            // get current process
            PetscMPIInt   rank;
            MPI_Comm_rank(mesh_p->MESH_COMM, &rank);

            if(!rank)
            {
                FILE *f;
                char filen[80];
                PetscInt width = -10;
                sprintf(filen, "precursor.xyz");

                // open a new file
                f = fopen(filen, "w");

                if(bz == DM_BOUNDARY_PERIODIC) PetscFPrintf(mesh_p->MESH_COMM, f, "-kPeriodicType 2\n");
                if(by == DM_BOUNDARY_PERIODIC) PetscFPrintf(mesh_p->MESH_COMM, f, "-jPeriodicType 2\n");
                if(bx == DM_BOUNDARY_PERIODIC) PetscFPrintf(mesh_p->MESH_COMM, f, "-iPeriodicType 2\n");
                PetscFPrintf(mesh_p->MESH_COMM, f, "%ld\t%ld\t%ld\n", mesh_p->KM-1, mesh_p->IM-1, mesh_p->JM-1);

                PetscReal zero = 0.0;

                for (k=0; k<mesh_p->KM-1; k++)
                {
                    PetscFPrintf(mesh_p->MESH_COMM, f, "%*.4f\t%*.4f\t%*.4f\n", width, gxpoints[k], width, zero, width, zero);
                }
                for (i=0; i<mesh_p->IM-1; i++)
                {
                    PetscFPrintf(mesh_p->MESH_COMM, f, "%*.4f\t%*.4f\t%*.4f\n", width, zero, width, gypoints[i], width, zero);
                }
                for (j=0; j<mesh_p->JM-1; j++)
                {
                    PetscFPrintf(mesh_p->MESH_COMM, f, "%*.4f\t%*.4f\t%*.4f\n", width, zero, width, zero, width, gzpoints[j]);
                }

                fclose(f);
            }

            // free memory
            std::vector<PetscReal> ().swap(gxpoints);
            std::vector<PetscReal> ().swap(gypoints);
            std::vector<PetscReal> ().swap(gzpoints);
            std::vector<PetscReal> ().swap(lxpoints);
            std::vector<PetscReal> ().swap(lypoints);
            std::vector<PetscReal> ().swap(lzpoints);

            // initialize mesh fields and evaliate metrics
            DMCreateGlobalVector(fda, &(mesh_p->Cent));
            DMCreateLocalVector (fda, &(mesh_p->lCent));

            VecDuplicate(mesh_p->lCent, &(mesh_p->lCsi));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lEta));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lZet));

            VecDuplicate(mesh_p->lCent, &(mesh_p->lICsi));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lIEta));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lIZet));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lJCsi));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lJEta));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lJZet));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lKCsi));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lKEta));
            VecDuplicate(mesh_p->lCent, &(mesh_p->lKZet));

            VecDuplicate(mesh_p->Cent, &(mesh_p->fluxLimiter));

            DMCreateLocalVector (da,  &(mesh_p->lAj));

            VecDuplicate(mesh_p->lAj, &(mesh_p->lIAj));
            VecDuplicate(mesh_p->lAj, &(mesh_p->lJAj));
            VecDuplicate(mesh_p->lAj, &(mesh_p->lKAj));
            VecDuplicate(mesh_p->lAj, &(mesh_p->lNvert));
            VecDuplicate(mesh_p->lAj, &(mesh_p->lNvert_o));

            DMCreateGlobalVector(da,  &(mesh_p->Nvert));
            DMCreateGlobalVector(da,  &(mesh_p->Nvert_o));

            // set to zero
            VecSet(mesh_p->Nvert, 0.0);
            VecSet(mesh_p->lNvert, 0.0);
            VecSet(mesh_p->Nvert_o, 0.0);
            VecSet(mesh_p->lNvert_o, 0.0);

            PetscPrintf(mesh_p->MESH_COMM, "done\n");

            // set curvilinear coordinates metrics
            PetscPrintf(mesh_p->MESH_COMM, "      -> setting curvilinear coordinates...");
            SetMeshMetrics(mesh_p);
            PetscPrintf(mesh_p->MESH_COMM, "done\n");

            // bounding box initialize
            SetBoundingBox(mesh_p);

            // set inflow functions
            SetInflowFunctionsPrecursor(mesh_p);
        }
    }

    MPI_Barrier(mesh_s->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode concurrentPrecursorSolve(abl_ *abl)
{
    // set pointer to precursor database
    precursor_ *precursor = abl->precursor;
    domain_    *domain    = precursor->domain;
    clock_     *clock     = domain->clock;
    flags_     *flags     = domain->access.flags;

    if(precursor->thisProcessorInFringe)
    {
        setRunTimeWrite(domain);

        // set initial field
        if(clock->it == clock->itStart)
        {
            SetInitialFieldPrecursor(abl);
        }

        if(domain->flags.isIBMActive)
        {
            UpdateIBM(domain->ibm);
        }

        // create old fields
        VecCopy(domain->ueqn->Ucont, domain->ueqn->Ucont_o);

        if(domain->flags.isTeqnActive)
        {
            VecCopy(domain->teqn->Tmprt, domain->teqn->Tmprt_o);
            UpdateWallModelsT(domain->teqn);

            if(domain->teqn->pTildeFormulation)
            {
                ghGradRhoK(domain->teqn);
            }

			Buoyancy(domain->ueqn, 1.0);
        }

        if(domain->ueqn->centralUpwindDiv || domain->flags.isTeqnActive)
        {
            UpdateFluxLimiter(domain->ueqn);
        }

        if(domain->flags.isLesActive)
        {
            UpdateCs (domain->les);
            UpdateNut(domain->les);
            UpdateWallModelsU(domain->ueqn);

            if(domain->flags.isIBMActive)
            {
                if(domain->ibm->wallShearOn)
                {
                    findIBMWallShearChester(domain->ibm);
                }
            }
        }

        if(domain->flags.isAblActive)
        {
			if(domain->abl->controllerActive)
			{
				CorrectSourceTerms(domain->ueqn, 1);
			}
			if(domain->abl->controllerActiveT)
			{
				CorrectSourceTermsT(domain->teqn, 1);
			}
        }

        if(domain->flags.isXDampingActive || domain->flags.isZDampingActive)
        {
            correctDampingSources(domain->ueqn);
        }

        // compute pressure gradient term
        GradP(domain->peqn);

        // Predictor Step
        SolveUEqn(domain->ueqn);

        // Pressure Correction
        SolvePEqn(domain->peqn);

        // transform contravariant to cartesian
        contravariantToCartesian(domain->ueqn);

        // temperature step
        if(domain->flags.isTeqnActive)
        {
			UpdateWallModelsT(domain->teqn);

            SolveTEqn(domain->teqn);
        }

        MPI_Barrier(domain->mesh->MESH_COMM);

        // print time step continuity errors (slower)
        // ContinuityErrors(domain->peqn);

        // print time step continuity errors (optimized)
        ContinuityErrorsOptimized(domain->peqn);

        // save momentum right hand side
        if(domain->ueqn->ddtScheme=="backwardEuler")
        {
            VecSet(domain->ueqn->Rhs_o, 0.0);

            if(domain->flags.isIBMActive)
            {

                if (domain->ibm->IBInterpolationModel == "CURVIB")
                {

                    if(domain->ibm->wallShearOn)
                    {
                        CurvibInterpolationInternalCell(domain->ibm);
                    }
                    else 
                    {
                        if(domain->ibm->curvibType == "CurvibTrilinear")
                        {
                            if(domain->ibm->curvibOrder == "linear")
                            {
                                CurvibInterpolation(domain->ibm);
                            }
                            else if(domain->ibm->curvibOrder == "quadratic")
                            {
                                CurvibInterpolationQuadratic(domain->ibm);
                            }
                            else
                            {
                                char error[512];
                                sprintf(error, "wrong interpolation order chosen. Available options are linear and quadratic\n");
                                fatalErrorInFunction("concurrentPrecursorSolve",  error);
                            }
                        }
                        else if(domain->ibm->curvibType == "CurvibTriangular")
                        {
                            CurvibInterpolationTriangular(domain->ibm);
                        }
                        else
                        {
                            char error[512];
                            sprintf(error, "wrong curvib interpolation type\n");
                            fatalErrorInFunction("concurrentPrecursorSolve", error);
                        }
                    }
                }

                if(domain->ibm->wallShearOn)
                {
                    findIBMWallShearChester(domain->ibm);
                }

                UpdateImmersedBCs(domain->ibm);
            }
            FormU (domain->ueqn, domain->ueqn->Rhs_o, 1.0);
        }

        if(domain->flags.isTeqnActive)
        {
            UpdateTemperatureBCs(domain->teqn);
        }

        // update cartesian BC
        UpdateCartesianBCs(domain->ueqn);

        // update contravariant BC
        UpdateContravariantBCs(domain->ueqn);

        WriteAcquisition(domain);
    }

    // sync processors
    MPI_Barrier(abl->access->mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetBoundaryConditionsPrecursor(mesh_ *mesh)
{
    // overwrite k-left and k-right BCs for concurrent precursor run
    // to interpolated mapped after reading

    flags_ *flags = mesh->access->flags;

    // always take the boundary conditions of the
    word location = "./boundary/" + mesh->meshName + "/";

    // read U boundary conditions
    readVectorBC(location, "U", &(mesh->boundaryU));

    if(flags->isPrecursorSpinUp)
    {
        mesh->boundaryU.kLeft    = "inletFunction";
        mesh->boundaryU.kRight   = "zeroGradient" ;
    }

    // read nut boundary conditions
    if (mesh->access->flags->isLesActive)
    {
        readScalarBC(location, "nut", &(mesh->boundaryNut));

        if(flags->isPrecursorSpinUp)
        {
            mesh->boundaryNut.kLeft  = "inletFunction";
            mesh->boundaryNut.kRight = "zeroGradient" ;
        }
    }

    // read T boundary conditions
    if (mesh->access->flags->isTeqnActive)
    {
        readScalarBC(location, "T", &(mesh->boundaryT));

        if(flags->isPrecursorSpinUp)
        {
            mesh->boundaryT.kLeft    = "inletFunction";
            mesh->boundaryT.kRight   = "zeroGradient" ;
        }
    }

    // check boundary conditions
    checkBoundaryConditions(mesh);

    // get name of precursor mesh file and read perodic connectivity info
    word meshFileName;
    if(mesh->meshName == ".") meshFileName ="mesh.xyz";
    else                      meshFileName = mesh->meshName + ".xyz";

    // read connectivity from successor mesh file
    SetPeriodicConnectivity(mesh, meshFileName);

    // overwrite periodic connectivity info
    if(flags->isPrecursorSpinUp)
    {
        mesh->k_periodic         = 0;
        mesh->kk_periodic        = 0;
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInflowFunctionsPrecursor(mesh_ *mesh)
{
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k, l;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***cent;
    PetscReal     ***aj;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    if(mesh->boundaryU.kLeft == "inletFunction")
    {
        // allocate memory for this patch inlet function data
        PetscMalloc(sizeof(inletFunctionTypes), &mesh->inletF.kLeft);

        // set pointer to this inlet function type
        inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

        // read inlet function type (only 4 (interpolatedMappedInflow) and 7 (periodic with i-shift) are available)
        readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "inletFunctionType",   &(ifPtr->typeU));

        if(ifPtr->typeU == 4)
        {
            // set tBarSelectionType the same as uBarSelectionType
            ifPtr->typeT   = 4; ifPtr->mapT   = 0;
            ifPtr->typeNut = 4; ifPtr->mapNut = 0;

            if(mesh->access->flags->isLesActive)  ifPtr->mapNut = 1;
            if(mesh->access->flags->isTeqnActive) ifPtr->mapT   = 1;

            // create inflow communicator to speed up the initialization
            PetscMPIInt rank;
            PetscInt    commColor = MPI_UNDEFINED;

            if(zs==0)
            {
                commColor = 1;
            }

            MPI_Comm_rank(mesh->MESH_COMM, &rank);
            MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(ifPtr->IFFCN_COMM));

            // read if source mesh is uniform or grading
            readSubDictWord("ABLProperties.dat", "xDampingProperties", "sourceType",  &(ifPtr->sourceType));

            // read interpolation method
            readSubDictWord("ABLProperties.dat", "xDampingProperties", "interpolation",&(ifPtr->interpMethod));

            readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Inflow",   &(ifPtr->n1));
            readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Inflow",   &(ifPtr->n2));
            readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Periods",  &(ifPtr->prds1));
            readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Periods",  &(ifPtr->prds2));
			readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Shift",    &(ifPtr->shift2));
            readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Merge",    &(ifPtr->merge1));

    		if(ifPtr->sourceType == "uniform")
    		{
    			PetscPrintf(mesh->MESH_COMM, "   -> using uniform source mesh type\n");

    			readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth1", &(ifPtr->width1));
    			readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth2", &(ifPtr->width2));

                // height of the inflow database (it is duplicate consider removing)
                ifPtr->inflowHeigth = ifPtr->n1*ifPtr->prds1*ifPtr->width1;

                // height of the inflow database
                ifPtr->avgTopLength = ifPtr->n1*ifPtr->prds1*ifPtr->width1;

                // width of the merging region
                ifPtr->avgTopDelta  = 10.0 * ifPtr->width1;

                // 5 top points coordinates
                PetscMalloc(10*sizeof(PetscReal), &(ifPtr->avgTopPointCoords));

                for (i=0; i<10; i++)
                {
                    ifPtr->avgTopPointCoords[i] = ifPtr->avgTopLength - (10-i)*ifPtr->width1 + 0.5*ifPtr->width1;
                }
            }
    		else if(ifPtr->sourceType == "grading")
    		{
    			PetscPrintf(mesh->MESH_COMM, "   -> using grading source mesh type\n");

    			std::vector<PetscReal>  Zcart;

    			word pointsFileName     = "./inflowDatabase/inflowMesh.xyz";
    			FILE *meshFileID        = fopen(pointsFileName.c_str(), "r");

    			if(!meshFileID)
    			{
    				char error[512];
    				sprintf(error, "cannot open inflow points file %s\n", pointsFileName.c_str());
    				fatalErrorInFunction("SetInflowWeights", error);
    			}
    			else
    			{
    				// read the source mesh file in .xyz format
    				PetscReal bufferDouble;
    				PetscInt  npx, npy, npz;
    				PetscInt  error = fscanf(meshFileID, "%ld %ld %ld\n", &npx, &npy, &npz);

    				if(ifPtr->n1 != npz - 1 || ifPtr->n2 != npy - 1)
    				{
    					char error[512];
    					sprintf(error, "source mesh given in %s and expected number of cells do not match\n", pointsFileName.c_str());
    					fatalErrorInFunction("SetInflowFunctionsPrecursor", error);
    				}

    				Zcart.resize(npz);

    				for (PetscInt k = 0; k < npx; k++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
    				for (PetscInt i = 0; i < npy; i++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
    				for (PetscInt j = 0; j < npz; j++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &Zcart[j]);

    				fclose(meshFileID);

    				// height of the inflow database
    				ifPtr->inflowHeigth = Zcart[npz-1] - Zcart[0];

                    // height of the inflow database
                    ifPtr->avgTopLength = Zcart[npz-1] - Zcart[0];

                    // width of the merging region
                    ifPtr->avgTopDelta  = Zcart[npz-1] - Zcart[npz-11];

                    // 5 top points coordinates
                    PetscMalloc(10*sizeof(PetscReal), &(ifPtr->avgTopPointCoords));

                    for (i=0; i<10; i++)
                    {
                        ifPtr->avgTopPointCoords[i] = 0.5 * (Zcart[npz - 11 + i] + Zcart[npz - 11 + i + 1]);
                    }

    				// wipe vectors
    				std::vector<PetscReal> ().swap(Zcart);
    			}
    		}
    		else
    		{
    			char error[512];
    			sprintf(error, "unknown sourceType in inletFunction type 4, available types are\n    1: uniform\n    2: grading\n");
    			fatalErrorInFunction("SetInflowFunctions",  error);
    		}

            // increase n1 and n2 accounting for side ghost cells
            ifPtr->n1wg = ifPtr->n1 + 2;
            ifPtr->n2wg = ifPtr->n2 + 2;

            // build interpolation weights and find periodized inflow cells indices
            SetInflowWeights(mesh, ifPtr);

            // initialize inflow data
            mappedInflowInitialize(ifPtr);

            if(ifPtr->merge1)
            {
                PetscPrintf(mesh->MESH_COMM, "   -> averaging inflow at 10 top cells...");

                // top average to avoid top oscillations
                PetscMalloc(10*sizeof(Cmpnts),    &(ifPtr->uBarAvgTopX));

                for(j=0; j<10; j++)
                {
                    mSetValue(ifPtr->uBarAvgTopX[j], 0.0);
                }

                if(mesh->access->flags->isTeqnActive)
                {
                    PetscMalloc(10*sizeof(PetscReal), &(ifPtr->tBarAvgTopX));

                    for(j=0; j<10; j++)
                    {
                        ifPtr->tBarAvgTopX[j] = 0.0;
                    }
                }

                // variable to store inflow function data
                std::vector<std::vector<Cmpnts>>    ucat_plane_tmp(ifPtr->n1wg);
                std::vector<std::vector<PetscReal>> t_plane_tmp(ifPtr->n1wg);

                // set to zero
                for(j=0; j<ifPtr->n1wg; j++)
                {
                    ucat_plane_tmp[j].resize(ifPtr->n2wg);
                    t_plane_tmp[j].resize(ifPtr->n2wg);

                    for(i=0; i<ifPtr->n2wg; i++)
                    {
                        mSetValue(ucat_plane_tmp[j][i], 0.0);
                        t_plane_tmp[j][i] = 0.0;
                    }
                }

                PetscInt  ti, nAvg;
                PetscReal ntimes;

                if(ifPtr->mapT)
                {
                    ntimes = std::min(ifPtr->inflowT.nInflowTimes, ifPtr->inflowU.nInflowTimes);
                }
                else
                {
                    ntimes = ifPtr->inflowU.nInflowTimes;
                }

                word      fname_U, fname_T;
                FILE      *fp_U, *fp_T;

                for(ti=0; ti<ntimes; ti++)
                {
                    fname_U = "inflowDatabase/U/" + getArbitraryTimeName(mesh->access->clock, ifPtr->inflowU.inflowTimes[ti]);

                    // open the inflow files and read
                    fp_U = fopen(fname_U.c_str(), "rb");

                    if(!fp_U)
                    {
                        char error[512];
                        sprintf(error, "cannot open file:\n    %s\n", fname_U.c_str());
                        fatalErrorInFunction("SetInflowFunctions",  error);
                    }

                    for(j=0; j<ifPtr->n1wg; j++)
                    {
                        PetscInt err1;
                        err1 = fread(&(ucat_plane_tmp[j][0]), sizeof(Cmpnts), ifPtr->n2wg, fp_U);
                    }

                    fclose(fp_U);

                    // now average the top 10 cells (exclude ghosts)
                    PetscInt jAvg = 0;
                    for(j=ifPtr->n1wg-11; j<ifPtr->n1wg-1; j++)
                    {
                        for(i=1; i<ifPtr->n2; i++)
                        {
                            mSum(ifPtr->uBarAvgTopX[jAvg], ucat_plane_tmp[j][i]);
                        }

                        jAvg++;
                    }

                    // do the same with temperature
                    if(mesh->access->flags->isTeqnActive)
                    {
                        fname_T = "inflowDatabase/T/" + getArbitraryTimeName(mesh->access->clock, ifPtr->inflowT.inflowTimes[ti]);
                        fp_T = fopen(fname_T.c_str(), "rb");

                        if(!fp_T)
                        {
                            char error[512];
                            sprintf(error, "cannot open file:\n    %s\n", fname_T.c_str());
                            fatalErrorInFunction("SetInflowFunctions",  error);
                        }

                        for(j=0; j<ifPtr->n1wg; j++)
                        {
                            PetscInt err2;
                            err2 = fread(&(t_plane_tmp[j][0]), sizeof(PetscReal), ifPtr->n2wg, fp_T);
                        }

                        fclose(fp_T);

                        // now average the top 10 cells (exclude ghosts)
                        PetscInt jAvg = 0;
                        for(j=ifPtr->n1wg-11; j<ifPtr->n1wg-1; j++)
                        {
                            for(i=1; i<ifPtr->n2; i++)
                            {
                                ifPtr->tBarAvgTopX[jAvg] += t_plane_tmp[j][i];
                            }

                            jAvg++;
                        }
                    }
                }

                // number of data summed per level (ntimes times n levels in direction 2)
                nAvg  = (ifPtr->n2-1) * ntimes;

                PetscPrintf(mesh->MESH_COMM, "done\n");

                // now average the top 10 cells (exclude ghosts)
                for(j=0; j<10; j++)
                {
                    mScale(1.0/nAvg, ifPtr->uBarAvgTopX[j]);
                    PetscPrintf(mesh->MESH_COMM, "   - Uavg at %.2f m = (%.2f %.2f %.2f) m/s", ifPtr->avgTopPointCoords[j], ifPtr->uBarAvgTopX[j].x, ifPtr->uBarAvgTopX[j].y, ifPtr->uBarAvgTopX[j].z);

                    if(mesh->access->flags->isTeqnActive)
                    {
                        ifPtr->tBarAvgTopX[j] /= nAvg;
                        PetscPrintf(mesh->MESH_COMM, ", thetaAvg = %.2f K", ifPtr->tBarAvgTopX[j]);
                    }

                    PetscPrintf(mesh->MESH_COMM, "\n");

                }

                // temporary (basically forces zero gradient at the top for velocity)
                mSet(ifPtr->uBarAvgTopX[9], ifPtr->uBarAvgTopX[8]);

                // wipe vectors
                for( j=0; j<ifPtr->n1wg; j++)
                {
                    std::vector<Cmpnts>    ().swap(ucat_plane_tmp[j]);
                    std::vector<PetscReal> ().swap(t_plane_tmp[j]);
                }

            }

			// see if must apply shift
			if (ifPtr->shift2 == 1)
			{
				readSubDictDouble("ABLProperties.dat", "xDampingProperties", "shiftSpeed", &(ifPtr->shiftSpeed));

                // compute y coordinates assuming that the mesh has straight z lines, this is required to speed up the
                // search to find the index from which the shifted velocity has to be sourced
                Cmpnts ***cent;
                DMDAVecGetArray(fda, mesh->lCent, &cent);

                std::vector<PetscReal> ycent(mx, 0.0);
                PetscMalloc(sizeof(PetscReal)*mx, &(ifPtr->ycent));
                PetscMalloc(sizeof(PetscReal)*mx, &(ifPtr->yWeights));
                PetscMalloc(sizeof(PetscInt )*mx, &(ifPtr->yIDs));

                if(zs==0 && ys==0)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        ycent[i] = cent[lzs][lys][i].y;
                    }
                }

                DMDAVecRestoreArray(fda, mesh->lCent, &cent);

                if(zs==0)
                {
                    MPI_Allreduce(&ycent[0], &ifPtr->ycent[0], mx, MPIU_REAL, MPIU_SUM, ifPtr->IFFCN_COMM);
                }

			}
        }
        else
        {
            char error[512];
            sprintf(error, "unknown inletFunctionType in xDampingProperties for concurrent precursor method, available types are\n    4: interpolatedMappedInflow\n");
            fatalErrorInFunction("SetInflowFunctionsPrecursor",  error);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode MapInitialConditionPrecursor(abl_ *abl)
{
    // get pointer to precursor database
    precursor_ *precursor  = abl->precursor;
    domain_    *domain     = precursor->domain;

    // get meshes
    mesh_      *mesh_p     = domain->mesh;
    mesh_      *mesh_s     = abl->access->mesh;

    // map cartesian velocity
    successorPrecursorMapVectorField(abl, abl->access->ueqn->Ucat,  domain->ueqn->Ucat, domain->ueqn->lUcat);

    // map contravariant velocity
    successorPrecursorMapVectorField(abl, abl->access->ueqn->Ucont, domain->ueqn->Ucont, domain->ueqn->lUcont);

    // map pressure
    successorPrecursorMapScalarField(abl, abl->access->peqn->P,  domain->peqn->P, domain->peqn->lP);

    // map temperature
    if(domain->flags.isTeqnActive)
    {
        successorPrecursorMapScalarField(abl, abl->access->teqn->Tmprt,  domain->teqn->Tmprt, domain->teqn->lTmprt);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode successorPrecursorMapVectorField(abl_ *abl, Vec &Source, Vec &Target, Vec &lTarget)
{
    // get pointer to precursor database
    precursor_ *precursor  = abl->precursor;
    domain_    *domain     = precursor->domain;

    // get meshes
    mesh_      *mesh_p     = domain->mesh;
    mesh_      *mesh_s     = abl->access->mesh;

    DMDALocalInfo info = mesh_p->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***source, ***target;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    // mapping info
    PetscInt kStart = abl->precursor->map.kStart;

    DMDAVecGetArray(mesh_s->fda, Source, &source);
    DMDAVecGetArray(mesh_p->fda, Target, &target);

    // read x coords - loop over k,i,j with global indexing
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                target[k][j][i].x = source[k+kStart][j][i].x;
                target[k][j][i].y = source[k+kStart][j][i].y;
                target[k][j][i].z = source[k+kStart][j][i].z;
            }
        }
    }

    DMDAVecRestoreArray(mesh_s->fda, Source, &source);
    DMDAVecRestoreArray(mesh_p->fda, Target, &target);

    DMGlobalToLocalBegin(mesh_p->fda, Target, INSERT_VALUES, lTarget);
    DMGlobalToLocalEnd  (mesh_p->fda, Target, INSERT_VALUES, lTarget);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode successorPrecursorMapScalarField(abl_ *abl, Vec &Source, Vec &Target, Vec &lTarget)
{
    // get pointer to precursor database
    precursor_ *precursor  = abl->precursor;
    domain_    *domain     = precursor->domain;

    // get meshes
    mesh_      *mesh_p     = domain->mesh;
    mesh_      *mesh_s     = abl->access->mesh;

    DMDALocalInfo info = mesh_p->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***source, ***target;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    // mapping info
    PetscInt kStart = abl->precursor->map.kStart;

    DMDAVecGetArray(mesh_s->da, Source, &source);
    DMDAVecGetArray(mesh_p->da, Target, &target);

    // read x coords - loop over k,i,j with global indexing
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                target[k][j][i] = source[k+kStart][j][i];
            }
        }
    }

    DMDAVecRestoreArray(mesh_s->da, Source, &source);
    DMDAVecRestoreArray(mesh_p->da, Target, &target);

    DMGlobalToLocalBegin(mesh_p->da, Target, INSERT_VALUES, lTarget);
    DMGlobalToLocalEnd  (mesh_p->da, Target, INSERT_VALUES, lTarget);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ABLInitializePrecursor(domain_ *domain)
{
    if(domain->abl != NULL)
    {
        mesh_         *mesh = domain->mesh;
        flags_        *flags = domain->access.flags;
        DM            da   = mesh->da, fda = mesh->fda;
        DMDALocalInfo info = mesh->info;

        PetscInt      xs = info.xs, xe = info.xs + info.xm;
        PetscInt      ys = info.ys, ye = info.ys + info.ym;
        PetscInt      zs = info.zs, ze = info.zs + info.zm;
        PetscInt      mx = info.mx, my = info.my, mz = info.mz;

        PetscMPIInt   rank, nProcs;
        PetscInt      i, j, k, l;
        PetscInt      lxs, lxe, lys, lye, lzs, lze;

        Vec           Coor;

        PetscReal     lxmin, lymin, lzmin,       // local min coordinates
                      lxmax, lymax, lzmax,       // local max coordinates
                      xmin, ymin, zmin,          // global min coordinates
                      xmax, ymax, zmax;          // global max coordinates

        PetscInt      nLevels = my-2;            // the vertical direction is the j direction in curvilinear coordinates

        Cmpnts        ***coor, ***cent;          // point and cell center coordinates
        PetscReal     ***aj;

        char          dataLoc[256], fileName[500];

        lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
        lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
        lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

        MPI_Comm_rank(mesh->MESH_COMM, &rank);
        MPI_Comm_size(mesh->MESH_COMM, &nProcs);

        if(mesh->meshFileType != "cartesian")
        {
            char error[512];
            sprintf(error, "ABL capabilites only available for cartesian meshes\n");
            fatalErrorInFunction("ABLInitializePrecursor",  error);
        }

        abl_ *abl =  domain->abl;

        readDictDouble("ABLProperties.dat", "hRough",                    &(abl->hRough));
        readDictDouble("ABLProperties.dat", "uRef",                      &(abl->uRef));
        readDictDouble("ABLProperties.dat", "hRef",                      &(abl->hRef));
        readDictDouble("ABLProperties.dat", "hInv",                      &(abl->hInv));
        readDictDouble("ABLProperties.dat", "dInv",                      &(abl->dInv));
        readDictDouble("ABLProperties.dat", "gInv",                      &(abl->gInv));
        readDictDouble("ABLProperties.dat", "tRef",                      &(abl->tRef));
        readDictDouble("ABLProperties.dat", "gTop",                      &(abl->gTop));
        readDictDouble("ABLProperties.dat", "vkConst",                   &(abl->vkConst));
        readDictDouble("ABLProperties.dat", "smearT",                    &(abl->smear));
        readDictInt   ("ABLProperties.dat", "coriolisActive",            &(abl->coriolisActive));
        readDictInt   ("ABLProperties.dat", "controllerActive",          &(abl->controllerActive));
		readDictInt   ("ABLProperties.dat", "controllerActivePrecursorT",&(abl->controllerActiveT));

        // find friction velocity based on neutral log law
        abl->uTau = abl->uRef * abl->vkConst / std::log(abl->hRef / abl->hRough);

        if(abl->coriolisActive)
        {
            readDictDouble("ABLProperties.dat", "fCoriolis", &(abl->fc));
        }

		if(abl->controllerActive)
		{
			readDictWord  ("ABLProperties.dat", "controllerType",   &(abl->controllerType));

			// overwrite controller type if precursor is periodic, use same as successor otherwise
			if(!flags->isPrecursorSpinUp)
			{
				// use the same as successor (average) because this may introduce
				// streamwise oscillations in the successor from the fringe
				// abl->controllerType = "pressure";
			}
		}

        // initialize some useful parameters used in fringe and velocity controller 'write'
        {
            PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->cellLevels));
            PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->totVolPerLevel));
            PetscMalloc(sizeof(PetscInt)  * nLevels, &(abl->totCelPerLevel));

            // initialize height levels for the velocity controller
            DMDAVecGetArray(fda, mesh->lCent, &cent);
            DMDAVecGetArray(da,  mesh->lAj,  &aj);

            std::vector<PetscReal> lLevels(nLevels);
            std::vector<PetscReal> gLevels(nLevels);
            std::vector<PetscReal> lVolumes(nLevels);
            std::vector<PetscReal> gVolumes(nLevels);
            std::vector<PetscInt>  lCells(nLevels);
            std::vector<PetscInt>  gCells(nLevels);

            for(l=0; l<nLevels; l++)
            {
                lLevels[l]  = 0.0;
                gLevels[l]  = 0.0;
                lVolumes[l] = 0.0;
                gVolumes[l] = 0.0;
                lCells[l]   = 0;
                gCells[l]   = 0;
            }

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        lLevels[j-1]  += (cent[k][j][i].z - mesh->grndLevel);
                        lVolumes[j-1] += 1.0 / aj[k][j][i];
                        lCells[j-1]++;
                    }
                }
            }

            MPI_Allreduce(&lLevels[0], &gLevels[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lVolumes[0], &gVolumes[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lCells[0], &gCells[0], nLevels, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            for(l=0; l<nLevels; l++)
            {
                gLevels[l] = gLevels[l] / gCells[l];
            }

            DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
            DMDAVecRestoreArray(fda, mesh->lCent, &cent);

            for(l=0; l<nLevels; l++)
            {
                abl->cellLevels[l]     = gLevels[l];
                abl->totVolPerLevel[l] = gVolumes[l];
                abl->totCelPerLevel[l] = gCells[l];
            }

            std::vector<PetscReal> ().swap(lLevels);
            std::vector<PetscReal> ().swap(gLevels);
            std::vector<PetscReal> ().swap(lVolumes);
            std::vector<PetscReal> ().swap(gVolumes);
            std::vector<PetscInt>  ().swap(lCells);
            std::vector<PetscInt>  ().swap(gCells);
        }

        if(abl->controllerActive)
        {
            readDictDouble("ABLProperties.dat", "controllerMaxHeight", &(abl->controllerMaxHeight));
            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "controllerAction",   &(abl->controllerAction));

            abl->cumulatedSource.x = 0.0;
            abl->cumulatedSource.y = 0.0;
            abl->cumulatedSource.z = 0.0;

            if(abl->controllerAction == "write")
            {
                if(abl->controllerType == "pressure" || abl->controllerType == "geostrophic")
                {

                    readDictDouble("ABLProperties.dat", "relaxPI",          &(abl->relax));
                    readDictDouble("ABLProperties.dat", "alphaPI",          &(abl->alpha));
                    readDictDouble("ABLProperties.dat", "timeWindowPI",     &(abl->timeWindow));

                    // calculating levels interpolation weights at reference height
                    {
                        PetscMalloc(sizeof(PetscReal) * 2,       &(abl->levelWeights));
                        PetscMalloc(sizeof(PetscInt)  * 2,       &(abl->closestLabels));
                        std::vector<PetscReal> absLevelDelta(nLevels);

                        for(l=0; l<nLevels; l++)
                        {
                            absLevelDelta[l] = std::fabs(abl->cellLevels[l] - abl->hRef);
                        }

                        for(PetscInt errI=0; errI<2; errI++)
                        {
                            PetscReal errMin   = 1e20;
                            PetscReal errValue = 0.0;
                            PetscInt  minLabel = 0;

                            for(PetscInt errJ = errI; errJ < nLevels; errJ++)
                            {
                                if(absLevelDelta[errJ] < errMin)
                                {
                                    errValue = absLevelDelta[errJ];
                                    minLabel = errJ;
                                    errMin   = errValue;
                                }
                            }

                            // exchange values so that elements are not ovwerwritten
                            absLevelDelta[minLabel] = absLevelDelta[errI];

                            // put the min value on the unchanged part at the last index of changed part
                            absLevelDelta[errI] = errValue;

                            // save the label adding one since DMDA labeling starts from physical ghost cells
                            abl->closestLabels[errI] = minLabel + 1;
                        }

                        abl->levelWeights[0] = (abl->cellLevels[abl->closestLabels[1]-1]-abl->hRef) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);
                        abl->levelWeights[1] = (abl->hRef-abl->cellLevels[abl->closestLabels[0]-1]) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);

                        std::vector<PetscReal> ().swap(absLevelDelta);

                        if(abl->controllerType == "pressure")
                        {
                            // read if geostrophic damping is active
                            readSubDictInt("ABLProperties.dat", "controllerProperties", "geostrophicDamping", &(abl->geostrophicDampingActive));

                            if(abl->geostrophicDampingActive)
                            {
                                if(!abl->coriolisActive)
                                {
                                    readDictDouble("ABLProperties.dat", "fCoriolis", &(abl->fc));
                                }

                                readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingAlpha",      &(abl->geoDampAlpha));
                                readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingStartTime",  &(abl->geoDampStart));
                                readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingTimeWindow", &(abl->geoDampWindow));

                                abl->geoDampH     = abl->hInv + 0.5 * abl->dInv;
                                abl->geoDampDelta = abl->dInv;
                                abl->geoDampC     = 2.0*(2.0*abl->fc);
                                abl->geoDampUBar  = nSetZero();
                                abl->geoDampAvgS  = nSetZero();
                                abl->geoDampAvgDT = mesh->access->clock->dt;

                                // allocate memory for filtered geostrophic velocity
                                PetscMalloc(sizeof(Cmpnts)*nLevels, &(abl->geoDampU));

                                for(j=0; j<nLevels; j++)
                                {
                                    abl->geoDampU[j] = nSetZero();
                                }

                                // see if must read the average
                                std::stringstream stream;
                                stream << std::fixed << std::setprecision(mesh->access->clock->timePrecision) << mesh->access->clock->startTime;
                                word location = "./fields/" + mesh->meshName + "/" + stream.str();
                                word fileName = location + "/geostrophicDampingInfo";

                                FILE *fp=fopen(fileName.c_str(), "r");

                                if(fp==NULL)
                                {
                                    // if start time > 0 should find the file
                                    if(mesh->access->clock->startTime != 0.0)
                                    {
                                        char error[512];
                                        sprintf(error, "cannot open file %s\n", fileName.c_str());
                                        fatalErrorInFunction("ABLInitialize",  error);
                                    }
                                }
                                else
                                {
                                    fclose(fp);
                                    readDictVector(fileName.c_str(), "filteredS",   &(abl->geoDampAvgS));
                                    readDictDouble(fileName.c_str(), "filteredDT",  &(abl->geoDampAvgDT));
                                    abl->geoDampUBar = nScale(1.0 / (2.0 * abl->fc * abl->geoDampAvgDT), nSetFromComponents(abl->geoDampAvgS.y, abl->geoDampAvgS.x, 0.0));
                                    PetscPrintf(mesh->MESH_COMM, "   -> reading filtered geostrophic wind: Ug = (%.3lf, %.3lf, 0.000)\n",abl->geoDampUBar.x, abl->geoDampUBar.y);
                                }
                            }
                        }
                    }

                    // calculating geostrophic wind and interpolation weight at geostrophic wind
                    if(abl->controllerType == "geostrophic")
                    {
                        // read geosptrophic height
                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "hGeo", &(abl->hGeo));
                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaGeo", &(abl->geoAngle));

                        // initial parameters (should be correct at ABL convergence)
                        abl->geoAngle = abl->geoAngle*M_PI/180;
                        abl->hubAngle = 0.0;
                        abl->omegaBar = 0.0;

                        // compute geostrophic speed
                        abl->uGeoBar  = nSetFromComponents(NieuwstadtGeostrophicWind(abl), 0.0, 0.0);

                        // rotate according to initial angle
                        Cmpnts uGeoBarTmp = nSetZero();
                        uGeoBarTmp.x = std::cos(abl->geoAngle) * abl->uGeoBar.x - std::sin(abl->geoAngle) * abl->uGeoBar.y;
                        uGeoBarTmp.y = std::sin(abl->geoAngle) * abl->uGeoBar.x + std::cos(abl->geoAngle) * abl->uGeoBar.y;
                        mSet(abl->uGeoBar, uGeoBarTmp);

                        // printf information
                        PetscPrintf(mesh->MESH_COMM, "   -> Ug = (%f, %f, %f) m/s, UgMag = %f m/s\n", abl->uGeoBar.x, abl->uGeoBar.y, abl->uGeoBar.z, nMag(abl->uGeoBar));

                        PetscMalloc(sizeof(PetscInt)  * 2,       &(abl->closestLabelsGeo));
                        PetscMalloc(sizeof(PetscReal) * 2,       &(abl->levelWeightsGeo));
                        std::vector<PetscReal> absLevelDelta(nLevels);

                        for(l=0; l<nLevels; l++)
                        {
                            absLevelDelta[l] = std::fabs(abl->cellLevels[l] - abl->hGeo);
                        }

                        for(PetscInt errI=0; errI<2; errI++)
                        {
                            PetscReal errMin   = 1e20;
                            PetscReal errValue = 0.0;
                            PetscInt  minLabel = 0;

                            for(PetscInt errJ = errI; errJ < nLevels; errJ++)
                            {
                                if(absLevelDelta[errJ] < errMin)
                                {
                                    errValue = absLevelDelta[errJ];
                                    minLabel = errJ;
                                    errMin   = errValue;
                                }
                            }

                            // exchange values so that elements are not ovwerwritten
                            absLevelDelta[minLabel] = absLevelDelta[errI];

                            // put the min value on the unchanged part at the last index of changed part
                            absLevelDelta[errI] = errValue;

                            // save the label adding one since DMDA labeling starts from physical ghost cells
                            abl->closestLabels[errI] = minLabel + 1;
                        }

                        abl->levelWeights[0] = (abl->cellLevels[abl->closestLabels[1]-1]-abl->hGeo) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);
                        abl->levelWeights[1] = (abl->hGeo-abl->cellLevels[abl->closestLabels[0]-1]) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);

                        PetscPrintf(mesh->MESH_COMM, "   -> l1 = %lf, l2 = %lf, hGeo = %lf\n", abl->cellLevels[abl->closestLabels[0]-1],abl->cellLevels[abl->closestLabels[1]-1], abl->hGeo);
                        PetscPrintf(mesh->MESH_COMM, "   -> sum of weights = %lf, w1 = %lf, w2 = %lf\n", abl->levelWeights[0]+abl->levelWeights[1], abl->levelWeights[0], abl->levelWeights[1]);

                        std::vector<PetscReal> ().swap(absLevelDelta);
                    }
                }
                else if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") )
                {
                    // read PI controller properties
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaPI",          &(abl->alpha));
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "timeWindowPI",     &(abl->timeWindow));
                    readSubDictInt   ("ABLProperties.dat", "controllerProperties", "avgSources",       &(abl->averageSource));
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "lowestSrcHeight",  &(abl->lowestSrcHt));
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "highestSrcHeight", &(abl->highestSrcHt));

                    // allocate memory for the cumulated sources and average velocity at all mesh heights 
                    PetscMalloc(sizeof(Cmpnts) * nLevels, &(abl->cumulatedSourceHt));

                    PetscMalloc(sizeof(Cmpnts) * nLevels, &(abl->avgVel));

                    for(PetscInt i = 0; i < nLevels; i++)
                    {
                        abl->cumulatedSourceHt[i] = nSetZero();
                        abl->avgVel[i] = nSetZero();
                    }

                    PetscPrintf(mesh->MESH_COMM, "   -> controller type: %s\n", abl->controllerType.c_str());
                    
                    readMesoScaleVelocityData(abl);

                    if(abl->averageSource)
                    {
                        PetscReal timeScaleMeso = 0.1 * abl->timeV[abl->numtV - 1] - abl->timeV[0];

                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "movingAvgWindow",       &(abl->tAvgWindow));

                        if(abl->tAvgWindow > timeScaleMeso)
                        {
                            char error[512];
                            sprintf(error, "moving average window for assimilation source averaging too high. reduce to less than %lf\n", timeScaleMeso);
                            fatalErrorInFunction("ABLInitialize",  error); 
                        }

                        PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->avgsrc));

                        for(PetscInt i = 0; i < my-2; i++)
                        {
                            abl->avgsrc[i] = nSetZero();
                        }

                        abl->currAvgtime = 0.0;
                    }

                    //find the interpolation points and weights for the velocity and temperature fields from the available heights 
                    findVelocityInterpolationWeights(abl);

                    // allocate memory for variables
                    PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->luMean));
                    PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->guMean));
                    PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->srcPA));

                    if(abl->controllerType=="indirectProfileAssimilation")
                    {
                        readSubDictInt   ("ABLProperties.dat", "controllerProperties", "polynomialOrder",   &(abl->polyOrder));

                        //precompute the polynomial coefficient matrix using least square regression method.
                        computeLSqPolynomialCoefficientMatrix(abl); 
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "unknown controllerType for controller action write, available types are:\n        1 : pressure\n        2 : geostrophic\n        3 : directProfileAssimilation\n");
                    fatalErrorInFunction("ABLInitializePrecursor",  error);
                }
            }
            else if(abl->controllerAction == "read")
            {
                if(abl->controllerType=="timeSeries" || abl->controllerType=="timeAverageSeries")
                {

                    // use the vector class to append data in the vector (simpler)
                    std::vector<std::vector<PetscReal>> preCompSourcesTmp;

                    // word by word read
                    char word[256];

                    // buffer for read data
                    PetscReal buffer;

                    // time counter
                    PetscInt ntimes;

                    // file is located into the main folder
                    std::ifstream indata;
                    indata.open("inflowDatabase/momentumSource");

                    if(!indata)
                    {
                    char error[512];
                        sprintf(error, "cannot open file inflowDatabase/momentumSource\n");
                        fatalErrorInFunction("ABLInitializePrecursor",  error);
                    }
                    else
                    {
                        std::string tmpStr;

                        // read lines and get number of saved times
                        for (ntimes = 0; std::getline(indata, tmpStr); ntimes++);

                        // first line is header
                        ntimes--;

                        // save the number of times
                        abl->nSourceTimes = ntimes;
                        abl->currentCloseIdx = 0;

                        // go back on top of file
                        indata.close();
                        indata.open("inflowDatabase/momentumSource");

                        // skip header line
                        std::getline(indata, tmpStr);

                        // resize the source table
                        preCompSourcesTmp.resize(ntimes);

                        for(PetscInt t=0; t<ntimes; t++)
                        {
                            // read along the line: time | sourceX | sourceY | sourceZ
                            for(PetscInt i=0; i<4; i++)
                            {
                                indata >> word;
                                std::sscanf(word, "%lf", &buffer);

                                preCompSourcesTmp[t].push_back(buffer);
                            }

                        }

                        indata.close();
                    }

                    // now store the source data into preCompSources and free the temporary variable
                    PetscMalloc(sizeof(PetscReal*) * ntimes, &(abl->preCompSources));
                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        PetscMalloc(sizeof(PetscReal) * 4, &(abl->preCompSources[t]));
                    }

                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        for(PetscInt i=0; i<4; i++)
                        {
                        abl->preCompSources[t][i] =  preCompSourcesTmp[t][i];
                        }
                    }

                    // clean the temporary variables
                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        std::vector<PetscReal> ().swap(preCompSourcesTmp[t]);
                    }

                    // if controllerType is timeAverageSeries then average the preCompSources
                    if(abl->controllerType=="timeAverageSeries")
                    {
                        readDictDouble("ABLProperties.dat", "controllerAvgStartTime", &(abl->sourceAvgStartTime));

                        // check that average start time is in the list
                        if(abl->sourceAvgStartTime < abl->preCompSources[0][0])
                        {
                        char error[512];
                            sprintf(error, "parameter 'controllerAvgStartTime' is lower than the first available time");
                            fatalErrorInFunction("ABLInitializePrecursor",  error);
                        }
                        // check that more than 100 s of history are used to average
                        else if(abl->sourceAvgStartTime > abl->preCompSources[ntimes-1][0] - 100.00)
                        {
                        char error[512];
                            sprintf(error, "Lower 'controllerAvgStartTime' parameter. Average is too poor (less than 100 s)");
                            fatalErrorInFunction("ABLInitializePrecursor",  error);
                        }
                        else
                        {
                            // warn if less then 1000 s of history are used to average
                            if(abl->sourceAvgStartTime > abl->preCompSources[ntimes-1][0] - 1000.00)
                            {
                            char error[512];
                                sprintf(error, "Lower 'controllerAvgStartTime' parameter. Average could be too poor (less than 1000 s)");
                                warningInFunction("ABLInitializePrecursor",  error);
                            }

                            // initialize average counter to zero
                            PetscInt  nAvgSources = 0;
                            PetscInt  timeOldSet  = 0;
                            PetscReal timeOld;

                            abl->avgTimeStep = 0.0;

                            // average source terms
                            for(PetscInt t=0; t<ntimes; t++)
                            {
                                if(abl->preCompSources[t][0] > abl->sourceAvgStartTime)
                                {
                                    if(!timeOldSet)
                                    {
                                        timeOld    = abl->preCompSources[t][0];
                                        timeOldSet = 1;
                                    }

                                    abl->cumulatedSource.x += abl->preCompSources[t][1];
                                    abl->cumulatedSource.y += abl->preCompSources[t][2];
                                    abl->cumulatedSource.z += abl->preCompSources[t][3];
                                    abl->avgTimeStep       += (abl->preCompSources[t][0] - timeOld);
                                    timeOld                =  abl->preCompSources[t][0];
                                    nAvgSources++;
                                }
                            }

                            // divide by total number of sources used for the average
                            abl->cumulatedSource.x = abl->cumulatedSource.x / nAvgSources;
                            abl->cumulatedSource.y = abl->cumulatedSource.y / nAvgSources;
                            abl->cumulatedSource.z = abl->cumulatedSource.z / nAvgSources;
                            abl->avgTimeStep       = abl->avgTimeStep       / nAvgSources;

                            PetscPrintf(mesh->MESH_COMM, "average driving sources = (%e %e %e), average time step = %lf\n\n", abl->cumulatedSource.x, abl->cumulatedSource.y, abl->cumulatedSource.z, abl->avgTimeStep);
                        }
                    }
                }

                else if(abl->controllerType=="timeHeightSeries")
                {
                    std::vector<std::vector<PetscReal>> preCompSourcesTmp;

                    // word by word read
                    char word[256];

                    // buffer for read data
                    PetscReal buffer;

                    // time counter
                    PetscInt ntimes, numLevels;

                    std::ifstream indata;
                    indata.open("inflowDatabase/momentumSource");

                    if(!indata)
                    {
                        char error[512];
                        sprintf(error, "cannot open file inflowDatabase/momentumSource\n");
                        fatalErrorInFunction("ABLInitialize",  error);
                    }
                    else 
                    {
                        std::string tmpStr;

                        ntimes = 0, numLevels = 0;
                        for (PetscInt t = 0; std::getline(indata, tmpStr); t++)
                        {
                            if (!tmpStr.empty())
                            {
                                ntimes++;
                            }
                        }  

                        //first three lines are header
                        ntimes = ntimes-3;
            
                        // save the number of times
                        abl->nSourceTimes = ntimes;
                        abl->currentCloseIdx = 0;

                        // go back on top of file
                        indata.close();
                        indata.open("inflowDatabase/momentumSource");
                                                
                        // skip first line - levels
                        std::getline(indata, tmpStr);

                        //read the number of levels 
                        indata >> word;
                        while(strcmp("time", word) !=0 )
                        {
                            numLevels++;
                            indata >> word;
                        }

                        //save number of levels in the original data
                        abl->numhV = numLevels;

                        // go back on top of file
                        indata.close();
                        indata.open("inflowDatabase/momentumSource");

                        // skip first line - levels
                        std::getline(indata, tmpStr);

                        PetscMalloc(sizeof(PetscReal) * numLevels, &(abl->hV));

                        //read each of the levels
                        for (PetscInt j = 0; j<numLevels; j++)
                        {
                            indata >> abl->hV[j];
                        }

                        // go back on top of file
                        indata.close();
                        indata.open("inflowDatabase/momentumSource");

                        //skip the first 3 lines 
                        for (PetscInt t = 0; t<3; t++)
                        {
                            std::getline(indata, tmpStr);
                        }      

                        // resize the source table
                        preCompSourcesTmp.resize(ntimes);

                        //x,y,z velocity per level + time
                        PetscInt wPerLine = numLevels * 3 + 1;

                        for (PetscInt t = 0; t<ntimes; t++)
                        {
                            for (PetscInt j = 0; j<wPerLine; j++)
                            {
                                indata >> word;
                                std::sscanf(word, "%lf", &buffer);

                                preCompSourcesTmp[t].push_back(buffer);
                            }
                        }                   
                        
                        indata.close();
                    }
                    
                    PetscMalloc(sizeof(PetscReal **) * ntimes, &(abl->timeHtSources));

                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        PetscMalloc(sizeof(PetscReal*) * (my-2), &(abl->timeHtSources[t]));
                    }

                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        for(PetscInt j=0; j<(my-2); j++)
                        {
                            PetscMalloc(sizeof(PetscReal) * 4, &(abl->timeHtSources[t][j]));
                        }
                    }
                    
                    //save the time for each height
                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        for(PetscInt j=0; j<(my-2); j++)
                        {
                            abl->timeHtSources[t][j][0] = preCompSourcesTmp[t][0];
                        }
                    }

                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        for(PetscInt j=0; j<(my-2); j++)
                        {
                            abl->timeHtSources[t][j][1] = preCompSourcesTmp[t][3*j + 1];
                            abl->timeHtSources[t][j][2] = preCompSourcesTmp[t][3*j + 2];
                            abl->timeHtSources[t][j][3] = preCompSourcesTmp[t][3*j + 3];
                        }
                    }     

                    // clean the temporary variables
                    for(PetscInt t=0; t<ntimes; t++)
                    {
                        std::vector<PetscReal> ().swap(preCompSourcesTmp[t]);
                    }

                    if(numLevels != my-2)
                    {
                        //find interpolation ids and weights to interpolate from the original precursor mesh source data to the current mesh
                        findTimeHeightSeriesInterpolationWts(abl);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "unknown controllerType for controller action read, available types are:\n        1 : timeSeries\n        2 : timeAverageSeries\n        3 : timeHeightSeries\n");
                    fatalErrorInFunction("ABLInitializePrecursor",  error);
                }
            
            }
            else 
            {
                    char error[512];
                    sprintf(error, "unknown controllerAction, available types are:\n        1 : write\n        2 : read\n");
                    fatalErrorInFunction("ABLInitializePrecursor",  error);
            }

            if(abl->controllerActiveT)
			{
				PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->tDes));

				for(l=0; l<nLevels; l++)
				{
					abl->tDes[l] = 0.0;
				}

				// read proportional controller relaxation factor (same as the velocity one)
				readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
                readDictWord     ("ABLProperties.dat", "controllerTypeT",    &(abl->controllerTypeT));

                if(abl->controllerTypeT=="indirectProfileAssimilation" || abl->controllerTypeT=="directProfileAssimilation")
                {   
                    readMesoScaleTemperatureData(abl);

                    findTemperatureInterpolationWeights(abl);
                }

                if(abl->controllerTypeT=="indirectProfileAssimilation")
                {
                    readSubDictInt   ("ABLProperties.dat", "controllerProperties", "polynomialOrder",   &(abl->polyOrder));

                    //precompute the polynomial coefficient matrix using least square regression method.
                    if(abl->controllerType !="indirectProfileAssimilation")
                    {
                        computeLSqPolynomialCoefficientMatrix(abl); 
                    }
                }
			}
        }

        // read the Rayleigh damping layer properties
        if(domain->flags.isZDampingActive)
        {
            readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingStart",   &(abl->zDampingStart));
            readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingEnd",     &(abl->zDampingEnd));
            readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingAlpha",   &(abl->zDampingAlpha));
            readSubDictInt   ("ABLProperties.dat", "zDampingProperties", "zDampingAlsoXY",  &(abl->zDampingAlsoXY));

            if(abl->zDampingAlsoXY)
            {
                // set average weight to zero
                abl->avgWeight = 0;

                // allocate memory for uBarMeanZ and set to zero
                PetscMalloc(my*sizeof(Cmpnts), &(abl->uBarMeanZ));

                for(j=0; j<my; j++)
                {
                    abl->uBarMeanZ[j].x = 0.0;
                    abl->uBarMeanZ[j].y = 0.0;
                    abl->uBarMeanZ[j].z = 0.0;
                }
            }
        }

        // read the recycling fringe region properties
        if(domain->flags.isXDampingActive)
        {
            // no fringe region in the concurrent domain
        }
    }

    return(0);
}
