//! \file  acquisition.c
//! \brief Acquisition functions definition

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

//***************************************************************************************************************//

PetscErrorCode InitializeAcquisition(domain_ *domain)
{
    // read from control file
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    // get number of domains and flags (same for all domains)
    PetscInt nDomains = domain[0].info.nDomains;
    flags_   *flags   = domain[0].access.flags;
    clock_   *clock   = domain[0].clock;

    // domain folder flag
    PetscInt createDomainFolder = 0;

    // test if some post processing is present
    if(flags->isAquisitionActive || flags->isWindFarmActive)
    {
        // create post processing folder
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        if(!rank)
        {
            errno = 0;
            PetscInt dirRes = mkdir("./postProcessing", 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create postProcessing directory\n");
                fatalErrorInFunction("InitializeAcquisition",  error);
            }
        }

        // read/allocate acquisition objects
        if(flags->isAquisitionActive)
        {
            for(PetscInt d=0; d<nDomains; d++)
            {
                acquisition_ *acquisition = domain[d].acquisition;

                // read acquisition flags
                acquisition->isProbesActive     = 0;
                acquisition->isSectionsActive   = 0;
                acquisition->isAverageABLActive = 0;
                acquisition->isAverage3LMActive = 0;

                PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-probes",        &(acquisition->isProbesActive),     PETSC_NULL);
                PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sections",      &(acquisition->isSectionsActive),   PETSC_NULL);
                PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averageABL",    &(acquisition->isAverageABLActive), PETSC_NULL);
                PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-average3LM",    &(acquisition->isAverage3LMActive), PETSC_NULL);

                // set domain folder flag to 1 is sections are active
                if(acquisition->isSectionsActive) createDomainFolder = 1;
            }
        }

        // set domain folder flag to 1 if wind farm is active
        if(flags->isWindFarmActive) createDomainFolder = 1;

        // create domain folders
        if(createDomainFolder)
        {
            for(PetscInt d=0; d<nDomains; d++)
            {
              if(!rank)
              {
                errno = 0;
                word domainFolderName = "./postProcessing/" + domain[d].mesh->meshName;
                PetscInt dirRes = mkdir(domainFolderName.c_str(), 0777);
                if(dirRes != 0 && errno != EEXIST)
                {
                   char error[512];
                    sprintf(error, "could not create %s directory", domainFolderName.c_str());
                    fatalErrorInFunction("InitializeAcquisition",  error);
                }
              }
            }
        }
    }

    if(flags->isAquisitionActive)
    {
        for(PetscInt d=0; d<nDomains; d++)
        {
            acquisition_ *acquisition = domain[d].acquisition;

            // initialize sections
            sectionsInitialize(acquisition);

            MPI_Barrier(domain[d].mesh->MESH_COMM);

            // initialize averages
            averageFieldsInitialize(acquisition);

        }

        // initialize probes
        ProbesInitialize(domain);

        // initialize 3LM averaging
        averaging3LMInitialize(domain);

        // initialize ABL averaging
        averagingABLInitialize(domain);
    }

    // sync processors
    MPI_Barrier(PETSC_COMM_WORLD);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeAcquisitionPrecursor(domain_ *domain)
{
    // get mesh communicator (precursor has always single domain)
    mesh_ *mesh = domain->mesh;

    // read from control file
    PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

    // get number of domains and flags (same for all domains)
    PetscInt nDomains = domain->info.nDomains;
    flags_ *flags     = domain->access.flags;
    clock_ *clock     = domain->clock;

    // test if some post processing is present
    if(flags->isAquisitionActive)
    {
        PetscMPIInt rank;
        MPI_Comm_rank(mesh->MESH_COMM, &rank);

        if(!rank)
        {
            errno = 0;
            PetscInt dirRes = mkdir("./postProcessing", 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create postProcessing directory\n");
                fatalErrorInFunction("InitializeAcquisition",  error);
            }
        }

        // read/allocate acquisition objects
        if(flags->isAquisitionActive)
        {
            acquisition_ *acquisition = domain->acquisition;

            // read acquisition flags
            acquisition->isProbesActive     = 0;
            acquisition->isSectionsActive   = 0;
            acquisition->isAverageABLActive = 1;
            acquisition->isAverage3LMActive = 0;
        }

        // create precursor folders
        if(!rank)
        {
            errno = 0;
            word domainFolderName = "./postProcessing/" + domain->mesh->meshName;
            PetscInt dirRes = mkdir(domainFolderName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", domainFolderName.c_str());
                fatalErrorInFunction("InitializeAcquisition",  error);
            }
        }
    }

    if(flags->isAquisitionActive)
    {
        acquisition_ *acquisition = domain->acquisition;

        // initialize sections
        sectionsInitialize(acquisition);

        // initialize averages
        averageFieldsInitialize(acquisition);

        // initialize probes
        ProbesInitialize(domain);

        // initialize 3LM averaging
        averaging3LMInitialize(domain);

        // initialize ABL averaging
        averagingABLInitialize(domain);
    }

    // sync processors
    MPI_Barrier(mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode WriteAcquisition(domain_ *domain)
{
    // get number of domains and flags (same for all domains)
    PetscInt nDomains = domain[0].info.nDomains;
    flags_   flags    = domain[0].flags;

    for(PetscInt d=0; d<nDomains; d++)
    {
        if(flags.isAquisitionActive)
        {
            acquisition_ *acquisition = domain[d].acquisition;

            // write sections
            writeSections(acquisition);

            // average fields
            averageFields(acquisition);
        }

        // write the fields
        writeFields(domain[d].io);
    }

    if(flags.isAquisitionActive)
    {
        // write probe data
        writeProbes(domain);

        // write 3LM data
        writeAveraging3LM(domain);

        // write ABL data
        writeAveragingABL(domain);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode averageFieldsInitialize(acquisition_ *acquisition)
{
    io_    *io    = acquisition->access->io;
    mesh_  *mesh  = acquisition->access->mesh;
    flags_ *flags = acquisition->access->flags;

    if(io->averaging || io->phaseAveraging || io->qCrit || io->l2Crit || io->sources)
    {
        PetscMalloc(sizeof(avgFields), &(acquisition->fields));
        avgFields *avg = acquisition->fields;

        if(io->qCrit)
        {
            VecDuplicate(mesh->Nvert, &(avg->Q));  VecSet(avg->Q,0.);
        }

        if(io->l2Crit)
        {
            VecDuplicate(mesh->Nvert, &(avg->L2));  VecSet(avg->L2,0.);
        }

        if(io->sources)
        {
            VecDuplicate(mesh->Cent, &(avg->Coriolis)); VecSet(avg->Coriolis, 0.);
            VecDuplicate(mesh->Cent, &(avg->Driving));  VecSet(avg->Driving,  0.);
            VecDuplicate(mesh->Cent, &(avg->xDamping)); VecSet(avg->xDamping, 0.);
            VecDuplicate(mesh->Cent, &(avg->SideForce));VecSet(avg->SideForce,0.);
        }

        // allocate averaging vectors
        if(io->averaging)
        {

            // allocate memory and set to zero
            VecDuplicate(mesh->Cent,  &(avg->avgU));        VecSet(avg->avgU,0.);
            VecDuplicate(mesh->Nvert, &(avg->avgP));        VecSet(avg->avgP,0.);
            DMCreateGlobalVector(mesh->sda, &(avg->avgUU)); VecSet(avg->avgUU,0.);

            if(io->averaging > 1)
            {
                VecDuplicate(mesh->Cent, &(avg->avgOmega)); VecSet(avg->avgOmega,0.);

                if(io->averaging > 2)
                {
                    VecDuplicate(mesh->Nvert, &(avg->avgP2));        VecSet(avg->avgP2,0.);
                    VecDuplicate(mesh->Nvert, &(avg->avgUdotGradP)); VecSet(avg->avgUdotGradP,0.);
                    VecDuplicate(mesh->Nvert, &(avg->avgMagGradU));  VecSet(avg->avgMagGradU,0.);
                    VecDuplicate(mesh->Cent,  &(avg->avgMagUU));     VecSet(avg->avgMagUU,0.);

                    DMCreateGlobalVector(mesh->sda, &(avg->avgOmegaOmega));  VecSet(avg->avgOmegaOmega,0.);
                }
            }

            if(flags->isLesActive)
            {
                VecDuplicate(mesh->Nvert, &(avg->avgNut)); VecSet(avg->avgNut, 0.);
                VecDuplicate(mesh->Nvert, &(avg->avgCs));  VecSet(avg->avgCs, 0.);
            }

        }

        // allocate phase averaging vectors
        if(io->phaseAveraging)
        {

            // allocate memory and set to zero
            VecDuplicate(mesh->Cent,  &(avg->pAvgU));        VecSet(avg->pAvgU,0.);
            VecDuplicate(mesh->Nvert, &(avg->pAvgP));        VecSet(avg->pAvgP,0.);
            DMCreateGlobalVector(mesh->sda, &(avg->pAvgUU)); VecSet(avg->pAvgUU,0.);

            if(io->phaseAveraging > 1)
            {
                VecDuplicate(mesh->Cent, &(avg->pAvgOmega)); VecSet(avg->pAvgOmega,0.);

                if(io->phaseAveraging > 2)
                {
                    VecDuplicate(mesh->Nvert, &(avg->pAvgP2));        VecSet(avg->pAvgP2,0.);
                    VecDuplicate(mesh->Nvert, &(avg->pAvgUdotGradP)); VecSet(avg->pAvgUdotGradP,0.);
                    VecDuplicate(mesh->Nvert, &(avg->pAvgMagGradU));  VecSet(avg->pAvgMagGradU,0.);
                    VecDuplicate(mesh->Cent,  &(avg->pAvgMagUU));     VecSet(avg->pAvgMagUU,0.);

                    DMCreateGlobalVector(mesh->sda, &(avg->pAvgOmegaOmega));  VecSet(avg->pAvgOmegaOmega,0.);
                }
            }

            if(flags->isLesActive)
            {
                VecDuplicate(mesh->Nvert, &(avg->pAvgNut)); VecSet(avg->pAvgNut, 0.);
                VecDuplicate(mesh->Nvert, &(avg->pAvgCs));  VecSet(avg->pAvgCs, 0.);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode averageFields(acquisition_ *acquisition)
{
    io_    *io    = acquisition->access->io;
    clock_ *clock  = acquisition->access->clock;

    if(io->averaging || io->phaseAveraging)
    {
        // accumulation flags for current time step
        PetscInt    accumulate         = 0;
        PetscInt    accumulateAvg      = 0;
        PetscInt    accumulatePhaseAvg = 0;
        PetscReal   epsilon            = 1e-8;

        if(io->averaging)
        {
            PetscReal startTimeAvg         = io->avgStartTime;
            PetscReal timeIntervalAvg      = io->avgPrd;
            // check if must accumulate averaged fields
            if
            (
                clock->time >= startTimeAvg &&
                (clock->time - startTimeAvg ) / timeIntervalAvg -
                std::floor((clock->time - startTimeAvg) / timeIntervalAvg + epsilon) < 1e-10
            )
            {
                accumulateAvg = 1;
            }
        }

        if(io->phaseAveraging)
        {
            PetscReal startTimePhaseAvg    = io->phAvgStartTime;
            PetscReal timeIntervalPhaseAvg = io->phAvgPrd;

            // check if must accumulate phase averaged fields
            if
            (
                clock->time >= startTimePhaseAvg &&
                (clock->time - startTimePhaseAvg ) / timeIntervalPhaseAvg -
                std::floor((clock->time - startTimePhaseAvg) / timeIntervalPhaseAvg + epsilon) < 1e-10
            )
            {
                accumulatePhaseAvg = 1;
            }
        }

        accumulate = accumulateAvg + accumulatePhaseAvg;

        if(accumulate)
        {
            mesh_  *mesh   = acquisition->access->mesh;
            flags_ *flags  = acquisition->access->flags;

            ueqn_  *ueqn   = acquisition->access->ueqn;
            peqn_  *peqn   = acquisition->access->peqn;
            les_   *les    = NULL;

            avgFields *avg = acquisition->fields;

            DMDALocalInfo info = mesh->info;
            DM            da = mesh->da, fda = mesh->fda, sda = mesh->sda;

            PetscInt       xs = info.xs, xe = info.xs + info.xm;
            PetscInt       ys = info.ys, ye = info.ys + info.ym;
            PetscInt       zs = info.zs, ze = info.zs + info.zm;
            PetscInt       mx = info.mx, my = info.my, mz = info.mz;

            PetscInt       i, j, k;
            PetscInt       lxs, lxe, lys, lye, lzs, lze;

            Cmpnts    ***ucat, ***csi, ***eta, ***zet;
            PetscReal ***p, ***nut, ***cs, ***nvert, ***aj;

            Cmpnts    ***u_mean, ***maguu_mean, ***w_mean,
                      ***u_phase, ***maguu_phase, ***w_phase;

            PetscReal ***p_mean, ***p2_mean, ***udgp_mean, ***maggradu_mean, ***nut_mean, ***cs_mean,
                      ***p_phase, ***p2_phase, ***udgp_phase, ***maggradu_phase, ***nut_phase, ***cs_phase;

            symmTensor ***uu_mean, ***uu_phase, ***ww_mean, ***ww_phase;

            PetscReal     ts, te;

            // averaging weights
            PetscReal        aN, pN;
            PetscReal        m1, m2, p1, p2;

            PetscTime(&ts);

            // indices for internal cells
            lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
            lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
            lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

            // get solution arrays
            DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
            DMDAVecGetArray(da,  mesh->lAj,    &aj);
            DMDAVecGetArray(da,  mesh->lNvert, &nvert);
            DMDAVecGetArray(fda, mesh->lCsi,   &csi);
            DMDAVecGetArray(fda, mesh->lEta,   &eta);
            DMDAVecGetArray(fda, mesh->lZet,   &zet);
            DMDAVecGetArray(da,  peqn->lP,     &p);

            if(flags->isLesActive)
            {
                les = acquisition->access->les;

                DMDAVecGetArray(da, les->lNu_t, &nut);
                DMDAVecGetArray(da, les->lCs, &cs);
            }

            // get averaged arrays
            if(accumulateAvg)
            {
                aN = (PetscReal)io->avgWeight;
                m1 = aN  / (aN + 1.0);
                m2 = 1.0 / (aN + 1.0);

                DMDAVecGetArray(fda, avg->avgU,  &u_mean);
                DMDAVecGetArray(da,  avg->avgP,  &p_mean);
                DMDAVecGetArray(sda, avg->avgUU, &uu_mean);

                if(io->averaging > 1)
                {
                    DMDAVecGetArray(fda, avg->avgOmega, &w_mean);

                    if(io->averaging > 2)
                    {
                        DMDAVecGetArray(da,  avg->avgP2,        &p2_mean);
                        DMDAVecGetArray(da,  avg->avgUdotGradP, &udgp_mean);
                        DMDAVecGetArray(da,  avg->avgMagGradU,  &maggradu_mean);
                        DMDAVecGetArray(fda, avg->avgMagUU,     &maguu_mean);

                        DMDAVecGetArray(sda, avg->avgOmegaOmega, &ww_mean);
                    }
                }

                if(flags->isLesActive)
                {
                    DMDAVecGetArray(da, avg->avgNut, &nut_mean);
                    DMDAVecGetArray(da, avg->avgCs,  &cs_mean);
                }
            }

            // get phase averaged arrays
            if (accumulatePhaseAvg)
            {
                pN = (PetscReal)io->pAvgWeight;
                p1 = pN  / (pN + 1.0);
                p2 = 1.0 / (pN + 1.0);

                DMDAVecGetArray(fda, avg->pAvgU,  &u_phase);
                DMDAVecGetArray(da,  avg->pAvgP,  &p_phase);
                DMDAVecGetArray(sda, avg->pAvgUU, &uu_phase);

                if(io->phaseAveraging > 1)
                {
                    DMDAVecGetArray(fda, avg->pAvgOmega, &w_phase);

                    if(io->phaseAveraging > 2)
                    {
                        DMDAVecGetArray(da,  avg->pAvgP2,        &p2_phase);
                        DMDAVecGetArray(da,  avg->pAvgUdotGradP, &udgp_phase);
                        DMDAVecGetArray(da,  avg->pAvgMagGradU,  &maggradu_phase);
                        DMDAVecGetArray(fda, avg->pAvgMagUU,     &maguu_phase);

                        DMDAVecGetArray(sda, avg->pAvgOmegaOmega, &ww_phase);
                    }
                }

                if(flags->isLesActive)
                {
                    DMDAVecGetArray(da, avg->pAvgNut, &nut_phase);
                    DMDAVecGetArray(da, avg->pAvgCs,  &cs_phase);
                }
            }

            for (k = lzs; k < lze; k++)
            {
                for (j = lys; j < lye; j++)
                {
                    for (i = lxs; i < lxe; i++)
                    {
                        // pre-set base variables for speed
                        PetscReal U = ucat[k][j][i].x,
                               V = ucat[k][j][i].y,
                               W = ucat[k][j][i].z,
                               P = p[k][j][i];

                        PetscReal dudc, dvdc, dwdc,
                               dude, dvde, dwde,
                               dudz, dvdz, dwdz;
                        PetscReal du_dx, du_dy, du_dz,
                               dv_dx, dv_dy, dv_dz,
                               dw_dx, dw_dy, dw_dz;
                        PetscReal dpdc, dpde, dpdz;
                        PetscReal dp_dx, dp_dy, dp_dz;

                        PetscReal csi0 = csi[k][j][i].x,
                               csi1 = csi[k][j][i].y,
                               csi2 = csi[k][j][i].z;
                        PetscReal eta0 = eta[k][j][i].x,
                               eta1 = eta[k][j][i].y,
                               eta2 = eta[k][j][i].z;
                        PetscReal zet0 = zet[k][j][i].x,
                               zet1 = zet[k][j][i].y,
                               zet2 = zet[k][j][i].z;
                        PetscReal ajc  = aj[k][j][i];

                        Compute_du_center
                        (
                            mesh,
                            i, j, k, mx, my, mz, ucat, nvert, &dudc,
                            &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                        );

                        Compute_dscalar_center
                        (
                            mesh,
                            i, j, k, mx, my, mz, p, nvert, &dpdc, &dpde, &dpdz
                        );

                        Compute_du_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                            zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                            dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                            &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                        );

                        Compute_dscalar_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2,
                            zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx,
                            &dp_dy, &dp_dz
                        );

                        PetscReal omega_x = dw_dy - dv_dz,
                               omega_y = du_dz - dw_dx,
                               omega_z = dv_dx - du_dy;

                        // cumulate the averaged fields
                        if(accumulateAvg)
                        {
                            // cumulate avgU
                            u_mean[k][j][i].x = m1 * u_mean[k][j][i].x + m2 * U;
                            u_mean[k][j][i].y = m1 * u_mean[k][j][i].y + m2 * V;
                            u_mean[k][j][i].z = m1 * u_mean[k][j][i].z + m2 * W;

                            // cumulate avgP
                            p_mean[k][j][i] = m1 * p_mean[k][j][i] + m2 * P;

                            // cumulate avgUU
                            uu_mean[k][j][i].xx = m1 * uu_mean[k][j][i].xx + m2 * U * U;
                            uu_mean[k][j][i].yy = m1 * uu_mean[k][j][i].yy + m2 * V * V;
                            uu_mean[k][j][i].zz = m1 * uu_mean[k][j][i].zz + m2 * W * W;

                            uu_mean[k][j][i].xy = m1 * uu_mean[k][j][i].xy + m2 * U * V;
                            uu_mean[k][j][i].xz = m1 * uu_mean[k][j][i].xz + m2 * U * W;
                            uu_mean[k][j][i].yz = m1 * uu_mean[k][j][i].yz + m2 * V * W;

                            if (flags->isLesActive)
                            {
                                // cumulate avgNut and avgCs
                                nut_mean[k][j][i] = m1 * nut_mean[k][j][i] + m2 * nut[k][j][i];
                                cs_mean[k][j][i]  = m1 * cs_mean[k][j][i]  + m2 * cs[k][j][i];
                            }

                            if(io->averaging > 1)
                            {
                                // cumulate avgOmega
                                w_mean[k][j][i].x = m1 * w_mean[k][j][i].x + m2 * omega_x;
                                w_mean[k][j][i].y = m1 * w_mean[k][j][i].y + m2 * omega_y;
                                w_mean[k][j][i].z = m1 * w_mean[k][j][i].z + m2 * omega_z;

                                if(io->averaging > 2)
                                {
                                    // cumulate avgOmegaOmega
                                    ww_mean[k][j][i].xx = m1 * ww_mean[k][j][i].xx + m2 * omega_x * omega_x;
                                    ww_mean[k][j][i].yy = m1 * ww_mean[k][j][i].yy + m2 * omega_y * omega_y;
                                    ww_mean[k][j][i].zz = m1 * ww_mean[k][j][i].zz + m2 * omega_z * omega_z;

                                    ww_mean[k][j][i].xy = m1 * ww_mean[k][j][i].xy + m2 * omega_x * omega_y;
                                    ww_mean[k][j][i].xz = m1 * ww_mean[k][j][i].xz + m2 * omega_x * omega_z;
                                    ww_mean[k][j][i].yz = m1 * ww_mean[k][j][i].yz + m2 * omega_y * omega_z;

                                    // cumulate avgP2
                                    p2_mean[k][j][i] = m1 * p2_mean[k][j][i] + m2 * P * P;

                                    // cumulate avgUdotGradP
                                    udgp_mean[k][j][i] = m1 * udgp_mean[k][j][i] + m2 * (U * dp_dx + V * dp_dy + W * dp_dz);

                                    // cumulate avgMagGradU
                                    maggradu_mean[k][j][i]
                                    =
                                    m1 * maggradu_mean[k][j][i] +
                                    m2 * sqrt
                                    (
                                        2.0*
                                        (
                                            pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.) +
                                            pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.) +
                                            pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.)
                                        )
                                    );

                                    // cumulate avgMagUU
                                    maguu_mean[k][j][i].x = m1 * maguu_mean[k][j][i].x + m2 * (U * U + V * V + W * W) * U;
                                    maguu_mean[k][j][i].y = m1 * maguu_mean[k][j][i].y + m2 * (U * U + V * V + W * W) * V;
                                    maguu_mean[k][j][i].z = m1 * maguu_mean[k][j][i].z + m2 * (U * U + V * V + W * W) * W;
                                }
                            }
                        }

                        // cumulate the phase averaged fields
                        if(accumulatePhaseAvg)
                        {
                            // cumulate avgU
                            u_phase[k][j][i].x = p1 * u_phase[k][j][i].x + p2 * U;
                            u_phase[k][j][i].y = p1 * u_phase[k][j][i].y + p2 * V;
                            u_phase[k][j][i].z = p1 * u_phase[k][j][i].z + p2 * W;

                            // cumulate avgP
                            p_phase[k][j][i] = p1 * p_phase[k][j][i] + p2 * P;

                            // cumulate pAvgUU
                            uu_phase[k][j][i].xx = p1 * uu_phase[k][j][i].xx + p2 * U * U;
                            uu_phase[k][j][i].yy = p1 * uu_phase[k][j][i].yy + p2 * V * V;
                            uu_phase[k][j][i].zz = p1 * uu_phase[k][j][i].zz + p2 * W * W;

                            uu_phase[k][j][i].xy = p1 * uu_phase[k][j][i].xy + p2 * U * V;
                            uu_phase[k][j][i].xz = p1 * uu_phase[k][j][i].xz + p2 * U * W;
                            uu_phase[k][j][i].yz = p1 * uu_phase[k][j][i].yz + p2 * V * W;

                            if(flags->isLesActive)
                            {
                                // cumulate pAvgNut and pAvgCs
                                nut_phase[k][j][i] = p1 * nut_phase[k][j][i] + p2 * nut[k][j][i];
                                cs_phase[k][j][i]  = p1 * cs_phase[k][j][i]  + p2 * cs[k][j][i];
                            }

                            if(io->phaseAveraging > 1)
                            {
                                // cumulate pAvgOmega
                                w_phase[k][j][i].x = p1 * w_phase[k][j][i].x + p2 * omega_x;
                                w_phase[k][j][i].y = p1 * w_phase[k][j][i].y + p2 * omega_y;
                                w_phase[k][j][i].z = p1 * w_phase[k][j][i].z + p2 * omega_z;

                                if(io->phaseAveraging > 2)
                                {
                                    // cumulate pAvgOmegaOmega
                                    ww_phase[k][j][i].xx = p1 * ww_phase[k][j][i].xx + p2 * omega_x * omega_x;
                                    ww_phase[k][j][i].yy = p1 * ww_phase[k][j][i].yy + p2 * omega_y * omega_y;
                                    ww_phase[k][j][i].zz = p1 * ww_phase[k][j][i].zz + p2 * omega_z * omega_z;

                                    ww_phase[k][j][i].xy = p1 * ww_phase[k][j][i].xy + p2 * omega_x * omega_y;
                                    ww_phase[k][j][i].xz = p1 * ww_phase[k][j][i].xz + p2 * omega_x * omega_z;
                                    ww_phase[k][j][i].yz = p1 * ww_phase[k][j][i].yz + p2 * omega_y * omega_z;

                                    // cumulate pAvgP2
                                    p2_phase[k][j][i] = p1 * p2_phase[k][j][i] + p2 * P * P;

                                    // cumulate pAvgUdotGradP
                                    udgp_phase[k][j][i] = p1 * udgp_phase[k][j][i] + p2 * (U * dp_dx + V * dp_dy + W * dp_dz);

                                    // cumulate pAvgMagGradU
                                    maggradu_phase[k][j][i]
                                    =
                                    p1 * maggradu_phase[k][j][i] +
                                    p2 * sqrt
                                    (
                                        2.0*
                                        (
                                            pow(du_dx, 2.) + pow(du_dy, 2.) + pow(du_dz, 2.) +
                                            pow(dv_dx, 2.) + pow(dv_dy, 2.) + pow(dv_dz, 2.) +
                                            pow(dw_dx, 2.) + pow(dw_dy, 2.) + pow(dw_dz, 2.)
                                        )
                                    );

                                    // cumulate pAvgMagUU
                                    maguu_phase[k][j][i].x = p1 * maguu_phase[k][j][i].x + p2 * (U * U + V * V + W * W) * U;
                                    maguu_phase[k][j][i].y = p1 * maguu_phase[k][j][i].y + p2 * (U * U + V * V + W * W) * V;
                                    maguu_phase[k][j][i].z = p1 * maguu_phase[k][j][i].z + p2 * (U * U + V * V + W * W) * W;
                                }
                            }
                        }
                    }
                }
            }

            // restore solution arrays
            DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
            DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
            DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
            DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
            DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
            DMDAVecRestoreArray(fda, mesh->lZet,   &zet);
            DMDAVecRestoreArray(da,  peqn->lP,     &p);

            if(flags->isLesActive)
            {
                DMDAVecRestoreArray(da, les->lNu_t, &nut);
                DMDAVecRestoreArray(da, les->lCs, &cs);
            }

            // restore averaged arrays
            if(accumulateAvg)
            {
                DMDAVecRestoreArray(fda, avg->avgU,  &u_mean);
                DMDAVecRestoreArray(da,  avg->avgP,  &p_mean);
                DMDAVecRestoreArray(sda, avg->avgUU, &uu_mean);

                if(io->averaging > 1)
                {
                    DMDAVecRestoreArray(fda, avg->avgOmega, &w_mean);

                    if(io->averaging > 2)
                    {
                        DMDAVecRestoreArray(da,  avg->avgP2,        &p2_mean);
                        DMDAVecRestoreArray(da,  avg->avgUdotGradP, &udgp_mean);
                        DMDAVecRestoreArray(da,  avg->avgMagGradU,  &maggradu_mean);
                        DMDAVecRestoreArray(fda, avg->avgMagUU,     &maguu_mean);

                        DMDAVecRestoreArray(sda, avg->avgOmegaOmega, &ww_mean);
                    }
                }

                if(flags->isLesActive)
                {
                    DMDAVecRestoreArray(da, avg->avgNut, &nut_mean);
                    DMDAVecRestoreArray(da, avg->avgCs,  &cs_mean);
                }

                // increase snapshot weighting
                io->avgWeight++;
            }

            // restore phase averaged arrays
            if(accumulatePhaseAvg)
            {
                DMDAVecRestoreArray(fda, avg->pAvgU,  &u_phase);
                DMDAVecRestoreArray(da,  avg->pAvgP,  &p_phase);
                DMDAVecRestoreArray(sda, avg->pAvgUU, &uu_phase);

                if(io->phaseAveraging > 1)
                {
                    DMDAVecRestoreArray(fda, avg->pAvgOmega, &w_phase);

                    if(io->phaseAveraging > 2)
                    {
                        DMDAVecRestoreArray(da,  avg->pAvgP2,        &p2_phase);
                        DMDAVecRestoreArray(da,  avg->pAvgUdotGradP, &udgp_phase);
                        DMDAVecRestoreArray(da,  avg->pAvgMagGradU,  &maggradu_phase);
                        DMDAVecRestoreArray(fda, avg->pAvgMagUU,     &maguu_phase);

                        DMDAVecRestoreArray(sda, avg->pAvgOmegaOmega, &ww_phase);
                    }
                }

                if(flags->isLesActive)
                {
                    DMDAVecRestoreArray(da, avg->pAvgNut, &nut_phase);
                    DMDAVecRestoreArray(da, avg->pAvgCs,  &cs_phase);
                }

                // increase snapshot weighting
                io->pAvgWeight++;
            }

            PetscTime(&te);
            PetscPrintf(mesh->MESH_COMM, "Averaged fields in %lf s\n", te-ts);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode sectionsInitialize(acquisition_ *acquisition)
{
    flags_        *flags = acquisition->access->flags;
    mesh_         *mesh = acquisition->access->mesh;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***cent;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;
    PetscMPIInt        rank, nprocs;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    MPI_Comm_size(mesh->MESH_COMM, &nprocs);

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation for search (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    // initialize sections to save and create directories
    if(acquisition->isSectionsActive)
    {
        // surface location is specified with the index of the generalized direction
        // normal to the plane. E.g. the surface normal to the k-direction storing the
        // velocity at the k-left face is specified with 1.
        // Files must be located inside "sampling/surfaces/kSections" for k-normal,
        // "sampling/surfaces/jSections" for j-normal, "sampling/surfaces/iSections" for i-normal.

        PetscPrintf(mesh->MESH_COMM, "Reading surfaces in sampling/surfaces/...");

        DMDAVecGetArray(mesh->fda, mesh->lCent, &cent);

        PetscInt iskavail = 0, isjavail = 0, isiavail = 0;

        word dataLoc, kSecName, jSecName, iSecName;

        dataLoc = "sampling/surfaces/";

        kSecName = dataLoc + "kSections";
        jSecName = dataLoc + "jSections";
        iSecName = dataLoc + "iSections";

        // see which sections are available
        FILE *fp;

        fp = fopen(kSecName.c_str(), "r");
        if(fp!=NULL)
        {
            iskavail = 1;
            fclose(fp);
        }

        fp = fopen(jSecName.c_str(), "r");
        if(fp!=NULL)
        {
            isjavail = 1;
            fclose(fp);
        }

        fp = fopen(iSecName.c_str(), "r");
        if(fp!=NULL)
        {
            isiavail = 1;
            fclose(fp);
        }

        // read k-section input file
        if(iskavail)
        {
            PetscInt atLeastOneVector = 0;
            PetscInt atLeastOneScalar = 0;

            // read number of sections
            PetscInt nkSec;
            readDictInt(kSecName.c_str(), "surfaceNumber" ,&nkSec);

            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->kSections));
            PetscMalloc(sizeof(PetscInt) * nkSec,  &(acquisition->kSections->indices));
            PetscMalloc(sizeof(PetscReal) * nkSec,  &(acquisition->kSections->coordinates));

            sections *kSections = acquisition->kSections;

            // store number of sections
            kSections->nSections = nkSec;

            // set available to 1
            kSections->available = iskavail;

            // read acquisition start time and type of interval
            readDictDouble(kSecName.c_str(), "timeStart", &(kSections->timeStart));
            readDictWord(kSecName.c_str(),   "intervalType", &(kSections->intervalType));
            readDictDouble(kSecName.c_str(), "timeInterval", &(kSections->timeInterval));

            // check if intervalType is known
            if(kSections->intervalType != "timeStep" && kSections->intervalType != "adjustableTime")
            {
               char error[512];
                sprintf(error, "unknown interval type %s. Known types are timeStep and adjustableTime\n", kSections->intervalType.c_str());
                fatalErrorInFunction("sectionsInitialize",  error);
            }

            // read section indices
            std::ifstream indata;

            indata.open(kSecName.c_str());

            char buffer[256];

            while(!indata.eof())
            {
                indata >> buffer;

                if
                (
                    strcmp
                    (
                        "coordinates",
                        buffer
                    ) == 0
                )
                {
                    for(PetscInt s=0; s<nkSec; s++)
                    {
                        indata >> buffer;
                        std::sscanf(buffer, "%lf", &(kSections->coordinates[s]));
                    }
                }
            }
            indata.close();

            // find the closest point, then take is k index
            for(PetscInt s=0; s<nkSec; s++)
            {
                Cmpnts surfacePoint;
                       surfacePoint.x = kSections->coordinates[s];
                       surfacePoint.y = (mesh->bounds.ymax + mesh->bounds.ymin) / 2.0;
                       surfacePoint.z = (mesh->bounds.zmax + mesh->bounds.zmin) / 2.0;

                PetscReal  lminDist = 1e30,  gminDist = 1e30;
                cellIds closestIds;
                PetscInt     lclosestK = 0;

                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            Cmpnts distVec = nSub(surfacePoint, cent[k][j][i]);
                            PetscReal distMag = nMag(distVec) + procContrib;

                            if(distMag < lminDist)
                            {
                                lminDist = distMag;
                                closestIds.i = i;
                                closestIds.j = j;
                                closestIds.k = k;
                            }

                        }
                    }
                }

                MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

                if(lminDist == gminDist)
                {
                    lclosestK = closestIds.k;
                }

                MPI_Allreduce(&lclosestK, &(kSections->indices[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
            }

            // make directories for saving
            errno = 0;

            // kSurfaces directory
            word kslicesFolder = "./postProcessing/" + mesh->meshName + "/kSurfaces";
            PetscInt dirRes = mkdir(kslicesFolder.c_str(), 0777);
            if (dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", kslicesFolder.c_str());
                fatalErrorInFunction("sectionsInitialize", error);
            }
            else
            {
                // each k-slice directory
                for(PetscInt s=0; s<kSections->nSections; s++)
                {
                    char ksliceName[256];
                    sprintf(ksliceName, "%s/%ld", kslicesFolder.c_str(), kSections->indices[s]);

                    errno = 0;
                    PetscInt dirRes = mkdir(ksliceName, 0777);
                    if (dirRes != 0 && errno != EEXIST)
                    {
                       char error[512];
                        sprintf(error, "could not create %s directory\n", ksliceName);
                        fatalErrorInFunction("sectionsInitialize", error);
                    }
                    else
                    {
                        // create U directory in which time snapshots are saved
                        PetscInt dirRes;
                        char ksliceNameU[260];
                        sprintf(ksliceNameU, "%s/U", ksliceName);

                        errno = 0;
                        dirRes = mkdir(ksliceNameU, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", ksliceNameU);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, ksliceNameU);
                            atLeastOneVector++;
                        }

                        char ksliceNameP[260];
                        sprintf(ksliceNameP, "%s/p", ksliceName);

                        errno = 0;
                        dirRes = mkdir(ksliceNameP, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", ksliceNameP);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, ksliceNameP);
                            atLeastOneScalar++;
                        }

                        // create T directory in which time snapshots are saved
                        if(flags->isTeqnActive)
                        {
                            char ksliceNameT[260];
                            sprintf(ksliceNameT, "%s/T", ksliceName);

                            errno = 0;
                            dirRes = mkdir(ksliceNameT, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", ksliceNameT);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, ksliceNameT);
                                atLeastOneScalar++;
                            }
                        }

                        // create nut directory in which time snapshots are saved
                        if(flags->isLesActive)
                        {
                            char ksliceNameNut[260];
                            sprintf(ksliceNameNut, "%s/nut", ksliceName);

                            errno = 0;
                            dirRes = mkdir(ksliceNameNut, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", ksliceNameNut);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, ksliceNameNut);
                                atLeastOneScalar++;
                            }
                        }
                    }
                }
            }

            // allocate variables where data are stored
            if(atLeastOneVector)
            {
                kSections->vectorSec = (Cmpnts **)malloc( sizeof(Cmpnts *) * my );

                for(j=0; j<my; j++)
                {
                    kSections->vectorSec[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * mx );
                }
            }

            if(atLeastOneScalar)
            {
                kSections->scalarSec = (PetscReal **)malloc( sizeof(PetscReal *) * my );

                for(j=0; j<my; j++)
                {
                    kSections->scalarSec[j] = (PetscReal *)malloc( sizeof(PetscReal) * mx );
                }
            }
        }
        else
        {
            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->kSections));

            // set available to 0
            acquisition->kSections->available = 0;
        }

        // read j-section input file
        if(isjavail)
        {
            PetscInt atLeastOneVector = 0;
            PetscInt atLeastOneScalar = 0;

            // read number of sections
            PetscInt njSec;
            readDictInt(jSecName.c_str(), "surfaceNumber" ,&njSec);

            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->jSections));
            PetscMalloc(sizeof(PetscInt) * njSec,  &(acquisition->jSections->indices));
            PetscMalloc(sizeof(PetscReal) * njSec,  &(acquisition->jSections->coordinates));

            sections *jSections = acquisition->jSections;

            // store number of sections
            jSections->nSections = njSec;

            // set available to 1
            jSections->available = isjavail;

            // read acquisition start time and type of interval
            readDictDouble(jSecName.c_str(), "timeStart", &(jSections->timeStart));
            readDictWord(jSecName.c_str(),   "intervalType", &(jSections->intervalType));
            readDictDouble(jSecName.c_str(), "timeInterval", &(jSections->timeInterval));

            // check if intervalType is known
            if(jSections->intervalType != "timeStep" && jSections->intervalType != "adjustableTime")
            {
               char error[512];
                sprintf(error, "unknown interval type %s. Known types are timeStep and adjustableTime\n", jSections->intervalType.c_str());
                fatalErrorInFunction("sectionsInitialize",  error);
            }

            // read section indices
            std::ifstream indata;

            indata.open(jSecName.c_str());

            char buffer[256];

            while(!indata.eof())
            {
                indata >> buffer;

                if
                (
                    strcmp
                    (
                        "coordinates",
                        buffer
                    ) == 0
                )
                {
                    for(PetscInt s=0; s<njSec; s++)
                    {
                        indata >> buffer;
                        std::sscanf(buffer, "%lf", &(jSections->coordinates[s]));
                    }
                }
            }
            indata.close();

            // find the closest point, then take is k index
            for(PetscInt s=0; s<njSec; s++)
            {
                Cmpnts surfacePoint;
                       surfacePoint.x = (mesh->bounds.xmax + mesh->bounds.xmin) / 2.0;
                       surfacePoint.y = (mesh->bounds.ymax + mesh->bounds.ymin) / 2.0;
                       surfacePoint.z = jSections->coordinates[s];

                PetscReal  lminDist = 1e30,  gminDist = 1e30;
                cellIds closestIds;
                PetscInt     lclosestJ = 0;

                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            Cmpnts distVec = nSub(surfacePoint, cent[k][j][i]);
                            PetscReal distMag = nMag(distVec) + procContrib;

                            if(distMag < lminDist)
                            {
                                lminDist = distMag;
                                closestIds.i = i;
                                closestIds.j = j;
                                closestIds.k = k;
                            }

                        }
                    }
                }

                MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

                if(lminDist == gminDist)
                {
                    lclosestJ = closestIds.j;
                }

                MPI_Allreduce(&lclosestJ, &(jSections->indices[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
            }

            // make directories for saving
            errno = 0;


            // kSurfaces directory
            word jslicesFolder = "./postProcessing/" + mesh->meshName + "/jSurfaces";
            PetscInt dirRes = mkdir(jslicesFolder.c_str(), 0777);
            if (dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", jslicesFolder.c_str());
                fatalErrorInFunction("sectionsInitialize", error);
            }
            else
            {
                // each k-slice directory
                for(PetscInt s=0; s<jSections->nSections; s++)
                {
                    char jsliceName[256];
                    sprintf(jsliceName, "%s/%ld", jslicesFolder.c_str(), jSections->indices[s]);

                    errno = 0;
                    PetscInt dirRes = mkdir(jsliceName, 0777);
                    if (dirRes != 0 && errno != EEXIST)
                    {
                       char error[512];
                        sprintf(error, "could not create %s directory\n", jsliceName);
                        fatalErrorInFunction("sectionsInitialize", error);
                    }
                    else
                    {
                        // create U directory in which time snapshots are saved
                        PetscInt dirRes;
                        char jsliceNameU[260];
                        sprintf(jsliceNameU, "%s/U", jsliceName);

                        errno = 0;
                        dirRes = mkdir(jsliceNameU, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", jsliceNameU);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, jsliceNameU);
                            atLeastOneVector++;
                        }

                        char jsliceNameP[260];
                        sprintf(jsliceNameP, "%s/p", jsliceName);

                        errno = 0;
                        dirRes = mkdir(jsliceNameP, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", jsliceNameP);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, jsliceNameP);
                            atLeastOneScalar++;
                        }

                        // create T directory in which time snapshots are saved
                        if(flags->isTeqnActive)
                        {
                            char jsliceNameT[260];
                            sprintf(jsliceNameT, "%s/T", jsliceName);

                            errno = 0;
                            dirRes = mkdir(jsliceNameT, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", jsliceNameT);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, jsliceNameT);
                                atLeastOneScalar++;
                            }
                        }

                        // create nut directory in which time snapshots are saved
                        if(flags->isLesActive)
                        {
                            char jsliceNameNut[260];
                            sprintf(jsliceNameNut, "%s/nut", jsliceName);

                            errno = 0;
                            dirRes = mkdir(jsliceNameNut, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", jsliceNameNut);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, jsliceNameNut);
                                atLeastOneScalar++;
                            }
                        }
                    }
                }
            }

            // allocate variables where data are stored
            if(atLeastOneVector)
            {
                jSections->vectorSec = (Cmpnts **)malloc( sizeof(Cmpnts *) * mz );

                for(k=0; k<mz; k++)
                {
                    jSections->vectorSec[k] = (Cmpnts *)malloc( sizeof(Cmpnts) * mx );
                }
            }

            if(atLeastOneScalar)
            {
                jSections->scalarSec = (PetscReal **)malloc( sizeof(PetscReal *) * mz );

                for(k=0; k<mz; k++)
                {
                    jSections->scalarSec[k] = (PetscReal *)malloc( sizeof(PetscReal) * mx );
                }
            }
        }
        else
        {
            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->jSections));

            // set available to 0
            acquisition->jSections->available = 0;
        }

        // read i-section input file
        if(isiavail)
        {
            PetscInt atLeastOneVector = 0;
            PetscInt atLeastOneScalar = 0;

            // read number of sections
            PetscInt niSec;
            readDictInt(iSecName.c_str(), "surfaceNumber" ,&niSec);

            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->iSections));
            PetscMalloc(sizeof(PetscInt) * niSec,  &(acquisition->iSections->indices));
            PetscMalloc(sizeof(PetscReal) * niSec,  &(acquisition->iSections->coordinates));

            sections *iSections = acquisition->iSections;

            // store number of sections
            iSections->nSections = niSec;

            // set available to 1
            iSections->available = isiavail;

            // read acquisition start time and type of interval
            readDictDouble(iSecName.c_str(), "timeStart", &(iSections->timeStart));
            readDictWord(iSecName.c_str(),   "intervalType", &(iSections->intervalType));
            readDictDouble(iSecName.c_str(), "timeInterval", &(iSections->timeInterval));

            // check if intervalType is known
            if(iSections->intervalType != "timeStep" && iSections->intervalType != "adjustableTime")
            {
               char error[512];
                sprintf(error, "unknown interval type %s. Known types are timeStep and adjustableTime\n", iSections->intervalType.c_str());
                fatalErrorInFunction("sectionsInitialize",  error);
            }

            // read section indices
            std::ifstream indata;

            indata.open(iSecName.c_str());

            char buffer[256];

            while(!indata.eof())
            {
                indata >> buffer;

                if
                (
                    strcmp
                    (
                        "coordinates",
                        buffer
                    ) == 0
                )
                {
                    for(PetscInt s=0; s<niSec; s++)
                    {
                        indata >> buffer;
                        std::sscanf(buffer, "%lf", &(iSections->coordinates[s]));
                    }
                }
            }
            indata.close();

            // find the closest point, then take is k index
            for(PetscInt s=0; s<niSec; s++)
            {
                Cmpnts surfacePoint;
                       surfacePoint.x = (mesh->bounds.xmax + mesh->bounds.xmin) / 2.0;
                       surfacePoint.y = iSections->coordinates[s];
                       surfacePoint.z = (mesh->bounds.zmax + mesh->bounds.zmin) / 2.0;

                PetscReal  lminDist = 1e30,  gminDist = 1e30;
                cellIds closestIds;
                PetscInt     lclosestI = 0;

                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            Cmpnts distVec = nSub(surfacePoint, cent[k][j][i]);
                            PetscReal distMag = nMag(distVec) + procContrib;

                            if(distMag < lminDist)
                            {
                                lminDist = distMag;
                                closestIds.i = i;
                                closestIds.j = j;
                                closestIds.k = k;
                            }

                        }
                    }
                }

                MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

                if(lminDist == gminDist)
                {
                    lclosestI = closestIds.i;
                }

                MPI_Allreduce(&lclosestI, &(iSections->indices[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
            }

            // make directories for saving
            errno = 0;

            // kSurfaces directory
            word islicesFolder = "./postProcessing/" + mesh->meshName + "/iSurfaces";
            PetscInt dirRes = mkdir(islicesFolder.c_str(), 0777);
            if (dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", islicesFolder.c_str());
                fatalErrorInFunction("sectionsInitialize", error);
            }
            else
            {
                // each k-slice directory
                for(PetscInt s=0; s<iSections->nSections; s++)
                {
                    char isliceName[256];
                    sprintf(isliceName, "%s/%ld", islicesFolder.c_str(), iSections->indices[s]);

                    errno = 0;
                    PetscInt dirRes = mkdir(isliceName, 0777);
                    if (dirRes != 0 && errno != EEXIST)
                    {
                       char error[512];
                        sprintf(error, "could not create %s directory\n", isliceName);
                        fatalErrorInFunction("sectionsInitialize", error);
                    }
                    else
                    {
                        // create U directory in which time snapshots are saved
                        PetscInt dirRes;
                        char isliceNameU[260];
                        sprintf(isliceNameU, "%s/U", isliceName);

                        errno = 0;
                        dirRes = mkdir(isliceNameU, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", isliceNameU);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, isliceNameU);
                            atLeastOneVector++;
                        }

                        char isliceNameP[260];
                        sprintf(isliceNameP, "%s/p", isliceName);

                        errno = 0;
                        dirRes = mkdir(isliceNameP, 0777);
                        if (dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", isliceNameP);
                            fatalErrorInFunction("sectionsInitialize", error);
                        }
                        else
                        {
                            //remove_subdirs(mesh->MESH_COMM, isliceNameP);
                            atLeastOneScalar++;
                        }

                        // create T directory in which time snapshots are saved
                        if(flags->isTeqnActive)
                        {
                            char isliceNameT[260];
                            sprintf(isliceNameT, "%s/T", isliceName);

                            errno = 0;
                            dirRes = mkdir(isliceNameT, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", isliceNameT);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, isliceNameT);
                                atLeastOneScalar++;
                            }
                        }

                        // create nut directory in which time snapshots are saved
                        if(flags->isLesActive)
                        {
                            char isliceNameNut[260];
                            sprintf(isliceNameNut, "%s/nut", isliceName);

                            errno = 0;
                            dirRes = mkdir(isliceNameNut, 0777);
                            if (dirRes != 0 && errno != EEXIST)
                            {
                               char error[512];
                                sprintf(error, "could not create %s directory\n", isliceNameNut);
                                fatalErrorInFunction("sectionsInitialize", error);
                            }
                            else
                            {
                                //remove_subdirs(mesh->MESH_COMM, isliceNameNut);
                                atLeastOneScalar++;
                            }
                        }
                    }
                }
            }

            // allocate variables where data are stored
            if(atLeastOneVector)
            {
                iSections->vectorSec = (Cmpnts **)malloc( sizeof(Cmpnts *) * mz );

                for(k=0; k<mz; k++)
                {
                    iSections->vectorSec[k] = (Cmpnts *)malloc( sizeof(Cmpnts) * my );
                }
            }

            if(atLeastOneScalar)
            {
                iSections->scalarSec = (PetscReal **)malloc( sizeof(PetscReal *) * mz );

                for(k=0; k<mz; k++)
                {
                    iSections->scalarSec[k] = (PetscReal *)malloc( sizeof(PetscReal) * my );
                }
            }
        }
        else
        {
            // allocate memory
            PetscMalloc(sizeof(sections), &(acquisition->iSections));

            // set available to 0
            acquisition->iSections->available = 0;
        }

        DMDAVecRestoreArray(mesh->fda, mesh->lCent, &cent);

        PetscPrintf(mesh->MESH_COMM, "done\n\n");
    }

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode writeSections(acquisition_ *acquisition)
{
    if(acquisition->isSectionsActive)
    {
        // get pointers
        clock_ *clock = acquisition->access->clock;
        mesh_  *mesh  = acquisition->access->mesh;
        flags_ *flags = acquisition->access->flags;

        ueqn_  *ueqn  = acquisition->access->ueqn;
        peqn_  *peqn  = acquisition->access->peqn;
        teqn_  *teqn;
        les_   *les;

        if(flags->isTeqnActive) teqn = acquisition->access->teqn;
        if(flags->isLesActive)  les  = acquisition->access->les;

        DMDALocalInfo info = mesh->info;
        PetscInt           mx = info.mx, my = info.my, mz = info.mz;

        PetscMPIInt           rank;
        MPI_Comm_rank(mesh->MESH_COMM, &rank);

        // write i sections to file
        if(acquisition->iSections->available)
        {
            sections *iSections = acquisition->iSections;

            PetscReal timeStart    = iSections->timeStart;
            PetscReal timeInterval = iSections->timeInterval;
            PetscReal epsilon      = 1e-8;
            // check the time and see if must write
            if
            (
                clock->time >= timeStart && // write only if past acquisition start
                (
                    // adjustableTime: write only if clock->time is multiple of time interval
                    (
                        iSections->intervalType == "adjustableTime" &&
                        (clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval + epsilon) < 1e-10
                    ) ||
                    // timeStep: write every timeInterval iterations
                    (
                        iSections->intervalType == "timeStep" &&
                        (clock->it / timeInterval - std::floor(clock->it / timeInterval + epsilon) < 1e-10)
                    )
                )

            )
            {
                PetscPrintf(mesh->MESH_COMM, "Saving i-sections:\n");

                for(PetscInt i=0; i<iSections->nSections; i++)
                {
                    PetscInt iplane = iSections->indices[i];

                    // exclude ghost nodes
                    if (iplane<1 || iplane>mx-2) continue;

                    // store data
                    iSectionSaveVector(mesh, iSections, iplane, ueqn->Ucat, "U");
                    iSectionSaveScalar(mesh, iSections, iplane, peqn->P, "p");

                    if(flags->isTeqnActive) iSectionSaveScalar(mesh, iSections, iplane, teqn->Tmprt, "T");
                    if(flags->isLesActive)  iSectionSaveScalar(mesh, iSections, iplane, les->lNu_t, "nut");
                }
            }
        }

        // write j sections to file
        if(acquisition->jSections->available)
        {
            sections *jSections = acquisition->jSections;

            PetscReal timeStart    = jSections->timeStart;
            PetscReal timeInterval = jSections->timeInterval;
            PetscReal epsilon      = 1e-8;
            // check the time and see if must write
            if
            (
                clock->time >= timeStart && // write only if past acquisition start
                (
                    // adjustableTime: write only if clock->time is multiple of time interval
                    (
                        jSections->intervalType == "adjustableTime" &&
                        (clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval + epsilon) < 1e-10
                    ) ||
                    // timeStep: write every timeInterval iterations
                    (
                        jSections->intervalType == "timeStep" &&
                        (clock->it / timeInterval - std::floor(clock->it / timeInterval + epsilon) < 1e-10)
                    )
                )

            )
            {
                PetscPrintf(mesh->MESH_COMM, "Saving j-sections:\n");

                for(PetscInt j=0; j<jSections->nSections; j++)
                {
                    PetscInt jplane = jSections->indices[j];

                    // exclude ghost nodes
                    if (jplane<1 || jplane>my-2) continue;

                    // store data
                    jSectionSaveVector(mesh, jSections, jplane, ueqn->Ucat, "U");
                    jSectionSaveScalar(mesh, jSections, jplane, peqn->P, "p");

                    if(flags->isTeqnActive) jSectionSaveScalar(mesh, jSections, jplane, teqn->Tmprt, "T");
                    if(flags->isLesActive)  jSectionSaveScalar(mesh, jSections, jplane, les->lNu_t, "nut");
                }
            }
        }

        // write k sections to file
        if(acquisition->kSections->available)
        {
            sections *kSections = acquisition->kSections;

            PetscReal timeStart    = kSections->timeStart;
            PetscReal timeInterval = kSections->timeInterval;
            PetscReal epsilon      = 1e-8;

            // check the time and see if must write
            if
            (
                clock->time >= timeStart && // write only if past acquisition start
                (
                    // adjustableTime: write only if clock->time is multiple of time interval
                    (
                        kSections->intervalType == "adjustableTime" &&
                        (clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval + epsilon) < 1e-10
                    ) ||
                    // timeStep: write every timeInterval iterations
                    (
                        kSections->intervalType == "timeStep" &&
                        (clock->it / timeInterval - std::floor(clock->it / timeInterval + epsilon) < 1e-10)
                    )
                )

            )
            {
                PetscPrintf(mesh->MESH_COMM, "Saving k-sections:\n");

                for(PetscInt k=0; k<kSections->nSections; k++)
                {
                    PetscInt kplane = kSections->indices[k];

                    // exclude ghost nodes
                    if (kplane<1 || kplane>mz-2) continue;

                    // store data
                    kSectionSaveVector(mesh, kSections, kplane, ueqn->Ucat, "U");
                    kSectionSaveScalar(mesh, kSections, kplane, peqn->P, "p");

                    if(flags->isTeqnActive) kSectionSaveScalar(mesh, kSections, kplane, teqn->Tmprt, "T");
                    if(flags->isLesActive)  kSectionSaveScalar(mesh, kSections, kplane, les->lNu_t, "nut");
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt iplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts        ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(k=0; k<mz; k++)
    {
        for(j=0; j<my; j++)
        {
            sec->vectorSec[k][j].x = 0;
            sec->vectorSec[k][j].y = 0;
            sec->vectorSec[k][j].z = 0;
        }
    }

    DMDAVecGetArray(mesh->fda, V, &v);

    if(iplane>=xs && iplane<xe)
    {
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                mSet(sec->vectorSec[k][j], v[k][j][iplane]);
            }
        }
    }

    DMDAVecRestoreArray(mesh->fda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->vectorSec[k], my*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->vectorSec[k], sec->vectorSec[k], my*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/iSurfaces/%ld/%s/%s", mesh->meshName.c_str(), iplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, iplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(k=0; k<mz; k++)
        {
            fwrite(&sec->vectorSec[k][0], sizeof(Cmpnts), my, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt jplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts        ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(k=0; k<mz; k++)
    {
        for(i=0; i<mx; i++)
        {
            sec->vectorSec[k][i].x = 0;
            sec->vectorSec[k][i].y = 0;
            sec->vectorSec[k][i].z = 0;
        }
    }

    DMDAVecGetArray(mesh->fda, V, &v);

    if(jplane>=ys && jplane<ye)
    {
        for (k=zs; k<ze; k++)
        {
            for (i=xs; i<xe; i++)
            {
                mSet(sec->vectorSec[k][i], v[k][jplane][i]);
            }
        }
    }

    DMDAVecRestoreArray(mesh->fda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->vectorSec[k], mx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->vectorSec[k], sec->vectorSec[k], mx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/jSurfaces/%ld/%s/%s", mesh->meshName.c_str(), jplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, jplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(k=0; k<mz; k++)
        {
            fwrite(&sec->vectorSec[k][0], sizeof(Cmpnts), mx, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt kplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts        ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(j=0; j<my; j++)
    {
        for(i=0; i<mx; i++)
        {
            sec->vectorSec[j][i].x = 0;
            sec->vectorSec[j][i].y = 0;
            sec->vectorSec[j][i].z = 0;
        }
    }

    DMDAVecGetArray(mesh->fda, V, &v);

    if(kplane>=zs && kplane<ze)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                mSet(sec->vectorSec[j][i], v[kplane][j][i]);
            }
        }
    }

    DMDAVecRestoreArray(mesh->fda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->vectorSec[j], mx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(sec->vectorSec[j], sec->vectorSec[j], mx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/kSurfaces/%ld/%s/%s", mesh->meshName.c_str(), kplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, kplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(j=0; j<my; j++)
        {
            fwrite(&sec->vectorSec[j][0], sizeof(Cmpnts), mx, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt iplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal     ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(k=0; k<mz; k++)
    {
        for(j=0; j<my; j++)
        {
            sec->scalarSec[k][j] = 0;
        }
    }

    DMDAVecGetArray(mesh->da, V, &v);

    if(iplane>=xs && iplane<xe)
    {
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                sec->scalarSec[k][j] = v[k][j][iplane];
            }
        }
    }

    DMDAVecRestoreArray(mesh->da, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if (!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->scalarSec[k], my, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->scalarSec[k], sec->scalarSec[k], my, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/iSurfaces/%ld/%s/%s",  mesh->meshName.c_str(), iplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, iplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(k=0; k<mz; k++)
        {
            fwrite(&sec->scalarSec[k][0], sizeof(PetscReal), my, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt jplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal     ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(k=0; k<mz; k++)
    {
        for(i=0; i<mx; i++)
        {
            sec->scalarSec[k][i] = 0;
        }
    }

    DMDAVecGetArray(mesh->da, V, &v);

    if(jplane>=ys && jplane<ye)
    {
        for (k=zs; k<ze; k++)
        {
            for (i=xs; i<xe; i++)
            {
                sec->scalarSec[k][i] = v[k][jplane][i];
            }
        }
    }

    DMDAVecRestoreArray(mesh->da, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->scalarSec[k], mx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->scalarSec[k], sec->scalarSec[k], mx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/jSurfaces/%ld/%s/%s", mesh->meshName.c_str(), jplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, jplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(k=0; k<mz; k++)
        {
            fwrite(&sec->scalarSec[k][0], sizeof(PetscReal), mx, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt kplane, Vec &V, const char* fieldName)
{
    clock_        *clock = mesh->access->clock;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal     ***v;

    // formatted print width
    PetscInt width1 = -5;
    PetscInt width2 = -50;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set velocity values to zero in all processes
    for(j=0; j<my; j++)
    {
        for(i=0; i<mx; i++)
        {
            sec->scalarSec[j][i] = 0;
        }
    }

    DMDAVecGetArray(mesh->da, V, &v);

    if(kplane>=zs && kplane<ze)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                sec->scalarSec[j][i] = v[kplane][j][i];
            }
        }
    }

    DMDAVecRestoreArray(mesh->da, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->scalarSec[j], mx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(sec->scalarSec[j], sec->scalarSec[j], mx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    if(!rank)
    {
        word timeName = getTimeName(clock);

        char fname[256];
        sprintf(fname, "./postProcessing/%s/kSurfaces/%ld/%s/%s", mesh->meshName.c_str(), kplane, fieldName, timeName.c_str());

        PetscPrintf(mesh->MESH_COMM, "    %*d/%s to %*s\n", width1, kplane, fieldName, width2, fname);

        FILE *fp=fopen(fname, "wb");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s", fname);
            fatalErrorInFunction("save_inflow_section",  error);
        }

        for(j=0; j<my; j++)
        {
            fwrite(&sec->scalarSec[j][0], sizeof(PetscReal), mx, fp);
        }
        fclose(fp);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ProbesInitialize(domain_ *domain)
{
    // indices for domain, rake and probes
    PetscInt d, r, p;

    // get number of domains
    PetscInt nDomains = domain[0].info.nDomains;

    // get clock anf flags from first domain (they are all the same)
    clock_        *clock  = domain[0].clock;
    flags_        *flags  = domain[0].access.flags;

    // get probes flag from first domain
    PetscInt isProbesActive = domain[0].acquisition->isProbesActive;

    if(isProbesActive)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Creating probes acquisition system...\n\n");

        // locally initialize probe rakes
        rakes         *probes = new rakes;

        PetscInt           nRakes;
        char          dataLoc[256], fileName[500];
        PetscInt           readError;

        // read inside sampling/probes and determine how many rakes are present
        sprintf(dataLoc, "sampling/probes/");
        nRakes = count_files(dataLoc);

        // allocate top level data
        probes->rakes  = new probeRake[nRakes];
        probes->nRakes = nRakes;

        // initialize rake counter to zero
        nRakes = 0;

        DIR *dir; struct dirent *diread;

        // loop through the rakes
        if((dir = opendir(dataLoc)) != nullptr)
        {
            while ((diread = readdir(dir)) != nullptr)
            {
                // get rake name
                char* rakeName = diread->d_name;
                if
                (
                    strcmp(rakeName, ".")  !=0 &&
                    strcmp(rakeName, "..") !=0
                )
                {
                    // allocate memory
                    char tmp[256];
                    char fields[256];

                    sprintf(fileName, "sampling/probes/%s", rakeName);
                    FILE *f = fopen(fileName, "r");

                    if(!f)
                    {
                       char error[530];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("ProbesInitialize",  error);
                    }

                    // get number of probes in the rake
                    readError = fscanf(f, "%s %ld\n", tmp, &(probes->rakes[nRakes].probesNumber));

                    // set the total number of probes in rakes
                    PetscMalloc(sizeof(PetscReal*) * probes->rakes[nRakes].probesNumber, &(probes->rakes[nRakes].locations));
                    PetscMalloc(sizeof(PetscInt*)    * probes->rakes[nRakes].probesNumber, &(probes->rakes[nRakes].cells));
                    PetscMalloc(sizeof(PetscInt)     * probes->rakes[nRakes].probesNumber, &(probes->rakes[nRakes].owner));
                    PetscMalloc(sizeof(PetscInt)     * probes->rakes[nRakes].probesNumber, &(probes->rakes[nRakes].domainID));

                    for(p=0; p<probes->rakes[nRakes].probesNumber; p++)
                    {
                        // set the dimension of data attributes (vector and indices set)
                        PetscMalloc(sizeof(PetscReal) * 3, &(probes->rakes[nRakes].locations[p]));
                        PetscMalloc(sizeof(PetscInt)    * 3, &(probes->rakes[nRakes].cells[p]));
                    }

                    // set rake name
                    probes->rakes[nRakes].rakeName = rakeName;

                    // set time name
                    probes->rakes[nRakes].timeName = "./postProcessing/" + probes->rakes[nRakes].rakeName + "/" + getTimeName(clock);

                    // get timeStart
                    readError = fscanf(f, "%s %lf\n", tmp, &(probes->rakes[nRakes].timeStart));

                    // get timeInterval
                    readError = fscanf(f, "%s %lf\n", tmp, &(probes->rakes[nRakes].timeInterval));

                    // get fields to gather
                    readError = fscanf(f, "%[^\n]", fields);
                    probes->rakes[nRakes].Uflag = foundInString(fields, "U");
                    probes->rakes[nRakes].Tflag = foundInString(fields, "T");

                    // set T flag to zero if T equation is not active
                    if(!flags->isTeqnActive) probes->rakes[nRakes].Tflag = 0;

                    readError = fscanf(f, "\n");
                    readError = fscanf(f, "%s\n", tmp);
                    if(strcmp(tmp, "locations") != 0)
                    {
                       char error[512];
                        sprintf(error, "probe rake %s has an incorrect file format, the correct format is\n\nprobesNumber 2\nrakeName     A1\ntimeStart    10.0\ntimeInterval 20.0\nfields       U,T\n\nlocations\n\n20.0 20.0 20.0\n23.1 45.1 45.9\n\n", probes->rakes[nRakes].rakeName.c_str());
                        fatalErrorInFunction("ProbesInitialize",  error);
                    }
                    readError = fscanf(f, "\n");

                    // initialize rake communicator color for this processor
                    PetscInt commColor = 0;

                    // loop over probes in the rake
                    for(p=0; p<probes->rakes[nRakes].probesNumber; p++)
                    {
                        PetscReal px, py, pz;
                        PetscInt    pi, pj, pk;

                        // read rake location
                        readError = fscanf(f, "%lf %lf %lf\n", &px, &py, &pz);
                        probes->rakes[nRakes].locations[p][0] = px;
                        probes->rakes[nRakes].locations[p][1] = py;
                        probes->rakes[nRakes].locations[p][2] = pz;

                        std::vector<cellIds> cells(nDomains);
                        std::vector<PetscInt>     owner(nDomains);
                        std::vector<PetscInt>     domainID(nDomains);

                        // loop over domains and see which domain controls this probe
                        for(d=0; d<nDomains; d++)
                        {
                            mesh_         *mesh   = domain[d].mesh;
                            DM            da = mesh->da, fda = mesh->fda;
                            DMDALocalInfo info = mesh->info;

                            PetscInt           xs = info.xs, xe = info.xs + info.xm;
                            PetscInt           ys = info.ys, ye = info.ys + info.ym;
                            PetscInt           zs = info.zs, ze = info.zs + info.zm;
                            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

                            PetscMPIInt        rank, nProcs;
                            PetscInt           i, j, k;
                            PetscInt           lxs, lxe, lys, lye, lzs, lze;

                            Cmpnts        ***cent;
                            PetscReal     ***aj;

                            lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
                            lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
                            lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

                            MPI_Comm_rank(mesh->MESH_COMM, &rank);
                            MPI_Comm_size(mesh->MESH_COMM, &nProcs);

                            DMDAVecGetArray(fda, mesh->lCent, &cent);
                            DMDAVecGetArray(da, mesh->lAj, &aj);

                            // max perturbation amplitude
                            PetscReal maxPerturb  = 1e-10;

                            // processor perturbation (changes between processors)
                            PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nProcs;

                            // find the closest cell
                            PetscReal lminDist = 1e20;
                            PetscReal gminDist = 1e20;
                            std::vector<PetscInt> indices{0, 0, 0};

                            for (k=lzs; k<lze; k++)
                            {
                                for (j=lys; j<lye; j++)
                                {
                                    for (i=lxs; i<lxe; i++)
                                    {
                                        PetscReal dist
                                        =
                                        std::sqrt
                                        (
                                            (cent[k][j][i].x - px - procContrib) * (cent[k][j][i].x - px - procContrib) +
                                            (cent[k][j][i].y - py - procContrib) * (cent[k][j][i].y - py - procContrib) +
                                            (cent[k][j][i].z - pz - procContrib) * (cent[k][j][i].z - pz - procContrib)
                                        );

                                        if(dist < lminDist)
                                        {
                                            // save distance value
                                            lminDist = dist;

                                            // save indices
                                            indices[0] = k;
                                            indices[1] = j;
                                            indices[2] = i;
                                        }
                                    }
                                }
                            }

                            // scatter the local distance to global using MIN operator
                            MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

                            std::vector<PetscInt> lindices{0, 0, 0};
                            std::vector<PetscInt> gindices{0, 0, 0};
                            PetscInt lownerRank = 0;
                            PetscInt gownerRank = 0;
                            PetscInt ldomainID  = 0;
                            PetscInt gdomainID  = 0;;

                            // now compare the distances: where they agree, this processor controls the probe
                            if(lminDist == gminDist)
                            {
                                // furthest possible distance as 2 * cellWidth
                                PetscReal cellWidth = 2*pow(1.0/aj[indices[0]][indices[1]][indices[2]], 1.0/3.0);

                                // probe is inside this domain
                                if(gminDist < cellWidth)
                                {
                                    // save indices
                                    lindices[0] = indices[0];
                                    lindices[1] = indices[1];
                                    lindices[2] = indices[2];

                                    // save owner rank
                                    lownerRank = rank;

                                    // set domain ID
                                    ldomainID = d;
                                }
                                // probe is outside this domain
                                else
                                {
                                    // set everything to -1
                                    lindices[0] = -1;
                                    lindices[1] = -1;
                                    lindices[2] = -1;
                                    lownerRank  = -1;
                                    ldomainID   = -1;
                                }
                            }

                            // scatter indices and owner
                            MPI_Allreduce(&lindices[0], &gindices[0], 3, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
                            MPI_Allreduce(&lownerRank, &gownerRank, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
                            MPI_Allreduce(&ldomainID, &gdomainID, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                            // save info for this domain
                            cells[d].k  = gindices[0]; // k index
                            cells[d].j  = gindices[1]; // j index
                            cells[d].i  = gindices[2]; // i index
                            owner[d]    = gownerRank;
                            domainID[d] = gdomainID;

                            DMDAVecRestoreArray(fda, mesh->lCent, &cent);
                            DMDAVecRestoreArray(da, mesh->lAj, &aj);

                        }

                        // choose the domain characterized by the finer mesh to be the
                        // one controlling this probe
                        PetscInt finerID = -1;

                        for(d=0; d<nDomains; d++)
                        {
                            finerID = PetscMax(finerID, domainID[d]);
                        }

                        probes->rakes[nRakes].cells[p][0] = cells[finerID].k; // k index
                        probes->rakes[nRakes].cells[p][1] = cells[finerID].j; // j index
                        probes->rakes[nRakes].cells[p][2] = cells[finerID].i; // i index
                        probes->rakes[nRakes].owner[p]    = owner[finerID];
                        probes->rakes[nRakes].domainID[p] = finerID;

                        if(finerID != -1)
                        {
                            // get this processor number from the communicator of the domain controlling this probe
                            PetscMPIInt rank; MPI_Comm_rank(domain[finerID].mesh->MESH_COMM, &rank);

                            // get the owner of this probe in the domain controlling this probe
                            PetscInt owner = probes->rakes[nRakes].owner[p];

                            // if this processor ownes this probe set color to 1
                            if(rank == owner)
                            {
                                commColor = 1;
                            }
                        }
                        else
                        {
                            char warning[256];
                            sprintf(warning, "Ignoring probe at location (%lf, %lf, %lf) in rake %s since it is outside the domain", probes->rakes[r].locations[p][0], probes->rakes[r].locations[p][1], probes->rakes[r].locations[p][2], probes->rakes[r].rakeName.c_str());
                            warningInFunction("ProbesInitialize",  warning);
                        }
                    }

                    // get this processor rank in the global communicator
                    PetscMPIInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

                    // create rake communicator
                    PetscSubcomm  psubcomm;
                    PetscSubcommCreate(PETSC_COMM_WORLD, &(psubcomm));
                    PetscSubcommSetTypeGeneral(psubcomm, commColor, rank);
                    probes->rakes[nRakes].RAKE_COMM  = PetscSubcommChild(psubcomm);
                    probes->rakes[nRakes].thisRakeControlled = commColor;

                    // increase rake count
                    nRakes++;

                    // close probe rake file
                    fclose(f);
                }
            }

            closedir(dir);
        }
        else
        {
            fatalErrorInFunction
            (
                "ProbesInitialize",
                "problem accessing acquisition data in sampling/probes/\n"
            );
        }

        // set domain probe pointers to the same memory location
        for(d=0; d<nDomains; d++)
        {
            domain[d].acquisition->probes = probes;
        }

        // set local probe rakes pointer (to be deleted, just for testing)
        probes = domain[0].acquisition->probes;

        PetscMPIInt rank;

        // create folders
        for(r=0; r<probes->nRakes; r++)
        {
            if(probes->rakes[r].thisRakeControlled)
            {
                MPI_Comm_rank(probes->rakes[r].RAKE_COMM, &rank);

                if(!rank)
                {
                    // create rake folder
                    errno = 0;
                    word rakeFolderName = "./postProcessing/" + probes->rakes[r].rakeName;
                    PetscInt dirRes = mkdir(rakeFolderName.c_str(), 0777);
                    if(dirRes != 0 && errno != EEXIST)
                    {
                       char error[512];
                        sprintf(error, "could not create %s directory\n", rakeFolderName.c_str());
                        fatalErrorInFunction("ProbesInitialize",  error);
                    }
                    else
                    {
                        // create time folder
                        errno = 0;
                        dirRes = mkdir(probes->rakes[r].timeName.c_str(), 0777);
                        if(dirRes != 0 && errno != EEXIST)
                        {
                           char error[512];
                            sprintf(error, "could not create %s directory\n", probes->rakes[r].timeName.c_str());
                            fatalErrorInFunction("ProbesInitialize",  error);
                        }
                        // if the time name exists remove everything inside
                        else if(errno == EEXIST)
                        {
                            remove_subdirs(probes->rakes[r].RAKE_COMM, probes->rakes[r].timeName.c_str());
                        }

                        // initialize velocity file
                        if(probes->rakes[r].Uflag) InitRakeFile(&(probes->rakes[r]), "U");

                        // initialize temperature file
                        if(probes->rakes[r].Tflag) InitRakeFile(&(probes->rakes[r]), "T");
                    }
                }

                // print acquisition system information
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   SAMPLING RAKE %ld\n",r);
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   number of probes                     : %ld\n", probes->rakes[r].probesNumber);
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   rake identification name             : %s\n", probes->rakes[r].rakeName.c_str());
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   path to rake data                    : %s\n", probes->rakes[r].timeName.c_str());
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   acq. starting time for this rake (s) : %lf\n", probes->rakes[r].timeStart);
                PetscPrintf(probes->rakes[r].RAKE_COMM, "   acq. time interval for this rake (s) : %lf\n", probes->rakes[r].timeInterval);
                if(probes->rakes[r].Uflag)
                {
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   get U                                : yes\n");
                }
                else
                {
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   get U                                : no\n");
                }
                if(probes->rakes[r].Tflag)
                {
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   get T                                : yes\n");
                }
                else
                {
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   get T                                : no\n");
                }

                for(p=0; p<probes->rakes[r].probesNumber; p++)
                {
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   probe %ld location                     : (%.2lf %.2lf %.2lf)\n", p, probes->rakes[r].locations[p][0], probes->rakes[r].locations[p][1], probes->rakes[r].locations[p][2]);
                    if(probes->rakes[r].domainID[p] != -1)
                    {
                        Cmpnts        ***cent;
                        mesh_         *mesh = domain[probes->rakes[r].domainID[p]].mesh;

                        DMDAVecGetArray(mesh->fda, mesh->lCent, &cent);

                        // scatter the closest cell center to the master processor for screen messaging purposes
                        std::vector<PetscReal> lc{0.0, 0.0, 0.0};
                        std::vector<PetscReal> gc{0.0, 0.0, 0.0};

                        MPI_Comm_rank(mesh->MESH_COMM, &rank);

                        if(rank==probes->rakes[r].owner[p])
                        {
                            lc[0] = cent[probes->rakes[r].cells[p][0]][probes->rakes[r].cells[p][1]][probes->rakes[r].cells[p][2]].x;
                            lc[1] = cent[probes->rakes[r].cells[p][0]][probes->rakes[r].cells[p][1]][probes->rakes[r].cells[p][2]].y;
                            lc[2] = cent[probes->rakes[r].cells[p][0]][probes->rakes[r].cells[p][1]][probes->rakes[r].cells[p][2]].z;
                        }

                        MPI_Allreduce(&lc[0], &gc[0], 3, MPIU_REAL, MPIU_SUM, probes->rakes[r].RAKE_COMM);

                        PetscPrintf(probes->rakes[r].RAKE_COMM, "   probe %ld closest cell center location : (%.2lf\t%.2lf\t%.2lf)\n", p, gc[0], gc[1], gc[2]);

                        DMDAVecRestoreArray(mesh->fda, mesh->lCent, &cent);
                    }

                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   probe %ld closest cell center indices  : (%ld\t%ld\t%ld)\n", p, probes->rakes[r].cells[p][0], probes->rakes[r].cells[p][1], probes->rakes[r].cells[p][2]);
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   probe %ld is owned by processor        : %ld\n", p, probes->rakes[r].owner[p]);
                    PetscPrintf(probes->rakes[r].RAKE_COMM, "   probe %ld is contained by domain       : %ld\n", p, probes->rakes[r].domainID[p]);
                }

                PetscPrintf(probes->rakes[r].RAKE_COMM, "\n");
            }

            MPI_Barrier(PETSC_COMM_WORLD);
        }

        PetscPrintf(PETSC_COMM_WORLD, "done\n\n");
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitRakeFile(probeRake *rake, const char *fieldName)
{
    PetscInt p;
    char fileName[500];
    sprintf(fileName, "%s/%s", rake->timeName.c_str(), fieldName);
    FILE *f = fopen(fileName, "w");

    if(!f)
    {
       char error[530];
        sprintf(error, "cannot open file %s\n", fileName);
        fatalErrorInFunction("ProbesInitialize",  error);
    }

    for(p=0; p<rake->probesNumber; p++)
    {
        if(rake->domainID[p] != -1)
        {
            PetscReal xp = rake->locations[p][0];
            PetscReal yp = rake->locations[p][1];
            PetscReal zp = rake->locations[p][2];
            fprintf(f, "# Probe %ld (%lf %lf %lf)\n", p, xp, yp, zp);
        }
    }

    fprintf(f, "#\tProbe\t\t\t\t\t\t");

    for(p=0; p<rake->probesNumber; p++)
    {
        if(rake->domainID[p] != -1)
        {
            fprintf(f, "%ld\t\t\t\t\t\t", p);
        }
    }

    fprintf(f, "\n");
    fprintf(f, "#\t Time\n");

    fclose(f);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeProbes(domain_ *domain)
{
    if(domain[0].acquisition->isProbesActive)
    {
        PetscInt r, p;
        rakes  *probes = domain[0].acquisition->probes;
        clock_ *clock  = domain[0].clock;

        for(r=0; r<probes->nRakes; r++)
        {
            // get ptr to this probe rake
            probeRake *rake = &(probes->rakes[r]);

            PetscMPIInt rakeRank;
            MPI_Comm_rank(rake->RAKE_COMM, &rakeRank);

            if(rake->thisRakeControlled)
            {
                double timeStart    = rake->timeStart;
                double timeInterval = rake->timeInterval;

                // check the time and see if must write
                if
                (
                    clock->time >= timeStart && // past acquisition start
                    (clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval) < 1e-10

                )
                {
                    // initialize local vectors
                    std::vector<Cmpnts> lprobeValuesU;
                    std::vector<Cmpnts> gprobeValuesU;

                    // initialize local vectors
                    std::vector<PetscReal> lprobeValuesT;
                    std::vector<PetscReal> gprobeValuesT;

                    if(rake->Uflag)
                    {
                        lprobeValuesU.resize(rake->probesNumber);
                        gprobeValuesU.resize(rake->probesNumber);
                    }

                    if(rake->Tflag)
                    {
                        lprobeValuesT.resize(rake->probesNumber);
                        gprobeValuesT.resize(rake->probesNumber);
                    }

                    for(p=0; p<rake->probesNumber; p++)
                    {
                        if(rake->domainID[p] != -1)
                        {
                            mesh_         *mesh = domain[rake->domainID[p]].mesh;
                            DMDALocalInfo info = mesh->info;
                            PetscInt           xs = info.xs, xe = info.xs + info.xm;
                            PetscInt           ys = info.ys, ye = info.ys + info.ym;
                            PetscInt           zs = info.zs, ze = info.zs + info.zm;
                            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

                            PetscInt           i, j, k;
                            PetscInt           lxs, lxe, lys, lye, lzs, lze;

                            lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
                            lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
                            lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

                            PetscMPIInt meshRank;
                            MPI_Comm_rank(mesh->MESH_COMM, &meshRank);

                            if(rake->Uflag)
                            {
                                Cmpnts ***ucat;

                                // get u equation
                                ueqn_ *ueqn = domain[rake->domainID[p]].ueqn;

                                DMDAVecGetArray(mesh->fda, ueqn->lUcat, &ucat);

                                if(meshRank == rake->owner[p])
                                {
                                    k = rake->cells[p][0];
                                    j = rake->cells[p][1];
                                    i = rake->cells[p][2];
                                    lprobeValuesU[p].x = ucat[k][j][i].x;
                                    lprobeValuesU[p].y = ucat[k][j][i].y;
                                    lprobeValuesU[p].z = ucat[k][j][i].z;
                                }

                                DMDAVecRestoreArray(mesh->fda, ueqn->lUcat, &ucat);
                            }

                            if(rake->Tflag)
                            {
                                PetscReal ***t;

                                // get u equation
                                teqn_ *teqn = domain[rake->domainID[p]].teqn;

                                DMDAVecGetArray(mesh->da, teqn->lTmprt, &t);

                                if(meshRank == rake->owner[p])
                                {
                                    k = rake->cells[p][0];
                                    j = rake->cells[p][1];
                                    i = rake->cells[p][2];
                                    lprobeValuesT[p] = t[k][j][i];
                                }

                                DMDAVecRestoreArray(mesh->da, teqn->lTmprt, &t);
                            }
                        }
                    }

                    if(rake->Uflag)
                    {
                        MPI_Reduce(&lprobeValuesU[0], &gprobeValuesU[0], 3*rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);

                        // only master process writes on the file
                        if(!rakeRank)
                        {
                            // write velocity
                            FILE *fu;
                            word fileName = rake->timeName + "/U";
                            fu = fopen(fileName.c_str(), "a");

                            fprintf(fu, "\t %.3lf\t\t\t\t\t\t", clock->time);
                            for(p=0; p<rake->probesNumber; p++)
                            {
                                if(rake->domainID[p] != -1)
                                {
                                    fprintf(fu, "(%.10lf  %.10lf  %.10lf)\t\t", gprobeValuesU[p].x, gprobeValuesU[p].y, gprobeValuesU[p].z);
                                }
                            }
                            fprintf(fu, "\n");
                            fclose(fu);
                        }
                    }

                    if(rake->Tflag)
                    {
                        MPI_Reduce(&lprobeValuesT[0], &gprobeValuesT[0], rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);

                        // only master process writes on the file
                        if(!rakeRank)
                        {
                            // write velocity
                            FILE *ft;
                            word fileName = rake->timeName + "/T";
                            ft = fopen(fileName.c_str(), "a");

                            fprintf(ft, "\t %.3lf\t\t\t\t\t\t", clock->time);
                            for(p=0; p<rake->probesNumber; p++)
                            {
                                if(rake->domainID[p] != -1)
                                {
                                    fprintf(ft, "%.10lf\t\t", gprobeValuesT[p]);
                                }
                            }
                            fprintf(ft, "\n");
                            fclose(ft);
                        }
                    }

                    std::vector<Cmpnts> ().swap(lprobeValuesU);
                    std::vector<Cmpnts> ().swap(gprobeValuesU);

                    std::vector<PetscReal> ().swap(lprobeValuesT);
                    std::vector<PetscReal> ().swap(gprobeValuesT);
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeQCritIO(acquisition_ *acquisition)
{
    mesh_         *mesh = acquisition->access->mesh;
    ueqn_         *ueqn = acquisition->access->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***ucat;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert, ***aj, ***q, ***ocode;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;


    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, acquisition->fields->Q, &q);

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        if( isIBMSolidCell(k, j, i, nvert))
        {
            continue;
        }

        PetscReal ajc = aj[k][j][i];
        PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;
        PetscReal dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

        Compute_du_center
        (
            mesh, i, j, k,
            mx, my, mz,
            ucat, nvert,
            &dudc, &dvdc, &dwdc,
            &dude, &dvde, &dwde,
            &dudz, &dvdz, &dwdz
        );

        Compute_du_dxyz
        (
            mesh,
            csi0, csi1, csi2,
            eta0, eta1, eta2,
            zet0, zet1, zet2,
            ajc,
            dudc, dvdc, dwdc,
            dude, dvde, dwde,
            dudz, dvdz, dwdz,
            &du_dx, &dv_dx, &dw_dx,
            &du_dy, &dv_dy, &dw_dy,
            &du_dz, &dv_dz, &dw_dz
        );

        PetscReal w11 = 0;
        PetscReal w12 = 0.5*(du_dy - dv_dx);
        PetscReal w13 = 0.5*(du_dz - dw_dx);
        PetscReal w21 = -w12;
        PetscReal w22 = 0.;
        PetscReal w23 = 0.5*(dv_dz - dw_dy);
        PetscReal w31 = -w13;
        PetscReal w32 = -w23;
        PetscReal w33 = 0.;

        q[k][j][i]
        =
        w11*w11 + w12*w12 + w13*w13 +
        w21*w21 + w22*w22 + w23*w23 +
        w31*w31 + w32*w32 + w33*w33;
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, acquisition->fields->Q, &q);

    MPI_Barrier ( mesh->MESH_COMM );

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode computeCoriolisIO(acquisition_ *acquisition)
{
    mesh_         *mesh = acquisition->access->mesh;
    ueqn_         *ueqn = acquisition->access->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***source, ***ucat, ***cent;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert;

    PetscReal     fc = ueqn->access->abl->fc; // coriolis parameter

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(acquisition->fields->Coriolis, 0.);

    // damping viscosity for fringe region exclusion
    double nu_fringe;

    // fringe region parameters (set only if active)
    double xS;
    double xE;
    double xD;

    if(ueqn->access->flags->isXDampingActive)
    {
        xS     = ueqn->access->abl->xDampingStart;
        xE     = ueqn->access->abl->xDampingEnd;
        xD     = ueqn->access->abl->xDampingDelta;
    }
    else
    {
        nu_fringe = 1.0;
    }

    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);
    DMDAVecGetArray(da,  mesh->lNvert,&nvert);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, acquisition->fields->Coriolis,  &source);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(ueqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k
                    double x = (cent[k][j][i].x   - mesh->bounds.xmin);

                    // compute Stipa viscosity at i,j,k,
                    nu_fringe = viscStipa(xS, xE, xD, x);
                }
                // might need an else to reset x_fringe to 1 at each iteration

                if
                (
                    isFluidCell(k, j, i, nvert)
                )
                {
                    source[k][j][i].x
                    +=
                    nu_fringe *
                    (
                        -2.0 *
                        (
                            - fc * ucat[k][j][i].y * csi[k][j][i].x +
                              fc * ucat[k][j][i].x * csi[k][j][i].y
                        )
                    );

                    source[k][j][i].y
                    +=
                    nu_fringe *
                    (
                        -2.0 *
                        (
                          - fc * ucat[k][j][i].y * eta[k][j][i].x +
                            fc * ucat[k][j][i].x * eta[k][j][i].y
                        )
                    );

                    source[k][j][i].z
                    +=
                    nu_fringe *
                    (
                        -2.0 *
                        (
                          - fc * ucat[k][j][i].y * zet[k][j][i].x +
                            fc * ucat[k][j][i].x * zet[k][j][i].y
                        )
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,  &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert,&nvert);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(fda, acquisition->fields->Coriolis,  &source);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeDrivingSourceIO(acquisition_ *acquisition)
{
    mesh_         *mesh = acquisition->access->mesh;
    ueqn_         *ueqn = acquisition->access->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***source, ***sourceu, ***ucat, ***cent;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(acquisition->fields->Driving, 0.);

    // damping viscosity for fringe region exclusion
    double nu_fringe;

    // fringe region parameters (set only if active)
    double xS;
    double xE;
    double xD;

    if(ueqn->access->flags->isXDampingActive)
    {
        xS     = ueqn->access->abl->xDampingStart;
        xE     = ueqn->access->abl->xDampingEnd;
        xD     = ueqn->access->abl->xDampingDelta;
    }
    else
    {
        nu_fringe = 1.0;
    }

    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);
    DMDAVecGetArray(da,  mesh->lNvert,&nvert);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, acquisition->fields->Driving,  &source);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(fda, ueqn->sourceU, &sourceu);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(ueqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k
                    double x = (cent[k][j][i].x   - mesh->bounds.xmin);

                    // compute Stipa viscosity at i,j,k,
                    nu_fringe = viscStipa(xS, xE, xD, x);
                }
                // might need an else to reset x_fringe to 1 at each iteration

                if
                (
                    isFluidCell(k, j, i, nvert)
                )
                {
                    source[k][j][i].x
                    +=
                    nu_fringe *
                    (
                        sourceu[k][j][i].x * csi[k][j][i].x +
                        sourceu[k][j][i].y * csi[k][j][i].y +
                        sourceu[k][j][i].z * csi[k][j][i].z
                    );


                    source[k][j][i].y
                    +=
                    nu_fringe *
                    (
                        sourceu[k][j][i].x * eta[k][j][i].x +
                        sourceu[k][j][i].y * eta[k][j][i].y +
                        sourceu[k][j][i].z * eta[k][j][i].z
                    );


                    source[k][j][i].z
                    +=
                    nu_fringe *
                    (
                        sourceu[k][j][i].x * zet[k][j][i].x +
                        sourceu[k][j][i].y * zet[k][j][i].y +
                        sourceu[k][j][i].z * zet[k][j][i].z
                    );

                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,  &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert,&nvert);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(fda, acquisition->fields->Driving,  &source);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, ueqn->sourceU, &sourceu);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeXDampingIO(acquisition_ *acquisition)
{
    abl_          *abl  = acquisition->access->abl;
    mesh_         *mesh = acquisition->access->mesh;
    ueqn_         *ueqn = acquisition->access->ueqn;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***source, ***ucat, ***ucatP, ***cent;
    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***ucont;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    precursor_    *precursor;
    domain_       *pdomain;
    PetscInt      kStart;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(acquisition->fields->xDamping, 0.);

    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, acquisition->fields->xDamping,  &source);

    if(ueqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->fda, pdomain->ueqn->lUcat,  &ucatP);
                kStart = precursor->map.kStart;
            }
        }
    }

    // z damping layer
    PetscReal alphaZ = abl->zDampingAlpha;
    PetscReal zS     = abl->zDampingStart;
    PetscReal zE     = abl->zDampingEnd;

    // x damping layer
    PetscReal alphaX = abl->xDampingAlpha;
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;
    PetscReal xD     = abl->xDampingDelta;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(ueqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k
                    PetscReal x = (cent[k][j][i].x   - mesh->bounds.xmin);

                    // compute Nordstrom viscosity at i,j,k
                    PetscReal nu_fringe = viscNordstrom(alphaX, xS, xE, xD, x);

                    // X DAMPING LAYER
                    // ---------------

                    if(abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2)
                    {
                        Cmpnts uBar  = nSet(abl->uBarInstX[j][i]);

                        // i-fluxes
                        source[k][j][i].x
                        +=
                        nu_fringe *
                        (
                            (
                                (uBar.x - ucat[k][j][i].x) * csi[k][j][i].x +
                                (uBar.y - ucat[k][j][i].y) * csi[k][j][i].y +
                                (uBar.z - ucat[k][j][i].z) * csi[k][j][i].z
                            )
                        );

                        // j-fluxes
                        source[k][j][i].y
                        +=
                        nu_fringe *
                        (
                            (
                                (uBar.x - ucat[k][j][i].x) * eta[k][j][i].x +
                                (uBar.y - ucat[k][j][i].y) * eta[k][j][i].y +
                                (uBar.z - ucat[k][j][i].z) * eta[k][j][i].z
                            )
                        );

                        // k-fluxes
                        source[k][j][i].z
                        +=
                        nu_fringe *
                        (
                            (
                                (uBar.x - ucat[k][j][i].x) * zet[k][j][i].x +
                                (uBar.y - ucat[k][j][i].y) * zet[k][j][i].y +
                                (uBar.z - ucat[k][j][i].z) * zet[k][j][i].z
                            )
                        );
                    }
                    else if(abl->xFringeUBarSelectionType == 3)
                    {
                        Cmpnts uBar;

                        if(precursor->thisProcessorInFringe)
                        {
                            uBar = nSet(ucatP[k+kStart][j][i]);
                        }
                        else
                        {
                            uBar = nSet(ucat[k][j][i]);
                        }

                        source[k][j][i].x
                        +=
                        nu_fringe *
                        (
                            (uBar.x - ucat[k][j][i].x) * csi[k][j][i].x +
                            (uBar.y - ucat[k][j][i].y) * csi[k][j][i].y +
                            (uBar.z - ucat[k][j][i].z) * csi[k][j][i].z
                        );

                        source[k][j][i].y
                        +=
                        nu_fringe *
                        (
                            (uBar.x - ucat[k][j][i].x) * eta[k][j][i].x +
                            (uBar.y - ucat[k][j][i].y) * eta[k][j][i].y +
                            (uBar.z - ucat[k][j][i].z) * eta[k][j][i].z
                        );

                        source[k][j][i].z
                        +=
                        nu_fringe *
                        (
                            (uBar.x - ucat[k][j][i].x) * zet[k][j][i].x +
                            (uBar.y - ucat[k][j][i].y) * zet[k][j][i].y +
                            (uBar.z - ucat[k][j][i].z) * zet[k][j][i].z
                        );
                    }
                }

                // Z DAMPING LAYER
                // ---------------
                if(ueqn->access->flags->isZDampingActive)
                {
                    // compute cell center z at i,j,k
                    PetscReal z = (cent[k][j][i].z   - mesh->bounds.zmin);

                    // compute Rayleigh viscosity at i,j,k and i,j+1,k points
                    PetscReal nud_rayleigh = viscRayleigh(alphaZ, zS, zE, z);

                    // damp also x and y components
                    if(abl->zDampingAlsoXY)
                    {
                        // i-fluxes: dampen w.r.t. uBarMean
                        source[k][j][i].x
                        +=
                        nud_rayleigh *
                        (
                            abl->uBarMeanZ[j].x -
                            (
                                ucat[k][j][i].x * csi[k][j][i].x +
                                ucat[k][j][i].y * csi[k][j][i].y +
                                ucat[k][j][i].z * csi[k][j][i].z
                            )
                        );

                        // k-fluxes: dampen w.r.t. uBarMean
                        source[k][j][i].z
                        +=
                        nud_rayleigh *
                        (
                            abl->uBarMeanZ[j].z -
                            (
                                ucat[k][j][i].x * zet[k][j][i].x +
                                ucat[k][j][i].y * zet[k][j][i].y +
                                ucat[k][j][i].z * zet[k][j][i].z
                            )
                        );
                    }

                    // j-fluxes: total damping to reach no penetration at jRight
                    source[k][j][i].y
                    +=
                    nud_rayleigh *
                    (
                        -1.0 *
                        (
                            ucat[k][j][i].x * eta[k][j][i].x +
                            ucat[k][j][i].y * eta[k][j][i].y +
                            ucat[k][j][i].z * eta[k][j][i].z
                        )
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,  &zet);

    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, acquisition->fields->xDamping,  &source);

    if(ueqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->fda, pdomain->ueqn->lUcat,  &ucatP);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeSideForceIO(acquisition_ *acquisition)
{
    mesh_         *mesh = acquisition->access->mesh;
    ueqn_         *ueqn = acquisition->access->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***source, ***ucat, ***cent;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert;

    PetscReal     fc      = ueqn->access->abl->fc; // coriolis parameter
    PetscReal     xStart  = ueqn->access->abl->xStartSideF,
                  zStart  = ueqn->access->abl->zStartSideF,
                  xEnd    = ueqn->access->abl->xEndSideF,
                  zEnd    = ueqn->access->abl->zEndSideF;

    double        K       = 5.0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(acquisition->fields->SideForce, 0.);

    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);
    DMDAVecGetArray(da,  mesh->lNvert,&nvert);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, acquisition->fields->SideForce,  &source);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                double coeff = 0.0;

                // compute cell center x at i,j,k
                double x    = cent[k][j][i].x - mesh->bounds.xmin;

                // compute cell center z at i,j,k
                double z    = cent[k][j][i].z - mesh->bounds.zmin;

                if(x < xEnd && x > xStart && z < zEnd && z > zStart) coeff = 1;

                if
                (
                    isFluidCell(k, j, i, nvert)
                )
                {
                    source[k][j][i].x
                    +=
                    coeff *
                    (
                        2.0 * K *
                        (
                            - fc * 10.0 * csi[k][j][i].x
                        )
                    );

                    source[k][j][i].y
                    +=
                    coeff *
                    (
                        2.0 * K *
                        (
                          0.0
                        )
                    );

                    source[k][j][i].z
                    +=
                    coeff *
                    (
                        2.0 * K *
                        (
                           0.0
                        )
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,  &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert,&nvert);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(fda, acquisition->fields->SideForce,  &source);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode averaging3LMInitialize(domain_ *domain)
{
    // 3LM averaging is done only on background domain
    acquisition_ *acquisition = domain[0].acquisition;

    if(acquisition->isAverage3LMActive)
    {
        mesh_ *mesh = acquisition->access->mesh;

        PetscPrintf(mesh->MESH_COMM, "Creating 3LM acquisition grid...\n");
        PetscPrintf(mesh->MESH_COMM, "   reading parameters...\n");

        if(mesh->meshFileType != "cartesian")
        {
           char error[512];
            sprintf(error, "3LM averaging only available for cartesian meshes\n");
            fatalErrorInFunction("averaging3LMInitialize",  error);
        }

        // allocate memory for the 3LM data
        PetscMalloc(sizeof(data3LM), &(acquisition->LM3));
        data3LM *lm3 = acquisition->LM3;

        char path2dict[256];
        sprintf(path2dict, "./sampling/3LM");

        // read input file located inside sampling/ and called 3LM
        readDictDouble(path2dict, "avgStartTime",  &(lm3->avgStartTime));
        readDictDouble(path2dict, "avgPeriod",     &(lm3->avgPrd));
        readDictInt(path2dict,    "nspw",          &(lm3->nspw));
        readDictInt(path2dict,    "nstw",          &(lm3->nstw));
        readDictVector(path2dict, "upDir",         &(lm3->upDir));
        readDictVector(path2dict, "streamDir",     &(lm3->streamDir));

        // make the direction vector unitary
        mUnit(lm3->upDir);
        mUnit(lm3->streamDir);

        // set the spanwise unit vector
        lm3->spanDir = nCross(lm3->upDir, lm3->streamDir);

        // initialize snapshot weighting to zero
        lm3->avgWeight = 0;

        // allocate memory for points
        lm3->points       = (Cmpnts**)malloc(lm3->nstw*sizeof(Cmpnts*));
        lm3->closestCells = (cellIds**)malloc(lm3->nstw*sizeof(cellIds*));

        for(PetscInt pk=0; pk<lm3->nstw; pk++)
        {
            lm3->points[pk]       = (Cmpnts*)malloc(lm3->nspw*sizeof(Cmpnts));
            lm3->closestCells[pk] = (cellIds*)malloc(lm3->nspw*sizeof(cellIds));
        }

        // allocate memory for levels
        lm3->levels = (level3LM**)malloc(3*sizeof(level3LM*));

        for(PetscInt li=0; li<3; li++)
        {
            lm3->levels[li] = (level3LM*)malloc(sizeof(level3LM));

            level3LM *lev = lm3->levels[li];

            // set level name
            PetscMalloc(sizeof(word), &(lev->levelName));
            char levName[256];
            sprintf(levName, "level%ld", li+1);
            lev->levelName = levName;

            // set start and end height
            PetscReal hStart, hEnd;
            PetscMalloc(sizeof(Cmpnts), &(lev->hStart));
            PetscMalloc(sizeof(Cmpnts), &(lev->hEnd));
            readSubDictDouble(path2dict, levName, "hStart",  &hStart);
            readSubDictDouble(path2dict, levName, "hEnd",    &hEnd);
            lev->hStart = nScale(hStart, lm3->upDir);
            lev->hEnd   = nScale(hEnd, lm3->upDir);

            // allocate memory for cells, velocity and pressure
            lev->U = (Cmpnts**)malloc(lm3->nstw*sizeof(Cmpnts*));
            lev->P = (PetscReal**)malloc(lm3->nstw*sizeof(PetscReal*));

            for(PetscInt pk=0; pk<lm3->nstw; pk++)
            {
                lev->U[pk] = (Cmpnts*)malloc(lm3->nspw*sizeof(Cmpnts));
                lev->P[pk] = (PetscReal*)malloc(lm3->nspw*sizeof(PetscReal));

                for(PetscInt pi=0; pi<lm3->nspw; pi++)
                {
                    lev->U[pk][pi].x = 0.0;
                    lev->U[pk][pi].y = 0.0;
                    lev->U[pk][pi].z = 0.0;

                    lev->P[pk][pi] = 0.0;
                }
            }
        }

        // build the point mesh
        Cmpnts max, min;
        min.x = mesh->bounds.xmin;
        max.x = mesh->bounds.xmax;
        min.y = mesh->bounds.ymin;
        max.y = mesh->bounds.ymax;
        min.z = mesh->bounds.zmin;
        max.z = mesh->bounds.zmax;

        Cmpnts maxDist       = nSub(max, min);

        PetscReal dStreamMag = nDot(maxDist, lm3->streamDir) / (lm3->nstw - 1);
        Cmpnts dStream       = nScale(dStreamMag, lm3->streamDir);

        PetscReal dSpanMag   = nDot(maxDist, lm3->spanDir) / (lm3->nspw - 1);
        Cmpnts dSpan         = nScale(dSpanMag, lm3->spanDir);

        PetscReal dVertMag   = nDot(min, lm3->upDir);
        Cmpnts dVert         = nScale(dVertMag, lm3->upDir);

        Cmpnts point         = min;

        // subtract vertical component so that up coord is always zero
        mSub(point, dVert);

        // loop over 3LM grid points
        for(PetscInt pk=0; pk<lm3->nstw; pk++)
        {
            // scaling in the streamwise direction
            Cmpnts dk = nScale((PetscReal)pk, dStream);

            for(PetscInt pi=0; pi<lm3->nspw; pi++)
            {
                // scaling in the spanwise direction
                Cmpnts di = nScale((PetscReal)pi, dSpan);

                // total scaling
                Cmpnts scale = nSum(dk, di);

                // set the point
                lm3->points[pk][pi] = nSum(point, scale);

                // initialize closest cells to zero
                lm3->closestCells[pk][pi].i = 0;
                lm3->closestCells[pk][pi].j = 0;
                lm3->closestCells[pk][pi].k = 0;
            }
        }

        write3LMPoints(acquisition);

        // now do the search and save the closest i,k indices to the vertical
        // averaging line in each processor. Here is where the assumptions of
        // cartesian mesh and vertical direction j are made. Return error if these
        // hypotesis are not verified.
        findAvgLineIds(acquisition);

        PetscPrintf(mesh->MESH_COMM, "done\n\n");

        MPI_Barrier(mesh->MESH_COMM);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeAveraging3LM(domain_ *domain)
{
    acquisition_ *acquisition = domain[0].acquisition;

    if(acquisition->isAverage3LMActive)
    {
        mesh_  *mesh  = acquisition->access->mesh;
        clock_ *clock = acquisition->access->clock;
        data3LM *lm3  = acquisition->LM3;

        ueqn_   *ueqn = acquisition->access->ueqn;
        peqn_   *peqn = acquisition->access->peqn;

        // check if must accumulate averaged fields
        PetscInt accumulateAvg = 0;

        PetscReal startTimeAvg         = lm3->avgStartTime;
        PetscReal timeIntervalAvg      = lm3->avgPrd;

        if
        (
            clock->time >= startTimeAvg &&
            (clock->time - startTimeAvg ) / timeIntervalAvg -
            std::floor((clock->time - startTimeAvg) / timeIntervalAvg) < 1e-10
        )
        {
            accumulateAvg = 1;
        }

        if(accumulateAvg)
        {
            DM               da = mesh->da, fda = mesh->fda;
            DMDALocalInfo    info = mesh->info;
            PetscInt         xs = info.xs, xe = info.xs + info.xm;
            PetscInt         ys = info.ys, ye = info.ys + info.ym;
            PetscInt         zs = info.zs, ze = info.zs + info.zm;
            PetscInt         mx = info.mx, my = info.my, mz = info.mz;

            PetscInt         lxs, lxe, lys, lye, lzs, lze;
            PetscInt         i, j, k, pi, pk, r, l;

            Cmpnts           ***cent, ***ucat;
            PetscReal        ***p;

            PetscReal        ts, te;

            PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
            PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

            lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
            lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
            lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

            PetscTime(&ts);

            DMDAVecGetArray(fda, mesh->lCent, &cent);
            DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
            DMDAVecGetArray(da, peqn->lP, &p);

            // set time averaging weights
            PetscReal mN = (PetscReal)lm3->avgWeight;
            PetscReal m1 = mN  / (mN + 1.0);
            PetscReal m2 = 1.0 / (mN + 1.0);

            // get up direction
            Cmpnts  upDir = lm3->upDir;

            // loop over levels
            for(l=0; l<3; l++)
            {
                std::vector<std::vector<Cmpnts>>    ldepthU(lm3->nstw);
                std::vector<std::vector<PetscReal>> ldepthP(lm3->nstw);
                std::vector<std::vector<PetscInt>>  ldepthN(lm3->nstw);
                std::vector<std::vector<Cmpnts>>    gdepthU(lm3->nstw);
                std::vector<std::vector<PetscReal>> gdepthP(lm3->nstw);
                std::vector<std::vector<PetscInt>>  gdepthN(lm3->nstw);

                // get level pointer for speed
                level3LM *lev = lm3->levels[l];

                PetscReal hStart = nDot(lev->hStart, upDir);
                PetscReal hEnd   = nDot(lev->hEnd, upDir);

                // loop over 3LM grid points
                for(pk=0; pk<lm3->nstw; pk++)
                {
                    // resize
                    ldepthU[pk].resize(lm3->nspw);
                    ldepthP[pk].resize(lm3->nspw);
                    ldepthN[pk].resize(lm3->nspw);
                    gdepthU[pk].resize(lm3->nspw);
                    gdepthP[pk].resize(lm3->nspw);
                    gdepthN[pk].resize(lm3->nspw);

                    for(pi=0; pi<lm3->nspw; pi++)
                    {
                        // set to zero
                        ldepthU[pk][pi].x = 0.0;
                        ldepthU[pk][pi].y = 0.0;
                        ldepthU[pk][pi].z = 0.0;
                        ldepthP[pk][pi]   = 0.0;
                        ldepthN[pk][pi]   = 0;
                        gdepthU[pk][pi].x = 0.0;
                        gdepthU[pk][pi].y = 0.0;
                        gdepthU[pk][pi].z = 0.0;
                        gdepthP[pk][pi]   = 0.0;
                        gdepthN[pk][pi]   = 0;

                        // number of depth averaging points

                        // get averaging line ids
                        cellIds lineIds = lm3->closestCells[pk][pi];

                        // test if this averaging line intercepts this processor
                        if
                        (
                            lineIds.k <= lze && lineIds.k >= lzs &&
                            lineIds.i <= lxe && lineIds.i >= lxs
                        )
                        {
                            k = lineIds.k;
                            i = lineIds.i;

                            // loop over line points
                            for(j=lys; j<lye; j++)
                            {
                                // get cell height
                                PetscReal hCell = nDot(cent[k][j][i], upDir);

                                // test if inside this level bounds
                                if(hCell < hEnd && hCell > hStart)
                                {
                                    // cumulate for average
                                    ldepthU[pk][pi].x += ucat[k][j][i].x;
                                    ldepthU[pk][pi].y += ucat[k][j][i].y;
                                    ldepthU[pk][pi].z += ucat[k][j][i].z;
                                    ldepthP[pk][pi]   += p[k][j][i];
                                    ldepthN[pk][pi]++;
                                }
                            }
                        }
                    }
                }

                // scatter this level depth averages on master node
                for(pk=0; pk<lm3->nstw; pk++)
                {
                    MPI_Reduce(&(ldepthU[pk][0]), &(gdepthU[pk][0]), 3*lm3->nspw, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
                    MPI_Reduce(&(ldepthP[pk][0]), &(gdepthP[pk][0]), lm3->nspw, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
                    MPI_Reduce(&(ldepthN[pk][0]), &(gdepthN[pk][0]), lm3->nspw, MPIU_INT, MPI_SUM, 0, mesh->MESH_COMM);
                }

                // time average
                for(pk=0; pk<lm3->nstw; pk++)
                {
                    for(pi=0; pi<lm3->nspw; pi++)
                    {
                        if(!rank)
                        {
                            lev->U[pk][pi].x = m1 * lev->U[pk][pi].x + m2 * gdepthU[pk][pi].x / gdepthN[pk][pi];
                            lev->U[pk][pi].y = m1 * lev->U[pk][pi].y + m2 * gdepthU[pk][pi].y / gdepthN[pk][pi];
                            lev->U[pk][pi].z = m1 * lev->U[pk][pi].z + m2 * gdepthU[pk][pi].z / gdepthN[pk][pi];
                            lev->P[pk][pi]   = m1 * lev->P[pk][pi]   + m2 * gdepthP[pk][pi]   / gdepthN[pk][pi];
                        }
                    }

                    // clean memory
                    std::vector<Cmpnts>    ().swap(ldepthU[pk]);
                    std::vector<Cmpnts>    ().swap(gdepthU[pk]);
                    std::vector<PetscReal> ().swap(ldepthP[pk]);
                    std::vector<PetscReal> ().swap(gdepthP[pk]);
                    std::vector<PetscInt>  ().swap(ldepthN[pk]);
                    std::vector<PetscInt>  ().swap(gdepthN[pk]);
                }
            }

            DMDAVecRestoreArray(fda, mesh->lCent, &cent);
            DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
            DMDAVecRestoreArray(da, peqn->lP, &p);

            lm3->avgWeight++;

            PetscTime(&te);

            PetscPrintf(mesh->MESH_COMM, "Averaged 3LM data in %lf s\n", te-ts);
        }

        // write fields
        write3LMFields(acquisition);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode write3LMPoints(acquisition_ *acquisition)
{
    mesh_  *mesh  = acquisition->access->mesh;
    clock_ *clock = acquisition->access->clock;
    data3LM *lm3  = acquisition->LM3;

    PetscPrintf(mesh->MESH_COMM, "   writing 3LM mesh...\n");

    PetscMPIInt rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // create directories
    if
    (
        clock->it == clock->itStart && !rank
    )
    {
        errno = 0;
        PetscInt dirRes = mkdir("./postProcessing/3LM/", 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create postProcessing/3LM directory\n");
            fatalErrorInFunction("write3LMPoints",  error);
        }
    }

    // write 3LM mesh points
    if(!rank)
    {
        FILE *fp = fopen("postProcessing/3LM/points", "w");
        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file postProcessing/3LM/points\n");
            fatalErrorInFunction("write3LMPoints",  error);
        }
        else
        {
            // write header lines

            PetscFPrintf(mesh->MESH_COMM, fp, "Points given as table of (xyz) coordinates where each line\n");
            PetscFPrintf(mesh->MESH_COMM, fp, "is an array of points in the spanwise direction and each\n");
            PetscFPrintf(mesh->MESH_COMM, fp, "column is an array of points in the streamwise direction\n\n");
            PetscFPrintf(mesh->MESH_COMM, fp, "x: streamwise direction at hub height\n");
            PetscFPrintf(mesh->MESH_COMM, fp, "y: spanwise direction at hub height\n");
            PetscFPrintf(mesh->MESH_COMM, fp, "z: non-meaningful direction (always zero)\n\n");

            // write mesh points
            for(PetscInt pk=0; pk<lm3->nstw; pk++)
            {
                for(PetscInt pi=0; pi<lm3->nspw; pi++)
                {
                    PetscFPrintf(mesh->MESH_COMM, fp, "(%*.2lf %*.2lf %*.2lf) ", -15, lm3->points[pk][pi].x, -10, lm3->points[pk][pi].y, 10, lm3->points[pk][pi].z);
                }
                // new line
                PetscFPrintf(mesh->MESH_COMM, fp, "\n");
            }

            fclose(fp);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode write3LMFields(acquisition_ *acquisition)
{
    mesh_  *mesh  = acquisition->access->mesh;
    data3LM *lm3  = acquisition->LM3;
    io_     *io   = acquisition->access->io;

    // check if must write
    if(io->runTimeWrite)
    {
        PetscMPIInt rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

        // write 3LM mesh points
        if(!rank)
        {
            for(PetscInt l=0; l<3; l++)
            {
                // get this level pointer
                level3LM *lev = lm3->levels[l];

                // write velocity
                {
                    char fileName[256];
                    sprintf(fileName, "postProcessing/3LM/U_L%ld",l+1);

                    FILE *fp = fopen(fileName, "w");
                    if(!fp)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("write3LMFields",  error);
                    }
                    else
                    {
                        // write header lines

                        PetscFPrintf(mesh->MESH_COMM, fp, "U given as table of (x y z) values for each point, where each\n");
                        PetscFPrintf(mesh->MESH_COMM, fp, "line is an array of points in the spanwise direction and each\n");
                        PetscFPrintf(mesh->MESH_COMM, fp, "column is an array of points in the streamwise direction\n\n");

                        // write mesh points
                        for(PetscInt pk=0; pk<lm3->nstw; pk++)
                        {
                            for(PetscInt pi=0; pi<lm3->nspw; pi++)
                            {
                                PetscFPrintf(mesh->MESH_COMM, fp, "(%*.4lf %*.4lf %*.4lf) ", -15, lev->U[pk][pi].x, -10, lev->U[pk][pi].y, 10, lev->U[pk][pi].z);
                            }
                            // new line
                            PetscFPrintf(mesh->MESH_COMM, fp, "\n");
                        }

                        fclose(fp);
                    }
                }

                // write pressure
                {
                    char fileName[256];
                    sprintf(fileName, "postProcessing/3LM/P_L%ld",l+1);

                    FILE *fp = fopen(fileName, "w");
                    if(!fp)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("write3LMFields",  error);
                    }
                    else
                    {
                        // write header lines

                        PetscFPrintf(mesh->MESH_COMM, fp, "P given as table of values for each point, where each\n");
                        PetscFPrintf(mesh->MESH_COMM, fp, "line is an array of points in the spanwise direction and each\n");
                        PetscFPrintf(mesh->MESH_COMM, fp, "column is an array of points in the streamwise direction\n\n");

                        // write mesh points
                        for(PetscInt pk=0; pk<lm3->nstw; pk++)
                        {
                            for(PetscInt pi=0; pi<lm3->nspw; pi++)
                            {
                                PetscFPrintf(mesh->MESH_COMM, fp, "%*.4lf ", -15, lev->P[pk][pi]);
                            }
                            // new line
                            PetscFPrintf(mesh->MESH_COMM, fp, "\n");
                        }

                        fclose(fp);
                    }
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findAvgLineIds(acquisition_ *acquisition)
{
    mesh_  *mesh  = acquisition->access->mesh;
    data3LM *lm3  = acquisition->LM3;

    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, pi, pk, r;

    Cmpnts           ***cent;

    PetscReal        ts, te;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // After this function every processor has information about the closest cell
    // to each point of the 3LM mesh. The z of the 3LM points is meaningless since
    // 3LM is a 2D model (x,y) for each one of the three layers. As a consequence,
    // only the i and k indices at each closest cell are meaningful, they define a
    // vertical line of cells in the mesh, on which the average must be done for
    // those cells whose cell center height is between hStart and hEnd.

    PetscPrintf(mesh->MESH_COMM, "   build avg lines - cell interceptions ");

    PetscTime(&ts);

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    std::vector<std::vector<PetscReal>>  gdist(lm3->nstw);
    std::vector<std::vector<PetscReal>>  ldist(lm3->nstw);
    std::vector<std::vector<cellIds>>    lclosestCells(lm3->nstw);

    for(pk=0; pk<lm3->nstw; pk++)
    {
        lclosestCells[pk].resize(lm3->nspw);
        ldist[pk].resize(lm3->nspw);
        gdist[pk].resize(lm3->nspw);

        for(pi=0; pi<lm3->nspw; pi++)
        {
            gdist[pk][pi] = 1e20;
            ldist[pk][pi] = 1e20;

            lclosestCells[pk][pi].i = 0;
            lclosestCells[pk][pi].j = 0;
            lclosestCells[pk][pi].k = 0;
        }
    }

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over 3LM grid points
    for(pk=0; pk<lm3->nstw; pk++)
    {
        for(pi=0; pi<lm3->nspw; pi++)
        {
            // find the cell closest to this 3LM mesh point on this processor
            PetscReal minDistMag = 1e20;
            cellIds   closestCell;

            Cmpnts perturbVec;
                   perturbVec.x = procContrib;
                   perturbVec.y = procContrib;
                   perturbVec.z = procContrib;
            Cmpnts point        = nSub(lm3->points[pk][pi], perturbVec);

            // loop over this processor mesh points
            for(k=lzs; k<lze; k++)
            for(j=lys; j<lye; j++)
            for(i=lxs; i<lxe; i++)
            {
                // compute distance from this cell center
                Cmpnts dist = nSub(point, cent[k][j][i]);

                // compute distance magnitue
                PetscReal distMag = nMag(dist);

                if(distMag < minDistMag)
                {
                    minDistMag = distMag;
                    lclosestCells[pk][pi].i = i;
                    lclosestCells[pk][pi].j = j;
                    lclosestCells[pk][pi].k = k;
                }
            }

            // save this processor minimum
            ldist[pk][pi] = minDistMag;
        }
    }

    // scatter the minimum dist to all processors
    for(pk=0; pk<lm3->nstw; pk++)
    {
        MPI_Allreduce(&ldist[pk][0], &gdist[pk][0], lm3->nspw, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);
    }

    // now compare the lists, where those agree this processor controls this 3LM line
    for(pk=0; pk<lm3->nstw; pk++)
    {
        for(pi=0; pi<lm3->nspw; pi++)
        {
            if(gdist[pk][pi] == ldist[pk][pi])
            {
                // do nothing: idx have been found on this processor
            }
            else
            {
                // zero ids back
                lclosestCells[pk][pi].i = 0;
                lclosestCells[pk][pi].j = 0;
                lclosestCells[pk][pi].k = 0;
            }
        }

        MPI_Allreduce(&(lclosestCells[pk][0]), &(lm3->closestCells[pk][0]), 3*lm3->nspw, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

        std::vector<PetscReal>  ().swap(gdist[pk]);
        std::vector<PetscReal>  ().swap(ldist[pk]);
        std::vector<cellIds>    ().swap(lclosestCells[pk]);
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    PetscTime(&te);

    PetscPrintf(mesh->MESH_COMM, "in %lf s\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode averagingABLInitialize(domain_ *domain)
{
    // ABL averaging is done on the background mesh only
    acquisition_ *acquisition = domain[0].acquisition;

    if(acquisition->isAverageABLActive)
    {

          // check if temperature transport is active
          if(!acquisition->access->flags->isTeqnActive)
          {
             char error[512];
              sprintf(error, "ABL averaging not available without temperature transport (set potentialT to true)");
              fatalErrorInFunction("averagingABLInitialize",  error);
          }

          // check if les is active
          if(!acquisition->access->flags->isLesActive)
          {
             char error[512];
              sprintf(error, "ABL averaging not available without les closure (set les to true)");
              fatalErrorInFunction("averagingABLInitialize",  error);
          }

        // allocate memory
        PetscMalloc(sizeof(dataABL), &(acquisition->statisticsABL));
        dataABL *ablStat = acquisition->statisticsABL;

        clock_  *clock   = acquisition->access->clock;
        mesh_   *mesh    = acquisition->access->mesh;

        ueqn_   *ueqn    = acquisition->access->ueqn;
        teqn_   *teqn    = acquisition->access->teqn;

        DM            da = mesh->da, fda = mesh->fda;
        DMDALocalInfo info = mesh->info;

        PetscInt           xs = info.xs, xe = info.xs + info.xm;
        PetscInt           ys = info.ys, ye = info.ys + info.ym;
        PetscInt           zs = info.zs, ze = info.zs + info.zm;
        PetscInt           mx = info.mx, my = info.my, mz = info.mz;

        PetscMPIInt        rank, nProcs;
        PetscInt           i, j, k, l;
        PetscInt           lxs, lxe, lys, lye, lzs, lze;

        Vec           Coor;

        PetscReal     lxmin, lymin, lzmin,       // local min coordinates
                      lxmax, lymax, lzmax,       // local max coordinates
                      xmin, ymin, zmin,          // global min coordinates
                      xmax, ymax, zmax;          // global max coordinates

        Cmpnts        ***coor, ***cent;          // point and cell center coordinates
        PetscReal     ***aj;

        char          dataLoc[256], fileName[500];

        lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
        lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
        lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

        MPI_Comm_rank(mesh->MESH_COMM, &rank);
        MPI_Comm_size(mesh->MESH_COMM, &nProcs);

        PetscPrintf(mesh->MESH_COMM, "Creating ABL statistics acquisition system...");

        // read averaging parameters
        readDictDouble("control.dat", "-avgABLPeriod", &(ablStat->avgPrd));
        readDictDouble("control.dat", "-avgABLStartTime", &(ablStat->avgStartTime));

        // create directories
        if
        (
            clock->it == clock->itStart && !rank
        )
        {
            word averageFolder = "./postProcessing/" + mesh->meshName + "/averaging";
            errno = 0;
            PetscInt dirRes = mkdir(averageFolder.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory",averageFolder.c_str());
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                errno = 0;
                word timeName = averageFolder + "/" + getTimeName(clock);

                // construct the timeName
                new(&(ablStat->timeName)) std::string{};
                ablStat->timeName = timeName;

                dirRes = mkdir(timeName.c_str(), 0777);
                if(dirRes != 0 && errno != EEXIST)
                {
                   char error[512];
                    sprintf(error, "could not create %s directory\n", timeName.c_str());
                    fatalErrorInFunction("averagingABLInitialize",  error);
                }
                else if(errno == EEXIST)
                {
                    remove_subdirs(mesh->MESH_COMM, timeName.c_str());
                }
            }
        }

        // find cell levels: the wall normal direction is assumed to be the j direction in curv. coords
        // obviously this kind of averaging only works for cartesian meshes, for which the flow
        // has invariant wall parallel statistics.

        DMDAVecGetArray(fda, mesh->lCent, &cent);
        DMDAVecGetArray(da,  mesh->lAj,   &aj);

        // the vertical direction is the j direction in curvilinear coordinates
        PetscInt nLevels = my-2;

        std::vector<PetscReal> lLevels(nLevels);
        std::vector<PetscReal> gLevels(nLevels);
        std::vector<PetscReal> lVolumes(nLevels);
        std::vector<PetscReal> gVolumes(nLevels);
        std::vector<PetscInt>    lCells(nLevels);
        std::vector<PetscInt>    gCells(nLevels);

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
                    lLevels[j-1]  += (cent[k][j][i].z - mesh->bounds.zmin);
                    lVolumes[j-1] += 1 / aj[k][j][i];
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

        // allocate data members
        ablStat->cellLevels     = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->totVolPerLevel = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->totCelPerLevel = (PetscInt* )malloc(sizeof(PetscInt ) * nLevels);
        ablStat->UMean          = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->VMean          = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->WMean          = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->TMean          = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->nutMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        ablStat->uuMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->uvMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->uwMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->vvMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->vwMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wwMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        ablStat->wuuMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wuvMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wuwMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wvvMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wvwMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->wwwMean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        ablStat->R11Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->R12Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->R13Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->R22Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->R23Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->R33Mean        = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        ablStat->TuMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->TvMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->TwMean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        ablStat->q1Mean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->q2Mean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);
        ablStat->q3Mean         = (PetscReal*)malloc(sizeof(PetscReal) * nLevels);

        VecDuplicate(ueqn->Ucat,  &(ablStat->UPrime)); VecSet(ablStat->UPrime, 0.);
        VecDuplicate(teqn->Tmprt, &(ablStat->TPrime)); VecSet(ablStat->TPrime, 0.);

        // set the already available data members
        for(l=0; l<nLevels; l++)
        {
            ablStat->cellLevels[l]     = gLevels[l];
            ablStat->totVolPerLevel[l] = gVolumes[l];
            ablStat->totCelPerLevel[l] = gCells[l];

            //PetscPrintf(mesh->MESH_COMM, "   level %ld: h = %lf, nCells = %ld, totVol = %lf\n", l+1, ablStat->cellLevels[l], ablStat->totCelPerLevel[l], ablStat->totVolPerLevel[l]);
        }

        FILE *f;

        // create the ABL averaging files
        if
        (
            clock->it == clock->itStart && !rank
        )
        {
            // hLevelsCell file
            sprintf(fileName, "%s/hLevelsCell", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                for(l=0; l<nLevels; l++)
                {
                    fprintf(f, "%.5lf\t", ablStat->cellLevels[l]);
                }
                fprintf(f, "\n");

                fclose(f);
            }

            // nu_SGS_mean file
            sprintf(fileName, "%s/nu_SGS_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // q1_mean file
            sprintf(fileName, "%s/q1_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // q2_mean file
            sprintf(fileName, "%s/q2_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // q3_mean file
            sprintf(fileName, "%s/q3_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R11_mean file
            sprintf(fileName, "%s/R11_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R12_mean file
            sprintf(fileName, "%s/R12_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R13_mean file
            sprintf(fileName, "%s/R13_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R22_mean file
            sprintf(fileName, "%s/R22_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R23_mean file
            sprintf(fileName, "%s/R23_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // R33_mean file
            sprintf(fileName, "%s/R33_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // T_mean file
            sprintf(fileName, "%s/T_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // Tu_mean file
            sprintf(fileName, "%s/Tu_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // Tv_mean file
            sprintf(fileName, "%s/Tv_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // Tw_mean file
            sprintf(fileName, "%s/Tw_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // U_mean file
            sprintf(fileName, "%s/U_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // uu_mean file
            sprintf(fileName, "%s/uu_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // uv_mean file
            sprintf(fileName, "%s/uv_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // uw_mean file
            sprintf(fileName, "%s/uw_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // V_mean file
            sprintf(fileName, "%s/V_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // vv_mean file
            sprintf(fileName, "%s/vv_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // vw_mean file
            sprintf(fileName, "%s/vw_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // W_mean file
            sprintf(fileName, "%s/W_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // wuu_mean file
            sprintf(fileName, "%s/wuu_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // wuv_mean file
            sprintf(fileName, "%s/wuv_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // wuw_mean file
            sprintf(fileName, "%s/wuw_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // wvv_mean file
            sprintf(fileName, "%s/wvv_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // wvw_mean file
            sprintf(fileName, "%s/wvw_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // ww_mean file
            sprintf(fileName, "%s/ww_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }

            // www_mean file
            sprintf(fileName, "%s/www_mean", ablStat->timeName.c_str());
            f = fopen(fileName, "w");

            if(!f)
            {
                char error[530];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("averagingABLInitialize",  error);
            }
            else
            {
                fclose(f);
            }
        }

        PetscPrintf(mesh->MESH_COMM, "done\n\n");
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeAveragingABL(domain_ *domain)
{
    // ABL averaging is done on the background mesh only
    acquisition_ *acquisition = domain[0].acquisition;

    if(acquisition->isAverageABLActive)
    {
        PetscInt    accumulateAvg   = 0;

        clock_  *clock   = acquisition->access->clock;
        dataABL *ablStat = acquisition->statisticsABL;

        PetscReal startTimeAvg    = ablStat->avgStartTime;
        PetscReal timeIntervalAvg = ablStat->avgPrd;

        // check if must accumulate averaged fields
        if
        (
            clock->time >= startTimeAvg &&
            (clock->time - startTimeAvg ) / timeIntervalAvg -
            std::floor((clock->time - startTimeAvg) / timeIntervalAvg) < 1e-10
        )
        {
            accumulateAvg = 1;
        }

        if(accumulateAvg)
        {
            mesh_   *mesh    = acquisition->access->mesh;
            ueqn_   *ueqn    = acquisition->access->ueqn;
            teqn_   *teqn    = acquisition->access->teqn;
            les_    *les     = acquisition->access->les;

            DM            da = mesh->da, fda = mesh->fda;
            DMDALocalInfo info = mesh->info;

            PetscInt           xs = info.xs, xe = info.xs + info.xm;
            PetscInt           ys = info.ys, ye = info.ys + info.ym;
            PetscInt           zs = info.zs, ze = info.zs + info.zm;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt        rank, nProcs;
            PetscInt           i, j, k, l;
            PetscInt           lxs, lxe, lys, lye, lzs, lze;

            Cmpnts        ***ucat, ***uprime;           // cartesian vel. and fluct. part
            Cmpnts        ***csi, ***eta, ***zet;
            PetscReal     ***tmprt, ***tprime, ***nut;  // potential temp. and fluct. part and turb. visc.
            PetscReal     ***aj, ***nvert;

            word          fileName;

            PetscReal        te, ts;

            PetscTime(&ts);

            lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
            lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
            lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

            MPI_Comm_rank(mesh->MESH_COMM, &rank);
            MPI_Comm_size(mesh->MESH_COMM, &nProcs);

            DMDAVecGetArray(da,  mesh->lAj,       &aj);
            DMDAVecGetArray(da,  mesh->lNvert,    &nvert);
            DMDAVecGetArray(fda, mesh->lCsi,      &csi);
            DMDAVecGetArray(fda, mesh->lEta,      &eta);
            DMDAVecGetArray(fda, mesh->lZet,      &zet);
            DMDAVecGetArray(fda, ueqn->lUcat,     &ucat);
            DMDAVecGetArray(da,  teqn->lTmprt,    &tmprt);
            DMDAVecGetArray(da,  les->lNu_t,      &nut);
            DMDAVecGetArray(fda, ablStat->UPrime, &uprime);
            DMDAVecGetArray(da,  ablStat->TPrime, &tprime);

            PetscInt nLevels = my-2;

            // compute the average velocity fields
            std::vector<Cmpnts> lVelocity(nLevels);
            std::vector<Cmpnts> gVelocity(nLevels);
            std::vector<PetscReal> lTemperature(nLevels);
            std::vector<PetscReal> gTemperature(nLevels);
            std::vector<PetscReal> lnut(nLevels);
            std::vector<PetscReal> gnut(nLevels);

            for(l=0; l<nLevels; l++)
            {
                lVelocity[l].x = lVelocity[l].y = lVelocity[l].z = 0.0;
                gVelocity[l].x = gVelocity[l].y = gVelocity[l].z = 0.0;
                lTemperature[l] = 0.0;
                gTemperature[l] = 0.0;
                lnut[l] = 0.0;
                gnut[l] = 0.0;
            }

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        lVelocity[j-1].x  += ucat[k][j][i].x / aj[k][j][i];
                        lVelocity[j-1].y  += ucat[k][j][i].y / aj[k][j][i];
                        lVelocity[j-1].z  += ucat[k][j][i].z / aj[k][j][i];
                        lTemperature[j-1] += tmprt[k][j][i]  / aj[k][j][i];
                        lnut[j-1]         += nut[k][j][i]    / aj[k][j][i];
                    }
                }
            }

            MPI_Allreduce(&lVelocity[0], &gVelocity[0], 3*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lTemperature[0], &gTemperature[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnut[0], &gnut[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

            for(l=0; l<nLevels; l++)
            {
                PetscReal totVolPerLevel = ablStat->totVolPerLevel[l];

                ablStat->UMean[l]   = gVelocity[l].x / totVolPerLevel;
                ablStat->VMean[l]   = gVelocity[l].y / totVolPerLevel;
                ablStat->WMean[l]   = gVelocity[l].z / totVolPerLevel;
                ablStat->TMean[l]   = gTemperature[l] / totVolPerLevel;
                ablStat->nutMean[l] = gnut[l] / totVolPerLevel;
            }

            // compute fluctuating fields
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        uprime[k][j][i].x = ucat[k][j][i].x - ablStat->UMean[j-1];
                        uprime[k][j][i].y = ucat[k][j][i].y - ablStat->VMean[j-1];
                        uprime[k][j][i].z = ucat[k][j][i].z - ablStat->WMean[j-1];
                        tprime[k][j][i]   = tmprt[k][j][i]  - ablStat->TMean[j-1];
                    }
                }
            }

            // compute the mean of the fluctuating fields to output
            std::vector<Cmpnts>     lTU(nLevels);
            std::vector<Cmpnts>     gTU(nLevels);
            std::vector<symmTensor> lS(nLevels);
            std::vector<symmTensor> gS(nLevels);
            std::vector<symmTensor> lSw(nLevels);
            std::vector<symmTensor> gSw(nLevels);
            std::vector<symmTensor> lR(nLevels);
            std::vector<symmTensor> gR(nLevels);

            for(l=0; l<nLevels; l++)
            {
                // set to zero
                lTU[l].x  = lTU[l].y  = lTU[l].z  = 0.0;
                gTU[l].x  = gTU[l].y  = gTU[l].z  = 0.0;
                lS[l].xx  = lS[l].xy  = lS[l].xz  = lS[l].yy = lS[l].yz = lS[l].zz = 0.0;
                gS[l].xx  = gS[l].xy  = gS[l].xz  = gS[l].yy = gS[l].yz = gS[l].zz = 0.0;
                lSw[l].xx = lSw[l].xy = lSw[l].xz = lSw[l].yy = lSw[l].yz = lSw[l].zz = 0.0;
                gSw[l].xx = gSw[l].xy = gSw[l].xz = gSw[l].yy = gSw[l].yz = gSw[l].zz = 0.0;
                lR[l].xx = lR[l].xy = lR[l].xz = lR[l].yy = lR[l].yz = lR[l].zz = 0.0;
                gR[l].xx = gR[l].xy = gR[l].xz = gR[l].yy = gR[l].yz = gR[l].zz = 0.0;
            }

            PetscReal volCell;
            PetscReal tprimeCell;
            Cmpnts uprimeCell;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        // pre-set base variables for speed
                        PetscReal dudc, dvdc, dwdc,
                               dude, dvde, dwde,
                               dudz, dvdz, dwdz;
                        PetscReal du_dx, du_dy, du_dz,
                               dv_dx, dv_dy, dv_dz,
                               dw_dx, dw_dy, dw_dz;
                        PetscReal ajc  = aj[k][j][i];

                        PetscReal csi0 = csi[k][j][i].x,
                               csi1 = csi[k][j][i].y,
                               csi2 = csi[k][j][i].z;
                        PetscReal eta0 = eta[k][j][i].x,
                               eta1 = eta[k][j][i].y,
                               eta2 = eta[k][j][i].z;
                        PetscReal zet0 = zet[k][j][i].x,
                               zet1 = zet[k][j][i].y,
                               zet2 = zet[k][j][i].z;

                        // pre-set variables for speed
                        volCell      = 1.0 / ajc;
                        tprimeCell   = tprime[k][j][i];
                        uprimeCell.x = uprime[k][j][i].x;
                        uprimeCell.y = uprime[k][j][i].y;
                        uprimeCell.z = uprime[k][j][i].z;

                        Compute_du_center
                        (
                            mesh,
                            i, j, k, mx, my, mz, ucat, nvert, &dudc,
                            &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                        );

                        Compute_du_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                            zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                            dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                            &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                        );

                        // compute statistics multipling by the volume
                        lS[j-1].xx  += uprimeCell.x * uprimeCell.x * volCell;
                        lS[j-1].xy  += uprimeCell.x * uprimeCell.y * volCell;
                        lS[j-1].xz  += uprimeCell.x * uprimeCell.z * volCell;
                        lS[j-1].yy  += uprimeCell.y * uprimeCell.y * volCell;
                        lS[j-1].yz  += uprimeCell.y * uprimeCell.z * volCell;
                        lS[j-1].zz  += uprimeCell.z * uprimeCell.z * volCell;

                        lSw[j-1].xx += uprimeCell.z * uprimeCell.x * uprimeCell.x * volCell;
                        lSw[j-1].xy += uprimeCell.z * uprimeCell.x * uprimeCell.y * volCell;
                        lSw[j-1].xz += uprimeCell.z * uprimeCell.x * uprimeCell.z * volCell;
                        lSw[j-1].yy += uprimeCell.z * uprimeCell.y * uprimeCell.y * volCell;
                        lSw[j-1].yz += uprimeCell.z * uprimeCell.y * uprimeCell.z * volCell;
                        lSw[j-1].zz += uprimeCell.z * uprimeCell.z * uprimeCell.z * volCell;

                        lTU[j-1].x  += tprimeCell * uprimeCell.x * volCell;
                        lTU[j-1].y  += tprimeCell * uprimeCell.y * volCell;
                        lTU[j-1].z  += tprimeCell * uprimeCell.z * volCell;

                        PetscReal tke = 0.5 * (uprimeCell.x * uprimeCell.x + uprimeCell.y * uprimeCell.y + uprimeCell.z * uprimeCell.z);

                        lR[j-1].xx += (2.0/3.0 * tke - 2.0 * nut[k][j][i] * du_dx) * volCell;
                        lR[j-1].xy += ( - nut[k][j][i] * (dv_dx + du_dy)) * volCell;
                        lR[j-1].xz += ( - nut[k][j][i] * (dw_dx + du_dz)) * volCell;
                        lR[j-1].yy += (2.0/3.0 * tke - 2.0 * nut[k][j][i] * dv_dy) * volCell;
                        lR[j-1].yz += ( - nut[k][j][i] * (dv_dz + dw_dy)) * volCell;
                        lR[j-1].zz += (2.0/3.0 * tke - 2.0 * nut[k][j][i] * du_dx) * volCell;
                    }
                }
            }

            // sum statistics over processors
            MPI_Allreduce(&lS[0],  &gS[0],  6*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lSw[0], &gSw[0], 6*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lR[0],  &gR[0],  6*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lTU[0], &gTU[0], 3*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

            // store statistics after dividing for the total volume per level
            for (l=0; l<nLevels; l++)
            {
                PetscReal totVolPerLevel = ablStat->totVolPerLevel[l];

                ablStat->uuMean[l] = gS[l].xx / totVolPerLevel;
                ablStat->uvMean[l] = gS[l].xy / totVolPerLevel;
                ablStat->uwMean[l] = gS[l].xz / totVolPerLevel;
                ablStat->vvMean[l] = gS[l].yy / totVolPerLevel;
                ablStat->vwMean[l] = gS[l].yz / totVolPerLevel;
                ablStat->wwMean[l] = gS[l].zz / totVolPerLevel;

                ablStat->wuuMean[l] = gSw[l].xx / totVolPerLevel;
                ablStat->wuvMean[l] = gSw[l].xy / totVolPerLevel;
                ablStat->wuwMean[l] = gSw[l].xz / totVolPerLevel;
                ablStat->wvvMean[l] = gSw[l].yy / totVolPerLevel;
                ablStat->wvwMean[l] = gSw[l].yz / totVolPerLevel;
                ablStat->wwwMean[l] = gSw[l].zz / totVolPerLevel;

                ablStat->TuMean[l] = gTU[l].x / totVolPerLevel;
                ablStat->TvMean[l] = gTU[l].y / totVolPerLevel;
                ablStat->TwMean[l] = gTU[l].z / totVolPerLevel;

                // set to zero remaining statistics (not yet implemented)
                ablStat->R11Mean[l] = gR[j].xx / totVolPerLevel;
                ablStat->R12Mean[l] = gR[j].xy / totVolPerLevel;
                ablStat->R13Mean[l] = gR[j].xz / totVolPerLevel;
                ablStat->R22Mean[l] = gR[j].yy / totVolPerLevel;
                ablStat->R23Mean[l] = gR[j].yz / totVolPerLevel;
                ablStat->R33Mean[l] = gR[j].zz / totVolPerLevel;

                ablStat->q1Mean[l] = 0.0 / totVolPerLevel;
                ablStat->q2Mean[l] = 0.0 / totVolPerLevel;
                ablStat->q3Mean[l] = 0.0 / totVolPerLevel;
            }

            DMDAVecRestoreArray(da,  mesh->lAj,       &aj);
            DMDAVecRestoreArray(da,  mesh->lNvert,    &nvert);
            DMDAVecRestoreArray(fda, mesh->lCsi,      &csi);
            DMDAVecRestoreArray(fda, mesh->lEta,      &eta);
            DMDAVecRestoreArray(fda, mesh->lZet,      &zet);
            DMDAVecRestoreArray(fda, ueqn->lUcat,     &ucat);
            DMDAVecRestoreArray(da,  teqn->lTmprt,    &tmprt);
            DMDAVecRestoreArray(da,  les->lNu_t,      &nut);
            DMDAVecRestoreArray(fda, ablStat->UPrime, &uprime);
            DMDAVecRestoreArray(da,  ablStat->TPrime, &tprime);

            // clean the allocated vectors
            std::vector<Cmpnts> ().swap(lVelocity);
            std::vector<Cmpnts> ().swap(gVelocity);
            std::vector<PetscReal> ().swap(lTemperature);
            std::vector<PetscReal> ().swap(gTemperature);
            std::vector<PetscReal> ().swap(lnut);
            std::vector<PetscReal> ().swap(gnut);

            std::vector<Cmpnts>     ().swap(lTU);
            std::vector<Cmpnts>     ().swap(gTU);
            std::vector<symmTensor> ().swap(lS);
            std::vector<symmTensor> ().swap(gS);
            std::vector<symmTensor> ().swap(lSw);
            std::vector<symmTensor> ().swap(gSw);
            std::vector<symmTensor> ().swap(lR);
            std::vector<symmTensor> ().swap(gR);

            // write statistics to files

            FILE *f;

            // create the ABL averaging files
            if
            (
                !rank
            )
            {
                // nu_SGS_mean file
                fileName = ablStat->timeName + "/nu_SGS_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->nutMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // q1_mean file
                fileName = ablStat->timeName + "/q1_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->q1Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // q2_mean file
                fileName = ablStat->timeName + "/q2_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->q2Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // q3_mean file
                fileName = ablStat->timeName + "/q3_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->q3Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R11_mean file
                fileName = ablStat->timeName + "/R11_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R11Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R12_mean file
                fileName = ablStat->timeName + "/R12_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R12Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R13_mean file
                fileName = ablStat->timeName + "/R13_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R13Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R22_mean file
                fileName = ablStat->timeName + "/R22_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R22Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R23_mean file
                fileName = ablStat->timeName + "/R23_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R23Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // R33_mean file
                fileName = ablStat->timeName + "/R33_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->R33Mean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // T_mean file
                fileName = ablStat->timeName + "/T_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->TMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // Tu_mean file
                fileName = ablStat->timeName + "/Tu_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->TuMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // Tv_mean file
                fileName = ablStat->timeName + "/Tv_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->TvMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // Tw_mean file
                fileName = ablStat->timeName + "/Tw_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->TwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // U_mean file
                fileName = ablStat->timeName + "/U_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->UMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // uu_mean file
                fileName = ablStat->timeName + "/uu_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->uuMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // uv_mean file
                fileName = ablStat->timeName + "/uv_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->uvMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // uw_mean file
                fileName = ablStat->timeName + "/uw_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->uwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // V_mean file
                fileName = ablStat->timeName + "/V_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->VMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // vv_mean file
                fileName = ablStat->timeName + "/vv_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->vvMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // vw_mean file
                fileName = ablStat->timeName + "/vw_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->vwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // W_mean file
                fileName = ablStat->timeName + "/W_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->WMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // wuu_mean file
                fileName = ablStat->timeName + "/wuu_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wuuMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // wuv_mean file
                fileName = ablStat->timeName + "/wuv_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wuvMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // wuw_mean file
                fileName = ablStat->timeName + "/wuw_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wuwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // wvv_mean file
                fileName = ablStat->timeName + "/wvv_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wvvMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // wvw_mean file
                fileName = ablStat->timeName + "/wvw_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wvwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // ww_mean file
                fileName = ablStat->timeName + "/ww_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }

                // www_mean file
                fileName = ablStat->timeName + "/www_mean";
                f = fopen(fileName.c_str(), "a");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName.c_str());
                    fatalErrorInFunction("writeAveragingABL",  error);
                }
                else
                {
                    fprintf(f, "%.5lf\t", clock->time);
                    fprintf(f, "%.5lf\t", clock->dt);
                    for(l=0; l<nLevels; l++)
                    {
                        fprintf(f, "%.5lf\t", ablStat->wwwMean[l]);
                    }
                    fprintf(f, "\n");
                    fclose(f);
                }
            }

            PetscTime(&te);

            PetscPrintf(mesh->MESH_COMM, "Averaged ABL data in %lf s\n", te-ts);
        }
    }

    return(0);
}
