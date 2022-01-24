//! \file  io.c
//! \brief Contains i/o operations function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"

//***************************************************************************************************************//

PetscErrorCode InitializeIO(io_ *io)
{
    mesh_          *mesh = io->access->mesh;
    flags_         *flags = io->access->flags;
    PetscMPIInt    rank, nProcs;

    char          dataLoc[256], fileName[500];

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    MPI_Comm_size(mesh->MESH_COMM, &nProcs);

    PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

    // read writing settings (mandatory entries)
    readDictWord  ("control.dat", "-intervalType", &(io->intervalType));
    readDictDouble("control.dat", "-timeInterval", &(io->timeInterval));

    // read purge write (optional)
    io->purgeWrite     = 0;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-purgeWrite", &(io->purgeWrite), PETSC_NULL);

    // read averaging settings (optional entries)
    io->averaging      = 0;
    io->phaseAveraging = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averaging", &(io->averaging), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-phaseAveraging", &(io->phaseAveraging), PETSC_NULL);

    // read q-criterion flag
    io->qCrit          = 0;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQ", &(io->qCrit), PETSC_NULL);

    // read l2 criterion flag
    io->l2Crit          = 0;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeL2", &(io->l2Crit), PETSC_NULL);

    // allocate memory and initialize averaged fields
    if(io->averaging)
    {
        PetscPrintf(mesh->MESH_COMM, "Averaging: yes\n\n");

        readDictDouble("control.dat", "-avgPeriod", &(io->avgPrd));
        readDictDouble("control.dat", "-avgStartTime", &(io->avgStartTime));

        // initialize snapshot weighting
        io->avgWeight = 0;

    }
    // allocate memory and initialize phase averaged fields
    if(io->phaseAveraging)
    {
        PetscPrintf(mesh->MESH_COMM, "Phase averaging: yes\n\n");

        readDictDouble("control.dat", "-phaseAvgPeriod", &(io->phAvgPrd));
        readDictDouble("control.dat", "-phaseAvgStartTime", &(io->phAvgStartTime));

        // initialize snapshot weighting
        io->pAvgWeight = 0;
    }

    MPI_Barrier(mesh->MESH_COMM);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readFields(domain_ *domain, PetscReal timeValue)
{
    PetscViewer viewer;
    PetscInt    N;
    ueqn_       *ueqn = domain->ueqn;
    mesh_       *mesh = domain->mesh;
    peqn_       *peqn = domain->peqn;
    teqn_       *teqn = domain->teqn;
    les_        *les  = domain->les;
    clock_      *clock= domain->clock;
    io_         *io   = domain->io;
    acquisition_ *acquisition = domain->acquisition;
    word        fileName, field, location;

    VecGetSize(ueqn->Ucont, &N);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << timeValue;

    location = "./fields/" + mesh->meshName + "/" + stream.str();

    // read contravariant flux
    PetscPrintf(mesh->MESH_COMM, "Reading V...\n");
    field = "/V";
    fileName = location + field;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(ueqn->Ucont, viewer);
    PetscViewerDestroy(&viewer);

    // read cartesian velocity
    PetscPrintf(mesh->MESH_COMM, "Reading U...\n");
    field = "/U";
    fileName = location + field;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(ueqn->Ucat,viewer);
    PetscViewerDestroy(&viewer);

    // read pressure
    PetscPrintf(mesh->MESH_COMM, "Reading p...\n");
    field = "/p";
    fileName = location + field;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(peqn->P,viewer);
    PetscViewerDestroy(&viewer);

    // read temperature
    if(domain->flags.isTeqnActive)
    {
        PetscPrintf(mesh->MESH_COMM, "Reading T...\n");
        field = "/T";
        fileName = location + field;
        PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
        VecLoad(teqn->Tmprt, viewer);
        PetscViewerDestroy(&viewer);
    }

    //read ibm field
    PetscPrintf(mesh->MESH_COMM, "Reading nv...\n");
    field = "/nv";
    fileName = location + field;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(mesh->Nvert,viewer);
    PetscViewerDestroy(&viewer);

    if(domain->flags.isLesActive)
    {
        Vec Cs, Nut;
        FILE *fp;
        VecDuplicate(peqn->P, &Cs);
        VecDuplicate(peqn->P, &Nut);

        PetscPrintf(mesh->MESH_COMM, "Reading cs...\n");
        field = "/cs";
        fileName = location + field;

        fp=fopen(fileName.c_str(), "r");

        if(fp==NULL)
        {
            char warning[256];
            sprintf(warning, "Cannot open %s, setting cs to 0 and continuing the computation...", fileName.c_str());
            warningInFunction("readFields",  warning);

            VecSet(Cs, 0.0);
        }
        else
        {
            fclose(fp);

            MPI_Barrier(mesh->MESH_COMM);
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(Cs,viewer);
            PetscViewerDestroy(&viewer);
        }

        DMGlobalToLocalBegin(mesh->da, Cs, INSERT_VALUES, les->lCs);
        DMGlobalToLocalEnd(mesh->da, Cs, INSERT_VALUES, les->lCs);

        PetscPrintf(mesh->MESH_COMM, "Reading nut...\n");
        field = "/nut";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp==NULL)
        {
            char warning[256];
            sprintf(warning, "Cannot open %s, setting nut to 0 and continuing the computation...", fileName.c_str());
            warningInFunction("readFields",  warning);

            VecSet(Nut, 0.0);
        }
        else
        {
            fclose(fp);

            MPI_Barrier(mesh->MESH_COMM);
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(Nut,viewer);
            PetscViewerDestroy(&viewer);
        }

        DMGlobalToLocalBegin(mesh->da, Nut, INSERT_VALUES, les->lNu_t);
        DMGlobalToLocalEnd(mesh->da, Nut, INSERT_VALUES, les->lNu_t);

        // destroy vectors
        VecDestroy(&Cs);
        VecDestroy(&Nut);
    }

    // scatter fields
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    VecCopy(ueqn->Ucont, ueqn->Ucont_o);

    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    DMGlobalToLocalBegin(mesh->da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd(mesh->da, peqn->P, INSERT_VALUES, peqn->lP);

    if(domain->flags.isTeqnActive)
    {
        DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        VecCopy(teqn->Tmprt, teqn->Tmprt_o);
    }

    DMGlobalToLocalBegin(mesh->da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(mesh->da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    DMGlobalToLocalBegin(mesh->da, mesh->Nvert_o, INSERT_VALUES, mesh->lNvert_o);
    DMGlobalToLocalEnd(mesh->da, mesh->Nvert_o, INSERT_VALUES, mesh->lNvert_o);

    // read averaged fields
    if(io->averaging)
    {
        // open file to check the existence, then read it with PETSc
        FILE *fp;

        // read avgU
        field = "/avgU";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading avgU...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->avgU,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        // read avgP
        field = "/avgP";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading avgP...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->avgP,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        // read avgUU
        field = "/avgUU";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading avgUU...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->avgUU,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        if(les)
        {
            // read avgNut
            field = "/avgNut";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading avgNut...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->avgNut,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);

            // read avgCs
            field = "/avgCs";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading avgCs...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->avgCs,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);
        }

        if(io->averaging > 1)
        {
            // read avgOmega
            field = "/avgOmega";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading avgOmega...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->avgOmega,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);

            if(io->averaging > 2)
            {
                // read avgP2
                field = "/avgP2";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading avgPsq...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->avgP2,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read avgOmegaOmega
                field = "/avgOmegaOmega";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading avgOmegaOmega...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->avgOmegaOmega,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read avgUdotGradP
                field = "/avgUdotGradP";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading avgUdotGradP...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->avgUdotGradP,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read avgMagGradU
                field = "/avgMagGradU";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading avgMagGradU...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->avgMagGradU,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read avgMagUU
                field = "/avgMagUU";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading avgMagUU...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->avgMagUU,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);
            }
        }
    }

    // read Q-Criterion
    if(io->qCrit)
    {
        // open file to check the existence, then read it with PETSc
        FILE *fp;

        // read pAvgU
        field = "/Q";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading Q...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->Q,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);
    }

    // read phase averaged fields
    if(io->phaseAveraging)
    {
        // open file to check the existence, then read it with PETSc
        FILE *fp;

        // read pAvgU
        field = "/phAvgU";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading phAvgU...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->pAvgU,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        // read pAvgP
        field = "/phAvgP";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading phAvgP...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->pAvgP,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        // read pAvgUU
        field = "/phAvgUU";
        fileName = location + field;
        fp=fopen(fileName.c_str(), "r");

        if(fp!=NULL)
        {
            fclose(fp);

            PetscPrintf(mesh->MESH_COMM, "Reading phAvgUU...\n");
            PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
            VecLoad(acquisition->fields->pAvgUU,viewer);
            PetscViewerDestroy(&viewer);
        }
        MPI_Barrier(mesh->MESH_COMM);

        if(domain->flags.isLesActive)
        {
            // read pAvgNut
            field = "/phAvgNut";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading phAvgNut...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->pAvgNut,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);

            // read pAvgCs
            field = "/phAvgCs";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading phAvgCs...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->pAvgCs,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);
        }

        if(io->phaseAveraging > 1)
        {
            // read pAvgOmega
            field = "/phAvgOmega";
            fileName = location + field;
            fp=fopen(fileName.c_str(), "r");

            if(fp!=NULL)
            {
                fclose(fp);

                PetscPrintf(mesh->MESH_COMM, "Reading phAvgOmega...\n");
                PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                VecLoad(acquisition->fields->pAvgOmega,viewer);
                PetscViewerDestroy(&viewer);
            }
            MPI_Barrier(mesh->MESH_COMM);

            if(io->phaseAveraging > 2)
            {
                // read pAvgP2
                field = "/phAvgP2";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading phAvgPsq...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->pAvgP2,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read pAvgOmegaOmega
                field = "/phAvgOmegaOmega";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading phAvgOmegaOmega...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->pAvgOmegaOmega,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read pAvgUdotGradP
                field = "/phAvgUdotGradP";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading phAvgUdotGradP...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->pAvgUdotGradP,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read pAvgMagGradU
                field = "/phAvgMagGradU";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading phAvgMagGradU...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->pAvgMagGradU,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);

                // read pAvgMagUU
                field = "/phAvgMagUU";
                fileName = location + field;
                fp=fopen(fileName.c_str(), "r");

                if(fp!=NULL)
                {
                    fclose(fp);

                    PetscPrintf(mesh->MESH_COMM, "Reading phAvgMagUU...\n");
                    PetscViewerBinaryOpen(mesh->MESH_COMM, fileName.c_str(), FILE_MODE_READ, &viewer);
                    VecLoad(acquisition->fields->pAvgMagUU,viewer);
                    PetscViewerDestroy(&viewer);
                }
                MPI_Barrier(mesh->MESH_COMM);
            }
        }
    }

    PetscPrintf(mesh->MESH_COMM, "\n");

    return(0);
}

//***************************************************************************************************************//

void writeBinaryField(MPI_Comm comm, Vec V, const char *file)
{
    PetscViewer viewer;
    char s[256];

    PetscMPIInt rank;  MPI_Comm_rank(comm, &rank);

    PetscViewerCreate(comm, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERBINARY);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerFileSetName(viewer, file);
    VecView(V, viewer);
    PetscViewerDestroy(&viewer);

    if(!rank)
    {
        sprintf(s, "%s.info", file);
        unlink(s);
    }
}

//***************************************************************************************************************//

PetscErrorCode writeFields(io_ *io)
{
    mesh_       *mesh  = io->access->mesh;
    ueqn_       *ueqn  = io->access->ueqn;
    peqn_       *peqn  = io->access->peqn;
    teqn_       *teqn  = io->access->teqn;
    les_        *les   = io->access->les;
    clock_      *clock = io->access->clock;
    flags_      *flags = io->access->flags;
    acquisition_ *acquisition = io->access->acquisition;

    PetscViewer viewer;
    word        timeName, fieldName;
    word        fieldsDir = "./fields/" + mesh->meshName;

    PetscMPIInt         rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // check and write
    if(io->runTimeWrite)
    {
        PetscPrintf(mesh->MESH_COMM, "Writing fields for time %lf\n", clock->time);

        // current time name path
        timeName = fieldsDir + "/"  + getTimeName(clock);

        // write contravariant velocity
        fieldName = timeName + "/V";
        writeBinaryField(mesh->MESH_COMM, ueqn->Ucont, fieldName.c_str());
        MPI_Barrier(mesh->MESH_COMM);

        // write cartesian velocity
        fieldName = timeName + "/U";
        writeBinaryField(mesh->MESH_COMM, ueqn->Ucat, fieldName.c_str());
        MPI_Barrier(mesh->MESH_COMM);

        // write pressure
        fieldName = timeName + "/p";
        writeBinaryField(mesh->MESH_COMM, peqn->P, fieldName.c_str());
        MPI_Barrier(mesh->MESH_COMM);

        // write nvert
        fieldName = timeName + "/nv";
        writeBinaryField(mesh->MESH_COMM, mesh->Nvert, fieldName.c_str());
        MPI_Barrier(mesh->MESH_COMM);

        // write temperature
        if(flags->isTeqnActive)
        {
            fieldName = timeName + "/T";
            writeBinaryField(mesh->MESH_COMM, teqn->Tmprt, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);
        }

        if(flags->isLesActive)
        {
            Vec Cs, Nut;

            VecDuplicate(peqn->P, &Cs);
            VecDuplicate(peqn->P, &Nut);
            VecSet(Cs, 0.);
            VecSet(Nut, 0.);

            // we cannot use INSERT_VALUE mode from local to global because it is not capable
            // of handling more than one periodic boundary, PETSc will throw an error.
            // So we set it explicitly

            DM            da = mesh->da, fda = mesh->fda;
            DMDALocalInfo info = mesh->info;
            PetscInt      xs = info.xs, xe = info.xs + info.xm;
            PetscInt      ys = info.ys, ye = info.ys + info.ym;
            PetscInt      zs = info.zs, ze = info.zs + info.zm;
            PetscInt      mx = info.mx, my = info.my, mz = info.mz;

            PetscInt      lxs, lxe, lys, lye, lzs, lze;
            PetscInt      i, j, k;

            PetscReal     ***cs, ***lcs, ***nut, ***lnut;

            lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
            lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
            lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

            DMDAVecGetArray(da, Cs, &cs);
            DMDAVecGetArray(da, les->lCs, &lcs);
            DMDAVecGetArray(da, Nut, &nut);
            DMDAVecGetArray(da, les->lNu_t, &lnut);

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        cs[k][j][i] = lcs[k][j][i];
                        nut[k][j][i] = lnut[k][j][i];
                    }
                }
            }

            DMDAVecRestoreArray(da, Cs, &cs);
            DMDAVecRestoreArray(da, les->lCs, &lcs);
            DMDAVecRestoreArray(da, Nut, &nut);
            DMDAVecRestoreArray(da, les->lNu_t, &lnut);

            fieldName = timeName + "/cs";
            writeBinaryField(mesh->MESH_COMM, Cs, fieldName.c_str());

            fieldName = timeName + "/nut";
            writeBinaryField(mesh->MESH_COMM, Nut, fieldName.c_str());

            MPI_Barrier(mesh->MESH_COMM);
            VecDestroy(&Cs);
            VecDestroy(&Nut);
        }

        // write Q-Criterion
        if(io->qCrit)
        {
            computeQ(acquisition);
            fieldName = timeName + "/Q";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->Q, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);
        }

        if(io->averaging)
        {
            // write avgU
            fieldName = timeName + "/avgU";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgU, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            // write avgP
            fieldName = timeName + "/avgP";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgP, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            // write avgUU
            fieldName = timeName + "/avgUU";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgUU, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            if(flags->isLesActive)
            {
                // write avgNut
                fieldName = timeName + "/avgNut";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgNut, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);

                // write avgCs
                fieldName = timeName + "/avgCs";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgCs, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);
            }

            if(io->averaging > 1)
            {
                // write avgOmega
                fieldName = timeName + "/avgOmega";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgOmega, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);

                if(io->averaging > 2)
                {
                    // write avgP2
                    fieldName = timeName + "/avgP2";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgP2, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write avgOmegaOmega
                    fieldName = timeName + "/avgOmegaOmega";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgOmegaOmega, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write avgUdotGradP
                    fieldName = timeName + "/avgUdotGradP";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgUdotGradP, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write avgMagGradU
                    fieldName = timeName + "/avgMagGradU";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgMagGradU, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write avgMagUU
                    fieldName = timeName + "/avgMagUU";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->avgMagUU, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);
                }
            }
        }

        if(io->phaseAveraging)
        {
            // write pAvgU
            fieldName = timeName + "/phAvgU";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgU, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            // write pAvgP
            fieldName = timeName + "/phAvgP";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgP, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            // write pAvgUU
            fieldName = timeName + "/phAvgUU";
            writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgUU, fieldName.c_str());
            MPI_Barrier(mesh->MESH_COMM);

            if(flags->isLesActive)
            {
                // write pAvgNut
                fieldName = timeName + "/phAvgNut";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgNut, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);

                // write pAvgCs
                fieldName = timeName + "/phAvgCs";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgCs, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);
            }

            if(io->phaseAveraging > 1)
            {
                // write pAvgOmega
                fieldName = timeName + "/phAvgOmega";
                writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgOmega, fieldName.c_str());
                MPI_Barrier(mesh->MESH_COMM);

                if(io->phaseAveraging > 2)
                {
                    // write pAvgP2
                    fieldName = timeName + "/phAvgP2";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgP2, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write pAvgOmegaOmega
                    fieldName = timeName + "/phAvgOmegaOmega";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgOmegaOmega, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write pAvgUdotGradP
                    fieldName = timeName + "/phAvgUdotGradP";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgUdotGradP, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write pAvgMagGradU
                    fieldName = timeName + "/phAvgMagGradU";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgMagGradU, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);

                    // write pAvgMagUU
                    fieldName = timeName + "/phAvgMagUU";
                    writeBinaryField(mesh->MESH_COMM, acquisition->fields->pAvgMagUU, fieldName.c_str());
                    MPI_Barrier(mesh->MESH_COMM);
                }
            }
        }

        // delete all other folders if purge is active (recommended for big cases)
        if(io->purgeWrite)
        {
            word keep, writeDir = "./fields/" + mesh->meshName;

            // current time name path
            if(clock->it == clock->itStart)
            {
                keep = getStartTimeName(clock);
            }
            else
            {
                keep = getTimeName(clock);
            }

            remove_subdirs_except3(mesh->MESH_COMM, writeDir.c_str(), keep, "turbines", "precursor");
        }
    }

    return(0);
}

//***************************************************************************************************************//

void fatalErrorInFunction(const char* functionName, const char* errorMsg)
{
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    printf("\n\n[%d] --> Fatal error in function ", rank);
    printf("%s\n", functionName);
    printf("[%d]     %s\n\n", rank, errorMsg);

    // exit on current process, will force other processes to exit.
    // MPIFinalize is not called since would not force all other processors to
    // exit after the call to this function.
    exit(0);
}

//***************************************************************************************************************//

void warningInFunction(const char* functionName, const char* wrngMsg)
{
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    printf("\n[%d] --> Warning in function ", rank);
    printf("%s\n", functionName);
    printf("[%d]     %s\n\n", rank, wrngMsg);

    return;
}

//***************************************************************************************************************//

void createDir(MPI_Comm comm, const char* path)
{
    PetscMPIInt  rank; MPI_Comm_rank(comm, &rank);

    errno = 0;

    if(!rank)
    {
      PetscInt dirRes = mkdir(path, 0777);
      if(dirRes != 0 && errno != EEXIST)
      {
         char error[512];
          sprintf(error, "could not create %s directory\n", path);
          fatalErrorInFunction("createDir",  error);
      }

      // if directory already exist remove everything inside
      if(errno == EEXIST)
      {
          remove_subdirs(comm, path);
      }

    }
    return;
}

//***************************************************************************************************************//

word thisCaseName()
{
    word caseName;

    char buff[FILENAME_MAX];
    char *ptr = getcwd(buff, FILENAME_MAX);

    if(ptr != NULL)
    {
        word s(buff);

        char sep = '/';

        size_t i = s.rfind(sep, s.length());

        // if there is at least one '/' inside the path get the last
        // word as the case name
        if (i != std::string::npos)
        {
            caseName = s.substr(i+1);
        }
        // if there are no '/' get the whole word
        else
        {
            caseName = s;
        }

        return(caseName);
    }
    else
    {
        printf("\n\n --> Fatal error in function thisCaseName: could not retrieve case name\n");
        exit(0);
    }
}

//***************************************************************************************************************//

PetscInt foundInString(const char *str, word keyword)
{
    PetscInt flag = 0;

    word inputstring  = str;
    word sentence     = " " + inputstring + " ";

    // matching expression: " keyword "
    word fullKeyword1 = " " + keyword + " ";
    // matching expression: " keyword,"
    word fullKeyword2 = " " + keyword + ",";
    // matching expression: ",keyword"
    word fullKeyword3 = "," + keyword + " ";
    // matching expression: ",keyword "
    word fullKeyword4 = "," + keyword + ",";

    std::size_t found1 = sentence.find(fullKeyword1);
    std::size_t found2 = sentence.find(fullKeyword2);
    std::size_t found3 = sentence.find(fullKeyword3);
    std::size_t found4 = sentence.find(fullKeyword4);

    if(found1!=std::string::npos || found2!=std::string::npos || found3!=std::string::npos || found4!=std::string::npos)
    {
        flag = 1;
    }

    return(flag);
}

//***************************************************************************************************************//

PetscInt file_exist(const char *str)
{
    PetscInt r = 0;

    FILE *fp=fopen(str, "r");
    if(!fp)
    {
        r=0;
        printf("\n\n--> Warning: ");
        printf("file %s does not exist !!!\n\n", str);
    }
    else
    {
        fclose(fp);
        r = 1;
    }

    return r;
}

//***************************************************************************************************************//

PetscInt dir_exist(const char *str)
{
    DIR* dir = opendir(str);

    if(dir)
    {
        // drecory exists
        closedir(dir);
        return(1);
    }
    else if (ENOENT == errno)
    {
        // directory does not exist
        return(0);
    }
    else
    {
       char error[512];
        sprintf(error, "could not open %s directory\n", str);
        fatalErrorInFunction("dir_exists",  error);
    }

  return(0);
}

//***************************************************************************************************************//

PetscInt count_files(const char* path)
{
    // get file names
    DIR *dir; struct dirent *diread;

    // number of files
    PetscInt nFiles = 0;

    if((dir = opendir(path)) != nullptr)
    {
        while ((diread = readdir(dir)) != nullptr)
        {
            char* fileName = diread->d_name;
            if
            (
                strcmp(fileName, ".") !=0 &&
                strcmp(fileName, "..") !=0
            )
            {
                // increase rake count
                nFiles++;
            }
        }
        closedir(dir);
    }
    else
    {
       char error[512];
        sprintf(error, "could not open %s directory", path);
        fatalErrorInFunction("count_files",  error);
    }

    return(nFiles);
}

//***************************************************************************************************************//

void remove_dir(MPI_Comm comm, const char *path2dir)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank)
    {
        struct dirent *entry = NULL;
        DIR *dir = NULL;
        dir = opendir(path2dir);
        while(entry = readdir(dir))
        {
            DIR *sub_dir = NULL;
            FILE *file = NULL;
            char abs_path[456] = {0};
            if(*(entry->d_name) != '.')
            {
                sprintf(abs_path, "%s/%s", path2dir, entry->d_name);
                if(sub_dir = opendir(abs_path))
                {
                    closedir(sub_dir);
                    remove_dir(comm, abs_path);
                }
                else
                {
                    if(file = fopen(abs_path, "r"))
                    {
                        fclose(file);
                        remove(abs_path);
                    }
                }
            }
        }
        closedir(dir);
        remove(path2dir);
    }

    return;
}

//***************************************************************************************************************//

void remove_subdirs(MPI_Comm comm, const char *path2dir)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank)
    {
        struct dirent *entry = NULL;
        DIR *dir = NULL;
        dir = opendir(path2dir);
        while(entry = readdir(dir))
        {
            DIR *sub_dir = NULL;
            FILE *file = NULL;
            char abs_path[456] = {0};
            if(*(entry->d_name) != '.')
            {
                sprintf(abs_path, "%s/%s", path2dir, entry->d_name);
                if(sub_dir = opendir(abs_path))
                {
                    closedir(sub_dir);
                    remove_dir(comm, abs_path);
                }
                else
                {
                    if(file = fopen(abs_path, "r"))
                    {
                        fclose(file);
                        remove(abs_path);
                    }
                }
            }
        }
        closedir(dir);
    }

    return;
}

//***************************************************************************************************************//

void remove_subdirs_except(MPI_Comm comm, const char *path2dir, const word name)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank)
    {
        struct dirent *entry = NULL;
        DIR *dir = NULL;
        dir = opendir(path2dir);
        while(entry = readdir(dir))
        {
            DIR *sub_dir = NULL;
            FILE *file = NULL;
            char abs_path[456] = {0};
            if
            (
                *(entry->d_name) != '.' &&
                (entry->d_name) != name
            )
            {
                sprintf(abs_path, "%s/%s", path2dir, entry->d_name);
                if(sub_dir = opendir(abs_path))
                {
                    closedir(sub_dir);
                    remove_dir(comm, abs_path);
                }
                else
                {
                    if(file = fopen(abs_path, "r"))
                    {
                        fclose(file);
                        remove(abs_path);
                    }
                }
            }
        }
        closedir(dir);
    }

    return;
}

//***************************************************************************************************************//

void remove_subdirs_except2(MPI_Comm comm, const char *path2dir, const word name1, const word name2)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank)
    {
        struct dirent *entry = NULL;
        DIR *dir = NULL;
        dir = opendir(path2dir);
        while(entry = readdir(dir))
        {
            DIR *sub_dir = NULL;
            FILE *file = NULL;
            char abs_path[456] = {0};
            if
            (
                *(entry->d_name) != '.' &&
                (entry->d_name) != name1 &&
                (entry->d_name) != name2
            )
            {
                sprintf(abs_path, "%s/%s", path2dir, entry->d_name);
                if(sub_dir = opendir(abs_path))
                {
                    closedir(sub_dir);
                    remove_dir(comm, abs_path);
                }
                else
                {
                    if(file = fopen(abs_path, "r"))
                    {
                        fclose(file);
                        remove(abs_path);
                    }
                }
            }
        }
        closedir(dir);
    }

    return;
}

//***************************************************************************************************************//

void remove_subdirs_except3(MPI_Comm comm, const char *path2dir, const word name1, const word name2, const word name3)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    if(!rank)
    {
        struct dirent *entry = NULL;
        DIR *dir = NULL;
        dir = opendir(path2dir);
        while(entry = readdir(dir))
        {
            DIR *sub_dir = NULL;
            FILE *file = NULL;
            char abs_path[456] = {0};
            if
            (
                *(entry->d_name) != '.' &&
                (entry->d_name) != name1 &&
                (entry->d_name) != name2 &&
                (entry->d_name) != name3
            )
            {
                sprintf(abs_path, "%s/%s", path2dir, entry->d_name);
                if(sub_dir = opendir(abs_path))
                {
                    closedir(sub_dir);
                    remove_dir(comm, abs_path);
                }
                else
                {
                    if(file = fopen(abs_path, "r"))
                    {
                        fclose(file);
                        remove(abs_path);
                    }
                }
            }
        }
        closedir(dir);
    }

    return;
}

//***************************************************************************************************************//

word getTimeName(clock_ *clock)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << clock->time;
    return(stream.str());
}

//***************************************************************************************************************//

word getArbitraryTimeName(clock_ *clock, double timeValue)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << timeValue;
    return(stream.str());
}

//***************************************************************************************************************//

word getStartTimeName(clock_ *clock)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << clock->startTime;
    return(stream.str());
}


//***************************************************************************************************************//

PetscErrorCode setRunTimeWrite(domain_ *domain)
{
    PetscInt        nDomains = domain[0].info.nDomains;
    clock_          *clock   = domain[0].clock;
    word            timeName;                                             // time directory name
    word            writeDir;                                             // write directory
    PetscInt        dirRes;
    PetscMPIInt     rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);         // non-blocking

    // create/initialize fields directory (at simulation start only)
    if
    (
      clock->it == clock->itStart && !rank
    )
    {
      errno = 0;
      dirRes = mkdir("./fields", 0777);
      if(dirRes != 0 && errno != EEXIST)
      {
         char error[512];
          sprintf(error, "could not create fields directory \n");
          fatalErrorInFunction("setRunTimeWrite",  error);
      }
    }

    for(PetscInt d=0; d < nDomains; d++)
    {
        io_         *io    = domain[d].io;
        mesh_       *mesh  = domain[d].mesh;

        writeDir = "fields/" + mesh->meshName;

        // all fields are written in fields/time_/fieldName.
        // Writing settings are specified in control.dat as
        // intervalType : adjustableTime or timeStep
        // timeInterval : in seconds if adjustableTime, in iterations if timeStep

        // this domain communicator master rank
        MPI_Comm_rank(mesh->MESH_COMM, &rank);

        // create/initialize fields directory (at simulation start only)
        if
        (
            clock->it == clock->itStart && !rank
        )
        {
            errno = 0;
            dirRes = mkdir(writeDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", writeDir.c_str());
                fatalErrorInFunction("setRunTimeWrite",  error);
            }

            // if directory already exist remove everything inside except start time and turbines folder
            if(errno == EEXIST)
            {
                word startTimeName = getStartTimeName(clock);
                remove_subdirs_except2(mesh->MESH_COMM, writeDir.c_str(), startTimeName.c_str(), "turbines");
            }
        }

        word      intervalType = io->intervalType;
        PetscReal timeInterval = io->timeInterval;

        // write every "timeInterval" seconds
        if
        (
            (intervalType == "adjustableTime") &&
            (
                (clock->time - clock->startTime) / timeInterval -
                std::floor
                (
                    (clock->time - clock->startTime) / timeInterval
                ) < 1e-10
            )
        )
        {
            io->runTimeWrite = 1;
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
            io->runTimeWrite = 1;
        }
        // check for unkwnown type
        else if
        (
            (intervalType != "timeStep") &&
            (intervalType != "adjustableTime")
        )
        {
           char error[512];
            sprintf(error, "unknown interval type %s. Known types are timeStep and adjustableTime\n", intervalType.c_str());
            fatalErrorInFunction("setRunTimeWrite",  error);
        }
        // don't write at this time step
        else
        {
            io->runTimeWrite = 0;
        }

        // first processor creates the time folder
        if
        (
            !rank && io->runTimeWrite
        )
        {
            // create time directory for writing fields at current time
            timeName = writeDir + "/" + getTimeName(clock);

            PetscInt dirRes = mkdir(timeName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory\n", timeName.c_str());
                fatalErrorInFunction("setRunTimeWrite",  error);
            }

            // if directory already exist remove everything inside
            else if(errno == EEXIST)
            {
                remove_subdirs(mesh->MESH_COMM, timeName.c_str());
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictDouble(const char *dictName, const char *keyword, PetscReal *value)
{
    std::ifstream indata;

    char word[256];

    // pointer for stdtod and stdtol
    char *eptr;

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictDouble",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
               indata >> word;
               *value = std::strtod(word, &eptr);

               indata.close();
               return(0);
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictDouble",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictVector(const char *dictName, const char *keyword, Cmpnts *value)
{
    std::ifstream indata;

    char word[256];

    // pointer for stdtod and stdtol
    char *eptr;

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictVector",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
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

                   value->x = std::strtod(word, &eptr);

                   // check if the first component is a PetscReal, throw error otherwise
                   std::string cmp1(word);

                   if(!isNumber(cmp1))
                   {
                      char error[512];
                       sprintf(error, "expected number after keyword '%s' in dictionary %s\n", keyword, dictName);
                       fatalErrorInFunction("readDictVector",  error);
                   }

                   // read the second component
                   indata >> word;
                   value->y = std::strtod(word, &eptr);

                   // check if the second component is a PetscReal, throw error otherwise
                   std::string cmp2(word);

                   if(!isNumber(cmp2))
                   {
                      char error[512];
                       sprintf(error, "expected number in vector defined by keyword '%s' in dictionary %s\n", keyword, dictName);
                       fatalErrorInFunction("readDictVector",  error);
                   }

                   // read the third component (contains ")" character)
                   indata >> word;

                   std::string last(word);
                   if (last.find (")") != std::string::npos)
                   {
                       // remove ") character from the last component and store
                       value->z = std::strtod(word, &eptr);

                       // check if the first component is a PetscReal, throw error otherwise
                       std::string cmp3(word);

                       if(!isNumber(cmp3))
                       {
                          char error[512];
                           sprintf(error, "expected number in vector defined by keyword '%s' in dictionary %s\n", keyword, dictName);
                           fatalErrorInFunction("readDictVector",  error);
                       }

                       // close the file
                       indata.close();

                       // exit
                       return(0);

                   }
                   else
                   {
                      char error[512];
                       sprintf(error, "expected <)>  after vector defined by keyword '%s' in dictionary %s\n", keyword, dictName);
                       fatalErrorInFunction("readDictVector",  error);
                   }
               }
               else
               {
                char error[512];
                 sprintf(error, "expected <(>  after keyword '%s' in dictionary %s\n", keyword, dictName);
                 fatalErrorInFunction("readDictVector",  error);
               }
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictVector",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictInt(const char *dictName, const char *keyword, PetscInt *value)
{
    std::ifstream indata;

    char word[256];

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictInt",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
               indata >> word;
               std::sscanf(word, "%ld", value);

               indata.close();
               return(0);
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictInt",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictWord(const char *dictName, const char *keyword, word *value)
{
    std::ifstream indata;

    char word[256];

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictWord",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
               indata >> word;

               // allocate memory
               PetscMalloc(sizeof(std::string), value);
               new(value) std::string{};

               // set the value
               *value = word;

               indata.close();
               return(0);
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictWord",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictWordAndDouble(const char *dictName, const char *keyword, word *value1, PetscReal *value2)
{
    std::ifstream indata;

    char word[256];

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictWordAndDouble",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
               // read the string
               indata >> word;
               *value1 = word;

               // if the word contains 'fixed' also read the following value
               std::string key(word);
               if (key.find ("fixed") != std::string::npos)
               {
                   // read the value
                   indata >> word;
                   std::sscanf(word, "%lf", value2);

                   // check if the value is a PetscReal, throw error otherwise
                   std::string str(word);

                   if(isNumber(str))
                   {
                       if (str.find ('.') != std::string::npos)
                       {
                           indata.close();
                           return(0);
                       }
                       else
                       {
                          char error[512];
                           sprintf(error, "expected <PetscReal>  after keyword '%s' in dict %s\n", keyword, dictName);
                           fatalErrorInFunction("readDictWordAndDouble",  error);
                       }
                   }
                   else
                   {
                      char error[512];
                       sprintf(error, "expected <PetscReal> after keyword '%s' in %s dictionary\n", keyword, dictName);
                       fatalErrorInFunction("readDictWordAndDouble",  error);
                   }
               }
               // if the string does not contain fixed, value is set to 0
               else
               {
                   *value2 = 0.0;
                   return(0);
               }
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictWordAndDouble",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readDictWordAndVector(const char *dictName, const char *keyword, word *value1, Cmpnts *value2)
{
    std::ifstream indata;

    char word[256];

    // pointer for stdtod and stdtol
    char *eptr;

    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readDictWordAndDouble",  error);
    }
    else
    {
        while(!indata.eof())
        {
            indata >> word;

            if
            (
                strcmp
                (
                    keyword,
                    word
                ) == 0
            )
            {
               // read the string
               indata >> word;
               *value1 = word;

               // if the word contains 'fixed' also read the following value
               std::string key(word);
               if (key.find ("fixed") != std::string::npos)
               {
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
                          word[i] = word[i+1];
                      }

                      // save the first component
                      value2->x = std::strtod(word, &eptr);

                      // check if the first component is a PetscReal, throw error otherwise
                      std::string str1(word);

                      if(isNumber(str1))
                      {

                          if (str1.find ('.') == std::string::npos)
                          {
                             char error[512];
                              sprintf(error, "expected <PetscReal>  after keyword '%s' in dict %s\n", keyword, dictName);
                              fatalErrorInFunction("readDictWordAndVector",  error);
                          }
                      }
                      else
                      {
                         char error[512];
                          sprintf(error, "expected <PetscReal> after keyword '%s' in %s dictionary\n", keyword, dictName);
                          fatalErrorInFunction("readDictWordAndVector",  error);
                      }

                      // read the second component
                      indata >> word;
                      value2->y = std::strtod(word, &eptr);

                      // check if the second component is a PetscReal, throw error otherwise
                      std::string str2(word);

                      if(isNumber(str2))
                      {
                          if (str2.find ('.') == std::string::npos)
                          {
                             char error[512];
                              sprintf(error, "expected <PetscReal>  in vector defined by keyword '%s' in dict %s\n", keyword, dictName);
                              fatalErrorInFunction("readDictWordAndVector",  error);
                          }
                      }
                      else
                      {
                         char error[512];
                          sprintf(error, "expected <PetscReal> in vector defined by keyword '%s' in %s dictionary\n", keyword, dictName);
                          fatalErrorInFunction("readDictWordAndVector",  error);
                      }

                      // read the third component (contains ")" character)
                      indata >> word;

                      std::string last(word);
                      if (last.find (")") != std::string::npos)
                      {
                          // remove ") character from the last component and store
                          value2->z = std::strtod(word, &eptr);

                          // check if the first component is a PetscReal, throw error otherwise
                          std::string str3(word);

                          if(isNumber(str3))
                          {
                              if (str3.find ('.') != std::string::npos)
                              {
                                  indata.close();
                                  return(0);
                              }
                              else
                              {
                                 char error[512];
                                  sprintf(error, "expected <PetscReal>  in vector defined by keyword '%s' in dict %s\n", keyword, dictName);
                                  fatalErrorInFunction("readDictWordAndVector",  error);
                              }
                          }
                          else
                          {
                             char error[512];
                              sprintf(error, "expected <PetscReal> in vector defined by keyword '%s' in %s dictionary\n", keyword, dictName);
                              fatalErrorInFunction("readDictWordAndVector",  error);
                          }

                      }
                      else
                      {
                         char error[512];
                          sprintf(error, "expected <)>  after vector defined by keyword '%s' in dict %s\n", keyword, dictName);
                          fatalErrorInFunction("readDictWordAndVector",  error);
                      }

                   }
                   else
                   {
                      char error[512];
                       sprintf(error, "expected <(>  after keyword '%s' in dict %s\n", keyword, dictName);
                       fatalErrorInFunction("readDictWordAndVector",  error);
                   }
               }
               // if the string does not contain fixed, value is set to 0
               else
               {
                   value2->x = 0.0;
                   value2->y = 0.0;
                   value2->z = 0.0;
                   return(0);
               }
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find keyword %s in dictionary %s\n", keyword, dictName);
        fatalErrorInFunction("readDictWordAndDouble",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readSubDictDouble(const char *dictName, const char *subdict, const char *keyword, PetscReal *value)
{
    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // pointer for strtod
    char *eptr;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readSubDictDouble",  error);
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
                    subdict,
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
                   // start reading inside the subdictionary
                   indata >> token2;

                   // read until "}" or end of dictionary is encountered
                   while(trim(token2)!="}" && !indata.eof())
                   {
                       // look for the keyword
                       if(token2==keyword)
                       {
                           indata >> token2;

                           // check if the value is a PetscReal or PetscInt
                           // (not really useful since can be casted on calling the function)

                           // chech if value contains character (throws errors if yes)
                           if(isNumber(token2))
                           {
                               *value = std::strtod(token2.c_str(), &eptr);

                               // now look for the terminating "}", trows error if not found
                               // if find another "{" means another subdict is entered: throws error

                               indata >> token2;

                               while(!indata.eof() && trim(token2)!="{" )
                               {
                                   if(trim(token2)=="}")
                                   {
                                       indata.close();
                                       return(0);
                                   }

                                   indata >> token2;
                               }

                              char error[512];
                               sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                               fatalErrorInFunction("readSubDictDouble",  error);
                           }
                           else
                           {
                              char error[512];
                               sprintf(error, "expected number after keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                               fatalErrorInFunction("readSubDictDouble",  error);
                           }

                       }

                       indata >> token2;
                   }

                  char error[512];
                   sprintf(error, "could not find keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                   fatalErrorInFunction("readSubDictDouble",  error);
               }
               else
               {
                   char error[512];
                   sprintf(error, "expected '{' token in subdict %s of %s dictionary, found '%s'\n", subdict, dictName, word);
                   fatalErrorInFunction("readSubDictDouble",  error);
               }
           }
       }

      char error[512];
       sprintf(error, "could not find subdictionary %s in dictionary %s\n", subdict, dictName);
       fatalErrorInFunction("readSubDictDouble",  error);
   }
   indata.close();

   return(0);
}

//***************************************************************************************************************//

PetscErrorCode readSubDictInt(const char *dictName, const char *subdict, const char *keyword, PetscInt *value)
{
    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // pointer for strtod
    char *eptr;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readSubDictInt",  error);
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
                    subdict,
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
                   // start reading inside the subdictionary
                   indata >> token2;

                   // read until "}" or end of dictionary is encountered
                   while(trim(token2)!="}" && !indata.eof())
                   {
                       // look for the keyword
                       if(token2==keyword)
                       {
                           indata >> token2;

                           // check if the value is a PetscReal or PetscInt
                           // (not really useful since can be casted on calling the function)

                           // chech if value contains character (throws errors if yes)
                           if(isNumber(token2))
                           {
                               *value = std::strtol(token2.c_str(), &eptr, 10);

                               // now look for the terminating "}", trows error if not found
                               // if find another "{" means another subdict is entered: throws error

                               indata >> token2;

                               while(!indata.eof() && trim(token2)!="{" )
                               {
                                   if(trim(token2)=="}")
                                   {
                                       indata.close();
                                       return(0);
                                   }

                                   indata >> token2;
                               }

                              char error[512];
                               sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                               fatalErrorInFunction("readSubDictInt",  error);
                           }
                           else
                           {
                              char error[512];
                               sprintf(error, "expected number after keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                               fatalErrorInFunction("readSubDictInt",  error);
                           }

                       }

                       indata >> token2;
                   }

                  char error[512];
                   sprintf(error, "could not find keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                   fatalErrorInFunction("readSubDictInt",  error);
               }
               else
               {
                  char error[512];
                   sprintf(error, "expected '{' token in subdict %s of %s dictionary, found '%s'\n", subdict, dictName, word);
                   fatalErrorInFunction("readSubDictInt",  error);
               }
           }
       }

      char error[512];
       sprintf(error, "could not find subdictionary %s in dictionary %s\n", subdict, dictName);
       fatalErrorInFunction("readSubDictInt",  error);
   }
   indata.close();

   return(0);
}

//***************************************************************************************************************//

PetscErrorCode readSubDictWord(const char *dictName, const char *subdict, const char *keyword, word *value)
{
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
        fatalErrorInFunction("readSubDictWord",  error);
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
                    subdict,
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
                   // start reading inside the subdictionary
                   indata >> token2;

                   // read until "}" or end of dictionary is encountered
                   while(trim(token2)!="}" && !indata.eof())
                   {
                       // look for the keyword
                       if(token2==keyword)
                       {
                           indata >> word;

                           // check if the value is a word
                           std::string str(word);

                           if(isNumber(str))
                           {
                              char error[512];
                               sprintf(error, "expected word after keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                               fatalErrorInFunction("readSubDictWord",  error);
                           }
                           else
                           {
                               // allocate memory
                               PetscMalloc(sizeof(std::string), value);
                               new(value) std::string{};

                               // store the value
                               *value = word;

                               // now look for the terminating "}", trows error if not found
                               // if find another "{" means another subdict is entered: throws error

                               indata >> token2;

                               while(!indata.eof() && trim(token2)!="{")
                               {
                                   if(trim(token2)=="}")
                                   {
                                       indata.close();
                                       return(0);
                                   }

                                   indata >> token2;
                               }

                              char error[512];
                               sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                               fatalErrorInFunction("readSubDictWord",  error);
                           }

                       }

                       indata >> token2;
                   }

                  char error[512];
                   sprintf(error, "could not find keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                   fatalErrorInFunction("readSubDictWord",  error);
               }
               else
               {
                  char error[512];
                   sprintf(error, "expected '{' token in subdict %s of %s dictionary, found '%s'\n", subdict, dictName, word);
                   fatalErrorInFunction("readSubDictWord",  error);
               }
           }
       }

      char error[512];
       sprintf(error, "could not find subdictionary %s in dictionary %s\n", subdict, dictName);
       fatalErrorInFunction("readSubDicWord",  error);
   }
   indata.close();

   return(0);
}

//***************************************************************************************************************//

PetscErrorCode readSubDictVector(const char *dictName, const char *subdict, const char *keyword, Cmpnts *value)
{
    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // pointer for strtod and strtol
    char *eptr;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readSubDictVector",  error);
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
                    subdict,
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
                   // start reading inside the subdictionary
                   indata >> token2;

                   // read until "}" or end of dictionary is encountered
                   while(trim(token2)!="}" && !indata.eof())
                   {
                       // look for the keyword
                       if(token2==keyword)
                       {
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

                              value->x = std::strtod(word, &eptr);

                              // check if the first component is a PetscReal, throw error otherwise
                              std::string cmp1(word);

                              if(isNumber(cmp1))
                              {

                                  if (cmp1.find ('.') == std::string::npos)
                                  {
                                     char error[512];
                                      sprintf(error, "expected <PetscReal>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                      fatalErrorInFunction("readSubDictVector",  error);
                                  }
                              }
                              else
                              {
                                 char error[512];
                                  sprintf(error, "expected <PetscReal> after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                  fatalErrorInFunction("readSubDictVector",  error);
                              }

                              // read the second component
                              indata >> word;
                              value->y = std::strtod(word, &eptr);

                              // check if the second component is a PetscReal, throw error otherwise
                              std::string cmp2(word);

                              if(isNumber(cmp2))
                              {
                                  if (cmp2.find ('.') == std::string::npos)
                                  {
                                     char error[512];
                                      sprintf(error, "expected <PetscReal>  in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                      fatalErrorInFunction("readSubDictVector",  error);
                                  }
                              }
                              else
                              {
                                 char error[512];
                                  sprintf(error, "expected <PetscReal> in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                  fatalErrorInFunction("readSubDictVector",  error);
                              }

                              // read the third component (contains ")" character)
                              indata >> word;

                              std::string last(word);
                              if (last.find (")") != std::string::npos)
                              {
                                  // remove ") character from the last component and store
                                  value->z = std::strtod(word, &eptr);

                                  // check if the first component is a PetscReal, throw error otherwise
                                  std::string cmp3(word);

                                  if(isNumber(cmp3))
                                  {
                                      if (cmp3.find ('.') == std::string::npos)
                                      {
                                         char error[512];
                                          sprintf(error, "expected <PetscReal>  in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                          fatalErrorInFunction("readSubDictVector",  error);
                                      }
                                  }
                                  else
                                  {
                                     char error[512];
                                      sprintf(error, "expected <PetscReal> in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                      fatalErrorInFunction("readSubDictVector",  error);
                                  }

                                  // now look for the terminating "}", trows error if not found
                                  // if find another "{" means another subdict is entered: throws error

                                  indata >> token2;

                                  while(!indata.eof() && trim(token2)!="{")
                                  {
                                      if(trim(token2)=="}")
                                      {
                                          indata.close();
                                          return(0);
                                      }

                                      indata >> token2;
                                  }

                                 char error[512];
                                  sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                                  fatalErrorInFunction("readSubDictWord",  error);
                              }
                              else
                              {
                                 char error[512];
                                  sprintf(error, "expected <)>  after vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                  fatalErrorInFunction("readSubDictVector",  error);
                              }
                          }
                          else
                          {
                               char error[512];
                                sprintf(error, "expected <(>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                fatalErrorInFunction("readSubDictVector",  error);
                           }
                       }

                       indata >> token2;
                   }

                  char error[512];
                   sprintf(error, "could not find keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                   fatalErrorInFunction("readSubDictVector",  error);
               }
               else
               {
                  char error[512];
                   sprintf(error, "expected '{' token in subdict %s of %s dictionary, found '%s'\n", subdict, dictName, word);
                   fatalErrorInFunction("readSubDictVector",  error);
               }
           }
       }
       indata.close();

      char error[512];
       sprintf(error, "could not find subdictionary %s in dictionary %s\n", subdict, dictName);
       fatalErrorInFunction("readSubDictVector",  error);
   }

   return(0);
}

//***************************************************************************************************************//

PetscErrorCode readSubDictIntArray(const char *dictName, const char *subdict, const char *keyword, labelList &value)
{
    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // pointer for strtod and strtol
    char *eptr;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readSubDictIntArray",  error);
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
                            subdict,
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
                    // start reading inside the subdictionary
                    indata >> token2;

                    // read until "}" or end of dictionary is encountered
                    while(trim(token2)!="}" && !indata.eof())
                    {
                        // look for the keyword
                        if(token2==keyword)
                        {
                            // read the first component (contains "(" character)
                            indata >> word;

                            std::string first(word);
                            if (first.find ("(") != std::string::npos)
                            {
                                // remove "("" character from the first component
                                PetscInt l1 = first.size();
                                PetscInt sum = 0;
                                for(PetscInt i=0;i<l1;i++)
                                {
                                    sum += isdigit(first[i]);
                                    // save the first component
                                    word[i] = word[i+1];

                                }

                                if(sum == 0)
                                {
                                   char error[512];
                                    sprintf(error, "empty array or expected <PetscInt>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                    fatalErrorInFunction("readSubDictIntArray",  error);
                                }

                                PetscInt *numPtr = new PetscInt;
                                *numPtr = std::strtod(word, &eptr);

                                std::string cmp1(word);

                                if(isNumber(cmp1))
                                {

                                    if (cmp1.find ('.') != std::string::npos)
                                    {
                                       char error[512];
                                        sprintf(error, "expected <PetscInt>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                        fatalErrorInFunction("readSubDictIntArray",  error);
                                    }
                                }
                                else
                                {
                                   char error[512];
                                    sprintf(error, "expected <PetscInt> after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                    fatalErrorInFunction("readSubDictIntArray",  error);
                                }

                                value.push_back(*numPtr);

                                // already at the end of the array
                                if(first.find (")") != std::string::npos){
                                    indata >> token2;

                                    while(!indata.eof() && trim(token2)!="{")
                                    {
                                        if(trim(token2)=="}")
                                        {
                                            indata.close();
                                            return(0);
                                        }

                                        indata >> token2;
                                    }

                                   char error[512];
                                    sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                                    fatalErrorInFunction("readSubDictWord",  error);
                                }
                                else
                                {
                                    indata >> word;

                                    std::string next(word);

                                    while (next.find (")") == std::string::npos)
                                    {
                                        std::string cmp1(word);

                                        PetscInt l1 = cmp1.size();
                                        for(PetscInt i=0;i<l1;i++)
                                        {
                                            if (!isdigit(cmp1[i]))
                                            {
                                               char error[512];
                                                sprintf(error, "expected <PetscInt> or expected <)>  after vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                                fatalErrorInFunction("readSubDictIntArray",  error);
                                            }
                                        }

                                        if(isNumber(cmp1))
                                        {

                                            if (cmp1.find ('.') != std::string::npos)
                                            {
                                               char error[512];
                                                sprintf(error, "expected <PetscInt>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                                fatalErrorInFunction("readSubDictIntArray",  error);
                                            }
                                        }
                                        else
                                        {
                                           char error[512];
                                            sprintf(error, "expected <PetscInt> after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                            fatalErrorInFunction("readSubDictIntArray",  error);
                                        }

                                        PetscInt num = std::strtod(word, &eptr);

                                        value.push_back(num);

                                        indata >> word;

                                        next = word;

                                    }

                                    if (next.find (")") != std::string::npos)
                                    {
                                        std::string cmp3(word);

                                        if(isNumber(cmp3))
                                        {
                                            if (cmp3.find ('.') != std::string::npos)
                                            {
                                               char error[512];
                                                sprintf(error, "expected <PetscInt>  in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                                fatalErrorInFunction("readSubDictIntArray",  error);
                                            }
                                        }
                                        else
                                        {
                                           char error[512];
                                            sprintf(error, "expected <PetscInt> in vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                            fatalErrorInFunction("readSubDictIntArray",  error);
                                        }

                                        PetscInt num = std::strtod(word, &eptr);

                                        value.push_back(num);

                                        indata >> token2;

                                        while(!indata.eof() && trim(token2)!="{")
                                        {
                                            if(trim(token2)=="}")
                                            {
                                                indata.close();
                                                return(0);
                                            }

                                            indata >> token2;
                                        }

                                       char error[512];
                                        sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                                        fatalErrorInFunction("readSubDictIntArray",  error);
                                    }
                                    else
                                    {
                                       char error[512];
                                        sprintf(error, "expected <)>  after vector defined by keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                        fatalErrorInFunction("readSubDictIntArray",  error);
                                    }
                                }

                               char error[512];
                                sprintf(error, "missing '}' token in subdict %s of %s dictionary\n", subdict, dictName);
                                fatalErrorInFunction("readSubDictIntArray",  error);
                            }
                            else
                            {
                                 char error[512];
                                  sprintf(error, "expected <(>  after keyword '%s' in subdict %s of %s dictionary\n", keyword, subdict, dictName);
                                  fatalErrorInFunction("readSubDictIntArray",  error);
                             }
                        }

                        indata >> token2;
                    }

                   char error[512];
                    sprintf(error, "could not find keyword '%s' in subdictionary %s of %s dictionary\n", keyword, subdict, dictName);
                    fatalErrorInFunction("readSubDictIntArray",  error);

                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token in subdict %s of %s dictionary, found '%s'\n", subdict, dictName, word);
                    fatalErrorInFunction("readSubDictIntArray",  error);
                }
            }
        }
        indata.close();

       char error[512];
        sprintf(error, "could not find subdictionary %s in dictionary %s\n", subdict, dictName);
        fatalErrorInFunction("readSubDictIntArray",  error);
    }

    return(0);
}

//***************************************************************************************************************//

word trim(const word& str)
{
    auto first = str.find_first_not_of(' ');
    if(first == std::string::npos) first = 0;

    auto last = str.find_last_not_of(' ');
    if(last == std::string::npos) last = str.size();

    return str.substr(first, last-first+1);
}

//***************************************************************************************************************//

bool isNumber(const word& str)
{
    PetscInt i;
    PetscInt nKeys = 48;
    char keyboard[nKeys];
    sprintf(keyboard, "qwertyuiopasdfghjklzxcvbnmè+ùàòé*§°ç,:;");

    for(i=0;i<nKeys; i++)

    if (str.find(keyboard[i]) != std::string::npos)
    {
        return false;
    }

    return true;
}
