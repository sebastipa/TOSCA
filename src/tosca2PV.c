
//#include </usr/include/hdf5/serial/hdf5.h>
#include "include/base.h"
#include "include/domain.h"
#include "include/initialization.h"
#include "include/inline.h"
#include <hdf5.h>
#include "include/tosca2PV.h"

static char head[] = "TOSCA Post Processor";

//***************************************************************************************************************//

int main(int argc, char **argv)
{

  // initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, head);

  // domains array
  domain_ *domain;

  // postProcess
  postProcess pp;

  // simulation clock
  clock_  clock;
  ReadTimeControls(&clock);

  // simulation flags
  flags_ flags;
  SetSimulationFlags(&flags);

  // simulation info
  simInfo_ info;
  SetSimulationInfo(&info);

  // initialize post-processing clock
  PetscReal ppTimeStart,ppTimeEnd; PetscTime(&ppTimeStart);

  PetscMPIInt rank;   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscMPIInt nProcs; MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);



  // initialize the flags
  pp.postProcessFields    = 0;
  pp.writeRaster          = 0;
  pp.samplingSections     = 0;
  pp.postProcessPrecursor = 0;

  // read from control file
  PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-postProcessFields", (PetscInt*)&(pp.postProcessFields), PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-writeRaster", (PetscInt*)&(pp.writeRaster), PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sections", (PetscInt*)&(pp.samplingSections), PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-postProcessPrecursor", (PetscInt*)&(pp.postProcessPrecursor), PETSC_NULL);

  // create XMF folder
  if (!rank) createDir(PETSC_COMM_WORLD, "./XMF");

  // initialize post processing
  if(pp.postProcessFields || pp.writeRaster || pp.samplingSections)
  {
      postProcessInitialize(&domain, &clock, &info, &flags);

      // turn-off field post processing if parallel
      if(nProcs > 1 && pp.postProcessFields == 1)
      {
          char warning[256];
          sprintf(warning, "-postProcessFields is not available in parallel: setting to 0");
          warningInFunction("main",  warning);
          pp.postProcessFields = 0;
      }

      // do probes acquisition
      if(flags.isAquisitionActive)
      {
          postProcessWriteProbes(domain);
      }
  }

  // write 3D fields into XMF
  if(pp.postProcessFields || pp.writeRaster)
  {
      binary3DToXMF(domain, &pp);
  }

  // initialize precursor post processing
  if(pp.postProcessPrecursor)
  {
      postProcessInitializePrecursor(&pp, &clock);

      binary3DToXMF(pp.precursor->domain, &pp);

      // read sections
      sectionsReadAndAllocate(pp.precursor->domain);

      // write k-sections into XMF
      binaryKSectionsToXMF(pp.precursor->domain);

      // write j-sections into XMF
      binaryJSectionsToXMF(pp.precursor->domain, &pp);

      // write i-sections into XMF
      binaryISectionsToXMF(pp.precursor->domain);

      // write perturbation data
      binaryKSectionsPerturbToXMF(pp.precursor->domain);
      binaryJSectionsPerturbToXMF(pp.precursor->domain, &pp);
      binaryISectionsPerturbToXMF(pp.precursor->domain);
  }

  if(pp.samplingSections)
  {
      // read sections
      sectionsReadAndAllocate(domain);

      // write k-sections into XMF
      binaryKSectionsToXMF(domain);

      // write j-sections into XMF
      binaryJSectionsToXMF(domain, &pp);

      // write i-sections into XMF
      binaryISectionsToXMF(domain);

      // write perturbation data
      binaryKSectionsPerturbToXMF(domain);
      binaryJSectionsPerturbToXMF(domain, &pp);
      binaryISectionsPerturbToXMF(domain);

      // on-the-fly read and cut sections from average fields
      fieldKSectionsToXMF(domain);
      fieldJSectionsToXMF(domain);
      fieldISectionsToXMF(domain);

      // on-the-fly user defined surface from average fields
      fieldUserDefinedPlaneToXMF(domain);

  }

    PetscTime(&ppTimeEnd);

    PetscPrintf(PETSC_COMM_WORLD, "Cpu Time = %lf s, Finalizing run\n", ppTimeEnd - ppTimeStart);

    PetscBarrier(PETSC_NULL);

    PetscPrintf(PETSC_COMM_WORLD, "\nEnd\n\n");

    PetscFinalize();
}

//***************************************************************************************************************//

PetscErrorCode postProcessInitialize(domain_ **domainAddr, clock_ *clock, simInfo_ *info, flags_ *flags)
{
  // print logo
  PrintOkWindLogo();

  // reads parent/child tree and allocates os memory
  SetDomainsAndAllocate(domainAddr, flags, info);

  domain_ *domain = *domainAddr;

  // set simulation start time
  SetStartTime(clock, domain,info);

  for(PetscInt d=0; d<info->nDomains; d++)
  {
      PetscPrintf(PETSC_COMM_WORLD, "\nDomain %ld\n", d);
      PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");

      // set pointer to time controls
      domain[d].clock = clock;

      // read physical constants
      ReadPhysicalConstants(&domain[d]);

      // set access database pointers
      SetAccessPointers(&domain[d]);

      // set boundary conditions
      SetBoundaryConditions(domain[d].mesh);

      // initialize mesh
      InitializeMesh(domain[d].mesh);

      // initialize i/o controls and initialization type
      InitializeIO(domain[d].io);

      // momentum equation initialize
      InitializeUEqn(domain[d].ueqn);

      // poisson equation initialize
      InitializePEqn(domain[d].peqn);

      // temperature equation initialize
      InitializeTEqn(domain[d].teqn);

      // LES model initialize
      InitializeLES(domain[d].les);

      // averages and kebudgets initialize
      if(flags->isAquisitionActive)
      {
          acquisition_ *acquisition = domain[d].acquisition;

          // read acquisition flags
          acquisition->isProbesActive     = 0;
          acquisition->isSectionsActive   = 0;
          acquisition->isAverageABLActive = 0;
          acquisition->isAverage3LMActive = 0;
          acquisition->isPerturbABLActive = 0;

          PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sections",      &(acquisition->isSectionsActive),   PETSC_NULL);
          PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-average3LM",    &(acquisition->isAverage3LMActive), PETSC_NULL);
          PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-perturbABL",    &(acquisition->isPerturbABLActive), PETSC_NULL);
          PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-probes",        &(acquisition->isProbesActive),     PETSC_NULL);

          averageFieldsInitialize(domain[d].acquisition);
          averageKEBudgetsInitialize(domain[d].acquisition);
          perturbationABLInitialize(domain[d].acquisition);

          // temporary (eliminate for 3LM fields)
          acquisition->isAverage3LMActive = 0;
      }

      PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");
  }

  // initialize probes
  if
  (
      domain[0].access.io->averaging ||
      domain[0].access.io->phaseAveraging
  )
  {
      // initialize probes with post processing flag activated
      ProbesInitialize(domain,1);
  }

  //averaging3LMInitialize(domain);

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode postProcessInitializePrecursor(postProcess *pp, clock_ *clock)
{
    PetscPrintf(PETSC_COMM_WORLD, "\n Precursor \n");
    PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");

    PetscMalloc(sizeof(precursor_), &(pp->precursor));
    precursor_ *precursor = pp->precursor;

    PetscMalloc(sizeof(domain_), &(precursor->domain));
    domain_ *domain = precursor->domain;

    // set solution flags
    SetSolutionFlagsPrecursor(domain);

    // set simulation info
    SetSimulationInfo(&(domain->info));

	// read physical constants
    ReadPhysicalConstants(domain);

    // allocate domain memory
    SetDomainMemory(domain);

    // set pointer to time controls
    domain->clock = clock;

    // set access database pointers
    SetAccessPointers(domain);

    // set mesh name for boundary conditions
    domain->mesh->meshName = ".";

    // set boundary conditions
    SetBoundaryConditions(domain->mesh);

    // set mesh name for fields
    domain->mesh->meshName = "precursor";

    // initialize mesh
    InitializeMesh(domain->mesh);

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

    // initialize equations
    InitializeUEqn(domain->ueqn);
    InitializePEqn(domain->peqn);
    InitializeTEqn(domain->teqn);
    InitializeLES(domain->les);

    // initialize acquisition
    InitializeAcquisitionPrecursor(domain);

    return(0);
}

PetscErrorCode binary3DToXMF(domain_ *domain, postProcess *pp)
{
    PetscMPIInt  rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscInt     nDomains = domain[0].info.nDomains;
    PetscInt     readError;
    word         fieldsPath, XMFDir, fieldsFileName;
    PetscInt     ntimes;


    for (PetscInt d=0; d < nDomains; d++)
    {
      fieldsPath = "./fields/" + domain[d].mesh->meshName;
      XMFDir     = "./XMF/" + domain[d].mesh->meshName;

      // create domain directory within XMF folder
      errno = 0;
      PetscInt dirRes = mkdir(XMFDir.c_str(), 0777);
      if(dirRes != 0 && errno != EEXIST)
      {
          char error[512];
          sprintf(error, "could not create directory %s for XMF\n", XMFDir.c_str());
          fatalErrorInFunction("binary3DToXMF",  error);
      }

      // get list of times from the fieldsPath of the current domain
      std::vector<PetscReal> timeSeries;
      getTimeList(fieldsPath.c_str(), timeSeries, ntimes);

      // create XMF file
      fieldsFileName = XMFDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + ".xmf";
      FILE *xmf;

      if(!rank)
      {
          // write XMF file intro
          xmf = fopen(fieldsFileName.c_str(), "w");
          fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
          fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
          fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
          fprintf(xmf, " <Domain>\n");
          fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
          fclose(xmf);
      }

      // loop over times
      for(PetscInt ti=0; ti<ntimes; ti++)
      {
          PetscPrintf(PETSC_COMM_WORLD, "\n\nTime: %lf\n\n", timeSeries[ti]);

          // read fields inside time folder
          readFields(&domain[d], timeSeries[ti]);

          // write fields inside XMF folder
          writeFieldsToXMF(&domain[d], fieldsFileName.c_str(), timeSeries[ti]);
      }

      if(!rank)
      {
        // write XMF file end
        xmf = fopen(fieldsFileName.c_str(), "a");
        fprintf(xmf, "   </Grid>\n");
        fprintf(xmf, " </Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
      }
    }

    return(0);
}

//***************************************************************************************************************//

// write fields inside XMF folder
PetscErrorCode writeFieldsToXMF(domain_ *domain, const char* filexmf, PetscReal time)
{
    mesh_  *mesh  = domain->mesh;
    clock_ *clock = domain->clock;
    io_    *io    = domain->io;
    flags_ *flags = &(domain->flags);
    acquisition_ *acquisition = domain->acquisition;

    // HDF5 file with path
    word fileName;

      // HDF5 file w/o path
  	word hdfileName;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;


    fileName   = "./XMF/" + mesh->meshName + "/" + thisCaseName() + "_" + stream.str();
    hdfileName = thisCaseName() + "_" + stream.str();

    DMDALocalInfo info = mesh->info;

    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    FILE          *xmf;

    // open this time section in the XMF file
    xmfWriteFileStartTimeSection(xmf, filexmf,mx-1, my-1, mz-1,"3DSMesh",time);

    // Write the data file.
    hid_t     dataspace_id;
    hsize_t   dims[3];
    herr_t    status;

    // write HDF file
    hid_t	file_id;
    file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // ************************** write mesh points **************************

    dims[0]	= mz - 1;
    dims[1]	= my - 1;
    dims[2]	= mx - 1;

    dataspace_id = H5Screate_simple(3, dims, NULL);

    writePointsToXMF
    (
        mesh,
        filexmf,
        hdfileName.c_str(),
        &file_id,
        &dataspace_id,
        time
    );

    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    // ************************** write cell fields **************************

    dims[0]	= mz - 2;
    dims[1]	= my - 2;
    dims[2]	= mx - 2;

    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    writeVectorToXMF
    (
        domain,
        filexmf,
        hdfileName.c_str(),
        &file_id,
        &dataspace_id,
        time,
        "U",
        domain->ueqn->Ucat
    );

    writeScalarToXMF
    (
        domain,
        filexmf,
        hdfileName.c_str(),
        &file_id,
        &dataspace_id,
        time,
        "p",
        domain->peqn->P
    );

    writeScalarToXMF
    (
        domain,
        filexmf,
        hdfileName.c_str(),
        &file_id,
        &dataspace_id,
        time,
        "nv",
        mesh->Nvert
    );

    writeScalarToXMF
    (
        domain,
        filexmf,
        hdfileName.c_str(),
        &file_id,
        &dataspace_id,
        time,
        "meshTag",
        mesh->meshTag
    );

    if(domain->flags.isLesActive)
    {
        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "nut",
            domain->les->lNu_t
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "cs",
            domain->les->lCs
        );

        if(domain->flags.isLesActive == 2)
        {
            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "sgsL",
                domain->les->L
            );
        }
    }

    if(domain->flags.isTeqnActive)
    {
        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "T",
            domain->teqn->Tmprt
        );
    }

    if(io->qCrit)
    {
        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "Q",
            acquisition->fields->Q
        );
    }

    if(io->windFarmForce)
    {
        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "bf",
            acquisition->fields->windFarmForce
        );
    }

    if(io->sources)
    {
        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "Coriolis",
            acquisition->fields->Coriolis
        );

        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "Driving",
            acquisition->fields->Driving
        );

        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "Damping",
            acquisition->fields->xDamping
        );

        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "CanopyForce",
            acquisition->fields->CanopyForce
        );

        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "yDampU",
            acquisition->fields->yDampMappedU
        );
    }

    if(io->continuity)
    {
        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "divU",
            acquisition->fields->divU
        );
    }

    if(io->buoyancy)
    {
        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "buoyancy",
            domain->ueqn->bTheta
        );
    }

    if(io->averaging)
    {
        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "avgU",
            acquisition->fields->avgU
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "avgP",
            acquisition->fields->avgP
        );

        writeSymmTensorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "avgUU",
            acquisition->fields->avgUU
        );

        if(domain->flags.isLesActive)
        {
            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "avgNut",
                acquisition->fields->avgNut
            );

            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "avgCs",
                acquisition->fields->avgCs
            );
        }

        if(io->averaging > 1)
        {
            writeVectorToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "avgOmega",
                acquisition->fields->avgOmega
            );

            if(io->averaging > 2)
            {
                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "avgPsq",
                    acquisition->fields->avgP2
                );

                writeSymmTensorToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "avgOmegaOmega",
                    acquisition->fields->avgOmegaOmega
                );

                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "avgUdotGradP",
                    acquisition->fields->avgUdotGradP
                );

                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "avgMagGradU",
                    acquisition->fields->avgMagGradU
                );

                writeVectorToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "avgMagUU",
                    acquisition->fields->avgMagUU
                );
            }
        }

    }

    if(io->phaseAveraging)
    {
        writeVectorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "phAvgU",
            acquisition->fields->pAvgU
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "phAvgP",
            acquisition->fields->pAvgP
        );

        writeSymmTensorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "phAvgUU",
            acquisition->fields->pAvgUU
        );

        if(domain->flags.isLesActive)
        {
            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "phAvgNut",
                acquisition->fields->pAvgNut
            );

            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "phAvgCs",
                acquisition->fields->pAvgCs
            );
        }

        if(io->phaseAveraging > 1)
        {
            writeVectorToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "phAvgOmega",
                acquisition->fields->pAvgOmega
            );

            if(io->phaseAveraging > 2)
            {
                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "phAvgPsq",
                    acquisition->fields->pAvgP2
                );

                writeSymmTensorToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "phAvgOmegaOmega",
                    acquisition->fields->pAvgOmegaOmega
                );

                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "phAvgUdotGradP",
                    acquisition->fields->pAvgUdotGradP
                );

                writeScalarToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "phAvgMagGradU",
                    acquisition->fields->pAvgMagGradU
                );

                writeVectorToXMF
                (
                    domain,
                    filexmf,
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    time,
                    "phAvgMagUU",
                    acquisition->fields->pAvgMagUU
                );
            }
        }
    }

    if(io->keBudgets)
    {
        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keErr",
            acquisition->keBudFields->Error
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keEm",
            acquisition->keBudFields->lEm
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keD",
            acquisition->keBudFields->D
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keF",
            acquisition->keBudFields->F
        );

        writeScalarToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keEps",
            acquisition->keBudFields->Eps
        );

        writeSymmTensorToXMF
        (
            domain,
            filexmf,
            hdfileName.c_str(),
            &file_id,
            &dataspace_id,
            time,
            "keUpUp",
            acquisition->keBudFields->lVpVp
        );

        if(flags->isWindFarmActive)
        {
            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "keFarm",
                acquisition->keBudFields->Pf
            );
        }
    }

    if(flags->isAquisitionActive)
    {
        if(domain->acquisition->isAverage3LMActive)
        {
            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "dTdz3LM",
                acquisition->LM3->avgdTdz
            );
        }

        if(domain->acquisition->isPerturbABLActive)
        {
            writeVectorToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "UpABL",
                acquisition->perturbABL->pertU
            );

            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "PpABL",
                acquisition->perturbABL->pertP
            );

            writeScalarToXMF
            (
                domain,
                filexmf,
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                time,
                "TpABL",
                acquisition->perturbABL->pertT
            );
        }
    }

    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    // close this time section in the XMF file
    xmfWriteFileEndTimeSection(xmf, filexmf);

    return(0);
}

//***************************************************************************************************************//

void setToZero(float *vec, PetscInt n)
{
    for(PetscInt p=0; p<n; p++)
    {
        vec[p] = 0.0;
    }

    return;
}

//***************************************************************************************************************//

void xmfWriteFileStartTimeSection
(
    FILE *Xmf, const char *FileXmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *Topology,
    PetscReal Time
)
{
    // Topology - mesh topology, 3DSMesh - 3d structured mesh
	Xmf = fopen(FileXmf, "a");
	fprintf(Xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
	fprintf(Xmf, "	   <Time Value=\"%03f\"/>\n", Time);
	fprintf(Xmf, "     <Topology TopologyType=\"%s\" NumberOfElements=\"%ld %ld %ld\"/>\n",Topology, Size_z, Size_y, Size_x);
    fclose(Xmf);

    return;
}

//***************************************************************************************************************//

void xmfWriteFileEndTimeSection
(
    FILE *Xmf, const char *FileXmf
)
{
	Xmf = fopen(FileXmf, "a");
    fprintf(Xmf, "   </Grid>\n");
    fclose(Xmf);

    return;
}

//***************************************************************************************************************//

void xmfWriteFileGeometry
(
    FILE *Xmf, const char *FileXmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *PathSave
)
{
	Xmf = fopen(FileXmf, "a");
	fprintf(Xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(Xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Size_z, Size_y, Size_x);
    fprintf(Xmf, "        %s:/X\n", PathSave);
    fprintf(Xmf, "       </DataItem>\n");
    fprintf(Xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Size_z, Size_y, Size_x);
    fprintf(Xmf, "        %s:/Y\n", PathSave);
    fprintf(Xmf, "       </DataItem>\n");
    fprintf(Xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Size_z, Size_y, Size_x);
    fprintf(Xmf, "        %s:/Z\n", PathSave);
    fprintf(Xmf, "       </DataItem>\n");
    fprintf(Xmf, "     </Geometry>\n");
	fclose(Xmf);

	return;
}

//***************************************************************************************************************//

void xmfWriteFileSymmTensor
(
    FILE *xmf, const char *filexmf,
    PetscInt size_x, PetscInt size_y, PetscInt size_z,
    const char *PathSave, const char *symmTensorName,
    const char * XX, const char * YY, const char *ZZ,
    const char * XY, const char * XZ, const char *YZ,
    const char *center
)
{

	// center  - Cell or Node
	// Vecname - Vector name as shown in paraview
	// XX,XY,XZ,YY,YZ,ZZ: tensor component names as saved in hdf5 file
	// pathsave = relative path with respect to the saved files

	xmf = fopen(filexmf, "a");
    fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Tensor6\" Center=\"%s\">\n", symmTensorName, center);
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 6\" Function=\"JOIN($0, $1, $2, $3, $4, $5)\" ItemType=\"Function\">\n", size_z, size_y, size_x);
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, XX);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, XY);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, XZ);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, YY);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, YZ);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, ZZ);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fclose(xmf);

    return;
}

//***************************************************************************************************************//

void xmfWriteFileVector
(
    FILE *xmf, const char *filexmf,
    PetscInt size_x, PetscInt size_y, PetscInt size_z,
    const char *PathSave, const char *Vecname,
    const char * V1, const char * V2, const char *V3,
    const char *center
)
{
	// center  - Cell or Node
	// Vecname - Vector name as shown in paraview
	// V1, V2 ,V3 vector component names as saved in hdf5 file
	// pathsave = relative path with respect to the saved files

	xmf = fopen(filexmf, "a");
    fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"%s\">\n", Vecname, center);
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 3\" Function=\"JOIN($0, $1, $2)\" ItemType=\"Function\">\n", size_z, size_y, size_x);
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, V1);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, V2);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%ld %ld %ld 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\" ItemType=\"Uniform\">\n", size_z, size_y, size_x);
	fprintf(xmf, "        %s:/%s\n", PathSave, V3);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fclose(xmf);

    return;
}

//***************************************************************************************************************//

void xmfWriteFileScalar
(
    FILE *Xmf, const char *Filexmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *PathSave, const char *ScalName,
    const char *Scal,
    const char *Center
)
{
	// Center - Cell or Node
	// ScalName - Scalar name as shown in paraview
	// Scal - Scalar name as saved in hdf5 file
	// PathSave = relative path with respect to the saved files

	Xmf = fopen(Filexmf, "a");
    fprintf(Xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"%s\">\n", ScalName, Center);
    fprintf(Xmf, "       <DataItem Dimensions=\"%ld %ld %ld\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", Size_z, Size_y, Size_x);
	fprintf(Xmf, "        %s:/%s\n", PathSave, Scal);
    fprintf(Xmf, "       </DataItem>\n");
    fprintf(Xmf, "     </Attribute>\n");
    fclose(Xmf);

    return;
}

//***************************************************************************************************************//

void hdfWriteDataset
(
    hid_t *file_id,
    hid_t *dataspace_id,
    char const *var,
    float *x
)
{

	hid_t  dataset_id = H5Dcreate(*file_id, var, H5T_NATIVE_FLOAT, *dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status     = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);

	status            = H5Dclose(dataset_id);

	return;

}

//***************************************************************************************************************//

PetscErrorCode writePointsToXMF
(
    mesh_      *mesh,        // domain mesh
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time          // time
)
{
    DM                fda  = mesh->fda, da = mesh->da;
    DMDALocalInfo     info = mesh->info;

    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    Vec		      Coor;

    Cmpnts        ***coor;

    FILE          *xmf;

    // create 1D vector which stores the 3D points
    float *x; x = new float[(mx-1)*(my-1)*(mz-1)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // write x component
    for (k=zs; k<ze-1; k++)
    for (j=ys; j<ye-1; j++)
    for (i=xs; i<xe-1; i++)
    {
        x[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].x;
    }
    hdfWriteDataset(file_id, dataspace_id, "/X", x);

    // write y component
    for (k=zs; k<ze-1; k++)
    for (j=ys; j<ye-1; j++)
    for (i=xs; i<xe-1; i++)
    {
        x[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].y;
    }
    hdfWriteDataset(file_id, dataspace_id, "/Y", x);

    // write z component
    for (k=zs; k<ze-1; k++)
    for (j=ys; j<ye-1; j++)
    for (i=xs; i<xe-1; i++)
    {
        x[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].z;
    }
    hdfWriteDataset(file_id, dataspace_id, "/Z", x);

    // write reference to HDF file in XMF file
    xmfWriteFileGeometry(xmf, filexmf, mx-1, my-1, mz-1, hdfilen);

    DMDAVecRestoreArray(fda, Coor, &coor);

    delete [] x;

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode writeUserSectionPointsToXMF
(
    mesh_      *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    uSections  *uSection    // pointer the user defined section
)
{
    FILE          *xmf;
    PetscInt      i, j;

    PetscMPIInt   rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    PetscInt nx = uSection->nx, ny = uSection->ny;

    // create 1D vectors which stores the 2D points
    float *lx; lx = new float[(nx)*(ny)];
    float *gx; gx = new float[(nx)*(ny)];

    // write x component
    setToZero(lx, (nx)*(ny));
    setToZero(gx, (nx)*(ny));

    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        if(uSection->hasCoor[j][i])
        {
            lx[j*(nx) + i] = uSection->coor[j][i].x;
        }
    }

    MPI_Reduce(lx, gx, (nx)*(ny), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/X", gx);

    // write y component
    setToZero(lx, (nx)*(ny));
    setToZero(gx, (nx)*(ny));

    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        if(uSection->hasCoor[j][i])
        {
            lx[j*(nx) + i] = uSection->coor[j][i].y;
        }
    }

    MPI_Reduce(lx, gx, (nx)*(ny), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Y", gx);

    // write z component
    setToZero(lx, (nx)*(ny));
    setToZero(gx, (nx)*(ny));

    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        if(uSection->hasCoor[j][i])
        {
            lx[j*(nx) + i] = uSection->coor[j][i].z;
        }
    }

    MPI_Reduce(lx, gx, (nx)*(ny), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Z", gx);

    // write reference to HDF file in XMF file
    if(!rank) xmfWriteFileGeometry(xmf, filexmf, nx, ny, 1, hdfilen);

    delete [] lx;
    delete [] gx;

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode writeKSectionPointsToXMF
(
    mesh_      *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    PetscInt   kIndex        // index of the k-section
)
{
    DM            fda = mesh->fda, da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      surfaceCuts = 0;

    Vec		        Coor;

    Cmpnts        ***coor;

    FILE          *xmf;

    PetscMPIInt   rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // test if the surface cuts this proc
    if
    (
        (zs <= kIndex-1 && ze > kIndex-1) ||
        (zs <= kIndex   && ze > kIndex)
    )
    {
        surfaceCuts = 1;
    }

    // create 1D vectors which stores the 2D points
    float *lx; lx = new float[(mx-1)*(my-1)];
    float *gx; gx = new float[(mx-1)*(my-1)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // write x component
    setToZero(lx, (mx-1)*(my-1));
    setToZero(gx, (mx-1)*(my-1));

    for (j=ys; j<lye; j++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[j*(mx-1) + i] = 0.5 * (coor[kIndex-1][j][i].x + coor[kIndex][j][i].x);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(my-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/X", gx);

    // write y component
	  setToZero(lx, (mx-1)*(my-1));
    setToZero(gx, (mx-1)*(my-1));

    for (j=ys; j<lye; j++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[j*(mx-1) + i] = 0.5 * (coor[kIndex-1][j][i].y + coor[kIndex][j][i].y);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(my-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Y", gx);

    // write z component
	  setToZero(lx, (mx-1)*(my-1));
    setToZero(gx, (mx-1)*(my-1));

    for (j=ys; j<lye; j++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[j*(mx-1) + i] = 0.5 * (coor[kIndex-1][j][i].z + coor[kIndex][j][i].z);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(my-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Z", gx);

    // write reference to HDF file in XMF file
    if(!rank) xmfWriteFileGeometry(xmf, filexmf, mx-1, my-1, 1, hdfilen);

    DMDAVecRestoreArray(fda, Coor, &coor);

    delete [] lx;
    delete [] gx;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeJSectionPointsToXMF
(
    mesh_      *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    PetscInt   jIndex        // index of the k-section
)
{
    DM            fda = mesh->fda, da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      surfaceCuts = 0;

    Vec		        Coor;

    Cmpnts        ***coor;

    FILE          *xmf;

    PetscMPIInt   rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // test if the surface cuts this proc
    if
    (
        (ys <= jIndex-1 && ye > jIndex-1) ||
        (ys <= jIndex   && ye > jIndex)
    )
    {
        surfaceCuts = 1;
    }

    // create 1D vectors which stores the 2D points
    float *lx; lx = new float[(mx-1)*(mz-1)];
    float *gx; gx = new float[(mx-1)*(mz-1)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // write x component
    setToZero(lx, (mx-1)*(mz-1));
    setToZero(gx, (mx-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[k*(mx-1) + i] = 0.5 * (coor[k][jIndex-1][i].x + coor[k][jIndex][i].x);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/X", gx);

    // write y component
    setToZero(lx, (mx-1)*(mz-1));
    setToZero(gx, (mx-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[k*(mx-1) + i] = 0.5 * (coor[k][jIndex-1][i].y + coor[k][jIndex][i].y);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Y", gx);

    // write z component
    setToZero(lx, (mx-1)*(mz-1));
    setToZero(gx, (mx-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (i=xs; i<lxe; i++)
    {
        if(surfaceCuts)
        {
            lx[k*(mx-1) + i] = 0.5 * (coor[k][jIndex-1][i].z + coor[k][jIndex][i].z);
        }
    }

    MPI_Reduce(lx, gx, (mx-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Z", gx);

    // write reference to HDF file in XMF file
    if(!rank) xmfWriteFileGeometry(xmf, filexmf, mx-1, 1, mz-1, hdfilen);

    DMDAVecRestoreArray(fda, Coor, &coor);

    delete [] lx;
    delete [] gx;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeISectionPointsToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    PetscInt        iIndex        // index of the k-section
)
{
    DM            fda = mesh->fda, da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      surfaceCuts = 0;

    Vec		        Coor;

    Cmpnts        ***coor;

    FILE          *xmf;

    PetscMPIInt   rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // test if the surface cuts this proc
    if
    (
        (xs <= iIndex-1 && xe > iIndex-1) ||
        (xs <= iIndex   && xe > iIndex)
    )
    {
        surfaceCuts = 1;
    }

    // create 1D vectors which stores the 2D points
	float *lx = (float*)malloc((my-1)*(mz-1)*sizeof(float));
	float *gx = (float*)malloc((my-1)*(mz-1)*sizeof(float));

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // write x component
    setToZero(lx, (my-1)*(mz-1));
    setToZero(gx, (my-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (j=ys; j<lye; j++)
    {
        if(surfaceCuts)
        {
            lx[k*(my-1) + j] = 0.5 * (coor[k][j][iIndex-1].x + coor[k][j][iIndex].x);
        }
    }

    MPI_Reduce(lx, gx, (my-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/X", gx);

    // write y component
    setToZero(lx, (my-1)*(mz-1));
    setToZero(gx, (my-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (j=ys; j<lye; j++)
    {
        if(surfaceCuts)
        {
            lx[k*(my-1) + j] = 0.5 * (coor[k][j][iIndex-1].y + coor[k][j][iIndex].y);
        }
    }

    MPI_Reduce(lx, gx, (my-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Y", gx);

    // write z component
    setToZero(lx, (my-1)*(mz-1));
    setToZero(gx, (my-1)*(mz-1));

    for (k=zs; k<lze; k++)
    for (j=ys; j<lye; j++)
    {
        if(surfaceCuts)
        {
            lx[k*(my-1) + j] = 0.5 * (coor[k][j][iIndex-1].z + coor[k][j][iIndex].z);
        }
    }

    MPI_Reduce(lx, gx, (my-1)*(mz-1), MPI_FLOAT, MPI_SUM, 0, PETSC_COMM_WORLD);

    if(!rank) hdfWriteDataset(file_id, dataspace_id, "/Z", gx);

    // write reference to HDF file in XMF file
    if(!rank) xmfWriteFileGeometry(xmf, filexmf, 1, my-1, mz-1, hdfilen);

    DMDAVecRestoreArray(fda, Coor, &coor);

    free(lx);
    free(gx);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeSymmTensorToXMF
(
    domain_    *domain,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Vec        V             // Petsc vector to write
)
{
    DM            sda = domain->mesh->sda;
    DMDALocalInfo info = domain->mesh->info;

    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    symmTensor    ***v;

    FILE          *xmf;

    float         *x; x = new float [(mx-2)*(my-2)*(mz-2)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1xx[256], str1yy[256], str1zz[256],
         str1xy[256], str1xz[256], str1yz[256],
         str2xx[256], str2yy[256], str2zz[256],
         str2xy[256], str2xz[256], str2yz[256];

    sprintf(str1xx, "/%s_xx",fieldName);
    sprintf(str1yy, "/%s_yy",fieldName);
    sprintf(str1zz, "/%s_zz",fieldName);
    sprintf(str1xy, "/%s_xy",fieldName);
    sprintf(str1xz, "/%s_xz",fieldName);
    sprintf(str1yz, "/%s_yz",fieldName);

    sprintf(str2xx, "%s_xx",fieldName);
    sprintf(str2yy, "%s_yy",fieldName);
    sprintf(str2zz, "%s_zz",fieldName);
    sprintf(str2xy, "%s_xy",fieldName);
    sprintf(str2xz, "%s_xz",fieldName);
    sprintf(str2yz, "%s_yz",fieldName);

    DMDAVecGetArray(sda, V, &v);

    // write xx component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].xx;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xx, x);

    // write yy component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].yy;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yy, x);

    // write zz component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].zz;
    }
    hdfWriteDataset(file_id, dataspace_id, str1zz, x);

    // write xy component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].xy;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xy, x);

    // write xz component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].xz;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xz, x);

    // write yz component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].yz;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yz, x);

    // write reference to HDF file in XMF file
    xmfWriteFileSymmTensor(xmf, filexmf, mx-2, my-2, mz-2, hdfilen, fieldName, str2xx, str2yy, str2zz, str2xy, str2xz, str2yz);

    DMDAVecRestoreArray(sda, V, &v);

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeVectorToXMF
(
    domain_    *domain,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Vec        V             // Petsc vector to write
)
{

    DM            fda = domain->mesh->fda;
    DMDALocalInfo info = domain->mesh->info;

    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***v;

    FILE          *xmf;

    float         *x; x = new float [(mx-2)*(my-2)*(mz-2)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1x[256], str1y[256], str1z[256],
         str2x[256], str2y[256], str2z[256];

    sprintf(str1x, "/%s_x",fieldName);
    sprintf(str1y, "/%s_y",fieldName);
    sprintf(str1z, "/%s_z",fieldName);

    sprintf(str2x, "%s_x",fieldName);
    sprintf(str2y, "%s_y",fieldName);
    sprintf(str2z, "%s_z",fieldName);

    DMDAVecGetArray(fda, V, &v);

    // write x component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].x;
    }
    hdfWriteDataset(file_id, dataspace_id, str1x, x);

    // write y component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].y;
    }
    hdfWriteDataset(file_id, dataspace_id, str1y, x);

    // write z component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1].z;
    }
    hdfWriteDataset(file_id, dataspace_id, str1z, x);

    // write reference to HDF file in XMF file
    xmfWriteFileVector(xmf, filexmf, mx-2, my-2, mz-2, hdfilen, fieldName, str2x, str2y, str2z);

    DMDAVecRestoreArray(fda, V, &v);

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeUserSectionVectorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    uSections  *uSection    // pointer to the user defined section
)
{
    PetscInt nx = uSection->nx, ny = uSection->ny;

    PetscInt      i, j;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(nx)*(ny)];
    Cmpnts **field = uSection->vectorSec;

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1x[256], str1y[256], str1z[256],
         str2x[256], str2y[256], str2z[256];

    sprintf(str1x, "/%s_x",fieldName);
    sprintf(str1y, "/%s_y",fieldName);
    sprintf(str1z, "/%s_z",fieldName);

    sprintf(str2x, "%s_x",fieldName);
    sprintf(str2y, "%s_y",fieldName);
    sprintf(str2z, "%s_z",fieldName);

    // write x component
    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        x[j * (nx) + i] = field[j][i].x;
    }

    hdfWriteDataset(file_id, dataspace_id, str1x, x);

    // write y component
    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        x[j * (nx) + i] = field[j][i].y;
    }
    hdfWriteDataset(file_id, dataspace_id, str1y, x);

    // write z component
    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        x[j * (nx) + i] = field[j][i].z;
    }
    hdfWriteDataset(file_id, dataspace_id, str1z, x);

    // write reference to HDF file in XMF file
    xmfWriteFileVector(xmf, filexmf, nx, ny, 1, hdfilen, fieldName, str2x, str2y, str2z, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode writeKSectionVectorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Cmpnts     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda  = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(my-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1x[256], str1y[256], str1z[256],
         str2x[256], str2y[256], str2z[256];

    sprintf(str1x, "/%s_x",fieldName);
    sprintf(str1y, "/%s_y",fieldName);
    sprintf(str1z, "/%s_z",fieldName);

    sprintf(str2x, "%s_x",fieldName);
    sprintf(str2y, "%s_y",fieldName);
    sprintf(str2z, "%s_z",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write x component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].x;
        else        VW = field[j+1][i  ].x;

        if(i==mx-2) VE = field[j+1][i  ].x;
        else        VE = field[j+1][i+1].x;

        if(j==0)    VS = field[j+1][i+1].x;
        else        VS = field[j  ][i+1].x;

        if(j==my-2) VN = field[j  ][i+1].x;
        else        VN = field[j+1][i+1].x;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1x, x);

    // write y component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].y;
        else        VW = field[j+1][i  ].y;

        if(i==mx-2) VE = field[j+1][i  ].y;
        else        VE = field[j+1][i+1].y;

        if(j==0)    VS = field[j+1][i+1].y;
        else        VS = field[j  ][i+1].y;

        if(j==my-2) VN = field[j  ][i+1].y;
        else        VN = field[j+1][i+1].y;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1y, x);

    // write z component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].z;
        else        VW = field[j+1][i  ].z;

        if(i==mx-2) VE = field[j+1][i  ].z;
        else        VE = field[j+1][i+1].z;

        if(j==0)    VS = field[j+1][i+1].z;
        else        VS = field[j  ][i+1].z;

        if(j==my-2) VN = field[j  ][i+1].z;
        else        VN = field[j+1][i+1].z;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1z, x);

    // write reference to HDF file in XMF file
    xmfWriteFileVector(xmf, filexmf, mx-1, my-1, 1, hdfilen, fieldName, str2x, str2y, str2z, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeJSectionVectorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Cmpnts     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda  = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(mz-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1x[256], str1y[256], str1z[256],
         str2x[256], str2y[256], str2z[256];

    sprintf(str1x, "/%s_x",fieldName);
    sprintf(str1y, "/%s_y",fieldName);
    sprintf(str1z, "/%s_z",fieldName);

    sprintf(str2x, "%s_x",fieldName);
    sprintf(str2y, "%s_y",fieldName);
    sprintf(str2z, "%s_z",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write x component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].x;
        else        VW = field[k+1][i  ].x;

        if(i==mx-2) VE = field[k+1][i  ].x;
        else        VE = field[k+1][i+1].x;

        if(k==0)    VS = field[k+1][i+1].x;
        else        VS = field[k  ][i+1].x;

        if(k==mz-2) VN = field[k  ][i+1].x;
        else        VN = field[k+1][i+1].x;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1x, x);

    // write y component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].y;
        else        VW = field[k+1][i  ].y;

        if(i==mx-2) VE = field[k+1][i  ].y;
        else        VE = field[k+1][i+1].y;

        if(k==0)    VS = field[k+1][i+1].y;
        else        VS = field[k  ][i+1].y;

        if(k==mz-2) VN = field[k  ][i+1].y;
        else        VN = field[k+1][i+1].y;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1y, x);

    // write z component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].z;
        else        VW = field[k+1][i  ].z;

        if(i==mx-2) VE = field[k+1][i  ].z;
        else        VE = field[k+1][i+1].z;

        if(k==0)    VS = field[k+1][i+1].z;
        else        VS = field[k  ][i+1].z;

        if(k==mz-2) VN = field[k  ][i+1].z;
        else        VN = field[k+1][i+1].z;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1z, x);

    // write reference to HDF file in XMF file
    xmfWriteFileVector(xmf, filexmf, mx-1, 1, mz-1, hdfilen, fieldName, str2x, str2y, str2z, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeISectionVectorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Cmpnts     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x = (float*)malloc((my-1)*(mz-1)*sizeof(float));

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1x[256], str1y[256], str1z[256],
         str2x[256], str2y[256], str2z[256];

    sprintf(str1x, "/%s_x",fieldName);
    sprintf(str1y, "/%s_y",fieldName);
    sprintf(str1z, "/%s_z",fieldName);

    sprintf(str2x, "%s_x",fieldName);
    sprintf(str2y, "%s_y",fieldName);
    sprintf(str2z, "%s_z",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write x component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].x;
        else        VW = field[k+1][j  ].x;

        if(j==my-2) VE = field[k+1][j  ].x;
        else        VE = field[k+1][j+1].x;

        if(k==0)    VS = field[k+1][j+1].x;
        else        VS = field[k  ][j+1].x;

        if(k==mz-2) VN = field[k  ][j+1].x;
        else        VN = field[k+1][j+1].x;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1x, x);

    // write y component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].y;
        else        VW = field[k+1][j  ].y;

        if(j==my-2) VE = field[k+1][j  ].y;
        else        VE = field[k+1][j+1].y;

        if(k==0)    VS = field[k+1][j+1].y;
        else        VS = field[k  ][j+1].y;

        if(k==mz-2) VN = field[k  ][j+1].y;
        else        VN = field[k+1][j+1].y;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1y, x);

    // write z component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].z;
        else        VW = field[k+1][j  ].z;

        if(j==my-2) VE = field[k+1][j  ].z;
        else        VE = field[k+1][j+1].z;

        if(k==0)    VS = field[k+1][j+1].z;
        else        VS = field[k  ][j+1].z;

        if(k==mz-2) VN = field[k  ][j+1].z;
        else        VN = field[k+1][j+1].z;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1z, x);

    // write reference to HDF file in XMF file
    xmfWriteFileVector(xmf, filexmf, 1, my-1, mz-1, hdfilen, fieldName, str2x, str2y, str2z, "Node");

    free(x);

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode writeKSectionSymmTensorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    symmTensor     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda  = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(my-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1xx[256], str1yy[256], str1zz[256],
         str1xy[256], str1xz[256], str1yz[256],
         str2xx[256], str2yy[256], str2zz[256],
         str2xy[256], str2xz[256], str2yz[256];

    sprintf(str1xx, "/%s_xx",fieldName);
    sprintf(str1yy, "/%s_yy",fieldName);
    sprintf(str1zz, "/%s_zz",fieldName);
    sprintf(str1xy, "/%s_xy",fieldName);
    sprintf(str1xz, "/%s_xz",fieldName);
    sprintf(str1yz, "/%s_yz",fieldName);

    sprintf(str2xx, "%s_xx",fieldName);
    sprintf(str2yy, "%s_yy",fieldName);
    sprintf(str2zz, "%s_zz",fieldName);
    sprintf(str2xy, "%s_xy",fieldName);
    sprintf(str2xz, "%s_xz",fieldName);
    sprintf(str2yz, "%s_yz",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write xx component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].xx;
        else        VW = field[j+1][i  ].xx;

        if(i==mx-2) VE = field[j+1][i  ].xx;
        else        VE = field[j+1][i+1].xx;

        if(j==0)    VS = field[j+1][i+1].xx;
        else        VS = field[j  ][i+1].xx;

        if(j==my-2) VN = field[j  ][i+1].xx;
        else        VN = field[j+1][i+1].xx;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xx, x);

    // write yy component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].yy;
        else        VW = field[j+1][i  ].yy;

        if(i==mx-2) VE = field[j+1][i  ].yy;
        else        VE = field[j+1][i+1].yy;

        if(j==0)    VS = field[j+1][i+1].yy;
        else        VS = field[j  ][i+1].yy;

        if(j==my-2) VN = field[j  ][i+1].yy;
        else        VN = field[j+1][i+1].yy;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yy, x);

    // write z component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].zz;
        else        VW = field[j+1][i  ].zz;

        if(i==mx-2) VE = field[j+1][i  ].zz;
        else        VE = field[j+1][i+1].zz;

        if(j==0)    VS = field[j+1][i+1].zz;
        else        VS = field[j  ][i+1].zz;

        if(j==my-2) VN = field[j  ][i+1].zz;
        else        VN = field[j+1][i+1].zz;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1zz, x);

    // write xx component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].xy;
        else        VW = field[j+1][i  ].xy;

        if(i==mx-2) VE = field[j+1][i  ].xy;
        else        VE = field[j+1][i+1].xy;

        if(j==0)    VS = field[j+1][i+1].xy;
        else        VS = field[j  ][i+1].xy;

        if(j==my-2) VN = field[j  ][i+1].xy;
        else        VN = field[j+1][i+1].xy;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xy, x);

    // write xx component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].xz;
        else        VW = field[j+1][i  ].xz;

        if(i==mx-2) VE = field[j+1][i  ].xz;
        else        VE = field[j+1][i+1].xz;

        if(j==0)    VS = field[j+1][i+1].xz;
        else        VS = field[j  ][i+1].xz;

        if(j==my-2) VN = field[j  ][i+1].xz;
        else        VN = field[j+1][i+1].xz;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xz, x);

    // write xx component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1].yz;
        else        VW = field[j+1][i  ].yz;

        if(i==mx-2) VE = field[j+1][i  ].yz;
        else        VE = field[j+1][i+1].yz;

        if(j==0)    VS = field[j+1][i+1].yz;
        else        VS = field[j  ][i+1].yz;

        if(j==my-2) VN = field[j  ][i+1].yz;
        else        VN = field[j+1][i+1].yz;

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yz, x);

    // write reference to HDF file in XMF file
    xmfWriteFileSymmTensor(xmf, filexmf, mx-1, my-1, 1, hdfilen, fieldName, str2xx, str2yy, str2zz, str2xy, str2xz, str2yz, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeJSectionSymmTensorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    symmTensor     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda  = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(mz-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1xx[256], str1yy[256], str1zz[256],
         str1xy[256], str1xz[256], str1yz[256],
         str2xx[256], str2yy[256], str2zz[256],
         str2xy[256], str2xz[256], str2yz[256];

    sprintf(str1xx, "/%s_xx",fieldName);
    sprintf(str1yy, "/%s_yy",fieldName);
    sprintf(str1zz, "/%s_zz",fieldName);
    sprintf(str1xy, "/%s_xy",fieldName);
    sprintf(str1xz, "/%s_xz",fieldName);
    sprintf(str1yz, "/%s_yz",fieldName);

    sprintf(str2xx, "%s_xx",fieldName);
    sprintf(str2yy, "%s_yy",fieldName);
    sprintf(str2zz, "%s_zz",fieldName);
    sprintf(str2xy, "%s_xy",fieldName);
    sprintf(str2xz, "%s_xz",fieldName);
    sprintf(str2yz, "%s_yz",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write xx component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].xx;
        else        VW = field[k+1][i  ].xx;

        if(i==mx-2) VE = field[k+1][i  ].xx;
        else        VE = field[k+1][i+1].xx;

        if(k==0)    VS = field[k+1][i+1].xx;
        else        VS = field[k  ][i+1].xx;

        if(k==mz-2) VN = field[k  ][i+1].xx;
        else        VN = field[k+1][i+1].xx;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xx, x);

    // write yy component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].yy;
        else        VW = field[k+1][i  ].yy;

        if(i==mx-2) VE = field[k+1][i  ].yy;
        else        VE = field[k+1][i+1].yy;

        if(k==0)    VS = field[k+1][i+1].yy;
        else        VS = field[k  ][i+1].yy;

        if(k==mz-2) VN = field[k  ][i+1].yy;
        else        VN = field[k+1][i+1].yy;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1yy, x);

    // write zz component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].zz;
        else        VW = field[k+1][i  ].zz;

        if(i==mx-2) VE = field[k+1][i  ].zz;
        else        VE = field[k+1][i+1].zz;

        if(k==0)    VS = field[k+1][i+1].zz;
        else        VS = field[k  ][i+1].zz;

        if(k==mz-2) VN = field[k  ][i+1].zz;
        else        VN = field[k+1][i+1].zz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1zz, x);

    // write xy component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].xy;
        else        VW = field[k+1][i  ].xy;

        if(i==mx-2) VE = field[k+1][i  ].xy;
        else        VE = field[k+1][i+1].xy;

        if(k==0)    VS = field[k+1][i+1].xy;
        else        VS = field[k  ][i+1].xy;

        if(k==mz-2) VN = field[k  ][i+1].xy;
        else        VN = field[k+1][i+1].xy;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xy, x);

    // write xz component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].xz;
        else        VW = field[k+1][i  ].xz;

        if(i==mx-2) VE = field[k+1][i  ].xz;
        else        VE = field[k+1][i+1].xz;

        if(k==0)    VS = field[k+1][i+1].xz;
        else        VS = field[k  ][i+1].xz;

        if(k==mz-2) VN = field[k  ][i+1].xz;
        else        VN = field[k+1][i+1].xz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1xz, x);

    // write yz component
    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1].yz;
        else        VW = field[k+1][i  ].yz;

        if(i==mx-2) VE = field[k+1][i  ].yz;
        else        VE = field[k+1][i+1].yz;

        if(k==0)    VS = field[k+1][i+1].yz;
        else        VS = field[k  ][i+1].yz;

        if(k==mz-2) VN = field[k  ][i+1].yz;
        else        VN = field[k+1][i+1].yz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }

    hdfWriteDataset(file_id, dataspace_id, str1yz, x);

    // write reference to HDF file in XMF file
    xmfWriteFileSymmTensor(xmf, filexmf, mx-1, 1, mz-1, hdfilen, fieldName, str2xx, str2yy, str2zz, str2xy, str2xz, str2yz, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeISectionSymmTensorToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    symmTensor     **field       // 2D array containing the field data
)
{
    // note: the field data is contained inside user->ucat_plane_ks
    //       and must be loaded prior to call this function
    DM            fda = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x = (float*)malloc((my-1)*(mz-1)*sizeof(float));

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1xx[256], str1yy[256], str1zz[256],
         str1xy[256], str1xz[256], str1yz[256],
         str2xx[256], str2yy[256], str2zz[256],
         str2xy[256], str2xz[256], str2yz[256];

    sprintf(str1xx, "/%s_xx",fieldName);
    sprintf(str1yy, "/%s_yy",fieldName);
    sprintf(str1zz, "/%s_zz",fieldName);
    sprintf(str1xy, "/%s_xy",fieldName);
    sprintf(str1xz, "/%s_xz",fieldName);
    sprintf(str1yz, "/%s_yz",fieldName);

    sprintf(str2xx, "%s_xx",fieldName);
    sprintf(str2yy, "%s_yy",fieldName);
    sprintf(str2zz, "%s_zz",fieldName);
    sprintf(str2xy, "%s_xy",fieldName);
    sprintf(str2xz, "%s_xz",fieldName);
    sprintf(str2yz, "%s_yz",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write xx component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].xx;
        else        VW = field[k+1][j  ].xx;

        if(j==my-2) VE = field[k+1][j  ].xx;
        else        VE = field[k+1][j+1].xx;

        if(k==0)    VS = field[k+1][j+1].xx;
        else        VS = field[k  ][j+1].xx;

        if(k==mz-2) VN = field[k  ][j+1].xx;
        else        VN = field[k+1][j+1].xx;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xx, x);

    // write yy component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].yy;
        else        VW = field[k+1][j  ].yy;

        if(j==my-2) VE = field[k+1][j  ].yy;
        else        VE = field[k+1][j+1].yy;

        if(k==0)    VS = field[k+1][j+1].yy;
        else        VS = field[k  ][j+1].yy;

        if(k==mz-2) VN = field[k  ][j+1].yy;
        else        VN = field[k+1][j+1].yy;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yy, x);

    // write zz component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].zz;
        else        VW = field[k+1][j  ].zz;

        if(j==my-2) VE = field[k+1][j  ].zz;
        else        VE = field[k+1][j+1].zz;

        if(k==0)    VS = field[k+1][j+1].zz;
        else        VS = field[k  ][j+1].zz;

        if(k==mz-2) VN = field[k  ][j+1].zz;
        else        VN = field[k+1][j+1].zz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1zz, x);

    // write xy component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].xy;
        else        VW = field[k+1][j  ].xy;

        if(j==my-2) VE = field[k+1][j  ].xy;
        else        VE = field[k+1][j+1].xy;

        if(k==0)    VS = field[k+1][j+1].xy;
        else        VS = field[k  ][j+1].xy;

        if(k==mz-2) VN = field[k  ][j+1].xy;
        else        VN = field[k+1][j+1].xy;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xy, x);

    // write xz component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].xz;
        else        VW = field[k+1][j  ].xz;

        if(j==my-2) VE = field[k+1][j  ].xz;
        else        VE = field[k+1][j+1].xz;

        if(k==0)    VS = field[k+1][j+1].yz;
        else        VS = field[k  ][j+1].yz;

        if(k==mz-2) VN = field[k  ][j+1].xz;
        else        VN = field[k+1][j+1].xz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1xz, x);

    // write yz component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1].yz;
        else        VW = field[k+1][j  ].yz;

        if(j==my-2) VE = field[k+1][j  ].yz;
        else        VE = field[k+1][j+1].yz;

        if(k==0)    VS = field[k+1][j+1].yz;
        else        VS = field[k  ][j+1].yz;

        if(k==mz-2) VN = field[k  ][j+1].yz;
        else        VN = field[k+1][j+1].yz;

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1yz, x);

    // write reference to HDF file in XMF file
    xmfWriteFileSymmTensor(xmf, filexmf, 1, my-1, mz-1, hdfilen, fieldName, str2xx, str2yy, str2zz, str2xy, str2xz, str2yz, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeScalarToXMF
(
    domain_    *domain,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	     *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    Vec        V             // Petsc vector to write
)
{
    DM            da = domain->mesh->da;
    DMDALocalInfo info = domain->mesh->info;

    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***v;

    FILE          *xmf;

    float         *x; x = new float [(mx-2)*(my-2)*(mz-2)];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1[256], str2[256];
    sprintf(str1, "/%s_",fieldName);
    sprintf(str2, "%s_",fieldName);

    DMDAVecGetArray(da, V, &v);

    // write x component
    for (k=zs; k<ze-2; k++)
    for (j=ys; j<ye-2; j++)
    for (i=xs; i<xe-2; i++)
    {
        x[ (k-zs) * (mx - 2)*(my - 2) + (j-ys) * (mx - 2) + (i-xs)] = v[k+1][j+1][i+1];
    }
    hdfWriteDataset(file_id, dataspace_id, str1, x);

    // write reference to HDF file in XMF file
    xmfWriteFileScalar(xmf, filexmf, mx-2, my-2, mz-2, hdfilen, fieldName, str2);

    DMDAVecRestoreArray(da, V, &v);

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeUserSectionScalarToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    uSections  *uSection    // pointer to the user defined section
)
{
    PetscInt nx = uSection->nx, ny = uSection->ny;

    PetscInt      i, j;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(nx)*(ny)];
    PetscReal **field = uSection->scalarSec;

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1[256], str2[256];
    sprintf(str1, "/%s_",fieldName);
    sprintf(str2, "%s_",fieldName);

    // write x component
    for (j=0; j<ny; j++)
    for (i=0; i<nx; i++)
    {
        x[j * (nx) + i] = field[j][i];
    }
    hdfWriteDataset(file_id, dataspace_id, str1, x);

    // write reference to HDF file in XMF file
    xmfWriteFileScalar(xmf, filexmf, nx, ny, 1, hdfilen, fieldName, str2, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode writeKSectionScalarToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    PetscReal     **field       // array containing field data
)
{
    DM            da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    PetscReal     ***v;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(my-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1[256], str2[256];
    sprintf(str1, "/%s_",fieldName);
    sprintf(str2, "%s_",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write x component
    for (j=0; j<my-1; j++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[j+1][i+1];
        else        VW = field[j+1][i  ];

        if(i==mx-2) VE = field[j+1][i  ];
        else        VE = field[j+1][i+1];

        if(j==0)    VS = field[j+1][i+1];
        else        VS = field[j  ][i+1];

        if(j==my-2) VN = field[j  ][i+1];
        else        VN = field[j+1][i+1];

        V = 0.25 * (VN + VS + VW + VW);

        x[j * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1, x);

    // write reference to HDF file in XMF file
    xmfWriteFileScalar(xmf, filexmf, mx-1, my-1, 1, hdfilen, fieldName, str2, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeJSectionScalarToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    PetscReal     **field       // array containing field data
)
{
    DM            da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    PetscReal     ***v;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x; x = new float[(mx-1)*(mz-1)];

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1[256], str2[256];
    sprintf(str1, "/%s_",fieldName);
    sprintf(str2, "%s_",fieldName);

    PetscReal VN, VS, VE, VW, V;

    for (k=0; k<mz-1; k++)
    for (i=0; i<mx-1; i++)
    {
        // build cell-to-node interpolation points
        if(i==0)    VW = field[k+1][i+1];
        else        VW = field[k+1][i  ];

        if(i==mx-2) VE = field[k+1][i  ];
        else        VE = field[k+1][i+1];

        if(k==0)    VS = field[k+1][i+1];
        else        VS = field[k  ][i+1];

        if(k==mz-2) VN = field[k  ][i+1];
        else        VN = field[k+1][i+1];

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (mx-1) + i] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1, x);

    // write reference to HDF file in XMF file
    xmfWriteFileScalar(xmf, filexmf, mx-1, 1, mz-1, hdfilen, fieldName, str2, "Node");

    delete [] x;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeISectionScalarToXMF
(
    mesh_    *mesh,        // user context
    const char *filexmf,     // name of the XMF file to append
    const char *hdfilen,     // name of the HDF file to refer
    hid_t	   *file_id,     // id of the HDF file to refer
    hid_t      *dataspace_id,// id of the HDF dataspace to refer
    PetscReal     time,         // time
    const char *fieldName,   // field name as it appears in ParaView
    PetscReal     **field       // array containing field data
)
{
    DM            da = mesh->da;
    DMDALocalInfo info = mesh->info;

    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    PetscReal     ***v;

    FILE          *xmf;

    // create 1D vector which stores the 2D points
    float *x = (float*)malloc((my-1)*(mz-1)*sizeof(float));

    // set internal variables used to recall the HDF5 field from the XMF file
    char str1[256], str2[256];
    sprintf(str1, "/%s_",fieldName);
    sprintf(str2, "%s_",fieldName);

    PetscReal VN, VS, VE, VW, V;

    // write x component
    for (k=0; k<mz-1; k++)
    for (j=0; j<my-1; j++)
    {
        // build cell-to-node interpolation points
        if(j==0)    VW = field[k+1][j+1];
        else        VW = field[k+1][j  ];

        if(j==my-2) VE = field[k+1][j  ];
        else        VE = field[k+1][j+1];

        if(k==0)    VS = field[k+1][j+1];
        else        VS = field[k  ][j+1];

        if(k==mz-2) VN = field[k  ][j+1];
        else        VN = field[k+1][j+1];

        V = 0.25 * (VN + VS + VW + VW);

        x[k * (my-1) + j] = V;
    }
    hdfWriteDataset(file_id, dataspace_id, str1, x);

    // write reference to HDF file in XMF file
    xmfWriteFileScalar(xmf, filexmf, 1, my-1, mz-1, hdfilen, fieldName, str2, "Node");

    free(x);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeJSectionToRaster
(
    mesh_    *mesh,        // user context
    PetscInt jIndex        // index of the j-section
)
{
    io_ *io = mesh->access->io;
    acquisition_ *acquisition = mesh->access->acquisition;
    ueqn_ *ueqn = mesh->access->ueqn;

    if(io->averaging)
    {
        DM            fda = mesh->fda, da = mesh->da;
        DMDALocalInfo info = mesh->info;

        PetscInt           xs = info.xs, xe = info.xs + info.xm;
        PetscInt           ys = info.ys, ye = info.ys + info.ym;
        PetscInt           zs = info.zs, ze = info.zs + info.zm;
        PetscInt           mx = info.mx, my = info.my, mz = info.mz;

        PetscInt           i, j, k;
        PetscInt           lxs, lxe, lys, lye, lzs, lze;
        PetscInt           surfaceCuts = 0;

        Vec		      Coor, lavgU;

        Cmpnts        ***umean, ***coor, ***cent;

        PetscMPIInt           rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

        lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
        lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
        lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

        // create the local average velocity
        VecDuplicate(ueqn->lUcat, &lavgU);

        // scatter the average velocity to local for the trilinear interpolation
        DMGlobalToLocalBegin(fda, acquisition->fields->avgU, INSERT_VALUES, lavgU);
        DMGlobalToLocalEnd(fda, acquisition->fields->avgU, INSERT_VALUES, lavgU);

        DMGetCoordinatesLocal(da, &Coor);
        DMDAVecGetArray(fda, Coor, &coor);
        DMDAVecGetArray(fda, lavgU, &umean);
        DMDAVecGetArray(fda, mesh->lCent, &cent);

        // raster format only works with uniform spaced mesh. So we have to
        // re-mesh in order things to work also with graded meshes.
        // Note: CARTESIAN MESH IS ASSUMED with direction correspondace as follows:
        //       x: k direction
        //       y: i direction
        //       z: j direction

        // set cell edge dimensions in i and k directions
        PetscReal        dC  = 50.0;

        // remesh the 2D slice
        PetscReal        Lx  = mesh->bounds.Lx;
        PetscReal        Ly  = mesh->bounds.Ly;
        PetscInt      Nx  = floor(Lx / dC) + 1;
        PetscInt      Ny  = floor(Ly / dC) + 1;

        // get the z of the raster plane (mesh can be graded but cartesian)
        PetscReal    rasterZ = cent[zs][jIndex][xs].z;

        std::vector<std::vector<Cpt2D>> rasterMesh(Nx);

        // set the coordinates
        for(k=0; k<Nx; k++)
        {
            rasterMesh[k].resize(Ny);

            for(i=0; i<Ny; i++)
            {
                rasterMesh[k][i].x = mesh->bounds.xmin + k*dC;
                rasterMesh[k][i].y = mesh->bounds.ymin + i*dC;
            }
        }

        // test if the surface cuts this proc
        if
        (
            (ys <= jIndex-1 && ye > jIndex-1) ||
            (ys <= jIndex   && ye > jIndex)
        )
        {
            surfaceCuts = 1;
        }

        // build a global vector and alocal vector of the field, which is zero only
        // in this processor's intersection with the plane.
        std::vector<std::vector<Cmpnts>> lrasterUmean(Nx);
        std::vector<std::vector<Cmpnts>> grasterUmean(Nx);

        // set the coordinates
        for(k=0; k<Nx; k++)
        {
            lrasterUmean[k].resize(Ny);
            grasterUmean[k].resize(Ny);

            for(i=0; i<Ny; i++)
            {
                lrasterUmean[k][i].x = 0.0;
                lrasterUmean[k][i].y = 0.0;
                lrasterUmean[k][i].z = 0.0;

                grasterUmean[k][i].x = 0.0;
                grasterUmean[k][i].y = 0.0;
                grasterUmean[k][i].z = 0.0;
            }
        }

        // now interpolate the values at the raster points
        if(surfaceCuts)
        {

            // this processor x,y bounds
            PetscReal thisPxmin = 0.5 * (coor[zs][jIndex-1][xs].x + coor[zs][jIndex][xs].x);
            PetscReal thisPxmax = 0.5 * (coor[lze-1][jIndex-1][lxe-1].x + coor[lze-1][jIndex][lxe-1].x);
            PetscReal thisPymin = 0.5 * (coor[zs][jIndex-1][xs].y + coor[zs][jIndex][xs].y);
            PetscReal thisPymax = 0.5 * (coor[lze-1][jIndex-1][lxe-1].y + coor[lze-1][jIndex][lxe-1].y);

            for(k=0; k<Nx; k++)
            {
                for(i=0; i<Ny; i++)
                {
                    Cmpnts rasterPt;
                    rasterPt.x = rasterMesh[k][i].x;
                    rasterPt.y = rasterMesh[k][i].y;
                    rasterPt.z = rasterZ;

                    // check that x and y are inside this processor bounds
                    if
                    (
                        rasterPt.x >= thisPxmin && rasterPt.x <= thisPxmax &&
                        rasterPt.y >= thisPymin && rasterPt.y <= thisPymax
                    )
                    {

                        PetscReal  minDist = 1.0e10;
                        cellIds closestLabels;

                        // find the closest point in the CFD mesh
                        for (PetscInt kk=lzs; kk<lze; kk++)
                        for (PetscInt ii=lxs; ii<lxe; ii++)
                        {
                           Cmpnts  distVec = nSub(cent[kk][jIndex][ii], rasterPt);
                           PetscReal  distMag = nMag(distVec);

                            if(distMag < minDist)
                            {
                                minDist = distMag;
                                closestLabels.i = ii;
                                closestLabels.j = jIndex;
                                closestLabels.k = kk;
                            }
                        }

                        // trilinear interpolate
                        Cmpnts rasterU;

                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            rasterPt.x, rasterPt.y, rasterPt.z,
                            closestLabels.i, closestLabels.j, closestLabels.k,
                            cent,
                            umean,
                            rasterU
                        );

                        lrasterUmean[k][i].x = rasterU.x;
                        lrasterUmean[k][i].y = rasterU.y;
                        lrasterUmean[k][i].z = rasterU.z;

                    }
                }
            }
        }

        // parallel sum the local raster fields on master node
        for(k=0; k<Nx; k++)
        {
            MPI_Reduce(&lrasterUmean[k][0], &grasterUmean[k][0], Ny*3, MPIU_REAL, MPI_SUM, 0, mesh->MESH_COMM);
        }

        // write to file
        if(!rank)
        {
            FILE *raster = fopen("uMean.asc", "w");
            fprintf(raster, "ncols %ld\n", Ny);
            fprintf(raster, "nrows %ld\n", Nx);
            fprintf(raster, "xllcorner %lf\n", mesh->bounds.xmin);
            fprintf(raster, "yllcorner %lf\n", mesh->bounds.ymin);
            fprintf(raster, "cellsize %lf\n", dC);
            fprintf(raster, "nodata_value -32768\n");
            for(k=0; k<Nx; k++)
            {
                for(i=0; i<Ny; i++)
                {
                    fprintf(raster, "%lf ", nMag(grasterUmean[k][i]));
                }

                fprintf(raster, "\n");
            }

            fclose(raster);
        }

        DMDAVecRestoreArray(fda, Coor, &coor);
        DMDAVecRestoreArray(fda, lavgU, &umean);
        DMDAVecRestoreArray(fda, mesh->lCent, &cent);

        VecDestroy(&lavgU);

        // clear memory
        for(k=0; k<Nx; k++)
        {
            std::vector<Cpt2D>  ().swap(rasterMesh[k]);
            std::vector<Cmpnts> ().swap(lrasterUmean[k]);
            std::vector<Cmpnts> ().swap(grasterUmean[k]);
        }
    }
    else
    {
        char warning[256];
        sprintf(warning, "cannot write averages on raster file, averaging is disabled\n");
        warningInFunction("writeJSectionToRaster",  warning);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( j=0; j<my; j++)
    {
        for(i=0; i<mx; i++)
        {
            sec->vectorSec[j][i].x = 0;
            sec->vectorSec[j][i].y = 0;
            sec->vectorSec[j][i].z = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/kSurfaces/" + std::to_string(kplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
        char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("kSectionLoadVector",  error);
    }

    for(j=0; j<my; j++)
    {
        PetscInt err;
        err = fread(&(sec->vectorSec[j][0]), sizeof(Cmpnts), mx, fp);
    }
    fclose(fp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( k=0; k<mz; k++)
    {
        for(i=0; i<mx; i++)
        {
          sec->vectorSec[k][i].x = 0;
          sec->vectorSec[k][i].y = 0;
          sec->vectorSec[k][i].z = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/jSurfaces/" + std::to_string(jplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
       char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("jSectionLoadVector",  error);
    }

    for(k=0; k<mz; k++)
    {
        PetscInt err;
        err = fread(&(sec->vectorSec[k][0]), sizeof(Cmpnts), mx, fp);
    }

    fclose(fp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( k=0; k<mz; k++)
    {
        for(j=0; j<my; j++)
        {
          sec->vectorSec[k][j].x = 0;
          sec->vectorSec[k][j].y = 0;
          sec->vectorSec[k][j].z = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/iSurfaces/" + std::to_string(iplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
       char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("iSectionLoadVector",  error);
    }

    for(k=0; k<mz; k++)
    {
        PetscInt err;
        err = fread(&(sec->vectorSec[k][0]), sizeof(Cmpnts), my, fp);
    }

    fclose(fp);

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode userSectionLoadVectorFromField(Vec &V, mesh_ *mesh, uSections *uSection, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    ueqn_              *ueqn  = mesh->access->ueqn;
    DM                 da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;
    Vec                lV;

    Cmpnts             ***v, ***cent;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // create the local field vector
    VecDuplicate(ueqn->lUcat, &lV);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // scatter from the global field to local for the trilinear interpolation
    DMGlobalToLocalBegin(fda, V, INSERT_VALUES, lV);
    DMGlobalToLocalEnd(fda, V, INSERT_VALUES, lV);

    PetscInt nx = uSection->nx, ny = uSection->ny;

    // set to zero
    for(j=0; j<ny; j++)
    {
        for(i=0; i<nx; i++)
        {
            uSection->vectorSec[j][i].x = 0;
            uSection->vectorSec[j][i].y = 0;
            uSection->vectorSec[j][i].z = 0;
        }
    }

    DMDAVecGetArray(fda, lV, &v);

    PetscInt ci, cj, ck;

    for(j=0; j<ny; j++)
    {
        for(i=0; i<nx; i++)
        {
            ci = uSection->closestId[j][i].i;
            cj = uSection->closestId[j][i].j;
            ck = uSection->closestId[j][i].k;

            Cmpnts cr = nSet(uSection->coor[j][i]);

            if(uSection->hasCoor[j][i])
            {
                vectorPointLocalVolumeInterpolation
                (
                        mesh,
                        cr.x, cr.y, cr.z,
                        ci, cj, ck,
                        cent,
                        v,
                        uSection->vectorSec[j][i]
                );
            }

        }
    }

    DMDAVecRestoreArray(fda, lV, &v);

    // reduce the values by storing only on the master node
    if(!rank)
    {
        for(j=0; j<ny; j++)
        {
            MPI_Reduce(MPI_IN_PLACE, uSection->vectorSec[j], nx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(j=0; j<ny; j++)
        {
            MPI_Reduce(uSection->vectorSec[j], uSection->vectorSec[j], nx*3, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    // apply periodic boundary conditions as average fields are only defined internally
    if(!rank)
    {
        for (j=0; j<ny; j++)
        {
            mSet(uSection->vectorSec[j][0],  uSection->vectorSec[j][1]);
            mSet(uSection->vectorSec[j][nx-1], uSection->vectorSec[j][nx-2]);
        }
        for (i=0; i<nx; i++)
        {
            mSet(uSection->vectorSec[0][i],  uSection->vectorSec[1][i]);
            mSet(uSection->vectorSec[ny-1][i], uSection->vectorSec[ny-2][i]);
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    VecDestroy(&lV);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts             ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (j=0; j<my; j++)
        {
            mSet(sec->vectorSec[j][0],    sec->vectorSec[j][1]);
            mSet(sec->vectorSec[j][mx-1], sec->vectorSec[j][mx-2]);
        }
        for (i=0; i<mx; i++)
        {
            mSet(sec->vectorSec[0][i],    sec->vectorSec[1][i]);
            mSet(sec->vectorSec[my-1][i], sec->vectorSec[my-2][i]);
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts             ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
    for( k=0; k<mz; k++)
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            mSet(sec->vectorSec[k][0],    sec->vectorSec[k][1]);
            mSet(sec->vectorSec[k][mx-1], sec->vectorSec[k][mx-2]);
        }
        for (i=0; i<mx; i++)
        {
            mSet(sec->vectorSec[0][i],    sec->vectorSec[1][i]);
            mSet(sec->vectorSec[mz-1][i], sec->vectorSec[mz-2][i]);
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    Cmpnts             ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            mSet(sec->vectorSec[k][0],    sec->vectorSec[k][1]);
            mSet(sec->vectorSec[k][my-1], sec->vectorSec[k][my-2]);
        }
        for (j=0; j<my; j++)
        {
            mSet(sec->vectorSec[0][j],    sec->vectorSec[1][j]);
            mSet(sec->vectorSec[mz-1][j], sec->vectorSec[mz-2][j]);
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionLoadSymmTensorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    symmTensor          ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
    for(j=0; j<my; j++)
    {
        for(i=0; i<mx; i++)
        {
            sec->symmTensorSec[j][i].xx = 0;
            sec->symmTensorSec[j][i].yy = 0;
            sec->symmTensorSec[j][i].zz = 0;
            sec->symmTensorSec[j][i].xy = 0;
            sec->symmTensorSec[j][i].xz = 0;
            sec->symmTensorSec[j][i].yz = 0;
        }
    }

    DMDAVecGetArray(mesh->sda, V, &v);

    if(kplane>=zs && kplane<ze)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                sec->symmTensorSec[j][i].xx = v[kplane][j][i].xx;
                sec->symmTensorSec[j][i].yy = v[kplane][j][i].yy;
                sec->symmTensorSec[j][i].zz = v[kplane][j][i].zz;
                sec->symmTensorSec[j][i].xy = v[kplane][j][i].xy;
                sec->symmTensorSec[j][i].xz = v[kplane][j][i].xz;
                sec->symmTensorSec[j][i].yz = v[kplane][j][i].yz;
            }
        }
    }

    DMDAVecRestoreArray(mesh->sda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->symmTensorSec[j], mx*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(j=0; j<my; j++)
        {
            MPI_Reduce(sec->symmTensorSec[j], sec->symmTensorSec[j], mx*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    // apply zeroGradient boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (j=0; j<my; j++)
        {
            sec->symmTensorSec[j][0].xx = sec->symmTensorSec[j][1].xx;
            sec->symmTensorSec[j][0].yy = sec->symmTensorSec[j][1].yy;
            sec->symmTensorSec[j][0].zz = sec->symmTensorSec[j][1].zz;
            sec->symmTensorSec[j][0].xy = sec->symmTensorSec[j][1].xy;
            sec->symmTensorSec[j][0].xz = sec->symmTensorSec[j][1].xz;
            sec->symmTensorSec[j][0].yz = sec->symmTensorSec[j][1].yz;

            sec->symmTensorSec[j][mx-1].xx = sec->symmTensorSec[j][mx-2].xx;
            sec->symmTensorSec[j][mx-1].yy = sec->symmTensorSec[j][mx-2].yy;
            sec->symmTensorSec[j][mx-1].zz = sec->symmTensorSec[j][mx-2].zz;
            sec->symmTensorSec[j][mx-1].xy = sec->symmTensorSec[j][mx-2].xy;
            sec->symmTensorSec[j][mx-1].xz = sec->symmTensorSec[j][mx-2].xz;
            sec->symmTensorSec[j][mx-1].yz = sec->symmTensorSec[j][mx-2].yz;
        }
        for (i=0; i<mx; i++)
        {
            sec->symmTensorSec[0][i].xx = sec->symmTensorSec[1][i].xx;
            sec->symmTensorSec[0][i].yy = sec->symmTensorSec[1][i].yy;
            sec->symmTensorSec[0][i].zz = sec->symmTensorSec[1][i].zz;
            sec->symmTensorSec[0][i].xy = sec->symmTensorSec[1][i].xy;
            sec->symmTensorSec[0][i].xz = sec->symmTensorSec[1][i].xz;
            sec->symmTensorSec[0][i].yz = sec->symmTensorSec[1][i].yz;

            sec->symmTensorSec[my-1][i].xx = sec->symmTensorSec[my-2][i].xx;
            sec->symmTensorSec[my-1][i].yy = sec->symmTensorSec[my-2][i].yy;
            sec->symmTensorSec[my-1][i].zz = sec->symmTensorSec[my-2][i].zz;
            sec->symmTensorSec[my-1][i].xy = sec->symmTensorSec[my-2][i].xy;
            sec->symmTensorSec[my-1][i].xz = sec->symmTensorSec[my-2][i].xz;
            sec->symmTensorSec[my-1][i].yz = sec->symmTensorSec[my-2][i].yz;
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionLoadSymmTensorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    symmTensor         ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
    for( k=0; k<mz; k++)
    {
        for(i=0; i<mx; i++)
        {
            sec->symmTensorSec[k][i].xx = 0;
            sec->symmTensorSec[k][i].yy = 0;
            sec->symmTensorSec[k][i].zz = 0;
            sec->symmTensorSec[k][i].xy = 0;
            sec->symmTensorSec[k][i].xz = 0;
            sec->symmTensorSec[k][i].yz = 0;
        }
    }

    DMDAVecGetArray(mesh->sda, V, &v);

    if(jplane>=ys && jplane<ye)
    {
        for (k=zs; k<ze; k++)
        {
            for (i=xs; i<xe; i++)
            {
                sec->symmTensorSec[k][i].xx = v[k][jplane][i].xx;
                sec->symmTensorSec[k][i].yy = v[k][jplane][i].yy;
                sec->symmTensorSec[k][i].zz = v[k][jplane][i].zz;
                sec->symmTensorSec[k][i].xy = v[k][jplane][i].xy;
                sec->symmTensorSec[k][i].xz = v[k][jplane][i].xz;
                sec->symmTensorSec[k][i].yz = v[k][jplane][i].yz;
            }
        }
    }

    DMDAVecRestoreArray(mesh->sda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->symmTensorSec[k], mx*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->symmTensorSec[k], sec->symmTensorSec[k], mx*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            sec->symmTensorSec[k][0].xx = sec->symmTensorSec[k][1].xx;
            sec->symmTensorSec[k][0].yy = sec->symmTensorSec[k][1].yy;
            sec->symmTensorSec[k][0].zz = sec->symmTensorSec[k][1].zz;
            sec->symmTensorSec[k][0].xy = sec->symmTensorSec[k][1].xy;
            sec->symmTensorSec[k][0].xz = sec->symmTensorSec[k][1].xz;
            sec->symmTensorSec[k][0].yz = sec->symmTensorSec[k][1].yz;

            sec->symmTensorSec[k][mx-1].xx = sec->symmTensorSec[k][mx-2].xx;
            sec->symmTensorSec[k][mx-1].yy = sec->symmTensorSec[k][mx-2].yy;
            sec->symmTensorSec[k][mx-1].zz = sec->symmTensorSec[k][mx-2].zz;
            sec->symmTensorSec[k][mx-1].xy = sec->symmTensorSec[k][mx-2].xy;
            sec->symmTensorSec[k][mx-1].xz = sec->symmTensorSec[k][mx-2].xz;
            sec->symmTensorSec[k][mx-1].yz = sec->symmTensorSec[k][mx-2].yz;
        }
        for (i=0; i<mx; i++)
        {
            sec->symmTensorSec[0][i].xx = sec->symmTensorSec[1][i].xx;
            sec->symmTensorSec[0][i].yy = sec->symmTensorSec[1][i].yy;
            sec->symmTensorSec[0][i].zz = sec->symmTensorSec[1][i].zz;
            sec->symmTensorSec[0][i].xy = sec->symmTensorSec[1][i].xy;
            sec->symmTensorSec[0][i].xz = sec->symmTensorSec[1][i].xz;
            sec->symmTensorSec[0][i].yz = sec->symmTensorSec[1][i].yz;

            sec->symmTensorSec[mz-1][i].xx = sec->symmTensorSec[mz-2][i].xx;
            sec->symmTensorSec[mz-1][i].yy = sec->symmTensorSec[mz-2][i].yy;
            sec->symmTensorSec[mz-1][i].zz = sec->symmTensorSec[mz-2][i].zz;
            sec->symmTensorSec[mz-1][i].xy = sec->symmTensorSec[mz-2][i].xy;
            sec->symmTensorSec[mz-1][i].xz = sec->symmTensorSec[mz-2][i].xz;
            sec->symmTensorSec[mz-1][i].yz = sec->symmTensorSec[mz-2][i].yz;
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionLoadSymmTensorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    symmTensor         ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
    for(k=0; k<mz; k++)
    {
        for(j=0; j<my; j++)
        {
            sec->symmTensorSec[k][j].xx = 0;
            sec->symmTensorSec[k][j].yy = 0;
            sec->symmTensorSec[k][j].zz = 0;
            sec->symmTensorSec[k][j].xy = 0;
            sec->symmTensorSec[k][j].xz = 0;
            sec->symmTensorSec[k][j].yz = 0;
        }
    }

    DMDAVecGetArray(mesh->sda, V, &v);

    if(iplane>=xs && iplane<xe)
    {
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                sec->symmTensorSec[k][j].xx = v[k][j][iplane].xx;
                sec->symmTensorSec[k][j].yy = v[k][j][iplane].yy;
                sec->symmTensorSec[k][j].zz = v[k][j][iplane].zz;
                sec->symmTensorSec[k][j].xy = v[k][j][iplane].xy;
                sec->symmTensorSec[k][j].xz = v[k][j][iplane].xz;
                sec->symmTensorSec[k][j].yz = v[k][j][iplane].yz;
            }
        }
    }

    DMDAVecRestoreArray(mesh->sda, V, &v);

    // reduce the values by storing only on the master node
    // mpi operation is sum because the other values are zero
    if(!rank)
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(MPI_IN_PLACE, sec->symmTensorSec[k], my*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(k=0; k<mz; k++)
        {
            MPI_Reduce(sec->symmTensorSec[k], sec->symmTensorSec[k], my*6, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            sec->symmTensorSec[k][0].xx = sec->symmTensorSec[k][1].xx;
            sec->symmTensorSec[k][0].yy = sec->symmTensorSec[k][1].yy;
            sec->symmTensorSec[k][0].zz = sec->symmTensorSec[k][1].zz;
            sec->symmTensorSec[k][0].xy = sec->symmTensorSec[k][1].xy;
            sec->symmTensorSec[k][0].xz = sec->symmTensorSec[k][1].xz;
            sec->symmTensorSec[k][0].yz = sec->symmTensorSec[k][1].yz;

            sec->symmTensorSec[k][my-1].xx = sec->symmTensorSec[k][my-2].xx;
            sec->symmTensorSec[k][my-1].yy = sec->symmTensorSec[k][my-2].yy;
            sec->symmTensorSec[k][my-1].zz = sec->symmTensorSec[k][my-2].zz;
            sec->symmTensorSec[k][my-1].xy = sec->symmTensorSec[k][my-2].xy;
            sec->symmTensorSec[k][my-1].xz = sec->symmTensorSec[k][my-2].xz;
            sec->symmTensorSec[k][my-1].yz = sec->symmTensorSec[k][my-2].yz;
        }
        for (j=0; j<my; j++)
        {
            sec->symmTensorSec[0][i].xx = sec->symmTensorSec[1][i].xx;
            sec->symmTensorSec[0][i].yy = sec->symmTensorSec[1][i].yy;
            sec->symmTensorSec[0][i].zz = sec->symmTensorSec[1][i].zz;
            sec->symmTensorSec[0][i].xy = sec->symmTensorSec[1][i].xy;
            sec->symmTensorSec[0][i].xz = sec->symmTensorSec[1][i].xz;
            sec->symmTensorSec[0][i].yz = sec->symmTensorSec[1][i].yz;

            sec->symmTensorSec[mz-1][j].xx = sec->symmTensorSec[mz-2][j].xx;
            sec->symmTensorSec[mz-1][j].yy = sec->symmTensorSec[mz-2][j].yy;
            sec->symmTensorSec[mz-1][j].zz = sec->symmTensorSec[mz-2][j].zz;
            sec->symmTensorSec[mz-1][j].xy = sec->symmTensorSec[mz-2][j].xy;
            sec->symmTensorSec[mz-1][j].xz = sec->symmTensorSec[mz-2][j].xz;
            sec->symmTensorSec[mz-1][j].yz = sec->symmTensorSec[mz-2][j].yz;
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( j=0; j<my; j++)
    {
        for(i=0; i<mx; i++)
        {
            sec->scalarSec[j][i] = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/kSurfaces/" + std::to_string(kplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
       char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("kSectionLoadScalar",  error);
    }

    for(j=0; j<my; j++)
    {
        PetscInt err;
        err = fread(&(sec->scalarSec[j][0]), sizeof(PetscReal), mx, fp);
    }
    fclose(fp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( k=0; k<mz; k++)
    {
        for(i=0; i<mx; i++)
        {
            sec->scalarSec[k][i] = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/jSurfaces/" + std::to_string(jplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
       char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("jSectionLoadScalar",  error);
    }

    for(k=0; k<mz; k++)
    {
        PetscInt err;
        err = fread(&(sec->scalarSec[k][0]), sizeof(PetscReal), mx, fp);
    }
    fclose(fp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // set to zero
    for( k=0; k<mz; k++)
    {
        for(j=0; j<my; j++)
        {
            sec->scalarSec[k][j] = 0;
        }
    }

    // get the file name to read
    word fname;
    fname = "./postProcessing/" + mesh->meshName + "/iSurfaces/" + std::to_string(iplane) + "/" + fieldName + "/" + stream.str();

    // open the file and read
    FILE *fp = fopen(fname.c_str(), "rb");

    if(!fp)
    {
       char error[512];
        sprintf(error, "cannot open file: %s\n", fname.c_str());
        fatalErrorInFunction("iSectionLoadScalar",  error);
    }

    for(k=0; k<mz; k++)
    {
        PetscInt err;
        err = fread(&(sec->scalarSec[k][0]), sizeof(PetscReal), my, fp);
    }
    fclose(fp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode userSectionLoadScalarFromField(Vec &V, mesh_ *mesh, uSections *uSection, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    peqn_              *peqn  = mesh->access->peqn;
    DMDALocalInfo      info = mesh->info;
    DM                 da   = mesh->da, fda = mesh->fda;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;
    Vec                lV;

    PetscReal          ***v;
    Cmpnts             ***cent;
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // create the local field vector
    VecDuplicate(peqn->lP, &lV);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // scatter from the global field to local for the trilinear interpolation
    DMGlobalToLocalBegin(da, V, INSERT_VALUES, lV);
    DMGlobalToLocalEnd(da, V, INSERT_VALUES, lV);

    PetscInt nx = uSection->nx, ny = uSection->ny;

    // set to zero
    for(j=0; j<ny; j++)
    {
        for(i=0; i<nx; i++)
        {
            uSection->scalarSec[j][i] = 0;
        }
    }

    DMDAVecGetArray(da, lV, &v);

    PetscInt ci, cj, ck;

    for(j=0; j<ny; j++)
    {
        for(i=0; i<nx; i++)
        {
            ci = uSection->closestId[j][i].i;
            cj = uSection->closestId[j][i].j;
            ck = uSection->closestId[j][i].k;

            Cmpnts cr = nSet(uSection->coor[j][i]);

            if(uSection->hasCoor[j][i])
            {
                scalarPointLocalVolumeInterpolation
                (
                        mesh,
                        cr.x, cr.y, cr.z,
                        ci, cj, ck,
                        cent,
                        v,
                        uSection->scalarSec[j][i]
                );
            }

        }
    }

    DMDAVecRestoreArray(da, lV, &v);

    // reduce the values by storing only on the master node
    if(!rank)
    {
        for(j=0; j<ny; j++)
        {
            MPI_Reduce(MPI_IN_PLACE, uSection->scalarSec[j], nx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }
    else
    {
        for(j=0; j<ny; j++)
        {
            MPI_Reduce(uSection->scalarSec[j], uSection->scalarSec[j], nx, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        }
    }

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (j=0; j<ny; j++)
        {
            uSection->scalarSec[j][0]    = uSection->scalarSec[j][1];
            uSection->scalarSec[j][nx-1] = uSection->scalarSec[j][nx-2];
        }
        for (i=0; i<mx; i++)
        {
            uSection->scalarSec[0][i]     = uSection->scalarSec[1][i];
            uSection->scalarSec[ny-1][i]  = uSection->scalarSec[ny-2][i];
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode kSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal          ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (j=0; j<my; j++)
        {
            sec->scalarSec[j][0]    = sec->scalarSec[j][1];
            sec->scalarSec[j][mx-1] = sec->scalarSec[j][mx-2];
        }
        for (i=0; i<mx; i++)
        {
            sec->scalarSec[0][i]     = sec->scalarSec[1][i];
            sec->scalarSec[my-1][i]  = sec->scalarSec[my-2][i];
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode jSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal          ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
    for( k=0; k<mz; k++)
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            sec->scalarSec[k][0]     = sec->scalarSec[k][1];
            sec->scalarSec[k][mx-1]  = sec->scalarSec[k][mx-2];
        }
        for (i=0; i<mx; i++)
        {
            sec->scalarSec[0][i]     = sec->scalarSec[1][i];
            sec->scalarSec[mz-1][i]  = sec->scalarSec[mz-2][i];
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode iSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time)
{
    clock_             *clock = mesh->access->clock;
    DMDALocalInfo      info = mesh->info;
    PetscInt           xs = info.xs, xe = info.xs + info.xm;
    PetscInt           ys = info.ys, ye = info.ys + info.ym;
    PetscInt           zs = info.zs, ze = info.zs + info.zm;
    PetscInt           mx = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscMPIInt        rank;

    PetscReal          ***v;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::stringstream stream;
    stream << std::fixed << std::setprecision(clock->timePrecision) << time;

    // get the file name to read
    word fname;
    fname = "./fields/" + mesh->meshName + "/" + stream.str() + "/" + fieldName;

    // read field
    PetscPrintf(mesh->MESH_COMM, " > %s\n", fieldName.c_str());

    PetscViewer  viewer;
    PetscViewerBinaryOpen(mesh->MESH_COMM, fname.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(V,viewer);
    PetscViewerDestroy(&viewer);

    // set to zero
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
    if(!rank)
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

    // apply periodic boundary conditions are average fields are only defined internally
    if(!rank)
    {
        for (k=0; k<mz; k++)
        {
            sec->scalarSec[k][0]    = sec->scalarSec[k][1];
            sec->scalarSec[k][my-1] = sec->scalarSec[k][my-2];
        }
        for (j=0; j<my; j++)
        {
            sec->scalarSec[0][j]    = sec->scalarSec[1][j];
            sec->scalarSec[mz-1][j] = sec->scalarSec[mz-2][j];
        }
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getIndexList(const char* dataLoc, std::vector<PetscInt> &indexSeries, PetscInt &ntimes)
{
    // get file names inside path
    DIR *dir; struct dirent *diread;

    // pointer for stdtod and stdtol
    char *eptr;

    ntimes = 0;

    if ((dir = opendir(dataLoc)) != nullptr)
    {
        while ((diread = readdir(dir)) != nullptr)
        {
            char* indexName = diread->d_name;
            if
            (
                strcmp(indexName, ".") !=0 &&
                strcmp(indexName, "..") !=0
            )
            {
                PetscInt indexValue;
                indexValue = std::strtol(indexName, &eptr, 10);

                // make sure the folder's name is a number by comparing char value after the name:
                // it should be the null character
                if(*eptr == '\0')
                {
                    indexSeries.push_back(indexValue);
                    ntimes++;
                }
            }
        }
        closedir (dir);
    }
    else
    {
       char error[512];
        sprintf(error, "could not access %s directory\n", dataLoc);
        fatalErrorInFunction("getIndexList", error);
    }

    // sort the timeSeries
    for(PetscInt i=0; i<ntimes; i++)
    {
        PetscInt min   = 1e10;
        PetscInt value = 0;
        PetscInt label    = 0;

        for(PetscInt s=i; s<ntimes; s++)
        {
            if(indexSeries[s] < min)
            {
                value = indexSeries[s];
                label = s;
                min   = value;
            }
        }
        // exchange values so that elements are not lost
        indexSeries[label] = indexSeries[i];
        // put the min value on the unchanged part at the last index of changed part
        indexSeries[i] = value;
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryKSectionsToXMF(domain_ *domain)
{
    PetscInt    nDomains = domain[0].info.nDomains;
    word        meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_      flags    = domain[d].flags;

        if(flags.isAquisitionActive)
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;
            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
                fatalErrorInFunction("binaryKSectionsToXMF",  error);
            }

            if(acquisition->kSections->available)
            {
              PetscPrintf(mesh->MESH_COMM, "Processing k-sections for mesh: %s...", mesh->meshName.c_str());
              sections *kSections = acquisition->kSections;

              // create kSections folder
              sectionDir = "./XMF/" + mesh->meshName + "/kSections";
              createDir(mesh->MESH_COMM, sectionDir.c_str());

              for(PetscInt k=0; k<kSections->nSections; k++)
              {
                PetscInt kplane = kSections->indices[k];

                // exclude ghost nodes
                if (kplane<1 || kplane>mz-2) continue;

                // get list of available times
                word path2times = "./postProcessing/" + mesh->meshName + "/kSurfaces/" + std::to_string(kplane) + "/U" ;

                // see if can access the directory (shouldn't give error
                // as samplig/kSection could have been created just for averages)
                DIR* dir;
                if ((dir = opendir(path2times.c_str())) == nullptr)
                {
                    PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
                    return(0);
                }

                std::vector<PetscReal>      timeSeries;
                PetscInt                      ntimes;
                getTimeList(path2times.c_str(), timeSeries, ntimes);

                if(!rank)
                {
                  //create the k index folder
                  indexDir = sectionDir + "/" + std::to_string(kplane);
                  createDir(mesh->MESH_COMM, indexDir.c_str());

                  // create XMF file
                  fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_kSec" + std::to_string(kplane) + ".xmf";

                  // write XMF file intro
                  xmf = fopen(fieldsFileName.c_str(), "w");
                  fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                  fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                  fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                  fprintf(xmf, " <Domain>\n");
                  fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                  fclose(xmf);
                }

                // loop over times
                for(PetscInt ti=0; ti<ntimes; ti++)
                {
                  // HDF5 file with path
                  word fileName;

                  // HDF5 file w/o path
                  word hdfileName;

                  std::stringstream stream;
                  stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

                  fileName   = indexDir + "/" + thisCaseName() + "_" + "kSec" + std::to_string(kplane) + "_" + stream.str();
                  hdfileName = thisCaseName() + "_" + "kSec" + std::to_string(kplane) + "_" + stream.str();

                  // open this time section in the XMF file
                  if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, my-1, 1,"3DSMesh", timeSeries[ti]);

                  // Write the data file.
                  hid_t     dataspace_id;
                  hsize_t   dims[3];
                  herr_t    status;

                  // write HDF file
                  hid_t	file_id;
                  if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                  // ************************** write mesh points **************************
                  dims[0] = 1;
                  dims[1]	= my - 1;
                  dims[2]	= mx - 1;

                  if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                  writeKSectionPointsToXMF
                  (
                      mesh,
                      fieldsFileName.c_str(),
                      hdfileName.c_str(),
                      &file_id,
                      &dataspace_id,
                      timeSeries[ti],
                      kplane
                  );

                  if(!rank) status = H5Sclose(dataspace_id);
                  if(!rank) status = H5Fclose(file_id);

                  // load velocity
                  kSectionLoadVector(mesh, kSections, kplane, "U", timeSeries[ti]);

                  if(!rank)
                  {
                    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeKSectionVectorToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        "U",
                        kSections->vectorSec
                    );

                    status = H5Sclose(dataspace_id);
                    status = H5Fclose(file_id);
                  }

                  // load pressure
                  kSectionLoadScalar(mesh, kSections, kplane, "p", timeSeries[ti]);

                  if(!rank)
                  {
                    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeKSectionScalarToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        "p",
                        kSections->scalarSec
                    );

                    status = H5Sclose(dataspace_id);
                    status = H5Fclose(file_id);
                  }

                  // load nut
                  if(flags.isLesActive)
                  {
                    kSectionLoadScalar(mesh, kSections, kplane, "nut", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeKSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "nut",
                          kSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                  }

                  // load temperature
                  if(flags.isTeqnActive)
                  {
                    kSectionLoadScalar(mesh, kSections, kplane, "T", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeKSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "T",
                          kSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                  }

                  // load temperature
                  if(flags.isIBMActive)
                  {
                    kSectionLoadScalar(mesh, kSections, kplane, "nv", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeKSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "nv",
                          kSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                  }

                  // close this time section in the XMF file
                  if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                  // wait all processes
                  MPI_Barrier(mesh->MESH_COMM);
                }

                if(!rank)
                {
                    // write XMF file end
                    xmf = fopen(fieldsFileName.c_str(), "a");
                    fprintf(xmf, "   </Grid>\n");
                    fprintf(xmf, " </Domain>\n");
                    fprintf(xmf, "</Xdmf>\n");
                    fclose(xmf);
                }

              }

              PetscPrintf(mesh->MESH_COMM, "done\n\n");

            }
        }
    }

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode fieldUserDefinedPlaneToXMF(domain_ *domain)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, surfaceDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;
        io_       *io      = domain[d].io;

        if
        (
            flags.isAquisitionActive &&
            (io->averaging || io->phaseAveraging)
        )
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;

            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt      mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
              char error[512];
              sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
              fatalErrorInFunction("fieldUserDefinedPlaneToXMF",  error);
            }

            if(acquisition->userSections->available)
            {
                PetscPrintf(mesh->MESH_COMM, "On-the-fly user defined sections extraction for mesh: %s...\n", mesh->meshName.c_str());

                for(PetscInt s=0; s<acquisition->userSections->nSections; s++)
                {
                    uSections *uSection = acquisition->userSections->uSection[s];

                    // create userdefined section folder
                    sectionDir = "./XMF/" + mesh->meshName + "/userSections";
                    createDirNoRemove(mesh->MESH_COMM, sectionDir.c_str());

                    // get list of available times
                    word path2times = "./fields/" + mesh->meshName;

                    std::vector<PetscReal>        timeSeries;
                    PetscInt                      ntimes;
                    getTimeList(path2times.c_str(), timeSeries, ntimes);

                    if(!rank)
                    {
                        //create the userSection folder
                        surfaceDir = sectionDir + "/" + uSection->sectionName + "_section";
                        createDirNoRemove(mesh->MESH_COMM, surfaceDir.c_str());

                        // create XMF file
                        fieldsFileName = surfaceDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_" + uSection->sectionName + "_averages.xmf";

                        // write XMF file intro
                        xmf = fopen(fieldsFileName.c_str(), "w");
                        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                        fprintf(xmf, " <Domain>\n");
                        fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                        fclose(xmf);
                    }

                    PetscPrintf(mesh->MESH_COMM, "\nSection %s:  \n", uSection->sectionName.c_str());

                    // HDF5 file with path
                    word fileName;

                    // HDF5 file w/o path
                    word hdfileName;

                    // set time precision
                    std::stringstream stream;
                    stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ntimes-1];

                    fileName   = surfaceDir + "/" + thisCaseName() + "_" +  uSection->sectionName  + "_averages_" + stream.str();
                    hdfileName = thisCaseName() + "_" + uSection->sectionName  + "_averages_" + stream.str();
                    
                    // open this time section in the XMF file
                    if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), uSection->nx, uSection->ny, 1, "3DSMesh", timeSeries[ntimes-1]);

                    // Write the data file.
                    hid_t     dataspace_id;
                    hsize_t   dims[3];
                    herr_t    status;

                    // write HDF file
                    hid_t	file_id;
                    if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                    dims[0] = 1;
                    dims[1]	= uSection->ny;
                    dims[2]	= uSection->nx;

                    if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeUserSectionPointsToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ntimes-1],
                        uSection
                    );

                    if(!rank) status = H5Sclose(dataspace_id);
                    if(!rank) status = H5Fclose(file_id);

                    if(io->averaging)
                    {
                        // load average velocity
                        userSectionLoadVectorFromField(acquisition->fields->avgU, mesh, uSection, "avgU", timeSeries[ntimes-1]);

                        if(!rank)
                        {
                            file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                            dataspace_id = H5Screate_simple(3, dims, NULL);

                            writeUserSectionVectorToXMF
                            (
                                mesh,
                                fieldsFileName.c_str(),
                                hdfileName.c_str(),
                                &file_id,
                                &dataspace_id,
                                timeSeries[ntimes-1],
                                "avgU",
                                uSection
                            );

                            status = H5Sclose(dataspace_id);
                            status = H5Fclose(file_id);
                        }

                        // load average pressure
                        userSectionLoadScalarFromField(acquisition->fields->avgP, mesh, uSection, "avgP", timeSeries[ntimes-1]);

                        if(!rank)
                        {
                            file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                            dataspace_id = H5Screate_simple(3, dims, NULL);

                            writeUserSectionScalarToXMF
                            (
                                mesh,
                                fieldsFileName.c_str(),
                                hdfileName.c_str(),
                                &file_id,
                                &dataspace_id,
                                timeSeries[ntimes-1],
                                "avgP",
                                uSection
                            );

                            status = H5Sclose(dataspace_id);
                            status = H5Fclose(file_id);
                        }

                        if(flags.isLesActive)
                        {
                            // load average nut
                            userSectionLoadScalarFromField(acquisition->fields->avgNut, mesh, uSection, "avgNut", timeSeries[ntimes-1]);

                            if(!rank)
                            {
                                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                                dataspace_id = H5Screate_simple(3, dims, NULL);

                                writeUserSectionScalarToXMF
                                (
                                    mesh,
                                    fieldsFileName.c_str(),
                                    hdfileName.c_str(),
                                    &file_id,
                                    &dataspace_id,
                                    timeSeries[ntimes-1],
                                    "avgNut",
                                    uSection
                                );

                                status = H5Sclose(dataspace_id);
                                status = H5Fclose(file_id);
                            }

                            // load average cs
                            userSectionLoadScalarFromField(acquisition->fields->avgCs, mesh, uSection, "avgCs", timeSeries[ntimes-1]);

                            if(!rank)
                            {
                                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                                dataspace_id = H5Screate_simple(3, dims, NULL);

                                writeUserSectionScalarToXMF
                                (
                                    mesh,
                                    fieldsFileName.c_str(),
                                    hdfileName.c_str(),
                                    &file_id,
                                    &dataspace_id,
                                    timeSeries[ntimes-1],
                                    "avgCs",
                                    uSection
                                );

                                status = H5Sclose(dataspace_id);
                                status = H5Fclose(file_id);
                            }

                        }

                        // load temperature
                        if(flags.isTeqnActive)
                        {
                            // no fields for now
                        }
                    }

                    if(io->phaseAveraging)
                    {
                        // load average velocity
                        userSectionLoadVectorFromField(acquisition->fields->pAvgU, mesh, uSection, "phAvgU", timeSeries[ntimes-1]);

                        if(!rank)
                        {
                            file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                            dataspace_id = H5Screate_simple(3, dims, NULL);

                            writeUserSectionVectorToXMF
                            (
                                mesh,
                                fieldsFileName.c_str(),
                                hdfileName.c_str(),
                                &file_id,
                                &dataspace_id,
                                timeSeries[ntimes-1],
                                "phAvgU",
                                uSection
                            );

                            status = H5Sclose(dataspace_id);
                            status = H5Fclose(file_id);
                        }

                        // load average pressure
                        userSectionLoadScalarFromField(acquisition->fields->pAvgP, mesh, uSection, "phAvgP", timeSeries[ntimes-1]);

                        if(!rank)
                        {
                            file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                            dataspace_id = H5Screate_simple(3, dims, NULL);

                            writeUserSectionScalarToXMF
                            (
                                mesh,
                                fieldsFileName.c_str(),
                                hdfileName.c_str(),
                                &file_id,
                                &dataspace_id,
                                timeSeries[ntimes-1],
                                "phAvgP",
                                uSection
                            );

                            status = H5Sclose(dataspace_id);
                            status = H5Fclose(file_id);
                        }

                        if(flags.isLesActive)
                        {
                            // load average nut
                            userSectionLoadScalarFromField(acquisition->fields->pAvgNut, mesh, uSection, "phAvgNut", timeSeries[ntimes-1]);

                            if(!rank)
                            {
                                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                                dataspace_id = H5Screate_simple(3, dims, NULL);

                                writeUserSectionScalarToXMF
                                (
                                    mesh,
                                    fieldsFileName.c_str(),
                                    hdfileName.c_str(),
                                    &file_id,
                                    &dataspace_id,
                                    timeSeries[ntimes-1],
                                    "phAvgNut",
                                    uSection
                                );

                                status = H5Sclose(dataspace_id);
                                status = H5Fclose(file_id);
                            }

                            // load average cs
                            userSectionLoadScalarFromField(acquisition->fields->pAvgCs, mesh, uSection, "phAvgCs", timeSeries[ntimes-1]);

                            if(!rank)
                            {
                                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                                dataspace_id = H5Screate_simple(3, dims, NULL);

                                writeUserSectionScalarToXMF
                                (
                                    mesh,
                                    fieldsFileName.c_str(),
                                    hdfileName.c_str(),
                                    &file_id,
                                    &dataspace_id,
                                    timeSeries[ntimes-1],
                                    "phAvgCs",
                                    uSection
                                );

                                status = H5Sclose(dataspace_id);
                                status = H5Fclose(file_id);
                            }

                        }

                        // load temperature
                        if(flags.isTeqnActive)
                        {
                            // no fields for now
                        }
                    }                    
                    // close this time section in the XMF file
                    if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                    // wait all processes
                    MPI_Barrier(mesh->MESH_COMM);

                    if(!rank)
                    {
                        // write XMF file end
                        xmf = fopen(fieldsFileName.c_str(), "a");
                        fprintf(xmf, "   </Grid>\n");
                        fprintf(xmf, " </Domain>\n");
                        fprintf(xmf, "</Xdmf>\n");
                        fclose(xmf);
                    }
                }
                
            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");         
        }        
    }

    return 0;
}

//***************************************************************************************************************//
PetscErrorCode fieldKSectionsToXMF(domain_ *domain)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;
        io_       *io      = domain[d].io;

        if
        (
            flags.isAquisitionActive &&
            (io->averaging || io->phaseAveraging)
        )
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;

            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt      mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
              char error[512];
              sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
              fatalErrorInFunction("fieldKSectionsToXMF",  error);
            }

            if(acquisition->kSections->available)
            {
                PetscPrintf(mesh->MESH_COMM, "On-the-fly k-sections extraction for mesh: %s...\n", mesh->meshName.c_str());
                sections *kSections = acquisition->kSections;

                // create jSections folder
                sectionDir = "./XMF/" + mesh->meshName + "/kSections";
                createDirNoRemove(mesh->MESH_COMM, sectionDir.c_str());

                for(PetscInt k=0; k<kSections->nSections; k++)
                {
                  PetscInt kplane = kSections->indices[k];

                  // exclude ghost nodes
                  if (kplane<1 || kplane>mz-2) continue;

                  // get list of available times
                  word path2times = "./fields/" + mesh->meshName;

                  std::vector<PetscReal>        timeSeries;
                  PetscInt                      ntimes;
                  getTimeList(path2times.c_str(), timeSeries, ntimes);

                  if(!rank)
                  {
                    //create the k index folder
                    indexDir = sectionDir + "/" + std::to_string(kplane);
                    createDirNoRemove(mesh->MESH_COMM, indexDir.c_str());

                    // create XMF file
                    fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_kSec" + std::to_string(kplane) + "_averages.xmf";

                    // write XMF file intro
                    xmf = fopen(fieldsFileName.c_str(), "w");
                    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                    fprintf(xmf, " <Domain>\n");
                    fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                    fclose(xmf);
                  }

                  // now we get only the last time (the meaningful one)

                  // HDF5 file with path
                  word fileName;

                  // HDF5 file w/o path
                  word hdfileName;

                  // set time precision
                  std::stringstream stream;
                  stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ntimes-1];

                  fileName   = indexDir + "/" + thisCaseName() + "_" + "kSec" + std::to_string(kplane) + "_averages_" + stream.str();
                  hdfileName = thisCaseName() + "_" + "kSec" + std::to_string(kplane) + "_averages_" + stream.str();

                  // open this time section in the XMF file
                  if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, my-1, 1, "3DSMesh", timeSeries[ntimes-1]);

                  // Write the data file.
                  hid_t     dataspace_id;
                  hsize_t   dims[3];
                  herr_t    status;

                  // write HDF file
                  hid_t	file_id;
                  if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                  // ************************** write mesh points **************************
                  dims[0] = 1;
                  dims[1]	= my - 1;
                  dims[2]	= mx - 1;

                  if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                  writeKSectionPointsToXMF
                  (
                      mesh,
                      fieldsFileName.c_str(),
                      hdfileName.c_str(),
                      &file_id,
                      &dataspace_id,
                      timeSeries[ntimes-1],
                      kplane
                  );

                  if(!rank) status = H5Sclose(dataspace_id);
                  if(!rank) status = H5Fclose(file_id);

                  if(io->averaging)
                  {
                      // load average velocity
                      kSectionLoadVectorFromField(acquisition->fields->avgU, mesh, kSections, kplane, "avgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgU",
                              kSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      kSectionLoadSymmTensorFromField(acquisition->fields->avgUU, mesh, kSections, kplane, "avgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgUU",
                              kSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      kSectionLoadScalarFromField(acquisition->fields->avgP, mesh, kSections, kplane, "avgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgP",
                              kSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          kSectionLoadScalarFromField(acquisition->fields->avgNut, mesh, kSections, kplane, "avgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeKSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgNut",
                                  kSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          kSectionLoadScalarFromField(acquisition->fields->avgCs, mesh, kSections, kplane, "avgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeKSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgCs",
                                  kSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  if(io->phaseAveraging)
                  {
                      // load average velocity
                      kSectionLoadVectorFromField(acquisition->fields->pAvgU, mesh, kSections, kplane, "phAvgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgU",
                              kSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      kSectionLoadSymmTensorFromField(acquisition->fields->pAvgUU, mesh, kSections, kplane, "phAvgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgUU",
                              kSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      kSectionLoadScalarFromField(acquisition->fields->pAvgP, mesh, kSections, kplane, "phAvgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeKSectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgP",
                              kSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          kSectionLoadScalarFromField(acquisition->fields->pAvgNut, mesh, kSections, kplane, "phAvgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeKSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgNut",
                                  kSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          kSectionLoadScalarFromField(acquisition->fields->pAvgCs, mesh, kSections, kplane, "phAvgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeKSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgCs",
                                  kSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  // close this time section in the XMF file
                  if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                  // wait all processes
                  MPI_Barrier(mesh->MESH_COMM);

                    if(!rank)
                    {
                        // write XMF file end
                        xmf = fopen(fieldsFileName.c_str(), "a");
                        fprintf(xmf, "   </Grid>\n");
                        fprintf(xmf, " </Domain>\n");
                        fprintf(xmf, "</Xdmf>\n");
                        fclose(xmf);
                    }
                }
            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryJSectionsToXMF(domain_ *domain, postProcess *pp)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;

        if(flags.isAquisitionActive)
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;
            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
                fatalErrorInFunction("binaryJSectionsToXMF",  error);
            }

            if(acquisition->jSections->available)
            {
                PetscPrintf(mesh->MESH_COMM, "Processing j-sections for mesh: %s...", mesh->meshName.c_str());
                sections *jSections = acquisition->jSections;

                // create jSections folder
                sectionDir = "./XMF/" + mesh->meshName + "/jSections";
                createDir(mesh->MESH_COMM, sectionDir.c_str());

                for(PetscInt j=0; j<jSections->nSections; j++)
                {
                  PetscInt jplane = jSections->indices[j];

                  // exclude ghost nodes
                  if (jplane<1 || jplane>my-2) continue;

                  // get list of available times
                  word path2times = "./postProcessing/" + mesh->meshName + "/jSurfaces/" + std::to_string(jplane) + "/U" ;

                  // see if can access the directory (shouldn't give error
                  // as samplig/jSection could have been created just for averages)
                  DIR* dir;
                  if ((dir = opendir(path2times.c_str())) == nullptr)
                  {
                      PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
                      return(0);
                  }

                  std::vector<PetscReal>      timeSeries;
                  PetscInt                      ntimes;
                  getTimeList(path2times.c_str(), timeSeries, ntimes);

                  if(!rank)
                  {
                    //create the j index folder
                    indexDir = sectionDir + "/" + std::to_string(jplane);
                    createDir(mesh->MESH_COMM, indexDir.c_str());

                    // create XMF file
                    fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_jSec" + std::to_string(jplane) + ".xmf";

                    // write XMF file intro
                    xmf = fopen(fieldsFileName.c_str(), "w");
                    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                    fprintf(xmf, " <Domain>\n");
                    fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                    fclose(xmf);
                  }

                  // loop over times
                  for(PetscInt ti=0; ti<ntimes; ti++)
                  {
                    // HDF5 file with path
                    word fileName;

                    // HDF5 file w/o path
                    word hdfileName;

                    // set time precision
                    std::stringstream stream;
                    stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

                    fileName   = indexDir + "/" + thisCaseName() + "_" + "jSec" + std::to_string(jplane) + "_" + stream.str();
                    hdfileName = thisCaseName() + "_" + "jSec" + std::to_string(jplane) + "_" + stream.str();

                    // open this time section in the XMF file
                    if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, 1, mz-1, "3DSMesh", timeSeries[ti]);

                    // Write the data file.
                    hid_t     dataspace_id;
                    hsize_t   dims[3];
                    herr_t    status;

                    // write HDF file
                    hid_t	file_id;
                    if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                    // ************************** write mesh points **************************
                    dims[0] = mz - 1;
                    dims[1]	= 1;
                    dims[2]	= mx - 1;

                    if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeJSectionPointsToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        jplane
                    );

                    if(!rank) status = H5Sclose(dataspace_id);
                    if(!rank) status = H5Fclose(file_id);

                    // load velocity
                    jSectionLoadVector(mesh, jSections, jplane, "U", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeJSectionVectorToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "U",
                          jSections->vectorSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                    // load pressure
                    jSectionLoadScalar(mesh, jSections, jplane, "p", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeJSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "p",
                          jSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                    // load nut
                    if(flags.isLesActive)
                    {
                      jSectionLoadScalar(mesh, jSections, jplane, "nut", timeSeries[ti]);

                      if(!rank)
                      {
                        file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                        dataspace_id = H5Screate_simple(3, dims, NULL);

                        writeJSectionScalarToXMF
                        (
                            mesh,
                            fieldsFileName.c_str(),
                            hdfileName.c_str(),
                            &file_id,
                            &dataspace_id,
                            timeSeries[ti],
                            "nut",
                            jSections->scalarSec
                        );

                        status = H5Sclose(dataspace_id);
                        status = H5Fclose(file_id);
                      }

                    }

                    // load temperature
                    if(flags.isTeqnActive)
                    {
                      jSectionLoadScalar(mesh, jSections, jplane, "T", timeSeries[ti]);

                      if(!rank)
                      {
                        file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                        dataspace_id = H5Screate_simple(3, dims, NULL);

                        writeJSectionScalarToXMF
                        (
                            mesh,
                            fieldsFileName.c_str(),
                            hdfileName.c_str(),
                            &file_id,
                            &dataspace_id,
                            timeSeries[ti],
                            "T",
                            jSections->scalarSec
                        );

                        status = H5Sclose(dataspace_id);
                        status = H5Fclose(file_id);
                      }

                    }

                    // load nut
                    if(flags.isIBMActive)
                    {
                      jSectionLoadScalar(mesh, jSections, jplane, "nv", timeSeries[ti]);

                      if(!rank)
                      {
                        file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                        dataspace_id = H5Screate_simple(3, dims, NULL);

                        writeJSectionScalarToXMF
                        (
                            mesh,
                            fieldsFileName.c_str(),
                            hdfileName.c_str(),
                            &file_id,
                            &dataspace_id,
                            timeSeries[ti],
                            "nv",
                            jSections->scalarSec
                        );

                        status = H5Sclose(dataspace_id);
                        status = H5Fclose(file_id);
                      }

                    }

                    // close this time section in the XMF file
                    if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                    // wait all processes
                    MPI_Barrier(mesh->MESH_COMM);
                  }

                  if(!rank)
                  {
                      // write XMF file end
                      xmf = fopen(fieldsFileName.c_str(), "a");
                      fprintf(xmf, "   </Grid>\n");
                      fprintf(xmf, " </Domain>\n");
                      fprintf(xmf, "</Xdmf>\n");
                      fclose(xmf);
                  }
                }

                PetscPrintf(mesh->MESH_COMM, "done\n\n");

                if(pp->writeRaster)
                {
                    for(PetscInt s=0; s<jSections->nSections; s++)
                    {
                        writeJSectionToRaster(mesh, jSections->indices[s]);
                    }
                }

            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode fieldJSectionsToXMF(domain_ *domain)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;
        io_       *io      = domain[d].io;

        if
        (
            flags.isAquisitionActive &&
            (io->averaging || io->phaseAveraging)
        )
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;

            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
              char error[512];
              sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
              fatalErrorInFunction("fieldJSectionsToXMF",  error);
            }

            if(acquisition->jSections->available)
            {
                PetscPrintf(mesh->MESH_COMM, "On-the-fly j-sections extraction for mesh: %s...\n", mesh->meshName.c_str());
                sections *jSections = acquisition->jSections;

                // create jSections folder
                sectionDir = "./XMF/" + mesh->meshName + "/jSections";
                createDirNoRemove(mesh->MESH_COMM, sectionDir.c_str());

                for(PetscInt j=0; j<jSections->nSections; j++)
                {
                  PetscInt jplane = jSections->indices[j];

                  // exclude ghost nodes
                  if (jplane<1 || jplane>my-2) continue;

                  // get list of available times
                  word path2times = "./fields/" + mesh->meshName;

                  std::vector<PetscReal>        timeSeries;
                  PetscInt                      ntimes;
                  getTimeList(path2times.c_str(), timeSeries, ntimes);

                  if(!rank)
                  {
                    //create the j index folder
                    indexDir = sectionDir + "/" + std::to_string(jplane);
                    createDirNoRemove(mesh->MESH_COMM, indexDir.c_str());

                    // create XMF file
                    fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_jSec" + std::to_string(jplane) + "_averages.xmf";

                    // write XMF file intro
                    xmf = fopen(fieldsFileName.c_str(), "w");
                    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                    fprintf(xmf, " <Domain>\n");
                    fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                    fclose(xmf);
                  }

                  // now we get only the last time (the meaningful one)

                  // HDF5 file with path
                  word fileName;

                  // HDF5 file w/o path
                  word hdfileName;

                  // set time precision
                  std::stringstream stream;
                  stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ntimes-1];

                  fileName   = indexDir + "/" + thisCaseName() + "_" + "jSec" + std::to_string(jplane) + "_averages_" + stream.str();
                  hdfileName = thisCaseName() + "_" + "jSec" + std::to_string(jplane) + "_averages_" + stream.str();

                  // open this time section in the XMF file
                  if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, 1, mz-1, "3DSMesh", timeSeries[ntimes-1]);

                  // Write the data file.
                  hid_t     dataspace_id;
                  hsize_t   dims[3];
                  herr_t    status;

                  // write HDF file
                  hid_t	file_id;
                  if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                  // ************************** write mesh points **************************
                  dims[0] = mz - 1;
                  dims[1] = 1;
                  dims[2] = mx - 1;

                  if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                  writeJSectionPointsToXMF
                  (
                      mesh,
                      fieldsFileName.c_str(),
                      hdfileName.c_str(),
                      &file_id,
                      &dataspace_id,
                      timeSeries[ntimes-1],
                      jplane
                  );

                  if(!rank) status = H5Sclose(dataspace_id);
                  if(!rank) status = H5Fclose(file_id);

                  if(io->averaging)
                  {
                      // load average velocity
                      jSectionLoadVectorFromField(acquisition->fields->avgU, mesh, jSections, jplane, "avgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgU",
                              jSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      jSectionLoadSymmTensorFromField(acquisition->fields->avgUU, mesh, jSections, jplane, "avgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgUU",
                              jSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      jSectionLoadScalarFromField(acquisition->fields->avgP, mesh, jSections, jplane, "avgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgP",
                              jSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          jSectionLoadScalarFromField(acquisition->fields->avgNut, mesh, jSections, jplane, "avgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeJSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgNut",
                                  jSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          jSectionLoadScalarFromField(acquisition->fields->avgCs, mesh, jSections, jplane, "avgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeJSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgCs",
                                  jSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  if(io->phaseAveraging)
                  {
                      // load average velocity
                      jSectionLoadVectorFromField(acquisition->fields->pAvgU, mesh, jSections, jplane, "phAvgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgU",
                              jSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      jSectionLoadSymmTensorFromField(acquisition->fields->pAvgUU, mesh, jSections, jplane, "phAvgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgUU",
                              jSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      jSectionLoadScalarFromField(acquisition->fields->pAvgP, mesh, jSections, jplane, "phAvgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeJSectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgP",
                              jSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          jSectionLoadScalarFromField(acquisition->fields->pAvgNut, mesh, jSections, jplane, "phAvgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeJSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgNut",
                                  jSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          jSectionLoadScalarFromField(acquisition->fields->pAvgCs, mesh, jSections, jplane, "phAvgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeJSectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgCs",
                                  jSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  // close this time section in the XMF file
                  if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                  // wait all processes
                  MPI_Barrier(mesh->MESH_COMM);

                    if(!rank)
                    {
                        // write XMF file end
                        xmf = fopen(fieldsFileName.c_str(), "a");
                        fprintf(xmf, "   </Grid>\n");
                        fprintf(xmf, " </Domain>\n");
                        fprintf(xmf, "</Xdmf>\n");
                        fclose(xmf);
                    }
                }

            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryISectionsToXMF(domain_ *domain)
{
  PetscInt  nDomains = domain[0].info.nDomains;
  word      meshDir, sectionDir, indexDir;

  for(PetscInt d=0; d<nDomains; d++)
  {
    flags_    flags    = domain[d].flags;

    if(flags.isAquisitionActive)
    {


      // get pointers
      clock_ *clock = domain[d].clock;
      mesh_  *mesh  = domain[d].mesh;

      acquisition_ *acquisition = domain[d].acquisition;
      ueqn_  *ueqn  = domain[d].ueqn;
      peqn_  *peqn  = domain[d].peqn;
      teqn_  *teqn;
      les_   *les;

      if(flags.isTeqnActive) teqn = domain[d].teqn;
      if(flags.isLesActive)  les  = domain[d].les;

      DMDALocalInfo info = mesh->info;
      PetscInt           mx = info.mx, my = info.my, mz = info.mz;

      PetscMPIInt           rank;
      MPI_Comm_rank(mesh->MESH_COMM, &rank);

      word fieldsFileName;
      FILE *xmf;

      meshDir     = "./XMF/" + mesh->meshName;

      // create domain directory within XMF folder
      errno = 0;
      PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
      if(dirRes != 0 && errno != EEXIST)
      {
         char error[512];
          sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
          fatalErrorInFunction("binaryISectionsToXMF",  error);
      }

      if(acquisition->iSections->available)
      {
        PetscPrintf(mesh->MESH_COMM, "Processing i-sections for mesh: %s...", mesh->meshName.c_str());
        sections *iSections = acquisition->iSections;

        // create iSections folder
        sectionDir = "./XMF/" + mesh->meshName + "/iSections";
        createDir(mesh->MESH_COMM, sectionDir.c_str());

        for(PetscInt i=0; i<iSections->nSections; i++)
        {
          PetscInt iplane = iSections->indices[i];

          // exclude ghost nodes
          if (iplane<1 || iplane>mx-2) continue;

          // get list of available times
          word path2times = "./postProcessing/" + mesh->meshName + "/iSurfaces/" + std::to_string(iplane) + "/U" ;

          // see if can access the directory (shouldn't give error
          // as samplig/iSection could have been created just for averages)
          DIR* dir;
          if ((dir = opendir(path2times.c_str())) == nullptr)
          {
              PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
              return(0);
          }

          std::vector<PetscReal>      timeSeries;
          PetscInt                      ntimes;
          getTimeList(path2times.c_str(), timeSeries, ntimes);

          if(!rank)
          {
            //create the i index folder
            indexDir = sectionDir + "/" + std::to_string(iplane);
            createDir(mesh->MESH_COMM, indexDir.c_str());

            // create XMF file
            fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_iSec" + std::to_string(iSections->indices[i]) + ".xmf";

            // write XMF file intro
            xmf = fopen(fieldsFileName.c_str(), "w");
            fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
            fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
            fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
            fprintf(xmf, " <Domain>\n");
            fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
            fclose(xmf);
          }

          // loop over times
          for(PetscInt ti=0; ti<ntimes; ti++)
          {
            // HDF5 file with path
            word fileName;

            // HDF5 file w/o path
            word hdfileName;

            // set time precision
            std::stringstream stream;
            stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

            fileName   = indexDir + "/" + thisCaseName() + "_" + "iSec" + std::to_string(iplane) + "_" + stream.str();
            hdfileName = thisCaseName() + "_" + "iSec" + std::to_string(iplane) + "_" + stream.str();

            // open this time section in the XMF file
            if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), 1, my-1, mz-1, "3DSMesh", timeSeries[ti]);

            // Write the data file.
            hid_t     dataspace_id;
            hsize_t   dims[3];
            herr_t    status;

            // write HDF file
            hid_t	file_id;
            if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            // ************************** write mesh points **************************
            dims[0] = mz - 1;
            dims[1]	= my - 1;
            dims[2]	= 1;

            if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

            writeISectionPointsToXMF
            (
                mesh,
                fieldsFileName.c_str(),
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                timeSeries[ti],
                iplane
            );

            if(!rank) status = H5Sclose(dataspace_id);
            if(!rank) status = H5Fclose(file_id);

            // load velocity
            iSectionLoadVector(mesh, iSections, iplane, "U", timeSeries[ti]);

            if(!rank)
            {
              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
              dataspace_id = H5Screate_simple(3, dims, NULL);

              writeISectionVectorToXMF
              (
                  mesh,
                  fieldsFileName.c_str(),
                  hdfileName.c_str(),
                  &file_id,
                  &dataspace_id,
                  timeSeries[ti],
                  "U",
                  iSections->vectorSec
              );

              status = H5Sclose(dataspace_id);
              status = H5Fclose(file_id);
            }

            // load pressure
            iSectionLoadScalar(mesh, iSections, iplane, "p", timeSeries[ti]);

            if(!rank)
            {
              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
              dataspace_id = H5Screate_simple(3, dims, NULL);

              writeISectionScalarToXMF
              (
                  mesh,
                  fieldsFileName.c_str(),
                  hdfileName.c_str(),
                  &file_id,
                  &dataspace_id,
                  timeSeries[ti],
                  "p",
                  iSections->scalarSec
              );

              status = H5Sclose(dataspace_id);
              status = H5Fclose(file_id);
            }

            // load nut
            if(flags.isLesActive)
            {
              iSectionLoadScalar(mesh, iSections, iplane, "nut", timeSeries[ti]);

              if(!rank)
              {
                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                dataspace_id = H5Screate_simple(3, dims, NULL);

                writeISectionScalarToXMF
                (
                    mesh,
                    fieldsFileName.c_str(),
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    timeSeries[ti],
                    "nut",
                    iSections->scalarSec
                );

                status = H5Sclose(dataspace_id);
                status = H5Fclose(file_id);
              }

            }

            // load temperature
            if(flags.isTeqnActive)
            {
              iSectionLoadScalar(mesh, iSections, iplane, "T", timeSeries[ti]);

              if(!rank)
              {
                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                dataspace_id = H5Screate_simple(3, dims, NULL);

                writeISectionScalarToXMF
                (
                    mesh,
                    fieldsFileName.c_str(),
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    timeSeries[ti],
                    "T",
                    iSections->scalarSec
                );

                status = H5Sclose(dataspace_id);
                status = H5Fclose(file_id);
              }

            }

            // load nut
            if(flags.isIBMActive)
            {
              iSectionLoadScalar(mesh, iSections, iplane, "nv", timeSeries[ti]);

              if(!rank)
              {
                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                dataspace_id = H5Screate_simple(3, dims, NULL);

                writeISectionScalarToXMF
                (
                    mesh,
                    fieldsFileName.c_str(),
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    timeSeries[ti],
                    "nv",
                    iSections->scalarSec
                );

                status = H5Sclose(dataspace_id);
                status = H5Fclose(file_id);
              }

            }

            // close this time section in the XMF file
            if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

            // wait all processes
            MPI_Barrier(mesh->MESH_COMM);
          }

          if(!rank)
          {
              // write XMF file end
              xmf = fopen(fieldsFileName.c_str(), "a");
              fprintf(xmf, "   </Grid>\n");
              fprintf(xmf, " </Domain>\n");
              fprintf(xmf, "</Xdmf>\n");
              fclose(xmf);
          }
        }

        PetscPrintf(mesh->MESH_COMM, "done\n\n");

      }
    }
  }
  return(0);
}

//***************************************************************************************************************//

PetscErrorCode fieldISectionsToXMF(domain_ *domain)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;
        io_       *io      = domain[d].io;

        if
        (
            flags.isAquisitionActive &&
            (io->averaging || io->phaseAveraging)
        )
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;

            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
              char error[512];
              sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
              fatalErrorInFunction("fieldISectionsToXMF",  error);
            }

            if(acquisition->iSections->available)
            {
                PetscPrintf(mesh->MESH_COMM, "On-the-fly i-sections extraction for mesh: %s...\n", mesh->meshName.c_str());
                sections *iSections = acquisition->iSections;

                // create jSections folder
                sectionDir = "./XMF/" + mesh->meshName + "/iSections";
                createDirNoRemove(mesh->MESH_COMM, sectionDir.c_str());

                for(PetscInt i=0; i<iSections->nSections; i++)
                {
                  PetscInt iplane = iSections->indices[i];

                  // exclude ghost nodes
                  if (iplane<1 || iplane>mx-2) continue;

                  // get list of available times
                  word path2times = "./fields/" + mesh->meshName;

                  std::vector<PetscReal>        timeSeries;
                  PetscInt                      ntimes;
                  getTimeList(path2times.c_str(), timeSeries, ntimes);

                  if(!rank)
                  {
                    //create the j index folder
                    indexDir = sectionDir + "/" + std::to_string(iplane);
                    createDirNoRemove(mesh->MESH_COMM, indexDir.c_str());

                    // create XMF file
                    fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_iSec" + std::to_string(iplane) + "_averages.xmf";

                    // write XMF file intro
                    xmf = fopen(fieldsFileName.c_str(), "w");
                    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                    fprintf(xmf, " <Domain>\n");
                    fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                    fclose(xmf);
                  }

                  // now we get only the last time (the meaningful one)

                  // HDF5 file with path
                  word fileName;

                  // HDF5 file w/o path
                  word hdfileName;

                  // set time precision
                  std::stringstream stream;
                  stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ntimes-1];

                  fileName   = indexDir + "/" + thisCaseName() + "_" + "iSec" + std::to_string(iplane) + "_averages_" + stream.str();
                  hdfileName = thisCaseName() + "_" + "iSec" + std::to_string(iplane) + "_averages_" + stream.str();

                  // open this time section in the XMF file
                  if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), 1, my-1, mz-1, "3DSMesh", timeSeries[ntimes-1]);

                  // Write the data file.
                  hid_t     dataspace_id;
                  hsize_t   dims[3];
                  herr_t    status;

                  // write HDF file
                  hid_t	file_id;
                  if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                  // ************************** write mesh points **************************
                  dims[0] = mz - 1;
                  dims[1] = my - 1;
                  dims[2] = 1;

                  if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                  writeISectionPointsToXMF
                  (
                      mesh,
                      fieldsFileName.c_str(),
                      hdfileName.c_str(),
                      &file_id,
                      &dataspace_id,
                      timeSeries[ntimes-1],
                      iplane
                  );

                  if(!rank) status = H5Sclose(dataspace_id);
                  if(!rank) status = H5Fclose(file_id);

                  if(io->averaging)
                  {
                      // load average velocity
                      iSectionLoadVectorFromField(acquisition->fields->avgU, mesh, iSections, iplane, "avgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgU",
                              iSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      iSectionLoadSymmTensorFromField(acquisition->fields->avgUU, mesh, iSections, iplane, "avgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgUU",
                              iSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      iSectionLoadScalarFromField(acquisition->fields->avgP, mesh, iSections, iplane, "avgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "avgP",
                              iSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          iSectionLoadScalarFromField(acquisition->fields->avgNut, mesh, iSections, iplane, "avgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeISectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgNut",
                                  iSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          iSectionLoadScalarFromField(acquisition->fields->avgCs, mesh, iSections, iplane, "avgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeISectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "avgCs",
                                  iSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }
                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  if(io->phaseAveraging)
                  {
                      // load average velocity
                      iSectionLoadVectorFromField(acquisition->fields->pAvgU, mesh, iSections, iplane, "phAvgU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionVectorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgU",
                              iSections->vectorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average Reynold stresses
                      iSectionLoadSymmTensorFromField(acquisition->fields->pAvgUU, mesh, iSections, iplane, "phAvgUU", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionSymmTensorToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgUU",
                              iSections->symmTensorSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      // load average pressure
                      iSectionLoadScalarFromField(acquisition->fields->pAvgP, mesh, iSections, iplane, "phAvgP", timeSeries[ntimes-1]);

                      if(!rank)
                      {
                          file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                          dataspace_id = H5Screate_simple(3, dims, NULL);

                          writeISectionScalarToXMF
                          (
                              mesh,
                              fieldsFileName.c_str(),
                              hdfileName.c_str(),
                              &file_id,
                              &dataspace_id,
                              timeSeries[ntimes-1],
                              "phAvgP",
                              iSections->scalarSec
                          );

                          status = H5Sclose(dataspace_id);
                          status = H5Fclose(file_id);
                      }

                      if(flags.isLesActive)
                      {
                          // load average nut
                          iSectionLoadScalarFromField(acquisition->fields->pAvgNut, mesh, iSections, iplane, "phAvgNut", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeISectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgNut",
                                  iSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }

                          // load average cs
                          iSectionLoadScalarFromField(acquisition->fields->pAvgCs, mesh, iSections, iplane, "phAvgCs", timeSeries[ntimes-1]);

                          if(!rank)
                          {
                              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                              dataspace_id = H5Screate_simple(3, dims, NULL);

                              writeISectionScalarToXMF
                              (
                                  mesh,
                                  fieldsFileName.c_str(),
                                  hdfileName.c_str(),
                                  &file_id,
                                  &dataspace_id,
                                  timeSeries[ntimes-1],
                                  "phAvgCs",
                                  iSections->scalarSec
                              );

                              status = H5Sclose(dataspace_id);
                              status = H5Fclose(file_id);
                          }
                      }

                      // load temperature
                      if(flags.isTeqnActive)
                      {
                          // no fields for now
                      }
                  }

                  // close this time section in the XMF file
                  if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                  // wait all processes
                  MPI_Barrier(mesh->MESH_COMM);

                    if(!rank)
                    {
                        // write XMF file end
                        xmf = fopen(fieldsFileName.c_str(), "a");
                        fprintf(xmf, "   </Grid>\n");
                        fprintf(xmf, " </Domain>\n");
                        fprintf(xmf, "</Xdmf>\n");
                        fclose(xmf);
                    }
                }
            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryKSectionsPerturbToXMF(domain_ *domain)
{
    PetscInt    nDomains = domain[0].info.nDomains;
    word        meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_      flags    = domain[d].flags;

        if(flags.isAquisitionActive)
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;
            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
                fatalErrorInFunction("binaryKSectionsPerturbToXMF",  error);
            }

            if(acquisition->kSections->available && acquisition->isPerturbABLActive)
            {
              PetscPrintf(mesh->MESH_COMM, "Processing k-sections (perturbation fields) for mesh: %s...", mesh->meshName.c_str());
              sections *kSections = acquisition->kSections;

              // create kSections folder
              sectionDir = "./XMF/" + mesh->meshName + "/kSections";

              for(PetscInt k=0; k<kSections->nSections; k++)
              {
                PetscInt kplane = kSections->indices[k];

                // exclude ghost nodes
                if (kplane<1 || kplane>mz-2) continue;

                // get list of available times
                word path2times = "./postProcessing/" + mesh->meshName + "/kSurfaces/" + std::to_string(kplane) + "/UpABL" ;

                // see if can access the directory (shouldn't give error
                // as samplig/kSection could have been created just for averages)
                DIR* dir;
                if ((dir = opendir(path2times.c_str())) == nullptr)
                {
                    PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
                    return(0);
                }

                std::vector<PetscReal>      timeSeries;
                PetscInt                      ntimes;
                getTimeList(path2times.c_str(), timeSeries, ntimes);

                if(!rank)
                {
                  //create the k index folder
                  indexDir = sectionDir + "/" + std::to_string(kplane);

                  // create XMF file
                  fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_kSecPerturb" + std::to_string(kplane) + ".xmf";

                  // write XMF file intro
                  xmf = fopen(fieldsFileName.c_str(), "w");
                  fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                  fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                  fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                  fprintf(xmf, " <Domain>\n");
                  fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                  fclose(xmf);
                }

                // loop over times
                for(PetscInt ti=0; ti<ntimes; ti++)
                {
                  // HDF5 file with path
                  word fileName;

                  // HDF5 file w/o path
                  word hdfileName;

                  std::stringstream stream;
                  stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

                  fileName   = indexDir + "/" + thisCaseName() + "_" + "kSecPerturb" + std::to_string(kplane) + "_" + stream.str();
                  hdfileName = thisCaseName() + "_" + "kSecPerturb" + std::to_string(kplane) + "_" + stream.str();

                  // open this time section in the XMF file
                  if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, my-1, 1,"3DSMesh", timeSeries[ti]);

                  // Write the data file.
                  hid_t     dataspace_id;
                  hsize_t   dims[3];
                  herr_t    status;

                  // write HDF file
                  hid_t	file_id;
                  if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                  // ************************** write mesh points **************************
                  dims[0] = 1;
                  dims[1]	= my - 1;
                  dims[2]	= mx - 1;

                  if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                  writeKSectionPointsToXMF
                  (
                      mesh,
                      fieldsFileName.c_str(),
                      hdfileName.c_str(),
                      &file_id,
                      &dataspace_id,
                      timeSeries[ti],
                      kplane
                  );

                  if(!rank) status = H5Sclose(dataspace_id);
                  if(!rank) status = H5Fclose(file_id);

                  // load velocity
                  kSectionLoadVector(mesh, kSections, kplane, "UpABL", timeSeries[ti]);

                  if(!rank)
                  {
                    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeKSectionVectorToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        "Up",
                        kSections->vectorSec
                    );

                    status = H5Sclose(dataspace_id);
                    status = H5Fclose(file_id);
                  }

                  // load pressure
                  kSectionLoadScalar(mesh, kSections, kplane, "PpABL", timeSeries[ti]);

                  if(!rank)
                  {
                    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeKSectionScalarToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        "Pp",
                        kSections->scalarSec
                    );

                    status = H5Sclose(dataspace_id);
                    status = H5Fclose(file_id);
                  }

                  // load Tp
                    kSectionLoadScalar(mesh, kSections, kplane, "TpABL", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeKSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "Tp",
                          kSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                  // close this time section in the XMF file
                  if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                  // wait all processes
                  MPI_Barrier(mesh->MESH_COMM);
                }

                if(!rank)
                {
                    // write XMF file end
                    xmf = fopen(fieldsFileName.c_str(), "a");
                    fprintf(xmf, "   </Grid>\n");
                    fprintf(xmf, " </Domain>\n");
                    fprintf(xmf, "</Xdmf>\n");
                    fclose(xmf);
                }

              }

              PetscPrintf(mesh->MESH_COMM, "done\n\n");

            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryJSectionsPerturbToXMF(domain_ *domain, postProcess *pp)
{
    PetscInt  nDomains = domain[0].info.nDomains;
    word      meshDir, sectionDir, indexDir;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_    flags    = domain[d].flags;

        if(flags.isAquisitionActive)
        {
            // get pointers
            clock_ *clock = domain[d].clock;
            mesh_  *mesh  = domain[d].mesh;

            acquisition_ *acquisition = domain[d].acquisition;
            ueqn_  *ueqn  = domain[d].ueqn;
            peqn_  *peqn  = domain[d].peqn;
            teqn_  *teqn;
            les_   *les;

            if(flags.isTeqnActive) teqn = domain[d].teqn;
            if(flags.isLesActive)  les  = domain[d].les;

            DMDALocalInfo info = mesh->info;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            PetscMPIInt           rank;
            MPI_Comm_rank(mesh->MESH_COMM, &rank);

            word fieldsFileName;
            FILE *xmf;

            meshDir     = "./XMF/" + mesh->meshName;

            // create domain directory within XMF folder
            errno = 0;
            PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
              char error[512];
              sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
              fatalErrorInFunction("binaryJSectionsPerturbToXMF",  error);
            }

            if(acquisition->jSections->available && acquisition->isPerturbABLActive)
            {
                PetscPrintf(mesh->MESH_COMM, "Processing j-sections (perturbation fields) for mesh: %s...", mesh->meshName.c_str());
                sections *jSections = acquisition->jSections;

                // create jSections folder
                sectionDir = "./XMF/" + mesh->meshName + "/jSections";

                for(PetscInt j=0; j<jSections->nSections; j++)
                {
                  PetscInt jplane = jSections->indices[j];

                  // exclude ghost nodes
                  if (jplane<1 || jplane>my-2) continue;

                  // get list of available times
                  word path2times = "./postProcessing/" + mesh->meshName + "/jSurfaces/" + std::to_string(jplane) + "/UpABL" ;

                  // see if can access the directory (shouldn't give error
                  // as samplig/jSection could have been created just for averages)
                  DIR* dir;
                  if ((dir = opendir(path2times.c_str())) == nullptr)
                  {
                      PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
                      return(0);
                  }

                  std::vector<PetscReal>      timeSeries;
                  PetscInt                      ntimes;
                  getTimeList(path2times.c_str(), timeSeries, ntimes);

                  if(!rank)
                  {
                    //create the j index folder
                    indexDir = sectionDir + "/" + std::to_string(jplane);

                    // create XMF file
                    fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_jSecPerturb" + std::to_string(jplane) + ".xmf";

                    // write XMF file intro
                    xmf = fopen(fieldsFileName.c_str(), "w");
                    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
                    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
                    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
                    fprintf(xmf, " <Domain>\n");
                    fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
                    fclose(xmf);
                  }

                  // loop over times
                  for(PetscInt ti=0; ti<ntimes; ti++)
                  {
                    // HDF5 file with path
                    word fileName;

                    // HDF5 file w/o path
                    word hdfileName;

                    // set time precision
                    std::stringstream stream;
                    stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

                    fileName   = indexDir + "/" + thisCaseName() + "_" + "jSecPerturb" + std::to_string(jplane) + "_" + stream.str();
                    hdfileName = thisCaseName() + "_" + "jSecPerturb" + std::to_string(jplane) + "_" + stream.str();

                    // open this time section in the XMF file
                    if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), mx-1, 1, mz-1, "3DSMesh", timeSeries[ti]);

                    // Write the data file.
                    hid_t     dataspace_id;
                    hsize_t   dims[3];
                    herr_t    status;

                    // write HDF file
                    hid_t	file_id;
                    if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                    // ************************** write mesh points **************************
                    dims[0] = mz - 1;
                    dims[1]	= 1;
                    dims[2]	= mx - 1;

                    if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeJSectionPointsToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        jplane
                    );

                    if(!rank) status = H5Sclose(dataspace_id);
                    if(!rank) status = H5Fclose(file_id);

                    // load velocity
                    jSectionLoadVector(mesh, jSections, jplane, "UpABL", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeJSectionVectorToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "Up",
                          jSections->vectorSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                    // load pressure
                    jSectionLoadScalar(mesh, jSections, jplane, "PpABL", timeSeries[ti]);

                    if(!rank)
                    {
                      file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                      dataspace_id = H5Screate_simple(3, dims, NULL);

                      writeJSectionScalarToXMF
                      (
                          mesh,
                          fieldsFileName.c_str(),
                          hdfileName.c_str(),
                          &file_id,
                          &dataspace_id,
                          timeSeries[ti],
                          "Pp",
                          jSections->scalarSec
                      );

                      status = H5Sclose(dataspace_id);
                      status = H5Fclose(file_id);
                    }

                    // load nut
                  jSectionLoadScalar(mesh, jSections, jplane, "TpABL", timeSeries[ti]);

                  if(!rank)
                  {
                    file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    dataspace_id = H5Screate_simple(3, dims, NULL);

                    writeJSectionScalarToXMF
                    (
                        mesh,
                        fieldsFileName.c_str(),
                        hdfileName.c_str(),
                        &file_id,
                        &dataspace_id,
                        timeSeries[ti],
                        "Tp",
                        jSections->scalarSec
                    );

                    status = H5Sclose(dataspace_id);
                    status = H5Fclose(file_id);
                  }

                    // close this time section in the XMF file
                    if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

                    // wait all processes
                    MPI_Barrier(mesh->MESH_COMM);
                  }

                  if(!rank)
                  {
                      // write XMF file end
                      xmf = fopen(fieldsFileName.c_str(), "a");
                      fprintf(xmf, "   </Grid>\n");
                      fprintf(xmf, " </Domain>\n");
                      fprintf(xmf, "</Xdmf>\n");
                      fclose(xmf);
                  }
                }

                PetscPrintf(mesh->MESH_COMM, "done\n\n");

                if(pp->writeRaster)
                {
                    for(PetscInt s=0; s<jSections->nSections; s++)
                    {
                        writeJSectionToRaster(mesh, jSections->indices[s]);
                    }
                }

            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode binaryISectionsPerturbToXMF(domain_ *domain)
{
  PetscInt  nDomains = domain[0].info.nDomains;
  word      meshDir, sectionDir, indexDir;

  for(PetscInt d=0; d<nDomains; d++)
  {
    flags_    flags    = domain[d].flags;

    if(flags.isAquisitionActive)
    {
      // get pointers
      clock_ *clock = domain[d].clock;
      mesh_  *mesh  = domain[d].mesh;

      acquisition_ *acquisition = domain[d].acquisition;
      ueqn_  *ueqn  = domain[d].ueqn;
      peqn_  *peqn  = domain[d].peqn;
      teqn_  *teqn;
      les_   *les;

      if(flags.isTeqnActive) teqn = domain[d].teqn;
      if(flags.isLesActive)  les  = domain[d].les;

      DMDALocalInfo info = mesh->info;
      PetscInt           mx = info.mx, my = info.my, mz = info.mz;

      PetscMPIInt           rank;
      MPI_Comm_rank(mesh->MESH_COMM, &rank);

      word fieldsFileName;
      FILE *xmf;

      meshDir     = "./XMF/" + mesh->meshName;

      // create domain directory within XMF folder
      errno = 0;
      PetscInt dirRes = mkdir(meshDir.c_str(), 0777);
      if(dirRes != 0 && errno != EEXIST)
      {
         char error[512];
          sprintf(error, "could not create mesh directory %s\n", meshDir.c_str());
          fatalErrorInFunction("binaryISectionsPerturbToXMF",  error);
      }

      if(acquisition->iSections->available && acquisition->isPerturbABLActive)
      {
        PetscPrintf(mesh->MESH_COMM, "Processing i-sections (perturbation fields) for mesh: %s...", mesh->meshName.c_str());
        sections *iSections = acquisition->iSections;

        // create iSections folder
        sectionDir = "./XMF/" + mesh->meshName + "/iSections";

        for(PetscInt i=0; i<iSections->nSections; i++)
        {
          PetscInt iplane = iSections->indices[i];

          // exclude ghost nodes
          if (iplane<1 || iplane>mx-2) continue;

          // get list of available times
          word path2times = "./postProcessing/" + mesh->meshName + "/iSurfaces/" + std::to_string(iplane) + "/UpABL" ;

          // see if can access the directory (shouldn't give error
          // as samplig/iSection could have been created just for averages)
          DIR* dir;
          if ((dir = opendir(path2times.c_str())) == nullptr)
          {
              PetscPrintf(mesh->MESH_COMM, "did not find any in postProcessing directory\n\n");
              return(0);
          }

          std::vector<PetscReal>      timeSeries;
          PetscInt                      ntimes;
          getTimeList(path2times.c_str(), timeSeries, ntimes);

          if(!rank)
          {
            //create the i index folder
            indexDir = sectionDir + "/" + std::to_string(iplane);

            // create XMF file
            fieldsFileName = indexDir + "/" + thisCaseName() + "_" + domain[d].mesh->meshName + "_iSecPerturb" + std::to_string(iSections->indices[i]) + ".xmf";

            // write XMF file intro
            xmf = fopen(fieldsFileName.c_str(), "w");
            fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
            fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
            fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
            fprintf(xmf, " <Domain>\n");
            fprintf(xmf, "   <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
            fclose(xmf);
          }

          // loop over times
          for(PetscInt ti=0; ti<ntimes; ti++)
          {
            // HDF5 file with path
            word fileName;

            // HDF5 file w/o path
            word hdfileName;

            // set time precision
            std::stringstream stream;
            stream << std::fixed << std::setprecision(clock->timePrecision) << timeSeries[ti];

            fileName   = indexDir + "/" + thisCaseName() + "_" + "iSecPerturb" + std::to_string(iplane) + "_" + stream.str();
            hdfileName = thisCaseName() + "_" + "iSecPerturb" + std::to_string(iplane) + "_" + stream.str();

            // open this time section in the XMF file
            if(!rank) xmfWriteFileStartTimeSection(xmf, fieldsFileName.c_str(), 1, my-1, mz-1, "3DSMesh", timeSeries[ti]);

            // Write the data file.
            hid_t     dataspace_id;
            hsize_t   dims[3];
            herr_t    status;

            // write HDF file
            hid_t	file_id;
            if(!rank) file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            // ************************** write mesh points **************************
            dims[0] = mz - 1;
            dims[1]	= my - 1;
            dims[2]	= 1;

            if(!rank) dataspace_id = H5Screate_simple(3, dims, NULL);

            writeISectionPointsToXMF
            (
                mesh,
                fieldsFileName.c_str(),
                hdfileName.c_str(),
                &file_id,
                &dataspace_id,
                timeSeries[ti],
                iplane
            );

            if(!rank) status = H5Sclose(dataspace_id);
            if(!rank) status = H5Fclose(file_id);

            // load velocity
            iSectionLoadVector(mesh, iSections, iplane, "UpABL", timeSeries[ti]);

            if(!rank)
            {
              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
              dataspace_id = H5Screate_simple(3, dims, NULL);

              writeISectionVectorToXMF
              (
                  mesh,
                  fieldsFileName.c_str(),
                  hdfileName.c_str(),
                  &file_id,
                  &dataspace_id,
                  timeSeries[ti],
                  "Up",
                  iSections->vectorSec
              );

              status = H5Sclose(dataspace_id);
              status = H5Fclose(file_id);
            }

            // load pressure
            iSectionLoadScalar(mesh, iSections, iplane, "PpABL", timeSeries[ti]);

            if(!rank)
            {
              file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
              dataspace_id = H5Screate_simple(3, dims, NULL);

              writeISectionScalarToXMF
              (
                  mesh,
                  fieldsFileName.c_str(),
                  hdfileName.c_str(),
                  &file_id,
                  &dataspace_id,
                  timeSeries[ti],
                  "Pp",
                  iSections->scalarSec
              );

              status = H5Sclose(dataspace_id);
              status = H5Fclose(file_id);
            }

              iSectionLoadScalar(mesh, iSections, iplane, "TpABL", timeSeries[ti]);

              if(!rank)
              {
                file_id      = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                dataspace_id = H5Screate_simple(3, dims, NULL);

                writeISectionScalarToXMF
                (
                    mesh,
                    fieldsFileName.c_str(),
                    hdfileName.c_str(),
                    &file_id,
                    &dataspace_id,
                    timeSeries[ti],
                    "Tp",
                    iSections->scalarSec
                );

                status = H5Sclose(dataspace_id);
                status = H5Fclose(file_id);
              }

            // close this time section in the XMF file
            if(!rank) xmfWriteFileEndTimeSection(xmf, fieldsFileName.c_str());

            // wait all processes
            MPI_Barrier(mesh->MESH_COMM);
          }

          if(!rank)
          {
              // write XMF file end
              xmf = fopen(fieldsFileName.c_str(), "a");
              fprintf(xmf, "   </Grid>\n");
              fprintf(xmf, " </Domain>\n");
              fprintf(xmf, "</Xdmf>\n");
              fclose(xmf);
          }
        }

        PetscPrintf(mesh->MESH_COMM, "done\n\n");

      }
    }
  }
  return(0);
}

//***************************************************************************************************************//

PetscErrorCode sectionsReadAndAllocate(domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_ flags = domain[d].flags;
        io_    *io   = domain[d].io;

        if(flags.isAquisitionActive)
        {
            acquisition_ *acquisition = domain[d].acquisition;
            mesh_         *mesh       = domain[d].mesh;

            DMDALocalInfo info = mesh->info;
            PetscInt           xs = info.xs, xe = info.xs + info.xm;
            PetscInt           ys = info.ys, ye = info.ys + info.ym;
            PetscInt           zs = info.zs, ze = info.zs + info.zm;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            Cmpnts             ***cent;
            PetscReal          ***aj;

            PetscInt           i, j, k;
            PetscInt           lxs, lxe, lys, lye, lzs, lze;
            PetscMPIInt        rank, nprocs;

            lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
            lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
            lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

            MPI_Comm_rank(mesh->MESH_COMM, &rank);
            MPI_Comm_size(mesh->MESH_COMM, &nprocs);


            if(d ==0)
            {
              PetscPrintf(mesh->MESH_COMM, "Reading surfaces in sampling/surfaces/...");
            }

            DMDAVecGetArray(mesh->fda, mesh->lCent, &cent);
            DMDAVecGetArray(mesh->da,  mesh->lAj, &aj);

            // max perturbation amplitude
            PetscReal maxPerturb  = 1e-10;
            // processor perturbation for search (changes between processors)
            PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

            PetscInt iskavail = 0, isjavail = 0, isiavail = 0;
            PetscInt nSurfaces = 0;

            word dataLoc, kSecName, jSecName, iSecName, userSecPath, userSecName;

            dataLoc = "sampling/surfaces/";

            kSecName = dataLoc + "kSections";
            jSecName = dataLoc + "jSections";
            iSecName = dataLoc + "iSections";

            userSecPath = dataLoc + "userSections";

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
                PetscInt atLeastOneVector     = 0;
                PetscInt atLeastOneScalar     = 0;
                PetscInt atLeastOneSymmTensor = 0;

                word path2index = "./postProcessing/" + mesh->meshName + "/kSurfaces/";
                std::vector<PetscInt>      kIndex;
                PetscInt                   nkSec;
                std::vector<PetscReal>     lcoor;

                // try to read from post procesing files
                DIR *dir;
                PetscInt indicesSet = 0;
                if ((dir = opendir(path2index.c_str())) != nullptr)
                {
                    getIndexList(path2index.c_str(), kIndex, nkSec);
                    indicesSet++;
                }
                // on-the-fly slices: calculate now
                else
                {
                    // read number of sections
                    readDictInt(kSecName.c_str(), "surfaceNumber" ,&nkSec);
                    kIndex.resize(nkSec);
                }

                lcoor.resize(nkSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->kSections));
                PetscMalloc(sizeof(PetscInt) * nkSec,  &(acquisition->kSections->indices));
                PetscMalloc(sizeof(PetscReal) * nkSec,  &(acquisition->kSections->coordinates));

                sections *kSections = acquisition->kSections;

                // set indices if didn't read from file
                if(!indicesSet)
                {
                    // read section indices
                    std::ifstream indata;
                    indata.open(kSecName.c_str());
                    char buffer[256];

                    while(!indata.eof())
                    {
                        indata >> buffer;
                        if(strcmp("coordinates",buffer) == 0)
                        {
                            for(PetscInt s=0; s<nkSec; s++)
                            {
                                indata >> buffer;
                                std::sscanf(buffer, "%lf", &(lcoor[s]));
                            }
                        }
                    }
                    indata.close();

                    // find the closest point, then take its k index
                    for(PetscInt s=0; s<nkSec; s++)
                    {
                        Cmpnts surfacePoint;
                               surfacePoint.x = lcoor[s];
                               surfacePoint.y = (mesh->bounds.ymax + mesh->bounds.ymin) / 2.0;
                               surfacePoint.z = (mesh->bounds.zmax + mesh->bounds.zmin) / 2.0;

                        PetscReal  lminDist = 1e30,  gminDist = 1e30;
                        cellIds    closestIds;
                        PetscInt   lclosestK = 0;

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

                            // exclude section if coordinate is greater than 5*cellWidth
                            PetscReal cellWidth = 5.0*pow(1.0/aj[closestIds.k][(PetscInt)(std::floor(0.5*(lys+lye)))][(PetscInt)(std::floor(0.5*(lxs+lxe)))], 1.0/3.0);

                            if(gminDist > cellWidth)
                            {
                                lclosestK = -1;
                            }
                        }

                        MPI_Allreduce(&lclosestK, &(kIndex[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                        lcoor[s] = 0.0;
                        kSections->coordinates[s] = 0.0;
                    }
                }
                else
                {
                    for(PetscInt s = 0; s< nkSec; s++)
                    {
                        lcoor[s] = 0.0;
                        kSections->coordinates[s] = 0.0;
                    }
                }

                // store number of sections
                kSections->nSections = nkSec;

                // set available to 1
                kSections->available = iskavail;

                for(PetscInt s=0; s<nkSec; s++)
                {
                    kSections->indices[s] = kIndex[s];

                    if((kIndex[s] >= lzs) && (kIndex[s] < lze))
                    {
                        lcoor[s] = cent[kIndex[s]][lys][lxs].x;
                    }
                }

                MPI_Allreduce(&lcoor[0], &kSections->coordinates[0], nkSec, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);

                atLeastOneVector++;
                if(flags.isTeqnActive) atLeastOneScalar++;
                if(flags.isLesActive)  atLeastOneScalar++;
                if(acquisition->isPerturbABLActive) {atLeastOneScalar++; atLeastOneVector++;}
                if(io->averaging || io->phaseAveraging) atLeastOneSymmTensor++;

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

                if(atLeastOneSymmTensor)
                {
                    kSections->symmTensorSec = (symmTensor **)malloc( sizeof(symmTensor *) * my );

                    for(j=0; j<my; j++)
                    {
                        kSections->symmTensorSec[j] = (symmTensor *)malloc( sizeof(symmTensor) * mx );
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
                PetscInt atLeastOneVector     = 0;
                PetscInt atLeastOneScalar     = 0;
                PetscInt atLeastOneSymmTensor = 0;

                // read number of sections
                word path2index = "./postProcessing/" + mesh->meshName + "/jSurfaces/";
                std::vector<PetscInt>      jIndex;
                PetscInt                   njSec;
                std::vector<PetscReal>     lcoor;

                // try to read from post procesing files
                DIR *dir;
                PetscInt indicesSet = 0;
                if ((dir = opendir(path2index.c_str())) != nullptr)
                {
                    getIndexList(path2index.c_str(), jIndex, njSec);
                    indicesSet++;
                }
                // on-the-fly slices: calculate now
                else
                {
                    // read number of sections
                    readDictInt(jSecName.c_str(), "surfaceNumber" ,&njSec);
                    jIndex.resize(njSec);
                }

                lcoor.resize(njSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->jSections));
                PetscMalloc(sizeof(PetscInt) * njSec,  &(acquisition->jSections->indices));
                PetscMalloc(sizeof(PetscReal) * njSec,  &(acquisition->jSections->coordinates));

                sections *jSections = acquisition->jSections;

                // set indices if didn't read from file
                if(!indicesSet)
                {
                    // read section indices
                    std::ifstream indata;
                    indata.open(jSecName.c_str());
                    char buffer[256];

                    while(!indata.eof())
                    {
                        indata >> buffer;

                        if(strcmp("coordinates",buffer) == 0)
                        {
                            for(PetscInt s=0; s<njSec; s++)
                            {
                                indata >> buffer;
                                std::sscanf(buffer, "%lf", &(lcoor[s]));
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
                               surfacePoint.z = lcoor[s];

                        PetscReal  lminDist = 1e30,  gminDist = 1e30;
                        cellIds    closestIds;
                        PetscInt   lclosestJ = 0;

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

                            // exclude section if coordinate is greater than 5*cellWidth
                            PetscReal cellWidth = 5.0*pow(1.0/aj[(PetscInt)(std::floor(0.5*(lzs+lze)))][closestIds.j][(PetscInt)(std::floor(0.5*(lxs+lxe)))], 1.0/3.0);

                            if(gminDist > cellWidth)
                            {
                                lclosestJ = -1;
                            }
                        }

                        MPI_Allreduce(&lclosestJ, &(jIndex[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                        lcoor[s] = 0.0;
                        jSections->coordinates[s] = 0.0;
                    }
                }
                else
                {
                    for(PetscInt s = 0; s< njSec; s++)
                    {
                        lcoor[s] = 0.0;
                        jSections->coordinates[s] = 0.0;
                    }
                }

                // store number of sections
                jSections->nSections = njSec;

                // set available to 1
                jSections->available = isjavail;

                for(PetscInt s=0; s<njSec; s++)
                {
                  jSections->indices[s] = jIndex[s];

                  if((jIndex[s] >= lys) && (jIndex[s] < lye))
                  {
                    lcoor[s] = cent[lzs][jIndex[s]][lxs].z;
                  }
                }

                MPI_Allreduce(&lcoor[0], &jSections->coordinates[0], njSec, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);

                atLeastOneVector++;
                if(flags.isTeqnActive) atLeastOneScalar++;
                if(flags.isLesActive)  atLeastOneScalar++;
                if(acquisition->isPerturbABLActive) {atLeastOneScalar++; atLeastOneVector++;}
                if(io->averaging || io->phaseAveraging) atLeastOneSymmTensor++;

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

                if(atLeastOneSymmTensor)
                {
                    jSections->symmTensorSec = (symmTensor **)malloc( sizeof(symmTensor *) * mz );

                    for(k=0; k<mz; k++)
                    {
                        jSections->symmTensorSec[k] = (symmTensor *)malloc( sizeof(symmTensor) * mx );
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
                PetscInt atLeastOneVector     = 0;
                PetscInt atLeastOneScalar     = 0;
                PetscInt atLeastOneSymmTensor = 0;

                // read number of sections
                word path2index = "./postProcessing/" + mesh->meshName + "/iSurfaces/";
                std::vector<PetscInt>      iIndex;
                PetscInt                   niSec;
                std::vector<PetscReal>   lcoor;

                // try to read from post procesing files
                DIR *dir;
                PetscInt indicesSet = 0;
                if ((dir = opendir(path2index.c_str())) != nullptr)
                {
                    getIndexList(path2index.c_str(), iIndex, niSec);
                    indicesSet++;
                }
                // on-the-fly slices: calculate now
                else
                {
                    // read number of sections
                    readDictInt(iSecName.c_str(), "surfaceNumber" ,&niSec);
                    iIndex.resize(niSec);
                }

                lcoor.resize(niSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->iSections));
                PetscMalloc(sizeof(PetscInt) * niSec,  &(acquisition->iSections->indices));
                PetscMalloc(sizeof(PetscReal) * niSec,  &(acquisition->iSections->coordinates));

                sections *iSections = acquisition->iSections;

                // set indices if didn't read from file
                if(!indicesSet)
                {
                    // read section indices
                    std::ifstream indata;
                    indata.open(iSecName.c_str());
                    char buffer[256];

                    while(!indata.eof())
                    {
                        indata >> buffer;
                        if(strcmp("coordinates",buffer) == 0)
                        {
                            for(PetscInt s=0; s<niSec; s++)
                            {
                                indata >> buffer;
                                std::sscanf(buffer, "%lf", &(lcoor[s]));
                            }
                        }
                    }
                    indata.close();

                    // find the closest point, then take is k index
                    for(PetscInt s=0; s<niSec; s++)
                    {
                        Cmpnts surfacePoint;
                               surfacePoint.x = (mesh->bounds.xmax + mesh->bounds.xmin) / 2.0;
                               surfacePoint.y = lcoor[s];
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

                            // exclude section if coordinate is greater than 5*cellWidth
                            PetscReal cellWidth = 5.0*pow(1.0/aj[(PetscInt)(std::floor(0.5*(lzs+lze)))][(PetscInt)(std::floor(0.5*(lys+lye)))][closestIds.i], 1.0/3.0);

                            if(gminDist > cellWidth)
                            {
                                lclosestI = -1;
                            }
                        }

                        MPI_Allreduce(&lclosestI, &(iIndex[s]), 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                        lcoor[s] = 0.0;
                        iSections->coordinates[s] = 0.0;
                    }
                }
                else
                {
                    for(PetscInt s = 0; s< niSec; s++)
                    {
                        lcoor[s] = 0.0;
                        iSections->coordinates[s] = 0.0;
                    }
                }

                // store number of sections
                iSections->nSections = niSec;

                // set available to 1
                iSections->available = isiavail;

                for(PetscInt s=0; s<niSec; s++)
                {
                    iSections->indices[s] = iIndex[s];

                    if((iIndex[s] >= lxs) && (iIndex[s] < lxe))
                    {
                        lcoor[s] = cent[lzs][lys][iIndex[s]].y;
                    }
                }

                MPI_Allreduce(&lcoor[0], &iSections->coordinates[0], niSec, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);

                atLeastOneVector++;
                if(flags.isTeqnActive) atLeastOneScalar++;
                if(flags.isLesActive)  atLeastOneScalar++;
                if(acquisition->isPerturbABLActive) {atLeastOneScalar++; atLeastOneVector++;}
                if(io->averaging || io->phaseAveraging) atLeastOneSymmTensor++;

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

                if(atLeastOneSymmTensor)
                {
                    iSections->symmTensorSec = (symmTensor **)malloc( sizeof(symmTensor *) * mz );

                    for(k=0; k<mz; k++)
                    {
                        iSections->symmTensorSec[k] = (symmTensor *)malloc( sizeof(symmTensor) * my );
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

            // get list of times from the fieldsPath of the current domain
            std::vector<word> surfaceSeries;
            getFileList(userSecPath.c_str(), surfaceSeries, nSurfaces);

            // read userDefined Sections input file
            if(nSurfaces)
            {
                acquisition->userSections = new udSections;

                // store number of sections
                acquisition->userSections->nSections = nSurfaces;

                acquisition->userSections->available = 1;

                acquisition->userSections->uSection = new uSections*[nSurfaces];

                for(PetscInt s=0; s<nSurfaces; s++)
                {

                    PetscInt atLeastOneVector     = 0;
                    PetscInt atLeastOneScalar     = 0;

                    //reverses the number of cells in the two surface directions - test feature needs to be removed
                    PetscInt flipIndexOrder;

                    // allocate memory
                    acquisition->userSections->uSection[s] = new uSections;
                    
                    uSections *uSection = acquisition->userSections->uSection[s];
                    
                    userSecName = userSecPath + "/" + surfaceSeries[s];

                    //set surface name 
                    uSection->sectionName = surfaceSeries[s]; 

                    // read acquisition start time and type of interval
                    readDictDouble(userSecName.c_str(), "timeStart", &(uSection->timeStart));
                    readDictWord(userSecName.c_str(),   "intervalType", &(uSection->intervalType));
                    readDictDouble(userSecName.c_str(), "timeInterval", &(uSection->timeInterval));
                    readDictInt(userSecName.c_str(), "flipIndexOrder", &(flipIndexOrder));

                    // check if intervalType is known
                    if(uSection->intervalType != "timeStep" && uSection->intervalType != "adjustableTime")
                    {
                        char error[512];
                        sprintf(error, "unknown interval type %s. Known types are timeStep and adjustableTime\n", uSection->intervalType.c_str());
                        fatalErrorInFunction("sectionsInitialize",  error);
                    }

                    // read section indices
                    std::ifstream indata;

                    indata.open(userSecName.c_str());

                    char buffer[256];

                    while(!indata.eof())
                    {
                        indata >> buffer;

                        if
                        (
                            strcmp
                            (
                                "meshPoints",
                                buffer
                            ) == 0
                        )
                        {
                            indata >> buffer;
                            std::sscanf(buffer, "%ld", &(uSection->ny));

                            indata >> buffer;
                            std::sscanf(buffer, "%ld", &(uSection->nx));

                            break;
                        }

                    }

                    if(flipIndexOrder == 1)
                    { 
                        PetscInt tempIndex;
                        tempIndex = uSection->ny;
                        uSection->ny = uSection->nx;
                        uSection->nx = tempIndex;
                    }

                    // allocate memory for the co-ordinates of the points in the user defined surface
                    PetscInt nx = uSection->nx, ny = uSection->ny;

                    uSection->coor = (Cmpnts **)malloc( sizeof(Cmpnts *) * (ny) );

                    for(j=0; j<ny; j++)
                    {
                        uSection->coor[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * (nx) );
                    }

                    for(j=0; j<ny; j++)
                    {
                        for(i=0; i<nx; i++)
                        {
                            indata >> buffer;
                            std::sscanf(buffer, "%le", &(uSection->coor[j][i].x));

                            indata >> buffer;
                            std::sscanf(buffer, "%le", &(uSection->coor[j][i].y));

                            indata >> buffer;
                            std::sscanf(buffer, "%le", &(uSection->coor[j][i].z));
                        }
                    }

                    indata.close();

                    //allocate memory 
                    uSection->closestId = (cellIds **)malloc( sizeof(cellIds *) * (ny) );

                    for(j=0; j<ny; j++)
                    {
                        uSection->closestId[j] = (cellIds *)malloc( sizeof(cellIds) * (nx) );
                    }

                    uSection->hasCoor = (PetscInt **)malloc( sizeof(PetscInt *) * (ny) );

                    for(j=0; j<ny; j++)
                    {
                        uSection->hasCoor[j] = (PetscInt *)malloc( sizeof(PetscInt) * (nx) );
                    }

                    cellIds **localId;
                    localId = (cellIds **)malloc( sizeof(cellIds *) * (ny) );

                    for(j=0; j<ny; j++)
                    {
                        localId[j] = (cellIds *)malloc( sizeof(cellIds) * (nx) );
                    }



                    for(j=0; j<ny; j++)
                    {
                        for(i=0; i<nx; i++)
                        {
                        localId[j][i].i = -100;
                        localId[j][i].j = -100;
                        localId[j][i].k = -100;

                        uSection->hasCoor[j][i] = 0;
                        }
                    }


                    for(PetscInt l=0; l<ny; l++)
                    {
                        for(PetscInt m=0; m<nx; m++)
                        {
                            Cmpnts surfacePoint = nSet(uSection->coor[l][m]);

                            PetscReal  lminDist = 1e30,  gminDist = 1e30;
                            cellIds closestIds;

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
                                localId[l][m].i = closestIds.i;
                                localId[l][m].j = closestIds.j;
                                localId[l][m].k = closestIds.k;

                                PetscReal cellWidth = 5.0*pow(1.0/aj[closestIds.k][closestIds.j][closestIds.i], 1.0/3.0);

                                if(gminDist > cellWidth)
                                {
                                    char warning[256];
                                    sprintf(warning, "User defined section has parts outside the domain, consider redefining the surface");
                                    warningInFunction("sectionsInitialize",  warning);
                                }

                                uSection->hasCoor[l][m] = 1;
                            }
                        }

                        MPI_Allreduce(&(localId[l][0]), &(uSection->closestId[l][0]), 3*nx, MPIU_INT, MPIU_MAX, mesh->MESH_COMM);
                    }

                    // free the allocated memory for coorPts
                    for(j = 0; j<ny; j++)
                    {
                        free(localId[j]);
                    }

                    free(localId);  

                    atLeastOneVector++;
                    if(flags.isTeqnActive) atLeastOneScalar++;
                    if(flags.isLesActive)  atLeastOneScalar++;
                    if(acquisition->isPerturbABLActive) {atLeastOneScalar++; atLeastOneVector++;}

                    if(atLeastOneVector)
                    {
                        // allocate memory for scalar and vector variables
                        uSection->vectorSec = (Cmpnts **)malloc( sizeof(Cmpnts *) * ny );

                        for(j=0; j<ny; j++)
                        {
                            uSection->vectorSec[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * nx );
                        }
                    }

                    if(atLeastOneScalar)
                    {
                        uSection->scalarSec = (PetscReal **)malloc( sizeof(PetscReal *) * ny );

                        for(j=0; j<ny; j++)
                        {
                            uSection->scalarSec[j] = (PetscReal *)malloc( sizeof(PetscReal) * nx );
                        }
                    }
                }
                
            }
            else
            {
                // allocate memory
                acquisition->userSections = new udSections;
                // set available to 0
                acquisition->userSections->available = 0;
            }

            DMDAVecRestoreArray(mesh->fda, mesh->lCent, &cent);
            DMDAVecRestoreArray(mesh->da,  mesh->lAj, &aj);

            if(d == 0)
            {
              PetscPrintf(mesh->MESH_COMM, "done\n\n");
            }

            MPI_Barrier(mesh->MESH_COMM);

        }

    }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode postProcessWriteProbes(domain_ *domain)
{
    if
    (
        domain[0].access.io->averaging ||
        domain[0].access.io->phaseAveraging
    )
    if(domain[0].acquisition->isProbesActive)
    {
        PetscInt  nDomains = domain[0].info.nDomains;

        // read the fields corresponding to the last time step on each domain
        for(PetscInt d=0; d<nDomains; d++)
        {
            PetscPrintf(domain[d].mesh->MESH_COMM, "On-the-fly probes extraction for mesh: %s...\n", domain[d].mesh->meshName.c_str());
            PetscInt     ntimes;
            word         fieldsPath = "./fields/" + domain[d].mesh->meshName;
            std::vector<PetscReal> timeSeries;
            getTimeList(fieldsPath.c_str(), timeSeries, ntimes);
            readFields(&domain[d], timeSeries[ntimes-1]);
        }

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
                // initialize local vectors
                std::vector<Cmpnts> lprobeValuesU(rake->probesNumber);
                std::vector<Cmpnts> gprobeValuesU(rake->probesNumber);

                // initialize local vectors
                std::vector<symmTensor> lprobeValuesUU(rake->probesNumber);
                std::vector<symmTensor> gprobeValuesUU(rake->probesNumber);

                // initialize local vectors
                std::vector<PetscReal> lprobeValuesP(rake->probesNumber);
                std::vector<PetscReal> gprobeValuesP(rake->probesNumber);

                // initialize local vectors
                std::vector<Cmpnts> lprobeValuesUph(rake->probesNumber);
                std::vector<Cmpnts> gprobeValuesUph(rake->probesNumber);

                // initialize local vectors
                std::vector<symmTensor> lprobeValuesUUph(rake->probesNumber);
                std::vector<symmTensor> gprobeValuesUUph(rake->probesNumber);

                // initialize local vectors
                std::vector<PetscReal> lprobeValuesPph(rake->probesNumber);
                std::vector<PetscReal> gprobeValuesPph(rake->probesNumber);

                for(p=0; p<rake->probesNumber; p++)
                {
                    if(rake->domainID[p] != -1)
                    {
                        io_               *io = domain[rake->domainID[p]].io;
                        acquisition_ *acquisition = domain[rake->domainID[p]].acquisition;
                        avgFields     *fields = acquisition->fields;
                        mesh_           *mesh = domain[rake->domainID[p]].mesh;
                        DMDALocalInfo    info = mesh->info;
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

                        if(io->averaging)
                        {
                            if(meshRank == rake->owner[p])
                            {
                                k = rake->cells[p][0];
                                j = rake->cells[p][1];
                                i = rake->cells[p][2];

                                Cmpnts ***uavg;
                                DMDAVecGetArray(mesh->fda, fields->avgU, &uavg);

                                lprobeValuesU[p].x = uavg[k][j][i].x;
                                lprobeValuesU[p].y = uavg[k][j][i].y;
                                lprobeValuesU[p].z = uavg[k][j][i].z;

                                DMDAVecRestoreArray(mesh->fda, fields->avgU, &uavg);

                                symmTensor ***uuavg;
                                DMDAVecGetArray(mesh->sda, fields->avgUU, &uuavg);

                                lprobeValuesUU[p].xx = uuavg[k][j][i].xx;
                                lprobeValuesUU[p].yy = uuavg[k][j][i].yy;
                                lprobeValuesUU[p].zz = uuavg[k][j][i].zz;
                                lprobeValuesUU[p].xy = uuavg[k][j][i].xy;
                                lprobeValuesUU[p].xz = uuavg[k][j][i].xz;
                                lprobeValuesUU[p].yz = uuavg[k][j][i].yz;

                                DMDAVecRestoreArray(mesh->sda, fields->avgUU, &uuavg);

                                PetscReal ***pavg;
                                DMDAVecGetArray(mesh->da, fields->avgP, &pavg);

                                lprobeValuesP[p] = pavg[k][j][i];

                                DMDAVecRestoreArray(mesh->da, fields->avgP, &pavg);
                            }
                        }

                        if(io->phaseAveraging)
                        {
                            if(meshRank == rake->owner[p])
                            {
                                k = rake->cells[p][0];
                                j = rake->cells[p][1];
                                i = rake->cells[p][2];

                                Cmpnts ***uavg;
                                DMDAVecGetArray(mesh->fda, fields->pAvgU, &uavg);

                                lprobeValuesUph[p].x = uavg[k][j][i].x;
                                lprobeValuesUph[p].y = uavg[k][j][i].y;
                                lprobeValuesUph[p].z = uavg[k][j][i].z;

                                DMDAVecRestoreArray(mesh->fda, fields->pAvgU, &uavg);

                                symmTensor ***uuavg;
                                DMDAVecGetArray(mesh->sda, fields->pAvgUU, &uuavg);

                                lprobeValuesUUph[p].xx = uuavg[k][j][i].xx;
                                lprobeValuesUUph[p].yy = uuavg[k][j][i].yy;
                                lprobeValuesUUph[p].zz = uuavg[k][j][i].zz;
                                lprobeValuesUUph[p].xy = uuavg[k][j][i].xy;
                                lprobeValuesUUph[p].xz = uuavg[k][j][i].xz;
                                lprobeValuesUUph[p].yz = uuavg[k][j][i].yz;

                                DMDAVecRestoreArray(mesh->sda, fields->pAvgUU, &uuavg);

                                PetscReal ***pavg;
                                DMDAVecGetArray(mesh->da, fields->pAvgP, &pavg);

                                lprobeValuesPph[p] = pavg[k][j][i];

                                DMDAVecRestoreArray(mesh->da, fields->pAvgP, &pavg);
                            }
                        }
                    }
                }

                if(domain[0].access.io->averaging)
                {
                    MPI_Reduce(&lprobeValuesU[0], &gprobeValuesU[0], 3*rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);
                    MPI_Reduce(&lprobeValuesUU[0], &gprobeValuesUU[0], 6*rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);
                    MPI_Reduce(&lprobeValuesP[0], &gprobeValuesP[0], rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);

                    PetscMPIInt rakeRank;
                    MPI_Comm_rank(rake->RAKE_COMM, &rakeRank);

                    // only master process writes on the file avgU
                    if(!rakeRank)
                    {
                        FILE *f;
                        word fileName;

                        // write velocity
                        fileName = rake->timeName + "/avgU";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "(%.10lf  %.10lf  %.10lf)\t\t", gprobeValuesU[p].x, gprobeValuesU[p].y, gprobeValuesU[p].z);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);

                        // write stresses
                        fileName = rake->timeName + "/avgUU";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "(%.10lf  %.10lf  %.10lf %.10lf  %.10lf  %.10lf)\t\t", gprobeValuesUU[p].xx, gprobeValuesUU[p].yy, gprobeValuesUU[p].zz, gprobeValuesUU[p].xy, gprobeValuesUU[p].xz, gprobeValuesUU[p].yz);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);

                        // write pressure
                        fileName = rake->timeName + "/avgP";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "%.10lf\t\t", gprobeValuesP[p]);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);
                    }
                }

                if(domain[0].access.io->phaseAveraging)
                {
                    MPI_Reduce(&lprobeValuesUph[0], &gprobeValuesUph[0], 3*rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);
                    MPI_Reduce(&lprobeValuesUUph[0], &gprobeValuesUUph[0], 6*rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);
                    MPI_Reduce(&lprobeValuesPph[0], &gprobeValuesPph[0], rake->probesNumber, MPIU_REAL, MPIU_SUM, 0, rake->RAKE_COMM);

                    PetscMPIInt rakeRank;
                    MPI_Comm_rank(rake->RAKE_COMM, &rakeRank);

                    // only master process writes on the file avgU
                    if(!rakeRank)
                    {
                        FILE *f;
                        word fileName;

                        // write velocity
                        fileName = rake->timeName + "/phAvgU";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "(%.10lf  %.10lf  %.10lf)\t\t", gprobeValuesUph[p].x, gprobeValuesUph[p].y, gprobeValuesUph[p].z);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);

                        // write stresses
                        fileName = rake->timeName + "/phAvgUU";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "(%.10lf  %.10lf  %.10lf %.10lf  %.10lf  %.10lf)\t\t", gprobeValuesUUph[p].xx, gprobeValuesUUph[p].yy, gprobeValuesUUph[p].zz, gprobeValuesUUph[p].xy, gprobeValuesUUph[p].xz, gprobeValuesUUph[p].yz);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);

                        // write pressure
                        fileName = rake->timeName + "/phAvgP";
                        f = fopen(fileName.c_str(), "a");

                        fprintf(f, "\t %.3lf\t\t\t\t\t\t", clock->time);
                        for(p=0; p<rake->probesNumber; p++)
                        {
                            if(rake->domainID[p] != -1)
                            {
                                fprintf(f, "%.10lf\t\t", gprobeValuesPph[p]);
                            }
                        }
                        fprintf(f, "\n");
                        fclose(f);
                    }
                }

                std::vector<Cmpnts> ().swap(lprobeValuesU);
                std::vector<Cmpnts> ().swap(gprobeValuesU);

                std::vector<symmTensor> ().swap(lprobeValuesUU);
                std::vector<symmTensor> ().swap(gprobeValuesUU);

                std::vector<PetscReal> ().swap(lprobeValuesP);
                std::vector<PetscReal> ().swap(gprobeValuesP);

                std::vector<Cmpnts> ().swap(lprobeValuesUph);
                std::vector<Cmpnts> ().swap(gprobeValuesUph);

                std::vector<symmTensor> ().swap(lprobeValuesUUph);
                std::vector<symmTensor> ().swap(gprobeValuesUUph);

                std::vector<PetscReal> ().swap(lprobeValuesPph);
                std::vector<PetscReal> ().swap(gprobeValuesPph);
            }
        }
    }

    return(0);
}
