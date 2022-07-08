
//#include </usr/include/hdf5/serial/hdf5.h>
#include <hdf5.h>
#include "include/base.h"
#include "include/domain.h"
#include "include/inline.h"
#include "include/initialization.h"
#include "include/windToPW.h"

static char head[] = "MARBLLES Post Processor";

//***************************************************************************************************************//

#undef __FUNCT__
#define __FUNCT__ "main"

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

  PetscMPIInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

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
  postProcessInitialize(&domain, &clock, &info, &flags);

  // write 3D fields into XMF
  binary3DToXMF(domain, &pp);

  // initialize precursor post processing
  if(pp.postProcessPrecursor)
  {
      postProcessInitializePrecursor(&pp, &clock);

      binary3DToXMF(pp.precursor->domain, &pp);
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

      // initialize ibm
      InitializeIBM(domain[d].ibm);

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

          averageFieldsInitialize(domain[d].acquisition);
          averageKEBudgetsInitialize(domain[d].acquisition);
          perturbationABLInitialize(domain[d].acquisition);
      }

      PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");
  }

  averaging3LMInitialize(domain);

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

    // read physical constants
    ReadPhysicalConstants(domain);

    // set solution flags
    SetSolutionFlagsPrecursor(domain);

    // set simulation info
    SetSimulationInfo(&(domain->info));

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
          if(pp->postProcessFields || pp->writeRaster)
          {
            readFields(&domain[d], timeSeries[ti]);

            if(pp->postProcessFields)
            {
              // write fields inside XMF folder
              writeFieldsToXMF(&domain[d], fieldsFileName.c_str(), timeSeries[ti]);
            }
          }
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
            "SideForce",
            acquisition->fields->SideForce
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
            "pAvgU",
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
            "pAvgP",
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
            "pAvgUU",
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
                "pAvgNut",
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
                "pAvgCs",
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
                "pAvgOmega",
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
                    "pAvgPsq",
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
                    "pAvgOmegaOmega",
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
                    "pAvgUdotGradP",
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
                    "pAvgMagGradU",
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
                    "pAvgMagUU",
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
        fatalErrorInFunction("getTimeList", error);
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
                fatalErrorInFunction("binary3DToXMF",  error);
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
              fatalErrorInFunction("binary3DToXMF",  error);
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
          fatalErrorInFunction("binary3DToXMF",  error);
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
                fatalErrorInFunction("binary3DToXMF",  error);
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
              fatalErrorInFunction("binary3DToXMF",  error);
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
          fatalErrorInFunction("binary3DToXMF",  error);
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

        if(flags.isAquisitionActive)
        {
            acquisition_ *acquisition = domain[d].acquisition;
            mesh_         *mesh       = domain[d].mesh;

            DMDALocalInfo info = mesh->info;
            PetscInt           xs = info.xs, xe = info.xs + info.xm;
            PetscInt           ys = info.ys, ye = info.ys + info.ym;
            PetscInt           zs = info.zs, ze = info.zs + info.zm;
            PetscInt           mx = info.mx, my = info.my, mz = info.mz;

            Cmpnts        ***cent;

            PetscInt           i, j, k;
            PetscInt           lxs, lxe, lys, lye, lzs, lze;
            PetscMPIInt           rank, nprocs;

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
                word path2index = "./postProcessing/" + mesh->meshName + "/kSurfaces/";
                std::vector<PetscInt>      kIndex;
                PetscInt                   nkSec;
                std::vector<PetscReal>   lcoor;
                getIndexList(path2index.c_str(), kIndex, nkSec);

                lcoor.resize(nkSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->kSections));
                PetscMalloc(sizeof(PetscInt) * nkSec,  &(acquisition->kSections->indices));
                PetscMalloc(sizeof(PetscReal) * nkSec,  &(acquisition->kSections->coordinates));

                sections *kSections = acquisition->kSections;

                for(PetscInt s = 0; s< nkSec; s++)
                {
                  lcoor[s] = 0.0;
                  kSections->coordinates[s] = 0.0;
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

                if(flags.isTeqnActive)
                {
                  atLeastOneScalar++;
                }

                if(flags.isLesActive)
                {
                  atLeastOneScalar++;
                }

                if(acquisition->isPerturbABLActive)
                {
                    atLeastOneScalar++;
                    atLeastOneVector++;
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
                word path2index = "./postProcessing/" + mesh->meshName + "/jSurfaces/";
                std::vector<PetscInt>      jIndex;
                PetscInt                   njSec;
                std::vector<PetscReal>   lcoor;
                getIndexList(path2index.c_str(), jIndex, njSec);

                lcoor.resize(njSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->jSections));
                PetscMalloc(sizeof(PetscInt) * njSec,  &(acquisition->jSections->indices));
                PetscMalloc(sizeof(PetscReal) * njSec,  &(acquisition->jSections->coordinates));

                sections *jSections = acquisition->jSections;

                for(PetscInt s = 0; s< njSec; s++)
                {
                  lcoor[s] = 0.0;
                  jSections->coordinates[s] = 0.0;
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

                if(flags.isTeqnActive)
                {
                  atLeastOneScalar++;
                }

                if(flags.isLesActive)
                {
                  atLeastOneScalar++;
                }

                if(acquisition->isPerturbABLActive)
                {
                    atLeastOneScalar++;
                    atLeastOneVector++;
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
                word path2index = "./postProcessing/" + mesh->meshName + "/iSurfaces/";
                std::vector<PetscInt>      iIndex;
                PetscInt                   niSec;
                std::vector<PetscReal>   lcoor;
                getIndexList(path2index.c_str(), iIndex, niSec);

                lcoor.resize(niSec);

                // allocate memory
                PetscMalloc(sizeof(sections), &(acquisition->iSections));
                PetscMalloc(sizeof(PetscInt) * niSec,  &(acquisition->iSections->indices));
                PetscMalloc(sizeof(PetscReal) * niSec,  &(acquisition->iSections->coordinates));

                sections *iSections = acquisition->iSections;

                for(PetscInt s = 0; s< niSec; s++)
                {
                  lcoor[s] = 0.0;
                  iSections->coordinates[s] = 0.0;
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

                if(flags.isTeqnActive)
                {
                  atLeastOneScalar++;
                }

                if(flags.isLesActive)
                {
                  atLeastOneScalar++;
                }

                if(acquisition->isPerturbABLActive)
                {
                    atLeastOneScalar++;
                    atLeastOneVector++;
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

            if(d == 0)
            {
              PetscPrintf(mesh->MESH_COMM, "done\n\n");
            }

            MPI_Barrier(mesh->MESH_COMM);

        }

    }

  return(0);
}
