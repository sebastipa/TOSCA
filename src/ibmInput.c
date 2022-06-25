//! \file  ibmInput.c
//! \brief Contains Immersed boundary method input and read function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/ibmInput.h"

//***************************************************************************************************************//

PetscErrorCode readIBMProperties(ibm_ *ibm)
{
  PetscMPIInt rank;

  mesh_ *mesh = ibm->access->mesh;
  io_ *io = ibm->access->io;

  MPI_Comm_rank(mesh->MESH_COMM, &rank);

  PetscPrintf(PETSC_COMM_WORLD, "\nIBM initialization...\n");

  // read debug switch
  readDictInt("./IBM/IBMProperties.dat", "debug", &(ibm->dbg));

  // read dynamic IBM switch
  readDictInt("./IBM/IBMProperties.dat", "dynamic", &(ibm->dynamic));

  // read compute force and moment
  readDictInt("./IBM/IBMProperties.dat", "computeForce", &(ibm->computeForce));

  // read check normals
  readDictInt("./IBM/IBMProperties.dat", "checkNormal", &(ibm->checkNormal));

  // write stl flag
  readDictInt("./IBM/IBMProperties.dat", "writeSTL", &(ibm->writeSTL));

  // read the number of ibm bodies
  readDictInt("./IBM/IBMProperties.dat", "NumberofBodies", &(ibm->numBodies));

  // read the interpolation method - MLS or CURVIB
  readDictWord("./IBM/IBMProperties.dat", "InterpolationMethod", &(ibm->IBInterpolationModel));

  // read the write settings
  readSubDictDouble("./IBM/IBMProperties.dat","writeSettings","timeStart", &(ibm->timeStart));
  readSubDictWord  ("./IBM/IBMProperties.dat","writeSettings","intervalType", &(ibm->intervalType));
  readSubDictDouble("./IBM/IBMProperties.dat","writeSettings","timeInterval", &(ibm->timeInterval));

  //counter for number of moving objects
  int movingObject = 0;

  // initialize pointers to NULL
  ibm->ibmFCells = NULL;
  ibm->sBox      = NULL;

  // allocate memory for each ibm object
  ibm->ibmBody = new ibmObject*[ibm->numBodies];

  // allocate memory for the search box Parameters
  ibm->sBox = new searchBox[ibm->numBodies];

  for (PetscInt i=0; i < ibm->numBodies; i++)
  {
    ibm->ibmBody[i] = new ibmObject;

    ibmObject   *ibmBody  = ibm->ibmBody[i];

    // set pointers to null
    ibmBody->bound          = NULL;
    ibmBody->searchCellList = NULL;
    ibmBody->ibMsh          = NULL;
    ibmBody->ibmRot         = NULL;

    // allocate memory for the IBM mesh of the object
    ibmBody->ibMsh = new ibmMesh;

    // allocate memory for the bounding box  of the object
    ibmBody->bound = new boundingBox;

    char objectName[256];
    sprintf(objectName, "object%ld", i);

    // read object name
    readSubDictWord("./IBM/IBMProperties.dat", objectName, "bodyName", &(ibmBody->bodyName));

    // read object type - surfaceBody or closedBody
    readSubDictWord("./IBM/IBMProperties.dat", objectName, "bodyType", &(ibmBody->bodyType));

    // read the body base location
    readSubDictVector("./IBM/IBMProperties.dat", objectName, "baseLocation", &(ibmBody->baseLocation));

    if(ibmBody->bodyType == "surfaceBody")
    {
        //read number of surfaces
        readSubDictInt("./IBM/IBMProperties.dat", objectName, "numSurfaces", &(ibmBody->numSurfaces));

        PetscInt numSurf = ibmBody->numSurfaces;

        // allocate memory for each surface object
        ibmBody->ibmSurface = new surface*[numSurf];

        // read ibm file type
        std::string* surfaceFileType = readSubDictWordArray("./IBM/IBMProperties.dat", objectName, "fileType", numSurf);

        //read the surfaces array
        std::string* surfaceNames = readSubDictWordArray("./IBM/IBMProperties.dat", objectName, "surfaceNames", numSurf);

        // loop through the surfaces and allocate memory
        for(PetscInt s = 0; s < numSurf; s++)
        {
            ibmBody->ibmSurface[s] = new surface;

            surface *ibmSurface = ibmBody->ibmSurface[s];

            ibmSurface->surfaceName = surfaceNames[s];
            ibmSurface->surfaceFileType = surfaceFileType[s];

            ibmSurface->surfaceId = s;

            ibmSurface->baseLocation = nSet(ibmBody->baseLocation);

            if(ibmSurface->surfaceFileType == "inp")
            {
                readSubDictWord("./IBM/IBMProperties.dat", objectName, "elementSet", &(ibmBody->elementSet));
                ibmSurface->elementSet = ibmBody->elementSet;
            }

            //allocate the mesh memory
            ibmSurface->ibMsh = new ibmMesh;
        }
    }
    else
    {
        // read ibm file type
        readSubDictWord("./IBM/IBMProperties.dat", objectName, "fileType", &(ibmBody->fileType));

        if(ibmBody->fileType == "inp")
        {
            // read name of the element set
            readSubDictWord("./IBM/IBMProperties.dat", objectName, "elementSet", &(ibmBody->elementSet));
        }

    }

    // read the ibm motion
    readSubDictWord("./IBM/IBMProperties.dat", objectName, "bodyMotion", &(ibmBody->bodyMotion));

    if(ibmBody->bodyMotion != "static")
    {
        movingObject ++;
    }

    if(ibm->dynamic)
    {
        if(ibmBody->bodyMotion == "rotation")
        {
            // allocate memory for ibm rotation
            ibmBody->ibmRot = new ibmRotation;

            ibmRotation *ibmRot = ibmBody->ibmRot;

            // set the initial rotation angle
            ibmRot->rotAngle = 0;

            readSubDictDouble("./IBM/IBMProperties.dat", objectName, "angularSpeed", &(ibmRot->angSpeed));
            readSubDictDouble("./IBM/IBMProperties.dat", objectName, "angularAcceleration", &(ibmRot->angAcc));
            readSubDictVector("./IBM/IBMProperties.dat", objectName, "rotationAxis", &(ibmRot->rotAxis));
            readSubDictVector("./IBM/IBMProperties.dat", objectName, "rotationCenter", &(ibmRot->rotCenter));
            readSubDictDouble("./IBM/IBMProperties.dat", objectName, "maxTipRadius", &(ibmRot->maxR));

            //transform angular speed and acceleration from rpm to rad/s
            ibmRot->angSpeed = ibmRot->angSpeed*2*M_PI/60.0;
            ibmRot->angAcc   = ibmRot->angAcc*2*M_PI/60.0;
        }
    }

    // read the search cell ratio wrt to the average cell size of the domain mesh
    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "searchCellRatio", &(ibmBody->searchCellRatio));

    //set the body index as the index of the loop
    ibm->ibmBody[i]->bodyID = i;

  }

  if(movingObject == 0 && ibm->dynamic == 1)
  {
     char error[512];
      sprintf(error, "ibm dynamic motion set to true but there are no moving objects. Set dynamic to 0 in IBMProperties\n");
      fatalErrorInFunction("readIBMProperties",  error);
  }

  if(movingObject > 0 && ibm->dynamic == 0)
  {
     char error[512];
      sprintf(error, "ibm dynamic motion set to false but there are moving objects. Set dynamic to 1 in IBMProperties\n");
      fatalErrorInFunction("readIBMProperties",  error);
  }

  for (PetscInt i=0; i < ibm->numBodies; i++)
  {
      // read the mesh file for the ibm object
      readIBMObjectMesh(ibm, i);
  }

  return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMObjectMesh(ibm_ *ibm, PetscInt b)
{
  ibmObject     *ibmBody = ibm->ibmBody[b];
  clock_        *clock   = ibm->access->clock;
  ibmMesh       *ibMesh  = ibmBody->ibMsh;

  if(ibmBody->bodyType == "surfaceBody")
  {
      //loop through the number of surfaces
      for (PetscInt i = 0; i < ibmBody->numSurfaces; i++)
      {
          surface *ibmSurface = ibmBody->ibmSurface[i];

          if(ibmSurface->surfaceFileType == "ucd")
          {
              readIBMSurfaceFileUCD(ibmSurface);
          }
          else if(ibmSurface->surfaceFileType == "ascii")
          {
              readIBMSurfaceFileASCIIRaster(ibmSurface);
          }
          else if(ibmSurface->surfaceFileType == "stl")
          {
              readIBMSurfaceFileSTL(ibmSurface);
          }
          else if(ibmSurface->surfaceFileType == "inp")
          {
              readIBMSurfaceFileAbaqusInp(ibmSurface);
          }
          else
          {
              char error[530];
              sprintf(error, "wrong ibm surface file type. Use ucd, ascii, inp or stl\n");
              fatalErrorInFunction("readIBMObjectMesh",  error);
          }
      }

      combineMesh(ibmBody);

  }
  else
  {
      if(ibmBody->fileType == "ucd")
      {
          readIBMBodyFileUCD(ibmBody);
      }
      else if(ibmBody->fileType == "ascii")
      {
          readIBMBodyFileASCIIRaster(ibmBody);
      }
      else if(ibmBody->fileType == "stl")
      {
          readIBMBodyFileSTL(ibmBody);
      }
      else if(ibmBody->fileType == "inp")
      {
          readIBMBodyFileAbaqusInp(ibmBody);
      }
      else
      {
          char error[530];
          sprintf(error, "wrong ibm file type. Use ucd, ascii, inp or stl\n");
          fatalErrorInFunction("readIBMObjectMesh",  error);
      }

  }

  // reset the co-ordinates if the angular position is read
  if(ibmBody->bodyMotion == "rotation")
  {
      word        timeName, dictName;
      ibmRotation   *ibmRot  = ibmBody->ibmRot;

      //read the angular position from file
      timeName = "./fields/" + ibm->access->mesh->meshName + "/ibm/" + getTimeName(clock);

      // check if this folder exists, if yes angular position is set, if no, angular position is 0
      if(dir_exist(timeName.c_str()))
      {
          PetscPrintf(ibm->access->mesh->MESH_COMM, "   reading IBM data for time: %s \n",timeName.c_str());

          dictName = timeName + "/" + ibmBody->bodyName;

          readDictDouble(dictName.c_str(), "AngularPosition", &(ibmRot->rotAngle));
          readDictDouble(dictName.c_str(), "AngularSpeed", &(ibmRot->angSpeed));
          readDictDouble(dictName.c_str(), "AngularAcceleration", &(ibmRot->angAcc));

          //transform angular speed and acceleration from rpm to rad/s
          ibmRot->angSpeed = ibmRot->angSpeed*2*M_PI/60.0;
          ibmRot->angAcc   = ibmRot->angAcc*2*M_PI/60.0;
      }
      else
      {
          ibmRot->rotAngle = 0.0;
      }

      for(PetscInt i = 0; i < ibMesh->nodes; i++)
      {
          Cmpnts rvec   = nSub(ibMesh->nCoor[i], ibmRot->rotCenter);
          Cmpnts angVel = nScale(ibmRot->angSpeed, ibmRot->rotAxis);

          // rotate the vector to the required angular position
          mRot(ibmRot->rotAxis, rvec, ibmRot->rotAngle);

          // find the new co-ordinate
          ibMesh->nCoor[i] = nSum(rvec, ibmRot->rotCenter);

          ibMesh->nU[i]  = nCross(angVel, rvec);
      }
  }

  // allocate memory for the element normal, area and center coordinate
  PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eN));
  PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eT1));
  PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eT2));

  PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->flipNormal));

  for (PetscInt i=0; i<ibMesh->elems; i++)
  {
      ibMesh->flipNormal[i] = 0;
  }

  PetscMalloc(ibMesh->elems * sizeof(PetscReal), &(ibMesh->eA));
  PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eCent));


  if(ibm->computeForce)
  {
      // allocate memory for pressure and surface stress
      PetscMalloc( ibMesh->elems * sizeof(Cmpnts), &(ibmBody->ibmPForce));

      // allocate memory for the closest cell id to the ibm mesh element
      PetscMalloc( ibMesh->elems * sizeof(cellIds), &(ibmBody->closestCells));
      PetscMalloc( ibMesh->elems * sizeof(PetscInt), &(ibmBody->thisPtControlled));
      PetscMalloc( ibMesh->elems * sizeof(PetscInt), &(ibmBody->thisPtControlTransfer));
  }

  //allocate memory for the smallest bounding sphere for each element - used for finding the nearest IBM Mesh element to IBM fluid cell
  PetscMalloc(ibMesh->elems * sizeof(PetscReal), &(ibMesh->eRVec));
  PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eQVec));

  return 0;
}

//***************************************************************************************************************//
PetscErrorCode combineMesh(ibmObject *ibmBody)
{
    // pointer to the mesh object of the current body
    ibmMesh *ibMesh = ibmBody->ibMsh;

    // total elements and nodes
    PetscInt elemTot = 0, nodeTot = 0;
    for (PetscInt s = 0; s < ibmBody->numSurfaces; s++)
    {
        surface *ibmSurface = ibmBody->ibmSurface[s];
        elemTot += ibmSurface->ibMsh->elems;
        nodeTot += ibmSurface->ibMsh->nodes;
    }

    ibMesh->nodes = nodeTot;
    ibMesh->elems = elemTot;

    // allocate memory for the x, y and z co-ordinates of the nodes
    ibMesh->nCoor = new Cmpnts[nodeTot];
    ibMesh->nU = new Cmpnts[nodeTot];

    PetscInt ctr = 0;
    for (PetscInt s = 0; s < ibmBody->numSurfaces; s++)
    {
        surface *ibmSurface = ibmBody->ibmSurface[s];

        for(PetscInt n=0; n < ibmSurface->ibMsh->nodes; n++)
        {
            ibMesh->nCoor[ctr] = nSet(ibmSurface->ibMsh->nCoor[n]);

            // initialize node velocity to 0
            mSetValue(ibMesh->nU[ctr], 0.0);

            // PetscPrintf(PETSC_COMM_WORLD, "node = %ld %ld, node coord = %lf %lf %lf, locally %lf %lf %lf\n", ctr, n, ibMesh->nCoor[ctr].x, ibMesh->nCoor[ctr].y, ibMesh->nCoor[ctr].z, ibmSurface->ibMsh->nCoor[n].x, ibmSurface->ibMsh->nCoor[n].y, ibmSurface->ibMsh->nCoor[n].z);

            ctr++;
        }
    }

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    ibMesh->nID1 = new PetscInt[elemTot];
    ibMesh->nID2 = new PetscInt[elemTot];
    ibMesh->nID3 = new PetscInt[elemTot];

    ctr = 0;
    PetscInt nTot = 0;  //total nodes in a surface body
    PetscInt eTot = 0;
    //allocate memory for the pointer to the surface body each element belongs to
    ibMesh->eSurface = new PetscInt[elemTot];

    //allocate memory for element mapping for .inp files
    ibmBody->elementMapping = new PetscInt[elemTot];

    for (PetscInt s = 0; s < ibmBody->numSurfaces; s++)
    {
        surface *ibmSurface = ibmBody->ibmSurface[s];

        for(PetscInt e=0; e < ibmSurface->ibMsh->elems; e++)
        {
            ibMesh->nID1[ctr] = ibmSurface->ibMsh->nID1[e] + nTot;
            ibMesh->nID2[ctr] = ibmSurface->ibMsh->nID2[e] + nTot;
            ibMesh->nID3[ctr] = ibmSurface->ibMsh->nID3[e] + nTot;

            ibMesh->eSurface[ctr] = s;

            if(ibmSurface->surfaceFileType == "inp")
            {
                ibmBody->elementMapping[ctr] = ibmSurface->elementMapping[e] + eTot;
            }
            else
            {
                ibmBody->elementMapping[ctr] = eTot;
            }
            ctr++;
        }

        nTot += ibmSurface->ibMsh->nodes;
        eTot += ibmSurface->ibMsh->elems;
    }

    return 0;
}
//***************************************************************************************************************//
PetscErrorCode readIBMBodyFileSTL(ibmObject *ibmBody)
{
    char      ibmFile[256];
    FILE      *fd;
    PetscInt  readError; char* charError;
    char      string[128];
    ibmMesh   *ibMesh = ibmBody->ibMsh;

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmBody->bodyName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmBody->bodyName.c_str());

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[512];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMSurfaceFileSTL(surface *ibmSurface)
{
    char      ibmFile[256];
    FILE      *fd;
    PetscInt  readError; char* charError;
    char      string[128];
    ibmMesh   *ibMesh = ibmSurface->ibMsh;

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmSurface->surfaceName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmSurface->surfaceName.c_str());

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[512];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMBodyFileAbaqusInp(ibmObject *ibmBody)
{
    // file stream
    std::ifstream indata;

    // ibm file
    char ibmFile[256];

    // word by word read
    char word[256];

    // pointer for strtod
    char *eptr;

    // total number of elements in the prescribed surface
    PetscInt totalElements = 0;

    // name of the element set to be read from the abaqus inp file
    std::string elementSet;

    // pointer to the mesh object of the current body
    ibmMesh   *ibMesh = ibmBody->ibMsh;

    // arrays to store the ibm mesh lists
    std::vector<std::vector<PetscInt>> elementStartList;

    std::vector<std::vector<PetscInt>> elemConnectivityList;

    std::vector<std::vector<PetscReal>> nodeList;

    // initialize the arrays
    elementStartList.resize(2);

    elemConnectivityList.resize(4);

    nodeList.resize(4);

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmBody->bodyName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmBody->bodyName.c_str());

    // open the ibm object file in read mode
    indata.open(ibmFile);

    // create the element start list for the given element surface
    if(!indata)
    {
       char error[512];
       sprintf(error, "could not open ibm file %s \n", ibmFile);
       fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *elset
            if
            (
                strcmp
                (
                    "*elset,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);
                std::string token2;

                token2 = token1.substr(6, token1.length() - 6);

                // test if the elementSet matches with the prescribed element set from IBMProperties.dat
                if
                (
                    strcmp
                    (
                        ibmBody->elementSet.c_str(),
                        token2.c_str()
                    ) == 0
                )
                {
                    PetscInt ctr = 0, eID, eIDLast = -1;
                    std::string elementID;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *elset or end of file
                        if
                        (
                            strcmp
                            (
                                "*elset,",
                                word
                            ) == 0
                        )
                        {
                            break;
                        }
                        else
                        {
                            std::string temp1(word);

                            // remove ending , from the read element
                            elementID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(elementID))
                            {

                                eID = std::strtol(elementID.c_str(), &eptr, 10);

                                if(fabs(eID-eIDLast) > 1)
                                {
                                    //save the start element
                                    elementStartList[0].push_back(eID);

                                    if(eIDLast != -1)
                                    {
                                        //save number of consecutive elements
                                        elementStartList[1].push_back(ctr);
                                    }

                                    //reset ctr to 1 and save last element
                                    ctr = 1;
                                    eIDLast = eID;
                                }
                                else
                                {
                                    ctr ++;
                                    eIDLast = eID;
                                }
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                        }

                    }

                    //save number of consecutive elements for the final set
                    elementStartList[1].push_back(ctr);
                }
            }
        }

    }

    indata.close();

    // open the ibm object file in read mode to read the entire element connectivity list
    indata.open(ibmFile);

    if(!indata)
    {
        char error[512];
        sprintf(error, "could not open ibm file %s \n", ibmFile);
        fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *element,
            if
            (
                strcmp
                (
                    "*element,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);

                if
                (
                    strcmp
                    (
                        "elset=S3,type=S3",
                        token1.c_str()
                    ) == 0
                )
                {
                    PetscInt eID, nID;
                    std::string elementID, nodeID;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *element or *elset or end of file
                        if
                        (
                            (   strcmp
                                (
                                    "*elset,",
                                    word
                                ) == 0
                            )

                            ||

                            (   strcmp
                                (
                                    "*element,",
                                    word
                                ) == 0
                            )
                        )
                        {
                            break;
                        }
                        else
                        {
                            //read the element number
                            std::string temp1(word);

                            // remove ending , from the read element
                            elementID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(elementID))
                            {
                                eID = std::strtol(elementID.c_str(), &eptr, 10);

                                elemConnectivityList[0].push_back(eID);
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                            //read the connectivity nodes
                            for (PetscInt i = 0; i < 3; i++)
                            {
                                indata >> word;

                                //read the node number
                                std::string temp2(word);

                                // remove ending , from the read element
                                nodeID = temp2.substr(0, temp2.length() - 1);

                                // check if nodeID contains character (throws errors if yes)
                                if(isNumber(nodeID))
                                {
                                    nID = std::strtol(nodeID.c_str(), &eptr, 10);

                                    elemConnectivityList[i+1].push_back(nID);
                                }
                                else
                                {
                                    char error[512];
                                     sprintf(error, "expected number in elementList\n");
                                     fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    indata.close();

    //read the nodes from the ibmfile
    indata.open(ibmFile);

    if(!indata)
    {
        char error[512];
        sprintf(error, "could not open ibm file %s \n", ibmFile);
        fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *element,
            if
            (
                strcmp
                (
                    "*node,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);

                if
                (
                    strcmp
                    (
                        "nset=Nall",
                        token1.c_str()
                    ) == 0
                )
                {
                    PetscInt nID;
                    std::string nodeID, nodeCoord;
                    PetscReal nCoord;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *element or *elset or end of file
                        if
                        (
                            strcmp
                            (
                                "*nset,",
                                word
                            ) == 0
                        )
                        {
                            break;
                        }
                        else
                        {
                            //read the node number
                            std::string temp1(word);

                            // remove ending , from the read element
                            nodeID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(nodeID))
                            {
                                nID = std::strtol(nodeID.c_str(), &eptr, 10);

                                nodeList[0].push_back((PetscReal)nID);
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                            //read the node coordinates - z coordinate - doesnt have comma
                            for (PetscInt i = 0; i < 3; i++)
                            {
                                indata >> word;

                                //read the node number
                                std::string temp2(word);

                                // remove ending , from the read element
                                if( (i == 0) || (i == 1))
                                {
                                    nodeCoord = temp2.substr(0, temp2.length() - 1);
                                }
                                else
                                {
                                    nodeCoord = temp2;
                                }

                                //convert to double
                                nCoord = std::strtod(nodeCoord.c_str(), &eptr);

                                nodeList[i+1].push_back(nCoord);

                            }

                        }
                    }
                }
            }
        }

    }

    //find total number of elements in the ibm body
    for (PetscInt i = 0; i < elementStartList[0].size(); i++)
    {
        totalElements += elementStartList[1][i];
    }

    ibMesh->nodes = nodeList[0].size();
    ibMesh->elems = totalElements;

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // read the node co-rdinates from the file and initialize the node velocity to 0
    for (PetscInt i = 0; i < nodeList[0].size(); i++)
    {
        ibMesh->nCoor[i].x = nodeList[1][i];
        ibMesh->nCoor[i].y = nodeList[2][i];
        ibMesh->nCoor[i].z = nodeList[3][i];

        // translate the body based on the base location
        mSum(ibMesh->nCoor[i], ibmBody->baseLocation);

        // initialize node velocity to 0
        mSetValue(ibMesh->nU[i], 0.0);
    }

    //mapping from the ibm element list to the global element list
    std::vector<PetscInt> local2GlobalConnectivity;
    local2GlobalConnectivity.resize(totalElements);

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID3));

    //local element counter
    PetscInt eCtr = 0;

    // create the IBM element list
    for (PetscInt i = 0; i < elemConnectivityList[0].size(); i++)
    {   //loop through the entire element list
        for (PetscInt j = 0; j < elementStartList[0].size(); j++)
        {
            if(elemConnectivityList[0][i] == elementStartList[0][j])
            {   // loop through the consequtive elements from the elementStartList
                for (PetscInt k = 0; k < elementStartList[1][j]; k++)
                {
                    local2GlobalConnectivity[eCtr] = elemConnectivityList[0][i+k];
                    ibMesh->nID1[eCtr] = elemConnectivityList[1][i+k];
                    ibMesh->nID2[eCtr] = elemConnectivityList[2][i+k];
                    ibMesh->nID3[eCtr] = elemConnectivityList[3][i+k];
                    eCtr ++;
                }
            }
        }
    }

    for (PetscInt i = 0; i < ibMesh->elems; i++)
    {
        // node numbers start from 0
        ibMesh->nID1[i] = ibMesh->nID1[i] - 1;
        ibMesh->nID2[i] = ibMesh->nID2[i] - 1;
        ibMesh->nID3[i] = ibMesh->nID3[i] - 1;
    }

    //save the element mapping vector
    PetscMalloc( ibMesh->elems * sizeof(PetscInt), &(ibmBody->elementMapping));

    for(PetscInt i = 0; i < ibMesh->elems; i++)
    {
        ibmBody->elementMapping[i] = local2GlobalConnectivity[i];
    }

    //clear temporary vector memory
    for (PetscInt i = 0; i < 2; i++)
    {
        std::vector<PetscInt> ().swap(elementStartList[i]);
    }

    for (PetscInt i = 0; i < 4; i++)
    {
        std::vector<PetscInt> ().swap(elemConnectivityList[i]);
        std::vector<PetscReal> ().swap(nodeList[i]);
    }

    std::vector<PetscInt> ().swap(local2GlobalConnectivity);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMSurfaceFileAbaqusInp(surface *ibmSurface)
{
    // file stream
    std::ifstream indata;

    // ibm file
    char ibmFile[256];

    // word by word read
    char word[256];

    // pointer for strtod
    char *eptr;

    // total number of elements in the prescribed surface
    PetscInt totalElements = 0;

    // pointer to the mesh object of the current body
    ibmMesh   *ibMesh = ibmSurface->ibMsh;

    // arrays to store the ibm mesh lists
    std::vector<std::vector<PetscInt>> elementStartList;

    std::vector<std::vector<PetscInt>> elemConnectivityList;

    std::vector<std::vector<PetscReal>> nodeList;

    // initialize the arrays
    elementStartList.resize(2);

    elemConnectivityList.resize(4);

    nodeList.resize(4);

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmSurface->surfaceName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmSurface->surfaceName.c_str());

    // open the ibm object file in read mode
    indata.open(ibmFile);

    // create the element start list for the given element surface
    if(!indata)
    {
       char error[512];
       sprintf(error, "could not open ibm file %s \n", ibmFile);
       fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *elset
            if
            (
                strcmp
                (
                    "*elset,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);
                std::string token2;

                token2 = token1.substr(6, token1.length() - 6);

                // test if the elementSet matches with the prescribed element set from IBMProperties.dat
                if
                (
                    strcmp
                    (
                        ibmSurface->elementSet.c_str(),
                        token2.c_str()
                    ) == 0
                )
                {
                    PetscInt ctr = 0, eID, eIDLast = -1;
                    std::string elementID;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *elset or end of file
                        if
                        (
                            strcmp
                            (
                                "*elset,",
                                word
                            ) == 0
                        )
                        {
                            break;
                        }
                        else
                        {
                            std::string temp1(word);

                            // remove ending , from the read element
                            elementID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(elementID))
                            {

                                eID = std::strtol(elementID.c_str(), &eptr, 10);

                                if(fabs(eID-eIDLast) > 1)
                                {
                                    //save the start element
                                    elementStartList[0].push_back(eID);

                                    if(eIDLast != -1)
                                    {
                                        //save number of consecutive elements
                                        elementStartList[1].push_back(ctr);
                                    }

                                    //reset ctr to 1 and save last element
                                    ctr = 1;
                                    eIDLast = eID;
                                }
                                else
                                {
                                    ctr ++;
                                    eIDLast = eID;
                                }
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                        }

                    }

                    //save number of consecutive elements for the final set
                    elementStartList[1].push_back(ctr);
                }
            }
        }

    }

    indata.close();

    // open the ibm object file in read mode to read the entire element connectivity list
    indata.open(ibmFile);

    if(!indata)
    {
        char error[512];
        sprintf(error, "could not open ibm file %s \n", ibmFile);
        fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *element,
            if
            (
                strcmp
                (
                    "*element,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);

                if
                (
                    strcmp
                    (
                        "elset=S3,type=S3",
                        token1.c_str()
                    ) == 0
                )
                {
                    PetscInt eID, nID;
                    std::string elementID, nodeID;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *element or *elset or end of file
                        if
                        (
                            (   strcmp
                                (
                                    "*elset,",
                                    word
                                ) == 0
                            )

                            ||

                            (   strcmp
                                (
                                    "*element,",
                                    word
                                ) == 0
                            )
                        )
                        {
                            break;
                        }
                        else
                        {
                            //read the element number
                            std::string temp1(word);

                            // remove ending , from the read element
                            elementID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(elementID))
                            {
                                eID = std::strtol(elementID.c_str(), &eptr, 10);

                                elemConnectivityList[0].push_back(eID);
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                            //read the connectivity nodes
                            for (PetscInt i = 0; i < 3; i++)
                            {
                                indata >> word;

                                //read the node number
                                std::string temp2(word);

                                // remove ending , from the read element
                                nodeID = temp2.substr(0, temp2.length() - 1);

                                // check if nodeID contains character (throws errors if yes)
                                if(isNumber(nodeID))
                                {
                                    nID = std::strtol(nodeID.c_str(), &eptr, 10);

                                    elemConnectivityList[i+1].push_back(nID);
                                }
                                else
                                {
                                    char error[512];
                                     sprintf(error, "expected number in elementList\n");
                                     fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    indata.close();

    //read the nodes from the ibmfile
    indata.open(ibmFile);

    if(!indata)
    {
        char error[512];
        sprintf(error, "could not open ibm file %s \n", ibmFile);
        fatalErrorInFunction("readIBMFileINP",  error);
    }
    else
    {
        // get word by word till end of ibmFile
        while(!indata.eof())
        {
            indata >> word;

            // test if found element set keyword *element,
            if
            (
                strcmp
                (
                    "*node,",
                    word
                ) == 0
            )
            {
                indata >> word;

                std::string token1(word);

                if
                (
                    strcmp
                    (
                        "nset=Nall",
                        token1.c_str()
                    ) == 0
                )
                {
                    PetscInt nID;
                    std::string nodeID, nodeCoord;
                    PetscReal nCoord;

                    while(!indata.eof())
                    {
                        indata >> word;

                        //read the elemet ids until encounters the next keyword *element or *elset or end of file
                        if
                        (
                            strcmp
                            (
                                "*nset,",
                                word
                            ) == 0
                        )
                        {
                            break;
                        }
                        else
                        {
                            //read the node number
                            std::string temp1(word);

                            // remove ending , from the read element
                            nodeID = temp1.substr(0, temp1.length() - 1);

                            // check if elementID contains character (throws errors if yes)
                            if(isNumber(nodeID))
                            {
                                nID = std::strtol(nodeID.c_str(), &eptr, 10);

                                nodeList[0].push_back((PetscReal)nID);
                            }
                            else
                            {
                                char error[512];
                                 sprintf(error, "expected number in elementList\n");
                                 fatalErrorInFunction("readIBMFileAbaqusInp",  error);
                            }

                            //read the node coordinates - z coordinate - doesnt have comma
                            for (PetscInt i = 0; i < 3; i++)
                            {
                                indata >> word;

                                //read the node number
                                std::string temp2(word);

                                // remove ending , from the read element
                                if( (i == 0) || (i == 1))
                                {
                                    nodeCoord = temp2.substr(0, temp2.length() - 1);
                                }
                                else
                                {
                                    nodeCoord = temp2;
                                }

                                //convert to double
                                nCoord = std::strtod(nodeCoord.c_str(), &eptr);

                                nodeList[i+1].push_back(nCoord);

                            }

                        }
                    }
                }
            }
        }

    }

    //find total number of elements in the ibm body
    for (PetscInt i = 0; i < elementStartList[0].size(); i++)
    {
        totalElements += elementStartList[1][i];
    }

    ibMesh->nodes = nodeList[0].size();
    ibMesh->elems = totalElements;

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // read the node co-rdinates from the file and initialize the node velocity to 0
    for (PetscInt i = 0; i < nodeList[0].size(); i++)
    {
        ibMesh->nCoor[i].x = nodeList[1][i];
        ibMesh->nCoor[i].y = nodeList[2][i];
        ibMesh->nCoor[i].z = nodeList[3][i];

        // translate the body based on the base location
        mSum(ibMesh->nCoor[i], ibmSurface->baseLocation);

        // initialize node velocity to 0
        mSetValue(ibMesh->nU[i], 0.0);
    }

    //mapping from the ibm element list to the global element list
    std::vector<PetscInt> local2GlobalConnectivity;
    local2GlobalConnectivity.resize(totalElements);

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(totalElements * sizeof(PetscInt), &(ibMesh->nID3));

    //local element counter
    PetscInt eCtr = 0;

    // create the IBM element list
    for (PetscInt i = 0; i < elemConnectivityList[0].size(); i++)
    {   //loop through the entire element list
        for (PetscInt j = 0; j < elementStartList[0].size(); j++)
        {
            if(elemConnectivityList[0][i] == elementStartList[0][j])
            {   // loop through the consequtive elements from the elementStartList
                for (PetscInt k = 0; k < elementStartList[1][j]; k++)
                {
                    local2GlobalConnectivity[eCtr] = elemConnectivityList[0][i+k];
                    ibMesh->nID1[eCtr] = elemConnectivityList[1][i+k];
                    ibMesh->nID2[eCtr] = elemConnectivityList[2][i+k];
                    ibMesh->nID3[eCtr] = elemConnectivityList[3][i+k];
                    eCtr ++;
                }
            }
        }
    }

    for (PetscInt i = 0; i < ibMesh->elems; i++)
    {
        // node numbers start from 0
        ibMesh->nID1[i] = ibMesh->nID1[i] - 1;
        ibMesh->nID2[i] = ibMesh->nID2[i] - 1;
        ibMesh->nID3[i] = ibMesh->nID3[i] - 1;
    }

    //save the element mapping vector
    PetscMalloc( ibMesh->elems * sizeof(PetscInt), &(ibmSurface->elementMapping));

    for(PetscInt i = 0; i < ibMesh->elems; i++)
    {
        ibmSurface->elementMapping[i] = local2GlobalConnectivity[i];
    }

    //clear temporary vector memory
    for (PetscInt i = 0; i < 2; i++)
    {
        std::vector<PetscInt> ().swap(elementStartList[i]);
    }

    for (PetscInt i = 0; i < 4; i++)
    {
        std::vector<PetscInt> ().swap(elemConnectivityList[i]);
        std::vector<PetscReal> ().swap(nodeList[i]);
    }

    std::vector<PetscInt> ().swap(local2GlobalConnectivity);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMBodyFileASCIIRaster(ibmObject *ibmBody)
{
    char      ibmFile[256];
    FILE      *fd;
    PetscInt  readError; char* charError;
    char      string[128];
    ibmMesh   *ibMesh = ibmBody->ibMsh;

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmBody->bodyName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmBody->bodyName.c_str());

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[512];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    //skip the first line
    charError = fgets(string, 128, fd);

    PetscInt    numX, numY;
    PetscReal   xMax, xMin, yMax, yMin, zMax, zMin;
    PetscReal   delX, delY;
    readError = fscanf(fd, "%ld %ld",  &numX, &numY);

    readError = fscanf(fd, "%le %le", &xMin, &xMax);
    readError = fscanf(fd, "%le %le", &yMin, &yMax);
    readError = fscanf(fd, "%le %le", &zMin, &zMax);

    delX = (xMax - xMin)/ (numX-1);
    delY = (yMax - yMin)/ (numY-1);

    ibMesh->nodes = numX * numY;
    ibMesh->elems = 2 * (numX-1) * (numY-1);

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // set the node co-rdinates from the file and initialize the node velocity to 0
    PetscInt n = 0;
    for(PetscInt j = 0; j < numY; j++)
    {
        for(PetscInt i = 0; i < numX; i++)
        {
            ibMesh->nCoor[n].x = xMin + i*delX;

            ibMesh->nCoor[n].y = yMin + j* delY;

            readError = fscanf(fd, "%le", &(ibMesh->nCoor[n].z));

            // translate the body based on the base location
            mSum(ibMesh->nCoor[n], ibmBody->baseLocation);

            // initialize node velocity to 0
            mSetValue(ibMesh->nU[n], 0.0);
            n++;
        }
    }

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID3));

    n = 0;
    for(PetscInt j = 0; j < numY-1; j++)
    {
        for(PetscInt i = 0; i < numX-1; i++)
        {
            ibMesh->nID1[n] = numX*j + i;
            ibMesh->nID2[n] = numX*(j) + i+1;
            ibMesh->nID3[n] = numX*(j+1) + i;

            n++;

            ibMesh->nID1[n] = numX*(j) + i+1;
            ibMesh->nID2[n] = numX*(j+1) + i+1;
            ibMesh->nID3[n] = numX*(j+1) + i;

            n++;
        }
    }

    if ( fd != NULL ) fclose(fd);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMSurfaceFileASCIIRaster(surface *ibmSurface)
{
    char      ibmFile[256];
    FILE      *fd;
    PetscInt  readError; char* charError;
    char      string[128];
    ibmMesh   *ibMesh = ibmSurface->ibMsh;

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmSurface->surfaceName.c_str());

    sprintf(ibmFile, "./IBM/%s", ibmSurface->surfaceName.c_str());

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[512];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    //skip the first line
    charError = fgets(string, 128, fd);

    PetscInt    numX, numY;
    PetscReal   xMax, xMin, yMax, yMin, zMax, zMin;
    PetscReal   delX, delY;
    readError = fscanf(fd, "%ld %ld",  &numX, &numY);

    readError = fscanf(fd, "%le %le", &xMin, &xMax);
    readError = fscanf(fd, "%le %le", &yMin, &yMax);
    readError = fscanf(fd, "%le %le", &zMin, &zMax);

    delX = (xMax - xMin)/ (numX-1);
    delY = (yMax - yMin)/ (numY-1);

    ibMesh->nodes = numX * numY;
    ibMesh->elems = 2 * (numX-1) * (numY-1);

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(ibMesh->nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // set the node co-rdinates from the file and initialize the node velocity to 0
    PetscInt n = 0;
    for(PetscInt j = 0; j < numY; j++)
    {
        for(PetscInt i = 0; i < numX; i++)
        {
            ibMesh->nCoor[n].x = xMin + i*delX;

            ibMesh->nCoor[n].y = yMin + j* delY;

            readError = fscanf(fd, "%le", &(ibMesh->nCoor[n].z));

            // translate the body based on the base location
            mSum(ibMesh->nCoor[n], ibmSurface->baseLocation);

            // initialize node velocity to 0
            mSetValue(ibMesh->nU[n], 0.0);
            n++;
        }
    }

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(ibMesh->elems * sizeof(PetscInt), &(ibMesh->nID3));

    n = 0;
    for(PetscInt j = 0; j < numY-1; j++)
    {
        for(PetscInt i = 0; i < numX-1; i++)
        {
            ibMesh->nID1[n] = numX*j + i;
            ibMesh->nID2[n] = numX*(j) + i+1;
            ibMesh->nID3[n] = numX*(j+1) + i;

            n++;

            ibMesh->nID1[n] = numX*(j) + i+1;
            ibMesh->nID2[n] = numX*(j+1) + i+1;
            ibMesh->nID3[n] = numX*(j+1) + i;

            n++;
        }
    }

    if ( fd != NULL ) fclose(fd);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode writeSTLFile(ibmObject *ibmBody)
{
    ibmMesh       *ibMesh  = ibmBody->ibMsh;

    PetscPrintf(PETSC_COMM_WORLD, "writing IBM mesh to STL...\n");

    PetscMPIInt rank;   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    PetscInt    n1, n2, n3;

    if(!rank)
    {
        FILE *f;
        char fileName[256];
        sprintf(fileName, "./IBM/%s.stl", ibmBody->bodyName.c_str());
        f = fopen(fileName, "w");

        if(!f)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName);
            fatalErrorInFunction("writeSTLFile",  error);
        }
        else
        {
            fprintf(f, "solid %s\n", ibmBody->bodyName.c_str());

            for(PetscInt i = 0; i < ibMesh->elems; i++)
            {
                n1 = ibMesh->nID1[i]; n2 = ibMesh->nID2[i]; n3 = ibMesh->nID3[i];

                fprintf(f, "   facet normal %le %le %le\n", ibMesh->eN[i].x, ibMesh->eN[i].y, ibMesh->eN[i].z);
                fprintf(f, "      outer loop\n");
                fprintf(f, "         vertex  %le %le %le\n", ibMesh->nCoor[n1].x, ibMesh->nCoor[n1].y, ibMesh->nCoor[n1].z);
                fprintf(f, "         vertex  %le %le %le\n", ibMesh->nCoor[n2].x, ibMesh->nCoor[n2].y, ibMesh->nCoor[n2].z);
                fprintf(f, "         vertex  %le %le %le\n", ibMesh->nCoor[n3].x, ibMesh->nCoor[n3].y, ibMesh->nCoor[n3].z);
                fprintf(f, "      endloop\n");
                fprintf(f, "   endfacet\n");


            }
            fprintf(f, "endsolid\n");

            fclose(f);
        }
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMBodyFileUCD(ibmObject *ibmBody)
{

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmBody->bodyName.c_str());

    char ibmFile[256];
    sprintf(ibmFile, "./IBM/%s", ibmBody->bodyName.c_str());

    FILE *fd;
    PetscInt  readError; char* charError;

    // pointer to the mesh object of the current body
    ibmMesh *ibMesh = ibmBody->ibMsh;

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[530];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    // local variables to read file variables
    PetscInt         nodes = 0;                          // number of ib nodes
    PetscInt         elem = 0;                           // numbe of elements
    PetscInt         tmp;                            // temporary access variables
    char             string[128], tmpS[128];                    // temporary string variable

    // skip the first 3 header lines of the file created by pointwise.
    charError = fgets(string, 128, fd);
    charError = fgets(string, 128, fd);
    charError = fgets(string, 128, fd);

    // scan the number of nodes and elements
    readError = fscanf(fd, "%ld %ld %ld %ld %ld",  &nodes, &elem, &tmp, &tmp, &tmp);

    ibMesh->nodes = nodes;
    ibMesh->elems = elem;

    if(nodes == 0)
    {
      char error[512];
      sprintf(error, "nodes read = 0. check the IBM file and make sure there are no header comments %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // read the node co-rdinates from the file and initialize the node velocity to 0
    for(PetscInt i = 0; i < nodes; i++)
    {
      readError = fscanf(fd, "%ld %le %le %le", &tmp, &ibMesh->nCoor[i].x, &ibMesh->nCoor[i].y, &ibMesh->nCoor[i].z);

      // translate the body based on the based location
      mSum(ibMesh->nCoor[i], ibmBody->baseLocation);

      // initialize node velocity to 0
      mSetValue(ibMesh->nU[i], 0.0);
    }

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID3));

    // read the nodes that form the 3 vertices of an element
    for(PetscInt i = 0; i < elem; i++)
    {
      readError = fscanf(fd, "%ld %ld %s %ld %ld %ld", &tmp, &tmp, tmpS, &ibMesh->nID1[i], &ibMesh->nID2[i], &ibMesh->nID3[i]);

      // node numbers start from 0
      ibMesh->nID1[i] = ibMesh->nID1[i] - 1;
      ibMesh->nID2[i] = ibMesh->nID2[i] - 1;
      ibMesh->nID3[i] = ibMesh->nID3[i] - 1;

    }

    if ( fd != NULL ) fclose(fd);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readIBMSurfaceFileUCD(surface *ibmSurface)
{

    PetscPrintf(PETSC_COMM_WORLD, "Reading IBM body: %s\n", ibmSurface->surfaceName.c_str());

    char ibmFile[256];
    sprintf(ibmFile, "./IBM/%s", ibmSurface->surfaceName.c_str());

    FILE *fd;
    PetscInt  readError; char* charError;

    // pointer to the mesh object of the current body
    ibmMesh *ibMesh = ibmSurface->ibMsh;

    // open the ibm object file in read mode
    fd = fopen(ibmFile, "r");

    if(fd == NULL)
    {
      char error[530];
      sprintf(error, "cannot open file %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    // local variables to read file variables
    PetscInt         nodes = 0;                          // number of ib nodes
    PetscInt         elem = 0;                           // numbe of elements
    PetscInt         tmp;                            // temporary access variables
    char             string[128], tmpS[128];                    // temporary string variable

    // skip the first 3 header lines of the file created by pointwise.
    charError = fgets(string, 128, fd);
    charError = fgets(string, 128, fd);
    charError = fgets(string, 128, fd);

    // scan the number of nodes and elements
    readError = fscanf(fd, "%ld %ld %ld %ld %ld",  &nodes, &elem, &tmp, &tmp, &tmp);

    ibMesh->nodes = nodes;
    ibMesh->elems = elem;

    if(nodes == 0)
    {
      char error[512];
      sprintf(error, "nodes read = 0. check the IBM file and make sure there are no header comments %s\n", ibmFile);
      fatalErrorInFunction("readIBMObjectMesh",  error);
    }

    // allocate memory for the x, y and z co-ordinates of the nodes
    PetscMalloc(nodes * sizeof(Cmpnts), &(ibMesh->nCoor));

    // allocate memory for the node velocity
    PetscMalloc(nodes * sizeof(Cmpnts), &(ibMesh->nU));

    // read the node co-rdinates from the file and initialize the node velocity to 0
    for(PetscInt i = 0; i < nodes; i++)
    {
      readError = fscanf(fd, "%ld %le %le %le", &tmp, &ibMesh->nCoor[i].x, &ibMesh->nCoor[i].y, &ibMesh->nCoor[i].z);

      // translate the body based on the based location
      mSum(ibMesh->nCoor[i], ibmSurface->baseLocation);

      // initialize node velocity to 0
      mSetValue(ibMesh->nU[i], 0.0);
    }

    // allocate memory for the pointer to nodes (n1, n2 and n3) of an element
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID1));
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID2));
    PetscMalloc(elem * sizeof(PetscInt), &(ibMesh->nID3));

    // read the nodes that form the 3 vertices of an element
    for(PetscInt i = 0; i < elem; i++)
    {
      readError = fscanf(fd, "%ld %ld %s %ld %ld %ld", &tmp, &tmp, tmpS, &ibMesh->nID1[i], &ibMesh->nID2[i], &ibMesh->nID3[i]);

      // node numbers start from 0
      ibMesh->nID1[i] = ibMesh->nID1[i] - 1;
      ibMesh->nID2[i] = ibMesh->nID2[i] - 1;
      ibMesh->nID3[i] = ibMesh->nID3[i] - 1;

    }

    if ( fd != NULL ) fclose(fd);

    return 0;
}

//***************************************************************************************************************//
