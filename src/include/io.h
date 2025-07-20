//! \file  io.h
//! \brief Contains i/o operations function headers

#ifndef _IO_H_
#define _IO_H_

//! \brief Struct defining io settings (simulation checkpointing)
struct io_
{
    word      intervalType;                   //!< timeStep or adjustableTime
    PetscReal timeInterval;                   //!< in iterations if intervalType = timeStep, in seconds if intervalType = adjustableTime

    PetscInt  runTimeWrite;                   //!< flag telling if must write at this iteration

    PetscInt  purgeWrite;                     //!< deletes all other files after writing the current

    PetscInt  averaging;                      //!< compute the time-averaged solution fields
    PetscReal avgPrd;                         //!< sampling period in seconds
    PetscReal avgStartTime;                   //!< start time of averaging procedure
    PetscInt  phaseAveraging;                 //!< compute the phase-averaged solution fields
    PetscReal phAvgPrd;                       //!< sampling period in seconds
    PetscReal phAvgStartTime;                 //!< start time of phase averaging procedure
    PetscInt  keBudgets;                      //!< compute kinetic energy budget terms
    PetscInt  writePForce;                    //!< compute pressure on individual elements

    PetscInt  avgWeight;                      //!< number of average snapshots (cumulated at runtime)
    PetscInt  pAvgWeight;                     //!< number of phase average snapshots (cumulated at runtime)
    PetscInt  keAvgWeight;                    //!< number of keBudget average snapshots (cumulated at runtime)

    PetscInt  qCrit;
    PetscInt  l2Crit;
    PetscInt  vgtQg;
    PetscInt  vgtRg;
    PetscInt  vgtQs;
    PetscInt  vgtRs;
    PetscInt  vgtQr;
    PetscInt  windFarmForce;
    PetscInt  sources;
    PetscInt  buoyancy;
    PetscInt  continuity;

    // runtime modifiable
    word      lastModification;

    // catalyst interface
    PetscReal     timeIntervalCatalyst;
    PetscReal     startTimeCatalyst;
    word          outputTypeCatalyst;
    word          ioTypeCatalyst;           //!< type of catalyst output (generic/script)
    word          scriptNameCatalyst;        //!< name of input script

    // access database
    access_   *access;
};


#endif

//! \brief Calls fatal error and exits
void fatalErrorInFunction(const char* functionName, const char* errorMsg);

//! \brief Calls warning
void warningInFunction(const char* functionName, const char* wrngMsg);

//! \brief Retrieves this case name
word thisCaseName();

//! \brief Returns the time list (from folder names) and size contained in a folder
PetscErrorCode getTimeList(const char* dataLoc, std::vector<PetscReal> &timeSeries, PetscInt &ntimes);

//! \brief Returns the file list (from folder names) and number of files contained in a folder
PetscErrorCode getFileList(const char* dataLoc, std::vector<word> &fileSeries, PetscInt &nfiles);

//! \brief Finds a word inside a string
PetscInt foundInString(const char *str, word keyword);

//! \brief Checks if file exists
PetscInt file_exist(const char *str);

//! \brief Checks if directory exists
PetscInt dir_exist(const char *str);

//! \brief Count number of files
PetscInt count_files(const char* path);

//! \brief Creat directory (remove stuff it exists)
void createDir(MPI_Comm comm, const char* path);

//! \brief Creat directory (leave stuff it exists)
void createDirNoRemove(MPI_Comm comm, const char* path);

//! \brief Removes directory
void remove_dir(MPI_Comm comm, const char *path2dir);

//! \brief Removes all subdirectory of a given directory
void remove_subdirs(MPI_Comm comm, const char *path2dir);

//! \brief Removes all subdirectory of a given directory except the name provided
void remove_subdirs_except(MPI_Comm comm, const char *path2dir, const word name);

//! \brief Removes all subdirectory except the latest n times inside a given directory and the name provided
void remove_subdirs_except_keep_n(MPI_Comm comm, const char *path2dir, const word name, PetscInt nKeep);

//! \brief Removes all subdirectory of a given directory except the 2 names provided
void remove_subdirs_except2(MPI_Comm comm, const char *path2dir, const word name1, const word name2);

//! \brief Removes all subdirectory of a given directory except the 3 names provided
void remove_subdirs_except3(MPI_Comm comm, const char *path2dir, const word name1, const word name2, const word name3);

//! \brief Removes all subdirectory of a given directory except the 4 names provided
void remove_subdirs_except4(MPI_Comm comm, const char *path2dir, const word name1, const word name2, const word name3, const word name4);

//! \brief Removes all subdirectory except the latest n times inside a given directory and the 4 names provided
void remove_subdirs_except4_keep_n(MPI_Comm comm, const char *path2dir, const word name1, const word name2, const word name3, const word name4, PetscInt nKeep);

//! \brief Set write directory for a given mesh
void SetWriteDir(mesh_ *mesh, const word str);

//! \brief Write binary vector file
void writeBinaryField(mesh_ *mesh, Vec &V, char *file);

//! \brief Initializes io data structure
PetscErrorCode InitializeIO(io_ *io);

//! \brief Re-read IO write parameters
PetscErrorCode RereadIO(domain_ *domain);

//! \brief Called by RereadIO, triggers updates in each domain
PetscErrorCode UpdateInput(io_ *io, word &modified);

//! \brief read and save the different fields stored in the fields/meshName/time folder
PetscErrorCode readFields(domain_ *domain, PetscReal timeValue);

//! \brief Get time name in string format with required digits
word getTimeName(clock_ *clock);

//! \brief Get start time name in string format with required digits
word getStartTimeName(clock_ *clock);

//! \brief Get time value in string format with required digits
word getArbitraryTimeName(clock_ *clock, double timeValue);

//! \brief Sets the runTimeWrite flag and creates initial output directory
PetscErrorCode setRunTimeWrite(domain_ *domain);

//! \brief Creates the checkpoint directory with a trailing _tmp (indicates that it is not finalized)
PetscErrorCode initializeTimeDir(io_ *io);

//! \brief Renames the temporary time directory to the final name (without _tmp) and deletes prev chekpoints if purge is active
PetscErrorCode finalizeTimeDir(io_ *io);

//! \brief Write output fields
PetscErrorCode writeFields(io_ *io);

//! \brief Read mandatory PetscReal from dictionary
PetscErrorCode readDictDouble(const char *dictName, const char *keyword, PetscReal *value);

//! \brief Read mandatory vector from dictionary
PetscErrorCode readDictVector(const char *dictName, const char *keyword, Cmpnts *value);

//! \brief Read mandatory 2D vector from dictionary
PetscErrorCode readDictVector2D(const char *dictName, const char *keyword, Cpt2D *value);

//! \brief Read mandatory PetscInt from dictionary
PetscErrorCode readDictInt(const char *dictName, const char *keyword, PetscInt *value);

//! \brief Read mandatory word from dictionary
PetscErrorCode readDictWord(const char *dictName, const char *keyword, word *value);

//! \brief Read mandatory word and PetscReal from dictionary
PetscErrorCode readDictWordAndDouble(const char *dictName, const char *keyword, word *value1, PetscReal *value2);

//! \brief Read mandatory word and vector from dictionary
PetscErrorCode readDictWordAndVector(const char *dictName, const char *keyword, word *value1, Cmpnts *value2);

//! \brief Read mandatory PetscReal from sub-dictionary
PetscErrorCode readSubDictDouble(const char *dictName, const char *subdict, const char *keyword, PetscReal *value);

//! \brief Read mandatory PetscInt from sub-dictionary
PetscErrorCode readSubDictInt(const char *dictName, const char *subdict, const char *keyword, PetscInt *value);

//! \brief Read mandatory vector from sub-dictionary
PetscErrorCode readSubDictVector(const char *dictName, const char *subdict, const char *keyword, Cmpnts *value);

//! \brief Read mandatory word from sub-dictionary
PetscErrorCode readSubDictWord(const char *dictName, const char *subdict, const char *keyword, word *value);

//! \brief Read mandatory array of PetscInt from sub-dictionary
PetscErrorCode readSubDictIntArray(const char *dictName, const char *subdict, const char *keyword, labelList &value);

std::string* readSubDictWordArray(const char *dictName, const char *subdict, const char *keyword, PetscInt numW);

//! \brief Trims a string removing pre and trail spaces
word trim(const word& str);

//! \brief Check if is a number
bool isNumber(const word& str);
