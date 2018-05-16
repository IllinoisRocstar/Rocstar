/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: Rocin.C,v 1.76 2009/04/08 18:01:57 mtcampbe Exp $

/** \file Rocin.C
 *  Implementation of Rocing for creation of Roccom window(s) from HDF files.
 */
/*  Author John Norris
 *  Initial date:   March 19, 2004
 */
/*  Modified by Phillip Alexander
 *  April 8, 2004
 *  -updated to Roccom3 compliance
 *  -added support for use as a Roccom module
 *
 * Modified by Eric Shaffer
 * Feb. 20, 2006
 * - fixes to read Boeing CGNS files
 */

#include <algorithm>
#include <set>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <strings.h>
#include <errno.h>
#include <iomanip>
#include <vector>
#include <list>

#ifndef _NO_GLOB_
#include <glob.h>
#endif
//#include "Matcher.h"
//#include "Pattern.h"
#include <fnmatch.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include "Directory.H"
//#endif

#include "Rocin.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
USE_COM_NAME_SPACE
#endif

// #define DEBUG_DUMP_PREFIX "."
// #define DEBUG_MSG(x) std::cout << "ROCIN DEBUG: " << __LINE__ << ": " << x << std::endl
#define DEBUG_MSG(x)

#ifdef DEBUG_DUMP_PREFIX
static void DebugDump(std::ofstream& fout, int n, int dType, void* data)
{
  int i;
  switch (dType) {
    case COM_CHAR:
    case COM_CHARACTER: {
      char* p = (char*)data;
      for (i=0; i<n; ++i)
        fout << i << ": " << (int)p[i] << '\n';
      break; }
    case COM_INT:
    case COM_INTEGER: {
      int* p = (int*)data;
      for (i=0; i<n; ++i)
        fout << i << ": " << p[i] << '\n';
      break; }
    case COM_FLOAT:
    case COM_REAL: {
      float* p = (float*)data;
      for (i=0; i<n; ++i)
        fout << i << ": " << p[i] << '\n';
      break; }
    case COM_DOUBLE:
    case COM_DOUBLE_PRECISION: {
      double* p = (double*)data;
      for (i=0; i<n; ++i)
        fout << i << ": " << p[i] << '\n';
      break; }
    default:
      fout << "Unhandled COM datatype " << dType << std::endl;
      break;
  }
}
#endif // DEBUG_DUMP_PREFIX

/**
 ** Error callback for glob.
 **/
#ifndef _NO_GLOB_
static int glob_error(const char* epath, int gerrno)
{
  if ( !COMMPI_Initialized())
    std::cerr << "Rocstar: Error (glob): error in pattern `" << epath << "': " 
              << strerror(gerrno) << std::endl;
  return 0;
}


/// Cast glob_error for portability (IBMSP has a different prototype)
template < class T>
inline T cast_err_func( int (*glob)(const char *, int, T, glob_t*),
                        int (*errfunc)(const char*, int)) 
{ return reinterpret_cast<T>(errfunc); }
#endif
/**
 ** Perform HDF4 error-checking.
 **
 ** Check the return value of the HDF4 function and print an error message
 ** if it is FAIL
 **
 ** \param routine The HDF4 function to call. (Input)
 ** \param args The argument list, in parentheses. (Input)
 ** \param val The return value. (Input)
 **/
#define HDF4_CHECK_RET(routine, args, retval) \
  {					      \
    if((routine args) == FAIL) {	      \
      sleep(1);				      \
      if((routine args) == FAIL) {	      \
	sleep(1);			      \
	if((routine args) == FAIL) {	      \
	  sleep(1);			      \
	  if((routine args) == FAIL) {			    \
	    std::cerr << "Rocstar: Error: " #routine " (line "		\
		      << __LINE__ << " in " << __FILE__ << ") failed: " \
		      << HDF4::error_msg() << std::endl;		\
	    retval;							\
	  } } } }							\
  }
/* 
#define HDF4_CHECK_RET(routine, args, retval)	\
{ \
  if ((routine args) == FAIL) {  \
    std::cerr << "Rocstar: Error: " #routine " (line " \
              << __LINE__ << " in " << __FILE__ << ") failed: " \
              << HDF4::error_msg() << std::endl; \
    retval; \
  } \
}
*/
#define HDF4_CHECK(routine, args) HDF4_CHECK_RET(routine, args, )

/**
 ** Automatically close open files, datasets, etc.
 **
 ** \param id The file id.
 ** \param f The function used to close the file.
 **/
template <typename T1, typename T2 = T1>
class AutoCloser
{
  public:
    inline AutoCloser(T1 id, T2 (*f)(T1)) : m_id(id), m_close(f) {}
    inline ~AutoCloser() { m_close(m_id); } 
  private:
    T1 m_id;
    T2 (*m_close)(T1);
};

#ifdef USE_CGNS
/**
 ** Perform CGNS error-checking.
 **
 ** Check the return value of the CGNS function and print an error message
 ** if it is non-zero.
 **
 ** \param routine The CGNS function to call. (Input)
 ** \param args The argument list, in parentheses. (Input)
 ** \param val The return value. (Input)
 **/
# define CG_CHECK_RET(routine, args, retval) \
{ \
  int ier = routine args; \
  if (ier != 0) { \
    std::cerr << "Rocstar: Error: " #routine " (line " \
              << __LINE__ << " in " << __FILE__ << ") failed: " \
              << cg_get_error() << std::endl; \
    retval; \
  } \
}
# define CG_CHECK(routine, args) CG_CHECK_RET(routine, args, )

/**
 ** Automatically save and restore the current working directory.
 **
 ** This class remembers the cwd and cds to it when it goes out of scope.
 **/
class AutoCDer
{
  public:
    inline AutoCDer()
    : m_cwd(getcwd(NULL, 0)) {}
    inline ~AutoCDer() { chdir(m_cwd); free(m_cwd); }

  private:
    char* m_cwd;
};

/** Boeing Fix
 ** Convert CGNS ElementType_t (defined in cgnslib.h) to a Roccom string 
 **/
bool elementType_to_string(ElementType_t etype, std::string & roccom_elem)
{

switch (etype)                  /* type of CGNS element */
    {

    case BAR_2: roccom_elem = ":b2:";
      break;
    case BAR_3: roccom_elem = ":b3:" ;
      break;
    case TRI_3: roccom_elem = ":t3:";
        break;
    case TRI_6: roccom_elem = ":t6:" ;
        break;
    case QUAD_4: roccom_elem = ":q4:" ;
        break;
    case QUAD_8: roccom_elem = ":q8:" ;
        break;
    case QUAD_9: roccom_elem = ":q9:" ;
      break;
    case TETRA_4: roccom_elem = ":T4:" ;
        break;
    case TETRA_10: roccom_elem = ":T10:";
      break;
    case PYRA_5: roccom_elem = ":P5:";
        break;
    case PYRA_14: roccom_elem = ":P14:" ;
        break;
    case PENTA_6: roccom_elem = ":P6:" ;
        break;
    case PENTA_15: roccom_elem = ":P15:" ;
        break;
    case PENTA_18: roccom_elem = ":P18:" ;
      break;
    case HEXA_8: roccom_elem = ":B8:" ;
        break;
    case HEXA_20: roccom_elem = ":B20:";
      break;
      //case 19: roccom_elem = ":B27:";
      //break;
    default: std::cerr << "Rocstar: Error: CGNS element type of " << etype << "has no Roccom equivalent\n";
      return false;
    }
 return true;
}
#endif //ifdef CGNS


/**
 ** Find the first geometry dataset for the given block.
 **
 ** The dataset may be right after the block header or in an external file.
 ** The caller is responsible for closing the sd_id.
 **/
static int32 FindFirstGeometryDataset(int32 sd_id, const std::string& geomFile,
                                      const std::string& matName,
                                      const std::string& label, int32 hIndex,
                                      int32* index)
{
  int32 geom_id;
  if (geomFile.empty()) {
    // We're not using a separate geometry file, so just return the dataset
    // immediately after the block header.
    geom_id = sd_id;
    *index = hIndex + 1;
    return geom_id;
  }

  // Open the geometry file.
  HDF4_CHECK_RET(geom_id = HDF4::SDstart, (geomFile.c_str(), DFACC_READ),
                 return FAIL);

  // Get the number of datasets.
  int32 dsCount = MAX_NC_VARS, nAttrs, status;
  HDF4_CHECK(HDF4::SDfileinfo, (geom_id, &dsCount, &nAttrs));

  int32 sds_id, rank, size[5], dType;
  char name[MAX_NC_NAME];
  // We're looking for the dataset for this material and section.
  std::string target = matName + '|' + label;

  // Scan through the datasets of the geometry file.
  for (*index=0; *index<dsCount; ++(*index)) {
    HDF4_CHECK_RET(sds_id = HDF4::Select,
                   (geom_id, *index, name, &rank, size, &dType, &nAttrs, dsCount),
                   HDF4::SDend(geom_id); return FAIL);

    HDF4_CHECK_RET(status = HDF4::SDgetdatastrs,
                   (sds_id, NULL, NULL, NULL, name, MAX_NC_NAME),
                   HDF4::SDendaccess(sds_id); continue);
    HDF4::SDendaccess(sds_id);

    // Check if this is the right dataset.
    if (strcasecmp(target.c_str(), name) == 0) { // We've found it!
      return geom_id;
    }
  }

  // No luck.
  HDF4::SDend(geom_id);
  return FAIL;
}

/**
 ** Extract metadata from the list of files.
 **
 ** Open each file, scan the datasets, identify the windows, panes, and
 ** attributes.
 **/
static void scan_files_HDF4( int pathc, char* pathv[],
                             BlockMM_HDF4& blocks, std::string& time,
                             std::map<int32, COM_Type>& HDF2COM)
{
  int32 sd_id, geom_id, sds_id, nDataSets, nAttrs, index, index3;
  int32 rank, dType, start[1];
  int32 size[5] = { 0, 0, 0, 0, 0 };
  int32 geomIdx[3];
  char name[MAX_NC_NAME], label[MAX_NC_NAME];
  char timeLevel[MAX_NC_NAME], units[MAX_NC_NAME];
  char matName[MAX_NC_NAME], varName[MAX_NC_NAME], varType[MAX_NC_NAME];
  int xyz;
  int numPoints=0, ghostPoints, gc, gridType=0;
  std::vector<int> ghostCells;
  char grid[32], ghost[32], location;
  std::string geomFile;
  char buffer[1024];
  Block_HDF4* block;
  int i;

  blocks.clear();

  // Iterate through each file.
  for (i=0; i<pathc; ++i) {
    // Make sure this is an HDF file.
    if (strcmp(&pathv[i][strlen(pathv[i])-4], ".hdf") != 0
        && strcmp(&pathv[i][strlen(pathv[i])-5], ".hdf4") != 0)
      continue;

    HDF4_CHECK_RET(sd_id = HDF4::SDstart, (pathv[i], DFACC_READ), continue);
    AutoCloser<int32, intn> auto0(sd_id, HDF4::SDend);

    // Determine the number of datasets.
    HDF4_CHECK_RET(HDF4::SDfileinfo, (sd_id, &nDataSets, &nAttrs),
                   continue);

    // Iterate through the datasets.
    for (index=0; index<nDataSets; ++index) {
      // Select each non-"fakeDim" dataset in the file.
      HDF4_CHECK_RET(sds_id = HDF4::Select,
                     (sd_id, index, name, &rank, size, &dType, &nAttrs, nDataSets),
                     continue);

      {
        AutoCloser<int32, intn> auto1(sds_id, HDF4::SDendaccess);

        // Get the paneId ('label'), time level ('timeLevel'), 
        // dataset type ('varName'), and window name ('matName'). 
        HDF4_CHECK_RET(HDF4::SDgetdatastrs,
                       (sds_id, label, timeLevel, varName, matName, MAX_NC_NAME),
                       continue);

        // If this isn't a block header, then it isn't interesting.
        if (strncasecmp(varName, "block header", 12) != 0) {
          continue;
        } else {
          // Otherwise, set the time level if not yet set, and skip
          // the blocks that have different time levels.
          if (time.empty()) 
            time = timeLevel;
          else if (time != timeLevel)
            continue;
        }

        // Block header datasets should have only one dimension.
        if (rank != 1) {
          std::cerr << "Rocstar: Warning: Dataset " << index << " in file " << pathv[i]
                    << " is marked as a block header, but its rank is not 1 "
                    << "(rank ==  " << rank << ')' << std::endl;
          continue;
        }

        ghostCells.clear();
        geomFile.erase();
        if (dType == DFNT_CHAR8) {
          // This data file is written for version 1.1.  Good!
          start[0] = 0;
          std::vector<char8> data(size[0]);
          HDF4_CHECK_RET(HDF4::SDreaddata,
                         (sds_id, start, NULL, size, (VOIDP)&data[0]),
                         continue);

          // Parse the data into grid type, ghost info, and external
          // geometry file.
          std::istringstream sin(&data[0]);
          sin.getline(grid, sizeof(grid), '|');
          sin.getline(ghost, sizeof(ghost), '|');
          buffer[0] = '\0';
          sin >> buffer;

          // The first token is the grid type.
          if (strcasecmp(grid, "structured") == 0 || strcmp(grid, "2") == 0) {
            gridType = 2;
          } else if (strcasecmp(grid, "unstructured") == 0
                     || strcmp(grid, "3") == 0) {
            gridType = 3;
          } else if (strcasecmp(grid, "discontinuous") == 0
                     || strcmp(grid, "4") == 0) {
            gridType = 4;
          } else if (strcasecmp(grid, "polydata") == 0
                     || strcmp(grid, "5") == 0) {
            gridType = 5;
          } else {
            std::cerr << "Rocstar: Warning: The block header (dataset " << index
                      << " in file " << pathv[i] << ") does not have a "
                      << "valid grid type (grid = " << grid << ')'
                      << std::endl;
          }

          // The second token is the number of ghost points and cells.  If
          // there is a single number, then this is (probably) the number of
          // layers of ghost cells on a structured grid.  Otherwise, it's the
          // number of ghost points followed by the number of ghost cells of
          // each connectivity table of an unstructured grid.
          // WARNING: When using external grid files, the ghost data will
          // probably not be accurate!

          ghostPoints = 0;
          /**
           ** istringstream is apparently fragile, so we'll use atoi instead.
           **
          std::istringstream sin2(ghost);

          sin2 >> gc;
          // Check if there are more data.
          if (!sin2.eof()) {
            if (sin2.get() == ',') {
              ghostPoints = gc;
              sin2 >> gc;
            }
          }
          if (gc < 0) std::cerr << "Rocstar: Error: Bad num ghosts " << gc << " from string: " << ghost << std::endl;
          ghostCells.push_back(gc);
          // Check if there are more data.
          while (!sin2.eof()) {
            if (sin2.get() == ',') {
              sin2 >> gc;
              if (gc < 0) std::cerr << "Rocstar: Error: Bad ghosts " << gc << std::endl;
              ghostCells.push_back(gc);
            }
          }
           **
           **
           **/
          gc = atoi(ghost);
          char* c = strchr(ghost, ',');
          if (c != NULL) {
            ++c;
            ghostPoints = gc;
            gc = atoi(c);
            c = strchr(c, ',');
          }
          if (gc < 0) std::cerr << "Rocstar: Warning: Bad num ghosts " << gc << " from string: " << ghost << std::endl;
          ghostCells.push_back(gc);
          // Check if there are more data.
          while (c != NULL) {
            ++c;
            gc = atoi(c);
            c = strchr(c, ',');
            if (gc < 0) std::cerr << "Rocstar: Warning: Bad ghosts " << gc << std::endl;
            ghostCells.push_back(gc);
          }

          // The third token (if present) is an external geometry file.
          if (buffer[0] != '\0') {
            if (buffer[0] == '/') {
              geomFile = buffer;
            } else {
              // Normalize the filename.
              geomFile = pathv[i];
              std::string::size_type pos = geomFile.find_last_of('/');
              if (pos != std::string::npos) {
                geomFile.erase(pos + 1);
                geomFile += buffer;

                // If normalized filename does not exist, then 
                // try to find the mesh file with the given name.
                struct stat sb;
                if (stat(geomFile.c_str(), &sb) != 0)
                  geomFile = buffer;
              }
              else
                geomFile = buffer;
            }

            // Make sure the file exists.
            struct stat sb;
            if (stat(geomFile.c_str(), &sb) != 0) {
              std::cerr << "Rocstar: Warning: The block header (dataset " << index
                        << " in file " << pathv[i] << ") has an invalid "
                        << "filename for the external geometry data.  The "
                        << "file does not exist: " << geomFile << std::endl;
              continue;
            }
          }
        } else {
          std::cerr << "Rocstar: Warning: Dataset " << index << " in file " << pathv[i]
                    << " is marked as a block header, but its datatype is "
                    << "not char8." << std::endl;
          continue;
        }
      }

      geom_id = FindFirstGeometryDataset(sd_id, geomFile, matName, label,
                                         index, &index3);
      if (geom_id == FAIL) {
        std::cerr << "Rocstar: Warning: Could not open geometry dataset for section "
                  << label << " of material " << matName << std::endl;
        continue;
      }

      // Collect data on the nodal coordinate datasets.
      for (xyz=0; xyz<6; xyz+=2,++index3) {
        HDF4_CHECK_RET(sds_id = HDF4::Select,
                       (geom_id, index3, name, &rank, size, &dType, &nAttrs),
                       continue);
        AutoCloser<int32, intn> auto2(sds_id, HDF4::SDendaccess);

        geomIdx[xyz/2] = index3;

        if (xyz == 0) {
          units[0] = '\0';
          if ((gridType == 2 && rank != 3) || (gridType != 2 && rank != 1)) {
            numPoints = size[0] = size[1] = size[2] = 0;
            geomIdx[0] = geomIdx[1] = geomIdx[2] = FAIL;
            break;
          }

          HDF4_CHECK(HDF4::SDgetdatastrs,
                     (sds_id, name, units, NULL, varType, MAX_NC_NAME));

          if (gridType == 2) {
            numPoints = size[0] * size[1] * size[2];
          } else {
            char* pos = strchr(name, '|');

            if (pos && *(pos+2)=='0')
              numPoints = ghostPoints = 0; // Empty pane.
            else {
              numPoints = size[0];
              // Get a trustworthy count of the number of ghost nodes.
              std::istringstream in(pos+3);
              in >> ghostPoints;
            }
          }
        }
      }

      // Determine the pane id.  Skip any non-numeric prefix such as "block".
      std::string lbl(label);
      std::string::size_type pos = lbl.find_first_of("0123456789");
      if (pos != std::string::npos)
        lbl.erase(0, pos);

      std::istringstream sin(lbl); 
      int paneId;
      sin >> paneId;

      // Create and initialize a Block object.
      block = new Block_HDF4(pathv[i], geomFile, geomIdx, paneId, timeLevel,
                        units, numPoints, ghostPoints);

      // Add grid info.
      if (gridType == 2) {
        // Structured grids only have dimensions.
        block->m_gridInfo.push_back(GridInfo_HDF4(size, ghostCells[0]));
      } else {
        // Get info for each connectivity table.  If we were given ghost
        // cell info, then we know how many connectivity tables we have.
        // Otherwise we have to keep going until we find a dataset that is
        // obviously not a connectivity table.
        std::vector<int>::iterator p;
        for (p=ghostCells.begin();
             ghostCells.empty() || p!=ghostCells.end();
             ++p,++index3) {
          sds_id = HDF4::Select(geom_id, index3, name, &rank, size,
                                &dType, &nAttrs);
          if (sds_id == FAIL) {
            if (!ghostCells.empty())
              std::cerr << "Rocstar: Warning: Select() failed on connectivity dataset "
                        << index3 << " in file "
                        << (geomFile.empty() ? pathv[i] : geomFile.c_str())
                        << std::endl;
            break;
          }
          AutoCloser<int32, intn> auto3(sds_id, HDF4::SDendaccess);

          // If the rank isn't 2, then it's not a connectivity table.
          if (rank != 2) {
            break;
          }

          // Get metadata on the connectivity table.
          HDF4_CHECK_RET(HDF4::SDgetdatastrs,
                         (sds_id, label, units, varName, varType, MAX_NC_NAME),
                         continue);

          int ghostCount = 0;
          if (!ghostCells.empty())
            ghostCount = *p;

          // A connectivity table must start with ':'
          if (varName[0] != ':') { /* Old, pre-roccom3 dataset */
            /* HIDEOUS HACK: should figure out real element type */
            if (gridType == 3 || gridType == 4) {
              if (size[0] == 4)
                strcpy(varName,":T4:"); 
              else if (size[0] == 5)
                strcpy(varName,":P5:"); 
              else if (size[0] == 6)
                strcpy(varName,":P6:"); 
              else if (size[0] == 8)
                strcpy(varName,":B8:"); 
              else if (size[0] == 10)
                strcpy(varName,":T10:"); 
              else if (size[0] == 20)
                strcpy(varName,":B20:"); 
              else
                strcpy(varName,":unknown");
            } else if (gridType == 5) {
              if (size[0] == 3)
                strcpy(varName,":t3:"); 
              else if (size[0] == 4)
                strcpy(varName,":q4:"); 
              else if (size[0] == 6)
                strcpy(varName,":t6:"); 
              else if (size[0] == 8)
                strcpy(varName,":q8:"); 
              else
                strcpy(varName,":unknown");
            } else {
              strcpy(varName,":unknown");
            }

            if (1) {
              std::cerr << "Rocstar: Warning: reading old-style connectivity dataset "
                      << varName << " from " << index3 << " in file "
                      << (geomFile.empty() ? pathv[i] : geomFile.c_str())
                      << std::endl;
            }
          } else {
            // If the table is empty, set number of elements to 0.
            char* pos = strchr(label, '|');
            if (pos && *(pos+2)=='0') {
              size[1] = ghostCount = 0;
            } else {
              // Get a trustworthy count of the number of ghost elements.
              std::istringstream in(pos+3);
              in >> ghostCount;
            }
          }

          if (size[1]<0 || ghostCount < 0) {
            std::cerr << "Rocstar: Error: Bad number of element or ghosts for attribute " << varName << " nElem=" << size[1] << " nGhost=" << ghostCount << std::endl;
            COM_assertion_msg(0,
                        "Bad number of element or ghosts for attribute.");
          }
          // Store the mesh info.
          block->m_gridInfo.push_back(GridInfo_HDF4(varName, size[1], ghostCount, index3));
        }
      }
      if (geom_id != sd_id) {
        HDF4::SDend(geom_id);
      } else {
        index = index3;
      }

      // Now search for datasets with variable data.
      for (++index; index<nDataSets; ++index) {
        HDF4_CHECK_RET(sds_id = HDF4::Select,
                       (sd_id, index, name, &rank, size, &dType, &nAttrs, nDataSets),
                       continue);
        AutoCloser<int32, intn> auto4(sds_id, HDF4::SDendaccess);

        // Get metadata on the variable.
        HDF4_CHECK_RET(HDF4::SDgetdatastrs,
                       (sds_id, label, units, varName, varType, MAX_NC_NAME),
                       continue);

        if (varName[0] == '\0' || varName[0] == ':') {
          // Hmm, this is a geometry dataset.
          continue;
        }

        if (strncasecmp(varName, "block header", 12) == 0) {
          // We've found another block header.  We're done.
          // Reverse index by one, as it will be incremented by the outer-loop
          --index;
          break;
        }

        // Determine the location. Format of label is (FIXME: complete?)
        //    "<component letter: 1,2,3,x,y,z>-<attribute>"
        //    "|<location letter: w,p,c,n,e>"
        //    " AT TIME <time SS.MMMUUU> "
        int ng=0; bool is_null=false;
        location = 'w';
        if ( label[0] != ':') {
          char* pos = strchr(label, '|');
          if (pos == NULL) {
            /* Old HDF file without nice labels. For backward compatability,
               guess the data location from the size of the data:
            */
            int sz=size[0]; // Should this be size[rank-1]?

            location='p'; /* guess a panel attribute */

            int num=block->m_numNodes, g=block->m_numGhostNodes;
            if (num==sz) {location='n'; ng=0; /* nodes */ }
            num+=g;
            if (num==sz) {location='n'; ng=g; /* nodes and ghosts */ }

            if (1) {
              num=block->m_gridInfo[0].m_numElements; g=block->m_gridInfo[0].m_numGhostElements;
              if (num==sz) {location='e'; ng=0; /* elements */ }
              num+=g;
              if (num==sz) {location='e'; ng=g; /* elements and ghosts */ }
            }

            if (0) { /* annoying compatability warning */
               std::cerr << "Rocstar: Warning: Compatability warning: guessing variable '"
                      << varName << "' in dataset " << index << " in file "
                      << pathv[i] << " with label " << label 
                      << " is of type " << location 
                      << "(" << sz << ")" << std::endl;
               /* display size of variable */
               std::cerr << "    size of variable: ";
               for (int i=0;i<rank;i++) std::cerr<<size[i]<<',';
               std::cerr << std::endl;
            }
          } 
          else { /* New HDF file: gives location of data directly. */
            location = *(pos + 1);
            // Also obtain the number of ghost items
            if ( *(pos+2) == '0')
              { size[0] = ng = 0; is_null=true; }
            else if ( (*(pos + 2)==',' || *(pos + 2)=='@'))
              { ng = std::atoi( pos+3); is_null=*(pos + 2)=='@'; }
          }
        }

        // Determine which component of the variable we're analyzing.
        std::istringstream str(varType);
        int vType = 0;
        str>> vType;

        int nc = 0;
        switch (vType) {
        case 0:
          // Scalar data.
          nc = 1;
          break;

        case 1:
          // First component of vector data.
          nc = 3;
          break;

        case 11:
          // First component of tensor data.
          nc = 9;
          break;

        case 2: case 3: {
          // Second or third component of vector data.
          VarInfo_HDF4 &vi = block->m_variables.back();
          if (vi.m_name != varName || vi.m_indices.size() != 3) {
            // We shouldn't be seeing 2 or 3 before we see 1.
            std::cerr << "Rocstar: Warning: Found a vector component named " << varName
                      << " for which there is no X component at "
                      << " dataset " << index << " (" << name
                      << ") in file " << pathv[i] << std::endl;
            continue;
          }
          vi.m_indices[vType-1] = index;
          vi.m_is_null[vType-1] = is_null;
          continue;
        }
        case 12: case 13:
        case 21: case 22: case 23:
        case 31: case 32: case 33: {
          // Second, third, etc component of tensor data.
          VarInfo_HDF4 &vi = block->m_variables.back();
          if (vi.m_name != varName || vi.m_indices.size() != 9) {
            // We shouldn't be seeing any of these before we see 11.
            std::cerr << "Rocstar: Warning: Found a tensor component named " << varName
                      << " for which there is no (0,0) component at "
                      << "dataset " << index << " (" << name
                      << ") in file " << pathv[i] << std::endl;
            continue;
          }
          vi.m_indices[3*(vType/10)+(vType%10)-4] = index;
           vi.m_is_null[3*(vType/10)+(vType%10)-4] = is_null;
          continue;
        }
        default:
          // I don't know what kind of variable this is supposed to be.
          // Assuming scalar.
          vType = 0;
          nc = 1;
          break;
        }

        int offset = 0; // Starting point for actual attribute name.
        if ( varName[0]>='0' && varName[0]<='9') {
          if ( std::strncmp( varName, "1-", 2)==0)
            offset = 2; // Skip the "1-" part
          else { // Append to the back
            VarInfo_HDF4 &vi = block->m_variables.back();
            unsigned int id = std::atoi( varName);
            const char *aName = std::strchr(varName, '-');

            if ( aName==NULL || vi.m_name != (aName+1) ||
                 vi.m_indices.size() != id-1) {
              // The individual components must be stored consecutively
              std::cerr << "Rocstar: Warning: Found a vector named " << vi.m_name
                        << " whose components are not stored together "
                        << " in file " << pathv[i] << std::endl;
              continue;
            }

            vi.m_indices.push_back( index);
            vi.m_is_null.push_back( is_null);
            DEBUG_MSG("Correction: variable '" << vi.m_name << "' actually has "
                      << vi.m_indices.size() << " components.");
            continue;
          }
        }

        DEBUG_MSG("Found variable '" << &varName[offset] << "', location == '"
                  << location << "', nComp == " << nc << ", nItems == "
                  << size[0] << ", nGhost == " << ng << ", isNull == "
                  << is_null);
        block->m_variables.push_back
          (VarInfo_HDF4(&varName[offset], location, HDF2COM[dType], 
                   units, nc, index, size[0], ng, is_null));
      }
      blocks.insert(BlockMM_HDF4::value_type(matName, block));
    }
  }
}

static void load_data_HDF4( BlockMM_HDF4::iterator p,
                            const BlockMM_HDF4::iterator& end,
                            const std::string& window, 
                            const MPI_Comm* comm, int rank, int nprocs)
{
  int local, i, ne = 0;
  Block_HDF4* block;
  std::string name;
  bool structured = 0;
  int32 sd_id = 0, sds_id, start[3], size[3];  
  void* data;
  std::vector<int32>::iterator r;
  std::vector<GridInfo_HDF4>::iterator s;
  int with_pane=0;

  for ( ; p!=end; ++p, ++with_pane) {
    block = (*p).second;

    // Panes are a local construct, so make sure that this pane is local.
    local = COM_get_status( window.c_str(), block->m_paneId)>=0;

    if (!local) continue;

    if (block->m_geomFile.empty()) {
      HDF4_CHECK_RET(sd_id = HDF4::SDstart,
                     (block->m_file.c_str(), DFACC_READ),
                     continue);
    } else {
      HDF4_CHECK_RET(sd_id = HDF4::SDstart,
                     (block->m_geomFile.c_str(), DFACC_READ),
                     continue);
    }

    structured = (block->m_gridInfo.size() && 
                  block->m_gridInfo.front().m_name.substr(0,3)== ":st");

    // Let Roccom manage memory for nodal coordinates
    // Moved allocation outside of dimension loop.  Not sure if this fully
    // accounts for single registration of "nc" versus "1-nc", "2-nc", etc.

    name = window + ".nc";
    DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == " << block->m_paneId << ", ptr == NULL, arg == 0 )");
    COM_resize_array(name.c_str(), block->m_paneId, NULL, 1);

    // Read the nodal coordinates.
    for (i=0; i<3; ++i) {
      // Read the mesh only if number of nodes is positive
      name = window + '.' + (char)('1' + i) + "-nc";
      COM_get_array(name.c_str(), block->m_paneId, &data);

      if (block->m_numNodes && data) {
        HDF4_CHECK_RET(sds_id = HDF4::SDselect, (sd_id, block->m_indices[i]),
                       continue);
        AutoCloser<int32, intn> auto1(sds_id, HDF4::SDendaccess);

        if (structured) {
          start[0] = start[1] = start[2] = 0;
          HDF4_CHECK(HDF4::SDreaddata,
                     (sds_id, start, NULL, block->m_gridInfo.front().m_size, data));
        } else {
          start[0] = 0;
          size[0] = block->m_numNodes;
          HDF4_CHECK(HDF4::SDreaddata, (sds_id, start, NULL, size, data));;
        }
#ifdef DEBUG_DUMP_PREFIX
        {
          std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.' + block->time_level + ".hdf").c_str());
          DebugDump(fout, block->m_numNodes, COM_DOUBLE, data);
        }
#endif // DEBUG_DUMP_PREFIX
      }
    }

    // Read the connectivity tables.
    if (!structured) {
      ne = 0;
      for (s=block->m_gridInfo.begin(); s!=block->m_gridInfo.end(); ++s) {
        name = window + "." + (*s).m_name;
        DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == " << block->m_paneId << ", ptr == " << &data << ", arg == 1 )");
        COM_resize_array(name.c_str(), block->m_paneId, &data, 1);

        start[0] = start[1] = 0;
        std::istringstream sin((*s).m_name);
        sin.get(); sin.get(); // Get rid if the leading two letter.
        sin >> size[0];
        size[1] = (*s).m_numElements;
        ne += (*s).m_numElements;

        if ( size[1] && data) {
          HDF4_CHECK_RET(sds_id = HDF4::SDselect, (sd_id, (*s).m_index),
                         continue);
          HDF4_CHECK(HDF4::SDreaddata, (sds_id, start, NULL, size, data));
#ifdef DEBUG_DUMP_PREFIX
        {
          std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.' + block->time_level + ".hdf").c_str());
          DebugDump(fout, size[0] * size[1], COM_INT, data);
        }
#endif // DEBUG_DUMP_PREFIX
          HDF4::SDendaccess(sds_id);
        }
      }
    }

    if (!block->m_geomFile.empty()) {
      HDF4::SDend(sd_id);
      HDF4_CHECK_RET(sd_id = HDF4::SDstart, (block->m_file.c_str(), DFACC_READ),
                     continue);
    }

    // Read the attribute data.
    std::vector<VarInfo_HDF4>::iterator q;
    for (q=block->m_variables.begin(); q!=block->m_variables.end(); ++q) {
      // Skip NULL attributes.
      if ( (*q).m_is_null[0]) continue;
      // Read in window attribute only for the first pane.
      if ( with_pane && (*q).m_position=='w') continue;

      int loop = 1;

      //edited
      name = window + '.' + (*q).m_name;
      if ( (*q).m_position!='w')  { // Allocate for nonwindow attribute
        DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == " << block->m_paneId << ", ptr == " << &data << ", arg == 1 )");
        COM_resize_array(name.c_str(), block->m_paneId, &data,1);
      } else {
        COM_get_array(name.c_str(), 0, &data);
      }

      for (r=(*q).m_indices.begin();
           r!=(*q).m_indices.end(); ++r, ++loop) {

        // If scalar, name = window.attribute
        // If vector, name = window.#-attribute
        if( (*q).m_indices.size()==1)
          name = window + '.' + (*q).m_name;
        else {
          std::ostringstream sout;
          sout << window << '.' << loop << '-' << (*q).m_name;
          name = sout.str();
        }

        COM_get_array(name.c_str(), block->m_paneId,&data);

        bool nonempty;
        if (structured && ((*q).m_position=='n' || (*q).m_position=='e')) {
          start[0] = start[1] = start[2] = 0;
          if ((*q).m_position == 'n') {
            size[0] = block->m_gridInfo.front().m_size[0];
            size[1] = block->m_gridInfo.front().m_size[1];
            size[2] = block->m_gridInfo.front().m_size[2];
          } else {
            size[0] = std::max(int32(1),
                               block->m_gridInfo.front().m_size[0]-1);
            size[1] = std::max(int32(1),
                               block->m_gridInfo.front().m_size[1]-1);
            size[2] = std::max(int32(1),
                               block->m_gridInfo.front().m_size[2]-1);
          }

          nonempty = block->m_numNodes && data;

        } else {
          start[0] = 0;
          size[0] = (*q).m_nitems;

          nonempty = size[0] && data;
        }

        if ( nonempty) {
          HDF4_CHECK_RET(sds_id = HDF4::SDselect, (sd_id, *r), continue);
          HDF4_CHECK(HDF4::SDreaddata, (sds_id, start, NULL, size, data));
#ifdef DEBUG_DUMP_PREFIX
          {
            std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.' + block->time_level + ".hdf").c_str());
            DebugDump(fout, (*q).m_nitems, (*q).m_dataType, data);
          }
#endif // DEBUG_DUMP_PREFIX
          HDF4::SDendaccess(sds_id);
        }
      }
    }
    HDF4::SDend(sd_id);
  }
}

#ifdef USE_CGNS
/**
 ** Extract metadata from the list of files.
 **
 ** Open each file, scan the datasets, identify the windows, panes, and
 ** attributes.
 **/
static void scan_files_CGNS( int pathc, char* pathv[],
                             BlockMM_CGNS& blocks, std::string& time,
                             std::map<DataType_t, COM_Type>& CGNS2COM)
{
  bool has_timesteps = true; //Boeing fix
  double t = -1.0;
  std::string timeLevel;
  Block_CGNS* block;


  if (!time.empty()) {
    std::istringstream in(time);
    in >> t;

    std::ostringstream sout;
    sout << t;
    timeLevel = sout.str();
  }

  blocks.clear();

  // Iterate through each file.
  int i;
  for (i=0; i<pathc; ++i) {
    // Make sure this is a CGNS file.
    if (strcmp(&pathv[i][strlen(pathv[i])-5], ".cgns") != 0)
      continue;

    AutoCDer autoCD;
    std::string fname(pathv[i]);
    std::string::size_type cloc = fname.rfind('/');
    if (cloc != std::string::npos) {
      chdir(fname.substr(0, cloc).c_str());
      fname.erase(0, cloc + 1);
    }

    int fn;
    CG_CHECK_RET(cg_open, (fname.c_str(), MODE_READ, &fn), continue);
    AutoCloser<int> auto0(fn, cg_close);

    // Determine the number of bases.
    int nBases;
    CG_CHECK_RET(cg_nbases, (fn, &nBases), continue);

    // Iterate through the bases.
    char buffer[33];
    int B;
    for (B=1; B<=nBases; ++B) {
      // Get the name of the base/material/window, as well as the dimensions.
      // Physical dimensions should always be 3 in Roccom, but cell dims will
      // be 2 for surface meshes.
      char baseName[33];
      int cellDim, physDim;
      CG_CHECK_RET(cg_base_read, (fn, B, baseName, &cellDim, &physDim),
                   continue);

      // If the baseName ends with "_ridges" then ignore it.
      int nameLen = std::strlen(baseName);
      if (nameLen > 7 && std::strcmp(&baseName[nameLen-7], "_ridges") == 0)
        continue;

      // If the no time was given, use the first time that we find.
      if (timeLevel.empty()) {
        int nSteps;

        int error_code = cg_biter_read(fn, B, buffer, &nSteps); 
        if (error_code == CG_NODE_NOT_FOUND) {
          std::cerr << "Rocstar: Warning: cg_biter_read() failed with CG_NODE_NOT_FOUND"
                    << " on Base " << B << " in file " << pathv[i] << std::endl;
          std::cerr << "Rocstar: Warning: No time steps given with Base "<< B  <<std::endl;
          has_timesteps = false; //Boeing fix
        } else if (error_code == CG_ERROR) {
          std::cerr << "Rocstar: Warning: cg_biter_read() failed  with CG_ERROR on Base " << B
                    << " in file " << pathv[i] << std::endl;
          continue;
        }
        if (has_timesteps) {
          CG_CHECK(cg_goto, (fn, B, "BaseIterativeData_t", 1, "end"));

          std::vector<double> times(nSteps);

          CG_CHECK(cg_array_read_as, (1, RealDouble, &times[0]));

          int ti;
          for (ti=0; ti<nSteps; ++ti)
            DEBUG_MSG("Time " << ti << " == " << times[ti]);

          std::ostringstream sout;
          sout << times[0];
          timeLevel = sout.str();
          DEBUG_MSG("timeLevel == " << timeLevel);

          //WHAT'S THE STANDARD FORMAT FOR TIME STRINGS????
          sprintf(buffer, "%9.6f", times[0]);
          time = buffer;
          DEBUG_MSG("time == " << time);
        }
      }

      int nZones;
      CG_CHECK_RET(cg_nzones, (fn, B, &nZones), continue);


      // Create the name for the GridCoordinates_t node.
      // HACK: Assume no timesteps means Boeing CGNS file, and use
      // the "GridCoordinates" name. Otherwise default to Rocout
      // standard "Grid"
      // FIX: Should be a better way to determine which name to use 

      std::string gridAttrName;
      if (has_timesteps)
        {
          gridAttrName = "Grid";
          if (timeLevel.length() + gridAttrName.length() > 32)
            gridAttrName.erase(32 - timeLevel.length());
          gridAttrName += timeLevel;
        }
      else
        {
          gridAttrName = "GridCoordinates";
        }


      // Create the name for the IntegralData_t node for window attributes.
      std::string winAttrName("WinData"); 
      if (timeLevel.length() + winAttrName.length() > 32)
        winAttrName.erase(32 - timeLevel.length());
      winAttrName += timeLevel;

      // Create the name for the IntegralData_t node for pane attributes.
      std::string paneAttrName("PaneData");
      if (timeLevel.length() + paneAttrName.length() > 32)
        paneAttrName.erase(32 - timeLevel.length());
      paneAttrName += timeLevel;

      // Create the name for the IntegralData_t node for connectivity attrs.
      std::string connAttrName("ConnData");
      if (timeLevel.length() + connAttrName.length() > 32)
        connAttrName.erase(32 - timeLevel.length());
      connAttrName += timeLevel;

      // Create the name for the FlowSolution_t node for element attributes.
      std::string elemAttrName("ElemData");
      if (timeLevel.length() + elemAttrName.length() > 32)
        elemAttrName.erase(32 - timeLevel.length());
      elemAttrName += timeLevel;

      // Create the name for the FlowSolution_t node for node attributes.
      std::string nodeAttrName("NodeData");
      if (timeLevel.length() + nodeAttrName.length() > 32)
        nodeAttrName.erase(32 - timeLevel.length());
      nodeAttrName += timeLevel;

      int Z;
      for (Z=1; Z<=nZones; ++Z) {
        int nPoints, nGhostPoints;
        std::string units;
        char zoneName[33];
        int gridSize[3], coreSize[9];
        CG_CHECK_RET(cg_zone_read, (fn, B, Z, zoneName, coreSize), continue);

        ZoneType_t zoneType;
        CG_CHECK_RET(cg_zone_type, (fn, B, Z, &zoneType), continue);

        int nGrids;
        CG_CHECK_RET(cg_ngrids, (fn, B, Z, &nGrids), continue);

        int G, rind[6] = { 0, 0, 0, 0, 0, 0 };
        for (G=1; G<=nGrids; ++G) {
          CG_CHECK_RET(cg_grid_read, (fn, B, Z, G, buffer), continue);

          if (gridAttrName != buffer)
            continue;

          CG_CHECK_RET(cg_goto,
                       (fn, B, "Zone_t", Z, "GridCoordinates_t", G, "end"),
                       break);

          if (cg_rind_read(rind) != 0)
            std::fill_n(rind, 6, 0);

          if (zoneType == Structured) {
            gridSize[0] = coreSize[0] + rind[0] + rind[1];
            gridSize[cellDim] = coreSize[cellDim] + rind[0] + rind[1];
            nPoints = gridSize[0];
            if (cellDim > 1) {
              gridSize[1] = coreSize[1] + rind[2] + rind[3];
              gridSize[cellDim+1] = coreSize[cellDim+1] + rind[2] + rind[3];
              nPoints *= gridSize[1];
              if (cellDim > 2) {
                gridSize[2] = coreSize[2] + rind[4] + rind[5];
                gridSize[5] = coreSize[5] + rind[5] + rind[5];
                nPoints *= gridSize[2];
              }
            }
          } else {
            nPoints = gridSize[0] = coreSize[0] + rind[0] + rind[1];
          }
          // Assuming that ROCCOM ghost layers are uniform for structured
          // meshes, and ghost nodes are at the end for unstructured meshes.
          nGhostPoints = rind[1];

          // Assume that CoordinateX is array 1.
          CG_CHECK_RET(cg_goto,
                       (fn, B, "Zone_t", Z, "GridCoordinates_t", G, "DataArray_t", 1, "end"),
                       break);

          // Also assume that if one of the coordinate arrays is null, then
          // they all are null.
          int nDescriptors;
          CG_CHECK_RET(cg_ndescriptors, (&nDescriptors), break);

          int D;
          char* text = NULL;
          for (D=1; D<=nDescriptors; ++D) {
            CG_CHECK_RET(cg_descriptor_read, (D, buffer, &text), continue);

            if (strcmp(buffer, "Range") == 0) {
              if (strcmp(text, "EMPTY") == 0)
                nPoints = nGhostPoints = 0;
            } else if (strcmp(buffer, "Units") == 0) {
              units = text;
            }
            free(text);
          }

          break;
        }

        // Determine the pane id.  Skip any non-numeric prefix such as "block".
        std::string lbl(zoneName);
        std::string::size_type pos = lbl.find_first_of("0123456789");
        if (pos != std::string::npos)
          lbl.erase(0, pos);

        int paneId;
        {
          std::istringstream in(lbl); 
          in >> paneId;
        }

        // Create and initialize a Block object.
        block = new Block_CGNS(pathv[i], B, Z, G, paneId, timeLevel, units,
                               nPoints, nGhostPoints);

        if (zoneType == Structured) {
          block->m_gridInfo.push_back(GridInfo_CGNS(gridSize, cellDim,
                                                    nGhostPoints));
        } else {
          // Get info on the connectivity tables.
          std::string rocElemType;
          int nSections;
          CG_CHECK_RET(cg_nsections, (fn, B, Z, &nSections),
                       delete block; continue);

          int totalGhostElems = 0, S;
          for (S=1; S<=nSections; ++S) {
            ElementType_t elemType;
            int iStart, iEnd, iBoundry, pFlag;
            CG_CHECK_RET(cg_section_read,
                         (fn, B, Z, S, buffer, &elemType, &iStart, &iEnd, &iBoundry, &pFlag),
                         delete block; continue);

            int nElems = iEnd - iStart + 1;
            int nGhostElems = 0;
            std::string label = buffer;
            if (label.substr(0, 4) == "Null") {
              nElems = 0;
              label.erase(0, 4);
            } else if (label.substr(0, 5) == "Empty") {
              nElems = 0;
              label.erase(0, 5);
            } else {
              CG_CHECK_RET(cg_goto,
                           (fn, B, "Zone_t", Z, "Elements_t", S, "end"),
                           delete block; continue);

              int rind[2];
              if (cg_rind_read(rind) == 0)
                nGhostElems = rind[1];
            }
            // Fix for the Boeing fix
            // If the conn table name doesn't begin with a ':', it's not
            // a Rocout-generated file.
            if (label[0] != ':') {
              elementType_to_string(elemType, rocElemType);//Boeing Fix
              label = rocElemType + label;
            }

            totalGhostElems += nGhostElems;
            block->m_gridInfo.push_back(GridInfo_CGNS(label, nElems,
                                                      nGhostElems));
          }
          gridSize[1] += totalGhostElems;
        }

        // Get the window attributes.
        if (cg_goto(fn, B, "end") != 0) {
          std::cerr << "Rocstar: Warning: cg_goto() failed to go to Base " << B << " in file "
                    << pathv[i] << std::endl;
        } else {
          int nIntegrals;
          if (cg_nintegrals(&nIntegrals) != 0) {
            std::cerr << "Rocstar: Warning: cg_nintegrals() failed in Base " << B << " in file "
                      << pathv[i] << std::endl;
          } else {
            int I;
            for (I=1; I<=nIntegrals; ++I) {
              CG_CHECK_RET(cg_integral_read, (I, buffer), continue);

              if (winAttrName != buffer)
                continue;

              block->m_W = I;

              CG_CHECK_RET(cg_goto, (fn, B, "IntegralData_t", I, "end"), break);

              int nArrays;
              CG_CHECK_RET(cg_narrays, (&nArrays), break);

              int A;
              for (A=1; A<=nArrays; ++A) {
                CG_CHECK_RET(cg_goto, (fn, B, "IntegralData_t", I, "end"),
                             break);

                DataType_t dataType;
                int rank, size[3];
                CG_CHECK_RET(cg_array_info, (A, buffer, &dataType, &rank, size),
                             continue);

                CG_CHECK_RET(cg_goto, (fn, B, "IntegralData_t", I, "DataArray_t", A, "end"),
                             continue);

                std::string label = buffer;
                int nDescriptors;
                CG_CHECK_RET(cg_ndescriptors, (&nDescriptors), continue);

                char* text = NULL;
                std::string units;
                bool isNull = false;
                int nGhostItems = 0, D;
                for (D=1; D<=nDescriptors; ++D) {
                  CG_CHECK_RET(cg_descriptor_read, (D, buffer, &text),
                               continue);

                  if (strcmp(buffer, "Range") == 0) {
                    if (strcmp(text, "EMPTY") == 0) {
                      nGhostItems = size[0] = 0;
                      isNull = true;
                    } else if (strcmp(text, "NULL") == 0)
                      isNull = true;
                  } else if (strcmp(buffer, "Ghost") == 0) {
                    if (size[0] > 0) {
                      std::istringstream in(text);
                      in >> nGhostItems;
                    }
                  } else if (strcmp(buffer, "Units") == 0) {
                    units = text;
                  }
                  free(text);
                }

                int nComps = 1, comp = 0;
                std::string::size_type i = label.rfind('#');
                if (i != std::string::npos) {
                  std::istringstream in(&label[i+1]);
                  in >> comp;
                  --comp;
                  in.get(); in.get();
                  in >> nComps;
                  label.erase(i);
                } else {
                  i = label.length() - 1;
                  if (label[i] >= 'X' && label[i] <= 'Z') {
                    comp = label[i] - 'X';
                    label.erase(i);
                    --i;
                    if (label[i] >= 'X' && label[i] <= 'Z') {
                      comp += 3 * (label[i] - 'X');
                      label.erase(i);
                      nComps = 9;
                    } else {
                      nComps = 3;
                    }
                  }
                }

                if (comp != 0) {
                  VarInfo_CGNS &vi = block->m_variables.back();
                  vi.m_indices[comp] = A;
                  vi.m_is_null[comp] = isNull;
                } else {
                  block->m_variables.push_back
                    (VarInfo_CGNS(label, 'w', CGNS2COM[dataType], units,
                                  nComps, A, size[0], nGhostItems, isNull));
                }
              }
              break;
            }
          }
        }

        // Get the pane and connectivity attributes.
        if (cg_goto(fn, B, "Zone_t", Z, "end") != 0) {
          std::cerr << "Rocstar: Warning: cg_goto() failed to go to Zone " << Z << " of Base "
                    << B << " in file " << pathv[i] << std::endl;
        } else {
          int nIntegrals;
          if (cg_nintegrals(&nIntegrals) != 0) {
            std::cerr << "Rocstar: Warning: cg_nintegrals() failed in Zone " << Z << " of Base "
                      << B << " in file " << pathv[i] << std::endl;
          } else {
            int I, found = 0;
            for (I=1; I<=nIntegrals && found<2; ++I) {
              CG_CHECK_RET(cg_integral_read, (I, buffer), continue);

              char attrType;
              if (paneAttrName == buffer) {
                attrType = 'p';
                block->m_P = I;
              } else if (connAttrName == buffer) {
                attrType = 'c';
                block->m_C = I;
              } else
                continue;

              ++found;

              CG_CHECK_RET(cg_goto,
                           (fn, B, "Zone_t", Z, "IntegralData_t", I, "end"),
                           break);

              int nArrays;
              CG_CHECK_RET(cg_narrays, (&nArrays), break);

              int A;
              for (A=1; A<=nArrays; ++A) {
                CG_CHECK_RET(cg_goto,
                             (fn, B, "Zone_t", Z, "IntegralData_t", I, "end"),
                             break);

                DataType_t dataType;
                int rank, size[3];
                CG_CHECK_RET(cg_array_info, (A, buffer, &dataType, &rank, size),
                             continue);

                CG_CHECK_RET(cg_goto,
                             (fn, B, "Zone_t", Z, "IntegralData_t", I, "DataArray_t", A, "end"),
                             continue);

                std::string label = buffer;
                int nDescriptors;
                CG_CHECK_RET(cg_ndescriptors, (&nDescriptors), continue);

                char* text = NULL;
                std::string units;
                bool isNull = false;
                int nGhostItems = 0, D;
                for (D=1; D<=nDescriptors; ++D) {
                  CG_CHECK_RET(cg_descriptor_read, (D, buffer, &text),
                               continue);

                  if (strcmp(buffer, "Range") == 0) {
                    if (strcmp(text, "EMPTY") == 0) {
                      nGhostItems = size[0] = 0;
                      isNull = true;
                    } else if (strcmp(text, "NULL") == 0)
                      isNull = true;
                  } else if (strcmp(buffer, "Ghost") == 0) {
                    if (size[0] > 0) {
                      std::istringstream in(text);
                      in >> nGhostItems;
                    }
                  } else if (strcmp(buffer, "Units") == 0) {
                    units = text;
                  }
                  free(text);
                }

                int nComps = 1, comp = 0;
                std::string::size_type i = label.rfind('#');
                if (i != std::string::npos) {
                  std::istringstream in(&label[i+1]);
                  in >> comp;
                  --comp;
                  in.get(); in.get();
                  in >> nComps;
                  label.erase(i);
                } else {
                  i = label.length() - 1;
                  if (label[i] >= 'X' && label[i] <= 'Z') {
                    comp = label[i] - 'X';
                    label.erase(i);
                    --i;
                    if (label[i] >= 'X' && label[i] <= 'Z') {
                      comp += 3 * (label[i] - 'X');
                      label.erase(i);
                      nComps = 9;
                    } else {
                      nComps = 3;
                    }
                  }
                }

                if (comp != 0) {
                  VarInfo_CGNS &vi = block->m_variables.back();
                  vi.m_indices[comp] = A;
                  vi.m_is_null[comp] = isNull;
                } else {
//std::cout << "Units now == " << units << std::endl;
                  block->m_variables.push_back
                    (VarInfo_CGNS(label, attrType, CGNS2COM[dataType], units,
                                  nComps, A, size[0], nGhostItems, isNull));
                }
              }
            }
          }
        }

        // Get the element and node attributes.
        int nSols;
        if (cg_nsols(fn, B, Z, &nSols) != 0) {
          std::cerr << "Rocstar: Warning: cg_nsols() failed in Zone " << Z << " of Base "
                    << B << " in file " << pathv[i] << std::endl;
        } else {
          GridLocation_t loc;
          int S, found = 0;
          for (S=1; S<=nSols && found<2; ++S) {
            CG_CHECK_RET(cg_sol_info, (fn, B, Z, S, buffer, &loc), continue);

            bool isNodal = false;
            if (buffer == elemAttrName ) {
              block->m_E = S;
            } else if (buffer == nodeAttrName ) {
              block->m_N = S;
              isNodal = true;
            } else
              continue;

            ++found;

            CG_CHECK_RET(cg_goto,
                         (fn, B, "Zone_t", Z, "FlowSolution_t", S, "end"),
                         break);

            int rind[6];
            if (cg_rind_read(rind) != 0)
              std::fill_n(rind, 6, 0);

            int nItems;
            int start = isNodal ? 0 : 1;
            if (zoneType == Structured) {
              start *= cellDim;
              // Note: the rind values should always be the same.
              nItems = coreSize[start] + rind[0] + rind[1];
              if (cellDim > 1) {
                nItems *= coreSize[start+1] + rind[2] + rind[3];
                if (cellDim > 2)
                  nItems *= coreSize[start+2] + rind[4] + rind[5];
              }
            } else {
              // Note: rind[0] should always be 0.
              nItems = coreSize[start] + rind[0] + rind[1];
            }

            int nFields;
            CG_CHECK_RET(cg_nfields, (fn, B, Z, S, &nFields), break);

            int F;
            for (F=1; F<=nFields; ++F) {
              DataType_t dataType;
              CG_CHECK_RET(cg_field_info, (fn, B, Z, S, F, &dataType, buffer),
                           continue);

              CG_CHECK_RET(cg_goto,
                           (fn, B, "Zone_t", Z, "FlowSolution_t", S, "DataArray_t", F, "end"),
                           continue);

              std::string label = buffer;
              int nDescriptors;
              CG_CHECK_RET(cg_ndescriptors, (&nDescriptors), continue);

              char* text = NULL;
              std::string units;
              bool isNull = false;
              int D;
              for (D=1; D<=nDescriptors; ++D) {
                CG_CHECK_RET(cg_descriptor_read, (D, buffer, &text), continue);

                if (strcmp(buffer, "Range") == 0) {
                  if (strcmp(text, "EMPTY") == 0) {
                    rind[1] = nItems = 0;
                    isNull = true;
                  } else if (strcmp(text, "NULL") == 0)
                    isNull = true;
                } else if (strcmp(buffer, "Units") == 0)
                  units = text;
                free(text);
              }

              int nComps = 1, comp = 0;
              std::string::size_type i = label.rfind('#');
              if (i != std::string::npos) {
                std::istringstream in(&label[i+1]);
                in >> comp;
                --comp;
                in.get(); in.get();
                in >> nComps;
                label.erase(i);
              } else {
                i = label.length() - 1;
                if (label[i] >= 'X' && label[i] <= 'Z') {
                  comp = label[i] - 'X';
                  label.erase(i);
                  --i;
                  if (label[i] >= 'X' && label[i] <= 'Z') {
                    comp += 3 * (label[i] - 'X');
                    label.erase(i);
                    nComps = 9;
                  } else {
                    nComps = 3;
                  }
                }
              }

              if (comp != 0) {
                // We assume that components are contiguous and in order.
                VarInfo_CGNS &vi = block->m_variables.back();
                if (vi.m_name != label)
                  std::cerr << "Rocstar: Warning: attributes for " << label
                            << " are not contiguous." << std::endl;
                else if (vi.m_indices.size() != nComps)
                  std::cerr << "Rocstar: Warning: attribute " << label
                            << " has inconsistent component count."
                            << std::endl;
                vi.m_indices[comp] = F;
                vi.m_is_null[comp] = isNull;
              } else {
                block->m_variables.push_back
                  (VarInfo_CGNS(label, isNodal ? 'n' : 'e', CGNS2COM[dataType],
                                units, nComps, F, nItems, rind[1], isNull));
              }
            }
          }
        }

        // Check for old node attributes.
        int nDiscrete;
        if (cg_ndiscrete(fn, B, Z, &nDiscrete) != 0) {
          std::cerr << "Rocstar:Warning: cg_ndiscrete() failed in Zone " << Z << " of Base "
                    << B << " in file " << pathv[i] << std::endl;
        } else {
          int S;
          for (S=1; S<=nDiscrete; ++S) {
            CG_CHECK_RET(cg_discrete_read, (fn, B, Z, S, buffer), continue);

            if (nodeAttrName != buffer)
              continue;

            block->m_N = S;

            CG_CHECK_RET(cg_goto,
                         (fn, B, "Zone_t", Z, "DiscreteData_t", S, "end"),
                         break);

            int rind[6];
            if (cg_rind_read(rind) != 0)
              std::fill_n(rind, 6, 0);

            int nItems;
            if (zoneType == Structured) {
              // Note: the rind values should always be the same.
              nItems = coreSize[0] + rind[0] + rind[1];
              if (cellDim > 1) {
                nItems *= coreSize[1] + rind[2] + rind[3];
                if (cellDim > 2)
                  nItems *= coreSize[2] + rind[4] + rind[5];
              }
            } else {
              // Note: rind[0] should always be 0.
              nItems = coreSize[0] + rind[0] + rind[1];
            }

            int nArrays;
            CG_CHECK_RET(cg_narrays, (&nArrays), break);

            int A;
            for (A=1; A<=nArrays; ++A) {
              CG_CHECK_RET(cg_goto,
                           (fn, B, "Zone_t", Z, "DiscreteData_t", S, "end"),
                           break);

              DataType_t dataType;
              int rank, size[3];
              CG_CHECK_RET(cg_array_info, (A, buffer, &dataType, &rank, size),
                           continue);
              std::string label = buffer;

              // 'rank' should equal the cell dimension for structured,
              // 1 for unstructured grids.  'size' should be the core size
              // plus rind layers.
              COM_assertion_msg(rank = (zoneType == Structured ? cellDim : 1),
                                ("Bad rank for nodal attribute '"
                                 + label + "'\n").c_str());
              // TODO: Verify the size, too.

              CG_CHECK_RET(cg_goto,
                           (fn, B, "Zone_t", Z, "DiscreteData_t", S, "DataArray_t", A, "end"),
                           continue);

              int nDescriptors;
              CG_CHECK_RET(cg_ndescriptors, (&nDescriptors), continue);

              char* text = NULL;
              std::string units;
              bool isNull = false;
              int D;
              for (D=1; D<=nDescriptors; ++D) {
                CG_CHECK_RET(cg_descriptor_read, (D, buffer, &text), continue);

                if (strcmp(buffer, "Range") == 0) {
                  if (strcmp(text, "EMPTY") == 0) {
                    rind[1] = nItems = 0;
                    isNull = true;
                  } else if (strcmp(text, "NULL") == 0)
                    isNull = true;
                } else if (strcmp(buffer, "Units") == 0)
                  units = text;
                free(text);
              }

              int nComps = 1, comp = 0;
              std::string::size_type i = label.rfind('#');
              if (i != std::string::npos) {
                std::istringstream in(&label[i+1]);
                in >> comp;
                --comp;
                in.get(); in.get();
                in >> nComps;
                label.erase(i);
              } else {
                i = label.length() - 1;
                if (label[i] >= 'X' && label[i] <= 'Z') {
                  comp = label[i] - 'X';
                  label.erase(i);
                  --i;
                  if (label[i] >= 'X' && label[i] <= 'Z') {
                    comp += 3 * (label[i] - 'X');
                    label.erase(i);
                    nComps = 9;
                  } else {
                    nComps = 3;
                  }
                }
              }

              if (comp != 0) {
                // We assume that components are contiguous and in order.
                VarInfo_CGNS &vi = block->m_variables.back();
                vi.m_indices[comp] = A;
                vi.m_is_null[comp] = isNull;
              } else {
                block->m_variables.push_back
                  (VarInfo_CGNS(label, 'n', CGNS2COM[dataType], units,
                                nComps, A, nItems, rind[1], isNull));
              }
            }

            break;
          }
        }

        blocks.insert(BlockMM_CGNS::value_type(baseName, block));
      }
    }
  }
}

static void load_data_CGNS( BlockMM_CGNS::iterator p,
                            const BlockMM_CGNS::iterator& end,
                            const std::string& window, 
                            const MPI_Comm* comm, int rank, int nprocs)
{
  int i;
  Block_CGNS* block;
  std::string name;
  void* data;
  int with_pane;

  for ( with_pane=0; p!=end; ++p, ++with_pane) {
    block = (*p).second;

    // Panes are a local construct, so make sure that this pane is local.
    bool local = COM_get_status( window.c_str(), block->m_paneId) >= 0;

    if (!local) continue;

    AutoCDer autoCD;
    std::string fname(block->m_file);
    std::string::size_type cloc = fname.rfind('/');
    if (cloc != std::string::npos) {
      chdir(fname.substr(0, cloc).c_str());
      fname.erase(0, cloc + 1);
    }

    int fn;
    CG_CHECK_RET(cg_open, (fname.c_str(), MODE_READ, &fn), continue);
    AutoCloser<int> auto0(fn, cg_close);

    bool structured = (!block->m_gridInfo.empty() && 
                       block->m_gridInfo.front().m_name.substr(0,3) == ":st");

    // Let Roccom manage memory for nodal coordinates
    // Moved allocation outside of dimension loop.  Not sure if this fully
    // accounts for single registration of "nc" versus "1-nc", "2-nc", etc.

    name = window + ".nc";
    DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == " << block->m_paneId << ", ptr == NULL, arg == 1 )");
    COM_resize_array(name.c_str(), block->m_paneId, NULL, 1);

    // Read the nodal coordinates.
    for (i=1; i<4; ++i) {
      // Read the mesh only if number of nodes is positive
      name = window + '.' + (char)('0' + i) + "-nc";
      COM_get_array(name.c_str(), block->m_paneId, &data);

      if (block->m_numNodes > 0 && data != NULL) {
        CG_CHECK_RET(cg_goto,
                     (fn, block->m_B, "Zone_t", block->m_Z,
                      "GridCoordinates_t", block->m_G, "end"),
                     continue);

        CG_CHECK(cg_array_read, (i, data));
#ifdef DEBUG_DUMP_PREFIX
        {
          std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.'
                              + block->time_level + ".cgns").c_str());
          DebugDump(fout, block->m_numNodes, COM_DOUBLE, data);
        }
#endif // DEBUG_DUMP_PREFIX
      }
    }

    // Read the connectivity tables.
    if (!structured) {
      std::vector<GridInfo_CGNS>::iterator s;
      for (i=1,s=block->m_gridInfo.begin();s!=block->m_gridInfo.end();++i,++s) {
        name = window + "." + (*s).m_name;
        DEBUG_MSG("Calling COM_resize_array( name == '" << name
                  << "', paneid == " << block->m_paneId << ", ptr == "
                  << &data << ", arg == 1 )");
        COM_resize_array(name.c_str(), block->m_paneId, &data, 1);

        if ((*s).m_numElements > 0 && data != NULL) {
          int nn;
          std::istringstream in(&(*s).m_name[2]);
          in >> nn;
          std::vector<int> cbuf((*s).m_numElements * nn);
          CG_CHECK_RET(cg_elements_read,
                       (fn, block->m_B, block->m_Z, i, &cbuf[0], NULL),
                       continue);

          // Scramble the conn table the way Roccom likes it.
          int elem, node, zz;
          for (elem=0,zz=0; elem<(*s).m_numElements; ++elem)
            for (node=0; node<nn; ++node,++zz)
              ((int*)data)[node*(*s).m_numElements+elem] = cbuf[zz];
#ifdef DEBUG_DUMP_PREFIX
          {
            std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.' + block->time_level + ".cgns").c_str());
            DebugDump(fout, (*s).m_numElements * nn, COM_INT, data);
          }
#endif // DEBUG_DUMP_PREFIX
        }
      }
    }

    // Read the attribute data.
    std::vector<VarInfo_CGNS>::iterator q;
    for (q=block->m_variables.begin(); q!=block->m_variables.end(); ++q) {
      // Skip NULL attributes.
      if ( (*q).m_is_null[0]) continue;
      // Read in window attribute only for the first pane.
      if ( with_pane && (*q).m_position=='w') continue;


      name = window + '.' + (*q).m_name;
      if ( (*q).m_position!='w')  { // Allocate for nonwindow attribute
        DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == " << block->m_paneId << ", ptr == " << &data << ", arg == 1 )");
        COM_resize_array(name.c_str(), block->m_paneId, &data, 1);
      } else {
        COM_get_array(name.c_str(), 0, &data);
      }

      switch((*q).m_position) {
        case 'w':
          CG_CHECK_RET(cg_goto,
                       (fn, block->m_B, "IntegralData_t", block->m_W, "end"),
                       continue);
          break;

        case 'p':
          CG_CHECK_RET(cg_goto,
                       (fn, block->m_B, "Zone_t", block->m_Z, "IntegralData_t", block->m_P, "end"),
                       continue);
          break;

        case 'c':
          CG_CHECK_RET(cg_goto,
                       (fn, block->m_B, "Zone_t", block->m_Z, "IntegralData_t", block->m_C, "end"),
                       continue);
          break;

        case 'e':
          CG_CHECK_RET(cg_goto,
                       (fn, block->m_B, "Zone_t", block->m_Z, "FlowSolution_t", block->m_E, "end"),
                       continue);
          break;

        case 'n':
          CG_CHECK_RET(cg_goto,
                       (fn, block->m_B, "Zone_t", block->m_Z, "FlowSolution_t", block->m_N, "end"),
                       continue);
          break;
      }

      int loop;
      std::vector<int>::iterator r;
      for (r=(*q).m_indices.begin(),loop=1; r!=(*q).m_indices.end();
           ++r,++loop) {

        // If scalar, name = window.attribute
        // If vector, name = window.#-attribute
        if( (*q).m_indices.size() == 1)
          name = window + '.' + (*q).m_name;
        else {
          std::ostringstream sout;
          sout << window << '.' << loop << '-' << (*q).m_name;
          name = sout.str();
        }

        COM_get_array(name.c_str(), block->m_paneId, &data);

        if ((*q).m_nitems == 0 || data == NULL) continue;

        CG_CHECK(cg_array_read, (*r, data));
#ifdef DEBUG_DUMP_PREFIX
        {
          std::ofstream fout((DEBUG_DUMP_PREFIX + name + '.' + block->time_level + ".cgns").c_str());
          DebugDump(fout, (*q).m_nitems, (*q).m_dataType, data);
        }
#endif // DEBUG_DUMP_PREFIX
      }
    }
  }
}
#endif // USE_CGNS

static void new_attributes(BlockMM_HDF4::iterator hdf4,
                           const BlockMM_HDF4::iterator& hdf4End,
#ifdef USE_CGNS
                           BlockMM_CGNS::iterator cgns,
                           const BlockMM_CGNS::iterator& cgnsEnd,
#endif // USE_CGNS
                           const std::string& window,
                           const MPI_Comm* comm, int rank, int nprocs)
{

  // Build an MPI message to synchronize variable names.
  std::ostringstream sout;
  std::set<std::string> added;

  while (hdf4 != hdf4End) {
    std::vector<VarInfo_HDF4> &vars = hdf4->second->m_variables;
    std::vector<VarInfo_HDF4>::const_iterator q;

    for (q=vars.begin(); q!=vars.end(); ++q){
      if (added.count((*q).m_name) == 0) {
        sout << (*q).m_name << '!' << (*q).m_position << '!'
             << (*q).m_dataType << '!' << (*q).m_indices.size()
             << '!' << (*q).m_units << '!';
        if ( (*q).m_position=='w')
          sout << (*q).m_nitems << '!' << (*q).m_ng;
        sout << '|';
        added.insert((*q).m_name);
      }
    }
    ++hdf4;
  }

#ifdef USE_CGNS
  while (cgns != cgnsEnd) {
    std::vector<VarInfo_CGNS> &vars = cgns->second->m_variables;
    std::vector<VarInfo_CGNS>::const_iterator q;

    for (q=vars.begin(); q!=vars.end(); ++q){
      if (added.count((*q).m_name) == 0) {
        sout << (*q).m_name << '!' << (*q).m_position << '!'
             << (*q).m_dataType << '!' << (*q).m_indices.size()
             << '!' << (*q).m_units << '!';
        if ( (*q).m_position=='w')
          sout << (*q).m_nitems << '!' << (*q).m_ng;
        sout << '|';
        added.insert((*q).m_name);
      }
    }
    ++cgns;
  }
#endif // USE_CGNS

  // Get the lengths of all of the messages.
  std::string msg = sout.str();
  int len = msg.length();
  std::vector<int> lengths(nprocs);
  if ( *comm != MPI_COMM_NULL)
    MPI_Allgather(&len, 1, MPI_INT, &lengths[0], 1, MPI_INT, *comm);
  else
    lengths[0] = len;

  // Compute the total length and displacements for each process's message.
  int globlen = lengths[0];
  std::vector<int> disp(nprocs, 0);
  for (int i=1; i<nprocs; ++i) {
    globlen += lengths[i];
    disp[i] = disp[i-1] + lengths[i-1];
  }

  // Give each process a huge glob with every attribute of every process.
  std::vector<char> buffer(len+1,'\0');
  std::vector<char> glob(globlen+1,'\0');
  std::strcpy( &buffer[0], msg.c_str());
  if ( *comm != MPI_COMM_NULL)
    MPI_Allgatherv(&buffer[0], len, MPI_CHAR, &glob[0], &lengths[0], &disp[0],
                   MPI_CHAR, *comm);
  else
    glob = buffer;

//std::cout << "glob == \"" << &glob[0] << '"' << std::endl;
  std::string name;
  char buf[256];
  char position;
  COM_Type dType;
  int numComponents;
  char* token = &glob[0];
  added.clear();
  // printf("parsing string: '%s'\n",token);
  while (1) {
    char *nextToken = strchr(token,'|');
    if (nextToken==0) break;
    *nextToken++=0; /* eliminate and advance over separator */
    // printf("parsed token: '%s'\n",token);
    std::istringstream sin(token);
    sin.getline(buf, sizeof(buf), '!');
    if (added.count(buf) == 0) {
      added.insert(buf);
      name = window + '.';
      name += buf;
      position = sin.get(); sin.get();
      sin >> dType; sin.get();
      sin >> numComponents; sin.get();
      sin.getline(buf, sizeof(buf), '!');
      DEBUG_MSG("Calling COM_new_attribute( name == '" << name << "', position == '" << position << "', datatype == " << dType << ", nComp == " << numComponents << ", units == '" << buf << "' )");
      COM_new_attribute(name.c_str(), position, dType, numComponents, buf);

      if ( position == 'w') { // Set size and allocate for window attributes
        int nitems, ng;
        sin >> nitems; sin.get();
        sin >> ng; sin.get();
        DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == 0, nitems == " << nitems << ", nghost == " << ng << " )");
        COM_set_size(name.c_str(), 0, nitems, ng);
        DEBUG_MSG("Calling COM_resize_array( name == '" << name << "', paneid == 0 )");
        COM_resize_array(name.c_str(), 0);
      }
    }
    token = nextToken;
  }
}

// Broadcast window attributes if a process does not have any pane.
static void broadcast_win_attributes( bool isEmpty, const std::string& window, 
                                      const MPI_Comm *comm, int rank,
                                      int nprocs)
{
  // Build an MPI message to synchronize variable names.
  int empty = isEmpty ? 1 : 0;
  std::vector<int> is_empty(nprocs);
  if ( *comm != MPI_COMM_NULL)
    MPI_Allgather(&empty, 1, MPI_INT, &is_empty[0], 1, MPI_INT, *comm);
  else
    is_empty[0] = empty;

  // Broadcast from the nonempty processor with highest rank.
  int last_empty=-1, root=-1;
  for (int i=0; i<nprocs; ++i) {
    if ( is_empty[i]) last_empty = i;
    else root = i;
  }

  if ( last_empty<0) return; // No need to broadcast.

  // Obtain the list of attributes
  int na; char *atts;
  COM_get_attributes( window.c_str(), &na, &atts);

  // Loop through the attributes
  std::istringstream is(atts);
  for ( int i=0; i<na; ++i) {
    // Obtain the attribute name
    std::string aname;  is >> aname; 
    char loc;
    int  type, ncomp;

    COM_get_attribute( (window+"."+aname).c_str(), &loc, &type, &ncomp, NULL);
    if ( loc !='w') continue; // Skip non-window attributes;

    void *addr;
    int count, strd;
    COM_get_array( (window+"."+aname).c_str(), 0, &addr, &strd, &count);

    int len = COM_get_sizeof( type, count*ncomp);
    if ( *comm != MPI_COMM_NULL)
      MPI_Bcast( addr, len, MPI_CHAR, root, *comm);
  }
  COM_free_buffer( &atts);
}

void Rocin::register_panes(BlockMM_HDF4::iterator hdf4,
                           const BlockMM_HDF4::iterator& hdf4End,
#ifdef USE_CGNS
                           BlockMM_CGNS::iterator cgns,
                           const BlockMM_CGNS::iterator& cgnsEnd,
#endif // USE_CGNS
                           const std::string& window, RulesPtr is_local,
                           const MPI_Comm* comm, int rank, int nprocs)
{
  int local;
  std::string name;
  bool is_first=true;
  for ( ; hdf4!=hdf4End; ++hdf4) {
    Block_HDF4* block = (*hdf4).second;

    // Panes are a local construct, so make sure that this pane is local.
    if ( m_is_local)
      (this->*m_is_local)(block->m_paneId, rank, nprocs, &local);
    else if (is_local)
      is_local(block->m_paneId, rank, nprocs, &local);
    else
      local = 1;

    if (!local) continue;

    // Register the mesh unit & number of nodes
    name = window + ".nc";
    if (is_first && !block->m_units.empty()) {

      DEBUG_MSG("Calling COM_new_attribute( name == '" << name << "', position == 'n', datatype == COM_DOUBLE, nComp == 3, units == '" << block->m_units << "' )");
      COM_new_attribute(name.c_str(),'n',COM_DOUBLE,3,block->m_units.c_str());
      is_first = false;
    }

    DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << block->m_numNodes << ", nghost == " << block->m_numGhostNodes << " )");
    COM_set_size(name.c_str(), block->m_paneId,
                 block->m_numNodes,block->m_numGhostNodes);      

    // Register the connectivity tables or the mesh dimensions
    // need to specify attribute names for multiple connectivity tables
    if (block->m_gridInfo.size() && 
        block->m_gridInfo.front().m_name.substr(0,3) == ":st") {
      int mysize[3], ndim = block->m_gridInfo.front().m_name[3] - '0';
      mysize[0] = block->m_gridInfo.front().m_size[2];
      mysize[1] = block->m_gridInfo.front().m_size[1];
      mysize[2] = block->m_gridInfo.front().m_size[0];
      name = window + '.' + block->m_gridInfo.front().m_name;

      DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << ndim << ", nghost == " << block->m_gridInfo.front().m_numGhostElements << " )");
      COM_set_size(name.c_str(), block->m_paneId, ndim,
                   block->m_gridInfo.front().m_numGhostElements);
      COM_set_array(name.c_str(),block->m_paneId,mysize);
    } else {
      std::vector<GridInfo_HDF4>::iterator q;
      for (q=block->m_gridInfo.begin(); q!=block->m_gridInfo.end(); ++q){
        name = window + "." + (*q).m_name;
        DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << (*q).m_numElements << ", nghost == " << (*q).m_numGhostElements << " )");
        COM_set_size(name.c_str(), block->m_paneId,
                     (*q).m_numElements, (*q).m_numGhostElements);
      }
    }

    // Set sizes for pane attributes.
    std::vector<VarInfo_HDF4> &vars = block->m_variables;
    std::vector<VarInfo_HDF4>::const_iterator q;

    for ( q=vars.begin(); q!=vars.end(); ++q) {
      if ( (*q).m_position == 'p' || (*q).m_position == 'c') {
        name = window + "." + (*q).m_name;
        DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << (*q).m_nitems << ", nghost == " << (*q).m_ng << " )");
        COM_set_size( name.c_str(), block->m_paneId, 
                      (*q).m_nitems, (*q).m_ng);
      }
    }
  }

#ifdef USE_CGNS
  for ( ; cgns!=cgnsEnd; ++cgns) {
    Block_CGNS* block = (*cgns).second;

    // Panes are a local construct, so make sure that this pane is local.
    if ( m_is_local)
      (this->*m_is_local)(block->m_paneId, rank, nprocs, &local);
    else if (is_local)
      is_local(block->m_paneId, rank, nprocs, &local);
    else
      local = 1;

    if (!local) continue;

    // Register the mesh unit & number of nodes
    name = window + ".nc";
    if (is_first && !block->m_units.empty()) {

      DEBUG_MSG("Calling COM_new_attribute( name == '" << name << "', position == 'n', datatype == COM_DOUBLE, nComp == 3, units == '" << block->m_units << "' )");
      COM_new_attribute(name.c_str(),'n',COM_DOUBLE,3,block->m_units.c_str());
      is_first = false;
    }

    DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << block->m_numNodes << ", nghost == " << block->m_numGhostNodes << " )");
    COM_set_size(name.c_str(), block->m_paneId,
                 block->m_numNodes,block->m_numGhostNodes);      

    // Register the connectivity tables or the mesh dimensions
    // need to specify attribute names for multiple connectivity tables
    if (block->m_gridInfo.size() && 
        block->m_gridInfo.front().m_name.substr(0,3) == ":st") {
      int ndim = block->m_gridInfo.front().m_name[3] - '0';
      name = window + '.' + block->m_gridInfo.front().m_name;

      DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << ndim << ", nghost == " << block->m_gridInfo.front().m_numGhostElements << " )");
      COM_set_size(name.c_str(), block->m_paneId, ndim,
                   block->m_gridInfo.front().m_numGhostElements);
      COM_set_array(name.c_str(), block->m_paneId,
                    block->m_gridInfo.front().m_size);
    } else {
      std::vector<GridInfo_CGNS>::iterator q;
      for (q=block->m_gridInfo.begin(); q!=block->m_gridInfo.end(); ++q){
        name = window + "." + (*q).m_name;
        DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << (*q).m_numElements << ", nghost == " << (*q).m_numGhostElements << " )");
        COM_set_size(name.c_str(), block->m_paneId,
                     (*q).m_numElements, (*q).m_numGhostElements);
      }
    }

    // Set sizes for pane attributes.
    std::vector<VarInfo_CGNS> &vars = block->m_variables;
    std::vector<VarInfo_CGNS>::const_iterator q;

    for ( q=vars.begin(); q!=vars.end(); ++q) {
      if ( (*q).m_position == 'p' || (*q).m_position == 'c') {
        name = window + "." + (*q).m_name;
        DEBUG_MSG("Calling COM_set_size( name == '" << name << "', paneid == " << block->m_paneId << ", nitems == " << (*q).m_nitems << ", nghost == " << (*q).m_ng << " )");
        COM_set_size( name.c_str(), block->m_paneId, 
                      (*q).m_nitems, (*q).m_ng);
      }
    }
  }
#endif // USE_CGNS
}

template<class BLOCK>
void free_blocks(BLOCK& blocks)
{
  typename BLOCK::iterator p;
  for (p=blocks.begin(); p!=blocks.end(); ++p)
    delete (*p).second;

  blocks.clear();
}

void Rocin::init(const std::string &mname) {
  HDF4::init();

  Rocin *rin = new Rocin();

  /// This data structure maps HDF4 data types to roccom data types.
  rin->m_HDF2COM[DFNT_CHAR8] = COM_CHAR;
  rin->m_HDF2COM[DFNT_INT32] = COM_INT;
  rin->m_HDF2COM[DFNT_FLOAT32] = COM_FLOAT;
  rin->m_HDF2COM[DFNT_FLOAT64] = COM_DOUBLE;

#ifdef USE_CGNS
  /// This data structure maps CGNS data types to roccom data types.
  rin->m_CGNS2COM[Character] = COM_CHAR;
  rin->m_CGNS2COM[Integer] = COM_INT;
  rin->m_CGNS2COM[RealSingle] = COM_FLOAT;
  rin->m_CGNS2COM[RealDouble] = COM_DOUBLE;
#endif // USE_CGNS

  COM_new_window( mname.c_str(), MPI_COMM_SELF);

  std::string glb=mname+".global";

  DEBUG_MSG("Calling COM_new_attribute( name == '" << glb << "', position == 'w', datatype == COM_OBJECT, nComp == 1, units == '' )");
  COM_new_attribute( glb.c_str(), 'w', COM_OBJECT, 1, "");
  COM_set_object( glb.c_str(), 0, rin);

  // Register the function read_windows
  COM_Type types[8] = { COM_RAWDATA, COM_STRING, COM_STRING, COM_STRING,
                        COM_MPI_COMM, COM_VOID, COM_STRING, COM_INT};
  COM_set_member_function( (mname+".read_windows").c_str(),
                           (Member_func_ptr)&Rocin::read_windows, 
                           glb.c_str(), "biiIIIBI", types);

  // Register the function read_window
  types[3] = COM_MPI_COMM; types[4] = COM_VOID;
  types[5] = COM_STRING; types[6] = COM_INT;
  COM_set_member_function( (mname+".read_window").c_str(),
                           (Member_func_ptr)&Rocin::read_window, 
                           glb.c_str(), "biiIIBI", types);

  // Register the function read_by_control_file
  types[4] = COM_STRING; types[5] = COM_INT;
  COM_set_member_function( (mname+".read_by_control_file").c_str(),
                           (Member_func_ptr)&Rocin::read_by_control_file, 
                           glb.c_str(), "biiIBI", types);

  // Register the function obtain_attribute
  types[1] = COM_METADATA;
  types[2] = COM_METADATA;
  types[3] = COM_INT;
  COM_set_member_function( (mname+".obtain_attribute").c_str(),
                           (Member_func_ptr)&Rocin::obtain_attribute, 
                           glb.c_str(), "bioI", types);

  // Regsister the function read_parameter_file
  types[1] = COM_STRING;
  types[2] = COM_STRING;
  types[3] = COM_MPI_COMM;
  COM_set_member_function( (mname+".read_parameter_file").c_str(),
                           (Member_func_ptr)&Rocin::read_parameter_file,
                           glb.c_str(), "biiI", types);

  COM_window_init_done( mname.c_str());
}

void Rocin::finalize(const std::string &mname) {
  // Retrieve Rocin object from Roccom
  Rocin *rin;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &rin);

  COM_delete_window( mname.c_str());

  // Delete the object
  delete rin;
  HDF4::finalize();
}

void Rocin::obtain_attribute( const COM::Attribute *attribute_in,
                              COM::Attribute *user_attribute,
                              int *pane_id)
{
  // obtain a valid attribute object from user_attribute's name
  COM_assertion_msg((attribute_in != NULL && user_attribute != NULL),
                    "Null attributes are not valid arguments to Rocin::obtain_attribute\n");

  if ( attribute_in != user_attribute) {
    COM::Window *win = user_attribute->window();
    win->inherit(const_cast<Attribute*>(attribute_in), user_attribute->name(), 
                 COM::Pane::INHERIT_COPY, true, NULL, pane_id?*pane_id:0);
  }
}

void Rocin::explicit_local( const int& pid, const int& comm_rank,
                            const int& comm_size, int* il)
{
  *il = m_pane_ids.count(pid);
}

void Rocin::cyclic_local(const int& pid, const int& comm_rank,
                         const int& comm_size, int* il)
{
  *il = (pid - m_offset) % comm_size == comm_rank;
}

void Rocin::blockcyclic_local(const int& pid, const int& comm_rank,
                              const int& comm_size, int* il)
{
  int proc = (pid-m_offset) / m_base;
  if ( proc>=comm_size) proc %= comm_size;
  *il = proc == comm_rank;
}

void Rocin::read_by_control_file( const char* control_file_name,
                                  const char* window_name,
                                  const MPI_Comm* comm,
                                  char* time_level,
                                  const int* str_len)
{
  const MPI_Comm default_comm=COM_get_default_communicator();
  const MPI_Comm comm_null=MPI_COMM_NULL;
  const MPI_Comm* myComm;
  int myRank, comm_size;

  // Determine the directory name of the control_file_name
  std::string dir_name = control_file_name;
  std::string::size_type pos = dir_name.find_last_of('/');
  if ( pos == std::string::npos)
    dir_name = "";
  else
    dir_name.erase( pos+1);

  if ( comm && *comm == MPI_COMM_NULL || !comm) 
    myComm = COMMPI_Initialized() ? &default_comm : &comm_null;
  else
    myComm = comm;

  // Set myRank to a wild card if MPI is not initialized
  if ( !COMMPI_Initialized() || *myComm == MPI_COMM_NULL) {
    myRank = -1;
    comm_size = 1;
  }
  else {
    MPI_Comm_rank(*myComm, &myRank);
    MPI_Comm_size(*myComm, &comm_size);
    if(comm_size == 1)
      myRank = -1;
  }

  std::ifstream fin(control_file_name);
  if (!fin.is_open()) {
    std::cerr << "Rocstar: Error: read_by_control_file unable to open " << control_file_name
              << " for reading." << std::endl;
    return;
  }

  std::vector<std::string> patterns;
  std::string buffer;
  int rank = -54321;
  m_pane_ids.clear();

  while (!fin.eof()) {
    std::vector<std::string> local_patterns;
    while (!fin.eof() && (myRank<0 || rank != myRank) ) {
      if (buffer != "@Proc:") {
        fin >> buffer;
        continue;
      }

      while (fin.peek() == ' ' || fin.peek() == '\t' || fin.peek() == '\n')
        fin.get();

      if (fin.peek() == '*') {
        fin.get();
        rank = myRank;
        break;
      } else {
        fin >> rank;
        // If myRank is a wild card, then any rank will match with it.
        if (myRank < 0) {
          rank = myRank;
          break;
        } else if (rank == myRank)
          break;
      }

      buffer=""; 
    }

    if (fin.eof()) {
      if ( rank != myRank)
        std::cerr << "Rocstar: Error (read_by_control_file): control file " 
                  << control_file_name
                  << " does not contain information for process " << myRank
                  << std::endl;
      break;
    }

    fin >> buffer;
    if (buffer != "@Files:") {
      std::cerr << "Rocstar: Error (read_by_control_file): in control file "
                << control_file_name
                << ": expected '@Files' but found '" << buffer << '\''
                << std::endl;
      return;
    }

    fin >> buffer;

    while (buffer != "@Panes:" && buffer != "@Proc:") {
      // The @Panes section may not be present
      if (fin.eof()) break;

      std::string::size_type pos = buffer.find("%t");
      if (pos != std::string::npos)
        buffer.replace(pos, 2, time_level);

      pos = buffer.find('%');
      if ( pos != std::string::npos) {
        std::string::size_type pos_key = 
          buffer.find_first_not_of("0123456789", pos+1);

        COM_assertion_msg( pos_key != std::string::npos, 
                           (std::string("Incomplete placeholder in file name ")+buffer).c_str());

        std::string width = buffer.substr(pos+1, pos_key-pos-1);

        int w = 4; // The default rankwidth is 4, if width is empty     
        if ( !width.empty()) {
          std::istringstream sin(width);
          sin >> w;
        }

        std::ostringstream sout;

        switch( buffer[pos_key]) {
        case 'p': {
          // If myRank<0, then match the rank with any integer of width w
          if ( myRank<0)
            for ( int i=0; i<w; ++i) sout << "[0-9]";
          else
            sout << std::setfill('0') << std::setw(w) << myRank;

          break;
        }
        case 'i':
        case 'b': {
          // TODO: Use the mapping rule on the next line to distribute files.
          // For now, have all processes read all files that match the number
          // of digits for the block or pane ID in the filenames.
          for ( int i=0; i<w; ++i) sout << "[0-9]";
          break;
        }
        default:
          COM_assertion_msg( false, (std::string("Unknown keyword in file name ")+buffer).c_str());
        }

        buffer.replace(pos, pos_key - pos + 1, sout.str());
      }

      // Determine whether buffer has directory part. If not, prepend that
      // of the control file. If so, then use it as is.
      if ( !dir_name.empty()) {
/*
        std::string::size_type pos = buffer.find_last_of('/');
        if ( pos == std::string::npos) buffer.insert(0, dir_name);
*/
          // when it is not absolute path and does not have same prefix
        if (buffer[0] != '/' &&
             strncmp(buffer.c_str(), dir_name.c_str(), std::min(buffer.size(), dir_name.size())))   buffer.insert(0, dir_name);
      }

      // Insert buffer into patterns.
      local_patterns.push_back(buffer);

      fin >> buffer;
    }

    // The @Panes section may not be present
    if (buffer == "@Proc:") {
      patterns.insert(patterns.end(), local_patterns.begin(),
                      local_patterns.end());
      continue;
    }

    m_base = 0;
    m_offset = 0;
    m_is_local = NULL;

    if ( !fin.eof()) 
      fin >> buffer;
    else
      buffer = "-1";

    // Read in panes.
    if (buffer == "@Cyclic") {
      fin >> m_offset;
      m_is_local = &Rocin::cyclic_local;
    } else if (buffer == "@BlockCyclic") {
      fin >> m_base; 
      fin >> m_offset;
      m_is_local = &Rocin::blockcyclic_local;
    } else if (buffer == "@Block" || buffer == "@BlockBlockCyclic") {
      int block;
      fin >> block;
      fin >> m_base; 
      fin >> m_offset;
      int quot = block/comm_size, rem = block-quot*comm_size;
      if ( rank < rem) { 
        m_base *= quot+1;
      }
      else {
        m_offset += rem;
        m_base *= quot;
      }
      m_is_local = &Rocin::blockcyclic_local;
    } else if (buffer == "@All" || buffer == "*") {
      m_is_local = NULL;
    } else if (buffer[0] == '@' && buffer != "@Panes:" ) {
      if ( buffer != "@Proc:") { // Skip empty @Panes section.
        std::cerr << "Rocstar: Error (read_by_control_file): in control file "
                  << control_file_name << ": expected pane info but found '"
                  << buffer << '\'' << std::endl;
        return;
      } else {
        // ignore file when no pane
        local_patterns.clear();
      }
    } else {
      m_is_local = &Rocin::explicit_local;

      for (;;) {
        if (buffer[0] < '0' || buffer[0] > '9')
          break;

        std::istringstream sin(buffer);
        int p=-1;

        sin >> p;

        if ( p>=0) {
          m_pane_ids.insert(p);
        } else {
          break;
        }

        if (!fin.eof()) fin >> buffer;
        else break;
      }
    }

    patterns.insert(patterns.end(), local_patterns.begin(),
                    local_patterns.end());

    if ( myRank >= 0 && rank == myRank) break;
  }
  fin.close();

  if ( myRank>=0 && rank != myRank) 
    std::cerr << "Rocstar: Warning: Did not find matching control blocks for process " << myRank << std::endl;

  // If myRank is a wild card and panes aren't specified explicitly, then
  // make all panes local.
  if (myRank < 0 && m_is_local != &Rocin::explicit_local)
    m_is_local = NULL;

  // Convert the patterns into a single string
  std::string files;
  std::vector<std::string>::const_iterator p;
  for (p=patterns.begin(); p!=patterns.end(); ++p)
    files = files+" "+(*p).c_str();

  // Invoke read_window
  //std::cout << __FILE__ << __LINE__ << " files = " << files << std::endl;
  read_window(files.c_str(), window_name, myComm, NULL, time_level, str_len);

  m_pane_ids.clear();
  m_offset = 0;
  m_base = 0;
  m_is_local = NULL;
}

//! Read in metadata from files, and optionally read in array data as well
/*!
 *  \param filename_patterns Specifies patterns (reg. expressions) of files 
 *         to be read in.
 *  \param window_name Specifies a name of the window to be created.
 *  \param comm The MPI communicator to use. If is NULL, the default 
 *         communicator of Roccom will be used.
 *  \param is_local A function pointer, which determines wheter a pane 
 *         should be read by a process.
 *  \param time_level the time stamp of the dataset to be read.  If time
 *         level is NULL (default) or "", then the first time_level in the
 *         file will be used. 
 *  \param str_len if present and positive, the time stamp of the dataset 
 *         will be copied to the time_level string. 
 */
void Rocin::read_window( const char* filename_patterns,
                         const char* window_name,
                         const MPI_Comm* comm,
                         RulesPtr is_local,
                         char* time_level,
                         const int* str_len)
{
  read_windows( filename_patterns, window_name, "", comm, is_local,
                time_level, str_len);
}

std::string 
CWD()
{
  char buf[1024];
  return(std::string(getcwd(buf,1024)));
}
//! Read in metadata from files, and optionally read in array data as well
//! Read in metadata from files, and optionally read in array data as well
/*!
 *  \param filename_patterns Specifies patterns (reg. expressions) of files 
 *         to be read in.
 *  \param window_prefix Specifies a prefix of the name(s) of the window(s) 
 *         to be created.
 *  \param material_names Spcifies a space-delimited list of materials 
 *         to be read from the files.
 *  \param comm The MPI communicator to use. If is NULL, the default 
 *         communicator of Roccom will be used.
 *  \param is_local A function pointer, which determines wheter a pane 
 *         should be read by a process.
 *  \param time_level the time stamp of the dataset to be read.  If time
 *         level is NULL (default) or "", then the first time_level in the
 *         file will be used. 
 *  \param str_len if present and positive, the time stamp of the dataset 
 *         will be copied to the time_level string. 
 */
void Rocin::read_windows(const char* filename_patterns,
                         const char* window_prefix,
                         const char* material_names,
                         const MPI_Comm* comm,
                         RulesPtr is_local,
                         char* time_level,
                         const int* str_len)
{
  //std::cout << __FILE__ << __LINE__ 
  //          << " filename_patterns = " << filename_patterns
  //          << std::endl;
  const MPI_Comm default_comm=COM_get_default_communicator();
  const MPI_Comm comm_null=MPI_COMM_NULL;
  const MPI_Comm* myComm = 
    COMMPI_Initialized()? (comm ? comm : &default_comm ) : &comm_null;

  int rank, nprocs;
  if ( *myComm != MPI_COMM_NULL) {
    MPI_Comm_rank(*myComm, &rank);
    MPI_Comm_size(*myComm, &nprocs);
  }
  else {
    rank = 0; nprocs=1;
  }

  DEBUG_MSG("time_level == " << (time_level == NULL ? "<Null>" : time_level));
  std::string time( time_level?time_level:"");
  DEBUG_MSG("time == " << time);

  std::set<std::string> materials;
  char* buffer;
  char* token;

  // Insert list of material names into materials.
  if (material_names != NULL) {
    buffer = new char[strlen(material_names)+1];
    strcpy(buffer, material_names);

    token = strtok(buffer, " \t\n");
    while (token != NULL) {
      materials.insert(token);
      token = strtok(NULL, " \t\n");
    }
    delete[] buffer;
  }

  buffer = new char[strlen(filename_patterns)+1];
  strcpy(buffer, filename_patterns);

  BlockMM_HDF4 blocks_HDF4;
#ifdef USE_CGNS
  BlockMM_CGNS blocks_CGNS;
#endif // USE_CGNS

  token = strtok(buffer, " \t\n");
  if (token != NULL) {
#ifndef _NO_GLOB_
    glob_t globbuf;
    
    glob(token, 0, cast_err_func(glob,glob_error), &globbuf);
    token = strtok(NULL, " \t\n");
    while (token != NULL) {
      glob(token, GLOB_APPEND, cast_err_func(glob,glob_error), &globbuf);
      token = strtok(NULL, " \t\n");
    }

    if ( globbuf.gl_pathc==0 && buffer[0]!='\0') 
      std::cerr << "Rocstar: Warning: Found no matching files for pattern "
                << buffer << std::endl;

    // Extracts metadata from a list of files.
    // Opens each file, scans dataset, identifies windows, panes, and attribute
    // Puts this information into blocks.
    // MS
    //std::cout << __FILE__ << __LINE__ 
    //          << "gl_pathv = " << globbuf.gl_pathv[0]
    //          << std::endl;  
    // search for hdf or cgns
    bool isHDF = false; 
    bool isCGNS = false;
    if ((std::string(globbuf.gl_pathv[0])).find(".hdf") != std::string::npos){
       isHDF = true;
       std::cout << "npos = " << (std::string(globbuf.gl_pathv[0])).find(".hdf") <<std::endl; 
    }
    isCGNS = isHDF ? false:true;
    //std::cout << "isHdf = " << isHDF << " isCGNS = " << isCGNS << std::endl;
    // MS End
    if (isHDF){
       scan_files_HDF4(globbuf.gl_pathc, globbuf.gl_pathv, blocks_HDF4, time,
                       m_HDF2COM);
       //std::cout << "Read HDF!" << std::endl;
    } else {
#ifdef USE_CGNS 
    scan_files_CGNS(globbuf.gl_pathc, globbuf.gl_pathv, blocks_CGNS, time,
                    m_CGNS2COM);
    //std::cout << "Read CGNS!" << std::endl;
#endif
    }
// Original
//#ifndef USE_CGNS
//    scan_files_HDF4(globbuf.gl_pathc, globbuf.gl_pathv, blocks_HDF4, time,
//                    m_HDF2COM);
//#else 
//    scan_files_CGNS(globbuf.gl_pathc, globbuf.gl_pathv, blocks_CGNS, time,
//                    m_CGNS2COM);
//#endif
// Original End
    globfree(&globbuf);
#else // No glob function on this system
    // Create a char** of n filenames matching patterns stored in buffer
    //  each token is a new pattern
    std::list<std::string> matching_filenames;
    while (token != NULL) {
      std::string dirname(CWD());
      std::string tstring(token);
      std::string::size_type x = tstring.find_last_of("/");
      if(x != std::string::npos){
	dirname += ("/" + tstring.substr(0,x));
	tstring.erase(0,x+1);
      }
      Directory directory(dirname);
      if(directory){
	Directory::iterator di = directory.begin();
	while(di != directory.end()){
	  if(!fnmatch(tstring.c_str(),di->c_str(),0))
	    matching_filenames.push_back(dirname + "/" + *di);
	  di++;
	}
      }
      token = strtok(NULL, " \t\n");
    }
    if(matching_filenames.empty() && buffer[0] != '\0')
      std::cerr << "Rocstar: Warning: Found no matching files for pattern "
		<< buffer << std::endl;
    unsigned int nmatch = matching_filenames.size();
    std::vector<char *> matches(nmatch);
    std::list<std::string>::iterator li = matching_filenames.begin();
    unsigned int ccount = 0;
    while(li != matching_filenames.end()){
      unsigned int lis = li->size();
      matches[ccount++] = new char [lis+1];
      strcpy(matches[ccount-1],li->c_str());
      matches[ccount-1][lis] = '\0';
      li++;
    }
#ifndef USE_CGNS
    scan_files_HDF4(nmatch, &matches[0], blocks_HDF4, time,
                    m_HDF2COM);
#else 
    scan_files_CGNS(nmatch, &matches[0], blocks_CGNS, time,
                    m_CGNS2COM);
#endif
    ccount = 0;
    while(ccount < nmatch){
      if(matches[ccount])
	delete [] matches[ccount];
      ccount++;
    }
#endif
  }

  delete[] buffer;

  // Copy out time level
  if ( time_level && str_len && *str_len) { 
    // TODO: Run MPI_Allgather to send time level to those with no data
    // and check whether all processes obtained the same non-empty time level.
    std::strncpy( time_level, time.c_str(), *str_len-1); 
    time_level[*str_len-1] = '\0';
  }

  std::string name;
  std::set<std::string>::iterator p = materials.begin();
  std::pair<BlockMM_HDF4::iterator, BlockMM_HDF4::iterator> range_HDF4;
#ifdef USE_CGNS
  std::pair<BlockMM_CGNS::iterator, BlockMM_CGNS::iterator> range_CGNS;
#endif // USE_CGNS

  if (materials.empty()) {
    // Default value of material_names is NULL
    // we assume the file only contains one type of material, and 
    // the window name is the window_prefix.
    range_HDF4.first = blocks_HDF4.begin();
    range_HDF4.second = blocks_HDF4.end();
#ifdef USE_CGNS
    range_CGNS.first = blocks_CGNS.begin();
    range_CGNS.second = blocks_CGNS.end();
#endif // USE_CGNS
    name = window_prefix;

    COM_new_window(name.c_str(), *myComm);

    new_attributes( range_HDF4.first, range_HDF4.second,
#ifdef USE_CGNS
                    range_CGNS.first, range_CGNS.second,
#endif // USE_CGNS
                    name, myComm, rank, nprocs);
    register_panes( range_HDF4.first, range_HDF4.second,
#ifdef USE_CGNS
                    range_CGNS.first, range_CGNS.second,
#endif // USE_CGNS
                    name, is_local, myComm, rank, nprocs);

    load_data_HDF4(range_HDF4.first, range_HDF4.second,
                   name, myComm, rank, nprocs);
#ifdef USE_CGNS
    load_data_CGNS(range_CGNS.first, range_CGNS.second,
                   name, myComm, rank, nprocs);
#endif // USE_CGNS

    broadcast_win_attributes( range_HDF4.first == range_HDF4.second
#ifdef USE_CGNS
                              && range_CGNS.first == range_CGNS.second
#endif // USE_CGNS
                              , name, myComm, rank, nprocs);

    COM_window_init_done(name.c_str());
  } else {

    // Iterate through each material (one window per material).
    while (p != materials.end()) {
      if (blocks_HDF4.count(*p) == 0
#ifdef USE_CGNS
          && blocks_CGNS.count(*p) == 0
#endif // USE_CGNS
          ) {
        std::cerr << "Rocstar: Warning (read_windows): could not find '" << *p << "'."
                  << std::endl;
        ++p; // Increment p.
        continue;
      }

      range_HDF4 = blocks_HDF4.equal_range(*p);
#ifdef USE_CGNS
      range_CGNS = blocks_CGNS.equal_range(*p);
#endif // USE_CGNS
      name = window_prefix + *p;

      COM_new_window(name.c_str(), *myComm);

      new_attributes( range_HDF4.first, range_HDF4.second,
#ifdef USE_CGNS
                      range_CGNS.first, range_CGNS.second,
#endif // USE_CGNS
                      name, myComm, rank, nprocs);
      register_panes( range_HDF4.first, range_HDF4.second,
#ifdef USE_CGNS
                      range_CGNS.first, range_CGNS.second,
#endif // USE_CGNS
                      name, is_local, myComm, rank, nprocs);

      load_data_HDF4(range_HDF4.first, range_HDF4.second,
                     name, myComm, rank, nprocs);
#ifdef USE_CGNS
      load_data_CGNS(range_CGNS.first, range_CGNS.second,
                     name, myComm, rank, nprocs);
#endif // USE_CGNS

      broadcast_win_attributes( range_HDF4.first == range_HDF4.second
#ifdef USE_CGNS
                                && range_CGNS.first == range_CGNS.second
#endif // USE_CGNS
                                , name, myComm, rank, nprocs);

      COM_window_init_done(name.c_str());

      ++p;
    }
  }

  // Free memory.
  free_blocks(blocks_HDF4);
#ifdef USE_CGNS
  free_blocks(blocks_CGNS);
#endif // USE_CGNS
}

extern "C" void Rocin_load_module( const char *name)
{ Rocin::init( std::string(name)); }

extern "C" void Rocin_unload_module( const char *name) 
{ Rocin::finalize( std::string(name)); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// All possible Fortran bindings

extern "C" void COM_F_FUNC2(rocin_load_module, ROCIN_LOAD_MODULE)( const char *name, long int length)
{ Rocin_load_module( std::string(name, length).c_str()); }
extern "C" void COM_F_FUNC2(rocin_unload_module, ROCIN_UNLOAD_MODULE)( const char *name, long int length) 
{ Rocin_unload_module( std::string(name, length).c_str()); }

#endif // DOXYGEN_SHOULD_SKIP_THIS







