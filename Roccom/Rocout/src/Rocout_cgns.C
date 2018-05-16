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
/** \file Rocout_cgns.C
 *  Implementation of Rocout CGNS routines.
 */
/*  Author John Norris
 *  Initial date:   Oct. 19, 2004
 */

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "cgnslib.h"
#include "Rocout.h"
#include "Element_accessors.h"
#include "Rocout_cgns.h"
#include <unistd.h>

// MS
char cwd[FILENAME_MAX];
#define GetCurrentDir getcwd
// MS End
#ifndef DOXYGEN_SHOULD_SKIP_THIS
USE_COM_NAME_SPACE
#endif

// #define DEBUG_DUMP_PREFIX "."
//#define DEBUG_MSG(x) std::cout << "ROCOUT DEBUG: " << __LINE__ << ": " << x << std::endl
#define DEBUG_MSG(x)

/**
 ** Perform CGNS error-checking.
 **
 ** Check the return value of the CGNS function and print an error message
 ** if it is non-zero.
 **
 ** \param routine The CGNS function to call. (Input)
 ** \param args The argument list, in parentheses. (Input)
 **/
#define CG_CHECK(routine, args) \
{ \
  int ier = routine args; \
  if (ier != 0 && errorhandle != "ignore") { \
    std::cerr << "Rocout::write_attribute: " #routine " (line " \
              << __LINE__ << " in " << __FILE__ << ") failed: " \
              << cg_get_error() << std::endl; \
    if (errorhandle == "abort") { \
      if (COMMPI_Initialized()) \
        MPI_Abort(MPI_COMM_WORLD, 0); \
      else \
        abort(); \
    } \
  } \
}

/**
 ** Call a templated function based on a Roccom data type.
 **
 ** A switch statement performs a typedef based on the data type provided,
 ** then calls the templated function.
 **
 ** \param dType The Roccom data type. (Input)
 ** \param funcCall The templated function, with arguments. (Input)
 **/
#define SwitchOnCOMDataType(dType, funcCall) \
  switch (dType) { \
    case COM_CHAR: \
    case COM_CHARACTER: { \
      typedef char COM_TT; \
      funcCall; \
      break; } \
    case COM_INT: \
    case COM_INTEGER: { \
      typedef int COM_TT; \
      funcCall; \
      break; } \
    case COM_FLOAT: \
    case COM_REAL: { \
      typedef float COM_TT; \
      funcCall; \
      break; } \
    case COM_DOUBLE: \
    case COM_DOUBLE_PRECISION: { \
      typedef double COM_TT; \
      funcCall; \
      break; } \
  }

/**
 ** Automatically close a CGNS file.
 **
 ** This class closes a CGNS file automatically when it goes out of scope.
 **/
class AutoCloser
{
  public:
    inline AutoCloser(int fn, const std::string& eh)
    : m_fn(fn), errorhandle(eh) {}
    inline ~AutoCloser() { CG_CHECK(cg_close, (m_fn)); }

  private:
    int m_fn;
    const std::string& errorhandle;
};

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

/**
 ** Convert a Roccom data type to a CGNS data type.
 **
 ** \param dType The Roccom data type. (Input)
 ** \returns The equivalent CGNS data type.
 **/
static DataType_t COM2CGNS(int dType)
{
  switch (dType) {
    case COM_CHAR:
    case COM_CHARACTER:
      return Character;

    case COM_INT:
    case COM_INTEGER:
      return Integer;

    case COM_FLOAT:
    case COM_REAL:
      return RealSingle;

    case COM_DOUBLE:
    case COM_DOUBLE_PRECISION:
      return RealDouble;

    default:
      std::cerr << "ROCOUT ERROR: Cannot convert COM type " << dType
                << " to a CGNS type." << std::endl;
  }

  return DataTypeNull;
}

/**
 ** Build and write out an Elements_t node.
 **
 ** Convert a Roccom connectivity table to a CGNS Elements_t node, and write
 ** it to file.
 **
 ** \param fn The CGNS file number. (Input)
 ** \param B The CGNS base. (Input)
 ** \param Z The CGNS zone. (Input)
 ** \param eType The Roccom element Type. (Input)
 ** \param nElems The number of elements to write out. (Input)
 ** \param elems The Roccom element data. (Input)
 ** \param isStaggered True is the table is staggered. (Input)
 ** \param length The table capacity. (Input)
 ** \param offset The index of the first element in this table. (Input)
 ** \param label The name for the Elements_t node. (Input)
 ** \returns The index of the newly-written Elements_t node.
 **/
static int WriteElements(int fn, int B, int Z, int eType, int nElems,
                         const int* elems, bool isStaggered, int length,
                         int offset, const std::string& label,
                         const std::string& errorhandle)
{
  ElementType_t elemType;
  int nNodes;

  // Determine the number of nodes per element and the CGNS element type.
  switch (eType) {
    case Connectivity::BAR2:
      elemType = BAR_2;
      nNodes = 2;
      break;

    case Connectivity::BAR3:
      elemType = BAR_3;
      nNodes = 3;
      break;

    case Connectivity::TRI3:
      elemType = TRI_3;
      nNodes = 3;
      break;

    case Connectivity::TRI6:
      elemType = TRI_6;
      nNodes = 6;
      break;

    case Connectivity::QUAD4:
      elemType = QUAD_4;
      nNodes = 4;
      break;

    case Connectivity::QUAD8:
      elemType = QUAD_8;
      nNodes = 8;
      break;

    case Connectivity::QUAD9:
      elemType = QUAD_9;
      nNodes = 9;
      break;

    case Connectivity::TET4:
      elemType = TETRA_4;
      nNodes = 4;
      break;

    case Connectivity::TET10:
      elemType = TETRA_10;
      nNodes = 10;
      break;

    case Connectivity::PYRIMID5:
      elemType = PYRA_5;
      nNodes = 5;
      break;

    case Connectivity::PYRIMID14:
      elemType = PYRA_14;
      nNodes = 14;
      break;

    case Connectivity::PRISM6:
      elemType = PENTA_6;
      nNodes = 6;
      break;

    case Connectivity::PRISM15:
      elemType = PENTA_15;
      nNodes = 15;
      break;

    case Connectivity::PRISM18:
      elemType = PENTA_18;
      nNodes = 18;
      break;

    case Connectivity::HEX8:
      elemType = HEXA_8;
      nNodes = 8;
      break;

    case Connectivity::HEX20:
      elemType = HEXA_20;
      nNodes = 20;
      break;

    case Connectivity::HEX27:
      elemType = HEXA_27;
      nNodes = 27;
      break;
  };

  // Allocate memory for the new connectivity table.
  std::vector<int> buffer;
  const int* pData;

  // Loop through the elements and get the indices of their nodes.
  if (elems != NULL) {
    buffer.resize(std::max(nElems * nNodes, 1));
    int i, j, k = 0;
    if (isStaggered) {
      for (i=0; i<nElems; ++i)
        for (j=0; j<nNodes; ++j,++k)
          buffer[k] = elems[j*length+i];
    pData = &buffer[0];
    } else {
      pData = elems;
    }
  } else {
    buffer.resize(1);
    pData = &buffer[0];
    ++nElems;
  }

  // Write out the Elements_t node.
  int S = 0;
  CG_CHECK(cg_section_write, (fn, B, Z, label.c_str(), elemType, offset,
                              offset + nElems - 1, 0, &(pData[0]), &S));

  return S;
}

/**
 ** Find the minimum and maximum core value.
 **
 ** Scan the given scalar data to find the minimum and maximum values,
 ** ignoring ghost data.  Use these values to build a string of the form
 ** "<minimum_value>, <maximum_value>" to be used in a Descriptor_t node.
 **
 ** \param rank The dimensionality of the data. (Input)
 ** \param size The size of the data in each dimension. (Input)
 ** \param rind The number of ghost values at the beginning and end of
 **             each dimension. (Input)
 ** \param pData The data to scan. (Input)
 ** \param range The string with the min and max values. (Output)
 **/
template <typename T>
void FindRange(int rank, const int* size, const int* rind, const T* pData,
               std::string& range)
{
  DEBUG_MSG("Entering FindRange<T>");
  int i, j, k;
  const T* t = pData;
  T min = std::numeric_limits<T>::max();
  T max = -min;
  for (k=0; k<(rank>2?size[2]-rind[5]:1); ++k) {
    for (j=0; j<(rank>1?size[1]:1); ++j) {
      for (i=0; i<size[0]; ++i,++t) {
        if (i >= std::min(rind[0], size[0] - 1)
            && i < std::max(1, size[0] - rind[1])
            && j >= std::min(rind[2], size[1] - 1)
            && j < std::max(1, size[1] - rind[3])
            && k >= std::min(rind[4], size[2] - 1)) {
          if (min > *t)
            min = *t;
          if (max < *t)
            max = *t;
        }
      }
    }
  }

  std::ostringstream sout;
  sout << min << ", " << max;
  range = sout.str();
}

void FindRange(int rank, const int* size, const int* rind, const char* pData,
               std::string& range)
{
  DEBUG_MSG("Entering FindRange<char>");
  int i, j, k;
  const char* t = pData;
  char min = std::numeric_limits<char>::max();
  char max = -min;
  for (k=0; k<(rank>2?size[2]-rind[5]:1); ++k) {
    for (j=0; j<(rank>1?size[1]:1); ++j) {
      for (i=0; i<size[0]; ++i,++t) {
        if (i >= std::min(rind[0], size[0] - 1)
            && i < std::max(1, size[0] - rind[1])
            && j >= std::min(rind[2], size[1] - 1)
            && j < std::max(1, size[1] - rind[3])
            && k >= std::min(rind[4], size[2] - 1)) {
          if (min > *t)
            min = *t;
          if (max < *t)
            max = *t;
        }
      }
    }
  }

  std::ostringstream sout;
  sout << (int)min << ", " << (int)max;
  range = sout.str();
}

/**
 ** Find the minimum and maximum core magnitude.
 **
 ** Scan the given vector data to find the minimum and maximum magnitudes,
 ** ignoring ghost data.  Use these values to build a string of the form
 ** "<minimum_value>, <maximum_value>" to be used in a Descriptor_t node.
 **
 ** \param physDim The dimensionality of the vectors. (Input)
 ** \param rank The dimensionality of the data. (Input)
 ** \param size The size of the data in each dimension. (Input)
 ** \param rind The number of ghost values at the beginning and end of
 **             each dimension. (Input)
 ** \param pData The data to scan. (Input)
 ** \param range The string with the min and max values. (Output)
 **/
template <typename T>
void FindMagnitudeRange(int physDim, int rank, const int* size,
                        const int* rind, const T** pData, std::string& range)
{
  int i, j, k;
  std::vector<T> buf;
  const T* t[3] = { pData[0], pData[1], pData[2] };
  if (physDim == 2) {
    int s = size[0];
    if (rank == 2)
      s *= size[1];
    buf.resize(s, (T)0);
    pData[2] = &(buf[0]);
  }

  T min = std::numeric_limits<T>::max();
  T max = (T)0;
  for (k=0; k<(rank>2?size[2]-rind[5]:1); ++k) {
    for (j=0; j<(rank>1?size[1]:1); ++j) {
      for (i=0; i<size[0]; ++i,++t[0],++t[1],++t[2]) {
        if (i >= std::min(rind[0], size[0] - 1)
            && i < std::max(1, size[0] - rind[1])
            && j >= std::min(rind[2], size[1] - 1)
            && j < std::max(1, size[1] - rind[3])
            && k >= std::min(rind[4], size[2] - 1)) {
          T val = *t[0] * *t[0] + *t[1] * *t[1] + *t[2] * *t[2];
          if (val < min)
            min = val;
          if (val > max)
            max = val;
        }
      }
    }
  }

  std::ostringstream sout;
  sout << std::sqrt((double)min) << ", " << std::sqrt((double)max);
  range = sout.str();
}

/**
 ** Find the minimum and maximum core trace.
 **
 ** Scan the given tensor data to find the minimum and maximum traces,
 ** ignoring ghost data.  Use these values to build a string of the form
 ** "<minimum_value>, <maximum_value>" to be used in a Descriptor_t node.
 **
 ** \param physDim The dimensionality of the tensors. (Input)
 ** \param rank The dimensionality of the data. (Input)
 ** \param size The size of the data in each dimension. (Input)
 ** \param rind The number of ghost values at the beginning and end of
 **             each dimension. (Input)
 ** \param pData The data to scan. (Input)
 ** \param range The string with the min and max values. (Output)
 **/
template <typename T>
void FindTraceRange(int physDim, int rank, const int* size, const int* rind,
                    const T** pData, std::string& range)
{
  int i, j, k;
  std::vector<T> buf;
  const T* t[3];
  t[0] = pData[0];
  if (physDim == 2) {
    t[1] = pData[3];
    int s = size[0];
    if (rank == 2)
      s *= size[1];
    buf.resize(s, (T)0);
    t[2] = &(buf[0]);
  } else {
    t[1] = pData[4];
    t[2] = pData[8];
  }

  T min = std::numeric_limits<T>::max();
  T max = -min;
  for (k=0; k<(rank>2?size[2]-rind[5]:1); ++k) {
    for (j=0; j<(rank>1?size[1]:1); ++j) {
      for (i=0; i<size[0]; ++i,++t[0],++t[1],++t[2]) {
        if (i >= std::min(rind[0], size[0] - 1)
            && i < std::max(1, size[0] - rind[1])
            && j >= std::min(rind[2], size[1] - 1)
            && j < std::max(1, size[1] - rind[3])
            && k >= std::min(rind[4], size[2] - 1)) {
          T val = *t[0] + *t[1] + *t[2];
          if (val < min)
            min = val;
          if (val > max)
            max = val;
        }
      }
    }
  }

  std::ostringstream sout;
  sout << min << ", " << max;
  range = sout.str();
}

//@{
/**
 ** Determine and write the exponents for the basic units of measurement.
 **
 ** These rountines converts a string (e.g. "m" or "(kg K)/(m^2 s)") into five
 ** exponents for mass (kg), length (m), temperature (K), time (s), and
 ** angle (degrees).
 **
 ** \note Currently recognises kg, m, s, K, J, N, and Pa.
 **/
float ExtractExponent(std::string& unit)
{
  if (unit[0] != '^')
    return 1.f;

  std::istringstream in(&unit[1]);
  float e;
  in >> e;

  if (in.eof())
    unit.erase();
  else
    in >> unit;

  return e;
}

void ModifyExponents(const std::string& unit, float* exponents, float dir)
{
  std::string u = unit;
  while (!u.empty()) {
    float exp;
    if (u.substr(0, 2) == "kg") { // Kilograms
      u.erase(0, 2);
      exp = ExtractExponent(u);
      exponents[0] += exp * dir;
    } else if (u[0] == 'm') { // Meters
      u.erase(0, 1);
      exp = ExtractExponent(u);
      exponents[1] += exp * dir;
    } else if (u[0] == 's') { // Seconds
      u.erase(0, 1);
      exp = ExtractExponent(u);
      exponents[2] += exp * dir;
    } else if (u[0] == 'K') { // Kelvin
      u.erase(0, 1);
      exp = ExtractExponent(u);
      exponents[3] += exp * dir;
    } else if (u[0] == 'J') { // Joules (kg m^2 / s^2)
      u.erase(0, 1);
      exp = ExtractExponent(u);
      exponents[0] += exp * dir;
      exponents[1] += 2 * exp * dir;
      exponents[2] -= 2 * exp * dir;
    } else if (u[0] == 'N') { // Newtons (kg m / s^2)
      u.erase(0, 1);
      exp = ExtractExponent(u);
      exponents[0] += exp * dir;
      exponents[1] += exp * dir;
      exponents[2] -= 2 * exp * dir;
    } else if (u.substr(0, 2) == "Pa") { // Pascals (kg / m s^2)
      u.erase(0, 2);
      exp = ExtractExponent(u);
      exponents[0] += exp * dir;
      exponents[1] -= exp * dir;
      exponents[2] -= 2 * exp * dir;
    } else {
      std::cerr << "ROCOUT ERROR: Unable to parse unit '" << unit
                << "': stuck at '" << u << "'." << std::endl;
      break;
    }
  }
}

template <typename T>
static int cg_array_core_write_internal(char const* name, DataType_t dType,
                                        int rank, int* rind, int const* size,
                                        T const* data)
{
  int i, j, k, x = 0, total = 1;
  int lMin[3] = { 0, 0, 0 }, lMax[3] = { 1, 1, 1 };
  int lSize[3] = { 1, 1, 1 }, newSize[3] = { 1, 1, 1 };
  for (i=0; i<rank; ++i) {
    lMin[i] = rind[2*i];
    lMax[i] = size[i] - rind[2*i+1];
    lSize[i] = size[i];
    newSize[i] = lMax[i] - lMin[i];
    total *= newSize[i];
  }
  std::vector<T> buffer(total);

  // Extract subset.
  for (k=lMin[2]; k<lMax[2]; ++k)
    for (j=lMin[1]; j<lMax[1]; ++j)
      for (i=lMin[0]; i<lMax[0]; ++i,++x)
        buffer[x] = data[i+j*lSize[0]+k*lSize[0]*lSize[1]];

  return cg_array_write(name, dType, rank, newSize, &buffer[0]);
}

static int cg_array_core_write(char const* name, DataType_t dType,
                               int rank, int* rind, int const* size,
                               void const* data)
{
  // Check for a full write.
  bool full = true;
  int i;
  for (i=0; i<rank && full; ++i)
    full = (rind[2*i] == 0 && rind[2*i+1] == 0);
  if (full)
    return cg_array_write(name, dType, rank, size, data);

  switch (dType) {
    case Character:
       return cg_array_core_write_internal(name, dType, rank, rind, size,
                                           static_cast<char const*>(data));
    case Integer:
       return cg_array_core_write_internal(name, dType, rank, rind, size,
                                           static_cast<int const*>(data));
    case RealSingle:
       return cg_array_core_write_internal(name, dType, rank, rind, size,
                                           static_cast<float const*>(data));
    case RealDouble:
       return cg_array_core_write_internal(name, dType, rank, rind, size,
                                           static_cast<double const*>(data));
  }

  return CG_ERROR;
}


static void cg_exponents_as_string_write(std::string unit,
                                         const std::string& errorhandle)
{
  // Eliminate spaces and parentheses.
  std::string::size_type i = 0;
  while (i < unit.length()) {
    if (unit[i] == ' ' || unit[i] == '(' || unit[i] == ')')
      unit.erase(i, 1);
    else
      ++i;
  }

  // Normalize exponent operator.
  i = unit.find("**");
  while (i != std::string::npos) {
    unit.replace(i, 2, "^");
    i = unit.find("**", i);
  }

  // Find the numerator and denominator.
  std::string top, bottom;
  i = unit.find('/');
  if (i == std::string::npos) {
    top = unit;
  } else {
    top = unit.substr(0, i);
    bottom = unit.substr(i + 1);
  }

  float exponents[5] = { 0.f, 0.f, 0.f, 0.f, 0.f };
  if (!top.empty() && top != "1")
    ModifyExponents(top, exponents, 1.f);
  if (!bottom.empty() && bottom != "1")
    ModifyExponents(bottom, exponents, -1.f);

  CG_CHECK(cg_exponents_write, (RealSingle, exponents));
}
//@}

/**
 ** Find the named Base_t node, or create one if it doesn't exist.
 **
 ** Search the existing Base_t nodes for one with the given name, physical
 ** dimensions and cell dimensions.  If no node of that name exists, then
 ** create one.
 **
 ** \param fn The CGNS file number. (Input)
 ** \param name The name of the Base_t node. (Input)
 ** \param cellDim The required cell dimension. (Input)
 ** \param physDim The required physical dimension. (Input)
 ** \param B The index of the Base_t node. (Output)
 ** \returns 0 on success, 1 otherwise.
 **/
static int cg_base_find_or_create(int fn, const char* name, int cellDim,
                                  int physDim, int* B,
                                  const std::string& errorhandle)
{
  char baseName[33];
  int nBases, cDim, pDim;
  CG_CHECK(cg_nbases, (fn, &nBases));

  // Check each existing Base_t node for a name match.
  for (*B=1; *B<=nBases; ++(*B)) {
    CG_CHECK(cg_base_read, (fn, *B, baseName, &cDim, &pDim));
    if (std::strncmp(name, baseName, 32) == 0) {
      // Confirm that the dimensions are correct.
      if (cDim == cellDim && pDim == physDim)
        return 0;
      DEBUG_MSG("ERROR: Dimensions don't match (cellDim == " << cDim
                << ", physDim == " << pDim << ")");
      return 1;
    }
  }

  CG_CHECK(cg_base_write, (fn, name, cellDim, physDim, B));
  return 0;
}

/**
 ** Find the named Zone_t node, or create one if it doesn't exist.
 **
 ** Search the existing Zone_t nodes under the given Base_t node for one
 ** with the given name, size and zone type.  If no node of that name
 ** exists, then create one.
 **
 ** \param fn The CGNS file number. (Input)
 ** \param B The CGNS Base_t node index. (Input)
 ** \param name The name of the Zone_t node. (Input)
 ** \param sizes The size of the zone. (Input)
 ** \param zType The type of zone (Structured or Unstructured). (Input)
 ** \param Z The index of the Zone_t node. (Output)
 ** \returns 0 on success, 1 otherwise.
 **/
static int cg_zone_find_or_create(int fn, int B, const char* name,
                                  const int* sizes, ZoneType_t zType, int* Z,
                                  const std::string& errorhandle)
{
  char zoneName[33];
  int nZones;
  std::vector<int> sz1(&(sizes[0]), &(sizes[9]));
  std::vector<int> sz2(9);
  CG_CHECK(cg_nzones, (fn, B, &nZones));
  for (*Z=1; *Z<=nZones; ++(*Z)) {
    CG_CHECK(cg_zone_read, (fn, B, *Z, zoneName, &(sz2[0])));
    if (std::strncmp(name, zoneName, 32) == 0) {
      if (sz1 != sz2) {
        // MS
        //DEBUG_MSG("ERROR: Zone dimensions don't match");
        //return 1;
        // MS End
      }

      ZoneType_t zt;
      CG_CHECK(cg_zone_type, (fn, B, *Z, &zt));
      if (zType == zt)
        return 0;

      DEBUG_MSG("ERROR: Zone type doesn't match");
      return 1;
    }
  }

  if (zType == Unstructured && sizes[2] > sizes[0])
    std::cerr << "ROCOUT ERROR: Creating unstructured zone \"" << name
              << "\" with invalid dimensions { " << sizes[0] << ", " << sizes[1]
              << ", " << sizes[2] << " }" << std::endl;

  CG_CHECK(cg_zone_write, (fn, B, name, sizes, zType, Z));
  return 0;
}

/**
 ** Find the named IntegralData_t node, or create one if it doesn't exist.
 **
 ** Search the IntegralData_t nodes under the current location for one with
 ** the given name.  If no node of that name exists, then create one.
 **
 ** \param name The name of the IntegralData_t node. (Input)
 ** \param I The index of the IntegralData_t node. (Output)
 ** \returns 0 on success, 1 otherwise.
 **/
static int cg_integral_find_or_create(const char* name, int* I,
                                      const std::string& errorhandle)
{
  char intName[33];
  int nIntegrals;
  CG_CHECK(cg_nintegrals, (&nIntegrals));
  for (*I=1; *I<=nIntegrals; ++(*I)) {
    CG_CHECK(cg_integral_read, (*I, intName));
    if (std::strncmp(name, intName, 32) == 0)
      return 0;
  }

  CG_CHECK(cg_integral_write, (name));
  return 0;
}

/**
 ** Find the named FlowSolution_t node, or create one if it doesn't exist.
 **
 ** Search the existing FlowSolution_t nodes under the given Zone_t node for
 ** one with the given name.  If no node of that name exists, then create one.
 **
 ** \param fn The CGNS file number. (Input)
 ** \param B The CGNS Base_t node index. (Input)
 ** \param Z The CGNS Zone_t node index. (Input)
 ** \param name The name of the FlowSolution_t node. (Input)
 ** \param location The location of the data. (Input)
 ** \param S The index of the FlowSolution_t node. (Output)
 ** \returns 0 on success, 1 otherwise.
 **/
static int cg_sol_find_or_create(int fn, int B, int Z, const char* name,
                                 GridLocation_t location, int* S,
                                 const std::string& errorhandle)
{
  char solName[33];
  GridLocation_t loc;
  int nSols;
  CG_CHECK(cg_nsols, (fn, B, Z, &nSols));
  for (*S=1; *S<=nSols; ++(*S)) {
    CG_CHECK(cg_sol_info, (fn, B, Z, *S, solName, &loc));
    if (std::strncmp(name, solName, 32) == 0) {
      return (loc == location ? 0 : 1);
    }
  }

  CG_CHECK(cg_sol_write, (fn, B, Z, name, location, S));
  return 0;
}

/**
 ** Return information on the named DataArray_t node.
 **
 ** Search the DataArray_t nodes under the current location for one with
 ** the given name, and return the index, data type, dimensionality and
 ** size.
 **
 ** \param name The name of the DataArray_t node. (Input)
 ** \param A The index of the DataArray_t node. (Output)
 ** \param dType The data type. (Output)
 ** \param rank The dimensionality of the data. (Output)
 ** \param size The size of each dimension of the data. (Output)
 ** \returns 0 on success, 1 otherwise.
 **/
static int cg_array_info_by_name(const char* name, int* A, DataType_t* dType,
                                 int* rank, int* size,
                                 const std::string& errorhandle)
{
  char arrayName[33];
  int nArrays;
  CG_CHECK(cg_narrays, (&nArrays));
  //std::cout << __FILE__<< __LINE__ << std::endl;
  //std::cout << " nArrays = " << nArrays << std::endl;
  for (*A=1; *A<=nArrays; ++(*A)) {
    CG_CHECK(cg_array_info, (*A, arrayName, dType, rank, size));
    //std::cout << " arrayName = " << arrayName << std::endl;
    if (std::strncmp(name, arrayName, 32) == 0)
      return 0;
  }
  return 1;
}

/**
 ** Write the data for the given attribute to file.
 **
 ** Write the given attribute to file using the CGNS format.  The attribute
 ** may be a "mesh", "all" or some other predefined attribute.
 **
 ** \param fname The name of the main datafile. (Input)
 ** \param mname The name of the optional mesh datafile. (Input)
 ** \param attr The attribute to write out. (Input)
 ** \param material The name of the material. (Input)
 ** \param timelevel The simulation time for this data. (Input)
 ** \param pane_id The id for the local pane. (Input)
 ** \param mode Write == 0, append == 1. (Input)
 **/
void write_attr_CGNS(const std::string& fname_in, const std::string& mfile,
                     const COM::Attribute* attr, const char* material,
                     const char* timelevel, int pane_id,
                     const std::string& ghosthandle,
                     const std::string& errorhandle, int mode)
{
  /*
  std::cout << " ------------------------------------------------------" << std::endl;
  std::cout << " Starting to write \n Data File = " << fname_in << std::endl;
  std::cout << " Mesh File = " << mfile << std::endl;
  std::cout << " timelevel = " << timelevel << std::endl;
  std::cout << " pane_id = " << pane_id << std::endl;
  std::cout << " mode = " << mode << std::endl;
  std::cout << " ------------------------------------------------------" << std::endl;
  */

  std::string fname(fname_in);
  bool writeGhost = (ghosthandle == "write");
  const Window* w = attr->window();
  COM_assertion(w != NULL);
  const Pane& pane = w->pane(pane_id);


  // Convert the timelevel from a string to a number.
  double timeValue;
  {
    std::istringstream sin(timelevel);
    sin >> timeValue;
  }

  // Covert the timeValue back to a string.  This "normalizes" the string.
  std::string timeLevel;
  {
    std::ostringstream sout;
    sout << timeValue;
    timeLevel = sout.str();
  }

  DEBUG_MSG("Writing to file '" << fname << "', mode == '"
            << (mode ? "append'" : "write'"));
  DEBUG_MSG("Using mesh file '" << mfile << "'");
  AutoCDer autoCD;
  std::string::size_type loc = fname.rfind('/');
  // MS
  GetCurrentDir(cwd, sizeof(cwd));
  //std::cout << __FILE__ << __LINE__ << "\n" << cwd << std::endl;
  std::string prefix = fname.substr(0, loc+1);
  // original
  /*
  if (loc != std::string::npos) {
    std::cout << __FILE__ << __LINE__ << " <<<<<<<<<<<<<<";
    chdir(fname.substr(0, loc).c_str());
    fname.erase(0, loc + 1);
  }
  */
  // original
  // MS

  // Open or create the file.
  int fn;
  if (mode == 0) {
    CG_CHECK(cg_open, (fname.c_str(), MODE_WRITE, &fn));
    CG_CHECK(cg_close, (fn));
  }
  CG_CHECK(cg_open, (fname.c_str(), MODE_MODIFY, &fn));

  // The file will be closed automagically when we exit this function.
  AutoCloser autoCloser(fn, errorhandle);

  // Find or create the base (corresponds to window/material).
  int i, B, nSteps = 0;
  std::vector<double> times;
  char buffer[33];
  std::string label;
  int cellDim = pane.dimension();
  int physDim = pane.attribute(COM::COM_NC)->size_of_components();
  //std::cout << __FILE__ << __LINE__ << std::endl;
  //std::cout << " cell dim = " << cellDim << std::endl;
  //std::cout << " phys dim = " << physDim << std::endl;
  // MS 
  if (cellDim == 0) 
      cellDim = 3;
  // MS End
  CG_CHECK(cg_base_find_or_create, (fn, material, cellDim, physDim, &B,
                                    errorhandle));
  //std::cout << __FILE__ << __LINE__ << std::endl;

  // Write the default Roccom units at the top.
  CG_CHECK(cg_goto, (fn, B, "end"));
  CG_CHECK(cg_units_write, (Kilogram, Meter, Second, Kelvin, Degree));

  // Read the time values in the BaseIterativeData_t node.
  // MS
  if (timeValue == 0){
     times.resize(1);
  } else {
     if (mode > 0 && cg_biter_read(fn, B, buffer, &nSteps) == 0) {
       CG_CHECK(cg_goto, (fn, B, "BaseIterativeData_t", 1, "end"));
       times.resize(nSteps + 1);
       CG_CHECK(cg_array_read_as, (1, RealDouble, &(times[0])));
     } else {
       times.resize(1);
     }
  }
  // MS End
  /* Original
  if (mode > 0 && cg_biter_read(fn, B, buffer, &nSteps) == 0) {
    CG_CHECK(cg_goto, (fn, B, "BaseIterativeData_t", 1, "end"));
    times.resize(nSteps + 1);
    CG_CHECK(cg_array_read_as, (1, RealDouble, &(times[0])));
  } else {
    times.resize(1);
  }
  */

  // Search for the current time.
  int timeIndex;
  for (timeIndex=0; timeIndex<nSteps; ++timeIndex)
    if (times[timeIndex] == timeValue) {
      //std::cout << " timeIndex = " << timeIndex << std::endl; 
      break;
    }

  // Update the time values in the BaseIterativeData_t node if necessary.
  if (timeIndex == nSteps) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    times[timeIndex] = timeValue;
    ++nSteps;

    // Write/rewrite the time values to the BaseIterativeData_t node.
    CG_CHECK(cg_biter_write, (fn, B, "TimeIterValues", nSteps));
    CG_CHECK(cg_goto, (fn, B, "BaseIterativeData_t", 1, "end"));
    CG_CHECK(cg_array_write, ("TimeValues", RealDouble, 1, &nSteps,
                              &(times[0])));
    --nSteps;
  }

  // Find or create the zone (corresponds to pane/section).  Use the pane
  // id as the zone name.
  std::string zoneName;
  {
    std::ostringstream sout;
    sout << std::setw(4) << std::setfill('0') << pane.id();
    zoneName = sout.str();
  }

  // Set up the sizes based on dimensionality and zone type.
  int Z, sizes[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  ZoneType_t zType;
  i = 0;
  if (pane.is_structured()) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    // Sizes should be core nodes/elements only.
    int ghost = pane.size_of_ghost_layers();
    sizes[i++] = pane.size_i() - 2 * ghost;
    sizes[i++] = pane.size_j() - 2 * ghost;
    if (physDim > 2 && cellDim > 2)
      sizes[i++] = pane.size_k() - 2 * ghost;
    sizes[i++] = sizes[0] - 1;
    sizes[i++] = sizes[1] - 1;
    if (physDim > 2 && cellDim > 2)
      sizes[i++] = sizes[2] - 1;
    zType = Structured;
  } else {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    sizes[0] = pane.size_of_real_nodes();
    sizes[1] = pane.size_of_real_elements();
    zType = Unstructured;
  }
  
  //std::cout << __FILE__ << __LINE__ << std::endl;
  // MS
  if (sizes[0] == 0){
     int meshfn;
     char zName[33];
     //std::vector<int> sz2(9);
     if (!mfile.empty()) {
      CG_CHECK(cg_open,
               ((prefix + mfile).c_str(), MODE_READ, &meshfn));
      cg_zone_read(meshfn, 1, 1, zName, &(sizes[0]));
      CG_CHECK(cg_close, (meshfn));
     } else {
     sizes[0] = 1;
     }
  }
  // MS End
  CG_CHECK(cg_zone_find_or_create, (fn, B, zoneName.c_str(), sizes,
                                    zType, &Z, errorhandle));
  //std::cout << __FILE__ << __LINE__ << std::endl;

  // Create the name for the GridCoordinates_t node.
  std::string gridName("Grid");
  if (timeLevel.length() + gridName.length() > 32)
    gridName.erase(32 - timeLevel.length());
  gridName += timeLevel;

  // Create the name for the IntegralData_t node for window attributes.
  std::string winName("WinData");
  if (timeLevel.length() + winName.length() > 32)
    winName.erase(32 - timeLevel.length());
  winName += timeLevel;

  // Create the name for the IntegralData_t node for pane attributes.
  std::string paneName("PaneData");
  if (timeLevel.length() + paneName.length() > 32)
    paneName.erase(32 - timeLevel.length());
  paneName += timeLevel;

  // Create the name for the IntegralData_t node for conn attributes.
  std::string connName("ConnData");
  if (timeLevel.length() + connName.length() > 32)
    connName.erase(32 - timeLevel.length());
  connName += timeLevel;

  // Create the name for the FlowSolution_t node for element attributes.
  std::string elemName("ElemData");
  if (timeLevel.length() + elemName.length() > 32)
    elemName.erase(32 - timeLevel.length());
  elemName += timeLevel;

  // Create the name for the FlowSolution_t node for node attributes.
  std::string nodeName("NodeData");
  if (timeLevel.length() + nodeName.length() > 32)
    nodeName.erase(32 - timeLevel.length());
  nodeName += timeLevel;

  // Read and update the grid coordinate and solution names in the
  // ZoneIterativeData_t node.
  char (*gridNames)[32] = NULL;
  char (*nodeNames)[32] = NULL;
  int rank, size[3];
  if (mode > 0 && cg_ziter_read(fn, B, Z, buffer) == 0) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    //std::cout << " Going to " << "Base = " << B << " Zone = "<< Z << std::endl;
    CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "ZoneIterativeData_t", 1, "end"));
    int numStr = nSteps;
    if (timeIndex >= nSteps)
      ++numStr;

    gridNames = (char(*)[32])new char[32*numStr];
    std::fill_n((char*)gridNames, 32 * numStr, '\0');
    nodeNames = (char(*)[32])new char[32*numStr];
    std::fill_n((char*)nodeNames, 32 * numStr, '\0');

    int A;
    DataType_t dataType;
    // Add the grid coordinates name to the GridCoordinatesPointers array,
    // unless it's already there.
    CG_CHECK(cg_array_info_by_name,
             ("GridCoordinatesPointers", &A, &dataType, &rank, size,
               errorhandle));
    COM_assertion_msg((dataType == Character && rank == 2
                       && size[1] == numStr),
                    "GridCoordinatesPointers in ZoneIterativeData is corrupt");
    CG_CHECK(cg_array_read, (A, gridNames));
    if (timeIndex < nSteps) {
      COM_assertion_msg(gridName == gridNames[timeIndex],
                    "GridCoordinatesPointers in ZoneIterativeData is corrupt");
    } else {
      std::strncpy(gridNames[timeIndex], gridName.c_str(), 32);
    }

    // Add the flow solutions name to the FlowSolutionsPointers array,
    // unless it's already there.
    CG_CHECK(cg_array_info_by_name,
             ("FlowSolutionsPointers", &A, &dataType, &rank, size,
              errorhandle));
    COM_assertion_msg((dataType == Character && rank == 2
                       && size[1] == numStr),
                      "FlowSolutionsPointers in ZoneIterativeData is corrupt");
    CG_CHECK(cg_array_read, (A, nodeNames));
    if (timeIndex < nSteps) {
      COM_assertion_msg(nodeName == nodeNames[timeIndex],
                      "FlowSolutionsPointers in ZoneIterativeData is corrupt");
    } else {
      std::strncpy(nodeNames[timeIndex], nodeName.c_str(), 32);
    }
    //std::cout << __FILE__ << __LINE__ << std::endl;
  } else {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    //std::cout << "Writing ZoneIterativeData" << std::endl;
    cg_ziter_write(fn, B, Z, "ZoneIterativeData");

    gridNames = (char(*)[32])new char[32];
    std::fill_n((char*)gridNames, 32, '\0');
    std::strncpy(gridNames[timeIndex], gridName.c_str(), 32);

    nodeNames = (char(*)[32])new char[32];
    std::fill_n((char*)nodeNames, 32, '\0');
    std::strncpy(nodeNames[timeIndex], nodeName.c_str(), 32);
  }
  //std::cout << __FILE__ << __LINE__ << std::endl;
  //std::cout << " timeIndex = " << timeIndex
  //          << " nSteps  = " << nSteps << std::endl;
  if (timeIndex == nSteps) {
    //std::cout << " Writing under ZoneIterativeData " 
    //          << " timeIndex = " << timeIndex
    //          << " nSteps = " << nSteps << std::endl;
    size[0] = 32;
    size[1] = ++nSteps;

    // Write/rewrite the grid coordinate and solution names in the
    // ZoneIterativeData_t node.
    CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "ZoneIterativeData_t", 1, "end"));
    //std::cout << __FILE__ << __LINE__ << std::endl;
    CG_CHECK(cg_array_write,
             ("GridCoordinatesPointers", Character, 2, size, gridNames));
    //std::cout << __FILE__ << __LINE__ << std::endl;
    CG_CHECK(cg_array_write,
             ("FlowSolutionsPointers", Character, 2, size, nodeNames));
  }
  //std::cout << __FILE__ << __LINE__ << std::endl;

  delete[] (char*)gridNames;
  delete[] (char*)nodeNames;

  std::vector<char> buf;
  int mfn = fn, mB = B, mZ = Z, nGC = 0, G = 0, numItems, dType;
  int rind[6] = { 0, 0, 0, 0, 0, 0 };
  // Build a CGNS path to the current zone.
  std::string path("/"), range;
  path += material;
  path += '/' + zoneName + '/';

  // Write the grid coordinates and the connectivity tables, if any.
  /*
  std::cout << " COM::COM_MESH = " << COM::COM_MESH << "\n"
            << " COM::COM_PMESH = " << COM::COM_PMESH << "\n"
            << " COM::COM_CONN = " << COM::COM_CONN << "\n"
            << " COM::COM_NC1 = " << COM::COM_NC1 << "\n"
            << " COM::COM_NC3 = " << COM::COM_NC3 << "\n"
            << " COM::COM_NC = " << COM::COM_NC << "\n"
            << " COM::COM_ALL = " << COM::COM_ALL << "\n"
            << " attr->id()= " << attr->id()  << "\n"
            << " attr->fullname()= " << attr->fullname()  << "\n" 
            << " attr->location()= " << attr->location()  << "\n" 
            << " attr->is_panel()= " << attr->is_panel()  << "\n" 
            << " attr->is_windowed()= " << attr->is_windowed()  << "\n" 
            << " attr->is_elemental()= " << attr->is_elemental()  << "\n" 
            << " attr->is_nodal()= " << attr->is_nodal()  << "\n" 
            << " attr->size_of_real_items()= " << attr->size_of_real_items()  << "\n"
            << " attr->size_of_ghost_items()= " << attr->size_of_ghost_items()  << std::endl;
  */
  if (attr->id() == COM::COM_MESH || attr->id() == COM::COM_PMESH
      || attr->id() == COM::COM_CONN
      || (attr->id() >= COM::COM_NC1 && attr->id() <= COM::COM_NC3)
      || attr->id() == COM::COM_NC || attr->id() == COM::COM_ALL) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    // Find or create the grid coordinates and element tables.  If they're in a
    // different file, then we should create links from the main file to the
    // mesh file nodes.
    if (!mfile.empty()) {
      DEBUG_MSG("Writing to mesh file '" << mfile << "', mode == '"
                << (mode ? "append'" : "write'"));
      GetCurrentDir(cwd, sizeof(cwd));
      //std::cout << __FILE__ << __LINE__ << std::endl;
      CG_CHECK(cg_open,
               ((prefix + mfile).c_str(), mode ? MODE_MODIFY : MODE_WRITE, &mfn));
      //std::cout << __FILE__ << __LINE__ << std::endl;

      //std::cout << __FILE__ << __LINE__ << std::endl;
      CG_CHECK(cg_base_find_or_create, (mfn, material, cellDim, physDim, &mB,
                                        errorhandle));
      //std::cout << __FILE__ << __LINE__ << std::endl;

      CG_CHECK(cg_zone_find_or_create, (mfn, mB, zoneName.c_str(),
                                        sizes, zType, &mZ, errorhandle));
    }
    //std::cout << __FILE__ << __LINE__ << std::endl;

    // The GridCoordinates node is referenced by the GridCoordinatesPointers,
    // so we should create it even if we don't put any data in it.
    // The first grid written must always be names "GridCoordinates".
    // Since we are naming our data based on timeLevel, we'll need to use
    // a link for the first one.
    CG_CHECK(cg_ngrids, (mfn, mB, mZ, &nGC));

    if (nGC == 0) {
      CG_CHECK(cg_grid_write, (mfn, mB, mZ, "GridCoordinates", &G));
      if (mfile.empty()) {
        CG_CHECK(cg_goto, (mfn, mB, "Zone_t", mZ, "end"));
        CG_CHECK(cg_link_write, (gridName.c_str(), "",
                                 (path + "GridCoordinates").c_str()));
      } else {
        // Create links from the data file to the mesh file.
        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));
        CG_CHECK(cg_ngrids, (fn, B, Z, &nGC));
        if (nGC == 0) {
          CG_CHECK(cg_link_write, ("GridCoordinates", (mfile).c_str(),
                                   (path + "GridCoordinates").c_str()));
        }
        CG_CHECK(cg_link_write, (gridName.c_str(), (mfile).c_str(),
                                 (path + "GridCoordinates").c_str()));
      }
    } else {
      if (mfile.empty()) {
        CG_CHECK(cg_grid_write, (mfn, mB, mZ, gridName.c_str(), &G));
      } else {
        // Create links from the data file to the mesh file.
        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));
        CG_CHECK(cg_ngrids, (fn, B, Z, &nGC));
        if (nGC == 0) {
          CG_CHECK(cg_link_write, ("GridCoordinates", (mfile).c_str(),
                                   (path + "GridCoordinates").c_str()));
        }
        CG_CHECK(cg_link_write, (gridName.c_str(), (mfile).c_str(),
                                 (path + "GridCoordinates").c_str()));
      }
    }
    //std::cout << __FILE__ << __LINE__ << std::endl;

    if (attr->id() == COM::COM_MESH || attr->id() == COM::COM_PMESH
        || attr->id() == COM::COM_NC || attr->id() == COM::COM_ALL
        || (attr->id() >= COM::COM_NC1 && attr->id() <= COM::COM_NC3)) {
      //std::cout << __FILE__ << __LINE__ << std::endl;
      // Write the actual mesh data.
      if (G > 0) {
        CG_CHECK(cg_goto, (mfn, mB, "Zone_t", mZ, "GridCoordinates_t", G,
                           "end"));

        if (pane.is_structured()) {
          rank = cellDim;
          size[0] = pane.size_i();
          size[1] = pane.size_j();
          size[2] = (physDim == 3 && cellDim == 3 ? pane.size_k(): 1);
          numItems = size[0] * size[1] * size[2];
          std::fill_n(rind, 2 * physDim, pane.size_of_ghost_layers());
        } else {
          rank = 1;
          size[0] = numItems = pane.size_of_nodes();
          std::fill_n(rind, 6, 0);
          rind[1] = pane.size_of_ghost_nodes();
        }

        // Write rind data.
        if (writeGhost)
          CG_CHECK(cg_rind_write, (rind));

        // std::string ranges[physDim];
        std::string ranges[3];
        int n, start = 0, finish = physDim;
        if (attr->id() >= COM::COM_NC1 && attr->id() <= COM::COM_NC3) {
          start = attr->id() - COM::COM_NC1;
          finish = start + 1;
        }
        for (n=start; n<finish; ++n) {
          const Attribute* nc = pane.attribute(COM::COM_NC1 + n);
          dType = nc->data_type();

          // Check for null/empty/strided data.
          const void* pData = nc->pointer();
          bool isNull = false;
          if (numItems <= 0) {
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData = &(buf[0]);
            size[0] = size[1] = size[2] = 1;
            ranges[n] = "EMPTY";
          } else if (pData == NULL) {
            buf.resize(COM::Attribute::get_sizeof(dType,
                                                  std::max(numItems, 1)), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData = &(buf[0]);
            // size[0] = size[1] = size[2] = 1;
            isNull = true;
            ranges[n] = "NULL";
          } else {
            int stride = nc->stride();
            if (stride > 1) {
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
              for (i=0; i<numItems; ++i)
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData)[i*stride*dSize]), dSize);
              pData = &(buf[0]);
            }

#ifdef DEBUG_DUMP_PREFIX
            {
              std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + ".nc" + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)pData)[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)pData)[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)pData)[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)pData)[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            SwitchOnCOMDataType(dType, FindRange(rank, size, rind,
                                                 (const COM_TT*)pData,
                                                 ranges[n]));
          }
          label = "Coordinate";
          label += (char)('X' + n);
          if (writeGhost) {
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType), rank,
                                      size, pData));
          } else {
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType), rank,
                                           rind, size, pData));
          }
          DEBUG_MSG("Writing attribute '" << (char)('x' + n) << "-nc', id == "
                    << nc->id() << ", components == "
                    << nc->size_of_components() << ", location == '"
                    << nc->location() << '\'');
        }

        // Write dimensional exponents and ranges.
        for (n=start; n<finish; ++n) {
          CG_CHECK(cg_goto, (mfn, mB, "Zone_t", mZ, "GridCoordinates_t", G,
                             "DataArray_t", n + 1, "end"));
          CG_CHECK(cg_descriptor_write, ("Range", ranges[n].c_str()));
          DEBUG_MSG("Writing descriptor 'Range': '" << ranges[n] << '\'');
          if (!attr->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            DEBUG_MSG("Writing descriptor 'Units': '" << attr->unit()
                      << '\'');
            CG_CHECK(cg_descriptor_write, ("Units", attr->unit().c_str()));
          } else { // Assume meters.
            cg_exponents_as_string_write("m", errorhandle);
            DEBUG_MSG("Writing descriptor 'Units': 'm'");
            CG_CHECK(cg_descriptor_write, ("Units", "m"));
          }
        }
      }
    }

    if (pane.is_unstructured() &&
        (attr->id() == COM::COM_MESH || attr->id() == COM::COM_PMESH
         || attr->id() == COM::COM_CONN || attr->id() == COM::COM_ALL)) {
      // Now do the element tables.
      int nS;
      CG_CHECK(cg_nsections, (mfn, mB, mZ, &nS));

      if (nS == 0) {
        if (!mfile.empty())
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));

        std::vector<const Connectivity*> conns;
        pane.elements(conns);

        int S = 0;
        i = 0;
        std::vector<const Connectivity*>::const_iterator c;
        for (c=conns.begin(); c!=conns.end(); ++c) {
          int nItems = (*c)->size_of_items();
          int nGhost = (*c)->size_of_ghost_items();

          DEBUG_MSG("Writing attribute '" << (*c)->name() << "', nItems == "
                    << nItems << ", elemType == " << (*c)->element_type());
          if (nItems == 0) {
            label = "Empty" + (*c)->name();
            S = WriteElements(mfn, mB, mZ, (*c)->element_type(), nItems,
                              NULL, false, 0, (*c)->index_offset() + 1, label,
                              errorhandle);
          } else if ((*c)->pointer() == NULL) {
            label = "Null" + (*c)->name();
            S = WriteElements(mfn, mB, mZ, (*c)->element_type(), nItems,
                              NULL, false, 0, (*c)->index_offset() + 1, label,
                              errorhandle);
          } else {
          // Write out a table of real nodes.
            label = (*c)->name();
            S = WriteElements(mfn, mB, mZ, (*c)->element_type(), nItems,
                              (*c)->pointer(), ((*c)->stride() == 1),
                              (*c)->capacity(), (*c)->index_offset() + 1,
                              label, errorhandle);
          }

          // Write out a table of ghost nodes.
          if (nGhost > 0) {
            //std::cout << __FILE__ << __LINE__ << " Writing Gost node data " << std::endl; 
            int elemRind[2] = { 0, nGhost };
            CG_CHECK(cg_goto, (mfn, mB, "Zone_t", mZ, "Elements_t", S, "end"));
            CG_CHECK(cg_rind_write, (elemRind));
          }
          
          // Link from the real file to the mesh file.
          if (!mfile.empty()) {
             /*
             std::cout << __FILE__ << __LINE__ 
                       << " label =" << label.c_str()
                       << std::endl;
             std::cout << __FILE__ << __LINE__ 
                       << " mfile =" << mfile.c_str()
                       << std::endl;
             std::cout << __FILE__ << __LINE__ 
                       << " path + label =" << (path + label).c_str()
                       << std::endl;
             */
             CG_CHECK(cg_link_write, (label.c_str(), (mfile).c_str(),
                                     (path + label).c_str()));
          }
          
        }
      }
    }
  } else if (!mfile.empty()) {
    // If we have a mesh file, we must link to the eternal GridCoordinates_t
    // and Elements_t nodes, even if attr->id() doesn't include mesh data.
    // First link the GridCoordinates_t node.
    //std::cout << __FILE__ << __LINE__ << std::endl; 
    CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));
    CG_CHECK(cg_ngrids, (fn, B, Z, &nGC));
       
    if (nGC == 0) {
      CG_CHECK(cg_link_write, ("GridCoordinates", (mfile).c_str(),
                               (path + "GridCoordinates").c_str()));
    }
    if (nGC == 0) {
    if (fname_in.compare(mfile)!=0){
       //std::cout << __FILE__ << __LINE__ << std::endl; 
       CG_CHECK(cg_link_write, (gridName.c_str(), (mfile).c_str(),
				(path + "GridCoordinates").c_str()));
    } else {
       //std::cout << __FILE__ << __LINE__ << std::endl; 
       CG_CHECK(cg_link_write, (gridName.c_str(), "",
				(path + "GridCoordinates").c_str()));
    }
    
    // Then link any Elements_t nodes.
    if (pane.is_unstructured()) {
      CG_CHECK(cg_open,
               ((prefix+mfile).c_str(), MODE_MODIFY, &mfn));

      CG_CHECK(cg_base_find_or_create, (mfn, material, cellDim, physDim, &mB,
                                        errorhandle));

      //std::cout << __FILE__ << __LINE__ << std::endl;
      CG_CHECK(cg_zone_find_or_create, (mfn, mB, zoneName.c_str(),
                                        sizes, zType, &mZ, errorhandle));

      CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));

      int S, nS;
      CG_CHECK(cg_nsections, (mfn, mB, mZ, &nS));

      char label[33];
      ElementType_t eType;
      int min, max, nBoundary, parent;
      for (S=1; S<=nS; ++S) {
        CG_CHECK(cg_section_read, (mfn, mB, mZ, S, label, &eType, &min, &max,
                                   &nBoundary, &parent));

        CG_CHECK(cg_link_write, (label, (mfile).c_str(), (path + label).c_str()));
      }
    }

    }

  }
  //std::cout << __FILE__ << __LINE__ << std::endl;

  if (mfn != fn){
    //std::cout << __FILE__ << __LINE__ << std::endl;
    CG_CHECK(cg_close, (mfn));
  }

  // This node is referenced in the FlowSolutionPointers, so we should make
  // sure it exists even if its left empty.
  int T;
  //std::cout << " nodeName = " << nodeName << std::endl;
  //std::cout << " Vertex = " << Vertex << std::endl;
  CG_CHECK(cg_sol_find_or_create, (fn, B, Z, nodeName.c_str(), Vertex, &T,
                                   errorhandle));
  //std::cout << __FILE__ << __LINE__ << std::endl;

  if (attr->id() == COM::COM_CONN || attr->id() == COM::COM_NC
      || (attr->id() >= COM::COM_NC1 && attr->id() <= COM::COM_NC3))
    return;
  //std::cout << __FILE__ << __LINE__ << std::endl;

  // Write out attributes.  Window attributes should be stored in IntegralData
  // nodes under the base node.  Pane, element, and node attributes should be
  // stored under the zone node, pane attrs in IntegralData nodes, element and
  // node attributes in FlowSolution nodes.

  std::vector<const Attribute*> attrs;
  if (attr->id() != COM::COM_MESH && attr->id() != COM::COM_PMESH
      && attr->id() != COM::COM_PCONN && attr->id() != COM::COM_RIDGES)
    pane.attributes(attrs);
  // MS Listing existing attributes in the pane
  /*
  std::cout << " //////////////////////////////////////////////////// " << std::endl;
  std::vector<const Attribute*> paneAttrs;
  pane.attributes(paneAttrs);
  std::vector<const Attribute*>::const_iterator ii;
  std::cout << " Listing existing attributes in the pane : " << std::endl;
  for (ii = paneAttrs.begin(); ii!=paneAttrs.end(); ii++)
     std::cout << " Attribute = " << (*ii)->fullname() << std::endl;
  std::cout << " //////////////////////////////////////////////////// " << std::endl;
  */
  // MS End
  //std::cout << __FILE__ << __LINE__ << std::endl;

  if (attr->id() == COM::COM_ALL || attr->id() == COM::COM_PMESH) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    attrs.insert(attrs.begin(), pane.attribute(COM::COM_RIDGES));
    attrs.insert(attrs.begin(), pane.attribute(COM::COM_PCONN));
  } else if (attr->id() == COM::COM_MESH) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    attrs.insert(attrs.begin(), pane.attribute(COM::COM_RIDGES));
  } else if (COM::COM_PCONN) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    attrs.insert(attrs.begin(), pane.attribute(COM::COM_PCONN));
  } else if (attr->id() == COM::COM_RIDGES) {
    //std::cout << __FILE__ << __LINE__ << std::endl;
    attrs.insert(attrs.begin(), pane.attribute(COM::COM_RIDGES));
  }
  //std::cout << "Number of attributes to write = " << attrs.size() << std::endl;
  //std::cout << __FILE__ << __LINE__ << std::endl;

  std::vector<const Attribute*>::const_iterator begin = attrs.begin();
  std::vector<const Attribute*>::const_iterator end = attrs.end();
  if ((attr->id() == COM::COM_MESH || attr->id() > COM::COM_ALL)  && attr->id() != COM::COM_PMESH
      && attr->id() != COM::COM_ATTS && attr->id() != COM::COM_ALL) {
    while (begin != end && (*begin)->id() != attr->id())
      ++begin;
    COM_assertion_msg(begin != end,
                      "Rocout::write_attribute: can't find attribute in pane");
    end = begin + 1;
  }

  std::vector<const Attribute*>::const_iterator a;
  // MS
  bool goToNextItem;
  char arrName[33];
  int tt2, tt1;
  DataType_t dT;
  // MS End
  for (a=begin; a!=end; ++a) {
    DEBUG_MSG("Writing attribute '" << (*a)->name() << "', id == "
              << (*a)->id() << ", components == " << (*a)->size_of_components()
              << ", location == '" << (*a)->location() << "\', datatype == "
              << (*a)->data_type());
    goToNextItem = false;
    int A, size[3], nComp = (*a)->size_of_components(), offset;
    const void* pData[9];
    dType = (*a)->data_type();
    switch ((*a)->location()) {
      case 'w':
        //std::cout << __FILE__ << __LINE__ << std::endl;
        size[0] = numItems = (*a)->size_of_items();
        size[1] = size[2] = 1;
        rind[0] = 0;
        rind[1] = (*a)->size_of_ghost_items();

        CG_CHECK(cg_goto, (fn, B, "end"));
        CG_CHECK(cg_integral_find_or_create, (winName.c_str(), &T,
                                              errorhandle));

        CG_CHECK(cg_goto, (fn, B, "IntegralData_t", T, "end"));
        CG_CHECK(cg_narrays, (&offset));
        //++offset;
        // MS
        // If the attribute we are trying to write is already written to
        // the file we just skip it. 
	//std::cout << __FILE__ << __LINE__ << std::endl;
	//std::cout << " Zone = " << Z << " IntegralData_t = " << T << std::endl;
	for (int ii=1; ii<=offset; ii++)
	{
	   cg_array_info(ii, arrName, &dT, &tt1, &tt2);
	   //std::cout << "offset = " << ii <<" arrayName = " << arrName << std::endl;
	   if ((std::string(arrName)).find(std::string((*a)->name())) != std::string::npos) {
	      //std::cout << "The request for duplicate dataitem is ignored." << std::endl;
	      goToNextItem = true;
	   }
	}
        if (goToNextItem)
           break;
        ++offset;

        for (A=0; A<nComp; ++A) {
          CG_CHECK(cg_goto, (fn, B, "IntegralData_t", T, "end"));

          label = (*a)->name();

          const Attribute* pa;
          if (nComp == 1)
            pa = pane.attribute((*a)->id());
          else {
            pa = pane.attribute((*a)->id() + A + 1);
            std::ostringstream sout;
            sout << label << '#' << A + 1 << "of" << nComp;
            label = sout.str();
          }

#ifdef DEBUG_DUMP_PREFIX
          std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + '.' + (*a)->name() + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
#endif // DEBUG_DUMP_PREFIX
          bool isNull = false;
          pData[A] = pa->pointer();
          if (numItems <= (writeGhost ? 0 : rind[1])) {
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            size[0] = 1;
            if (!writeGhost)
              rind[1] = 0;
            range = "EMPTY";
          } else if (pData[A] == NULL) {
            buf.resize(COM::Attribute::get_sizeof(dType,
                                                  std::max(size[0], 1)), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            isNull = true;
            size[0] = 1;
            range = "NULL";
          } else {
            int stride = pa->stride();
            if (stride > 1) {
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
              for (i=0; i<numItems; ++i)
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData[A])[i*stride*dSize]), dSize);
              pData[A] = &(buf[0]);
            }

#ifdef DEBUG_DUMP_PREFIX
            {
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)(pData[A]))[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)(pData[A]))[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)(pData[A]))[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)(pData[A]))[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            // SwitchOnCOMDataType(pa->data_type(),
            SwitchOnCOMDataType(dType,
                                FindRange(1, size, rind,
                                          (const COM_TT*)pData[A], range));
          }
          if (writeGhost) {
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType),
                                      1, size, pData[A]));
          } else {
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType),
                                           1, rind, size, pData[A]));
          }

          CG_CHECK(cg_goto, (fn, B, "IntegralData_t", T, "DataArray_t",
                             A + offset, "end"));

          CG_CHECK(cg_descriptor_write, ("Range", range.c_str()));
          DEBUG_MSG("Writing descriptor 'Range': '" << range << '\'');
          if (!pa->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            CG_CHECK(cg_descriptor_write, ("Units", pa->unit().c_str()));
            DEBUG_MSG("Writing descriptor 'Units': '" << pa->unit() << '\'');
          }

          // We can't put a Rind_t node under a DataArray_t node.
          if (writeGhost) {
            std::ostringstream sout;
            sout << rind[1];
            CG_CHECK(cg_descriptor_write, ("Ghost", sout.str().c_str()));
            DEBUG_MSG("Writing descriptor 'Ghost': '" << sout.str() << '\'');
          }
        }
        break;

      case 'p':
        size[0] = numItems = (*a)->size_of_items();
        size[1] = size[2] = 1;
        rind[0] = 0;
        rind[1] = (*a)->size_of_ghost_items();
        //std::cout << __FILE__ << __LINE__ << std::endl;
        DEBUG_MSG("numItems == " << numItems << ", numGhostItems == "
                  << rind[1]);

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));
        CG_CHECK(cg_integral_find_or_create, (paneName.c_str(), &T,
                                              errorhandle));

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T, "end"));
        CG_CHECK(cg_narrays, (&offset));
        //++offset;
        // MS
        // If the attribute we are trying to write is already written to
        // the file we just skip it. 
	//std::cout << __FILE__ << __LINE__ << std::endl;
	//std::cout << " Zone = " << Z << " IntegralData_t = " << T << std::endl;
	for (int ii=1; ii<=offset; ii++)
	{
	   cg_array_info(ii, arrName, &dT, &tt1, &tt2);
	   //std::cout << "offset = " << ii <<" arrayName = " << arrName << std::endl;
	   if ((std::string(arrName)).find(std::string((*a)->name())) != std::string::npos) {
	      //std::cout << "The request for duplicate dataitem is ignored." << std::endl;
	      goToNextItem = true;
	   }
	}
        if (goToNextItem)
           break;
        ++offset;
        //if (offset == 4)
        //   return;
        // MS End

        for (A=0; A<nComp; ++A) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T, "end"));

          label = (*a)->name();

          const Attribute* pa;
          if (nComp == 1)
            pa = pane.attribute((*a)->id());
          else {
            pa = pane.attribute((*a)->id() + A + 1);
            std::ostringstream sout;
            sout << label << '#' << A + 1 << "of" << nComp;
            label = sout.str();
          }

#ifdef DEBUG_DUMP_PREFIX
          std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + '.' + (*a)->name() + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
#endif // DEBUG_DUMP_PREFIX
          bool isNull = false;
          pData[A] = pa->pointer();
          if (numItems <= (writeGhost ? 0 : rind[1])) {
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            size[0] = 1;
            if (!writeGhost)
              rind[1] = 0;
            range = "EMPTY";
          } else if (pData[A] == NULL) {
            buf.resize(COM::Attribute::get_sizeof(dType,
                                                  std::max(size[0], 1)), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            isNull = true;
            size[0] = 1;
            range = "NULL";
          } else {
            int stride = pa->stride();
            if (stride > 1) {
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
              for (i=0; i<numItems; ++i)
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData[A])[i*stride*dSize]), dSize);
              pData[A] = &(buf[0]);
            }

#ifdef DEBUG_DUMP_PREFIX
            {
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)(pData[A]))[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)(pData[A]))[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)(pData[A]))[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)(pData[A]))[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            SwitchOnCOMDataType(dType,
                                FindRange(1, size, rind,
                                          (const COM_TT*)pData[A], range));
          }
          if (writeGhost) {
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType),
                                      1, size, pData[A]));
          } else {
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType),
                                           1, rind, size, pData[A]));
          }
          //std::cout << __FILE__ << __LINE__ 
          //          << " B = " << B << " Z = " << Z 
          //          << " T = " << T << " A = " << A
          //          << " offset = " << offset
          //          << std::endl;
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T,
                             "DataArray_t", A + offset, "end"));

          CG_CHECK(cg_descriptor_write, ("Range", range.c_str()));
          DEBUG_MSG("Writing descriptor 'Range' == '" << range << '\'');
          if (!pa->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            DEBUG_MSG("Writing descriptor 'Units' == '" << pa->unit() << '\'');
            CG_CHECK(cg_descriptor_write, ("Units", pa->unit().c_str()));
          }

          // We can't put a Rind_t node under a DataArray_t node.
          if (writeGhost) {
            std::ostringstream sout;
            sout << rind[1];
            CG_CHECK(cg_descriptor_write, ("Ghost", sout.str().c_str()));
            DEBUG_MSG("Writing descriptor 'Ghost' == '" << sout.str() << '\'');
          }
        }

        if ((*a)->id() == COM::COM_PCONN) {
          // Create CGNS connectivity node.
        } else if ((*a)->id() == COM::COM_RIDGES) {
          DEBUG_MSG("Special COM_RIDGES handling: zType == "
                    << (zType == Structured ? "Structured" : "Unstructured")
                    << ", size_of_components() == "
                    << (*a)->size_of_components() << ", size_of_items() == "
                    << (*a)->size_of_items());
          if (zType == Unstructured && (*a)->size_of_components() == 2
              && (*a)->size_of_items() > 0) {
            // Create parallel base and zone with BAR_2 elements.
            std::string newName(material);
            newName += "_ridges";
            int RB;
            DEBUG_MSG("Creating base \"" << newName
                      << "\", cellDim == 1, physDim == " << physDim);
            CG_CHECK(cg_base_find_or_create, (fn, newName.c_str(), 1, physDim,
                                              &RB, errorhandle));

            // Write the default Roccom units at the top.
            CG_CHECK(cg_goto, (fn, RB, "end"));
            CG_CHECK(cg_units_write, (Kilogram, Meter, Second, Kelvin, Degree));

            std::string rpath("/");
            rpath += material;
            rpath += '/';
            // Check to see if a BaseIterativeData node exists.
            if (cg_goto(fn, RB, "BaseIterativeData_t", 1, "end") != 0) {
              // Link to the main BaseIterativeData_t node.
              DEBUG_MSG("Linking TimeIterValues  to "
                        << rpath + "TimeIterValues");
              CG_CHECK(cg_goto, (fn, RB, "end"));
              CG_CHECK(cg_link_write, ("TimeIterValues", "",
                                       (rpath + "TimeIterValues").c_str()));
            }

            // Find or create the zone (corresponds to pane/section).  Use the
            // pane id as the zone name.  Set up the sizes based on
            // dimensionality and zone type.
            int RZ, rsizes[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            rsizes[0] = sizes[0];
            rsizes[1] = (*a)->size_of_real_items();
            DEBUG_MSG("Creating zone \"" << zoneName << "\", size == [ "
                      << rsizes[0] << ", " << rsizes[1] << " ]");
            //std::cout << __FILE__ << " line : " << __LINE__ << std::endl;
            CG_CHECK(cg_zone_find_or_create, (fn, RB, zoneName.c_str(), rsizes,
                                              Unstructured, &RZ, errorhandle));

            // Check to see if a ZoneIterativeData node exists.
            if (cg_goto(fn, RB, "Zone_t", RZ,
                        "ZoneIterativeData_t", 1, "end") != 0) {
              CG_CHECK(cg_goto, (fn, RB, "Zone_t", RZ, "end"));
              rpath += zoneName;
              rpath += '/';

              // Link to the main ZoneIterativeData_t node.
              DEBUG_MSG("Linking ZoneIterativeData  to "
                        << rpath + "ZoneIterativeData");
              CG_CHECK(cg_link_write, ("ZoneIterativeData", "",
                                       (rpath + "ZoneIterativeData").c_str()));

              // Link to the main GridCoordinates_t node.
              DEBUG_MSG("Linking GridCoordinates  to "
                        << rpath + "GridCoordinates");
              CG_CHECK(cg_link_write, ("GridCoordinates", "",
                                       (rpath + "GridCoordinates").c_str()));
              DEBUG_MSG("Linking " << gridName << "  to " << rpath + gridName);
              CG_CHECK(cg_link_write, (gridName.c_str(), "",
                                       (rpath + gridName).c_str()));

              // Link to the FlowSolution_t node with nodal data.
              DEBUG_MSG("Linking " << nodeName << "  to "
                        << rpath + nodeName);
              CG_CHECK(cg_link_write, (nodeName.c_str(), "",
                                       (rpath + nodeName).c_str()));
            }

            const Attribute* pa[2];
            pa[0] = pane.attribute((*a)->id() + 1);
            pa[1] = pane.attribute((*a)->id() + 2);
            DEBUG_MSG("Retreived the 'ridges' components: '" << pa[0]->name()
                      << "' and '" << pa[1]->name() << '\'');

            // Now do the BAR_2 element tables.
            int RS = 0;
            int nItems = pa[0]->size_of_items();
            int nGhost = pa[0]->size_of_ghost_items();
            DEBUG_MSG("Writing " << nItems << " BAR_2 elements (with "
                      << nGhost << " ghost elements)");
            COM_assertion(nItems == pa[1]->size_of_items());
            COM_assertion(nGhost == pa[1]->size_of_ghost_items());
            COM_assertion((pa[0]->pointer()==NULL) == (pa[1]->pointer()==NULL));

            if (!writeGhost)
              nItems = pa[0]->size_of_real_items();

            // Write out a table of real nodes.
            label = (*a)->name();
            std::vector<int> bar(2 * nItems);
            for (i=0; i<2; ++i) {
              int* ptr = (int*)pa[i]->pointer();
              int j, stride = std::max(pa[i]->stride(), 1);
              for (j=0; j<nItems; ++j)
                bar[2*j+i] = ptr[j*stride];
            }
            RS = WriteElements(mfn, RB, RZ, Connectivity::BAR2, nItems,
                               &bar[0], false, 0, 1, label, errorhandle);

            // Write out a table of ghost elements.
            if (writeGhost) {
              int elemRind[2] = { 0, nGhost };
              CG_CHECK(cg_goto,
                       (fn, RB, "Zone_t", RZ, "Elements_t", RS, "end"));
              CG_CHECK(cg_rind_write, (elemRind));
            }
          }
        } else if ((*a)->name() == "bcflag") {
          // Create CGNS boundary conditions node.
          // Boundary conditions might have other names, too.
        }
        break;

      case 'c':
        size[0] = numItems = (*a)->size_of_items();
        size[1] = size[2] = 1;
        rind[0] = 0;
        rind[1] = (*a)->size_of_ghost_items();

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "end"));
        CG_CHECK(cg_integral_find_or_create, (connName.c_str(), &T,
                                              errorhandle));

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T, "end"));
        CG_CHECK(cg_narrays, (&offset));
        ++offset;

        for (A=0; A<nComp; ++A) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T, "end"));

          label = (*a)->name();

          const Attribute* pa;
          if (nComp == 1)
            pa = pane.attribute((*a)->id());
          else {
            pa = pane.attribute((*a)->id() + A + 1);
            std::ostringstream sout;
            sout << label << '#' << A + 1 << "of" << nComp;
            label = sout.str();
          }

#ifdef DEBUG_DUMP_PREFIX
          std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + '.' + (*a)->name() + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
#endif // DEBUG_DUMP_PREFIX
          bool isNull = false;
          pData[A] = pa->pointer();
          if (numItems <= (writeGhost ? 0 : rind[1])) {
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            size[0] = 1;
            if (!writeGhost)
              rind[1] = 0;
            range = "EMPTY";
          } else if (pData[A] == NULL) {
            buf.resize(COM::Attribute::get_sizeof(dType,
                                                  std::max(size[0], 1)), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            isNull = true;
            size[0] = 1;
            range = "NULL";
          } else {
            int stride = pa->stride();
            if (stride > 1) {
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
              for (i=0; i<numItems; ++i)
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData[A])[i*stride*dSize]), dSize);
              pData[A] = &(buf[0]);
            }

#ifdef DEBUG_DUMP_PREFIX
            {
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)(pData[A]))[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)(pData[A]))[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)(pData[A]))[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)(pData[A]))[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            SwitchOnCOMDataType(dType,
                                FindRange(1, size, rind,
                                          (const COM_TT*)pData[A], range));
          }
          if (writeGhost) {
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType),
                                      1, size, pData[A]));
          } else {
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType),
                                           1, rind, size, pData[A]));
          }

          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "IntegralData_t", T,
                             "DataArray_t", A + offset, "end"));

          CG_CHECK(cg_descriptor_write, ("Range", range.c_str()));
          DEBUG_MSG("Writing descriptor 'Range': '" << range << '\'');
          if (!pa->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            CG_CHECK(cg_descriptor_write, ("Units", pa->unit().c_str()));
            DEBUG_MSG("Writing descriptor 'Units': '" << pa->unit() << '\'');
          }

          // We can't put a Rind_t node under a DataArray_t node.
          if (writeGhost) {
            std::ostringstream sout;
            sout << rind[1];
            CG_CHECK(cg_descriptor_write, ("Ghost", sout.str().c_str()));
            DEBUG_MSG("Writing descriptor 'Ghost': '" << sout.str() << '\'');
          }
        }
        break;

      case 'n':
        if (pane.is_structured()) {
          rank = cellDim;
          size[0] = pane.size_i();
          if (cellDim > 1) {
            size[1] = pane.size_j();
            size[2] = cellDim > 2 ? pane.size_k() : 1;
          } else {
            size[1] = size[2] = 1;
          }
          numItems = size[0] * size[1] * size[2];
          std::fill_n(rind, 2 * physDim, pane.size_of_ghost_layers());
        } else {
          rank = 1;
          size[0] = numItems = (*a)->size_of_items();
          size[1] = size[2] = 1;
          std::fill_n(rind, 6, 0);
          rind[1] = (*a)->size_of_ghost_items();
        }

        CG_CHECK(cg_sol_find_or_create, (fn, B, Z, nodeName.c_str(), Vertex, &T,
                                         errorhandle));

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T, "end"));
        CG_CHECK(cg_gridlocation_write, (Vertex));

        // Write rind data.
        if (writeGhost) {
          CG_CHECK(cg_rind_write, (rind));
          DEBUG_MSG("Ghost == " << rind[1]);
        }

        CG_CHECK(cg_narrays, (&offset));
        ++offset;

        for (A=0; A<nComp; ++A) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T, "end"));

          label = (*a)->name();

          const Attribute* pa;
          if (nComp == 1)
            pa = pane.attribute((*a)->id());
          else {
            pa = pane.attribute((*a)->id() + A + 1);
            if (nComp == physDim) {
              label += (char)('X' + (char)A);
            } else if (nComp == physDim * physDim) {
              label += (char)('X' + (char)(A / physDim));
              label += (char)('X' + (char)(A % physDim));
            } else {
              std::ostringstream sout;
              sout << label << '#' << A + 1 << "of" << nComp;
              label = sout.str();
            }
          }

#ifdef DEBUG_DUMP_PREFIX
          std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + '.' + (*a)->name() + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
#endif // DEBUG_DUMP_PREFIX
          bool isNull = false;
          pData[A] = pa->pointer();
          if (numItems <= 0) {
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            size[0] = size[1] = size[2] = 1;
            range = "EMPTY";
          } else if (pData[A] == NULL) {
            buf.resize(COM::Attribute::get_sizeof(dType, std::max(numItems, 1)),
                                                  '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            isNull = true;
            // size[0] = size[1] = size[2] = 1;
            range = "NULL";
          } else {
            int stride = pa->stride();
            if (stride > 1) {
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
              for (i=0; i<numItems; ++i) {
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData[A])[i*stride*dSize]), dSize);
              }
              pData[A] = &(buf[0]);
            }
#ifdef DEBUG_DUMP_PREFIX
            {
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)(pData[A]))[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)(pData[A]))[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)(pData[A]))[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)(pData[A]))[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            SwitchOnCOMDataType(dType,
                                FindRange(rank, size, rind,
                                          (const COM_TT*)pData[A], range));
          }
          if (writeGhost) {
            DEBUG_MSG("Calling cg_array_write( name == '" << label
                      << "', dataType == "
                      << (COM2CGNS(dType) == RealSingle ? "float"
                         : (COM2CGNS(dType) == RealDouble ? "double"
                         : (COM2CGNS(dType) == Character ? "char" : "int?")))
                      << ", rank == " << rank << ", size[] = { " << size[0]
                      << ", " << size[1] << ", " << size[2] << " } )");
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType),
                                      rank, size, pData[A]));
          } else {
            DEBUG_MSG("Calling cg_array_core_write( name == '" << label
                      << "', dataType == "
                      << (COM2CGNS(dType) == RealSingle ? "float"
                         : (COM2CGNS(dType) == RealDouble ? "double"
                         : (COM2CGNS(dType) == Character ? "char" : "int?")))
                      << ", rank == " << rank << ", rind[] = { " << rind[0] 
                      << ", " << rind[1] << ", " << rind[2] << ", " << rind[3]
                      << ", " << rind[4] << ", " << rind[5] << " }, size[] = { "
                      << size[0] << ", " << size[1] << ", " << size[2]
                      << " } )");
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType),
                                           rank, rind, size, pData[A]));
          }

          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", A + offset, "end"));

          DEBUG_MSG("Calling cg_descriptor_write( name == 'Range', "
                    << "value == '" << range << "' )");
          CG_CHECK(cg_descriptor_write, ("Range", range.c_str()));
          if (!pa->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            DEBUG_MSG("Calling cg_descriptor_write( name == 'Units', "
                      << "value == '" << pa->unit() << "' )");
            CG_CHECK(cg_descriptor_write, ("Units", pa->unit().c_str()));
          }
        }

        // Vectors and tensors need a MagnitudeRange or TraceRange
        // descriptor under the first DataArray_t node.
        if (nComp == 3) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", offset, "end"));
          SwitchOnCOMDataType(dType,
                              FindMagnitudeRange(physDim, rank, size, rind,
                                                 (const COM_TT**)pData,
                                                 range));
          CG_CHECK(cg_descriptor_write, ("MagnitudeRange", range.c_str()));
        } else if (nComp == 9) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", offset, "end"));
          SwitchOnCOMDataType(dType,
                              FindTraceRange(physDim, rank, size, rind,
                                             (const COM_TT**)pData, range));
          CG_CHECK(cg_descriptor_write, ("TraceRange", range.c_str()));
        }
        break;

      case 'e':
        if (pane.is_structured()) {
          rank = cellDim;
          size[0] = pane.size_i() - 1;
          if (cellDim > 1) {
            size[1] = pane.size_j() - 1;
            size[2] = cellDim > 2 ? pane.size_k() - 1 : 1;
          } else {
            size[1] = size[2] = 1;
          }
          numItems = size[0] * size[1] * size[2];
          std::fill_n(rind, 2 * physDim, pane.size_of_ghost_layers());
        } else {
          rank = 1;
          size[0] = numItems = (*a)->size_of_items();
          size[1] = size[2] = 1;
          std::fill_n(rind, 6, 0);
          rind[1] = (*a)->size_of_ghost_items();
        }

        CG_CHECK(cg_sol_find_or_create, (fn, B, Z, elemName.c_str(), CellCenter,
                                         &T, errorhandle));

        CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T, "end"));
        CG_CHECK(cg_gridlocation_write, (CellCenter));

        // Write rind data.
        if (writeGhost) {
          CG_CHECK(cg_rind_write, (rind));
          DEBUG_MSG("Ghost == " << rind[1]);
        }

        CG_CHECK(cg_narrays, (&offset));
        //++offset;
        // MS
        // If the attribute we are trying to write is already written to
        // the file we just skip it. 
	//std::cout << __FILE__ << __LINE__ << std::endl;
	//std::cout << " Zone = " << Z << " IntegralData_t = " << T << std::endl;
	//std::cout << " Read offset = " << offset << std::endl;
	for (int ii=1; ii<=offset; ii++)
	{
	   cg_array_info(ii, arrName, &dT, &tt1, &tt2);
	   //std::cout << "offset = " << ii <<" arrayName = " << arrName << std::endl;
	   if ((std::string(arrName)).find(std::string((*a)->name())) != std::string::npos) {
	      //std::cout << "The request for duplicate dataitem is ignored." << std::endl;
	      goToNextItem = true;
	   }
	}
        if (goToNextItem)
           break;
        ++offset;
        //if (offset == 4)
        //   return;
        // MS End
	//std::cout << __FILE__ << __LINE__ << std::endl;

        for (A=0; A<nComp; ++A) {
          //std::cout << __FILE__ << __LINE__ << " Component = " << A << std::endl;
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T, "end"));
	  //std::cout << __FILE__ << __LINE__ << std::endl;

          label = (*a)->name();
	  //std::cout << __FILE__ << __LINE__ << std::endl;

          const Attribute* pa;
          if (nComp == 1)
            pa = pane.attribute((*a)->id());
          else {
            pa = pane.attribute((*a)->id() + A + 1);
            if (nComp == physDim) {
              label += (char)('X' + (char)A);
            } else if (nComp == physDim * physDim) {
              label += (char)('X' + (char)(A / physDim));
              label += (char)('X' + (char)(A % physDim));
            } else {
              std::ostringstream sout;
              sout << label << '#' << A + 1 << "of" << nComp;
              label = sout.str();
            }
          }
	  //std::cout << __FILE__ << __LINE__ << std::endl;

#ifdef DEBUG_DUMP_PREFIX
          std::ofstream fout((DEBUG_DUMP_PREFIX + std::string(material) + '.' + (*a)->name() + '_' + timeLevel + ".cgns").c_str(), std::ios_base::app);
#endif // DEBUG_DUMP_PREFIX
          bool isNull = false;
          pData[A] = pa->pointer();
	  //std::cout << __FILE__ << __LINE__ << " numItems = " << numItems << std::endl;
          if (numItems <= 0) {
	    //std::cout << __FILE__ << __LINE__ << std::endl;
            buf.resize(COM::Attribute::get_sizeof(dType, 1), '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            size[0] = size[1] = size[2] = 1;
            range = "EMPTY";
          } else if (pData[A] == NULL) {
	    //std::cout << __FILE__ << __LINE__ << std::endl;
            buf.resize(COM::Attribute::get_sizeof(dType, std::max(numItems, 1)),
                                                  '\0');
            std::fill(buf.begin(), buf.end(), '\0');
            pData[A] = &(buf[0]);
            isNull = true;
            // size[0] = size[1] = size[2] = 1;
            range = "NULL";
          } else {
	    //std::cout << __FILE__ << __LINE__ << std::endl;
            int stride = pa->stride();
	    //std::cout << __FILE__ << __LINE__ << std::endl;
            if (stride > 1) {
	      //std::cout << __FILE__ << __LINE__ << std::endl;
              int dSize = COM::Attribute::get_sizeof(dType, 1);
              buf.resize(dSize * std::max(numItems, 1));
	      //std::cout << __FILE__ << __LINE__ << std::endl;
              for (i=0; i<numItems; ++i) {
                std::memcpy(&(buf[i*dSize]),
                            &(((const char*)pData[A])[i*stride*dSize]), dSize);
              }
	      //std::cout << __FILE__ << __LINE__ << std::endl;
              pData[A] = &(buf[0]);
            }
#ifdef DEBUG_DUMP_PREFIX
            {
              switch (COM2CGNS(dType)) {
                case Character:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << (int)((const char*)(pData[A]))[i]
                         << '\n';
                    break;
                case Integer:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const int*)(pData[A]))[i] << '\n';
                  break;
                case RealSingle:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const float*)(pData[A]))[i] << '\n';
                  break;
                case RealDouble:
                  for (i=0; i<numItems; ++i)
                    fout << i << " : " << ((const double*)(pData[A]))[i]
                         << '\n';
                  break;
                default:
                  break;
              }
              fout << "###########################################\n";
            }
#endif // DEBUG_DUMP_PREFIX

            SwitchOnCOMDataType(dType,
                                FindRange(rank, size, rind,
                                          (const COM_TT*)pData[A], range));
          }
          if (writeGhost) {
            DEBUG_MSG("Calling cg_array_write( name == '" << label
                      << "', dataType == "
                      << (COM2CGNS(dType) == RealSingle ? "float"
                         : (COM2CGNS(dType) == RealDouble ? "double"
                         : (COM2CGNS(dType) == Character ? "char" : "int?")))
                      << ", rank == " << rank << ", size[] = { " << size[0]
                      << ", " << size[1] << ", " << size[2] << " } )");
            CG_CHECK(cg_array_write, (label.c_str(), COM2CGNS(dType),
                                      rank, size, pData[A]));
          } else {
            DEBUG_MSG("Calling cg_array_core_write( name == '" << label
                      << "', dataType == "
                      << (COM2CGNS(dType) == RealSingle ? "float"
                         : (COM2CGNS(dType) == RealDouble ? "double"
                         : (COM2CGNS(dType) == Character ? "char" : "int?")))
                      << ", rank == " << rank << ", rind[] = { " << rind[0] 
                      << ", " << rind[1] << ", " << rind[2] << ", " << rind[3]
                      << ", " << rind[4] << ", " << rind[5] << " }, size[] = { "
                      << size[0] << ", " << size[1] << ", " << size[2]
                      << " } )");
            CG_CHECK(cg_array_core_write, (label.c_str(), COM2CGNS(dType),
                                           rank, rind, size, pData[A]));
          }

          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", A + offset, "end"));

          DEBUG_MSG("Calling cg_descriptor_write( name == 'Range', "
                    << "value == '" << range << "' )");
          CG_CHECK(cg_descriptor_write, ("Range", range.c_str()));
          if (!pa->unit().empty()) {
            cg_exponents_as_string_write(attr->unit().c_str(), errorhandle);
            DEBUG_MSG("Calling cg_descriptor_write( name == 'Units', "
                      << "value == '" << pa->unit() << "' )");
            CG_CHECK(cg_descriptor_write, ("Units", pa->unit().c_str()));
          }
        }
	//std::cout << __FILE__ << __LINE__ << std::endl;

        // Vectors and tensors need a MagnitudeRange or TraceRange
        // descriptor under the first DataArray_t node.
        if (nComp == 3) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", offset, "end"));
          SwitchOnCOMDataType(dType,
                              FindMagnitudeRange(physDim, rank, size, rind,
                                                 (const COM_TT**)pData,
                                                 range));
          CG_CHECK(cg_descriptor_write, ("MagnitudeRange", range.c_str()));
        } else if (nComp == 9) {
          CG_CHECK(cg_goto, (fn, B, "Zone_t", Z, "FlowSolution_t", T,
                             "DataArray_t", offset, "end"));
          SwitchOnCOMDataType(dType,
                              FindTraceRange(physDim, rank, size, rind,
                                             (const COM_TT**)pData, range));
          CG_CHECK(cg_descriptor_write, ("TraceRange", range.c_str()));
        }
        break;
    }
  }
}






