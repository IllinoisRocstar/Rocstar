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

/*******************************************************************
 *
 *  RocMeshData.h
 *
 *  Declaration for class RocMeshData. RocMeshData contains mesh  
 *  data and ghost information from Rocstar's HDF file. 
 *
 *  COM_init() needed to be called before using this class.
 *
 *
 *  Pornput Suriyamongkol 
 *  10/10/06
 *
 *
 *  Modification
 *  ------------
 *  Add field information - Pornput 11/02/06
 *
 *******************************************************************/

#ifndef _ROCMESHDATA_
#define _ROCMESHDATA_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::cerr;
using std::setw;
using std::string;
using std::vector;
using std::ofstream;
using std::ios;

#include "com.h"

COM_EXTERN_MODULE(SimIN);


class RocMeshData
{
  /****************************************************************
   *
   * Member variables
   *
   ***************************************************************/

 private:

  // Mesh Info
  string   _filename;  
  int      _indexBase;
  int      _numVer;
  int      _numVerGhost;
  int      _numElem;
  int      _elemDataSize;
  int      _pconnSize;
  int *    _elemData;
  int *    _elemType;
  int *    _pconnArray;
  double * _coords;

  // Field Info
  int               _numFields;
  vector <string>   _fieldNames;
  vector <string>   _fieldUnits;
  vector <char>     _fieldTypes;
  vector <int>      _fieldComp;
  vector <double *> _fieldVals;


  /****************************************************************
   *
   *  Public Methods
   *
   ***************************************************************/

 public:
  //
  // Default constructor
  //
  RocMeshData(); 

  //
  // Destructor
  //
  ~RocMeshData(); 

  //
  // Constructor
  //
  RocMeshData(string const & filename); 

  //
  // Constructor (overloaded)
  //
  RocMeshData(char const * filename); 

  //
  // Copy constructor
  //
  RocMeshData(RocMeshData const & origVal); 

  //
  // Assignment operator
  //
  RocMeshData & operator=(RocMeshData const & origVal);

  //
  // Show all RocMeshData info
  //
  void showAllInfo() const; 

  //
  // Check pconn. Return 0 if ok
  //
  int pconnCheck(string pltFile = "") const;

  //
  // Set base index to 0
  //
  void set0IndexBase();

  //
  // Set base index to 1
  //
  void set1IndexBase();

  //
  // Return index base
  //
  int getIndexBase() const;

  //
  // Return number of vertices
  //
  int numVer() const;

  //
  // Return number of ghost vertices
  //
  int numVerGhost() const;

  //
  // Return number of elements
  //
  int numElem() const;

  // 
  // Return pointer to coordinates
  //
  double * getCoords() const;

  // 
  // Return pointer to element types
  //
  int * getElemType() const;

  //
  // Return pointer to element connectivity
  //
  int * getElemData() const;

  //
  // Return PCONN array
  //
  int * getPCONN() const;

  //
  // Return number of fields
  //
  int numField() const;

  //
  // Return field types
  //
  vector <char> getFieldTypes() const;

  //
  // Return field names
  //
  vector <string> getFieldNames() const;

  //
  // Return field units
  //
  vector <string> getFieldUnits() const;

  //
  // Return field numbers of components
  //
  vector <int> getFieldNumComp() const;

  //
  // Return field value
  //
  vector <double *> getFieldData() const;

  /****************************************************************
   *
   * Helper functions
   *
   ***************************************************************/
 private:

  //
  // Free all dynamic memory
  //
  void free(); 

  //
  // Copy all member variables
  //
  void copy(RocMeshData const & origVal); 

  //
  // Initialize RocMeshData in constructors
  //
  void initRocMeshData(char const * filename);

  //
  // Load data from Rocstar's HDF file to window
  //
  void loadData(char const * hdfFilename, const char * winName);

  //
  // Get pane id and connectivity table info
  //
  void getPaneInfo(const char * winName, 
		   int & paneId, int & numCon, char * & conNames);

  //
  // Get coordinates
  //
  void getCoord(const char * winName, int paneId);

  //
  // Get total number of elements
  //
  void getElemInfo(const char * winName, 
		   int paneId, int numCon, char * conNames, 
		   int * & types, int * & numElemArray);

  //
  // Allocate memory for elemType and elemData
  //
  void allocElemMem(int numCon, int * numElemArray, int * types);

  //
  // Get element connectivity
  //
  void getElemConn(const char * winName, int paneId, int numCon, 
		   char * conNames, int * types, int * numElemArray);

  //
  // Get PCONN array and its size
  //
  void getPconnArray(const char * winName, int paneId);

  //
  // Get field information
  //
  void getFieldInfo(const char * winName, int paneId);

  //
  // Get field attributes and allocate field data memory
  //
  void getFieldAttributes(const char * winName);

  //
  // Get field data 
  //
  void getFieldData(const char * winName, int paneId);

  //
  // Get Simmetrix element type number
  //
  int getSimElemTypeNum(string type);

  //
  // Get #nodes / element given element type in Simmetrix format
  //
  int getNodePerElem(int num) const;
};

#endif /* _ROCMESHDATA_ */






