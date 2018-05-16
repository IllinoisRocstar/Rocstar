/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Version: Sandia Evaluation Distribution                           *
 * Licensed To: Sandia National Laboratories                         *
 * License Type: Evaluation                                          *
 * License Expiration: March 13, 2013                                *
 *********************************************************************/
/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2012, IllinoisRocstar LLC. All rights reserved.         *
 *                                                                   *
 * The Rocstar Simulation Suite is the property of IllinoisRocstar   *
 * LLC. No use or distribution of this version of the Rocstar        *
 * Simulation Suite beyond the license provided through separate     *
 * contract is permitted.                                            *
 *                                                                   *
 * IllinoisRocstar LLC                                               *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *********************************************************************/
/* *******************************************************************
 *  Initial open source Rocstar software developed by                *
 *     Center for Simulation of Advanced Rockets                     *
 *     University of Illinois at Urbana-Champaign                    *
 *     Urbana, IL                                                    *
 *     www.csar.uiuc.edu                                             *
 *                                                                   *
 * Copyright@2008, University of Illinois.  All rights reserved.     *
 *                                                                   *
 * @ Redistributions of source code must retain the above copyright  * 
 *   notice, this list of conditions and the following disclaimers.  *
 *                                                                   * 
 * @ Redistributions in binary form must reproduce the above         *
 *   copyright notice, this list of conditions and the following     *
 *   disclaimers in the documentation and/or other materials         *
 *   provided with the distribution.                                 *
 *                                                                   *
 * @ Neither the names of the Center for Simulation of Advanced      *
 *   Rockets, the University of Illinois, nor the names of its       *
 *   contributors may be used to endorse or promote products derived * 
 *   from this Software without specific prior written permission.   *
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
 *  RocMeshData.cc
 *
 *  Definitions for class RocMeshData. RocMeshData contains mesh  
 *  data and ghost information from Rocstar's HDF file. 
 *
 *  COM_init() need to be called before using this class.
 *
 *  Pornput Suriyamongkol 
 *  10/10/06
 *
 *
 *  Modification
 *  ------------
 *  Add field information - Pornput 11/02/06
 *  Add pconnCheck()      - Pornput 04/20/07
 *
 *******************************************************************/

#include "RocMeshData.h"
#include <cstring>
#include <cstdlib>


/*******************************************************************
 *
 *  RocMeshData:
 *
 *    Constructor with no argument
 *
 *******************************************************************/

RocMeshData::RocMeshData()
{
  _filename = "";
  _indexBase = 0;
  _numVer = 0;
  _numVerGhost = 0;
  _numElem = 0;
  _elemDataSize = 0;
  _pconnSize = 0;
  _elemData = NULL;
  _elemType = NULL;
  _pconnArray = NULL;
  _coords = NULL;

  _numFields = 0;
  _fieldNames.resize(0);
  _fieldUnits.resize(0);
  _fieldTypes.resize(0);
  _fieldComp.resize(0);
  _fieldVals.resize(0);
}


/*******************************************************************
 *
 *  ~RocMeshData:
 *
 *    Destructor
 *
 *******************************************************************/

RocMeshData::~RocMeshData()
{
  free();
}


/*******************************************************************
 *
 *  RocMeshData:
 *
 *    Constructor with a filename as an argument
 *
 *******************************************************************/

RocMeshData::RocMeshData(string const & filename)
{
  // filename in chars
  char hdfFilename [40];
  strcpy(hdfFilename, filename.c_str());

  // initialize RocMeshData
  initRocMeshData(hdfFilename);
}


/*******************************************************************
 *
 *  RocMeshData: (Overloaded)
 *
 *    Constructor with a filename as an argument
 *
 *******************************************************************/

RocMeshData::RocMeshData(char const * filename)
{
  initRocMeshData(filename);
}


/*******************************************************************
 *
 *  RocMeshData:
 *
 *    Copy constructor
 *
 *******************************************************************/

RocMeshData::RocMeshData(RocMeshData const & origVal)
{
  copy(origVal);
}


/*******************************************************************
 *
 *  Operator=:
 *
 *    Assignment operator overloaded
 *
 *******************************************************************/

RocMeshData & RocMeshData::operator=(RocMeshData const & origVal)
{
  if (this != &origVal)
  {
    free();                // free all memory
    copy(origVal);         // copy
  }
  return *this;
}


/*******************************************************************
 *
 *  showAllInfo:
 *
 *    Show all Info of RocMeshData
 *
 *******************************************************************/

void RocMeshData::showAllInfo() const
{
  cout << endl << endl;

  cout << "Filename: " << _filename << endl;
  cout << "Index base: " << _indexBase << endl;
  cout << setw(10) << _numVer << " vertices" << endl;
  cout << setw(10) << _numElem << " elements" << endl;

  cout << endl << "Coordinates" << endl;
  for (int i = 0; i < (_numVer + _numVerGhost); i++)
  {
    cout << "node " << setw(5) << i << ": "
	 << setw(10) << _coords[3*i + 0] << ", "
	 << setw(10) << _coords[3*i + 1] << ", "
	 << setw(10) << _coords[3*i + 2] << endl;
  }

  cout << endl << "Connectivity" << endl;
  int counter = 0;
  for (int i = 0; i < _numElem; i++)
  {
    cout << "elem " << setw(5) << i << ": "
	 << "type " << _elemType[i] << ": ";
    
    int nodePerElem = getNodePerElem(_elemType[i]);
    for (int j = 0; j < nodePerElem; j++)
    {
      cout << setw(10) << _elemData[counter] << ", ";
      counter++;
    }
    cout << endl;
  }

  cout << endl << "Partition Connectivity" << endl;
  for (int i = 0; i < _pconnSize; i++)
    cout << setw(10) << _pconnArray[i] << endl;

  cout << endl << endl;
}


/*******************************************************************
 *
 *  pconnCheck:
 *
 *    Check PCONN block 3.  Return 0 if ok.
 *
 *******************************************************************/

int RocMeshData::pconnCheck(string pltFile) const
{
  int result = 0;
  bool OORGhostNodes  = false;    // found out of range ghost node
  bool missGhostNodes = false;    // ghost nodes missed from pconn
  int minGhostNodeID = 0;
  int maxGhostNodeID = 0;

  minGhostNodeID = _numVer + 1;
  maxGhostNodeID = _numVer + _numVerGhost;

  // Allocate marking array
  int * isMarked = new int [_numVerGhost];
  for (int i = 0; i < _numVerGhost; i++)
    isMarked[i] = 0;

  // Jump to pconn block 3
  int index = 0;
  int numNbrs = 0;
  int numItems = 0;

  for (int i = 0; i < 2; i++)
  {
    numNbrs = _pconnArray[index];
    index++;

    for (int j = 0; j < numNbrs; j++)
    {
      index++;                         // step over proc id
      numItems = _pconnArray[index];
      index++;
      index += numItems;
    }
  }

  // Begin marking
  numNbrs = _pconnArray[index];
  index++;

  for (int j = 0; j < numNbrs; j++)
  {
    index++;                         // step over proc id
    numItems = _pconnArray[index];
    index++;
  
    for (int k = 0; k < numItems; k++)
    {
      int lid = _pconnArray[index];
      index++;

      if (lid < minGhostNodeID || lid > maxGhostNodeID)
      {
	cout << "ERROR: Ghost node " << lid 
	     << " is out of ghost node range.\n";
	OORGhostNodes = true;
      }
      else if (isMarked[lid - minGhostNodeID] == 1)
      {
	cout << "WARNING: Ghost node " << lid 
	     << " is listed in pconn block 3"
	     << " more than once.\n";	
      }
      else
	isMarked[lid - minGhostNodeID] = 1;
    }
  }

  // Check if any node is not marked. 
  cout << "\nGhost nodes not listed in pconn block 3.\n";
  cout << "-----------------------------------------------\n\n";
  int numGNodesMissed = 0;
  for (int i = 0; i < _numVerGhost; i++)
    if (isMarked[i] == 0)
    {
      missGhostNodes = true;
      numGNodesMissed++;
      int lid = (i + minGhostNodeID);
      cout << "ERROR: Ghost node " << lid
	   << ": " << setw(10) << *(_coords + 3 * (lid - 1)) 
	   << ", " << setw(10) << *(_coords + 3 * (lid - 1) + 1)
	   << ", " << setw(10) << *(_coords + 3 * (lid - 1) + 2)
	   << endl;
    }

  // Write out to Tecplot file if needed.
  if ( (pltFile.compare("") != 0) && (numGNodesMissed > 0) )
  {
    ofstream outFile(pltFile.c_str(), ios::out);
    outFile << "TITLE=" << pltFile << endl;
    outFile << "VARIABLES=\"x\", \"y\", \"z\"" << endl;
    outFile << "ZONE I=" << numGNodesMissed << " F=POINT\n";
    for (int i = 0; i < _numVerGhost; i++)
      if (isMarked[i] == 0)
      {
	int lid = (i + minGhostNodeID);
	outFile << setw(10) << *(_coords + 3 * (lid - 1)) 
		<< " " << setw(10) << *(_coords + 3 * (lid - 1) + 1)
		<< " " << setw(10) << *(_coords + 3 * (lid - 1) + 2)
		<< endl;
    }
    outFile.close();
  }

  if (numGNodesMissed == 0)
    cout << "None.\n\n";
  else
    cout << endl;

  // Free up memory
  delete isMarked;

  // Return the result
  if (OORGhostNodes == false && missGhostNodes == false)
    return 0;
  else
    return 1;
}


/*******************************************************************
 *
 *  set0IndexBase:
 *
 *    Set index base to 0
 *
 *******************************************************************/

void RocMeshData::set0IndexBase()
{
  if (_indexBase == 0)         // return if already 0
    return;

  if (_indexBase != 1)         // return if base is not 1
  {
    cerr << "\n\nWARNING: Unexpected index base "
	 << _indexBase << ", Nothing Done!\n\n";
    return;
  }

  _indexBase = 0;              // set new index base

  for (int i = 0; i < _elemDataSize; i++)
    _elemData[i]--;
  
  int iter = 0;                // iterator of _pconnArray
  while (iter < _pconnSize)
  {
    int numPar = _pconnArray[iter];
    iter++;
    
    for (int i = 0; i < numPar; i++)
    {
      iter++;
      int numItem = _pconnArray[iter];
      iter++;

      for (int j = 0; j < numItem; j++)
      {
	_pconnArray[iter]--;
	iter++;
      }
    }
  }
}


/*******************************************************************
 *
 *  set1IndexBase:
 *
 *    Set index base to 1
 *
 *******************************************************************/

void RocMeshData::set1IndexBase()
{
  if (_indexBase == 1)         // return if already 1
    return;

  if (_indexBase != 0)         // return if base is not 0
  {
    cerr << "\n\nWARNING: Unexpected index base "
	 << _indexBase << ", Nothing Done!\n\n";
    return;
  }

  _indexBase = 1;              // set new index base

  for (int i = 0; i < _elemDataSize; i++)
    _elemData[i]++;
  
  int iter = 0;                // iterator of _pconnArray
  while (iter < _pconnSize)
  {
    int numPar = _pconnArray[iter];
    iter++;
    
    for (int i = 0; i < numPar; i++)
    {
      iter++;
      int numItem = _pconnArray[iter];
      iter++;

      for (int j = 0; j < numItem; j++)
      {
	_pconnArray[iter]++;
	iter++;
      }
    }
  }
}


/*******************************************************************
 *
 *  getIndexBase:
 *
 *    Return index base
 *
 *******************************************************************/

int RocMeshData::getIndexBase() const
{
  return _indexBase;
}


/*******************************************************************
 *
 *  numVer:
 *
 *    Return number of vertices
 *
 *******************************************************************/

int RocMeshData::numVer() const
{
  return _numVer;
}

/*******************************************************************
 *
 *  numVerGhost:
 *
 *    Return number of ghost vertices
 *
 *******************************************************************/

int RocMeshData::numVerGhost() const
{
  return _numVerGhost;
}


/*******************************************************************
 *
 *  numElem:
 *
 *    Return number of elements
 *
 *******************************************************************/
int RocMeshData::numElem() const
{
  return _numElem;
}


/*******************************************************************
 *
 *  getCoords:
 *
 *    Return pointer to coordinates
 *
 *******************************************************************/

double * RocMeshData::getCoords() const
{
  return _coords;
}


/*******************************************************************
 *
 *  getElemType:
 *
 *    Return pointer to element types
 *
 *******************************************************************/

int * RocMeshData::getElemType() const
{
  return _elemType;
}


/*******************************************************************
 *
 *  getElemData:
 *
 *    Return pointer to element connectivity
 *
 *******************************************************************/

int * RocMeshData::getElemData() const
{
  return _elemData;
}

/*******************************************************************
 *
 *  getPCONN:
 *
 *    Return PCONN array
 *
 *******************************************************************/

int * RocMeshData::getPCONN() const
{
  return _pconnArray;
}


/*******************************************************************
 *
 *  numField:
 *
 *    Return number of fields
 *
 *******************************************************************/

int RocMeshData::numField() const
{
  return _numFields;
}


/*******************************************************************
 *
 *  getFieldTypes:
 *
 *    Return field types
 *
 *******************************************************************/

vector <char> RocMeshData::getFieldTypes() const
{
  return _fieldTypes;
}


/*******************************************************************
 *
 *  getFieldNames:
 *
 *    Return field names
 *
 *******************************************************************/

vector <string> RocMeshData::getFieldNames() const
{
  return _fieldNames;
}


/*******************************************************************
 *
 *  getFieldUnits:
 *
 *    Return field units
 *
 *******************************************************************/

vector <string> RocMeshData::getFieldUnits() const
{
  return _fieldUnits;
}


/*******************************************************************
 *
 *  getFieldNumComp:
 *
 *    Return field numbers of components
 *
 *******************************************************************/

vector <int> RocMeshData::getFieldNumComp() const
{
  return _fieldComp;
}


/*******************************************************************
 *
 *  getFieldData:
 *
 *    Return field value
 *
 *******************************************************************/

vector <double *> RocMeshData::getFieldData() const
{
  return _fieldVals;
}


/*******************************************************************
 *                                                                 *
 *                                                                 *
 *                     HELPER FUNCTIONS                            *
 *                                                                 *
 *                                                                 *
 *******************************************************************/


/*******************************************************************
 *
 *  free:
 *
 *    Free all dynamic memory
 *
 *******************************************************************/

void RocMeshData::free()
{
  if (_elemData != NULL)   delete [] _elemData;
  if (_elemType != NULL)   delete [] _elemType;
  if (_pconnArray != NULL) delete [] _pconnArray;
  if (_coords != NULL)     delete [] _coords;

  _fieldNames.clear();
  _fieldUnits.clear();
  _fieldTypes.clear();
  _fieldComp.clear();
  
  for (int i = 0; i < _fieldVals.size(); i++)
    if (_fieldVals[i] != NULL) delete [] _fieldVals[i];
  _fieldVals.clear();  
}


/*******************************************************************
 *
 *  copy:
 *
 *    Copy all member variables
 *
 *******************************************************************/

void RocMeshData::copy(RocMeshData const & origVal)
{
  _filename     = origVal._filename;
  _indexBase    = origVal._indexBase;
  _numVer       = origVal._numVer;
  _numVerGhost  = origVal._numVerGhost;
  _numElem      = origVal._numElem;
  _elemDataSize = origVal._elemDataSize;
  _pconnSize    = origVal._pconnSize;

  _elemData = new int [_elemDataSize];
  for (int i = 0; i < _elemDataSize; i++)
    _elemData[i] = origVal._elemData[i];

  _elemType = new int [_numElem];
  for (int i = 0; i < _numElem; i++)
    _elemType[i] = origVal._elemType[i];

  _pconnArray = new int [_pconnSize];
  for (int i = 0; i < _pconnSize; i++)
    _pconnArray[i] = origVal._pconnArray[i];

  _coords = new double [(_numVer + _numVerGhost) * 3];
  for (int i = 0; i < 3*(_numVer + _numVerGhost); i++)
    _coords[i] = origVal._coords[i];

  _numFields    = origVal._numFields;
  _fieldNames   = origVal._fieldNames;
  _fieldUnits   = origVal._fieldUnits;
  _fieldTypes   = origVal._fieldTypes;
  _fieldComp    = origVal._fieldComp;
  
  for (int i = 0; i < origVal._fieldVals.size(); i++)
  {
    int size;
    if (_fieldTypes[i] == 'n')
      size = _numVer * _fieldComp[i];
    else
      size = _numElem * _fieldComp[i];

    double * newFieldVal = new double [size];
    for (int j = 0; j < size; j++)
      newFieldVal[j] = origVal._fieldVals[i][j];

    _fieldVals.push_back(newFieldVal);
  }
}


/*******************************************************************
 *
 *  initRocMeshData:
 *
 *    Initialize RocMeshData in constructors
 *
 *******************************************************************/

void RocMeshData::initRocMeshData(char const * filename)
{
  // helper variables
  const char * winName = "random";
  int paneId;          
  int numCon;            // number of connec. tables 
  char * conNames;       // names of connec. tables
  int * types;           // elem type for each connec. table
  int * numElemArray;    // num  elem for each connec. table 

  // load Rocin handle and load file to window

  loadData(filename, winName);

  // set all class variable member

  _filename = filename;
  _indexBase = 1;

  getPaneInfo(winName, paneId, numCon, conNames);

  getCoord(winName, paneId);

  getElemInfo(winName, paneId, numCon,         // _elemDataSize is 
	      conNames, types, numElemArray);  // set in this function

  allocElemMem(numCon, numElemArray, types);

  getElemConn(winName, paneId, numCon, 
	      conNames, types, numElemArray);

  getPconnArray(winName, paneId);

  getFieldInfo(winName, paneId);

  // free memory, buffer and window
  delete [] types;
  delete [] numElemArray;
  COM_free_buffer(&conNames);
  COM_delete_window(winName);
  COM_UNLOAD_MODULE_STATIC_DYNAMIC(Rocin, "IN");
}


/*******************************************************************
 *
 *  loadData:
 *
 *    Load data from Rocstar's HDF file to window
 *
 *******************************************************************/

void RocMeshData::loadData(char const * hdfFilename, const char * winName)
{
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocin, "IN");
  int IN_read = COM_get_function_handle("IN.read_windows");
  COM_set_verbose(0);
  COM_set_profiling(0);
  COM_call_function(IN_read, hdfFilename, winName, "");
}


/*******************************************************************
 *
 *  getPaneInfo:
 *
 *    Get pane id and connectivity table info
 *
 *******************************************************************/

void RocMeshData::getPaneInfo(const char * winName, int & paneId, 
			      int & numCon, char * & conNames)
{
  // Obtain the list of panes
  int numPane;
  int* paneIds;
  COM_get_panes(winName, &numPane, &paneIds);
  if (numPane > 1)
  {
    cerr << "\n\nError: currently not support multiple panes ...\n\n";
    exit(1);
  }

  // Obtain the connectivity tables
  COM_get_connectivities(winName, paneIds[0], &numCon, &conNames);

  // exit if structured mesh (table name starts with :st)
  if (numCon == 1 && strncmp(conNames,":st",3) == 0)
  {
    cerr << "\n\nError: currently not support structured mesh ...\n\n";
    exit(1);
  }

  // get pane id
  paneId = paneIds[0];

  // free buffer
  COM_free_buffer(&paneIds);
}


/*******************************************************************
 *
 *  getCoord:
 *
 *    Get coordinates
 *
 *******************************************************************/

void RocMeshData::getCoord(const char * winName, int paneId)
{
  // Obtain the size of nodes
  int  numNodes;                            // total number of nodes
  int  numGhostNodes;                       // Number of ghost nodes
  string w(winName); w+=".";
  COM_get_size((w+"nc").c_str(), paneId, &numNodes, &numGhostNodes);

  // get coordinates
  double *coor;
  int strd;
  COM_get_array_const((w+"nc").c_str(), paneId,
                      (const void **)&coor, &strd);

  _numVer      = numNodes - numGhostNodes;       // get only real nodes.
  _numVerGhost = numGhostNodes;
  _coords = new double [numNodes * 3];

  // put coordinates in coords
  int cjump =1, rjump =1;
  int m, n;
  if(strd == 1)
    cjump = numNodes;
  else
    rjump = strd;

  for (m = 0; m < numNodes; m++)
    for (n = 0; n < 3; n++)
      _coords[m * 3 + n] = coor[m * rjump + n * cjump];

  delete [] coor;
}


/*******************************************************************
 *
 *  getElemInfo:
 *
 *    Get total number of elements
 *
 *******************************************************************/

void RocMeshData::getElemInfo(const char * winName, int paneId, 
			      int numCon, char * conNames, 
			      int * & types, int * & numElemArray)
{
  // allocate memory

  types = new int [numCon];       
  numElemArray = new int [numCon];

  std::istringstream conNamesStream(conNames);
  for (int c=0; c < numCon; ++c)
  {
    string conName;
    conNamesStream >> conName;

    int numElems, numGhostElems;
    string w(winName); 
    w+=".";
    COM_get_size((w+conName).c_str(), paneId, &numElems, &numGhostElems);

    // assuming no mix between ghost and real elements
    if (numGhostElems == 0)
      numElemArray[c] = numElems;
    else
      numElemArray[c] = 0;

    string tableType = conName.substr(0,3);
    int simmetrixNumber = getSimElemTypeNum(tableType);

    if (simmetrixNumber == -1)
    {
      cerr << "\n\nError: unknown mesh type ...\n\n";
      exit(1);
    }
    else
      types[c] = simmetrixNumber;

    if (numGhostElems != 0)
      types[c] = -1;
  }
}


/*******************************************************************
 *
 *  allocElemMem:
 *
 *    Allocate memory for elemType and elemData
 *
 *******************************************************************/

void RocMeshData::allocElemMem(int numCon, int * numElemArray, int * types)
{
  // initialize
  _numElem = 0;                      
  _elemDataSize = 0;
  int sum = 0;

  // allocate _elemType
  for (int c = 0; c < numCon; ++c)     
    _numElem += numElemArray[c];
  _elemType = new int [_numElem];

  // allocate _elemData
  for (int c = 0; c < numCon; ++c)
  {
    for (int i = sum; i < sum + numElemArray[c]; i++)
      _elemType[i] = types[c];

    sum += numElemArray[c];
    int nodePerElem = getNodePerElem(types[c]);
    _elemDataSize += numElemArray[c] * nodePerElem;
  }
  _elemData = new int [_elemDataSize];
}

/*******************************************************************
 *
 *  getElemConn:
 *
 *    Get element connectivity
 *
 *******************************************************************/

void RocMeshData::getElemConn(const char * winName, int paneId, int numCon, 
			      char * conNames, int * types, int * numElemArray)
{
  std::istringstream conNamesStream(conNames);
  int counter = 0;
  string w(winName);
  w += ".";
  string conName;
  int numElems, numGhostElems;

  for ( int c = 0; c < numCon; ++c)
  {
    conNamesStream >> conName;
    COM_get_size((w+conName).c_str(), paneId, &numElems, &numGhostElems);

    int * connec;
    int stride;
    COM_get_array_const((w+conName).c_str(), paneId,
                        (const void **)&connec, &stride);

    int nrow =1; 
    int ncol =1;
    int row, col;
    if(stride == 1)
      nrow = numElems;
    else
      ncol = stride;

    int nodePerElem = getNodePerElem(types[c]);
    for (row = 0; row < numElemArray[c]; row++)
    {
      for (col = 0; col < nodePerElem; col++)
        _elemData[counter + row * nodePerElem + col] 
	  = connec[row * ncol + col * nrow];

      // Make mesh orientation positive for Simmetrix.
      // This assumes Tetrahedra in HDF have negative orientation.
      int temp = _elemData[counter + row * nodePerElem + 1];

      _elemData[counter + row * nodePerElem + 1]
	= _elemData[counter + row * nodePerElem];

      _elemData[counter + row * nodePerElem] = temp;
    }
    counter += numElemArray[c] * nodePerElem;

    COM_free_buffer(&connec);
  }
}


/******************************************************************
 *
 *  getPconnArray:
 *
 *    Get PCONN array and its size
 *
 *****************************************************************/

void RocMeshData::getPconnArray(const char * winName, int paneId)
{
  string w(winName);
  int numGhost;
  COM_get_size((w+".pconn").c_str(), paneId,
               &_pconnSize, &numGhost);

  COM_get_array_const((w+".pconn").c_str(), paneId,
                      (const void **)&_pconnArray);
}

/******************************************************************
 *
 *  getFieldInfo:
 *
 *    Get field information
 *
 *****************************************************************/

void RocMeshData::getFieldInfo(const char * winName, int paneId)
{

  // get field attributes
  getFieldAttributes(winName);

  // get field data
  getFieldData(winName, paneId);
}


/******************************************************************
 *
 *  getFieldAttributes:
 *
 *    Get field names, type, number of components and units.
 *    Allocate memory for field values.
 *
 *****************************************************************/

void RocMeshData::getFieldAttributes(const char * winName)
{
  // get raw field info
  int numAllFields;
  string allFieldNames;
  COM_get_dataitems(winName, &numAllFields, allFieldNames); 

  // populate field info
  std::istringstream fieldStream(allFieldNames);
  string fieldName;

  while (fieldStream >> fieldName)
  {
    // ignore ridges and grid speed
    if (fieldName == "ridges" || fieldName == "gs")
      continue;
    
    // process this field info except data
    string w(winName);
    w += ".";
    char loc;
    COM_Type comType;
    int ncomp;
    string unit;
    COM_get_dataitem((w + fieldName).c_str(), &loc, 
		      &comType, &ncomp, &unit);

    _fieldNames.push_back(fieldName);        
    _fieldUnits.push_back(unit);
    _fieldTypes.push_back(loc);
    _fieldComp.push_back(ncomp);

    // allocate data to _FieldVals
    int numItem;             // number of items (numVer/numElem)   
    if (loc == 'n')
      numItem = _numVer;
    else
      numItem = _numElem;
    
    double * fieldVal = new double [ncomp * numItem];
    _fieldVals.push_back(fieldVal);
  }
  _numFields = _fieldNames.size();
}


/******************************************************************
 *
 *  getFieldData:
 *
 *    Get field data from Roccom
 *
 *****************************************************************/

void RocMeshData::getFieldData(const char * winName, int paneId)
{
  string w(winName);         // window name string
  w += ".";

  for (int i = 0; i < _numFields; i++)
  {
    int numItem;             // number of items (numVer/numElem)
    if (_fieldTypes[i] == 'n')
      numItem = _numVer;
    else
      numItem = _numElem;

    // obtain field data from Roccom
    int stride;
    double * tempFld;        // pointer to data in Roccom
    if (_fieldComp[i] == 1)
    {
      COM_get_array_const((w + _fieldNames[i]).c_str(), paneId,
			  (const void **)&tempFld, &stride);

      for (int j = 0; j < numItem; j++)
	_fieldVals[i][j] = tempFld[j];
    }
    else                    // more than 1 component
    {
      for (int j = 0; j < _fieldComp[i]; j++)
      {
	std::stringstream ss;
	ss << (j + 1);       // index start from 1
	string prefix;
	prefix = ss.str() + "-";
	COM_get_array_const((w + prefix + _fieldNames[i]).c_str(), 
			    paneId, (const void **)&tempFld, &stride);

	// stride here = number of components
	for (int k = 0; k < numItem; k++)
	  _fieldVals[i][_fieldComp[i] * k + j] = tempFld[k]; 
      }
    }
  }
}


/*******************************************************************
 *
 *  getSimElemTypeNum:
 *
 *    Number of element type in Simmetrix format
 *
 *******************************************************************/

int RocMeshData::getSimElemTypeNum(string type)
{
  int simType;
  if (type == ":b2") simType = 0;
  else if (type == ":b3") simType = 1;
  else if (type == ":t3") simType = 5;
  else if (type == ":t4") simType = 6;
  else if (type == ":t6") simType = 7;
  else if (type == ":t8") simType = 8;
  else if (type == ":T4") simType = 10;
  else if (type == ":P5") simType = 11;
  else if (type == ":W6") simType = 12;
  else simType = -1;
  return simType;
}


/*******************************************************************
 *
 *  getNodePerElem:
 *
 *    get number of nodes per element given element type 
 *    (in Simmetrix format)
 *
 *******************************************************************/

int RocMeshData::getNodePerElem(int num) const
{
  int nodePerElem;
  switch (num)
  {
  case -1: nodePerElem = 0; break;
  case 0: nodePerElem = 2; break;
  case 1: nodePerElem = 3; break;
  case 5: nodePerElem = 3; break;
  case 6: nodePerElem = 4; break;
  case 7: nodePerElem = 6; break;
  case 8: nodePerElem = 8; break;
  case 10: nodePerElem = 4; break;
  case 11: nodePerElem = 5; break;
  case 12: nodePerElem = 6; break;
  }
  return nodePerElem;
}



