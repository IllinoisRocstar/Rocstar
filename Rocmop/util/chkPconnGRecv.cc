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
/********************************************************************************
 *                                                                              *
 *   chkPconnGRecv.cc                                                           *
 *                                                                              *
 *   Take HDF files and checks PCONN ghost receive (block 3) if complete.       *
 *   Write out missing ghost nodes to Tecplot if requested.                     *
 *                                                                              *
 *   Pornput Suriyamongkol                                                      *
 *   4/26/07                                                                    *
 *                                                                              *
 *******************************************************************************/


#include "RocMeshData.h"
#include <sstream>
#include <vector>
#include <unistd.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ostringstream;
using std::vector;

// function prototype
void writeUsage();


int main(int argc, char** argv)
{
  char* sampleFName;
  int numPars;
  int parBase;
  int option;
  vector <string> errFiles;

  // check argument
  if(argc != 5) 
    writeUsage();
  else
  {
    sampleFName = argv[1];
    sscanf(argv[2], "%d", &numPars);
    sscanf(argv[3], "%d", &parBase);
    sscanf(argv[4], "%d", &option);
  }

  // set up
  COM_init(&argc, &argv);
  int numDigits = 0;
  int locLast_ = 0;     // location of last under score
  string baseFName(sampleFName);
  baseFName = baseFName.substr(0, baseFName.size() - 4);
  locLast_ = baseFName.rfind("_", baseFName.size());
  numDigits = baseFName.size() - locLast_ - 1;
  baseFName = baseFName.substr(0, baseFName.size() - numDigits);

  // check each mesh
  for (int i = 0; i < numPars; i++)
  {
    // get filename and print it to screen
    string filename;
    int index;
    string indexStr;
    index = i + parBase;
    ostringstream indexStream;
    indexStream << index;
    indexStr = indexStream.str();
    int digit2pad = numDigits - indexStr.size();
    for (int j = 0; j < digit2pad; j++)
      indexStr.insert(0, "0");
    
    filename = baseFName + indexStr + ".hdf";
    cout << "\n\n---------------------------------------------------\n";
    cout << "Reading HDF file:   " << filename << endl;
    cout << "---------------------------------------------------\n\n";

    // check pconn block 3
    RocMeshData rocMesh(filename.c_str());

    cout << "... PCONN Checking ...\n\n";
    int pconnErr = 0;

    filename = filename.substr(0, filename.size() - 4);
    if (option == 1)
      pconnErr = rocMesh.pconnCheck((filename + "_missingGNodes.plt").c_str());
    else
      pconnErr = rocMesh.pconnCheck();

    if (pconnErr)
    {
      cout << "... PCONN ERROR ...\n\n";
      errFiles.push_back(filename);
    }
    else
      cout << "... PCONN OK ...\n\n";
  }

  // Write out Summary
  cout << "==============================\n";
  cout << "     chkPconnGRecv SUMMARY\n";
  cout << "==============================\n\n";
  if (errFiles.size() == 0)
    cout << "NO ERROR.\n";
  else
  {
    cout << "There are " << errFiles.size()
	 << " files with incomplete PCONN.\n";
    for (int i = 0; i < errFiles.size(); i++)
    {
      cout << errFiles[i] << ".hdf: ";
      if (option != 1)
	cout << "No file is written.\n";
      else
	cout << "Missing ghost nodes are written out in "
	     << errFiles[i] << "_missingGNodes.plt\n";
    }
  }
  cout << endl;

  // Return
  COM_finalize();
  return 0;
}


void writeUsage()
{
    cout << "\n\n-----------------------------------"
	 << "------------------------------------\n";
    cout << "Usage: chkPconnGRecv <sample_filename> "
	 << "<num_pars> <par_base> <option>\n";
    cout << "------------------------------------" 
	 << "-----------------------------------\n";
    cout << "<sammple_filename>: Sample Filename\n";
    cout << "<num_pars>: Number of partitions\n";
    cout << "<par_base>: Partition number starts from 0 or 1\n";
    cout << "<option>: Set to 1 to write out missing "
	 << "ghost nodes to Tecplot\n\n";
    exit(1);
}

