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
// $Id: Rocin.h,v 1.13 2008/12/06 08:43:19 mtcampbe Exp $

/** \file Rocin.h
 *  Rocin creates a series of Roccom windows by reading in a list of files.
 *  Rocin can also copy Roccom attributes from window to window.
 *  Only HDF4 files are supported in the current implementation.
 */
/*  Author John Norris
 *  Initial date:   March 19, 2004
 */
/*  Modified by Phillip Alexander
 *  April 8, 2004
 *  -updated to Roccom3 compliance
 *  -added support for use as a Roccom module
 */


#ifndef _ROCIN_H_
#define _ROCIN_H_

#include "roccom.h"
#include "HDF4.h"
#include "rocin_block.h"
#include "roccom_devel.h"

#ifdef USE_CGNS
#include "cgnslib.h"
#define MODE_READ CG_MODE_READ
#endif // USE_CGNS

/** \name Module loading and unloading
 * \{
 */
/// Load the module Rocin into Roccom using the given module name. 
/// This module provides three subroutines: "read_windows", 
/// "read_by_control_file", and "obtain_attribute".
extern "C" void Rocin_load_module( const char *name);
/// Unload the module Rocin from Roccom.
extern "C" void Rocin_unload_module( const char *name);
//\}

class Rocin : public COM_Object {
public:
  /// Default constructor
  Rocin() : m_is_local( NULL), m_base(0), m_offset(0) {}

  /// Pointer to a function to determine locality of a pane.
  typedef void (*RulesPtr)(const int &pane_id, const int &comm_rank,
                           const int &comm_size, int* is_local);

  /** \name User interface
   *  \{
   */
  /** Create a single Roccom window by reading in the files specified by
   *  the given control file.
   *
   * \param control_file_name specifies the control file.
   * \param window_name the name of the window to be created. 
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of Roccom will be used.
   * \param time_level the time stamp of the dataset to be read.  If time
   *        level is NULL (default) or "", then the first time_level in the
   *        file will be assumed. 
   * \param str_len if present and positive, the time stamp of the dataset 
   *        will be copied to the time_level string. 
   */
  void read_by_control_file( const char* control_file_name,
                             const char* window_name,
                             const MPI_Comm* comm=NULL,
                             char* time_level=NULL,
                             const int* str_len=NULL);

  /** Create a single Roccom window by reading in a list of files.
   *
   * \param filename_patterns specifies the patterns of the files to be
   *        read in.  Patterns take the form of regular expressions in 
   *        the same format as the shell commands use.  Separate multiple
   *        patterns with empty space.
   * \param window_name the name of the window to be created. 
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of Roccom will be used.
   * \param is_local a function pointer which determines wheter a pane
   *        should be read by a process.
   * \param time_level the time stamp of the dataset to be read.  If time
   *        level is NULL (default) or "", then the first time_level in the
   *        file will be assumed. 
   * \param str_len if present and positive, the time stamp of the dataset 
   *        will be copied to the time_level string. 
   */
  void read_window( const char* filename_patterns,
                    const char* window_name,
                    const MPI_Comm* comm=NULL,
                    RulesPtr is_local=NULL,
                    char* time_level=NULL,
                    const int *str_len=NULL);
  
  /** Create a series of Roccom windows by reading in a list of files.
   *
   * \param filename_patterns specifies the patterns of the files to be
   *        read in.  Patterns take the form of regular expressions in 
   *        the same format as the shell commands use.  Separate multiple
   *        patterns with empty space.
   * \param window_prefix a prefix of the name(s) of the window(s) to be
   *        created.  Actual window name is appended by a material name
   *        and there can be multiple materials per file.
   * \param material_names a list of space-separated materials to be read
   *        from files.  If material name is NULL (default), then "" will 
   *        be used and it is assumed that files contain only one type of 
   *        material.
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of Roccom will be used.
   * \param is_local a function pointer which determines wheter a pane should
   *        be read by a process.
   * \param time_level the time stamp of the dataset to be read.  If time
   *        level is NULL (default) or "", then the first time_level in the
   *        file will be assumed. 
   * \param str_len if present and positive, the time stamp of the dataset 
   *        will be copied to the time_level string. 
   */
  void read_windows(const char* filename_patterns,
                    const char* window_prefix,
                    const char* material_names=NULL,
                    const MPI_Comm* comm=NULL,
                    RulesPtr is_local=NULL,
                    char* time_level=NULL,
                    const int *str_len=NULL);
  
  /** Fill the destination (second) attribute from files using the data 
   *  corresponding to the source (first) attribute. Note that the 
   *  destination and the source attributes can be the same. 
   *
   * \param attribute_in   source attribute
   * \param user_attribute destination attribute
   * \param pane_id if present, then copy only the specific pane
   *
   * \see obtain_attribute()
   */
  
  void obtain_attribute( const COM::Attribute *attribute_in,
                         COM::Attribute *user_attribute,
                         int *pane_id=NULL);

  /** Read in parameters from a given file into the given window. 
   *  If the window exists, then only read in the values for the attributes 
   *  which exist in the window. If not, then create the window
   *  and define all options as window attributes of base type of character 
   *  strings. Process 0 of the communicator should read in the parameters, 
   *  and then broadcast to all the processes. If comm is not specified, then 
   *  use the communicator of the window as default. If an option is listed
   *  more than once in the parameter_file, then the last value for that
   *  option will overwrite the others.
   *
   * \param file_name the file containing the parameters to be read in.
   * \param window_name window where parameters will be placed
   * \param comm The MPI communicator for parameter broadcast. If NULL,
   *        then the communicator of window_name will be used.
   */

void read_parameter_file( const char* file_name, 
                          const char* window_name, 
                          const MPI_Comm *comm = NULL);

  //\}

protected:
  /** \name Pane distribution
   *  Assign panes onto different processes.
   */
  typedef void (Rocin::*MemberRulePtr)(const int&, const int&,
                                        const int&, int*);

  void explicit_local( const int& pid, const int& comm_rank,
                       const int& comm_size, int* il);

  void cyclic_local(const int& pid, const int& comm_rank,
                    const int& comm_size, int* il);
  
  void blockcyclic_local(const int& pid, const int& comm_rank,
                         const int& comm_size, int* il);

  void register_panes(BlockMM_HDF4::iterator hdf4,
                      const BlockMM_HDF4::iterator& hdf4End,
#ifdef USE_CGNS
                      BlockMM_CGNS::iterator cgns,
                      const BlockMM_CGNS::iterator& cgnsEnd,
#endif // USE_CGNS
                      const std::string& window, RulesPtr is_local,
                      const MPI_Comm* comm, int rank, int nprocs);

  //\}

protected:  
  /** \name Initialization and finalization
   * \{
   */
  /** Initialize the module by registering it to Roccom with 
   *  the given module name. This function is called 
   *  Rocin_load_module.
   *  \see Rocin_load_module()
   */
  static 
  void init(const std::string &mname);
  
  /** Finalize the module by deregistering it from Roccom.
   *  \see Rochdf_unload_module()
   */
  static 
  void finalize(const std::string &mname);
  //\}

  /** \name Module loading and unloading
   * \{
   */
  friend void Rocin_load_module( const char *name);
  friend void Rocin_unload_module( const char *name);
  //\}

protected:
  MemberRulePtr m_is_local;
  std::set<int> m_pane_ids;
  int m_base;
  int m_offset;

  std::map<int32, COM_Type> m_HDF2COM;
#ifdef USE_CGNS
  std::map<DataType_t, COM_Type> m_CGNS2COM;
#endif // USE_CGNS
};

#endif






