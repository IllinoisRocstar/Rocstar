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
/** \file Rocout.h
 *  Rocout creates a series of Roccom windows by reading in a list of files.
 *  Rocout can also copy Roccom attributes from window to window.
 *  Only HDF4 files are supported in the current implementation.
 */
/*  Author John Norris
 *  Initial date:   May 17, 2004
 */


#ifndef _ROCOUT_H_
#define _ROCOUT_H_

#include <map>
#include "roccom.h"
#include "HDF4.h"
#include "roccom_devel.h"

#ifdef USE_CGNS
#include "Rocout_cgns.h"
#define MODE_READ CG_MODE_READ 
#define MODE_WRITE CG_MODE_WRITE 
#define MODE_MODIFY CG_MODE_MODIFY
#endif // USE_CGNS

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
 ** \param ghosthandle "ignore" or "write" on ghost data.
 ** \param errorhandle "ignore", "warn", or "abort" on errors.
 ** \param mode Write == 0, append == 1. (Input)
 **/
void write_attr_CGNS(const std::string& fname, const std::string& mfile, 
                     const COM::Attribute* attr, const char* material, 
                     const char* timelevel, int pane_id,
                     const std::string& ghosthandle,
                     const std::string& errorhandle, int mode);

/** \name Module loading and unloading
 * \{
 */
/// Load the module Rocout into Roccom using the given module name. 
/// This module provides one subroutine: "write_attribute".
extern "C" void Rocout_load_module( const char *name);
/// Unload the module Rocout from Roccom.
extern "C" void Rocout_unload_module( const char *name);
//\}


class Rocout : public COM_Object {
 public:
  /** \name User interface
   *  \{
   */
  /** Write a (possibly aggregate) attribute to a file.
   *
   * \param filename_pre the prefix of the file name.
   * \param attr a reference to the attribute to be written.
   * \param material the name of the material (usually the window name).
   * \param time_level the time stamp of the dataset.
   * \param mfile_pre the prefix of the file that contains the mesh data.
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of the Roccom window will be used.
   * \param pane_id the pane to be written.
   */
  void write_attribute( const char* filename_pre,
                        const COM::Attribute* attr,
                        const char* material,
                        const char* timelevel,
                        const char* mfile_pre=NULL,
                        const MPI_Comm* comm=NULL,
                        const int* pane_id=NULL);

  /** Write a (possibly aggregate) attribute to a new file.
   *
   * \param filename_pre the prefix of the file name.
   * \param attr a reference to the attribute to be written.
   * \param material the name of the material (usually the window name).
   * \param time_level the time stamp of the dataset.
   * \param mfile_pre the prefix of the file that contains the mesh data.
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of the Roccom window will be used.
   * \param pane_id the pane to be written.
   */
  void put_attribute( const char* filename_pre,
                      const COM::Attribute* attr,
                      const char* material,
                      const char* timelevel,
                      const char* mfile_pre=NULL,
                      const MPI_Comm* comm=NULL,
                      const int* pane_id=NULL);

  /** Append a (possibly aggregate) attribute to a file.
   *
   * \param filename_pre the prefix of the file name.
   * \param attr a reference to the attribute to be written.
   * \param material the name of the material (usually the window name).
   * \param time_level the time stamp of the dataset.
   * \param mfile_pre the prefix of the file that contains the mesh data.
   * \param comm The MPI communicator to use. If is NULL, the default 
   *        communicator of the Roccom window will be used.
   * \param pane_id the pane to be written.
   */
  void add_attribute( const char* filename_pre,
                      const COM::Attribute* attr,
                      const char* material,
                      const char* timelevel,
                      const char* mfile_pre=NULL,
                      const MPI_Comm* comm=NULL,
                      const int* pane_id=NULL);

  /** Wait for the completion of an asychronous write operation.
   */
  void sync();

  /** Generate a control file for Rocin.
   *
   * \param window_name The name of the Roccom window.
   * \param file_prfixes The prefixes of the data files.
   * \param control_file_name The name of the file to be written.
   */
  void write_rocin_control_file( const char* window_name,
                                 const char* file_prefixes,
                                 const char* control_file_name);

  /** Set an option for Rocout, such as controlling the output format.
   *
   * \param option_name the option name: "format", "async", "mode",
   *        "localdir", "rankwidth", "pnidwidth", "separator" or "errorhandle".
   * \param option_val the option value.
   */
  void set_option( const char* option_name,
                   const char* option_val);

  /** Write out the parameters defined in the fiven window into
   * a parameter file. Only process 0 of the communicator writes
   * the parameter file.
   *
   * \param file_name the file to which the parameters will be
   * written.
   * \param window_name the window to read parameters from
   * \param comm the MPI_communicator to use. If NULL, then the
   * communicator of the window associated with window_name
   * is used.
   */
  void write_parameter_file( const char* file_name, 
			     const char* window_name, 
			     const MPI_Comm *comm = NULL);

  /** Set options for Rocout via a control file.
   *
   * \param filename the path to the control file.
   */
  void read_control_file( const char* filename);
    //\}

protected:  
  /** \name Initialization and finalization
   * \{
   */
  /** Initialize the module by registering it to Roccom with 
   *  the given module name. This function is called 
   *  Rocout_load_module.
   *  \see Rocout_load_module()
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
  friend void Rocout_load_module( const char *name);
  friend void Rocout_unload_module( const char *name);
  //\}

  /** \name Implementation
   * \{
   */
  /** Does the actual writing to file.
   *
   * \param attrInfo Information on what to write and where to write it.
   */
  static
  void* write_attr_internal(void* attrInfo);

  /** Builds a filename from the given prefix and rank.
   *
   * \param pre Filename prefix.
   * \param rank The rank of this MPI process.
   * \param pPaneId A pointer to the pane id.
   * \param check Check for an existing file.
   */
  std::string get_fname(const std::string& pre, int rank = -1,
                        const int paneId = 0, bool check = false);
  //\}

  std::map<std::string, std::string> _options;
# ifdef USE_PTHREADS
  std::vector<pthread_t> _writers;
  static Semaphore _writesem;
# endif // USE_PTHREADS
};

#endif






