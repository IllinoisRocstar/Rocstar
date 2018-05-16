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
/** \file Rocout.C
 *  Implementation of Rocout for creation of a data file from a Roccom window.
 */
/*  Author John Norris
 *  Initial date:   May 17, 2004
 */

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <errno.h>

#include "Rocout.h"
#include "Rocout_hdf4.h"
#ifdef USE_CGNS
#include "Rocout_cgns.h"
#endif // USE_CGNS

#ifndef DOXYGEN_SHOULD_SKIP_THIS
USE_COM_NAME_SPACE
#endif

static const int MAX_ASYNC_WRITES = 10;

//! Pass write_attribute arguments to the background worker thread.
struct WriteAttrInfo {
  WriteAttrInfo( Rocout *rout, const char* filename_pre, const Attribute* attr,
                 const char* material, const char* timelevel,
                 const char* mfile_pre = NULL, const MPI_Comm* pComm = NULL,
                 const int* pane_id = NULL, int append = -1,
                 bool cloned = false)
  : m_rout(rout), m_prefix(filename_pre), m_attr(attr), m_material(material),
    m_timelevel(timelevel), m_meshPrefix(mfile_pre != NULL ? mfile_pre : ""),
    m_pComm(pComm), m_pPaneId(pane_id), m_append(append), m_cloned(cloned) {}

  Rocout *m_rout;
  const std::string m_prefix;
  const Attribute* m_attr;
  const std::string m_material;
  const std::string m_timelevel;
  const std::string m_meshPrefix;
  const MPI_Comm* m_pComm;
  const int* m_pPaneId;
  const int m_append;
  const bool m_cloned;
};

#define SwitchOnDataType(dType, funcCall) \
   switch (dType) { \
      case COM_CHAR: \
      case COM_BYTE: \
         { typedef char TT; \
           funcCall; } \
         break; \
      case COM_UNSIGNED_CHAR: \
         { typedef unsigned char TT; \
           funcCall; } \
         break; \
      case COM_SHORT: \
         { typedef short TT; \
           funcCall; } \
         break; \
      case COM_UNSIGNED_SHORT: \
         { typedef unsigned short TT; \
           funcCall; } \
         break; \
      case COM_INT: \
         { typedef int TT; \
           funcCall; } \
         break; \
      case COM_UNSIGNED: \
         { typedef unsigned int TT; \
           funcCall; } \
         break; \
      case COM_LONG: \
         { typedef long TT; \
           funcCall; } \
         break; \
      case COM_UNSIGNED_LONG: \
         { typedef unsigned long TT; \
           funcCall; } \
         break; \
      case COM_FLOAT: \
         { typedef float TT; \
           funcCall; } \
         break; \
      case COM_DOUBLE: \
         { typedef double TT; \
           funcCall; } \
         break; \
      case COM_LONG_DOUBLE: \
         { typedef long double TT; \
           funcCall; } \
         break; \
   }

#define ERROR_MSG(msg) \
{ if (_options["errorhandle"] != "ignore") { \
    std::cerr << msg << std::endl; \
    if (_options["errorhandle"] == "abort") { \
      if (COMMPI_Initialized()) \
        MPI_Abort(MPI_COMM_WORLD, 0); \
      else \
        abort(); \
    } \
  } \
}

#ifdef USE_PTHREADS
Semaphore Rocout::_writesem(MAX_ASYNC_WRITES, MAX_ASYNC_WRITES);
#endif // USE_PTHREADS

void Rocout::init(const std::string &mname) {

  HDF4::init();

  Rocout *rout = new Rocout();

  // Masoud, changing default to CGNS
  //rout->_options["format"] = "HDF4";
  rout->_options["format"] = "CGNS";
  //
  rout->_options["async"] = "off";
  rout->_options["mode"] = "w";
  rout->_options["localdir"] = "";
  rout->_options["rankwidth"] = "4";
  rout->_options["pnidwidth"] = "0";
  rout->_options["separator"] = "_";
  rout->_options["errorhandle"] = "abort";
  rout->_options["rankdir"] = "off";
  rout->_options["ghosthandle"] = "write";

  COM_new_window( mname.c_str(), MPI_COMM_SELF);

  std::string glb=mname+".global";

  COM_new_attribute( glb.c_str(), 'w', COM_VOID, 1, "");
  COM_set_object( glb.c_str(), 0, rout);


  // Register the function write_attribute
  COM_Type types[8] = { COM_RAWDATA, COM_STRING, COM_METADATA, COM_STRING, 
			COM_STRING, COM_STRING, COM_MPI_COMM, COM_VOID};
  COM_set_member_function( (mname+".write_attribute").c_str(),
			   (Member_func_ptr)&Rocout::write_attribute, 
			   glb.c_str(), "biiiiIII", types);
  COM_set_member_function( (mname+".put_attribute").c_str(),
			   (Member_func_ptr)&Rocout::put_attribute, 
			   glb.c_str(), "biiiiIII", types);
  COM_set_member_function( (mname+".add_attribute").c_str(),
			   (Member_func_ptr)&Rocout::add_attribute, 
			   glb.c_str(), "biiiiIII", types);
  COM_set_member_function( (mname+".write_window").c_str(),
			   (Member_func_ptr)&Rocout::write_attribute,
			   glb.c_str(), "biiiiIII", types);

  // Register the function write_rocin_control_file
  types[2] = types[3] = COM_STRING;
  COM_set_member_function( (mname+".write_rocin_control_file").c_str(),
			   (Member_func_ptr)&Rocout::write_rocin_control_file,
			   glb.c_str(), "biii", types);

  // Register the function sync
  COM_set_member_function( (mname+".sync").c_str(),
			   (Member_func_ptr)&Rocout::sync, 
			   glb.c_str(), "b", types);

  // Register the function set_option
  COM_set_member_function( (mname+".set_option").c_str(), 
			   (Member_func_ptr)&Rocout::set_option, 
			   glb.c_str(), "bii", types);

  // Register the function write_parameter_file
  types[3] = COM_MPI_COMM;
  COM_set_member_function( (mname+".write_parameter_file").c_str(),
			   (Member_func_ptr)&Rocout::write_parameter_file,
			   glb.c_str(), "biiI", types);

  // Register the function read_control_file
  COM_set_member_function( (mname+".read_control_file").c_str(), 
			   (Member_func_ptr)&Rocout::read_control_file, 
			   glb.c_str(), "bi", types);

  COM_window_init_done( mname.c_str());
}

void Rocout::finalize(const std::string &mname) {
  Rocout *rout;
  std::string glb=mname+".global";

  COM_get_object( glb.c_str(), 0, &rout);
  
  COM_delete_window( mname.c_str());

#ifdef USE_PTHREADS
  // Wait for any writer threads to finish.
  void* ret;
  std::vector<pthread_t>::iterator p = rout->_writers.begin();
  while (p != rout->_writers.end()) {
    pthread_join(*p, &ret);
    ++p;
  }
# endif // USE_PTHREADS

  delete rout;
  HDF4::finalize();
}

//! Write an attribute to file.
/*!
 * \param filename_pre the prefix of the file name.
 * \param attr a reference to the attribute to be written.
 * \param material the name of the material (usually the window name).
 * \param time_level the time stamp of the dataset.
 * \param mfile_pre the prefix of the file that contains the mesh data.
 * \param pComm The MPI communicator to use. If is NULL, the default 
 *        communicator of the Roccom window will be used.
 * \param pane_id the pane to be written.
 */
void Rocout::write_attribute( const char* filename_pre,
                              const Attribute* attr, const char* material,
                              const char* timelevel, const char* mfile_pre,
                              const MPI_Comm* pComm, const int* pane_id)
{
  WriteAttrInfo* pWAI;

#ifdef USE_PTHREADS
  if (_options["async"] == "off") {
#endif // USE_PTHREADS
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id);
    write_attr_internal(pWAI);
#ifdef USE_PTHREADS
  } else {
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id, -1, true);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    pthread_t id;
    pthread_create(&id, &attr, write_attr_internal, pWAI);
    pthread_attr_destroy(&attr);
    _writers.push_back(id);
  }
#endif // USE_PTHREADS
}

//! Write an attribute to a new file.
/*!
 * \param filename_pre the prefix of the file name.
 * \param attr a reference to the attribute to be written.
 * \param material the name of the material (usually the window name).
 * \param time_level the time stamp of the dataset.
 * \param mfile_pre the prefix of the file that contains the mesh data.
 * \param pComm The MPI communicator to use. If is NULL, the default 
 *        communicator of the Roccom window will be used.
 * \param pane_id the pane to be written.
 */
void Rocout::put_attribute( const char* filename_pre,
                            const Attribute* attr, const char* material,
                            const char* timelevel, const char* mfile_pre,
                            const MPI_Comm* pComm, const int* pane_id)
{
  //std::cout << __FILE__ << __LINE__ << std::endl;
  //std::cout << "Attribute = " << (attr)->fullname() << std::endl;
  //std::cout << "Size = " << (attr)->size_of_items() << std::endl;
  WriteAttrInfo* pWAI;

#ifdef USE_PTHREADS
  if (_options["async"] == "off") {
#endif // USE_PTHREADS
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id, 0);
    write_attr_internal(pWAI);
#ifdef USE_PTHREADS
  } else {
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id, 0, true);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    pthread_t id;
    pthread_create(&id, &attr, write_attr_internal, pWAI);
    pthread_attr_destroy(&attr);
    _writers.push_back(id);
  }
#endif // USE_PTHREADS
}

//! Append an attribute to a file.
/*!
 * \param filename_pre the prefix of the file name.
 * \param attr a reference to the attribute to be written.
 * \param material the name of the material (usually the window name).
 * \param time_level the time stamp of the dataset.
 * \param mfile_pre the prefix of the file that contains the mesh data.
 * \param pComm The MPI communicator to use. If is NULL, the default 
 *        communicator of the Roccom window will be used.
 * \param pane_id the pane to be written.
 */
void Rocout::add_attribute( const char* filename_pre,
                            const Attribute* attr, const char* material,
                            const char* timelevel, const char* mfile_pre,
                            const MPI_Comm* pComm, const int* pane_id)
{
  //std::cout << __FILE__ << __LINE__ << std::endl;
  //std::cout << "Attribute = " << (attr)->fullname() << std::endl;
  //std::cout << "Size = " << (attr)->size_of_items() << std::endl;
  WriteAttrInfo* pWAI;

#ifdef USE_PTHREADS
  if (_options["async"] == "off") {
#endif // USE_PTHREADS
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id, 1);
    write_attr_internal(pWAI);
#ifdef USE_PTHREADS
  } else {
    pWAI = new WriteAttrInfo(this, filename_pre, attr, material, timelevel, 
			     mfile_pre, pComm, pane_id, 1, true);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    pthread_t id;
    pthread_create(&id, &attr, write_attr_internal, pWAI);
    pthread_attr_destroy(&attr);
    _writers.push_back(id);
  }
#endif // USE_PTHREADS
}

void Rocout::write_rocin_control_file(const char* window_name,
                                      const char* file_prefixes,
                                      const char* control_file_name)
{
  const MPI_Comm default_comm = COM_get_default_communicator();
  const MPI_Comm comm_null = MPI_COMM_NULL;
  const MPI_Comm* myComm = COMMPI_Initialized() ? &default_comm : &comm_null;


  // Obtain process rank
  int flag = 0, rank = 0, size = 1;
  MPI_Initialized(&flag);

  if (flag) {
    MPI_Comm_rank(*myComm, &rank);
    MPI_Comm_size(*myComm, &size);
  }
 
  if (rank == 0) {
    std::vector<std::vector<int> > paneIds(size);

    /* Get list of panes local each process */
    int i;
    for (i=0; i<size; ++i)
      COM_get_panes(window_name, paneIds[i], i);
 
    int rw, pw;
    {
      std::istringstream sin(_options["rankwidth"]);
      sin >> rw;
    }
    {
      std::istringstream sin(_options["pnidwidth"]);
      sin >> pw;
    }

/*
    std::ostringstream sout;
    sout << "@Files:";

    // Write out file name pattern only if number of panes is nonzero.
    std::istringstream sin(file_prefixes);
    while (!sin.eof()) {
      std::string prefix;
      sin >> prefix;
      if (prefix.empty())
        continue;

      sout << ' ' << prefix;
      if (rw > 0)
        sout << "%0" << rw << 'p';
      if (pw > 0) {
        if (rw > 0)
	  sout << _options["separator"];
        sout << "%0" << pw << 'i';
      }

      const std::string fmt = _options["format"];
      if ( fmt == "HDF4" || fmt == "HDF")
        sout << ".hdf";
      else if (fmt == "HDF5")
        sout << ".hdf5";
      else if (fmt == "CGNS")
        sout << ".cgns";
    }
*/

    /* Now write data into control file */
    std::ofstream fout(control_file_name);
    COM_assertion_msg( fout.is_open(),
		       (std::string("Rocout cannot open control file:")+control_file_name+" for writing.\n").c_str());
    
    const std::string fmt = _options["format"];
    for (i=0; i<size; i++) {
      fout << "@Proc: " << i << std::endl;
      if ( !paneIds[i].empty()) {
        std::ostringstream sout;
        sout << "@Files:";

        // Write out file name pattern only if number of panes is nonzero.
        std::istringstream sin(file_prefixes);
        while (!sin.eof()) {
          std::string prefix;
          sin >> prefix;
          if (prefix.empty())
            continue;

          sout << ' ';
            // write output file in <rank> dir
          if (_options["rankdir"] == "on") {
            std::ostringstream rank_prefix; 
            rank_prefix << i << "/";
            sout << rank_prefix.str();
          }
          sout << prefix;
          if (rw > 0)
            sout << "%0" << rw << 'p';
          if (pw > 0) {
            if (rw > 0)
	      sout << _options["separator"];
            sout << "%0" << pw << 'i';
          }

          if ( fmt == "HDF4" || fmt == "HDF")
            sout << ".hdf";
          else if (fmt == "HDF5")
            sout << ".hdf5";
          else if (fmt == "CGNS")
            sout << ".cgns";
        }
	fout << sout.str() << std::endl;
      }
      else // Write out an empty block
	fout << "@Files: " << std::endl;

      fout << "@Panes:";
      
      /* Write all the panes of this process to control file */
      std::vector<int>::iterator p;
      for (p=paneIds[i].begin(); p!=paneIds[i].end(); ++p)
	fout << ' ' << *p;
      fout << '\n' << std::endl;
    }

    fout.close();
  }
}

/** Wait for the completion of an asychronous write operation.
 */
void Rocout::sync()
{
#ifdef USE_PTHREADS
  void* ret;
  std::vector<pthread_t>::iterator p = _writers.begin();
  while (p != _writers.end()) {
    pthread_join(*p, &ret);
    ++p;
  }
  _writers.clear();
#endif // USE_PTHREADS
}

/** Return true if the given string is the name of a Rocout option.
 */
static bool is_option_name(const std::string& name)
{
  return (name == "format" || name == "async" || name == "mode"
          || name == "localdir" || name == "rankwidth" || name == "pnidwidth"
          || name == "separator" || name == "errorhandle" || name == "rankdir"
          || name == "ghosthandle");
}

// Return true if the given string is a whole number.
static bool is_whole(const std::string& s)
{
  if (s.empty())
    return false;

  std::istringstream sin(s);
  int i = 0;
  sin >> i;
  return (i >= 0 && sin.eof());
}

/** Return true if the given Rocout option name/value pair is valid.
 */
static bool is_option_value(const std::string& name, const std::string& val)
{
  return ((name == "format"
           && (val == "HDF" || val == "HDF4" || val == "HDF5" || val == "CGNS"))
          || (name == "async"
              && (val == "on" || val == "off"))
          || (name == "mode"
              && (val == "w" || val == "a"))
          || (name == "localdir" /* && is_valid_path(val) */ )
          || ((name == "rankwidth" || name == "pnidwidth") && is_whole(val))
          || (name == "rankdir" 
              && (val == "on" || val == "off"))
          || (name == "errorhandle"
              && (val == "abort" || val == "ignore" || val == "warn"))
          || (name == "ghosthandle"
              && (val == "write" || val == "ignore")));
}

/** Set an option for Rocout, such as controlling the output format.
 *
 * \param option_name the option name: "format", "async", "mode", "localdir",
 *        "rankdir", "rankwidth", "pnidwidth", "errorhandle" or "ghosthandle".
 * \param option_val the option value.
 */
void Rocout::set_option( const char* option_name, const char* option_val)
{
  const std::string name(option_name);
  const std::string val(option_val);

  if (!is_option_name(name)) {
    ERROR_MSG("Rocout::set_option(): unknown option name \"" << name << "\".");
    return;
  }

  if (!is_option_value(name, val)) {
    ERROR_MSG("Rocout::set_option(): invalid value \"" << val
              << "\" for option \"" << name << "\".");
    return;
  }

#ifndef USE_CGNS
  COM_assertion_msg(name != "format" || val != "CGNS",
                    "Roccom not built with option CGNS=1.\n");
#endif // USE_CGNS

  _options[name] = val;
}

/** Set options for Rocout via a control file.
 *
 * \param filename the path to the control file.
 */
void Rocout::read_control_file( const char* filename)
{
  // Attempt to open the control file.
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    ERROR_MSG("Rocout::read_control_file(): could not open file \""
              << filename << "\".");
    return;
  }

  // Each line of the control file should be of the form "option_name=value",
  // e.g. "format=HDF4".  White space at the beginning, end, and around the
  // '=' is ignored, as are empty lines and comments (lines that start with
  // '#').
  int lineNum = 0;
  char buffer[256];
  std::string line, name, val;
  std::string::size_type eq;
  while (!fin.eof()) {
    ++lineNum;
    // Read a line.
    fin.getline(buffer, sizeof(buffer));
    line = buffer;

    // Skip empty lines and comments.
    if (line.empty() || line.find_first_not_of(" \t\n") == std::string::npos
        || line[line.find_first_not_of(" \t\n")] == '#')
      continue;

    // Find the '='.  It's mandatory.
    eq = line.find('=');
    if (eq == std::string::npos) {
      ERROR_MSG("Rocout::read_control_file(): option without value at line "
                << lineNum << " of " << filename << '.');
      continue;
    }

    // Split the line into an option name and a value.
    name = line.substr(0, eq);
    val = line.substr(eq + 1);

    // Remove leading whitespace from the option name.
    if (name.find_first_not_of(" \t\n") != 0)
      name.erase(0, name.find_first_not_of(" \t\n"));

    // Check for an empty option name.
    if (name.empty()) {
      ERROR_MSG("Rocout::read_control_file(): missing option name at line "
                << lineNum << " of " << filename << '.');
      continue;
    }

    // Remove trailing whitespace from the option name.
    if (name.find_last_not_of(" \t\n") != name.length() - 1)
      name.erase(name.find_last_not_of(" \t\n") + 1);

    // Make sure it is a valid option name.
    if (!is_option_name(name)) {
      ERROR_MSG("Rocout::read_control_file(): unknown option name \"" << name
                << "\" at line " << lineNum << " of " << filename << '.');
      continue;
    }

    // Remove leading whitespace from the option value.
    if (val.find_first_not_of(" \t\n") != 0)
      val.erase(0, val.find_first_not_of(" \t\n"));

    // Check for an empty option value.
    if (val.empty()) {
      ERROR_MSG("Rocout::read_control_file(): option without value at line " 
                << lineNum << " of " << filename << '.');
      continue;
    }

    // Remove trailing whitespace from the option value.
    if (val.find_last_not_of(" \t\n") != val.length() - 1)
      val.erase(val.find_last_not_of(" \t\n") + 1);

    // Make sure it is a valid value for the named option.
    if (!is_option_value(name, val)) {
      ERROR_MSG("Rocout::read_control_file(): invalid value \"" << val
                << "\" for option \"" << name << "\" at line " << lineNum
                << " of " << filename << '.');
      continue;
    }

    set_option(name.c_str(), val.c_str());
  }
}

void* Rocout::write_attr_internal(void* attrInfo)
{
  WriteAttrInfo* ai = static_cast<WriteAttrInfo*>(attrInfo);
  const Attribute* attr = ai->m_attr;
  Window* tempWin = NULL;

  if (ai->m_cloned) {
#ifdef USE_PTHREADS
    _writesem.Wait();
#endif // USE_PTHREADS
    tempWin = new Window(attr->window()->name(),
                         attr->window()->get_communicator());
    attr = tempWin->inherit(const_cast<Attribute*>(attr), attr->name(),
                            Pane::INHERIT_CLONE, true, NULL, 0);
  }

  int flag = 0; MPI_Initialized(&flag);

  // Obtain process rank
  int rank;
  if ( flag) {
    MPI_Comm wcomm;
    const MPI_Comm* pComm;

    if ( ai->m_pComm) 
      pComm = ai->m_pComm;
    else { 
      wcomm = attr->window()->get_communicator();
      pComm = &wcomm;
    }

    if (*pComm != MPI_COMM_NULL)
      MPI_Comm_rank(*pComm, &rank); 
    else
      rank = 0;
  }
  else 
    rank = 0;

  int append = ai->m_append;
  //std::cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
  //std::cout << __FILE__ << __LINE__
  //          << " append = " << append << std::endl;
  if (append < 0) {
    if (ai->m_rout->_options["mode"] == "w")
      append = 0;
    else
      append = 1;
  }

  std::vector<int> paneIds;
  COM_get_panes(attr->window()->name().c_str(), paneIds);

  std::vector<int>::iterator begin = paneIds.begin(), end = paneIds.end(), p;
  if (ai->m_pPaneId != NULL && *(ai->m_pPaneId) > 0) {
    begin = std::find(begin, end, *(ai->m_pPaneId));
    if (begin != end)
      end = begin + 1;
  }

  // MS
  /*
  std::cout << __FILE__ << __LINE__ << std::endl;
  std::cout << "Attribute requested = " << attr->fullname() << std::endl;
  std::cout << "Location = " << attr->location() << std::endl;
  std::cout << "Size of items = " << attr->size_of_items() << std::endl;
  std::cout << "Size of components = " << attr->size_of_components() << std::endl;
  std::cout << "Window information :" << std::endl;
  std::string wname = attr->window()->name();
  int nNde; 
  for (p=begin; p!=end; ++p) {
    std::cout << " pane = " << *p << std::endl;
    COM_get_size((wname+".nc").c_str(), *p, &nNde);
    std::cout << " number of nodes = " << nNde << std::endl;
    std::string stringNames;
    int numConn;
    COM_get_connectivities(wname.c_str(), *p, &numConn, stringNames);
    std::istringstream ConnISS(stringNames);
    std::vector<std::string> connNames;
    std::cout << "\t # \t Type \t #Elm \t #ElmNde \t Loc \n";
    std::cout << "\t---\t------\t------\t---------\t-----\n";
    for (int i=0; i<numConn; ++i) {
       std::string name;
       ConnISS >> name;
       connNames.push_back(name);
       std::string fullConnName(wname+"."+name);
       int nElm;
       COM_get_size(fullConnName, *p, &nElm);
       int nElmNde;
       std::string dataItemUnits;
       char dataItemLoc;
       COM_Type dataItemType;
       COM_get_attribute(fullConnName, &dataItemLoc, &dataItemType,
		      &nElmNde, &dataItemUnits);
       std::cout << "\t " << i+1
  		 << "\t " << name
		 << "\t  " << nElm
		 << "\t    " << nElmNde
		 << "\t          " << dataItemLoc
		 << std::endl;
    }
  }
  std::cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
  */
  // MS End

  std::set<std::string> written;
  for (p=begin; p!=end; ++p) {
    //std::cout << __FILE__ << __LINE__
    //          << " pane = " << *p << std::endl;
    std::string fname, mfile;
    fname = ai->m_rout->get_fname(ai->m_prefix, rank, *p, true);
    //std::cout << " fname = " << fname << std::endl;
    if (!ai->m_meshPrefix.empty())
      mfile = ai->m_rout->get_fname(ai->m_meshPrefix, rank, *p);

    int ap = append + written.count(fname);
    written.insert(fname);

    const std::string fmt = ai->m_rout->_options["format"];
    if ( fmt == "HDF4" || fmt == "HDF") {

      write_attr_HDF4(fname, mfile, attr, ai->m_material.c_str(),
                      ai->m_timelevel.c_str(), *p,
                      ai->m_rout->_options["errorhandle"], ap);

    } else if ( fmt == "CGNS") {
#ifdef USE_CGNS
     write_attr_CGNS(fname, mfile, attr, ai->m_material.c_str(),
                     ai->m_timelevel.c_str(), *p,
                     ai->m_rout->_options["ghosthandle"],
                     ai->m_rout->_options["errorhandle"], ap);
#endif // USE_CGNS

    }
  }

  if (ai->m_cloned) {
#ifdef USE_PTHREADS
    _writesem.Post();
#endif // USE_PTHREADS
    delete tempWin;
  }

  delete ai;

  return NULL;
}

/** Build a filename.
 *
 * Get a file name by appending an underscore, a 4-digit rank id,
 * and an extension to the "localdir" and given prefix.
 * If the pre contains .hdf or .cgns, then use it as the file name.
 */
std::string Rocout::get_fname(const std::string& prefix,
                              int rank /* = -1 */,
                              int paneId /* = 0 */,
                              bool check /* = false */)
{
  // Modify the prefix using the "localdir" option.
  std::string pre(_options["localdir"]);
  if (!pre.empty()) {
    // Make sure there's exactly one '/' between the localdir and given prefix.
    if (prefix[0] == '/') {
      if (pre[pre.length()-1] == '/')
        pre.erase(pre.length() - 1);
    } else {
      if (pre[pre.length()-1] != '/')
        pre += '/';
    }
  }
  pre += prefix;

  if (_options["rankdir"] == "on") {   // write output file in <rank> dir
    std::ostringstream rank_prefix; 
    rank_prefix << "/" << rank;
    std::string::size_type s = pre.find_last_of('/');
    if (s != std::string::npos)
      pre.insert(s, rank_prefix.str());
  }

  // Make sure the directory exists.  Ignore any errors except for the last.
  int result = 0;
  std::string::size_type s = pre.find('/', 1);
  while (s != std::string::npos) {
    result = mkdir(pre.substr(0, s).c_str(), 0755);
    s = pre.find('/', s + 1);
  }
  if (result < 0 && errno != EEXIST) {
    ERROR_MSG("Rocout::write_attribute(): could not create directory '" 
              << pre.substr(0, pre.rfind('/') + 1) << "'.");
  }

  if ( pre.find(".hdf") == pre.size()-4) {
    _options["format"] = "HDF4";
    return pre;
  }
  else if (pre.find(".cgns") == pre.size()-5) {
    _options["format"] = "CGNS";
    return pre;
  }
  
  if (rank < 0) {
    int flag = 0;
    MPI_Initialized(&flag);
    if (flag)
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    else
      rank = 0; 
  }             
                
  std::string name;
    
  if (!pre.empty()) {
    int rw, pw;
    {
      std::istringstream sin(_options["rankwidth"]);
      sin >> rw;
    }
    {
      std::istringstream sin(_options["pnidwidth"]);
      sin >> pw;
    }
    
    std::ostringstream sout;
    sout << pre;
    if (rw > 0)
      sout << std::setw(rw) << std::setfill('0') << rank;
    if (pw > 0 && paneId > 0) {
      if (rw > 0)
        sout << _options["separator"];
      sout << std::setw(pw) << std::setfill('0') << paneId;
    }
  
    const std::string fmt = _options["format"];
    if ( fmt == "HDF" || fmt == "HDF4")
      sout << ".hdf";
    else if (fmt == "HDF5")
      sout << ".hdf5";
    else if (fmt == "CGNS")
      sout << ".cgns";
    
    name = sout.str();
    
    if (check) {
      std::FILE* f;
      if ((f = std::fopen(name.c_str(), "r")) != NULL) 
        std::fclose(f);
      else if ((f = std::fopen(name.c_str(), "w")) != NULL) {
        std::fclose(f);
        std::remove(name.c_str());
      } else {
        // Change the directory to current directory
        std::string::size_type n = name.rfind('/');
        if (n != std::string::npos)
          name.erase(0, n + 1);
      }
    }
  }

  return name;
}

extern "C" void Rocout_load_module( const char *name)
{ Rocout::init( std::string(name)); }

extern "C" void Rocout_unload_module( const char *name) 
{ Rocout::finalize( std::string(name)); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// All possible Fortran bindings

extern "C" void COM_F_FUNC2(rocout_load_module, ROCOUT_LOAD_MODULE)( const char *name, long int length)
{ Rocout_load_module( std::string(name, length).c_str()); }
extern "C" void COM_F_FUNC2(rocout_unload_module, ROCOUT_UNLOAD_MODULE)( const char *name, long int length) 
{ Rocout_unload_module( std::string(name, length).c_str()); }

#endif // DOXYGEN_SHOULD_SKIP_THIS







