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
// $Id:

#include "Roccom_base.h"
#include "roccom.h"
#include "roccom_devel.h"


#include <iostream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cassert>
#include <map>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <cgnslib.h>
#include "Rocout.h"

using namespace std;

#define CG_CHECK(routine, args) \
{ \
  int ier = routine args; \
  if (ier != 0) { \
    cerr << "ghostbuster: " #routine " (line " << __LINE__ << " in " \
         << __FILE__ << ") failed: " << cg_get_error() << std::endl; \
    return -1; \
  } \
}

class AutoCloser
{
  public:
    inline AutoCloser(int fn)
    : m_fn(fn) {}
    inline ~AutoCloser() { cg_close(m_fn); }
  private:
    int m_fn;
};

COM_EXTERN_MODULE( Rocin);
COM_EXTERN_MODULE( Rocout);

int main(int argc, char *argv[]) {
  if ( argc < 2) {
    cout << "Usage: " << argv[0] << " <inputfile>" << endl;
    return -1;
  }

  // Determine true material names.  There's no clean ROCCOM way to do this.
  vector<string> wins;
  string materials;
  string file_in(argv[1]);
  string::size_type c = file_in.length();
  if (file_in.substr(c - 5) == ".cgns") {
    int fn;
    CG_CHECK(cg_open, (file_in.c_str(), MODE_READ, &fn));
    AutoCloser auto1(fn);

    int nBases;
    CG_CHECK(cg_nbases, (fn, &nBases));

    int B, cellDim, physDim;
    char baseName[33];
    for (B=1; B<=nBases; ++B) {
      CG_CHECK(cg_base_read, (fn, B, baseName, &cellDim, &physDim));
      std::string name(baseName);
      c = name.length();
      if (c <= 7 || name.substr(c - 7) != "_ridges") {
        wins.push_back(name);
        materials += name + ' ';
      }
    }
    if (!materials.empty())
      materials.erase(materials.size() - 1);
  } else {
    cout << "Error: " << argv[0] << " only handles CGNS files at this time."
         << endl;
    return -1;
  }

  std::cout << "Found materials: " << materials << std::endl;

  string file_out(file_in);
  c = file_out.rfind('/');
  if (c != string::npos)
    file_out.replace(0, c + 1, "gb_");
  else
    file_out.insert(0, "gb_");

  COM_init( &argc, &argv);

  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocin, "IN");
  COM_LOAD_MODULE_STATIC_DYNAMIC(Rocout, "OUT");

  int IN_read = COM_get_function_handle( "IN.read_windows");
  int OUT_set = COM_get_function_handle( "OUT.set_option");
  int OUT_write = COM_get_function_handle( "OUT.write_attribute");

  COM_call_function( OUT_set, "format", "CGNS");
  COM_call_function( OUT_set, "ghosthandle", "ignore");

  char time_level[33] = "";
  int length = 32;
  vector<string>::iterator w;
  for (w=wins.begin(); w!=wins.end(); ++w) {
    COM_call_function( IN_read, file_in.c_str(), wins[0].c_str(), NULL, NULL,
                       NULL, time_level, &length);

    int IN_all = COM_get_attribute_handle((*w+".all").c_str());

    COM_call_function( OUT_write, file_out.c_str(), &IN_all, (*w).c_str(),
                       time_level);

    COM_call_function( OUT_set, "mode", "a");
  }

  COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocin, "IN");
  COM_UNLOAD_MODULE_STATIC_DYNAMIC( Rocout, "OUT");

  COM_finalize();

  return 0;
}






