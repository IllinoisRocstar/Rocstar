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
#include <iostream>
#include <vector>
#include <list>
#include <sstream>

using namespace std;

#include "clop.H"

extern "C" {
  void RFLUMP(const char *,int *,long int);
}

void
Usage(const string &pn)
{
  cout << endl << "Usage: " << pn << " -c <casename> [-v 0-2]" << endl << endl
       << "   -c | --casename     :     Specifies the casename." << endl
       << "   -v | --verbosity    :     0 - Quiet" << endl
       << "                             1 - Moderately verbose" << endl
       << "                             2 - Ridiculously verbose" << endl
       << endl;
}

int
main(int argc,char *argv[])
{
  // Get the commandline into a string vector - it's easier
  // to deal with that way.
  vector<string> args = Vectize((const char **)argv,argc);

  // Get the name of the executable by stripping off any leading
  // directory names
  string program_name(stripdir(args[0]));

  // Specify the allowable options to the program
  AddOp("casename",'c');
  AddOp("verbosity",'v');
  AddOp("help",'h');
  
  // Declare some variables for command line argument handling
  string casename;
  string sverb;
  int verbosity;
  bool help;
  bool isset;
  
  // See if the help option is specified, if so give'm the usage text
  if(help = GetOp("help",args)){
    Usage(program_name);
    exit(0);
  }

  // Process casename option, if it's not set then fail
  if(GetOp("casename",casename,args)){
    if(casename.empty()){ // casename was empty
      cerr << program_name 
	   << ": Expected casename after casename option."
	   << " Use -h for usage instructions."
	   << endl;
      exit(1);
    }
  }
  else{ // option not specified (but it's required!)
    cerr << program_name
	 << ": Missing required casename option."
	 << "  Use -h for usage instructions."
	 << endl;
    exit(1);
  }      

  // Process verbosity option
  if(GetOp("verbosity",sverb,args)){
    if(sverb.empty()){
      cerr << program_name
	   << ": Expected verbosity level. "
	   << "Use -h for usage instructions." << endl;
      exit(1);
    }
    istringstream Istr(sverb);
    Istr >> verbosity;
    if(verbosity < 0 || verbosity > 2){ // Some jerk specified a non numeric or negative
      cerr << program_name
	   << ": Invalid verbosity value.  Use -h for usage "
	   << "instructions." << endl;
      exit(1);
    }
  }
  else{ // Default verbosity
    verbosity = 1;
  }

  RFLUMP(casename.c_str(),&verbosity,casename.length());
}






