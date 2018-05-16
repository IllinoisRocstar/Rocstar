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
#ifndef _TRAIL_REMESH_H_
#define _TRAIL_REMESH_H_
bool TRAIL_RemeshInitFluSurfData(GEM_Partition &gp);
bool TRAIL_RemeshInitFluSurfData(GEM_DomainBoundary &db,int src_index=0);
bool TRAIL_RemeshInitFluVolData(GEM_Partition &gp,int src_index=0);
bool TRAIL_RemeshWrite(GEM_Partition &gp,const std::string &path,double t=0.0,
		      MPI_Comm comm = MPI_COMM_NULL,bool=false);
void TRAIL_RemeshRunDirSetup(const std::string &path,double t,MPI_Comm comm,
			    bool shuffle = false);
void TRAIL_RemeshShutdown(MPI_Comm comm);
void TRAIL_RemeshShuffle(const std::string &solver,double t=0,bool debug=false);
bool TRAIL_RemeshWrite(GEM_Partition &gp,
		      const std::string &solver,
		      const std::string &source_directory,
		      const std::string &target_directory,
		      double t,MPI_Comm comm,
		      bool surftrans,bool restartshuffle);
#endif






