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
#ifndef _TRAIL_FLU_H_
#define _TRAIL_FLU_H_

bool TRAIL_FluInitSolver(GEM_Partition &gp,const std::string &prefix = "./");
bool TRAIL_FluReadControlFile(GEM_Partition &gp,
			     const std::string &prefix = "./");
bool TRAIL_FluWriteCOM(GEM_Partition &gp,const std::string &path);
bool TRAIL_FluReadCOM(GEM_Partition &gp);
bool TRAIL_FluWriteMAP(GEM_Partition &gp,unsigned int nproc,
		      unsigned int nregions,const std::string &path);
bool TRAIL_FluWriteDIM(GEM_Partition &gp,double t=0.0,bool unsteady=true,
		      const std::string &path="./");

bool TRAIL_FluReadDIM(GEM_Partition &gp,double t=0.0,bool unsteady=true);

bool TRAIL_FluPopRemBordIndFILE(GEM_Partition &gp,double t=0.0,
			       bool unsteady=true,
			       const std::string &path="./");
bool TRAIL_FluWriteCMP(GEM_Partition &gp,const std::string &path);
void TRAIL_FluResizeVolSoln(GEM_Partition &gp);
bool TRAIL_FluInitVolSoln(GEM_Partition &gp);
void TRAIL_FluResizeSurfSoln(GEM_DomainBoundary &db);
bool TRAIL_FluInitSurfSoln(GEM_DomainBoundary &db,const std::string &,
			  bool all=true);
bool TRAIL_FluInitSurfSoln(GEM_Partition &gp);
int  TRAIL_FluNumPatches(GEM_Partition &gp);
bool TRAIL_FluRegisterVolSoln(GEM_Partition &gp,bool all=true);
bool TRAIL_FluRegisterSurfSoln(GEM_DomainBoundary &db,const std::string &wname,
			      bool all=true);
bool TRAIL_FluRegisterSurfSoln(GEM_Partition &gp,bool all=true);
bool TRAIL_FluRegisterSurfMesh(GEM_Partition &gp);
bool TRAIL_FluBuildPatchMapping(std::map<unsigned int,unsigned int> &pmap,
			       const std::string &prefix);
bool TRAIL_FluWriteNative(GEM_Partition &gp,const std::string &path);
bool TRAIL_FluPopulatePatches(GEM_Partition &gp);
bool TRAIL_FluReadGridASCII(GEM_Partition &gp,double t=0.0,bool unsteady=true);
void TRAIL_FluCopyCaseFiles(GEM_Partition &gp,const std::string &path);


#endif


















