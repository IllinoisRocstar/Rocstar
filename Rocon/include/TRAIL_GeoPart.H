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
#ifndef _TRAIL_GEOPART_
#define _TRAIL_GEOPART_

bool TRAIL_GeoPartReadPatchASCII(GEM_DomainBoundary &db,std::ifstream &Inf);
bool TRAIL_GeoPartReadPatchesASCII(GEM_Partition &gp,std::ifstream &Inf);
bool TRAIL_GeoPartReadPartitionBoundaryASCII(GEM_Partition &gp,
					    unsigned int pb_index,
					    std::ifstream &Inf);
bool TRAIL_GeoPartReadPartitionBoundariesASCII(GEM_Partition &gp,
					      std::ifstream &Inf);
bool TRAIL_GeoPartReadASCII(GEM_Partition &gp,std::ifstream &Inf);

#endif














