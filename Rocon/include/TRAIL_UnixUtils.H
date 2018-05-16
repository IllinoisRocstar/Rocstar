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
#ifndef _UNIX_UTIL_H_
#define _UNIX_UTIL_H_

#include <ostream>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>


std::string TRAIL_CWD(void);
void TRAIL_SafeRemove(const std::string &fname,const std::string &ext);
int TRAIL_CD(const std::string &path,std::ostream * = NULL);
bool TRAIL_FILEEXISTS(const std::string &fname);
bool TRAIL_ISDIR(const std::string &fname);
bool TRAIL_ISLINK(const std::string &fname);
int TRAIL_CreateDirectory(const std::string &fname);
std::string TRAIL_ResolveLink(const std::string &path);
std::string ResolveLink(const std::string &path);
class Directory : public std::vector<std::string>
{
protected:
  std::string _path;
  bool   _good;
  DIR    *_dir;
public:
  Directory(const std::string &s = "");
  ~Directory();
  int open(const std::string &s = "");
  void close();
  operator void* ();
  bool operator ! ();
};

#endif






