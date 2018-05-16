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
#include <string>
#include <vector>
#include <cstdio>
#include <sstream>

#include "TRAIL_UnixUtils.H"

void
TRAIL_SafeRemove(const std::string &fname,const std::string &ext)
{
  if(!TRAIL_FILEEXISTS(fname))
    return;
  if(TRAIL_ISLINK(fname))
    unlink(fname.c_str());
  std::string savename_base(fname+"."+ext);
  std::string savename(savename_base);
  unsigned int n = 1;
  while(TRAIL_FILEEXISTS(savename)){
    std::ostringstream Ostr;
    Ostr << savename_base << n++;
    savename.assign(Ostr.str());
  }
  rename(fname.c_str(),savename.c_str());
}

bool
TRAIL_FILEEXISTS(const std::string &fname)
{
  struct stat fstat;
  if(lstat(fname.c_str(),&fstat))
    return false;
  return true;
}

bool
TRAIL_ISDIR(const std::string &fname)
{
  struct stat fstat;
  if(stat(fname.c_str(),&fstat))
    return false;
  if(S_ISDIR(fstat.st_mode))
    return true;
  return false;
}

bool
TRAIL_ISLINK(const std::string &fname)
{
  struct stat fstat;
  if(lstat(fname.c_str(),&fstat))
    return false;
  if(S_ISLNK(fstat.st_mode))
    return true;
  return(false);
}

int
TRAIL_CreateDirectory(const std::string &fname)
{
  return(mkdir(fname.c_str(),S_IRGRP | S_IXGRP  | S_IRWXU));  
}

std::string
TRAIL_ResolveLink(const std::string &path)
{
  std::string retVal;
  char buf[1024];
  size_t s = readlink(path.c_str(),buf,1024);
  if(!(s <= 0)){
    buf[s] = '\0';
    retVal.assign(buf);
  }
  std::string::size_type x = retVal.find_last_of("/");
  if(x != std::string::npos)
    retVal.erase(x);
  return (retVal);
}

std::string
ResolveLink(const std::string &path)
{
  std::string retVal;
  char buf[1024];
  size_t s = readlink(path.c_str(),buf,1024);
  if(!(s <= 0)){
    buf[s] = '\0';
    retVal.assign(buf);
  }
  std::string::size_type x = retVal.find_last_of("/");
  if(x != std::string::npos)
    retVal.erase(x);
  return (retVal);
}

Directory::Directory(const std::string &path)
{
  _good = false;
  _dir = NULL;
  _path.assign(path);
  if(open(path))
    std::cerr << "Directory::Error: Could not open " << path 
	      << " as a directory." << std::endl;
}

Directory::~Directory()
{
  if (_good)
    closedir(_dir);
}

Directory::operator void* ()
{
  return(_good ? this : NULL);
}

bool 
Directory::operator ! ()
{
  return(!_good);
}

void 
Directory::close()
{
  if(_good)
    closedir(_dir);
}

int 
Directory::open(const std::string &path)
{
  if(_good){
    this->close();
    _path = path;
  }
  if(path.empty())
    return(1);
  if(!(_dir = opendir(path.c_str())))
    return(1);
  _path = path;
  _good = true;
  struct dirent *entry;
  // Skip . and ..
  entry = readdir(_dir);
  entry = readdir(_dir);
  while((entry = readdir(_dir)) != NULL)
    this->push_back(entry->d_name);
  return(0);
}

std::string TRAIL_CWD(void)
{
  char buf[1024];
  return(std::string(getcwd(buf,1024)));
}

int
TRAIL_CD(const std::string &path,std::ostream *ouf)
{
  if(ouf)
    *ouf << "TRAIL_CD: Switching directory from " 
	 << TRAIL_CWD() << " to " << path << std::endl;
  int result = chdir(path.c_str());
  if(ouf)
    *ouf << "TRAIL_CD: CWD(" << TRAIL_CWD() << ")" << std::endl;
  return(result);

}






