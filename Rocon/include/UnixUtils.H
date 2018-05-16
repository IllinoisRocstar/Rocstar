///
/// @file
/// @ingroup irad_group
/// @brief Unix System Tools interface
///
#ifndef _UNIX_UTIL_H_
#define _UNIX_UTIL_H_

#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <istream>
#include <ostream>
#include <streambuf>
#include <cstdio>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <errno.h>

namespace IRAD {

  /// Unix system utilities
  namespace Sys {
    std::string LogTime();
    void TokenizePath(std::vector<std::string> rv,const std::string &path);
    std::string TempFileName(const std::string &stub);
    int OpenTemp(std::string &stub);
    const std::string Hostname(bool longname = false);
    const std::string StripDirs(const std::string &);
    const std::string CWD();
    int SymLink(const std::string &source,const std::string &target);
    void SafeRemove(const std::string &fname,const std::string &ext);
    int ChDir(const std::string &path);
    bool FILEEXISTS(const std::string &fname);
    bool ISDIR(const std::string &fname);
    bool ISLINK(const std::string &fname);
    int CreateDirectory(const std::string &fname);
    const std::string ResolveLink(const std::string &path);
    int Remove(const std::string &fname);
    int Rename(const std::string &source_file,const std::string &target_file);
    class Directory : public std::vector<std::string>
    {
    private:
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


    class Environment : public std::vector< std::pair<std::string,std::string> >
    {
    public:
      Environment();
      int SetEnv(const std::string &,const std::string &,bool);
      void UnSetEnv(const std::string &);
#ifndef DARWIN
      int ClearEnv();
#endif
      const std::string GetEnv(const std::string &) const;
      std::string &GetEnv(const std::string &);
      int PutEnv(char *);
      void Refresh();
      char **GetRawEnv();
    private:
      void init();
      std::string empty_string;
    };
  }; 
};
#endif
