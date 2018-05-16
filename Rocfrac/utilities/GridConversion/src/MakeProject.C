/// 
/// @file
/// @ingroup ir_group
/// @brief Project name changing utility
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>

#include "UnixUtils.H"
#include "FDUtils.H"

namespace GridConversion 
{
  ///
  /// Exclude certain template files from conversion in project creation.
  ///
  bool Excluded(const std::string &filename)
  {
    std::string::size_type x = filename.find("IR.dox");
    if(x != std::string::npos)
      return(true);
    x = filename.find("MakeProject.dox");
    if(x != std::string::npos)
      return(true);
    x = filename.find("ExistingProject.dox");
    if(x != std::string::npos)
      return(true);
    return(false);
  }

  
  ///
  /// Creates a new project from an GridConversion base
  ///
  /// This program is designed to make a new IllinoisRocstar
  /// software project starting from an intallation of the
  /// GridConversion project template.
  ///
  /// Usage of this program should go like the following:
  /// 1. Get and Install the GridConversion template.
  /// <blockquote>
  ///  svn co $IRREPO/GridConversion/trunk GridConversion\n
  ///  mkdir gridconversion_build\n
  ///  cd gridconversion_build\n
  ///  cmake ../GridConversion\n
  ///  make\n
  ///  cd ../\n
  /// </blockquote>
  /// 2. Create the new project, <NewProject>,  with the "make_project" command.
  /// <blockquote>
  ///  gridconversion_build/bin/make_project GridConversion <NewProject>\n
  /// </blockquote>
  ///
  /// @note Be sure to replace <NewProject> with your new project name.
  ///
  /// The usage for the make_project command is: 
  /// <blockquote>
  ///       make_project usage:\n
  /// 
  ///       make_project \<template name\> \<new project\> [verb level]\n
  ///
  ///       This program will read the project template from the \<template name\>\n 
  ///       directory and create a new "blank" IllinoisRocstar project\n
  ///       named \<new name\>, in a directory named \<new name\>.\n
  ///       An optional verblevel of 1 or 2 can be given to make the process\n
  ///       more verbose.\n
  /// </blockquote>
  ///
  int MakeProject(int argc,char *argv[])
  {
    if(argc < 3){
      std::cout << "Usage:" << std::endl << std::endl
                << argv[0] << " <template name> <new name> [verb level]" << std::endl
                << std::endl
                << "This program will read the project template from the <template name> "
                << std::endl
                << "directory and create a new \"blank\" IllinoisRocstar project"
                << std::endl
                << "named <new name>, in a directory named <new name>."
                << std::endl
                << "An optional verblevel of 1 or 2 can be given to make the process"
                << std::endl << "more verbose." << std::endl;
      return(0);
    }
    std::string OriginalName(argv[1]);
    std::string NewName(argv[2]);
    int verb = 0;
    if(argv[3]){
      verb = 1;
      int v = atoi(argv[3]);
      if(v > 0)
        verb = v;
    }
    std::vector<std::string> ProtectedFiles;
    ProtectedFiles.push_back("AUTHORS");
    ProtectedFiles.push_back("CTestConfig.cmake");
    ProtectedFiles.push_back("LICENSE");
    ProtectedFiles.push_back(".svn");
    std::vector<std::string> CommonFiles;
    CommonFiles.push_back("CMakeLists.txt");
    int syserr = 0;
    std::string olower(OriginalName);
    std::string oupper(OriginalName);
    std::string nlower(NewName);
    std::string nupper(NewName);
    std::transform(olower.begin(),olower.end(),olower.begin(),tolower);
    std::transform(oupper.begin(),oupper.end(),oupper.begin(),toupper);
    std::transform(nlower.begin(),nlower.end(),nlower.begin(),tolower);
    std::transform(nupper.begin(),nupper.end(),nupper.begin(),toupper);
  
    if(verb)
      std::cout << "Creating a new project (" << NewName 
                << ") from project template (" << OriginalName
                << ")." << std::endl << std::endl
                << "Creating top level project directories...";
    std::string dirname(NewName);
    if(!IRAD::Sys::FILEEXISTS(NewName)){
      if(verb > 1)
        std::cout << "   Creating directory " << dirname << "...";
      syserr = IRAD::Sys::CreateDirectory(dirname);
      if(syserr){
        std::cout << "Unable to create directory " << dirname << "."
                  << std::endl;
        return(1);
      }
      if(verb > 1)
        std::cout << "done." << std::endl;
    }
    dirname += "/branches";
    if(!IRAD::Sys::FILEEXISTS(dirname)){
      if(verb > 1)
        std::cout << "   Creating directory " << dirname << "...";
      syserr = IRAD::Sys::CreateDirectory(dirname);
      if(syserr){
        std::cout << "Unable to create directory " << dirname << "."
                  << std::endl;
        return(1);
      }
      if(verb > 1)
        std::cout << "done." << std::endl;
    }
    dirname = NewName+"/tags";
    if(!IRAD::Sys::FILEEXISTS(dirname)){
      if(verb > 1)
        std::cout << "   Creating directory " << dirname << "...";
      syserr = IRAD::Sys::CreateDirectory(dirname);
      if(syserr){
        std::cout << "Unable to create directory " << dirname << "."
                  << std::endl;
        return(1);
      }
      if(verb > 1)
        std::cout << "done." << std::endl;
    }
    dirname.assign(NewName+"/examples");
    if(!IRAD::Sys::FILEEXISTS(dirname)){
      if(verb > 1)
        std::cout << "   Creating directory " << dirname << "...";
      syserr = IRAD::Sys::CreateDirectory(dirname);
      if(syserr){
        std::cout << "Unable to create directory " << dirname << "."
                  << std::endl;
        return(1);
      }
      if(verb > 1)
        std::cout << "done." << std::endl;
    }
    bool protect_svn = false;
    dirname = NewName+"/trunk";
    std::vector<std::string>::iterator pfi = ProtectedFiles.begin();
    while(pfi != ProtectedFiles.end()){
      std::string ProtectThisFile(dirname+"/"+*pfi++);
      std::string ProtectedFile(ProtectThisFile+".backup");
      if(IRAD::Sys::FILEEXISTS(ProtectThisFile))
        IRAD::Sys::Rename(ProtectThisFile,ProtectedFile);
    }
    pfi = CommonFiles.begin();
    while(pfi != CommonFiles.end()){
      std::string ProtectThisFile(dirname+"/"+*pfi++);
      std::string ProtectedFile(ProtectThisFile+".backup");
      if(IRAD::Sys::FILEEXISTS(ProtectThisFile))
        IRAD::Sys::Rename(ProtectThisFile,ProtectedFile);
    }
    std::ostringstream ComStr;
    if(!IRAD::Sys::FILEEXISTS(dirname)){
      if(verb > 1)
        std::cout << "   Creating directory " << dirname << "...";
      ComStr << "cp -r " << OriginalName << " " << dirname;
    } else {
      if(verb > 1)
        std::cout <<   "   Making project files from template ...";
      ComStr << "cp -r " << OriginalName << "/* " << dirname;
    }
    IRAD::Sys::InProcess System(ComStr.str());
    std::string comline;
    while(std::getline(System,comline)){
      if(verb > 1)
        std::cout << comline << std::endl;
    }
    if(verb)
      std::cout << "done." << std::endl;
    if(verb)
      std::cout << "Cleaning up ...";
    ComStr.str("");
    ComStr << "rm -rf " << dirname << "/.svn";
    System.Execute(ComStr.str());
    int n = 0;
    while(std::getline(System,comline))
      n++;
    pfi = ProtectedFiles.begin();
    while(pfi != ProtectedFiles.end()){
      std::string ProtectThisFile(dirname+"/"+*pfi++);
      std::string ProtectedFile(ProtectThisFile+".backup");
      if(IRAD::Sys::FILEEXISTS(ProtectedFile)){
        if(IRAD::Sys::FILEEXISTS(ProtectThisFile))
          IRAD::Sys::Remove(ProtectThisFile);
        IRAD::Sys::Rename(ProtectedFile,ProtectThisFile);
      }
    }
    pfi = CommonFiles.begin();
    while(pfi != CommonFiles.end()){
      std::string ProtectThisFile(dirname+"/"+*pfi++);
      std::string ProtectedFile(ProtectThisFile+".backup");
      std::string CommonFileTemplate(ProtectThisFile+".template");
      if(IRAD::Sys::FILEEXISTS(ProtectedFile)){
        if(IRAD::Sys::FILEEXISTS(ProtectThisFile))
          IRAD::Sys::Rename(ProtectThisFile,CommonFileTemplate);
        IRAD::Sys::Rename(ProtectedFile,ProtectThisFile);
      }
    }
    if(verb)
      std::cout << "done." << std::endl;
    if(verb > 1)
      std::cout << "Done creating new project files." 
                << std::endl;
    if(verb)
      std::cout << "Renaming project...";
    if(IRAD::Sys::ChDir(dirname)){
      std::cout << "Something went wrong, cannot find new project directory." 
                << std::endl;
      return(1);
    }
    ComStr.str("");
    ComStr << "grep -i " << OriginalName << " -r * | cut -d \":\" -f 1 | sort | uniq";
    if(verb > 1)
      std::cout << "   " << ComStr.str() << std::endl;
    System.Execute(ComStr.str());
    std::vector<std::string> filenames;
    if(verb > 1)
      std::cout << "   Files to change:" << std::endl;
    while(std::getline(System,comline)){
      if(!Excluded(comline)){
        filenames.push_back(comline);
        if(verb > 1)
          std::cout << "     " << comline << std::endl;
      }
    }
    std::vector<std::string>::iterator fni = filenames.begin();
    if(verb > 1)
      std::cout << "   Processing files...."; 
    while(fni != filenames.end()){
      std::string filename(*fni++);
      if(verb > 1)
        std::cout << "     File: " << filename << std::endl;
      ComStr.str("");
      ComStr << "sed -i 's/" << OriginalName << "/" << NewName << "/g' " << filename;
      if(verb > 1)
        std::cout << "       " << ComStr.str() << std::endl;
      System.Execute(ComStr.str());
      int n = 0;
      while(std::getline(System,comline))
        n++;
      ComStr.str("");
      ComStr << "sed -i 's/" << oupper << "/" << nupper << "/g' " << filename;
      if(verb > 1)
        std::cout << "       " << ComStr.str() << std::endl;
      //    std::cout << ComStr.str() << std::endl;
      System.Execute(ComStr.str());
      while(std::getline(System,comline))
        n++;
      ComStr.str("");
      ComStr << "sed -i 's/" << olower << "/" << nlower << "/g' " << filename;
      if(verb > 1)
        std::cout << "       " << ComStr.str() << std::endl;
      //    std::cout << ComStr.str() << std::endl;
      System.Execute(ComStr.str());
      while(std::getline(System,comline))
        n++;
    }
    if(verb > 1)
      std::cout << "   Done processing file contents." << std::endl
                << "   Renaming files..." << std::endl;
    ComStr.str("");
    // Now the inside of all files is fixed, need to fix filenames
    ComStr << "find . -name \"*" << OriginalName << "*\"";
    System.Execute(ComStr.str());
    std::string::size_type olen = OriginalName.length();
    std::string::size_type nlen = NewName.length();
    while(std::getline(System,comline)){
      if(verb > 1) 
        std::cout << "     Renaming " << comline << " to ";
      std::string newname(comline);
      std::string::size_type x = newname.find(OriginalName);
      while(x != std::string::npos){
        newname.replace(x,olen,NewName);
        x = newname.find(OriginalName);
      }
      ComStr.str("");
      if(verb > 1)
        std::cout << newname << std::endl;
      ComStr << "mv " << comline << " " << newname;
      //    std::cout << ComStr.str() << std::endl;
      IRAD::Sys::InProcess MV(ComStr.str());
    }
    if(verb > 1)
      std::cout << "done." << std::endl;
    return(0);  
  }
};

int main(int argc,char *argv[])
{
  return(GridConversion::MakeProject(argc,argv));
}
