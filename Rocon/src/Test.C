#include <iostream>
#include <string>
#include <sstream>

#include "Global.H"
#include "Profiler.H"
#include "COMM.H"

///
/// \brief ComLineObject for testing app
///
class TestComLine : public ComLineObject
{
public:
  TestComLine()
    : ComLineObject()
  {};
  TestComLine(const char *args[])
    : ComLineObject(args)
  {};
  void Initialize(){
    AddOption('h',"help");
    AddOption('p',"partition",2,"number");
    AddOption('v',"verb",1,"level");
    AddOption('c',"clone",2,"number");
    AddOption('a',"assembly");
    AddOption('m',"mesh");
    AddOption('d',"debug");
    AddOption('k',"checking");
    AddOption('r',"reorient");
    AddOption('g',"generate");
    AddOption('s',"sparse");
    AddOption('t',"metis");
    AddOption('n',"renumber");
    AddArgument("input",1);
    AddHelp("metis","Metis testing stub.");
    AddHelp("sparse","Write out the sparse matrix for visualization.");
    AddHelp("clone","Generate <number> of partitions identical to the input mesh.");
    AddHelp("generate","Generate a uniform mesh with N nodes and quit.");
    AddHelp("checking","Paranoid and insanely verbose dumping of all important arrays to Log.");
    AddHelp("partition","Performs Metis partitioning of input mesh into <number> partitions.");
    AddHelp("assembly","Performs assembly test.");
    AddHelp("mesh","Performs mesh tests.");
    AddHelp("help","Prints this long version of help.");
    AddHelp("verb","Makes the test more verbose. Default level is 1.");
    AddHelp("config","Specifies the name of the configuration file.");
    AddHelp("out","Specifies the name of the output file.");
    AddHelp("renumber","Uses ParMETIS to do optimal graph reordering.");
    AddArgHelp("input","Mode dependent arguments");
    std::ostringstream Ostr;
    //    Ostr << "Use fixed problem size in scalability analysis.  Only makes"
    //	 << "\n\t\tsense when scalability mode is enabled.";
    //    Ostr.str("");
    Ostr << "Test tool for exercising the mesh library.";
    _description.assign(Ostr.str());
  };
};
  
typedef Global::ParallelGlobalObj<Comm::CommunicatorObject,std::string,int,Profiler::ProfilerObj> TestGlobal;

class TestProgram : public Global::Program<TestGlobal,TestComLine>
{
protected:
  std::string sverb, spart, sclone;
  bool do_part,do_meshtest,do_assem,do_orient,debug;
  bool do_check, do_gen, do_clone, do_dump, do_renum;
  
public:
  TestProgram() :
    Global::Program<TestGlobal,TestComLine>()
  {};
  TestProgram(int nargs,char **args) :
    Global::Program<TestGlobal,TestComLine>(nargs,args)
  {};
  virtual int Initialize()
  {
    int retval = Global::Program<TestGlobal,TestComLine>::Initialize();
    sverb       =  _command_line.GetOption("verb");
    spart       =  _command_line.GetOption("partition");
    sclone      =  _command_line.GetOption("clone");
    do_part     = !_command_line.GetOption("partition").empty();
    do_meshtest = !_command_line.GetOption("mesh").empty();
    do_assem    = !_command_line.GetOption("assembly").empty();
    do_orient   = !_command_line.GetOption("reorient").empty();
    debug       = !_command_line.GetOption("debug").empty();
    do_check    = !_command_line.GetOption("checking").empty();
    do_gen      = !_command_line.GetOption("generate").empty();
    do_clone    = !_command_line.GetOption("clone").empty();
    do_renum    = !_command_line.GetOption("renumber").empty();
    do_dump     = !_command_line.GetOption("sparse").empty();
    if(!_command_line.GetOption("help").empty()){
      if(_OutStream) 
	*_OutStream << _command_line.LongUsage() << std::endl;
      _communicator.SetExit(1);
    }
    if(_communicator.Check())
      return(1);
    if(retval){
      if(_ErrStream) 
	*_ErrStream << _command_line.ErrorReport() << std::endl
		    << std::endl << _command_line.ShortUsage() << std::endl;
      _communicator.SetExit(1);
    }
    if(_communicator.Check())
      return(1);
    return(0);
  };
};

 
int main(int argc,char *argv[])
{
  TestProgram MyProgram(argc,argv);
  if(MyProgram.Initialize())
    return(1);
  MyProgram.OutStream() << "Calling RUN." << std::endl;
  MyProgram.StdOut("Testing\n");
  if(MyProgram.Run())
    return(1);
  std::cout << "Calling Finalize." << std::endl;
  if(MyProgram.Finalize())
    return(1);
  std::cout << "All done." << std::endl;
  return(0);
}
   
