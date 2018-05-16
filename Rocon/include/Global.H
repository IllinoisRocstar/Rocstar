///
/// \file
/// \ingroup irad_group
/// \brief Base global object header
///
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "ComLine.H"
#include "UnixUtils.H"
#include "primitive_utilities.H"

namespace IRAD {

  /// Base global namespace
  namespace Global {
    ///
    /// \brief Base global object
    ///
    /// Provides an interface for tracking errors, profiling,
    /// and user defined stack information.  The idea is that this object or a
    /// derived object will be available to every \e rank  \e local method.
    ///
    template<typename StackType,typename ErrorType,typename ProfilerType>
    class GlobalObj {
    private:

      /// Debugging level 
      unsigned char _debug_level;

      /// Logging level 
      unsigned char _log_level;

      /// Verbosity level
      unsigned char _verb_level;

      /// Stack provides a list for storing the user defined "stack".  
      /// It's not a real call stack, but a stack of user defined code constructs.
      std::list< StackType > _Stack;

      /// Errors provides a list for storing caught execution errors for later processing.
      std::list< ErrorType > _Errors;
    
      /// has been initialized?
      bool _initd;

    protected:
      /// DebugStream is useful for parallel programs that cannot stream 
      /// debugging information to the screen.  Each rank can use this to 
      /// stream debugging output to it's own file for separate analysis.
      std::ostream          *_DebugStream;
      std::ostream          *_OutStream;
      std::ostream          *_LogStream;
      std::ostream          *_ErrStream;
      /// This is a placeholder for a generic profiling tool
      ProfilerType           *_Profiler;
    private:
      bool                   _profiling_enabled;
      bool                   _profiler_mine;
      std::string            _name;
      std::ostream          *_OutStreamSave;
      std::ostream          *_ErrStreamSave;
      std::ostream          *_DebugStreamSave;
    public:

      /// Default constructor.  Profiling is OFF by default
      GlobalObj() : _debug_level(0), _log_level(0), _verb_level(0),
		    _DebugStream(NULL), 
		    _OutStream(&std::cout), _LogStream(NULL),
		    _ErrStream(&std::cerr),_Profiler(NULL),
		    _profiling_enabled(false),_profiler_mine(false), 
		    _name(""), _OutStreamSave(NULL), _ErrStreamSave(NULL),
		    _DebugStreamSave(NULL)
      {
	//      _Stack.push_back("main");
	//      std::cout << "default constr" << std::endl;
      };
    
      /// Default constructor.  If no args, then profiling is OFF
      /// by default.
      GlobalObj(int narg,char **args) :
	_debug_level(0), _log_level(0), _verb_level(0),
	_DebugStream(NULL), 
	_OutStream(&std::cout), _LogStream(NULL),
	_ErrStream(&std::cerr),_Profiler(NULL),
	_profiling_enabled(true),_profiler_mine(false),
	_name(""), _OutStreamSave(NULL), _ErrStreamSave(NULL),
	_DebugStreamSave(NULL)
      {
	//      std::cout << "std args constr" << std::endl;
	Init(Util::stripdirs(args[0]));
      };

      /// Constructor with args.   If this constructor is used, Init is
      /// invoked which turns profiling ON by default.
      GlobalObj(GlobalObj<StackType,ErrorType,ProfilerType> &globin) :
	_debug_level(globin._debug_level), _log_level(globin._log_level),_verb_level(globin._verb_level),
	_DebugStream(globin._DebugStream), _OutStream(globin._OutStream), 
	_LogStream(globin._LogStream), _ErrStream(globin._ErrStream),
	_Profiler(globin._Profiler),_profiling_enabled(globin._profiling_enabled),
	_profiler_mine(false), _name(globin._name), _OutStreamSave(globin._OutStreamSave),
	_ErrStreamSave(globin._ErrStreamSave), _DebugStreamSave(globin._DebugStreamSave)
      {
	//      std::cout << "copy constr" << std::endl;
	if(!globin._Stack.empty()){
	  typename std::list<StackType>::iterator si = globin._Stack.begin();
	  while(si != globin._Stack.end())
	    _Stack.push_back(*si++);
	}
	_initd = false;
      };
    

    
      /// Constructor with args.   If this constructor is used, Init is
      /// invoked which turns profiling ON by default.
      GlobalObj(const StackType &name) : _debug_level(0) , _log_level(0), _verb_level(0)
      {
	//      std::cout << "named constr" << std::endl;

	Init(name);
      };

      /// Constructor with args.   If this constructor is used, Init is
      /// invoked which turns profiling ON by default.
      GlobalObj(const StackType &name,int id) : _debug_level(0), _log_level(0), _verb_level(0)
      {
	//      std::cout << "named&&numbered constr" << std::endl;
	Init(name,id);
      };

      /// Initializes the global object and it's profiler object    /// Initializes the global object and it's profiler object
      /// Profiling is ON by default if Init is invoked
      virtual int Init(const StackType &name)
      {
	unsigned int id = 0;
	_debug_level = 0;
	_log_level = 0;
	_verb_level = 0;
	_name = name;
	//      _profiling_enabled = true;
	if(_Profiler && _profiler_mine)
	  delete _Profiler;
	_Profiler = new ProfilerType;
	_Profiler->Init(name,id);
	_profiler_mine = true;
	_profiling_enabled = true;
	//      _DebugStream = &std::cout;
	_DebugStream = NULL;
	_ErrStream = &std::cerr;
	_OutStream = &std::cout;
	_LogStream = NULL;
	_OutStreamSave =   _OutStream;
	_ErrStreamSave =   _ErrStream;
	_DebugStreamSave = _DebugStream;
	//      _Stack.clear();
	_Stack.push_front(name);
	_Errors.clear();
	_initd = true;
	return(0);
      };

      /// Initializes the global object and it's profiler object
      /// Profiling is ON by default if Init is invoked
      virtual int Init(const StackType &name,unsigned int id)
      {
	_name = name;
	_verb_level = 0;
	_debug_level = 0;
	_log_level = 0;
	_profiling_enabled = true;
	if(_Profiler && _profiler_mine)
	  delete _Profiler;
	_Profiler = new ProfilerType;
	_Profiler->Init(name,id);
	_profiler_mine = true;
	//      _DebugStream = &std::cout;
	_DebugStream = NULL;
	if(id == 0){
	  _DebugStream = &std::cout;
	  _ErrStream   = &std::cerr;
	  _OutStream   = &std::cout;
	}
	else{
	  _OutStream = NULL;
	  _ErrStream = NULL;
	}
	_LogStream = NULL;
	_OutStreamSave =   _OutStream;
	_ErrStreamSave =   _ErrStream;
	_DebugStreamSave = _DebugStream;
	//      _Stack.clear();
	_Stack.push_front(name);
	_Errors.clear();
	_initd = true;
	return(0);
      };

      /// Finalizes the global object, and it's profiler object
      virtual int Finalize()
      {
	int errcode = 0;
	if(_initd){
	  if(_profiler_mine){
	    if(_Profiler->FinalizeReady() && _Profiler)
	      _Profiler->Finalize();
	    else 
	      errcode = 1;
	  }
	  _profiling_enabled = false;
	  _initd = false;
	}
	return(errcode);
      };

      /// Get profiling state
      virtual bool Profiling(){
	return(_profiling_enabled);
      };

      /// Enable or disable profiling
      virtual bool Profiling(bool yn){
	_profiling_enabled = yn;
	return(_profiling_enabled);
      };

      /// Creates the final profiling report on the specified stream
      virtual void Report(std::ostream &ReportStream)
      {
	if(_Profiler){
	  _Profiler->Finalize();
	  _Profiler->SummarizeSerialExecution(ReportStream);
	}
      };

      /// Register will push to the Stack only.  Can be called upon entry 
      /// of user defined code constructs where profiling is not desired.
      virtual void Register(const StackType &stackentry){
	_Stack.push_front(FName()+stackentry);
	if(_debug_level > 1 && _DebugStream)
	  *_DebugStream << "GlobalDebugging:Register(" << FName()+stackentry << ")" << std::endl;
      }

      /// DeRegister will pop from the stack only.  Should be called upon exit 
      /// of user defined code constructs where a previous call to Register was 
      /// made upon entry.
      virtual void DeRegister(const StackType &stackentry){
	_Stack.pop_front();
	if(_debug_level > 1 && _DebugStream)
	  *_DebugStream << "GlobalDebugging:DeRegister(" << FName() + stackentry << ")" << std::endl;
      }

      /// FunctionEntry updates the Stack as well as the Profiler
      virtual void FunctionEntry(const StackType &stackentry){
	_Stack.push_front(FName()+stackentry);
	if(_profiling_enabled && _Profiler)
	  _Profiler->FunctionEntry(FName()+stackentry);
	if(_debug_level > 1 && _DebugStream)
	  *_DebugStream << "GlobalDebugging:FunctionEntry(" << FName() + stackentry << ")" << std::endl;
      };

      /// FunctionExit updates the Stack as well as the Profiler
      virtual void FunctionExit(const StackType &stackentry){
	_Stack.pop_front();
	if(_profiling_enabled && _Profiler)
	  _Profiler->FunctionExit(FName()+stackentry);
	if(_debug_level > 1 && _DebugStream)
	  *_DebugStream << "GlobalDebugging:FunctionExit(" << FName() + stackentry << ")" << std::endl;
      };

      /// Pushes an error onto the error stack
      virtual void Error(const ErrorType &error){
	_Errors.push_front(error);
      };

      /// Dumps the errors themselves to the indicated stream
      virtual void DumpErrors(std::ostream &Ostr){
	Util::DumpContents(Ostr,_Errors);
      };

      /// Dumps the code construct statck to the indicated stream
      virtual void DumpStack(std::ostream &Ostr){
	Util::DumpContents(Ostr,_Stack);
      };

      /// Dumps the timing profile to the indicated stream
      virtual void DumpProfile(std::ostream &Ostr){
	if(_Profiler)
	  _Profiler->Dump(Ostr);
      };
    
      /// Checks debug stream 
      bool DebugStreamReady() { return (_DebugStream != NULL); };

      /// Sets the debug stream object
      void SetDebugStream(std::ostream &dbstream){
	_DebugStream = &dbstream;
      };

      /// Gets the debug stream object
      std::ostream &DebugStream(){ assert(_DebugStream); return(*_DebugStream); };

      /// Gets the debug stream pointer
      std::ostream *DebugStreamPtr(){ return(_DebugStream); };

      /// Checks error stream 
      bool ErrStreamReady() { return (_ErrStream != NULL); };

      /// Sets the error stream object
      void SetErrStream(std::ostream &dbstream){
	_ErrStream = &dbstream;
      };

      /// Gets the error stream object
      std::ostream &ErrStream(){ assert(_ErrStream); return(*_ErrStream); };

      /// Gets the error stream pointer
      std::ostream *ErrStreamPtr(){ return(_ErrStream); };

      /// Checks out stream 
      bool OutStreamReady() { return (_OutStream != NULL); };

      /// Sets the stdout stream object
      void SetOutStream(std::ostream &outstream){
	_OutStream = &outstream;
      };

      /// Gets the stdout stream object
      std::ostream &OutStream(){ assert(_OutStream); return(*_OutStream); };

      /// Gets the stdout stream pointer
      std::ostream *OutStreamPtr(){ return(_OutStream); };

      /// Checks log stream 
      bool LogStreamReady() { return (_LogStream != NULL); };

      /// Sets the stdout stream object
      void SetLogStream(std::ostream &logstream){
	_LogStream = &logstream;
      };

      /// Gets the log stream object
      std::ostream &LogStream(){ assert(_LogStream); return(*_LogStream); };

      /// Gets the log stream pointer
      std::ostream *LogStreamPtr(){ return(_LogStream); };

      /// Sets the Profiler object
      void SetProfiler(ProfilerType &profiler){
	if(_profiler_mine){
	  if(_profiling_enabled && _Profiler && _Profiler->FinalizeReady())
	    _Profiler->Finalize();
	  if(_Profiler)
	    delete _Profiler;
	}
	_Profiler = &profiler;
      };
    
      virtual void SetName(const std::string &name) { _name = FName()+name; };
      virtual void ReName(const std::string &name) { _name = name; };
      virtual const std::string &Name(){return(_name);};
      /// Gets the debug stream object
      ProfilerType &Profiler(){ assert(_Profiler); return(*_Profiler); };
    
      /// Sets the level of debugging
      void SetDebugLevel(unsigned char l){  _debug_level = l; };

      /// Returns the debugging level
      unsigned char DebugLevel(){ return(_debug_level); };

      void SetVerbLevel(unsigned char l) { _verb_level = l; };

      /// Returns the verbosity level
      unsigned char VerbLevel(){ return (_verb_level); };

      void SetLogLevel(unsigned char l) { _log_level = l; };

      /// Returns the verbosity level
      unsigned char LogLevel(){ return (_log_level); };
    
      /// DeSync IO - relinquish control over IO, revert to standard
      void DeSyncIO()
      {
	_OutStreamSave = _OutStream;
	_ErrStreamSave = _ErrStream;
	_DebugStreamSave = _DebugStream;
	_OutStream = &std::cout;
	_ErrStream = &std::cerr;
	_DebugStream = &std::cout;
      };

      /// Syncronize IO - assume control over IO, using established streams
      void SyncIO()
      {
	_OutStream   = _OutStreamSave;      
	_ErrStream   = _ErrStreamSave;
	_DebugStream = _DebugStreamSave;
      };

      virtual int StdOut(const std::string &outstr,unsigned char inlev = 0,bool stack = false)
      {
        if(_OutStream && _verb_level >= inlev){
	  std::string line;
	  std::istringstream Istr(outstr);
	  while(std::getline(Istr,line)){
            if(stack){
              *_OutStream << _Stack.front() << ": " << line;
            }
            else{
	      *_OutStream << (Name().empty() ? "" : std::string(Name()+std::string(": "))) 
	     		  << line;
            }
	    *_OutStream << (Istr.good() ? "\n" : "");
	  }
	}
	else
	  return(1);
	return(0);
      };

      virtual int ErrOut(const std::string &outstr)
      {
	if(_ErrStream){
	  std::string line;
	  std::istringstream Istr(outstr);
	  while(std::getline(Istr,line)){
	    *_ErrStream << (Name().empty() ? "" : std::string(Name()+std::string(":"))) 
			<< "Error: " << line;
	    *_ErrStream << (Istr.good() ? "\n" : "");
	  }
	}
	else
	  return(1);
	return(0);
      };
    
      virtual int DebugOut(const std::string &outstr,int inlev=1)
      {
	if(_DebugStream && _debug_level >= inlev){
	  std::string line;
	  std::istringstream Istr(outstr);
	  while(std::getline(Istr,line)){
	    *_DebugStream << (Name().empty() ? "" : std::string(Name()+std::string(":")))
			  << "db: " << line;
	    *_DebugStream << (Istr.good() ? "\n" : "");
	  }
	}
	else
	  return(1);
	return(0);
      };

      virtual int LogOut(const std::string &outstr,int inlev=0)
      {
	std::string timestring(Sys::LogTime());
	if((_LogStream != NULL) && (_log_level >= inlev)){
	  std::string line;
	  std::istringstream Istr(outstr);
	  while(std::getline(Istr,line)){
	    *_LogStream << timestring 
			<< (Name().empty() ? "" : std::string(Name()+std::string(": ")))
			<< line;
	    *_LogStream << (Istr.good() ? "\n" : "");
	  }
	}
	else
	  return(1);
	return(0);
      };

      /// Destructor
      virtual ~GlobalObj()
      { 
	if(_initd) 
	  Finalize(); 
	if(_profiler_mine && _Profiler){
	  delete _Profiler;
	}
      };
    private:
      std::string FName() { return((Name().empty() ? 
				    std::string("") :
				    std::string(Name()+std::string(":")))); };
    };
  
    template<typename CommunicatorType,typename StackType,
	     typename ErrorType,typename ProfilerType>
    class ParallelGlobalObj : public GlobalObj<StackType,ErrorType,ProfilerType> 
    {
      typedef GlobalObj<StackType,ErrorType,ProfilerType> GlobalType;
    protected:
      CommunicatorType _communicator;
    public:
      ParallelGlobalObj() :
	GlobalType()
      {};
      ParallelGlobalObj(GlobalType &inglob) :
	GlobalType(inglob)
      {};
      ParallelGlobalObj(ParallelGlobalObj<CommunicatorType,StackType,ErrorType,ProfilerType> &pglobin) :
	GlobalObj<StackType,ErrorType,ProfilerType>(pglobin)
      {
	_communicator.Initialize(pglobin._communicator);
      };
      ParallelGlobalObj(const std::string &name) :
	GlobalObj<StackType,ErrorType,ProfilerType>(name)
      {
	//      Init(name);
      };
      ParallelGlobalObj(const std::string &name,unsigned int id) :
	GlobalObj<StackType,ErrorType,ProfilerType>(name,id)
      {
	//      Init(name);
      };
      ParallelGlobalObj(int narg,char **args) :
	GlobalObj<StackType,ErrorType,ProfilerType>(narg,args)
      {
	assert(Init(narg,args)==0);
      };
      virtual int Init(const std::string &name,CommunicatorType &incomm)
      {
	int retval = _communicator.Initialize(incomm);
	if(retval)
	  return(retval);
	return(GlobalObj<StackType,ErrorType,ProfilerType>::Init(name,_communicator.Rank()));
      };
      virtual int Init(const std::string &name)
      {
	return(GlobalObj<StackType,ErrorType,ProfilerType>::Init(name));
      };
      virtual int Init(const std::string &name,unsigned int id) 
      {
	return(GlobalObj<StackType,ErrorType,ProfilerType>::Init(name,id));
      };
      virtual int Init(int narg,char **args)
      {
	int retval = _communicator.Initialize(&narg,&args);
	if(retval)
	  return(retval);
	return(GlobalObj<StackType,ErrorType,ProfilerType>::Init(Util::stripdirs(args[0]),_communicator.Rank()));
      };
      virtual int Rank()
      {
	return(_communicator.Rank());
      };
      virtual int NProc()
      {
	return(_communicator.Size());
      };
      virtual int Finalize()
      {
	int retval = GlobalObj<StackType,ErrorType,ProfilerType>::Finalize();
	if(retval)
	  return(retval);
	return(_communicator.Finalize());
      };
      virtual CommunicatorType &Communicator(){return _communicator; };
      virtual ~ParallelGlobalObj(){ Finalize();};
    };

    template<typename GlobalType,typename ComLineType>
    class Program : public GlobalType
    {
    protected:
      ComLineType _command_line;
    public:
      Program() :
	GlobalType()
      {};
      Program(GlobalType &inglob) :
	GlobalType(inglob){};
      Program(ComLineType &incom) :
	GlobalType(),_command_line(incom)
      {};
      Program(ComLineType &incom,GlobalType &inglob) :
	GlobalType(inglob), _command_line(incom) {};
      Program(int narg,char **args) :
	GlobalType(narg,args)
      {
	_command_line.Record((const char **)args);
	//      SetName(_command_line.ProgramName());
      };
      virtual int Initialize()
      {
	_command_line.Initialize();
	return(_command_line.ProcessOptions());
      };
      ComLineType &CommandLine() {return _command_line; };
      virtual int Run(){return(0);};
      virtual int Finalize(){return(GlobalType::Finalize());};
      virtual ~Program(){};
    };
  };
};
#endif
