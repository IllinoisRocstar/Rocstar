///
/// \file
/// \ingroup irad_group
/// \brief ComLineObject header.
/// 
/// Specification and interface for the ComLineObject command line
/// processing utility class.
///
#ifndef _COMLINE_H_
#define _COMLINE_H_
#include <string>
#include <vector>
#include <map>
#include <istream>
#include <ostream>
#include <utility>


namespace IRAD { 

  namespace Util {

    /// 
    /// Command line processing.
    ///
    /// Provides a quick and easy way to deal with the arguments
    /// coming into an application from the command line.
    /// Several useful methods for parsing the command line tokens and  
    /// generation of detailed usage information are provided for the 
    /// application.
    ///
    /// Command line tokens are separated into three kinds.
    /// The zeroth kind is the first command line token which is the 
    /// name of the program as it was invoked on the command line. This
    /// token is stripped of any absolute path information and stored.
    /// The first kind are _options_.
    /// Any command line token which begins with a "-" is interpreted as
    /// an option. There are several types of options:
    /// * Type 0) A simple flag.
    /// * Type 1) May accept an _option argument_. An option argument is the 
    /// command line token immediatly following the option token.
    /// * Type 2) Requires an option argument.
    /// * Type 3) Non-optional option.  Application requires this option to 
    /// function.
    /// 
    /// Option usage is summarized by the following table.
    ///
    /// | Type  |     Description    |      Usage     |    Value        |
    /// | ----: | :----------------: | :------------: | :-------------- |
    /// |  0    | simple flag        | -t             | .true.          |
    /// |  1    | argument optional  | -t or -t [arg] | .true. or [arg] |
    /// |  2    | argument required  | -t <arg>       | <arg>           |
    /// |  3    | required w/arg     | -t <arg>       | <arg>           |
    ///
    /// The second kind of command line tokens are _arguments_ to the application.
    /// These kind of tokens are simply all the tokens not used up in the 
    /// processing of the zeroth and first kinds of command line tokens. A 
    /// mechanism for specifying whether these kind of tokens are required or 
    /// optional is also provided.
    ///
    class ComLineObject : public std::vector< std::pair<char,std::string> >
    {
      friend std::ostream &operator<<(std::ostream &Out,const ComLineObject &cl);
      friend std::istream &operator>>(std::istream &In,ComLineObject &cl);
    protected:
      ///           
      /// application description.
      ///
      /// This description will be dumped as the first section 
      /// of the output of LongUsage().
      ///
      std::string _description; 
      /// Notes to be displayed at the end of LongUsage().
      std::string _notes;
      /// the name of the program
      std::string _program_name; 
      /// unformatted command line
      std::string _line;         
      /// stores error messages
      std::vector<std::string> _error_messages; 
      /// stores non-option arguments
      std::vector<std::string> _nonops; 
      /// stores the value of each option
      std::map<char,std::string> _options;
      /// stores the help string for each op
      std::map<char,std::string> _help;     
      /// stores a name for arguments
      std::map<char,std::string> _argname;
      /// stores the type of option   
      std::map<char,int> _type;       
      /// application arguments;
      std::vector<std::pair<std::string,int> > _args;
      /// help string for args
      std::map<std::string,std::string> _arghelp;
      /// commandline tokens
      std::vector<std::string> _toks;
    public:
      ///
      /// Default constructor.
      ///
      ComLineObject(){};
      ///
      /// Constructor.
      ///
      /// Typically, args comes directly from \e argv.  One must be careful about
      /// environments in which extra arguments are stored in argv.  A good example
      /// is MPI.  MPI_Init(&argc,&argv) should be invoked before any ComLineObject is 
      /// constructed.
      ///
      ComLineObject(const char *args[])
      {
	Record(args);
      };
      ///
      /// Copy constructor.
      ///
      ComLineObject(ComLineObject &incom);
      ///
      /// Copy method.
      ///
      void Copy(ComLineObject &incom);
      ///
      /// One fell swoop processing of command line.
      ///
      int ProcessCommandLine(const char *args[])
      {
	Record(args);
	return(ProcessOptions());
      };
      ///
      /// \brief Raw Command Line Access.
      ///
      std::string GetRawComLine() const { return (_program_name+std::string(" ")+_line);};
      ///
      /// \brief Raw Command Line Access.
      ///
      void SetRawComLine(const std::string &incl);
      ///
      /// \brief Processes all command line tokens.
      ///
      /// This function examines all the command line tokens, sorts them into
      /// their respective types, and assigns their values.  The number of 
      /// errors encountered during processing is returned.
      ///
      int ProcessOptions();
      ///
      /// \brief Minimal recording of command line.
      ///
      /// This function simply parses the command line to extract the executble
      /// name and record the full command line in a string.
      ///
      void Record(const char *args[]);
      /// 
      /// \brief User interface to describe simple option.
      ///
      /// This function is provided for the specification of a simple option.
      /// The arguments are the short, single letter option, the long "word"
      /// version, and the type of the option.  If the option is of the type
      /// which may have or requires an argument itself, then the argument
      /// will be given the default name, "arg".
      ///
      void AddOption(char s,const std::string &l,int = 0);
      /// 
      /// \brief User interface to describe option with argument.
      ///
      /// This function is provided for the specification of an option along
      /// with a name for it's argument.  Argument names are useful for the
      /// generation of usage information to the user.
      /// The arguments are the short, single letter option, the long "word"
      /// version, the type, and the argument name.
      ///
      void AddOption(char s,const std::string &l,int,const std::string argname);
      /// 
      /// \brief User interface to describe an application argument.
      ///
      /// This function is provided for the specification of an argument
      /// to the application (i.e. a non-option argument), and whether the
      /// argument is required for the application to function.
      /// The function arguments are the name of the argument, and an 
      /// integer which is > 0 for required application arguments. 
      ///
      void AddArgument(const std::string &a,int reqd = 0)
      {
	//_args.push_back(std::make_pair<std::string,int>(a,reqd));
	_args.push_back(std::make_pair(a,reqd));
      };
      /// 
      /// \brief Specify the usage info for application argument.
      ///
      /// This function is provided for the specification of the usage
      /// information for application arguments.  This is called after
      /// the argument has been set by AddArgument(...).  The string
      /// specified in "help" string is displayed in the output of 
      /// LongUsage().
      ///
      void AddArgHelp(const std::string &a,const std::string &help)
      {
	_arghelp[a] = help;
      };
      /// 
      /// \brief Specify name of an option argument.
      ///
      /// Function is used to set the argument name for an option 
      /// if the option can use an argument.
      ///
      void SetArgName(const std::string opstr,const std::string argname)
      {
	_argname[GetShort(opstr)] = argname;
      };
      /// 
      /// \brief Specify name of an option argument.
      ///
      /// Function is used to set the argument name for an option 
      /// if the option can use an argument.
      ///
      void SetArgName(char s,const std::string &argname)
      {
	_argname[s] = argname;
      };
      /// 
      /// \brief Specify usage for an option.
      ///
      /// Function is used to set the usage string for an option 
      /// which is displayed in the output of LongUsage().
      ///
      void AddHelp(char s,const std::string &help){ _help[s] = help;};
      /// 
      /// \brief Specify usage for an option.
      ///
      /// Function is used to set the usage string for an option 
      /// which is displayed in the output of LongUsage().
      ///
      void AddHelp(const std::string &l,const std::string &help)
      {
	_help[GetShort(l)] = help;
      };
      ///
      /// \brief Obtains option strings by type.
      ///
      /// Internal utility function for building a string of 
      /// option characters which are of a certain type.
      ///
      std::string GetOpStringByType(int mintype,int maxtype); 
      ///
      /// \brief Obtain the usage string for an option.
      ///
      std::string GetHelp(char s){return(_help[s]);};
      ///
      /// \brief Obtain the usage string for an option.
      ///
      std::string GetHelp(const std::string &l){return(_help[GetShort(l)]);};
      ///
      /// \brief Obtain the short one char option from the long word version.
      ///
      char GetShort(const std::string &l);
      ///
      /// \brief Obtain the long word option from the short char version.
      ///
      std::string GetLong(const char &s);
      ///
      /// \brief Generate short usage string.
      ///
      /// Generates a short usage string similar to what one would see
      /// in the typical short command line usage output for a UNIX 
      /// command.  Note that the format of this usage string follows
      /// the convention that optional quantities are enclosed in []'s
      /// while required quantities are enclosed in <>'s.
      ///
      std::string ShortUsage();
      ///
      /// \brief Generate long usage string.
      ///
      /// Generates a long usage string similar to what one would see
      /// in the typical UNIX man-like description.  Each option is 
      /// shown with it's short and long version, along with usage 
      /// information for each option and argument if such information 
      /// has been specified with AddHelp, or AddArgHelp, respectively. 
      /// Note that the format of this usage string follows
      /// the convention that optional quantities are enclosed in []'s
      /// while required quantities are enclosed in <>'s.
      ///
      std::string LongUsage();
      ///
      /// \brief Get the value of an option.
      ///
      /// Returns the value of the specified option.  If the option 
      /// has not been set on the command line, then the value will
      /// be an empty string.
      ///
      std::string GetOption(const char &s){return(_options[s]);};
      ///
      /// \brief Get the value of an option.
      ///
      /// Returns the value of the specified option.  If the option 
      /// has not been set on the command line, then the value will
      /// be an empty string.
      ///
      std::string GetOption(const std::string &l){return(_options[GetShort(l)]);};
      ///
      /// \brief Error reporting.
      ///
      /// Returns a report of the errors encountered in the ProcessOptions() 
      /// function.
      ///
      std::string ErrorReport();
      ///
      /// \brief Program name access.
      ///
      std::string ProgramName() const { return _program_name;};
      ///
      /// \brief Argument access.
      ///
      std::vector<std::string> GetArgs() const { return _nonops;};
      ///
      /// \brief Set description string.
      ///
      void SetDescription(const std::string &desc){ _description.assign(desc); };
      ///
      /// \brief Set notes string.
      ///
      void SetNotes(const std::string &notes){ _notes.assign(notes); };
      ///
      /// \brief Write an RC string that can be used for config.
      ///
      void WriteRC(std::ostream &Ostr) const;
      ///
      /// \brief Read a config from RC string.
      ///
      void ReadRC(const std::string &RCString);
      ///
      /// \brief Process a range in the format "t1:t2".
      ///
      template<typename NumType>
      void ProcessRange(NumType &t1,NumType &t2,const std::string stinter)
      {
	if(!stinter.empty()){
	  std::string::size_type x = stinter.find(":");
	  std::string t1s = stinter;
	  std::string t2s = stinter;
	  if(x != std::string::npos){
	    t1s = stinter.substr(0,x);
	    t2s = stinter.substr(x+1);
	    if(!t2s.empty()){
	      std::istringstream Istr(t2s);
	      Istr >> t2;
	    }
	  }
	  else {
	    t1s.erase();
	  }
	  if(!t1s.empty()){
	    std::istringstream Istr(t1s);
	    Istr >> t1;
	  }
	  if(!t2s.empty()){
	    std::istringstream Istr(t2s);
	    Istr >> t2;
	  }
  
	}
      };
      ///
      /// \brief Resolve an option (i.e. obtain it's "value").
      ///
      template<typename NumType>
      void ResolveOption(NumType &t2,const std::string stinter)
      {
	if(stinter == ".true." || stinter == ".TRUE."){
	  t2 = 1;
	  return;
	}
	if(stinter == ".false." || stinter == ".FALSE."){
	  t2 = 0;
	  return;
	}
	if(!stinter.empty()){
	  std::istringstream Istr(stinter);
	  Istr >> t2;
	}
      };
      ///
      /// \brief virtual function for program specific Initialization.
      ///
      virtual void Initialize(void){};
      virtual ~ComLineObject(){};
    };
  };
};
#endif


