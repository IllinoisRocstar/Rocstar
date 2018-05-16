

class FlotsamComLine : public IRAD::Util::ComLineObject
{
public:
  FlotsamComLine(const char *args[])
    : ComLineObject(args)
  {};
  void Initialize(){
    AddOption('h',"help");
    AddOption('v',"verb",1,"level");
    AddOption('d',"debug");
    AddOption('k',"checking");
    AddOption('c',"config");
    AddArgument("source",1);
    AddArgument("timestamp",1);
    AddArgument("target",1);
    AddHelp("checking","Paranoid and insanely verbose dumping of all important arrays to Log.");
    AddHelp("help","Prints this long version of help.");
    AddHelp("verb","Makes the test more verbose. Default level is 1.");
    AddHelp("config","Specifies the name of the configuration file.");
    AddArgHelp("source","Specifies the name of the source data set.");
    AddArgHelp("target","Specifies the name of the target grid file.");
    AddArgHelp("timestamp","Specifies the timestamp on the source solution.");
    std::ostringstream Ostr;
    //    Ostr << "Use fixed problem size in scalability analysis.  Only makes"
    //	 << "\n\t\tsense when scalability mode is enabled.";
    //    Ostr.str("");
    Ostr << "Rocflo tool for transferring a solution to another mesh.";
    _description.assign(Ostr.str());
  };
};

