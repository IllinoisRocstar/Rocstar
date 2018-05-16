/// 
/// @file
/// @ingroup gridconversion_group
/// @brief GridConversion interface
/// @author Jessica Kress (jkress@illinoisrocstar.com)
/// 
#ifndef __EXAMPLE_PROGRAM_H__
#define __EXAMPLE_PROGRAM_H__
#ifdef _GRIDCONVERSION_PARALLEL_
#include "COMM.H"
#endif
#include "Global.H"
#include "Profiler.H"
#include "Mesh.H"

namespace GridConversion {
  
  ///
  /// Encapsulate GridConversion-specific code constructs.
  ///
  /// Based off the example program from IRProject
  /// For Example Program:
  /// This namespace contains all of the base code constructs
  /// that are used in the example programs.  These example
  /// programs have several features which make them very 
  /// easy to integrate in composite software systems. 
  ///
  /// * Program is entirely encapsulated in program object
  /// * All screen output is automatically tagged with the program name
  /// * Common look-and-feel interface
  /// * Interface functionality allows menu-driven integration
  ///
  /// Other convenient features of this API are:
  ///
  /// * Built-in stack tracing for debugging
  /// * Built-in profiling utilities
  /// * Advanced command-line handling
  /// * Built-in parallelism 
  ///
  namespace DriverProgram {
    ///
    /// Convenience type definition for Profiler object.
    ///
    typedef IRAD::Profiler::ProfilerObj ProfilerType;
    ///
    /// Convenience type definition for program stack.
    ///
    typedef std::string                 StackType;
    ///
    /// Convenience type definition for the serial global object.
    ///
    /// The global object encapsulates all of the built-in functionalities
    /// for a program.  Its interface and implementation are defined in 
    /// IRAD/include/Global.H 
    ///
    typedef IRAD::Global::GlobalObj<StackType,int,ProfilerType> GlobalType;
#ifdef _GRIDCONVERSION_PARALLEL_
    ///
    /// Convenience type definition for parallel communication object.
    ///
    typedef IRAD::Comm::CommunicatorObject CommType;
    ///
    /// Convenience type definition for parallel global object
    ///
    typedef IRAD::Global::ParallelGlobalObj<CommType,StackType,int,ProfilerType> PGlobalType;
#endif
    ///
    /// The ComLineObject for the example programs.
    ///
    /// This comline object implements a basic command-line interface
    /// handler for the example programs.  These command-line objects
    /// handle all of the generation of help messages, command line 
    /// argument validation, and encapsulation.
    ///
    class DriverProgramComLine : public IRAD::Util::ComLineObject
    {
    public:
      DriverProgramComLine()
        : IRAD::Util::ComLineObject()
      {};
      DriverProgramComLine(const char *args[])
        : IRAD::Util::ComLineObject(args)
      {};
      ///
      /// This is where the actual options and arguments are described.
      ///
      /// For details about what can be done in this function, see the
      /// documentation for the IRAD::Util::ComLineObject.
      ///
      void Initialize(){
        AddOption('h',"help");
        //AddOption('q',"quadratic");
        AddOption('v',"verbosity",1);
        AddOption('o',"output",2,"filename");
        AddArgument("grid input> <bc input",1);
        AddHelp("help","Prints this long version of help.");
        AddHelp("output","Specifies the name of the output file.");
        std::ostringstream Ostr;
        //Ostr << "If using a mesh with tetrahedral elements, specifying" << std::endl
        //     << "this flag will cause the application to change the 4 node" << std::endl
        //     << " tetrahedra from input to 10 node quadratic tetrahedra";
        //AddHelp("quadratic",Ostr.str());
        Ostr.str("");
        Ostr << "GridConversion requires two arguments in order:" << std::endl
             << "                "
             << "su2 grid file> <bc file> where the su2 is a grid" << std::endl
             << "                 "
             << "analysis file output from Pointwise in Stanford format." << std::endl
             << "                "
             << "The bc file indicates a Rocfrac boundary condition for "
             << "each" << std::endl
             << "                domain." << std::endl
             << "                See the Grid Conversion User Manual for more "
             << "information on" << std::endl
             << "                the formats of these two files" << std::endl;
        AddArgHelp("grid input> <bc input",Ostr.str());
        // The following commented string formation is kept around to 
        // illustrate the required tabs and newlines to make the output 
        // look proper:
        //
        //    Ostr << "Use fixed problem size in scalability analysis.  Only makes"
        //         << "\n\t\tsense when scalability mode is enabled.";
        //    Ostr.str("");
        Ostr.str("");
        Ostr << "Main GridConversion program.";
        _description.assign(Ostr.str());

      };
    };
    
    ///
    /// Convenience type definition for the Main comline object.
    ///
    typedef DriverProgramComLine ComLineType;
    ///
    /// Convenience type definition for the serial program.
    /// 
    /// For a detailed description of the program type, see IRAD::Global::Program
    /// documentation.
    ///
    typedef IRAD::Global::Program<GlobalType,ComLineType> SerialProgramType;
#ifdef _GRIDCONVERSION_PARALLEL_
    ///
    /// Convenience type definition for the parallel program.
    /// 
    /// For a detailed description of the program type, see IRAD::Global::Program
    /// documentation.
    ///
    typedef IRAD::Global::Program<PGlobalType,ComLineType> ParallelProgramType;
#endif
    ///
    /// Implementation of the basic parts of the serial program example.
    ///
    /// This object encapsulates the example serial program.  It inherits
    /// from the IRAD::Global::Program type.  
    ///
    /// The program itself just copies a specified input file to the 
    /// specified output file (or stdout if no file is given). The 
    /// command line usage goes:
    /// 
    ///
    ///            sep Usage:
    /// 
    ///            sep [-h] [-v [arg] -o <filename> ] <input> 
    ///
    ///            -h,--help
    ///               Prints this long version of help.
    ///
    ///                   -v,--verbosity [arg]
    /// 
    ///            -o,--output <filename>
    ///               Specifies the name of the output file.
    ///
    ///                   <input>
    ///               Mode-dependent input: input <filename> for serial example
    ///               program, or <number of divisions> for parallel example.
    ///
    class SerialProgram : public SerialProgramType
    {
    protected:
      /// Name of file for output.
      std::string output_name;
      /// Name of input file.
      std::string input_name;
      /// Name of bc input file
      std::string bc_input_name;
      /// Verbosity level. 
      int verblevel;
      /// Outfile stream for output.
      std::ofstream Ouf;
      /// Infile stream for input.
      std::ifstream Inf;
      /// boolean for changing the elements to higher order
      bool quadratic;
      /// number of nodes in mesh
      int numNodes;
      /// vector for holding nodal coordinates read from input
      std::vector<double> nodes;
      // vector of vectors to hold nodal bc flags
      std::vector< std::vector<int> > nodeBCs;
      /// vector to hold list of node bc values (these are different for
      /// each bc type - see the Rocfrac documentation for more clarification)
      std::vector< std::vector< std::vector<double> > > nodeBCValues;
      // vector to hold the list of nodes to domains
      SolverUtils::Mesh::Connectivity nodesToDomains;
      /// number of element edges in mesh
      int numElemEdges;
      /// vector for holding element edges
      std::vector<unsigned int> elemEdges;
      /// number of nodes per element edge
      int numNodesPerElemEdge;
      /// map of elements to element edges
      std::vector<unsigned int> elemToElemEdges;
      //number of edges per element
      int numEdgesPerElem;
      /// number of elements in mesh
      int numElems;
      /// vector for holding node to node map
      std::vector< std::vector<unsigned int> > nodeToNode;
      /// vector for holding element connectivies read from input
      std::vector<unsigned int> elems;
      /// vector to hold list of element bc values 
      std::vector< std::vector<int> > elemBCs;
      /// vector to hold list of element bc values (these are different for
      /// each bc type - see the Rocfrac documentation for more clarification)
      std::vector< std::vector< std::vector<double> > > elemBCValues;
      /// number of nodes per element
      int numNodesPerElem;
      /// vector to hold the domains for each element
      std::vector< std::vector<int> >elemsToDomains;
      /// integer to denote the shape of elements used in the mesh
      int elemShape;
      /// vector to hold the element shapes (only used for Patran input format)
      std::vector<std::string> shapes;
      /// number of domains in mesh
      int numDomains;
      /// vector of vectors to hold the domains
      std::vector< std::vector<unsigned int> > domains;
      /// vector to hold list of domain bc flags
      std::vector< std::vector<int> > domainBCs;
      /// vector to hold list of domain bc values (these are different for
      /// each bc type - see the Rocfrac documentation for more clarification)
      std::vector< std::vector< std::vector<double> > > domainBCValues;
      /// vector of vectors to hold the edges
      std::vector< std::vector<int> > edges;
      /// vector to hold list of edge bc flags
      std::vector< std::vector<int> > edgeBCs;
      /// vector to hold list of edge bc values (these are different for
      /// each bc type - see the Rocfrac documentation for more clarification)
      std::vector< std::vector< std::vector<double> > > edgeBCValues;
    public:
      ///
      /// Default constructor.
      ///
      SerialProgram() :
        SerialProgramType()
      {};
      ///
      /// Constructor designed to take the commandline args.
      /// 
      SerialProgram(int nargs,char **args) :
        SerialProgramType(nargs,args)
      {};
      ///
      /// Initializes native data structures from commandline args.
      /// 
      virtual int Initialize()
      {
        int retval = SerialProgramType::Initialize();
        if(!_command_line.GetOption("help").empty()){
          std::ostringstream Ostr;
          Ostr << _command_line.LongUsage() << std::endl;
          StdOut(Ostr.str());
          return(-1);
        }
        if(retval){
          std::ostringstream Ostr;
          Ostr << _command_line.ErrorReport() << std::endl
               << std::endl << _command_line.ShortUsage() << std::endl;
          ErrOut(Ostr.str());
          return(retval);
        }

        //Initialize mesh values to 0
        numNodes = 0;
        numElems = 0;
        numDomains = 0;
        numNodesPerElem = 0;
        elemShape = 0;

        //Set element shape vector (shape numbers are those used in Patran grid files)
        shapes.resize(10);
        shapes[0] = shapes[1] = shapes[6] = "none";
        shapes[2]= "bar";
        shapes[3] = "triangle";
        shapes[4] = "quadrilateral";
        shapes[5] = "tetrahedron";
        shapes[7] = "triangular prism";
        shapes[8] = "hexahedron";
        shapes[9] = "pyramid";

        // Check if we should create higher order elements
        quadratic = false;
        if(!_command_line.GetOption("quadratic").empty())
          quadratic = true;        

        // Check if output file is specified
        output_name = _command_line.GetOption("output");

        // Set input file to first argument
        std::vector<std::string> args(_command_line.GetArgs());
        input_name = args[0];
        bc_input_name = args[1];

        // Check the verbosity level
        std::string sverb(_command_line.GetOption("verbosity"));
        if(sverb.empty() || sverb == ".true.")
          verblevel = 1;
        else {
          std::istringstream Vin(sverb);
          Vin >> verblevel;
          if(verblevel < 0)
            verblevel = 1;
        }

        // If high verbosity, stick a configuration blurb to stdout.
        if(verblevel > 1){
          std::ostringstream Ostr;
          Ostr << "Configuration:" << std::endl
               << "verbosity           = " << verblevel << std::endl
               << "input_name          = " << input_name << std::endl
               << "bc_input_name       = " << bc_input_name << std::endl;
          if(!output_name.empty())
            Ostr << "output file         = " << output_name << std::endl;
          Ostr << std::endl;
          StdOut(Ostr.str());
        } 

        return(0);
      };
      ///
      /// Returns verbosity level.
      ///
      int VerbLevel() const { return verblevel;};
      ///
      /// Destructor
      ///
      virtual ~SerialProgram() {
        if(Ouf){
          Ouf.close();
          SetOutStream(std::cout);
        }
      };
      ///
      /// This function implements the main function executed by the program.
      /// 
      virtual int Run();
      ///
      /// This function parses a Patran format input file.
      /// 
      virtual int ReadPatranInput();
      ///
      /// This function parses a Stanford format input file.
      /// 
      virtual int ReadStanfordInput();
      ///
      /// This function creates higher order tets from linear tets.
      /// (For now it only takes 4 node tets and creates only 10 node tets.
      /// 
      virtual int HigherOrderTets(SolverUtils::Mesh::Connectivity elemEdgeToElems,
                                  SolverUtils::Mesh::Connectivity nodesToDomains);
      ///
      /// This function creates a few different connectivity maps
      /// that are needed.
      /// 
      virtual int ConnectivityMaps(SolverUtils::Mesh::Connectivity nodesToElems);
      ///
      /// This function writes the output grid file
      /// 
      virtual int WriteOutput();
    };
    

#ifdef _GRIDCONVERSION_PARALLEL_
    ///
    /// Implementation of the basic parts of the parallel program example.
    ///
    /// This object encapsulates the example parallel program.  It inherits
    /// from the IRAD::Global::Program type. 
    /// 
    /// The program itself calculates \f$\pi\f$ by integrating 
    /// \f$\frac{d}{dx}(4\,\arctan(x))\f$ over the interval [0,1]:
    ///
    /// \f$\pi = \int_{0}^{1}\frac{4}{(1+x^2)}{dx}\f$
    ///
    /// This method is MPI parallel and uses domain decomposition to 
    /// calculate the total integral. On _N_ processors, the 
    /// _n_ partitions of the domain, [0,1], are partitioned as follows:
    ///
    /// First, the domain is broken into _n_ equal pieces of size _h_:
    ///
    /// \f$h = \frac{1}{n}\f$
    ///
    /// Each processor gets a roughly equal number of pieces, \f$n_i\f$. 
    /// The domain on each processor,\f$\mathcal{D}_i\f$, is defined as:
    ///
    /// \f$\mathcal{D}_i \equiv [x_i,x_{i+1}],\;\forall\,i \in [0,N-1]\f$
    ///
    /// With \f$x_0 = 0\f$ and \f$x_i = x_{i-1} + hn_{i-1},\;\forall\, i > 0\f$, 
    /// each processor performs the integration on its own portion of the domain:
    ///
    /// \f$\mathcal{I}_i = \int_{\mathcal{D}_i}\frac{4}{(1+x^2)}{dx}\f$ 
    /// 
    /// The approximate values for the integrals on each processor,
    /// \f$\bar{\mathcal{I}_i}\f$, are determined by numerical quadrature
    /// using two different methods:\n
    ///
    /// 1. GridConversion::TrapezoidQuadrature 
    /// 2. GridConversion::MidPointQuadrature
    ///
    /// The results of these quadrature methods from each processors 
    /// are then summed to recover a total numerical 
    /// estimate for \f$\pi\f$:
    ///
    /// \f$\bar{\pi} = \sum_{i=0}^{N-1}\bar{\mathcal{I}_i}\f$
    ///
    /// The approximate value, \f$\bar{\pi}\f$ and the error, 
    /// \f$|\pi-\bar{\pi}|\f$ are reported to stdout or, if specified, to
    /// a file.
    ///
    /// @note If the number of partitions, _n_, is less than the number
    /// of processors, _N_, then the program behaves as if it were serial,
    /// performing all calculations on processor 0.
    ///
    /// The usage for this program goes like:
    ///               
    ///           pepi Usage: 
    /// 
    ///           pepi [-h] [-v [arg] -o <filename> ] <input> 
    /// 
    ///           -h,--help
    ///              Prints this long version of help.
    ///        
    ///                  -v,--verbosity [arg]
    ///        
    ///                  -o,--output <filename>
    ///                     Specifies the name of the output file.
    ///
    ///           <input>
    ///                     Mode-dependent input: Input <filename> for serial example
    ///                     program, or <number of divisions> for parallel example.
    ///
    class ParallelProgram : public ParallelProgramType
    {
    protected:
      /// Name of file for output.
      std::string output_name;
      /// Verbosity level.
      int verblevel;           
      /// Output file stream object for output.
      std::ofstream Ouf;      
      /// Number of partitions for domain.
      int ndiv;
    public:
      /// 
      /// Default constructor.
      ///
      ParallelProgram() :
        ParallelProgramType()
      {verblevel = 0;};
      /// 
      /// Constructor designed to take all commandline args.
      ///
      ParallelProgram(int nargs,char **args) :
        ParallelProgramType(nargs,args)
      {verblevel = 0;};
      /// 
      /// Constructor designed to take all commandline args.
      ///
      ParallelProgram(ComLineType &comline,CommType &incomm) 
      {
        this->_command_line.Copy(comline);
        this->Init(_command_line.ProgramName(),incomm);
        verblevel = 0;
      };
      ///
      /// Populates native data structures from commandline.
      ///
      virtual int Initialize()
      {
        int retval = ParallelProgramType::Initialize();
        // Check if the user just entered -h or --help
        if(!_command_line.GetOption("help").empty()){
          std::ostringstream Ostr;
          Ostr << _command_line.LongUsage() << std::endl;
          StdOut(Ostr.str());
          return(-1);
        }
        if(retval){ // then there were errors on the command line
          std::ostringstream Ostr;
          Ostr << _command_line.ErrorReport() << std::endl
               << std::endl << _command_line.ShortUsage() << std::endl;
          ErrOut(Ostr.str());
          return(retval);
        }

        // Everything is OK up to here, so check for an output file.
        output_name = _command_line.GetOption("output");

        // Set number of divisions to first argument
        std::vector<std::string> args(_command_line.GetArgs());
        std::istringstream Istr(args[0]);
        Istr >> ndiv;

        // Check the verbosity level
        std::string sverb(_command_line.GetOption("verbosity"));
        if(sverb.empty() || sverb == ".true.")
          verblevel = 1;
        else {
          std::istringstream Vin(sverb);
          Vin >> verblevel;
          if(verblevel < 0)
            verblevel = 1;
        }

        // If high verbosity, stick a configuration blurb to stdout.
        if(verblevel > 1){
          std::ostringstream Ostr;
          Ostr << "Configuration:" << std::endl
               << "verbosity           = " << verblevel << std::endl
               << "number of divisions = " << ndiv << std::endl;
          if(!output_name.empty())
            Ostr << "output file         = " << output_name << std::endl;
          Ostr << std::endl;
          StdOut(Ostr.str());
        } 

        return(0);
      };
      ///
      /// Returns verbosity level.
      ///
      int VerbLevel() const { return verblevel;};
      ///
      /// Destructor.
      ///
      virtual ~ParallelProgram() {
        if(Ouf){
          Ouf.close();
          SetOutStream(std::cout);
        }
      };
      ///
      /// This function implements the main function executed by the program.
      /// 
      virtual int Run();
    };
#endif    


    ///
    /// Unified driver for the example programs
    ///
    /// This driver is exactly the same for both 
    /// the serial and parallel example programs.
    /// It is a direct replacement for main and
    /// is responsible for creating the "Program object"
    /// and calling its Initialization, Run, and 
    /// Finalize methods.
    ///
    template<typename ProgramType>
    int Driver(int argc,char *argv[])
    {
      ProgramType MyProgram(argc,argv);
      int retval = MyProgram.Initialize();
      if(retval){
        if(retval > 0){ // then there was an error
          MyProgram.ErrOut("Error in Program Initialization.\n");
          return(1);
        }
        else return(0); // no error, but the program should not continue
      }
      if(MyProgram.VerbLevel() > 1)
        MyProgram.StdOut("Calling RUN.\n");
      if(MyProgram.Run()){
        MyProgram.ErrOut("Error in Program RUN method.\n");
        return(1);
      }
      if(MyProgram.VerbLevel() > 1){
        std::ostringstream Ostr;
        MyProgram.Report(Ostr);
        MyProgram.StdOut(Ostr.str());
        MyProgram.StdOut("Calling Finalize.\n");
      }
      if(MyProgram.Finalize()){
        MyProgram.ErrOut("Error in Program Finalization.\n");
        return(1);
      }
      if(MyProgram.VerbLevel() > 1)
        MyProgram.StdOut("All done.\n");
      return(0);
    }


    ///
    /// Convenience definition for serial program type
    ///
    typedef GridConversion::DriverProgram::SerialProgram SEProgramType;
#ifdef _GRIDCONVERSION_PARALLEL_
    /// Support function for parallel example. 
    ///
    /// @param x The x-coordinate.
    /// @returns The value of the function at x
    ///
    /// This function returns the derivative of
    /// (4 arctan(x)) at x.  Specifically, it implements:
    ///
    /// \f$f(x) = \frac{d}{dx}({4}\arctan(x)) = \frac{4}{(1+x^2)}\f$.
    /// 
    inline double f(double x)
    {
      return (4.0 / (1.0 + x*x));
    }
    ///
    /// Convenience definition for parallel program type.
    ///
    typedef GridConversion::DriverProgram::ParallelProgram PEProgramType;
#endif
  };
};
#endif
