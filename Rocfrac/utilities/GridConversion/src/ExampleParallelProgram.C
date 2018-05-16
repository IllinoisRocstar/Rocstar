/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Example parallel program.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 

#include "ExampleProgram.H"
#include "ExampleHeader.H"
#include <cmath>
#include <iomanip>

namespace GridConversion {

  namespace ExampleProgram {

    int ParallelProgram::Run()
    {
      // FunctionEntry("NAME"): Updates the user-defined stack and 
      // the program profiles with timing information.  The placement
      // of FunctionEntry and FunctionExit calls is at the developer's
      // discretion.  
      FunctionEntry("Run");

      int myid  = _communicator.Rank();
      int nproc = _communicator.Size();
      bool use_file = !output_name.empty();

      // Put out a quick blurb about number of procs just
      // to give the user confidence we are actually running
      // in parallel.  Note that when using the "StdOut" method
      // that it automatically does output only on rank 0. If 
      // the developer wants asyncrhonous output, then they 
      // have to revert to using standard streams.
      std::ostringstream RepStr;
      if(verblevel > 1)
        RepStr << "Running on " << nproc << " processors." << std::endl;
      StdOut(RepStr.str());
      _communicator.Barrier();
      if(verblevel > 1)
        StdOut("All procesors ready.\n");

      // Open the specified output file for writing on rank 0
      if(use_file && !myid){
        Ouf.open(output_name.c_str(),std::ios::app);
        if(!Ouf){
          // If the output file failed to open, notify
          // to error stream and return non-zero
          std::ostringstream Ostr;
          Ostr << "Error: Unable to open output file, " << output_name << ".";
          ErrOut(Ostr.str());
          // In parallel, we don't return right away since 
          // this part of the code is only done on proc 0.
          // Instead, the error value is set in the communicator
          // which will indicate to all processors that there 
          // has been some error.
          _communicator.SetExit(1);
        }
      }
      // Check to see if an error condition was set 
      // in the file open block above.  If so, then
      // return with an error code.
      if(_communicator.Check()){
          // don't forget to tell the profiler/stacker the
          // function is exiting.
          FunctionExit("Run");
          return(1);
      }
        
      // Correct value of PI to several places
      double PIVAL = 3.141592653589793238462643383;

      // All processors should already have this as it 
      // came from the command line - but just in case,
      // we broadcast it here to make sure.
      FunctionEntry("Broadcast");
      _communicator.BroadCast(ndiv,0);
      FunctionExit("Broadcast");
      
      // This block partitions the 
      // domain [0,1] with ndiv 
      // intervals among nproc processors
      int nper = ndiv/nproc;
      int nvol = nper*nproc;
      int leftover = ndiv - nvol;
      int myn = nper;
      if(myid < leftover)
        myn++;
      int mystart = 0;
      for(int i = 0; i < myid;i++){
        mystart += nper;
        if(i < leftover)
          mystart++;
      }   
      // If there are less divisions than processors, then
      // just forget about domain decomposition and do everything
      // on processor 0.
      if(nper == 0){
        if(myid == 0) myn = ndiv;
        else myn = 0;
      }
        
      // Use the results from domain decomposition to 
      // set the domain [a,b] for this processor and get the 
      // stepsize.
      double h   = 1.0/(static_cast<double>(ndiv));
      double a = h*mystart;
      double b = a + h*myn;

      // Integrate f on this processor's subdomain using trapezoid quadrature
      FunctionEntry("TrapezoidMethod");
      double my_tpi = 0.0;
      if(myn > 0){ // if there are points in this processor's domain
        try {
          FunctionEntry("TrapezoidQuadrature");
          my_tpi = GridConversion::TrapezoidQuadrature(f,a,b,myn);
          FunctionExit("TrapezoidQuadrature");
        } catch (...){
          StdOut("Numerical limits of GridConversion::TrapezoidQuadrature exceeded.");
          my_tpi = 0.0;
        }
      }
      double tpi = 0.0;
      // Sum up the contributions from each processor and store the results on
      // processor 0 in "tpi".
      FunctionEntry("Reduce"); // time the communication
      _communicator.Reduce(my_tpi, tpi,IRAD::Comm::DTDOUBLE, IRAD::Comm::SUMOP, 0);
      FunctionExit("Reduce");  // end of the communication
      FunctionExit("TrapezoidMethod");

      // Integrate f on this processor's subdomain using midpoint quadrature
      FunctionEntry("MidPointMethod");
      double my_mppi = 0.0;
      if(myn > 0) {  // if there are points in this processor's domain
        try{
          FunctionEntry("MidPointQuadrature");
          my_mppi = GridConversion::MidPointQuadrature(f,a,b,myn);
          FunctionExit("MidPointQuadrature");
        } catch (...){
          StdOut("Numerical limits of GridConversion::MidPointQuadrature exceeded.");
          my_mppi = 0.0;
        }
      }
      double mppi = 0.0;
      // Sum up the contributions from each processor and store the results on
      // processor 0 in "mppi".
      FunctionEntry("Reduce"); // time the communication
      _communicator.Reduce(my_mppi, mppi,IRAD::Comm::DTDOUBLE, IRAD::Comm::SUMOP, 0);
      FunctionExit("Reduce");  // end of the communication
      FunctionExit("MidPointMethod");
      
      // Report the results to stdout or to file if one was specified.
      FunctionEntry("IO");
      if (!myid) {
        if(use_file){
          Ouf << ndiv << " " << std::setprecision(16) << tpi << " " 
              << std::setprecision(4)  << std::fabs(tpi-PIVAL) << " "
              << std::setprecision(16) << mppi << " " 
              << std::setprecision(4)  << std::fabs(mppi-PIVAL) << " "
              << std::endl;
          Ouf.close();
        }
        else if(verblevel) {
          std::ostringstream Ostr;
          Ostr << "With " << ndiv << " divisions, PI was calculated:" << std::endl
               << "MidPointQuadrature:  "
               << std::setprecision(16) << mppi << "\t\t" 
               << std::setprecision(4) << std::fabs(mppi - PIVAL) << std::endl
               << "TrapezoidQuadrature: " 
               << std::setprecision(16) << tpi << "\t\t"
               << std::setprecision(4)  << std::fabs(tpi - PIVAL) << std::endl;
          StdOut(Ostr.str());
        }
      }
      FunctionExit("IO");
      //
      // ---------- Program End ----------------  
      
      // Update the stacker/profiler that we are exiting 
      // this function.
      FunctionExit("Run");
      // return 0 for success
      return(0);
    };
  };
};
   
