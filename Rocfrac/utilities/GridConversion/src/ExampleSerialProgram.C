/// 
/// @file
/// @ingroup gridconversion_group
/// @brief Example serial program.
/// @author Mike Campbell (mtcampbe@illinois.edu)
/// @date 
/// 

#include "ExampleProgram.H"


namespace GridConversion {

  namespace ExampleProgram {

    int SerialProgram::Run()
    {
      // FunctionEntry("NAME"): Updates the user-defined stack and 
      // the program profiles with timing information. The placement
      // of FunctionEntry and FunctionExit calls is at the developer's
      // discretion.  
      FunctionEntry("Run");

      // ---------- The Program -----------------
      // This program just "copies" the input file
      // to the output file.
      //

      // Open the specified input file for reading
      Inf.open(input_name.c_str());
      if(!Inf){
        // If the input file failed to open, notify to
        // the error stream and return non-zero
        std::ostringstream Ostr;
        Ostr << "Error: Could not open input file, '" 
             << input_name << "'.\n";
        ErrOut(Ostr.str());
        // don't forget to tell the profiler/stacker the
        // function is exiting.
        FunctionExit("Run");
        return(1);
      }

      // Open the specified output file for writing
      bool use_outfile = false;
      if(!output_name.empty()){
        use_outfile = true;
        Ouf.open(output_name.c_str());
        if(!Ouf){
          // If the output file failed to open, notify
          // to error stream and return non-zero
          std::ostringstream Ostr;
          Ostr << "Error: Unable to open output file, " << output_name << ".";
          ErrOut(Ostr.str());
          // don't forget to tell the profiler/stacker the
          // function is exiting.
          FunctionExit("Run");
          return(1);
        }
      }

      // Read lines from the input file and repeat them 
      // to the output file.
      std::string line;
      while(std::getline(Inf,line)){
        if(use_outfile)
          Ouf << line << std::endl;
        else
          StdOut(line+"\n");
      }
      // Close both files
      Ouf.close();
      Inf.close();
      
      //
      // ---------- Program End -----------------


      // Update the stacker/profiler that we are exiting 
      // this function.
      FunctionExit("Run");
      // return 0 for success
      return(0);
    };
  };
};
