#The RocfaceProcessor package inherits from ModuleProcessor.pm. 

package RocfaceProcessor;

require "ModuleProcessor.pm";
@RocfaceProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor
use strict;


#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocfaceProcessor";
}

#************************************************************************************
#************************************************************************************
sub preProcess {
    my ($log, $props);
    my ($fluidsmodule, $solidsmodule, $targetDir, $currentDir, $binDir);
    my ($error);
    my ($errorflag) = 0;
    my ($fluidsfile) = "Rocin/ifluid_in_00.000000.txt";
    my ($solidsfile) = "Rocin/isolid_in_00.000000.txt";

    $log = ModuleProcessor::getLog();
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR();
    $binDir = $props->getBINDIR();
    $currentDir = `pwd`;
    chomp($currentDir);
    chdir($targetDir);

    #Process all the necessary code combinations with Surfdiver (e.g., RocfloRocfrac). 
    foreach $fluidsmodule ("Rocflo", "Rocflu") {
       foreach $solidsmodule("Rocfrac", "Rocsolid") {
	   my ($combo) = $fluidsmodule.$solidsmodule;
	   my ($logfile) = 'surf_'.lc($combo).'.log';
	   $logfile =~ s/roc//g;
	   if ($props->processRocfaceCombo(uc($combo))){
	       my ($fluid) = $fluidsmodule."/".$fluidsfile;
	       my ($solid) = $solidsmodule."/".$solidsfile;
	       my ($commandline) = $binDir."surfdiver $fluid $solid Rocman/$combo/ > $logfile 2>&1";
	       $error = system($commandline);
	       if ($error) {
		   $errorflag++;
		   $log->processErrorCode(11, $commandline);
	       }
	   }	
	}
     }

    chdir($currentDir);
    return($errorflag);
}

#************************************************************************************
# Subroutine verifies that surfdiver executed on the correct combinations of fluid &
# solid interface meshes. Uses a very primative method to check if surfdiver was
# successful: list the resulting output files with a wildcard for each fluid and solid,
# assuming that if one or more files are present for each surfdiver succeeded.
#************************************************************************************
sub checkResults {
    my ($log, $props);
    my ($fluidsmodule, $solidsmodule, $targetDir);
    $log = ModuleProcessor::getLog();
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR()."Rocman/";

    #check all code combinations were successfully surfdived (e.g., RocfloRocfrac). 
    foreach $fluidsmodule ("Rocflo", "Rocflu") {
       foreach $solidsmodule("Rocfrac", "Rocsolid") {
	   my ($combo) = $fluidsmodule.$solidsmodule;
	   if ($props->processRocfaceCombo(uc($combo))){
	      unless (glob($targetDir.$combo."/ifluid*.hdf") or glob($targetDir.$combo."/ifluid*.cgns")) {
		 $log->processErrorCode(27, "$combo fluid interface missing?");
		 return(1);
	      }
	      unless (glob($targetDir.$combo."/isolid*.hdf") or glob($targetDir.$combo."/isolid*.cgns")) {
		 $log->processErrorCode(27, "$combo solid interface missing?");
		 return(1);
	      }
           }
       }
    }
    return(0);
}


1;
