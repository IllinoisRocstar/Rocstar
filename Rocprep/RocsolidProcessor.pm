#The RocsolidProcessor package inherits from ModuleProcessor.pm. 

package RocsolidProcessor;

require "ModuleProcessor.pm";
@RocsolidProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor

use strict;

my (%NDAfiles);

#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocsolidProcessor";
}

#************************************************************************************
#************************************************************************************
sub getRocstarModuleName {
	return "Rocsolid";
}

#************************************************************************************
#************************************************************************************
sub checkNDA {
    my ($props, $log, $badFileArray, $numBadFiles, $filename);
    my ($datadir, $griddir, $tempName, $i);
    $props = ModuleProcessor::getProps();
    $log   = ModuleProcessor::getLog  ();

    $datadir = $props->getSOURCEDIR()."Rocsolid/".$props->getDATADIR("ROCSOLID")."/";
    $griddir = $props->getSOURCEDIR()."Rocsolid/".$props->getGRIDDIR("ROCSOLID")."/";

    #check to see if the Data and Grid directories exist before we make the file list
    #if not, then return an empty array of filenames.  

    if (!ModuleProcessor::checkDataGridExists($datadir, $griddir)) {
	return(1);
    }
    getNDAFiles($datadir, $griddir);
    my (@tempArray) = keys(%NDAfiles);
    $badFileArray = ModuleProcessor::checkFilesExist(\@tempArray);
    $numBadFiles = @$badFileArray;

    #The following is the expected behavior.  Since there can be one of two filenames
    #for the grid file, and we're checking for both of them, it is expected that one
    #of the two won't exist.  If numBadFiles = 2 or 3, then neither were found, and an
    #error should occur. 

    if ($numBadFiles == 1) {
	delete $NDAfiles{@$badFileArray[0]};
	$numBadFiles--;
    }
    
    if ($numBadFiles) {
	foreach $filename (@$badFileArray) {
	    $log->processErrorCode(23, $filename);
	}
	return(1);
    }

    return(0);
}

#************************************************************************************
#
#************************************************************************************
sub extract {
    my ($props);
    $props = ModuleProcessor::getProps();

    return(ModuleProcessor::performExtract(\%NDAfiles));
}

#************************************************************************************
#************************************************************************************
sub preProcess {
    my ($log, $props, $targetDir,$currentDir, $errorCode, $caseName, $numProcs);
    my ($rsolidprep, $binDir, $workingDir, @lines);
    $log = ModuleProcessor::getLog();
    $props=ModuleProcessor::getProps();

    $targetDir = $props->getTARGETDIR();
    $workingDir = $targetDir."Rocsolid/";
    $binDir = $props->getBINDIR();
    $caseName = getCaseName($workingDir.'Modin/');
    $currentDir = `pwd`;
    chomp($currentDir);
    $numProcs = $props->getNumProcs();
    
    chdir($workingDir);
    # open the RocsolidPrep input file and read its contents into an array
    if (!open(INPFILE, "Modin/RocsolidPrep.inp"))
    {
      $log->processErrorCode(23, "RocsolidPrep.inp");
      $log->printMessage("RocsolidPrep.inp is needed to specify number of partitions.");
      return(1);
    }

    @lines = <INPFILE>;
    close(INPFILE);
    if (@lines < 6)
    {
      $log->processErrorCode(24, "RocsolidPrep.inp");
      $log->printMessage("RocsolidPrep.inp truncated, there is no input for number of partitions.");
      return(1);
    }

    # dynamically edit and replace the number of partitions at line #6
    $lines[5] = "$numProcs            ! NumberOfPartitions\n";
    
    # reopen the file clobbering it, and write out the modified array contents
    if (!open(INPFILE, ">Modin/RocsolidPrep.inp"))
    {
      $log->processErrorCode(23, "RocsolidPrep.inp");
      $log->printMessage("RocsolidPrep.inp cannot be overwritten to modify number of partitions.");
      return(1);
    }
    print INPFILE @lines;
    close(INPFILE);

    $rsolidprep = $binDir."rsolidprep > $targetDir/rsolidprep.log 2>&1";

    $errorCode = system($rsolidprep);

    if ($errorCode) {
	$log->processErrorCode(11, $rsolidprep);
	return(1);
    }
    chdir($currentDir);
    
    return(0);
}

#************************************************************************************
#************************************************************************************
sub checkResults {
    my ($log, $props, $badFileArray, $numBadFiles, $filename);
    $log = ModuleProcessor::getLog();
    $props = ModuleProcessor::getProps();

    my (@tempArray) = getRuntimeFiles();
    $badFileArray = ModuleProcessor::checkFilesExist(\@tempArray);
    $numBadFiles = @$badFileArray;

    if ($numBadFiles) {
	foreach $filename (@$badFileArray) {
	    $log->processErrorCode(23, $filename);
	}
	return(1);
    }

    return(0);
}

#************************************************************************************
#  Subroutine that develops a list of files that should be in the NDA given user input
#  of a specific Data and Grid combination. 
#************************************************************************************
sub getNDAFiles {
    my ($datadir,$griddir) = @_;
    my ($props, $log, $casename);
    $log = ModuleProcessor::getLog();
    $props = ModuleProcessor::getProps();

    $casename = getCaseName($griddir);

    if ($casename eq "NOT FOUND") {return;}

    $NDAfiles{$datadir."RocsolidControl.txt"} = $props->getTARGETDIR()."Rocsolid/" ;
    $NDAfiles{$griddir."RocsolidPrep.inp"} =  $props->getTARGETDIR()."Rocsolid/Modin/";
    $NDAfiles{$griddir.$casename.'.trugrd'} =  $props->getTARGETDIR()."Rocsolid/Modin/";
}

#************************************************************************************
#  Subroutine that develops a list of files and directories that should be in the
#  dataset so that it is considered ready to run. 
#************************************************************************************
sub getRuntimeFiles {
    my ($props);
    my ($targetDir, @runtimeFiles);
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR()."Rocsolid/";

    push(@runtimeFiles, $targetDir."RocsolidControl.txt");
    push(@runtimeFiles, $targetDir."Modin/RocsolidPrep.inp");

# FIX THIS!!!!   Once we understand what the expected .hdf files are from Rocsolid
# NOTE:  the RocsolidPrep.inp contains a prefix for all hdf files that will be needed here

    #parse volume Rocin Control file
#   push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/solid_in_00.000000.txt"));

    #parse surface Rocin Control file
#   push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/isolid_in_00.000000.txt"));

    return(@runtimeFiles);
}


#************************************************************************************
#
#************************************************************************************
sub getCaseName {
    my ($dataDir) = shift;
    my ($controlFile, $log, $title, $caseName);
    $controlFile = "RocsolidPrep.inp";
    $controlFile = $dataDir.$controlFile;

    $log = ModuleProcessor::getLog();

    if (!open(CONTROLFILE, $controlFile))
        {
	    $log->processErrorCode(23, $controlFile);
	    $log->printMessage("$controlFile must exist to supply casename for grid file.");
	    return("NOT FOUND");
	}

    $title = <CONTROLFILE>;
    $caseName = <CONTROLFILE>;

    close (CONTROLFILE);

    if (($caseName eq undef)||($caseName eq "")) {
	    $log->processErrorCode(24, "In file: ".$controlFile.", can't read case name.");
	    return("NOT FOUND");
    }
    
    chomp($caseName);
    
    #the following strips spaces from the beginning and end of the casename. 
    for ($caseName) {
	s/^\s+//;
	s/\s+$//;
    }  


    return ($caseName);
}
    

#************************************************************************************
#   Returns Boolean whether the module should be processed or not. 
#************************************************************************************
sub isActive {
    my ($props);
    $props = ModuleProcessor::getProps();
    return($props->processModule("ROCSOLID"));
}

1;
