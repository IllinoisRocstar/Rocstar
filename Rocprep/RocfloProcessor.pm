#The RocfloProcessor package inherits from ModuleProcessor.pm. 

package RocfloProcessor;

require "ModuleProcessor.pm";
@RocfloProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor

use strict;

my (%NDAfiles);

#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocfloProcessor";
}

#************************************************************************************
#************************************************************************************
sub getRocstarModuleName {
	return "Rocflo";
}

#************************************************************************************
#************************************************************************************
sub checkNDA {
    my ($props,$log, $badFileArray, $numBadFiles, $filename);
    my ($datadir, $griddir);
    $props = ModuleProcessor::getProps();
    $log   = ModuleProcessor::getLog  ();

    $datadir = $props->getSOURCEDIR()."Rocflo/".$props->getDATADIR("ROCFLO")."/";
    $griddir = $props->getSOURCEDIR()."Rocflo/".$props->getGRIDDIR("ROCFLO")."/";

    #check to see if the Data and Grid directories exist before we make the file list
    #if not, then return an empty array of filenames.  

    if (!ModuleProcessor::checkDataGridExists($datadir, $griddir)) {
	return(1);
    }
    getNDAFiles($datadir, $griddir);
    my (@tempArray) = keys(%NDAfiles);
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
    my ($makeflo, $rfloprep, $plagprep, $binDir, $axis, $axisflag, $plagbc);
    $log = ModuleProcessor::getLog();
    $props=ModuleProcessor::getProps();

    $targetDir = $props->getTARGETDIR();
    $binDir = $props->getBINDIR();
    $caseName = getCaseName($targetDir."Rocflo/");
    $currentDir = `pwd`;
    chomp($currentDir);
    $numProcs = $props->getNumProcs();
    if ($props->propExists("SPLITAXIS")) {
       $axis = $props->getSplitAxis();
       $axisflag = "-splitaxis $axis";
    }

    chdir($targetDir."Rocflo/Modin");

    $makeflo = $binDir."makeflo $axisflag $caseName"."-PLOT3D.grd "."$numProcs $caseName.top $caseName.grda > $targetDir/makeflo.log 2>&1";

#DEBUG LINE
    $log->printMessage($makeflo);

    $errorCode = system($makeflo);

    if ($errorCode) {
	    $log->processErrorCode(11, $makeflo);
	    return(1);
    }

    if( $props->propExists("PLAG") ) {
      
      $plagbc   = $props->getPLAG();
      $plagprep = $binDir."plagprep -p $caseName -bc $plagbc > $targetDir/plagprep.log 2>&1";

      $log->printMessage( $plagprep );
      $errorCode = system($plagprep);

      if( $errorCode ) {
        $log->processErrorCode(11, $plagprep);
        return(1);
      }


    }

    $rfloprep = $binDir."rfloprep $caseName 1 0 1 > $targetDir/rfloprep.log 2>&1";
    $errorCode = system($rfloprep);

    if ($errorCode) {
	$log->processErrorCode(11, $rfloprep);
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
    $props=ModuleProcessor::getProps();

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
#  Subroutine that develops a list of files and directories that should be in the
#  NDA given user input of a specific Data and Grid combination. 
#************************************************************************************
sub getNDAFiles {
    my ($datadir,$griddir) = @_;
    my ($props, $casename);
    $props = ModuleProcessor::getProps();

    $casename = getCaseName($datadir);

    if ($casename eq "NOT FOUND") {return;}

    $NDAfiles{$datadir."RocfloControl.txt"} = $props->getTARGETDIR()."Rocflo/" ;
    $NDAfiles{$datadir.$casename.'.inp'} =  $props->getTARGETDIR()."Rocflo/Modin/";
    $NDAfiles{$datadir.$casename.'.bc'} =  $props->getTARGETDIR()."Rocflo/Modin/";
    $NDAfiles{$griddir.$casename.'-PLOT3D.bcmp'} =  $props->getTARGETDIR()."Rocflo/Modin/";
    $NDAfiles{$griddir.$casename.'-PLOT3D.inp'} =  $props->getTARGETDIR()."Rocflo/Modin/";
    $NDAfiles{$griddir.$casename.'-PLOT3D.grd'} =  $props->getTARGETDIR()."Rocflo/Modin/";
}

#************************************************************************************
#  Subroutine that develops a list of files and directories that should be in the
#  dataset so that it is considered ready to run. 
#************************************************************************************
sub getRuntimeFiles {
    my ($props, $casename);
    my ($targetDir, @runtimeFiles);
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR()."Rocflo/";
    $casename = getCaseName($targetDir);

    if ($casename eq "NOT FOUND") {return;}

    push(@runtimeFiles, $targetDir."RocfloControl.txt");
    push(@runtimeFiles, $targetDir."Modin/".$casename.".inp");
    push(@runtimeFiles, $targetDir."Modin/".$casename.".bc");
    push(@runtimeFiles, $targetDir."Modin/".$casename.".top");

    #parse volume Rocin Control file
    push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/fluid_in_00.000000.txt"));

    return(@runtimeFiles);
}


#************************************************************************************
#
#************************************************************************************
sub getCaseName {
    my ($dataDir) = shift;
    my ($controlFile, $log, $caseName);
    $controlFile = "RocfloControl.txt";
    $controlFile = $dataDir.$controlFile;

    $log = ModuleProcessor::getLog();

    if (!open(CONTROLFILE, $controlFile))
        {
	    $log->processErrorCode(23, $controlFile);
	    $log->printMessage("$controlFile must exist to supply casename for run files.");
	    return("NOT FOUND");
	}

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
    return($props->processModule("ROCFLO"));
}

1;
