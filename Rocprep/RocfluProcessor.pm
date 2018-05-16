#The RocfluProcessor package inherits from ModuleProcessor.pm. 

package RocfluProcessor;

require "ModuleProcessor.pm";
@RocfluProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor

use strict;

my (%NDAfiles);

#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocfluProcessor";
}

#************************************************************************************
#************************************************************************************
sub getRocstarModuleName {
	return "Rocflu";
}

#************************************************************************************
#************************************************************************************
sub checkNDA {
    my ($props, $log, $badFileArray, $numBadFiles, $filename);
    my ($datadir, $griddir);
    $props = ModuleProcessor::getProps();
    $log   = ModuleProcessor::getLog  ();

    $datadir = $props->getSOURCEDIR()."Rocflu/".$props->getDATADIR("ROCFLU")."/";
    $griddir = $props->getSOURCEDIR()."Rocflu/".$props->getGRIDDIR("ROCFLU")."/";

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
    my ($log, $props, $targetDir,$currentDir, $errorCode, $caseName, $numProcs, $numRegions);
    my ($rflumap, $rflupart, $rfluinit, $binDir);
    $log = ModuleProcessor::getLog();
    $props=ModuleProcessor::getProps();

    $targetDir = $props->getTARGETDIR();
    $binDir = $props->getBINDIR();
    $caseName = getCaseName($targetDir."Rocflu/");
    $currentDir = `pwd`;
    chomp($currentDir);
    $numProcs = $props->getNumProcs();
    if ($props->propExists("NUMREGIONS")) {
       $numRegions = $props->getNumRegions();
    }
    else {
       $numRegions = $numProcs;
    }

    chdir($targetDir);

    # since rflumap is interactive, this uses an obscure shell feature called
    # a "here document" to redirect shell input for multiple lines of responses
    #$rflumap = $binDir."rflumap $caseName 1 > rflumap1.log 2>&1 << ENDINPUT
    #			  1
    #			  $numRegions
    #			  $numProcs
    #			ENDINPUT";
    
    #print "rflumap -m 1 -v 1 -c $caseName -r $numRegions -p $numProcs > rflumap1.log 2>&1";
    $rflumap = $binDir."rflumap -m 1 -v 1 -c $caseName -r $numRegions -p $numProcs > rflumap1.log 2>&1";	

    $errorCode = system($rflumap);

    if ($errorCode) {
	$log->processErrorCode(11, $rflumap);
	return(1);
    }

    #5/2/05 Addition: MDB
    #broke up rfluprep call into rflupart and rfluinit per new MPI version
    #of rocflu code and preprocessors. 

    $rflupart = $binDir."rflupart -c $caseName -v 1 > rflupart.log 2>&1";
    $errorCode = system($rflupart);

    if ($errorCode) {
	$log->processErrorCode(11, $rflupart);
	return(1);
    } 

    $rfluinit = $binDir."rfluinit -c $caseName -v 1 > rfluinit.log 2>&1";
    #print "$rfluinit";
    $errorCode = system($rfluinit);

    if ($errorCode) {
	$log->processErrorCode(11, $rfluinit);
	return(1);
    }

    # since rflumap is interactive, this uses an obscure shell feature called
    # a "here document" to redirect shell input for multiple lines of responses
    #$rflumap = $binDir."rflumap $caseName 1 > rflumap2.log 2>&1 << ENDINPUT
    #			  2
    #			ENDINPUT";


    $rflumap = $binDir."rflumap -m 2 -v 1 -c $caseName > rflumap2.log 2>&1";

    $errorCode = system($rflumap);

    if ($errorCode) {
	$log->processErrorCode(11, $rflumap);
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
#  Subroutine that develops a list of files that should be in the NDA given user input
#  of a specific Data and Grid combination. 
#************************************************************************************
sub getNDAFiles {
    my ($datadir,$griddir) = @_;
    my ($props, $log, $casename);
    $log = ModuleProcessor::getLog();
    $props = ModuleProcessor::getProps();

    $casename = getCaseName($datadir);

    if ($casename eq "NOT FOUND") {return;}

    $NDAfiles{$datadir."RocfluControl.txt"} = $props->getTARGETDIR()."Rocflu/" ;
    $NDAfiles{$datadir.$casename.'.inp'} =  $props->getTARGETDIR()."Rocflu/Modin/";
    $NDAfiles{$datadir.$casename.'.bc'} =  $props->getTARGETDIR()."Rocflu/Modin/";
    $NDAfiles{$griddir.$casename.'.cgi'} =  $props->getTARGETDIR()."Rocflu/Modin/";
    $NDAfiles{$griddir.$casename.'-COBALT.inp'} =  $props->getTARGETDIR()."Rocflu/Modin/$casename.cgr";
    $NDAfiles{$griddir.$casename.'-COBALT.bc'} =  $props->getTARGETDIR()."Rocflu/Modin/";
}

#************************************************************************************
#  Subroutine that develops a list of files and directories that should be in the
#  dataset so that it is considered ready to run. 
#************************************************************************************
sub getRuntimeFiles {
    my ($props, $casename);
    my ($targetDir, @runtimeFiles);
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR()."Rocflu/";
    $casename = getCaseName($targetDir);

    if ($casename eq "NOT FOUND") {return;}

    push(@runtimeFiles, $targetDir."RocfluControl.txt");
    push(@runtimeFiles, $targetDir."Modin/".$casename.".inp");
    push(@runtimeFiles, $targetDir."Modin/".$casename.".bc");

    #TBD: There may be more files to be checked - we don't actually know what rfluprep will put out yet.

    #parse volume Rocin Control file
    push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/fluid_in_00.000000.txt"));

    #parse surface Rocin Control file
    push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/ifluid_in_00.000000.txt"));

    #TBD - when particles are supported, this pattern will probably need to be checked
    #parse particle Rocin Control file
    #push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/fluid_plag_in_00.000000.txt"));

    return(@runtimeFiles);
}


#************************************************************************************
#
#************************************************************************************
sub getCaseName {
    my ($dataDir) = shift;
    my ($controlFile, $log, $caseName);
    $controlFile = "RocfluControl.txt";
    $controlFile = $dataDir.$controlFile;

#TBD decide whether to grab casename from filename and make control file. 

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
    return($props->processModule("ROCFLU"));
}

1;
