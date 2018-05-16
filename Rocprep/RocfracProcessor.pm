#The RocfracProcessor package inherits from ModuleProcessor.pm. 

package RocfracProcessor;

require "ModuleProcessor.pm";
@RocfracProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor

use strict;

my (%NDAfiles);

#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocfracProcessor";
}

#************************************************************************************
#************************************************************************************
sub getRocstarModuleName {
	return "Rocfrac";
}

#************************************************************************************
#************************************************************************************
sub checkNDA {
    my ($props, $log, $badFileArray, $numBadFiles, $filename);
    my ($datadir, $griddir, $tempName, $i);
    my($casename, $form1Name, $form2Name, $suffix, $file);
    $props = ModuleProcessor::getProps();
    $log   = ModuleProcessor::getLog  ();

    $datadir = $props->getSOURCEDIR()."Rocfrac/".$props->getDATADIR("ROCFRAC")."/";
    $griddir = $props->getSOURCEDIR()."Rocfrac/".$props->getGRIDDIR("ROCFRAC")."/";

    #check to see if the Data and Grid directories exist before we make the file list
    #if not, then return an empty array of filenames.  

    if (!ModuleProcessor::checkDataGridExists($datadir, $griddir)) {
	return(1);
    }
    getNDAFiles($datadir, $griddir);
    my (@tempArray) = keys(%NDAfiles);
    $badFileArray = ModuleProcessor::checkFilesExist(\@tempArray);
    $numBadFiles = @$badFileArray;

    #The following is the expected behavior.  Since there can be one of four filenames
    #for the grid file, and we're checking for all of them, it is expected that
    #three of the four won't exist.  If numBadFiles = 4 or 5, then none were found, and an
    #error should occur. 

    if ($numBadFiles == 3) {
	delete $NDAfiles{@$badFileArray[0]};       
        delete $NDAfiles{@$badFileArray[1]};
        delete $NDAfiles{@$badFileArray[2]};

	$numBadFiles=$numBadFiles-3;
    }
    
    if ($numBadFiles) {
	foreach $filename (@$badFileArray) {
	    $log->processErrorCode(23, $filename);
	}
	return(1);
    }

    #Now, here comes the hack. Since Rocfrac allows it's patran grid files
    #to come in multiple pieces, we have to find them all and put them on
    #this $NDAfiles list. Otherwise, they won't be copied.  Rocfrac
    #unfortunately can't detect this problem, and will merrily go on and
    #preprocess whatever is there, making an incomplete grid, which will
    #then fail to surfdive. Worse, if this is Rocfrac alone, it will
    #complete the preprocessing, and act like everything is OK. The
    #filenames look like <casename>.1.out or .1.pat, and then will have 
    #.2.out, .3.out, .4.out, etc...with no way to know how many except to
    #search. Fun. 

    $casename = getCaseName($datadir);
    if ($casename eq "NOT FOUND") {
        print "COULD NOT FIND CASE NAME IN RocfracProcessor.checkNDA!";
	return;
    }

    $form1Name = $griddir.$casename.'.1.pat';
    $form2Name = $griddir.$casename.'.1.out';
    $suffix = "";

    foreach $file (keys(%NDAfiles)) {
        #look for either of the form1Name or form2Name.
        #should never actually see both... 
        if ($file =~ m/$form1Name/) {$suffix = ".pat";}
	if ($file =~ m/$form2Name/) {$suffix = ".out";} 
    }

    #Non-empty strings are true in Perl...
    if ($suffix) {
       #if suffix is non-empty, assume we have to find more parts of the
       #grid. Otherwise, assume there's a single part with no .1 or .2, etc.
   
       $log->printMessage("Ready to look for more grid parts. Multipart Rocfrac Grid format found.");

       addGridPartsToNDAFiles($griddir.$casename, $suffix);
    }

    return(0);
}



#************************************************************************************
#
#************************************************************************************
sub addGridPartsToNDAFiles {
    my ($prefix,$suffix) = @_;
    my ($testFile, $fileNumber,$props, $log);

   
    $props = ModuleProcessor::getProps();
    $log = ModuleProcessor::getLog();
    $fileNumber=2;

    $log->printMessage("prefix = $prefix, and suffix = $suffix");    

    while(1) {
        $testFile = $prefix.'.'.$fileNumber.$suffix;
        if (-e $testFile) {
	    $log->printMessage("Found Grid filename: $testFile");
            $NDAfiles{$testFile} = $props->getTARGETDIR()."Rocfrac/Modin/";
	    $fileNumber++;
        }
        else {
            $log->printMessage("Assume that all parts of the Rocfrac Grid have been found. There is no way to validate this.\n");
            $fileNumber++;
            last;
        }
        #Safety for coding screwups: 
        if ($fileNumber > 50) {
           $log->printMessage("Coding error in RocfracProcessor.addGridPartsToNDAFiles:more than 50 files tried!");
           last;
        }
    }
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
    my ($rfracprep, $binDir, $workingDir, $units, $unitflag);
    $log = ModuleProcessor::getLog();
    $props=ModuleProcessor::getProps();

    $targetDir = $props->getTARGETDIR();
    $workingDir = $targetDir."Rocfrac/";
    $binDir = $props->getBINDIR();
    $caseName = getCaseName($workingDir);
    $currentDir = `pwd`;
    chomp($currentDir);
    $numProcs = $props->getNumProcs();
    if ($props->propExists("ROCFRACUNITS")) {
       $units = $props->getRocfracUnits();
       $unitflag = "-un $units";
    }
    
    chdir($workingDir);
    
    $rfracprep = $binDir."rfracprep -np $numProcs $unitflag > $targetDir/rfracprep.log 2>&1";

    $errorCode = system($rfracprep);

    if ($errorCode) {
	$log->processErrorCode(11, $rfracprep);
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

    $casename = getCaseName($datadir);

    if ($casename eq "NOT FOUND") {return;}

    #Note: The grid file will one of the following four, but we don't know which
    #Also note that the ones with .1.pat or .1.out will also be followed by
    #some unknown number of .2.out, .3.out, etc. files. We won't check
    #those.

    $NDAfiles{$datadir."RocfracControl.txt"} = $props->getTARGETDIR()."Rocfrac/" ;
    $NDAfiles{$griddir.$casename.'.pat'} =  $props->getTARGETDIR()."Rocfrac/Modin/";
    $NDAfiles{$griddir.$casename.'.out'} =  $props->getTARGETDIR()."Rocfrac/Modin/";
    $NDAfiles{$griddir.$casename.'.1.pat'} = $props->getTARGETDIR()."Rocfrac/Modin/";
    $NDAfiles{$griddir.$casename.'.1.out'} = $props->getTARGETDIR()."Rocfrac/Modin/";
}

#************************************************************************************
#  Subroutine that develops a list of files and directories that should be in the
#  dataset so that it is considered ready to run. 
#************************************************************************************
sub getRuntimeFiles {
    my ($props);
    my ($targetDir, @runtimeFiles);
    $props = ModuleProcessor::getProps();
    $targetDir = $props->getTARGETDIR()."Rocfrac/";

    push(@runtimeFiles, $targetDir."RocfracControl.txt");

    #parse volume Rocin Control file
    push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/solid_in_00.000000.txt"));

    #parse surface Rocin Control file
    push(@runtimeFiles,ModuleProcessor::getRocinFiles($targetDir."Rocin/isolid_in_00.000000.txt"));

    return(@runtimeFiles);
}


#************************************************************************************
#
#************************************************************************************
sub getCaseName {
    my ($dataDir) = shift;
    my ($controlFile, $log, $caseName);
    $controlFile = "RocfracControl.txt";
    $controlFile = $dataDir.$controlFile;

    $log = ModuleProcessor::getLog();

    if (!open(CONTROLFILE, $controlFile))
        {
	    $log->processErrorCode(23, $controlFile);
	    $log->printMessage("$controlFile must exist to supply casename for grid file.");
	    return("NOT FOUND");
	}

    while(<CONTROLFILE>) {
        chomp;
	if ($_ =~ /^\*PREFIX/) {
	    $caseName = <CONTROLFILE>;
	    last;
	}
    }

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
    return($props->processModule("ROCFRAC"));
}

1;
