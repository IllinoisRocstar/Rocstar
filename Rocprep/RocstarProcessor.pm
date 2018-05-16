#The RocstarProcessor package inherits from ModuleProcessor.pm. 

package RocstarProcessor;

require "ModuleProcessor.pm";
@RocstarProcessor::ISA = qw(ModuleProcessor); #Inherits from ModuleProcessor
use strict;

use RocfaceProcessor;

my (%NDAfiles);

my ($rface);
$rface=0;

my ($TRUE, $FALSE, $timestamp);
$TRUE=1;
$FALSE=0;

#************************************************************************************
#************************************************************************************
sub getModuleName {
	return "RocstarProcessor";
}


#************************************************************************************
#************************************************************************************
sub checkNDA {
    my ($thePackage, $processorObjs) = @_;
    my ($log, $props, $obj, @error);
    my ($badFileArray, $numBadFiles, $badFile, $baseDir);

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();

    $timestamp = localtime;
    $log->printMessage("\n$timestamp: Checking NDA files\n");

    foreach $obj (@$processorObjs) {
	if ($obj->isActive()) {
	    push(@error, $obj->checkNDA());
	    $log->printEndPhase("Check NDA Files", ref($obj));
        }
    }
    
    $baseDir = $props->getSOURCEDIR()."Rocstar/";

    getNDAFiles();
    my (@tempArray) = keys(%NDAfiles);
    $badFileArray = ModuleProcessor::checkFilesExist(\@tempArray);

    $numBadFiles = @$badFileArray;
    if ($numBadFiles) {
       push(@error, 23);
       foreach $badFile (@$badFileArray) {
          $log->processErrorCode(23, $badFile);
       }
    }

    # check Rocburn files
    # It is an error to have zero Rocburn models if the Rocburn flag is true. 
    # As of 9/16/04, the ROCBURN flag is forced to be true in Rocprep.pm. 
    if ($props->processModule("ROCBURN")) {
        if (!analyzeRocburn($baseDir)) {
            $log->processErrorCode(301, $baseDir);
            push(@error, 301);
        };
    }

    # check Rocmop files
    if ($props->processModule("ROCMOP")) {
        if (!analyzeRocburn($baseDir)) {
            $log->processErrorCode(405, $baseDir);
            push(@error, 405);
        };
    }

    $log->printSeparator();
    return (hasError(\@error));
}

#************************************************************************************
# This routine checks to see if any of the three Rocburn flavors are available.  That
# means that the appropriate directory and control file are both "there".  The routine
# is called during the extract phase, with $basedir = the NDA Rocstar directory.  It
# sets the source and extract path in $NDAfiles in this case.  During the checkResults
# phase, this routine is called with the Rocstar3.0 dataset directory as $baseDir. 
# In the checkResults phase, the $NDAfiles hash is still set, but is never used.  The
# goal is to get one or more of the three keywords set in the properties, so from then
# on we know which Rocburn models are available. 
#************************************************************************************
sub analyzeRocburn {
    my ($basedir) = @_;
    my ($log, $props, $badFileArray, $numBadFiles, $rbDir);
    my ($numFound, @fileList, $targetDir);
    $numFound = 0;

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();
    $targetDir = $props->getTARGETDIR();

    $rbDir = $basedir."RocburnAPN/";
    @fileList = ($rbDir, $rbDir."RocburnAPNControl.txt");
    $badFileArray = ModuleProcessor::checkFilesExist(\@fileList);
    $numBadFiles = @$badFileArray;
    if (!$numBadFiles) {
	$props->setKeyValuePair("ROCBURNAPN", $TRUE);
	$NDAfiles{$rbDir."RocburnAPNControl.txt"} = $targetDir."RocburnAPN/";
	$numFound++;
    }

    $rbDir = $basedir."RocburnPY/";
    @fileList = ($rbDir, $rbDir."RocburnPYControl.txt");
    $badFileArray = ModuleProcessor::checkFilesExist(\@fileList);
    $numBadFiles = @$badFileArray;
    if (!$numBadFiles) {
	$props->setKeyValuePair("ROCBURNPY", $TRUE);
	$NDAfiles{$rbDir."RocburnPYControl.txt"} = $targetDir."RocburnPY/";
	$numFound++;
    }

    $rbDir = $basedir."RocburnZN/";
    @fileList = ($rbDir, $rbDir."RocburnZNControl.txt");
    $badFileArray = ModuleProcessor::checkFilesExist(\@fileList);
    $numBadFiles = @$badFileArray;
    if (!$numBadFiles) {
	$props->setKeyValuePair("ROCBURNZN", $TRUE);
	$NDAfiles{$rbDir."RocburnZNControl.txt"} = $targetDir."RocburnZN/";
	$numFound++;
    }

    return ($numFound);

}

sub analyzeRocmop {
    my ($basedir) = @_;
    my ($log, $props, $badFileArray, $numBadFiles, $rbDir);
    my ($numFound, @fileList, $targetDir);
    $numFound = 0;

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();
    $targetDir = $props->getTARGETDIR();

    $rbDir = $basedir."Rocmop/";
    @fileList = ($rbDir, $rbDir."RocmopControl.txt");
    $badFileArray = ModuleProcessor::checkFilesExist(\@fileList);
    $numBadFiles = @$badFileArray;
    if (!$numBadFiles) {
    $props->setKeyValuePair("ROCMOP", $TRUE);
    $NDAfiles{$rbDir."RocmopControl.txt"} = $targetDir."Rocmop/";
    $numFound++;
    }

    return ($numFound);

}
#************************************************************************************
#************************************************************************************
sub extract {
    my ($thePackage, $processorObjs) = @_;
    my ($log, $props, $obj, @error);

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();

    if (buildDirectoryStructure()) {return(1);} #if non-zero returned, problem making target directories

    $timestamp = localtime;
    $log->printMessage("\n$timestamp: Extracting NDA files to rocstar dataset\n");

    push(@error, ModuleProcessor::performExtract(\%NDAfiles));

    foreach $obj (@$processorObjs) {
	if ($obj->isActive()) {
	    push(@error, $obj->extract());
	    $log->printEndPhase("Extract NDA Files", ref($obj));
        }
    }
    
    $log->printSeparator();
    return (hasError(\@error));
}

#************************************************************************************
#************************************************************************************
 sub buildDirectoryStructure(){
    my ($log, $props);

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();

    if ($props->getROCSTARVERS() eq "3.0") {
	return(makeRS3Structure());
    }
    else {
	$log->printMessage("Unknown Rocstar version detected.  Only Rocstar 3.0 currently supported\n");
	return(1);
    }
 }

#************************************************************************************
#************************************************************************************
 sub makeRS3Structure(){
    my ($log, $props);
    my ($target, $module, $fluidsmodule, $solidsmodule);

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();

    $target = $props->getTARGETDIR();

    if (ModuleProcessor::makeDirIfNeeded($target)) {return(1);}  #if failed to make it

    foreach $module ("Rocflo", "Rocflu", "Rocfrac", "Rocsolid", "RocburnAPN", "RocburnPY", "RocburnZN") {
	if ($props->processModule(uc($module))){
	    if (ModuleProcessor::makeDirIfNeeded($target.$module)) {
		return(1);  #if failed to make it
	    } 
	    ModuleProcessor::makeModuleDirectories($target.$module."/");
	}
    }

    #Make the Rocmop Directory
    if(ModuleProcessor::makeDirIfNeeded($target."Rocmop")) {
    	return(1); #if failed to make it
    }

    #Make the Rocman Directories
    if(ModuleProcessor::makeDirIfNeeded($target."Rocman")) {
	return(1); #if failed to make it
    }
    if(ModuleProcessor::makeDirIfNeeded($target."Rocman/Modout")) {
	return(1); #if failed to make it
    }
    if(ModuleProcessor::makeDirIfNeeded($target."Rocman/Profiles")) {
	return(1); #if failed to make it
    }

    #Make all the necessary module-combination directories under Rocman (e.g., RocfloRocfrac). 
    foreach $fluidsmodule ("Rocflo", "Rocflu") {
	if ($props->processModule(uc($fluidsmodule))){
	    foreach $solidsmodule("Rocfrac", "Rocsolid") {
		if ($props->processModule(uc($solidsmodule))){
		    if(ModuleProcessor::makeDirIfNeeded($target."Rocman/".$fluidsmodule.$solidsmodule)) {
			return(1);
		    }	
		}
	    }
	}
    } 
 }


#************************************************************************************
#************************************************************************************
sub preProcess {
    my ($thePackage, $processorObjs) = @_;
    my ($log, $props, $obj, @error);

    $props = ModuleProcessor::getProps();
    $log = ModuleProcessor::getLog();

    if (checkRS3PreprocDirs()) {return(1);}

    $timestamp = localtime;
    $log->printMessage("\n$timestamp: Running preprocessor codes to make rocstar dataset\n");

    if (!$rface) {$rface = RocfaceProcessor->new($props, $log);}

    foreach $obj (@$processorObjs) {
	if ($obj->isActive()) {
	    push(@error, $obj->preProcess());
	    $log->printEndPhase("Run Preprocessors", ref($obj));
        }
    }
    
    #Now, call the preProcess stage for Rocface, so surfdiver gets called
    push(@error, $rface->preProcess());
    $log->printEndPhase("Run Preprocessors", ref($rface));

    $log->printSeparator();
    return (hasError(\@error));
}

#************************************************************************************
#************************************************************************************
 sub checkRS3PreprocDirs(){
    my ($log, $props);
    my ($target, $module, $fluidsmodule, $solidsmodule);
    my ($error, $endflag);
    $error = 0;

    $props = ModuleProcessor::getProps();
    $log  = ModuleProcessor::getLog();

    $target = $props->getTARGETDIR();

    if (!ModuleProcessor::checkDirExists($target)) {return(1);}  #if not there

    foreach $module ("Rocflo", "Rocflu", "Rocfrac", "Rocsolid") {
	if ($props->processModule(uc($module))){
	    if (!ModuleProcessor::checkDirExists($target.$module)) {
		$log->processErrorCode(22, $target.$module);
		$error++;
	    }
	    if (!ModuleProcessor::checkDirExists($target.$module."/Modin")) {
		$log->processErrorCode(22, $target.$module."Modin");
		$error++;
	    }
	    if (!ModuleProcessor::checkDirExists($target.$module."/Rocin")) {
		$log->processErrorCode(22, $target.$module."Rocin");
		$error++;
	    }
	}
    }
    if ($error) {
	return(1);
    }

    #Check the Rocman Directory
    if(!ModuleProcessor::checkDirExists($target."Rocman")) {
	$log->processErrorCode(22, $target."Rocman");
	return(1); #if not there
    }

    #Check all the necessary module-combination directories under Rocman (e.g., RocfloRocfrac). 
    foreach $fluidsmodule ("Rocflo", "Rocflu") {
	if ($props->processModule(uc($fluidsmodule))){
	    foreach $solidsmodule("Rocfrac", "Rocsolid") {
		if ($props->processModule(uc($solidsmodule))){
		    if (!ModuleProcessor::checkDirExists($target."Rocman/".$fluidsmodule.$solidsmodule)) {
			$log->processErrorCode(22,$target."Rocman/".$fluidsmodule.$solidsmodule);
			$error++;
		    }	
		}
	    }
	}
    }
    return($error);
 }


#************************************************************************************
#************************************************************************************
sub checkResults {
    my ($thePackage, $processorObjs) = @_;
    my ($log, $props, $obj, @error);

    $props = ModuleProcessor::getProps();
    $log = ModuleProcessor::getLog();

    if (!$rface) {$rface = RocfaceProcessor->new($props, $log);}

    $timestamp = localtime;
    $log->printMessage("\n$timestamp: Checking rocstar dataset files for consistency\n");

    push(@error, checkRocstarFiles());

    foreach $obj (@$processorObjs) {
	if ($obj->isActive()) {
	    push(@error, $obj->checkResults());
	    $log->printEndPhase("Check Rocstar Dataset Files", ref($obj));
        }
    }

    #Now, call the checkResults stage for Rocface
    push(@error, $rface->checkResults());
    $log->printEndPhase("Check Rocstar Dataset Files", ref($rface));

    $log->printSeparator();
    return (hasError(\@error));
}

#************************************************************************************
#************************************************************************************
sub checkRocstarFiles {
    my ($log, $props, $module);
    my ($error)=0;

    $props = ModuleProcessor::getProps();
    $log = ModuleProcessor::getLog();

    my ($targetDir) = $props->getTARGETDIR();

    if ($props->processModule("ROCBURN")) {
    analyzeRocburn($targetDir);
    #If any of the three Rocburn models is there, then it's OK
    unless ($props->processModule("ROCBURNAPN")||
            $props->processModule("ROCBURNPY")||
            $props->processModule("ROCBURNZN")) {
	$error +=1;
	$log->processErrorCode(27, "No Rocburn Modules Available");
    }
    }

    if ($props->processModule("ROCMOP")) {
        analyzeRocmop($targetDir);
        #If any of the three Rocburn models is there, then it's OK
        unless ($props->processModule("ROCMOP")){
            $error +=1;
            $log->processErrorCode(27, "No Rocmop module available");
        }
    }
#    if (checkRocstarControlKey("Rocpanda")) {
#	unless (ModuleProcessor::checkFileExists($targetDir."Rocman/RocpandaControl.txt")){
#	    $error +=1;
#	    $log->processErrorCode(27, "RocstarControl contains Rocpanda, but RocpandaControl.txt doesn't exist.");
#	};
#    }

    #The following check considers the current physics flags versus how RocstarControl.txt is set up.
    #Right now, if you are doing "-A" with all 4 physics codes, you will get two warnings.  This check
    #is really more suited to a checkonly run where you want to check the dataset to see if it is 
    #completely ready to run.  Right now this check may cause an error to be returned from this
    #routine.  Currently, nothing follows checkResults, so we don't short-circuit because of this, but
    #if we add any more stages, this might become a problem...

     foreach $module ("Rocflo", "Rocflu", "Rocfrac", "Rocsolid") {
	if ($props->processModule(uc($module))){
	    unless (checkRocstarControlKey($module)) {
		$log->processErrorCode(27, "RocstarControl.txt doesn't contain $module.");
		$error+=1;
	    }
	}
    }
    
    return($error);
}

#************************************************************************************
#************************************************************************************
sub checkRocstarControlKey {
    my ($keyWord) = @_;
    my ($log, $props, $controlFile);

    $props = ModuleProcessor::getProps();
    $log = ModuleProcessor::getLog();
    
    $controlFile = $props->getTARGETDIR()."RocstarControl.txt";

    if (!open(CONTROLFILE, $controlFile))
    {
	$log->processErrorCode(27, $controlFile);
	$log->printMessage("$controlFile must exist in problem set.");
        #false - key doesn't exist because file is missing
	return(0); 
    }    
    
    while (<CONTROLFILE>) {
	$_ =~ $keyWord;
	return(1); 
    }

    return(0);
}


#************************************************************************************
#  Subroutine that develops a list of files that should be in the NDA given user input
#  of a specific Data and Grid combination. 
#************************************************************************************
sub getNDAFiles {
    my ($sourcedir, $targetdir);
    my ($props);

    $props = ModuleProcessor::getProps();

    $sourcedir = $props->getSOURCEDIR()."Rocstar/";
    $targetdir = $props->getTARGETDIR();

    $NDAfiles{$sourcedir."RocstarControl.txt"} = $targetdir ;
 #   $NDAfiles{$sourcedir."Rocman/RocpandaControl.txt"} = $targetdir."Rocman/";
    if($props->processModule("ROCMOP")){
        $NDAfiles{$sourcedir."Rocmop/RocmopControl.txt"} = $targetdir."Rocmop";
    }
    $NDAfiles{$sourcedir."Rocman/RocmanControl.txt"} = $targetdir."Rocman";
}


#************************************************************************************
#
#************************************************************************************
sub hasError {
    my ($errors) = @_;
    my ($error);

    foreach $error (@$errors) {
	if (($error) > 0) {return(1)};
    }
    return(0);
}


1;
