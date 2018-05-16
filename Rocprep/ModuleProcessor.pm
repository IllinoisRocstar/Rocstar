#************************************************************************************
#  Class ModuleProcessor

#  Superclass to all other Processor classes. 
#************************************************************************************

package ModuleProcessor;

use strict;
use File::Spec;

my ($props, $log);

#************************************************************************************
#  Constructor for new ModuleProcessor objects.
#  Calls allocate to do the real work.
#************************************************************************************
sub new {
	my($pkg, $myprops, $mylog) = @_;
	my ($prop);

	$props = $myprops;
	$log = $mylog;

	$prop = $pkg->allocate();
	
	return $prop;
}

#************************************************************************************
# Accessor methods for class attributes - properties list and logfile
#************************************************************************************
sub getProps {
    return ($props);
}

sub getLog {
    return ($log);
}

#************************************************************************************
#  Then echos out the properties read and "blesses" the
#  hash with the appropriate class name (passed in through the $pkg variable).  
#************************************************************************************
sub allocate {
    my ($pkg) = @_;
    my ($moduleName);

    $moduleName = \getModuleName();

    bless ($moduleName, $pkg);
	
    return $moduleName;
	
}

#************************************************************************************
#  Must be overridden by subclasses.  Returns the name of this type of properties
#  object to be used in messages. 
#************************************************************************************
sub getModuleName {
	return "ModuleProcessor";
}


#************************************************************************************
#  Must be overridden by subclasses. 
#  Subroutine for each module processor subclass to check it's NDA hierarchy. 
#************************************************************************************
sub checkNDA {
	print "Must Override in subclass!!!! (ModuleProcessor::checkNDA)\n";
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses. 
#  Subroutine for each module processor subclass to check existance of files in
#  the Rocstar 3.0 directory before preprocessing. 
#************************************************************************************
sub checkPreProcessorInput {
	print "Must Override in subclass!!!! (ModuleProcessor::extract)\n";
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses. 
#  Subroutine for each module processor subclass to run the extract steps. 
#************************************************************************************
sub extract {
	print "Must Override in subclass!!!! (ModuleProcessor::extract)\n";
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses. 
#  Subroutine for each module processor subclass to run the preprocessing steps. 
#************************************************************************************
sub preProcess {
	print "Must Override in subclass!!!! (ModuleProcessor::preProcess)\n";
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses. 
#  Subroutine for each module processor subclass to check the preprocessor output. 
#************************************************************************************
sub checkResults {
	print "Must Override in subclass!!!! (ModuleProcessor::checkResults)\n";
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses.  
#  Subroutine that develops a list of files that should be in the NDA given user input
#  of a specific Data and Grid combination. 
#************************************************************************************
sub getNDAFiles {
	print "Must Override in subclass!!!! (ModuleProcessor::getNDAFiles)\n";
	#real function should return an array of file path/names. 
	return (512);
}

#************************************************************************************
#  Must be overridden by subclasses.  
#  Subroutine that develops a list of files that should be in the dataset after
#  preprocessing, so that the dataset is considered ready to run.  Also used for
#  checking an existing dataset using the -C flag. 
#************************************************************************************
sub getRuntimeFiles {
	print "Must Override in subclass!!!! (ModuleProcessor::getRuntimeFiles)\n";
	#real function should return an array of file path/names. 
	return (512);
}



#************************************************************************************
# Takes a list of files and directories (files with relative or full paths), and 
# checks for their existance and read status (also executability for directories).  
# Returns a reference to an array containing a list of any that fail the test. 
#************************************************************************************
sub checkFilesExist {
    my ($filesToCheck) = @_;
    my ($filename);
    my (@missingFiles);
    
    foreach $filename (@$filesToCheck) {
	if (!((-e $filename)&&(-r $filename))){
	    push(@missingFiles, $filename);
	}
	if (-d $filename) {
	    if (!(-x $filename)) {
		push(@missingFiles, $filename);
	    }
	}
    }
    return (\@missingFiles);
}

#************************************************************************************
#
#************************************************************************************
sub checkFileExists {
    my ($file) = @_;
    my ($badFile, $numBad);
    my (@files);

    @files = ($file);

    $badFile = checkFilesExist(\@files);
    $numBad = @$badFile;

    if ($numBad) {return (0);}
    else {return(1);}
}

#************************************************************************************
#
#************************************************************************************
sub checkDirExists {
    my ($dir) = @_;
    my ($badDir, $numBad);
    my (@dirs);

    @dirs = ($dir);

    $badDir = checkFilesExist(\@dirs);
    $numBad = @$badDir;

    if ($numBad) {return (0);}
    else {return(1);}
}

#************************************************************************************
#
#************************************************************************************
sub checkDataGridExists {
    my ($datadir, $griddir) = @_;
    my ($result);
    $result = 1;
    
    if (!checkDirExists($datadir)) {
	$log->processErrorCode(22, $datadir);
	$result = 0;
    }
    if (!checkDirExists($griddir)) {
	$log->processErrorCode(22, $griddir);
	$result = 0;
    }
    
    return($result);
}


#************************************************************************************
#
#************************************************************************************
sub makeModuleDirectories {
    my ($basedir) = @_;
    my ($directory);

    foreach $directory ("Rocin", "Rocout", "Modin", "Modout") {

	if(makeDirIfNeeded($basedir.$directory)) {
	    return(1); #couldn't successfully make directory
	}
    }

    return(0);
}

#************************************************************************************
#************************************************************************************
sub makeDirIfNeeded {
    my ($dirToMake) = @_;
    
    if (!ModuleProcessor::checkDirExists($dirToMake)) { #if directory doesn't already exist...
	if(!mkdir($dirToMake, 0750)) {
	    $log->processErrorCode(25, $dirToMake);
	    return(1); #directory doesn't exist, and error making it.
	}
    }
    return(0); #directory already exists (OK), or was make successfully
}

#************************************************************************************
#************************************************************************************
sub performExtract {
    my ($fileHash) = @_;
    my ($fileName, $copyTo, $result, $error);
    $result = 0;

    foreach $fileName (keys(%$fileHash)) { #fileArray has the paths on the files already
	$copyTo = $$fileHash{$fileName};
        $error = system("cp $fileName $copyTo");
        if ($error) {
	    $log->processErrorCode(26, $fileName);
	    $result = $error;
        }
    }
    return($result);
}

#************************************************************************************
# returns an array of .hdf files needed by Rocin
#************************************************************************************
sub getRocinFiles {
    my ($rocinFile) = @_;
    my ($rocinFH, $filePattern, $keyword, $workingDir, @fileList);
    my ($volume, $rocinDir, $filename) = File::Spec->splitpath($rocinFile);

    if (!open($rocinFH, $rocinFile))
    {
	$log->processErrorCode(23, $rocinFile);
	$log->printMessage("$rocinFile must exist to supply names for Rocin .hdf run files.");
	return;
    }

    $keyword = '^\s*@Files:\s*';

    while(<$rocinFH>) {
       chomp;
       if ($_ =~ /$keyword/) {
	   $filePattern = $_;
	   last;
       }
    }

    # strip @Files keyword and spaces from the beginning and end of filename pattern
    $filePattern =~ s/$keyword//;
    $filePattern =~ s/\s+$//;

    # control file pattern must contain one valid character (a-z,A-Z,0-9,_)
    if ($filePattern =~ /\w+/) {
       # control file pattern can be relative to target dir or control file
       if ($filePattern =~ '/') {
          $workingDir = $props->getTARGETDIR();
       }
       else {
          $workingDir = $rocinDir;
       }
       @fileList = buildPatternFiles($filePattern, $workingDir, $rocinFH);
    }
    close($rocinFH);
    return @fileList;
}


#************************************************************************************
#
#************************************************************************************
sub buildPatternFiles {
    my ($pattern, $directory, $rocinFH) = @_;
    my ($i, $numProcs, $numBlocks, @fileList);
    my ($fileName, $fieldWidth, $proc, $block, $keyword);

    # first, find and replace any timestamp patterns %t
    my ($timestamp) = '%t';
    if ($pattern =~ $timestamp) {
	$pattern =~ s/$timestamp/00.000000/g;
    }

    # find any process rank placeholder %p and replace with proc number i = 0..N
    if ($pattern =~ /%(\d*)p/) {
        $numProcs = $props->getNumProcs();
	$fieldWidth = $1;
        unless ($fieldWidth) {
           $fieldWidth = 4;
        }
	for ($i= 0; $i< $numProcs; $i++) {
	    $fileName = $pattern;
	    $proc = sprintf("%0$fieldWidth"."d", $i);  # formatted print pads with zeroes
	    $fileName =~ s/%\d*p/$proc/g;
	    push(@fileList, $directory.$fileName);
	}
    }

    # find any block ID placeholder %b and replace with block number i = 0..N
    if ($pattern =~ /%(\d*)b/) {
	$fieldWidth = $1;
        unless ($fieldWidth) {
           $fieldWidth = 4;
        }

        $keyword = '^\s*@Panes:\s*@Block\s*(\d+)';
        while(<$rocinFH>) {
           if ($_ =~ /$keyword/) {
	       $numBlocks = $1;
	       last;
           }
        }

        unless ($numBlocks) {
	   $log->processErrorCode(24, $numBlocks);
	   $log->printMessage("Rocin control file has bad block number, cannot verify .hdf files.");
	}

	for ($i=1; $i<= $numBlocks; $i++) {
	    $fileName = $pattern;
	    $block = sprintf("%0$fieldWidth"."d", $i);  # formatted print pads with zeroes
	    $fileName =~ s/%\d*b/$block/g;
	    push(@fileList, $directory.$fileName);
	}
    }

    return(@fileList);
}


1;
