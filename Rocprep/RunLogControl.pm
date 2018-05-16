#************************************************************************************
#  Class RunLogControl
#************************************************************************************

package RunLogControl;

use strict;

my ($props);

#************************************************************************************
#  Constructor for new RunLogControl objects.  Takes the Log filename as it's
#  parameter.  Calls allocate to do the real work.
#************************************************************************************
sub new {
	my($pkg, $filename, $theProps) = @_;
	my ($file);
	$props = $theProps;
	
	$file = $pkg->allocate($filename);

	printBanner();
	
	return $file;
}

#************************************************************************************
#  Opens the log file supplied in $filename and  "blesses" the
#  modulename string  with the appropriate class name (passed in through the $pkg variable).  
#************************************************************************************
sub allocate {
    my ($pkg, $filename) = @_;
    my ($moduleName);
	
    open (LOGFILE, "| tee $filename") ||
        warn "Could not open logfile $filename!!!\n";

    $moduleName = \getModuleName();

    bless ($moduleName, $pkg);
	
    return $moduleName;
}

sub closeLog {
    close (LOGFILE);
}

#************************************************************************************
#  Must be overridden by subclasses.  Returns the name of this type of properties
#  object to be used in messages. 
#************************************************************************************
sub getModuleName {
	return "RunLogControl";
}


#************************************************************************************
# Examine an array of error codes and process them
#************************************************************************************
sub processErrorList {
    my ($obj, $errorList, $msg) = @_;
    my ($error);

    foreach $error (@$errorList) {
	processErrorCode ($obj, $error, $msg);
    }
}

#************************************************************************************
# Examine the error code and process it appropriately.

# Code-specific errors:
# Range 10 - 100; generic errors (file or dir not found, etc.)
# Range 100-150; Rocflo specific errors
# Range 151-200; Rocflu specific errors
# Range 201-250; Rocfrac specific errors
# Range 251-300; Rocsolid specific errors
# Range 301-350; Rocburn specific errors
# Range 351-400; Rocstar specific errors
#
# Code 22 is for a bad directory
# Code 23 is for a bad file
# Code 24 is for a bad parameter read from a file
# Code 404 is an internal error for an unimplemented function
# Code 512 is an internal error for an unknown error. 
#************************************************************************************
sub processErrorCode {
    my ($obj, $errorCode, $msg) = @_;
    my ($theCaller);
    
    $theCaller = caller();  #introspectively find the calling package

    SWITCH:  {
	if ($errorCode == 0) {last SWITCH;}
	if ($errorCode == 10) {
	    print LOGFILE "\n--->ERROR in $theCaller: Bad Input Parameter Value: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 11) {
	    print LOGFILE "\n--->ERROR in $theCaller: Preprocessor Code Run failure:\n$msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 22) {
	    print LOGFILE "\n--->ERROR in $theCaller: directory not accessible: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 23) { 
	    print LOGFILE "\n--->ERROR in $theCaller: file not accessible: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 24) {
	    print LOGFILE "\n--->ERROR in $theCaller: bad parameter read: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 25) {
	    print LOGFILE "\n--->ERROR in $theCaller: cannot create directory: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 26) {
	    print LOGFILE "\n--->ERROR in $theCaller: Could not copy NDA source file: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 27) {
	    print LOGFILE "\n--->ERROR in $theCaller: Dataset incomplete: $msg\n\n";
	    last SWITCH;
	}
	if ($errorCode == 301) {
	    print LOGFILE "\n--->ERROR in $theCaller: No valid Rocburn Models found at: $msg\n\n";
	    last SWITCH;
	}
    if ($errorCode == 404) {
        print LOGFILE "\n--->ERROR in $theCaller: unimplemented function: $msg\n\n";
        last SWITCH;
    }
    if ($errorCode == 405) {
        print LOGFILE "\n--->ERROR in $theCaller: No Rocmop config found at: $msg\n\n";
        last SWITCH;
    }
    if ($errorCode == 512) {
	    print LOGFILE "\n\nMajor Internal Program Error: Error code.  Contact Code Authors\n\n";
	    printContactInfo();
	    last SWITCH;
        }
    }
}

#************************************************************************************
# Prints banner information and contact info
#************************************************************************************
sub printBanner {
    my ($prepvers, $starvers, $timestamp);
    $prepvers = $props->getROCPREPVERS();
    $starvers = $props->getROCSTARVERS();
    $timestamp = localtime;
    print LOGFILE "***************************************************************************\n";
    print LOGFILE "                  Rocprep Tool Version $prepvers\n";
    print LOGFILE "                  For Rocstar Version $starvers File formats\n\n";
    print LOGFILE "                  Center for Simulation of Advanced Rockets\n";
    print LOGFILE "                  University of Illinois, Urbana, IL  61801\n";
    print LOGFILE "                  www.csar.uiuc.edu\n\n";
    printContactInfo();
    print LOGFILE "***************************************************************************\n\n";
    print LOGFILE "$timestamp: Rocprep Initialized\n\n";
}


#************************************************************************************
# Prints banner information and contact info
#************************************************************************************
sub printContactInfo {
    print LOGFILE "                  Code Authors:\n";
    print LOGFILE "                  Mark Brandyberry (mdbrandy\@uiuc.edu)\n";
    print LOGFILE "                  Court McLay (cmclay\@uiuc.edu)\n";
}


#************************************************************************************
# Prints end of run message
#************************************************************************************
sub printRunTerminated {
    my ($obj, $message) = @_;
    print LOGFILE "\n***************************************************************************\n";
    print LOGFILE "\nRun terminated with error: $message\n\n";
    print LOGFILE "***************************************************************************\n";
}


#************************************************************************************
# Prints ending of phase message
#************************************************************************************
sub printEndPhase {
    my ($obj, $phase, $theCaller) = @_;

    print LOGFILE "Ending phase: $phase for module $theCaller. \n";
}

#************************************************************************************
# Prints supplied message to the Log
#************************************************************************************
sub printMessage {
    my ($obj, $msg) = @_;
    
    print LOGFILE $msg."\n";
}

#************************************************************************************
#************************************************************************************
sub printMessageList {
    my ($obj, @msgList) = @_;
    my ($msg);

    foreach $msg (@msgList) {
       print LOGFILE $msg."\n";
    }
}

#************************************************************************************
#************************************************************************************
sub printSeparator {
    print LOGFILE "***************************************************************************\n";
}

1;
