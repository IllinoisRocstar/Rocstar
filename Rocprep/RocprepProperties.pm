#The RocprepProperties package inherits from Properties.pm.  It adds access methods
#for all propeties in the main RocprepControl.txt Configuration file. 

package RocprepProperties;

require "Properties.pm";
@RocprepProperties::ISA = qw(Properties); #Inherits from Properties
use strict;

#************************************************************************************
# If the result of a property lookup is an empty string, then the property was not
# found - print error message.  Don't call this function if the property is not
# required to be in the properites object. 
#************************************************************************************
sub processNotFound {
	my ($result, $property) = @_;
	
	if ($result eq "") {
		print "\n$property variable is unavailable\n";
		$result="NOT FOUND";
	}
	
	return $result;
}

#************************************************************************************
# Get the value in the ROCFLO property
#************************************************************************************
sub getROCFLO {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCFLO");

	return 	(processNotFound($result, "ROCFLO"));
}

#************************************************************************************
# Get the value in the ROCFLU property
#************************************************************************************
sub getROCFLU {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCFLU");
	
	return 	(processNotFound($result, "ROCFLU"));
}

#************************************************************************************
# Get the value in the ROCFRAC property
#************************************************************************************
sub getROCFRAC {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCFRAC");
	
	return 	(processNotFound($result, "ROCFRAC"));
}

#************************************************************************************
# Get the value in the ROCSOLID property
#************************************************************************************
sub getROCSOLID {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCSOLID");
	
	return 	(processNotFound($result, "ROCSOLID"));
}


#************************************************************************************
# Get the value in the *DATA parameter for the module indicated
#************************************************************************************
sub getDATADIR {
	my ($props, $module) = @_;
	my ($data) = "";
	my ($grid) = "";
	my ($result) = "";
	
	$data = $props->getProp($module."DATA");
	$grid = $props->getProp($module."GRID");
	
	$result = "$grid/$data";
	
	return 	(processNotFound($result, $module."DATA"));
}


#************************************************************************************
# Get the value of the BC corresponding to an injecting surface
#************************************************************************************
sub getPLAG {
  my ($props) = shift;
  my ($result) = "";

  $result = $props->getProp("PLAG");

  return( processNotFound($result,"PLAG") );
}

#************************************************************************************
# Get the value in the *GRID parameter for the module indicated
#************************************************************************************
sub getGRIDDIR {
	my ($props, $module) = @_;
	my ($result) = "";
	
	$result = $props->getProp($module."GRID");
	
	return 	(processNotFound($result, $module."GRID"));
}


#************************************************************************************
# Get the value in the ROCPREPVERS property
#************************************************************************************
sub getROCPREPVERS {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCPREPVERS");
	
	return 	(processNotFound($result, "ROCPREPVERS"));
}


#************************************************************************************
# Get the value in the ROCSTARVERS property
#************************************************************************************
sub getROCSTARVERS {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCSTARVERS");
	
	return 	(processNotFound($result, "ROCSTARVERS"));
}

#************************************************************************************
# Get the value in the TARGETDIR  property
#************************************************************************************
sub getTARGETDIR {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("TARGETDIR");
	
	return 	(processNotFound($result, "TARGETDIR"));
}

#************************************************************************************
# Get the value in the SOURCEDIR property
#************************************************************************************
sub getSOURCEDIR {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("SOURCEDIR");
	
	return 	(processNotFound($result, "SOURCEDIR"));
}

#************************************************************************************
# Get the value in the BINDIR property
#************************************************************************************
sub getBINDIR {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("BINDIR");
	
	# don't use processNotFound() because null-string is the default for BINDIR!
	# using processNotFound() will return "NOT FOUND" instead of null, which causes
	# big problems when you assemble the commandline path for the preptools
	return 	($result);
}

#************************************************************************************
# Get the value in the NUMPROCS property
#************************************************************************************
sub getNumProcs {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("NUMPROCS");
	
	return 	(processNotFound($result, "NUMPROCS"));
}

#************************************************************************************
# Get the value in the NUMREGIONS property
#************************************************************************************
sub getNumRegions {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("NUMREGIONS");
	
	return 	(processNotFound($result, "NUMREGIONS"));
}

#************************************************************************************
# Get the value in the ROCFRACUNITS property
#************************************************************************************
sub getRocfracUnits {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ROCFRACUNITS");
	
	return 	(processNotFound($result, "ROCFRACUNITS"));
}

#************************************************************************************
# Get the value in the SPLITAXIS property
#************************************************************************************
sub getSplitAxis {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("SPLITAXIS");
	
	return 	(processNotFound($result, "SPLITAXIS"));
}

#************************************************************************************
# Get the value in the IGNOREFILE property
#************************************************************************************
sub getIgnoreFile {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("IGNOREFILE");
	
	return 	(processNotFound($result, "IGNOREFILE"));
}

#************************************************************************************
# If the -C flag was submitted, then the user wants to only check an existing problem
# dataset - return true from this function if CHECKONLY is set.  Otherwise, false. 
#************************************************************************************
sub runCheckOnly {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("CHECKONLY");

	return decideBoolean($result);
	
	return $result;
}

#************************************************************************************
# If the -E flag was submitted, then the user wants to only extract files to a target
# dataset - return true from this function if EXTRACTONLY is set.  Otherwise, false. 
#************************************************************************************
sub runExtractOnly {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("EXTRACTONLY");

	return decideBoolean($result);
	
	return $result;
}

#************************************************************************************
# If the -P flag was submitted, then the user wants to only extract files to a target
# dataset - return true from this function if EXTRACTONLY is set.  Otherwise, false. 
#************************************************************************************
sub runPreprocessOnly {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("PREPROCESS");

	return decideBoolean($result);
	
	return $result;
}

#************************************************************************************
# If the -E flag was submitted, then the user wants to only extract files to a target
# dataset - return true from this function if EXTRACTONLY is set.  Otherwise, false. 
#************************************************************************************
sub runAll {
	my ($props) = shift;
	my ($result) = "";
	
	$result = $props->getProp("ALL");

	return decideBoolean($result);
	
	return $result;
}


#************************************************************************************
# 
#************************************************************************************
sub processRocfaceCombo {
    my ($props, $combination) = @_;

    return(processModule($props, $combination));
}

#************************************************************************************
# 
#************************************************************************************
sub processModule {
	my ($props, $module) = @_;
	my ($result) = "";
	
	$result = Properties::getProp($props, $module);
	
	return decideBoolean($result);
}

#************************************************************************************
# Takes an input value, checks if it is a valid 0 or 1, and returns that value if so.
# Otherwise, it returns zero (false). 
#************************************************************************************
sub decideBoolean {
    my ($result) = @_;
    
    #if parameter is blank string, then the property was not found. 
    if (($result eq "NOT FOUND")||($result eq "")) {$result = 0};

    #if it isn't zero or one, then it is invalid
    if (!(($result == 0) || ($result == 1))) {$result = 0};

    #otherwise, return the original value. 
    return($result);

}

#************************************************************************************
# Function returns a list of properties that must be available in the properties file
# After it is constructed by reading in the file.  For this program, there are none. 
#************************************************************************************
sub getPropsToCheck {
	my (@propList) = ();
	return \@propList;
}

#************************************************************************************
# Function writes properties file for restarting preprocessing after extract-only from
# the NDA.  Only a subset of the properties will be dumped. The file to be written is
# passed in as a parameter with the full path of the file.  
#************************************************************************************
sub writePropFile {
    my($props, $fileName) = @_;
    my ($success, $key, @keys);

    $success = open(PROPFILE, ">$fileName");

    #TBD - need to write the following message to the log instead. 
    
    if (!$success) {
	print "Cannot open file $fileName for writing\n\n";
    }

    @keys = qw(ROCFLO ROCFLU ROCFRAC ROCSOLID NUMPROCS NUMREGIONS ROCFLOROCFRAC ROCFLOROCSOLID ROCFLUROCFRAC ROCFLUROCSOLID);
    
    foreach $key (@keys) {
	if (Properties::propExists($props, $key)) {
	    print PROPFILE $key."->".$$props{$key}."\n";
	}
    }
}

sub getPropfileName {
	return "";
}

1;
