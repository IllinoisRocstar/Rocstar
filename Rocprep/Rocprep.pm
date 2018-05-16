#!/usr/bin/perl

#************************************************************************************
#CSAR Rocstar Dataset Preparation Program
#Mark Brandyberry, Court McLay:  9/10/2004
#************************************************************************************

package Rocprep;

use strict;

# need File::Spec for handling relative to absolute path conversions
use File::Spec;

# FindBin locates the full path to the executable. Next, we dynamically
# add that path to the search path for supporting perl modules that follow
use FindBin qw($RealBin);
use lib "$RealBin";

use RocprepProperties;
use RocfluProcessor;
use RocfloProcessor;
use RocfracProcessor;
use RocsolidProcessor;
use RocstarProcessor;
use RunLogControl;


#************************************************************************************
#  The Rocprep.pm perl module contains the main program and supporting subroutines
#  for the Dataset Preparation Program.  It depends on a properties file (see package
#  Properties in file Properties.pm and package RocprepProperties.pm), as well as a number
#  of other packages. 
# 
#  The main configuration file must be named RocprepControl.txt and be placed in the root
#  problem directory, OR command line parameters used OR both (command line parameters override 
#  configuration file information.  A sample RocprepControl.txt file has:

#  TBD

#************************************************************************************


my ($RocProps, $PrepLog, $TRUE, $FALSE );

$TRUE=1;
$FALSE=0;

main();

#************************************************************************************
#************************************************************************************
sub main {

	# check perl version before continuing
	if ($] < 5.006) {
	   print "Perl version $] is too old, exiting\n";
	   print "Please install Perl 5.6.0 or greater\n\n";
	   exit;
	}

	my ($rflo, $rflu, $rfrac, $rsolid, $rstar);
	my (@processorObjs);
	
	#Read RocprepControl.txt, and process any command line arguments
	processARGs(@ARGV );

	#start log file
	$PrepLog = RunLogControl->new("rocprep.log", $RocProps);
        $PrepLog->printMessageList($RocProps->dumpProperties());

	#instantiate Processor Objects
	$rstar = RocstarProcessor->new($RocProps, $PrepLog);
	$rflo = RocfloProcessor->new($RocProps, $PrepLog);
	$rflu = RocfluProcessor->new($RocProps, $PrepLog);
	$rfrac = RocfracProcessor->new($RocProps, $PrepLog);
	$rsolid = RocsolidProcessor->new($RocProps, $PrepLog);

	@processorObjs = ($rflo, $rflu, $rfrac, $rsolid);

	if (($RocProps->runExtractOnly())||($RocProps->runAll())) {

	     if ($rstar->checkNDA(\@processorObjs)) {
	        $PrepLog->printRunTerminated("Error(s) were found in the NDA. Correct and rerun $0.");
	        exit;
	     }

	     if ($rstar->extract(\@processorObjs)) {
	        $PrepLog->printRunTerminated("Error(s) were found during file extraction from the NDA.");
	        exit;
	     }

	     if ($RocProps->runExtractOnly()) {
		 $RocProps->writePropFile($RocProps->getTARGETDIR()."RocprepControl.txt");
	         $PrepLog->printRunTerminated("NO ERRORS");
	         exit;
     	     }
        }

	unless ($RocProps->runCheckOnly()) {
	   if ($rstar->preProcess(\@processorObjs)) {
	      $PrepLog->printRunTerminated("Errors occurred during preprocessing.");
	      exit;
	   }
	}

	if ($rstar->checkResults(\@processorObjs)) {
	   $PrepLog->printRunTerminated("Error(s) were found in the rocstar dataset.");
	   exit;
	}

	$PrepLog->printRunTerminated("NO ERRORS");
}


#************************************************************************************
#************************************************************************************
sub processARGs {
    my (@cline) = @_;
    my ($arg);
    my ($surfdone)=$FALSE;
	
    $RocProps = RocprepProperties->new();
    
    #Check the first switch, which determines Rocprep's major mode
    $arg = $cline[0];
    if (($arg !~ /^-[ACEPU]$/)&&($arg ne '-h')&&($arg ne '--help')&&($arg ne '--all')&&
        ($arg ne '--check')&&($arg ne '--extract')&&($arg ne '--preprocess') ) {
      print "First switch must be mode switch -A|C|E|P|U, not: $arg\n\n";
      printUsage();
      exit;
    }

    while (@cline) {
      $arg = shift @cline;
      SWITCH: {
	  if (($arg eq "-A")||($arg eq "--all")) { 
	      A_command(\@cline);
	      last SWITCH;
	  }
	  if (($arg eq "-C")||($arg eq "--check")) { 
	      C_command(\@cline);
	      last SWITCH;
	  }
	  if (($arg eq "-E")||($arg eq "--extract")) { 
	      E_command(\@cline);
	      last SWITCH;
	  }
	  if (($arg eq "-P")||($arg eq "--preprocess")) { 
	      P_command(\@cline);
	      last SWITCH;
	  }
	  if ($arg eq "-o") { 
	      o_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-u") {
	      u_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-f") {
	      f_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-s") { 
	      s_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-i") { 
	      i_command(\@cline);
	      $surfdone=$TRUE;
	      last SWITCH;
          }
	  if ($arg eq "-b") { 
	      b_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-d") {  
	      d_command(\@cline);
	      last SWITCH;
          }
	  if (($arg eq "-h") || ($arg eq "--help")) { 
	      h_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-m") { 
	      m_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-n") { 
	      n_command(\@cline);
	      last SWITCH;
          }
	  if ($arg eq "-t") { 
	      t_command(\@cline);
	      last SWITCH;
          } 
	  if ($arg eq "-p") { 
	      p_command(\@cline);
	      last SWITCH;
          } 
	  if ($arg eq "-r") { 
	      r_command(\@cline);
	      last SWITCH;
          } 
	  if ($arg eq "-un") { 
	      un_command(\@cline);
	      last SWITCH;
          } 
	  if ($arg eq "-splitaxis") { 
	      split_command(\@cline);
	      last SWITCH;
          } 
	  if (($arg eq "-x")||($arg eq "--ignore")) { 
	      x_command(\@cline);
	      last SWITCH;
	  }
    if (($arg eq "-plag")||($arg eq "--particles")) {
      plag_command(\@cline);
      last SWITCH;
    }
	  bad_command($arg);
        }
    }

    setCriticalDefaults();

    # In every mode read the RocprepControl.txt file without overwriting values already
    # set by the commandline parsing. Advanced switch (-x, --ignore) will not read this file.
    if (!$RocProps->getIgnoreFile()) {
       $RocProps->readFile($RocProps->getSOURCEDIR()."RocprepControl.txt", $FALSE);
    }

    # set fallback properties not already read in from commandline or .txt file
    setDefaultProps();

    checkInputConsistency($surfdone);
}

#************************************************************************************
# set up default key-value pairs that must be available in Rocprops.  The
# default values are generally false (i.e, don't process), but don't overwrite any
# pre-existing settings read from the commandline or control file.
#************************************************************************************
sub setDefaultProps {
   my ($keyword);

   $RocProps->setKeyValuePair("ROCPREPVERS", "1.0");
   $RocProps->setKeyValuePair("ROCSTARVERS", "3.0");

   # physics modules
   foreach $keyword (qw(ROCFLO ROCFLU ROCFRAC ROCSOLID ROCBURNAPN ROCBURNPY ROCBURNZN ROCMOP)) {
      unless ($RocProps->propExists($keyword)) {
         $RocProps->setKeyValuePair($keyword, $FALSE);
      }
   }

   # When genx no longer requires a rocburn module to be present, CHANGE THE DEFAULT TO $FALSE:
   unless ($RocProps->propExists("ROCBURN")) {
      $RocProps->setKeyValuePair("ROCBURN", $FALSE);
   }

   # Rocmop
   unless ($RocProps->propExists("ROCMOP")) {
     $RocProps->setKeyValuePair("ROCMOP", $FALSE);
   }

   # surfdiver combinations
   foreach $keyword (qw(ROCFLOROCFRAC ROCFLOROCSOLID ROCFLUROCFRAC ROCFLUROCSOLID)) {
      unless ($RocProps->propExists($keyword)) {
         $RocProps->setKeyValuePair($keyword, $FALSE);
      }
   }

   # default path to preptools will search user's shell path
   unless ($RocProps->propExists("BINDIR")) {
      $RocProps->setKeyValuePair("BINDIR", "");
   }
}

#************************************************************************************
# build a default source directory path if none has been given. Also, do not ignore
# the RocprepControl.txt file by default.
#************************************************************************************
sub setCriticalDefaults {
   my ($currentDir);

   unless ($RocProps->propExists("SOURCEDIR")) {
      $currentDir = File::Spec->curdir();
      $currentDir = File::Spec->rel2abs($currentDir);
      $currentDir = checkPathSlash($currentDir);
      $RocProps->setKeyValuePair("SOURCEDIR", $currentDir);
   }

   unless ($RocProps->propExists("IGNOREFILE")) {
      $RocProps->setKeyValuePair("IGNOREFILE", $FALSE);
   }
}

#************************************************************************************
# sanity check properties after the commandline & properties file have been processed
#************************************************************************************
sub checkInputConsistency {
    my ($surfdone) = @_;

#   It makes no sense to continue if nothing is to be extracted, processed, or checked. 
#   So we check to see if one or more of the Physics code flags has been specified, after
#   which we check for any surfdiver flags.
    my($keyword);
    my($doPrep)=0;
    foreach $keyword (qw(ROCFLO ROCFLU ROCFRAC ROCSOLID)) {
	if ($RocProps->processModule($keyword)){
	    $doPrep++;
	    last;
	}
    }
    foreach $keyword (qw(ROCFLOROCFRAC ROCFLOROCSOLID ROCFLUROCFRAC ROCFLUROCSOLID)) {
	if ($RocProps->propExists($keyword)) {
	    $doPrep++;
	    last;
	}
    }
    if (!$doPrep) {
	print "You have not specified any module to extract/process/check.\n";
	print "Please use the -o, -u, -f, -s, or -i flags to specify which module(s) to use.\n";
	printUsage();
	exit;
    }

    # Overwrite the target directory with the source directory for the 2 modes below.
    if (($RocProps->runCheckOnly())||($RocProps->runPreprocessOnly()))
    {
	my $target = $RocProps->getSOURCEDIR();
	if (!ModuleProcessor::checkDirExists($target)) {
	    print "The directory you specified: $target does not exist, or is not writable!\n";
	    print "Please specify a valid directory with the -d switch.\n";
	    printUsage();
	    exit;
	}
	#Set target to source, because preprocessing code works on target directory.  This is
	#done because in "All" mode, the NDA is the source, and the Rocstar directory is the
	#target.  In runCheckOnly, and preProcessOnly the -d (source) flag is used to specify
	#the single directory to be acted upon, but the same code is used for preprocess or 
	#checkResults, and that code expects targetdir. Sourcedir is used to find RocprepControl.txt,
	#which is read in automatically unless -x or --ignore is specified.
	$RocProps->setKeyValuePair("TARGETDIR", $RocProps->getSOURCEDIR());
    }

    if (($RocProps->runExtractOnly() || $RocProps->runAll())&&(!$RocProps->propExists("TARGETDIR")))
    {
	print "You MUST specify a target directory for $0 to extract to!\n"; 
	print "Please specify a target directory using the -t switch.\n\n";
	printUsage();
	exit;
    }

    unless ($RocProps->runExtractOnly()) {
	unless ($RocProps->propExists("NUMPROCS")) {
	    print "You have not specified the number of processors to use for partitioning.\n";
	    print "Please use the -n switch to provide the number of processors.\n";
	    printUsage();
	    exit;
	}
    }
    
    #now, if we have received no command line switch to tell us what to do with surfdiver, then 
    #build the combinations required given the physics codes specified. 
    my ($fluidsmodule, $solidsmodule);
    if (!$surfdone) {
        foreach $fluidsmodule ("Rocflo", "Rocflu") {
	    if ($RocProps->processModule(uc($fluidsmodule))){
		foreach $solidsmodule("Rocfrac", "Rocsolid") {
		    if ($RocProps->processModule(uc($solidsmodule))){
			$RocProps->setKeyValuePair(uc($fluidsmodule).uc($solidsmodule), $TRUE);
		    }
		}
	    }
	}
    }

}

#************************************************************************************
#************************************************************************************
sub o_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -o [m] [n]     Rocflo preprocessing, optional NDA Data<m> & Grid<n> dirs\n";
      return;
    }

    if ($RocProps->runExtractOnly() || $RocProps->runAll()) {
       #Rocflo: format -o 1 2  where 1 is the dataset and 2 is the grid
       $RocProps->setKeyValuePair("ROCFLO", $TRUE);
       $RocProps->setKeyValuePair("ROCFLODATA", "Data".shift(@$cline)); #integer input
       $RocProps->setKeyValuePair("ROCFLOGRID", "Grid".shift(@$cline)); #integer input
    }
    elsif (($RocProps->runPreprocessOnly())||($RocProps->runCheckOnly())) {
       $RocProps->setKeyValuePair("ROCFLO", $TRUE);
    }
}

sub u_command {
    my ($cline) = @_;
 
    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -u [m] [n]     Rocflu preprocessing, optional NDA Data<m> & Grid<n> dirs\n";
      return;
    }

    if ($RocProps->runExtractOnly() || $RocProps->runAll()) {
       #Rocflu: format -u 1 2  where 1 is the dataset and 2 is the grid
       $RocProps->setKeyValuePair("ROCFLU", $TRUE);
       $RocProps->setKeyValuePair("ROCFLUDATA", "Data".shift(@$cline)); #integer input
       $RocProps->setKeyValuePair("ROCFLUGRID", "Grid".shift(@$cline)); #integer input
    }
    elsif  (($RocProps->runPreprocessOnly())||($RocProps->runCheckOnly())) {
       $RocProps->setKeyValuePair("ROCFLU", $TRUE);
    }
}

sub f_command {
    my ($cline) = @_;
 
    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -f [m] [n]     Rocfrac preprocessing, optional NDA Data<m> & Grid<n> dirs\n";
      return;
    }

    if ($RocProps->runExtractOnly() || $RocProps->runAll()) {
       #Rocfrac: format -f 1 2  where 1 is the dataset and 2 is the grid
       $RocProps->setKeyValuePair("ROCFRAC", $TRUE);
       $RocProps->setKeyValuePair("ROCFRACDATA", "Data".shift(@$cline)); #integer input
       $RocProps->setKeyValuePair("ROCFRACGRID", "Grid".shift(@$cline)); #integer input
    }
    elsif  (($RocProps->runPreprocessOnly())||($RocProps->runCheckOnly())) {
       $RocProps->setKeyValuePair("ROCFRAC", $TRUE);
    }
}

sub s_command {
    my ($cline) = @_;
 
    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -s [m] [n]     Rocsolid preprocessing, optional NDA Data<m> & Grid<n> dirs\n";
      return;
    }

    if ($RocProps->runExtractOnly() || $RocProps->runAll()) {
	#Rocsolid: format -s 1 2  where 1 is the dataset and 2 is the grid
	$RocProps->setKeyValuePair("ROCSOLID", $TRUE);
	$RocProps->setKeyValuePair("ROCSOLIDDATA", "Data".shift(@$cline)); #integer input
	$RocProps->setKeyValuePair("ROCSOLIDGRID", "Grid".shift(@$cline)); #integer input
    }
    elsif  (($RocProps->runPreprocessOnly())||($RocProps->runCheckOnly())) {
        $RocProps->setKeyValuePair("ROCSOLID", $TRUE);
    }
}

sub i_command {
    my ($cline) = @_;
    my ($args);

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -i <o|u|f|s>   surfdive interface meshes, default infers from physics options\n";
      return;
    }

    #string does not contain only combinations of [oufs] min length 2, max length 4
    if ($$cline[0] !~ /^[oufs]{2,4}$/) {
       print "-i switch must be followed by valid combination of module letters [oufs]\n";
       printUsage();
       exit;
    }
    
    $args = shift(@$cline);
    if ($args =~ /[ufs]*o[ufs]*/) {
	if ($args =~ /[ous]*f[ous]*/) {
	    $RocProps->setKeyValuePair("ROCFLOROCFRAC",$TRUE);
	}
	if ($args =~ /[ouf]*s[ouf]*/) {
	    $RocProps->setKeyValuePair("ROCFLOROCSOLID",$TRUE);
	}
    }
    if ($args =~ /[ofs]*u[ofs]*/) {
	if ($args =~ /[ous]*f[ous]*/) {
	    $RocProps->setKeyValuePair("ROCFLUROCFRAC",$TRUE);
	}
	if ($args =~ /[ouf]*s[ouf]*/) {
	    $RocProps->setKeyValuePair("ROCFLUROCSOLID",$TRUE);
	}
    }
}


sub b_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -b             Rocburn preprocessing\n";
      return;
    }

    #process Rocburn
    $RocProps->setKeyValuePair("ROCBURN", $TRUE);
}

sub m_command {
  my ($cline) = @_;

  # if nothing is passed in, print the usage message for this option
  unless (defined($cline)) {
    print "  -m             Rocmop preprocessing\n";
    return;
  }

  #process Rocmop
  $RocProps->setKeyValuePair("ROCMOP", $TRUE);
}

sub A_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -A, --all         extract and preprocess\n";
      return;
    }

    $RocProps->setKeyValuePair("ALL", $TRUE);
    $RocProps->setKeyValuePair("ROCBURN", $FALSE);
}

sub C_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -C, --check       check an existing dataset at -d <path>\n";
      return;
    }

    $RocProps->setKeyValuePair("CHECKONLY", $TRUE);
}

sub E_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -E, --extract     copy NDA files to target at -t <path>\n";
      return;
    }

    $RocProps->setKeyValuePair("EXTRACTONLY", $TRUE);
}

sub P_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -P, --preprocess  run module preptools on data at -d <path>\n";
      return;
    }

    $RocProps->setKeyValuePair("PREPROCESS", $TRUE);

}


sub x_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -x, --ignore   ignore RocprepControl.txt control file\n";
      return;
    }

    $RocProps->setKeyValuePair("IGNOREFILE", $TRUE);
}


sub d_command {
    my ($cline) = @_; 

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -d <path>      path to source data, default is current working directory\n";
      return;
    }

    # make sure the source directory is a valid absolute path that's rwx
    my ($path) = File::Spec->rel2abs(shift(@$cline));
    $path = checkPathSlash($path);
    if (ModuleProcessor::checkDirExists($path)) {
       $RocProps->setKeyValuePair("SOURCEDIR", $path);
    }
    else {
       print "Directory is invalid or lacks rwx permission: $path\n";
       print "Please specify a valid directory with the -d switch.\n";
       printUsage();
       exit;
    }
}

sub h_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -h, --help     print this help message and terminate\n";
      return;
    }

    printUsage();
    exit;
}

sub n_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -n <m>         specify <m> processors/partitions\n";
      return;
    }

    my ($numprocs) = shift(@$cline);
    if (($numprocs !~ /^\d+$/) || ($numprocs <= 0)) {
	print "Number of processors must be an valid positive integer: $numprocs\n";
	print "Please specify an integer value with the -n switch\n";
        printUsage();
	exit;
    }

    #number of processors - integer follows switch
    $RocProps->setKeyValuePair("NUMPROCS", $numprocs); #integer input
}

sub plag_command {
  
  my ($cline) = @_;
  
  # if nothing is passed in, print the usage message for this option
  unless (defined($cline)) {
    print " -plag <bc> specify the injecting surface boundary condition\n";
    return;
  }

  my ($bc) = shift(@$cline);

  $RocProps->setKeyValuePair("PLAG",$bc);

}

sub t_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -t <path>      target path for new rocstar dataset\n";
      return;
    }

    # make sure the target directory is a valid absolute path
    my ($path) = File::Spec->rel2abs(shift(@$cline));
    $path = checkPathSlash($path);
    $RocProps->setKeyValuePair("TARGETDIR", $path);
}

sub p_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -p <path>      path to preptool binaries, default will use shell path\n";
      return;
    }

    # make sure the target directory is a valid absolute path
    my ($path) = File::Spec->rel2abs(shift(@$cline));
    $path = checkPathSlash($path);
    if (ModuleProcessor::checkDirExists($path)) {
       $RocProps->setKeyValuePair("BINDIR", $path);
    }
    else {
       print "Directory is invalid or lacks rwx permission: $path\n";
       print "Please specify a valid directory with the -p switch.\n";
       printUsage();
       exit;
    }
}

sub r_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -r <m>         specify <m> regions (rocflu only), default is -n value\n";
      return;
    }

    my ($numregions) = shift(@$cline);
    if (($numregions !~ /^\d+$/) || ($numregions <= 0)) {
	print "Number of regions must be an valid positive integer: $numregions\n";
	print "Please specify an integer value with the -r switch\n";
        printUsage();
	exit;
    }

    #number of regions - integer follows switch
    $RocProps->setKeyValuePair("NUMREGIONS", $numregions); #integer input
}

sub un_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -un <units>    convert model units to meters (rocfrac only)\n";
      return;
    }

    my ($units) = shift(@$cline);
    if (($units !~ /^\d*\.\d+$/) || ($units <= 0)) {
	print "Conversion to meters must be a valid positive number: $units\n";
	print "Please specify an new value with the -un switch\n";
        printUsage();
	exit;
    }

    #number of regions - integer follows switch
    $RocProps->setKeyValuePair("ROCFRACUNITS", $units);
}

sub split_command {
    my ($cline) = @_;

    # if nothing is passed in, print the usage message for this option
    unless (defined($cline)) {
      print "  -splitaxis <n> force split along n=0,1, or 2 axis (rocflo only)\n";
      return;
    }

    my ($axis) = shift(@$cline);
    if (($axis !~ /^\d+$/) && ($axis != 0) && ($axis != 1) && ($axis != 2)) {
	print "Makeflo split axis identifier must be 0,1, or 2: $axis\n";
	print "Please specify an new value with the -splitaxis switch\n";
        printUsage();
	exit;
    }

    #number of regions - integer follows switch
    $RocProps->setKeyValuePair("SPLITAXIS", $axis);
}

sub bad_command(){
    my ($arg) = shift;
    print "Unknown command line argument encountered: $arg\n\n";
}

#************************************************************************************
# print a useful usage message
#************************************************************************************
sub printUsage {
    print "****************************************************************************\n";
    print "Usage: $0 -A|C|E|P [OPTION]...\n";

    print "\nMajor modes of operation:\n";
    A_command();
    C_command();
    E_command();
    P_command();

    print "\nPhysics & service module selection:\n";
    o_command();
    u_command();
    f_command();
    s_command();
    b_command();
    m_command();

    print "\nModule-specific flags:\n";
    r_command();
    split_command();
    un_command();

    print "\nGeneral options:\n";
    i_command();
    d_command();
    h_command();
    n_command();
    t_command();
    p_command();
    x_command();

    print "\n";
    print "Example: $0 -A -o 1 1 -f 2 4 -d archiveDir/ -t newDataset/ -n 8\n";
    print "****************************************************************************\n\n";
}

#************************************************************************************
# place a trailing / on directory path if necessary, and trim any whitespace
#************************************************************************************
sub checkPathSlash {
    my ($path) = @_;

    if ($path !~ /\/\s*$/) {
        $path =~ s/\s*$/\//;
    }

    return ($path);
}
