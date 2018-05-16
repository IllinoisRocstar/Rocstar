#************************************************************************************
#  Class Properties

#  Superclass to all other properties classes. May be used as a generic 
#  properties class by using the getProps accessor method.  The only way to
#  add properties is through reading the properties file that is passed in 
#  as the second parameter to new. The new method then uses the allocate method
#  to read the file and fill a hash table with the key/value pairs. The properties
#  file has the form:

#  #->This configuration file consists of param name=param value pairs.
#  #->Comment lines must start with '#->'.  Parameter name/value pairs
#  #->must have a (->) between the name and value, with 
#  #->NO SPACES AROUND THE -> SIGN!
#  #->Example:
#  #->      paramname->/usr/bin/perl

#  The -> characters delimit the key and the value, as the last line shows. 
#************************************************************************************

package Properties;

use strict;

#************************************************************************************
#  Constructor for new properties objects.  Takes the properties filename as it's
#  parameter.  Calls allocate to do the real work. Then calls checkProps to make sure
#  that all the required properties are present. 
#************************************************************************************
sub new {
	my($pkg, $filename) = @_;
	my ($prop);
	
	$prop = $pkg->allocate($filename);
	
#	$pkg->checkProps($prop);
	
	return $prop;

}

#************************************************************************************
#  Reads the properties file supplied in $filename and fills an anonymous hash table
#  with the key/value pairs.  Then echos out the properties read and "blesses" the
#  hash with the appropriate class name (passed in through the $pkg variable).  
#************************************************************************************
sub allocate {
    my ($pkg, $fileName) = @_;
	my (%tempHash, @param, $key);
	
    #The following should probably place a warning in the Log, but we don't have access to
    #the log object here.... TBD

    if (defined($fileName) && ($fileName ne "")) {readFile(\%tempHash, $fileName, 1);}
	
    my $paramlist = \%tempHash;

    bless ($paramlist, $pkg);
	
    return $paramlist;
}


#************************************************************************************
#  Read Properties file upon explicit request. 
#************************************************************************************
sub readFile {
    my($props, $fileName, $replace) = @_;
    my($keyword, $value);
	
    open (JOBDATA, $fileName);

    # Fill hash table with parameter name/value pairs from a file.
    # if $replace=TRUE, values will always be overwritten, if $replace=FALSE then
    # hash values will only be stored if it is a fresh key name
    while(<JOBDATA>) {
        chomp;
        ($keyword, $value)= split(/->/,$_);
	# check that its not a comment or empty line
        if (($keyword ne '')&&($keyword !~ /^\s*\#/)) {
           if ($replace || !exists($$props{$keyword})) {
    	      $$props{$keyword}=$value;
           }
        }
    }
    close (JOBDATA);
}


#************************************************************************************
#  Returns the property requested by the $key parameter.  Returns null if property
#  doesn't exist. 
#************************************************************************************
sub getProp {
	my($properties, $key) = @_;
	
	return $$properties{$key};
}

#************************************************************************************
# write out all properties and their values to a message array and return it
#************************************************************************************
sub dumpProperties {
    my ($properties) = @_;
    my ($key, @messages);

    foreach $key (keys(%$properties)) {
        push(@messages, sprintf("%-18s= $$properties{$key}",$key));
    }

    return(sort @messages);
}

#************************************************************************************
#  Checks to see if the property name passed in as $key exists in this properties object. 
#************************************************************************************
sub propExists {
	my($properties, $key) = @_;
	#print "propExists:  properties = ". $properties . "\n";
	#print "propExists:  key = ". $key . "\n";
	
	return(exists($$properties{$key}));
}


#************************************************************************************
#  For debugging.  Prints the list of keys/values from a hash along with debug markers
#  and the name of the calling routine if supplied. Prints the hash reference too, so
#  you can see if you have the properties object you expected. 
#************************************************************************************
sub debugProperties {
	my ($prop, $caller) = @_;
	my ($key);
	print "\n***********DEBUG**************\n";
	print "Properties Hash supplied: " . $prop . "\n";
	print "Called from: " . $caller . "\n";
	foreach $key (keys(%$prop)) {
	    print "key = $key with value = $$prop{$key}\n";
	}
	print "***********DEBUG**************\n\n";

}


#************************************************************************************
#  Checks for missing properties in the properties file.  Uses a list of required
#  properties provided by the getPropsToCheck() method in the subclass. Note that 
#  the $pkg variable is a reference to the calling class, which will be a subclass
#  normally, and thus the getPropsToCheck() method from the subclass will be called.
#************************************************************************************
sub checkProps {
	my ($pkg, $props) = @_;
	my (@missing, $propList);
	my ($index, $key, $result);

	$index=0;
	$propList = $pkg->getPropsToCheck();
	foreach $key (@$propList) {
		$result = $props->propExists($key);
		if (!$result) {
			$missing[$index] = $key;
			$index=$index + 1;
		}
	}
	if ($index > 0) {
		printMissingProps($pkg, \@missing, $pkg->getPropfileName());
	}
	
	print "\n".$pkg->getPropfileName(). " Properties file OK. Continuing.\n\n";
}

#************************************************************************************
#  Prints a list of any missing properties in the properties file.  Uses a list of
#  missing properties passed in from the checkProps method. 
#************************************************************************************
sub printMissingProps {
	my ($pkg, $missing, $filename) = @_;
	my ($parameter);
	
	print "\nMISSING PROPERTIES IN $filename PROPERTIES FILE:\n";
	foreach $parameter (@$missing) {
        print "property = $parameter\n";
    }
	die "CAN'T CONTINUE WITHOUT THESE PROPERTIES!!!\n";
}

#************************************************************************************
#  Must be overridden by subclasses.  Returns the name of this type of properties
#  object to be used in messages. 
#************************************************************************************
sub getPropfileName {
	return "Generic Properties";
}

#************************************************************************************
#  Must be overridden by subclasses.  No checking performed if Properties class used
#  alone, since we don't know what properties will be there!
#************************************************************************************
sub getPropsToCheck {
	my (@empty) = ();
	return \@empty;
}


#************************************************************************************
#  Sets a key-value pair in the properties hash.  Adds if not there, changes if it is.
#  Send in key first, value second. 
#************************************************************************************
sub setKeyValuePair {
	my ($props, $theKey, $theValue) =  @_;

	$$props{$theKey} = $theValue;
}

1;
