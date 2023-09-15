#!/usr/bin/perl

#----------------------------------------------------------------
# File:     mergeOutput.pl
# Module:   tRIBS
# Author:   Sue Mniszewski 
# Created:  November 23, 2007
# Description: This script merges the individual output files
# generated from the parallel tRIBS runs.
#
# Input args are:
#   1. The input parameter file (ex. peach_run.in)
#
# Output:
#   1. Merged *.cntrl file containing reach information.
#   2. Merged *.time_00d file per selected interval containing dynamic 
#      data per node.
#   3. Merged  *.time-OOi containing integrated data at beginning and end.
#   4. Merged *_width file.
#   5. Merged *_area file.
#   6. Merged *_voi file. 
#
# To run:
#
#   mergeOutput.pl <parameter_file> 
#
# Example:
#
#   mergeOutput.pl peach_run.in
#
#----------------------------------------------------------------

if ($#ARGV != 0) {
    print "Incorrect # of arguments.\n";
    print "usage: mergeOutput.pl <parameter_file>\n";
    exit;
}

# Read args
$paramFile = shift(@ARGV);

print "Parameter file = $paramFile\n";

# Read in needed parameters for file names
sub readParameterFile() {
    open(PF, $paramFile)
        || die "readParameterFile: Parameterfile $paramFile is not found! \n";

    # Read for specific parameters
    while (<PF>) {
        chomp; # chop off <return>

        # Runtime or endtime
        if (/RUNTIME/) {
            chomp($runtime = <PF>);
             print "RUNTIME = $runtime\n";
        }

        # Dynamic file output interval
        if (/SPOPINTRVL/) {
            chomp($spopintrvl = <PF>);
            print "SPOPINTRVL = $spopintrvl\n";
        }

        # Output file name
        if (/OUTFILENAME/) {
           chomp( $outfilename = <PF>);
           print "OUTFILENAME = $outfilename\n";
        }

        # Output hydro file name
        if (/OUTHYDROFILENAME/) {
            chomp($outhydrofilename = <PF>);
            print "OUTHYDROFILENAME = $outhydrofilename\n";
        }

    }
    close(PF);
}


#
# Merge voronoi files
sub mergeVoronoiFiles() {
    # Dynamic files
#    $dsuffix = "_00d";
#    $dtime = 0;
#    while ($dtime <= $runtime) {
#        $otime = "$dtime";
#        while (length($otime)  < 4) {
#            $otime = "0$otime";
#        }
#        $dynamicFile = "$outfilename.$otime$dsuffix";
#
#        # Merge
#        print "Merging $dynamicFile.* ...\n";
#        sortMergeFile($dynamicFile, 1);
#
#        # Next time
#        $dtime += $spopintrvl;
#    }

    # Integrated files
    $dtime = 0;
    $isuffix = "_00i";
    while ($dtime <= $runtime) {
        $otime = "$dtime";
        while (length($otime)  < 4) {
            $otime = "0$otime";
        }
        $intFile = "$outfilename.$otime$isuffix";

        # Merge
        print "Merging $intFile.* ...\n";
        sortMergeFile($intFile, 1);

        # Next time
        $dtime += $spopintrvl;
    }

    # Width files
    $wsuffix = "_width";
    $widthFile = "$outfilename$wsuffix";
    print "Merging $widthFile.* ...\n";
    sortMergeFile($widthFile, 0);

    # Area files
    $asuffix = "_area";
    $areaFile = "$outfilename$asuffix";
    print "Merging $areaFile.* ...\n";
    sortMergeFile($areaFile, 0);

    # Voi files
    $vsuffix = "_voi";
    $voiFile = "$outfilename$vsuffix";
    $tsuffix = ".tmp";
    $tempfile = "$voiFile$tsuffix";
    $msuffix = ".merge";
    $mergefile = "$voiFile$msuffix";
    $first = 1;
    print "Merging $voiFile.* ...\n";
    # Merge multi-line records 
    while (defined ($nextf = <$voiFile.*>)) {
        #print "Merging $nextf\n";
        # First file gets copied to temp file
        if ($first == 1) {
            system("cp $nextf $tempfile");
            $first = 0;
        }
        else {
            # Open last merged file (temp file),
            # next file, and new merged file
            open(TF, $tempfile)
                || die "mergeVoronoiFiles: Last merged file $tempfile is not found! \n";         
            open(NF, $nextf)
                || die "mergeVoronoiFiles: Next file $nextf is not found! \n";
            open(MF, ">$mergefile")
                || die "mergeVoronoiFiles: Next merged file $mergefile is not created! \n";

            # Merge multi-line records from the last merged file and
            # the next file
            $tend = 0;
            $nend = 0;
            chomp($tbuffer = <TF>);
            @tfields = split(',',$tbuffer);
            chomp($nbuffer = <NF>);
            @nfields = split(',',$nbuffer);

            # Do merging till the end of one of the files has been 
            # encountered
            while ($tend != 1 && $nend != 1) {
                if ($tfields[0] <= $nfields[0]) {
                    # Write temp file record
                    print MF "$tbuffer\n";
                    while ($tbuffer ne "END") {
                        chomp($tbuffer = <TF>);
                        print MF "$tbuffer\n";
                    }
                    # Read first line of next record
                    if (not eof(TF)) {
                        chomp($tbuffer = <TF>);
                        @tfields = split(',',$tbuffer);
                        # Check for last 'end'
                        if ($tbuffer eq "END") {
                            $tend = 1;
                        }
                    }
                    # If end of file, set flag
                    else {
                        $tend = 1;
                    }
                }
                else {
                    # Write next file record
                    print MF "$nbuffer\n";
                    while ($nbuffer ne "END") {
                        chomp($nbuffer = <NF>);
                        print MF "$nbuffer\n";
                    }
                    # Read first line of next record
                    if (not eof(NF)) {
                        chomp($nbuffer = <NF>);
                        @nfields = split(',',$nbuffer);
                        # Check for last 'end'
                        if ($nbuffer eq "END") {
                            $nend = 1;
                        }
                    }
                    # If end of file, set flag
                    else {
                        $nend = 1;
                    }
                }
            }

            # One of the files has ended, write out the rest
            # of its lines to the merged file
            $endt = 0;
            if ($tend == 1){
                # It is the temp file
                # Finish writing out next file
                # Write out last line read
                print MF "$nbuffer\n";
                # Transfer rest of lines
                # Skip the last 'end'
                while (<NF>) {
                    print MF $_;
                }
            }
            else {
                # It is the next file
                # Finish writing out temp file
                # Write out last line read
                print MF "$tbuffer\n";
                # Transfer rest of lines
                # Skip the last end
                while (<TF>) {
                    print MF $_;
                }
            }
           
            # Make the merged file the temp file for the next time
            close(NF);
            close(MF);
            close(TF);
            system("mv $mergefile $tempfile");
        }
    } 

    # File is merged
    system("mv $tempfile $voiFile");
}

#
# Merge hyd files
sub mergeHydFiles() {

    $csuffix = ".cntrl";
    $rsuffix = ".reach";
    $cntrlFile = "$outhydrofilename$csuffix";
    $reachFile = "$outhydrofilename$rsuffix";
    print "Merging $cntrlFile.* ...\n";
    $rmin = 1; # Minimum reach number
    $rmax = 0; # Maximum reach number

    # Merge multi-line records 
    while (defined ($nextf = <$cntrlFile.*>)) {
       # Open each file, write out a reach's info to it's own file
       open(CN, $nextf)
        || die "mergeHydFiles: File $nextf could not be opened! \n";
       while (<CN>) {
          if (/REACH/) {
              @fields = split(" ", $_); 
              $reachNum = $fields[4];
              # Check for max/min
              if ($reachNum < $rmin) { $rmin = $reachNum; }
              if ($reachNum > $rmax) { $rmax = $reachNum; }
              # Open output file
              open(RE, ">$reachFile.$reachNum")        
                  || die "mergeHydFiles: Output file $reachFile.$reachNum could not be created! \n";  
              print RE $_; # Reach line
              # Read width, length, roughbness, slope, C, Y1, Y2, Y3  
              for ($i = 0; $i < 16; $i++) {
                   $buffer = <CN>;
                   print RE $buffer;
              }
              #print "Created $reachFile.$reachNum\n";
              close(RE);
           }
        }       
        close(CN);
    }

    # Concatenate all reach files together    
    system("mv $reachFile.$rmin $cntrlFile");
    for ($i = ($rmin+1); $i <= $rmax; $i++) {
        system("cat $reachFile.$i >> $cntrlFile");
    } 
    system ("rm $reachFile.*");
}

#
# Merge files and remove extra header copies
sub sortMergeFile() {
    my $mfile = shift;
    my $hascomma = shift;
    $tempfile = "$mfile-xxx";
    # Check if there is a comma delimiter
    $sortMergeLine = "";
    if ($hascomma == 1) {
        system("sort -mg -k 1 --field-separator=, $mfile.* > $tempfile");
    }
    # White space as delimiter
    else {
        system("sort -mg -k 1 $mfile.* > $tempfile");
    }

    # Remove extra headers
    removeExtraHeaders($mfile);
}

#
# Remove extra header copies
sub removeExtraHeaders() {
    my $mfile = shift;

    # Open output file
    open(MF, ">$mfile")
        || die "removeExtraHeaders: Output file $mfile could not be created! \n";
    # Open file with extra headers
    open(XF, "$mfile-xxx")
        || die "removeExtraHeaders: File $mfile.xxx could not be opened! \n";

    # Write first line
    $firstid = 1;
    $firstblank = 1;
    # Rewrite file without extra headers
    while (<XF>) {
        # Skip extra header lines
        if (/ID/) {
            # Only write first header line
            if ($firstid == 1) {
                print MF $_;
                $firstid = 0;
            }
        } 
        # Skip blank line
        elsif (/^\n/) {
            # Only write first blank line
            if ($firstblank == 1) {
                print MF $_;
                $firstblank = 0;
            } 
        }
        # Copy over other lines
        else {
             print MF $_;
        }
    }
    close(XF);
    close(MF);

    # Remove temporary file with extra headers
    system("rm $mfile-xxx");
}

# Read needed parameters
print "\nReading input parameters ...\n";
readParameterFile();

# Merge vononoi files
print "\nMerging voronoi files ...\n";
mergeVoronoiFiles();

# Merge hydfiles
print "\nMerging hyd files ...\n";
mergeHydFiles();
