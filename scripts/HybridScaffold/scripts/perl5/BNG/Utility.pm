# $Id: Utility.pm 4052 2015-08-17 22:37:01Z apang $

package BNG::Utility;

use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
sub dieLog;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
			readCMap
			writeCMapFile
			getContigStat
			getSubsetCMap
			getCMapIds
			shiftCMapIds
			countOverangLables
			readXMap
			writeXMapFile
			logMessage
			find_cmap_list_w_prefix
			parsingMrgStdOut
			parsingFastMrgStdOut
			writeAllMrgPairs
			writeAllFastMrgPairs
			getArrayFileNames
			uniqueIDs
			dieLog
			split_xmap_byRefContig
               );
our @EXPORT_OK = qw();
our $VERSION = 1.00;

##
# read a CMAP file
##
sub readCMap	{
	my ($cmap_file) = @_;
	# read cmap file (first time to see if there is the keyword "ChimQuality")
	open(IN, "$cmap_file") or dieLog ("ERROR: Unable to read in file $cmap_file: $!\n");
	my $foundChimQuality = 0;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /#h\s+/){
				$foundChimQuality = 1 if ($line =~ /ChimQuality/);
				last;	# done
			}
		} # if line
	} # while line
	close IN;
	# then read it again, will read in according whether ChimQuality is there or not
	
	my $cmap={};
	my $c_cmapId = 0;
	my $numcontig = 0;
	my @contigLength = ();
	$cmap->{"FileName"} = abs_path($cmap_file);
	my @header = ();
	my @data_name = ();
	my @data_type;
	my $cmap_version = "";
	open(IN, "$cmap_file") or dieLog ("ERROR: Unable to read in file $cmap_file: $!\n");
	my $NUM_C_LIMITED = ($foundChimQuality == 1) ? (10) : (9);	# assumption: column 10 has ChimQuality, if it is present
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /^#h\s+/)	{
				# data name:
				my $headTag = "";	# the #h
				($headTag, @data_name) = split(/\s+/, $line);
				splice @data_name, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_name[0]", @data_name[1..$#data_name]));
			} elsif ($line =~ /^#f\s+/)	{
				# data type:
				my $headTag = "";       # the #f
				($headTag, @data_type) = split(/\s+/, $line);
				splice @data_type, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_type[0]", @data_type[1..$#data_type]));
			} elsif ($line =~ /^#\s+CMAP File Version:/)	{
				my (@tmp) = split(/\s+/, $line);
				$cmap_version = $tmp[-1];
				push(@header, $line);
			} else	{
				push(@header, $line);
			} # if 
			next;
		} # if header line
		
		# data lines
		if (! exists $cmap->{"headers"})	{
			# first data line, populate the relevant hash values
			$cmap->{"headers"} = \@header;
			$cmap->{"dataName"} = \@data_name;
			$cmap->{"dataType"} = \@data_type;
			$cmap->{"MAPFileVersion"} = $cmap_version;
		} # if exists
		
		my $numc = $NUM_C_LIMITED;
		
		# now store information of that data line
		$line =~ s/^\s*//;	$line =~ s/\s+/\t/g;
		my @d_items = split(/\t/, $line);
		if ($d_items[0] != $c_cmapId)	{
			# a new contig:
			$numcontig += 1;
			$c_cmapId=$d_items[0];
			# ContigLength:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[1]} = $d_items[1];
			push(@contigLength, $d_items[1]);
			# NumSites:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[2]} = $d_items[2];
			
			# the rest of the columns
			for (my $i = 3; $i < $numc; $i += 1)	{
				my $ad = [];
				push(@$ad, $d_items[$i]);
				$cmap->{contigs}->{$c_cmapId}->{$data_name[$i]} = $ad;
			} # for i
		} else	{
			# another site of the same contig:
			for (my $i = 3; $i < $numc; $i += 1)	{
				my $ad = $cmap->{contigs}->{$c_cmapId}->{$data_name[$i]};
				push(@$ad, $d_items[$i]);
				$cmap->{contigs}->{$c_cmapId}->{$data_name[$i]} = $ad;
			} # for i
		} # if same cmap id
	} # while line
	close IN;
	$cmap->{nContigs} = $numcontig;
	return ($cmap, $numcontig, \@contigLength);
} # readCMap

# get all CMapIds in an array 
sub getCMapIds {
	my $cmap = shift;
	my @ids = keys %{ $cmap->{contigs} };
	return (\@ids);
}

sub writeCMapFile 
{
	my $cmap = shift;
	my $cmap_file = shift;
	my $want_full_header = shift;
	print "   write a cmap file $cmap_file ...\n";
	
	open CMF, ">$cmap_file" or dieLog ("ERROR: Unable to write to file " + $cmap_file + " - " + $!);	
	
	## headers:
	my @headers = @{ $cmap->{headers} };
	if ($want_full_header) { 
		for (my $i=0; $i<=$#headers; $i++) { 
			print CMF $headers[$i] . "\n";
		}
	} else { 
		# only header and format:
		print CMF "$headers[-2]\n";
		print CMF "$headers[-1]\n";
	}
	## contigs:
	my @data_name = @{ $cmap->{dataName} };
	#print "@data_name\n";
	my @cmapIds = sort keys %{ $cmap->{contigs} };
#	print "\n\n cmapid: @cmapIds\n\n";
#	print $cmap->{contigs}->{1}->{SiteID}->[2] . "\n";
#	exit;
	for (my $i=0; $i <= $#cmapIds; $i++) {
		my $c = $cmapIds[$i];
		my $numsites = $cmap->{contigs}->{$c}->{$data_name[2]};
		for (my $j=0; $j <= $numsites; $j++) { 
			print CMF $c . "\t" . $cmap->{contigs}->{$c}->{$data_name[1]} . "\t" . $numsites;
			for (my $m=3; $m <= $#data_name; $m++) {
				print CMF  "\t" . $cmap->{contigs}->{$c}->{$data_name[$m]}->[$j];
			}
			print CMF "\n";
		}
	}
	
	close CMF;
}

##
# create a new cmap with some cmap ids are excluded:
##
sub getSubsetCMap
{
	my $cmap = shift;
	my $givenCMapID_href = shift;
	my $flag = shift;
	if (!defined $flag) { 
		$flag = "exclude";
	}
	# also deal with list as an array ref
	if (ref $givenCMapID_href eq "ARRAY") { 
		my %givenCMapID = map { $_ => 1 } @$givenCMapID_href;
		$givenCMapID_href = \%givenCMapID;
	}
    my $allEID = join(",", keys %{$givenCMapID_href});	
	my $sCMap={};
	# - no header - $sCMap->{headers}
	print "   getting subset of cmap data from $cmap->{FileName}...\n";
    $sCMap->{dataName} = $cmap->{dataName};
    $sCMap->{dataType} = $cmap->{dataType};
    $sCMap->{MAPFileVersion} = $cmap->{MAPFileVersion};

    my @comment=("# This CMAP file was subset of $cmap->{FileName} with CMapID " . $flag .": " . $allEID);
    push(@comment, @{$cmap->{headers}});
    $sCMap->{headers} = \@comment;
#contigs:
    my $numc = 0;
    for my $contig ( keys %{ $cmap->{contigs} })   { 
    	#print "contig - $contig ";
    	if ($flag =~ /excl/) { 
    		# exclude given ids:
	    	if ($givenCMapID_href->{$contig}) { 
	    		# print "$contig excluded *******\n";
	    	} else { 
	    		#print "included\n";
	    		$sCMap->{contigs}->{$contig} = $cmap->{contigs}->{$contig};
	    		$numc++;
	    	}
    	} else { 
    		# include given ids only:
	    	if ($givenCMapID_href->{$contig}) {
	    	  # print "included\n";
	    		$sCMap->{contigs}->{$contig} = $cmap->{contigs}->{$contig};
	    		$numc++; 
	    	} else { 
	    		# print "$contig excluded *******\n";
	    	}    		
    	}
    }
    $sCMap->{nContigs} = $numc;  
    return $sCMap;
}

sub getContigStat
{
 # give a contig length array, caculate its statistical data:
    my $len = @_;
    if ($len == 0) {
    	return ("NA", "NA", "NA", "NA", "NA", "NA");
    }
    my @vals = sort {$b <=> $a} @_;
    my $median;
    my $r = $len%2;
    my $m = ($len-$r)/2;
    if($r != 0) #odd?
    {
        $median = $vals[$m];
        #print "\n$r, $m\n";
    }
    else #even
    {
        $median = ($vals[$m-1] + $vals[$m])/2;
        #print "\n$r, $m\n";
    }
    my $min = $vals[-1];
    my $max = $vals[0];
    my $total = 0.0;
    for (my $i=0; $i<$len; $i++) { 
    	$total = $total + $vals[$i];
    }
    my $mean = $total/($len);

#	N50 length is $n50: 
#   N50 value is $n50value: 
#   L50 is $L50:
    my $n50=0; 
    my $L50=0;
    my $n50value=0.0;
    if ($len==1) {
    	$n50value = $max;
    	$n50=$max;
    	$L50=1;
    } else {
	    for (my $i=0; $i<$len; $i++) {
	     $n50=$n50+$vals[$i];
	     $L50++;
	     if($n50 == $total/2.0){
	     	$n50value=$vals[$i];
	     	last; 
	     }
	     if($n50 > $total/2.0){
	     	$n50value=($vals[$i-1]+$vals[$i])/2.0;
	        last; 
	     }
	    }
    }
   	# print "N50 length is $n50 and N50 value is: $n50value and L50 is $L50\n"; 
    return ($min, $max, $mean, $median, $n50value, $total);
}

sub shiftCMapIds { 
	my $cmap = shift;
	my $nshift = shift;

	my @cmap_ids = sort {$b <=> $a} keys %{$cmap->{contigs}};
	for (my $i = 0; $i <= $#cmap_ids; $i++) { 
		my $contig_k = $cmap->{contigs}->{$cmap_ids[$i]};
		#delete the hash element:
		delete $cmap->{contigs}->{$cmap_ids[$i]};
		$cmap->{contigs}->{ $cmap_ids[$i] + $nshift } = $contig_k;
	}
}
##
# read xmap file:
##
sub readXMap {

	my $xmap_file =shift;
	print "   reading xmap file $xmap_file ...\n";
	open XMF, "<$xmap_file" or dieLog ("ERROR: Unable to read in file " + $xmap_file + " - " + $!);
	my $xmap={};
	$xmap->{"FileName"} = abs_path($xmap_file);
	my @header=();
	my @data_name;
	my @data_type;
	my $xmap_version;
	#processing header/comment lines:
	my $line=<XMF>;
	while ($line =~ /^#/) {
		chomp($line);
        	# bug fix in case sometime col-name is RefcontigID in some earlier verisons (should be RefContigID):
        	if ($line =~ /^#h\s+/) {
        		if ($line =~ /RefcontigID/) { 
        			$line =~ s/RefcontigID/RefContigID/;
        		}
        	}
	    	push (@header, $line);
		if ($line =~ /^#h\s+/) { 
			# data name:
			(undef, @data_name) = split(/\s+/, $line);
		}
		if ($line =~ /^#f\s+/) { 
			# data type:
			(undef, @data_type) = split(/\s+/, $line);
			last;
		}
		if ($line =~ /^#\s+XMAP File Version:/) {
			my (@tmp) = split(/\s+/, $line);
			$xmap_version = $tmp[-1];
		}
		$line=<XMF>;		
	}
	$xmap->{"headers"} = \@header;
	$xmap->{"dataName"} = \@data_name;
	$xmap->{"dataType"} = \@data_type;
	$xmap->{"MAPFileVersion"} = $xmap_version;

	#processing data:
	#XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum
	my $numc = @data_name;
	for (my $i=0; $i<$numc; $i++) { 
		$xmap->{hits}->{$data_name[$i]} = [];
	}
	$line=<XMF>;
	my $c_xmapId=0;
    my $numm = 0;
	#while ($line =~ /^\d/) {
	while (defined $line) {
		if($line =~ /^#/){
			$line = <XMF>;
			next;
		}
		# print $line;
		chomp($line);
		my (@d_items) = split(/\s+/, $line);
		$numm++;
        for (my $i=0; $i<$numc; $i++) {
        	my $a = $xmap->{hits}->{$data_name[$i]};
        	push(@$a, $d_items[$i]);
        	$xmap->{hits}->{$data_name[$i]} = $a;
        }
		$line=<XMF>;
	}
	$xmap->{totalHits} = $numm;
	close XMF;

	return ($xmap);	
}

##
# write a xmap data to a file
##
sub writeXMapFile {
	my $xmap_file = shift;
	my $xmap = shift;
	my $wantHeader = shift;
	my $rCmapFileName = shift;
	my $qCmapFileName = shift;

	print "   write xmap file $xmap_file ...\n";
	open XMF, ">$xmap_file" or dieLog ("ERROR: Unable to write in file " + $xmap_file + " - " + $!);

	if ($wantHeader == 1)	{
		# print the header lines
                my $tmp = $xmap->{headers};
		my $numLines = @$tmp;
		for (my $i = 0; $i < $numLines; $i++)	{
			my $line = $xmap->{headers}->[$i];
			if ($line =~ /^#\s+Reference/ && $rCmapFileName ne "")	{
				print XMF "# Reference Maps From:\t$rCmapFileName\n";
			} elsif ($line =~ /^#\s+Query/ && $qCmapFileName ne "")	{
				print XMF "# Query Maps From:\t$qCmapFileName\n";
			} else	{
				print XMF "$line\n";
			} # if line
		} # for i
	} # if wantHeader
	
    	my @data_name = @{ $xmap->{dataName} };
    	my $numc = @data_name;
	for (my $i=0; $i < $xmap->{totalHits}; $i++) { 
		for (my $j=0; $j < $numc; $j++) { 
	  		print XMF $xmap->{hits}->{$data_name[$j]}->[$i];
	  		if ($j < $numc - 1) { 
	  			print XMF "\t";
	  		} else { 
	  			print XMF "\n";
	  		}
	  	}
	}
	#print XMF "\n";
	close XMF;
}

##
# count number of lables extending outside a region:
##
sub countOverangLables { 
	my $cmap = shift;
	my $cmapId = shift;
	my $leftPosition = shift;
	my $rightPosition = shift;
	my @data_name = @{ $cmap->{dataName} };
	my $Position = $cmap->{contigs}->{$cmapId}->{Position};
	
#	if ($cmapId == 11 && $leftPosition== 271459.7) { 
#		print "\n\n\n\n @{ $Position } \n\n\n";
#	}
	# my $LabelChannel = $cmap->{contigs}->{$cmapId}->{LabelChannel};
	my $totalNS = $cmap->{contigs}->{$cmapId}->{NumSites};

	my $numL = 0;
	my $numR = 0;
	my $i;
	for ($i=0; $i<$totalNS; $i++) { 
		if ( ($Position->[$i]) < $leftPosition ) { 
			$numL++;
		};
		if ( ($Position->[$i]) > $rightPosition ) { 
			$numR++;
		};
	}
#	if ($cmapId == 11){print "\ncohl: $cmapId, $leftPosition,$rightPosition, $numL, $numR\n"};
	#exit;
	return ($numL, $numR);
}

sub logMessage{
	my $fh = shift;
	my $msg1 = shift;
	my $msg2 = shift;
	print $fh "$msg1 " . "$msg2\n";	
}

sub find_cmap_list_w_prefix {
	my $prefix = shift;
	my @cmap_files = glob("$prefix" . "_contig[0-9]*.cmap" );
	return (@cmap_files);
}

sub parsingMrgStdOut { 
	my $stdout_file = shift; 
	open SOFH, "<$stdout_file" or dieLog ("ERROR: Unable to read in file " + $stdout_file + " - " + $!);
	my $merge_pairs = {};
	my @ContigID1=();
	my @ContigID2=();
	my @ResultContigID=();
	my $numMrgP = 0;
	while (my $line = <SOFH>) {
		chomp($line);
#       parsing out three id numbers from these lines:
#          "Creating merged based on alignment between ..." 		
        if ($line =~ m/^Creating merged map based on alignment between .*\(id=(\d+),len=.*\(id=(\d+),len=.*\:id=(\d+),len=.*/) { 
	        push(@ContigID1, $1);
	        push(@ContigID2, $2);
	        push(@ResultContigID, $3);
	        $numMrgP++;
        }
	}
	close(SOFH);
	$merge_pairs->{ContigID1} = \@ContigID1;
	$merge_pairs->{ContigID2} = \@ContigID2;
	$merge_pairs->{ResultContigID} = \@ResultContigID;
	$merge_pairs->{numMrgPairs} = $numMrgP;
	return($numMrgP, $merge_pairs);
}

sub parsingFastMrgStdOut	{
	my $stdout_file = shift;
	open (IN, $stdout_file) or dieLog ("ERROR: Unable to read in file " + $stdout_file + " - " + $!);
	my $merge_pairs = {};
	my @contigId1 = ();
	my @contigId2 = ();
	my @resultsContigId = ();
	my @mergeRounds = ();
	my $numMrgP = 0;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		my ($theStage, $firstId, $firstSize, $secondId, $secondSize) = (-1, -1, -1, -1, -1);
		my $hybridId = -1;
		# parsing out from these lines (\d+):Created merged map based on alignment between .*\(id=(\d+),len=.*\(id=(\d+),len=.*\:id=(\d+),len=.*
		if ($line =~ /^(\d+).+Creat.* merged map based on alignment between .*\(id=(\d+),len=(\d+\.\d*)\).*\(id=(\d+),len=(\d+\.\d*)\).+(Just keeping Xmap)?/)	{
			($theStage, $firstId, $firstSize, $secondId, $secondSize) = ($1, $2, $3, $4, $5);
			if ($line =~ /Just keeping Xmap/)       {
				# one of the maps is completely encompassed by another map, then the larger map's id is the final hybrid id
				$hybridId = ($firstSize < $secondSize) ? ($secondId) : ($firstId);
			} else	{
				# hybrid Id is the smaller of the two ids
				$hybridId = ($firstId < $secondId) ? ($firstId) : ($secondId);
			} # if line has that string
			push(@mergeRounds, $theStage);
			push(@contigId1, $firstId);
			push(@contigId2, $secondId);
			push(@resultsContigId, $hybridId);
			$numMrgP += 1;
		} # if line
	} # while line
	close IN;
	$merge_pairs->{ContigID1} = \@contigId1;
	$merge_pairs->{ContigID2} = \@contigId2;
	$merge_pairs->{ResultContigID} = \@resultsContigId;
	$merge_pairs->{numMrgPairs} = $numMrgP;
	$merge_pairs->{mergeRounds} = \@mergeRounds;
	return ($numMrgP, $merge_pairs);
} # parsingFastMrgStdOut

# writeAllMrgPairs($all_step1_mrg_pairs, $scratchDir . "step1.merge.pairs.txt", $rounds, \@mrg_rounds_ids);
sub writeAllMrgPairs { 
	my $all_step1_mrg_pairs = shift;
	my $fn = shift;
	my $final_round_num = shift;
	my $mrgRoundIds_ref = shift;
    open TMPH, ">$fn" or dieLog ("ERROR: Unable to write to file " + $fn + " - " + $!);
    print TMPH '"RowID" "ContigID1" "ContigID2" "ResultContigID" "MrgRound"' . "\n";
    my $n =0;
    for (my $i=1; $i<=$final_round_num; $i++) {
      #print  "$final_round_num, $i, mergeCnt=$all_step1_mrg_pairs->{$i}->{mergeCnt}\n";
      for (my $j=0; $j < $all_step1_mrg_pairs->{$i}->{mergeCnt}; $j++) {
      	$n++;      	
    	print TMPH '"' . $n . '" "' . $all_step1_mrg_pairs->{$i}->{mergePairs}->{ContigID1}->[$j] . '" "' . 
    	    $all_step1_mrg_pairs->{$i}->{mergePairs}->{ContigID2}->[$j] . '" "' . 
    	    $all_step1_mrg_pairs->{$i}->{mergePairs}->{ResultContigID}->[$j] . '" "Mrg' .
    	    $mrgRoundIds_ref->[$i] . "\"\n";
      }
    }
    close(TMPH);	
}

sub writeAllFastMrgPairs	{
	my ($fn, $mrg_pairs, $numMrgP, $mrgRoundsIdsRef) = @_;
	open(OUT, ">$fn") or dieLog("ERROR: Unable to write to file " + $fn + " - " + $!);
	print OUT '"RowID" "ContigID1" "ContigID2" "ResultContigID" "MrgRound"'."\n";
	for (my $i = 0; $i < $numMrgP; $i += 1)	{
		my $letterMrgRound = ($mrg_pairs->{mergeRounds}[$i] >= scalar(@$mrgRoundsIdsRef)) ? ($mrg_pairs->{mergeRounds}[$i]) : ($mrgRoundsIdsRef->[$mrg_pairs->{mergeRounds}[$i] - 1]);	
		print OUT '"'.($i + 1).'" "'.$mrg_pairs->{ContigID1}[$i].'" "'.
			$mrg_pairs->{ContigID2}[$i].'" "'.
			$mrg_pairs->{ResultContigID}[$i].'" "Mrg'.
			$letterMrgRound."\"\n";	
	} # for i
	close OUT;
} # writeAllFastMrgPairs

sub getArrayFileNames { 
	my $prefix = shift;
	my $ids_ref = shift;
	my $suffix = shift;
	my @fns=();
	foreach (@{$ids_ref} ) { 
		my $fn = $prefix . $_ . $suffix;
		push(@fns, $fn);
	}
	return \@fns;
}

sub uniqueIDs{ 
	my $firstA = shift;
	my $secondA = shift;
	my %second = map {$_=>1} @{$secondA};
    my @only_in_firstA = grep { !$second{$_} } @{$firstA};
    return (\@only_in_firstA);
}

sub dieLog{
	my $string = shift;
	print "$string";
	die "$string";
	exit;
}

sub split_xmap_byRefContig {
	my ($xmap, $rcmap, $qcmap, $contig_prefix, $molPath) = @_;
	 	
	my @xmap_in = read_file($xmap);	
	my @rcmap_in = read_file($rcmap);
	my @qcmap_in = read_file($qcmap);
	
	my %rcmap;
	my $rcmap_header="";
	my %qcmap;
	my $qcmap_header="";
	
	my @xmapheader;
	
	print "\nSplitting $xmap...\n";
	
	for (my $i=0; $i <= 4; $i++) {
		push (@xmapheader, $xmap_in[$i]); }
		
	eval { mkpath($molPath) };
	if ($@) {
		dieLog ("ERROR: Couldn't create $molPath: $@"); }
		
	
	foreach my $rline (@rcmap_in) {
	if ($rline =~ /^#/ ) {
		$rcmap_header = $rcmap_header . $rline; }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$Position,$StdDev,$Coverage,$Occurrence,$GmeanSNR,$lnSNRsd) = split("\t",$rline);
		if ($CMapId) {
			push @{ $rcmap{$CMapId} }, $rline; }}}
			
	foreach my $qline (@qcmap_in) {
	if ($qline =~ /^#/ ) {
		$qcmap_header = $qcmap_header . $qline;  }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$PositionStdDev,$Coverage,$Occurrence) = split("\t",$qline);
		if ($CMapId) {
			push @{ $qcmap{$CMapId} }, $qline; }}}
	
	
	my $prevContigID = 0;
	foreach my $xline (@xmap_in) {
		if ($xline !~ /^#/ ) {			
			my ($XmapEntryID, $QryContigID, $RefContigID, $QryStartPos, $QryEndPos, $RefStartPos, $RefEndPos, $Orientation, $Confidence, $HitEnum) = split("\t",$xline);  
			if ($RefContigID != $prevContigID) {
				
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				my $rcmapName = $contig_prefix.$RefContigID."_r.cmap";
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (RCMAP, ">>$molPath/$rcmapName");
				open (QCMAP, ">>$molPath/$qcmapName");
				
				foreach my $s (@xmapheader) {
					$s =~ s/^\s+//;
					print XMAP $s; }
				print XMAP "# Reference Maps From:	".$rcmapName."\n";
				print XMAP "# Query Maps From:	".$qcmapName."\n";
				print XMAP "#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum\n";
				print XMAP "#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string\n"; 
				print XMAP $xline;
				

				
				# RCMAP output
				print RCMAP $rcmap_header; 
				if($rcmap{$RefContigID}){
					print RCMAP join("", @{ $rcmap{$RefContigID} }); }
				
				#QCMAP output
				print QCMAP $qcmap_header;
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }
				
				$prevContigID = $RefContigID;
					
				close(QCMAP);
				close(RCMAP);
				close(XMAP);
				}				
				
			elsif ($RefContigID == $prevContigID) {
				
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (QCMAP, ">>$molPath/$qcmapName");				

				print XMAP $xline; 
				
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }

				close(QCMAP);
				close(XMAP);				
				
				} } }
	
	print "Splitting $xmap complete.\n"; 	
}

1;
