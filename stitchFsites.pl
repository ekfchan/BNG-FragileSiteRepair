#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long qw(:config bundling); 
use Cwd; 
use File::Basename;
use File::Copy;
use File::Find; 
use Sys::Info;
# use Sys::Info::Constants qw( :device_cpu );
# use Sys::MemInfo qw(totalmem freemem totalswap);
use Sys::MemInfo qw(totalmem);
use List::Util qw( max );	#qw( min max )
use FindBin;	#locate this script 
use lib "$FindBin::Bin";
use stitchFsites qw(mergeContigs doalignment appendStitched);
use summarise qw(getCmapIds calcCmapStats getFsiteStats); 

sub makeNewID; 
sub unique;

# my $scriptspath = Cwd::abs_path(dirname($0));#library path 
my $scriptspath = $FindBin::Bin; 

print "\n";
print qx/ps -o args $$/;
print "\n";

## << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'xmap|x=s', 'qcmap|q=s', 'rcmap|r=s', 'errbin|e=s', 'output|o=s', 'bed|b=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'n:i', 'round:i', 'verbose'); 

if( !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{errbin} | !exists $inputs{output} | !exists $inputs{bed} ) {
	print "Usage: perl stitchFsites.pl --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output_folder> --bed <.bed fragile sites file> [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>] [--n <CPU cores to use>]\n"; 
	exit 0; 
}
if( exists $inputs{verbose} ) {
	print "Scripts and modules:  $scriptspath \n"; 
	foreach my $key ("xmap","qcmap","rcmap","errbin","output","bed") { 
		$inputs{$key} = Cwd::abs_path($inputs{$key}); 
		if( !-e $inputs{$key} ) { die "Input file $inputs{$key} does not exist: $$!\n"; }
		print "--",$key, ": ", $inputs{$key}, "\n"; 
	}
}

## Default values
my $maxBp = 10000; 
if( exists $inputs{maxfill} ) { $maxBp = $inputs{maxfill}; } else { $inputs{maxfill} = $maxBp; }
my $wobble = 250; 
if( exists $inputs{wobble} ) { $wobble = $inputs{wobble}; } else { $inputs{wobble} = $wobble; }
my $maxlab = 0; 
if( exists $inputs{maxlab} ) { $maxlab = $inputs{maxlab}; } else { $inputs{maxlab} = $maxlab; }
foreach my $opt ("maxfill", "maxlab", "wobble") { print "--",$opt," : ",$inputs{$opt},"\n"; }
if( !exists $inputs{round} ) { $inputs{round} = 1; }

# my $idOffset = 100000;

#get number of CPU cores
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );	#Number of available CPUs
my $cpuCount = $cpu->count;
if( exists $inputs{n} ) { $cpuCount = $inputs{n}; }
my $mem = (((&totalmem / 1024) / 1024)) / 1024; #Available RAM

# Define prefix
my $prefix = basename($inputs{rcmap}); 	#the reference never changes
$prefix =~ s/_r.cmap//; 
if( exists $inputs{verbose} ) { print "Output file prefix: $prefix \n"; }

## << read input reference CMAP >> 
print "\n\nReading in reference CMAP..\n"; 
my $r_sites=0;
my %r_cmap;
open RCMAP, $inputs{rcmap} or die $!;
while (my $line = <RCMAP>) {
	chomp($line);
	if ($line =~ /^#/) {	#skip header
		next; }
	else {
		$r_sites++;
		my @s = split(/\t/,$line);
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	
		$r_cmap{"$s[0]"}{"ContigLength"} = "$s[1]";
		$r_cmap{"$s[0]"}{"NumSites"} = "$s[2]";
		$r_cmap{"$s[0]"}{"$s[5]"}{"SiteID"} = "$s[3]";
		# $r_cmap{"$s[0]"}{"$s[5]"}{"LabelChannel"} = "$s[4]";
		# $r_cmap{"$s[0]"}{"$s[5]"}{"TheRest"} = join("\t",@s[6..$#s]);
	}
}
close RCMAP; 
print "\tRead in $r_sites sites from $inputs{rcmap}\n";
print "\tRead in ".scalar(keys %r_cmap)." maps from $inputs{rcmap}\n";
if (scalar(keys %r_cmap) > 1) {
	die "\tInput _r.cmap should have only 1 map\n"; }
my $RefContigID = (keys %r_cmap)[0]; 

## << read input XMAP >> 
## Theoretically, allow XMAP to have >1 RefContigID, but will only keep alignments where RefContigID matching that one in RCMAP. 
print "\n\nReading in reference XMAP..\n"; 
my @xmap; 
open XMAP, $inputs{xmap} or die $!;
while (my $line = <XMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum
		#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 
		# 1	14839	17	344.2	999946.6	42640026.3	43628302.0	+	138.45	2M1D4M1D8M1D10M1D7M1D7M1D2M1D6M1D2M1D2M1I10M2D3M1D13M2D6M1D10M1D10M1D6M5I4M
		if( $s[2] eq $RefContigID ) {
			my %xmap_line = (
					"XmapEntryID"  => "$s[0]", 
					"QryContigID" => "$s[1]", 
					# "RefContigID"  => "$s[2]",	#should be identical (one of the assumptions)
					"QryStartPos"  => "$s[3]", 
					"QryEndPos"  => "$s[4]", 
					"RefStartPos"  => "$s[5]", 
					"RefEndPos" => "$s[6]", 
					"Orientation" => "$s[7]"
					# "Orientation" => "$s[7]", 
					# "Confidence" => "$s[8]", 
					# "HitEnum" => "$s[9]"
				); 
			push @xmap, \%xmap_line; }
		}
}
close XMAP; 
print "\tRead in ", scalar(@xmap), " alignments from $inputs{xmap}\n";
if( scalar(@xmap) == 0 ) { die "\tThere are no alignments matching $RefContigID (in $inputs{rcmap})\n"; } 
@xmap =  sort { $a->{'RefStartPos'} <=> $b->{'RefStartPos'} } @xmap;

## << read input query CMAP >> 
print "\n\nReading in query CMAP..\n"; 
my $q_sites=0;
my %q_cmap; 
my $qmapHeader; 
open QCMAP, $inputs{qcmap} or die $!;
while (my $line = <QCMAP>) {
	chomp($line);
	if ($line =~ /^#/) {
		if( $line =~ /h\sCMapId\sContigLength/ ) { $qmapHeader = $line; }
		next; }	#skip header
	else {
		$q_sites++;
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	GmeanSNR	lnSNRsd	SNR
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	12.4650	0.4814	0
		$q_cmap{"$s[0]"}{"ContigLength"} = "$s[1]";
		$q_cmap{"$s[0]"}{"NumSites"} = "$s[2]";
		$q_cmap{"$s[0]"}{"$s[5]"}{"SiteID"} = "$s[3]";
		$q_cmap{"$s[0]"}{"$s[5]"}{"LabelChannel"} = "$s[4]";
		$q_cmap{"$s[0]"}{"$s[5]"}{"TheRest"} = join("\t",@s[6..$#s]);
	}
}
close QCMAP; 
print "\tRead in $q_sites sites from $inputs{qcmap}\n";
print "\tRead in ".scalar(keys %q_cmap)." maps from $inputs{qcmap}\n";
my $idOffset = "1".("0" x (length(max (keys %q_cmap))-1)); 

## << fsite bed >> 
# Filter fsites.bed and keep only those corresp to contig in _r.cmap
# Assumption:  FSITEs must be provided and only gaps co-locating with fsites are merged. 
my @fsites;
my $fragileSites = 0;
print "\n\nFiltering fragile site BED list..\n";

open FSITES, $inputs{bed} or die $!;
while (my $line = <FSITES>) {
	chomp($line);
	if ($line =~ /^#/) {#if header then skip
		next; }
	else {
		my @s = split(/\t/,$line);
		if ($s[0] eq $RefContigID) {
			$fragileSites++;
			#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
			my %fsites_line = (
				# "CMapId"  => "$s[0]",	#CMapId must equal $r_mapIds[0] (and scalar($r_mapIds) should be 1)
				"Start" => "$s[1]", 
				"End"  => "$s[2]",
				"Type"  => "$s[3]",
				# "Score" => "$s[4]",
				# "Strand" => "$s[5]",
				# "ThickStart" => "$s[6]",
				# "ThickEnd" => "$s[7]",
				# "ItemRgba" => "$s[8]",
				"Sequence" => "$s[9]",		#??? DO WE NEED THIS HERE??  
				"line" => $line
			);
			push @fsites, \%fsites_line;
		}
	}
}
close FSITES; 
@fsites =  sort { $a->{'Start'} <=> $b->{'Start'} } @fsites;	# @fsites sorted by start position, assuming end position is always >= start positions
print "\t", scalar(@fsites)." fragile sites identified in $inputs{bed} \n";


## << Initial Stitched BED file >> 
# For completenes, we want a stitched file created even if there are no contigs can be merged (the file will only have a header) 
my $outbed = $inputs{output}."/".$prefix."_fragileSiteRepaired_stitchPositions.bed";
print "\n\nStitched fragile sites will be record in $outbed\n"; 
if (-e $outbed) { print "\tWARNING:: File already exists and will be overwritten\n"; }
#print header 
if( $inputs{round} == 1 && !-e $outbed ) {
	open (OUT, ">", $outbed) or die "\tERROR: Cannot open $outbed for writing! $!\n";
	print OUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n"; 
	close OUT; 
} 

## << Merge Contigs >> 

my @mergedcmaps; #maps that have been merged and will be removed
my $didmerge = 0; 

print "\n +++++++++++++++++ Round $inputs{round} +++++++++++++++++ \n"; 
MAIN: for (my $i=0; $i < scalar(@xmap); $i++) {
	my $id1 = $xmap[$i]->{'XmapEntryID'};
	my $firstRefStartPos = $xmap[$i]->{'RefStartPos'};
	my $firstRefEndPos = $xmap[$i]->{'RefEndPos'};
	my $firstQryContigID = $xmap[$i]->{'QryContigID'};
	my $firstOrientation = $xmap[$i]->{'Orientation'};
	my $firstQryStartPos = $xmap[$i]->{'QryStartPos'};
	my $firstQryEndPos = $xmap[$i]->{'QryEndPos'};

	if (($i+1)<scalar(@xmap)) {	
		my $id2 = $xmap[$i+1]->{'XmapEntryID'};
		my $secondRefStartPos = $xmap[$i+1]->{'RefStartPos'};
		my $secondRefEndPos = $xmap[$i+1]->{'RefEndPos'};
		my $secondQryContigID = $xmap[$i+1]->{'QryContigID'}; 
		my $secondOrientation = $xmap[$i+1]->{'Orientation'};
		my $secondQryStartPos = $xmap[$i+1]->{'QryStartPos'};
		my $secondQryEndPos = $xmap[$i+1]->{'QryEndPos'};

		print "\nConsidering merge between XmapEntryID: $id1 and XmapEntryID: $id2 \n"; 
		print "\t1st Qry: $firstQryContigID\t2nd Qry: $secondQryContigID\n"; 
		if( grep($_ eq $firstQryContigID, @mergedcmaps) ) { 
			print "\tContig $firstQryContigID already merged\n"; 
			next MAIN; 
		}
		print "\tDistance: ", abs($secondRefStartPos - $firstRefEndPos), " bp \n";
		
		# Contig pair must not overlap on RefContig
		if (($secondRefStartPos >= $firstRefEndPos) && ($secondRefEndPos >= $firstRefEndPos) && ($secondRefStartPos >= $firstRefStartPos)) {
			# NOTE::: Assuming StartPos <= EndPost in XMAPs
			print "\tOverlap filter: PASS\n";
			
			# Gap distance must not excees $maxBp
			# NOTE:: Is the abs() necessary??? ASSUMPTION:: RefContig in same orientation in all alignments (i.e. throughput XMAP)
			if (abs($secondRefStartPos - $firstRefEndPos) <= $maxBp) {
				print "\tDistance filter: PASS\n";
				
				# Alignment must extend to start/end of contig at junction site
				if ( (($firstOrientation eq "+" && $secondOrientation eq "+") && ($q_cmap{$firstQryContigID}{$firstQryEndPos}{"SiteID"} eq $q_cmap{$firstQryContigID}{"NumSites"}) && ($q_cmap{$secondQryContigID}{$secondQryStartPos}{"SiteID"} eq 1)) || (($firstOrientation eq "+" && $secondOrientation eq "-") && ($q_cmap{$firstQryContigID}{$firstQryEndPos}{"SiteID"} eq $q_cmap{$firstQryContigID}{"NumSites"}) && ($q_cmap{$secondQryContigID}{$secondQryStartPos}{"SiteID"} eq $q_cmap{$secondQryContigID}{"NumSites"})) || (($firstOrientation eq "-" && $secondOrientation eq "+") && ($q_cmap{$firstQryContigID}{$firstQryEndPos}{"SiteID"} eq 1) && ($q_cmap{$secondQryContigID}{$secondQryStartPos}{"SiteID"} eq 1)) || (($firstOrientation eq "-" && $secondOrientation eq "-") && ($q_cmap{$firstQryContigID}{$firstQryEndPos}{"SiteID"} eq 1) && ($q_cmap{$secondQryContigID}{$secondQryStartPos}{"SiteID"} eq $q_cmap{$secondQryContigID}{"NumSites"})) ) {
					print "\tAlignment filter: PASS\n";
					
					# Number of labels on RefContigID in gap must be <= $maxlab
					# ASSUMPTION::  All RefContigs in same orientation for all alignments
					my $firstRefEndPosSite = $r_cmap{$RefContigID}{$firstRefEndPos}{"SiteID"} or die "Can't find site of RefContig $RefContigID} at $firstRefEndPos \n";
					my $secondRefStartPosSite = $r_cmap{$RefContigID}{$firstRefEndPos}{"SiteID"} or die "Can't find site of RefContig $RefContigID at $firstRefEndPos \n";
					my $labelsDistance = $secondRefStartPosSite - $firstRefEndPosSite - 1;	#-1 means start and stop at same label
					print "\tLabels: $labelsDistance\n";

					if ($labelsDistance <= $maxlab) {
						print "\t\tLabel filter: PASS\n";				

						# Only stitch contigs if they overlap with fsite. 
						my $loBound = $firstRefEndPos - $wobble; 
						my $upBound = $secondRefStartPos + $wobble; 
						my $fsiteFound; 
						print "\tLooking for fragile site within window corrected for wobble: $loBound to $upBound, Range: ".($upBound - $loBound)."\n";
						foreach my $hash (@fsites) {
							#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
							my $start = $hash->{'Start'};
							my $end = $hash->{'End'};
							my $type = $hash->{'Type'};
							my $seq = $hash->{'Sequence'};
							my $fsiteline = $hash->{'line'};
							if( $start >= $loBound && $end <= $upBound ) {
								print "\t\tFsite found: Start: $start End: $end Range: ".abs($end - $start)."\n";
								$fsiteFound = $fsiteline; 
								last; 
							}
						}
						if ( (!defined $fsiteFound) or ($fsiteFound eq 0) ) { print "\t\tNo fragile site found in Window\n"; }
						else {
							## Merge contig pair
							my $newCMapId = makeNewID($firstQryContigID, $idOffset); 
							print "\tNew Query ID: $newCMapId\n"; 
							if ( exists $inputs{verbose} ) { 
								print "\t\tFirst alignment: $id1 Orientation: $firstOrientation Query Contig:: $firstQryContigID align from $firstQryStartPos (site # $q_cmap{$firstQryContigID}{$firstQryStartPos}{'SiteID'}) to $firstQryEndPos (site # $q_cmap{$firstQryContigID}{$firstQryEndPos}{'SiteID'}) of $q_cmap{$firstQryContigID}{'NumSites'} total sites\n";
								print "\t\tSecond alignment: $id2 Orientation: $secondOrientation Query Contig:: $secondQryContigID align from $secondQryStartPos (site # $q_cmap{$secondQryContigID}{$secondQryStartPos}{'SiteID'}) to $secondQryEndPos (site # $q_cmap{$secondQryContigID}{$secondQryEndPos}{'SiteID'}) of $q_cmap{$secondQryContigID}{'NumSites'} total sites\n";
							}
							my $newCMapRef = mergeContigs(\%{$q_cmap{$firstQryContigID}}, \%{$q_cmap{$secondQryContigID}}, $firstOrientation, $secondOrientation, abs($secondRefStartPos - $firstRefEndPos), (exists $inputs{verbose})); 
							
							# add newly merged map to %q_cmap
							foreach my $key (keys $newCMapRef) { 
								if( $key eq 'ContigLength' || $key eq 'NumSites' ) { $q_cmap{$newCMapId}{$key} = $newCMapRef->{$key}; } 
								else {
									$q_cmap{$newCMapId}{$key}{'SiteID'} = $newCMapRef->{$key}{'SiteID'}; 
									$q_cmap{$newCMapId}{$key}{'LabelChannel'} = $newCMapRef->{$key}{'LabelChannel'}; 
									$q_cmap{$newCMapId}{$key}{'TheRest'} = $newCMapRef->{$key}{'TheRest'}; 
								}
							}
							push @mergedcmaps, $firstQryContigID;
							push @mergedcmaps, $secondQryContigID; 

							$didmerge = 1; 
							# print "\t\tMerging successful!!\n"; 

							# print stitched position 
							print "\n\tAppending stitched fragile site to $outbed ...\n"; 
							appendStitched($outbed, $RefContigID, int($firstRefEndPos-1), int($secondRefStartPos+1), $fsiteFound); 
							# printstitched($outbed, $RefContigID, int($firstRefEndPos-1), int($secondRefStartPos+1), $fsiteFound); 
						}
					} else { print "\t\tLabel filter: FAIL\n"; }
				}
				else {
					print "\tAlignment filter: FAIL\n";
					print "\t\tFirst alignment: $id1 Orientation: $firstOrientation Query Contig:: $firstQryContigID align from $firstQryStartPos (site # $q_cmap{$firstQryContigID}{$firstQryStartPos}{'SiteID'}) to $firstQryEndPos (site # $q_cmap{$firstQryContigID}{$firstQryEndPos}{'SiteID'}) of $q_cmap{$firstQryContigID}{'NumSites'} total sites\n";
					print "\t\tSecond alignment: $id2 Orientation: $secondOrientation Query Contig:: $secondQryContigID align from $secondQryStartPos (site # $q_cmap{$secondQryContigID}{$secondQryStartPos}{'SiteID'}) to $secondQryEndPos (site # $q_cmap{$secondQryContigID}{$secondQryEndPos}{'SiteID'}) of $q_cmap{$secondQryContigID}{'NumSites'} total sites\n";
				}
			} else { print "\tDistance filter: FAIL\n"; } 
		} else { print "\tOverlap filter: FAIL\n"; } 
	}
}

# << write current q_cmap >> 
my @cmapstoprint = grep ! ${{map{$_,1} @mergedcmaps}}{$_}, keys %q_cmap;
@cmapstoprint =  sort {$a <=> $b} @cmapstoprint;
if ( exists $inputs{verbose} ) { 
	print "\n***\n"; 
	print "All cmaps now: ", join("\t", keys %q_cmap), "\n"; 
	print "\nMerged maps: ", join("\t", @mergedcmaps), "\n"; 
	print "\nFinal cmaps: ", join("\t", @cmapstoprint), "\n"; 
}

my $tmpprefix = $inputs{output}."/".$prefix."_tmp".$inputs{round}; 
open( OUT, ">$tmpprefix.cmap" ) or die $!; 
print OUT "#Round $inputs{round} of fragile site repair\n"; 
print OUT $qmapHeader, "\n"; 
##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	GmeanSNR	lnSNRsd	SNR
foreach my $id (@cmapstoprint) {
	my @positions = (keys $q_cmap{$id}); 
	@positions = grep ! /ContigLength/, @positions;
	@positions = grep ! /NumSites/, @positions;
	@positions =  sort {$a <=> $b} @positions;
	foreach my $pos (@positions) {
		print OUT "$id\t",sprintf("%.1f",$q_cmap{$id}{'ContigLength'}),"\t",$q_cmap{$id}{'NumSites'},"\t",$q_cmap{$id}{$pos}{'SiteID'},"\t",$q_cmap{$id}{$pos}{'LabelChannel'},"\t",sprintf("%.1f",$pos),"\t",$q_cmap{$id}{$pos}{'TheRest'},"\n"; 
	}
}
close OUT; 



if( $didmerge eq 1 ) {

	# << Align new _q.cmap with merged contigs back to _r.cmap >>
	print "\nAligning merged query maps to $inputs{rcmap}...\n"; 
	# Usage: doalignment($xmap, $rcmap, $qcmap, $out, $errbin, $mem, $cpuCount, $verbose)  
	doalignment($inputs{xmap}, $inputs{rcmap}, $tmpprefix.".cmap", $tmpprefix, $inputs{errbin}, $mem, $cpuCount, (exists $inputs{verbose}));
	if ( exists $inputs{verbose} ) { 
		print "\n--- Santity Check ---\n"; 
		my @tmpids = getCmapIds($tmpprefix."_q.cmap"); 
		print "After alignment, ", $tmpprefix."_q.cmap", " has ", scalar(@tmpids), " maps and they are: \n"; 
		print join( "\n", @tmpids ); 
		print "\n"; 
	}


	## << Call stitchFsites again >> 
	print "\n\n-----------------------------\n\n"; 
	print "\nTrying another round of stitchFsites...\n\n"; 
	# Usage: perl stitchFsites.pl --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output_folder> --bed <.bed fragile sites file> [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>] [--n <CPU cores to use>]
	my %args = %inputs; 
	$args{xmap} = $tmpprefix.".xmap";  
	$args{qcmap} = $tmpprefix."_q.cmap";
	$args{round} = $args{round}+1; 
	
	my @args;
    foreach my $opt (keys %args)  { push @args, "--$opt $args{$opt}"; }
	
	my $syscall = $^X." $0 ". join(" ",@args);
    system($syscall);
}
elsif( $didmerge eq 0 ) {
	# No more merges possible
	print "\n *** No more merges possible *** \n"; 

	# Write final results for this reference contig
	my $finalmap = $inputs{output}."/".$prefix."_fragileSiteRepaired_q.cmap"; 
	print "Writing final results for RefContig $RefContigID to $finalmap...\n";  
	if( exists $inputs{verbose} ) { print "\t...copying from $inputs{qcmap}\n"; }
	copy($inputs{qcmap} , $finalmap) or die "Copy failed: $!";
	# Clean up
	unlink glob $inputs{output}."/".$prefix."_tmp*";

	# Summary for this reference contig
	if ( exists $inputs{verbose} ) { 
		print "\n--- Santity Check ---\n"; 
		my @tmpids = getCmapIds($finalmap); 
		print "Final map $finalmap has ", scalar(@tmpids), " maps and they are: \n"; 
		print join( "\n", @tmpids ); 
		print "\n"; 
	}
	print "\n\n*** CMap Summary: $finalmap ***\n"; 
	calcCmapStats($finalmap); 
	print "\n\n*** Stitched Sites Summary: $outbed ***\n"; 
	getFsiteStats($outbed); 

	# Exit 
	print "\n\n"; 
	exit 0; 	
}




## *************** SUBROUTINES *************** ##

sub makeNewID {
	# example useage: $newCMapId = makeNewID($firstQryContigID, $idOffset); 
	my ($oldID, $offset) = @_;
	my $newID = $oldID + $offset; 
	while( exists $q_cmap{$newID} ) { $newID = $newID + $offset; }
	return($newID);
}

sub wanted { 
	m/_tmp/ and do { unlink $_ or warn "Could not unlink file $_\n"; }; 
} 

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}



