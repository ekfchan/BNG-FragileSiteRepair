#!/usr/bin/perl

# A wrapper script to gap-fill: stitch together genome maps that are sufficiently close to each other and (optionally) overlapping fragile sites as predicted from reference genome map.

# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -e <input.errbin> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]

# Details: 
# * Assumption: that contigs on XMAP is being read from left to right and is sorted by RefStartPos
# * Assumption: input _r.cmap has only one contig
# * If fragile sites BED file is provided then only contigs spanning fragile sites will be merged. 
# * Designed to work with XMAP 0.1

use strict;
use warnings;
use IPC::System::Simple qw(system capture);
sub unique;
use Cwd;
use File::Find; 
use Data::Dumper;
use File::Basename;
use File::Copy;
use Getopt::Long qw(:config bundling); 
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Sys::MemInfo qw(totalmem freemem totalswap);

print "\n";
print qx/ps -o args $$/;
print "\n";

## << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'x|xmap=s', 'q|qcmap=s', 'r|rcmap=s', 'e|errbin=s', 'o|output|prefix=s', 'bed|b:s', 'round:i', 'maxlab:i', 'maxfill:i', 'wobble:i'); 

if( !exists $inputs{x} | !exists $inputs{q} | !exists $inputs{r} | !exists $inputs{e} | !exists $inputs{o} ) {
	print "usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -e <input.errbin> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round=1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 10000>] [--wobble <fragile site wobble in bp = 250>]\n"; 
	exit 0; 
}

my $scriptspath = Cwd::abs_path(dirname($0));	#gapFill.pl relies on two additional scripts (extractContigs.pl and mergeContigs.pl) and assumes thay are available from the same directory as gapFill.pl

## Default values
if( !exists $inputs{round} ) { $inputs{round} = 1; }
my $maxBp = 10000; 
if( exists $inputs{maxfill} ) { $maxBp = $inputs{maxfill}; }
my $wobble = 250; 
if( exists $inputs{wobble} ) { $wobble = $inputs{wobble}; }
my $maxlab = 0; 
if( exists $inputs{maxlab} ) { $maxlab = $inputs{maxlab}; }

open XMAP, $inputs{x} or die $!;
open QCMAP, $inputs{q} or die $!;
open RCMAP, $inputs{r} or die $!;

my $outName = $inputs{q}; 
$outName =~ s/_q.cmap/_merged_contigs_q.cmap/g;

my $outName2 = $inputs{q};
$outName2 =~ s/_q.cmap/_temp/g;

my $outName3 = $inputs{o};

my @xmap;
my @q_cmap;
my @r_cmap;

my $alignments=0;	#number of alignments
my @q_mapIds;
my @r_mapIds;	#RefContigID should always be the same 

my $q_sites=0;	#total number of sites across all query maps in _q.cmap
my $r_sites=0;	#total sites on single contig of _r.cmap

my $mergeCount=0;
my @firstContigList = ();
my @secondContigList= ();
#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
my @fsiteTypeList = ();
my @scoreList = ();
my @strandList = ();
my @thickStartList = ();
my @thickEndList = ();
my @itemRgbaList = ();
my @seqList = ();

my $idOffset = 100000;

#get number of CPU cores
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
#printf "There are %d CPUs\n"  , $cpu->count || 1;
my $cpuCount = $cpu->count;
#get system RAM
my $mem = (((&totalmem / 1024) / 1024)) / 1024; 

print "\n";

## << read input XMAP >> 
while (my $line = <XMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$alignments++;
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum
		#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 
		# 1	14839	17	344.2	999946.6	42640026.3	43628302.0	+	138.45	2M1D4M1D8M1D10M1D7M1D7M1D2M1D6M1D2M1D2M1I10M2D3M1D13M2D6M1D10M1D10M1D6M5I4M
		my %xmap_line = (
				"XmapEntryID"  => "$s[0]", 
				"QryContigID" => "$s[1]", 
				"RefContigID"  => "$s[2]",	#should be identical (one of the assumptions)
				"QryStartPos"  => "$s[3]", 
				"QryEndPos"  => "$s[4]", 
				"RefStartPos"  => "$s[5]", 
				"RefEndPos" => "$s[6]", 
				"Orientation" => "$s[7]", 
				"Confidence" => "$s[8]", 
				"HitEnum" => "$s[9]"
			); 
		push @xmap, \%xmap_line; }
}
print "Read in $alignments alignments from $inputs{x}\n";
print "\n";

## << read input query CMAP >> 
while (my $line = <QCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$q_sites++;
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	GmeanSNR	lnSNRsd	SNR
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	12.4650	0.4814	0
		push @q_mapIds, $s[0];
		my %cmap_line = (
			"CMapId"  => "$s[0]",
			"ContigLength" => "$s[1]",
			"NumSites"  => "$s[2]",
			"SiteID"  => "$s[3]",
			"LabelChannel"  => "$s[4]",
			"Position"  => "$s[5]",
			"StdDev" => "$s[6]",
			"Coverage" => "$s[7]",
			"Occurrence" => "$s[8]"
			);
		push @q_cmap, \%cmap_line;	
	}
}
print "Read in $q_sites sites from $inputs{q}\n";
@q_mapIds = unique(@q_mapIds);
print "Read in ".scalar(@q_mapIds)." maps from $inputs{q}\n";
print "\n";

## << read input reference CMAP >> 
while (my $line = <RCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$r_sites++;
		my @s = split(/\t/,$line);
		#load first contig into hash
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	
		push @r_mapIds, $s[0];
		my %cmap_line = (
			"CMapId"  => "$s[0]",	#should be identical (only single contig in _r.cmap)
			"ContigLength" => "$s[1]",
			"NumSites"  => "$s[2]",
			"SiteID"  => "$s[3]",
			"LabelChannel"  => "$s[4]",
			"Position"  => "$s[5]",
			"StdDev" => "$s[6]",
			"Coverage" => "$s[7]",
			"Occurrence" => "$s[8]"
			);
		push @r_cmap, \%cmap_line;	}
}
print "Read in $r_sites sites from $inputs{r}\n";
@r_mapIds = unique(@r_mapIds);
print "Read in ".scalar(@r_mapIds)." maps from $inputs{r}\n";
if (scalar(@r_mapIds) > 1) {
	die "Input _r.cmap should have only 1 map\n"; }
print "\n";

## << fsite bed >> 
# IF fsites BED file provided, only merge contigs that span fragile sites

my $hasfsites = 0; 
my @fsites;
my $fragileSites = 0;

if( exists $inputs{bed} ) {
	$hasfsites = 1; 
	my $fsitesFile = Cwd::abs_path($inputs{bed});

	if (-e $fsitesFile && defined $wobble) {
		print "Fragile sites BED file $fsitesFile provided\n";
		print "Filtering merge list based on fragile sites BED file: $fsitesFile\n";
		
		#read in .fsites file
		open FSITES, "$fsitesFile" or die $!;
		while (my $line = <FSITES>) {
			chomp($line);
			#if header then skip
			if ($line =~ /^#/) {
				next; }
			else {
				my @s = split(/\t/,$line);
				if ($s[0] eq $r_mapIds[0]) {
					$fragileSites++;
					#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
					my %fsites_line = (
						"CMapId"  => "$s[0]",
						"Start" => "$s[1]", 
						"End"  => "$s[2]",
						"Type"  => "$s[3]",
						"Score" => "$s[4]",
						"Strand" => "$s[5]",
						"ThickStart" => "$s[6]",
						"ThickEnd" => "$s[7]",
						"ItemRgba" => "$s[8]",
						"Sequence" => "$s[9]",
						"line" => $line
					);
					push @fsites, \%fsites_line;
				}
			}
		}
		my @fsites =  sort { $a->{'Start'} <=> $b->{'Start'} } @fsites;
		print scalar(@fsites)." fragile sites identified in $fsitesFile \n";
		close FSITES; 
	}
	else {
		print ".fsites and/or wobble undefined. Skipping .fsites filter\n\n"
	}
}

## << calculate distance in labels between each adjacent alignments >> 
# Generate list of contigs to be merged

print "\n";

my $previous = 0;

for (my $i=0; $i < scalar(@xmap); $i++) {
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
		print "\nConsidering merge between XmapEntryID: $id1 and XmapEntryID: $id2 Distance: ".abs($secondRefStartPos - $firstRefEndPos)."\n";
		if (($secondRefStartPos >= $firstRefEndPos) && ($secondRefEndPos >= $firstRefEndPos) && ($secondRefStartPos >= $firstRefStartPos)) {
		print "\tOverlap filter: PASS\n";
		if (abs($secondRefStartPos - $firstRefEndPos) <= $maxBp) {
			print "\tDistance filter: PASS\n";
			
			#check to make sure that the alignment extends to start/end of contig
			my $firstQryStartPosSite;
			my $firstQryEndPosSite;
			my $secondQryStartPosSite;
			my $secondQryEndPosSite;
			my $firstNumSites;
			my $secondNumSites;
			foreach my $hash (@q_cmap) {
				#print Dumper($hash);
				if (($hash->{'CMapId'} eq $firstQryContigID) && ($hash->{'Position'} eq $firstQryStartPos)) {
					$firstQryStartPosSite = $hash->{'SiteID'}; 
					$firstNumSites = $hash->{'NumSites'};
				}
				elsif (($hash->{'CMapId'} eq $firstQryContigID) && ($hash->{'Position'} eq $firstQryEndPos)) {
					$firstQryEndPosSite = $hash->{'SiteID'}; 
					$firstNumSites = $hash->{'NumSites'};
				}
				elsif (($hash->{'CMapId'} eq $secondQryContigID) && ($hash->{'Position'} eq $secondQryStartPos)) {
					$secondQryStartPosSite = $hash->{'SiteID'}; 
					$secondNumSites = $hash->{'NumSites'};
				}
				elsif (($hash->{'CMapId'} eq $secondQryContigID) && ($hash->{'Position'} eq $secondQryEndPos)) {
					$secondQryEndPosSite = $hash->{'SiteID'}; 
					$secondNumSites = $hash->{'NumSites'};
				}
			}
			if ( (($firstOrientation eq "+" && $secondOrientation eq "+") && ($firstQryEndPosSite eq $firstNumSites) && ($secondQryStartPosSite eq 1)) || (($firstOrientation eq "+" && $secondOrientation eq "-") && ($firstQryEndPosSite eq $firstNumSites) && ($secondQryStartPosSite eq $secondNumSites)) || (($firstOrientation eq "-" && $secondOrientation eq "+") && ($firstQryEndPosSite eq 1) && ($secondQryStartPosSite eq 1)) || (($firstOrientation eq "-" && $secondOrientation eq "-") && ($firstQryEndPosSite eq 1) && ($secondQryStartPosSite eq $secondNumSites)) ) {
			print "\tAlignment filter: PASS\n";
			my $firstRefEndPosSite = 0;
			my $secondRefStartPosSite = 0;
			foreach my $hash (@r_cmap) {
				if ($hash->{'Position'} eq $firstRefEndPos) {
					$firstRefEndPosSite = $hash->{'SiteID'}; 
				}
				if ($hash->{'Position'} eq $secondRefStartPos) {
					$secondRefStartPosSite = $hash->{'SiteID'}; 
				}
			}
			my $labelsDistance = $secondRefStartPosSite - $firstRefEndPosSite - 1;
			print "\tLabels: $labelsDistance\n";
			if ($labelsDistance <= $maxlab) {
				print "\tLabel filter: PASS\n";				
				if ( (grep {$_ eq $firstQryContigID} @firstContigList) || (grep {$_ eq $firstQryContigID} @secondContigList) || (grep {$_ eq $secondQryContigID} @secondContigList) || (grep {$_ eq $secondQryContigID} @firstContigList)) {
					print "\tOne of the contigs already ID'ed to merge...\n";
					next; }
				else {
					print "\tLooking for fragile sites...\n";				
					#IF fsites defined, only output contigs that span fragile sites
					if ($hasfsites!=0) {
						#previous and next reference label as borders for wobble
						my $prevLabelId = $firstRefEndPosSite - 1;
						my $nextLabelId = $secondRefStartPosSite + 1;
						my $prevLabelPos = 0;
						my $nextLabelPos = 0;
						foreach my $hash (@r_cmap) {
							if ($hash->{'SiteID'} eq $prevLabelId) {
								$prevLabelPos = $hash->{'Position'}; 
							}
							if ($hash->{'SiteID'} eq $nextLabelId) {
								$nextLabelPos = $hash->{'Position'}; 
							}
						}
						
						#print "\tPrevLabel Id: $prevLabelId Pos: $prevLabelPos\tNextLabel Id: $nextLabelId Pos: $nextLabelPos\n";
						my $globalMin=0; 
						my $globalMax=0;
						if ($labelsDistance != -1) {
							$globalMin = $firstRefEndPos;
							$globalMax = $secondRefStartPos;
						}
						else {
							$globalMin = $firstRefEndPos-1000;
							$globalMax = $secondRefStartPos+1000;
						}
						print "\tLooking for fsites in global window Min: $globalMin Max: $globalMax Range: ".abs($globalMax - $globalMin)."\n";
						#print "Working on $inputs{bed} with ".$fragileSites." fsites\n";
						
						my $fsiteFound = 0;
						foreach my $hash (@fsites) {
							#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
							my $start = $hash->{'Start'};
							my $end = $hash->{'End'};
							my $type = $hash->{'Type'};
							my $score = $hash->{'Score'};
							my $strand = $hash->{'Strand'};
							my $thickStart = $hash->{'ThickStart'};
							my $thickEnd = $hash->{'ThickEnd'};
							my $itemRgba = $hash->{'ItemRgba'};
							my $seq = $hash->{'Sequence'};
							my $min=0;
							my $max=0;
							if ( abs($start - $globalMin) > $wobble) {
								#print "\tDifference: ".abs($start - $globalMin)." Changing min from $globalMin to ".($start - $wobble)."\n";
								$min = $start - $wobble;	 
							}
							else {
								$min = $globalMin; 
							}
							if ( abs($globalMax - $end) > $wobble) {
								#print "\tDifference: ".abs($globalMax - $end)." Changing max from $globalMax to ".($end + $wobble)."\n";
								$max = $end + $wobble;
							}
							else {
								$max = $globalMax;
							}
							
							
							if ( ($firstRefEndPos >= $min) && ($firstRefEndPos <= $max) && ($secondRefStartPos <= $max) && ($secondRefStartPos >= $min) ) {		
								print "\tActual Fsite Start: $start End: $end Range: ".abs($end - $start)."\n";
								print "\tFragile site located within Window corrected for max wobble fsite window Min: $min Max: $max Range: ".abs($max - $min)."\n";
								
								#print "\tLeft wobble: $leftWobble Right wobble: $rightWobble\n";			
								next if ($firstQryContigID eq $previous);
								push @firstContigList,$firstQryContigID;
								push @secondContigList,$secondQryContigID;
								#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
								push @fsiteTypeList, $type;
								push @scoreList, $score;
								push @strandList, $strand;
								push @thickStartList, $thickStart;
								push @thickEndList, $thickEnd;
								push @itemRgbaList, $itemRgba;
								push @seqList, $seq;
								$previous = $firstQryContigID;
								$fsiteFound = 1;
							}
						}
						if ($fsiteFound eq 0) { print "\tNo fragile site found in Window\n"; }
					}
					else {
						#print $inputs{bed}, " undefined'."\n";
						push @firstContigList,$firstQryContigID;
						push @secondContigList,$secondQryContigID;
					}
				}				
			}
			else {
				print "\tLabel filter: FAIL\n";
			}
			}
			else {
				print "\tAlignment filter: FAIL\n";
				#if ( (($firstOrientation eq "+" && $secondOrientation eq "+") && ($firstQryEndPosSite eq $firstNumSites) && ($secondQryStartPosSite eq 1)) || (($firstOrientation eq "+" && $secondOrientation eq "-") && ($firstQryEndPosSite eq $firstNumSites) && ($secondQryStartPosSite eq $secondNumSites)) || (($firstOrientation eq "-" && $secondOrientation eq "+") && ($firstQryEndPosSite eq 1) && ($secondQryStartPosSite eq 1)) || (($firstOrientation eq "-" && $secondOrientation eq "-") && ($firstQryEndPosSite eq 1) && ($secondQryStartPosSite eq $secondNumSites)) ) {
				print "\tFirst contig: $id1 Orientation: $firstOrientation Start: $firstQryStartPos $firstQryStartPosSite End: $firstQryEndPos $firstQryEndPosSite NumSites: $firstNumSites\n";
				print "\tSecond contig: $id2 Orientation: $secondOrientation Start: $secondQryStartPos $secondQryStartPosSite End: $secondQryEndPos $secondQryEndPosSite NumSites: $secondNumSites\n";
			
			
			
			}
		}
		else {
			print "\tDistance filter: FAIL\n";
		}
		}
		else {
			print "\tOverlap filter: FAIL\n";
		}
	}
	else {
		print "All possible merges identified.\n\n";
	}	
}

# create AoH of pairs of contigs 
my @contigs;
for (my $i=0; $i<(scalar(@secondContigList)); $i++) {
	my %merge = (
		"firstContig" => $firstContigList[$i],
		"secondContig" => $secondContigList[$i]
		);
	push @contigs, \%merge;
}



## << Merge Contigs >> 
# Iterate over identified array list and merge contigs

my $nomerge = 0; 
if (scalar(@secondContigList) > 0) {
	if (scalar(@firstContigList) == scalar(@secondContigList)) {
		print "Gap filling between ".scalar(@secondContigList)." contigs\n\n";
		
		my $extractScript = $scriptspath."/extractContigs.pl";
		
		# Perform first merge
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		my @ARGS = ($inputs{x}, $inputs{q}, $inputs{r}, $firstContigList[0], $secondContigList[0],$fsiteTypeList[0],$scoreList[0],$strandList[0],$thickStartList[0],$thickEndList[0],$itemRgbaList[0],$seqList[0]);
		my $cwd = cwd();
		print "Running command: ".$^X." $extractScript ". join(" ",@ARGS)."\n";
		system($^X, "$extractScript", @ARGS);
		print "QryContigID $firstContigList[0] merged with QryContigID $secondContigList[0] into QueryContigID ".($firstContigList[0]+$idOffset)."\n";
		print "\n";
		
		my $temp = $secondContigList[0];
		
		# Perform next merges
		for (my $i=1; $i<(scalar(@secondContigList)); $i++) {
			next if ( ($secondContigList[$i] eq $secondContigList[0]) ) ;
			#for (my $i=1; $i<(3); $i++) {
			@ARGS = ($inputs{x}, $outName, $inputs{r}, $firstContigList[$i], $secondContigList[$i], $fsiteTypeList[$i],$scoreList[$i],$strandList[$i],$thickStartList[$i],$thickEndList[$i],$itemRgbaList[$i],$seqList[$i]);
			print "Running command: ".$^X." $extractScript ". join(" ",@ARGS)."\n";
			system($^X, "$extractScript", @ARGS);
			print "QryContigID $firstContigList[$i] merged with QryContigID $secondContigList[$i] into QueryContigID ".($firstContigList[$i]+$idOffset)."\n";
			print "\n";
		}	
	}
	else { print "ERROR: uneven number of contigs in firstQueryList and secondQueryList\n"; }
		
	print "DONE: merging ".scalar(@secondContigList)." contigs complete.\n\n";
		
	# Align new _q.cmap with merged contigs back to _r.cmap

	my $veto = q/-output-veto-filter '(_intervals.txt|.err|.maprate|[a-z|A-Z].map)$'/;
	#$veto = $veto." -output-veto-filter .err";
	#print "Veto: $veto\n";
	if ($inputs{round} == 1) {
		# Perform first alignment round
		print "=====  Performing round $inputs{round} alignment =====\n"; 
		# system("~/tools/RefAligner -ref $inputs{r} -i $outName -o $outName2 -maxthreads $cpuCount -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.10 -extend 1 -outlier 0.0001 -endoutlier 0.001 -deltaX 12 -deltaY 12 -xmapchim 14 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 2.9 -insertThreads 4 -nosplit 2 -biaswt 0 -f -maxmem $mem -T 1e-12 -BestRef 1 $veto -readparameters $inputs{e} -stdout -stderr");
		system("~/tools/RefAligner -ref $inputs{r} -i $outName -o $outName2 -maxthreads $cpuCount -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.10 -extend 1 -outlier 0.0001 -endoutlier 0.001 -deltaX 12 -deltaY 12 -xmapchim 14 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 2.9 -insertThreads 4 -nosplit 2 -biaswt 0 -f -maxmem $mem -T 1e-12 -BestRef 1 $veto -readparameters $inputs{e} -stdout");	
		print "\nFIRST ROUND COMPLETE.\n\n";
		
		# Launch second round
		# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round=1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		print "******   Launching round 2 ******\n"; 
		$outName2 =~ s/.xmap//g;
		my %args = %inputs; 
		$args{x} = $outName2.".xmap"; 
		$args{q} = $outName2."_q.cmap"; 
		$args{round} = 2; 	
		my @args;
        	foreach my $opt (keys %args)  { push @args, "--$opt $args{$opt}"; }
		
		my $syscall = $^X." $0 ". join(" ",@args);
        	print "Running command: $syscall \n";
        	system($syscall);
	}
	elsif ($inputs{round} == 2) {
		# If second round, perform 2nd round of merge
		print "======   Performing round $inputs{round} alignment ======= \n"; 
		# system("~/tools/RefAligner -ref $inputs{r} -i $outName -o $outName3 -maxthreads $cpuCount -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.10 -extend 1 -outlier 0.0001 -endoutlier 0.001 -deltaX 12 -deltaY 12 -xmapchim 14 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 2.9 -insertThreads 4 -nosplit 2 -biaswt 0 -f -maxmem $mem -T 1e-12 -BestRef 1 $veto -readparameters $inputs{e} -stdout -stderr");
		system("~/tools/RefAligner -ref $inputs{r} -i $outName -o $outName3 -maxthreads $cpuCount -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.10 -extend 1 -outlier 0.0001 -endoutlier 0.001 -deltaX 12 -deltaY 12 -xmapchim 14 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -mres 2.9 -insertThreads 4 -nosplit 2 -biaswt 0 -f -maxmem $mem -T 1e-12 -BestRef 1 $veto -readparameters $inputs{e} -stdout");	
		print "\nSECOND ROUND COMPLETE.\n\n";
		#$inputs{round} = 0;
		exit 0; #If second round, exit script gracefully to return back to first round script to do more work
	}
#	else { exit 0; }

	# If second round, exit script gracefully to return back to first round script to do more work
	#if ($inputs{round} == 2) { exit 0; }

	# Clean up temp files
	my $cwd = cwd();
	#find(\&wanted, $cwd); 
	#find(\&wanted2, $cwd); 
	
} else {
	print "No merges possible. Copying input data unchanged.\n";
	copy($inputs{x},$outName3.".xmap") or die "Copy failed: $!";
	copy($inputs{q},$outName3."_q.cmap") or die "Copy failed: $!";
	copy($inputs{r},$outName3."_r.cmap") or die "Copy failed: $!";
	$nomerge = 1; 
}


## << Calculate stats on original and new CMAPs >> 
# usage: calc_cmap_stats.pl <CMAP_File>
my $dir = glob("~/scripts/HybridScaffold/scripts");
my $script = "calc_cmap_stats.pl";
if (-e "$dir/$script") {
	print "Original query CMAP $inputs{q} stats:\n";
	#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
	my $cmd = "perl $dir/$script $inputs{q}";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	
	if( !$nomerge ) {
		my $file = $outName3."_q.cmap";
		print "Gap-filled query CMAP $file stats:\n";
		#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
		$cmd = "perl $dir/$script $file";
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";
	}
	
}
else {
	print "perl script calc_cmap_stats.pl not found at $script\n"; }





sub wanted { 
	m/_temp/ and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; } 

sub wanted2 { 
	m/_merged_contigs_q/ and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; }



sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
