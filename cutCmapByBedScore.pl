#!/usr/bin/perl

# Usage: perl cutCmapByBedScore.pl <--bed stitchPositions_scored.bed> <--stats stitchPositions_scoreStats.csv> <--cmap fragileSiteRepaired_merged.cmap> <--output output folder> [<--prefix output prefix>] [--threshold <score threshold below which to cut default=1.0>] [--maxfill <minlen filter bp default=30000>] [--breakNGSonly <break maps if stitch has no BioNano support default=off>]
# Example:

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long qw(:config bundling); 
use Scalar::Util qw(looks_like_number);

print "\n";
print qx/ps -o args $$/;
print "\n";

## << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'stats=s', 'bed=s', 'cmap=s', 'prefix=s', 'threshold=s', 'maxfill:i', 'output:s', 'breakNGSonly'); 

if( !exists $inputs{bed} | !exists $inputs{stats} | !exists $inputs{cmap} | !exists $inputs{output} ) {
	print "Usage: perl cutCmapByBedScore.pl <--bed stitchPositions_scored.bed> <--stats stitchPositions_scoreStats.csv> <--cmap fragileSiteRepaired_merged.cmap> <--output output folder> [<--prefix output prefix>] [--threshold <score threshold below which to cut default=1.0>] [--maxfill <minlen filter bp default=30000>] [--breakNGSonly <break maps if stitch has no BioNano support>]\n"; 
	exit 0; 
}
foreach my $key ("bed","stats","cmap", "output") {
	if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
		$inputs{$key} = abs_path($inputs{$key});
	}
}


my $scoreThreshold = 1.0;
my $maxfill = 300000;
if (exists $inputs{threshold} && looks_like_number($inputs{threshold})) { $scoreThreshold = $inputs{threshold}; }
if (exists $inputs{maxfill}) { $maxfill = $inputs{maxfill} / 1000; }

print "\nBreaking maps at stitchPositions with score <$scoreThreshold...\n\n";


# Load SCORED input BED file (stitchPositions_scored.bed)
my $bedFileIn = $inputs{bed};
my @bedIn;
open BED, "<$bedFileIn" or die "ERROR: Could not open $bedFileIn: $!\n";
my $hasSeq=0;

while (my $line = <BED>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		#if ($s[0] == $refId) {
			my %bedLine = (
				"CMapId"  => "$s[0]",
				"Start" => "$s[1]", 
				"End"  => "$s[2]",
				"Type"  => "$s[3]",
				"Score" => "$s[4]",
				"Strand" => "$s[5]",
				"ThickStart" => "$s[6]",
				"ThickEnd" => "$s[7]",
				"ItemRgba" => "$s[8]",
				#"Sequence" => "$s[9]",
				"line" => $line
			);
			# deal with BED file that has sequences
			if (scalar(@s) == 10) {
				$hasSeq=1;
				$bedLine{"Sequence"} = "$s[9]";		
			}			
			push @bedIn, \%bedLine;
		#}
	}
}
@bedIn =  sort { $a->{'CMapId'} <=> $b->{'CMapId'} || $a->{'Start'} <=> $b->{'Start'} || $a->{'End'} <=> $b->{'End'} } @bedIn;
print scalar(@bedIn)." stitch positions identified in $bedFileIn\n";
close BED; 


# Load SCORED input STATS file (stitchPositions_scoreStats.csv)
my $statsFileIn = $inputs{stats};
my @statsIn;
open STATS, "<$statsFileIn" or die "ERROR: Could not open $statsFileIn: $!\n";

while (my $line = <STATS>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(",",$line);
		#BedLine,Type,FsiteSize,RefId,RefStart-End,RefLabelStart-End,MapId,MapStart-End,MapLabelStart-End,AvgMapCov,TotalHits,UniqueMolecules,TotalNGSHits,TotalConf,AvgConfPerHit,AvgConfPerMol,Score,AltScore,MoleculeMatches,NGSMatches
			my %statsLine = (
				"BedLine"  => "$s[0]",
				#"Type" => "$s[1]", 
				#"FsiteSize" => "$s[2]",
				#"RefId"  => "$s[3]",
				#"RefStart-End"  => "$s[4]",
				#"RefLabelStart-End" => "$s[5]",
				"MapId" => "$s[6]",
				"MapStart-End" => "$s[7]",
				#"MapLabelStart-End" => "$s[8]",
				#"AvgMapCov" => "$s[9],
				#"TotalHits" => "$s[10]",
				"UniqueMolecules" => "$s[11]",
				"TotalNGSHits" => "$s[12]",
				#"TotalConf" => "$s[13]",
				#"AvgConfPerHit" => "$s[14]",
				#"AvgConfPerMol" => "$s[15]",
				#"Score" => "$s[16]",
				#"AltScore" => "$s[17]",
				#"MoleculeMatches" => "$s[18]",
				#"NGSMatches" => "$s[19]"
			);
						
			push @statsIn, \%statsLine;
		
	}
}
@statsIn =  sort { $a->{'BedLine'} <=> $b->{'BedLine'} } @statsIn;
print scalar(@statsIn)." stitch positions identified in $statsFileIn\n";
close STATS; 


# Load final merged fragileSiteRepaired CMAP
my $cmapFileIn = $inputs{cmap};
my @cmapIn;
open CMAP, "<$cmapFileIn" or die "ERROR: Could not open $cmapFileIn: $!\n";

while (my $line = <CMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	
		#if ($s[0] == $refId) {
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
			push @cmapIn, \%cmap_line;	
		#}
	}
}
close CMAP; 

## Load alignref_final XMAP
#my $xmapFileIn = $ARGV[2];
#my @xmapIn;
#my $alignments;
#open XMAP, "$xmapFileIn" or die "ERROR: Could not open $xmapFileIn: $!\n";
#while (my $line = <XMAP>) {
	#chomp($line);
	##if header then skip
	#if ($line =~ /^#/) {
		#next; }
	#else {
		#$alignments++;
		#my @s = split(/\t/,$line);
		#for (my $i=0; $i<scalar(@s); $i++) {
			#$s[$i] =~ s/^\s+|\s+$//g; }
		##h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment
		##f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 	float 	float 	int         	string   
		#my %xmap_line = (
				#"XmapEntryID"  => "$s[0]", 
				#"QryContigID" => "$s[1]", 
				#"RefContigID"  => "$s[2]",	
				#"QryStartPos"  => "$s[3]", 
				#"QryEndPos"  => "$s[4]", 
				#"RefStartPos"  => "$s[5]", 
				#"RefEndPos" => "$s[6]", 
				#"Orientation" => "$s[7]", 
				#"Confidence" => "$s[8]", 
				#"HitEnum" => "$s[9]",
				#"QryLen" => "$s[10]",
				#"RefLen" => "$s[11]",
				#"LabelChannel" => "$s[12]",
				#"Alignment" => "$s[13]"
			#); 
		#push @xmapIn, \%xmap_line; }
#}
#print "Read in $alignments alignments from $xmapFileIn\n";
#print "\n";

#if (scalar(@statsIn) != scalar(@bedIn)) {
	#die "ERROR: Unequal number of entries in BED and stats file!\n";
#}

# open output BED file
my $prefix = "";
if (exists $inputs{prefix}) { $prefix = $inputs{prefix}; }
else { $prefix = basename("$bedFileIn",".bed"); }
my $bedPrefix = basename("$bedFileIn",".bed");
my $bedFileOut = $inputs{output}."/".$bedPrefix."_final.bed";
open BEDOUT, ">$bedFileOut" or die "ERROR: Could not open $bedFileOut: $!\n";
my $outputPrefix = $prefix;
# open output breakpoints file
my $brkptFile = $inputs{output}."/".$bedPrefix."_breakpoints.bed";
open BREAK, ">$brkptFile" or die "ERROR: Could not open $brkptFile: $!\n";

if ($hasSeq==1) {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n";
	print BREAK "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n";
}
else {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\n";
	print BREAK "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\n";
}


### DO WORK ###

my %breakPoints;
my $undoCount=0;
my $leftover=0;

for (my $i=0; $i<scalar(@bedIn); $i++) {
	my $bedLine_ref = $bedIn[$i];
	my %bedLine = %$bedLine_ref;
	#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
	my $refCmapId = $bedLine{'CMapId'};
	my $bedStart = $bedLine{'Start'};
	my $bedEnd = $bedLine{'End'};
	my $score = $bedLine{'Score'};
	
	my $statsLine_ref = $statsIn[$i];
	my %statsLine = %$statsLine_ref; 			#BedLine,Type,FsiteSize,RefId,RefStart-End,RefLabelStart-End,MapId,MapStart-End,MapLabelStart-End,AvgMapCov,TotalHits,UniqueMolecules,TotalNGSHits,TotalConf,AvgConfPerHit,AvgConfPerMol,Score,AltScore,MoleculeMatches,NGSMatches
	my $qryId = $statsLine{'MapId'};
	my $ngsOnly=0;
	my $uniqMols = $statsLine{'UniqueMolecules'};
	my $ngsHits = $statsLine{'TotalNGSHits'};
	if ((int($uniqMols)==0) && (int($ngsHits)>0)) { $ngsOnly = 1; }
	my $qryStartEnd = $statsLine{'MapStart-End'};
	my @s = split("-",$qryStartEnd);
	my $qryStart = $s[0]; $qryStart = $qryStart/1000;
	my $qryEnd = $s[1];	$qryEnd = $qryEnd/1000;
	
	##### SCORE THRESHOLD FILTERING AND/OR BREAK-NGS-ONLY#####
	if ( (defined($inputs{breakNGSonly}) && $ngsOnly) || $score < $scoreThreshold) { ## get IDs of maps where the stitch score is less than threshold
		$undoCount++;
				
		#my $breakpoint = 0;
		my @breakpoints;
		
		# get breakpoint locations
		if (($qryStart == $qryEnd)) {
			#$breakpoint = $qryStart;
			push @breakpoints, $qryStart;
		}
		elsif ($qryEnd == 0) {
			#$breakpoint = $qryStart;
			push @breakpoints, $qryStart;
		}
		elsif ($qryStart == 0) {
			#$breakpoint = $qryEnd;
			push @breakpoints, $qryEnd;
		}		
		else {
			#$breakpoint = "$qryStart $qryEnd";
			push @breakpoints, $qryStart;
			push @breakpoints, $qryEnd;
		}
		
		# push breakpoints to hash	
		my $line = "";	
		if ( !exists $breakPoints{$qryId} ) {
			foreach (sort {$a<=>$b} @breakpoints) {
				$line = $line." ".$_;
			}
			$breakPoints{$qryId} = $line;
		}
		else {
			my $temp = $breakPoints{$qryId};
			my @old = split(" ",$temp);
			push @breakpoints, @old;
			foreach (sort {$a<=>$b} @breakpoints) {
				$line = $line." ".$_;
			}
			$breakPoints{$qryId} = $line;
		}	
		
		#print breakpoints to BED
		print BREAK $bedLine{'line'}."\n";
		
	}
	else {
		$leftover++;
		print BEDOUT $bedLine{'line'}."\n";
	}
}

my $mapCount=0;
# loop over identified maps and locations
my $breakLine="";
for my $id (sort keys %breakPoints) {
	$breakPoints{$id} =~ s/^\s+|\s+$//g;
	my $out = "$id $breakPoints{$id}";
	#print "\t$out\n";
	if ($id != 0) {
		$breakLine = $breakLine. " ".$out;
		$mapCount++;
	}
}
$breakLine =~ s/^\s+|\s+$//g;
print "\n";

print "Undo stitch count: $undoCount affecting total maps: $mapCount\n";
print "$leftover stitches with >$scoreThreshold score output in $bedFileOut\n\n";

# setup RefAligner command to break the input CMAP at specificed locations and output
my $cmd;
my $RefAligner = glob("~/tools/RefAligner");
if ($undoCount > 0) {
	$cmd = "cd $inputs{output}; $RefAligner -f -i $cmapFileIn -o $inputs{output}/$outputPrefix -merge -minlen $maxfill -minsites 2 -break 100000 ".$breakLine." -stdout -stderr";
}
else {
	$cmd = "cd $inputs{output}; $RefAligner -f -i $cmapFileIn -o $inputs{output}/$outputPrefix -merge -minlen $maxfill -minsites 2 -stdout -stderr";
}
print "Running command: $cmd\n";
system($cmd);
print "\n";
$cmd = "cd $inputs{output}; $RefAligner -f -i $inputs{output}/$outputPrefix".".cmap -o $outputPrefix -merge -minlen $maxfill -minsites 2 -stdout -stderr";
print "Running command: $cmd\n";
system($cmd);
print "\n";

