#!/usr/bin/perl

#Usage: perl scoreStitchPositionsBED_v4.pl <--bed stitchPositions.bed> <--xmap alignref_final.xmap> <--qcmap alignref_final_q.cmap> <--rcmap alignref_final_r.cmap> <--alignmolxmap alignmol.xmap> <--alignmolrcmap alignmol_r.cmap> <--fasta reference.fasta> <--key key file from fa2cmap> <--bam NGSalignments.bam to ref fasta> <--output output directory> <--ngsBuffer bp +/- stitch location to require NGS alignment support> [<--n cpu cores>] [--ngsBonus <raw score value give supporting NGS alignemnts =10>]
#Example: 

use strict;
use warnings;
use Cwd qw(abs_path);
use Data::Dumper;
use File::Basename;
use DateTime;
use DateTime::Format::Human::Duration;
use threads;
use threads::shared;
use Getopt::Long qw(:config bundling); 
use JSON::XS;
use Bio::DB::Sam;

my $sharedmem = 0; ### $sharedmem = 1 means less mem usage but slower
$| = 1; # set perl to flush buffer immediately after write out command

my $startDateTime = DateTime->now;
my $endDateTime = DateTime->now;
my $formatDateTime = DateTime::Format::Human::Duration->new();

print "\n";
print qx/ps -o args $$/;
print "\n";

## << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'xmap=s', 'qcmap=s', 'rcmap=s', 'cmap=s', 'alignmolxmap=s', 'alignmolrcmap=s','bed:s', 'n:i','output:s', 'bam:s', 'fasta:s', 'key:s', 'ngsBonus:i', 'ngsBuffer:i'); 

if( !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{alignmolxmap} | !exists $inputs{output} | !exists $inputs{alignmolrcmap} ) {
	print "Usage: perl scoreStitchPositionsBED_v4.pl <--bed stitchPositions.bed> <--xmap alignref_final.xmap> <--qcmap alignref_final_q.cmap> <--rcmap alignref_final_r.cmap> <--alignmolxmap alignmol.xmap> <--alignmolrcmap alignmol_r.cmap> <--fasta reference.fasta> <--key key file from fa2cmap> <--bam NGSalignments.bam to ref fasta> <--output output directory> <--ngsBuffer bp +/- stitch location to require NGS alignment support> [<--n cpu cores>] [--ngsBonus <raw score value give supporting NGS alignemnts =10>]\n"; 
	exit 0; 
}

my $scoreBuffer = 500;
if (exists $inputs{ngsBuffer}) {
	$scoreBuffer = $inputs{ngsBuffer};
}

my $scriptspath = Cwd::abs_path(dirname($0));

my $hostname = `hostname`;
chomp($hostname);
print "Running on host: $hostname\n\n";

print "\nStart time : "; print join ' ', $startDateTime->ymd, $startDateTime->hms; print "\n\n";

#get number of CPU cores
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
#printf "There are %d CPUs\n"  , $cpu->count || 1;
my $cpuCount = $cpu->count;
if (exists $inputs{n}) { $cpuCount = int($inputs{n}); }
my $ngsBonus = 10;
if (exists $inputs{ngsBonus}) { $ngsBonus = int($inputs{ngsBonus}); }

foreach my $key ("xmap","qcmap","rcmap","cmap","alignmolxmap", "alignmolrcmap", "bed", "output", "bam", "fasta", "key") {
	if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
		$inputs{$key} = abs_path($inputs{$key});
	}
}

# Load input BED file with stitchPositions 
my $bedFileIn = $inputs{bed};
my $prefix = basename("$bedFileIn",".bed");
my @bedIn;
share(@bedIn) if $sharedmem;
open BED, "<$bedFileIn" or die "ERROR: Could not open $bedFileIn: $!\n";
my $hasSeq=0;
my @bed_format = qw{CMapId   Start   End   Type   Score   Strand   ThickStart   ThickEnd   ItemRgba   Sequence};

while (my $line = <BED>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		push @bedIn, get_json_by_format($line, \@bed_format);
	}
}
@bedIn =  sort { decode_json($a)->{'CMapId'} <=> decode_json($b)->{'CMapId'} || decode_json($a)->{'Start'} <=> decode_json($b)->{'Start'} || decode_json($a)->{'End'} <=> decode_json($b)->{'End'} } @bedIn;
print scalar(@bedIn)." stitch positions identified in $bedFileIn\n";
close BED; 


# Load alignref XMAP
my $xmapFileIn = $inputs{xmap};
my @xmap;
share(@xmap) if $sharedmem;
my $alignments;
my @xmap_format = qw{XmapEntryID  QryContigID  RefContigID  QryStartPos  QryEndPos  RefStartPos  RefEndPos  Orientation  Confidence  HitEnum  QryLen  RefLen  LabelChannel  Alignment};
open XMAP, "$xmapFileIn" or die "ERROR: Could not open $xmapFileIn: $!\n";
while (my $line = <XMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$alignments++;
		my ($arr_json, $ref_config_id) = get_arr_json_by_format($line, \@xmap_format);
		if ($sharedmem && !defined($xmap[$ref_config_id])) {
			my @a :shared = ();
			$xmap[$ref_config_id] = \@a;
		}
		push @{$xmap[$ref_config_id]}, $arr_json; 
	}
}
print "Read in $alignments alignments from $xmapFileIn\n";
close XMAP;
		

# Load alignref_final _q.cmap
my $qcmapFileIn = $inputs{qcmap};
my @qcmap;
share(@qcmap) if $sharedmem;
my @cmap_format = qw{CMapId      ContigLength    NumSites        SiteID  LabelChannel    Position        StdDev  Coverage        Occurrence};
open QCMAP, "$qcmapFileIn" or die "ERROR: Could not open $qcmapFileIn: $!\n";
while (my $line = <QCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my ($arr_json, $cmap_id) = get_arr_json_by_format($line, \@cmap_format, 0);
		if ($sharedmem && !defined($qcmap[$cmap_id])) {
			my @a :shared = ();
			$qcmap[$cmap_id] = \@a;
		}
		push @{$qcmap[$cmap_id]}, $arr_json;
	}
}
close QCMAP;
		

# Load alignref_final _r.cmap
my $rcmapFileIn = $inputs{rcmap};
my @rcmap;
share(@rcmap) if $sharedmem;
open RCMAP, "$rcmapFileIn" or die "ERROR: Could not open $rcmapFileIn: $!\n";
while (my $line = <RCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my ($arr_json, $cmap_id) = get_arr_json_by_format($line, \@cmap_format, 0);
		if ($sharedmem && !defined($rcmap[$cmap_id])) {
			my @a :shared = ();
			$rcmap[$cmap_id] = \@a;
		}
		push @{$rcmap[$cmap_id]}, $arr_json;
	}
}		
close RCMAP; 

# Load merged alignmol XMAP
my $molxmapFileIn = $inputs{alignmolxmap};
my @molxmap;
share(@molxmap) if $sharedmem;
my $molalignments;
open MOLXMAP, "$molxmapFileIn" or die "ERROR: Could not open $molxmapFileIn: $!\n";
while (my $line = <MOLXMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$molalignments++;
		my ($arr_json, $ref_config_id) = get_arr_json_by_format($line, \@xmap_format);
                if ($sharedmem && !defined($molxmap[$ref_config_id])) {
                        my @a :shared = ();
                        $molxmap[$ref_config_id] = \@a;
                }
		push @{$molxmap[$ref_config_id]}, $arr_json; 
	}
}
print "Read in $molalignments molecule alignments from $molxmapFileIn\n";

# Load alignmol_r.cmap
my $molrcmapFileIn = $inputs{alignmolrcmap};
my @molrcmap;
share(@molrcmap) if $sharedmem;
open MOLRCMAP, "$molrcmapFileIn" or die "ERROR: Could not open $molrcmapFileIn: $!\n";
while (my $line = <MOLRCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my ($arr_json, $cmap_id) = get_arr_json_by_format($line, \@cmap_format, 0);
		if ($sharedmem && !defined($molrcmap[$cmap_id])) {
			my @a :shared = ();
			$molrcmap[$cmap_id] = \@a;
		}
		push @{$molrcmap[$cmap_id]}, $arr_json;
	}
}		
close MOLRCMAP; 


# if input BAM and fasta and key files specified, load
my $sam; my $bai;
my @keys; my $fai;
my $ngsScoring=0;
my $singleThread=0;
my @key_format = qw{CompntId   CompntName   CompntLength};
if ( defined $inputs{fasta} && exists $inputs{fasta} && defined $inputs{bam} && exists $inputs{bam} && defined $inputs{key} && exists $inputs{key}) {
	# force index of fasta file
	$bai = $inputs{bam} . ".bai"; share($bai) if $sharedmem;
	$fai = $inputs{fasta} . ".fai"; share($fai) if $sharedmem;
	$singleThread = 1; $cpuCount = 1;
	if (-e $bai && -e $fai) {
		$sam = Bio::DB::Sam->new( -fasta=>"$inputs{fasta}", -fai=>"$fai", -bam=>"$inputs{bam}", -bai=>"$bai" );
	}
	else {
		my $status_code = Bio::DB::Bam->index_build($inputs{bam});
		my $index = Bio::DB::Sam::Fai->load($inputs{fasta});
		$sam = Bio::DB::Sam->new( -fasta=>"$inputs{fasta}", -bam=>"$inputs{bam}" );
	}
	#my $bai = $sam->bam_index;
	#share($sam) if $sharedmem;
	$ngsScoring=1;
	print "\n";
	print "Input reference FASTA: $inputs{fasta}\n";
	print "Input reference FASTA key: $inputs{key}\n";
	print "Input NGS alignment BAM: $inputs{bam}\n";
	print "\tReducing maximum cores to 1 ...\n";
	print "\tRequire $scoreBuffer bp +/- stitchPosition start/end to count NGS alignment as support\n";
	#print "\n";
	
	#load key file
	open KEY, "$inputs{key}" or die "ERROR: Could not open $inputs{key}: $!\n";
	while (my $line = <KEY>) {
		chomp($line);
		#if header then skip
		if ($line =~ /^#/) {
			next; }
		elsif ($line =~ /^CompntId/) {
			next; }
		else {
			my @s = split("\t",$line);
			#CompntId	CompntName	CompntLength
			my %keyLine = (
				"CompntId"  => "$s[0]", 
				"CompntName" => "$s[1]", 
				"CompntLength"  => "$s[2]"
			);
			push @keys, \%keyLine; 
		}	
	}
	share(@keys) if $sharedmem;	
	print "Read in ".scalar(@keys)." keys from $inputs{key}\n";
	print "\n";		
}


$endDateTime = DateTime->now;
$formatDateTime = DateTime::Format::Human::Duration->new();
print 'Spent time reading input files : ', $formatDateTime->format_duration_between($endDateTime, $startDateTime); print "\n\n";

print "Maximum CPU cores to use: $cpuCount\n\n";

# set various locks for shared data
my $output_lock :shared; #output lock
my $genomeLabelsLock :shared;

# configure shared data variables
my %genomeMapIds :shared;
my %genomeMapStartLabelIds :shared;
my %genomeMapEndLabelIds :shared;
my %genomeMapStartPoss :shared;
my %genomeMapEndPoss :shared;
my $bedEntryCount :shared = 0;
my $bedEntries :shared = scalar(@bedIn);

my %bedStartLabelIds :shared;
my %bedEndLabelIds :shared;

my %scoredBeds :shared;
my %scoreStats :shared;

### DO WORK ###

if (!$singleThread) {
	for (my $thread=0; $thread<$cpuCount; $thread++) {
		my $t = threads->create( sub {
			while(1) {
				my $curBed;
				{
					lock($bedEntryCount);
					$curBed = $bedEntryCount++;
				}
				last if ($curBed >= $bedEntries);
				#$sam->clone;
				process_bed($curBed);
			}
			threads->detach();
		});
	}

	my $thread_wait_cnt = 0;
	while (threads->list(threads::all)) {
		sleep(2);
	}
}
else {
	for (my $i=0; $i<$bedEntries; $i++) {
		process_bed($i);
	}
}

### OUTPUT RESULTS ###

# Open output stats file
my $statsFileOut = $inputs{output}."/".$prefix."_scoreStats.csv";
open STATS, ">$statsFileOut" or die "ERROR: Could not open $statsFileOut: $!\n";
close MOLXMAP;

# open scored BED output file 
my $bedFileOut = $inputs{output}."/".$prefix."_scored.bed";
open BEDOUT, ">$bedFileOut" or die "ERROR: Could not open $bedFileOut: $!\n";

if ($hasSeq==1) {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n";
}
else {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\n";
}

print STATS "#BedLine,Type,FsiteSize,RefId,RefStart-End,RefLabelStart-End,MapId,MapStart-End,MapLabelStart-End,AvgMapCov,TotalHits,UniqueMolecules,TotalNGSHits,TotalConf,AvgConfPerHit,AvgConfPerMol,Score,AltScore,MoleculeMatches,NGSMatches\n";

my $scoreCount=0;
my $statCount=0;

foreach my $i (sort {$a <=> $b} keys %scoredBeds) {
  print BEDOUT $scoredBeds{$i};
  $scoreCount++;
}

foreach my $i (sort {$a <=> $b} keys %scoreStats) {
  print STATS $scoreStats{$i};
  $statCount++;
}
close BEDOUT;
close STATS;

print "\n";
print "Output STATS file: $statsFileOut with $statCount entries\n";
print "Output BED file: $bedFileOut with $scoreCount entries\n";
print "\n";

$endDateTime = DateTime->now;
$formatDateTime = DateTime::Format::Human::Duration->new();
print "Total time spent scoring BED file: $inputs{bed} : ", $formatDateTime->format_duration_between($endDateTime, $startDateTime); print "\n\n";



### SUBFUNCTIONS ###

sub process_bed {
	my ($bedEntryCount) = @_;
	my %bedEntry = %{decode_json($bedIn[$bedEntryCount])};
	#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
	my $bedCmapId = $bedEntry{'CMapId'};
	my $bedStart = int($bedEntry{'Start'});
	my $bedEnd = int($bedEntry{'End'});
	my $bedStartLabelId=0;
	my $bedEndLabelId=0;
	my $curPos=0;

	# Get alignref _r.cmap label IDs that correspond to nearest labels to the left of the start and right of the end of the BED stitch locations 
	foreach my $rcmapEntry_ref (defined($rcmap[$bedCmapId])?@{$rcmap[$bedCmapId]}:()) {
		my %rcmapEntry = %{get_hash_from_arr_json_by_format($rcmapEntry_ref, \@cmap_format)};
		if ($rcmapEntry{'CMapId'} eq $bedCmapId) {
			$curPos = $rcmapEntry{'Position'};
			if ($curPos<=($bedStart+10)) {
				#print "\t$curPos\n";
				$bedStartLabelId = $rcmapEntry{'SiteID'};
			}
			if ($curPos>=($bedEnd-10)) {
				#print "\t$curPos\n";
				$bedEndLabelId = $rcmapEntry{'SiteID'};
				last;
			}
		}
	}
	
	### DEBUG ###
	#if ($bedEntryCount == 14) {
		#lock($output_lock);
		#print "BedLine: $bedEntryCount RefLabel Start-End: $bedStartLabelId-$bedEndLabelId\n";
	#}
	
	# Load bedStart and bedEnd label IDs into hash
	{
		lock($genomeLabelsLock);
		$bedStartLabelIds{$bedEntryCount} = $bedStartLabelId;
		$bedEndLabelIds{$bedEntryCount} = $bedEndLabelId;	
	}
	
	# get the cooresponding genome map label IDs for each start/end position
	my $genomeMapId=0;
	my %matchPairs;
	my $genomeMapStartPos=0;
	my $genomeMapEndPos=0;
	my $genomeMapStart = 0;
	my $genomeMapEnd = 0;
	my $genomeMapAvgCov=0;
	
	foreach my $xmapEntry_ref (defined($xmap[$bedCmapId])?@{$xmap[$bedCmapId]}:()) {
		my %xmapEntry = %{get_hash_from_arr_json_by_format($xmapEntry_ref, \@xmap_format)};		
		if ($xmapEntry{'RefContigID'} eq $bedCmapId) {
			my $line = $xmapEntry{'Alignment'};
			my @s = split(/(?<=\))/, $line);
			my $refLabels = "";
			my $qryLabels = "";			
			# process alignment doublet string
			foreach my $pair (@s) {
				#print "Doublet: $pair\n";
				my @match = $pair =~ /(\d+)/g;
				#print "Ref: $match[0] Qry: $match[1]\n";
				$refLabels = $refLabels."$match[0] ";
				$qryLabels = $qryLabels."$match[1] ";
				$matchPairs{$match[0]} = $match[1];
			}
			
			# check to see if genome map alignment contains both reference labels cooresponding to BED start and end
			if ( ($refLabels =~ /\b$bedStartLabelId\b/) && ($refLabels =~ /\b$bedEndLabelId\b/) ) {
				
				my @qryLabelsArray = split(" ",$qryLabels);
				$genomeMapStart = $qryLabelsArray[0];
				$genomeMapEnd = $qryLabelsArray[-1];
				
				$genomeMapId = $xmapEntry{'QryContigID'};
				last;
			}
		}
	}
	
	# find corresponding genome map if no alignemnts containing exact flanking labels of genome map
	if ($genomeMapId == 0) {
		foreach my $xmapEntry_ref (defined($xmap[$bedCmapId])?@{$xmap[$bedCmapId]}:()) {
			my %xmapEntry = %{get_hash_from_arr_json_by_format($xmapEntry_ref, \@xmap_format)};		
			if ($xmapEntry{'RefContigID'} eq $bedCmapId) {
				my $line = $xmapEntry{'Alignment'};
				my @s = split(/(?<=\))/, $line);
				my $refLabels = "";
				my $qryLabels = "";			
				# process alignment doublet string
				foreach my $pair (@s) {
					#print "Doublet: $pair\n";
					my @match = $pair =~ /(\d+)/g;
					#print "Ref: $match[0] Qry: $match[1]\n";
					$refLabels = $refLabels."$match[0] ";
					$qryLabels = $qryLabels."$match[1] ";
					$matchPairs{$match[0]} = $match[1];
				}
				
				# check to see if genome map alignment contains both reference labels cooresponding to BED start and end
				my $start = $bedStartLabelId - 1; my $start2 = $bedStartLabelId - 2;
				my $end = $bedEndLabelId + 1; my $end2 = $bedEndLabelId + 2;
				#if ( (($refLabels =~ /\b$bedStartLabelId\b/) || ($refLabels =~ /\b$start\b/) || ($refLabels =~ /\b$start2\b/)) && (($refLabels =~ /\b$bedEndLabelId\b/) || ($refLabels =~ /\b$end\b/) || ($refLabels =~ /\b$end2\b/)) ) {
				if ( ($refLabels =~ /\b$bedStartLabelId\b/) || ($refLabels =~ /\b$bedEndLabelId\b/) ) {
					
					my @qryLabelsArray = split(" ",$qryLabels);
					$genomeMapStart = $qryLabelsArray[0];
					$genomeMapEnd = $qryLabelsArray[-1];
					
					$genomeMapId = $xmapEntry{'QryContigID'};
					last;
				}
			}
		}
	}	
				
	# loop thru _q.cmap to find coordinates of identified genome map labels
	foreach my $qcmapEntry_ref (defined($qcmap[$genomeMapId])?@{$qcmap[$genomeMapId]}:()) {
		my %qcmapEntry = %{get_hash_from_arr_json_by_format($qcmapEntry_ref, \@cmap_format)};
		if ($qcmapEntry{'CMapId'} eq $genomeMapId) {
			#if (exists $matchPairs{$bedStartLabelId} && exists $matchPairs{$bedEndLabelId}) {
				if (exists $matchPairs{$bedStartLabelId}) {
					if ($qcmapEntry{'SiteID'} eq $matchPairs{$bedStartLabelId}) {
						#$genomeMapStartPoss{$bedEntryCount} = $qcmapEntry{'Position'};
						$genomeMapStartPos = $qcmapEntry{'Position'};
					}
				}
				else {
					$genomeMapStartPos = 0;
				}
				
				if (exists $matchPairs{$bedEndLabelId}) {
					if ($qcmapEntry{'SiteID'} eq $matchPairs{$bedEndLabelId}) {
						#$genomeMapEndPoss{$bedEntryCount} = $qcmapEntry{'Position'};
						$genomeMapEndPos = $qcmapEntry{'Position'};
					}
				}
				else {
					$genomeMapEndPos = 0;
				}
			#}
		}
	}
	
	# loop thru molrcmap to get average coverage of genome map by single molecules
	my $count = 0;
	foreach my $molrcmapEntry_ref (defined($molrcmap[$genomeMapId])?@{$molrcmap[$genomeMapId]}:()) {
		my %molrcmapEntry = %{get_hash_from_arr_json_by_format($molrcmapEntry_ref, \@cmap_format)};
		if (($molrcmapEntry{'CMapId'} eq $genomeMapId) && (defined $molrcmapEntry{'SiteID'}) && (defined $matchPairs{$bedStartLabelId} && defined $matchPairs{$bedEndLabelId}) )  {
			if ( ($molrcmapEntry{'SiteID'} < $matchPairs{$bedStartLabelId}) || ($molrcmapEntry{'SiteID'} > $matchPairs{$bedEndLabelId}) ) {
				$genomeMapAvgCov = $genomeMapAvgCov + $molrcmapEntry{'Coverage'};	
				$count++;			
			}
		}
	}
	if ($count != 0) {
		$genomeMapAvgCov = ($genomeMapAvgCov/$count);		
	}
	
	# load stuff into hashes
	{
		lock($genomeLabelsLock);
		if ($genomeMapId != 0) {
			$genomeMapIds{$bedEntryCount} = $genomeMapId;
		}
		
		if (exists $matchPairs{$bedStartLabelId}) {
			$genomeMapStartLabelIds{$bedEntryCount} = $matchPairs{$bedStartLabelId};
		}
		else {
			$genomeMapStartLabelIds{$bedEntryCount} = 0;
		}
		
		if (exists $matchPairs{$bedEndLabelId}) {
			$genomeMapEndLabelIds{$bedEntryCount} = $matchPairs{$bedEndLabelId};
		}
		else {
			$genomeMapEndLabelIds{$bedEntryCount} = 0;
		}
		
		$genomeMapStartPoss{$bedEntryCount} = $genomeMapStartPos;
		$genomeMapEndPoss{$bedEntryCount} = $genomeMapEndPos;
		
	}
				
	my @allIds;
	my %confHash;
	
	# scan alignmol xmap to find molecules that align to identified genome map label IDs
	MOLLOOP: foreach my $molxmapEntry_ref (defined($molxmap[$genomeMapId])?@{$molxmap[$genomeMapId]}:()) {
		#print "Alignmol Loop:\n";
		my %molxmapEntry = %{get_hash_from_arr_json_by_format($molxmapEntry_ref, \@xmap_format)};
		if ($molxmapEntry{'RefContigID'} eq $genomeMapId ) {
			#print "Looking at moleculeID $molxmapEntry{'QryContigID'} for $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount}\n";
			my $molline = $molxmapEntry{'Alignment'};
			my @mols = split(/(?<=\))/, $molline);
			my $molrefLabels = "";
			my $molqryLabels = "";
			my %molmatchPairs;
		
			# process alignment doublet string
			foreach my $molpair (@mols) {
				#print "Doublet: $pair\n";
				my @molmatch = $molpair =~ /(\d+)/g;
				#print "Ref: $match[0] Qry: $match[1]\n";
				$molrefLabels = $molrefLabels."$molmatch[0] ";
				$molqryLabels = $molqryLabels."$molmatch[1] ";
				$molmatchPairs{$molmatch[0]} = $molmatch[1];
			}
						
			#print "\tmolRefLabels: $molrefLabels\n";
			my @molRefLabs = split(" ",$molrefLabels);
			my $firstLabId = $molRefLabs[0];
			my $lastLabId = $molRefLabs[-1];
				
			my $low = $genomeMapStart - 1;
			my $high = $genomeMapEnd + 1;
				
			#my ($genomeMapStartLabelId, $genomeMapEndLabelId, $bedStartLabelId, $bedEndLabelId, $genomeMapStartPos, $genomeMapEndPos);
			my ($genomeMapStartLabelId, $genomeMapEndLabelId, $bedStartLabelId, $bedEndLabelId);
			{
				lock($genomeLabelsLock);
				$genomeMapStartLabelId = $genomeMapStartLabelIds{$bedEntryCount};
				$genomeMapEndLabelId = $genomeMapEndLabelIds{$bedEntryCount};
				$bedStartLabelId = $bedStartLabelIds{$bedEntryCount};
				$bedEndLabelId = $bedEndLabelIds{$bedEntryCount};
				#$genomeMapStartPos = $genomeMapStartPoss{$bedEntryCount};
				#$genomeMapEndPos = $genomeMapEndPoss{$bedEntryCount};
			}
						
			if ($genomeMapStart < ($genomeMapStartLabelId - 25)) {
				$low = $genomeMapStartLabelId - 25;
			}
			if ($genomeMapEnd > ($genomeMapEndLabelId + 25)) {
				$high = $genomeMapEndLabelId + 25;
			}
											                         
			if ( !(($genomeMapStartLabelId != $firstLabId) && ($genomeMapStartLabelId != $lastLabId) && ($genomeMapEndLabelId != $firstLabId) && ($genomeMapEndLabelId != $lastLabId)) ) {
				next MOLLOOP;
			}
						
			ILOOP: for (my $i=$genomeMapStartLabelId; $i>=$low; $i--) {
				JLOOP: for (my $j=$genomeMapEndLabelId; $j<=$high; $j++) {								
					if ($j != $i) {
						if (($molrefLabels =~ /\b$i\b/) && ($molrefLabels =~ /\b$j\b/)) {	
							if ( ($genomeMapStartLabelId != $firstLabId) && ($genomeMapStartLabelId != $lastLabId) && ($genomeMapEndLabelId != $firstLabId) && ($genomeMapEndLabelId != $lastLabId) ) {
								my $leftLabelCount = 0;
								my $rightLabelCount = 0;
								my $conf = 0;
								for (my $a=0; $a<scalar(@molRefLabs); $a++) {
									if ($molRefLabs[$a] <= $genomeMapStartLabelId) {
										if (exists $molmatchPairs{$molRefLabs[$a]}) {
											$leftLabelCount++;
										}
									}
									elsif (($molRefLabs[$a] >= $genomeMapEndLabelIds{$bedEntryCount})) {
										if (exists $molmatchPairs{$molRefLabs[$a]}) {
											$rightLabelCount++;
										}
									}
								}
												
								if ($leftLabelCount <= $rightLabelCount) {
									$conf = $leftLabelCount;
								}
								else {
									$conf = $rightLabelCount;
								}
												
								if ($conf > 10) { $conf = 10; }
								if ($conf < 2) { next MOLLOOP; }
								
								my %allIdsHash = map { $_ => 1 } @allIds;
								if( !exists($allIdsHash{$molxmapEntry{'QryContigID'}}) ) {								
									push @allIds, $molxmapEntry{'QryContigID'};
								}
												
								if (!exists $confHash{$molxmapEntry{'QryContigID'}}) {										
									$confHash{$molxmapEntry{'QryContigID'}} = $conf;
								}
								else {
									next MOLLOOP;
								}																			
							}
						}
					}
				}
			}
		}
	} # end foreach molecule xmap entry
	
	my $totalConf=0;
	my @ngsIds;
	
	# if NGS scoring enabled
	my $ngsCount=0;
	if ($ngsScoring==1) {
		#$sam->clone;
		#my $sam = Bio::DB::Sam->new( -fasta=>"$inputs{fasta}", -fai=>"$fai", -bam=>"$inputs{bam}", -bai=>"$bai" );
		#my $sam = Bio::DB::Sam->new( -fasta=>"$inputs{fasta}", -bam=>"$inputs{bam}" );
		my $seqName;
		my $seqLen=0;
		foreach my $keyEntry_ref (@keys) {
			my %keyEntry = %$keyEntry_ref;
			#print "bedCmapId: $bedCmapId CompntId: $keyEntry{CompntId} CompntName: $keyEntry{CompntName}\n";
			if ($keyEntry{CompntId} eq $bedCmapId) {
				$seqName = $keyEntry{CompntName};
				$seqLen = $keyEntry{CompntLength};
				#print "Looking at $seqName with length $seqLen\n";
			}
		}
		my $bamStart = $bedStart - $scoreBuffer - 1000000; if ($bamStart < 1) { $bamStart = 1; }
		my $bamEnd = $bedEnd + $scoreBuffer + 100000; if ($bamEnd > $seqLen) { $bamEnd = $seqLen; }
		
		#print "Looking for seqName: $seqName\n";
		
		my @targets = $sam->seq_ids;
		my $seqId="";
		foreach my $i (@targets) {
			#print "\tNGS SeqId: $i\n";
			if ( ($i =~ /$seqName/i) || ($seqName =~ /$i/i) ) {
				$seqId = $i;
				#print "\tFound NGS SeqId: $i\n";
			}
		}
		
		#print "Querying for NGS $seqId start: $bamStart end: $bamEnd\n";
		my @alignments = $sam->get_features_by_location(-seq_id => "$seqId", -start => "$bamStart", -end => "$bamEnd");
		#print "\tFound ".scalar(@alignments)." NGS alignments\n";
		foreach my $a (@alignments) {
			#my $seqid  = $a->seq_id;
			my $start  = $a->start;
			my $end    = $a->end;
			#my $strand = $a->strand;
			#my $ref_dna= $a->dna;

			#my $query_start  = $a->query->start;
			#my $query_end    = $a->query->end;
			#my $query_strand = $a->query->strand;
			#my $query_dna    = $a->query->dna;
			my $read_name = $a->query->name;
		   
			#my $cigar     = $a->cigar_str;
			#my @scores    = $a->qscore;     # per-base quality scores
			#my $score     = $a->qstring;    # TAM-style quality string
			my $match_qual= $a->qual;       # quality of the match

			#my $paired = $a->get_tag_values('PAIRED');
			
			#print "\t\t$start $end $match_qual\n";
			
			if (($start < ($bedStart-$scoreBuffer)) && ($end > ($bedEnd+$scoreBuffer)) && ($match_qual >= 30)) {
				$ngsCount++;
				$totalConf = $totalConf + $ngsBonus;
				push @ngsIds, $read_name;
				#print "\t\tFOUND: $start $end $match_qual\n";
			}
		}
	}
				
	my $totalHits = scalar(@allIds);
	my $totalMol = scalar(unique(@allIds));

	
    foreach my $id (keys %confHash) {
		$totalConf = $totalConf + $confHash{$id};
	}

	my $avgConfPerHit = 0;
	if ($totalConf != 0 && $totalHits != 0) {
		$avgConfPerHit = $totalConf/($totalHits+$ngsCount);
	}
	else { $avgConfPerHit = 0; }
				
	my $avgConfPerMol = 0;
	if ($totalConf != 0 && $totalMol != 0) {
		$avgConfPerMol = $totalConf/($totalMol+$ngsCount);
	}
	else { $avgConfPerMol = 0; }
				
	my $score=0;
	#if ($totalConf != 0) { $score = $totalConf; }
	my $altScore = $totalConf;
	my $log10Score = log10($score) * 10;
	my $logScore = 0;
	if ($score != 0) { $logScore = log($score) * 10; }
	#$score = (sprintf "%.2f", $log10Score);
	$score = (sprintf "%.2f", ($totalConf/10) );
	{
		lock($output_lock);
		lock($genomeLabelsLock);
		print "\nBedLine: $bedEntryCount RefId: $bedCmapId Start-End: $bedStart-$bedEnd GenomeMapId: $genomeMapIds{$bedEntryCount} Start-End: $genomeMapStartPoss{$bedEntryCount}-$genomeMapEndPoss{$bedEntryCount}\n\tRefLabels: $bedStartLabelId $bedEndLabelId\tQryLabels: $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount}\n\tTotalHits: $totalHits UniqueMolecules: $totalMol TotalNGSHits: $ngsCount TotalConf: $totalConf AvgConfPerHit: $avgConfPerHit AvgConfPerMol $avgConfPerMol\n";
		print "\t\tRecommended score: $score TotalConf: $totalConf log10(score): $log10Score log(score): $logScore\n\n";
	}
				
	my $bedLineOut = "$bedCmapId\t$bedStart\t$bedEnd\t$bedEntry{Type}\t$score\t$bedEntry{Strand}\t$bedEntry{ThickStart}\t$bedEntry{ThickEnd}\t$bedEntry{ItemRgba}";
	if (defined $bedEntry{Sequence}) {
		$bedLineOut = $bedLineOut."\t$bedEntry{Sequence}";
	}
	$scoredBeds{$bedEntryCount} .= $bedLineOut."\n";
				
	$scoreStats{$bedEntryCount} .= "$bedEntryCount,$bedEntry{Type},".abs($bedEntry{ThickStart} - $bedEntry{ThickEnd}).",$bedCmapId,$bedStart-$bedEnd,$bedStartLabelId-$bedEndLabelId,$genomeMapIds{$bedEntryCount},$genomeMapStartPoss{$bedEntryCount}-$genomeMapEndPoss{$bedEntryCount},$genomeMapStartLabelIds{$bedEntryCount}-$genomeMapEndLabelIds{$bedEntryCount},$genomeMapAvgCov,$totalHits,$totalMol,$ngsCount,$totalConf,$avgConfPerHit,$avgConfPerMol,$score,$altScore,".join(" ",unique(@allIds)).",".join(" ",unique(@ngsIds))."\n";
} 

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

sub log10 {
    my $n = shift;
	my $out=0;
    if ($n != 0) {
		$out = log($n)/log(10);
	}
	return $out;
}

sub get_json_by_format {
  my ($str, $format) = @_;
  my @str = split(/\t/, $str);
  my %res;
  for (my $i=0; $i<scalar(@str); $i++) {
    $res{$$format[$i]} = $str[$i];
  }
  return encode_json(\%res);
}

sub get_arr_json_by_format {
  my ($str, $format, $key_index) = @_;
  $key_index = defined($key_index)?$key_index:2;
  my @str = split(/\t/, $str);
  return (encode_json(\@str), $str[$key_index]);
}

sub get_hash_from_arr_json_by_format {
	my ($json, $format) = @_;
	my @str = @{decode_json($json)};
	my @form = @{ $format };
	my %res;
	#print "\tScalar(format): ".scalar(@form)." Scalar(str): ".scalar(@str)."\n";
	for (my $i=0; $i<scalar(@form); $i++) {
		$res{$$format[$i]} = $str[$i];
	}
	
	return \%res;
}
