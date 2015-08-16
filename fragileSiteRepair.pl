#!/usr/bin/perl

# A wrapper script to run fragile site repair on assembly alignments to a reference using an input FASTA and alignref_final XMAP and CMAPs
# If .bed file of fragile sites is not provided, an input (reference) fasta is expected, from which potential fragiles sites will be calculated. 
# Assembly maps are stitched together based on alignments if the maps start and stop overlap a fragile site

# Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bed>] [--cmap <assembly.cmap>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> [--ref <reference.cmap>] --output <output folder> [--maxlab <max_label_gap_tolerence = 1>] [--maxfill <max basepairs to fill between contigs = 30000>] [--wobble <fragile site wobble in bp = 30000>] [--seq <sequence bp +/- fragile site to print in final BED = 50>] [--force <overwrite output folder>] [--n <number of CPU cores>] [--j <number of parallel jobs>] [--optArgs <optArguments_human.xml>] [--runSV <enable run SV detection>] [--bnx <input.bnx to run alignmol>]

use strict;
use warnings;
use Getopt::Long; 
use Cwd qw(abs_path);
use File::Basename;
use File::Path qw(rmtree mkpath);
use File::Copy;
use File::Find;
use IPC::System::Simple qw(system capture);
#use Parallel::Iterator qw(iterate_as_array);
#use Parallel::Loops;
use Parallel::ForkManager;
use DateTime;
use DateTime::Format::Human::Duration;
#use Data::Dumper;

my $dtStart = DateTime->now;
my $stime = DateTime->now;
my $etime = DateTime->now;
my $dtime = DateTime::Format::Human::Duration->new();
my $cpuCount=4;
my $jobs=1;
my $mem = 32;

print "\n";
print qx/ps -o args $$/;
print "\n";

# << usage statement and variable initialisation >>
my %inputs = (); 
$inputs{'force'}=0;
GetOptions( \%inputs, 'fasta:s', 'xmap=s', 'qcmap=s', 'rcmap=s', 'errbin=s', 'output=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'force', 'bed:s', 'seq:i', 'cmap:s', 'ref:s', 'n:i','j:i','optArgs:s', 'runSV', 'bnx:s', 'step:i'); 

my $hasbed = 0; 
my $getSeq = 0;

if ( (!exists $inputs{fasta} & !exists $inputs{bed}) | !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{errbin} | !exists $inputs{output} ) {
	print "Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bed>] [--cmap <assembly.cmap>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> [--ref <reference.cmap>] --output <output folder> [--maxlab <max_label_gap_tolerence = 1>] [--maxfill <max basepairs to fill between contigs = 30000>] [--wobble <fragile site wobble in bp = 30000>] [--seq <sequence bp +/- fragile site to print in final BED = 50>] [--force <overwrite output folder>] [--n <number of CPU cores>] [--j <number of parallel jobs>] [--optArgs <optArguments_human.xml>] [--runSV <enable run SV detection>] [--bnx <input.bnx to run alignmol>]\n"; 
	exit 0; 
}
else {
	print "\nStart time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";
	
	my $hostname = `hostname`;
	chomp($hostname);
	print "Running on host: $hostname\n";

	#get number of CPU cores
	use Sys::Info;
	use Sys::Info::Constants qw( :device_cpu );
	my $info = Sys::Info->new;
	my $cpu  = $info->device( CPU => my %options );
	#printf "There are %d CPUs\n"  , $cpu->count || 1;
	$cpuCount = $cpu->count;
	$jobs = int($cpuCount/6);
	#print "CPU count: $cpuCount\n";
	# get memory info
	use Sys::MemInfo qw(totalmem freemem totalswap);
	$mem = (((&totalmem / 1024) / 1024)) / 1024; 
	$mem = int($mem);
	if (exists $inputs{n}) { $cpuCount = int($inputs{n}); }
	if (exists $inputs{j}) { $jobs = int($inputs{j}); }

	print "Maximum CPU cores to use: $cpuCount\n";
	print "Maximum parallel jobs: $jobs\n";
	print "Maximum memory to use: ", int($mem)," GB\n";
	print "\n";
	
	foreach my $key ("xmap","qcmap","rcmap","errbin","output","bed","fasta", "cmap", "ref", "optArgs", "bnx") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
	if( !exists $inputs{maxfill} ) { $inputs{maxfill} = 30000; }
	if( !exists $inputs{wobble} ) { $inputs{wobble} = 30000; }
	if( !exists $inputs{maxlab} ) { $inputs{maxlab} = 1; }
		else {
			if ($inputs{maxlab} > 1) {
				die "ERROR: Maximum recommended --maxlab = 1\n";
			}
		}
	
	if( exists $inputs{seq} ) { 
		$getSeq = 1; 
	}

	if( exists $inputs{bed} ) {
		$hasbed = 1; 
		if( exists $inputs{fasta} ) {
			print "$inputs{bed} fragile sites file already provided.\n$inputs{fasta} will be ignored.\n";
		}
	} 
}
my $splitprefix = basename(abs_path($inputs{xmap}), ".xmap");

#print out input variables
if( $hasbed==0 ) {
	print "Input FASTA: $inputs{fasta}\n"; 
	if ( $getSeq == 1) {
		print "\tSequence bp +/- start/end of fragile to print to BED: $inputs{seq}\n";
	}
}
else { print "Input BED: $inputs{bed}\n";  }

print "Input XMAP: $inputs{xmap}\n";
print "Input QCMAP: $inputs{qcmap}\n";
print "Input RMCAP: $inputs{rcmap}\n";
print "Input ERRBIN: $inputs{errbin}\n";
if (exists $inputs{cmap}) {print "Input CMAP: $inputs{cmap}\n";} else {print "Input original assembly CMAP: not provided\n";}
if (exists $inputs{ref}) {print "Input reference CMAP: $inputs{ref}\n";} else {print "Input reference CMAP: not provided\n";}
if (exists $inputs{bnx}) {print "Input BNX: $inputs{bnx}\n";} else {print "Input BNX: not provided\n";}
print "Output folder: $inputs{output}\n\n";
print "Maximum labels between maps: $inputs{maxlab}\n";
print "Maximum basepairs to fill between maps: $inputs{maxfill}\n";
print "Maximum fragile site wobble: $inputs{wobble}\n";
print "\n";

if (exists $inputs{runSV} && exists $inputs{ref}) { print "SV module enabled\n\n"; }

# Make sure dependency scripts exist
my $scriptspath = abs_path(dirname($0));
if (!-e "$scriptspath/calcFragileSites.pl" | !-e "$scriptspath/split_xmap_standalone.pl" | !-e "$scriptspath/gapFill.pl" | !-e "$scriptspath/extractFASTA_fromBED.pl" | !-e "$scriptspath/runCharacterizeFinal.py" | !-e "$scriptspath/callSV.pl" | !-e "$scriptspath/scoreStitchPositionsBED_v3.pl" | !-e "$scriptspath/cutCmapByBedScore.pl") {
	die "ERROR: Dependency scripts not found at $scriptspath\n"; 
}
if( !-e $ENV{"HOME"}."/tools/RefAligner" | !-e $ENV{"HOME"}."/scripts/HybridScaffold/scripts/calc_cmap_stats.pl" | !-e $ENV{"HOME"}."/scripts/optArguments_human.xml" | !-e $ENV{"HOME"}."/scripts/runSV.py") {
	# die "ERROR: Pipeline tools and scripts not found. Please ensure ~/tools and ~/scripts are properly set up.\n"; 
	die "abs_path('~/tools/RefAligner') not found: $!\n"; 
}
if( !exists $inputs{optArgs} ) { $inputs{optArgs} = abs_path($ENV{"HOME"}."/scripts/optArguments_human.xml"); }


# check output folder
if (-d $inputs{output} and $inputs{force} eq 0) {
	die "ERROR: Output directory $inputs{output} exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $inputs{output} and $inputs{force} eq 1) {
	print "WARNING: Output directory $inputs{output} exists and will be overwritten\n";
	rmtree($inputs{output}) or die "ERROR: Cannot remove existing output directory $inputs{output}: $!\n";
}
mkpath($inputs{output}) or die "ERROR: Cannot create output directory $inputs{output}: $!\n";
print "\n";



# Step 1: Split input XMAP into individual anchor maps
print "===Step 1: Split merged input XMAP into individual anchor maps===\n";
# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
my $cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $splitprefix"."_contig $inputs{output}/contigs";
print "Running command: $cmd\n";
system($cmd);
print "\n";


# Step 2: Calculate fragile sites for input FASTA
print "===Step 2: Calculate fragile sites for input FASTA===\n";
my $bed; 
if ( $hasbed==1 and -e $inputs{bed}) {
	print "Fragiles sites BED file $inputs{bed} provided.\nSkipping fragile sites calculation...\n\n";
	$bed = $inputs{bed};
	copy("$inputs{bed}","$inputs{output}") or print "Copy of input BED $inputs{bed} failed: $!\n";
}
else {
	$stime = DateTime->now;
	# Usage: perl calcFragilesSites.pl <input FASTA> [output.bed] [sequence bp +/- fsite to print out]
	$cmd = "perl $scriptspath/calcFragileSites.pl $inputs{fasta} $inputs{output}";
	if ($getSeq == 1) {
		$cmd = $cmd." $inputs{seq}";
	}
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	# my $bed = findBED($inputs{output});
	$bed = findBED($inputs{output});
	$bed = abs_path($inputs{output}."/".$bed);
	#print "BED file: $bed\n";
	
	$etime = DateTime->now;
	$dtime = DateTime::Format::Human::Duration->new();
	print 'Spent time at this step: ', $dtime->format_duration_between($etime, $stime); print "\n\n";
}

# Step 3: Run gapFill.pl for each anchor map
print "===Step 3: Run gapFill.pl for each anchor map===\n";
print "Processing anchor maps in parallel fashion...\n";
print "\n";
$stime = DateTime->now;

my @xmaps = findXMAPs($inputs{output}."/contigs");
@xmaps = sort @xmaps;

# normal sequential processing

#foreach my $xmap (@xmaps) {
	#$xmap = abs_path($inputs{output}."/contigs/$xmap");
	#my $base = $xmap; $base =~ s/.xmap//i;
	##print "Base: $base\n";
	#my $qcmap = $base."_q.cmap";
	#my $rcmap = $base."_r.cmap";
	#if (-e $xmap and -e $qcmap and -e $rcmap) {
		#print "XMAP: $xmap\n";
		#print "QCMAP: $qcmap\n";
		#print "RCMAP: $rcmap\n";
		#print "\n";
		
		## usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		#$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -e $inputs{errbin} -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble}";
		#print "\tRunning command: $cmd\n";
		#print "\n";
		##system($cmd) or die "ERROR: $cmd failed: $!\n";
		#my @result = capture($cmd);
		#my $log = "$base"."_fragileSiteRepaired_log.txt";
		##print join("\t", @result);
		#open LOG, ">$log" or die "ERROR: Cannot open $log for writing! $!\n";
		#print LOG join("\t", @result);
		#close LOG;
	#}	
#}


# using Parallel::ForkManager

#my $pm =  new Parallel::ForkManager($cpuCount);
my $pm =  new Parallel::ForkManager($jobs);
foreach my $xmap (@xmaps) {
	$pm->start and next; # do the fork
	# do work
	$xmap = abs_path($inputs{output}."/contigs/$xmap");
	my $base = $xmap; $base =~ s/.xmap//i;
	#print "Base: $base\n";
	my $qcmap = $base."_q.cmap";
	my $rcmap = $base."_r.cmap";
	if (-e $xmap and -e $qcmap and -e $rcmap) {
		print "XMAP: $xmap\n";
		print "QCMAP: $qcmap\n";
		print "RCMAP: $rcmap\n";
		print "\n";
		
		# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -e <errbin> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>] [--n CPU cores to use]
		$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -e $inputs{errbin} -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble} --n $cpuCount";
		print "\tRunning command: $cmd\n";
		print "\n";
		#system($cmd) or die "ERROR: $cmd failed: $!\n";
		my @result = capture($cmd);
		my $log = "$base"."_fragileSiteRepaired_log.txt";
		#print join("\t", @result);
		open LOG, ">$log" or die "ERROR: Cannot open $log for writing! $!\n";
		print LOG join("\t", @result);
		close LOG;
		print "\n";
	}
	$pm->finish; # do the exit in the child process
}
print "Waiting for all processes to finish...\n";
print "\n";
$pm->wait_all_children;
# clean up temp files
find(\&wanted, "$inputs{output}/contigs"); 
find(\&wanted2, "$inputs{output}/contigs"); 
find(\&wanted3, "$inputs{output}/contigs"); 

$etime = DateTime->now;
$dtime = DateTime::Format::Human::Duration->new();
print 'Spent time at this step: ', $dtime->format_duration_between($etime, $stime); print "\n\n";

# using Parallel::Loops

#my $pl = Parallel::Loops->new($cpuCount);

#$pl->foreach( \@xmaps, sub {
	#my $xmap = $_;
	#$xmap = abs_path($inputs{output}."/contigs/$xmap");
	#my $base = $xmap; $base =~ s/.xmap//i;
	##$xmap = "$inputs{output}/contigs/$xmap";
	#my $qcmap = $base."_q.cmap";
	#my $rcmap = $base."_r.cmap";
	#if (-e $xmap and -e $qcmap and -e $rcmap) {
		#print "XMAP: $xmap\n";
		#print "QCMAP: $qcmap\n";
		#print "RCMAP: $rcmap\n";
		#print "\n";
		
		## usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		#$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble}";
		#print "\tRunning command: $cmd\n";
		#print "\n";
		##system($cmd) or die "ERROR: $cmd failed: $!\n";
		#my @result = capture($cmd);
		#my $log = "$base"."_fragileSiteRepaired_log.txt";
		##print join("\t", @result);
		#open LOG, ">$log" or die "ERROR: Cannot open $log for writing! $!\n";
		#print LOG join("\t", @result);
		#close LOG;
	#}	
#});


# Step 4: Merge individual fragileSiteRepaired anchor maps and stitched BEDs
print "\n===Step 4: Merge individual fragileSiteRepaired anchor maps===\n\n";
my @qcmaps = findQCMAPs($inputs{output}."/contigs");
@qcmaps = sort @qcmaps;

#generate file with list of maps to be merged
my $mergeFile = "$inputs{output}/contigs/mergeList.txt";
open (LIST, ">$mergeFile") or die "ERROR: Could not open $mergeFile: $!\n";
print "Merging ".scalar(@qcmaps)." fragileSiteRepaired anchor maps...\n\n";
foreach (@qcmaps) {
	#print "QCMAP: $_\n";
	$_ = abs_path($inputs{output}."/contigs/$_");
	print LIST "$_\n";
}
close LIST;

#my $input = join(" -i ",@qcmaps);
#$cmd = "cd $inputs{output}; ~/tools/RefAligner -i $input -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0";
$cmd = "cd $inputs{output}; ~/tools/RefAligner -f -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0";
print "Running command: $cmd\n";
print "\n";
system($cmd);
my $finalmap = "$inputs{output}/$splitprefix"."_fragileSiteRepaired.cmap";
print "\n";

print "Generating merged stitchPositions BED...";
my $mergedBED = ""; 
if ($getSeq == 1) {
	$mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.withSeqs.bed";
}
else {
	$mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.bed";
}
open BEDOUT, ">$mergedBED" or die "ERROR: Could not open $mergedBED: $!\n";
my @beds = findBEDs("$inputs{output}/contigs");
my @bedOut;
#print BEDOUT "#CMapId\tStart\tEnd\tType\tSequence\n";
if ($getSeq == 1) {
	print BEDOUT "#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence\n";
}
else {
	print BEDOUT "#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba\n";
}
foreach my $bedFile (@beds) {
	open BEDFILE, "<$bedFile" or die "ERROR: Could not open $bedFile: $!\n";
	while (<BEDFILE>) {
		if ($_ =~ /^#/) {
			next;
		}
		chomp($_);
		#print BEDOUT "$_\n";
		push @bedOut, $_;
	}
	close BEDFILE; 
}
my ($bedOut_ref, $typeCount_ref) = sortBEDandCount(\@bedOut);
@bedOut = @$bedOut_ref;
my %typeCount = %$typeCount_ref;
my $stitchCount = scalar(@bedOut);
print BEDOUT join("\n",@bedOut);
close BEDOUT; 
print "done\n\n";


# Step 5: Generate new merged cmap with unmapped maps
print "\n\n===Step 5: Process original assembly CMAP and generate new merged CMAP ===\n\n";
	if (exists $inputs{cmap}) {
	my @origCmapIds = getCmapIds($inputs{cmap});
	my @qryCmapIds = getCmapIds($inputs{qcmap});
	my %in_qryCmapIds = map {$_ => 1} @qryCmapIds;
	my @diff  = grep {not $in_qryCmapIds{$_}} @origCmapIds;
	# my $oldMap = abs_path($inputs{output}."/".$splitprefix."_fragileSiteRepaired.cmap");
	my $oldMap = $finalmap; 
	my $tempMap = $inputs{output}."/temp";
	my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged";
	#print "\tcmapIds in orig but not ref: ".join("\n\t",@diff)."\n";
	if( scalar(@diff) == 0 ) { 
		print "No unmapped cmaps were found in assembly.cmap. No merges necessary. \n$oldMap will be copied as is to $newMap.cmap\n"; 
		copy( $oldMap, $newMap.".cmap" ) or print "WARNING: Unable to copy $oldMap to $newMap.cmap. Step 7 may fail\n"; 
	}
	else {
		$cmd = "cd $inputs{output}; ~/tools/RefAligner -f -merge -i $inputs{cmap} -selectid ".join(" ",@diff)." -o $tempMap";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n";
		$cmd = "cd $inputs{output}; ~/tools/RefAligner -f -merge -i $tempMap".".cmap -i $oldMap -o $newMap; rm -f $tempMap".".cmap";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n";
	}
	$finalmap = $newMap.".cmap";
}
else {
	print "Skipping Step 5. Original assembly CMAP not provided\n";
	print "\n";
}


# Step 6: run CharacterizeFinal on newly merged cmap
print "\n\n===Step 6: Align fragileSiteRepaired merged CMAP to reference and get stats ===\n\n";
if (-e $finalmap && (exists $inputs{ref})) {
	copy("$inputs{ref}","$inputs{output}") or print "WARNING: Copy of input reference $inputs{ref} failed: $!\n"; 
	# my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged.cmap";
	#usage: runCharacterizeFinal.py [-h] [-t REFALIGNER] [-r REFERENCEMAP] [-q QUERYMAP] [-x XMAP] [-p PIPELINEDIR] [-a OPTARGUMENTS] [-n NUMTHREADS] 
	# $cmd = "cp runCharacterizeFinal.py ~/scripts/; python ~/scripts/runCharacterizeFinal.py -t ~/tools/RefAligner -r $inputs{ref} -q $newMap -p ~/scripts/ -n $cpuCount";
	$cmd = "cp $scriptspath/runCharacterizeFinal.py ~/scripts/; python ~/scripts/runCharacterizeFinal.py -t ~/tools/RefAligner -r $inputs{ref} -q $finalmap -p ~/scripts/ -n $cpuCount";
	print "\tRunning command: $cmd\n\n";
	system($cmd);
}
else {
	print "Skipping Step 6. Original assembly CMAP or reference CMAP not provided\n";
	print "\n";
}
my $alignreffinalXmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged.xmap";
my $alignreffinalQcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_q.cmap";
my $alignreffinalRcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_r.cmap";


# Step 7: align molecules to merged fragileSiteRepaired CMAP
print "\n\n===Step 7: Run single molecule alignments against fragileSiteRepaired merged cmap ===\n\n";
my $alignmolXmap="";
if (defined $inputs{bnx} && exists $inputs{bnx}) {
	print "Running alignment of $inputs{bnx} to $finalmap\n\n";
	mkpath("$inputs{output}/alignmol") or die "ERROR: Cannot create output directory $inputs{output}/alignmol: $!\n";
	my $alignmolDir = abs_path("$inputs{output}/alignmol");
	# $ cd /home/palak; /home/palak/tools/RefAligner.mic -ref /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -o /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/contigs/exp_refineFinal1/alignmol/EXP_REFINEFINAL1_1 -i /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/all_1_of_30.bnx -f -stdout -stderr -maxthreads 240 -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -maxmem 7
	$cmd = "cd $alignmolDir; ~/tools/RefAligner -f -ref $finalmap -o $splitprefix"."_fragileSiteRepaired_merged_alignmol -i $inputs{bnx} -maxthreads $cpuCount -maxmem ".($mem)." -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt\$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -stdout -stderr";
	print "\tRunning command: $cmd\n\n";
	system($cmd);
	sleep(5);
	
	#split xmap into _contig maps
	$alignmolXmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol.xmap";
	my $alignmolQmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol_q.cmap";
	my $alignmolRmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol_r.cmap";
	my $alignmolPrefix = "$splitprefix"."_fragileSiteRepaired_merged_alignmol_contig";
	#mkpath("$alignmolDir/contigs") or die "ERROR: Cannot create output directory $inputs{output}/alignmol/contigs: $!\n";
	# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
	$cmd = "perl $scriptspath/split_xmap_standalone.pl $alignmolXmap $alignmolQmap $alignmolRmap $alignmolPrefix $alignmolDir/contigs";
	print "\n\tRunning command: $cmd\n\n";
	system($cmd);
	print "\n";
}
else {
	print "Skipping Step 7 Single molecule alignments. Input BNX not provided or does not exist\n";
	print "\n";
}

print "\n";


# Step 8: score BED file
print "\n\n===Step 8: Score stitchPositions BED based on single molecule alignments ===\n\n";

my $scoredBED = basename($mergedBED,".bed");
my $scoredSTATS = $scoredBED;
$scoredBED = $scoredBED."_scored.bed";
$scoredSTATS = $scoredSTATS."_scoreStats.csv";

if (defined $inputs{bnx} && exists $inputs{bnx} && -e $finalmap && exists $inputs{ref}) {
	print "Scoring BED file $mergedBED based on molecule alignments in $alignmolXmap\n\n";
	mkpath("$inputs{output}/scoreBED") or warn "WARNING: Cannot create output directory $inputs{output}/scoreBED: $!\nNext steps may fail...\n";
	$scoredBED = $inputs{output}."/scoreBED/$scoredBED";
	
	$cmd = "perl $scriptspath/scoreStitchPositionsBED_v3.pl --bed $mergedBED --xmap $alignreffinalXmap --qcmap $alignreffinalQcmap --rcmap $alignreffinalRcmap --alignmol $alignmolXmap --n $cpuCount --output $inputs{output}/scoreBED/ > $inputs{output}/scoreBED/".basename($mergedBED,".bed")."_scored_log.txt 2>&1";
	print "\n\tRunning command: $cmd\n\n";
	system($cmd);
	print "\n"; 
	$scoredSTATS = $inputs{output}."/scoreBED/$scoredSTATS";
}
else {
	print "Skipping Step 8 BED scoring. Input BNX, reference, or fragileSiteRepaired merged CMAP not provided or does not exist\n";
	print "\n";
}


# Step 9: cut fragileSiteRepaired CMAP where score is too low
print "\n\n=== Step 9: Breaking fragileSiteRepaired merged CMAP based on scores in BED ===\n\n";
my $pre = basename($finalmap, ".cmap")."_final";
#$pre = $pre."_final";
my $newfinalmap = "";
my $newBED = "$inputs{output}/scoreBED/".basename($scoredBED, ".bed");
if (defined $inputs{bnx} && exists $inputs{bnx} && -e $finalmap && exists $inputs{ref}) {
	if (-e $scoredBED && -e $scoredSTATS) {
		# Usage: perl cutCmapByBedScore.pl <--bed stitchPositions_scored.bed> <--stats stitchPositions_scoreStats.csv> <--cmap fragileSiteRepaired_merged.cmap> <--prefix output prefix> [--threshold score threshold] [--maxfill maxfill <minlen filter>]
		$cmd = "perl $scriptspath/cutCmapByBedScore.pl --bed $scoredBED --stats $scoredSTATS --cmap $finalmap --prefix $pre --maxfill $inputs{maxfill} --output $inputs{output}/scoreBED/";
		print "\n\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n"; 
		
		$newfinalmap = $inputs{output}."/scoreBED/".$pre.".cmap";
		$newBED = $newBED."_final.bed";
		
		copy("$newfinalmap","$inputs{output}") or print "Copy of final merged fragileSiteRepaired CMAP $newfinalmap failed: $!\n";
		$newfinalmap = abs_path($inputs{output}."/".$pre.".cmap");
		copy("$newBED","$inputs{output}") or print "Copy of final merged fragileSiteRepaired CMAP $newBED failed: $!\n";
		#$newBED = "$inputs{output}/".basename($scoredBED, ".bed")."_final.bed";
	}
	else {
		print "Skipping Step 9 cutting CMAP based on scores. Scored BED and/or stats file not provided or does not exist\n";
	}
}
else {
	print "Skipping Step 9 cutting CMAP based on scores. Input BNX, reference, or fragileSiteRepaired merged CMAP not provided or does not exist\n";
	print "\n";
}


# Step 10: run CharacterizeFinal on new final merged CMAP if needed
if (-e $newfinalmap) {
	# move/rename old alignref_final folder
	$cmd = "mv $inputs{output}/alignref_final $inputs{output}/alignref";
	system($cmd);
	
	print "\n\n===Step 10: Align final fragileSiteRepaired merged CMAP to reference and get stats ===\n\n";
	if (-e $newfinalmap && (exists $inputs{ref})) {
		#copy("$inputs{ref}","$inputs{output}") or print "WARNING: Copy of input reference $inputs{ref} failed: $!\n"; 
		$cmd = "cp $scriptspath/runCharacterizeFinal.py ~/scripts/; python ~/scripts/runCharacterizeFinal.py -t ~/tools/RefAligner -r $inputs{ref} -q $newfinalmap -p ~/scripts/ -n $cpuCount";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		
		$alignreffinalXmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final.xmap";
		$alignreffinalQcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final_q.cmap";
		$alignreffinalRcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final_r.cmap";
	}
	else {
		print "Skipping Step 10 final characterizeFinal. Final merged fragileSiteRepaired CMAP or reference CMAP not found\n";
		print "\n";
	}
}
else {
	print "Skipping Step 10 final characterizeFinal. Final merged fragileSiteRepaired CMAP not found and/or input BNX not provided\n";
}


# Step 11: run alignmol on new final merged CMAP if needed
if (-e $newfinalmap) {
	$alignmolXmap = "";
	mkpath("$inputs{output}/alignmol_final") or warn "WARNING: Cannot create output directory $inputs{output}/alignmol_final : $!\nNext steps may fail...\n"; 
	print "\n\n===Step 11: Run final single molecule alignments against final fragileSiteRepaired merged cmap ===\n\n";
	#$alignmolXmap="";
	if (defined $inputs{bnx} && exists $inputs{bnx}) {
		print "Running alignment of $inputs{bnx} to $newfinalmap\n\n";
		my $alignmolDir = abs_path("$inputs{output}/alignmol_final");
		$cmd = "cd $alignmolDir; ~/tools/RefAligner -f -ref $newfinalmap -o $splitprefix"."_fragileSiteRepaired_merged_final_alignmol -i $inputs{bnx} -maxthreads $cpuCount -maxmem ".($mem)." -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt\$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -stdout -stderr";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		sleep(5);
	
		#split xmap into _contig maps
		$alignmolXmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol.xmap";
		my $alignmolQmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_q.cmap";
		my $alignmolRmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_r.cmap";
		my $alignmolPrefix = "$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_contig";
		#mkpath("$alignmolDir/contigs") or die "ERROR: Cannot create output directory $inputs{output}/alignmol/contigs: $!\n";
		# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
		$cmd = "perl $scriptspath/split_xmap_standalone.pl $alignmolXmap $alignmolQmap $alignmolRmap $alignmolPrefix $alignmolDir/contigs";
		print "\n\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n";
	}
	else {
		print "Skipping Step 11 Single molecule alignments. Input BNX not provided or does not exist\n";
		print "\n";
	}
}
else {
	print "Skipping Step 11 final single molecule alignments. Input BNX and/or Final merged fragileSiteRepaired CMAP not found\n";
}

# Step 12: re-score final BED file if needed
if (-e $newfinalmap && -e $newBED) {
	print "\n\n===Step 12: Re-score final stitchPositions BED based on final single molecule alignments ===\n\n";

	my $scoredBED = basename($newBED,".bed");
	my $scoredSTATS = $scoredBED;
	$scoredBED = $scoredBED."_scored.bed";
	$scoredSTATS = $scoredSTATS."_scoreStats.csv";

	if (defined $inputs{bnx} && exists $inputs{bnx} && -e $newfinalmap && exists $inputs{ref}) {
		print "Scoring final BED file $newBED based on molecule alignments in $alignmolXmap\n\n";
		#mkpath("$inputs{output}/scoreBED") or warn "WARNING: Cannot create output directory $inputs{output}/scoreBED: $!\nNext steps may fail...\n";
		$scoredBED = $inputs{output}."/scoreBED/$scoredBED";
		my $oldBED = $newBED;
		$newBED = "$inputs{output}/".basename($scoredBED);
		
		$cmd = "perl $scriptspath/scoreStitchPositionsBED_v3.pl --bed $oldBED --xmap $alignreffinalXmap --qcmap $alignreffinalQcmap --rcmap $alignreffinalRcmap --alignmol $alignmolXmap --n $cpuCount --output $inputs{output}/scoreBED/ > $inputs{output}/scoreBED/".basename($newBED,".bed")."_log.txt 2>&1";
		print "\n\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n"; 
		$scoredSTATS = $inputs{output}."/scoreBED/$scoredSTATS";
		#$scoredBED = $inputs{output}."/scoreBED/$scoredBED";
		
		copy("$scoredBED","$inputs{output}") or print "WARNING: Copy of re-scored BED $scoredBED failed: $!\n"; 	
		
	}
	else {
		print "Skipping Step 12 BED scoring. Input BNX, reference, or final fragileSiteRepaired merged CMAP not provided or does not exist\n";
		print "\n";
	}
}
else {
	print "Skipping Step 12. Final merged fragileSiteRepaired CMAP and/or final BED not found\n";
}


# Step 13: Calculate stats for merged fragileSiteRepaired BED and CMAP
print "=== Step 13: Calculate stats for final merged fragileSiteRepaired CMAP ===\n";
print "\n"; 

#generate FASTA of stitchPositions sequences
if ($getSeq == 1) {
	if (-e $newBED) {
		$cmd = "perl $scriptspath/extractFASTA_fromBED.pl --bed $newBED";
		print "Generating FASTA from final scored merged stitchPositions BED...\n";
	}
	else {
		$cmd = "perl $scriptspath/extractFASTA_fromBED.pl --bed $mergedBED";
		print "Generating FASTA from merged stitchPositions BED...\n";
	}
	#print "Running command $cmd\n";	
	system($cmd);
}
print "\n";

#get stats of original input BED
my @origBed;
open ORIGBEDFILE, "<$bed" or die "ERROR: Could not open $bed: $!\n";
while (<ORIGBEDFILE>) {
		if ($_ =~ /^#/) {
			next;
		}
		chomp($_);
		push @origBed, $_;
}
close ORIGBEDFILE; 
my ($origBedOut_ref, $origTypeCount_ref) = sortBEDandCount(\@origBed);
my %origTypeCount = %$origTypeCount_ref;
my $origCount = scalar(@origBed);
print "Total potential fragile sites: $origCount\n";
foreach my $type (sort keys %origTypeCount) {
	print "\t$type: $origTypeCount{$type}\n";
}

#get stats of scored BED
if (-e $newBED) {
	my @newBed;
	open NEWBEDFILE, "<$newBED" or warn "WARNING: Could not open $newBED: $!\n";
	while (<NEWBEDFILE>) {
		if ($_ =~ /^#/) {
			next;
		}
		chomp($_);
		push @newBed, $_;
	}
	close NEWBEDFILE; 
	my ($newBedOut_ref, $newTypeCount_ref) = sortBEDandCount(\@newBed);
	my %newTypeCount = %$newTypeCount_ref;
	my $newCount = scalar(@newBed);
	
	print "Total fragile sites repaired (stitched) and scored: $newCount\n";
	foreach my $type (sort keys %newTypeCount) {
		print "\t$type: $typeCount{$type}\n";
	}
}
else {
	print "Total fragile sites repaired (stitched): $stitchCount\n";
	foreach my $type (sort keys %typeCount) {
		print "\t$type: $typeCount{$type}\n";
	}
}
print "\n";

# usage: calc_cmap_stats.pl <CMAP_File>
my $dir = glob("~/scripts/HybridScaffold/scripts");
my $script = "calc_cmap_stats.pl";
$cmd = "";
if (-e "$dir/$script") {
	if (exists $inputs{cmap}) {
		print "Original assembly CMAP $inputs{cmap} stats:\n";
		$cmd = "perl $dir/$script $inputs{cmap}";
	}
	else {
		print "Original query CMAP $inputs{qcmap} stats:\n";
		$cmd = "perl $dir/$script $inputs{qcmap}";
	}
	#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	# my $file = "$inputs{output}/$splitprefix"."_fragileSiteRepaired.cmap";
	#print "Fragile site repaired merged CMAP $finalmap stats:\n";
	#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
	if (-e $newfinalmap) {
		print "Fragile site repaired final merged CMAP $newfinalmap stats:\n";
		$cmd = "perl $dir/$script $newfinalmap";
	}
	else {
		print "Fragile site repaired merged CMAP $finalmap stats:\n";
		$cmd = "perl $dir/$script $finalmap";
	}
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	
}
else {
	print "Skipping Step 13: perl script calc_cmap_stats.pl not found at $dir/$script\n"; }
print "\n";


# Step 14: run SV detection on final cmap
print "\n\n===Step 14: Run SV detection on new merged cmap ===\n\n";
if (defined $inputs{runSV}) {
	if ((exists $inputs{cmap}) && (exists $inputs{ref})) {
		my $errbin = "";
		my $map = $finalmap;
		if (-e $newfinalmap) {
			$errbin = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final.errbin";
			$map = $newfinalmap;
		}
		else {
			$errbin = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged.errbin";
		}
		my $mres = "2.0";		
		
		if (-e $map && -e $errbin) {
			# Usage: callSV.pl --ref <reference.cmap> --cmap <input.cmap> --errbin <input.errbin> --optarg <optArguments.xml> --outdir <output_folder> [--prefix <contig_prefix>] [--mres <pixels_to_reduce_ref> [--force]
			$cmd = "perl callSV.pl --ref $inputs{ref} --cmap $map --errbin $errbin --optarg $inputs{optArgs} --outdir $inputs{output}/alignref_final_sv --mres $mres --n $cpuCount --prefix $splitprefix"."_fragileSiteRepaired_merged_final";
			if ($inputs{force} eq 1) {$cmd = $cmd." --force";}
			print "\tRunning command: $cmd\n\n";
			system($cmd);
			
			#print SV stats summary
			print "==SV Summary ==\n\n";
			#$cmd = "grep -A 25 -i \"SV type\" $inputs{output}/alignref_final_sv/sv_log.txt";
			$cmd = "cat $inputs{output}/alignref_final_sv/sv_log.txt | sed -n -e \'/SV type/,\$p\'";
			my @SVresults = capture($cmd); 
			print join("",@SVresults);
		}
		else { print "WARNING: Merged CMAP $map and/or alignref_final errbin $errbin not found! Skipping Step 8 runSV...\n"; }
	}
	else {
		print "Skipping Step 14. Original assembly CMAP and/or reference CMAP not provided\n";
		print "\n";
	}
}
else {
	print "Skipping Step 14 SV detection. Add --runSV to enable SV detection\n";
	print "\n";
}
print "\n";











# << Parting words >> 
print "\n";
my $dtEnd = DateTime->now;
print "End time: "; print join ' ', $dtEnd->ymd, $dtEnd->hms; print "\n\n";
my $span = DateTime::Format::Human::Duration->new();
print 'Total elapsed time: ', $span->format_duration_between($dtEnd, $dtStart); print "\n\n";






sub sortBEDandCount {
	my @bedIn = @{$_[0]};
	#create array of arrays
	my @AoA;
	my %typeCount;
	foreach my $line (@bedIn) {
		chomp($line);
		#print "$line\n";
		my @s = split("\t",$line);
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		if (scalar(@s) == 10) {
			push (@AoA, [$s[0],$s[1],$s[2],$s[3], $s[4],$s[5],$s[6],$s[7],$s[8],$s[9]]);
		}
		else {
			push (@AoA, [$s[0],$s[1],$s[2],$s[3], $s[4],$s[5],$s[6],$s[7],$s[8],$s[9]]);
		}
		if ($s[3] =~ m/Type/i) {
			$typeCount{$s[3]}++;
		}
	}
	# my @sortedAoA = map  { $_->[0] }
					# sort { $a->[0] <=> $b->[0] }
					# map  { [ $_, $_->[0] ] }
					# @AoA;
	my @sortedAoA = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @AoA;
	#print Dumper(@sortedAoA)."\n";
	
	my @bedOut;
	for my $aref ( @sortedAoA ) {
		#print "\t [ @$aref ],\n";
		my $lineOut = "";
		if (defined @$aref[9]) {
			$lineOut = "@$aref[0]\t@$aref[1]\t@$aref[2]\t@$aref[3]\t@$aref[4]\t@$aref[5]\t@$aref[6]\t@$aref[7]\t@$aref[8]\t@$aref[9]";
		}
		else {
			$lineOut = "@$aref[0]\t@$aref[1]\t@$aref[2]\t@$aref[3]\t@$aref[4]\t@$aref[5]\t@$aref[6]\t@$aref[7]\t@$aref[8]";
		}
		push @bedOut, $lineOut;		
	}	
	
	# foreach my $line (@sortedAoA) {
		# foreach my $value (@{$line}) {
			# my $lineOut = "$value->[0]\t$value->[1]\t$value->[2]\n";
			# push @bedOut, $lineOut;
		# }
	# }
	return (\@bedOut,\%typeCount);
}

sub findBED {
	my $in = shift;
	opendir(DIR, $in);
	my @bed = grep(/\.bed$/,readdir(DIR));
	closedir(DIR);
	
	return $bed[0];
}

sub findBEDs {
	my $in = shift;
	opendir(DIR, $in);
	my @beds = grep(/\.bed$/,readdir(DIR));
	closedir(DIR);
	
	#@beds = map { $_->[1] }
	#		sort { $a->[0] <=> $b->[0] }
	#		map { [ ($_ =~ /(\d+)/)[0] || 0, $_ ] } @beds;
		
	#@beds = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } @beds;
	foreach (@beds) {
		$_ = $in."/".$_;
	}
	#@beds = sort { lc($a) cmp lc($b) } @beds;
	#@beds = sort { getnum($a) <=> getnum($b) } @beds;
	
	return @beds;
}

sub getnum {
	my $v = shift;
	return( ($v =~ /(\d+)/)[0] || 0);
}

sub findXMAPs {
	my $in = shift;
	opendir(DIR, $in);
	my @xmaps = grep(/\.xmap$/,readdir(DIR));
	closedir(DIR);
	
	return @xmaps;
}

sub findQCMAPs {
	my $in = shift;
	opendir(DIR, $in);
	my @qcmaps = grep(/fragileSiteRepaired_q.cmap$/,readdir(DIR));
	closedir(DIR);
	
	return @qcmaps;
}

sub wanted { 
	m/_temp/g and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; } 

sub wanted2 { 
	m/_merged_contigs_q/ and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; }
		
sub wanted3 { 
	m/_temp_/g and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; } 		
		
	
sub getCmapIds {
	my $cmapFile = shift;
	my @cmapIds;
	open CMAP, "<$cmapFile" or die "ERROR: Could not open $cmapFile: $!\n";
	while (<CMAP>) {
		my $line = $_;
		chomp($line);
		#if header then skip
		if ($line =~ /^#/) {
			next; }
		else {
			my @s = split(/\t/,$line);
			for (my $i=0; $i<scalar(@s); $i++) {
				$s[$i] =~ s/^\s+|\s+$//g; 
			}
			push @cmapIds, $s[0];
		}
	}
	@cmapIds = unique(@cmapIds);
	
	return @cmapIds;
}

sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
	
	
