#!/usr/bin/perl

# A wrapper script to run fragile site repair on assembly alignments to a reference using an input FASTA and alignref_final XMAP and CMAPs
# If .bed file of fragile sites is not provided, an input (reference) fasta is expected, from which potential fragiles sites will be calculated. 
# Assembly maps are stitched together based on alignments if the maps start and stop overlap a fragile site

# Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bed>] [--cmap <assembly.cmap>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> [--ref <reference.cmap>] --output <output folder> --maxlab <max_label_gap_tolerence=0> --maxfill <max basepairs to fill between contigs = 35000> --wobble <fragile site wobble in bp = 0> --seq <sequence bp +/- fragile site to print in final BED> --force <overwrite output folder>

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

print "\n";
print qx/ps -o args $$/;
print "\n";

my $hostname = `hostname`;
chomp($hostname);
print "Running on host: $hostname\n";

#get number of CPU cores
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
#printf "There are %d CPUs\n"  , $cpu->count || 1;
my $cpuCount = $cpu->count;
#print "CPU count: $cpuCount\n";
# get memory info
use Sys::MemInfo qw(totalmem freemem totalswap);
my $mem = (((&totalmem / 1024) / 1024)) / 1024; 

print "Maximum CPU cores to use: $cpuCount\n";
print "Maximum memory to use: $mem GB\n";
print "\n";

# << usage statement and variable initialisation >>
my %inputs = (); 
my $prefix;
$inputs{'force'}=0;
GetOptions( \%inputs, 'fasta:s', 'xmap=s', 'qcmap=s', 'rcmap=s', 'errbin=s', 'output=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'force', 'bed:s', 'seq:i', 'cmap:s', 'ref:s'); 

my $hasbed = 0; 
if ( (!exists $inputs{fasta} & !exists $inputs{bed}) | !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{errbin} | !exists $inputs{output} ) {
	print "Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bed>] [--cmap <assembly.cmap>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> [--ref <reference.cmap>] --output <output folder> --maxlab <max_label_gap_tolerence=0> --maxfill <max basepairs to fill between contigs = 35000> --wobble <fragile site wobble in bp = 0> --seq <sequence bp +/- fragile site to print in final BED> --force <overwrite output folder>\n"; 
	exit 0; 
}
else {
	print "\nStart time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";
	
	foreach my $key ("xmap","qcmap","rcmap","errbin","output","bed","fasta", "cmap", "ref") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
	$prefix = basename(abs_path($inputs{xmap}), ".xmap");
	if( !exists $inputs{maxfill} ) { $inputs{maxfill} = 10000; }
	if( !exists $inputs{wobble} ) { $inputs{wobble} = 250; }
	if( !exists $inputs{maxlab} ) { $inputs{maxlab} = 0; }
	if( !exists $inputs{seq} ) { $inputs{seq} = 50; }

	if( exists $inputs{bed} ) {
		$hasbed = 1; 
		if( exists $inputs{fasta} ) {
			print "$inputs{bed} fragile sites file already provided.\n$inputs{fasta} will be ignored.\n";
		}
	} 
}

#print out input variables
if( $hasbed==0 ) { print "Input FASTA: $inputs{fasta}\n"; } else { print "Input BED: $inputs{bed}\n"; }
print "Input XMAP: $inputs{xmap}\n";
print "Input QCMAP: $inputs{qcmap}\n";
print "Input RMCAP: $inputs{rcmap}\n";
print "Input ERRBIN: $inputs{errbin}\n";
if (exists $inputs{cmap}) {print "Input CMAP: $inputs{cmap}\n";}
if (exists $inputs{ref}) {print "Input reference CMAP: $inputs{ref}\n";}
print "Output folder: $inputs{output}\n";
print "Maximum labels between maps: $inputs{maxlab}\n";
print "Maximum basepairs to fill between maps: $inputs{maxfill}\n";
print "Maximum fragile site wobble: $inputs{wobble}\n";
print "Sequence bp +/- start/end of fragile to print to BED: $inputs{seq}\n";
print "\n";


# Make sure dependency scripts exist
my $scriptspath = abs_path(dirname($0));
if (!-e "$scriptspath/calcFragileSites.pl" | !-e "$scriptspath/split_xmap_standalone.pl" | !-e "$scriptspath/gapFill.pl") {
	die "ERROR: Dependency scripts not found at $scriptspath\n"; 
}

# check output folder
if (-d $inputs{output} and $inputs{force} eq 0) {
	die "ERROR: Output directory $inputs{output} exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $inputs{output} and $inputs{force} eq 1) {
	print "WARNING: Output directory $inputs{output} exists and will be overwritten\n";
	rmtree($inputs{output}) or die "ERROR: Cannot remove existing output directory $inputs{output}: $!\n";
	mkpath($inputs{output}) or die "ERROR: Cannot create output directory $inputs{output}: $!\n";
}
else {
	mkpath($inputs{output}) or die "ERROR: Cannot create output directory $inputs{output}: $!\n";
}
print "\n";



# Step 1: Split input XMAP into individual anchor maps
print "===Step 1: Split merged input XMAP into individual anchor maps===\n";
# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
my $splitprefix = basename($prefix);
my $cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $splitprefix"."_contig $inputs{output}/contigs";
print "Running command: $cmd\n";
system($cmd);
print "\n";


# Step 2: Calculate fragile sites for input FASTA
print "===Step 2: Calculate fragile sites for input FASTA===\n";
my $bed; 
if ( $hasbed==1 and -e $inputs{bed}) {
	print "Fragiles sites BED file $inputs{bed} provided.\nSkipping fragile sites calculation...\n\n";
	copy("$inputs{bed}","$inputs{output}/") or warn "Copy of input BED $inputs{bed} failed: $!";
	$bed = $inputs{bed};
}
else {
	my $stime = DateTime->now;
	# Usage: perl calcFragilesSites.pl <input FASTA> [output.bed] [sequence bp +/- fsite to print out]
	$cmd = "perl $scriptspath/calcFragileSites.pl $inputs{fasta} $inputs{output} $inputs{seq}";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	# my $bed = findBED($inputs{output});
	$bed = findBED($inputs{output});
	$bed = abs_path($inputs{output}."/".$bed);
	#print "BED file: $bed\n";
	
	my $etime = DateTime->now;
	my $dtime = DateTime::Format::Human::Duration->new();
	print 'Spent time at this step: ', $dtime->format_duration_between($etime, $stime); print "\n\n";
}

# Step 3: Run gapFill.pl for each anchor map
print "===Step 3: Run gapFill.pl for each anchor map===\n";
print "Processing anchor maps in parallel fashion...\n";
print "\n";
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
my $pm =  new Parallel::ForkManager(8);
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
		
		# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -e <errbin> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -e $inputs{errbin} -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble}";
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


# Step 4: Merge individual fragileSiteRepaired anchor maps
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
$cmd = "cd $inputs{output}; ~/tools/RefAligner -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0";
print "Running command: $cmd\n";
print "\n";
system($cmd);
print "\n";

print "Generating merged stitchPositions BED...";
my $mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.withSeqs.bed";
open BEDOUT, ">$mergedBED" or die "ERROR: Could not open $mergedBED: $!\n";
my @beds = findBEDs("$inputs{output}/contigs");
my @bedOut;
#print BEDOUT "#CMapId\tStart\tEnd\tType\tSequence\n";
print BEDOUT "#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence\n";
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
}
my ($bedOut_ref, $typeCount_ref) = sortBEDandCount(\@bedOut);
@bedOut = @$bedOut_ref;
my %typeCount = %$typeCount_ref;
my $stitchCount = scalar(@bedOut);
print BEDOUT join("\n",@bedOut);
print "done\n\n";

#generate FASTA of stitchPositions sequences
$cmd = "perl $scriptspath/extractFASTA_fromBED.pl --bed $mergedBED";
#print "Running command $cmd\n";
print "Generating FASTA from merged stitchPositions BED...\n";
system($cmd);
print "\n";

# Step 5: Calculate stats for merged fragileSiteRepaired CMAP
print "===Step 5: Calculate stats for merged fragileSiteRepaired CMAP===\n";
print "\n"; 
#get stats of original input BED
my @origBed;
open ORIGBEDFILE, "<$bed" or die "ERROR: Could not open $bed: $!\n";
while (<ORIGBEDFILE>) {
		if ($_ =~ /^#/) {
			next;
		}
		chomp($_);
		#print BEDOUT "$_\n";
		push @origBed, $_;
}
my ($origBedOut_ref, $origTypeCount_ref) = sortBEDandCount(\@origBed);
my %origTypeCount = %$origTypeCount_ref;
my $origCount = scalar(@origBed);
print "Total potential fragile sites: $origCount\n";
foreach my $type (sort keys %origTypeCount) {
	print "\t$type: $origTypeCount{$type}\n";
}
print "Total fragile sites repaired (stitched): $stitchCount\n";
foreach my $type (sort keys %typeCount) {
	print "\t$type: $typeCount{$type}\n";
}

print "\n";

# usage: calc_cmap_stats.pl <CMAP_File>
my $dir = glob("~/scripts/HybridScaffold/scripts");
my $script = "calc_cmap_stats.pl";
if (-e "$dir/$script") {
	print "Original query CMAP $inputs{qcmap} stats:\n";
	#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
	my $cmd = "perl $dir/$script $inputs{qcmap}";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	
	my $file = "$inputs{output}/$splitprefix"."_fragileSiteRepaired.cmap";
	print "Fragile site repaired merged CMAP $file stats:\n";
	#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
	$cmd = "perl $dir/$script $file";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	
}
else {
	print "perl script calc_cmap_stats.pl not found at $script\n"; }
print "\n";

# Step 6: Generate new merged cmap with unmapped maps
if (exists $inputs{cmap}) {
	print "===Step 6: Process original assembly CMAP and generate new merged CMAP ===\n";
	print "\n"; 
	my @origCmapIds = getCmapIds($inputs{cmap});
	my @qryCmapIds = getCmapIds($inputs{qcmap});
	my %in_qryCmapIds = map {$_ => 1} @qryCmapIds;
	my @diff  = grep {not $in_qryCmapIds{$_}} @origCmapIds;
	my $oldMap = abs_path($inputs{output}."/".$splitprefix."_fragileSiteRepaired.cmap");
	my $tempMap = $inputs{output}."/temp";
	my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged";
	#print "\tcmapIds in orig but not ref: ".join("\n\t",@diff)."\n";
	$cmd = "cd $inputs{output}; ~/tools/RefAligner -merge -i $inputs{cmap} -selectid ".join(" ",@diff)." -o $tempMap";
	print "\tRunning command: $cmd\n\n";
	system($cmd);
	print "\n";
	$cmd = "cd $inputs{output}; ~/tools/RefAligner -merge -i $tempMap".".cmap -i $oldMap -o $newMap; rm -f $tempMap".".cmap";
	print "\tRunning command: $cmd\n\n";
	system($cmd);
	print "\n";
	$newMap = $newMap.".cmap";
}
else {
	print "Skipping Step 6. Original assembly CMAP not provided\n";
	print "\n";
}

# Step 7: run CharacterizeFinal on newly merged cmap
if ((exists $inputs{cmap}) && (exists $inputs{ref})) {
	my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged.cmap";
	print "===Step 7: Align new merged CMAP to reference and get stats ===\n";
	print "\n"; 
	#usage: runCharacterizeFinal.py [-h] [-t REFALIGNER] [-r REFERENCEMAP] [-q QUERYMAP] [-x XMAP] [-p PIPELINEDIR] [-a OPTARGUMENTS] [-n NUMTHREADS] 
	$cmd = "cp runCharacterizeFinal.py ~/scripts/; python ~/scripts/runCharacterizeFinal.py -t ~/tools/RefAligner -r $inputs{ref} -q $newMap -p ~/scripts/ -n $cpuCount";
	print "\tRunning command: $cmd\n\n";
	system($cmd);
}
else {
	print "Skipping Step 7. Original assembly CMAP or reference CMAP not provided\n";
	print "\n";
}















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
		push (@AoA, [$s[0],$s[1],$s[2],$s[3], $s[4],$s[5],$s[6],$s[7],$s[8],$s[9]]);
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
		my $lineOut = "@$aref[0]\t@$aref[1]\t@$aref[2]\t@$aref[3]\t@$aref[4]\t@$aref[5]\t@$aref[6]\t@$aref[7]\t@$aref[8]\t@$aref[9]";
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
	m/_temp/ and do { 
		unlink $_ or warn "Could not unlink file $_\n"; }; } 

sub wanted2 { 
	m/_merged_contigs_q/ and do { 
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
	
	
