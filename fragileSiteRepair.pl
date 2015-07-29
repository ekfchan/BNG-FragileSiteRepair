#!/usr/bin/perl

# A wrapper script to run fragile site repair on assembly alignments to a reference using an input FASTA and alignref_final XMAP and CMAPs
# If .bed file of fragile sites is not provided, an input (reference) fasta is expected, from which potential fragiles sites will be calculated. 
# Assembly maps are stitched together based on alignments if the maps start and stop overlap a fragile site

# Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bed>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output folder> --maxlab <max_label_gap_tolerence=0> --maxfill <max basepairs to fill between contigs = 35000> --wobble <fragile site wobble in bp = 0> --force <overwrite output folder>

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

#get number of CPU cores
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
#printf "There are %d CPUs\n"  , $cpu->count || 1;
my $cpuCount = $cpu->count;
#print "CPU count: $cpuCount\n";


# << usage statement and variable initialisation >>
my %inputs = (); 
my $prefix;
$inputs{'force'}=0;
GetOptions( \%inputs, 'fasta:s', 'xmap=s', 'qcmap=s', 'rcmap=s', 'errbin=s', 'output=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'force', 'bed:s'); 

my $hasbed = 0; 
if ( (!exists $inputs{fasta} & !exists $inputs{bed}) | !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{errbin} | !exists $inputs{output} ) {
	print "Usage: perl fragileSiteRepair.pl [--fasta <reference.fasta>] [--bed <reference_fsites.bd>] --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output folder> --maxlab <max_label_gap_tolerence> --maxfill <max basepairs to fill between contigs> --wobble <fragile site wobble in bp> --force <overwrite output folder>\n"; 
	exit 0; 
}
else {
	print "\n\nStart time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";
	
	foreach my $key ("xmap","qcmap","rcmap","errbin","output","bed","fasta") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
	$prefix = basename(abs_path($inputs{xmap}), ".xmap");
	if( !exists $inputs{maxfill} ) { $inputs{maxfill} = 10000; }
	if( !exists $inputs{wobble} ) { $inputs{wobble} = 250; }
	if( !exists $inputs{maxlab} ) { $inputs{maxlab} = 0; }

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
print "Output folder: $inputs{output}\n";
print "Maximum labels between maps: $inputs{maxlab}\n";
print "Maximum basepairs to fill between maps: $inputs{maxfill}\n";
print "Maximum fragile site wobble: $inputs{wobble}\n";
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
	# Usage: perl calcFragilesSites.pl <input FASTA> [output.bed]
	$cmd = "perl $scriptspath/calcFragileSites.pl $inputs{fasta} $inputs{output}";
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
my $mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.bed";
open BEDOUT, ">$mergedBED" or die "ERROR: Could not open $mergedBED: $!\n";
my @beds = findBEDs("$inputs{output}/contigs");
my @bedOut;
print BEDOUT "#CMapId\tStart\tEnd\tType\n";
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

# Step 5: Calculate stats for merged fragileSiteRepaired CMAP
print "===Step 5: Calculate stats for merged fragileSiteRepaired CMAP===\n";
print "\n"; 
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
		push (@AoA, [$s[0],$s[1],$s[2],$s[3]]);
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
		my $lineOut = "@$aref[0]\t@$aref[1]\t@$aref[2]\t@$aref[3]";
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
	

