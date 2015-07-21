#!/usr/bin/perl

# A wrapper script to run fragile site repair on assembly alignments to a reference using an input FASTA and alignref_final XMAP and CMAPs
# we will first calculate potential fragiles sites based on an input fasta, then stitch assembly maps together based on alignments if the maps start and stop overlapping a fragile site

# Usage: perl fragileSiteRepair.pl --fasta <reference.fasta> --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output folder> --maxlab <max_label_gap_tolerence=0> --maxfill <max basepairs to fill between contigs = 35000> --wobble <fragile site wobble in bp = 0> --force <overwrite output folder>

use strict;
use warnings;
use Getopt::Long; 
use Cwd qw(abs_path);
use File::Basename;
use File::Path qw(rmtree mkpath);
use IPC::System::Simple qw(system capture);
#use Parallel::Iterator qw(iterate_as_array);
use Parallel::Loops;
use DateTime;
use DateTime::Format::Human::Duration;

print "\n";
my $dtStart = DateTime->now;
print "Start time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n";

#get number of CPU cores
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
#printf "There are %d CPUs\n"  , $cpu->count || 1;
my $cpuCount = $cpu->count;
#print "CPU count: $cpuCount\n";

print "\n";

# << usage statement and variable initialisation >>
my %inputs = (); 
my $prefix;
$inputs{'force'}=0;
GetOptions( \%inputs, 'fasta=s', 'xmap=s', 'qcmap=s', 'rcmap=s', 'errbin=s', 'output=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'force'); 

if ( !exists $inputs{fasta} | !exists $inputs{xmap} | !exists $inputs{qcmap} | !exists $inputs{rcmap} | !exists $inputs{errbin} |!exists $inputs{output} | !exists $inputs{maxlab} | !exists $inputs{maxfill} | !exists $inputs{wobble} ) {
	print "Usage: perl fragileSiteRepair.pl --fasta <reference.fasta> --xmap <input.xmap> --qcmap <input_q.cmap> --rcmap <input_r.cmap> --errbin <input.errbin> --output <output folder> --maxlab <max_label_gap_tolerence> --maxfill <max basepairs to fill between contigs> --wobble <fragile site wobble in bp> --force <overwrite output folder>\n"; 
	exit 0; 
}
else {
	$prefix = $inputs{xmap}; $prefix =~ s/.xmap//i;
	foreach my $key (keys %inputs) {
		#print "Key: $key Value: $inputs{$key}\n";
		if ($inputs{$key} =~ /[a-z]/i) {
			#print "Key: $key Value: $inputs{$key}\n";
			$inputs{$key} = abs_path($inputs{$key});
			#print "$inputs{$key}\n";
		}
	}
}

#print out input variables
print "Input FASTA: $inputs{fasta}\n";
print "Input XMAP: $inputs{xmap}\n";
print "Input QCMAP: $inputs{qcmap}\n";
print "Input RMCAP: $inputs{rcmap}\n";
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
my $cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $prefix"."_contig $inputs{output}/contigs";
print "Running command: $cmd\n";
system($cmd);
print "\n";


# Step 2: Calculate fragile sites for input FASTA
print "===Step 2: Calculate fragile sites for input FASTA===\n";
# Usage: perl calcFragilesSites.pl <input FASTA> [output.bed]
$cmd = "perl $scriptspath/calcFragileSites.pl $inputs{fasta} $inputs{output}";
print "Running command: $cmd\n";
system($cmd);
print "\n";
my $bed = findBED($inputs{output});
$bed = abs_path($inputs{output}."/".$bed);
#print "BED file: $bed\n";

# Step 3: Run gapFill.pl for each anchor map
print "===Step 3: Run gapFill.pl for each anchor map===\n";
#print "\n";
my @xmaps = findXMAPs($inputs{output}."/contigs");
@xmaps = sort @xmaps;

foreach my $xmap (@xmaps) {
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
		
		# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
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
	}	
}

#my $pl = Parallel::Loops->new($cpuCount);
#my %out;
#$pl->share(\%out);

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
print "===Step 4: Merge individual fragileSiteRepaired anchor maps===\n";
my @qcmaps = findQCMAPs($inputs{output}."/contigs");
@qcmaps = sort @qcmaps;
print "Merging ".scalar(@qcmaps)." fragileSiteRepaired anchor maps...\n";
foreach (@qcmaps) {
	#print "QCMAP: $_\n";
	$_ = abs_path($inputs{output}."/contigs/$_");
}
my $input = join(" -i ",@qcmaps);
$cmd = "cd $inputs{output}; ~/tools/RefAligner -i $input -merge -o $prefix"."_fragileSiteRepaired -minsites 0";
print "Running command: $cmd\n";
print "\n";
system($cmd);
print "\n";


# Step 5: Calculate stats for merged fragileSiteRepaired CMAP
print "===Step 5: Calculate stats for merged fragileSiteRepaired CMAP===\n";
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
	
	my $file = "$inputs{output}/$prefix"."_fragileSiteRepaired.cmap";
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







sub findBED {
	my $in = shift;
	opendir(DIR, $in);
	my @bed = grep(/\.bed$/,readdir(DIR));
	closedir(DIR);
	
	return $bed[0];
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

	
