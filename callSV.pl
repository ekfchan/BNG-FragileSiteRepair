#!/usr/bin/perl

# Script to merge call structural variations on FSITEREPAIRED_FINAL.cmap 
#
# Usage: callSV.pl --ref <reference.cmap> --cmap <input.cmap> --errbin <input.errbin> --optarg <optArguments.xml> --outdir <output_folder> [--gaps <N-base gaps BED file>] [--prefix <contig_prefix>] [--mres <pixels_to_reduce_ref> [--force] [--n <CPU cores to use>] [--stitchBED <BED file of stitch locations to filter final smap>]

use strict; 
use warnings; 
use Cwd; use Cwd qw(abs_path);
use File::Basename;
use File::Path qw(rmtree mkpath);
use File::Copy;
use Getopt::Long; 
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Sys::MemInfo qw(totalmem freemem totalswap);
use IPC::System::Simple qw(system capture);
use DateTime;
use DateTime::Format::Human::Duration;

my $dtStart = DateTime->now;

# Make sure dependency scripts exist
#my $splitCmapByIds = abs_path(dirname($0))."/splitCmapByIds.sh";
# print "Split script: $splitCmapByIds\n";
my $runSV = glob("~/scripts/runSV.py");
# print "Calc Stats script: $calc_cmap_stats\n";
if (!-e $runSV) {
	die "ERROR: Dependency scripts not found: $!\n"; 
}

# << compute requirement >> 
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
my $cpuCount = $cpu->count;
my $mem = int((((&totalmem / 1024) / 1024)) / 1024);


# << usage statement and variable initialisation >>
my %inputs = (); 
$inputs{'force'}=0;
GetOptions( \%inputs, 'ref=s', 'cmap=s', 'errbin=s', 'outdir=s', 'prefix:s', 'mres:f', 'optarg=s', 'force','n:i', 'gaps=s','stitchBED=s'); 

if ( (!exists $inputs{ref} & !exists $inputs{cmap}) | !exists $inputs{outdir} | !exists $inputs{errbin} | !exists $inputs{optarg}) {
	print "Usage: callSV.pl --ref <reference.cmap> --cmap <input.cmap> --errbin <input.errbin> --optarg <optArguments.xml> --outdir <output_folder> [--gaps <N-base gaps BED file>] [--prefix <contig_prefix>] [--mres <pixels_to_reduce_ref> [--force] [--n <CPU cores to use>] [--stitchBED <BED file of stitch locations to filter final smap>]\n"; 
	exit 0; 
} else {
	print "\nStart time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n";
}
foreach my $key ("ref","cmap","errbin","optarg") {
	$inputs{$key} = abs_path($inputs{$key});
	if ( ! -e $inputs{$key} ) { die "Input $inputs{$key} does not exists: $!\n"; }
}
if( exists $inputs{n} ) { $cpuCount = $inputs{n}; }
if( exists $inputs{gaps} ) { $inputs{gaps} = abs_path($inputs{gaps}); }
if( exists $inputs{stitchBED} ) { $inputs{stitchBED} = abs_path($inputs{stitchBED}); }

# create output folders
$inputs{outdir} = abs_path($inputs{outdir});
# if (-d $inputs{outdir} and $inputs{force} eq 0) {
# 	die "ERROR: Output directory $inputs{outdir} exists and force disabled. Use --force to overwrite.\n"
# }
# elsif (-d $inputs{outdir} and $inputs{force} eq 1) {
# 	print "WARNING: Output directory $inputs{outdir} exists and will be overwritten\n";
# 	rmtree($inputs{outdir}) or die "ERROR: Cannot remove existing output directory $inputs{outdir}: $!\n";
# }
# mkpath($inputs{outdir}) or die "ERROR: Cannot create output directory $inputs{outdir}: $!\n";

# my $svdir = $inputs{outdir}."/fsiterepaired_final_sv"; 
my $svdir = $inputs{outdir}; 
if (-d $svdir and $inputs{force} eq 0) {
	die "ERROR: Output directory $svdir exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $svdir and $inputs{force} eq 1) {
	print "WARNING: Output directory $svdir exists and will be overwritten\n";
	rmtree($svdir) or die "ERROR: Cannot remove existing output directory $svdir: $!\n";
}

# my $qcmapdir = $inputs{outdir}."/fsiterepaired_final"; 
my $qcmapdir = "$inputs{outdir}/contigs"; 
if (-d $qcmapdir and $inputs{force} eq 0) {
	die "ERROR: Output directory $qcmapdir exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $qcmapdir and $inputs{force} eq 1) {
	print "WARNING: Output directory $qcmapdir exists and will be overwritten\n";
	rmtree($qcmapdir) or die "ERROR: Cannot remove existing output directory $qcmapdir: $!\n";
}
# mkpath($qcmapdir) or die "ERROR: Cannot create output directory $qcmapdir: $!\n";
mkpath($qcmapdir) or die "ERROR: Cannot create output directory $qcmapdir: $!\n\n";

# create prefix 
if( !exists $inputs{prefix} ) { $inputs{prefix} = "FSITEREPAIRED_FINAL"; }


## << Split FSITEREPAIRED_FINAL.cmap into individual contigs >> 
# splitCmapByIds.sh <cmap to split> <output_folder> <output_prefix>
my $dtNow = DateTime->now;
print "\n\n===== Split $inputs{cmap} =====\n\n";

my $cmd = "cd $qcmapdir; ~/tools/RefAligner -i $inputs{cmap} -o $inputs{prefix} -split -stdout -stderr"; 
#my $cmd = "$splitCmapByIds $inputs{cmap} $qcmapdir $inputs{prefix}"; 
print "Running command: \n\n $cmd \n\n";
system($cmd); 

## Wait until qsub is done. 
#open( CMAPS, "$qcmapdir/cmapIdsToSplit.txt" ); 
#my $numcmaps=0; 
#while (<CMAPS>) { $numcmaps++ if  !/^\s+?$/;}
#close CMAPS;
#print "\nThere are $numcmaps CMAP IDs\n"; 
#my @numfiles;
#while( $numcmaps != scalar(@numfiles)-1 ) {
	##@numfiles = capture("ls -1 $qcmapdir/$inputs{prefix}"."*.cmap"); 
	#@numfiles = capture("ls -1 $qcmapdir/");
	#sleep(5);
	#my $left = ($numcmaps) - (scalar(@numfiles)-1);
	#print "\tWaiting...$left still to go...\n";
#}
#print "\nFinished. There are ",scalar(@numfiles)-1," cmaps in $qcmapdir\n"; 
#print "Waiting 10s for files to propogate...\n";
#sleep (9);

#copy($inputs{cmap},$qcmapdir."/") or warn "Copy failed: $!";

my $dtEnd = DateTime->now;
my $span = DateTime::Format::Human::Duration->new();
print 'Time spent on Split: ', $span->format_duration_between($dtEnd, $dtNow); print "\n\n";


## << De-res reference if --mres is provided  >> 
$dtNow = DateTime->now;
print "\n\n===== De-res $inputs{ref} =====\n\n";
my $ref; 
if( exists $inputs{mres} ) {
	if( $inputs{mres} =~ /^[+-]?\d+(\.\d+)?$/ ) {
		my $simplifiedres = $inputs{mres}; 
		$simplifiedres =~ s/[^\d]//g; 
		$ref = basename($inputs{ref}, ".cmap");
		$ref = $ref."_res".$simplifiedres; 
		my $cmd = "~/tools/RefAligner -i ".$inputs{ref}." -o $inputs{outdir}/".$ref." -f -merge -mres ".$inputs{mres}." -maxthreads ".$cpuCount." -maxmem ".$mem ; 
		print "Running command: \n\n $cmd \n\n";
		system($cmd); 
		$ref = $ref.".cmap"; 
		print "de-res cmap created: $ref \n"; 
	} else {
		print "Invalid --mres value: $inputs{mres}\n"; 
		exit 0; 
	}
} else {
	print "--mres was not provided. Assuming resolution of $inputs{ref} has already been reduced for runSV\n"; 
	$ref = $inputs{ref}; 
}

$dtEnd = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
print 'Time spent on Reducing Resolution of reference: ', $span->format_duration_between($dtEnd, $dtNow); print "\n\n";


## << Detect SV >> 
# runSV.py [-h] [-t REFALIGNER] [-r REFERENCEMAP] [-q QUERYDIR]
#                 [-o OUTPUTDIR] [-p PIPELINEDIR] [-a OPTARGUMENTS]
#                 [-T NUMTHREADS] [-j MAXTHREADS] [-b BEDFILE] [-e ERRFILE]
#                 [-E ERRBINFILE]
$dtNow = DateTime->now;
print "\n\n===== Detect SV =====\n\n";
#$cmd="python $runSV -t ~/tools/RefAligner -r $inputs{outdir}/$ref -q $qcmapdir -o $svdir -p ~/scripts/ -a $inputs{optarg} -T $cpuCount -j ".int($cpuCount/4)." -E $inputs{errbin}";
$cmd="python $runSV -t ~/tools/RefAligner -r $inputs{outdir}/$ref -q $qcmapdir -o $svdir -p ~/scripts/ -a $inputs{optarg} -T $cpuCount -E $inputs{errbin}";
if (defined $inputs{gaps}) {
	if (-e $inputs{gaps}) { $cmd = $cmd." -b $inputs{gaps}"; }
}
print "Running command: \n\n$cmd\n\n"; 
my @log = capture($cmd); 
my $logfile = $svdir."/runSV.log"; 
open LOG, ">$logfile" or die "ERROR: Cannot open $logfile for writing! $!\n";
print LOG join("\t", @log);
close LOG;
print "\nLog to runSV.py is in $logfile.\n\n";
print "\nSV summary is in $svdir/sv_log.txt\n\n"; 

$dtEnd = DateTime->now;
print "End time: "; print join ' ', $dtEnd->ymd, $dtEnd->hms; print "\n\n";
$span = DateTime::Format::Human::Duration->new();
print 'Time spent Calling SV: ', $span->format_duration_between($dtEnd, $dtNow); print "\n\n";


$dtEnd = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
print 'Total time elapsed: ', $span->format_duration_between($dtEnd, $dtStart); print "\n\n";

