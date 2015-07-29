#!/usr/bin/perl

# Script to merge call structural variations on FSITEREPAIRED_FINAL.cmap 
#
# Usage: callSV.pl --ref <reference.cmap> --cmap <input.cmap> --errbin <input.errbin> --optarg <optArguments.xml> --outdir <output_folder> [--prefix <contig_prefix>] [--mres <pixels_to_reduce_ref> [--force]

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
my $splitCmapByIds = abs_path(dirname($0))."/splitCmapByIds.sh";
# print "Split script: $splitCmapByIds\n";
my $runSV = glob("~/scripts/runSV.py");
# print "Calc Stats script: $calc_cmap_stats\n";
if (!-e $splitCmapByIds | !-e $runSV) {
	die "ERROR: Dependency scripts not found: $!\n"; 
}

# << compute requirement >> 
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
my $cpuCount = $cpu->count;
my $mem = (((&totalmem / 1024) / 1024)) / 1024; 


# << usage statement and variable initialisation >>
my %inputs = (); 
$inputs{'force'}=0;
GetOptions( \%inputs, 'ref=s', 'cmap=s', 'errbin=s', 'outdir=s', 'prefix:s', 'mres:f', 'optarg=s', 'force'); 

if ( (!exists $inputs{ref} & !exists $inputs{cmap}) | !exists $inputs{outdir} | !exists $inputs{errbin} | !exists $inputs{optarg}) {
	print "Usage: callSV.pl --ref <reference.cmap> --cmap <input.cmap> --errbin <input.errbin> --optarg <optArguments.xml> --outdir <output_folder> [--prefix <contig_prefix>] [--mres <pixels_to_reduce_ref> [--force]\n"; 
	exit 0; 
} else {
	print "\nStart time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n";
}
foreach my $key ("ref","cmap","errbin","optarg") {
	$inputs{$key} = abs_path($inputs{$key});
	if ( ! -e $inputs{$key} ) { die "Input $inputs{$key} does not exists: $!\n"; }
}

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

# my $qcmapdir = $inputs{outdir}."/fsiterepaired_final"; 
my $qcmapdir = $inputs{outdir}; 
if (-d $qcmapdir and $inputs{force} eq 0) {
	die "ERROR: Output directory $qcmapdir exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $qcmapdir and $inputs{force} eq 1) {
	print "WARNING: Output directory $qcmapdir exists and will be overwritten\n";
	rmtree($qcmapdir) or die "ERROR: Cannot remove existing output directory $qcmapdir: $!\n";
}
# mkpath($qcmapdir) or die "ERROR: Cannot create output directory $qcmapdir: $!\n";

# my $svdir = $inputs{outdir}."/fsiterepaired_final_sv"; 
my $svdir = $inputs{outdir}."_sv"; 
if (-d $svdir and $inputs{force} eq 0) {
	die "ERROR: Output directory $svdir exists and force disabled. Use --force to overwrite.\n"
}
elsif (-d $svdir and $inputs{force} eq 1) {
	print "WARNING: Output directory $svdir exists and will be overwritten\n";
	rmtree($svdir) or die "ERROR: Cannot remove existing output directory $svdir: $!\n";
}
mkpath($svdir) or die "ERROR: Cannot create output directory $svdir: $!\n\n";

# create prefix 
if( !exists $inputs{prefix} ) { $inputs{prefix} = "FSITEREPAIRED_FINAL"; }


## << Split FSITEREPAIRED_FINAL.cmap into individual contigs >> 
# splitCmapByIds.sh <cmap to split> <output_folder> <output_prefix>
my $dtNow = DateTime->now;
print "\n\n===== Split $inputs{cmap} =====\n\n";
my $cmd = "$splitCmapByIds $inputs{cmap} $qcmapdir FSITEREPAIRED_FINAL"; 
print "Running command: \n\n $cmd \n\n";
system($cmd); 

# Wait until qsub is done. 
open( CMAPS, "$qcmapdir/cmapIdsToSplit.txt" ); 
my $numcmaps=0; 
while (<CMAPS>) { $numcmaps++ if  !/^\s+?$/;}
close CMAPS;
print "\nThere are $numcmaps CMAP IDs\n"; 
my @numfiles;
while( $numcmaps != scalar(@numfiles) ) {
	@numfiles = capture("ls -1 $qcmapdir/FSITEREPAIRED_FINAL*.cmap"); 
	sleep(1);
	print "waiting...\n";
}
print "\nThere are ",scalar(@numfiles)," cmaps in $qcmapdir\n"; 

copy($inputs{cmap},$qcmapdir."/FSITEREPAIRED_FINAL.cmap") or warn "Copy failed: $!";

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
		$ref = $inputs{ref}."_res".$simplifiedres; 
		my $cmd = "~/tools/RefAligner -i ".$inputs{ref}." -o ".$ref." -f -merge -mres ".$inputs{mres}." -stdout -stderr -maxthreads ".$cpuCount." -maxmem ".$mem ; 
		print "Running command: \n\n $cmd \n\n";
		system($cmd); 
		$ref = $ref.".cmap"; 
		print "de-res cmap created: $ref \n"; 
	} else {
		print "Invalid --res value: $inputs{mres}\n"; 
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
$cmd="python $runSV -t ~/tools/RefAligner -r $ref -q $qcmapdir -o $svdir -p ~/scripts/ -a $inputs{optarg} -T $cpuCount -E $inputs{errbin}";
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

