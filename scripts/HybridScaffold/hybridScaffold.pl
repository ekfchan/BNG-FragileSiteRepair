# $Id: hybridScaffold.pl 4304 2015-11-24 22:43:43Z apang $

#!/usr/bin/perl -w
#################################################################################################################################################################################
# File: hybridScaffold.pl                                                                  					
# Date: 07/24/2014                                                                         				
# Purpose: Merge NGS scaffolds with BioNano CMAP                                           				
#                                                                                          			
#                                                                                          		
# Usage:                                                                                   	
#  hybridScaffold.pl <-h> <-n ngs_file> <-b bng_cmap_file> <-c hybrid_config_xml> <-o output_folder> <-B conflict_filter_level> <-N conflict_filter_level> <-f> <-m molecules_bnx> <-p de_novo_pipeline> <-q de_novo_xml> <-v>	
#     -h    : This help message                                                            												
#     -n    : Input NGS FASTA or CMAP file [required]                                      												
#     -b    : Input BioNano CMAP [required]                                                												
#     -c    : Merge configuration file [required]                                          												
#     -o    : Output folder [required]                                                     												 
#     -B    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contigs [required if not using -M option]
#     -N    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contigs [required if not using -M option]
#     -r    : RefAligner program [required]
#     -f    : Force output and overwrite any existing files                                												
#     -m    : Input BioNano molecules BNX for molecules to hybrid scaffold alignment [optional]    										
#     -p    : Input de novo assembly pipeline script [required for -m option]
#     -q    : Input de novo assembly pipeline optArguments XML script [required for -m option]
#     -v    : Print pipeline version information																
#     -M    : Input a conflict resolution file indicating which NGS and BioNano conflicting contigs to be cut [optional]
#
#################################################################################################################################################################################

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/scripts/perl5/" direcory at run time to the @INC array
# This script sould sit two levels above the additional Perl modules directory
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path = abs_path($script_path . "/scripts/perl5");
	unshift @INC, $module_path;
	my $lib3;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path = $module_path."/5.10.1"; 
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path = $module_path."/5.14.4"; 
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.016000 && $] <= 5.017000) {
        	$module_path = $module_path."/5.16.3";
        	$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.018000 && $] <= 5.019000) {
		$module_path = $module_path."/5.18.2";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n"; 
		exit; }
	unshift @INC, $module_path;
	unshift @INC, $lib3;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use IO::Handle;
use PerlIO::Util;
use Getopt::Std;
use Config::Simple;
use File::Path qw(mkpath rmtree);
use File::Slurp;
use File::Copy qw(copy move);
use File::Copy::Recursive qw(fcopy);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use XML::Simple;
use DateTime;
use DateTime::Format::Human::Duration;
use File::Basename;
use threads;
use IO::Select;
use IPC::Open3;
use File::Spec;
use File::Find;
sub Init;
sub Usage;
sub CHECK_Config;
sub split_xmap;
sub Version;
sub find_refaligner;
sub log10;
use List::MoreUtils qw(uniq);
use BNG::Utility;

my %opts;
my @cmd = ();	
my $cmdRef = \@cmd;
my ($outResults, $errResults) = ("", "");
my $plDir=abs_path(dirname($0));

#read comand line args
Init();

# this file indicates where the final results should be

my $hybridScaffoldFinalResultsDir = "hybrid_scaffolds";	
my $modifyNum = getCurMaxModifyNum($opts{o});	# to be used to distinctly identify the number of times a modified conflict status file was fed (ignored if no such file was inputted)
$hybridScaffoldFinalResultsDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
$hybridScaffoldFinalResultsDir = abs_path("$opts{o}/$hybridScaffoldFinalResultsDir");
eval {mkpath($hybridScaffoldFinalResultsDir)};
if ($@)	{
	print("ERROR: Could not create final output directory $hybridScaffoldFinalResultsDir: $@");
	print("ERROR: Output folder invalid or does not exist!\n");
	Usage();
}
open(OUT_DIR_POINTER, ">$opts{o}/cur_results.txt") or dieLog "ERROR: Cannot write a file to indicate which hybrid scaffold final results dir is: $!\n";
print OUT_DIR_POINTER "$hybridScaffoldFinalResultsDir\n";
if (defined($opts{m}))	{
	print OUT_DIR_POINTER "$opts{o}/alignmol_bionano";	print OUT_DIR_POINTER "_M$modifyNum" if (defined($opts{M}));	print OUT_DIR_POINTER "\n";
	print OUT_DIR_POINTER "$opts{o}/alignmol_hybrid";	print OUT_DIR_POINTER "_M$modifyNum" if (defined($opts{M}));	print OUT_DIR_POINTER "\n";
}
close OUT_DIR_POINTER;

my @s = split("/", $opts{o});

my @NGSpath = split(/\//, $opts{n});
my $NGSfile = $NGSpath[$#NGSpath];
$NGSfile =~ s/\./_/g;	# replaces all . to underscores

my @BNGpath = split(/\//, $opts{b});
my $BNGfile = $BNGpath[$#BNGpath];	# get the actual file name
$BNGfile =~ s/\./_/g;	# replaces all . to underscores

my $temp = $BNGfile;	$temp =~ s/_cmap$//;	# rid the file extension
my $modBNGfile = $temp."_bppAdjust_cmap";
my $log_file2 = "$hybridScaffoldFinalResultsDir/$modBNGfile"."_$NGSfile"."_HYBRID_SCAFFOLD_log.txt";
for (*STDOUT, *STDERR)	{
	# captures the stdout and stderr 
	$_->autoflush;	$_->push_layer(tee=>$log_file2)
} # for	
	
#print header
print "**************************\n";
print "*****BioNano Genomics*****\n";
print "******BNG-NGS Merge*******\n";
print "**************************\n\n";

Version();
print "\n";

my $hostname = `hostname`;
print "Running on host: $hostname\n";

my $dtStart = DateTime->now;
print "Start time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";

print qx/ps -o args $$/;

#make sure all input files exist 
print "\n";
if (-e $opts{n} && -e $opts{b} && -e $opts{c}) {
	
	print "NGS file: $opts{n} \n" ;
	print "BNG file: $opts{b} \n";
	print "Configuration file: $opts{c} \n"; 
	print "Output folder: $opts{o}\n"; 
	print "Molecules file: $opts{m}\n" if (defined($opts{m}));
	print "User-defined conflict status file: $opts{M}\n" if (defined($opts{M}));
} else {
	print("One or more input files do not exist. \n");
	Usage(); 
}

#make sure input BNG CMAP file contains at least one contig
my (undef, $numContig, undef) = readCMap($opts{b});
if (!($numContig > 0)) {
	dieLog( "ERROR: Input BNG file $opts{b} does not contain any contigs!\n\n");
	exit;
}

### RefAligner section
my $refaligner = $opts{r};
if ($refaligner =~ /~/)	{
	$refaligner = glob($refaligner);
} else	{
	$refaligner = abs_path($refaligner) if (-e $refaligner);
} # if refaligner

# make sure RefAligner exists
if (! -e $refaligner)	{
	print "WARNING: RefAligner binary does not exist at $refaligner as defined in config file. Trying to find it...\n";
	my $script_path = abs_path(dirname($0));
	my @s = split('/',$script_path);
	my $val = pop(@s); $val = pop(@s);
	my $home_path = join('/',@s);
	$refaligner = $home_path."/tools/RefAligner";
	
	if (! -e $refaligner)	{
		$refaligner = glob("~/tools/RefAligner");
		if (! -e $refaligner)	{
			dieLog ("ERROR: RefAligner binary cannot be found!\n");
		} # if -e refaligner (final try)
	} # if -e refaligner (second try)
} # if -e refaligner (first try)
print "RefAligner binary: $refaligner\n";
### end of RefAligner section

### config file section
#load config file
my $XML = new XML::Simple(KeyAttr=>[]);
my $config = $XML->XMLin($opts{c});

# parse config XML file
my %config = Parse($config);

#check config file
CHECK_Config(\%config);

copy "$opts{c}", "$hybridScaffoldFinalResultsDir" ;
### end of config file section
	
my $maxmem=$config{'global.maxmem'};
my $maxthreads=$config{'global.maxthreads'};

# define some variables that will be needed regardless of whether a conflict status file is fed in

my @enzymes = split(/\s+/, $config{'fasta2cmap.enzyme'});	
for (my $i = 0; $i < scalar(@enzymes); $i += 1)	{
	$enzymes[$i] = uc($enzymes[$i]);	# capitalize the entire enzyme name
} # foreach enzyme
my $channelNumber = $config{'fasta2cmap.channelNum'};
my $minLabels = $config{'fasta2cmap.minLabels'};
my $minLength = $config{'fasta2cmap.minLength'};
my @file = split(/\./, $opts{n});	my $filename = $opts{n};	$filename =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"};
my $filename_key = $filename."_".$minLength."kb_".$minLabels."labels_key.txt";
my $ngs_cmap = $filename."_".$minLength."kb_".$minLabels."labels.cmap";
my $bppAdjust = "";					# adjust the bpp of BioNano genome maps after align to sequence in align0
my $readparameters = "";				# to store the noise parameters in align1
my ($assignAlign_r, $assignAlign_q) = ("", "");		# to store those contigs that passed the contig removal process
my ($cutBNFile, $cutSeqFile) = ("", "");
my $outputDir = "";
my ($dtStartStage, $dtEndStage) = (DateTime->now, DateTime->now);
my $span = DateTime::Format::Human::Duration->new();

if (! defined($opts{M}))	{
	# if the user did not incorporate a conflict status file manually, then the pipeline starts normally from the beginning

	# first check if chimeric quality scores are needed
	if ($opts{B} == 2 || $opts{N} == 2)	{
		# if cut was indicated to be done, then the presence of chimeric quality score is checked, else reset the option to remove conflict
		my %qScores = ();	my $qScoresRef = \%qScores;	my $noQScoreFlag = 0;
		($qScoresRef, $noQScoreFlag) = getQScores($opts{b}, $qScoresRef);
		if ($noQScoreFlag == 1)	{
			# print a warning; switch opts{N} to 3 if opts{N} was 2; switch opts{B} to 3 if opts{B} was 2
			print "WARNING: cut at conflict option was selected, but no chimeric quality score is available in the genome map\n";
			if ($opts{B} == 2)	{
				$opts{B} = 3;
				print "WARNING: when conflict occurs, BioNano map is thrown out instead of cut\n";
			} # if opts{B}
			if ($opts{N} == 2)	{
				$opts{N} = 3;
				print "WARNING: when conflict occurs, sequence is thrown out instead of cut\n";
			} # if opts{N}
		} 
		print "B option is now $opts{B} and N option is now $opts{N}\n";
	} # if opts{B} or opts{N}

	#run fasta to cmap conversion 

	$dtStartStage = DateTime->now;
	eval{ mkpath("$opts{o}/fa2cmap") };	print "Couldn't create $opts{o}/fa2cmap: $@" if ($@);
	print "\nBeginning FASTA to CMAP conversion...\n";
	print "Using Enzyme: ".(join(" ", @enzymes))."\n";
	print "Minimum Length: $minLength Kb\nMinimum Labels: $minLabels\n";

	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	@$cmdRef = ($^X, "scripts/fa2cmap_multi_color.pl", "-i", $opts{n}, "-m", $minLabels, "-M", $minLength, "-o", "$opts{o}/fa2cmap", "-e");
	foreach my $enzyme (@enzymes)	{
		push(@$cmdRef, $enzyme, $channelNumber);
	} # foreach enzyme
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	open (OUT, ">$opts{o}/fa2cmap/fa2cmap.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/fa2cmap.errlog"); print ERR $errResults."\n"; close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "FASTA to CMAP conversion complete in $span.", "ERROR: FASTA to CMAP conversion cannot be completed.");


	# convert FASTA header to cmap id
	$dtStartStage = DateTime->now;
	print "\nBeginning FASTA header conversion...\n";
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");

	$filename_key = basename($opts{n});
	$filename_key =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename_key .= "_$enzyme"};
	$filename_key = $filename_key."_".$minLength."kb_".$minLabels."labels_key.txt";
	$filename_key = $opts{o}."/fa2cmap/".$filename_key;
	$filename = basename($opts{n});	$filename =~ s/(\S+)\.\w+$/$1/; foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"};
	$filename = $filename."_".$minLength."kb_".$minLabels."labels.cmap";
	$filename = $opts{o}."/fa2cmap/".$filename;

	@$cmdRef = ($^X, "scripts/fa_key_convert.pl", $opts{n}, $filename_key);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	open (OUT, ">$opts{o}/fa2cmap/faHeader_to_cmapId.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/faHeader_to_cmapId.errlog"); print ERR $errResults."\n"; close ERR;	
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check	
	errCheck($outResults, $errResults, "ERROR:", "FASTA header conversion complete in $span.", "ERROR: FASTA header conversion cannot be completed.");
	$filename =~ s/(\S+)\.\w+$/$1/;
	$temp = $filename."_CmapIdHeaders.fa";	
	print "\nNew FASTA with CMAP Ids as headers: $temp\n";

	# copy the fast file to the designated output directory fa2cmap sub-directory
	copy "$opts{n}", "$opts{o}/fa2cmap" ;
	$ngs_cmap = $filename.".cmap";

	$ngs_cmap = abs_path($ngs_cmap);
	print "NGS map path: $ngs_cmap\n";

	#make sure input NGS CMAP file contains at least one contit
	(undef, $numContig, undef) = readCMap($ngs_cmap);
	if (!($numContig > 0)) {
		dieLog( "ERROR: Input NGS file $ngs_cmap does not contain any contigs!\n\n");
		exit;
	}

	#perform initial alignment 
	$dtStartStage = DateTime->now;
	print "\nBeginning initial NGS CMAP to BioNano CMAP alignment...\n";
	$outputDir = $opts{o}."/align0";
	eval { mkpath($outputDir) };
	if ($@) {
		print "Couldn't create $outputDir: $@"; }

	copy "$opts{b}" , "$outputDir";
	
	my $endoutlier=$config{'align1.endoutlier'};
	my $outlier=$config{'align1.outlier'};
	my $extend=$config{'align1.extend'};
	my $FN=$config{'align1.FN'};
	my $FP=$config{'align1.FP'};
	my $sf=$config{'align1.sf'};
	my $sd=$config{'align1.sd'};
	my $sr=$config{'align1.sr'};
	my $res=$config{'align1.res'};
	my $resSD=$config{'align1.resSD'};
	my $mres=$config{'align1.mres'};
	my $A=$config{'align1.A'};
	my $biaswt=$config{'align1.biaswt'};
	my $M=$config{'align1.M'};
	my $Mfast=$config{'align1.Mfast'};
	my $deltaX=$config{'align1.deltaX'};
	my $deltaY=$config{'align1.deltaY'};
	my $xmapchim=$config{'align1.xmapchim'};
	my $xmapUnique=$config{'align1.xmapUnique'};
	my $RepeatMask=$config{'align1.RepeatMask'};  
	my $RepeatRec=$config{'align1.RepeatRec'}; 
	my $T=$config{'align1.T'};
	my $hash = $config{'align1.hash'};

	chdir $outputDir or dieLog( "ERROR: Cannot change directory to $outputDir: $!\n");
	@$cmdRef = ($refaligner, "-f", "-ref", $ngs_cmap, "-i", $opts{b}, "-o", "align0", "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-xmapchim", $xmapchim, "-xmapUnique", $xmapUnique,"-RepeatMask", split(/\s+/,$RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, split(/\s+/, $hash), "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheckRefAligner("align0.stdout", "END of output", "Initial alignment complete in $span.", "ERROR: Initial alignment cannot be completed.");

	my $align0_errbin = $outputDir."/align0.errbin";

	#rescale BNG input cmap
	$dtStartStage = DateTime->now;
	print "\nRescaling BioNano CMAP...\n";
	chdir $outputDir or dieLog( "ERROR: Cannot change directory to $outputDir: $!\n");
	my @suffixes = {".cmap",".CMAP"};
	$bppAdjust = $modBNGfile;	$bppAdjust =~ s/_cmap$//;	# remove the suffix
	$bppAdjust = $opts{o}."/align0/".$bppAdjust;
	@$cmdRef = ($refaligner, "-merge", "-i", $opts{b}, "-o", $bppAdjust, "-readparameters", $align0_errbin, "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheckRefAligner($bppAdjust.".stdout", "END of output", "Rescaling complete in $span.", "ERROR: Rescaling cannot be completed.");
	$bppAdjust = $bppAdjust.".cmap";


	#perform initial alignment using rescaled BNG
	$dtStartStage = DateTime->now;
	print "\nBeginning initial NGS CMAP to rescaled BioNano CMAP alignment...\n";
	$outputDir = $opts{o}."/align1";
	eval { mkpath($outputDir) };
	if ($@) {
	  print "Couldn't create $outputDir: $@"; }

	$endoutlier=$config{'align1.endoutlier'};
	$outlier=$config{'align1.outlier'};
	$extend=$config{'align1.extend'};
	$FN=$config{'align1.FN'};
	$FP=$config{'align1.FP'};
	$sf=$config{'align1.sf'};
	$sd=$config{'align1.sd'};
	$sr=$config{'align1.sr'};
	$res=$config{'align1.res'};
	$resSD=$config{'align1.resSD'};
	$mres=$config{'align1.mres'};
	$A=$config{'align1.A'};
	$biaswt=$config{'align1.biaswt'};
	$M=$config{'align1.M'};
	$Mfast=$config{'align1.Mfast'};
	$maxmem=$config{'global.maxmem'};
	$maxthreads=$config{'global.maxthreads'};
	$deltaX=$config{'align1.deltaX'};
	$deltaY=$config{'align1.deltaY'};
	$xmapchim=$config{'align1.xmapchim'};
	$xmapUnique=$config{'align1.xmapUnique'};
	$RepeatMask=$config{'align1.RepeatMask'};  
	$RepeatRec=$config{'align1.RepeatRec'}; 
	$T=$config{'align1.T'};
	$hash=$config{'align1.hash'};
	my $align1Reference = $ngs_cmap;
	my $align1Query = $bppAdjust;

	chdir $outputDir or dieLog( "ERROR: Cannot change directory to $outputDir: $!\n");
	@$cmdRef = ($refaligner, "-f", "-ref", $ngs_cmap, "-i", $bppAdjust, "-o", "align1", "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-xmapchim", $xmapchim, "-xmapUnique", $xmapUnique,"-RepeatMask", split(/\s+/,$RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, split(/\s+/, $hash), "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheckRefAligner("align1.stdout", "END of output", "Initial rescaled alignment complete in $span.", "ERROR: Initial rescaled alignment cannot be completed.");

	my $align1_xmap = $outputDir."/align1.xmap";
	my $align1_r_cmap = $outputDir."/align1_r.cmap";
	my $align1_q_cmap = $outputDir."/align1_q.cmap";
	$readparameters = abs_path("$opts{o}/align1/align1.errbin");


	#check to make sure that intial alignment XMAP contains alignments
	my @int_xmap = read_file($align1_xmap);
	my $headLines=0;
	my $alignLines=0;
	foreach (@int_xmap) {
		if ($_ =~ /^#/ ) {
			$headLines++; }
		else {
			$alignLines++; } }
	if ($alignLines < 1) {
		dieLog ("\nERROR: No intial alignments found between $opts{n} and $bppAdjust\n"); }
	else {
		print "\n$alignLines alignments found between $opts{n} and $bppAdjust\n"; }

	#run AssignAlignType script
	$dtStartStage = DateTime->now;
	print "\nBeginning AssignAlignType...\n";
	$outputDir = $opts{o}."/assignAlignType";
	eval { mkpath($outputDir) };
	if ($@) {
	  print "Couldn't create $outputDir: $@"; }
	my $assignAlignType_xmap = $outputDir."/assignAlignType.xmap";
	my $assignAlignType_r_cmap = $outputDir."/assignAlignType_r.cmap";
	my $assignAlignType_q_cmap = $outputDir."/assignAlignType_q.cmap";
	my $conflict_r_cmap = $outputDir."/conflict_ngs_r.cmap";
	my $conflict_q_cmap = $outputDir."/conflict_bng_q.cmap";
	my $breakPoint_txt = $outputDir."/conflicts.txt";
	my $T_cutoff=$config{'assignAlignType.T_cutoff'};
	$T_cutoff = -log10($T_cutoff);
	my $max_overhang=$config{'assignAlignType.max_overhang'};

	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	@$cmdRef = ("$^X", "./scripts/AssignAlignType.pl", $align1_xmap, $align1_r_cmap, $align1_q_cmap, $assignAlignType_xmap, $assignAlignType_r_cmap, $assignAlignType_q_cmap, $T_cutoff, $max_overhang, $ngs_cmap, $bppAdjust, $breakPoint_txt);

	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	open(OUT, ">$outputDir/assignAlignType.log") or dieLog ("ERROR: Cannot write to $outputDir/assignAlignType.log: $!\n");	print OUT "$outResults";	close OUT;
	open(ERR, ">$outputDir/assignAlignType.errlog") or dieLog ("ERROR: Cannot write to $outputDir/assignAlignType.errlog: $!\n");	print ERR "$errResults";	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "AssignAlignType complete in $span.", "ERROR: AssignAlignType cannot be completed.");

	#print number of BioNano and NGS contigs that have been flagged as conflicting
	my @bn; my @ngs;
	my @sticky = read_file($assignAlignType_xmap);
	for (my $i = 0; $i < scalar(@sticky); $i += 1)	{
		next if ($sticky[$i] =~ /^#/);
		my @s = split("\t", $sticky[$i]); 
		push (@bn, $s[1]);
		push (@ngs, $s[2]);
	}
	@bn = uniq @bn;
	@ngs = uniq @ngs;
	print scalar(@bn)." BNG contigs have been flagged as conflicting\n";
	print scalar(@ngs)." NGS contigs have been flagged as conflicting\n";

	#check AssignAlign outputs and Merge inputs contain at least one contig
	$assignAlign_r = abs_path("$opts{o}/assignAlignType/assignAlignType_r.cmap");
	(undef, $numContig, undef) = readCMap($assignAlign_r);
	if (!($numContig > 0)) {;
		print "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
		warn "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
	}
	$assignAlign_q = abs_path("$opts{o}/assignAlignType/assignAlignType_q.cmap");
	(undef, $numContig, undef) = readCMap($assignAlign_q);
	if (!($numContig > 0)) {
		print "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
		warn "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
	}

	# now perform cut at conflicts; run cut_conflicts.pl
	$dtStartStage = DateTime->now;
	print "\nBeginning cut_conflicts...\n";
	my $windowSize = $config{'cut_conflicts.window_size'};
	my $qScoreThreshold = $config{'cut_conflicts.min_quality_score_threshold'};
	my $covThreshold = $config{'cut_conflicts.min_coverage_threshold'};

	my $oriBNFile = $bppAdjust;	
	my @theOriBNFileContent = split(/\//, $oriBNFile);	my $theOriBNFileName = $theOriBNFileContent[$#theOriBNFileContent];	my $theOutBNFileName = $theOriBNFileName;	$theOutBNFileName =~ s/\.cmap$/_cut.cmap/;
	my $oriSeqFile = $ngs_cmap;	
	my @theOriSeqFileContent = split(/\//, $oriSeqFile);	my $theOriSeqFileName = $theOriSeqFileContent[$#theOriSeqFileContent];	my $theOutSeqFileName = $theOriSeqFileName;	$theOutSeqFileName =~ s/\.cmap$/_cut.cmap/;
	my $conflictFile = $breakPoint_txt;
	my $theOutputDir = "$outputDir/cut_conflicts";	mkdir $theOutputDir if (! -e $theOutputDir);
	$cutSeqFile = abs_path("$opts{o}/assignAlignType/cut_conflicts/$theOutSeqFileName");
	$cutBNFile = abs_path("$opts{o}/assignAlignType/cut_conflicts/$theOutBNFileName");

	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	@$cmdRef = ("$^X", "./scripts/cut_conflicts.pl", "-oriGMFile", $oriBNFile, "-oriSeqFile", $oriSeqFile, "-conflictFile", $conflictFile, "-outDir", $theOutputDir, "-outGMFile", $cutBNFile, "-outSeqFile", $cutSeqFile, "-windowSize", $windowSize, "-qScoreThreshold", $qScoreThreshold, "-covThreshold", $covThreshold, "-refAligner", $refaligner);

	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	open(OUT, ">$theOutputDir/cut_conflicts.log") or dieLog("ERROR: Cannot write to $theOutputDir/cut_conflicts.log: $!\n");	print OUT $outResults;	close OUT;
	open(ERR, ">$theOutputDir/cut_conflicts.errlog") or dieLog("ERROR: Cannot write to $theOutputDir/cut_conflicts.errlog: $!\n");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "cut_conflicts complete in $span.", "ERROR: cut_conflicts cannot be completed.");

	# print the number of cut contigs
	(undef, $numContig, undef) = readCMap($cutBNFile);
	print "$numContig BNG contigs are found after the cut-conflict step\n";
	(undef, $numContig, undef) = readCMap($cutSeqFile);
	print "$numContig NGS contigs are found after the cut-conflict step\n";
} else	{
	# if the user defined a conflict status file manually, then the pipeline starts from this point on
	
	# check for the existence of some important files
	$bppAdjust = $modBNGfile;	$bppAdjust =~ s/_cmap$/.cmap/;	$bppAdjust = "$opts{o}/align0/$bppAdjust"; dieLog ("ERROR: file does not exist: $bppAdjust\n") if (! -e $bppAdjust);
	$filename_key = basename($opts{n});	$filename_key =~ s/(\S+)\.\w+$/$1/; foreach my $enzyme (@enzymes) {$filename_key .= "_$enzyme"}; $filename_key = $filename_key."_".$minLength."kb_".$minLabels."labels_key.txt"; $filename_key = $opts{o}."/fa2cmap/".$filename_key; 
	dieLog("ERROR: fasta to cmap key file does not exist: $filename_key\n") if (! -e $filename_key);
	$filename = basename($opts{n});	$filename =~ s/(\S+)\.\w+$/$1/;	foreach my $enzyme (@enzymes) {$filename .= "_$enzyme"}; $filename = $filename."_".$minLength."kb_".$minLabels."labels.cmap"; $filename = $opts{o}."/fa2cmap/".$filename;
	$ngs_cmap = $filename;
	dieLog("ERROR: sequence cmap file does not exist: $ngs_cmap\n") if (! -e $ngs_cmap);
	my $fastaName = $opts{n};	my @fastaNameContent = split(/\//, $fastaName);	$fastaName = $fastaNameContent[$#fastaNameContent];	$fastaName = "$opts{o}/fa2cmap/$fastaName";
	dieLog("ERROR: the sequence fasta file has not been copied to the output fa2cmap directory: $fastaName\n") if (! -e $fastaName);
		
	# check for the presence of align0 error bin file
	$readparameters = abs_path("$opts{o}/align1/align1.errbin"); dieLog ("ERROR: file does not exist: $readparameters\n") if (! -e $readparameters);
	
	# now perform cut at conflicts; run cut_conflicts.pl
	print "\nBeginning cut_conflicts based on input a conflict file...\n";

	my $oriBNFile = $bppAdjust;	
	my @theOriBNFileContent = split(/\//, $oriBNFile);	my $theOriBNFileName = $theOriBNFileContent[$#theOriBNFileContent];	my $theOutBNFileName = $theOriBNFileName;	$theOutBNFileName =~ s/\.cmap$/_cut.cmap/;
	my $oriSeqFile = $ngs_cmap;	
	my @theOriSeqFileContent = split(/\//, $oriSeqFile);	my $theOriSeqFileName = $theOriSeqFileContent[$#theOriSeqFileContent];	my $theOutSeqFileName = $theOriSeqFileName;	$theOutSeqFileName =~ s/\.cmap$/_cut.cmap/;
	my $theOutputDir = "$opts{o}/assignAlignType/cut_conflicts_M$modifyNum";	mkdir $theOutputDir if (! -e $theOutputDir);
	$cutSeqFile = abs_path("$theOutputDir/$theOutSeqFileName");
	$cutBNFile = abs_path("$theOutputDir/$theOutBNFileName");
	
	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
	@$cmdRef = ("$^X", "./scripts/cut_conflicts.pl", "-oriGMFile", $oriBNFile, "-oriSeqFile", $oriSeqFile, "-outDir", $theOutputDir, "-outGMFile", $cutBNFile, "-outSeqFile", $cutSeqFile, "-modBkptStatusFile", $opts{M}, "-refAligner", $refaligner);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	open(OUT, ">$theOutputDir/cut_conflicts.log") or dieLog("ERROR: Cannot write to $theOutputDir/cut_conflicts.log: $!\n");	print OUT $outResults;	close OUT;
	open(ERR, ">$theOutputDir/cut_conflicts.errlog") or dieLog("ERROR: Cannot write to $theOutputDir/cut_conflicts.errlog: $!\n");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "cut_conflicts complete in $span.", "ERROR: cut_conflicts cannot be completed.");	

	# print the number of cut contigs
	(undef, $numContig, undef) = readCMap($cutBNFile);
	print "$numContig BNG contigs are found after the cut-conflict step\n";
	(undef, $numContig, undef) = readCMap($cutSeqFile);
	print "$numContig NGS contigs are found after the cut-conflict step\n";
} # if defined a conflict status file 

# run MergeNGS_BN
$dtStartStage = DateTime->now;
print "\nBeginning MergeNGS_BN...\n";
$outputDir = $opts{o}."/mergeNGS_BN";	$outputDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
eval { mkpath($outputDir) };
if ($@) {
	print "Couldn't create $outputDir: $@"; }
	
my $merge_Tvalue=$config{'mergeNGS_BN.merge_Tvalue'};
my $id_shift=$config{'mergeNGS_BN.id_shift'};
my $max_merge_rounds=$config{'mergeNGS_BN.max_merge_rounds'};
$maxthreads=$config{'global.maxthreads'};
my $endoutlier=$config{'mergeNGS_BN.endoutlier'};
my $outlier=$config{'mergeNGS_BN.outlier'};
my $biaswt=$config{'mergeNGS_BN.biaswt'};
my $sd=$config{'mergeNGS_BN.sd'};
my $res=$config{'mergeNGS_BN.res'};
my $sf=$config{'mergeNGS_BN.sf'};
my $RepeatMask=$config{'mergeNGS_BN.RepeatMask'};
my $RepeatRec=$config{'mergeNGS_BN.RepeatRec'};
my $pairmerge=$config{'mergeNGS_BN.pairmerge'};
my $minsites = 0;       # force the minsites to 0 to facilitate accurate accounting of the number of entries inputed/outputed
$maxmem = (! defined($config{'mergeNGS_BN.maxmem'})) ? ($config{'global.maxmem'}) : ($config{'mergeNGS_BN.maxmem'});	

### but first, we need to define a variable that stores the input files to the merge step
my ($mergeBN, $mergeNgs) = ("", "");
if (defined($opts{M}))	{
	# user defined a conflict file
	($mergeBN, $mergeNgs) = ($cutBNFile, $cutSeqFile);
	print "*Using user-defined conflict-cut BioNano and sequence CMAP files*\n";
} else	{
	if ($opts{B} == 1 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $ngs_cmap);
		print "*Using all BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 1 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $cutSeqFile);
		print "*Using all BioNano and conflict-cut sequence CMAP*\n";
	} elsif ($opts{B} == 1 && $opts{N} == 3)	{
		($mergeBN, $mergeNgs) = ($bppAdjust, $assignAlign_r);
		print "*Using all BioNano and non-conflicting-only sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $ngs_cmap);
		print "*Using conflict-cut BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $cutSeqFile);
		print "*Using conflict-cut BioNano and conflict-cut sequence CMAP*\n";
	} elsif ($opts{B} == 2 && $opts{N} == 3)	{
		($mergeBN, $mergeNgs) = ($cutBNFile, $assignAlign_r);
		print "*Using conflict-cut BioNano and non-conflicting-only sequence CMAP*\n";
	} elsif ($opts{B} == 3 && $opts{N} == 1)	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $ngs_cmap);
		print "*Using non-conflicting-only BioNano and all sequence CMAP*\n";
	} elsif ($opts{B} == 3 && $opts{N} == 2)	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $cutSeqFile);
		print "*Using non-conflicting-only BioNano and conflict-cut sequence CMAP*\n";
	} else	{
		($mergeBN, $mergeNgs) = ($assignAlign_q, $assignAlign_r);
		print "*Using non-conflicting-only BioNano and non-conflicting-only sequence CMAP*\n";
	} # if cut
} # if a modified conflict file is specified

# make sure that the input files contain at least one contig
(undef, $numContig, undef) = readCMap($mergeNgs);
if ($numContig < 1)	{
	dieLog("ERROR: before merge, sequence file $mergeNgs does not contain any contigs! You can try using all (no filter) sequence option.\n\n");
	exit;
} # if numContig
(undef, $numContig, undef) = readCMap($mergeBN);
if ($numContig < 1)	{
	dieLog("ERROR: before merge, BioNano file $mergeBN does not contain any contigs! You can try using all (no filter) BioNano option.\n\n");
	exit;
} # if numContig

# we need to have the id shift value to be bigger than the largest ngs cmap id value
my ($maxSeqId, undef) = readCmapIdLength($mergeNgs);
while ($id_shift < $maxSeqId)	{
	$id_shift = $id_shift * 10;
} # while id_shift

chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "-outputDir", $outputDir, "-refaligner", $refaligner, "-merge_Tvalue", $merge_Tvalue, "-scratchDir", "./", "-logFile", "mergeNGS_BN_script.log", "-ngs_cmap_fn", $mergeNgs, "-bng_cmap_fn", $mergeBN, "-id_shift", $id_shift, "-max_merge_rounds", $max_merge_rounds, "-maxthreads", $maxthreads, "-endoutlier", $endoutlier, "-outlier", $outlier, "-biaswt", $biaswt, "-sd", $sd, "-res", $res, "-sf", $sf, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-pairmerge", split(/\s+/, $pairmerge), "-readparameters", $readparameters, "-minsites", $minsites, "-maxmem", $maxmem);


print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/MergeNGS_BN.log") or dieLog ("ERROR: Cannot write to $outputDir/MergeNGS_BN.log: $!\n");	print OUT "$outResults";	close OUT;
open(ERR, ">$outputDir/MergeNGS_BN.errlog") or dieLog ("ERROR: Cannot write to $outputDir/MergeNGS_BN.errlog: $!\n");	print ERR "$errResults";	close ERR;

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "MergeNGS_BN complete in $span.", "ERROR: MergeNGS_BN cannot be completed");


#Align original NGS contigs to Hybrid map 

#find NGS contigs that went into hybrid assembly and make new NGS cmap
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
my $final_outputDir = "$opts{o}/align_final";	$final_outputDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");	
eval { mkpath($final_outputDir) };
if ($@) {
	print "Couldn't create $final_outputDir: $@"; }
my $hybridScaffoldFileName = "step2.hybrid.cmap";
copy "$outputDir/$hybridScaffoldFileName", "$final_outputDir/$hybridScaffoldFileName" ;

$dtStartStage = DateTime->now;
print "\nBeginning extraction of used and not NGS contigs in hybrid assembly...\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
if (defined($opts{M}) || $opts{N} == 2)	{
	# if the user indicated to cut
	@$cmdRef = ($^X, "scripts/find_used_not_used_ngs.pl", $outputDir, $id_shift, $cutSeqFile, $refaligner);
} else	{
	# if the user did not indicate cut, take the original input NGS 
	@$cmdRef = ($^X, "scripts/find_used_not_used_ngs.pl", $outputDir, $id_shift, $ngs_cmap, $refaligner);
} # if optsN
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/find_used_not_used_ngs.log") or dieLog ("ERROR: Cannot write to $outputDir/find_used_not_used_ngs.log: $!\n"); print OUT $outResults;  close OUT;
open(ERR, ">$outputDir/find_used_not_used_ngs.errlog") or dieLog ("ERROR: Cannot write to $outputDir/find_used_not_used_ngs.errlog: $!\n"); print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "Extraction of NGS contigs complete in $span.", "ERROR: Extraction of NGS contigs cannot be completed.");

copy "$outputDir/all_used_ngs.cmap", "$final_outputDir/filtered_NGS.cmap";
copy "$outputDir/all_not_used_ngs.cmap", "$final_outputDir/filtered_not_used_NGS.cmap";


#align just the used NGS contigs in hybrid back to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of filtered NGS cmap to Hybrid CMAP...\n";

my $prefixOrig = $BNGfile."_".$NGSfile."_HYBRID_SCAFFOLD";
my $prefix = $modBNGfile."_".$NGSfile."_NGScontigs_HYBRID_SCAFFOLD";

# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
my $filteredNGS_file = "$final_outputDir/filtered_NGS.cmap";
my $hybrid_cmap = "$final_outputDir/$hybridScaffoldFileName";

$endoutlier=$config{'align_final.endoutlier'};
$outlier=$config{'align_final.outlier'};
my $extend=$config{'align_final.extend'};
my $FN=$config{'align_final.FN'};
my $FP=$config{'align_final.FP'};
$sf=$config{'align_final.sf'};
$sd=$config{'align_final.sd'};
my $sr=$config{'align_final.sr'};
$res=$config{'align_final.res'};
my $resSD=$config{'align_final.resSD'};
my $mres=$config{'align_final.mres'};
my $A=$config{'align_final.A'};
$biaswt=$config{'align_final.biaswt'};
my $M=$config{'align_final.M'};
my $Mfast=$config{'align_final.Mfast'};
$maxmem=$config{'global.maxmem'};
$maxthreads=$config{'global.maxthreads'};
my $deltaX=$config{'align_final.deltaX'};
my $deltaY=$config{'align_final.deltaY'};
$RepeatMask=$config{'align_final.RepeatMask'}; 
$RepeatRec=$config{'align_final.RepeatRec'}; 
my $T=$config{'align_final.T'};
my $BestRef=$config{'align_final.BestRef'};
my $nosplit=$config{'align_final.nosplit'};
my $hash=$config{'align_final.hash'};

chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredNGS_file, "-o", $prefix, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix.stats", split(/\s+/, $hash), "-stdout", "-stderr");
print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
errCheckRefAligner("$prefix.stdout", "END of output", "Alignment of filtered NGS CMAP to Hybrid CMAP complete in $span.", "ERROR: alignment of filtered NGS CMAP to Hybrid CMAP cannot be completed.");


#Align original BNG contigs to Hybrid map 
#find BNG contigs that went into hybrid assembly and make new BNG cmap
$dtStartStage = DateTime->now;
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
print "\nBeginning extraction of used and not used BNG contigs in hybrid assembly...\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
if (defined($opts{M}) || $opts{B} == 2)	{
	# the user indicated to cut
	@$cmdRef = ($^X, "scripts/find_used_not_used_bn.pl", $outputDir, $id_shift, $cutBNFile, $refaligner);
} else	{
	# the user indicated not to cut (put in the original bpp-adjusted BioNano genome maps)
	@$cmdRef = ($^X, "scripts/find_used_not_used_bn.pl", $outputDir, $id_shift, $bppAdjust, $refaligner);
} # if optsB
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/find_used_not_used_bn.log") or dieLog ("ERROR: Cannot write to $outputDir/find_used_not_used_bn.log: $!\n");	print OUT $outResults;	close OUT;
open(ERR, ">$outputDir/find_used_not_used_bn.errlog") or dieLog ("ERROR: Cannot write to $outputDir/find_used_not_used_bn.errlog: $!\n");	print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "Extraction of BNG contigs complete in $span.", "ERROR: Extraction of BNG contigs cannot be completed.");

copy "$outputDir/all_used_bng.cmap", "$final_outputDir/filtered_BNG.cmap";
copy "$outputDir/all_not_used_bng.cmap", "$final_outputDir/filtered_not_used_BNG.cmap";

#align just the BNG contigs in hybrid back to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of filtered BNG cmap to Hybrid CMAP...\n";

my $prefixBNG = $modBNGfile."_".$NGSfile."_BNGcontigs_HYBRID_SCAFFOLD";

my $filteredBNG_file = "$final_outputDir/filtered_BNG.cmap";
$hybrid_cmap = "$final_outputDir/$hybridScaffoldFileName";

$endoutlier=$config{'align_final.endoutlier'};
$outlier=$config{'align_final.outlier'};
$extend=$config{'align_final.extend'};
$FN=$config{'align_final.FN'};
$FP=$config{'align_final.FP'};
$sf=$config{'align_final.sf'};
$sd=$config{'align_final.sd'};
$sr=$config{'align_final.sr'};
$res=$config{'align_final.res'};
$resSD=$config{'align_final.resSD'};
$mres=$config{'align_final.mres'};
$A=$config{'align_final.A'};
$biaswt=$config{'align_final.biaswt'};
$M=$config{'align_final.M'};
$Mfast=$config{'align_final.Mfast'};
$maxmem=$config{'global.maxmem'};
$maxthreads=$config{'global.maxthreads'};
$deltaX=$config{'align_final.deltaX'};
$deltaY=$config{'align_final.deltaY'};
$RepeatMask=$config{'align_final.RepeatMask'}; 
$RepeatRec=$config{'align_final.RepeatRec'}; 
$T=$config{'align_final.T'};
$BestRef=$config{'align_final.BestRef'};
$nosplit=$config{'align_final.nosplit'};
$hash=$config{'align_final.hash'};

chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredBNG_file, "-o", $prefixBNG, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefixBNG.stats", split(/\s+/, $hash), "-stdout", "-stderr");
print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
errCheckRefAligner("$prefixBNG.stdout", "END of output", "Alignment of filtered BNG CMAP to Hybrid CMAP complete in $span.", "ERROR: alignment of filtered BNG CMAP to Hybrid CMAP cannot be completed.");

#merge hybrid cmap with all NGS not participated in scaffolding process
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with NGS not participated in the hybrid scaffold...\n";
my $notUsedNGSCmap = "$outputDir/all_not_used_ngs.cmap";
my $hybridNotUsedNGS = "HYBRID_SCAFFOLD_notUsedNGS_merged.cmap";
chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $notUsedNGSCmap, "-o", "HYBRID_SCAFFOLD_notUsedNGS_merged", "-merge", "-minsites", 0, "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheckRefAligner("HYBRID_SCAFFOLD_notUsedNGS_merged.stdout", "END of output", "Merging Hybrid CMAP with naive NGS CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive NGS CMAP cannot be completed.");

#merge hybrid cmap with all BNG not participated in scaffolding process
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with naive BNG CMAP...\n";
my $notUsedBNGCmap = "$outputDir/all_not_used_bng.cmap";
my $hybridNotUsedBNG = "HYBRID_SCAFFOLD_notUsedBNG_merged.cmap";
chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $notUsedBNGCmap, "-o", "HYBRID_SCAFFOLD_notUsedBNG_merged", "-merge", "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheckRefAligner("HYBRID_SCAFFOLD_notUsedBNG_merged.stdout", "END of output", "Merging Hybrid CMAP with naive BNG CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive BNG CMAP cannot be completed.");

#generate AGP and FASTA file
# $outputDir is pointing to $opts{o}/mergeNGS_BN(_M$modifyNum)
# $final_outputDir is pointing to $opts{o}/align_final(_M$modifyNum)
$dtStartStage = DateTime->now;
my $agpFastaDir = "$opts{o}/agp_fasta";	$agpFastaDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");	
mkdir $agpFastaDir if (! -e $agpFastaDir);
my $fastaName = $opts{n};       my @fastaNameContent = split(/\//, $fastaName); $fastaName = $fastaNameContent[$#fastaNameContent];     $fastaName = "$opts{o}/fa2cmap/$fastaName";
print "\nBeginnning construction of AGP and FASTA file of the scaffolded and unscaffolded sequences...\n";
chdir $plDir or dieLog ("ERROR: Cannot change to $plDir: $!\n");
my $enzymeString = "";
foreach my $enzyme (@enzymes)	{
	$enzymeString .= "$enzyme $channelNumber ";	
} # foreach enzyme
$enzymeString =~ s/\s$//;
if (defined($opts{M}))	{
	# the user indicated to use manually cut results
	# enzyme argument here is somewhat weird
	my $cutNGSCoordFile = "$opts{o}/assignAlignType/cut_conflicts_M$modifyNum/auto_cut_NGS_coord_translation.txt"; 
	dieLog ("ERROR: file indicating NGS cut coordinates does not exists: $cutNGSCoordFile\n") if (! -e $cutNGSCoordFile);
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString, "-r", $cutNGSCoordFile);
} elsif ($opts{N} == 2)	{
	# no manual cut, auto cut was done instead
	my $cutNGSCoordFile = "$opts{o}/assignAlignType/cut_conflicts/auto_cut_NGS_coord_translation.txt";	
	dieLog ("ERROR: file indicating NGS cut coordinates does not exists: $cutNGSCoordFile\n") if (! -e $cutNGSCoordFile);
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString, "-r", $cutNGSCoordFile);
} else	{
	# the user indicated no cut to the results
	# enzyme argument here is somewhat weird
	@$cmdRef = ($^X, "scripts/ExportAGP.pl", "-i", "$final_outputDir/$prefix.xmap", "-c", "$final_outputDir/$prefix"."_r.cmap", "-o", $agpFastaDir, "-m", $filename_key, "-s", $fastaName, "-e", $enzymeString);	
} # if opts
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$agpFastaDir/xmap2agp.log") or dieLog ("ERROR: Cannot write to $agpFastaDir/xmap2agp.log: $!\n");	print OUT $outResults;	close OUT;
open(ERR, ">$agpFastaDir/xmap2agp.errlog") or dieLog ("ERROR: Cannot write to $agpFastaDir/xmap2agp.errlog: $!\n"); print ERR $errResults; close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "AGP and FASTA generation complete in $span.", "ERROR: AGP and FASTA generation cannot be completed.");


#Calculate some stats about CMAPs

$dtStartStage = DateTime->now;
print "\nCalculating statistics...\n\n";

print "Original BioNano Genome Map statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $opts{b});
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# errorCheck
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for the original BioNano CMAP file");
# print statistics
print "$outResults\n";

print "Bpp-adjusted BioNano Genome Map statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $bppAdjust);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# errorCheck
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for the bpp-adjusted BioNano CMAP file");
# print statistics
print "$outResults\n";

print "Original NGS Genome Map statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $ngs_cmap);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for original NGS CMAP file");
# print statistics
print "$outResults\n";

print "Before merge: BioNano Genome Map statistics:\n";
chdir $plDir or dieLog("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $mergeBN);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for BioNano Genome Map right before the Merge");
# print statistics
print "$outResults\n";

print "Before merge: NGS Genome Map statistics:\n";
chdir $plDir or dieLog("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $mergeNgs);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for NGS Genome Map right before the Merge");
# print statistics
print "$outResults\n";

print "BNG contigs in hybrid Genome Map statistics:\n";
$filteredBNG_file = abs_path($filteredBNG_file);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $filteredBNG_file);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for CMAP of BNG contigs used in hybrid scaffolding");
# print statistics
print "$outResults\n";

print "NGS contigs in hybrid Genome Map statistics:\n";
$filteredNGS_file = abs_path($filteredNGS_file);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $filteredNGS_file);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for CMAP of NGS contigs used in hybrid scaffolding");
# print statistics
print "$outResults\n";

print "Hybrid Genome Map statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/$hybridScaffoldFileName");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Hybrid CMAP");
# print statistics
print "$outResults\n";

print "The statistics of hybrid scaffold plus not scaffolded BNG:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/HYBRID_SCAFFOLD_notUsedBNG_merged.cmap");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for hybrid and not-scaffolded BNG");
# print statistics
print "$outResults\n";

print "The statistics of hybrid scaffold plus not scaffolded NGS:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$final_outputDir/HYBRID_SCAFFOLD_notUsedNGS_merged.cmap");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for hybrid and not-scaffolded NGS");
# print statistics
print "$outResults\n";

#calculate XMAP stats and output to .stats file

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "\nCalculating CMAP statistics complete in $span.\n";

$dtStartStage = DateTime->now;
print "\nCalculating XMAP statistics...\n";

chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $final_outputDir, "$prefix.xmap", $hybridScaffoldFileName, "$prefix.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$final_outputDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.log: $!\n");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix.xmap");
# print statistics
print "$outResults\n";

chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $final_outputDir, "$prefixBNG.xmap", $hybridScaffoldFileName, "$prefixBNG.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$final_outputDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.log: $!\n");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefixBNG.xmap");
# print statistics
print "$outResults\n";

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "XMAP statistics calculation complete in $span.\n";

### now copy the important results files to the final output directory ###
$hybridScaffoldFileName = $modBNGfile."_$NGSfile"."_HYBRID_SCAFFOLD.cmap";
copy "$final_outputDir/step2.hybrid.cmap", "$hybridScaffoldFinalResultsDir/$hybridScaffoldFileName";

# prefix is  $prefix = $modBNGfile."_".$NGSfile."_NGScontigs_HYBRID_SCAFFOLD";
copy "$final_outputDir/$prefix.xmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefix"."_q.cmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefix"."_r.cmap", "$hybridScaffoldFinalResultsDir";

# prefixBNG is $prefixBNG = $modBNGfile."_".$NGSfile."_BNGcontigs_HYBRID_SCAFFOLD";
copy "$final_outputDir/$prefixBNG.xmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefixBNG"."_q.cmap", "$hybridScaffoldFinalResultsDir";
copy "$final_outputDir/$prefixBNG"."_r.cmap", "$hybridScaffoldFinalResultsDir";

# copy the cut status text file, only if users specified cut was to be done (i.e. N 2 or B 2 or M)
if (defined($opts{M}) || $opts{N} == 2 || $opts{B} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";	$sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/conflicts_cut_status.txt", "$hybridScaffoldFinalResultsDir";
} # if to copy

# copy the bed files indicating whether the new fragments were cut or not
if (defined($opts{M}) || $opts{N} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";	$sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/ngs_pre_cut_annotations.bed", "$hybridScaffoldFinalResultsDir";
} # if to copy
if (defined($opts{M}) || $opts{B} == 2)	{
	my $sourceDir = "$opts{o}/assignAlignType/cut_conflicts";       $sourceDir .= "_M$modifyNum" if (defined($opts{M}));
	copy "$sourceDir/bn_pre_cut_projected_ngs_coord_annotations.bed", "$hybridScaffoldFinalResultsDir";
} # if to copy

# copy the AGP/FASTA/trim files 
my $theSourceDir = "$opts{o}/agp_fasta";	$theSourceDir .= "_M$modifyNum" if (defined($opts{M}));
copy "$theSourceDir/$prefix.agp", "$hybridScaffoldFinalResultsDir";
copy "$theSourceDir/$prefix.fasta", "$hybridScaffoldFinalResultsDir";
copy "$theSourceDir/$prefix"."_NOT_SCAFFOLDED.fasta", "$hybridScaffoldFinalResultsDir";
copy "$theSourceDir/$prefix"."_trimHeadTailGap.coord", "$hybridScaffoldFinalResultsDir";


# copy breakpoint file (delinearing conflicts between the sequence and genome maps), and align1 alignment files
$theSourceDir = "$opts{o}/assignAlignType";
copy "$theSourceDir/conflicts.txt", "$hybridScaffoldFinalResultsDir";	
$theSourceDir = "$opts{o}/align1"; 
my $prefixNGS_BNG = "$modBNGfile"."_$NGSfile"."_BNGcontigs_NGScontigs";
copy "$theSourceDir/align1.xmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG.xmap";
editReferenceQueryMapsHeader("$hybridScaffoldFinalResultsDir/$prefixNGS_BNG.xmap", "$prefixNGS_BNG");
copy "$theSourceDir/align1_r.cmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG"."_r.cmap";
copy "$theSourceDir/align1_q.cmap", "$hybridScaffoldFinalResultsDir/$prefixNGS_BNG"."_q.cmap";

###

# align molecules BNX to BioNano and hyrbid scaffold [optional]
if (defined($opts{m}))	{
	my $theRefAligner = $refaligner;
	my $globalMaxThreads = $config{'global.maxthreads'};
	my $deNovoPipelineDir = $opts{p};	# opts{p} should be absolute path
	my $deNovoOpArgsFile = $opts{q};	# opts{q} should be absolute path

	my ($moleculesDir, $moleculesFile) = ("", "");
	my @moleculesPath = split(/\//, $opts{m});
	$moleculesFile = $moleculesPath[$#moleculesPath];
	$moleculesDir = $opts{m};	$moleculesDir =~ s/\/$moleculesFile$//;

	my ($bioNanoMapDir, $bioNanoMapFile) = ("", "");
	my @bioNanoMapPath = split(/\//, $mergeBN);
	$bioNanoMapFile = $bioNanoMapPath[$#bioNanoMapPath];
	$bioNanoMapDir = $mergeBN;	$bioNanoMapDir =~ s/\/$bioNanoMapFile$//;
	my ($hybridMapDir, $hybridMapFile) = ($hybridScaffoldFinalResultsDir, $hybridScaffoldFileName);

	my $autoNoiseOutDir = "$opts{o}/auto_noise";	$autoNoiseOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
	my $alignMolBioNanoOutDir = "$opts{o}/alignmol_bionano";	$alignMolBioNanoOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");
	my $alignMolHybridOutDir  = "$opts{o}/alignmol_hybrid";	$alignMolHybridOutDir .= (defined($opts{M})) ? ("_M$modifyNum") : ("");

	@$cmdRef = ();
	$dtStartStage = DateTime->now;
	print "\nBeginning molecules alignment to BioNano genome maps and hybrid scaffolds...\n";
	chdir $plDir or dieLog ("ERROR: Cannot change to $plDir: $!\n");
	@$cmdRef = ($^X, "scripts/align_molecules.pl", "-refAligner", $refaligner, "-globalMaxThreads", $globalMaxThreads, "-deNovoPipelineDir", $deNovoPipelineDir, "-deNovoOpArgsFile", $deNovoOpArgsFile, "-moleculesDir", $moleculesDir, "-moleculesFile", $moleculesFile, "-bioNanoMapDir", $bioNanoMapDir, "-bioNanoMapFile", $bioNanoMapFile, "-hybridMapDir", $hybridMapDir, "-hybridMapFile", $hybridMapFile, "-autoNoiseOutDir", $autoNoiseOutDir, "-alignMolBioNanoOutDir", $alignMolBioNanoOutDir, "-alignMolHybridOutDir", $alignMolHybridOutDir);
	
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	open(OUT, ">$autoNoiseOutDir/align_molecules.log") or dieLog ("ERROR: Cannot write to $autoNoiseOutDir/align_molecules.log: $!\n");	print OUT $outResults;	close OUT;
	open(ERR, ">$autoNoiseOutDir/align_molecules.errlog") or dieLog ("ERROR: Cannot write to $autoNoiseOutDir/align_molecules.errlog: $!\n");	print ERR $errResults;	close ERR;
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "Alignment of molecules to BioNano genome maps and hybrid scaffolds complete in $span.", "ERROR: Alignment of molecules to BioNano genome maps and hybrid scaffolds cannot be completed.");
} # if defined m


print "\n";
my $dtEnd = DateTime->now;
print "End time: "; print join ' ', $dtEnd->ymd, $dtEnd->hms; print "\n\n";

$span = DateTime::Format::Human::Duration->new();
print 'Total elapsed time: ', $span->format_duration_between($dtEnd, $dtStart); print "\n";

print "\n\nMerging of $opts{b} with $opts{n} is complete.\n\n";


print "END of output\n";





######################################################################
#                           Subroutines                              #
######################################################################

# this subroutine is used to check RefAligner processes to see if it finishes successfully
# it assumes that -stderr and -stdout was used to run RefAligner
sub errCheckRefAligner	{
	my ($file, $completeString, $completeMsg, $dieMsg) = @_;
	open(IN, "$file") or dieLog ("ERROR: Cannot open $file: $!\n");
	my $presence = 0;
	while (my $line = <IN>)	{
		if ($line =~ /$completeString/)	{
			# if the line contains the string that indicates successful completion
			$presence = 1;
		} # if line
	} # while line
	if ($presence == 0)	{
		dieLog ("ERROR: $dieMsg\n");
	} else	{
		print "$completeMsg\n";
	} # if presence
	close IN;
} # errCheckRefAligner

# this subrountine is used to call non-RefAligner processes to see if there is error in the child
# it assumes that the child process prints out "ERROR:" tag to all messages
sub errCheck	{
	my ($outResults, $errResults, $errorString, $completeMsg, $dieMsg) = @_;
	if ($outResults =~ /$errorString/)	{
		dieLog ("ERROR: $dieMsg\n"); } 
	elsif ($errResults =~ /$errorString/)	{
		dieLog ("ERROR: $dieMsg\n"); } 
	else {
		print "$completeMsg\n"; } 
}

sub runCommand  {
        my ($argsRef) = @_;
#print "runCommand: before open3\n";
        my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef);
#print "runCommand: after open3\n";

        close($CMD_IN);

        my $outResults = "";
        my $errResults = "";
        my $sel = new IO::Select;
        $sel->add($out, $err);
        while(my @fhs = $sel->can_read) {
                foreach my $fh (@fhs) {
                        my $line = <$fh>;
                        unless(defined $line) {
                                $sel->remove($fh);
                                next;
                        } # unless line
                        if($fh == $out) {
                                $outResults .= "$line";
                                #print "$line";
                        }elsif($fh == $err) {
                                $errResults .= "$line";
                                #print "$line";
                        }else{
                                dieLog ("ERROR: This should never execute!");
                        } # if fh
                } # foreach fh
        } # while
#print "runCommand: before waitpid\n";
        my $ret=waitpid ($pid, 0); # reap the exit code
#print "runCommand: after waitpid\n";
        return ($outResults, "$errResults");
} # runCommand

sub Init{
	my $opt_string = 'hn:b:c:o:r:B:N:fm:vM:p:q:';
	if(!getopts("$opt_string", \%opts)){
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
	if ($opts{v}) {
		Version();
		exit; 
	} # if opts v
	if (! defined $opts{b} || $opts{b} !~ /\w+/ || ! -e $opts{b} )	{
		 print ("\nERROR: Please provide a valid BioNano genome map file, option -b\n");
		Usage();
	}
	if (! defined $opts{n} || $opts{n} !~ /\w+/ || ! -e $opts{n})	{
		print ("\nERROR: Please provide a valid sequence fasta file, option -n\n");
		Usage();
	} # if 
	if (! defined $opts{c} || $opts{c} !~ /\w+/ || ! -e $opts{c})	{
		print ("\nERROR: Please provide a valid XML configuration file for the hybrid scaffold pipeline, option -c\n");
		Usage();
	} # if 
	if (! defined $opts{r} || $opts{r} !~ /\w+/)	{
		print ("\nERROR: Please enter the RefAligner binary file, option -r\n");
		Usage();
	} # if 
	
	$opts{b} = abs_path($opts{b});
	$opts{n} = abs_path($opts{n}); 
	$opts{c} = abs_path($opts{c});
	$opts{r} = abs_path($opts{r});

	# output directory	
	if(defined($opts{o})) {
		$opts{o} = abs_path($opts{o}); 
	} # if 
	
	# user defined conflict resolution file, which can affect output directory, B and N options, and f overwrite option
	if (! defined($opts{M}))	{
		# if the user annotated conflict resolution file is not provided, then, must specify how to resolve conflict
		if ((! defined($opts{N}) || ! defined($opts{B})) || ($opts{N} !~ /^\d$/ || $opts{B} !~ /^\d$/ || !(1 <= $opts{B} && $opts{B} <= 3) || !(1 <= $opts{N} && $opts{N} <= 3)))	{
			print ("\nERROR: Please input in a conflict-resolution value between 1 and 3 for -N and -B to handle conflicts.\n");
			Usage();
		} # if optsB

		# a conflict status file was not given, then start anew
		if ($opts{f})	{
			rmtree($opts{o}) if (-d $opts{o});	# if exists an output directory
		} else	{
			dieLog( "\nERROR: Output folder: $opts{o} already exists! Use -f option to overwrite.\n") if (-d $opts{o});
		} # if force overwrite
		eval { mkpath($opts{o}) };
		if ($@) {
			print("\nERROR: Couldn't create $opts{o}: $@");
			print("\nERROR: Output folder invalid or does not exist! \n");
                	Usage();
        	} # if
	} else	{
		# check if the user annotated conflict resolution file exists
		if (! -e $opts{M})	{
			dieLog ("\nERROR: the conflict status file does not exist: $opts{M}\n");
		} # if not exists

		# a conflict status file was given, then there must already be data from a previous hybrid scaffold run (ignore -f)
		if (! -d $opts{o})	{
			dieLog("\nERROR: Output folder: $opts{o} does not exist! Please point to the output of a previous hybrid scaffold run.\n");
		} # if exists an output direcotry
	} # defined M

	# user wants to perform alignment of molecules to BioNano genome maps and Hybrid scaffolds
	if (defined($opts{m}))	{
		if ($opts{m} !~ /\w+/ || ! -e $opts{m})	{
			print ("\nERROR:  Please provide a valid molecules BNX file for alignment to BioNano genome maps and hybrid scaffolds, option -m\n");
			Usage();
		} # check molecules file
		if (! defined($opts{p}) || $opts{p} !~ /\w+/ || ! -d $opts{p})	{
			print ("\nERROR: alignment of molecules option has been indicated, but a valid de novo assembly pipeline script directory is not provided, option -p\n");
			Usage();
		} # check pipeline script directory
		if (! defined($opts{q}) || $opts{q} !~ /\w+/ || ! -e $opts{q})	{
			print ("\nERROR: alignment of molecules option has been indicated, but a valid de novo assembly pipeline parameter file (optArguments.xml) is not provided, option -q\n");
			Usage();
		} # check pipeline xml file
		
		$opts{m} = abs_path($opts{m});
		$opts{p} = abs_path($opts{p});
		$opts{q} = abs_path($opts{q});
	} # if defined m
} # Init

sub Version{
	my $dtStartStage = DateTime->now;

	my @revs;
	###my @pipelineFiles = process_files($plDir);
	# only process a shorter list of files to find the version number
	my @pipelineFiles = process_files_shortList($plDir);
	foreach (@pipelineFiles) {
		if (!-d $_) {
			#if ($_ !~ /svn/ && $_ !~ /perl5/) {
				open (FILE,$_);
				my @file = <FILE>; 
				close FILE;
				my $start = '$Id:'; my $end = '$';
				#print "$_\n";
				foreach (@file) {
					if ($_ =~ m/\$Id/) { 
						chomp($_);
						my ($wanted) = $_ =~ /\$Id:(.*)\$/;
						$wanted =~ s/^\s+|\s+$//g;
						push @revs, $wanted;
						#print $wanted."\n"; 
						last; 
						} } } } 
						#}
	
	my $revLine; my $revNum=0;
	@revs = uniq @revs;
	foreach (@revs) {
		my @split = split(/\s+/);
		#print "Version: $split[2]  Line: $_\n";
		if ($split[1] >= $revNum) {
			$revNum = $split[1];
			$revLine = $_; } }
	
	print "$revLine\n";  

	my $dtEndStage = DateTime->now;
	my $span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	#print "Total time: $span\n";
} # Version

sub process_files_shortList	{
	my ($thePath) = @_;
	my @theFiles = ();
	opendir(DIR, "$thePath") or die "process_files_shortList: cannot open dir $thePath: $!\n";
	my @topFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;	# add hybridScaffold.pl
	closedir DIR;
	foreach my $topFile (@topFiles)	{
		push(@theFiles, "$thePath/$topFile");
	} # foreach topFile

	# scripts sub-directory
	opendir(DIR, "$thePath/scripts") or die "process_files_shortList: cannot open dir $thePath/scripts: $!\n";
	my @secLevelFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;
	closedir DIR;	
	foreach my $secLevelFile (@secLevelFiles)	{
		push(@theFiles, "$thePath/scripts/$secLevelFile");
	} # foreach secLevelFile

	# BNG sub directory
	opendir(DIR, "$thePath/scripts/perl5/BNG") or die "process_files_shortList: cannot open dir $thePath/scripts/perl5/BNG: $!\n";
	my @thirdLevelFiles = grep {$_ =~ /\.pl$/ || $_ =~ /\.pm$/} readdir DIR;
	closedir DIR;
	foreach my $thirdLevelFile (@thirdLevelFiles)	{
		push(@theFiles, "$thePath/scripts/perl5/BNG/$thirdLevelFile");
	} # foreach thirdLevelFile

	return @theFiles;
} # process_files_shortList

sub process_files {
	my $path = shift;
    opendir (DIR, $path)
        or dieLog ("ERROR: Unable to open $path: $!");

    # We are just chaining the grep and map
    # This is the same as: LIST = map(EXP, grep(EXP, readdir()))
    my @files =
        # Third: Prepend the full path
        map { $path . '/' . $_ }
        # Second: take out '.' and '..'
        grep { !/^\.{1,2}$/ }
        # First: get all files
        readdir (DIR);

    closedir (DIR);

    for (@files) {
        if (-d $_) {
            # Add all of the new files from this directory
            # (and its subdirectories, and so on... if any)
            if (abs_path($_) !~ /5.10.1/ && abs_path($_) !~ /5.14.4/ && abs_path($_) !~ /svn/ && abs_path($_) =~ /scripts/) {
				push @files, process_files ($_); } }

        else {
            # Do whatever you want here =) .. if anything.
        }
    }
    # NOTE: we're returning the list of files
    return @files;
}

sub Usage{
	print << "EOF";
	
Usage: perl hybridScaffold.pl <-h> <-n ngs_file> <-b bng_cmap_file> <-c hybrid_config_xml> <-o output_folder> <-B conflict_filter_level> <-N conflict_filter_level> <-f> <-m molecules_bnx> <-p de_novo_pipeline> <-q de_novo_xml> <-v> 
      -h    : This help message         
      -n    : Input NGS FASTA or CMAP file [required]
      -b    : Input BioNano CMAP  [required]
      -c    : Merge configuration file [required]
      -o    : Output folder [required]
      -r    : RefAligner program [required]
      -B    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
      -N    : conflict filter level: 1 no filter, 2 cut contig at conflict, 3 exclude conflicting contig [required if not using -M option]
      -f    : Force output and overwrite any existing files
      -m    : Input BioNano molecules BNX for molecules to hybrid scaffold alignment [optional]
      -p    : Input de novo assembly pipeline script [required for -m option]
      -q    : Input de novo assembly pipeline optArguments XML script [required for -m option]
      -v    : Print pipeline version information
      -M    : Input a conflict resolution file indicating which NGS and BioNano conflicting contigs to be cut [optional] 
      
EOF
	exit;
} # Usage

sub numeric	{	$a	<=>	$b	}

sub Parse	{
	my ($config) = @_;
	my $tmp_key;
	my %tmpHash = ();

	my %hash;

my $toPrint = 0;
	foreach my $node (keys %$config)	{
		foreach my $item (keys %{$config->{$node}})	{
			%tmpHash = ();
			if (ref(($config->{$node}{$item})) eq 'ARRAY')	{
				# if there are multiple entries in each node
				foreach my $attr (@{$config->{$node}{$item}})	{
					%tmpHash = ();
					foreach my $value (keys %$attr)	{
						if ($value eq 'attr')	{
							$tmp_key = $node.".".$attr->{$value};
						} else	{
							if ($value =~ /val(\d+)/)	{
								$tmpHash{$1} = $attr->{$value};	
							} # if value
						} # if value eq attr
					} # foreach value

					# the end of one line, gather all the "val" pairs and put them in order as a string, and then store in the tmp_key bucket 
					# (e.g. store "0.7 0.6 1.4" in the "align1.RepeatRec" bucket)
					$hash{$tmp_key} = "";
					foreach my $tmpHashKey (sort numeric keys %tmpHash)	{
						$hash{$tmp_key} .= ($tmpHash{$tmpHashKey}." ");
					} # foreach tmpHashKey
					$hash{$tmp_key} =~ s/\s+$//;
				} # foreach attr
			} elsif (ref($config->{$node}{$item}) eq 'HASH')	{
				# the version node goes here, as it only has one entry
				foreach my $attr (keys %{$config->{$node}{$item}})	{
					if ($attr eq 'attr')	{
						$tmp_key = $node.".".$config->{$node}{$item}{$attr};
					} else	{
						if ($attr =~ /val(\d+)/)	{
							$tmpHash{$1} = $config->{$node}{$item}{$attr};
						} # if attr
					} # if attr eq attr
				} # foreach attr

				# the end of one line, gather all the "val" pairs, and put them in order as a string, and then store in the tmp_key bucket 
				# (e.g. store "$Id: hybridScaffold.pl 4304 2015-11-24 22:43:43Z apang $" in "version.version" bucket)
				$hash{$tmp_key} = "";
				foreach my $tmpHashKey (sort numeric keys %tmpHash)	{
					$hash{$tmp_key} .= ($tmpHash{$tmpHashKey}." ");
				} # foreach tmpHashKey
				$hash{$tmp_key} =~ s/\s+$//;
			} else	{
				dieLog("ERROR: in parsing XML configuration file!\n");
			} # if ref array
		} # foreach item
	} # foreach node

	return (%hash);
} # Parse

sub CHECK_Config{
	my ($config) = @_;
	foreach my $item (keys %$config){
		if($item =~ 'RepeatMask' || $item =~ 'RepeatRec' || $item =~ 'subset' || $item =~ 'resbias'){
			my @data = split(" ", $config->{$item});

			for (my $d = 0; $d < scalar(@data); $d += 1)	{
				if (!(defined($data[$d]) && length($data[$d])))	{
					dieLog("ERROR: Configuration file error 1, please make sure that $item is defined and has a value\n");
				} # if undefined or if empty string
			} # for d
		}
		else{
			# Test if the option exists and the value is defined
			if(!defined($config->{$item}) || !length($config->{$item})){
				dieLog ("ERROR: Configuration file error 2, please make sure that $item has a value.\n");
			}
		}
		
		# Test if the required numerical values are numerical
		my @non_num_params = qw(fasta2cmap.enzyme version hash pairmerge);
		push(@non_num_params, "global.refaligner");				# should be deprecated, as new XML no longer has refaligner field in it
		my $found = 0;
		for (my $i = 0; $i < scalar(@non_num_params); $i += 1)	{
			$found = 1 if ($item =~ /$non_num_params[$i]/);
		} # for i
		if($found == 0){
			if($item =~ 'RepeatMask' || $item =~ 'RepeatRec' || $item =~ 'subset' || $item =~ 'resbias'){
				my @data = split(" ", $config->{$item});
				for (my $d = 0; $d < scalar(@data); $d += 1)	{
					if(!looks_like_number($data[$d]) ){
						dieLog ("ERROR: Configuration file error 3, please make sure that $item is numerical.\n");
					} 
				} # for d
			}
			else{
				if(!looks_like_number($config->{$item})){
					dieLog ("ERROR: Configuration file error 4, please make sure that $item is numerical.\n");
				}
			}
		}
	}
}


sub split_xmap {
	my ($xmap, $rcmap, $qcmap, $contig_prefix, $molPath) = @_;
	 	
	my @xmap_in = read_file($xmap);	
	my @rcmap_in = read_file($rcmap);
	my @qcmap_in = read_file($qcmap);
	
	my %rcmap;
	my $rcmap_header="";
	my %qcmap;
	my $qcmap_header="";
	
	my @xmapheader;
	
	print "\nSplitting $xmap...\n";
	
	for (my $i=0; $i <= 4; $i++) {
		push (@xmapheader, $xmap_in[$i]); }
		
	eval { mkpath($molPath) };
	if ($@) {
		dieLog ("ERROR: Couldn't create $molPath: $@"); }
		
	
	foreach my $rline (@rcmap_in) {
	if ($rline =~ /^#/ ) {
		$rcmap_header = $rcmap_header . $rline; }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$Position,$StdDev,$Coverage,$Occurrence,$GmeanSNR,$lnSNRsd) = split("\t",$rline);
		if ($CMapId) {
			push @{ $rcmap{$CMapId} }, $rline; }}}
			
	foreach my $qline (@qcmap_in) {
	if ($qline =~ /^#/ ) {
		$qcmap_header = $qcmap_header . $qline;  }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$PositionStdDev,$Coverage,$Occurrence) = split("\t",$qline);
		if ($CMapId) {
			push @{ $qcmap{$CMapId} }, $qline; }}}
	
	
	my $prevContigID = 0;
	foreach my $xline (@xmap_in) {
		if ($xline !~ /^#/ ) {			
			my ($XmapEntryID, $QryContigID, $RefContigID, $QryStartPos, $QryEndPos, $RefStartPos, $RefEndPos, $Orientation, $Confidence, $HitEnum) = split("\t",$xline);  
			if ($RefContigID != $prevContigID) {
				
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				my $rcmapName = $contig_prefix.$RefContigID."_r.cmap";
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (RCMAP, ">>$molPath/$rcmapName");
				open (QCMAP, ">>$molPath/$qcmapName");
				
				foreach my $s (@xmapheader) {
					$s =~ s/^\s+//;
					print XMAP $s; }
				print XMAP "# Reference Maps From:	".$rcmapName."\n";
				print XMAP "# Query Maps From:	".$qcmapName."\n";
				print XMAP "#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum\n";
				print XMAP "#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string\n"; 
				print XMAP $xline;
				

				
				# RCMAP output
				print RCMAP $rcmap_header; 
				if($rcmap{$RefContigID}){
					print RCMAP join("", @{ $rcmap{$RefContigID} }); }
				
				#QCMAP output
				print QCMAP $qcmap_header;
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }
				
				$prevContigID = $RefContigID;
					
				close(QCMAP);
				close(RCMAP);
				close(XMAP);
				}				
				
			elsif ($RefContigID == $prevContigID) {
				
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (QCMAP, ">>$molPath/$qcmapName");				

				print XMAP $xline; 
				
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }

				close(QCMAP);
				close(XMAP);				
				
				} } }
	
	print "Splitting $xmap complete.\n"; 	
}

sub find_refaligner {
	if (! -d $_) {
		if ($_ =~ /RefAligner/i) {
			$refaligner = $File::Find::name; }}
}

sub log10 {
	my $n = shift;
    return log($n)/log(10);
}

# subroutine to determine the modify number for this current run (by looking at the largest modify number, if exists, in the "hybrid_scaffolds" output directory)
sub getCurMaxModifyNum	{
	my ($aDir) =@_;
	my $modifyNum = 1;
	opendir(DIR, $aDir) or dieLog ("ERROR: cannot access output directory $aDir: $!\n");
	my @alignFinalDirs = grep {$_ =~ /^hybrid_scaffolds/i} readdir DIR;
	closedir DIR;

	my @alignFinalModifyNum = ();
	for (my $i = 0; $i < scalar(@alignFinalDirs); $i += 1)	{
		if ($alignFinalDirs[$i] =~ /.+_M(\d+)$/)	{
			push(@alignFinalModifyNum, $1);
		} # if alignFinalDirs
	} # for i

	if (scalar(@alignFinalModifyNum) > 0)	{
		@alignFinalModifyNum = sort(@alignFinalModifyNum);
		$modifyNum = $alignFinalModifyNum[$#alignFinalModifyNum] + 1;
	} # if scalar

	return $modifyNum;
} # getCurMaxModifyNum

# this subroutine determines if chimeric quality scores exist in a given file
sub getQScores	{
	my ($file, $qScoresRef) = @_;
	# read in the genome map cmap file, ASSUMING that the chimQuality is on column 10
	my $noQScoreFlag = 0;
	open(IN, "$file") or dieLog "getQScores: cannot open $file: $!\n";
	my $foundChimQuality = -1;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# assumes that data line comes after header lines
			if ($line =~ /^#h\s+/)	{
				# find that column storing the ChimQuality
				$line =~ s/^#h\s+//;
				my @headerContent = split(/\t/, $line);
				for (my $i = 0; $i < scalar(@headerContent); $i += 1)	{
					$foundChimQuality = $i if ($headerContent[$i] =~ /ChimQuality/i);
				} # for i
			} # if line
			next;
		} # if line
		$line =~ s/^\s+/\t/;    $line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		if ($foundChimQuality == -1 || $content[$foundChimQuality] !~ /\d+/)	{
			warn "WARNING: getQScores cmap file=$file does not have chimeric quality score\n";
			$noQScoreFlag = 1;
			%$qScoresRef = ();	# empty any qScores currently stored, and just flag this file as having no qScores
			return ($qScoresRef, $noQScoreFlag);
		} # if scalar
		my ($id, $labelChannel, $start, $coverage, $chimQuality) = ($content[0], $content[4], $content[5], $content[7], $content[$foundChimQuality]);
		# skip if label channel is 0    (line indicating contig length)
		next if ($labelChannel == 0);
		push(@{$qScoresRef->{$id}}, {start => $start, coverage => $coverage, chimQuality => $chimQuality});
	} # while line
	close IN;
	return ($qScoresRef, $noQScoreFlag);
} # getQScores

sub editReferenceQueryMapsHeader	{
	my ($file, $headerSubString) = @_;
	# the sole purpose of this subroutine is to edit the # Reference Maps From: and # Query Maps From: header lines (so that IrysView will not complaining)
	my @lines = ();
	open(IN, $file) or dieLog("ERROR: editReferenceQueryMapsHeader: cannot read in $file: $!\n");
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^(#\s+Reference\s+Maps\s+From:)/)	{
			# _r.cmap
			$line = "$1\t$headerSubString"."_r.cmap";
		} elsif ($line =~ /^(#\s+Query\s+Maps\s+From:)/)	{
			# _q.cmap
			$line = "$1\t$headerSubString"."_q.cmap";
		} else	{
			# no op
		} # if line
		push(@lines, $line);
	} # while line
	close IN;

	open(OUT, ">$file") or dielog("ERROR: editReferenceQueryMapsHeader: cannot write to $file: $!\n");
	foreach my $line (@lines)	{
		print OUT "$line\n";
	} # foreach line
	close OUT;
} # editReferenceQueryMapsHeader

# this subrountine reads in a cmap file, and determines its largest cmap id
sub readCmapIdLength	{
	my ($file) = @_;
	# read in a cmap file, and determine the largest sequence id
	my $maxSeqId = -1;
	my %theLengths = ();
	open(IN, "$file") or dieLog "ERROR: readCmapId: reading in $file: $!\n";
	while (my $line = <IN>) {
		chomp $line;    $line =~ s/\r//g;
		next if ($line =~ /^#/);        # skip header lines
		my @content = split(/\t/, $line);
		my ($id, $aLength) = @content[0..1];
		$aLength = int($aLength);       # again change to integer
		$theLengths{$id} = $aLength;
		$maxSeqId = ($maxSeqId < $id) ? ($id) : ($maxSeqId);
	} # while line
	close IN;
	return ($maxSeqId, \%theLengths);
} # readCmapIdLength

