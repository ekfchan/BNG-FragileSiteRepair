#!/usr/bin/perl

# A wrapper script to run fragile site repair on a BioNano assembly as it aligns to a reference/NGS using an input FASTA and assembly CMAP
# A wrapper script to run fragile site repair on a BioNano single molecules as they align to a reference/NGS using an input FASTA and molecules BNX
# Maps are stitched together based on alignments if the maps start and stop overlap a fragile site

# Usage: perl fragileSiteRepair.pl --help

use strict;
use warnings;
use Getopt::Long; 
use Cwd qw(abs_path cwd);
use File::Basename;
use File::Path qw(rmtree mkpath);
use File::Copy;
use File::Find;
use File::Tee qw(tee);
use IPC::System::Simple qw(system capture);
#use Parallel::Iterator qw(iterate_as_array);
#use Parallel::Loops;
use Scalar::Util qw(looks_like_number);
use Parallel::ForkManager;
use DateTime;
use DateTime::Format::Human::Duration;
#use Data::Dumper;
use PerlIO::Util;
use Fcntl qw(:flock SEEK_END);
use Statistics::Lite qw(:all);
use subs qw(die);

my $dtStart = DateTime->now;
my $stime = DateTime->now;
my $etime = DateTime->now;
my $dtime = DateTime::Format::Human::Duration->new();
my $cpuCount=4;
my $jobs=1;
my $mem = 32;

# << usage statement and variable initialisation >>
my %inputs = (); 
$inputs{'force'}=0;
GetOptions( \%inputs, 'fasta:s', 'output=s', 'maxlab:i', 'maxfill:i', 'wobble:i', 'force', 'seq:i', 'cmap:s', 'n:i','j:i','optArgs:s', 'runSV', 'bnx:s', 'threshold:s', 'random', 'alignmolDir:s', 'break', 'h|help', 'endoutlier:s', 'minRatio:s','maxOverlap:i','aggressive', 'enzyme:s', 'bam:s', 'ngsBuffer:i', 'ngsBonus:i', 'breakNGSonly', 'gaps=s','maxOverlapLabels:i', 'pvalue:s', 'skipAlignMol', 'MapsMode','MoleculesMode', 'RefAligner:s','scripts','maxiter:i'); 

my $getSeq = 0;
my $threshold = 1.0;
my $ngsBonus = 10;

if ($inputs{h} || $inputs{help}) {
	print "\n";
	print qx/ps -o args $$/;
	print "\n";
	Usage();
	exit 0;
}

if ( (!$inputs{MapsMode} && !$inputs{MoleculesMode}) || ($inputs{MapsMode} && $inputs{MoleculesMode}) ) {
	print "ERROR: --MapsMode OR --MoleculesMode must be specified!\n";
	print qx/ps -o args $$/;
	print "\n";
	Usage();
	exit 1;
}

my $log_file = cwd()."/".basename($inputs{output})."_fragileSiteRepair_stdout_log.txt";
my $log_err_file = cwd()."/".basename($inputs{output})."_fragileSiteRepair_stderr_log.txt";
#for (*STDOUT, *STDERR)	{
#for (*STDOUT)	{
#	# captures the stdout
#	$_->autoflush;	$_->push_layer(tee=>$log_file)
#} 
#for (*STDERR)	{
#	# captures the stderr 
#	$_->autoflush;	$_->push_layer(tee=>$log_err_file)
#}
tee(STDOUT, '>', "$log_file");
tee(STDERR, '>', "$log_err_file");

print "\n";
print qx/ps -o args $$/;
print "\n";

if ( $inputs{MapsMode} && ( !exists $inputs{fasta} | !exists $inputs{cmap} | !exists $inputs{output} ) ) {
	print "ERROR: Required inputs for --MapsMode must be specified!\n";
	Usage();
	print "\n";
	exit 1; 
}
elsif ( $inputs{MoleculesMode} && ( !exists $inputs{fasta} | !exists $inputs{bnx} | !exists $inputs{output} ) ) {
	print "ERROR: Required inputs for --MoleculesMode must be specified!\n";
	Usage();
	print "\n";
	exit 1; 
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
	if ($inputs{MoleculesMode}) {
		$jobs = $cpuCount;
	}

	print "Maximum CPU cores to use: $cpuCount\n";
	print "Maximum parallel jobs: $jobs\n";
	print "Maximum memory to use: ", int($mem)," GB\n";
	print "\n";
	
	foreach my $key ("cmap","alignmolDir","output","fasta", "optArgs", "bnx", "bam", "gaps") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z|0-9]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
	if( !exists $inputs{maxfill} ) { $inputs{maxfill} = 30000; }
	if( !exists $inputs{wobble} ) { if ($inputs{MapsMode}) { $inputs{wobble} = 30000; } }
	if( !exists $inputs{maxlab} ) { $inputs{maxlab} = 1; }
		else {
			if ($inputs{maxlab} > 1 && $inputs{MapsMode}) {
				die "ERROR: Maximum recommended --maxlab = 1\n";
			}
			elsif ($inputs{maxlab}>0 && $inputs{MoleculesMode}) {
				die "ERROR: Maximum recommended --maxlab = 0\n";
			}
		}
	
	if( exists $inputs{seq} ) { 
		$getSeq = 1; 
	}
	if ($getSeq==1 && ($inputs{seq} < 1)) {
		$inputs{seq} = 50;
	}

	if( exists $inputs{threshold} && looks_like_number($inputs{threshold}) ) { $threshold = $inputs{threshold}; } 
	
	if( !exists $inputs{endoutlier} ) { $inputs{endoutlier} = "1e-5"; }
	
	if( !exists $inputs{minRatio} ) { $inputs{minRatio} = 0.70; }

	if( !exists $inputs{maxOverlap} ||  $inputs{maxOverlap}<0) { $inputs{maxOverlap} = 20000; }

	if( !exists $inputs{maxOverlapLabels} ||  $inputs{maxOverlapLabels}<0) { $inputs{maxOverlapLabels} = 5; }
	
	if ( !exists $inputs{enzyme} ) { $inputs{enzyme} = "GCTCTTC"; }

	if ( !exists $inputs{pvalue} ) { $inputs{pvalue} = "1e-12"; }
	
	if ( !exists $inputs{maxiter} ) { $inputs{maxiter} = 20; }

	if ( exists $inputs{ngsBonus} ) { $ngsBonus = int($inputs{ngsBonus}); }

	if ($inputs{MoleculesMode}) {
		if( !exists $inputs{maxOverlap} ||  $inputs{maxOverlap}<0) { $inputs{maxOverlap} = 10000; }
		if( !exists $inputs{maxOverlapLabels} ||  $inputs{maxOverlapLabels}<0) { $inputs{maxOverlapLabels} = 3; }
		if( !exists $inputs{maxfill} ) { $inputs{maxfill} = 10000; }
		if( !exists $inputs{wobble} ) { $inputs{wobble} = $inputs{maxfill}; }
		$inputs{maxlab} = 0; 
	}
	
}

my $splitprefix = "";
if ($inputs{MoleculesMode}) {
	$splitprefix = basename(abs_path($inputs{bnx}), ".bnx");
}
elsif ($inputs{MapsMode}) {
	$splitprefix = basename(abs_path($inputs{cmap}), ".cmap");
}
my @suffixlist = (".fa",".fasta",".fna");
my $splitfastaprefix = basename(abs_path($inputs{fasta}), @suffixlist);

# Make sure dependency scripts exist
my $scriptspath = abs_path(dirname($0));

if ( !exists $inputs{RefAligner} ) { $inputs{RefAligner} = $scriptspath."/tools/RefAligner"; } else { $inputs{RefAligner} = abs_path($inputs{RefAligner}); }

if ( !exists $inputs{scripts} ) { $inputs{scripts} = $scriptspath."/scripts/"; } else { $inputs{scripts} = abs_path($inputs{scripts}); }

if (!-e "$scriptspath/calcFragileSites_v2.pl" | !-e "$scriptspath/split_xmap_standalone.pl" | !-e "$scriptspath/gapFill_v3.pl" | !-e "$scriptspath/extractFASTA_fromBED.pl" | !-e "$scriptspath/runCharacterizeFinal.py" | !-e "$scriptspath/callSV.pl" | !-e "$scriptspath/scoreStitchPositionsBED_v5.pl" | !-e "$scriptspath/cutCmapByBedScore.pl" | !-e "$scriptspath/fa2cmap_multi_color.pl" | !-e "$scriptspath/alignmolAnalysis.pl") {
	die "ERROR: Dependency scripts not found at $scriptspath\n"; 
}
if( !-e "$inputs{RefAligner}" | !-e "$inputs{scripts}/HybridScaffold/scripts/calc_cmap_stats.pl" | !-e "$inputs{scripts}/optArguments_human.xml" | !-e "$inputs{scripts}/runSV.py" | !-e "$inputs{scripts}/runAlignMol.py") {
	die "ERROR: Pipeline scripts not found at $inputs{scripts} OR RefAligner not found at $inputs{RefAligner} : Please ensure ".$scriptspath."/tools and ".$scriptspath."/scripts are properly set up per IrysView Install Guide.\n\n"; 
	#die "abs_path('~/tools/RefAligner') not found: $!\n"; 
}

#print out input variables
if ($inputs{MoleculesMode}) {
	print "--MoleculesMode enabled\n\n";
}
else {
	print "--MapsMode enabled\n\n";
}

print "Input FASTA: $inputs{fasta}\n"; 
print "\tUsing enzyme: $inputs{enzyme}\n";
if ( $getSeq == 1) {
	print "\tSequence bp +/- start/end of fragile to print to BED: $inputs{seq}\n";
}
if (exists $inputs{aggressive}) {
	print "\tAggressive fragile site calculation enabled. Calculating TypeIII and TypeIV ...\n";
}

if ($inputs{MapsMode}) {
	print "Input assembly CMAP: $inputs{cmap}\n";
}
elsif ($inputs{MoleculesMode}) {
	print "Input molecules BNX: $inputs{bnx}\n";
	if (!-e $inputs{bnx}) {
		die "Inpux molecules BNX $inputs{bnx} supplied but does not exist!\n\n";
	}
}
print "Pvalue for alignments: $inputs{pvalue}\n";

if( !exists $inputs{optArgs} ) { $inputs{optArgs} = $inputs{scripts}."/optArguments_human.xml"; }
print "Using optArguments: $inputs{optArgs}\n";
copy("$inputs{optArgs}","$inputs{output}/");

if (exists $inputs{alignmolDir}) {
	if (-d $inputs{alignmolDir} && -e $inputs{alignmolDir}) {
		opendir(DIR, $inputs{alignmolDir});
		my @errList = grep(/\.err$/,readdir(DIR));
		closedir(DIR);
		$inputs{err} = "$inputs{alignmolDir}/$errList[0]";
		print "\nInput alignmolDir: $inputs{alignmolDir}\n";
		print "\tminRatio: $inputs{minRatio} for single molecule alignments at genome map ends\n";
	}
}
elsif (-e $inputs{bnx} && $inputs{MapsMode} && !$inputs{skipAlignMol} ) {
	print "\nInput alignmolDir not provided or does not exist!\n";
	print "\tSingle molecules will be performed using $inputs{bnx}\n";
	print "\tminRatio: $inputs{minRatio} for single molecule alignments at genome map ends\n";
}

print "\n";

#if (exists $inputs{bnx} && exists $inputs{err} && $inputs{MapsMode}) {
if (exists $inputs{bnx} && $inputs{MapsMode}) {
	print "Input molecules BNX used for scoring: $inputs{bnx}\n";
	if (-e $inputs{err}) {
		print "Input molecules ERR file: $inputs{err}\n";
	}
	else {
		print "Molecules ERR file not found!\n";
		print "\tSingle molecules will be performed using $inputs{bnx}\n";
	}
	print "\tUsing -endoutlier $inputs{endoutlier} for single molecule alignments\n"; }
else {
	print "\nInput BNX not provided and/or --MapsMode not enabled! Skipping single molecule scoring alignments...\n\n";
}

if (exists $inputs{bnx} && exists $inputs{err} && defined $inputs{break}) {print "\tBreaking maps at stitchPostions with score less than $threshold\n";}
if (exists $inputs{bnx} && exists $inputs{err} && defined $inputs{bam}) {print "\tUsing NGS alignments from $inputs{bam} to supplement single molecule scoring...\n";}
print "Output folder: $inputs{output}\n\n";
print "Maximum labels between maps: $inputs{maxlab}\n";
print "Maximum basepairs to fill between maps: $inputs{maxfill}\n";
print "Maximum fragile site wobble: $inputs{wobble}\n";
print "Maximum overlap bp between maps: $inputs{maxOverlap}\n";
print "Maximum overlap labels between maps: $inputs{maxOverlapLabels}\n";
print "\n";
if ($inputs{skipAlignMol}) { print "Skipping alignmol --skipAlignMol enabled!\n\n"; }
print "RefAligner path: $inputs{RefAligner}\n";
print "scripts path: $inputs{scripts}\n\n";

if (exists $inputs{runSV}) { print "SV module enabled\n\n"; }

# check output folder
if (-d $inputs{output} and $inputs{force} eq 0) {
	print "ERROR: Output directory $inputs{output} exists and force disabled. Use --force to overwrite.\n\n";
	die "ERROR: Output directory $inputs{output} exists and force disabled. Use --force to overwrite.\n";
}
elsif (-d $inputs{output} and $inputs{force} eq 1) {
	print "WARNING: Output directory $inputs{output} exists and will be overwritten\n";
	rmtree($inputs{output}) or die "ERROR: Cannot remove existing output directory $inputs{output}: $!\n";
}
mkpath($inputs{output}) or die "ERROR: Cannot create output directory $inputs{output}: $!\n";
print "\n";





# Step 1: Generate CMAP from input FASTA
print "=== Step 1: In-silico digest input FASTA into CMAP ===\n\n";
#Usage: fa2cmap_multi_color.pl [options] <Args>
#Options:
  #-h    : This help message
  #-v    : Verbose output  (Default: OFF)
  #-i    : Input fasta file  (Required)
  #-o    : Output folder  (Default: the same as the input file)
  #-e    : Names/sequences of the enzymes followed by channel # (Can be multiple)
  #-m    : Filter: Minimum labels  (Default: 0)
  #-M    : Filter: Minimum size (Kb)  (Default: 0)
my $cmd = "perl $scriptspath/fa2cmap_multi_color.pl -i $inputs{fasta} -o $inputs{output}/ -e $inputs{enzyme} 1";
print "Running command: $cmd\n";
system($cmd);
print "\n";

$inputs{ref} = abs_path("$inputs{output}/$splitfastaprefix"."_".uc($inputs{enzyme}).".cmap");
$inputs{key} = abs_path("$inputs{output}/$splitfastaprefix"."_".uc($inputs{enzyme})."_key.txt");

if ($inputs{MapsMode}) {
	copy("$inputs{cmap}","$inputs{output}") or print "WARNING: Copy of input assembly CMAP $inputs{cmap} failed: $!\n";
	$inputs{cmap} = "$inputs{output}/".$splitprefix.".cmap";
	print "\n";
}
#elsif ($inputs{MoleculesMode}) {
#	print "=== Convert BNX to CMAP ===\n";
#	$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix -merge -minsites 0 -minlen 0 -maxthreads $cpuCount -maxEnd 0.020";
#	print "Running command: $cmd\n";
#	system($cmd);
#	print "\n";
#	print "\n";
#	$inputs{cmap} = "$inputs{output}/$splitprefix".".cmap";
#}






if ($inputs{MapsMode}) {
	# Step 2: Align input assembly CMAP to input FASTA CMAP
	print "=== Step 2: Align input assembly CMAP to input reference FASTA CMAP ===\n\n";
	if (exists $inputs{ref} && defined $inputs{ref} && -e $inputs{ref}) {
		$cmd = "cp $scriptspath/runCharacterizeFinal.py $inputs{scripts}/; cd $inputs{output}; python $inputs{scripts}/runCharacterizeFinal.py -t $inputs{RefAligner} -r $inputs{ref} -q $inputs{cmap} -p $inputs{scripts}/ -n $cpuCount -a $inputs{optArgs} -v $inputs{pvalue}";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
	}
	move("$inputs{output}/alignref_final","$inputs{output}/alignref_initial");
	move("$inputs{output}/molecule_stats.txt", "$inputs{output}/alignref_initial/");
	$inputs{xmap} = abs_path("$inputs{output}/alignref_initial/$splitprefix".".xmap");
	$inputs{qcmap} = abs_path("$inputs{output}/alignref_initial/$splitprefix"."_q.cmap");
	$inputs{rcmap} = abs_path("$inputs{output}/alignref_initial/$splitprefix"."_r.cmap");
	$inputs{errbin} = abs_path("$inputs{output}/alignref_initial/$splitprefix".".errbin");
	$inputs{err} = abs_path("$inputs{output}/alignref_initial/$splitprefix".".err");
	$inputs{molStats} = abs_path("$inputs{output}/alignref_initial/molecule_stats.txt");

	#print "\n\t$inputs{xmap}\n\t$inputs{qcmap}\n\t$inputs{rcmap}\n\t$inputs{errbin}\n\n";
	print "\n";
}
elsif ($inputs{MoleculesMode}) {
	# Step 2: Align input molecules BNX to input FASTA CMAP
	print "=== Step 2: Align input molecules BNX to input FASTA CMAP ===\n\n";
	if (exists $inputs{ref} && defined $inputs{ref} && -e $inputs{ref}) {

		# first run autoNoise then use rescaled BNX (if exists) to perfrom stringent alignment to input FASTA CMAP
        	mkpath("$inputs{output}"."/autoNoise") or die "ERROR: Cannot create output directory $inputs{output}"."/autoNoise: $!\n";
       		print "\n";
        
		# autoNoise0
		$cmd = "cd $inputs{output}"."/autoNoise/; $inputs{RefAligner} -f -i $inputs{bnx} -ref $inputs{ref} -o autoNoise0 -stdout -stderr -MapRate 0.6 -usecolor 1 -A 7 -T 1e-4 -L 130 -BestRef 1 -BestRefPV 1 -mres 0.9 -res 3.3 -extend 1 -outlier 1e-5 -endoutlier 1e-4 -deltaX 4 -deltaY 6 -nosplit 2 -biaswt 0 -S -10000 -PVres 2 -AlignRes 1.5 -resEstimate -resbias 4.0 64 -f -MaxSE 0.5 -rres 1.2 -maptype 0 -outlierMax 40. -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -PVendoutlier -maxthreads $cpuCount -maxmem $mem -FP 2 -FN 0.2 -sf 0.2 -sd 0 -sr 0.02 -M 7 -maxEnd 0.020 -resbias 5.0 64 -relerr 0.001 -randomize 1 -subset 1 10000 -hashgen 5 3 3 1.7 0.05 5.0 1 1 2 -hash -hashdelta 10 -insertThreads 8";
        	print "\tRunning command: $cmd\n\n";
		system($cmd);
		if ( (-e "$inputs{output}"."/autoNoise/autoNoise0.err") && (-e "$inputs{output}"."/autoNoise/autoNoise0.errbin")) {
			$inputs{err} = "$inputs{output}"."/autoNoise/autoNoise0.err";
			$inputs{errbin} = "$inputs{output}"."/autoNoise/autoNoise0.errbin";
			#$inputs{molStats} = "$inputs{output}"."/autoNoise/molecule_stats.txt";
		}
                print "\n";

		# autoNoise1
		$cmd = "cd $inputs{output}"."/autoNoise/; $inputs{RefAligner} -f -i $inputs{bnx} -ref $inputs{ref} -o autoNoise1 -stdout -stderr -readparameters $inputs{errbin} -MapRate 0.6 -usecolor 1 -A 7 -T 1e-4 -L 130 -BestRef 1 -BestRefPV 1 -mres 0.9 -res 3.3 -extend 1 -outlier 1e-5 -endoutlier 1e-4 -deltaX 4 -deltaY 6 -nosplit 2 -biaswt 0 -S -10000 -PVres 2 -AlignRes 1.5 -resEstimate -resbias 4.0 64 -f -MaxSE 0.5 -rres 1.2 -maptype 0 -outlierMax 40. -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -PVendoutlier -maxmem $mem -maxthreads $cpuCount -M 3 3 -maxEnd 0.020 -resbias 5.0 64 -randomize 1 -ScanScaling 2 -subset 1 100000 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 2 -hash -hashdelta 10 -insertThreads 8 -XmapStatWrite molecule_stats.txt";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		if ( (-e "$inputs{output}"."/autoNoise/autoNoise1.err") && (-e "$inputs{output}"."/autoNoise/autoNoise1.errbin")) {
			$inputs{err} = "$inputs{output}"."/autoNoise/autoNoise1.err";
			$inputs{errbin} = "$inputs{output}"."/autoNoise/autoNoise1.errbin";
			$inputs{molStats} = "$inputs{output}"."/autoNoise/molecule_stats.txt";
		}
		if (-e "$inputs{output}"."/autoNoise/autoNoise1_rescaled.bnx") {
			$inputs{bnx} = "$inputs{output}"."/autoNoise/autoNoise1_rescaled.bnx";
		}
                print "\n";	

		my $errHash_ref = extractErr($inputs{err});
		my %errHash = %$errHash_ref;	

		# convert resclaed bnx to rescaled cmap
		if (-e "$inputs{output}"."/autoNoise/autoNoise1_rescaled.bnx") {
			print "=== Convert rescaled BNX to rescaled CMAP ===\n";

			$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_rescaled -merge -bpp $errHash{'bpp'} -minsites 0 -minlen 0 -maxthreads $cpuCount -bnx -stdout -stderr";
			print "Running command: $cmd\n";
			system($cmd);
			print "\n";
			$inputs{bnx} = "$inputs{output}/$splitprefix"."_rescaled.bnx";

			$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_rescaled -merge -minsites 0 -minlen 0 -maxthreads $cpuCount -maxEnd 0.20 -stdout -stderr";
			print "Running command: $cmd\n";
			system($cmd);
			print "\n";
			#unlink $inputs{cmap};
			$inputs{cmap} = "$inputs{output}/$splitprefix"."_rescaled.cmap";
			print "\n";
		}
		else {
			print "=== Convert BNX to rescaled CMAP ===\n";

			$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_rescaled -merge -bpp $errHash{'bpp'} -minsites 0 -minlen 0 -maxthreads $cpuCount -bnx -stdout -stderr";
			print "Running command: $cmd\n";
			system($cmd);
			print "\n";
			$inputs{bnx} = "$inputs{output}/$splitprefix"."_rescaled.bnx";

			$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_rescaled -merge -minsites 0 -minlen 0 -maxthreads $cpuCount -maxEnd 0.020 -stdout -stderr";
			print "Running command: $cmd\n";
			system($cmd);
			print "\n";
			#unlink $inputs{cmap};
			$inputs{cmap} = "$inputs{output}/$splitprefix"."_rescaled.cmap";
			print "\n";
		}

		# run alignmol 
		mkpath($inputs{output}."/alignmol") or die "ERROR: Cannot create output directory $inputs{output}/alignmol: $!\n";
		print "\n";

		# without hashing
		$cmd = "cd $inputs{output}/alignmol; $inputs{RefAligner} -ref $inputs{ref} -o $splitfastaprefix"."_$inputs{enzyme} -i $inputs{cmap} -f -stdout -stderr -maxthreads $cpuCount -usecolor 1 -FP 1.5 -FN 0.15 -sd 0.0 -sf 0.2 -sr 0.03 -res 3.3 -readparameters $inputs{errbin} -T $inputs{pvalue} -L 100 -nosplit 2 -biaswt 0 -extend 1 -BestRef 1 -maptype 0 -PVres 2 -HSDrange 1.0 -f -deltaX 4 -deltaY 6 -outlierExtend 2 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -resEstimate -outlier 1e-5 -I 20 -endoutlier 0.9 -E 0 -outlierMax 20. -maxEnd 0.020 -mres 0.9 -resbias 5.0 64 -MaxSE 0.5 -insertThreads 8 -rres 1.2 -maxmem $mem -XmapStatRead $inputs{molStats} -mapped $splitfastaprefix"."_$inputs{enzyme}"."_mapped -unmapped $splitfastaprefix"."_$inputs{enzyme}"."_unmapped -NoBpp";

		# with hashing
		$cmd = $cmd ." -hashoffset 1 -hashMultiMatch 20 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10";

		print "\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n";
		print "\n";
		
		
		#move("$inputs{output}/alignmol","$inputs{output}/alignmol_initial");
	    	$inputs{xmap} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}.xmap");
		$inputs{qcmap} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}_q.cmap");
		#$inputs{qcmap} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}_mapped.cmap");
		$inputs{rcmap} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}_r.cmap");
		#$inputs{err} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}.err");
		#$inputs{errbin} = abs_path("$inputs{output}/alignmol/$splitfastaprefix"."_$inputs{enzyme}.errbin");
		#$inputs{molStats} = abs_path("$inputs{output}/alignmol/molecule_stats.txt");

		#print "\n\t$inputs{xmap}\n\t$inputs{qcmap}\n\t$inputs{rcmap}\n\t$inputs{errbin}\n\n";
		print "\n";
	}
}




#if ($inputs{MapsMode}) {
	# Step 3: Split input XMAP into individual anchor maps
	print "=== Step 3: Split initial alignment XMAP into individual anchor XMAPs ===\n";
	# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
	$cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $splitprefix"."_contig $inputs{output}/contigs";
	if ($inputs{MoleculesMode}) {
		$cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $splitfastaprefix"."_$inputs{enzyme}"."_contig $inputs{output}/alignmol/merge";
	}
	elsif ($inputs{MapsMode}) {
		$cmd = "perl $scriptspath/split_xmap_standalone.pl $inputs{xmap} $inputs{qcmap} $inputs{rcmap} $splitprefix"."_contig $inputs{output}/contigs";
	}

	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	print "\n";
#}
#elsif ($inputs{MoleculesMode}) {
#	print "=== Step 3: Convert BNX to CMAP ===\n";
#	$cmd = "cd $inputs{output}; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix -merge -minsites 0 -minlen 0 -maxthreads $cpuCount";
#	print "Running command: $cmd\n";
#	system($cmd);
#	print "\n";
#	print "\n";
#	$inputs{cmap} = "$inputs{output}/$splitprefix".".cmap";
#}









## Step 4: Calculate potential fragile sites for input FASTA
print "=== Step 4: Calculate potential fragile sites for input FASTA ===\n";
my $bed; 

$stime = DateTime->now;
# Usage: perl calcFragileSites.pl --fasta <input FASTA> [--output <output.bed>] [--buffer <sequence buffer in bp>] [--random] [--enzyme <sequence to use to calculate fragile sites> [--agressive <calculate TypeIII and TypeIV fragile sites>]]
$cmd = "perl $scriptspath/calcFragileSites_v2.pl --fasta $inputs{fasta} --output $inputs{output} --enzyme $inputs{enzyme}";
if ($inputs{aggressive}) {
	$cmd = $cmd." --aggressive";
}
if ($getSeq == 1) {
	$cmd = $cmd." --buffer $inputs{seq}";
}
if (exists $inputs{random}) {
	$cmd = $cmd." --random";
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
print "\n";









print "=== Step 5: Run single molecule alignment analysis for each genome map ===\n";
if ($inputs{MapsMode}) {
	# Step 5: Run single molecule alignment analysis for each genome map
	#print "=== Step 5: Run single molecule alignment analysis for each genome map ===\n";
	if (-d $inputs{alignmolDir} && -e $inputs{alignmolDir}) {
		$stime = DateTime->now;
		my @xmaps = findXMAPs($inputs{alignmolDir});
		@xmaps = sort @xmaps;
		print "\n";
		print "Found ".scalar(@xmaps)." maps in $inputs{alignmolDir}\n\n";
		
		## normal sequential processing

		#foreach my $xmap (@xmaps) {
			#$xmap = abs_path("$inputs{alignmolDir}/$xmap");
			#my $base = $xmap; $base =~ s/.xmap//i;
			#my $qcmap = $base."_q.cmap";
			#my $rcmap = $base."_r.cmap";
			#if (-e $xmap and -e $qcmap and -e $rcmap) {
				##print "XMAP: $xmap\n";
				##print "QCMAP: $qcmap\n";
				##print "RCMAP: $rcmap\n";
				##print "\n";
				
				#$cmd = "perl $scriptspath/alignmolAnalysis.pl -x $xmap -q $qcmap -r $rcmap >> $inputs{output}/alignmolAnalysisOut.txt";
				#print "Running command: $cmd\n";
				#print "\n";
				#system($cmd);
			#}	
		#}

		# using Parallel::ForkManager
		print "Analyzing and processing genome maps in parallel fashion...\n";
		print "\n";
		my $pm =  new Parallel::ForkManager($cpuCount);
		#my $pm =  new Parallel::ForkManager($jobs);
		foreach my $xmap (@xmaps) {
			$pm->start and next; # do the fork
			# do work
			$xmap = abs_path("$inputs{alignmolDir}/$xmap");
			my $base = $xmap; $base =~ s/.xmap//i;
			#print "Base: $base\n";
			my $qcmap = $base."_q.cmap";
			my $rcmap = $base."_r.cmap";
			if (-e $xmap and -e $qcmap and -e $rcmap) {
				#print "XMAP: $xmap\n";
				#print "QCMAP: $qcmap\n";
				#print "RCMAP: $rcmap\n";
				#print "\n";
				
				$cmd = "perl $scriptspath/alignmolAnalysis.pl -x $xmap -q $qcmap -r $rcmap -t 10";
				#print "\tRunning command: $cmd\n";
				#print "\n";
				#system($cmd) or die "ERROR: $cmd failed: $!\n";
				my @result = capture($cmd);
				my $out = "$inputs{output}/alignmolAnalysisOut.txt";
				open OUT, "+>>$out" or die "ERROR: Cannot open $out for writing! $!\n";
				flock (OUT, LOCK_EX) or die "ERROR: Cannot open $out for locking! $!\n";
				seek (OUT, 0, 2);
				print OUT join("\n", @result);
				flock(OUT, LOCK_UN) or die "ERROR: Cannot unlock $out! $!\n";
				close OUT;
				#print "\n";
			}
			$pm->finish; # do the exit in the child process
		}
		print "Waiting for all processes to finish...\n";
		print "\n";
		$pm->wait_all_children;
		
		
		
		$etime = DateTime->now;
		$dtime = DateTime::Format::Human::Duration->new();
		print 'Spent time at this step: ', $dtime->format_duration_between($etime, $stime); print "\n\n";
	}
	else {
		print "--alignmolDir not specified. Skipping Step 5...\n";
	}
}
elsif ($inputs{MoleculesMode}) {
	print "\nStep 5 skipped due to --MoleculesMode\n\n";
}









# Step 6: Run gapFill.pl for each anchor map
print "=== Step 6: Run gapFill.pl for each anchor XMAP ===\n";
print "Processing anchor maps in parallel fashion...\n";
print "\n";
$stime = DateTime->now;

my @xmaps = "";
if ($inputs{MoleculesMode}) {
	@xmaps = findXMAPs($inputs{output}."/alignmol/merge");
}
elsif ($inputs{MapsMode}) {
	@xmaps = findXMAPs($inputs{output}."/contigs");
}
@xmaps = sort @xmaps;

## normal sequential processing

##foreach my $xmap (@xmaps) {
	##$xmap = abs_path($inputs{output}."/contigs/$xmap");
	##my $base = $xmap; $base =~ s/.xmap//i;
	###print "Base: $base\n";
	##my $qcmap = $base."_q.cmap";
	##my $rcmap = $base."_r.cmap";
	##if (-e $xmap and -e $qcmap and -e $rcmap) {
		##print "XMAP: $xmap\n";
		##print "QCMAP: $qcmap\n";
		##print "RCMAP: $rcmap\n";
		##print "\n";
		
		### usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		##$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -e $inputs{errbin} -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble} --minRatio $inputs{minRatio}";
		##print "\tRunning command: $cmd\n";
		##print "\n";
		###system($cmd) or die "ERROR: $cmd failed: $!\n";
		##my @result = capture($cmd);
		##my $log = "$base"."_fragileSiteRepaired_log.txt";
		###print join("\t", @result);
		##open LOG, ">$log" or die "ERROR: Cannot open $log for writing! $!\n";
		##print LOG join("\t", @result);
		##close LOG;
	##}	
##}


# using Parallel::ForkManager

#my $pm =  new Parallel::ForkManager($cpuCount);
my $pm =  new Parallel::ForkManager($jobs);
foreach my $xmap (@xmaps) {
	$pm->start and next; # do the fork
	# do work
	if ($inputs{MoleculesMode}) {
		$xmap = abs_path($inputs{output}."/alignmol/merge/$xmap");
	}
	elsif ($inputs{MapsMode}) {
		$xmap = abs_path($inputs{output}."/contigs/$xmap");
	}
	my $base = $xmap; $base =~ s/.xmap//i;
	#print "Base: $base\n";
	my $qcmap = $base."_q.cmap";
	my $rcmap = $base."_r.cmap";
	if (-e $xmap and -e $qcmap and -e $rcmap) {
		print "XMAP: $xmap\n";
		print "QCMAP: $qcmap\n";
		print "RCMAP: $rcmap\n";
		print "\n";
		
		# usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -e <errbin> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>] [--n CPU cores to use] [--alignmolAnalysis <alignmolAnalysisOut.txt>]
		$cmd = "perl $scriptspath/gapFill_v3.pl -x $xmap -q $qcmap -r $rcmap -e $inputs{errbin} --err $inputs{err} -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble} --n $cpuCount --minRatio $inputs{minRatio} --maxOverlap $inputs{maxOverlap} --maxOverlapLabels $inputs{maxOverlapLabels} --pvalue $inputs{pvalue} --RefAligner $inputs{RefAligner} --maxiter $inputs{maxiter}";
		if (exists $inputs{alignmolDir}) {
			$cmd = $cmd . " --alignmolAnalysis $inputs{output}/alignmolAnalysisOut.txt"; 
		}
		if ($inputs{MoleculesMode}) {
			$cmd = $cmd . " --MoleculesMode";
			$cmd = $cmd . " --molStats $inputs{molStats}"; 
		}
		elsif ($inputs{MapsMode}) {
			$cmd = $cmd . " --MapsMode"; 
			$cmd = $cmd . " --molStats $inputs{molStats}";
		}
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
if ($inputs{MapsMode}) {
	find(\&wanted, "$inputs{output}/contigs"); 
	find(\&wanted2, "$inputs{output}/contigs"); 
	find(\&wanted3, "$inputs{output}/contigs"); 
}
elsif ($inputs{MoleculesMode}) {
	find(\&wanted, "$inputs{output}/alignmol/merge"); 
	find(\&wanted2, "$inputs{output}/alignmol/merge"); 
	find(\&wanted3, "$inputs{output}/alignmol/merge"); 
}

$etime = DateTime->now;
$dtime = DateTime::Format::Human::Duration->new();
print 'Spent time at this step: ', $dtime->format_duration_between($etime, $stime); print "\n\n";

## using Parallel::Loops

##my $pl = Parallel::Loops->new($cpuCount);

##$pl->foreach( \@xmaps, sub {
	##my $xmap = $_;
	##$xmap = abs_path($inputs{output}."/contigs/$xmap");
	##my $base = $xmap; $base =~ s/.xmap//i;
	###$xmap = "$inputs{output}/contigs/$xmap";
	##my $qcmap = $base."_q.cmap";
	##my $rcmap = $base."_r.cmap";
	##if (-e $xmap and -e $qcmap and -e $rcmap) {
		##print "XMAP: $xmap\n";
		##print "QCMAP: $qcmap\n";
		##print "RCMAP: $rcmap\n";
		##print "\n";
		
		### usage: perl gapFill.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -o <output_prefix> [--bed <.bed fragile sites file>] [--round <start_round    =1>] [--maxlab <max_label_gap_tolerence=0>] [--maxfill <max basepairs to fill between contigs = 35000>] [--wobble <fragile site wobble in bp = 0>]
		##$cmd = "perl $scriptspath/gapFill.pl -x $xmap -q $qcmap -r $rcmap -o $base"."_fragileSiteRepaired --bed $bed --maxlab $inputs{maxlab} --maxfill $inputs{maxfill} --wobble $inputs{wobble}";
		##print "\tRunning command: $cmd\n";
		##print "\n";
		###system($cmd) or die "ERROR: $cmd failed: $!\n";
		##my @result = capture($cmd);
		##my $log = "$base"."_fragileSiteRepaired_log.txt";
		###print join("\t", @result);
		##open LOG, ">$log" or die "ERROR: Cannot open $log for writing! $!\n";
		##print LOG join("\t", @result);
		##close LOG;
	##}	
##});


# Step 7: Merge individual fragileSiteRepaired anchor maps and stitched BEDs
print "\n=== Step 7: Merge individual fragileSiteRepaired anchor maps ===\n\n";
my @qcmaps = "";
if ($inputs{MoleculesMode}) {
	@qcmaps = findQCMAPs($inputs{output}."/alignmol/merge");
}
elsif ($inputs{MapsMode}) {
	@qcmaps = findQCMAPs($inputs{output}."/contigs");
}
@qcmaps = sort @qcmaps;

#generate file with list of maps to be merged
my $mergeFile = "$inputs{output}/contigs/mergeList.txt";
if ($inputs{MoleculesMode}) {
	$mergeFile = "$inputs{output}/alignmol/merge/mergeList.txt";
}
open (LIST, ">$mergeFile") or die "ERROR: Could not open $mergeFile: $!\n";
print "Merging ".scalar(@qcmaps)." fragileSiteRepaired anchor maps...\n\n";
foreach (@qcmaps) {
	#print "QCMAP: $_\n";
	if ($inputs{MoleculesMode}) {
		$_ = abs_path($inputs{output}."/alignmol/merge/$_");
	}
	elsif ($inputs{MapsMode}) {
		$_ = abs_path($inputs{output}."/contigs/$_");
	}
	print LIST "$_\n";
}
close LIST;

#my $input = join(" -i ",@qcmaps);
#$cmd = "cd $inputs{output}; ~/tools/RefAligner -i $input -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0";

mkpath($inputs{output}."/merged_cmaps") or die "ERROR: Cannot create output directory $inputs{output}/merged_cmaps: $!\n";
#print "\n";
my $finalmap = "$inputs{output}/merged_cmaps/$splitprefix"."_fragileSiteRepaired.cmap";
if ($inputs{MapsMode}) {
	if ($inputs{NGSMode}) {
		$cmd = "cd $inputs{output}/merged_cmaps; $inputs{RefAligner} -f -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired_unmerged -minsites 0 -minlen 0 -stdout -stderr -maxthreads $cpuCount";
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";

		$cmd = "cd $inputs{output}/merged_cmaps; $inputs{RefAligner} -f -i $splitprefix"."_fragileSiteRepaired_unmerged.cmap -pairmerge 100 0.1 -T 1e-15 -mres 1e-3 -FP 0.2 -FN 0.02 -sr 0.01 -sd 0.0 -sf 0.2 -outlier 0 -endoutlier 1e-4 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -biaswt 0 -A 10 -HSDrange 1.0 -pairmergeRepeat -f -maxthreads $cpuCount -o $splitprefix"."_fragileSiteRepaired -minsites 0 -stdout -stderr";
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";

		my @pairmergemaps = findCMAPs($inputs{output}."/merged_cmaps","$splitprefix"."_fragileSiteRepaired_contig");
		@pairmergemaps = sort @pairmergemaps;
		$mergeFile = "$inputs{output}/merged_cmaps/mergeList.txt";
		open (LIST, ">$mergeFile") or die "ERROR: Could not open $mergeFile: $!\n";
		print "Merging ".scalar(@pairmergemaps)." fragileSiteRepaired genome maps...\n\n";
		foreach (@pairmergemaps) {
			#print "QCMAP: $_\n";
			$_ = abs_path($inputs{output}."/merged_cmaps/$_");
			print LIST "$_\n";
		}
		close LIST;

		$cmd = "cd $inputs{output}/merged_cmaps; $inputs{RefAligner} -f -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0 -minlen 0 -stdout -stderr -maxthreads $cpuCount";
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";
	}
	else {
		$cmd = "cd $inputs{output}/merged_cmaps; $inputs{RefAligner} -f -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0 -minlen 0 -stdout -stderr -maxthreads $cpuCount";
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";
	}
	
	$finalmap = "$inputs{output}/merged_cmaps/$splitprefix"."_fragileSiteRepaired.cmap";
}
elsif ($inputs{MoleculesMode}) {
	$cmd = "cd $inputs{output}/merged_cmaps; $inputs{RefAligner} -f -if $mergeFile -merge -o $splitprefix"."_fragileSiteRepaired -minsites 0 -minlen 0 -stdout -stderr -maxthreads $cpuCount";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n";
	
	$finalmap = "$inputs{output}/merged_cmaps/$splitprefix"."_fragileSiteRepaired.cmap";
}



print "\n\n";


print "Generating merged stitchPositions BED...";
my $mergedBED = ""; 
if ($getSeq == 1) {
	$mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.withSeqs.bed";
}
else {
	$mergedBED = "$inputs{output}/$splitprefix"."_fragileSiteRepaired_stitchPositions.bed";
}
open BEDOUT, ">$mergedBED" or die "ERROR: Could not open $mergedBED: $!\n";
my @beds = "";
if ($inputs{MoleculesMode}) {
	@beds = findBEDs("$inputs{output}/alignmol/merge");
}
elsif ($inputs{MapsMode}) {
	@beds = findBEDs("$inputs{output}/contigs");
}
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
my ($bedOut_ref, $typeCount_ref, $typeConf_ref, $typeConfStdDev_ref) = sortBEDandCount(\@bedOut);
@bedOut = @$bedOut_ref;
my %typeCount = %$typeCount_ref;
my %typeConf = %$typeConf_ref;
my %typeConfStdDev = %$typeConfStdDev_ref;
my $stitchCount = scalar(@bedOut);
print BEDOUT join("\n",@bedOut);
close BEDOUT; 
print "done\n\n";
print "\n";








# Step 8: Generate new merged cmap with unmapped maps
print "=== Step 8: Process original assembly/molecules CMAP and generate new merged CMAP ===\n\n";
if (exists $inputs{cmap}) {
	my @origCmapIds = getCmapIds($inputs{cmap});
	my @qryCmapIds = getCmapIds($inputs{qcmap});
	#if ($inputs{MoleculesMode}) {
	#	@qryCmapIds = getCmapIds($);
	#}
	my %in_qryCmapIds = map {$_ => 1} @qryCmapIds;
	my @diff  = grep {not $in_qryCmapIds{$_}} @origCmapIds;
	# my $oldMap = abs_path($inputs{output}."/".$splitprefix."_fragileSiteRepaired.cmap");
	
	@diff = sort @diff;
	$mergeFile = "$inputs{output}/merged_cmaps/idList.txt";
	open (LIST, ">$mergeFile") or die "ERROR: Could not open $mergeFile: $!\n";
	#print "Merging ".scalar(@pairmergemaps)." fragileSiteRepaired genome maps...\n\n";
	foreach (@diff) {
		#print "QCMAP: $_\n";
		#$_ = abs_path($inputs{output}."/merged_cmaps/$_");
		print LIST "$_\n";
	}
	close LIST;

		my $oldMap = $finalmap; 
		my $tempMap = $inputs{output}."/temp";
		my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged";
		#print "\tcmapIds in orig but not ref: ".join("\n\t",@diff)."\n";
		if( scalar(@diff) == 0 ) { 
			print "No unmapped cmaps were found in assembly.cmap. No merges necessary. \n$oldMap will be copied as is to $newMap.cmap\n"; 
			copy( $oldMap, $newMap.".cmap" ) or print "WARNING: Unable to copy $oldMap to $newMap.cmap. Step 7 may fail\n"; 
		}
		else {
			$cmd = "cd $inputs{output}; $inputs{RefAligner} -f -merge -i $inputs{cmap} -selectidf $mergeFile -o $tempMap -maxthreads $cpuCount -minsites 0 -minlen 0";
			print "\tRunning command: $cmd\n\n";
			system($cmd);
			print "\n";
			$cmd = "cd $inputs{output}; $inputs{RefAligner} -f -merge -i $tempMap".".cmap -i $oldMap -o $newMap -minsites 0 -minlen 0 -maxthreads $cpuCount; rm -f $tempMap".".cmap";
			print "\tRunning command: $cmd\n\n";
			system($cmd);
			print "\n";
		}
		$finalmap = $newMap.".cmap";
	}
	print "\n";
	
	
	
	
	
if ($inputs{MapsMode}) {

	# Step 9: run CharacterizeFinal on newly merged cmap
	print "=== Step 9: Align fragileSiteRepaired merged CMAP to reference and get stats ===\n\n";
	if (-e $finalmap && (exists $inputs{ref})) {
		#copy("$inputs{ref}","$inputs{output}") or print "WARNING: Copy of input reference $inputs{ref} failed: $!\n"; 
		# my $newMap = "$inputs{output}"."/"."$splitprefix"."_fragileSiteRepaired_merged.cmap";
		#usage: runCharacterizeFinal.py [-h] [-t REFALIGNER] [-r REFERENCEMAP] [-q QUERYMAP] [-x XMAP] [-p PIPELINEDIR] [-a OPTARGUMENTS] [-n NUMTHREADS] 
		# $cmd = "cp runCharacterizeFinal.py ~/scripts/; python ~/scripts/runCharacterizeFinal.py -t ~/tools/RefAligner -r $inputs{ref} -q $newMap -p ~/scripts/ -n $cpuCount";
		$cmd = "cp $scriptspath/runCharacterizeFinal.py $inputs{scripts}/; cd $inputs{output}; python $inputs{scripts}/runCharacterizeFinal.py -t $inputs{RefAligner} -r $inputs{ref} -q $finalmap -p $inputs{scripts} -n $cpuCount -a $inputs{optArgs} -v $inputs{pvalue}";
		print "\tRunning command: $cmd\n\n";
		system($cmd);
	}
	print "\n";
	my $alignreffinalXmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged.xmap";
	my $alignreffinalQcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_q.cmap";
	my $alignreffinalRcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_r.cmap";


	# Step 10: align molecules to merged fragileSiteRepaired CMAP
	my $alignmolDir="";
	my $alignmolErrbin="";
	print "=== Step 10: Run single molecule alignments against fragileSiteRepaired merged cmap ===\n\n";
	my $alignmolXmap="";
	my $alignmolRmap="";
	if (!$inputs{skipAlignMol}) {
	if (defined $inputs{bnx} && exists $inputs{bnx} && defined $inputs{err} && exists $inputs{err}) {
		print "Running alignment of $inputs{bnx} to $finalmap\n\n";
		mkpath("$inputs{output}/alignmol") or die "ERROR: Cannot create output directory $inputs{output}/alignmol: $!\n";
		$alignmolDir = abs_path("$inputs{output}/alignmol");
		my $errHash_ref = extractErr($inputs{err});
		my %errHash = %$errHash_ref;	
		# $ cd /home/palak; /home/palak/tools/RefAligner.mic -ref /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -o /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/contigs/exp_refineFinal1/alignmol/EXP_REFINEFINAL1_1 -i /home/palak/data/HuRef_fragileSiteRepair_IrysView/HuRef_105x_hg19/output/all_1_of_30.bnx -f -stdout -stderr -maxthreads 240 -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -maxmem 7
		#$cmd = "cd $alignmolDir; ~/tools/RefAligner -f -ref $finalmap -o $splitprefix"."_fragileSiteRepaired_merged_alignmol -i $inputs{bnx} -maxthreads $cpuCount -maxmem ".($mem)." -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt\$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -stdout -stderr";
		
	$cmd = "cd $alignmolDir; $inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_fragileSiteRepaired_merged_alignmol -ref $finalmap -f -stdout -stderr -maxthreads $cpuCount -usecolor 1 -FP $errHash{'FP(/100kb)'} -FN $errHash{'FN(rate)'} -sd $errHash{'sd'} -sf $errHash{'sf'} -sr $errHash{'sr'} -res $errHash{'res(pixels)'} -bpp $errHash{'bpp'} -output-veto-filter intervals.txt\$ -T 1e-10 -L 130 -nosplit 2 -biaswt 0 -resSD 0.75 -extend 1 -BestRef 1 -maptype 0 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -deltaX 4 -deltaY 6 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -outlier 0.00001 -endoutlier $inputs{endoutlier} -outlierMax 40. -cres 5.4 3 0.1 -MaxSE 0.5 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 -M 3 3 -MaxSE 0.5 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -maxmem $mem";
		
		print "\tRunning command: $cmd\n\n";
		system($cmd);
		sleep(5);
		
		#split xmap into _contig maps
		$alignmolXmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol.xmap";
		my $alignmolQmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol_q.cmap";
		$alignmolRmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol_r.cmap";
		$alignmolErrbin = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol.errbin";
			$inputs{err} = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_alignmol.err";
		my $alignmolPrefix = "$splitprefix"."_fragileSiteRepaired_merged_alignmol_contig";
		#mkpath("$alignmolDir/contigs") or die "ERROR: Cannot create output directory $inputs{output}/alignmol/contigs: $!\n";
		# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
		$cmd = "perl $scriptspath/split_xmap_standalone.pl $alignmolXmap $alignmolQmap $alignmolRmap $alignmolPrefix $alignmolDir/contigs";
		print "\n\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n";
		copy("$inputs{err}","$alignmolDir/contigs/") or die "\nERROR: Copy of $inputs{err} failed: $!\n";
	}
	else {
		print "Skipping Step 10 Single molecule alignments. Input BNX and/or molecules ERR not provided or does not exist! $!\n";
		print "\n";
	}
	}
	else {
		print "Skipping Step 10 Single molecule alignments. --skipAlignMol enabled! $!\n";
		print "\n";
	}

	print "\n";


	# Step 11: score BED file
	print "=== Step 11: Score stitchPositions BED based on single molecule alignments ===\n\n";

	my $scoredBED = basename($mergedBED,".bed");
	my $scoredSTATS = $scoredBED;
	$scoredBED = $scoredBED."_scored.bed";
	$scoredSTATS = $scoredSTATS."_scoreStats.csv";

	if (!$inputs{skipAlignMol}) {
	if (defined $inputs{bnx} && exists $inputs{bnx} && -e $finalmap && exists $inputs{ref} && -e $alignmolErrbin && -e $alignmolXmap && -e $alignmolRmap) {
		print "Scoring BED file $mergedBED based on molecule alignments in $alignmolXmap\n\n";
		mkpath("$inputs{output}/scoreBED") or warn "WARNING: Cannot create output directory $inputs{output}/scoreBED: $!\nNext steps may fail...\n";
		my $scoredBED2 = $inputs{output}."/scoreBED/$scoredBED";
		#Usage: perl scoreStitchPositionsBED_v4.pl <--bed stitchPositions.bed> <--xmap alignref_final.xmap> <--qcmap alignref_final_q.cmap> <--rcmap alignref_final_r.cmap> <--alignmolxmap alignmol.xmap> <--alignmolrcmap alignmol_r.cmap> <--fasta reference.fasta> <--key key file from fa2cmap> <--bam NGSalignments.bam to ref fasta> <--output output directory> <--ngsBuffer bp +/- stitch location to require NGS alignment support =1000> [<--n cpu cores>]
		$cmd = "perl $scriptspath/scoreStitchPositionsBED_v5.pl --bed $mergedBED --xmap $alignreffinalXmap --qcmap $alignreffinalQcmap --rcmap $alignreffinalRcmap --alignmolxmap $alignmolXmap --alignmolrcmap $alignmolRmap --n $cpuCount --output $inputs{output}/scoreBED/ --fasta $inputs{fasta} --key $inputs{key}";
		if ($inputs{bam}) {
			$cmd = $cmd." --bam $inputs{bam} --ngsBonus $ngsBonus";
		} 
		if ($inputs{ngsBuffer}) {
			$cmd = $cmd." --ngsBuffer $inputs{ngsBuffer}";
		}	
		$cmd = $cmd." > $inputs{output}/scoreBED/".basename($mergedBED,".bed")."_scored_log.txt 2>&1";
			
		print "\n\tRunning command: $cmd\n\n";
		system($cmd);
		print "\n"; 
		$scoredSTATS = $inputs{output}."/scoreBED/$scoredSTATS";
		copy("$scoredBED2","$inputs{output}/$scoredBED") or die "\nERROR: Copy of $scoredBED2 failed: $!\n";
		$scoredBED = $scoredBED2;
	}
	else {
		print "Skipping Step 11 BED scoring. Input BNX, reference, or fragileSiteRepaired merged CMAP not provided or does not exist! $!\n";
		print "\n";
	}
	}
	else {
		print "Skipping Step 11 BED scoring. --skipAlignMol enabled! $!\n";
		print "\n";
	}

	print "\n";


	# Step 12: cut fragileSiteRepaired CMAP where score is too low
	my $newfinalmap;
	my $newBED;
	print "=== Step 12: Break fragileSiteRepaired merged CMAP at stitchPositions with score less than $threshold ===\n\n";
	if (defined $inputs{break}) {	
		my $pre = basename($finalmap, ".cmap")."_final";
		#$pre = $pre."_final";
		$newBED = "$inputs{output}/scoreBED/".basename($scoredBED, ".bed");
		if (defined $inputs{bnx} && exists $inputs{bnx} && -e $finalmap && exists $inputs{ref}) {
			if (-e $scoredBED && -e $scoredSTATS) {
				# Usage: perl cutCmapByBedScore.pl <--bed stitchPositions_scored.bed> <--stats stitchPositions_scoreStats.csv> <--cmap fragileSiteRepaired_merged.cmap> <--prefix output prefix> [--threshold score threshold] [--maxfill maxfill <minlen filter>]
				$cmd = "perl $scriptspath/cutCmapByBedScore.pl --bed $scoredBED --stats $scoredSTATS --cmap $finalmap --prefix $pre --threshold $threshold --maxfill $inputs{maxfill} --output $inputs{output}/scoreBED/";
				if ($inputs{breakNGSonly}) {
					$cmd = $cmd." --breakNGSonly";
				}
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
				print "Skipping Step 12 Breaking fragileSiteRepaired CMAP based on scores. Scored BED and/or stats file not provided or does not exist! $!\n\n";
			}
		}
		else {
			print "Skipping Step 12 Breaking fragileSiteRepaired CMAP based on scores. Input BNX, reference, or fragileSiteRepaired merged CMAP not provided or does not exist! $!\n";
			print "\n";
		}
	}
	else {
		$newBED = $scoredBED;
		print "Disabled Step 12 Breaking fragileSiteRepaired CMAP based on scores. Use --break to enable\n";
		print "\n";
	}
	print "\n";


	# Step 13: run CharacterizeFinal on new final merged CMAP if needed
	if (-e $newfinalmap) {
		print "=== Step 13: Align final scored fragileSiteRepaired merged CMAP to reference and get stats ===\n\n";
		# move/rename old alignref_final folder
		$cmd = "mv $inputs{output}/alignref_final $inputs{output}/alignref";
		print "\tRunning command: $cmd\n\n";
		system($cmd);	
		
		if (-e $newfinalmap && (exists $inputs{ref})) {
			#copy("$inputs{ref}","$inputs{output}") or print "WARNING: Copy of input reference $inputs{ref} failed: $!\n"; 
			$cmd = "cp $scriptspath/runCharacterizeFinal.py $inputs{scripts}; python $inputs{scripts}/runCharacterizeFinal.py -t $inputs{RefAligner} -r $inputs{ref} -q $newfinalmap -p $inputs{scripts} -n $cpuCount -a $inputs{optArgs}";
			print "\tRunning command: $cmd\n\n";
			system($cmd);
			
			$alignreffinalXmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final.xmap";
			$alignreffinalQcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final_q.cmap";
			$alignreffinalRcmap = "$inputs{output}/alignref_final/$splitprefix"."_fragileSiteRepaired_merged_final_r.cmap";
		}
		else {
			print "Skipping Step 13 final characterizeFinal. Final merged fragileSiteRepaired CMAP or reference CMAP not found! $!\n";
			print "\n";
		}
	}
	else {
		print "Skipping Step 13 final characterizeFinal because it is not needed. Input BNX not provided and/or --break not enabled. $!\n\n";
	}
	print "\n";

	print "\n";
	# Step 14: run alignmol on new final merged CMAP if needed
	if (-e $newfinalmap) {
		print "=== Step 14: Run final single molecule alignments against final scored fragileSiteRepaired merged cmap ===\n\n";
		$alignmolXmap = "";
		#$alignmolXmap="";
		if (defined $inputs{bnx} && exists $inputs{bnx} && defined $inputs{err} && exists $inputs{err}) {
			mkpath("$inputs{output}/alignmol_final") or warn "WARNING: Cannot create output directory $inputs{output}/alignmol_final : $!\nNext steps may fail...\n"; 
			print "Running alignment of $inputs{bnx} to $newfinalmap\n\n";
			my $errHash_ref = extractErr($inputs{err});
			#my $errHash_ref = extractErr($alignmolErr);
			my %errHash = %$errHash_ref;	
			$alignmolDir = abs_path("$inputs{output}/alignmol_final");
			#$cmd = "cd $alignmolDir; ~/tools/RefAligner -f -ref $newfinalmap -o $splitprefix"."_fragileSiteRepaired_merged_final_alignmol -i $inputs{bnx} -maxthreads $cpuCount -maxmem ".($mem)." -usecolor 1 -FP 1.5 -FN 0.15 -sd 0. -sf 0.2 -sr 0.03 -res 3.3 -output-veto-filter intervals.txt\$ -T 1e-9 -usecolor 1 -S -1000 -biaswt 0 -res 3.3 -resSD 0.75 -outlier 0.0001 -extend 1 -BestRef 1 -maptype 1 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -stdout -stderr";
			$cmd = "cd $alignmolDir; ~$inputs{RefAligner} -i $inputs{bnx} -o $splitprefix"."_fragileSiteRepaired_merged_final_alignmol -ref $newfinalmap -f -stdout -stderr -maxthreads $cpuCount -usecolor 1 -FP $errHash{'FP(/100kb)'} -FN $errHash{'FN(rate)'} -sd $errHash{'sd'} -sf $errHash{'sf'} -sr $errHash{'sr'} -res $errHash{'res(pixels)'} -bpp $errHash{'bpp'} -output-veto-filter intervals.txt\$ -T 1e-10 -L 130 -nosplit 2 -biaswt 0 -resSD 0.75 -extend 1 -BestRef 1 -maptype 0 -PVres 2 -HSDrange 1.0 -hashoffset 1 -hashMultiMatch 20 -f -deltaX 4 -deltaY 6 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 1.4 -outlier 0.00001 -endoutlier $inputs{endoutlier} -outlierMax 40. -cres 5.4 3 0.1 -MaxSE 0.5 -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 3 -M 3 3 -MaxSE 0.5 -hash -hashdelta 10 -mres 0.9 -insertThreads 4 -rres 1.2 -maxmem $mem -readparameters $alignmolErrbin";
			
			print "\tRunning command: $cmd\n\n";
			system($cmd);
			sleep(5);
		
			#split xmap into _contig maps
			$alignmolXmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol.xmap";
			my $alignmolQmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_q.cmap";
			$alignmolRmap = "$alignmolDir/$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_r.cmap";
			my $alignmolPrefix = "$splitprefix"."_fragileSiteRepaired_merged_final_alignmol_contig";
			#mkpath("$alignmolDir/contigs") or die "ERROR: Cannot create output directory $inputs{output}/alignmol/contigs: $!\n";
			# Usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
			$cmd = "perl $scriptspath/split_xmap_standalone.pl $alignmolXmap $alignmolQmap $alignmolRmap $alignmolPrefix $alignmolDir/contigs";
			print "\n\tRunning command: $cmd\n\n";
			system($cmd);
			print "\n";
		}
		else {
			print "Skipping Step 14 Single molecule alignments. Input BNX not provided or does not exist. $!\n";
			print "\n";
		}
	}
	else {
		print "Skipping Step 14 final single molecule alignments. because it is not needed. Input BNX not provided and/or --break not enabled. $!\n\n";
	}
	print "\n\n";

	# Step 15: re-score final BED file if needed
	print "=== Step 15: Re-score final stitchPositions BED based on final single molecule alignments ===\n\n";
	if (-e $newfinalmap && -e $newBED) {
		
		my $scoredBED = basename($newBED,".bed");
		my $scoredSTATS = $scoredBED;
		$scoredBED = $scoredBED."_scored.bed";
		$scoredSTATS = $scoredSTATS."_scoreStats.csv";

		if (defined $inputs{bnx} && exists $inputs{bnx} && -e $newfinalmap && exists $inputs{ref} && exists $inputs{break}) {
			print "Scoring final BED file $newBED based on molecule alignments in $alignmolXmap\n\n";
			#mkpath("$inputs{output}/scoreBED") or warn "WARNING: Cannot create output directory $inputs{output}/scoreBED: $!\nNext steps may fail...\n";
			$scoredBED = $inputs{output}."/scoreBED/$scoredBED";
			my $oldBED = $newBED;
			$newBED = "$inputs{output}/".basename($scoredBED);
			
			#$cmd = "perl $scriptspath/scoreStitchPositionsBED_v5.pl --bed $oldBED --xmap $alignreffinalXmap --qcmap $alignreffinalQcmap --rcmap $alignreffinalRcmap --alignmolxmap $alignmolXmap --alignmolrcmap $alignmolRmap --n $cpuCount --output $inputs{output}/scoreBED/ > $inputs{output}/scoreBED/".basename($newBED,".bed")."_log.txt 2>&1";
			
			$cmd = "perl $scriptspath/scoreStitchPositionsBED_v5.pl --bed $oldBED --xmap $alignreffinalXmap --qcmap $alignreffinalQcmap --rcmap $alignreffinalRcmap --alignmolxmap $alignmolXmap --alignmolrcmap $alignmolRmap --n $cpuCount --output $inputs{output}/scoreBED/ --fasta $inputs{fasta} --key $inputs{key}";
			if ($inputs{bam}) {
				$cmd = $cmd." --bam $inputs{bam} --ngsBonus $ngsBonus";
			} 
			if ($inputs{ngsBuffer}) {
				$cmd = $cmd." --ngsBuffer $inputs{ngsBuffer}";
			}	
			$cmd = $cmd." > $inputs{output}/scoreBED/".basename($newBED,".bed")."_log.txt 2>&1";		
			
			print "\n\tRunning command: $cmd\n\n";
			system($cmd);
			print "\n"; 
			$scoredSTATS = $inputs{output}."/scoreBED/$scoredSTATS";
			#$scoredBED = $inputs{output}."/scoreBED/$scoredBED";
			
			copy("$scoredBED","$inputs{output}") or print "WARNING: Copy of re-scored BED $scoredBED failed: $!\n"; 	
			
		}
		else {
			print "Skipping Step 15 re-score BED \. Input BNX, reference, or final fragileSiteRepaired merged CMAP not provided or does not exist! $!\n";
			print "\n";
		}
	}
	else {
		print "Skipping Step 15 re-score BED because it is not needed. Input BNX not provided and/or --break not enabled. $!\n\n";
	}
	print "\n\n";


	# Step 16: Calculate stats for merged fragileSiteRepaired BED and CMAP
	print "=== Step 16: Calculate stats for final merged fragileSiteRepaired CMAP ===\n";
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
	my ($origBedOut_ref, $origTypeCount_ref, $origTypeConf_ref, $origtypeConfStdDev_ref) = sortBEDandCount(\@origBed);
	my %origTypeCount = %$origTypeCount_ref;
	my %origTypeConf = %$origTypeConf_ref;
	my %origTypeConfStdDev = %$origtypeConfStdDev_ref;
	my $origCount = scalar(@origBed);
	print "Total potential predicted fragile sites: $origCount\n";
	foreach my $type (sort keys %origTypeCount) {
		print "\t$type: $origTypeCount{$type} AvgConf: $origTypeConf{$type} StdDevConf: $origTypeConfStdDev{$type}\n";
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
		my ($newBedOut_ref, $newTypeCount_ref, $newTypeConf_ref, $newTypeConfStdDev_ref) = sortBEDandCount(\@newBed);
		my %newTypeCount = %$newTypeCount_ref;
		my %newTypeConf = %$newTypeConf_ref;
		my %newTypeConfStdDev = %$newTypeConfStdDev_ref;
		my $newCount = scalar(@newBed);
		
		print "Total fragile sites repaired (stitched) and scored: $newCount\n";
		foreach my $type (sort keys %newTypeCount) {
			print "\t$type: $newTypeCount{$type} AvgConf: $newTypeConf{$type} StdDevConf: $newTypeConfStdDev{$type}\n";
		}
	}
	else {
		print "Total fragile sites repaired (stitched): $stitchCount\n";
		foreach my $type (sort keys %typeCount) {
			print "\t$type: $typeCount{$type} AvgConf: $typeConf{$type} StdDevConf: $typeConfStdDev{$type}\n";
		}
	}
	print "\n";

	# usage: calc_cmap_stats.pl <CMAP_File>
	my $dir = glob("$inputs{scripts}/HybridScaffold/scripts");
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
		print "WARNING: Skipping Step 15: perl script calc_cmap_stats.pl not found at $dir/$script\n\n"; }
	print "\n\n";


	# Step 17: run SV detection on final cmap
	print "=== Step 17: Run SV detection on fragileSiteRepaired CMAP ===\n\n";
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
				$cmd = "perl $scriptspath/callSV.pl --ref $inputs{ref} --cmap $map --errbin $errbin --optarg $inputs{optArgs} --outdir $inputs{output}/alignref_final_sv --mres $mres --n $cpuCount --prefix $splitprefix"."_fragileSiteRepaired_merged_final";
				if ($inputs{force} eq 1) {$cmd = $cmd." --force";}
				if (-e $inputs{gaps}) {$cmd = $cmd." --gaps $inputs{gaps}";}
				print "\tRunning command: $cmd\n\n";
				system($cmd);
				
				#print SV stats summary
				print "==SV Summary ==\n\n";
				#$cmd = "grep -A 25 -i \"SV type\" $inputs{output}/alignref_final_sv/sv_log.txt";
				$cmd = "cat $inputs{output}/alignref_final_sv/sv_log.txt | sed -n -e \'/SV type/,\$p\'";
				my @SVresults = capture($cmd); 
				print join("",@SVresults);
			}
			else { print "WARNING: Merged CMAP $map and/or alignref_final errbin $errbin not found! Skipping Step 16 runSV...\n"; }
		}
		else {
			print "Skipping Step 17. Original assembly CMAP and/or reference CMAP not provided\n";
			print "\n";
		}
	}
	else {
		print "Disabled Step 17 SV detection. Add --runSV to enable SV detection\n";
		print "\n";
	}
	print "\n";



} # end if --MapsMode
elsif ($inputs{MoleculesMode}) {

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
	my ($origBedOut_ref, $origTypeCount_ref, $origTypeConf_ref, $origtypeConfStdDev_ref) = sortBEDandCount(\@origBed);
	my %origTypeCount = %$origTypeCount_ref;
	my %origTypeConf = %$origTypeConf_ref;
	my %origTypeConfStdDev = %$origtypeConfStdDev_ref;
	my $origCount = scalar(@origBed);
	print " ===== Total potential predicted fragile sites: $origCount ===== \n";
	foreach my $type (sort keys %origTypeCount) {
		print "\t$type: $origTypeCount{$type} AvgConf: $origTypeConf{$type} StdDevConf: $origTypeConfStdDev{$type}\n";
	}
	print "\n\n";

	#get stats of stitchPositions BED
	print "===== Total fragile sites repaired (stitched): $stitchCount =====\n";
	foreach my $type (sort keys %typeCount) {
		print "\t$type: $typeCount{$type} AvgConf: $typeConf{$type} StdDevConf: $typeConfStdDev{$type}\n";
	}
	print "\n";
	print "\n";
	print " ===== Create merged fragileSiteRepaired BNX ===== \n";
	$cmd = "cd $inputs{output}; $inputs{RefAligner} -merge -bnx -minsites 0 -minlen 0 -i $finalmap -o $splitprefix"."_fragileSiteRepaired_merged -stdout -stderr";
	print "Running command: $cmd\n";
	system($cmd);
	print "\n\n";

	$cmd = "cd $inputs{output}; $scriptspath/convertBNX0.1_to_1.0.pl $splitprefix"."_fragileSiteRepaired_merged.bnx";
	system($cmd);
	print "\n\n";
	
	
	my $dir = glob("$inputs{scripts}/HybridScaffold/scripts");
	my $script = "calc_cmap_stats.pl";
	$cmd = "";
	if (-e "$dir/$script") {
		if (exists $inputs{cmap}) {
			print " ===== Original molecules BNX $inputs{cmap} stats ===== \n";
			$cmd = "perl $dir/$script $inputs{cmap}";
		}
		#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";
		# my $file = "$inputs{output}/$splitprefix"."_fragileSiteRepaired.cmap";
		#print "Fragile site repaired merged CMAP $finalmap stats:\n";
		#chdir $dir or die "ERROR: Cannot change directory to $dir: $!\n";	
		if (-e $finalmap) {
			print " ===== Fragile site repaired merged BNX $finalmap stats =====\n";
			$cmd = "perl $dir/$script $finalmap";
		}
		print "Running command: $cmd\n";
		system($cmd);
		print "\n";
	
	
		}
	}

























copy("$log_file","$inputs{output}");
copy("$log_err_file","$inputs{output}");








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
	my %typeConf;
	my %confHashArray;
	my %typeAvgConf;
	my %typeStdDevConf;
	foreach my $line (@bedIn) {
		chomp($line);
		#print "$line\n";
		my @s = split("\t",$line);
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		if (scalar(@s) == 10) {
			push (@AoA, [$s[0],int($s[1]),int($s[2]),$s[3], $s[4],$s[5],int($s[6]),int($s[7]),$s[8],$s[9]]);
		}
		else {
			push (@AoA, [$s[0],int($s[1]),int($s[2]),$s[3],$s[4],$s[5],int($s[6]),int($s[7]),$s[8]]);
		}
		if ($s[3] =~ m/Type/i) {
			$typeCount{$s[3]}++;
		}
		elsif ($s[3] =~ m/Observed/i) {
			$typeCount{$s[3]}++;
		}
		if (!exists $typeConf{$s[3]}) {
			$typeConf{$s[3]} = $s[4];
			push( @{$confHashArray{$s[3]} }, $s[4]);  
			
		}
		else {
			$typeConf{$s[3]} += $s[4];
			push( @{$confHashArray{$s[3]} }, $s[4]);  
		}
	}
	foreach my $type (sort keys %typeConf) {
		$typeAvgConf{$type} = $typeConf{$type} / $typeCount{$type};
		$typeAvgConf{$type} = (sprintf "%.2f", $typeAvgConf{$type});
		
		my $confArray_ref = $confHashArray{$type};
		my @confArray = @{$confArray_ref};
		$typeStdDevConf{$type} = stddev(@confArray);
		$typeStdDevConf{$type} = (sprintf "%.2f", $typeStdDevConf{$type});

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
	return (\@bedOut,\%typeCount,\%typeAvgConf, \%typeStdDevConf);
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

sub findCMAPs {
	my $in = shift;
	my $str = shift;
	opendir(DIR, $in);
	my @cmaps = grep(/$str/,readdir(DIR));
	closedir(DIR);
	
	return @cmaps;
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

sub wanted4 { 
	m/.bnx$/g and do { 
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


sub extractErr {
	my $fileIn = shift;
	my @errArr;
	open FILE, "<$fileIn" or die "ERROR: Could not open $fileIn: $!\n";
	while (<FILE>) {
		if ($_ =~ /^#/) {
			next;
		}
		chomp($_);
		push @errArr, $_;
	}
	close FILE;
	
	my %errHash;
	#Iteration	FP(/100kb)	FN(rate)	sf	sd	bpp	res(pixels)	Maps	Log10LR(/Maps)	GoodMaps	log10LR(/GoodMaps)	bppSD	FPrate	sr	se	LabelDensity(/100kb)	resSD(pixels)	mres(x500bp)	mresSD(x500bp)
	#OLD: Iteration	FP(/100kb)	FNrate	SiteSD(Kb)	ScalingSD(Kb^1/2)	bpp	res(pixels)	Maps	Log10LR(/Maps)	GoodMaps	log10LR(/GoodMaps)	bppSD	FPrate	RelativeSD	ResolutionSD	LabelDensity(/100kb)	resSD	mres	mresSD
	my @header = split("\t",$errArr[0]);
	my @params = split("\t",$errArr[-1]);
	for (my $i=0; $i<scalar(@params); $i++) {
		my $head = $header[$i];
		my $param = $params[$i];
		if (defined $head) { chomp($head); } 
		if (defined $param) { chomp($param); }
		if ($head =~ /FNrate/i) { $head = "FN(rate)"; }
		elsif ($head =~ m/SiteSD/i) { $head = "sf"; }
		elsif ($head =~ m/ScalingSD/i) { $head = "sd"; }
		elsif ($head =~ m/RelativeSD/i) { $head = "sr"; }
		#load into hash
		$errHash{$head} = $param;
	}
	
	return \%errHash;
} 
	
sub Usage {
	#Usage: perl fragileSiteRepair.pl --fasta <reference.fasta> --cmap <assembly.cmap> --output <output folder> [--bnx <input.bnx to run alignmol>] [--err <original assembly alignmol_merge.err or autonoise1.err>] [--enzyme <sequence of enzyme to use =GCTCTTC>] [--bam <.bam alignment of NGS reads/contigs to input ref FASTA>] [--ngsBuffer <bp +/- stitchPosition to require NGS alignment =500>] [--threshold <minimum score threshold below which to cut stitched maps =1>] [--maxlab <max_label_gap_tolerence =1>] [--maxfill <max basepairs to fill between contigs =30000>] [--wobble <fragile site wobble in bp =30000>] [--seq <sequence bp +/- fragile site to print in final BED =off>] [--force <overwrite output folder =off>] [--n <number of CPU cores =nproc>] [--j <number of parallel jobs =nproc/6>] [--optArgs <optArguments.xml =optArguments_human.xml>] [--runSV <enable run SV detection =off>]  [--break <break maps at stitch locations with score below threshold =off>] [--endoutlier <Pvalue penalty for endoutlier =1e-3>] [--aggressive <calculate TypeIII and TypeIV fragile sites =off>] [--ngsBonus <raw score bonus for supporting NGS alignments =10>] [--breakNGSonly <break maps if stitch has no BioNano support default=off>]
	print "Usage: perl fragileSiteRepair.pl [OPTIONS] --fasta <ref.fa> --cmap <assembly.cmap> --output <output folder>\n\n";
	print "\t=== REQUIRED ===\n";
	print "\t  OPERATION MODE\n";
	print "\t--MapsMode : Operate in \"MapsMode\" where primary input is BioNano assembly CMAP and reference/NGS fasta. If NGS fasta used, recommended optArguments_NGS.xml\n";
	print "\t--MoleculesMode : Operate in \"MoleculesMode\" where primary input is BioNano single molecules and reference/NGS fasta. Recommended to use filtered BNX to be used in assembly\n";
	print "\n\t  INPUTS\n";
	print "\t--fasta <.fasta> : input reference FASTA file to be used for fragile site prediction and alignment. REQUIRED\n";
	print "\t--output <path/to/output> : output folder. REQUIRED\n";
	print "\t--cmap <.cmap> : input BioNano assembly CMAP to perform fragileSiteRepair on. REQUIRED with --MapsMode\n";
	print "\t--bnx <.bnx> : input BioNano single molecules BNX. REQUIRED with --MoleculesMode. OPTIONAL for --MapsMode. In --MapsMode, BNX used to perform single molecule alignments and stitch scoring\n";	
	print "\n";
	print "\t=== KEY OPTIONS ===\n";
	
	print "\t--RefAligner <path/to/RefAligner> : Path to RefAligner binary to use. Default: ./tools/RefAligner\n";
	print "\t--scripts <path/to/scripts> : Path to IrysSolve pipeline scripts folder. Default: ./scripts\n";
	print "\t--alignmolDir <path/to/assembly/output/contigs/exp_refineFinal1/alignmol/merge> : directory containing molecules aligned to assembly maps XMAPs and ERR file. Used to stitch across \"observed\" fragile sites. Requires --bnx to perform single molecule alignment scoring. Default: NONE\n";
	print "\t--bam <.bam> : BAM file of NGS scaffolds/contigs/reads aligned to --fasta. Used to supplement scoring of stitchPositions. Default: OFF\n";
	print "\n";
	print "\t--break : Flag to enable breaking of maps at stitchPositions with score less than threshold value. Requires --bnx and --err. Default: OFF\n";
	print "\t--runSV : Flag to run SV module on final cmap. Default: OFF\n";
	print "\t--gaps  : N-base gaps BED file\n";
	print "\t--force : Flag to overwrite --output folder. Default: OFF\n";
	print "\t--skipAlignMol : Flag to enable skipping single molecule alignments and scoring Default: OFF\n";
	print "\t--aggressive : Flag to calculate TypeIII and TypeIV fragile sites in addition. Default: OFF\n";
	print "\n";
	print "\t=== OTHER OPTIONS ===\n";
	print "\t--enzyme <nickase sequence> : Nickase enzyme for in-silico digestion and fragile site prediction. Default: GCTCTTC\n";
	print "\t--ngsBuffer <basepairs> : Number of basepairs that a single NGS alignment must extend past the ends of a stitchPosition to supplement score. Default: 500\n";
	print "\t--ngsBonus <raw score value> : Raw score bonus for each NGS alignment supporting a fragile site repair. Default: 10\n";
	print "\t--breakNGSonly : Flag to break maps at stitchPositions that have only NGS alignment support and no BioNano single molecule support. Default: OFF\n";
	print "\t--threshold <scaled score> : Minimum stitchPoisitons scaled score below which to break maps. Default: 1.0\n";
	print "\t--maxlab <label count> : Maximum number of reference labels to allow between adjacent maps. Default: 1\n";
	print "\t--maxfill <basepairs> : Maximum number of basepairs to allow between adjacent maps. Default: 30000\n";
	print "\t--wobble <basepairs> : Maximum number of basepairs to allow the fragile site to vary. Default: 30000\n";
	print "\t--seq <basepairs> : Number of basepairs of reference sequence +/- fragile site to output into BED file. Default: OFF\n";
	print "\t--optArgs <optArguments.xml> : optArguments.xml to use for alignment. Default: ./scripts/optArguments_human.xml\n";
	print "\t--endoutlier <pvalue> : endoutlier penalty for single molecule alignments. Default: 1e-5\n";
	print "\t--minRatio <ratio> : minimum ratio of single molecule alignments that must end at genome maps ends to be classified as a potential fragile site. Requires --alignmolDir. Default: 0.70\n";
	print "\t--maxOverlap <bp> : maximum number of basepairs overlap between maps to allow merge. Default: 20000\n";
	print "\t--maxOverlapLabels <int labels> : maximum number of labels overlap between maps to allow merge. Default: 5\n";
	print "\t--pvalue <pvalue> : Minimum threshold pvalue used for alignments. Recommended to match pvalue from optArguments characterizeFinal. See RefAligner help for more info. Default: 1e-12\n";
	print "\t--maxiter <int> : Maximum number of rounds to perform fragile site repair Default: 20\n";
	print "\n";
	print "\t--n <CPU cores> : Maximum number of CPU cores/threads to use. Default: nproc\n";
	print "\t--j <number jobs> : Maximum number of parallel jobs. Default: nproc/6\n";
	print "\n";
	print "\t--h : Display this help menu\n";	

	print "\n";
}

sub die {
	my $txt = shift;
	print "$txt\n";
	exit 1;
}
	
