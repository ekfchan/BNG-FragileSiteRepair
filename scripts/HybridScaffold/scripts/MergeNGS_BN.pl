# $Id: MergeNGS_BN.pl 4290 2015-11-20 00:13:03Z apang $

#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Getopt::Long;


# This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
# This script sould sit one level above the additional Perl modules directory.
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path2 = abs_path($script_path . "/perl5");
	unshift @INC, $module_path2;
	my $lib4;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.016000 && $] <= 5.017000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.018000 && $] <= 5.019000) {
		$module_path2 = $module_path2."/5.18.2";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n"; 
		exit; }
	unshift @INC, $module_path2;
	unshift @INC, $lib4;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

print "\nInfo: Running the command: $0 @ARGV\n";

use File::Path qw(make_path);
use File::Copy;
use BNG::Utility;
use BNG::refAlignerRun;

# declare and read in command line args
my $outputDir = "";
my $refaligner = "";
my $merge_Tvalue = 1e-13;
my $scratchDir = "";
my $logFile = "";
my $ngs_cmap_fn = "";
my $bng_cmap_fn = "";
my $id_shift=100000;
my $max_merge_rounds=40;			# probably useless now
my $maxthreads=128;
my $endoutlier=0;
my $outlier=1e-4;
my $biaswt=0;
my $sd=0.1;
my $res=2.9;
my $sf=0.2;
my @RepeatMaskArray = ();       my $RepeatMask="";
my @RepeatRecArray = ();        my $RepeatRec="";
my @pairmergeArray = ();                my $pairmerge="";
my $readparameters = "";
my $minsites = 0;
my $maxmem = 128;
my $unMergeableMaxSites = 4;
my $mergeableMinSites = $unMergeableMaxSites + 1;

GetOptions      (
	"outputDir=s"		=>	\$outputDir,
	"refaligner=s"		=>	\$refaligner,
	"merge_Tvalue=f"	=>	\$merge_Tvalue,
	"scratchDir=s"		=>	\$scratchDir,
	"logFile=s"		=>	\$logFile,
	"ngs_cmap_fn=s"		=>	\$ngs_cmap_fn,
	"bng_cmap_fn=s"		=>	\$bng_cmap_fn,
	"id_shift=i"		=>	\$id_shift,
	"max_merge_rounds=i"	=>	\$max_merge_rounds,
	"maxthreads=i"		=>	\$maxthreads,
	"endoutlier=f"		=>	\$endoutlier,
	"outlier=f"		=>	\$outlier,
	"biaswt=i"		=>	\$biaswt,
	"sd=f"			=>	\$sd,
	"res=f"			=>	\$res,
	"sf=f"			=>	\$sf,
	"RepeatMask=f{,}"	=>	\@RepeatMaskArray,
	"RepeatRec=f{,}"	=>	\@RepeatRecArray,
	"pairmerge=f{,}"	=>	\@pairmergeArray,
	"readparameters=s"	=>	\$readparameters,
	"minsites=i"		=>	\$minsites,
	"maxmem=i"		=>	\$maxmem) or dieLog ("ERROR: MergeNGS_BN: error in command line arguments.\n");

$RepeatMask = join(" ", @RepeatMaskArray);      
$RepeatRec = join(" ", @RepeatRecArray);
$pairmerge = join(" ", @pairmergeArray);

#####################
# Step 0. Initiation.
#####################
chdir $outputDir;		# change current directory to $outputDir
print "cwd=$outputDir\n";
make_path($scratchDir); 	# make sure $scratchDir exist
open(LOG, ">$logFile");
my $logFH = *LOG;

# 0.1. Read in the BioNanoAssembled Cmap files
my ($bng_cmap, $nc_bng_cmap, undef) = readCMap($bng_cmap_fn);
logMessage($logFH, "Read BioNano contig '".$bng_cmap_fn."' completed with ", $nc_bng_cmap, " cmaps.");

# 0.2 Shift BioNano ContigId by offset to prevent ContigID collision,
shiftCMapIds($bng_cmap, $id_shift);
my $bionano_idshift_cmap_fn = $scratchDir."bionano.idshift.cmap";
writeCMapFile($bng_cmap, $bionano_idshift_cmap_fn, 0);
logMessage($logFH, "Output BioNano contig with ID shift of '" . $id_shift . "' completed and output to '" . $bionano_idshift_cmap_fn, "'.");
logMessage($logFH, "Step 0. Initiation Completed");


# 0.3 Separate the sequence entries into those with too few sites and those with normal number of sites
my $runContigExtract =
	new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		i			=>	$ngs_cmap_fn,
		o			=>	$scratchDir."ngs_contig_minsites_$minsites"."_maxsites_$unMergeableMaxSites",
		merge			=>	1,
		minsites		=>	$minsites,
		maxsites		=>	$unMergeableMaxSites,
		f			=>	1,
		stdout			=>	1,
		stderr			=>	1
});
logMessage($logFH, $runContigExtract->getCMD());
my ($outResults1, $errResults1, $jobStatus1) = $runContigExtract->runCMD();
exit $jobStatus1 if ($jobStatus1 != 0);

$runContigExtract =
	new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		i			=>	$ngs_cmap_fn,
		o			=>	$scratchDir."ngs_contig_minsites_$mergeableMinSites",
		merge			=>	1,
		minsites		=>	$mergeableMinSites,
		f			=>	1,
		stdout			=>	1,
		stderr			=>	1
});
logMessage($logFH, $runContigExtract->getCMD());
($outResults1, $errResults1, $jobStatus1) = $runContigExtract->runCMD();
exit $jobStatus1 if ($jobStatus1 != 0);
my $mergeable_ngs_cmap_fn = $scratchDir . "ngs_contig_minsites_$mergeableMinSites.cmap";        # only look at those NGS that are potentially mergeable

my @mrg_rounds_ids = ("A".."Z", "AA".."AZ", "BA".."BZ", "CA".."CZ", "DA".."DZ");

#####################
# Step 1. Merge
#####################
logMessage($logFH, "Step 1. Pair Merge Repeat between NGS\/BioNano started.");
my $mrg_output_prefix = $scratchDir."Mrg";
my $pairmerge_single_round_2files_run = 
	new BNG::refAlignerRun({	
		binPath			=>	$refaligner,
		i			=>	[$bionano_idshift_cmap_fn, $mergeable_ngs_cmap_fn],
		o			=>	$mrg_output_prefix,
		T			=> 	$merge_Tvalue,
		preferNGS		=>	1,
		pairmergeRepeat 	=>	1,
		endoutlier		=>	$endoutlier,
		outlier			=>	$outlier,
		biaswt			=>	$biaswt,
		maxthreads		=>	$maxthreads,
		sd			=>	$sd,                 #noise.params
		res			=>	$res,
		sf			=>	$sf,
		f			=>	1,                   # overwrite
		stdout			=>	1,
		stderr			=>	1,
		first			=>	-1,
		RepeatMask		=>	"$RepeatMask",
		RepeatRec		=>	"$RepeatRec",          
		pairmerge		=>	"$pairmerge",
		readparameters		=>	"$readparameters",
		minsites		=>	$minsites,
		maxmem			=>	$maxmem,
		NoBpp			=>	1 
});
my $cmd0 = $pairmerge_single_round_2files_run->getCMD();
logMessage($logFH, $cmd0);
my ($outResults0, $errResults0, $job_status0) = $pairmerge_single_round_2files_run->runCMD();
exit $job_status0 if ($job_status0 != 0);

my ($numMrg, $this_mrg_pairs) = parsingFastMrgStdOut("$mrg_output_prefix.stdout");

if ($numMrg == 0)	{
	logMessage($logFH, "1.1 No merge was possible. Stop the merge program.");
	logMessage($logFH, "ERROR: 1.1. There was no merge possible between NGS and BioNano map.");
	dieLog("ERROR: 1.1. There are no merge possible between NGS and BioNano map. Please loosen parameters (e.g. -T or -pairmerge) and check input files.");
} # if numMrg

# keep record of which pairs of fragments were merged
writeAllFastMrgPairs($scratchDir."step1.merge.pairs.txt", $this_mrg_pairs, $numMrg, \@mrg_rounds_ids);
logMessage($logFH, "1. We successfully merged $numMrg pairs.");

# now figure out which contigs are hybrid scaffolds, and sequences, and BioNano
my ($usedBioNanoRef, $usedSeqRef, $hybridRef) = determineParticipants($this_mrg_pairs, $numMrg, $id_shift);

# now iterates the output directory to figure out which file is hybrid, sequence leftover or bionano leftover
my ($hybridCmapsRef, $seqLeftOverCmapsRef, $bioNanoLeftOverCmapsRef) = categorizeCmapFiles($scratchDir, "Mrg_contig", $id_shift, $usedBioNanoRef, $usedSeqRef, $hybridRef);
my ($numHybridFiles, $numSeqLeftOverFiles, $numBioNanoLeftOverFiles) = (scalar(keys %$hybridCmapsRef), scalar(keys %$seqLeftOverCmapsRef), scalar(keys %$bioNanoLeftOverCmapsRef));
logMessage($logFH, "1.1 There are $numHybridFiles hybrid scaffold cmap files, $numSeqLeftOverFiles sequence left over cmap files, $numBioNanoLeftOverFiles BioNano left over cmap files");

# now make a single hybrid cmap, a single sequence left over and a single BioNano left over cmap files
logMessage($logFH, "1.2 Concatenating into single hybrid, and single left over cmap files now");
# print a test file
printIdFile($scratchDir."step1.hybrid.id.txt", $hybridCmapsRef);
printIdFile($scratchDir."step1.BN.id.txt", $bioNanoLeftOverCmapsRef);
printIdFile($scratchDir."step1.NGS.id.txt", $seqLeftOverCmapsRef);

my $quickConcatCmap = new BNG::refAlignerRun({
	binPath			=>	$refaligner,
	"if"			=>	$scratchDir."step1.hybrid.id.txt",
	o			=>	$scratchDir."step1.hybrid",
	merge			=>	1,
	stdout			=>	1,
	f			=>	1
});
my $cmd = $quickConcatCmap->getCMD();
logMessage($logFH, $cmd);
my ($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
exit $job_status if ($job_status != 0);

$quickConcatCmap = new BNG::refAlignerRun({
	binPath 		=>	$refaligner,
	"if"			=>	$scratchDir."step1.hybrid.id.txt",
	o			=>	$scratchDir."step2.hybrid",
	merge			=>	1,
	stdout			=>	1,
	f			=>	1
});
$cmd = $quickConcatCmap->getCMD();
logMessage($logFH, $cmd);
($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
exit $job_status if ($job_status != 0);

if ($numBioNanoLeftOverFiles != 0)	{
	$quickConcatCmap = new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		"if"			=>	$scratchDir."step1.BN.id.txt",
		o			=>	$scratchDir."step1.BN.naive",
		merge			=>	1,
		stdout			=>	1,
		f			=>	1
	});
	$cmd = $quickConcatCmap->getCMD();
	logMessage($logFH, $cmd);
	($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
	exit $job_status if ($job_status != 0);
} else	{
	# no genome map left over, just create an empty naive file
	open(OUT, ">$scratchDir"."step1.BN.naive.cmap") or dieLog ("ERROR: Cannot write to $scratchDir"."step1.BN.naive.cmap: $!\n");
	close OUT;
} # if numBioNanoLeftOverFiles

if ($numSeqLeftOverFiles != 0)	{
	$quickConcatCmap = new BNG::refAlignerRun({
		binPath			=>	$refaligner,
		"if"			=>	$scratchDir."step1.NGS.id.txt",
		o			=>	$scratchDir."step1.NGS.naive",
		merge			=>	1,
		stdout			=>	1,
		f			=>	1
	});
	$cmd = $quickConcatCmap->getCMD();
	logMessage($logFH, $cmd);
	($outResults, $errResults, $job_status) = $quickConcatCmap->runCMD();
	exit $job_status if ($job_status != 0);
} else	{
	# no sequence left over, just create an empty naive file
	open(OUT, ">$scratchDir"."step1.NGS.naive.cmap") or dieLog ("ERROR: Cannot write to $scratchDir"."step1.NGS.naive.cmap: $!\n");
	close OUT;	
} # if numSeqLeftOverFiles

logMessage($logFH, "Step 1.2 file concatenation completed successfully");

close($logFH);
print "$0 finished successfully.";

sub printIdFile	{
	my ($file, $idsRef) = @_;
	open(OUT, ">$file") or dieLog("ERROR: printIdFile: cannot write to $file: $!\n");
	foreach my $id (sort numeric keys %$idsRef)	{
		print OUT "$idsRef->{$id}\n";	# print the file name
	} # foreach id
	close OUT;
} # printIdFile

sub numeric	{	$a	<=>	$b	}

sub categorizeCmapFiles	{
	my ($dir, $filePrefix, $idShift, $usedBioNanoRef, $usedSeqRef, $hybridRef) = @_;
	my %hybridCmaps = ();
	my %seqLeftOverCmaps = ();
	my %bioNanoLeftOverCmaps = ();
	opendir(DIR, $dir) or dieLog "ERROR: categorizeCmapFiles: cannot open dir $dir: $!\n";
	my @cmapFiles = grep {$_ =~ /^$filePrefix\d+\.cmap$/i} readdir DIR;
	closedir DIR;

	for (my $i = 0; $i < scalar(@cmapFiles); $i += 1)	{
		my $theId = $cmapFiles[$i];	$theId =~ s/^$filePrefix//;	$theId =~ s/\.cmap$//;
		if (exists $hybridRef->{$theId})	{
			# this file is a hybrid
			$hybridCmaps{$theId} = $dir."$cmapFiles[$i]";
		} else	{
			# this file is a left over, check whether it is a sequence or a BioNano
			if ($theId < $idShift)	{
				$seqLeftOverCmaps{$theId} = $dir."$cmapFiles[$i]";
			} else	{
				$bioNanoLeftOverCmaps{$theId} = $dir."$cmapFiles[$i]";
			} # if theId
		} # if exists
	} # for i	
	return (\%hybridCmaps, \%seqLeftOverCmaps, \%bioNanoLeftOverCmaps);
} # categorizeCmapFiles

sub determineParticipants	{
	my ($mrgPairsRef, $numMrg, $idShift) = @_;
	my %usedBioNano = ();
	my %usedSeq = ();
	my %hybrid = ();	# stores hybrid ids, but would not store those that are not present in the output directory in the end, BECAUSE they have been merged with an entity either with a lower id or was totally encompassed

	for (my $i = 0; $i < $numMrg; $i += 1)	{
		my ($contig1, $contig2, $theHybrid) = ($mrgPairsRef->{ContigID1}[$i], $mrgPairsRef->{ContigID2}[$i], $mrgPairsRef->{ResultContigID}[$i]);
		if ($contig1 < $idShift)	{
			# a sequence
			$usedSeq{$contig1} = 1;
		} else	{
			# a genome map
			$usedBioNano{$contig1} = 1;
		} # if mrgPairsRef

		if ($contig2 < $idShift)	{
			# a sequence
			$usedSeq{$contig2} = 1;
		} else	{
			# a genome map
			$usedBioNano{$contig2} = 1;
		} # if mrgPairsRef

		$hybrid{$theHybrid} = 1;

		# now figure out if the particpant was a hybrid already AND that this merge resulted in an hybrid with a different id (either be a lower id, or that the participant was completely encompassed)
		delete $hybrid{$contig1} if (exists $hybrid{$contig1} && $theHybrid != $contig1);
		delete $hybrid{$contig2} if (exists $hybrid{$contig2} && $theHybrid != $contig2);
	} # for i
	return (\%usedBioNano, \%usedSeq, \%hybrid);
} # determineParticipants

