# $Id: align_molecules.pl 4311 2015-11-30 21:30:32Z apang $
#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use IPC::Open3;
use IO::Select;

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
use BNG::Utility;
print "\nInfo: Running the command: $0 @ARGV\n";

# global section
my $refAligner = "";
my ($globalMaxThreads, $numThreads) = (32, 32);
# common files between auto noise and align mol
my ($moleculesDir, $moleculesFile) = ("", "");
my ($bioNanoMapDir, $bioNanoMapFile) = ("", "");
my ($hybridMapDir, $hybridMapFile) = ("", "");

# auto noise part
my $autoNoisePythonScript = "pipelineCL.py";
my $autoNoiseOutDir = "autoNoise";	# should be passed in by hybridScaffold.pl (needs to know which _M\d+ subdirectory to output to)
my $deNovoOpArgsFile = "optArguments.xml";

# align molecule part
my $alignMolPythonScript = "runAlignMol.py";
my $deNovoPipelineDir = "Pipeline";
# should be passed in by hybridScaffold.pl (needs to know which _M\d+ subdirectory to output to)
my $alignMolBioNanoOutDir = "alignmol_bionano";	
my $alignMolHybridOutDir = "alignmol_hybrid";	 


# read in arguments section
GetOptions	(
	"refAligner=s"			=>	\$refAligner,
	"globalMaxThreads=i"		=>	\$globalMaxThreads,
	"deNovoPipelineDir=s"		=>	\$deNovoPipelineDir,
	"deNovoOpArgsFile=s"		=>	\$deNovoOpArgsFile,
	"moleculesDir=s"		=>	\$moleculesDir,
	"moleculesFile=s"		=>	\$moleculesFile,
	"bioNanoMapDir=s"		=>	\$bioNanoMapDir,
	"bioNanoMapFile=s"		=>	\$bioNanoMapFile,
	"hybridMapDir=s"		=>	\$hybridMapDir,
	"hybridMapFile=s"		=>	\$hybridMapFile,
	"autoNoiseOutDir=s"		=>	\$autoNoiseOutDir,
	"alignMolBioNanoOutDir=s"	=>	\$alignMolBioNanoOutDir,
	"alignMolHybridOutDir=s"	=>	\$alignMolHybridOutDir
) or dieLog("ERROR: align_molecules.pl: error in command line arguments\n");

my $refAlignerDir = $refAligner;	$refAlignerDir =~ s/RefAligner$//;	# just give a directory
$numThreads = $globalMaxThreads;	# since we are just running one job, then numThreads and maxThreads per job should be the same
$autoNoisePythonScript = "$deNovoPipelineDir/pipelineCL.py";
$alignMolPythonScript = "$deNovoPipelineDir/runAlignMol.py";
mkdir $autoNoiseOutDir if (! -e $autoNoiseOutDir);
mkdir $alignMolBioNanoOutDir if (! -e $alignMolBioNanoOutDir);
mkdir $alignMolHybridOutDir if (! -e $alignMolHybridOutDir);
my $bioNanoMapFilePrefix = $bioNanoMapFile;	$bioNanoMapFilePrefix =~ s/\.cmap$//;
my $hybridMapFilePrefix = $hybridMapFile;	$hybridMapFilePrefix =~ s/\.cmap$//;

my @cmd = ();	my $cmdRef = \@cmd;
my ($outResults, $errResults) = ("", "");
my $jobStatus = 0;

############### auto noise ######################
# check if auto noise results are already present
# if yes, then skip this step, else perform autonoise
print "\nBeginning alignment molecules to BioNano genome maps (auto noise)\n";
dieLog ("ERROR: align_molecules.pl: cannot perform auto noise because $autoNoisePythonScript does not exist\n") if (! -e $autoNoisePythonScript);
@$cmdRef = ("/usr/bin/python", $autoNoisePythonScript, "-x", "-y", "-T", $numThreads, "-j", $globalMaxThreads, "-l", $autoNoiseOutDir, "-t", $refAlignerDir, "-b", "$moleculesDir/$moleculesFile", "-r", "$bioNanoMapDir/$bioNanoMapFile", "-a", $deNovoOpArgsFile);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
if ($jobStatus != 0)	{
	# print out error
	$errResults = "ERROR: in calling auto noise ".(join(" ", @$cmdRef))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
	dieLog("$errResults\n");
} # if jobStatus
# check the stdout of the alignment for "END of output" keyword
errCheckRefAligner("$autoNoiseOutDir/contigs/auto_noise/autoNoise1.stdout", "END of output", "Alignment of molecules to BioNano genome map (auto noise) complete.", "ERROR: Alignment of molecules to BioNano genome map (auto noise) cannot be completed.");
	
# further check to make sure that auto noise errbin file is present, just in case there were errors in auto noise step
if (! -e "$autoNoiseOutDir/contigs/auto_noise/autoNoise1.errbin")	{
	dieLog("ERROR: after auto noise is called, errbin still does not exist $autoNoiseOutDir/contigs/auto_noise/autoNoise1.errbin\n");
} # if e 2

#################### align molecules to BioNano ####################
print "\nBeginning aligning molecules to BioNano genome maps\n";
dieLog ("ERROR: align_molecules.pl: cannot perform alignmol because $alignMolPythonScript does not exist\n") if (! -e $alignMolPythonScript);
@$cmdRef = ("/usr/bin/python", $alignMolPythonScript, "-T", $numThreads, "-j", $globalMaxThreads, "-o", $alignMolBioNanoOutDir, "-t", $refAlignerDir, "-b", "$moleculesDir/$moleculesFile", "-q", "$bioNanoMapDir/$bioNanoMapFile", "-a", $deNovoOpArgsFile, "-E", "$autoNoiseOutDir/contigs/auto_noise/autoNoise1.errbin", "-p", $deNovoPipelineDir); 
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
if ($jobStatus != 0)	{
	# print out error
	$errResults = "ERROR: in calling molecule alignment to BioNano genome map ".(join(" ", @$cmdRef))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
	dieLog("$errResults\n");
} # if jobStatus
# check the stdout of the alignment for "END of output" keyword
errCheckRefAligner("$alignMolBioNanoOutDir/$bioNanoMapFilePrefix.stdout", "END of output", "Alignment of molecules to BioNano genome map completed.", "ERROR: Alignment of molecules to BioNano genome map cannot be completed.");

#################### align molecules to BioNano ####################
print "\nBeginning aligning molecules to hybrid scaffolds\n";
@$cmdRef = ("/usr/bin/python", $alignMolPythonScript, "-T", $numThreads, "-j", $globalMaxThreads, "-o", $alignMolHybridOutDir, "-t", $refAlignerDir, "-b", "$moleculesDir/$moleculesFile", "-q", "$hybridMapDir/$hybridMapFile", "-a", $deNovoOpArgsFile, "-E", "$autoNoiseOutDir/contigs/auto_noise/autoNoise1.errbin", "-p", $deNovoPipelineDir);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
if ($jobStatus != 0)	{
	# print out error
	$errResults = "ERROR: in calling molecule alignment to hybrid map ".(join(" ", @$cmdRef))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
	dieLog("$errResults\n");
} # if jobStatus
errCheckRefAligner("$alignMolHybridOutDir/$hybridMapFilePrefix.stdout", "END of output", "Alignment of molecules to BioNano genome map completed.", "ERROR: Alignment of molecules to BioNano genome map cannot be completed.");

sub errCheckRefAligner  {
	my ($file, $completeString, $completeMsg, $dieMsg) = @_;
	open(IN, "$file") or dieLog ("ERROR: Cannot open $file: $!\n");
	my $presence = 0;
	while (my $line = <IN>) {
		if ($line =~ /$completeString/) {
			# if the line contains the string that indicates successful completion
			$presence = 1;
		} # if line
	} # while line
	if ($presence == 0)     {
		dieLog ("ERROR: $dieMsg\n");
	} else  {
		print "$completeMsg\n";
	} # if presence
	close IN;
} # errCheckRefAligner

sub runCommand  {
	my ($argsRef) = @_;

	my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef);


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
			}elsif($fh == $err) {
				$errResults .= "$line";
			}else{
				dieLog ("ERROR: This should never execute!");
			} # if fh
		} # foreach fh
	} # while

	my $ret=waitpid ($pid, 0); # reap the exit code
	return ($outResults, "$errResults");
} # runCommand

