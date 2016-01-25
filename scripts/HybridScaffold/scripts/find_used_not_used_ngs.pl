# $Id: find_used_not_used_ngs.pl 4116 2015-09-16 21:23:20Z apang $

#!/usr/bin/perl -w

use strict;
use warnings;
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

# this script identifies which sequence entries (id < $idShift) went into the hybrid scaffold, it then looks at the original NGS and separates those that did and did not paricipate

my $homeDir = $ARGV[0];
my $inFile = "step1.merge.pairs.txt";
my $idShift = $ARGV[1];	# this is the value used during the merge step to shift the BioNano cmapIds in order to distinguish them from sequence entries

my $oriCmapFile = $ARGV[2];
my $oriCmapDir = $oriCmapFile;  my @oriCmapDirContent = split(/\//, $oriCmapDir);       my $oriCmapFileName = $oriCmapDirContent[$#oriCmapDirContent];  $oriCmapDir =~ s/\/*$oriCmapFileName$//; # removed the file name
$oriCmapFile = $oriCmapFileName;

my $outIdDir = $homeDir;
my $outUsedIdFile = "all_used_ngs_id.txt";
my $outNotUsedIdFile = "all_not_used_ngs_id.txt";

my $outCmapDir = $homeDir;
my $outUsedCmapPrefix = "all_used_ngs";
my $outNotUsedCmapPrefix = "all_not_used_ngs";

my $refaligner = $ARGV[3];

my $usedIdsRef = getUsedIds($homeDir, $inFile, $idShift);
# now read in the original cmap file and record those cmap ids that are not in the recorded hash
my $notUsedIdsRef = getNotUsed($oriCmapDir, $oriCmapFile, $usedIdsRef);
# print a list file that contains those used ids
printIds($outIdDir, $outUsedIdFile, $usedIdsRef);
# print a list file that contains those not used ids
printIds($outIdDir, $outNotUsedIdFile, $notUsedIdsRef);

my @cmd = ();	my $cmdRef = \@cmd;
my ($outResults, $errResults, $jobStatus) = ("", "", 0);
if (scalar(keys %$usedIdsRef) > 0)	{
	# print a cmap file of those used ids 
	@$cmdRef = ($refaligner, "-f", "-merge", "-selectidf", "$outIdDir/$outUsedIdFile", "-i", "$oriCmapDir/$oriCmapFile", "-o", "$outCmapDir/$outUsedCmapPrefix", "-minsites", 0, "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	my ($outResults, $errResults, $jobStatus) = runCommand($cmdRef);
	if ($jobStatus != 0)    {
        	# print out error
        	$errResults = "ERROR: in calling ".(join(" ", @$cmdRef))."with exit code $jobStatus and info: $outResults; ".$errResults;
        	dieLog("$errResults\n");
	} # if jobStatus
} else	{
	# if there is no usedIds, something is wrong
	my $errResults = "ERROR: the number of used ids is zero. That is, nothing was used in the pairmerge step in MergeNGS_BN";
	dieLog("$errResults\n");
} # if scalar

@$cmdRef = ();
($outResults, $errResults, $jobStatus) = ("", "", 0);
if (scalar(keys %$notUsedIdsRef) > 0)	{
	# print a cmap file of those not used ids
	@$cmdRef = ($refaligner, "-f", "-merge", "-selectidf", "$outIdDir/$outNotUsedIdFile", "-i", "$oriCmapDir/$oriCmapFile", "-o", "$outCmapDir/$outNotUsedCmapPrefix", "-minsites", 0, "-stdout", "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults, $jobStatus) = runCommand($cmdRef);
	if ($jobStatus != 0)    {
        	# print out error
        	$errResults = "ERROR: in calling ".(join(" ", @$cmdRef))."with exit code $jobStatus and info: $outResults; ".$errResults;
        	dieLog("$errResults\n");
	} # if jobStatus
} else	{
	# if there is no NotUsedIds (i.e. everything was used in the pairmerge step in MergeNGS_BN)
	# then create an empty file
	open(OUT, ">$outCmapDir/$outNotUsedCmapPrefix.cmap") or dieLog "ERROR: Cannot write to $outCmapDir/$outNotUsedCmapPrefix.cmap: $!\n";
	close OUT;
} # if scalar




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
        my $childExitStatus = $? >> 8;
        return ($outResults, $errResults, $childExitStatus);
} # runCommand

sub numeric     {$a <=> $b}

sub printIds    {
        my ($dir, $file, $idsRef) = @_;
        open(OUT, ">$dir/$file") or dieLog "ERROR: printIds: cannot write to $dir/$file: $!\n";
my $count = 0;  print "printIds: writing to $file\n";
        foreach my $id (sort numeric keys %$idsRef)      {
                print OUT "$id\n";
$count += 1;
        } # foreach id
print "\twritten $count records\n";
        close OUT;
} # printIds

sub getNotUsed	{
	my ($dir, $file, $usedIdsRef) = @_;
	my %notUsedIds = ();
	open(IN, "$dir/$file") or die "getNotUsed: cannot open $dir/$file: $!\n";
print "reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		next if ($line =~ /^#/);
		my @content = split(/\t/, $line);
		my $id = $content[0];
		next if (exists $usedIdsRef->{$id});
		$notUsedIds{$id} = 1;	# return those ids of those entries that did not participate in the scaffolding process
	} # while line
	close IN;
print "\tread in ".(keys %notUsedIds)." records\n";
	return \%notUsedIds;
} # getNotUsed

# extracts ids of entries used in the scaffolding process
sub getUsedIds	{
	my ($dir, $file, $idShift) = @_;
	my %usedIds = ();
	open(IN, "$dir/$file") or dieLog "ERROR: getUsedIds: cannot open $dir/$file: $!\n";
my $count = 0;	print "reading in $file\n";
	my $skip = <IN>;
	while (my $line = <IN>) {
                chomp $line;
                $line =~ s/\r//g;
                $line =~ s/"//g;
                $line =~ s/\s+/\t/g;

                my @content = split(/\t/, $line);
                my ($id1, $id2, $hybridId, $stage) = @content[1..3];
                
                if ($id1 < $idShift)    {
                        $usedIds{$id1} = 1;
$count += 1;
                } # if id1

                if ($id2 < $idShift)    {
                        $usedIds{$id2} = 1;
$count += 1;
                } # if id2
        } # while line
print "\tread in $count entries\n";
print "\tbut there are ".(scalar(keys %usedIds)). " entries\n";
	close IN;
	return \%usedIds;
} # getUsedIds


