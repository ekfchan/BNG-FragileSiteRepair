# $Id: cut_conflicts.pl 4323 2015-12-04 19:35:49Z apang $

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

my ($oriGMFile) = ("");	# assumes that the oriGMFile has chimeric quality score
my ($oriSeqFile) = ("");
my ($conflictFile) = ("", "break_point.txt");
my ($outDir) = ("");	# sequence sub dir is $outDir/NGS, and GM sub dir is $outDir/BN
my ($outGMFile) = ("");
my ($outSeqFile) = ("");
my ($cutGmFile, $cutSeqFile) = ("", "");
my ($idOffSet) = 50000;
my ($windowSize, $qScoreThreshold, $covThreshold) = (10000, 35, 10);
my $refAligner = "";
my $modBkptStatusFile = "";

GetOptions	(
	"oriGMFile=s"		=>	\$oriGMFile,
	"oriSeqFile=s"		=>	\$oriSeqFile,
	"conflictFile=s"	=>	\$conflictFile,
	"outDir=s"		=>	\$outDir,
	"outGMFile=s"		=>	\$outGMFile,
	"outSeqFile=s"		=>	\$outSeqFile,
	"windowSize=i"		=>	\$windowSize,
	"qScoreThreshold=i"	=>	\$qScoreThreshold,
	"covThreshold=i"	=>	\$covThreshold,
	"refAligner=s"		=>	\$refAligner,
	"modBkptStatusFile:s"	=>	\$modBkptStatusFile		# optional
) or dieLog ("ERROR: cut_conflicts: error in command line arguments.\n");

# define some of the outputs (except, outGMFile and outSeqFile)
mkdir $outDir if (! -e $outDir);
my ($outBkptFile, $outSeqCoordFile, $outGMCoordFile) = ("$outDir/conflicts_cut_status.txt", "$outDir/auto_cut_NGS_coord_translation.txt", "$outDir/auto_cut_BN_coord_translation.txt");
my ($outSeqPreCutAnnotBedFile, $outGMPreCutAnnotBedFile, $outGMPreCutProjectedNGSCoordBedFile) = ("$outDir/ngs_pre_cut_annotations.bed", "$outDir/bn_pre_cut_annotations.bed", "$outDir/bn_pre_cut_projected_ngs_coord_annotations.bed");
my ($outSeqCutPreDiscardCmap, $outGMCutPreDiscardCmap) = ("$outDir/ngs_cut_pre_exclude.cmap", "$outDir/bn_cut_pre_exclude.cmap");
my ($outDiscardSeqFile, $outDiscardGMFile) = ("$outDir/ngs_exclude.cmap", "$outDir/bn_exclude.cmap");
my ($outKeepSeqIdList, $outDiscardSeqIdList) = ("$outDir/keep_ngs_ids.txt", "$outDir/exclude_ngs_ids.txt");
my ($outKeepGMIdList, $outDiscardGMIdList) = ("$outDir/keep_bn_ids.txt", "$outDir/exclude_bn_ids.txt");
# bed file annotation window size, a small size to annotate where cuts were made
my $bedAnnotLength = 1000;	


my %conflictsAll = ();	my $conflictsAllRef = \%conflictsAll;
my $noQScoreFlag = 0;	# remember if there is chimeric quality score in the BioNano assembly
if ($modBkptStatusFile eq "")	{
	# no manually modified break point status file was given, so
	# 1) read in break point file
	# 2) check the BioNano chimeric scores at the conflict loci
	# 3) for each conflict locus, determine whether it is the BioNano or the sequence assembly that is needed to be cut
	# 4) write an updated break point status file

	# extract the breakpoints (the query is the BioNano gm)
	my @headerLines = ();	my %conflictsQuery = ();	my $headerLinesRef = \@headerLines;	my $conflictsQueryRef = \%conflictsQuery;
	($headerLinesRef, $conflictsQueryRef, $conflictsAllRef) = getConflicts($conflictFile, $headerLinesRef, $conflictsQueryRef, $conflictsAllRef);
	$conflictsQueryRef = extendByWindowSize($conflictsQueryRef, $windowSize);
	$conflictsQueryRef = sortByCoord($conflictsQueryRef, "start");
	# now extract the chimeric quality score
	my %qScores = ();	my $qScoresRef = \%qScores;
	($qScoresRef, $noQScoreFlag) = getQScores($oriGMFile, $qScoresRef);
	print "noQScore = $noQScoreFlag\n";
	if ($noQScoreFlag == 0)	{
		# assumes that the cmap file is already sorted by the label position in the gm
		$conflictsQueryRef = findOverlap($conflictsQueryRef, $qScoresRef);
		$conflictsQueryRef = flagCut($conflictsQueryRef, $qScoreThreshold, $covThreshold);
		
		# update the conflictsAllRef hash
		$conflictsAllRef = updateConflictsAll($conflictsQueryRef, $conflictsAllRef);
	} else	{
		### remember, if there is no quality column, print warning message
		print "WARNING: In the cut conflict stage, but there were no quality scores in the BioNano assembly\n";
	} # if noQScoreFlag

	# print an updated breakpoint file
	printUpdatedBreakPointFile($outBkptFile, $headerLinesRef, $conflictsAllRef);
} else	{
	# a manually modified break point status file was given, so
	# read in the break point status file
	my @headerLines = ();	my $headerLinesRef = \@headerLines;
	($headerLinesRef, $conflictsAllRef) = getManualBreakPointFile($modBkptStatusFile, $headerLinesRef, $conflictsAllRef);
	# print an updated breakpoint file
	printUpdatedBreakPointFile($outBkptFile, $headerLinesRef, $conflictsAllRef);
} # if modBkptStatusFile

# now read in the sequence file and the gm file to figure out the largest id in each file
my ($maxSeqId, $seqLengthsRef) = readCmapIdLength($oriSeqFile);
my ($maxGMId, $gMLengthsRef) = readCmapIdLength($oriGMFile);

# print coordinate translations, unless the fragment has no conflict, is correct whenever there is conflict, and  not to be discarded
my ($seqCutIdsRef, $seqDiscardIdsRef, $seqToCutIdFragments) = printCoordTranslation($outSeqCoordFile, $conflictsAllRef, "ref", $seqLengthsRef, $maxSeqId);	
my ($gMCutIdsRef, $gMDiscardIdsRef, $gMToCutIdFragments) = printCoordTranslation($outGMCoordFile, $conflictsAllRef, "qry", $gMLengthsRef, $maxGMId);
# print bed file co-ordinates (ngs pre-cut ngs coordinates; BN pre-cut BN coordinates; BN pre-cut BN coordinate projected onto NGS coordinates)
printPreCutBedRespectiveCoords($outSeqPreCutAnnotBedFile, $seqCutIdsRef, $bedAnnotLength, "NGS", "255,102,0");
printPreCutBedRespectiveCoords($outGMPreCutAnnotBedFile, $gMCutIdsRef, $bedAnnotLength, "BN", "102,0,255");
printGMPreCutBedRefCoords($outGMPreCutProjectedNGSCoordBedFile, $conflictsAllRef, $seqDiscardIdsRef, $gMDiscardIdsRef, $bedAnnotLength, "102,0,255");

# now tell RefAligner to do the cut (remember that the co-ordinates fed into RefAligner MUST be in kilobase not bp)
callRefAlignerBreak($refAligner, $oriSeqFile, $outSeqCutPreDiscardCmap, $maxSeqId, $seqCutIdsRef, $noQScoreFlag, "ref");
callRefAlignerBreak($refAligner, $oriGMFile, $outGMCutPreDiscardCmap, $maxGMId, $gMCutIdsRef, $noQScoreFlag, "qry");

### now this part deals with the fragments that need to be discarded
my $toKeepSeqIdsRef = distinguishKeepDiscards($seqLengthsRef, $seqToCutIdFragments, $seqDiscardIdsRef);	# toKeepSeqIdsRef are the ones to keep, seqDiscardIdsRef are the ones to discard	
my $toKeepGMIdsRef = distinguishKeepDiscards($gMLengthsRef, $gMToCutIdFragments, $gMDiscardIdsRef);	# toKeepGMIdsRef are the ones to keep, gMDiscardIdsRef are the ones to discard
# call RefAligner to make two set of files
callRefAlignerKeepDiscard($refAligner, $outSeqCutPreDiscardCmap, $outSeqFile,  $toKeepSeqIdsRef, $outKeepSeqIdList);	# the input is just what RefAligner results from the break command above
callRefAlignerKeepDiscard($refAligner, $outSeqCutPreDiscardCmap, $outDiscardSeqFile, $seqDiscardIdsRef, $outDiscardSeqIdList);
callRefAlignerKeepDiscard($refAligner, $outGMCutPreDiscardCmap, $outGMFile, $toKeepGMIdsRef, $outKeepGMIdList);
callRefAlignerKeepDiscard($refAligner, $outGMCutPreDiscardCmap, $outDiscardGMFile, $gMDiscardIdsRef, $outDiscardGMIdList);

sub callRefAlignerKeepDiscard	{
	my ($refAligner, $inFile, $outFile, $idsRef, $idListFile) = @_;
	my $outFilePrefix = $outFile;	$outFilePrefix =~ s/\.cmap$//;

	# now make a list txt file, whose purpose is to feed into RefAligner with selectidf
	printIdListFile($idListFile, $idsRef);

	my @cmd = ();	my $cmdRef = \@cmd;
	@$cmdRef = ($refAligner, "-f", "-minsites", 0, "-i", "$inFile", "-o", "$outFilePrefix", "-merge", "-selectidf", "$idListFile");
	push(@$cmdRef, "-stdout");	push(@$cmdRef, "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	my ($outResults, $errResults, $jobStatus) = runCommand($cmdRef);
	if ($jobStatus != 0)	{
		# print out error
		$errResults = "ERROR: callRefAlignerKeepDiscard: in calling ".(join(" ", @$cmdRef))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
		dieLog("$errResults\n");
	} # if jobStatus
} # callRefAlignerKeepDiscard

sub printIdListFile	{
	my ($file, $idsRef) = @_;
	# this subroutine just output to a file with id on each line
	open(OUT, ">$file") or dieLog "ERROR: printIdListFile: cannot write to $file: $!\n";
	foreach my $id (keys %$idsRef)	{
		print OUT "$id\n";
	} # foreach id
	close OUT;
} # printIdListFile

sub distinguishKeepDiscards	{
	my ($allIdsRef, $cuttedIdsRef, $discardIdsRef) = @_;
	# remember allIdsRef has the lengths of all the pre-cut fragments
	# cuttedIdsRef has ids of those that went through cut, its old id and all new ids
	# discardIdsRef has ids of those need to be thrown away (MUST be mutually exclusive from the cuttedIdsRef)
	my %toKeepIds = ();
	
	# now deal with those fragments where no cut was involved
	foreach my $id (keys %$allIdsRef)	{
		next if (exists $discardIdsRef->{$id});	# skip those that have been flagged to be discard
		$toKeepIds{$id} = 1;			
	} # foreach id

	# add in the new ids generated from performing cuts
	foreach my $id (keys %$cuttedIdsRef)	{
		for (my $i = 0; $i < scalar(@{$cuttedIdsRef->{$id}}); $i += 1)	{
			$toKeepIds{$cuttedIdsRef->{$id}[$i]{newId}} = 1;		
		} # for i
	} # foreach id

	return \%toKeepIds;
} # distinguishKeepDiscards

sub printPreCutBedRespectiveCoords	{
	my ($file, $toCutIdsRef, $bedAnnotLength, $dataType, $colourString) = @_;
	# this subroutine looks at an array of cut coordinates and determine the bed co-ordinates (in reference co-ordinates)
	my %toCutBedFragments = ();
	foreach my $id (keys %$toCutIdsRef)	{
		for (my $i = 0; $i < scalar(@{$toCutIdsRef->{$id}}); $i += 1)	{
			my $cRef = $toCutIdsRef->{$id}[$i];
			push(@{$toCutBedFragments{$id}}, {start => $cRef->{start} - $bedAnnotLength, end => $cRef->{start} + $bedAnnotLength, oriPosition => $cRef->{start},
				descString => "$dataType"."_proposed_cut_$dataType"."Id_$id"."_$dataType"."Position_$cRef->{start}"});
		} # for i
	} # foreach id

	# print to file
	open(OUT, ">$file") or dieLog "ERROR: printNGSCutBedRefCoords: cannot write to $file: $!\n";
	my $count = 1;
	foreach my $id (sort numeric keys %toCutBedFragments)	{
		for (my $i = 0; $i < scalar(@{$toCutBedFragments{$id}}); $i += 1)	{
			my $aRef = $toCutBedFragments{$id}[$i];
			print OUT join("\t", $id, $aRef->{start}, $aRef->{end}, $aRef->{descString}, $count, "+", $aRef->{start}, $aRef->{end}, $colourString)."\n";
			$count += 1;
		} # for i
	} # foreach id
	close OUT;
	
} # printPreCutBedRespectiveCoords

sub printGMPreCutBedRefCoords	{
	my ($file, $conflictsAllRef, $toDiscardRefIdsRef, $toDiscardQryIdsRef, $bedAnnotLength, $colourString) = @_;
	my %toCutRefIds = ();		# it is okay for this array to store redundant reference coordinates, as multiple queries can conflict with the same reference at the same reference position
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);	
		# if either refId or qryId is -1, then skip
		next if ($refId =~ /^-1$/ || $qryId =~ /^-1$/);
		# if either ref or qry needs to be discard, then skip
		next if (exists $toDiscardRefIdsRef->{$refId} || exists $toDiscardQryIdsRef->{$qryId});
		
		# now build an array of cut co-ordinates on the reference, but only if the qry needs to be cut
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			
			if ($aRef->{qryLeftBkptStatus} =~ /^cut$/)	{
				push(@{$toCutRefIds{$refId}}, {refCoord => $aRef->{refLeftBkpt}, 
					refStart => $aRef->{refLeftBkpt} - $bedAnnotLength, refEnd => $aRef->{refLeftBkpt} + $bedAnnotLength,
					qryId => $qryId, qryCoord => $aRef->{qryLeftBkpt},
					descString => "BN_proposed_cut_BNId_$qryId"."_BNPosition_$aRef->{qryLeftBkpt}"."_projectedNGSId_$refId"."_projectedNGSPosition_$aRef->{refLeftBkpt}"});
			} # if left query breakpoint
			
			if ($aRef->{qryRightBkptStatus} =~ /^cut$/)	{
				push(@{$toCutRefIds{$refId}}, {refCoord => $aRef->{refRightBkpt}, 
					refStart => $aRef->{refRightBkpt} - $bedAnnotLength, refEnd => $aRef->{refRightBkpt} + $bedAnnotLength, 
					qryId => $qryId, qryCoord => $aRef->{qryRightBkpt},
					descString => "BN_proposed_cut_BNId_$qryId"."_BNPosition_$aRef->{qryRightBkpt}"."_projectedNGSId_$refId"."_projectedNGSPosition_$aRef->{refRightBkpt}"});
			} # if right query breakpoint

		} # for i
	} # foreach theKey
	
	# now print a bed file
	open(OUT, ">$file") or dieLog "ERROR: printGMPreCutBedRefCoords: cannot write to $file: $!\n";
	my $count = 1;
	foreach my $refId (sort numeric keys %toCutRefIds)	{
		for (my $i = 0; $i < scalar(@{$toCutRefIds{$refId}}); $i += 1)	{
			my $aRef = $toCutRefIds{$refId}[$i];
			print OUT join("\t", $refId, $aRef->{refStart}, $aRef->{refEnd}, $aRef->{descString}, $count, "+", $aRef->{refStart}, $aRef->{refEnd}, $colourString)."\n";
			$count += 1;
		} # for i
	} # foreach refId
	close OUT;
} # printGMPreCutBedRefCoords


sub callRefAlignerBreak	{
	my ($refAligner, $inFile, $outFile, $maxId, $cutIdsRef, $noQScoreFlag, $assemblyType) = @_;
	# actual call to RefAligner to do cut (remember that the co-ordinates must be in kilobase)
	# Aside: remember that AssignAlignType prints the breakpoints such that the conflict label is always kept with the aligned match group
	my $outFilePrefix = $outFile;	$outFilePrefix =~ s/\.cmap$//;
	
	my @cmd = (); my $cmdRef = \@cmd;
	@$cmdRef = ($refAligner, "-f", "-minsites", 0,"-i", "$inFile", "-o", $outFilePrefix, "-merge");
	#if ($assemblyType eq "ref")	{
		# if the assembly is sequence: make sure we do not pad extra 20bp to the front (or back) of the contig [Remember we don't change sequence length/label UNLESS we are doing just cut]
	#	push(@$cmdRef, "-minEnd", "0.0");
	#} # if assemblyType
	if ($noQScoreFlag == 0 && scalar(keys %$cutIdsRef) > 0)	{
		# only break if there was chimeric score provided 
		push(@$cmdRef, "-break");
		push(@$cmdRef, $maxId);
	} # if noQScoreFlag
	foreach my $id (keys %$cutIdsRef)	{
		push(@$cmdRef, $id);
		for (my $i = 0; $i < scalar(@{$cutIdsRef->{$id}}); $i += 1)	{
			my $theCoord = $cutIdsRef->{$id}[$i]{start} / 1000;
			my $theCoordString = sprintf("%.4f", $theCoord);
			push(@$cmdRef, $theCoordString);
		} # for i
	} # foreach id
	push(@$cmdRef, "-stdout");	push(@$cmdRef, "-stderr");
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	my ($outResults, $errResults, $jobStatus) = runCommand($cmdRef);
	if ($jobStatus != 0)	{
		# print out error
		$errResults = "ERROR: callRefAlignerBreak: in calling ".(join(" ", @$cmdRef))." with exit code $jobStatus; out info: $outResults; error info: $errResults";
		dieLog("$errResults\n");
	} # if jobStatus
} # callRefAlignerBreak

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

sub printCoordTranslation	{
	my ($file, $conflictsAllRef, $assemblyType, $lengthsRef, $maxId) = @_;
	# this subroutine will print out the pre-cut and post-cut ids and co-ordinates in a file
	# and return a list of ids and the co-ordinates to be cut
	# but of course, do NOT mark a contig to be cut, if it is to be discarded 

	my $toDiscardIdsRef = determineDiscardFragments($conflictsAllRef, $assemblyType);

	my $toCutIdsRef = determineCutCoords($conflictsAllRef, $assemblyType, $toDiscardIdsRef);

	# sort by coord
	$toCutIdsRef = sortByCoord($toCutIdsRef, "start");	
	
	# determine the fragments
	my $toCutIdFragmentsRef = determineFragments($toCutIdsRef, $lengthsRef, $maxId);

	open(OUT, ">$file") or dieLog ("ERROR printCoordTranslation: cannot write to $file");
	print OUT join("\t", "oldId", "oldStart", "oldEnd", "newId", "newStart", "newEnd")."\n";
	# first print out those that will be cut
	foreach my $id (sort numeric keys %$toCutIdFragmentsRef)	{
		for (my $i = 0; $i < scalar(@{$toCutIdFragmentsRef->{$id}}); $i += 1)	{
			my $aRef = $toCutIdFragmentsRef->{$id}[$i];
			print OUT join("\t", $id, $aRef->{oldStart}, $aRef->{oldEnd}, $aRef->{newId}, $aRef->{newStart}, $aRef->{newEnd})."\n";
		} # for i
	} # foreach id

	# then print out the ones that do not need to be cut
	foreach my $id (sort numeric keys %$lengthsRef)	{
		next if (exists $toCutIdFragmentsRef->{$id});	# skip those that needed to be cut
		my $endCoord = $lengthsRef->{$id} - 1;
		print OUT join("\t", $id, 0, $endCoord, $id, 0, $endCoord)."\n";
	} # foreach id
	close OUT;	

	return ($toCutIdsRef, $toDiscardIdsRef, $toCutIdFragmentsRef);
} # printCoordTranslation

#sub printCutAnnotationsBed	{
#	my ($file, $toCutIdsRef, $lengthsRef, $maxId) = @_;
#	# this subroutine looks at an array to determine the cut coordinates and determine the fragment coordinates (same as sub determineFragments), but also annotates 1000 bp to the beginning (or end) of the new fragment
#	# this is purely for visualization purpose for users to know whether a new fragment was generated from a cut, and whether the cut was at its head or tail
#	my $bedAnnotLength = 1000;
#	open(OUT, ">$file") or dieLog "ERROR: printCutAnnotationsBed: cannot write to $file: $!\n";
#	foreach my $id (keys %$toCutIdsRef)	{
#		dieLog "ERROR: printCutAnnotationsBed: cannot find the original length of the contig whose id=$id\n" if (! exists $lengthsRef->{$id});
#		for (my $i = 0; $i < scalar(@{$toCutIdsRef->{$id}}); $i += 1)	{
#			my $cRef = $toCutIdsRef->{$id}[$i];
#			if ($i == 0)	{
#				# first fragment
#				if (scalar(@{$toCutIdsRef->{$id}}) == 1)	{
#					# also the last fragment
#					my ($newStart, $newEnd, $newId) = (0, $cRef->{start}, $maxId * $i + $id);	my ($newBedStart, $newBedEnd) = ($newEnd - $bedAnnotLength, $newEnd);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";
#					($newStart, $newEnd, $newId) = (0, ($lengthsRef->{$id} - 1) - ($cRef->{start} + 1), $maxId * ($i + 1) + $id);	($newBedStart, $newBedEnd) = ($newStart, $newStart + $bedAnnotLength);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";
#				} else	{
#					my ($newStart, $newEnd, $newId) = (0, $cRef->{start}, $maxId * $i + $id);	my ($newBedStart, $newBedEnd) = ($newEnd - $bedAnnotLength, $newEnd);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";
#					($newStart, $newEnd, $newId) = (0, $toCutIdsRef->{$id}[$i + 1]{start} - $cRef->{start} - 1, $maxId * ($i + 1) + $id);	($newBedStart, $newBedEnd) = ($newStart, $newStart + $bedAnnotLength);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";
#				} # if scalar
#			} elsif ($i == $#{$toCutIdsRef->{$id}})	{
#				# last fragment
#					my ($newStart, $newEnd, $newId) = (0, ($lengthsRef->{$id} - 1) - ($cRef->{start} + 1), $maxId * ($i + 1) + $id);	my ($newBedStart, $newBedEnd) = ($newStart, $newStart + $bedAnnotLength);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";	
#			} else	{
#				# middle fragment
#					my ($newStart, $newEnd, $newId) = (0, $toCutIdsRef->{$id}[$i + 1]{start} - $cRef->{start} - 1, $maxId * ($i + 1) + $id);	my ($newBedStart, $newBedEnd) = ($newStart, $newStart + $bedAnnotLength);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";	
#					($newBedStart, $newBedEnd) = ($newEnd - $bedAnnotLength, $newEnd);
#					print OUT join("\t", $newId, $newBedStart, $newBedEnd, "cut_annotation", -1, "+", $newBedStart, $newBedEnd, "255,102,0")."\n";	
#			} # if i
#		} # for i
#	} # foreach id
#	close OUT;
#} # printCutAnnotationsBed

sub numeric	{	$a	<=>	$b	}

sub determineFragments	{
	my ($toCutIdsRef, $lengthsRef, $maxId) = @_;
	# this subroutine looks at an array of cut coordinates and determine the fragment coordinates
	my %toCutIdFragments = ();
	foreach my $id (keys %$toCutIdsRef)	{
		dieLog "ERROR: determineFragments: cannot find the original length of the contig whose id=$id\n" if (! exists $lengthsRef->{$id});
		for (my $i = 0; $i < scalar(@{$toCutIdsRef->{$id}}); $i += 1)	{
			my $cRef = $toCutIdsRef->{$id}[$i];
			if ($i == 0)	{
				# first fragment
				if (scalar(@{$toCutIdsRef->{$id}}) == 1)	{
					# also the last fragment
					push(@{$toCutIdFragments{$id}}, {oldStart => 0, oldEnd => $cRef->{start}, newStart => 0, newEnd => $cRef->{start}, newId => $maxId * $i + $id});
					push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $lengthsRef->{$id} - 1, newStart => 0, newEnd => ($lengthsRef->{$id} - 1) - ($cRef->{start} + 1), newId => $maxId * ($i + 1) + $id});
				} else	{
					push(@{$toCutIdFragments{$id}}, {oldStart => 0, oldEnd => $cRef->{start}, newStart => 0, newEnd => $cRef->{start}, newId => $maxId * $i + $id});
					push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $toCutIdsRef->{$id}[$i + 1]{start}, newStart => 0, newEnd => $toCutIdsRef->{$id}[$i + 1]{start} - $cRef->{start} - 1, newId => $maxId * ($i + 1) + $id});
				} # if scalar
			} elsif ($i == $#{$toCutIdsRef->{$id}})	{
				# last fragment
					push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $lengthsRef->{$id} - 1, newStart => 0, newEnd => ($lengthsRef->{$id} - 1) - ($cRef->{start} + 1), newId => $maxId * ($i + 1) + $id});
			} else	{
				# middle fragment
					push(@{$toCutIdFragments{$id}}, {oldStart => $cRef->{start} + 1, oldEnd => $toCutIdsRef->{$id}[$i + 1]{start}, newStart => 0, newEnd => $toCutIdsRef->{$id}[$i + 1]{start} - $cRef->{start} - 1, newId => $maxId * ($i + 1) + $id});
			} # if i
		} # for i
	} # foreach id
	return \%toCutIdFragments;
} # determineFragments

sub determineCutCoords	{
	my ($conflictsAllRef, $assemblyType, $toDiscardIdsRef) = @_;
	# this subroutines puts the cut coordinates into an array (after making a non-redundant list of co-ordinates first)
	# again, ignore those fragments that have been flagged as discard
	my %nrCutIds = ();
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			if ($assemblyType eq "ref")	{
				# examine the sequence bkpt status
				$nrCutIds{$refId}{$aRef->{refLeftBkpt}} = 1 if ($aRef->{refLeftBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$refId} && $refId !~ /^-1$/ && $aRef->{refLeftBkpt} != -1);
				$nrCutIds{$refId}{$aRef->{refRightBkpt}} = 1 if ($aRef->{refRightBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$refId} && $refId !~ /^-1$/ && $aRef->{refRightBkpt} != -1);
			} else	{
				# examine the genome map bkpt status
				$nrCutIds{$qryId}{$aRef->{qryLeftBkpt}} = 1 if ($aRef->{qryLeftBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$qryId} && $qryId !~ /^-1$/ && $aRef->{qryLeftBkpt} != -1);
				$nrCutIds{$qryId}{$aRef->{qryRightBkpt}} = 1 if ($aRef->{qryRightBkptStatus} =~ /^cut$/ && ! exists $toDiscardIdsRef->{$qryId} && $qryId !~ /^-1$/ && $aRef->{qryRightBkpt} != -1);
			} # if assemblyType
		} # for i
	} # foreach theKey

	my %toCutIds = ();
	foreach my $theId (keys %nrCutIds)	{
		foreach my $theStart (keys %{$nrCutIds{$theId}})	{
			push(@{$toCutIds{$theId}}, {start => $theStart});
		} # foreach theStart
	} # foreach theId

	return \%toCutIds;
} # determineCutCoords

sub determineDiscardFragments	{
	my ($conflictsAllRef, $assemblyType) = @_;
	my %toDiscardIds = ();
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $cRef = $conflictsAllRef->{$theKey}[$i];
			$toDiscardIds{$refId} = 1 if ($cRef->{refStatus} =~ /^exclude$/ && $refId !~ /^-1$/ && $assemblyType =~ /^ref$/i);
			$toDiscardIds{$qryId} = 1 if ($cRef->{qryStatus} =~ /^exclude$/ && $qryId !~ /^-1$/ && $assemblyType =~ /^qry$/i);
		} # for i
	} # foreach theKey
	return \%toDiscardIds;
} # determineDiscardFragments

sub readCmapIdLength	{
	my ($file) = @_;
	# read in a cmap file, and determine the largest sequence id
	my $maxSeqId = -1;
	my %theLengths = ();
	open(IN, "$file") or dieLog "ERROR: readCmapId: reading in $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;		
		next if ($line =~ /^#/);	# skip header lines
		my @content = split(/\t/, $line);
		my ($id, $aLength) = @content[0..1];
		$aLength = int($aLength);	# again change to integer
		$theLengths{$id} = $aLength;
		$maxSeqId = ($maxSeqId < $id) ? ($id) : ($maxSeqId);
	} # while line
	close IN;
	return ($maxSeqId, \%theLengths);	
} # readCmapIdLength

sub getManualBreakPointFile	{
	# this subroutine reads in manually modified break point status file, and return a hash containing that information
	my ($file, $headerLinesRef, $conflictsAllRef) = @_;
	open(IN, $file) or dieLog "ERROR: getManualBreakPointFile: cannot read in $file\n: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^(xMapId|id)/i)	{
			push(@$headerLinesRef, $line);		# capture the header line
			next;
		} # if line
		# store the conflict information into the conflictsAll hash
		$line =~ s/\s+/\t/g;	$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		my $xId = $content[0];
		my ($refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus) = @content[2..8];
		my ($qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus) = @content[10..16]; 
		
		breakPointFileCheck($line, $xId, $refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus,
						$qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus);
		# change the bkpt to integer, and status to lower case
		($refLeftBkpt, $refRightBkpt, $qryLeftBkpt, $qryRightBkpt) = (int($refLeftBkpt), int($refRightBkpt), int($qryLeftBkpt), int($qryRightBkpt));
		($refLeftBkptStatus, $refRightBkptStatus, $qryLeftBkptStatus, $qryRightBkptStatus) = (lc($refLeftBkptStatus), lc($refRightBkptStatus), lc($qryLeftBkptStatus), lc($qryRightBkptStatus));
		($refStatus, $qryStatus) = (lc($refStatus), lc($qryStatus));
		# record
		push(@{$conflictsAllRef->{"$xId\t$refId\t$qryId"}}, {	refLeftBkpt => $refLeftBkpt, qryLeftBkpt => $qryLeftBkpt,
									refRightBkpt => $refRightBkpt, qryRightBkpt => $qryRightBkpt,
									refAlignOrientation => $refAlignOrientation, qryAlignOrientation => $qryAlignOrientation,
									refLeftBkptStatus => $refLeftBkptStatus, qryLeftBkptStatus => $qryLeftBkptStatus, 
									refRightBkptStatus => $refRightBkptStatus, qryRightBkptStatus => $qryRightBkptStatus,
									refStatus => $refStatus, qryStatus => $qryStatus
		});
	} # while line
	close IN;
	return ($headerLinesRef, $conflictsAllRef);
} # getManualBreakPointFile

sub breakPointFileCheck	{
	my ($line, $xId, $refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $refLeftBkptStatus, $refRightBkptStatus, $refStatus,
			$qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation, $qryLeftBkptStatus, $qryRightBkptStatus, $qryStatus) = @_;
	
	dieLog ("ERROR: getManualBreakPointFile: Incorrect XMapId, must be a positive integer or -1. Line = $line\n") if ($xId !~ /^\d+$/ && $xId !~ /^-1$/);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect refId, must be a positive integer or -1. Line = $line\n") if ($refId !~ /^\d+$/ && $refId !~ /^-1$/);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect qryId, must be a positive integer or -1. Line = $line\n") if ($qryId !~ /^\d+$/ && $qryId !~ /^-1$/);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect conflict reference co-ordinate. Line = $line\n") if (!($refRightBkpt =~ /^\d+\.?\d*$/ || $refRightBkpt == -1) || !($refLeftBkpt =~ /^\d+\.?\d*$/ || $refLeftBkpt == -1));
	dieLog ("ERROR: getManualBreakPointFile: Incorrect conflict query co-ordinate. Line = $line\n") if (!($qryRightBkpt =~ /^\d+\.?\d*$/ || $qryRightBkpt == -1) || !($qryLeftBkpt =~ /^\d+\.?\d*$/ || $qryLeftBkpt == -1));
	dieLog ("ERROR: getManualBreakPointFile: Incorrect conflict reference conflict point status, must be okay, cut or -. Line = $line\n") if ($refRightBkptStatus !~ /^(okay|cut|-)$/i || $refLeftBkptStatus !~ /^(okay|cut|-)$/i);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect conflict query conflict point status, must be okay, cut or -. Line = $line\n") if ($qryRightBkptStatus !~ /^(okay|cut|-)$/i || $qryLeftBkptStatus !~ /^(okay|cut|-)$/i);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect reference status, must be okay or exclude. Line = $line\n") if ($refStatus !~ /^(okay|exclude|-)$/i);
	dieLog ("ERROR: getManualBreakPointFile: Incorrect query status, must be okay or exclude. Line = $line\n") if ($qryStatus !~ /^(okay|exclude|-)$/i);
} # breakPointFileCheck

sub printUpdatedBreakPointFile	{
	my ($file, $headerLinesRef, $conflictsAllRef) = @_;
	open(OUT, ">$file") or dieLog "ERROR printUpdatedBreakPointFile: cannot write to $file: $!\n";
	my @headerContent = split(/\t/, $headerLinesRef->[0]);
	print OUT join("\t", @headerContent[0..5])."\t".join("\t", "ref_leftBkpt_toCut", "ref_rightBkpt_toCut", "ref_toDiscard")."\t";
	print OUT join("\t", @headerContent[6..10])."\t".join("\t", "qry_leftBkpt_toCut", "qry_rightBkpt_toCut", "qry_toDiscard")."\n";

	print OUT join("\t", "id/-1", 	"ref", "id/-1", "position/-1", "position/-1", "+/-", "okay/cut/-", "okay/cut/-", "okay/exclude/-")."\t";
	print OUT join("\t", 		"qry", "id/-1", "position/-1", "position/-1", "+/-", "okay/cut/-", "okay/cut/-", "okay/exclude/-")."\n";
	
	foreach my $theKey (keys %$conflictsAllRef)	{
		my ($xId, $refId, $qryId) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$conflictsAllRef->{$theKey}}); $i += 1)	{
			my $aRef = $conflictsAllRef->{$theKey}[$i];
			print OUT join("\t", $xId, 	"ref", $refId, $aRef->{refLeftBkpt}, $aRef->{refRightBkpt}, $aRef->{refAlignOrientation}, $aRef->{refLeftBkptStatus}, $aRef->{refRightBkptStatus}, $aRef->{refStatus})."\t";
			print OUT join("\t", 		"qry", $qryId, $aRef->{qryLeftBkpt}, $aRef->{qryRightBkpt}, $aRef->{qryAlignOrientation}, $aRef->{qryLeftBkptStatus}, $aRef->{qryRightBkptStatus}, $aRef->{qryStatus})."\n";
		} # for i
	} # foreach theKey
	close OUT;
} # printUpdatedBreakPointFile

sub updateConflictsAll	{
	my ($conflictsQueryRef, $conflictsAllRef) = @_;
	foreach my $id (keys %$conflictsQueryRef)	{
		for (my $i = 0; $i < scalar(@{$conflictsQueryRef->{$id}}); $i += 1)	{
			my $cRef = $conflictsQueryRef->{$id}[$i];
			my ($xId, $refId, $qryId) = ($cRef->{xId}, $cRef->{refId}, $id);
			# now loop through conflictsAll 
			dieLog "ERROR: updateConflictsAll: cannot find entry for xId=$xId; refId=$refId; qryId=$qryId in original breakpoint information\n" if (! exists $conflictsAllRef->{"$xId\t$refId\t$qryId"});
			for (my $a = 0; $a < scalar(@{$conflictsAllRef->{"$xId\t$refId\t$qryId"}}); $a += 1)	{
				my $aRef = $conflictsAllRef->{"$xId\t$refId\t$qryId"}[$a];
				# check if the query breakpoints are the same
				if (equal($aRef->{qryLeftBkpt}, $cRef->{oriStart}, 3))	{
					# left breakpoint
					if ($cRef->{toCut} == 1)	{
						# the left breakpoint was genome map at fault; label the leftBkpt Status with genome map being toCut
						$aRef->{qryLeftBkptStatus} = "cut";
					} else	{
						# the left breakpoint was sequence at fault; label the leftBkpt status with sequence being toCut
						$aRef->{refLeftBkptStatus} = "cut";
					} # if qRef
				} # if equal	
				if (equal($aRef->{qryRightBkpt}, $cRef->{oriStart}, 3))	{
					# right breakpoint
					if ($cRef->{toCut} == 1)	{
						$aRef->{qryRightBkptStatus} = "cut";	
					} else	{
						$aRef->{refRightBkptStatus} = "cut";
					} # if qRef
				} # if equal
			} # for a	
		} # for i
	} # foreach id
	return $conflictsAllRef;
} # updateConflictsAll

sub equal	{
	# this subrountine checks whether two decimal numbers are the same (up to $decimalPlace precision)
	my ($num1, $num2, $decimalPlace) = @_;
	return sprintf("%.${decimalPlace}g", $num1) eq sprintf("%.${decimalPlace}g", $num2);
}

sub flagCut	{
	# this subroutine loops through all breakpoints and determine if there is any label whose support is below the cutThreshold percentage value
	my ($conflictsRef, $qScoreThreshold, $covThreshold) = @_;
	foreach my $id (keys %$conflictsRef)	{
		for (my $i = 0; $i < scalar(@{$conflictsRef->{$id}}); $i += 1)	{
			my $toCut = 0;
			for (my $j = 0; $j < scalar(@{$conflictsRef->{$id}[$i]{neighbourQScoreRef}}) && $toCut == 0; $j += 1)	{
				if ($conflictsRef->{$id}[$i]{neighbourCovRef}[$j] < $covThreshold || $conflictsRef->{$id}[$i]{neighbourQScoreRef}[$j] < $qScoreThreshold)       {
					$toCut = 1;
				} # if qScore
			} # for j
			$conflictsRef->{$id}[$i]{toCut} = $toCut;
		} # for i
	} # foreach id
	return $conflictsRef;
} # flagCut

sub findOverlap {
        my ($data1Ref, $data2Ref) = @_;
        foreach my $id (keys %$data1Ref)        {
                dieLog "ERROR: findOverlap: cannot find quality scores information for genome map = $id\n" if (! exists $data2Ref->{$id});
                for (my $i = 0; $i < scalar(@{$data1Ref->{$id}}); $i += 1)      {
                        my $d1Ref = $data1Ref->{$id}[$i];
                        my $oIndeciesRef = searchForOverlap($d1Ref->{start}, $d1Ref->{end}, $data2Ref->{$id});  
                        for (my $o = 0; $o < scalar(@$oIndeciesRef); $o += 1)   {
                                # record the coverages and chimeric quality scores of the region
                                my $d2Ref = $data2Ref->{$id}[$oIndeciesRef->[$o]];
                                push(@{$d1Ref->{neighbourCovRef}}, $d2Ref->{coverage});
                                push(@{$d1Ref->{neighbourQScoreRef}}, $d2Ref->{chimQuality});
                        } # for o
                } # for i
        } # foreach id
        return $data1Ref;
} # findOverlap

sub searchForOverlap    {
        my ($qStart, $qEnd, $dataChromRef) = @_;
        my @oIndecies = ();
        my $top = scalar(@$dataChromRef);
        for (my $bot = -1; $top - $bot > 1;)    {
                my $mid = int(($top + $bot) / 2);
                if ($dataChromRef->[$mid]{start} < $qStart)       {
                        $bot = $mid;
                } else  {
                        $top = $mid;
                } # if dataChromRef
        } # for bot

        for (my $i = $top; 0 <= $i && $i < scalar(@$dataChromRef); $i += 1)     {
                my $dRef = $dataChromRef->[$i];
                # any overlap
                if ( $qStart <= $dRef->{start} && $dRef->{start} <= $qEnd )   {
                        push(@oIndecies, $i);
                } # if any overlap
                last if ($qEnd + 3000000 < $dRef->{start});
        } # for i
        return \@oIndecies;
} # searchForOverlap

sub getQScores	{
	my ($file, $qScoresRef) = @_;
	# read in the genome map cmap file, ASSUMING that the chimQuality is on column 10
	my $noQScoreFlag = 0;
	open(IN, "$file") or dieLog "ERROR: getQScores: cannot open $file: $!\n";
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
			warn "cut_conflicts.pl: getQScores: cmap file=$file does not have chimeric quality score\n";
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

sub sortByCoord {
        my ($dataRef, $coord) = @_;
        foreach my $id (keys %$dataRef) {
                @{$dataRef->{$id}} = sort       {
                        $a->{$coord}    <=>     $b->{$coord}
                } @{$dataRef->{$id}}
        } # foreach id
        return $dataRef;
} # sortByCoord

sub extendByWindowSize  {
        my ($conflictsRef, $windowSize) = @_;
        foreach my $id (keys %$conflictsRef)    {
                for (my $i = 0; $i < scalar(@{$conflictsRef->{$id}}); $i += 1)  {
                        my $hRef = $conflictsRef->{$id}[$i];
                        $hRef->{start} = ($hRef->{oriStart} - $windowSize < 1) ? (1) : ($hRef->{oriStart} - $windowSize);
                        $hRef->{end} = $hRef->{oriStart} + $windowSize;
                } # for i
        } # foreach id
        return $conflictsRef;
} # extendByWindowSize

sub getConflicts	{
	# this subroutine extracts the break point information from a given file
	# it stores conflicts (wrt the qId) into conflictsQuery to then cross check with chimeric quality scores
	# it stores conflicts (wrt to xMapId\trId\tqId triplet)
	my ($file, $headerLinesRef, $conflictsQueryRef, $conflictsAllRef) = @_;
	open(IN, "$file") or dieLog "ERROR: cut_conflicts, getConflicts: cannot open file $file: $!\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^xMapId/)	{
			push(@$headerLinesRef, $line);
			next;
		} # if line
		# store the conflict information into the conflictsQuery hash and conflictsAll hash
		my @content = split(/\t/, $line);
		my ($xId, $refId, $refLeftBkpt, $refRightBkpt, $refAlignOrientation, $qryId, $qryLeftBkpt, $qryRightBkpt, $qryAlignOrientation) = ($content[0], $content[2], $content[3], $content[4], $content[5], $content[7], $content[8], $content[9], $content[10]);
		# change the bkpt to integers
		($refLeftBkpt, $refRightBkpt, $qryLeftBkpt, $qryRightBkpt) = (int($refLeftBkpt), int($refRightBkpt), int($qryLeftBkpt), int($qryRightBkpt));
		push(@{$conflictsAllRef->{"$xId\t$refId\t$qryId"}}, {refLeftBkpt => $refLeftBkpt, qryLeftBkpt => $qryLeftBkpt, 
			refRightBkpt => $refRightBkpt, qryRightBkpt => $qryRightBkpt,
			refLeftBkptStatus => "okay", refRightBkptStatus => "okay", refStatus => "okay",
			qryLeftBkptStatus => "okay", qryRightBkptStatus => "okay", qryStatus => "okay",
			refAlignOrientation => $refAlignOrientation, qryAlignOrientation => $qryAlignOrientation});
		
		# only store conflicts in conflictsQuery
		if ($qryLeftBkpt != -1)	{
			push(@{$conflictsQueryRef->{$qryId}}, {oriStart => $qryLeftBkpt, start => -1, end => -1, xId => $xId, refId => $refId, line => $line, neighbourCovRef => [], neighbourQScoreRef => [], toCut => 0});
		} # if leftBkpt
		if ($qryRightBkpt != -1)	{
			push(@{$conflictsQueryRef->{$qryId}}, {oriStart => $qryRightBkpt, start => -1, end => -1, xId => $xId, refId => $refId, line => $line, neighbourCovRef => [], neighbourQScoreRef => [], toCut => 0});
		} # if rightBkpt
	} # while line
	close IN;
	return ($headerLinesRef, $conflictsQueryRef, $conflictsAllRef);
} # getConflicts
