#!/usr/bin/perl -w

use strict;
use warnings;

# Usage:
# /home/users/xzhou/tools/SV/MergeInversions/Merge_Inversions.pl <Merged_SMAP_DIR> <Merged_SMAP_FILE> <XMAP_DIR> <XMAP_FILE_PREFIX> <Overlap> <Flank_Align_Conf_Cutoff> <Inv_Align_Conf_Cutoff> <OUTPUT_DIR>
#
# Example:
# /home/users/xzhou/tools/SV/MergeInversions/Merge_Inversions.pl /home/users/xzhou/tools/Pipeline_2799_X/test/gasi/contigs/gasi_rFinal_sv/merged_smaps gasi_rFinal_merged.smap /home/users/xzhou/tools/Pipeline_2799_X/test/gasi/contigs/gasi_rFinal_sv gasi_rFinal_contig 30000 12 9 /home/users/xzhou/tools/Pipeline_2799_X/test/gasi/contigs/Merge_Inversions

my $homeDir = $ARGV[0];

my $allInvBkptDir = $homeDir;
my $allInvBkptFile = $ARGV[1];

my $invBkptDir = $homeDir;
my $invBkptFile = $allInvBkptFile;

my $xmapDir = $ARGV[2];
my $xmapFilePrefix = $ARGV[3];	my $xmapFileSuffix = "xmap";

my $outDir = $ARGV[7];	mkdir $outDir if (! -e $outDir);
my $outFile_same = "EXP_RFINAL_ALL.sameContig.inv";
my $outFile_diff = "EXP_RFINAL_ALL.diffContig.inv";
my $outFile = "EXP_RFINAL_ALL.inv";

my $offSetValue = $ARGV[4];	# a little bit of overlap allowed, give some off sets for alignments' positioning
my $flankAlignConfCutoff = $ARGV[5];
my $invAlignConfCutoff = $ARGV[6];

# read in all the inversion bkpt
my $invBkptRef = getInvBkpt_same($allInvBkptDir, $allInvBkptFile);

# put pairs of bkpt from the same contig
my ($invPairsRef, $nrMatchGroupIdsRef) = getInvPairsSameContig($invBkptRef);
$invPairsRef = calcValidInvPairs($invPairsRef, $nrMatchGroupIdsRef, $xmapDir, $xmapFilePrefix, $xmapFileSuffix, $offSetValue, $flankAlignConfCutoff, $invAlignConfCutoff);
printInv($outDir, $outFile_same, $invPairsRef);

# put pairs of bkpt from different contigs
my $invBkptsRef = getInvBkpts_diff($invBkptDir, $invBkptFile);
my ($invPairBkptsRef, $nrMatchIdsRef) = calcInvPairBkpts($invBkptsRef);
$invPairBkptsRef = findValidInvBkptPairs($invPairBkptsRef, $nrMatchIdsRef, $xmapDir, $xmapFilePrefix, $xmapFileSuffix, $offSetValue, $flankAlignConfCutoff, $invAlignConfCutoff);
printInvBkptResults($outDir, $outFile_diff, $invPairBkptsRef);

# Merge the output files
my $command = "type $outDir\\$outFile_same | findstr /b \"Ref\" > $outDir\\$outFile";
system($command);

$command = "type $outDir\\$outFile_same | findstr /v /b \"Ref\" >> $outDir\\$outFile";
system($command);

$command = "type $outDir\\$outFile_diff | findstr /v /b \"Ref\" >> $outDir\\$outFile";
system($command);


##############################################
#               Subroutines                  #
##############################################
sub printInv	{
	my ($dir, $file, $invPairsRef) = @_;
	open(OUT, ">$dir\\$file") or die "printInv: cannot write to $dir\\$file: $!\n";
print  "printInv: writing to $file\n";	my $count = 0;
	print OUT "Ref\tiContig\tiMatchId1\tiMatchId2\tiRefInvBkptStart\tiRefInvBkptEnd\tjContig\tjMatchId1\tjMatchId2\tjRefInvBkptStart\tjRefInvBkptEnd\tvalidPair\n";
	foreach my $rIdQId (keys %$invPairsRef)	{
		my ($rId, $qId) = split(/\t/, $rIdQId);
		for (my $i = 0; $i < scalar(@{$invPairsRef->{$rIdQId}}); $i += 1)	{
			my $aRef = $invPairsRef->{$rIdQId}[$i];
			print OUT join("\t", $rId, $qId, $aRef->{iMatchId1}, $aRef->{iMatchId2}, $aRef->{iInvBkptStart}, $aRef->{iInvBkptEnd}, $qId, $aRef->{jMatchId1}, $aRef->{jMatchId2}, $aRef->{jInvBkptStart}, $aRef->{jInvBkptEnd}, $aRef->{validPair})."\n";	
$count += 1;	
		} # for i
	} # foreach rIdQId
print  "\twritten $count records\n";
	close OUT;
} # printInv

sub calcValidInvPairs	{
	my ($invPairsRef, $nrMatchGroupIdsRef, $dir, $filePrefix, $fileSuffix, $offSetValue, $flankAlignConfCutoff, $invAlignConfCutoff) = @_;

	foreach my $rIdQId (keys %$invPairsRef)	{

		my %alignments = ();

		my ($rId, $qId) = split(/\t/, $rIdQId);
		my $alignmentsRef = getAlignments_same($dir, "$filePrefix$qId.$fileSuffix", $nrMatchGroupIdsRef, $rIdQId);

		for (my $r = 0; $r < scalar(@{$invPairsRef->{$rIdQId}}); $r += 1)	{
			my $pairsRef = $invPairsRef->{$rIdQId}[$r];
			my $bkptIAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId1}};
			my $bkptIAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId2}};
			my $bkptJAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId1}};
			my $bkptJAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId2}};

			# make sure that iAlign1Ref < iAlign2Ref, based on chromosome co-ordinates
			if ($bkptIAlign2Ref->{start} < $bkptIAlign1Ref->{start})	{
				# swap
				$bkptIAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId2}};
				$bkptIAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId1}};
				# assume bkptIAlign1Ref->{end} <= bkptJAlign2Ref->{start}
			} # if bkptIAlign2Ref
			# make sure that jAlign1Ref < jAlign2Ref, based on chromosome co-ordinates
			if ($bkptJAlign2Ref->{start} < $bkptJAlign1Ref->{start})	{
				# swap
				my $bkptJAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId2}};
				my $bkptJAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId1}};
			} # if bkptJAlign2Ref

			my $validPair = 1;
			# check to make sure that the config is correct
			if ($bkptIAlign1Ref->{start} < $bkptJAlign1Ref->{start})	{
				# allowed config: bkptIAlign1Ref < bkptJAlign1Ref and bkptIAlign2Ref < bkptJAlign2Ref
				if (!( $bkptIAlign1Ref->{end} - $offSetValue <= $bkptJAlign1Ref->{start} && $bkptIAlign2Ref->{end} <= $bkptJAlign2Ref->{start} + $offSetValue))	{
					$validPair = 0;
				} # if not valid
				# bkptJAlign1Ref must be in the inversion area of iAlignments AND bkptIAlign2Ref must be in the inversion area of jAlignments
				if (! ( ($pairsRef->{iInvBkptStart} - $offSetValue <= $bkptJAlign1Ref->{start} && $bkptJAlign1Ref->{start} <= $pairsRef->{iInvBkptEnd}) && ($pairsRef->{jInvBkptStart} <= $bkptIAlign2Ref->{end} && $bkptIAlign2Ref->{end} <= $pairsRef->{jInvBkptEnd} + $offSetValue) ))	{
					$validPair = 0
				} # if not valid
				# finally cut the inverted region alignments to have some slack
				if ( !( $flankAlignConfCutoff <= $bkptIAlign1Ref->{confidence} && $invAlignConfCutoff <= $bkptIAlign2Ref->{confidence} && $invAlignConfCutoff <= $bkptJAlign1Ref->{confidence} && $flankAlignConfCutoff <= $bkptJAlign2Ref->{confidence} ))	{
					$validPair = 0;
				} # if violated allowed alignment confidence
			} else	{
				# allowed config: bkptJAlign1Ref < bkptIAlign1Ref and bkptJAlign2Ref < bkptIAlign2Ref
				if (! ( $bkptJAlign1Ref->{end} - $offSetValue <= $bkptIAlign1Ref->{start} && $bkptJAlign2Ref->{end} <= $bkptIAlign2Ref->{start} + $offSetValue) )	{
					$validPair = 0;
				} # if not valid
				# bkptIAlign1Ref must be in the inversion area of jAlignments AND bkptJAlign2Ref must be in the inversion area of iAlignments
				if (! (($pairsRef->{jInvBkptStart} - $offSetValue <= $bkptIAlign1Ref->{start} && $bkptIAlign1Ref->{start} <= $pairsRef->{jInvBkptEnd}) && ($pairsRef->{iInvBkptStart} <= $bkptJAlign2Ref->{end} && $bkptJAlign2Ref->{end} <= $pairsRef->{iInvBkptEnd} + $offSetValue)))	{
					$validPair = 0;
				} # if not valid
				if ( !( $flankAlignConfCutoff <= $bkptJAlign1Ref->{confidence} && $invAlignConfCutoff <= $bkptJAlign2Ref->{confidence} && $invAlignConfCutoff <= $bkptIAlign1Ref->{confidence} && $flankAlignConfCutoff <= 
$bkptIAlign2Ref->{confidence}  ))	{
					$validPair = 0;
				} # if violated allowed alignment confidence
			} # if bkptIAlign1Ref


			### now make sure that the match groups positions on the query contig make sense
			### determine the order based on query co-ordinate
			if ($bkptIAlign2Ref->{qStart} < $bkptIAlign1Ref->{qStart})	{
				# swap
				$bkptIAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId2}};
				$bkptIAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId1}};
				# assume bkptIAlign1Ref->{qEnd} <= bkptIAlign2Ref->{qStart}
			} # bkptIAlign2Ref
			if ($bkptJAlign2Ref->{qStart} < $bkptJAlign1Ref->{qStart})	{
				 # swap
                                $bkptJAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId2}};
                                $bkptJAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId1}};
                                # assume that bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart}
			} # bkptJ

			if ($bkptIAlign1Ref->{qStart} < $bkptJAlign1Ref->{qStart})	{
				# allowed config bkptIAlign1Ref < bkptIAlign1Ref <= bkptJAlign1Ref < bkptJAlign2Ref
				if ( ! ($bkptIAlign1Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} && $bkptIAlign2Ref->{qStart} - $offSetValue <= $bkptJAlign1Ref->{qStart} && $bkptIAlign2Ref->{qEnd} <= $bkptJAlign2Ref->{qStart} + $offSetValue && $bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart}) ) {
                                        $validPair = 0;
                                } # if violated allowed config
			} else	{
				# allowed config bkptJAlign1Ref < bkptJAlign2Ref <= bkptIAlign1Ref < bkptIalign2Ref
				if ( ! ($bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart} && $bkptJAlign2Ref->{qStart} - $offSetValue <= $bkptIAlign1Ref->{qStart} && $bkptJAlign2Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} + $offSetValue && $bkptIAlign1Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} ) )        {
                                        $validPair = 0;
                                } # if violated allowed config
			} # if bkptIAlign1Ref
			$pairsRef->{validPair} = $validPair;

		} # for r
=head
		# iterate all inv bkpt pairs which are called by mapping the same contig to the same chromosome
		for (my $i = 0; $i < scalar(@{$invPairsRef->{$rIdQId}}); $i += 1)	{
			my $pairsRef = $invPairsRef->{$rIdQId}[$i];
			my $bkptIAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId1}};
			my $bkptIAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId2}};
			my $bkptJAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId1}};
			my $bkptJAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId2}};

			my $validPair = 1;

			### determine the order based on query co-ordinate
			if ($bkptIAlign2Ref->{qStart} < $bkptIAlign1Ref->{qStart})	{
				# swap
				$bkptIAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId2}};
				$bkptIAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{iMatchId1}};	
				# assume bkptIAlign1Ref->{qEnd} <= bkptIAlign2Ref->{qStart}
			} # bkptIAlign2Ref
			if ($bkptJAlign2Ref->{qStart} < $bkptJAlign1Ref->{qStart})	{
				# swap
				$bkptJAlign1Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId2}};
				$bkptJAlign2Ref = $alignmentsRef->{$rIdQId}{$pairsRef->{jMatchId1}};
				# assume that bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart}
			} # if bkptJ
			
			if ($bkptIAlign1Ref->{qStart} < $bkptJAlign1Ref->{qStart})	{
				# allowed config bkptIAlign1Ref < bkptIAlign1Ref <= bkptJAlign1Ref < bkptJAlign2Ref
				if ( ! ($bkptIAlign1Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} && $bkptIAlign2Ref->{qStart} - $offSetValue <= $bkptJAlign1Ref->{qStart} && $bkptIAlign2Ref->{qEnd} <= $bkptJAlign2Ref->{qStart} + $offSetValue && $bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart}) )	{
					$validPair = 0;
				} # if violated allowed config

				# now check the reference side, allowed config bkptJAlign1Ref->{start} contains within iInvBkpt area, and bkptIAlign2Ref->{end} contains within jInvBkpt area
				if ( ! ( ($pairsRef->{iInvBkptStart} - $offSetValue <= $bkptJAlign1Ref->{start} && $bkptJAlign1Ref->{start} <= $pairsRef->{iInvBkptEnd}) && ($pairsRef->{jInvBkptStart} <= $bkptIAlign2Ref->{end} && $bkptIAlign2Ref->{end} <= $pairsRef->{jInvBkptEnd} + $offSetValue) ) )	{
					$validPair = 0;
				} # if violated allowed config
			
				# finally cut the inverted region alignments to have some slack
				if ( !( $flankAlignConfCutoff <= $bkptIAlign1Ref->{confidence} && $invAlignConfCutoff <= $bkptIAlign2Ref->{confidence} && $invAlignConfCutoff <= $bkptJAlign1Ref->{confidence} && $flankAlignConfCutoff <= $bkptJAlign2Ref->{confidence} ) )	{
					$validPair = 0;
				} # if violated allowed alignment confidence
			} else	{
				# allowed config bkptJAlign1Ref < bkptJAlign2Ref <= bkptIAlign1Ref < bkptIAlign2Ref
				if ( ! ($bkptJAlign1Ref->{qEnd} <= $bkptJAlign2Ref->{qStart} && $bkptJAlign2Ref->{qStart} - $offSetValue <= $bkptIAlign1Ref->{qStart} && $bkptJAlign2Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} + $offSetValue && $bkptIAlign1Ref->{qEnd} <= $bkptIAlign2Ref->{qStart} ) )	{
					$validPair = 0;
				} # if violated allowed config

				# now check the reference side, allowed config bkptIAlignRef->{start} is within jInvBkpt area, and bkptJAlignRef->{end} is within IInvBkpt area
				if ( ! ( ($pairsRef->{jInvBkptStart} - $offSetValue <= $bkptIAlign1Ref->{start} && $bkptIAlign1Ref->{start} <= $pairsRef->{jInvBkptEnd}) && ($pairsRef->{iInvBkptStart} <= $bkptJAlign2Ref->{end} && $bkptJAlign2Ref->{end} <= $pairsRef->{iInvBkptEnd} + $offSetValue) ) )	{
print  "bkptIAlign1Ref->{qStart}=$bkptIAlign1Ref->{qStart};bkptIAlign1Ref->{qEnd}=$bkptIAlign1Ref->{qEnd};\nbkptIAlign1Ref->{start}=$bkptIAlign1Ref->{start};bkptIAlign1Ref->{end}=$bkptIAlign1Ref->{end};\n";
print  "bkptIAlign2Ref->{qStart}=$bkptIAlign2Ref->{qStart};bkptIAlign2Ref->{qEnd}=$bkptIAlign2Ref->{qEnd};\nbkptIAlign2Ref->{start}=$bkptIAlign2Ref->{start};bkptIAlign1Ref->{end}=$bkptIAlign2Ref->{end};\n";
print  "bkptJAlign1Ref->{qStart}=$bkptJAlign1Ref->{qStart};bkptJAlign1Ref->{qEnd}=$bkptJAlign1Ref->{qEnd};\nbkptJAlign1Ref->{start}=$bkptJAlign1Ref->{start};bkptJAlign1Ref->{end}=$bkptJAlign1Ref->{end};\n";
print  "bkptJAlign2Ref->{qStart}=$bkptJAlign2Ref->{qStart};bkptJAlign2Ref->{qEnd}=$bkptJAlign2Ref->{qEnd};\nbkptJAlign2Ref->{start}=$bkptJAlign2Ref->{start};bkptJAlign1Ref->{end}=$bkptJAlign2Ref->{end};\n";

print  "pairsRef->{jInvBkptStart}=$pairsRef->{jInvBkptStart};bkptIAlign1Ref->{start}=$bkptIAlign1Ref->{start};pairsRef->{jInvBkptEnd}=$pairsRef->{jInvBkptEnd};\npairsRef->{iInvBkptStart}=$pairsRef->{iInvBkptStart};bkptJAlign2Ref->{end}=$bkptJAlign2Ref->{end};pairsRef->{iInvBkptEnd}=$pairsRef->{iInvBkptEnd};\noffSetValue=$offSetValue;$rIdQId;\n";
					#$validPair = 0;
				} # if violated allowed config
				# finally cut the inverted region of alignments some slack
				if ( !( $flankAlignConfCutoff <= $bkptJAlign1Ref->{confidence} && $invAlignConfCutoff <= $bkptJAlign2Ref->{confidence} && $invAlignConfCutoff <= $bkptIAlign1Ref->{confidence} && $flankAlignConfCutoff <= $bkptIAlign2Ref->{confidence} ) )	{
					$validPair = 0;
				}
			} # bkptIAlign1Ref
		
		} # for i
=cut
	} # foreach rIdQId
	return $invPairsRef;
} # calcValidInvPairs

sub getAlignments_same	{
	my ($dir, $file, $nrMatchGroupIdsRef, $rIdQId) = @_;
	my %alignments = ();
	open(IN, "$dir\\$file") or die "getAlignments_same: cannot open $dir\\$file: $!\n";
my $count = 0;	print  "\t\tgetAlignments_same: reading in $file\n";
	while (my $line = <IN>)         {
		chomp $line;
		$line =~ s/\r//g; 
		next if ($line =~ /^#/);
		$line =~ s/^\s+//;      
		$line =~ s/\s+/\t/g;

		my @content = split(/\t/, $line);
		my ($xmapId, $aQId, $aRId, $qStart, $qEnd, $start, $end, $orientation, $confidence) = @content[0..8];

		next if (! exists $nrMatchGroupIdsRef->{$rIdQId}{$xmapId});    # dont read those alignments that do not contribute to any inversion breakpoint

		($qStart, $qEnd, $start, $end) = (int($qStart), int($qEnd), int($start), int($end));
		if ($qEnd < $qStart)    {
			# swap if qEnd is less than qStart
			my $temp = $qEnd;
			$qEnd = $qStart;
			$qStart = $temp;
		} # if qEnd

		if ($end < $start)      {
			my $temp = $end;
			$end = $start;
			$start = $temp;
		} # if end

		$alignments{$rIdQId}{$xmapId} = {start => $start, end => $end, qStart => $qStart, qEnd => $qEnd, orientation => $orientation, confidence => $confidence};
$count += 1;		
                } # while line
print  "\t\t\tread in $count records\n";
	close IN;
	return \%alignments;
} # getAlignments_same


sub getAlignments_diff	{
	my ($dir, $filePrefix, $fileSuffix, $nrMatchIdsRef) = @_;
	my %alignments = ();
	foreach my $qId (keys %$nrMatchIdsRef)	{
		# open the relevant contig alignments
		open(IN, "$dir/$filePrefix$qId.$fileSuffix") or die "getalignments_diff: cannot open $dir/$filePrefix$qId.$fileSuffix: $!\n";
my $count = 0;	print  "\tgetAlignments_diff: reading in $filePrefix$qId.$fileSuffix\n";
		while (my $line = <IN>)	{
			chomp $line;	
			$line =~ s/\r//g;
			next if ($line =~ /^#/);
			# remove leading white space
			$line =~ s/^\s+//;
			$line =~ s/\s+/\t/g;
			my @content = split(/\t/, $line);

			my ($xmapId, $qId, $rId, $qStart, $qEnd, $start, $end, $orientation, $confidence) = @content[0..8];
			
			# check if this match group is an inversion bkpt candidate
			next if (! exists $nrMatchIdsRef->{$qId}{$xmapId});
			# double check to make sure that the chromosome and contig id are right
			die "getAlignments_diff: alignment information is different from inv bkpt information\nxmap rId=$rId; qId=$qId\tinvBkpt alignment pair=$nrMatchIdsRef->{$qId}{$xmapId};\n" if ("$rId\t$qId" !~ /^$nrMatchIdsRef->{$qId}{$xmapId}$/);
			($qStart, $qEnd) = (int($qStart), int($qEnd));
			($start, $end) = (int($start), int($end));

			# make sure qStart and start are smaller than qEnd and end
			if ($qEnd < $qStart)	{
				my $temp = $qStart;
				$qStart = $qEnd;
				$qEnd = $temp;
			} # if qEnd
			if ($end < $start)	{
				my $temp = $start;
				$start = $end;
				$end = $temp;
			} # if end
			
			# store this alignment information
			$alignments{$qId}{$xmapId} = {chrom => $rId, start => $start, end => $end, qStart => $qStart, qEnd => $qEnd, orientation => $orientation, confidence => $confidence};
$count += 1;
		} # while line
print  "\tread in $count records\n";
		close IN;
	} # foreach qId
	return \%alignments;
} # getAlignments_diff

sub getInvPairsSameContig	{
	my ($invBkptRef) = @_;
	my %invPairs = ();	my %nrMatchGroupIds = ();

	my %invPairsSameContigMatchIds = ();
	foreach my $chrom (keys %$invBkptRef)	{
		my $rId = $chrom;
		for (my $i = 0; $i < scalar(@{$invBkptRef->{$chrom}}); $i += 1)	{
			my $iRef = $invBkptRef->{$chrom}[$i];
			my $qId = $iRef->{qId};
			# group them by the same queryId
			push(@{$invPairsSameContigMatchIds{$rId}{$qId}}, {matchId1 => $iRef->{matchId1}, invBkptStart => $iRef->{start}, invBkptEnd => $iRef->{end}, matchId2 => $iRef->{matchId2}});	
		} # for i
	} # foreach chrom	

	# now iterates every pair of inversion bkpt pairs
	foreach my $chrom (keys %invPairsSameContigMatchIds)	{
		foreach my $qId (keys %{$invPairsSameContigMatchIds{$chrom}})	{
			# pairing
			for (my $i = 0; $i < scalar(@{$invPairsSameContigMatchIds{$chrom}{$qId}}) - 1; $i += 1)	{
				# record all the matchgroup ids
				$nrMatchGroupIds{"$chrom\t$qId"}{$invPairsSameContigMatchIds{$chrom}{$qId}[$i]{matchId1}} = 1;
				$nrMatchGroupIds{"$chrom\t$qId"}{$invPairsSameContigMatchIds{$chrom}{$qId}[$i]{matchId2}} = 1;

				for (my $j = $i + 1; $j < scalar(@{$invPairsSameContigMatchIds{$chrom}{$qId}}); $j += 1)	{
					# record all the matchgroup ids
					$nrMatchGroupIds{"$chrom\t$qId"}{$invPairsSameContigMatchIds{$chrom}{$qId}[$j]{matchId1}} = 1;
					$nrMatchGroupIds{"$chrom\t$qId"}{$invPairsSameContigMatchIds{$chrom}{$qId}[$j]{matchId2}} = 1;

					# get all pairs
					push(@{$invPairs{"$chrom\t$qId"}}, {iMatchId1 => $invPairsSameContigMatchIds{$chrom}{$qId}[$i]{matchId1}, iMatchId2 => $invPairsSameContigMatchIds{$chrom}{$qId}[$i]{matchId2}, iInvBkptStart => $invPairsSameContigMatchIds{$chrom}{$qId}[$i]{invBkptStart}, iInvBkptEnd => $invPairsSameContigMatchIds{$chrom}{$qId}[$i]{invBkptEnd} , jMatchId1 => $invPairsSameContigMatchIds{$chrom}{$qId}[$j]{matchId1}, jMatchId2 => $invPairsSameContigMatchIds{$chrom}{$qId}[$j]{matchId2}, jInvBkptStart => $invPairsSameContigMatchIds{$chrom}{$qId}[$j]{invBkptStart}, jInvBkptEnd => $invPairsSameContigMatchIds{$chrom}{$qId}[$j]{invBkptEnd}, validPair => 0});

				} # for j
			} # for i
		} # foreach qId
	} # foreach chrom

	return (\%invPairs, \%nrMatchGroupIds);
} # getInvPairsSameContig

sub getInvBkpt_same	{
	my ($dir, $file) = @_;
	my %invBkpt = ();
	open(IN, "$dir\\$file") or die "getInvBkpt_same: cannot read in $dir\\$file: $!\n";
my $count = 0;	print  "getInvBkpt_same: reading in $file\n";
	my $all_count = 0;
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r//g;
		next if ($line =~ /^#/);
		$all_count += 1;
		$line =~ s/^\s+//;	$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		my ($qId, $rId, $qStart, $qEnd, $start, $end, $type, $matchId1, $matchId2) = ($content[1], $content[2], int($content[4]), int($content[5]), int($content[6]), int($content[7]), $content[9], $content[10], $content[11]);
		next if ($type !~ /inv/i);

		# make sure start < end
		if ($end < $start)	{
			my $temp = $start;
			$start = $end ;
			$end = $temp;
		} # end
		push(@{$invBkpt{$rId}}, {start => $start, end => $end, qStart => $qStart, qEnd => $qEnd, qId => $qId, type => $type, matchId1 => $matchId1, matchId2 => $matchId2});
$count += 1;
	} # whiel line
print  "\tread in $all_count total records\n";
print  "\tread in $count inversion records\n";
	close IN;
	return \%invBkpt;
} # getInvBkpt_same

##########################################

sub printInvBkptResults	{
	my ($dir, $file, $invPairBkptsRef) = @_;

	open(OUT, ">$dir\\$file") or die "printInvBkptResults: cannot write to $dir\\$file: $!\n";
print  "printInvBkptResults: writing to $file\n";	my $count = 0;
	print OUT join("\t", "Ref", "iContig", "iMatchId1", "iMatchId2", "iRefInvBkptStart", "iRefInvBkptEnd", "jContig", "jMatchId1", "jMatchId2", "jRefInvBkptStart", "jRefInvBkptEnd", "validPair")."\n";
	foreach my $chromQId (keys %$invPairBkptsRef)	{
		for (my $i = 0; $i < scalar(@{$invPairBkptsRef->{$chromQId}}); $i += 1)	{
			my $aRef = $invPairBkptsRef->{$chromQId}[$i];
			my ($chrom, $iQId, $jQId) = split(/\t/, $chromQId);
			print OUT join("\t", $chrom, $iQId, $aRef->{iMatchId1}, $aRef->{iMatchId2}, $aRef->{iInvBkptStart}, $aRef->{iInvBkptEnd}, $jQId, $aRef->{jMatchId1}, $aRef->{jMatchId2}, $aRef->{jInvBkptStart}, $aRef->{jInvBkptEnd}, $aRef->{validPair})."\n";
$count += 1;
		} # for i
	} # foreach chromQId
print  "\twritten $count records\n";
	close OUT;
} # printInvBkptResults

sub findValidInvBkptPairs	{
	my ($invPairBkptsRef, $nrMatchIdsRef, $dir, $filePrefix, $fileSuffix, $offSetValue, $flankAlignConfCutoff, $invAlignConfCutoff) = @_;

	my $alignmentsRef = getAlignments_diff($dir, $filePrefix, $fileSuffix, $nrMatchIdsRef);

	foreach my $chromQId (keys %$invPairBkptsRef)	{
		my ($chrom, $iQId, $jQId) = split(/\t/, $chromQId);
		for (my $r = 0; $r < scalar(@{$invPairBkptsRef->{$chromQId}}); $r += 1)	{
			my $iAlign1Ref = $alignmentsRef->{$iQId}{$invPairBkptsRef->{$chromQId}[$r]{iMatchId1}};
                	my $iAlign2Ref = $alignmentsRef->{$iQId}{$invPairBkptsRef->{$chromQId}[$r]{iMatchId2}};
                	my $jAlign1Ref = $alignmentsRef->{$jQId}{$invPairBkptsRef->{$chromQId}[$r]{jMatchId1}};
                	my $jAlign2Ref = $alignmentsRef->{$jQId}{$invPairBkptsRef->{$chromQId}[$r]{jMatchId2}};

                	# make sure that iAlign1Ref < iAlign2Ref, based on chromosome co-ordinates
                	if ($iAlign2Ref->{start} < $iAlign1Ref->{start})        {
                        	# swap
                        	$iAlign1Ref = $alignmentsRef->{$iQId}{$invPairBkptsRef->{$chromQId}[$r]{iMatchId2}};;
                        	$iAlign2Ref = $alignmentsRef->{$iQId}{$invPairBkptsRef->{$chromQId}[$r]{iMatchId1}};
                        	# assume that iAlign1Ref->{end} <= iAlign2Ref->{start}
                	} # if iAlign2Ref
                	# make sure that jAlign1Ref < jAlign2Ref, based on chromosome co-ordinates
                	if ($jAlign2Ref->{start} < $jAlign1Ref->{start})        {
                        	# swap
                        	$jAlign1Ref = $alignmentsRef->{$jQId}{$invPairBkptsRef->{$chromQId}[$r]{jMatchId2}};
                        	$jAlign2Ref = $alignmentsRef->{$jQId}{$invPairBkptsRef->{$chromQId}[$r]{jMatchId1}};
                        	# assume that jAlign1Ref->{end} <= jAlign2Ref->{start}
                	} # if jAlign2Ref

			my $validPair = 1;
                	# check to make sure that the config is correct
                	if ($iAlign1Ref->{start} < $jAlign1Ref->{start})        {
                        	# allowed config: iAlign1Ref < jAlign1Ref and iAlign2Ref < jAlign2Ref
				if ( ! (  $iAlign1Ref->{end} - $offSetValue <= $jAlign1Ref->{start} && $iAlign2Ref->{end} <= $jAlign2Ref->{start} + $offSetValue) )	{
                                	$validPair = 0;
                        	} # if not valid
                        	# jAlign1Ref must be in the inversion area of iAlignments AND iAlign2Ref must be in the inversion area of jAlignments
				if (! (($invPairBkptsRef->{$chromQId}[$r]{iInvBkptStart} - $offSetValue <= $jAlign1Ref->{start} && $jAlign1Ref->{start} <= $invPairBkptsRef->{$chromQId}[$r]{iInvBkptEnd}) && ($invPairBkptsRef->{$chromQId}[$r]{jInvBkptStart} <= $iAlign2Ref->{end} && $iAlign2Ref->{end} <= $invPairBkptsRef->{$chromQId}[$r]{jInvBkptEnd} + $offSetValue)) )	{
                                	$validPair = 0;
                        	} # if not valid
				# finally require the flanking region alignment to  be of higher confidence
				if (! ($flankAlignConfCutoff <= $iAlign1Ref->{confidence} && $invAlignConfCutoff <= $iAlign1Ref->{confidence} && $invAlignConfCutoff <= $jAlign1Ref->{confidence} && $flankAlignConfCutoff <= $jAlign2Ref->{confidence}) )	{
					$validPair = 0;
				} # if not valid
                	} else  {
                        	# allowed config: jAlign1Ref < iAlign2Ref and jAlign2Ref < iAlign2Ref
				if (! ($jAlign1Ref->{end} - $offSetValue <= $iAlign1Ref->{start} && $jAlign2Ref->{end} <= $iAlign2Ref->{start} + $offSetValue ) )	{
                                	$validPair = 0;
                        	} # if not valid
                        	# iAlign2Ref must be in the inversion area of jAlignment AND jAlign2Ref must be in the inversion area of iAlignments
				if ( ! (($invPairBkptsRef->{$chromQId}[$r]{jInvBkptStart} - $offSetValue <= $iAlign1Ref->{start} && $iAlign1Ref->{start} <= $invPairBkptsRef->{$chromQId}[$r]{jInvBkptEnd}) && ($invPairBkptsRef->{$chromQId}[$r]{iInvBkptStart} <= $jAlign2Ref->{end} && $jAlign2Ref->{end} <= $invPairBkptsRef->{$chromQId}[$r]{iInvBkptEnd} + $offSetValue) ))	{
                                	$validPair = 0;
                        	} # if not valid
				# finally require the flanking region alignments to be of high confidence
				if (! ($flankAlignConfCutoff <= $jAlign1Ref->{confidence} && $invAlignConfCutoff <= $jAlign2Ref->{confidence} && $invAlignConfCutoff <= $iAlign1Ref->{confidence} && $flankAlignConfCutoff <= $iAlign2Ref->{confidence}) )	{
					$validPair = 0;
				} # if not valid
                	} # if iAlign1Ref

                	$invPairBkptsRef->{$chromQId}[$r]{validPair} = $validPair;
		} # for r
	} # foreach chromQId

	return $invPairBkptsRef;
} # findValidInvBkptPairs


sub calcInvPairBkpts	{
	my ($invBkptsRef) = @_;
	my %invPairBkpts = ();
	my %nrMatchIds = ();

	foreach my $chrom (keys %$invBkptsRef)	{
		for (my $i = 0; $i < scalar(@{$invBkptsRef->{$chrom}}) - 1; $i += 1)	{
			my $iRef = $invBkptsRef->{$chrom}[$i];

			for (my $j = $i + 1; $j < scalar(@{$invBkptsRef->{$chrom}}); $j += 1)	{
				my $jRef = $invBkptsRef->{$chrom}[$j];

				next if ($iRef->{qId} =~ /^$jRef->{qId}$/);	# skip those inversion breakpoint pairs that are called by the same contig aligning to the same chromosome (this is done in previous script)

				my $theKey = "$chrom\t$iRef->{qId}\t$jRef->{qId}";

				# store the pairs
				push(@{$invPairBkpts{$theKey}}, {iMatchId1 => $iRef->{matchId1}, iMatchId2 => $iRef->{matchId2}, jMatchId1 => $jRef->{matchId1}, jMatchId2 => $jRef->{matchId2}, iInvBkptStart => $iRef->{start}, iInvBkptEnd => $iRef->{end}, jInvBkptStart => $jRef->{start}, jInvBkptEnd => $jRef->{end}, validPair => 0});

				# store all the necessary match ids
				$nrMatchIds{$iRef->{qId}}{$iRef->{matchId1}} = "$chrom\t$iRef->{qId}";
				$nrMatchIds{$iRef->{qId}}{$iRef->{matchId2}} = "$chrom\t$iRef->{qId}";
				$nrMatchIds{$jRef->{qId}}{$jRef->{matchId1}} = "$chrom\t$jRef->{qId}";
				$nrMatchIds{$jRef->{qId}}{$jRef->{matchId2}} = "$chrom\t$jRef->{qId}";
			} # for j
		} # for i
	} # foreach chrom

	return (\%invPairBkpts, \%nrMatchIds);
} # calcInvPairBkpts

sub getInvBkpts_diff	{
	my ($dir, $file) = @_;
	my %invBkpts = ();
	open(IN, "$dir\\$file") or die "getInvBkpts_diff: cannot read in $dir\\$file: $!\n";
my $count = 0;	print  "getInvBkpts_diff: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r//g;
		next if ($line =~ /^#/);

		$line =~ s/^\s+//;
		$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);

		my ($qId, $rId, $qStart, $qEnd, $start, $end, $type, $matchId1, $matchId2) = ($content[1], $content[2], int($content[4]), int($content[5]), int($content[6]), int($content[7]), $content[9], $content[10], $content[11]); 

		next if ($type !~ /inv/i);

		if ($qEnd < $qStart)	{
			# swap
			my $temp = $qStart;
			$qStart = $qEnd;
			$qEnd = $temp;
		} # qEnd

		if ($end < $start)	{
			# swap
			my $temp = $start;
			$start = $end;
			$end = $temp;
		} # if end

		push(@{$invBkpts{$rId}}, {start => $start, end => $end, qStart => $qStart, qEnd => $qEnd, type => $type, qId => $qId, matchId1 => $matchId1, matchId2 => $matchId2});
$count += 1;
	} # while line
print  "\tread in $count records\n";
	close IN;
	return \%invBkpts;
} # getInvBkpts_diff


