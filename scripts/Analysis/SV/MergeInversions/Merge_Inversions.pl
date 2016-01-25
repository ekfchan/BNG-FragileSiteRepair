#!/usr/bin/perl -w

use strict;
use warnings;

# Usage:
# /home/users/xzhou/tools/SV/MergeInversions/Merge_Inversions.pl <Merged_SMAP_DIR> <Merged_SMAP_FILE> <Overlap> <OUTPUT_DIR>
#
# Example:
# /home/users/xzhou/tools/SV/MergeInversions/Merge_Inversions.pl /home/users/xzhou/PERL/SV/merge_inversions/TEST/merged_smaps exp_refineFinal1_merged.smap 30000 /home/users/xzhou/PERL/SV/merge_inversions/TEST_out

die "Please enter in the appropriate input parameters: Merge_Inversions.pl <Merged_SMAP_DIR> <Merged_SMAP_FILE> <Overlap> <OUTPUT_DIR>\n" if (scalar(@ARGV) != 4);

my $homeDir = $ARGV[0];
my ($allInvBkptDir, $allInvBkptFile) = ($homeDir, $ARGV[1]);
my $outDir = $ARGV[3];	mkdir $outDir if (! -e $outDir);
my $outFile = $allInvBkptFile;	$outFile =~ s/\.smap$/_inversions.smap/;
my $logFile = $outFile;	$logFile =~ s/\.smap$/.log/;
my $offSetValue = $ARGV[2];	# a little bit of overlap allowed, give some offsets for alignments' positioning

open(LOG_HANDLE, ">$outDir/$logFile") or die "Cannot write to $outDir/$logFile: $!\n";
my ($headerLinesRef, $invBkptRef, $otherVarLinesRef, $relevantInvLinesRef) = readSmapInv($allInvBkptDir, $allInvBkptFile, \*LOG_HANDLE);

# listing all pair-able breakpoints on each chromosome
my ($invPairBkptsRef) = getInvPairs($invBkptRef);
# validate the pairings
$invPairBkptsRef = findValidInvPairs($invPairBkptsRef, $offSetValue);

# now deal with those that are not mergeable (and also those that are merged)
open(OUT, ">$outDir/$outFile") or die "Cannot write to $outDir/$outFile: $!\n";
my $count = 0;	my $countRef = \$count; print LOG_HANDLE "Main: printing to $outFile\n";

printHeaderLines(\*OUT, $headerLinesRef);
$countRef = printPairableBkpt(\*OUT, $invPairBkptsRef, $countRef);
$countRef = printNonPairableBkpt(\*OUT, $invPairBkptsRef, $relevantInvLinesRef, $countRef);
$countRef = printOtherVarLines(\*OUT, $otherVarLinesRef, $countRef);
print LOG_HANDLE "\twritten $$countRef lines\n";
close OUT;
close LOG_HANDLE;

sub printOtherVarLines	{
	my ($outRef, $otherVarLinesRef, $countRef) = @_;
	foreach my $otherVarLines (@$otherVarLinesRef)	{
		print $outRef "$otherVarLines\n";
$$countRef += 1;	
	} # foreach otherVarLines
	return $countRef;
} # printOtherVarLines

sub printNonPairableBkpt	{
	my ($outRef, $invPairBkptsRef, $relevantInvLinesRef, $countRef) = @_;
	# first determine among all possible pairings, which entries never pairs
	my %nrPairedEntries = ();
	foreach my $theKey (keys %$invPairBkptsRef)	{
		my ($rId, $iQId, $jQId, $iMatchId1, $iMatchId2, $jMatchId1, $jMatchId2) = split(/\t/, $theKey);
		for (my $i = 0; $i < scalar(@{$invPairBkptsRef->{$theKey}}); $i += 1)	{
			if ($invPairBkptsRef->{$theKey}[$i]{validPair} == 1)	{
				my $aRef = $invPairBkptsRef->{$theKey}[$i];
				my @iContent = split(/\t/, $aRef->{iLine});
				my @jContent = split(/\t/, $aRef->{jLine});
				my ($iSmapId, $iLinkId) = ($iContent[0], $iContent[12]);
				my ($jSmapId, $jLinkId) = ($jContent[0], $jContent[12]);
				my $iUniqueKey = ($iSmapId < $iLinkId) ? ("$iSmapId\t$iLinkId\t$rId\t$iQId") : ("$iLinkId\t$iSmapId\t$rId\t$iQId");	# assume smap id and linkId are integers, and put the smaller value first
				my $jUniqueKey = ($jSmapId < $jLinkId) ? ("$jSmapId\t$jLinkId\t$rId\t$jQId") : ("$jLinkId\t$jSmapId\t$rId\t$jQId");
	
				# record those that can be paired
				$nrPairedEntries{$iUniqueKey} = 1;
				$nrPairedEntries{$jUniqueKey} = 1;
			} # if validPair
		} # for i
	} # foreach theKey
	
	foreach my $entryKey (keys %$relevantInvLinesRef)	{
		# now iterates those relevant inversion lines from the original input file, and see which have no pairing, and hence be printed
		next if (exists $nrPairedEntries{$entryKey});
		foreach my $type (sort keys %{$relevantInvLinesRef->{$entryKey}})	{
			# assuming the type can be "inversion_XXX" and "inversion_XXX_partial"
			print $outRef "$relevantInvLinesRef->{$entryKey}{$type}\n";
$$countRef += 1;			
		} # foreach type
	} # foreach entryKey
	return $countRef;
} # sub printNonPairableBkpt

# pairedIInvBkptStartIndex => -1, pairedIInvBkptEndIndex => -1,

sub printPairableBkpt	{
	my ($outRef, $invPairBkptsRef, $countRef) = @_;
	foreach my $theKey (keys %$invPairBkptsRef)	{
		my @keyContent = split(/\t/, $theKey);
		my ($iQId, $jQId) = @keyContent[1..2];
		for (my $i = 0; $i < scalar(@{$invPairBkptsRef->{$theKey}}); $i += 1)	{
			my $aRef = $invPairBkptsRef->{$theKey}[$i];
			next if ($aRef->{validPair} == 0);	# not a valid pairing
			
			# print the line for the iEntry
			# only change is the breakpoint co-ordinates on the reference, and the linkId (assuming it's the column 13)
			my @iContent = split(/\t/, $aRef->{iLine});
			my @jContent = split(/\t/, $aRef->{jLine});
			
			my ($newIType, $newJType) = ("$iContent[9]_paired", "$jContent[9]_paired");
			print $outRef join("\t", @iContent[0..3], $aRef->{pairedIInvBkptQueryStart}, $aRef->{pairedIInvBkptQueryEnd}, $aRef->{pairedIInvBkptStart}, $aRef->{pairedIInvBkptEnd}, $iContent[8], $newIType, @iContent[10..11], $jContent[0])."\t";
			print $outRef join("\t", $aRef->{pairedIInvBkptQueryStartIndex}, $aRef->{pairedIInvBkptQueryEndIndex}, $aRef->{pairedIInvBkptStartIndex}, $aRef->{pairedIInvBkptEndIndex})."\t";	
			print $outRef join("\t", @iContent[17..$#iContent])."\n";
			# print the line for the jEntry
			print $outRef join("\t", @jContent[0..3], $aRef->{pairedJInvBkptQueryStart}, $aRef->{pairedJInvBkptQueryEnd}, $aRef->{pairedJInvBkptStart}, $aRef->{pairedJInvBkptEnd}, $jContent[8], $newJType, @jContent[10..11], $iContent[0])."\t";
			print $outRef join("\t", $aRef->{pairedJInvBkptQueryStartIndex}, $aRef->{pairedJInvBkptQueryEndIndex}, $aRef->{pairedJInvBkptStartIndex}, $aRef->{pairedJInvBkptEndIndex})."\t";
			print $outRef join("\t", @jContent[17..$#jContent])."\n";
$$countRef += 2;			
		} # for i
	} # foreach theKey
	return $countRef;
} # printPairableBkpt

sub printHeaderLines	{
	my ($outRef, $headerLinesRef) = @_;
	foreach my $headerLine (@$headerLinesRef)	{
		print $outRef "$headerLine\n";
	} # foreach headerLine
} # printHeaderLines

sub findValidInvPairs	{
	my ($invPairsRef, $offSetValue) = @_;
	foreach my $theKey (keys %$invPairsRef)	{
		my ($rId, $iQId, $jQId, $iMatchId1, $iMatchId2, $jMatchId1, $jMatchId2) = split(/\t/, $theKey);
		for (my $r = 0; $r < scalar(@{$invPairsRef->{$theKey}}); $r += 1)	{
			my $validPair = 0;

			my $pairsRef = $invPairsRef->{$theKey}[$r];
			### check reference configuration ###
			if ($pairsRef->{iOuterRefCoord} < $pairsRef->{iFlippedRefStart} && $pairsRef->{iFlippedRefStart} < $pairsRef->{iFlippedRefEnd})	{
				if ($pairsRef->{jFlippedRefStart} < $pairsRef->{jFlippedRefEnd} && $pairsRef->{jFlippedRefEnd} < $pairsRef->{jOuterRefCoord})	{
					# iOuterRefCoord, then inversion, then jOuterRefCoord
					# allowed config: iOuterRefCoord < jFlippedRefStart < iFlippedRefStart AND jFlippedRefEnd < iFlippedEnd < jOuterRefCoord
					if ( ($pairsRef->{iOuterRefCoord} - $offSetValue <= $pairsRef->{jFlippedRefStart} && $pairsRef->{jFlippedRefStart} <= $pairsRef->{iFlippedRefStart}) && 
						($pairsRef->{jFlippedRefEnd} <= $pairsRef->{iFlippedRefEnd} && $pairsRef->{iFlippedRefEnd} <= $pairsRef->{jOuterRefCoord} + $offSetValue) )	{

						### check query configuration ###
						$validPair = checkQueryValidPair($iQId, $pairsRef->{iFlippedQueryStart}, $pairsRef->{iFlippedQueryEnd}, $pairsRef->{iOuterQueryCoord}, 
							$jQId, $pairsRef->{jFlippedQueryStart}, $pairsRef->{jFlippedQueryEnd}, $pairsRef->{jOuterQueryCoord});	
					
						if ($validPair == 1)	{	
							# reference coordinates
							$pairsRef->{pairedIInvBkptStart} = $pairsRef->{iOuterRefCoord};
							$pairsRef->{pairedIInvBkptEnd} = $pairsRef->{jFlippedRefStart};
							$pairsRef->{pairedJInvBkptStart} = $pairsRef->{iFlippedRefEnd};
							$pairsRef->{pairedJInvBkptEnd} = $pairsRef->{jOuterRefCoord};

							# query coordinates - always match the outer-coordinates of the reference with the outer-coordinates of the query                                
							$pairsRef->{pairedIInvBkptQueryStart} =	$pairsRef->{iOuterQueryCoord};
							$pairsRef->{pairedIInvBkptQueryEnd} = (abs($pairsRef->{iOuterQueryCoord} - $pairsRef->{iFlippedQueryStart}) < abs($pairsRef->{iOuterQueryCoord} - $pairsRef->{iFlippedQueryEnd})) ? ($pairsRef->{iFlippedQueryStart}) : ($pairsRef->{iFlippedQueryEnd});
							$pairsRef->{pairedJInvBkptQueryStart} = (abs($pairsRef->{jOuterQueryCoord} - $pairsRef->{jFlippedQueryStart}) < abs($pairsRef->{jOuterQueryCoord} - $pairsRef->{jFlippedQueryEnd})) ? ($pairsRef->{jFlippedQueryStart}) : ($pairsRef->{jFlippedQueryEnd});
							$pairsRef->{pairedJInvBkptQueryEnd} = $pairsRef->{jOuterQueryCoord};

							# reference indecies
							$pairsRef->{pairedIInvBkptStartIndex} = $pairsRef->{iOuterRefCoordLabel};
							$pairsRef->{pairedIInvBkptEndIndex} = $pairsRef->{jFlippedRefStartLabel};
							$pairsRef->{pairedJInvBkptStartIndex} = $pairsRef->{iFlippedRefEndLabel};
							$pairsRef->{pairedJInvBkptEndIndex} = $pairsRef->{jOuterRefCoordLabel};

							# query indecies - always match the outer label index of the reference with the outer label index of the query
							$pairsRef->{pairedIInvBkptQueryStartIndex} = $pairsRef->{iOuterQueryCoordLabel};
							$pairsRef->{pairedIInvBkptQueryEndIndex} = (abs($pairsRef->{iOuterQueryCoordLabel} - $pairsRef->{iFlippedQueryStartLabel}) < abs($pairsRef->{iOuterQueryCoordLabel} - $pairsRef->{iFlippedQueryEndLabel})) ? ($pairsRef->{iFlippedQueryStartLabel}) : ($pairsRef->{iFlippedQueryEndLabel});
							$pairsRef->{pairedJInvBkptQueryStartIndex} = (abs($pairsRef->{jOuterQueryCoordLabel} - $pairsRef->{jFlippedQueryStartLabel}) < abs($pairsRef->{jOuterQueryCoordLabel} - $pairsRef->{jFlippedQueryEndLabel})) ? ($pairsRef->{jFlippedQueryStartLabel}) : ($pairsRef->{jFlippedQueryEndLabel});
							$pairsRef->{pairedJInvBkptQueryEndIndex} = $pairsRef->{jOuterQueryCoordLabel};

							
							$pairsRef->{validPair} = 1;
						} # if validPair
					} # if valid
				} # if pairsRef
			} # if pairsRef

			if ($pairsRef->{jOuterRefCoord} < $pairsRef->{jFlippedRefStart} && $pairsRef->{jFlippedRefStart} < $pairsRef->{jFlippedRefEnd})	{
				if ($pairsRef->{iFlippedRefStart} < $pairsRef->{iFlippedRefEnd} && $pairsRef->{iFlippedRefEnd} < $pairsRef->{iOuterRefCoord})	{
					# jOuterRefCoord, then inversion, then iOuterRefCoord
					# allowed config: jOuterRefCoord < iFlippedRefStart < jFlippedRefStart AND iFlippedRefEnd < jFlippedRefEnd < iOuterRefCoord
					if ( ($pairsRef->{jOuterRefCoord} - $offSetValue <= $pairsRef->{iFlippedRefStart} && $pairsRef->{iFlippedRefStart} <= $pairsRef->{jFlippedRefStart}) &&
						($pairsRef->{iFlippedRefEnd} <= $pairsRef->{jFlippedRefEnd} && $pairsRef->{jFlippedRefEnd} <= $pairsRef->{iOuterRefCoord} + $offSetValue) )	{

						### check query configuration ###
						$validPair = checkQueryValidPair($iQId, $pairsRef->{iFlippedQueryStart}, $pairsRef->{iFlippedQueryEnd}, $pairsRef->{iOuterQueryCoord}, 
							$jQId, $pairsRef->{jFlippedQueryStart}, $pairsRef->{jFlippedQueryEnd}, $pairsRef->{jOuterQueryCoord});	
					
						if ($validPair == 1)	{
							# reference coordinates
							$pairsRef->{pairedJInvBkptStart} = $pairsRef->{jOuterRefCoord};
							$pairsRef->{pairedJInvBkptEnd} = $pairsRef->{iFlippedRefStart};
							$pairsRef->{pairedIInvBkptStart} = $pairsRef->{jFlippedRefEnd};	
							$pairsRef->{pairedIInvBkptEnd} = $pairsRef->{iOuterRefCoord};
	
							# query coordinates - always match the outer-coordinates of the reference with the outer-coordinates of the query
							$pairsRef->{pairedJInvBkptQueryStart} = $pairsRef->{jOuterQueryCoord};
							$pairsRef->{pairedJInvBkptQueryEnd} = (abs($pairsRef->{jOuterQueryCoord} - $pairsRef->{jFlippedQueryStart}) < abs($pairsRef->{jOuterQueryCoord} - $pairsRef->{jFlippedQueryEnd})) ? ($pairsRef->{jFlippedQueryStart}) : ($pairsRef->{jFlippedQueryEnd});
							$pairsRef->{pairedIInvBkptQueryStart} = (abs($pairsRef->{iOuterQueryCoord} - $pairsRef->{iFlippedQueryStart}) < abs($pairsRef->{iOuterQueryCoord} - $pairsRef->{iFlippedQueryEnd})) ? ($pairsRef->{iFlippedQueryStart}) : ($pairsRef->{iFlippedQueryEnd});
							$pairsRef->{pairedIInvBkptQueryEnd} = $pairsRef->{iOuterQueryCoord};
			
							# reference indecies
							$pairsRef->{pairedJInvBkptStartIndex} = $pairsRef->{jOuterRefCoordLabel};
							$pairsRef->{pairedJInvBkptEndIndex} = $pairsRef->{iFlippedRefStartLabel};
							$pairsRef->{pairedIInvBkptStartIndex} = $pairsRef->{jFlippedRefEndLabel};
							$pairsRef->{pairedIInvBkptEndIndex} = $pairsRef->{iOuterRefCoordLabel};	

							# query indecies - always match the outer label index of the reference with the outer label index of the query
							$pairsRef->{pairedJInvBkptQueryStartIndex} = $pairsRef->{jOuterQueryCoordLabel};
							$pairsRef->{pairedJInvBkptQueryEndIndex} = (abs($pairsRef->{jOuterQueryCoordLabel} - $pairsRef->{jFlippedQueryStartLabel}) < abs($pairsRef->{jOuterQueryCoordLabel} - $pairsRef->{jFlippedQueryEndLabel})) ? ($pairsRef->{jFlippedQueryStartLabel}) : ($pairsRef->{jFlippedQueryEndLabel});
							$pairsRef->{pairedIInvBkptQueryStartIndex} = (abs($pairsRef->{iOuterQueryCoordLabel} - $pairsRef->{iFlippedQueryStartLabel}) < abs($pairsRef->{iOuterQueryCoordLabel} - $pairsRef->{iFlippedQueryEndLabel})) ? ($pairsRef->{iFlippedQueryStartLabel}) : ($pairsRef->{iFlippedQueryEndLabel});
							$pairsRef->{pairedIInvBkptQueryEndIndex} = $pairsRef->{iOuterQueryCoordLabel};
							
							$pairsRef->{validPair} = 1;
						} # if validPair
					} # if valid
				} # if pairsRef
			} # if pairsRef
		} # for r
	} # foreach theKey
	return $invPairsRef;
} # findValidInvPairs

sub checkQueryValidPair	{
	my ($iQId, $iFlippedQueryStart, $iFlippedQueryEnd, $iOuterQueryCoord, $jQId, $jFlippedQueryStart, $jFlippedQueryEnd, $jOuterQueryCoord) = @_;
	my $queryValidPair = 0;
	if ($iQId != $jQId)	{
		# if not the same query
		$queryValidPair = 1;
		return $queryValidPair;
	} # if qId

	# now one contig captures both breakpoints, must check the coordinates on the contig
	if ($iOuterQueryCoord < $iFlippedQueryStart && $iFlippedQueryStart < $iFlippedQueryEnd)	{
		if ($jFlippedQueryStart < $jFlippedQueryEnd && $jFlippedQueryEnd < $jOuterQueryCoord)	{
			# from left to right: iOuterQueryCoord, inverted region, jOuterQueryCoord
			if ($iFlippedQueryStart < $jFlippedQueryEnd)	{
				$queryValidPair = 1;
			} # if iFlippedQueryStart
		} # if jFlippedQueryStart
	} # if iOuterQueryCoord

	if ($jOuterQueryCoord < $jFlippedQueryStart && $jFlippedQueryStart < $jFlippedQueryEnd)	{
		if ($iFlippedQueryStart < $iFlippedQueryEnd && $iFlippedQueryEnd < $iOuterQueryCoord)	{
			# from left to right: jOuterQueryCoord, inverted region, iOuterQueryCoord
			if ($jFlippedQueryStart < $iFlippedQueryEnd)	{
				$queryValidPair = 1;
			} # if jFlippedQueryStart
		} # if iFlippedQueryStart
	} # if jOuterQuery

	return $queryValidPair;
} # checkQueryValidPair

sub getInvPairs	{
	my ($invBkptRef) = @_;
	my %invPairBkpts = ();

	foreach my $rId (keys %$invBkptRef)     {	
		# enumerate all possible pairing of inversion breakpoint along the chromosome, irrespective of whether the two entries have the same or different qId
		for (my $i = 0; $i < scalar(@{$invBkptRef->{$rId}}) - 1; $i += 1)       {
			my $iRef = $invBkptRef->{$rId}[$i];
			for (my $j = $i + 1; $j < scalar(@{$invBkptRef->{$rId}}); $j += 1)      {
				my $jRef = $invBkptRef->{$rId}[$j];
				my $theKey = "$rId\t$iRef->{qId}\t$jRef->{qId}\t$iRef->{matchId1}\t$iRef->{matchId2}\t$jRef->{matchId1}\t$jRef->{matchId2}";
				# store the pairs; there really should not be more than 1 entry in the array for each $theKey
				push(@{$invPairBkpts{$theKey}}, {iMatchId1 => $iRef->{matchId1}, iMatchId2 => $iRef->{matchId2},
					jMatchId1 => $jRef->{matchId1}, jMatchId2 => $jRef->{matchId2},
					
					iOuterRefCoord => $iRef->{outerRefCoord}, iFlippedRefStart => $iRef->{flippedRefStart}, iFlippedRefEnd => $iRef->{flippedRefEnd},
					jOuterRefCoord => $jRef->{outerRefCoord}, jFlippedRefStart => $jRef->{flippedRefStart}, jFlippedRefEnd => $jRef->{flippedRefEnd},
					iOuterQueryCoord => $iRef->{outerQueryCoord}, iFlippedQueryStart => $iRef->{flippedQueryStart}, iFlippedQueryEnd => $iRef->{flippedQueryEnd},
					jOuterQueryCoord => $jRef->{outerQueryCoord}, jFlippedQueryStart => $jRef->{flippedQueryStart}, jFlippedQueryEnd => $jRef->{flippedQueryEnd},
					
					iOuterRefCoordLabel => $iRef->{outerRefCoordLabel}, iFlippedRefStartLabel => $iRef->{flippedRefStartLabel}, iFlippedRefEndLabel => $iRef->{flippedRefEndLabel},
					jOuterRefCoordLabel => $jRef->{outerRefCoordLabel}, jFlippedRefStartLabel => $jRef->{flippedRefStartLabel}, jFlippedRefEndLabel => $jRef->{flippedRefEndLabel},
					iOuterQueryCoordLabel => $iRef->{outerQueryCoordLabel}, iFlippedQueryStartLabel => $iRef->{flippedQueryStartLabel}, iFlippedQueryEndLabel => $iRef->{flippedQueryEndLabel},
					jOuterQueryCoordLabel => $jRef->{outerQueryCoordLabel}, jFlippedQueryStartLabel => $jRef->{flippedQueryStartLabel}, jFlippedQueryEndLabel => $jRef->{flippedQueryEndLabel},

					pairedIInvBkptStart => -1, pairedIInvBkptEnd => -1,
					pairedJInvBkptStart => -1, pairedJInvBkptEnd => -1,
					pairedIInvBkptStartIndex => -1, pairedIInvBkptEndIndex => -1,
					pairedJInvBkptStartIndex => -1, pairedJInvBkptEndIndex => -1,
					iLine => $iRef->{line}, jLine => $jRef->{line},
					validPair => 0});	# the paired co-ordineates store potentially refined inversion breakpoint

			} # for j
		} # for i
	} # foreach rId
	return (\%invPairBkpts);
} # getInvPairs

sub readSmapInv	{
	my ($dir, $file, $logRef) = @_;
	my @headerLines = ();
	my %invBkpts = ();
	my @otherVarLines = ();
	my %relevantInvLines = ();

	open(IN, "$dir/$file") or die "readSmapInv: cannot open $dir/$file: $!\n";
my $count = 0;	print $logRef "readSmapInv: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;
		$line =~ s/\r//g;	
		if ($line =~ /^#/)	{
			# header line
			push(@headerLines, $line);
			next;
		} # if line

		$line =~ s/^\s+//;	$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		
		my ($smapId, $qId, $rId, $qStart, $qEnd, $start, $end, $type, $matchId1, $matchId2, $linkId, $qStartLabel, $qEndLabel, $startLabel, $endLabel) = ($content[0], $content[1], $content[2], $content[4], $content[5], $content[6], $content[7], $content[9], $content[10], $content[11], $content[12], $content[13], $content[14], $content[15], $content[16]);

		if ($type =~ /inversion/i)	{
			# inversion (various sub types)
			if ( ! ($type =~ /duplicate/i || $type =~ /overlap/i || $type =~ /repeat/i))	{
				# assuming the smapId and linkId are always positive integers (put the smaller of the two first)
				my $theKey = ($smapId < $linkId) ? ("$smapId\t$linkId\t$rId\t$qId") : ("$linkId\t$smapId\t$rId\t$qId");
				# record the line
				$relevantInvLines{$theKey}{$type} = $line;

				if ($type !~ /partial/i)	{
					# not the inferred line, record for pairing
					if ($qEnd < $qStart)	{
						# swap, but this should not occur for inversion entires
						my $temp = $qStart;	$qStart = $qEnd;	$qEnd = $temp;
					} # if qEnd
					if ($end < $start)	{
						# swap, but this should not occur for inversion entries
						my $temp = $start;	$start = $end;	$end = $temp;
					} # if end
					$invBkpts{$theKey} = {start => $start, end => $end, qStart => $qStart, qEnd => $qEnd, type => $type, qId => $qId, rId => $rId, matchId1 => $matchId1, matchId2 => $matchId2, smapId => $smapId, inferredRBkpt => -1, inferredQBkpt => -1, startLabel => $startLabel, endLabel => $endLabel, qStartLabel => $qStartLabel, qEndLabel => $qEndLabel, inferredRBkptLabel => -1, inferredQBkptLabel => -1, line => $line};
$count += 1;
				} else	{
					# inferred breakpoint line, take the qStart and start records for the inferred breakpoints
					if (! exists $invBkpts{$theKey})	{
						# if the main line has not been encountered yet
						$invBkpts{$theKey} = {start => -1, end => -1, qStart => -1, qEnd => -1, type => "", qId => $qId, rId => $rId, matchId1 => -1, matchId2 => -1, smapId => $smapId, inferredRBkpt => $start, inferredQBkpt => $qStart, startLabel => -1, endLabel => -1, qStartLabel => -1, qEndLabel => -1, inferredRBkptLabel => $startLabel, inferredQBkptLabel => $qStartLabel,  line => $line};
					} else	{
						# if the main line has already been encountered 
						($invBkpts{$theKey}{inferredRBkpt}, $invBkpts{$theKey}{inferredQBkpt}, $invBkpts{$theKey}{inferredRBkptLabel}, $invBkpts{$theKey}{inferredQBkptLabel}) = ($start, $qStart, $startLabel, $qStartLabel);
					} # if exists invBkpts
				} # if type partial
			} else	{
				# the other types (i.e. duplicate, overlap, etc) just put with the rest of SV, as they will not be participating in pairing breakpoints
				push(@otherVarLines, $line);
			} # if type duplicate
		} else	{
			# non-inversion variants
			push(@otherVarLines, $line);
		} # if type
	} # while line
	close IN;
print $logRef "\tread in $count relevant inversion records\n";

	my $refinedInvBkptsRef = refineInvBkpts(\%invBkpts, $logRef); 

	return (\@headerLines, $refinedInvBkptsRef, \@otherVarLines, \%relevantInvLines);
} # readSmapInv

sub refineInvBkpts	{
	my ($invBkptsRef, $logRef) = @_;
	my %refinedInvBkpts = ();
	# this subrountine figures out on the reference and query a) the non-flipped breakpoint, b) the flipped breakpoints
	foreach my $theKey (keys %$invBkptsRef)	{
		my $hRef = $invBkptsRef->{$theKey};
		if ($hRef->{start} == -1)	{
			print $logRef "ERROR: refineInvBkpts: inversion entry with smapId,linkId,refId,queryId=$theKey does not have the main inversion line\n";	next;
		} # if hRef start
		if ($hRef->{inferredRBkpt} == -1)	{
			print $logRef "ERROR: refineInvBkpts: inversion entry with smapId,linkId,refId,queryId=$theKey does not have partial inversion line\n";	next;
		} # if hRef inferredRBkpt
		
		# figure out which coordinates (on the reference and query) are the closest to the inferred bkpt, as these and the inferred bkpts denote the flipped region
		my ($outerRefCoord, $flippedRefStart, $flippedRefEnd, $outerRefCoordLabel, $flippedRefStartLabel, $flippedRefEndLabel, $refErrorFlag) = getFlipCoords($hRef->{start}, $hRef->{end}, $hRef->{inferredRBkpt}, $hRef->{startLabel}, $hRef->{endLabel}, $hRef->{inferredRBkptLabel}, $logRef);
		my ($outerQueryCoord, $flippedQueryStart, $flippedQueryEnd, $outerQueryCoordLabel, $flippedQueryStartLabel, $flippedQueryEndLabel, $queryErrorFlag) = getFlipCoords($hRef->{qStart}, $hRef->{qEnd}, $hRef->{inferredQBkpt}, $hRef->{qStartLabel}, $hRef->{qEndLabel}, $hRef->{inferredQBkptLabel}, $logRef);
		next if ($refErrorFlag == 1 || $queryErrorFlag == 1);

		# do not need to carry over the label information for the query, as that will not change
		push(@{$refinedInvBkpts{$hRef->{rId}}}, {
			outerRefCoord => $outerRefCoord, flippedRefStart => $flippedRefStart, flippedRefEnd => $flippedRefEnd, 
			outerQueryCoord => $outerQueryCoord, flippedQueryStart => $flippedQueryStart, flippedQueryEnd => $flippedQueryEnd, 
			type => $hRef->{type}, qId => $hRef->{qId}, matchId1 => $hRef->{matchId1}, matchId2 => $hRef->{matchId2}, smapId => $hRef->{smapId}, 
			outerRefCoordLabel => $outerRefCoordLabel, flippedRefStartLabel => $flippedRefStartLabel, flippedRefEndLabel => $flippedRefEndLabel, 
			outerQueryCoordLabel => $outerQueryCoordLabel, flippedQueryStartLabel => $flippedQueryStartLabel, flippedQueryEndLabel => $flippedQueryEndLabel,
			line => $hRef->{line}});
	} # foreach theKey
	return \%refinedInvBkpts;
} # refineInvBkpts

sub getFlipCoords	{
	my ($start, $end, $inferredBkpt, $startLabel, $endLabel, $inferredBkptLabel, $logRef) = @_;
	my ($outerCoord, $flippedStart, $flippedEnd) = (-1, -1, -1);    my ($outerCoordLabel, $flippedStartLabel, $flippedEndLabel) = (-1, -1, -1);
	my $errorFlag = 0;
	# note: the inferredBkpt is ALWAYS an coordinate of the inversion, and its closest neighbour is ALWAYS the other coordinate of the inversion
	if ($inferredBkpt < $start && $inferredBkpt < $end)	{
		# inferredBkpt is the smallest
		($flippedStart, $flippedEnd, $outerCoord) = ($start < $end) ? ($inferredBkpt, $start, $end) : ($inferredBkpt, $end, $start);
	} # inferredBkpt 
	if ($inferredBkpt > $start && $inferredBkpt > $end)	{
		# inferredBkpt is the largest
		($outerCoord, $flippedStart, $flippedEnd) = ($start < $end) ? ($start, $end, $inferredBkpt) : ($end, $start, $inferredBkpt);
	} # if inferredBkpt
	
	if ($inferredBkptLabel < $startLabel && $inferredBkptLabel < $endLabel)	{
		# inferredBkptLabel is the smallest
		($flippedStartLabel, $flippedEndLabel, $outerCoordLabel) = ($startLabel < $endLabel) ? ($inferredBkptLabel, $startLabel, $endLabel) : ($inferredBkptLabel, $endLabel, $startLabel);		
	} # inferredBkptLabel
	if ($inferredBkptLabel > $startLabel && $inferredBkptLabel > $endLabel)	{
		# inferredBkptLabel is the largest
		($outerCoordLabel, $flippedStartLabel, $flippedEndLabel) = ($startLabel < $endLabel) ? ($startLabel, $endLabel, $inferredBkptLabel) : ($endLabel, $startLabel, $inferredBkptLabel);
	} # if inferredBkptLabel

	if ($outerCoord == -1)	{
		print $logRef "ERROR: getFlipCoords: for start=$start, end=$end, inferredBkpt=$inferredBkpt triplet, it is impossible to determine which coordinates denote the flipped region, and which denotes the outer non-flip boundary\n";	$errorFlag = 1;
	} # if outerCoord
	if ($outerCoordLabel == -1)	{
		print $logRef "ERROR: getFlipCoords: for startLabel=$startLabel, endLabel=$endLabel, inferredBkptLabel=$inferredBkptLabel label triplet, it is impossible to determine which label denote the flipped region, and which denotes the outer non-flip boundary\n";	$errorFlag = 1;
	} # if outerCoordLabel
	
	return ($outerCoord, $flippedStart, $flippedEnd, $outerCoordLabel, $flippedStartLabel, $flippedEndLabel, $errorFlag);
} # getFlipCoords
