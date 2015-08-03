#!/usr/bin/perl

#Usage: perl scoreStitchPositionsBED.pl <stitchPositions.bed> <alignref_final_r.cmap> <alignref_final.xmap> <alignmol.xmap>
#Example: perl scoreStitchPositionsBED.pl EXP_REFINEFINAL1_fragileSiteRepaired_merged_fragileSiteRepaired_stitchPositions.withSeqs.bed EXP_REFINEFINAL1_fragileSiteRepaired_merged_fragileSiteRepaired_merged_r.cmap EXP_REFINEFINAL1_fragileSiteRepaired_merged_fragileSiteRepaired_merged.xmap EXP_REFINEFINAL1_fragileSiteRepaired_merged_fragileSiteRepaired_merged_alignmol.xmap

use strict;
use warnings;
use Cwd qw(abs_path);
use Data::Dumper;
use File::Basename;

#my $refId=$ARGV[2];

# Load input BED file with stitchPositions 
my $bedFileIn = $ARGV[0];
my $prefix = basename("$bedFileIn",".bed");
my $bedFileOut = $prefix."_scored.bed";
open BEDOUT, ">$bedFileOut" or die "ERROR: Could not open $bedFileOut: $!\n";
#print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n";

my @bedIn;
open BED, "<$bedFileIn" or die "ERROR: Could not open $bedFileIn: $!\n";
my $hasSeq=0;
while (my $line = <BED>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		#if ($s[0] == $refId) {
			my %bedLine = (
				"CMapId"  => "$s[0]",
				"Start" => "$s[1]", 
				"End"  => "$s[2]",
				"Type"  => "$s[3]",
				"Score" => "$s[4]",
				"Strand" => "$s[5]",
				"ThickStart" => "$s[6]",
				"ThickEnd" => "$s[7]",
				"ItemRgba" => "$s[8]",
				#"Sequence" => "$s[9]",
				#"line" => $line
			);
			# deal with BED file that has sequences
			if (scalar(@s) == 10) {
				$hasSeq=1;
				$bedLine{"Sequence"} = "$s[9]";		
			}			
			push @bedIn, \%bedLine;
		#}
	}
}
@bedIn =  sort { $a->{'CMapId'} <=> $b->{'CMapId'} || $a->{'Start'} <=> $b->{'Start'} || $a->{'End'} <=> $b->{'End'} } @bedIn;
print scalar(@bedIn)." stitch positions identified in $bedFileIn\n";
close BED; 


# Load alignref _r.cmap
my $rcmapFileIn = $ARGV[1];
my @rcmap;
open RCMAP, "$rcmapFileIn" or die "ERROR: Could not open $rcmapFileIn: $!\n";
while (my $line = <RCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	
		#if ($s[0] == $refId) {
			my %cmap_line = (
				"CMapId"  => "$s[0]",	#should be identical (only single contig in _r.cmap)
				"ContigLength" => "$s[1]",
				"NumSites"  => "$s[2]",
				"SiteID"  => "$s[3]",
				"LabelChannel"  => "$s[4]",
				"Position"  => "$s[5]",
				"StdDev" => "$s[6]",
				"Coverage" => "$s[7]",
				"Occurrence" => "$s[8]"
				);
			push @rcmap, \%cmap_line;	
		#}
	}
}
close RCMAP; 

# Load alignref XMAP
my $xmapFileIn = $ARGV[2];
my @xmap;
my $alignments;
open XMAP, "$xmapFileIn" or die "ERROR: Could not open $xmapFileIn: $!\n";
while (my $line = <XMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$alignments++;
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment
		#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 	float 	float 	int         	string   
		my %xmap_line = (
				"XmapEntryID"  => "$s[0]", 
				"QryContigID" => "$s[1]", 
				"RefContigID"  => "$s[2]",	
				"QryStartPos"  => "$s[3]", 
				"QryEndPos"  => "$s[4]", 
				"RefStartPos"  => "$s[5]", 
				"RefEndPos" => "$s[6]", 
				"Orientation" => "$s[7]", 
				"Confidence" => "$s[8]", 
				"HitEnum" => "$s[9]",
				"QryLen" => "$s[10]",
				"RefLen" => "$s[11]",
				"LabelChannel" => "$s[12]",
				"Alignment" => "$s[13]"
			); 
		push @xmap, \%xmap_line; }
}
print "Read in $alignments alignments from $xmapFileIn\n";

# Load merged alignmol XMAP
my $molxmapFileIn = $ARGV[3];
my @molxmap;
my $molalignments;
open MOLXMAP, "$molxmapFileIn" or die "ERROR: Could not open $molxmapFileIn: $!\n";
while (my $line = <MOLXMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		$molalignments++;
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment
		#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 	float 	float 	int         	string   
		my %xmap_line = (
				"XmapEntryID"  => "$s[0]", 
				"QryContigID" => "$s[1]", 
				"RefContigID"  => "$s[2]",	
				"QryStartPos"  => "$s[3]", 
				"QryEndPos"  => "$s[4]", 
				"RefStartPos"  => "$s[5]", 
				"RefEndPos" => "$s[6]", 
				"Orientation" => "$s[7]", 
				"Confidence" => "$s[8]", 
				"HitEnum" => "$s[9]",
				"QryLen" => "$s[10]",
				"RefLen" => "$s[11]",
				"LabelChannel" => "$s[12]",
				"Alignment" => "$s[13]"
			); 
		push @molxmap, \%xmap_line; }
}
print "Read in $molalignments alignments from $molxmapFileIn\n";

my $bedEntryCount=0;
my %bedStartLabelIds;
my %bedEndLabelIds;
my %genomeMapIds;
my %genomeMapStartLabelIds;
my %genomeMapEndLabelIds;

if ($hasSeq==1) {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n";
}
else {
	print BEDOUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\n";
}

foreach my $bedEntry_ref (@bedIn) {
	my %bedEntry = %$bedEntry_ref;
	#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
	my $bedCmapId = $bedEntry{'CMapId'};
	my $bedStart = int($bedEntry{'Start'});
	my $bedEnd = int($bedEntry{'End'});
	#print "$bedCmapId $bedStart $bedEnd\n";
	
	
	
	my $bedStartLabelId=0;
	my $bedEndLabelId=0;
	my $curPos=0;
	# Get alignref _r.cmap label IDs that correspond to nearst labels to the left of the start and right of the end  
	foreach my $rcmapEntry_ref (@rcmap) {
		my %rcmapEntry = %$rcmapEntry_ref;
		if ($rcmapEntry{'CMapId'} eq $bedCmapId) {
			$curPos = $rcmapEntry{'Position'};
			#print "\t$curPos\n";
			if ($curPos<=($bedStart+2)) {
				#print "\t$curPos\n";
				$bedStartLabelId = $rcmapEntry{'SiteID'};
			}
			elsif ($curPos>=($bedEnd-2)) {
				#print "\t$curPos\n";
				$bedEndLabelId = $rcmapEntry{'SiteID'};
				last;
			}
		}
	}
	#print "$bedCmapId $bedStart $bedStartLabelId $bedEnd $bedEndLabelId\n";
	$bedStartLabelIds{$bedEntryCount} = $bedStartLabelId;
	$bedEndLabelIds{$bedEntryCount} = $bedEndLabelId;
		
	
	# get the cooresponding genome map label IDs for each start/end position
	my $genomeMapId=0;
	foreach my $xmapEntry_ref (@xmap) {
	#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment
	#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 	float 	float 	int         	string   
		my %xmapEntry = %$xmapEntry_ref;
		if ($xmapEntry{'RefContigID'} eq $bedCmapId) {
			my $line = $xmapEntry{'Alignment'};
			my @s = split(/(?<=\))/, $line);
			my $refLabels = "";
			my $qryLabels = "";
			my %matchPairs;
			# process alignment doublet string
			foreach my $pair (@s) {
				#print "Doublet: $pair\n";
				my @match = $pair =~ /(\d+)/g;
				#print "Ref: $match[0] Qry: $match[1]\n";
				$refLabels = $refLabels."$match[0] ";
				$qryLabels = $qryLabels."$match[1] ";
				$matchPairs{$match[0]} = $match[1];
			}
			#print "RefLabels: $refLabels\n";
			#print "QryLabels: $refLabels\n";
			#print "";
			#print Dumper(\%matchPairs);
			# check to see if alignment contains both reference labels cooresponding to BED start and end
			if (($refLabels =~ /\b$bedStartLabelId\b/)) {
			if (($refLabels =~ /\b$bedEndLabelId\b/)) {
				$genomeMapId = $xmapEntry{'QryContigID'};
				#print "XmapEntryID:$xmapEntry{'RefContigID'} BED line: $bedEntryCount Matching genome map ID: $genomeMapId\n";
				$genomeMapStartLabelIds{$bedEntryCount} = $matchPairs{$bedStartLabelId};
				$genomeMapEndLabelIds{$bedEntryCount} = $matchPairs{$bedEndLabelId};
				
				#print "RefId: $bedCmapId QryId: $genomeMapId BEDstart: $bedStart RefStartLabel: $bedStartLabelId QryStartLabel: $genomeMapStartLabelIds{$bedEntryCount} BEDend: $bedEnd RefEndLabel: $bedEndLabelId QryEndLabel: $genomeMapEndLabelIds{$bedEntryCount}\n";
				
			
				my %molMatchHash;
				my @molMatchIds;
				my @molMatchConfs;
				my @perfectMatchIds;
				my %perfectMatchHash;
				#my $molMatchConfCulm=0;
				#my $molMatchCountCulm=0;
				my @allIds;
				my %confHash;
				foreach my $molxmapEntry_ref (@molxmap) {
				#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen	RefLen	LabelChannel	Alignment
				#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 	float 	float 	int         	string   
					my %molxmapEntry = %$molxmapEntry_ref;
					if ($molxmapEntry{'RefContigID'} eq $genomeMapId ) {
						#print "Looking at moleculeID $molxmapEntry{'QryContigID'} for $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount}\n";
						my $molline = $molxmapEntry{'Alignment'};
						my @mols = split(/(?<=\))/, $molline);
						my $molrefLabels = "";
						my $molqryLabels = "";
						my %molmatchPairs;
						
						# process alignment doublet string
						foreach my $molpair (@mols) {
							#print "Doublet: $pair\n";
							my @molmatch = $molpair =~ /(\d+)/g;
							#print "Ref: $match[0] Qry: $match[1]\n";
							$molrefLabels = $molrefLabels."$molmatch[0] ";
							$molqryLabels = $molqryLabels."$molmatch[1] ";
							$molmatchPairs{$molmatch[0]} = $molmatch[1];
						}
						#print "\tmolRefLabels: $molrefLabels\n";
						
						my $low = $genomeMapStartLabelIds{$bedEntryCount} - 50;
						my $high = $genomeMapEndLabelIds{$bedEntryCount} + 50;
						#for (my $i=$low; $i<$genomeMapEndLabelIds{$bedEntryCount}; $i++) {
						ILOOP: for (my $i=$low; $i<=$genomeMapStartLabelIds{$bedEntryCount}; $i++) {
							#print "STARTING NEW ILOOP\n";
							#for (my $j=$genomeMapStartLabelIds{$bedEntryCount}; $j<=$high; $j++) {
							my $repeatIcount=0;
								JLOOP: for (my $j=$genomeMapEndLabelIds{$bedEntryCount}; $j<=$high; $j++) {								
								if ($j != $i) {
								#print "\t\tLooking for labels $i and $j\n";
									if (($molrefLabels =~ /\b$i\b/) && ($molrefLabels =~ /\b$j\b/)) {									
										
										my $conf=1; #confidence match value
										#my $conf = $molxmapEntry{'Confidence'};
										
										if (($i ==  $genomeMapStartLabelIds{$bedEntryCount}) && ($j == $genomeMapEndLabelIds{$bedEntryCount})) {
											#push @perfectMatchIds, $molxmapEntry{'QryContigID'};
											#$conf = $molxmapEntry{'Confidence'};
											#$conf = $conf * 2;
											#$perfectMatchHash{$molxmapEntry{'QryContigID'}} = $conf * 10;
											if (exists $confHash{$molxmapEntry{'QryContigID'}}) {
												$confHash{$molxmapEntry{'QryContigID'}} = $confHash{$molxmapEntry{'QryContigID'}} + ($conf*2);
											} 
											else {$confHash{$molxmapEntry{'QryContigID'}} = ($conf*2);}
											push @allIds, $molxmapEntry{'QryContigID'};
											#$conf = 1;
											#$repeatIcount++;
										}
										else {
											#$conf = $molxmapEntry{'Confidence'};
											# my $multiplier = 1;
											my $iDist = abs($genomeMapStartLabelIds{$bedEntryCount} - $i);
											my $jDist = abs($j - $genomeMapEndLabelIds{$bedEntryCount});
											# if (! ($iDist==0 && $jDist==0)) {
												# $multiplier = $iDist + $jDist;
											# }
											$conf = $conf + $iDist + $jDist;				
											#if ($iDist< 1) { $iDist = 1;}
											#if ($jDist< 1) { $jDist = 1;}
											#my $multiplier=1;
											#my $multiplier = $iDist + $jDist;
											#if ($multiplier < 1) {$multiplier = 1;}
											#$conf = $conf * ($multiplier);
											if (exists $confHash{$molxmapEntry{'QryContigID'}}) {
												#$multiplier=1;
												#$conf = $conf * ($multiplier);
												#if ($repeatIcount < 99999999999) {
													#$molMatchHash{$molxmapEntry{'QryContigID'}} = $molMatchHash{$molxmapEntry{'QryContigID'}} + $conf;
													$confHash{$molxmapEntry{'QryContigID'}} = $confHash{$molxmapEntry{'QryContigID'}} + ($conf);
													$repeatIcount++;
													#print "\t\tMolecule Match repeatCounter: $repeatIcount\n";
												#}
												#else {
												#	#print "Maximum Matches for $molxmapEntry{'QryContigID'} i: $i $j reached repeatCounter: $repeatIcount\n";	
												#	next ILOOP;
												#}
											}
											else {
												#$multiplier=5;
												#$conf = $conf * ($multiplier);
												#$molMatchHash{$molxmapEntry{'QryContigID'}} = $conf;
												$confHash{$molxmapEntry{'QryContigID'}} = ($conf);
												#$repeatIcount++;
												#print "\tUniqueMolecule Match repeatCounter: $repeatIcount\n";
											}
											push @allIds, $molxmapEntry{'QryContigID'};
											#$molMatchHash{$molxmapEntry{'QryContigID'}} = $conf;
											#$molMatchConfCulm = $molMatchConfCulm + $conf;
											#$molMatchCountCulm++;
										}
										
										#print "BedLine: $bedEntryCount MATCH: $molxmapEntry{QryContigID} want: $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount} found: $i $j\n";
										
										#print "\t\tMATCH FOUND! $i $j\n";
										#$molMatchCount++;
										#print "\nLooking for: $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount} Found: $i $j\n";
										#print "molRefLabels: $molrefLabels\n";
										#push @molMatchIds, $molxmapEntry{'QryContigID'};
										#push @molMatchConfs, $molxmapEntry{'Confidence'};
										#$molMatchHash{$molxmapEntry{'QryContigID'}} = $conf;
										
										#$molMatchConf = $molMatchConf + $molxmapEntry{'Confidence'};
										#print "\tMolXmapEntry $molxmapEntry{'Alignment'}\t\tmolRefLabels: $molrefLabels\n";	
									}
									#else { print "\t\tMATCH NOT FOUND! $i $j\n"; }
								}
							}
						}
					}
				}

				
				my $molMatchCount=0;
				my $molMatchConf=0;
				foreach my $id (keys %molMatchHash) {
					#print "\t\tMolId: $id Conf: $molMatchHash{$id}\n";
					$molMatchCount++;
					$molMatchConf = $molMatchConf + $molMatchHash{$id};
				}
				
				# final method to rank matches to be determined
				#$molMatchConf = $molMatchConfCulm;
				
				my $perfectMatchCount=0;
				my $perfectMatchConf=0;
				foreach my $id (keys %perfectMatchHash) {
					#print "\t\tPerfectMolId: $id PerfectConf: $perfectMatchHash{$id}\n";
					$perfectMatchCount++;
					$perfectMatchConf = $perfectMatchConf + $perfectMatchHash{$id};
				}
				
				
				my $molMatchConfAvg=0;
				#$molMatchConf = scalar(unique(@molMatchConfs));				
				if (($molMatchConf != 0) && ($molMatchCount != 0)) {
					$molMatchConfAvg = $molMatchConf/$molMatchCount;
				}
				else { $molMatchConfAvg = 0;}
				
				my $perfectMatchConfAvg=0;
				if (($perfectMatchConf != 0) && ($perfectMatchCount != 0)) {
					$perfectMatchConfAvg = $perfectMatchConf/$perfectMatchCount;
				}
				else { $perfectMatchConfAvg = 0;}
				
				#my $totalHits = $molMatchCount + $perfectMatchCount;
				my $totalHits = scalar(@allIds);
				my $totalMol = scalar(unique(@allIds));
				#my $totalConf = $molMatchConf + $perfectMatchConf;
				
				my %ids;
				for (@allIds) {
					#print "Matched Molecule Ids:\n".join("\n\t",@allIds)."\n";
					$ids{$_}++;
				}
				
				# debug
				# if ($bedEntryCount == 16) {
					# for (unique(@allIds)) {
						# print "\n\tMatched Molecule Id: $_ TotalConf: $confHash{$_}\n";
					# }
				# }
				
				
				my $totalConf=0;
				foreach my $id (keys %confHash) {
					$totalConf = $totalConf + $confHash{$id};
				}
				
				#print "Count Hash:\n";
				#print "\n\t".Dumper(\%ids)."\n";
				#print "Confidence Hash:\n";
				#print "\n\t".Dumper(\%confHash)."\n";
				
				my $totalConfAvg=0;
				if (($totalConf != 0) && ($totalMol != 0)) {
					$totalConfAvg = $totalConf/$totalMol;
				}
				else { $totalConfAvg = 0;}
				
				#test
				#$totalMol = $totalHits;
				
				#my $perfectMatchs = scalar(unique(@perfectMatchIds));
				#my $logTotalConf = log10($totalConf);
				#my $logTotalConfAvg = log10($totalConfAvg);
				my $logTotalConf = log($totalConf);
				my $logTotalConfAvg = log($totalConfAvg);
				#my $logTotalConf = $totalConf;
				#my $logTotalConfAvg = $totalConfAvg;
				
				#print "\nBedLine: $bedEntryCount RefId: $bedCmapId QryId: $genomeMapId BEDstart: $bedStart RefStartLabel: $bedStartLabelId QryStartLabel: $genomeMapStartLabelIds{$bedEntryCount} BEDend: $bedEnd RefEndLabel: $bedEndLabelId QryEndLabel: $genomeMapEndLabelIds{$bedEntryCount} UniqueMoleculeCount: $molMatchCount MolCulmCount: $molMatchCountCulm PerfectMatches: $perfectMatchCount ConfCount: $molMatchConf PerfectMatchConf: $perfectMatchConf AvgConf: $molMatchConfAvg PerfectMatchConfAvg: $perfectMatchConfAvg TotalCount: $totalMol TotalConf: $totalConf log(TotalConf): $logTotalConf TotalConfAvg: $totalConfAvg log(TotalConfAvg): $logTotalConfAvg\n";
				
				print "\nBedLine: $bedEntryCount RefId: $bedCmapId Start-End: $bedStart-$bedEnd\n\tQryLabels: $genomeMapStartLabelIds{$bedEntryCount} $genomeMapEndLabelIds{$bedEntryCount} RefLabels: $bedStartLabelId $bedEndLabelId\n\tTotalHits: $totalHits UniqueMolecules: $totalMol TotalConf: $totalConf AvgConf: ".($totalConf/$totalMol)." AvgConf2 ".($totalConf/$totalHits)."\n";
				
				#my $score = (sprintf "%.3f", $logTotalConf);
				#my $score = (sprintf "%.2f", $totalConfAvg);
				my $score = int($totalConf);
				#my $score = int($totalConfAvg);
				my $log10Score = log10($score);
				my $logScore = log($score);
				my $score2 = log10($totalConf);
				my $score3 = $score*($totalMol/$totalHits);
				$score = (sprintf "%.3f", $logScore);
				#$score = (sprintf "%.2f", ($totalConf/$totalHits)); 
				print "\t\tRecommended score: $score log10(score): $log10Score log(score): $logScore AltScore: $score3\n";
				#" AltScore: ".int(($totalConf/$totalHits))." AltLog10Score: $score2 AltLogScore: $logTotalConf\n";	
				
				#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
				#my %bedEntry = %$bedEntry_ref;
				my $bedLineOut = "$bedCmapId\t$bedStart\t$bedEnd\t$bedEntry{Type}\t$score\t$bedEntry{Strand}\t$bedEntry{ThickStart}\t$bedEntry{ThickEnd}\t$bedEntry{ItemRgba}";
				if (defined $bedEntry{Sequence}) {
					$bedLineOut = $bedLineOut."\t$bedEntry{Sequence}";
				}
				print BEDOUT "$bedLineOut\n";
			}			
			}
		}		
	}	
	$bedEntryCount++;
	# if ($bedEntryCount>17) {
		# last; 
	# }
}


	



sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

sub log10 {
    my $n = shift;
	my $out=0;
    if ($n != 0) {
		$out = log($n)/log(10);
	}
	return $out;
}