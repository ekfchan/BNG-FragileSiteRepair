#!/usr/bin/perl

#Usage: perl mergeContigs.pl <input CMAP> <first_contig_ID> <+/- orientation> <second_contig_ID> <+/- orientation> <N bases to add/delete in-between> <bp from end of first contig to insert an extra label> <labels distance [-1 = start/end on same label maxValue=1]> <firstQryConf> <secondQryConf>
#Designed to be used with non-overlapping AND overlapping contigs

# 0=====>N0=====>N == + +
# N<=====0N<=====0 == - -
# 0=====>NN<=====0 == + -
# N<=====00=====>N == - +

#Designed to work with CMAP 0.1

use strict;
use warnings;
use Data::Dumper;
sub getMissedLabelHash;

#print "\n";
#print qx/ps -o args $$/;
#print "\n";

open FILE, "$ARGV[0]" or die $!;
my $outName = $ARGV[0]; 
$outName =~ s/_q.cmap//g;
$outName =~ s/_merged_contigs//g; 
$outName =~ s/_merged_contigs_q.cmap//g; 
 
$outName = $outName."_merged_contigs_q.cmap";
my $outName_orig = $outName;
#$outName = $outName."_q.cmap";
$outName = $outName."_.temp";
open OUT, ">$outName" or die $!;

my $doOutput = 1;

my $idOffset = 1000000000;
my @cmap_out;		#references of hash table of cmap to print out
my @firstContig;	#reference of hash table of cmaps whose contig is <first_contig_ID>
my @secondContig;	#reference of hash table of cmaps whose contig is <second_contig_ID>
my @mergedContig;

my $missedLabelPadding=0;
if (defined $ARGV[6] && $ARGV[6]>0) {
	$missedLabelPadding = $ARGV[6];
}

my $labelsDistance = $ARGV[7] + 1;
if ($labelsDistance<0 && $ARGV[5]<0) {
	$labelsDistance = $ARGV[7];
}

my $firstQryConf = $ARGV[8];
my $secondQryConf = $ARGV[9];

my $excludeFile = $ARGV[10];
open (EXCLUDE, '>>', $excludeFile) or die "ERROR; $!\n";

#read input CMAP
while (my $line = <FILE>) {
	chomp($line);
	#if header then print directly out to new CMAP
	if ($line =~ /^#/) {
		print OUT "$line\n"; }
	else {
		$line =~ s/\r//g;
		my @s = split(/\t/,$line);
		if ($s[0] eq $ARGV[1]) {
			#print "Adding contig $s[0] to firstContig\n";
			#load first contig into hash
			##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	GmeanSNR	lnSNRsd	SNR
			#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	12.4650	0.4814	0
			my %cmap_line = (
				"CMapId"  => shift @s, 
				"ContigLength" => shift @s, 
				"NumSites"  => shift @s, 
				"SiteID"  => shift @s, 
				"LabelChannel"  => shift @s, 
				"Position"  => shift @s, 
				"StdDev" => shift @s, 
				"Coverage" => shift @s, 
				"Occurrence" => shift @s, 
				"TheRest" => join("\t", @s)
			);
			push @firstContig, \%cmap_line;	next; }
		elsif ($s[0] eq $ARGV[3]) {
			my %cmap_line = (
				"CMapId"  => shift @s, 
				"ContigLength" => shift @s, 
				"NumSites"  => shift @s, 
				"SiteID"  => shift @s, 
				"LabelChannel"  => shift @s, 
				"Position"  => shift @s, 
				"StdDev" => shift @s, 
				"Coverage" => shift @s, 
				"Occurrence" => shift @s, 
				"TheRest" => join("\t", @s)
			);
			push @secondContig, \%cmap_line; next; }
		else {
			my %cmap_line = (
				"CMapId"  => shift @s, 
				"ContigLength" => shift @s, 
				"NumSites"  => shift @s, 
				"SiteID"  => shift @s, 
				"LabelChannel"  => shift @s, 
				"Position"  => shift @s, 
				"StdDev" => shift @s, 
				"Coverage" => shift @s, 
				"Occurrence" => shift @s, 
				"TheRest" => join("\t", @s)
			);
			push @cmap_out, \%cmap_line; next; }
		}
}

#merge firstContig and secondContig taking into account orientation
	
	my $mergedId = $ARGV[1] + $idOffset;

	my $firstContigSites = scalar(@firstContig) - 1;	#number of first contigs (elements in @firstContig)
	my $firstContigSitesOrig = $firstContig[$firstContigSites]->{'NumSites'};
	my $firstContigSitesStartIdx = 0; 
	my $firstContigSitesEndIdx = $firstContig[$firstContigSites]->{'NumSites'};
	#my $firstContigStart = $firstContig[0]->{'Position'} - 20;
	#my $firstContigStartOrig = $firstContig[0]->{'Position'} - 20;
	my $firstContigStart = 0;
	my $firstContigStartOrig = 0;
	my $firstContigEnd = $firstContig[$firstContigSites]->{'ContigLength'};	#assume the contigs are sorted
	my $firstContigEndOrig = $firstContig[$firstContigSites]->{'ContigLength'};
	
	
	my $secondContigSites = scalar(@secondContig) - 1;	#number of second contigs (elements in @secondContig)
	my $secondContigSitesOrig = $secondContig[$secondContigSites]->{'NumSites'};
	my $secondContigSitesStartIdx = 0;
	my $secondContigSitesEndIdx = $secondContig[$secondContigSites]->{'NumSites'};
	#my $secondContigStart = $secondContig[0]->{'Position'} - 20;
	#my $secondContigStartOrig = $secondContig[0]->{'Position'} - 20;
	my $secondContigStart = 0;
	my $secondContigStartOrig = 0;
	my $secondContigEnd = $secondContig[$secondContigSites]->{'ContigLength'};
	my $secondContigEndOrig = $secondContig[$secondContigSites]->{'ContigLength'};

	my $positionOffset = ($firstContigEnd - $firstContigStart) + abs($ARGV[5]) - 1;

	# deal with overlapping contigs
	my $sitesOffset = 0; 
	my $posOffset = 0;
	my $reverseOffset = 0;
	my $reverseOffset2 = 0;
	my $forwardOffset = 0;
	my $forwardPosOffset = 0;
	my $reversePosOffset = 0;
	my $forwardOffset2 = 0;
	if ($labelsDistance < 0 && $ARGV[5]<0) {
		#$mergedSites = $firstContigSites + $secondContigSites + $labelsDistance;
		$sitesOffset = abs($labelsDistance);
		$posOffset = abs($ARGV[5]);

		
		# if first contig conf < second contig conf, trim first contig
		if ($firstQryConf < $secondQryConf) { 
			$firstContigSites = $firstContigSites + $labelsDistance;
			if ($ARGV[2] eq "+") {
				$firstContigEnd = $firstContigEnd - $posOffset;
				$firstContigSitesEndIdx = $firstContigSitesEndIdx - $sitesOffset;
				$forwardOffset2 = abs($ARGV[7]);
			}
			else {
				$firstContigEnd = $firstContigEnd - $posOffset;
				#$firstContigStart = $firstContigStart + $posOffset;
				$firstContigSitesEndIdx = $firstContigSitesEndIdx - $sitesOffset;
				#$reversePosOffset = abs($ARGV[5]);
				
			}
		}
		else {
			$secondContigSites = $secondContigSites + $labelsDistance;
			if ($ARGV[4] eq "+") {
				$secondContigStart = $secondContigStart + $posOffset;
				$secondContigSitesStartIdx = $secondContigSitesStartIdx + $sitesOffset;	
				$forwardOffset = abs($ARGV[7]);
				$forwardPosOffset = abs($ARGV[5]);

						
			}
			else {
				#$secondContigEnd = $secondContigEnd - $posOffset;
				$secondContigStart = $secondContigStart + $posOffset;
				$secondContigSitesStartIdx = $secondContigSitesStartIdx + $sitesOffset;
				$reverseOffset = abs($ARGV[7]);
			}
		}
		
		$positionOffset = ($firstContigEnd - $firstContigStart);
	}


	my $mergedLength = ($firstContigEnd - $firstContigStart) + ($secondContigEnd - $secondContigStart);
	if ($ARGV[5] > 1) {
		$mergedLength = $mergedLength + abs($ARGV[5]) + 1;
	}
	
	my $mergedSites = $firstContigSites + $secondContigSites;
		
	
			
	if ($missedLabelPadding != 0) {
		#$positionOffset = $firstContigEnd + ($ARGV[5] - $missedLabelPadding) + 1;
		$mergedSites += 1;
	}

	my $theRestCount=0; 
	my $lastSite=0;
	my $lastPos=0;

	
	
#IF orientation is +/+
if (($ARGV[2] eq '+' && $ARGV[4] eq '+')) {	

	#output first contig	
	for (my $i=$firstContigSitesStartIdx; $i < $firstContigSitesEndIdx; $i++ ) {
		my $hash = $firstContig[$i];	#reference to the ith value in anonymous hash (containing first contig cmap) 
		
		if ($hash->{'Position'} > $firstContigEnd) {
			$doOutput = 0;
			last;
		}
		
		my $theRest = $hash->{'TheRest'};
		my @s = split("\t",$theRest);
		$theRestCount = scalar(@s);		
		$lastSite = $hash->{'SiteID'};
		$lastPos = $hash->{'Position'};
		
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => "$hash->{'SiteID'}", # 1 (hash element, scalar)
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1 (hash element, scalar)
				"Position"  => "$hash->{'Position'}", # 20.0 (hash element, scalar)
				"StdDev" => "$hash->{'StdDev'}", # 81.9 (hash element, scalar)
				"Coverage" => "$hash->{'Coverage'}", # 14.0 (hash element, scalar)
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0 (hash element, scalar)
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond (hash element, array)
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);
		push @mergedContig, \%new_hash_ref; }
		
	my $offset = $firstContigSites;
	if ($missedLabelPadding != 0) {
		my ($missedSite_ref, $missedLabelPos) = getMissedLabelHash($theRestCount, $missedLabelPadding, $lastPos, $lastSite);
		my %missedSite = %$missedSite_ref;
		push @mergedContig, \%missedSite; 
		$offset = $offset + 1;
	}
			

	#output second contig
	#my $positionOffset = $firstContigEnd + $ARGV[5];
	#my $positionOffset = $firstContigEnd + $missedLabelPadding + 1;
	 
	for (my $i=$secondContigSitesStartIdx; $i < ($secondContigSitesEndIdx+1); $i++ ) {
		my $hash = $secondContig[$i];

		my $pos = $hash->{'Position'} + $positionOffset - $forwardPosOffset;
		if ($i == $secondContigSitesEndIdx) {
			$pos = $mergedLength;
		}

		if ($pos < $lastPos) {
			$doOutput = 0;
			last;
		}
		
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => $hash->{'SiteID'} + $offset - $forwardOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => $pos, # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref; }
}


#IF orientation is +/- 
elsif (($ARGV[2] eq '+' && $ARGV[4] eq '-')) {
	
	#output first contig 
	for (my $i=$firstContigSitesStartIdx; $i < $firstContigSitesEndIdx; $i++ ) {
		my $hash = $firstContig[$i];

		if ($hash->{'Position'} > $firstContigEnd) {
			$doOutput = 0;
			last;
		}
		
		my $theRest = $hash->{'TheRest'};
		my @s = split("\t",$theRest);
		$theRestCount = scalar(@s);		
		$lastSite = $hash->{'SiteID'};
		$lastPos = $hash->{'Position'};
		
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => "$hash->{'SiteID'}", # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => "$hash->{'Position'}", # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref; }
		
	my $offset = $firstContigSites; 
	if ($missedLabelPadding != 0) {
		my ($missedSite_ref, $missedLabelPos) = getMissedLabelHash($theRestCount, $missedLabelPadding, $lastPos, $lastSite);
		my %missedSite = %$missedSite_ref;
		push @mergedContig, \%missedSite; 
		$offset = $offset + 1;
	}
	
	#output REVERSED second contig
	#my $positionOffset = $firstContigEnd + $ARGV[5];
	#my $positionOffset = $firstContigEnd + $missedLabelPadding + 1;

	@secondContig = reverse(@secondContig);
	if ($doOutput==0) {
		@secondContig = reverse(@secondContig);
	}

	
	for (my $i=(1+$secondContigSitesStartIdx); $i < (1+$secondContigSitesEndIdx); $i++ ) {
		my $hash = $secondContig[$i];
		
		my $pos = ($secondContigEndOrig - $hash->{'Position'}) + $positionOffset - $secondContigStart;
		if ($pos < $lastPos) {
			$doOutput = 0;
			@secondContig = reverse(@secondContig);
			last;
		}
	
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($secondContigSitesOrig - $hash->{'SiteID'}) + $offset - $reverseOffset + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => $pos, # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
		my $hash = $secondContig[0];
		my %new_hash_ref = (
			"CMapId"  => "$mergedId", # 2020
			"ContigLength" => "$mergedLength", # 718132.6
			"NumSites"  => "$mergedSites", # 74
			"SiteID"  => $mergedSites + 1, # 1
			"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
			"Position"  => $mergedLength, # 20.0
			"StdDev" => "$hash->{'StdDev'}", # 81.9
			"Coverage" => "$hash->{'Coverage'}", # 14.0
			"Occurrence" => "$hash->{'Occurrence'}", # 14.0
			"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
			#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
			#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
			#"SNR" => "$hash->{'SNR'}" # 0
		);
		push @mergedContig, \%new_hash_ref;		
}

#IF orientation is -/+ 
if (($ARGV[2] eq '-' && $ARGV[4] eq '+')) {

	#output REVERSED first contig
	my @reverse_firstContig = reverse(@firstContig);
	#print scalar(@reverse_firstContig)."\n";
	for (my $i=(1+$firstContigSitesStartIdx); $i < (1+$firstContigSitesEndIdx); $i++ ) {
		my $hash = $reverse_firstContig[$i];
		
		if ( ($firstContigEndOrig - $hash->{'Position'}) > ($firstContigEnd)) {
			$doOutput = 0;
			@firstContig = reverse(@reverse_firstContig);
			last;
		}
		
		my $theRest = $hash->{'TheRest'};
		my @s = split("\t",$theRest);
		$theRestCount = scalar(@s);		
		$lastSite = ($firstContigSitesOrig - $hash->{'SiteID'}) + 1;
		$lastPos = ($firstContigEndOrig - $hash->{'Position'});

		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($firstContigSitesOrig - $hash->{'SiteID'}) + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($firstContigEndOrig - $hash->{'Position'}), # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
		
	my $offset = $firstContigSites;
	if ($missedLabelPadding != 0) {
		my ($missedSite_ref, $missedLabelPos) = getMissedLabelHash($theRestCount, $missedLabelPadding, $lastPos, $lastSite);
		my %missedSite = %$missedSite_ref;
		push @mergedContig, \%missedSite; 
		$offset = $offset + 1;
	} 
	
	#output second contig
	#my $positionOffset = $firstContigEnd + $missedLabelPadding + 1;
	
	for (my $i=($secondContigSitesStartIdx); $i < ($secondContigSitesEndIdx+1); $i++ ) {
		my $hash = $secondContig[$i];

		my $pos = $hash->{'Position'} + $positionOffset - $forwardPosOffset + $reversePosOffset;
		if ($i == $secondContigSitesEndIdx) {
			$pos = $mergedLength;
		}

		if ($pos < $lastPos) {
			@firstContig = reverse(@reverse_firstContig);
			$doOutput = 0;
			last;
		}

		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => $hash->{'SiteID'} + $offset - $forwardOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				#"Position"  => $hash->{'Position'} + $positionOffset, # 20.0
				"Position"  => $pos,
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref; }
}



#IF orientation is -/-
elsif (($ARGV[2] eq "-" && $ARGV[4] eq "-")) {
	
	#output REVERSED first contig
	my @firstContig_reverse = reverse(@firstContig);
	
	for (my $i=(1+$firstContigSitesStartIdx); $i < (1+$firstContigSitesEndIdx); $i++ ) {
		my $hash = $firstContig_reverse[$i];

		if ( ($firstContigEndOrig - $hash->{'Position'}) > ($firstContigEnd)) {
			@firstContig = reverse(@firstContig_reverse);
			$doOutput = 0;
			last;
		}
		
		my $theRest = $hash->{'TheRest'};
		my @s = split("\t",$theRest);
		$theRestCount = scalar(@s);		
		$lastSite = ($firstContigSitesOrig - $hash->{'SiteID'}) + 1;
		$lastPos = ($firstContigEndOrig - $hash->{'Position'});

		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($firstContigSitesOrig - $hash->{'SiteID'}) + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($firstContigEndOrig - $hash->{'Position'}), # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
	
	my $offset = $firstContigSites;
	if ($missedLabelPadding != 0) {
		my ($missedSite_ref, $missedLabelPos) = getMissedLabelHash($theRestCount, $missedLabelPadding, $lastPos, $lastSite);
		my %missedSite = %$missedSite_ref;
		push @mergedContig, \%missedSite; 
		$offset = $offset + 1;
	} 
	
	
	#output REVERSED second contig
	#my $positionOffset = $firstContigEnd + $missedLabelPadding;
		
	@secondContig = reverse(@secondContig);	
	if ($doOutput==0) {
		@secondContig = reverse(@secondContig);
	}

	
	for (my $i=(1+$secondContigSitesStartIdx); $i < (1+$secondContigSitesEndIdx); $i++ ) {
		my $hash = $secondContig[$i];
		my $pos = ($secondContigEndOrig - $hash->{'Position'}) + $positionOffset - $secondContigStart;
		
		if ($pos < $lastPos) {
			$doOutput = 0;
			@firstContig = reverse(@firstContig);
			@secondContig = reverse(@secondContig);
			last;
		}
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($secondContigSitesOrig - $hash->{'SiteID'}) + $offset - $reverseOffset + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => $pos, # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
		my $hash = $secondContig[0];
		my %new_hash_ref = (
			"CMapId"  => "$mergedId", # 2020
			"ContigLength" => "$mergedLength", # 718132.6
			"NumSites"  => "$mergedSites", # 74
			"SiteID"  => $mergedSites + 1, # 1
			"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
			"Position"  => $mergedLength, # 20.0
			"StdDev" => "$hash->{'StdDev'}", # 81.9
			"Coverage" => "$hash->{'Coverage'}", # 14.0
			"Occurrence" => "$hash->{'Occurrence'}", # 14.0
			"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
			#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
			#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
			#"SNR" => "$hash->{'SNR'}" # 0
		);
		push @mergedContig, \%new_hash_ref;		
}

	


#print out new CMAP with merged contigs
#sort mergedContig by SiteID
#my @sorted_mergedContig =  sort { $a->{SiteID} <=> $b->{SiteID} } @mergedContig;

if ($doOutput == 1) {
	foreach (@mergedContig) { #change to @mergedContig
		push @cmap_out, $_; 
	print "TRUE\n";
	}
}
else {
	#if ($firstQryConf >= $secondQryConf) {
		push @cmap_out, @firstContig;
	#}
	#else {
		push @cmap_out, @secondContig;
	#}
	print EXCLUDE "$ARGV[1]\t$ARGV[3]\n"
}
	
foreach my $hash_ref (@cmap_out) {
	# mod: eva chan, 9 july 2015
	# my $str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}\t$hash_ref->{'GmeanSNR'}\t$hash_ref->{'lnSNRsd'}\t$hash_ref->{'SNR'}";
	my $str = "";
	if(defined($hash_ref->{'TheRest'}) && $hash_ref->{'TheRest'} ne "") { 
		$str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}\t$hash_ref->{'TheRest'}";
	} 
	else { 
		$str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}";
	}
	print OUT "$str\n";	
}


#rename output file
rename $outName, $outName_orig;
	
		

#print "END OF OUTPUT mergeContigs.pl\n";

close FILE;
close OUT;
close EXCLUDE;


sub getMissedLabelHash {
	my ($theRestCount, $missedLabelPadding, $lastPos, $lastSite) = @_;
	my @r;
	for (my $i=0; $i<$theRestCount; $i++) {
		push @r, "0";
	}
	my $theRestOut = join("\t", @r);
	$lastSite += 1;
	my $missedLabelPos = $missedLabelPadding + $lastPos + 1;
	my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => "$lastSite", # 1 (hash element, scalar)
				"LabelChannel"  => "1", # 1 (hash element, scalar)
				"Position"  => "$missedLabelPos", # 20.0 (hash element, scalar)
				"StdDev" => "1", # 81.9 (hash element, scalar)
				"Coverage" => "1", # 14.0 (hash element, scalar)
				"Occurrence" => "1", # 14.0 (hash element, scalar)
				"TheRest" => $theRestOut #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond (hash element, array)
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);
	return \%new_hash_ref, $missedLabelPos; 	
}
