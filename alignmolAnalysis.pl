#!/usr/bin/perl

# a standalone script to interogate single molecule alignments vs genome maps (alignmol)

# Usage: perl alignmolAnalysis.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -t <coverage threshold default=10>

# Details: 
# * Assumption: that contigs on XMAP is being read from left to right and is sorted by RefStartPos
# * Assumption: input _r.cmap has only one contig
# * Designed to work with XMAP 0.1

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Find; 
use File::Basename;
use File::Copy;
use Getopt::Long qw(:config bundling); 

#print "\n";
#print qx/ps -o args $$/;
#print "\n";

## << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'x|xmap=s', 'q|qcmap=s', 'r|rcmap=s', 't|threshold=i'); 

if( !exists $inputs{x} | !exists $inputs{q} | !exists $inputs{r} ) {
	print "Usage: perl alignmolAnalysis.pl -x <input.xmap> -q <input_q.cmap> -r <input_r.cmap> -t <coverage threshold default=10>\n"; 
	exit 0; 
}
else {
	foreach my $key ("x","q","r") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z|0-9]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
}

if( !exists $inputs{t} || $inputs{t} < 0 ) {
	$inputs{t} = 10;
}


open XMAP, $inputs{x} or die "ERROR: $!\n";
open QCMAP, $inputs{q} or die "ERROR: $!\n";
open RCMAP, $inputs{r} or die "ERROR: $!\n";
my @xmap;
my @qcmap;
my @rcmap;

## << read input XMAP >> 
while (my $line = <XMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
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
#print "\n";
close XMAP;

## << read input query CMAP >> 
while (my $line = <QCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		for (my $i=0; $i<scalar(@s); $i++) {
			$s[$i] =~ s/^\s+|\s+$//g; }
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	GmeanSNR	lnSNRsd	SNR
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	12.4650	0.4814	0
		my %cmap_line = (
			"CMapId"  => "$s[0]",
			"ContigLength" => "$s[1]",
			"NumSites"  => "$s[2]",
			"SiteID"  => "$s[3]",
			"LabelChannel"  => "$s[4]",
			"Position"  => "$s[5]",
			"StdDev" => "$s[6]",
			"Coverage" => "$s[7]",
			"Occurrence" => "$s[8]"
			);
		push @qcmap, \%cmap_line;	
	}
}
#print "\n";
close QCMAP;

## << read input reference CMAP >> 
while (my $line = <RCMAP>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		#load first contig into hash
		##h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence	
		#  2020	 718132.6	74	1	1	20.0	81.9	14.0	14.0	
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
		push @rcmap, \%cmap_line;	}
}
#print "\n";
close RCMAP;

# parse _r.cmap to get first and last labelID
my $rcmapEntry_ref = $rcmap[0]; 
my %rcmapEntry = %$rcmapEntry_ref;
my $firstRefId = 1;
my $firstRefIdOcc = $rcmapEntry{'Occurrence'};
my $lastRefId = $rcmapEntry{'NumSites'};
$rcmapEntry_ref = $rcmap[-2]; 
%rcmapEntry = %$rcmapEntry_ref;
my $lastRefIdOcc = $rcmapEntry{'Occurrence'};
my $cmapId = $rcmapEntry{'CMapId'};
$cmapId =~ s/^\s+|\s+$//g;

#print "RefMapId: $cmapId RefFirstLabel: $firstRefId RefFirstOcc: $firstRefIdOcc RefLastLabel: $lastRefId RefLastOcc: $lastRefIdOcc\n";

my $mapStartEnds=0;
my $mapEndEnds=0;

foreach my $xmapEntry_ref (@xmap) {
	my %xmapEntry = %$xmapEntry_ref;
	my %matchPairs;
	if ($xmapEntry{'RefContigID'} eq $cmapId) {
		my $line = $xmapEntry{'Alignment'};
		my @s = split(/(?<=\))/, $line);
		my $refLabels = "";
		my $qryLabels = "";			
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
		if ($refLabels =~ /\b$firstRefId\b/) {
			#print "firstRefLabel xmapId: $xmapEntry{'XmapEntryID'}\n";
			#print "\tRefLabels: $refLabels\n";
			my $queryLabelId = $matchPairs{$firstRefId};
			my $queryId = $xmapEntry{'QryContigID'};
			#print "FIRST RefLabel: $firstRefId QueryLabel: $queryLabelId QueryId: $queryId\n";
			foreach my $qcmapEntry_ref (@qcmap) {
				my %qcmapEntry = %$qcmapEntry_ref;
				if ($qcmapEntry{'CMapId'} eq $queryId) {		
					my $lastSite = $qcmapEntry{'NumSites'};
					#print "FIRST RefLabel: $firstRefId QueryLabel: $queryLabelId QueryId: $queryId NumSites: $lastSite\n";
					if ($queryLabelId==1) {
						$mapStartEnds++; 
					}
					elsif ($queryLabelId==$lastSite) {
						$mapStartEnds++;
					}
				last;
				}
			}
		}
		if ($refLabels =~ /\b$lastRefId\b/) {
			#print "lastRefLabel xmapId: $xmapEntry{'XmapEntryID'}\n";
			#print "\tRefLabels: $refLabels\n";
			my $queryLabelId = $matchPairs{$lastRefId};
			my $queryId = $xmapEntry{'QryContigID'};
			#print "LAST RefLabel: $lastRefId QueryLabel: $queryLabelId QueryId: $queryId\n";
			foreach my $qcmapEntry_ref (@qcmap) {
				my %qcmapEntry = %$qcmapEntry_ref;
				if ($qcmapEntry{'CMapId'} eq $queryId) {		
					my $lastSite = $qcmapEntry{'NumSites'};
					#print "LAST RefLabel: $lastRefId QueryLabel: $queryLabelId QueryId: $queryId NumSites: $lastSite\n";
					if ($queryLabelId==1) {
						$mapEndEnds++; 
					}
					elsif ($queryLabelId==$lastSite) {
						$mapEndEnds++;
					}
				last;
				}
			}
		}
		else {
			next;
		}
	}	
}

#print "RefMapId: $cmapId MapStartEnds: $mapStartEnds MapEndEnds: $mapEndEnds\n";
#print "RefMapId: $cmapId MapStart fragility: $mapStartEnds/$firstRefIdOcc MapEnd fragility: $mapEndEnds/$lastRefIdOcc\n";
my $startRatio = 0;
if ($firstRefIdOcc >= $inputs{t}) {
	$startRatio = $mapStartEnds/$firstRefIdOcc;
}
my $endRatio = 0;
if ($lastRefIdOcc >= $inputs{t}) {
	$endRatio = $mapEndEnds/$lastRefIdOcc;
}
print "$cmapId\t$startRatio\t$endRatio\n";
