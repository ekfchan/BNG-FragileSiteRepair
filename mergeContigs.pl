#!/usr/bin/perl

#Usage: perl mergeContigs.pl <input CMAP> <first_contig_ID> <+/- orientation> <second_contig_ID> <+/- orientation> <N bases to add in-between>
#Designed to be used with non-overlapping contigs

# 0=====>N0=====>N == + +
# N<=====0N<=====0 == - -
# 0=====>NN<=====0 == + -
# N<=====00=====>N == - +

#Designed to work with CMAP 0.1

use strict;
use warnings;
use Data::Dumper;

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

my $idOffset = 100000;
my @cmap_out;		#references of hash table of cmap to print out
my @firstContig;	#reference of hash table of cmaps whose contig is <first_contig_ID>
my @secondContig;	#reference of hash table of cmaps whose contig is <second_contig_ID>
my @mergedContig;

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

	
	my $firstSites = scalar(@firstContig) - 1;	#number of first contigs (elements in @firstContig)
	#print "Sites in first contig: $firstSites\n";
	#my $firstContigStart = $firstContig[0]->{'Position'};
	my $firstContigEnd = $firstContig[$firstSites]->{'Position'};	#assume the contigs are sorted
	my $firstContigSites = $firstContig[$firstSites]->{'NumSites'}; 
	
	#my $secondContigMolStart = ($secondContig[0]->{'Position'}) - 20;
	my $secondSites = scalar(@secondContig) - 1;	#number of second contigs (elements in @secondContig)
	#my $secondContigStart = $secondContig[0]->{'Position'};
	my $secondContigEnd = $secondContig[$secondSites]->{'Position'};
	my $secondContigSites = $secondContig[$secondSites]->{'NumSites'};
	
	my $mergedLength = $firstContigEnd + $secondContigEnd + $ARGV[5];
	my $mergedSites = $firstContigSites + $secondContigSites;
	
#IF orientation is +/+
if (($ARGV[2] eq '+' && $ARGV[4] eq '+')) {	
	
	#output first contig 
	for (my $i=0; $i < $firstSites; $i++ ) {
		my $hash = $firstContig[$i];	#reference to the ith value in anonymous hash (containing first contig cmap) 
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
	
	#output second contig
	my $positionOffset = $firstContigEnd + $ARGV[5];
	my $sitesOffset = $firstSites; 
	for (my $i=0; $i < ($secondSites+1); $i++ ) {
		my $hash = $secondContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => $hash->{'SiteID'} + $sitesOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => $hash->{'Position'} + $positionOffset, # 20.0
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
	for (my $i=0; $i < $firstSites; $i++ ) {
		my $hash = $firstContig[$i];
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
	
	#output REVERSED second contig
	my $positionOffset = $firstContigEnd + $ARGV[5];
	my $sitesOffset = $firstSites +1; 
	@secondContig = reverse(@secondContig);
	
	for (my $i=1; $i < ($secondSites+1); $i++ ) {
		my $hash = $secondContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($secondSites - $hash->{'SiteID'}) + $sitesOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($secondContigEnd - $hash->{'Position'}) + $positionOffset, # 20.0
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
			"Position"  => $secondContigEnd + $positionOffset, # 20.0
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
	for (my $i=1; $i < ($firstSites+1); $i++ ) {
		my $hash = $reverse_firstContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($firstSites - $hash->{'SiteID'}) + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($firstContigEnd - $hash->{'Position'}), # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
	
	#output second contig
	my $positionOffset = $firstContigEnd + $ARGV[5];
	my $sitesOffset = $firstSites; 
	for (my $i=0; $i < ($secondSites+1); $i++ ) {
		my $hash = $secondContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => $hash->{'SiteID'} + $sitesOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => $hash->{'Position'} + $positionOffset, # 20.0
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
	@firstContig = reverse(@firstContig);
	for (my $i=1; $i < ($firstSites+1); $i++ ) {
		my $hash = $firstContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($firstSites - $hash->{'SiteID'}) + 1, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($firstContigEnd - $hash->{'Position'}), # 20.0
				"StdDev" => "$hash->{'StdDev'}", # 81.9
				"Coverage" => "$hash->{'Coverage'}", # 14.0
				"Occurrence" => "$hash->{'Occurrence'}", # 14.0
				"TheRest" => $hash->{'TheRest'} #mod: eva chan, 8 july 2015, some _q.cmap have no columns beyond 
				#"GmeanSNR" => "$hash->{'GmeanSNR'}", # 12.4650
				#"lnSNRsd" => "$hash->{'lnSNRsd'}", # 0.4814				
				#"SNR" => "$hash->{'SNR'}" # 0
			);		
		push @mergedContig, \%new_hash_ref;	}
	
	#output REVERSED second contig
	my $positionOffset = $firstContigEnd + $ARGV[5];
	my $sitesOffset = $firstSites +1; 
	@secondContig = reverse(@secondContig);
	
	for (my $i=1; $i < ($secondSites+1); $i++ ) {
		my $hash = $secondContig[$i];
		my %new_hash_ref = (
				"CMapId"  => "$mergedId", # 2020
				"ContigLength" => "$mergedLength", # 718132.6
				"NumSites"  => "$mergedSites", # 74
				"SiteID"  => ($secondSites - $hash->{'SiteID'}) + $sitesOffset, # 1
				"LabelChannel"  => "$hash->{'LabelChannel'}", # 1
				"Position"  => ($secondContigEnd - $hash->{'Position'}) + $positionOffset, # 20.0
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
			"Position"  => $secondContigEnd + $positionOffset, # 20.0
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

foreach (@mergedContig) { #change to @mergedContig
	push @cmap_out, $_; 
	}
	
foreach my $hash_ref (@cmap_out) {
	# mod: eva chan, 9 july 2015
	# my $str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}\t$hash_ref->{'GmeanSNR'}\t$hash_ref->{'lnSNRsd'}\t$hash_ref->{'SNR'}";
	my $str = "";
	if(defined($hash_ref->{'TheRest'}) && $hash_ref->{'TheRest'} ne "") { 
		$str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}\t$hash_ref->{'TheRest'}";
	} else { 
		$str = "$hash_ref->{'CMapId'}\t$hash_ref->{'ContigLength'}\t$hash_ref->{'NumSites'}\t$hash_ref->{'SiteID'}\t$hash_ref->{'LabelChannel'}\t$hash_ref->{'Position'}\t$hash_ref->{'StdDev'}\t$hash_ref->{'Coverage'}\t$hash_ref->{'Occurrence'}";
	}
	print OUT "$str\n";	
}


#rename output file
rename $outName, $outName_orig;
	
		

print "END OF OUTPUT mergeContigs.pl\n";

close FILE;
close OUT;
