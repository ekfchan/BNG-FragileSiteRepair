#!/usr/bin/perl

#Usage: perl extractContigs.pl <input XMAP> <_q.cmap> <_r.cmap> <first_contig_ID> <second_contig_ID>
#assumption is that contigs on XMAP is being read from left to right and is sorted by RefStartPos
#if .fsites file is provided, 

#Designed to work with XMAP 0.1

use strict;
use warnings;
use IPC::System::Simple qw(system capture);
use Cwd;
use File::Basename; 	#mod: eva chan, 8 july 2015, to allow trace path of scripts

my $scriptspath = Cwd::abs_path(dirname($0));	#mod: eva chan, 8 july 2015

open FILE, "$ARGV[0]" or die $!;

#open output BED file to print stitch locations
my $out;
if (-e $ARGV[0]) { 
	$out = Cwd::abs_path($ARGV[0]); }
else {
	die "ERROR: Input $ARGV[0] does not exist! $!\n";
}
$out =~ s/.xmap//i;
$out =~ s/_temp//i;
$out = $out."_fragileSiteRepaired_stitchPositions.bed";
if (!-e $out) {
	open OUT, ">$out" or die "ERROR: Cannot open $out for writing! $!\n";
	print OUT "#CMapId\tStart\tEnd\n"; }
else {
	open OUT, ">>$out" or die "ERROR: Cannot open $out for writing! $!\n";
}

my %firstContigAlignment;
my %secondContigAlignment;

#read input XMAP
while (my $line = <FILE>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		if ($s[1] eq $ARGV[3]) {
			#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum
			#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string 
			# 1	14839	17	344.2	999946.6	42640026.3	43628302.0	+	138.45	2M1D4M1D8M1D10M1D7M1D7M1D2M1D6M1D2M1D2M1I10M2D3M1D13M2D6M1D10M1D10M1D6M5I4M
			%firstContigAlignment = (
				"XmapEntryID"  => "$s[0]", 
				"QryContigID" => "$s[1]", 
				"RefContigID"  => "$s[2]",
				"QryStartPos"  => "$s[3]", 
				"QryEndPos"  => "$s[4]", 
				"RefStartPos"  => "$s[5]", 
				"RefEndPos" => "$s[6]", 
				"Orientation" => "$s[7]", 
				"Confidence" => "$s[8]", 
				"HitEnum" => "$s[9]"
			); }
		elsif ($s[1] eq $ARGV[4]) {
			%secondContigAlignment = (
				"XmapEntryID"  => "$s[0]", 
				"QryContigID" => "$s[1]", 
				"RefContigID"  => "$s[2]",
				"QryStartPos"  => "$s[3]", 
				"QryEndPos"  => "$s[4]", 
				"RefStartPos"  => "$s[5]", 
				"RefEndPos" => "$s[6]", 
				"Orientation" => "$s[7]", 
				"Confidence" => "$s[8]", 
				"HitEnum" => "$s[9]"
			); 
		}
	}
}

#calculate padding between two contigs
my $firstContigEnd = $firstContigAlignment{'RefEndPos'};
my $secondContigStart = $secondContigAlignment{'RefStartPos'};
my $padding = ($secondContigStart - $firstContigEnd) + 1;

#print out stitch locations to BED file
#print "\tPrinting out stitch locations to $out\n";
#print "\t$firstContigAlignment{'RefContigID'}\t$firstContigAlignment{'RefEndPos'}\t$secondContigAlignment{'RefStartPos'}\n";
print OUT "$firstContigAlignment{'RefContigID'}\t".int(($firstContigAlignment{'RefEndPos'}-1))."\t".int(($secondContigAlignment{'RefStartPos'}+1))."\n";

#run script to merge two contigs
# mod: eva chan, 8 july 2015, assuming path of extractContigs.pl to be same as gapFill.pl (rather than in ../)
#my $mergeScript = Cwd::abs_path("../mergeContigs.pl");
my $mergeScript = $scriptspath."/mergeContigs.pl";

my @ARGS = ($ARGV[1], $ARGV[3], $firstContigAlignment{'Orientation'}, $ARGV[4], $secondContigAlignment{'Orientation'}, $padding);
print "Running command: ".$^X." $mergeScript ". join(" ",@ARGS)."\n";
system($^X, "$mergeScript", @ARGS);



print "END OF OUTPUT extractContigs.pl\n";


close OUT;
close FILE;