#!/usr/bin/perl

# a standalone script to validate stitchPositions based on orthogonal hybrid scaffold XMAP aligned to same reference

# Usage: perl validate_stitchPositionsBED.pl -x <hybridScaffold_vs_hg19> -b <stitchPositions_scored BED file> [-s NGS BAM file vs hg19 requires hg19 key file and fasta] [-a hg19 fasta] [-k reference key file from fa2cmap] [-f <minimum required number bases flanking stitchPosition in XMAP alignment default=1000>]

# Details: 
# * Assumption: that contigs on XMAP is being read from left to right and is sorted by RefStartPos
# * Assumption: stitchPositions_scored BED file is from fragileSiteRepair pipeline
# * Designed to work with XMAP 0.1

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Find; 
use File::Basename;
use File::Copy;
use Getopt::Long qw(:config bundling); 
use Bio::DB::Sam;

my %inputs = (); 
GetOptions( \%inputs, 'x|xmap=s', 'b|bed=s', 'f|flank=i', 's|bam=s', 'a|fasta=s', 'k|key=s'); 

if( !exists $inputs{x} | !exists $inputs{b} ) {
	print "Usage: perl validate_stitchPositionsBED.pl -x <hybridScaffold_vs_hg19> -b <stitchPositions_scored BED file> [-s NGS BAM file vs hg19 requires hg19 key file and fasta] [-a hg19 fasta] [-k reference key file from fa2cmap] [-f <minimum required number bp flanking stitchPosition in XMAP alignment default=1000>]\n"; 
	exit 0; 
}
else {
	foreach my $key ("x","b","s","k", "a") {
		if (exists $inputs{$key} and $inputs{$key} =~ /[a-z|0-9]/i) {
			$inputs{$key} = abs_path($inputs{$key});
		}
	}
}
if ( !exists $inputs{f} ) {
	$inputs{f} = 1000;
}


my $bedBase = basename($inputs{b},".bed");
my $bedValidated = $bedBase."_validated.bed";
my $bedUnvalidated = $bedBase."_unvalidated.bed";
open VALIDATED, ">$bedValidated" or die "ERROR: Could not open $bedValidated: $!\n";
open UNVALIDATED, ">$bedUnvalidated" or die "ERROR: Could not open $bedUnvalidated: $!\n";


open XMAP, $inputs{x} or die "ERROR: $!\n";
my @xmap;
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

open BED, $inputs{b} or die "ERROR: $!\n";
my @bed;
## << stitchPositions_scored bed >> 
while (my $line = <BED>) {
	chomp($line);
	#if header then skip
	if ($line =~ /^#/) {
		next; }
	else {
		my @s = split(/\t/,$line);
		#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
		my %bed_line = (
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
			"line" => $line
		);
		push @bed, \%bed_line;
	}
}

# if input BAM and fasta and key files specified, load
my $sam; my $bai;
my @keys; my $fai;
my $hasBam=0;
if ( defined $inputs{a} && exists $inputs{a} && defined $inputs{s} && exists $inputs{s} && defined $inputs{k} && exists $inputs{k}) {
	# force index of fasta file
	$bai = $inputs{s} . ".bai";
	$fai = $inputs{a} . ".fai"; 
	if (-e $bai && -e $fai) {
		$sam = Bio::DB::Sam->new( -fasta=>"$inputs{a}", -fai=>"$fai", -bam=>"$inputs{s}", -bai=>"$bai" );
	}
	else {
		my $status_code = Bio::DB::Bam->index_build($inputs{s});
		my $index = Bio::DB::Sam::Fai->load($inputs{a});
		$sam = Bio::DB::Sam->new( -fasta=>"$inputs{a}", -bam=>"$inputs{s}" );
	}
	$hasBam=1;
	#print "\n";
	#print "Input reference FASTA: $inputs{fasta}\n";
	#print "Input reference FASTA key: $inputs{key}\n";
	#print "Input NGS alignment BAM: $inputs{bam}\n";
	#load key file
	
	open KEY, "$inputs{k}" or die "ERROR: Could not open $inputs{k}: $!\n";
	while (my $line = <KEY>) {
		chomp($line);
		#if header then skip
		if ($line =~ /^#/) {
			next; }
		elsif ($line =~ /^CompntId/) {
			next; }
		else {
			my @s = split("\t",$line);
			#CompntId	CompntName	CompntLength
			my %keyLine = (
				"CompntId"  => "$s[0]", 
				"CompntName" => "$s[1]", 
				"CompntLength"  => "$s[2]"
			);
			push @keys, \%keyLine; 
		}	
	}
	#print "Read in ".scalar(@keys)." keys from $inputs{key}\n";
	#print "\n";		
}


# DO WORK
my $validatedCount=0;
my $unvalidatedCount=0;
foreach my $stitchPosition_ref (@bed) {
	my %stitchPosition = %$stitchPosition_ref;
	# increase window to including $inputs{f} flanking bases
	$stitchPosition{Start} = $stitchPosition{Start} - $inputs{f};
	if ($stitchPosition{Start} < 1) { $stitchPosition{Start} = 1; } 
	$stitchPosition{End} = $stitchPosition{End} + $inputs{f};
	
	my $validated = 0;
	
	# IF BAM file provided
	if ($hasBam==1 && $validated==0) {
		my $seqName;
		my $seqLen=0;
		foreach my $keyEntry_ref (@keys) {
			my %keyEntry = %$keyEntry_ref;
			#print "bedCmapId: $bedCmapId CompntId: $keyEntry{CompntId} CompntName: $keyEntry{CompntName}\n";
			if ($keyEntry{CompntId} eq $stitchPosition{CMapId}) {
				$seqName = $keyEntry{CompntName};
				$seqLen = $keyEntry{CompntLength};
				#print "Looking at $seqName with length $seqLen\n";
			}
		}
		my $bamStart = $stitchPosition{Start} - 1000000; if ($bamStart < 1) { $bamStart = 1; }
		my $bamEnd = $stitchPosition{End} + 1000000; if ($bamEnd > $seqLen) { $bamEnd = $seqLen; }
		
		my @targets = $sam->seq_ids;
		my $seqId="";
		foreach my $i (@targets) {
			#print "\tNGS SeqId: $i\n";
			if ( ($i =~ /$seqName/i) || ($seqName =~ /$i/i) ) {
				$seqId = $i;
				#print "\tFound NGS SeqId: $i\n";
			}
		}
		
		my @alignments = $sam->get_features_by_location(-seq_id => "$seqId", -start => "$bamStart", -end => "$bamEnd");
		
		foreach my $a (@alignments) {
			#my $seqid  = $a->seq_id;
			my $start  = $a->start;
			my $end    = $a->end;
			#my $strand = $a->strand;
			#my $ref_dna= $a->dna;

			#my $query_start  = $a->query->start;
			#my $query_end    = $a->query->end;
			#my $query_strand = $a->query->strand;
			#my $query_dna    = $a->query->dna;
			my $read_name = $a->query->name;
		   
			#my $cigar     = $a->cigar_str;
			#my @scores    = $a->qscore;     # per-base quality scores
			#my $score     = $a->qstring;    # TAM-style quality string
			my $match_qual= $a->qual;       # quality of the match

			#my $paired = $a->get_tag_values('PAIRED');
			
			#print "\t\t$start $end $match_qual\n";
			
			if (($start < $stitchPosition{Start}) && ($end > $stitchPosition{End}) && ($match_qual >= 30)) {
				print VALIDATED "$stitchPosition{line}\tBAM\n";
				$validatedCount++;
				$validated=1;
				last;
			}
		}
	}	
	
	#IF XMAP provided
	if ($validated==0) {
		foreach my $xmapLine_ref (@xmap) {
			my %xmapLine = %$xmapLine_ref;
			if ($xmapLine{RefContigID} eq $stitchPosition{CMapId}) {
				if ($stitchPosition{End} > $xmapLine{RefLen}) {
					$stitchPosition{End} = $xmapLine{RefLen};
				}
			
				if ( ($xmapLine{RefStartPos} <= $stitchPosition{Start}) &&  ($xmapLine{RefEndPos} >= $stitchPosition{End}) ) {
					print VALIDATED "$stitchPosition{line}\tXMAP\n";
					$validatedCount++;
					$validated=1;
					last;
				}
			}
		}
	}
	
	
	
	if ($validated eq 0) {
		print UNVALIDATED "$stitchPosition{line}\n";		
		$unvalidatedCount++;
	}			
}

print "FlankingBases: $inputs{f} Validated: $validatedCount Unvalidated: $unvalidatedCount\n";
