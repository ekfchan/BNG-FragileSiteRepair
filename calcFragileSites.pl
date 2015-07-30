#!/usr/bin/perl

# Usage: perl calcFragileSites.pl <input FASTA> [output.bed] [sequence buffer in bp]
# Example: perl calcFragileSites.pl hg19_chromosome.fa hg19_chromosome_bbspqi.bed 500

use strict;
use warnings;
use Cwd qw(realpath);
use File::Basename;
#use File::Copy;
use File::Path qw(make_path);
use Bio::Seq;
use Bio::SeqIO;
#use List::MoreUtils qw(uniq);
use Data::Dumper;

# usage statement
if( scalar(@ARGV)<1 or scalar(@ARGV)>3 ) {
	print "usage: $0 <input.fa> [output.bed] [sequence buffer in bp]\n"; 
	exit 0; 
}

# distance in basepairs between nicks to consider a potential fragile site
my $bp_t1 = 500; #maximum distance between nicks on opposite strands (moving towards each other) to count as TypeI
my $bp_t2 = 500; #maximum distance between nicks on opposite strands (moving away from each other) to count as TypeII
my $bp_t3 = 500; #maximum distance between nicks on same strands (moving towards each other) to count as TypeIII
my $bp_t4 = 500; #maximum distance between nicks on same strands (moving away from each other) to count as TypeIV

my $enzyme1 = "GCTCTTC"; #BspQI... motif is actually "GCTCTTCN" with nicking occuring after the "N"
#my $enzyme2 = "CACGAG"; #BssSI

my $outputfile; 
my $basename;
my $fasta;
my $fsitesCountTotal=0;

print "\n";

# check fasta input variable
if (defined $ARGV[0] && -e $ARGV[0]) {
	$fasta = realpath($ARGV[0]);
	#print "Calculating fragiles sites for FASTA file: $fasta\n";
	$basename = basename($fasta);
	$basename =~ s/\.[^.]+$//;
	#print "\tFile filename = $basename\n";
	}
else {
	die "ERROR: Input FASTA file $ARGV[0] not provided or does not exist! $!\n";
} 
print "Input FASTA: $fasta\n";	

# check optional output variable
if( defined $ARGV[1] and !-d $ARGV[1]) { 
	$outputfile = realpath($ARGV[1]); 
} 
elsif ( defined $ARGV[1] and -d $ARGV[1]) {
	#$outputfile = $ARGV[1]."/".$basename."_".uc($enzyme1)."_".uc($enzyme2).".fsites.withSeqs.bed";
	$outputfile = $ARGV[1]."/".$basename."_".uc($enzyme1)."_fsites.withSeqs.bed";
}
else {
	#$outputfile = realpath("./$basename.".uc($enzyme1)."_".uc($enzyme2).".fsites.withSeqs.bed");  #if output file not given use $basename as prefix 
	$outputfile = realpath("./$basename.".uc($enzyme1)."_fsites.withSeqs.bed");
}
if( -e $outputfile ) { 
	die "ERROR: Output file ", $outputfile, " already exists. Please specify alternate filename. $!\n"; 
}

#check input sequence buffer variable
my $buffer=50;	
if (defined $ARGV[2] && int($ARGV[2])>0) {
		$buffer = int($ARGV[2]);
}

print "Output BED: $outputfile\n\n";

print "\tTypeI fragile site threashold (bp): $bp_t1\n";
print "\tTypeII fragile site threashold (bp): $bp_t2\n";
print "\tTypeIII fragile site threashold (bp): $bp_t3\n";
print "\tTypeIV fragile site threashold (bp): $bp_t4\n\n";
print "\tExtracting sequences +/- $buffer"."bp from fragile site start/end positions\n";

print "\n";
## Read input fasta and find nicks and fragiles sites
#my $seqin = Bio::SeqIO->new( -format => 'Fasta' , -file => "$fasta");
my $seqin = Bio::SeqIO->new( -file => "$fasta");
my $count=1;
my @fsitesBed;
#my $fsitesHeaderBed = "#Nickase 1: $enzyme1\n#Nickase 2: $enzyme2\n#CMapId\tStart\tEnd\tType\tSequence";
my $fsitesHeaderBed = "#Nickase 1: $enzyme1\n#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence";
while((my $seqobj = $seqin->next_seq())) {   #for each sequence in FASTA (e.g. each contig)
	my $seq = $seqobj->seq();
	my $id  = $seqobj->display_id();
	my $length = $seqobj->length();

	# find nick sites
	print "\tHeader: $id Length: $length\n";
	#my %result = find_nick_sites($seq,$enzyme1,$enzyme2);
	my %result = find_nick_sites($seq,$enzyme1);	
	
	# sort nick sites by genomic position, creating an ordered array of hash references
	my @nick_sites;
	foreach my $key (sort{$a <=> $b} keys %result){
		my %hash = ($key, $result{$key});
		push @nick_sites, \%hash;
	}
	print "\t\tTotal nick sites found: ", scalar(@nick_sites), "\n"; 

	#calc fragile sites for each seq contig
	my (@fsiteBed) = calcFragilesSitesBED(@nick_sites);
	
	#extract sequence +/- $ARGV[2] or 500bp
	my @fsiteBedSeq;
	foreach my $line (@fsiteBed) {
		#print "\t\tBED Line: $line\n";
		#my $lineBed = "$count\t$pos1\t$pos2\t$curtype";
		my @s = split("\t",$line);
		my $start = $s[1]-$buffer;
		if ($start <1) {
			$start = 1;
		}
		my $end = $s[2]+$buffer;
		if ($end >$length) {
			$end = $length-1;
		}
		my $seq = $seqobj->subseq($start,$end);
		
		#set fsite type RGBA codes
		my $itemRgb;
		if ($s[3] =~ m/TypeI/i) {
			$itemRgb = "255,0,0,255"; }
		elsif ($s[3] =~ m/TypeII/i) {
			$itemRgb = "255,255,0,255"; }
		elsif ($s[3] =~ m/TypeIII/i) {
			$itemRgb = "255,0,255,255"; }
		elsif ($s[3] =~ m/TypeIV/i) {
			$itemRgb = "0,255,255,255"; }
		
		#build and push BED line
		#CMapId\tStart\tEnd\tType\tScore\tStrand\thickStart\thickEnd\titemRgb\tSequence";
		my $lineOut = "$s[0]\t$s[1]\t$s[2]\t$s[3]\t0\t0\t$s[1]\t$s[2]\t$itemRgb\t$seq";
		push @fsiteBedSeq, $lineOut;
	}
	
	#push @fsitesBed, @fsiteBed;	
	push @fsitesBed, @fsiteBedSeq;
	$count++;
}

## output merged fsites file
unshift @fsitesBed, $fsitesHeaderBed;
open (OUT2, ">$outputfile") or die "ERROR: Could not open $outputfile $!\n";
print OUT2 join("\n",@fsitesBed);
print "\n";

print "Total potential fragile sites found in $fasta: $fsitesCountTotal\n";
#print "Merged fragile sites file: $outputfile\n";
#print "Finished calculating fragile sites for FASTA file: $fasta\n";
print "\n";
	

	
	
	
	
	
	
	
	
	
	
	
	
sub find_nick_sites{
	# Find position of the "N" base in the enzyme motif "GCTCTTCN"
	#my ($seq, $enzyme1, $enzyme2) = @_;
	my ($seq, $enzyme1) = @_;
	my %result;
	my $slength = length($seq);
	my $elength1 = length($enzyme1);
	#my $elength2 = length($enzyme2);
	
	# Find the first enzyme in the forward strand, starting from the first nucleotide!!!
	my $current_loc = index($seq, $enzyme1, 0);
	while ($current_loc != -1){
		if($current_loc + $elength1 < $slength){
			# $result{$current_loc + $elength + 1} = 1;	#records position of base after "N" (not part of motif)
			#$result{$current_loc + $elength} = 1;	#records position at base "N" GCTCTTCN
			$result{$current_loc} = 1;	#records position at base "G" GCTCTTCN
		}
		$current_loc = index($seq, $enzyme1, $current_loc + 1);
	}
	
	# Find the second enzyme in the forward strand, starting from the first nucleotide!!!
	# $current_loc = index($seq, $enzyme2, 0);
	# while ($current_loc != -1){
		# if($current_loc + $elength2 < $slength){
			# # $result{$current_loc + $elength + 1} = 1;	#records position of base after "N" (not part of motif)
			# #$result{$current_loc + $elength} = 1;	#records position at base "N" GCTCTTCN
			# $result{$current_loc} = 1;	#records position at base "G" GCTCTTCN
		# }
		# $current_loc = index($seq, $enzyme2, $current_loc + 1);
	# }
	
	# Find the reverse(first enzyme) in the forward strand, starting from the first nucleotide!!!
	my $enzyme_rc1 = reverse($enzyme1);
	$current_loc = index($seq, $enzyme_rc1, 0);
	while ($current_loc != -1){
		if($current_loc + $elength1 < $slength){
			# $result{$current_loc} = 2;	#records position of base just after "N"
			#$result{$current_loc-1} = 2;	#records position at base "N" NCTTCTCG
			$result{$current_loc + $elength1} = 2;	#records position at base "G" NCTTCTCG
		}
		$current_loc = index($seq, $enzyme_rc1, $current_loc + 1);
	}
	
	# Find the reverse(second enzyme) in the forward strand, starting from the first nucleotide!!!
	# my $enzyme_rc2 = reverse($enzyme2);
	# $current_loc = index($seq, $enzyme_rc2, 0);
	# while ($current_loc != -1){
		# if($current_loc + $elength2 < $slength){
			# # $result{$current_loc} = 2;	#records position of base just after "N"
			# #$result{$current_loc-1} = 2;	#records position at base "N" NCTTCTCG
			# $result{$current_loc + $elength2} = 2;	#records position at base "G" NCTTCTCG
		# }
		# $current_loc = index($seq, $enzyme_rc2, $current_loc + 1);
	# }
	
	# Find the rc(first enzyme) in the forward strand, staring from the first nucleotide!!!
	$enzyme_rc1 =~ tr/ACGTUN/TGCAAN/;
	$current_loc = index($seq, $enzyme_rc1, 0);
	while ($current_loc != -1){
		# if($current_loc - 1 >= 0){	#.. we're not searching backwards... 
		if($current_loc + $elength1 < $slength){
			# $result{$current_loc} = -1;	#records position of base just before "N"
			#$result{$current_loc-1} = -1;	#records position at base "N" NGAAGAGC
			$result{$current_loc + $elength1} = -1;	#records position at base "C" NGAAGAGC
		}
		$current_loc = index($seq, $enzyme_rc1, $current_loc + 1);
	}
	
	# Find the rc(enzymes) in the forward strand, staring from the first nucleotide!!!
	# $enzyme_rc2 =~ tr/ACGTUN/TGCAAN/;
	# $current_loc = index($seq, $enzyme_rc2, 0);
	# while ($current_loc != -1){
		# # if($current_loc - 1 >= 0){	#.. we're not searching backwards... 
		# if($current_loc + $elength2 < $slength){
			# # $result{$current_loc} = -1;	#records position of base just before "N"
			# #$result{$current_loc-1} = -1;	#records position at base "N" NGAAGAGC
			# $result{$current_loc + $elength2} = -1;	#records position at base "C" NGAAGAGC
		# }
		# $current_loc = index($seq, $enzyme_rc2, $current_loc + 1);
	# }
	
			
	return %result;
}

sub complement {
	my $seq = shift;
	#$seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

sub calcFragilesSitesBED {
	# takes in AoH with nick sites sorted by position and the bp threshold between nick sites to count as a potential fragile site
	my (@nickSites) = @_;
	my @fsitesBed;
	my $fsiteCount=0;
	#my $elength = length($enzyme); 	#enzyme length

	for(my $i = 0; $i < scalar(@nickSites) - 1; $i++){
		my $firstNick = $nickSites[$i];
		my $secondNick = $nickSites[($i+1)];
			
		Dumper($firstNick);
		Dumper($secondNick);
	
		my ($pos1,$strand1) = each %$firstNick;	
		my ($pos2,$strand2) = each %$secondNick;		
				
		my $curtype; 
		# Want pos1 and pos2 to be the rightmost and leftmost positions of the adjacent nick sites
		if(($strand1 == 1) && ($strand2 == -1) && ($pos2 - $pos1 <= $bp_t1)) {
			# Type I
			# eva chan:  $pos1 and $pos2 are already cloest to each other
			# $pos1 = $pos1 - $elength;
			# $pos2 = $pos2 + $elength;
			$curtype = "TypeI";
		}
		elsif(($strand1 == -1) && ($strand2 == 1) && ($pos2 - $pos1 <= $bp_t2)) {
			# Type II
			# eva chan:: Because the nick sites have been sorted, we would expect the motif on the reverse strand to be encontered first (firstNick) in Type II	
			# $pos1 = $pos1 - $elength;
			# $pos2 = $pos2 + $elength;
			#$pos1 = $pos1 + $elength; 	#reverse strand
			#$pos2 = $pos2 - $elength; 	#forward strand
			$curtype = "TypeII";
		}
		elsif(($strand1 == 1) && ($strand2 == 2) && ($pos2 - $pos1 <= $bp_t3)) {
			# Type III
			# eva chan:  $pos1 and $pos2 are already cloest to each other
			# $pos1 = $pos1 - $elength;
			# $pos2 = $pos2 + $elength;
			$curtype = "TypeIII";
		}
		elsif(($strand1 == 2) && ($strand2 == 1) && ($pos2 - $pos1 <= $bp_t4)) {
			# Type IV
			# $pos1 = $pos1 - length($enzyme);
			# $pos2 = $pos2 + length($enzyme);
			#$pos1 = $pos1 + $elength;
			#$pos2 = $pos2 - $elength;
			$curtype = "TypeIV";
		} else {
			next;
		}
		$pos1 = int($pos1);
		$pos2 = int($pos2);
		my $lineBed = "$count\t$pos1\t$pos2\t$curtype";
		push @fsitesBed, $lineBed;
		$fsiteCount++;
	}
	print "\t\tTotal potential fragile sites found: $fsiteCount\n";
	$fsitesCountTotal = $fsitesCountTotal + $fsiteCount;
	return (@fsitesBed);
}






