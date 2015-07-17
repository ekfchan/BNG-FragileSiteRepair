#!/usr/bin/perl

# Usage: perl calcFragilesSites.pl <input FASTA> [output.bed]
# Example: perl calcFragileSites.pl hg19_chromosome.fa hg19_chromosome_bbspqi.bed

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
if( scalar(@ARGV)<1 or scalar(@ARGV)>2 ) {
	print "usage: $0 <input.fa> [output.bed]\n"; 
	exit 0; 
}

# distance in basepairs between nicks to consider a potential fragile site
my $bp_t1 = 500; #maximum distance between nicks on opposite strands (moving towards each other) to count as TypeI
my $bp_t2 = 500; #maximum distance between nicks on opposite strands (moving away from each other) to count as TypeII
my $bp_t3 = 500; #maximum distance between nicks on same strands (moving towards each other) to count as TypeIII
my $bp_t4 = 500; #maximum distance between nicks on same strands (moving away from each other) to count as TypeIV

my $enzyme = "GCTCTTC"; #BspQI... motif is actually "GCTCTTCN" with nicking occuring after the "N"

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
if( defined $ARGV[1] ) { 
	$outputfile = realpath($ARGV[1]); 
} else {
	$outputfile = realpath("./$basename.fsites.bed");  #if output file not given use $basename as prefix 
}
if( -e $outputfile ) { 
	die "ERROR: Output file ", $outputfile, " already exists. Please specify alternate filename. $!\n"; 
}
print "Output BED: $outputfile\n\n";

print "\tTypeI fragile site threashold (bp): $bp_t1\n";
print "\tTypeII fragile site threashold (bp): $bp_t2\n";
print "\tTypeIII fragile site threashold (bp): $bp_t3\n";
print "\tTypeIV fragile site threashold (bp): $bp_t4\n\n";


## Read input fasta and find nicks and fragiles sites
#my $seqin = Bio::SeqIO->new( -format => 'Fasta' , -file => "$fasta");
my $seqin = Bio::SeqIO->new( -file => "$fasta");
my $count=1;
my @fsitesBed;
my $fsitesHeaderBed = "#Nickase: $enzyme\n#CMapId\tStart\tEnd\tType";
while((my $seqobj = $seqin->next_seq())) {   #for each sequence in FASTA (e.g. each contig)
	my $seq = $seqobj->seq();
	my $id  = $seqobj->display_id();
	my $length = $seqobj->length();

	# find nick sites
	print "\tHeader: $id Length: $length\n";
	my %result = find_nick_sites($seq,$enzyme);	
	
	# sort nick sites by genomic position, creating an ordered array of hash references
	my @nick_sites;
	foreach my $key (sort{$a <=> $b} keys %result){
		my %hash = ($key, $result{$key});
		push @nick_sites, \%hash;
	}
	print "\t\tTotal nick sites found: ", scalar(@nick_sites), "\n"; 

	#calc fragile sites for each seq contig
	my (@fsiteBed) = calcFragilesSitesBED(@nick_sites);
	push @fsitesBed, @fsiteBed;	
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
	my ($seq, $enzyme) = @_;
	my %result;
	my $slength = length($seq);
	my $elength = length($enzyme);
	
	# Find the enzymes in the forward strand, starting from the first nucleotide!!!
	my $current_loc = index($seq, $enzyme, 0);
	while ($current_loc != -1){
		if($current_loc + $elength < $slength){
			# $result{$current_loc + $elength + 1} = 1;	#records position of base after "N" (not part of motif)
			$result{$current_loc + $elength} = 1;	#records position at base "N"
		}
		$current_loc = index($seq, $enzyme, $current_loc + 1);
	}
	
	# Find the reverse(enzymes) in the forward strand, starting from the first nucleotide!!!
	my $enzyme_rc = reverse($enzyme);
	$current_loc = index($seq, $enzyme_rc, 0);
	while ($current_loc != -1){
		if($current_loc + $elength < $slength){
			# $result{$current_loc} = 2;	#records position of base just before "N"
			$result{$current_loc-1} = 2;	#records position at base "N"
		}
		$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
	}
	
	# Find the rc(enzymes) in the forward strand, staring from the first nucleotide!!!
	$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
	$current_loc = index($seq, $enzyme_rc, 0);
	while ($current_loc != -1){
		# if($current_loc - 1 >= 0){	#.. we're not searching backwards... 
		if($current_loc + $elength < $slength){
			# $result{$current_loc} = -1;	#records position of base just before "N"
			$result{$current_loc-1} = -1;	#records position at base "N"
		}
		$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
	}
			
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
	my $elength = length($enzyme); 	#enzyme length

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
			$pos1 = $pos1 + $elength; 	#reverse strand
			$pos2 = $pos2 - $elength; 	#forward strand
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
			$pos1 = $pos1 + $elength;
			$pos2 = $pos2 - $elength;
			$curtype = "TypeIV";
		} else {
			next;
		}
		my $lineBed = "$count\t$pos1\t$pos2\t$curtype";
		push @fsitesBed, $lineBed;
		$fsiteCount++;
	}
	print "\t\tTotal potential fragile sites found: $fsiteCount\n";
	$fsitesCountTotal = $fsitesCountTotal + $fsiteCount;
	return (@fsitesBed);
}






