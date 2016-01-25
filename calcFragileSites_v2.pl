#!/usr/bin/perl

# Usage: perl calcFragileSites.pl --fasta <input FASTA> [--output <output.bed>] [--buffer <sequence buffer in bp>] [--random] [--enzyme <sequence to use to calculate fragile sites> [--agressive <calculate TypeIII, TypeIV, TypeV fragile sites>]]

use strict;
use warnings;
use Cwd qw(realpath abs_path);
use File::Basename;
#use File::Copy;
use File::Path qw(make_path);
use Bio::Seq;
use Bio::SeqIO;
#use List::MoreUtils qw(uniq);
use Data::Dumper;
use Getopt::Long qw(:config bundling); 
#use Number::Bytes::Human qw(format_bytes);

# usage statement
my %inputs = (); 
GetOptions( \%inputs, 'fasta=s', 'output=s', 'buffer=i', 'random', 'enzyme:s', 'aggressive'); 
if ( !exists $inputs{fasta} ) {
	print "Usage: perl calcFragileSites.pl --fasta <input FASTA> [--output <output.bed>] [--buffer <sequence buffer in bp>] [--random] [--enzyme <sequence to use to calculate fragile sites> [--agressive <calculate TypeIII, TypeIV, TypeV fragile sites>]]\n\n";
	exit 0;
}
foreach my $key ("fasta","output") {
	if (exists $inputs{$key} and $inputs{$key} =~ /[a-z]/i) {
		$inputs{$key} = abs_path($inputs{$key});
	}
}


#if( scalar(@ARGV)<1 or scalar(@ARGV)>3 ) {
	#print "usage: $0 <input.fa> [output.bed] [sequence buffer in bp]\n"; 
	#exit 0; 
#}

print "\n";
print qx/ps -o args $$/;
print "\n";

# distance in basepairs between nicks to consider a potential fragile site
my $bp_t1 = 1000; #maximum distance between nicks on opposite strands (moving towards each other) to count as TypeI
my $bp_t2 = 1000; #maximum distance between nicks on opposite strands (moving away from each other) to count as TypeII
my $bp_t3 = 1000; #maximum distance between nicks on same strands (moving towards each other) to count as TypeIII
my $bp_t4 = 1000; #maximum distance between nicks on same strands (moving away from each other) to count as TypeIV
my $bp_t5 = 100; #maximum distance between nicks on same strands (moving in the same direction) to count as TypeV
my $bp_t6 = 1000; #maximum distance between nicks of more than 1 strand type to count as TypeVI

my @enzymes;
#testing purposes for randomness
if (exists $inputs{random}) { 
	my @chars = ("A","C","T","G");
	$enzymes[0] .= $chars[rand @chars] for 1..7;
}
else {
	push @enzymes, $inputs{enzyme};
	#@enzymes = ("GCTCTTC"); #BspQI
	#@enzymes = ("CACGAG"); #BssSI
	#@enzymes = ("GCTCTTC","CACGAG"); #BspQI and BssSI dual-nick
}

#my $enzyme1 = "GCTCTTC"; #BspQI... motif is actually "GCTCTTCN" with nicking occuring after the "N"
#my $enzyme2 = "CACGAG"; #BssSI

my $outputfile; 
my $basename;
my $fasta;
my $fsitesCountTotal=0;

print "\n";

# check fasta input variable
if (defined $inputs{fasta} && -e $inputs{fasta}) {
	$fasta = realpath($inputs{fasta});
	#print "Calculating fragiles sites for FASTA file: $fasta\n";
	$basename = basename($fasta);
	$basename =~ s/\.[^.]+$//;
	#print "\tFile filename = $basename\n";
	}
else {
	die "ERROR: Input FASTA file $inputs{fasta} not provided or does not exist! $!\n";
} 
print "Input FASTA: $fasta\n";	

# check optional output variable
if( defined $inputs{output} and !-d $inputs{output}) { 
	$outputfile = realpath($inputs{output}); 
} 
elsif ( defined $inputs{output} and -d $inputs{output}) {
	#$outputfile = $inputs{output}."/".$basename."_".uc($enzyme1)."_".uc($enzyme2).".fsites.withSeqs.bed";
	$outputfile = $inputs{output}."/".$basename."_".uc(join("_",@enzymes));
}
else {
	#$outputfile = realpath("./$basename.".uc($enzyme1)."_".uc($enzyme2).".fsites.withSeqs.bed");  #if output file not given use $basename as prefix 
	$outputfile = realpath("./$basename.".uc(join("_",@enzymes)));
}
if( -e $outputfile ) { 
	die "ERROR: Output file ", $outputfile, " already exists. Please specify alternate filename. $!\n"; 
}

#check optional input sequence buffer variable
my $getSeq=0;
my $buffer=50;	
if (defined $inputs{buffer}) {
	$getSeq=1;
	if (int($inputs{buffer})>0) {
		$buffer = int($inputs{buffer}); 
	}
	$outputfile = $outputfile."_fsites.withSeqs".$buffer.".bed";
}
else {
	$outputfile = $outputfile."_fsites.bed";
}

print "Calculating potential fragile sites for $inputs{fasta}\n";
print "Using enzyme(s): ".join(" ",@enzymes)."\n\n"; 

print "Output BED: $outputfile\n\n";

if ($inputs{aggressive}) {
	print "Aggressive fragile site calculation enabled. Calculating TypeIII and TypeIV...\n";

	print "\tTypeI fragile site threashold (bp): $bp_t1\n";
	print "\tTypeII fragile site threashold (bp): $bp_t2\n";
	print "\tTypeIII fragile site threashold (bp): $bp_t3\n";
	print "\tTypeIV fragile site threashold (bp): $bp_t4\n";
	print "\tTypeV fragile site threashold (bp): $bp_t5\n";
	print "\tTypeVI fragile site threashold (bp): $bp_t6\n";
	print "\n";
}
else {
	print "\tTypeI fragile site threashold (bp): $bp_t1\n";
	print "\tTypeII fragile site threashold (bp): $bp_t2\n";
	#print "\tTypeIII fragile site threashold (bp): $bp_t3\n";
	#print "\tTypeIV fragile site threashold (bp): $bp_t4\n";
	#print "\tTypeV fragile site threashold (bp): $bp_t5\n";
	#print "\tTypeVI fragile site threashold (bp): $bp_t6\n";
	print "\n";
}

if ($getSeq==1) {
	print "\tExtracting sequences +/- $buffer"."bp from fragile site start/end positions\n";
	print "\n";
}

## Read input fasta and find nicks and fragiles sites
#my $seqin = Bio::SeqIO->new( -format => 'Fasta' , -file => "$fasta");
my $seqin = Bio::SeqIO->new( -file => "$fasta");
my $count=1;
my @fsitesBed;
my $fsitesHeaderBed="";
for (my $i=0; $i<scalar(@enzymes); $i++) {
	$fsitesHeaderBed = $fsitesHeaderBed."#Nickase ".($i+1).": $enzymes[$i]\n";
}
if ($getSeq==1) {
	$fsitesHeaderBed = $fsitesHeaderBed."#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence";
}
else {
	$fsitesHeaderBed = $fsitesHeaderBed."#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba";
}

#my $fsitesHeaderBed = "#Nickase 1: $enzyme1\n#Nickase 2: $enzyme2\n#CMapId\tStart\tEnd\tType\tSequence";
#my $fsitesHeaderBed = "#Nickase 1: $enzyme1\n#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence";
while((my $seqobj = $seqin->next_seq())) {   #for each sequence in FASTA (e.g. each contig)
	my $seq = $seqobj->seq();
	my $id  = $seqobj->display_id();
	my $length = $seqobj->length();

	# find nick sites
	print "\tHeader: $id Length: $length\n";
	
	#my %result = find_nick_sites($seq,$enzyme1,$enzyme2);
	#my %result = find_nick_sites($seq,$enzyme1);	
	my %result;
	if ($inputs{aggressive}) { 	
		%result = find_nick_sites_aggressive($seq,\@enzymes);
	}
	else {
		%result = find_nick_sites($seq,\@enzymes);
	}
	
	# sort nick sites by genomic position, creating an ordered array of hash references
	my @nick_sites;
	
	###DEBUG###
	my %nickStrandCounts;
	
	foreach my $key (sort{$a <=> $b} keys %result){
		my %hash = ($key, $result{$key});
		
		###DEBUG###
		$nickStrandCounts{$result{$key}}++;
		
		push @nick_sites, \%hash;
	}
	print "\t\tTotal potential nick sites found: ", scalar(@nick_sites), "\n"; 
	
	###DEBUG###
	#foreach my $key (sort keys %nickStrandCounts) {
	#	print "\t\t\tClass: $key Count: $nickStrandCounts{$key}\n";
	#} 

	#calc fragile sites for each seq contig
	my (@fsiteBed) = calcFragilesSitesBED(@nick_sites);
	
	my @fsiteBedSeq;
	foreach my $line (@fsiteBed) {
		my $lineOut="";
		my @s = split("\t",$line);

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
		elsif ($s[3] =~ m/TypeV/i) {
			$itemRgb = "0,255,0,255"; }
		elsif ($s[3] =~ m/TypeVI/i) {
			$itemRgb = "0,0,255,255"; }
		else {
			$itemRgb = "255,255,255,255"; }
			
		#build BED lineOut
		$lineOut = "$s[0]\t$s[1]\t$s[2]\t$s[3]\t0\t0\t$s[1]\t$s[2]\t$itemRgb";
		
		#if getting fasta sequences, extract sequence
		if ($getSeq==1) {
			#extract sequence +/- $inputs{buffer} or 500bp			
			my $start = $s[1]-$buffer;
			if ($start <1) {
				$start = 1;
			}
			my $end = $s[2]+$buffer;
			if ($end >$length) {
				$end = $length-1;
			}
			my $seq = $seqobj->subseq($start,$end);
			
			$lineOut = $lineOut."\t$seq";
		}
		
		push @fsiteBedSeq, $lineOut;
		
		
				
		#push @fsitesBed, @fsiteBed;	
		
	}
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
	
	my ($seq, $enzymes_ref) = @_;
	$seq = uc($seq);
	my $seq_comp = complement($seq);
	my @enzymes = @{$enzymes_ref};
	my %result;
	my $slength = length($seq);	
	
	foreach my $enzyme (@enzymes) {
		$enzyme = uc($enzyme);
		my $elength = length($enzyme);
		my %classCount;
		
		# Find the first enzyme in the forward strand, starting from the first nucleotide!!!
		my $runningCount=0;
		my $current_loc = index( $seq, $enzyme, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc}) {
					$result{$current_loc} = 1;	#records position at base "G" GCTCTTCN
					$runningCount++;
				}
				else {
					$result{$current_loc} = 1;
					#print "\t\tNickPos: $current_loc\n";
					$runningCount++;
				}

				#my $loc = format_bytes($current_loc, bs => 1000);
				#print "\t\tNickPos: $loc\n";
				#if ($runningCount % 10 == 0) {
				#	print "\t\tRunning nick count: $runningCount NickPos: $current_loc\n";
				#}
			}
		$current_loc = index( $seq, $enzyme, $current_loc + 1);
		}	
		
		## find the reverse(enzyme)) in the forward strand 5' -> 3' 
		my $enzyme_rc = reverse($enzyme);
		#$current_loc = index($seq, $enzyme_rc, 0);
		#while ($current_loc != -1){
			#if($current_loc + $elength < $slength){
				#if (!exists $result{$current_loc + $elength}) {
					#$result{$current_loc + $elength} = 2;	#records position at base "G" NCTTCTCG
				#}
				#else {
					#$result{$current_loc + $elength} = 3;
				#}
			#}
			#$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
		#}
		
		# Find the rc(first enzyme) in the forward strand, staring from the first nucleotide!!!
		#$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
		$enzyme_rc = complement($enzyme_rc);
		$current_loc = index( $seq, $enzyme_rc, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc + $elength}) {
					$result{$current_loc + $elength} = -1;	#records position at base "C" NGAAGAGC
					#print "\t\tRevCompNickPos: $current_loc\n";
				}
				else {
					$result{$current_loc + $elength} = -1;
					#print "\t\tRevCompNickPos: $current_loc\n";
				}
			}
			$current_loc = index( $seq, $enzyme_rc, $current_loc + 1);
		}
		
		## Find the complement(enzyme) in the forward strand 5' -> 3', starting from the first nucleotide!!!		
		#$enzyme_rc = $enzyme;
		#$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
		#$current_loc = index($seq, $enzyme_rc, 0);
		#while ($current_loc != -1){
			#if($current_loc + $elength < $slength){
				#if (!exists $result{$current_loc}) {
					#$result{$current_loc} = -2;	#records position at base "C" CGAGAAG
				#}
				#else {
					#$result{$current_loc} = 3;
				#}
			#}
			#$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
		#}
		
		#print "\t\tNickase Enzyme: $enzyme\n";
		#for my $key (sort{$b <=> $a} keys %classCount) {
			#print "\t\t\tNick class: $key count: $classCount{$key}\n";
		#}
	}

	return %result;
}

sub find_nick_sites_aggressive{
	
	my ($seq, $enzymes_ref) = @_;
	$seq = uc($seq);
	my $seq_comp = complement($seq);
	my @enzymes = @{$enzymes_ref};
	my %result;
	my $slength = length($seq);	
	
	foreach my $enzyme (@enzymes) {
		$enzyme = uc($enzyme);
		my $elength = length($enzyme);
		my %classCount;
		
		# Find the first enzyme in the forward strand, starting from the first nucleotide!!!
		my $current_loc = index( $seq, $enzyme, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc}) {
					$result{$current_loc} = 1;	#records position at base "G" GCTCTTCN
				}
				else {
					$result{$current_loc} = 3;
				}
			}
		$current_loc = index( $seq, $enzyme, $current_loc + 1);
		}	
		
		## find the reverse(enzyme)) in the forward strand 5' -> 3' 
		my $enzyme_rc = reverse($enzyme);
		$current_loc = index( $seq, $enzyme_rc, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc + $elength}) {
					$result{$current_loc + $elength} = 2;	#records position at base "G" NCTTCTCG
				}
				else {
					$result{$current_loc + $elength} = 3;
				}
			}
			$current_loc = index( $seq, $enzyme_rc, $current_loc + 1);
		}
		
		# Find the rc(first enzyme) in the forward strand, staring from the first nucleotide!!!
		#$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
		$enzyme_rc = complement($enzyme_rc);
		$current_loc = index( $seq, $enzyme_rc, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc + $elength}) {
					$result{$current_loc + $elength} = -1;	#records position at base "C" NGAAGAGC
				}
				else {
					$result{$current_loc + $elength} = 3;
				}
			}
			$current_loc = index( $seq, $enzyme_rc, $current_loc + 1);
		}
		
		# Find the complement(enzyme) in the forward strand 5' -> 3', starting from the first nucleotide!!!		
		$enzyme_rc = $enzyme;
		#$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
		$enzyme_rc = complement($enzyme_rc);
		$current_loc = index( $seq, $enzyme_rc, 0);
		while ($current_loc != -1){
			if($current_loc + $elength < $slength){
				if (!exists $result{$current_loc}) {
					$result{$current_loc} = -2;	#records position at base "C" CGAGAAG
				}
				else {
					$result{$current_loc} = 3;
				}
			}
			$current_loc = index( $seq, $enzyme_rc, $current_loc + 1);
		}
		
		#print "\t\tNickase Enzyme: $enzyme\n";
		#for my $key (sort{$b <=> $a} keys %classCount) {
			#print "\t\t\tNick class: $key count: $classCount{$key}\n";
		#}
	}

	return %result;
}

sub complement {
	my $seq = shift;
	$seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	#$seq =~ tr/ACGTacgt/TGCAtgca/;
	#$seq =~ tr/ACGTUN/TGCAAN/;
	return $seq;
}

sub calcFragilesSitesBED {
	# takes in AoH with nick sites sorted by position and the bp threshold between nick sites to count as a potential fragile site
	my (@nickSites) = @_;
	my @fsitesBed;
	my $fsiteCount=0;
	#my $elength = length($enzyme); 	#enzyme length
	my %fsiteTypeCount;

	for(my $i = 0; $i < scalar(@nickSites) - 1; $i++){
		my $firstNick = $nickSites[$i];
		my $secondNick = $nickSites[($i+1)];
			
		Dumper($firstNick);
		Dumper($secondNick);
	
		my ($pos1,$strand1) = each %$firstNick;	
		my ($pos2,$strand2) = each %$secondNick;		
				
		my $curtype; 

		#maximum distance between nicks on opposite strands (moving towards each other) to count as TypeI
		if(($strand1 == 1) && ($strand2 == -1) && ($pos2 - $pos1 <= $bp_t1)) {
			$curtype = "TypeI";
		}
		#maximum distance between nicks on opposite strands (moving away from each other) to count as TypeII
		elsif(($strand1 == -1) && ($strand2 == 1) && ($pos2 - $pos1 <= $bp_t2)) {
			$curtype = "TypeII";
		}
		
		#maximum distance between nicks on same strands (moving towards each other) to count as TypeIII
		elsif(($strand1 == 1) && ($strand2 == 2) && ($pos2 - $pos1 <= $bp_t3)) {
			$curtype = "TypeIII";
		}
		elsif(($strand1 == -2) && ($strand2 == -1)  && ($pos2 - $pos1 <= $bp_t3)) {
			$curtype = "TypeIII";
		}
		
		#maximum distance between nicks on same strands (moving away from each other) to count as TypeIV
		elsif(($strand1 == 2) && ($strand2 == 1) && ($pos2 - $pos1 <= $bp_t4)) {
			$curtype = "TypeIV";
		}
		elsif(($strand1 == -1) && ($strand2 == -2) && ($pos2 - $pos1 <= $bp_t4)) {
			$curtype = "TypeIV";
		}
		
		#maximum distance between nicks on same strands (moving in the same direction) to count as TypeV
		elsif ($inputs{aggressive}) {
			if(($strand1 == 1) && ($strand2 == 1) && ($pos2 - $pos1 <= $bp_t5)) {
				$curtype = "TypeV";
			}
			elsif(($strand1 == -1) && ($strand2 == -1) && ($pos2 - $pos1 <= $bp_t5)) {
				$curtype = "TypeV";
			}
			elsif(($strand1 == -2) && ($strand2 == -2) && ($pos2 - $pos1 <= $bp_t5)) {
				$curtype = "TypeV";
			}
			elsif(($strand1 == 2) && ($strand2 == 2) && ($pos2 - $pos1 <= $bp_t5)) {
				$curtype = "TypeV";
			}
		}

		#maximum distance between nicks of more than 1 strand type to count as TypeVI
		elsif((($strand1 == 3 || $strand2 == 3)) && ($pos2 - $pos1 <= $bp_t6)) {
			$curtype = "TypeVI";
		} 		
		else {
			next;
		}
		$pos1 = int($pos1);
		$pos2 = int($pos2);
		my $lineBed = "$count\t$pos1\t$pos2\t$curtype";
		push @fsitesBed, $lineBed;
		$fsiteCount++;		
		$fsiteTypeCount{$curtype}++;
	}
	print "\t\tTotal potential fragile sites found: $fsiteCount\n";
	foreach my $type (sort keys %fsiteTypeCount) {
		print "\t\t\t$type: $fsiteTypeCount{$type}\n";
	}
	$fsitesCountTotal = $fsitesCountTotal + $fsiteCount;
	return (@fsitesBed);
}






