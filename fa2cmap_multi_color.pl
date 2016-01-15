#!/usr/bin/perl -w
############################################################################
# File: fa2cmap_multi_color.pl                                             #
# Date: 06/05/2015                                                         #
# Purpose: Transform fasta file to BioNano cmap file format (Color-aware)  #
#                                                                          #
# Author: Xiang Zhou, Computational Biologist                              #
# Email : xzhou@bionanogenomics.com                                        #
# Affiliation: Research Department, BioNano Genomics Inc.                  #
#                                                                          #
# Usage:                                                                   #
#   fa2cmap_multi_color.pl [options] <Args>                                #
# Options:                                                                 #
#   -h    : This help message                                              #
#   -v    : Verbose output  (Default: OFF)                                 #
#   -i    : Input fasta file  (Required)                                   #
#   -o    : Output folder  (Default: the same as the input file)           #
#   -e    : Names/sequences of the enzymes followed by channel # (multiple)#
#   -m    : Filter: Minimum labels  (Default: 0)                           #
#   -M    : Filter: Minimum size (Kb)  (Default: 0)                        #
#                                                                          #
# NOTE: CMAP index is 1-based, and is color-aware.                         #
############################################################################

use strict;
use warnings;
use POSIX;
use File::Spec;
use File::Basename;
use File::Spec::Functions;
#use List::MoreUtils qw(uniq);
#use Getopt::Std;
use Getopt::Long qw(:config no_ignore_case);
#use File::Slurp qw(edit_file_lines);
#$Getopt::Std::STANDARD_HELP_VERSION = 1;

sub Init;
sub Usage;
sub Find_enzymes;
sub Print_cmap_header;
sub Generate_cmap;
sub Get_and_set_enzyme_sequence;
sub is_enzyme_name;
sub is_nt_sequence;
sub Uniq;

my %enzyme = (
	"BspQI" => "GCTCTTC",
	"BbvCI" => "CCTCAGC",
	"BsmI"  => "GAATGC",
	"BsrDI" => "GCAATG",
	"bseCI" => "ATCGAT"
	
	# You can add more enzymes here ...
	
);
my (@enzymes, @enzyme_sequences, @color, @color_uniq);
my ($min_labels, $min_length) = (0, 0);
my ($FASTA, $CMAP, $KEY);
my ($filename, $filename_key);

#my %opts;
my ($help, $verbose, $keyfile, $input, $output);
my (@loci, @loci_tmp);
my %color;
my $count = 0;
my @num_nicks;  # For each channel
my @num_nicks_global;  # For each channel
my @density;  # For each channel
my $density_global;
my ($total_length_original, $total_length_filtered) = (0, 0);
my $total_nicks = 0;
my $num_cmaps = 0;

my ($command, $seq, $tmp);
my ($fastaHeader, $fastaHeader_previous);
my ($A, $C, $G, $T, $N, $global_N, $global_GC, $global_ACGT, $global_ACGTN);
my ($N_percentage, $GC_percentage, $global_N_percentage, $global_GC_percentage);

Init();

for(my $i = 0; $i < @color_uniq; $i++){
	@num_nicks_global[$color_uniq[$i]] = 0;
}

open($FASTA, $input) || die ("ERROR: Can't open $input: $!\n");
$filename = $input;
$filename_key = $input;

for(my $i = 0; $i < @enzymes; $i++){
	$tmp .= "_$enzymes[$i]";
}

# If output folder is not defined, use the input folder
if($output){
	#if(substr($output, -1) eq "/"){
	#	$filename = $output . (split("/", $input))[-1];
	#	$filename_key = $output . (split("/", $input))[-1];
	#}
	#else{
	#	$filename = $output . "/" . (split("/", $input))[-1];
	#	$filename_key = $output . "/" . (split("/", $input))[-1];
	#}
	
	# TODO: Convert ".." to absolute path!
	my ($nothing, $out_dir) = fileparse($output);
	#print "\nOUTPUT_DIR\n", $out_dir, "\n\n";
	
	$filename = basename($filename);
	$filename = catfile($out_dir, $filename);
	#print "\nFILENAME\n", $filename, "\n\n";
	#print "\nFULL_PATH_FILENAME\n", $filename, "\n\n";
	
	$filename_key = basename($filename_key);
	$filename_key = catfile($out_dir, $filename_key);
	#print "\nFILENAME_KEY\n", $filename_key, "\n\n";
	#print "\nFULL_PATH_FILENAME_KEY\n", $filename_key, "\n\n";
}
$filename =~ s/(\S+)\.\w+$/$1$tmp\.cmap/;
$filename_key =~ s/(\S+)\.\w+$/$1$tmp\_key.txt/;

open($KEY, ">$filename_key") || die ("ERROR: Can't open $filename_key: $!\n");
print $KEY("# CMAP = ", File::Spec -> rel2abs($filename), "\n");
print $KEY("# filter: Minimum Labels = $min_labels\n");
print $KEY("# filter: Minimum Size (Kb) = $min_length\n");
print $KEY("CompntId\tCompntName\tCompntLength\n");

# Verbose output?
if($verbose){
	print("Input file:\t", File::Spec -> rel2abs($input), "\n");
	print("Output file:\t", File::Spec -> rel2abs($filename), "\n");
	print("ID keyfile:\t", File::Spec -> rel2abs($filename_key), "\n");
	for(my $i = 0; $i < @enzymes; $i++){
		print("Enzyme ", $i+1, " (Channel $color[$i]):\t$enzymes[$i]($enzyme_sequences[$i])\n");
	}
	print("\n");
	
	print("CMapId\tChannel\tLength(MB)\tFreq(nicks/100KB)\tGC%\tN%\n");
	print("------------------------------------------------------------------\n");
}

Print_cmap_header($filename);

while(my $line = <$FASTA>){
	chomp $line;
	$line =~ s/\r//g;
	
	if($line =~ /^>/){
		$fastaHeader_previous = $fastaHeader;
		$fastaHeader = substr($line, 1);
		
		if($count != 0){
			$seq = uc($seq);
			$A = ($seq =~ tr/A/A/);
			$C = ($seq =~ tr/C/C/);
			$G = ($seq =~ tr/G/G/);
			$T = ($seq =~ tr/T/T/);
			$N = ($seq =~ tr/N/N/);
			$global_N += $N;
			$global_GC += ($C+$G);
			$global_ACGT += ($A+$C+$G+$T);
			$global_ACGTN += ($A+$C+$G+$T+$N);
			
			$total_length_original += length($seq);
			
			@loci_tmp = ();
			%color = ();
			for(my $i = 0; $i < @color_uniq; $i++){
				@num_nicks[$color_uniq[$i]] = 0;
			}
			for(my $i = 0; $i < @enzyme_sequences; $i++){
				my @nicks_tmp = Find_enzymes($seq, $enzyme_sequences[$i], $color[$i]);
				@loci_tmp = (@loci_tmp, @nicks_tmp);
			}
			
			# Remove duplicated values!
			@loci = Uniq(@loci_tmp);
			
			if(scalar(@loci) >= $min_labels && length($seq) >= $min_length * 1000){
				Generate_cmap($filename, $count, \@loci, length($seq));
				#if(scalar(@loci)){
				print $KEY("$count\t$fastaHeader_previous\t", length($seq), "\n");
				#}
				
				my $length_tmp = sprintf("%11.3f", length($seq)/1000000);
				
				if($A+$C+$G+$T == 0){
					$GC_percentage = sprintf("%5.2f", 0);
				}
				else{
					$GC_percentage = sprintf("%5.2f", ($C+$G)/($A+$C+$G+$T)*100);
				}
				if($N == 0){
					$N_percentage = sprintf("%6.2f", 0);
				}
				else{
					$N_percentage = sprintf("%6.2f", $N/length($seq)*100);
				}
				
				# Verbose output?
				if($verbose){
					for(my $i = 0; $i < @color_uniq; $i++){
						$density[$count] = sprintf("%7.3f", $num_nicks[$color_uniq[$i]]/length($seq) * 100000);
						print("$count\t$color_uniq[$i]\t$length_tmp\t$density[$count]\t$GC_percentage\t$N_percentage\n");
					}
				}
				
				# NOTE: No duplicate sites!!!
				$total_nicks += @loci;
				$total_length_filtered += length($seq);
			}
		}
		
		$seq = "";
		$count++;
	}
	else{
		$seq .= uc($line);
	}
}
$seq = uc($seq);
$A = ($seq =~ tr/A/A/);
$C = ($seq =~ tr/C/C/);
$G = ($seq =~ tr/G/G/);
$T = ($seq =~ tr/T/T/);
$N = ($seq =~ tr/N/N/);
$global_N += $N;
$global_GC += ($C+$G);
$global_ACGT += ($A+$C+$G+$T);
$global_ACGTN += ($A+$C+$G+$T+$N);

$total_length_original += length($seq);

@loci_tmp = ();
%color = ();
for(my $i = 0; $i < @color_uniq; $i++){
	@num_nicks[$color_uniq[$i]] = 0;
}
for(my $i = 0; $i < @enzyme_sequences; $i++){
	my @nicks_tmp = Find_enzymes($seq, $enzyme_sequences[$i], $color[$i]);
	@loci_tmp = (@loci_tmp, @nicks_tmp);
}

# Remove duplicated values!
@loci = Uniq(@loci_tmp);

if(scalar(@loci) >= $min_labels && length($seq) >= $min_length * 1000){
	Generate_cmap($filename, $count, \@loci, length($seq));
	#if(scalar(@loci)){
	print $KEY("$count\t$fastaHeader\t", length($seq), "\n");
	#}
	
	my $length_tmp = sprintf("%11.3f", length($seq)/1000000);
	
	if($A+$C+$G+$T == 0){
		$GC_percentage = sprintf("%5.2f", 0);
	}
	else{
		$GC_percentage = sprintf("%5.2f", ($C+$G)/($A+$C+$G+$T)*100);
	}
	if($N == 0){
		$N_percentage = sprintf("%6.2f", 0);
	}
	else{
		$N_percentage = sprintf("%6.2f", $N/length($seq)*100);
	}
	# Verbose output?
	if($verbose){
		for(my $i = 0; $i < @color_uniq; $i++){
			$density[$count] = sprintf("%7.3f", $num_nicks[$color_uniq[$i]]/length($seq) * 100000);
			print("$count\t$color_uniq[$i]\t$length_tmp\t$density[$count]\t$GC_percentage\t$N_percentage\n");
		}
	}
	
	# NOTE: No duplicate sites!!!
	$total_nicks += @loci;
	$total_length_filtered += length($seq);
}
close($FASTA);
close($KEY);

#$command = "perl -pi -e 's|N/A|$num_cmaps| if " . '$.' . " == 7' $filename";
$command = "perl -pi -e 's|N/A|$num_cmaps|' $filename";
#print("$command\n");
system("$command");

#edit_file_lines {s |N/A|$num_cmaps| } $filename;

# Verbose output?
if($verbose){
	print("\n================================Summary================================\n");
	print("Total contigs processed: $count\n");
	print("Total contigs generated: $num_cmaps\n");
	
	print("Total length of the original contigs: $total_length_original\n");
	print("Total length of the filtered contigs: $total_length_filtered\n");
	
	for(my $i = 0; $i < @color_uniq; $i++){
		my $density_tmp;
		if($total_length_filtered){
			$density_tmp = sprintf("%6.3f", $num_nicks_global[$color_uniq[$i]]/$total_length_filtered * 100000);
		}
		else{
			$density_tmp = sprintf("%6.3f", 0);
		}
		print("Global nick frequency for channel $color_uniq[$i]:\t$density_tmp nicks/100KB\n");
	}
	
	if($total_length_filtered){
		$density_global = sprintf("%6.3f", $total_nicks/$total_length_filtered * 100000);
	}
	else{
		$density_global = sprintf("%6.3f", 0);
	}
	print("Global nick frequency:\t$density_global nicks/100KB\n");
	
	$global_GC_percentage = sprintf("%5.3f", $global_GC/$global_ACGT*100);
	$global_N_percentage = sprintf("%5.3f", $global_N/$global_ACGTN*100);
	print("Global GC percentage:\t$global_GC_percentage%\n");
	print("Global N percentage:\t$global_N_percentage%\n");
	print("=======================================================================\n");
}

######################################################################
#                           Subroutines                              #
######################################################################
sub Init{
=COMMENT
	my $opt_string = 'hvi:n:s:m:M:';
	if(!getopts("$opt_string", \%opts)){
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
=cut
	my @tmp;
	my $ret = GetOptions(
		'help|h|?'        => \$help,
		'verbose|v'       => \$verbose,
		'input|i=s'       => \$input,
		'output|o=s'      => \$output,
		'enzyme|e=s{1,8}' => \@tmp,
		'min_labels|m:i'  => \$min_labels,
		'min_length|M:i'  => \$min_length,
	);
	
	if(!$ret){
		print("ERROR: Missing or invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	
	Usage() if $help;
	
	if(!$input){
		print("ERROR: Missing parameter(s)! Try -h for more information.\n");
		Usage();
	}
	
	if(!@tmp){
		print("ERROR: Missing parameter(s)! Try -h for more information.\n");
		Usage();
	}
	
	# Total number of values must be even
	if(@tmp%2 != 0){
		print("ERROR: There must be even numbers of -e parameters! Try -h for more information.\n");
		Usage();
	}
	
	for(my $i = 1; $i < @tmp; $i += 2){
		# The values at odd positions must be positive integers
		if($tmp[$i] =~ /^\d+$/){
			$color[($i-1)/2] = $tmp[$i];
		}
		else{
			print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
			Usage();
		}
		
		# The values at even positions must be either name or sequence
		if(is_enzyme_name(\%enzyme, $tmp[$i-1])){
			$enzyme_sequences[($i-1)/2] = Get_and_set_enzyme_sequence(\%enzyme, \$tmp[$i-1]);
			$enzymes[($i-1)/2] = uc($tmp[$i-1]);
		}
		elsif(is_nt_sequence($tmp[$i-1])){
			$enzyme_sequences[($i-1)/2] = uc($tmp[$i-1]);
			$enzymes[($i-1)/2] = uc($tmp[$i-1]);
		}
		else{
			print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
			Usage();
		}
	}
	
	@color_uniq = sort {$a <=> $b} (Uniq(@color));
	#print "[@enzymes]", "\n";
	#print "[@enzyme_sequences]", "\n";
	#print "[@color]", "\n";
	#print "[@color_uniq]", "\n";
}

sub Usage{
	print << "EOF";

Usage: $0 [options] <Args>
Options:
  -h    : This help message
  -v    : Verbose output  (Default: OFF)
  -i    : Input fasta file  (Required)
  -o    : Output folder  (Default: the same as the input file)
  -e    : Names/sequences of the enzymes followed by channel # (Can be multiple)
  -m    : Filter: Minimum labels  (Default: 0)
  -M    : Filter: Minimum size (Kb)  (Default: 0)

NOTE: CMAP index is 1-based, and is color-aware.
EOF
	exit;
}

sub Find_enzymes{
	my ($seq, $enzyme, $color) = @_;
	my @result;
	$seq = uc($seq);
	
	# Find the enzymes in the forward strand, staring from the first nucleotide!!!
	my $current_loc = index($seq, $enzyme, 0);
	while ($current_loc != -1){
		
		# Add the current location into @result IFF
		# 1) it is not in @result, or
		# 2) it is in @result, and the current color is not in the color array.
		if( !defined($color{$current_loc+1}[0]) || !($color ~~ @{$color{$current_loc+1}}[1 .. $color{$current_loc+1}[0]]) ){
			# This is color-aware!
			push(@result, $current_loc+1);
			
			if( !defined($color{$current_loc+1}) ){
				$color{$current_loc+1}[0] = 1;
			}
			else{
				$color{$current_loc+1}[0]++;
			}
			$color{$current_loc+1}[ $color{$current_loc+1}[0] ] = $color;
		}
		
		$current_loc = index( $seq, $enzyme, $current_loc + 1);
	}

	my $enzyme_rc = reverse($enzyme);
	$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
	
	# Find the rc(enzymes) in the forward strand, staring from the first nucleotide!!!
	$current_loc = index( $seq, $enzyme_rc, 0);
	while ($current_loc != -1){
		
		# Add the current location into @result IFF
		# 1) it is not in @result, or
		# 2) it is in @result, and the current color is not in the color array.
		if( !defined($color{$current_loc+1}[0]) || !($color ~~ @{$color{$current_loc+1}}[1 .. $color{$current_loc+1}[0]]) ){
			# This is color-aware!
			push(@result, $current_loc+1);
			
			if( !defined($color{$current_loc+1}) ){
				$color{$current_loc+1}[0] = 1;
			}
			else{
				$color{$current_loc+1}[0]++;
			}
			$color{$current_loc+1}[ $color{$current_loc+1}[0] ] = $color;
		}
		
		$current_loc = index( $seq, $enzyme_rc, $current_loc + 1);
	}
	
	# Remove duplicated values!
	return Uniq(@result);
}

sub Print_cmap_header{
	my ($filename) = @_;
	my $OUT;
	my ($tmp_1, $tmp_2);
	
	open($OUT, ">$filename") || die ("ERROR: Can't open $filename: $!\n");
	
	for(my $i = 0; $i < @enzyme_sequences; $i++){
		$tmp_1 .= sprintf("# Nickase Recognition Site %d:	$enzyme_sequences[$i]\n", $i+1);
	}
	for(my $i = 0; $i < @enzymes; $i++){
		$tmp_2 .= sprintf("# Enzyme %d:	$enzymes[$i]\n", $i+1);
	}
	
	chomp($tmp_1);
	chomp($tmp_2);
	
	#my $total_enzymes = @enzyme_sequences;
	my $color_uniq = scalar(@color_uniq);
	
	my $str = << "EOF";
# CMAP File Version:	0.1
# Label Channels:	$color_uniq
$tmp_1
$tmp_2
# Number of Consensus Nanomaps:	N/A
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
#f int	float	int	int	int	float	float	int	int
EOF
	print $OUT($str);
	close($OUT);
}

sub Generate_cmap{
	my ($filename, $ID, $loci_ref, $length) = @_;
	my $i;
	my $siteID = 1;
	my $total_loci = 0;
	my $OUT;
	my $length_float = sprintf("%.1f", $length);
	my @sorted_loci = sort {$a <=> $b} @$loci_ref;
	
	open($OUT, ">>$filename") || die ("ERROR: Can't open $filename: $!\n");
	
	for($i = 0; $i < @sorted_loci; $i++){
		$total_loci += $color{$sorted_loci[$i]}[0];
	}
	
	for($i = 0; $i < @sorted_loci; $i++){
		my $loci_float = sprintf("%.1f", $sorted_loci[$i]);
		
		for(my $j = 1; $j <= $color{$sorted_loci[$i]}[0]; $j++){
			print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t", $color{$sorted_loci[$i]}[$j], "\t$loci_float\t1.0\t1\t1\n");
			
			$num_nicks[$color{$sorted_loci[$i]}[$j]]++;
			$num_nicks_global[$color{$sorted_loci[$i]}[$j]]++;
			
			$siteID++;
		}
	}
	
	#if(scalar(@sorted_loci)){
	print $OUT("$ID\t$length_float\t$total_loci\t$siteID\t0\t$length_float\t0.0\t1\t0\n");
	$num_cmaps++;
	#}
	
	close($OUT);
}

sub Get_and_set_enzyme_sequence{
	my ($hash_ref, $str_ref) = @_;
	
	foreach my $item (keys %$hash_ref){
		if(uc(substr($item, 0, 3)) eq uc(substr($$str_ref, 0, 3))){
			$$str_ref = $item;
			return uc($hash_ref -> {$item});
		}
	}
	
	print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
	Usage();
}

sub is_enzyme_name{
	my ($hash_ref, $str) = @_;
	
	my @array = map { uc(substr($_, 0, 3)) } keys %$hash_ref;
	
	if(uc(substr($str, 0, 3)) ~~ @array){
		return 1;
	}
	else{
		return 0;
	}
}

sub is_nt_sequence{
	my ($str) = @_;
	
	for(my $i = 0; $i < length($str); $i++){
		if("ACGTacgt" !~ substr($str, $i, 1)){
			return 0;
		}
	}
	return 1;
}

sub Uniq{
	my %seen;
	grep {!$seen{$_}++} @_;
}

sub Uniq_BAK{
	my %seen;
	my @arr_uniq;
	
	foreach (@_){
		if(!exists($seen{$_})){
			$seen{$_} = undef;
			push(@arr_uniq, $_);
		}
	}
	
	return @arr_uniq;
}

__END__

/home/users/xzhou/PERL/fa2cmap_multi_color.pl -i XXX.fa -e bspq 2 CACTTAAA 1 -v -o ./

=COMMENT
use File::Spec;
use File::Spec::Functions;
File::Spec -> rel2abs($f)

use File::Basename;
my $dir = dirname($f_in_abs);
my $filename = basename($f_in_abs);
print "dir = $dir\n";
print "filename = $filename\n";

use File::Spec::Functions;
$f_log = catfile($dir, $filename);
=cut


