#!/usr/bin/perl

#Extract sequences from stitchPositions BED file with sequences and output FASTA
#Usage: perl extractFASTA_fromBED.pl --bed <input BED> [--output <output FASTA filename>]

use strict;
use warnings;

use Cwd; use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long; 

my %inputs = (); 
my $cwd;
my $basename;
my $outputfile;

GetOptions( \%inputs, 'bed=s', 'output=s');
if ( (!exists $inputs{bed}) ) {
	print "Usage: perl extractFASTA_fromBED.pl --bed <input BED> [--output <output FASTA>]\n";
	exit 0;
}
else {
	$cwd = abs_path($inputs{bed});
	$cwd = dirname($cwd);
	$basename = basename($inputs{bed});
	$basename =~ s/\.[^.]+$//;
}

if( defined $inputs{output} and !-e $inputs{output}) { 
	$outputfile = abs_path($inputs{output}); 
} 
elsif ( defined $inputs{output} and -e $inputs{output}) {
	print "Output file $inputs{output} already exists! Using default output filename...\n";
	$outputfile = $cwd."/".$basename.".fasta";
}
else {
	$outputfile = $cwd."/".$basename.".fasta";  #if output file not given use $basename as prefix 
}
if( -e $outputfile ) { 
	die "ERROR: Output file ", $outputfile, " already exists. Please specify alternate filename. $!\n"; 
}

foreach my $key ("bed") {
	if (exists $inputs{$key}) { $inputs{$key} = abs_path($inputs{$key}); }
	if ( ! -e $inputs{$key} ) { die "Input file $inputs{$key} does not exists: $!\n"; }
}

print "\n";
print "\tInput BED file: $inputs{bed}\n";
print "\tOutput FASTA file: $outputfile\n";
print "\n";

open OUT, ">$outputfile" or die "ERROR: Could not open output file $outputfile: $!\n";

open BED, "<$inputs{bed}" or die "ERROR: Could not open input BED: $!\n";
while (<BED>) {
	my $line = $_;
	chomp($line);
	if ($line =~ /^#/) {
		next;
	}
	else {
		my @s = split("\t",$line);
		my $fastaHeader = ">$s[0]_$s[1]_$s[2]_$s[3]";
		print OUT "$fastaHeader\n$s[9]\n";
	}
}
		
