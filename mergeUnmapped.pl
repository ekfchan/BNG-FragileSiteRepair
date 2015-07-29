#!/usr/bin/perl

# Script to merge fragile-site-repaired cmap with unaligned cmap from original consensus genome map set. 
#
# Usage: mergeUnmapped.pl --ref <reference.cmap> --oricmap <original.cmap> --xmap <original.xmap> --errbin <original.errbin> --fsrcmap <fsiteRepaired.cmap> [--outdir <output_folder>]

use strict; 
use warnings; 
use Cwd; use Cwd qw(abs_path);
use File::Basename;
use Getopt::Long; 
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Sys::MemInfo qw(totalmem freemem totalswap);

# << compute requirement >> 
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => my %options );
my $cpuCount = $cpu->count;
my $mem = (((&totalmem / 1024) / 1024)) / 1024; 

# << usage statement and variable initialisation >>
my %inputs = (); 
GetOptions( \%inputs, 'ref=s', 'oricmap=s', 'fsrcmap=s', 'errbin=s', 'outdir:s', 'xmap=s'); 

if ( (!exists $inputs{ref} & !exists $inputs{oricmap}) | !exists $inputs{fsrcmap} | !exists $inputs{errbin} | !exists $inputs{xmap}) {
	print "Usage: mergeUnmapped.pl --ref <reference.cmap> --oricmap <original.cmap> --xmap <original.xmap> --errbin <original.errbin> --fsrcmap <fsiteRepaired.cmap> [--outdir <output_folder>]\n"; 
	exit 0; 
}
foreach my $key ("ref","oricmap","fsrcmap","errbin","xmap") {
	if (exists $inputs{$key}) { $inputs{$key} = abs_path($inputs{$key}); }
	if ( ! -e $inputs{$key} ) { die "Input file $inputs{$key} does not exists: $!\n"; }
}

if( ! exists $inputs{outdir} ) { 
	$inputs{outdir} = cwd(); 
	$inputs{outdir} = $inputs{outdir}."/unmapped";
}
my $oriprefix = basename(abs_path($inputs{oricmap}), ".cmap");

# check output folder
if (-d $inputs{outdir}) {
	die "ERROR: Output directory $inputs{outdir} exists.\n";
} else {
	mkdir($inputs{outdir}) or die "ERROR: Cannot create output directory $inputs{outdir}: $!\n";
}
print "\n";


# << get RefAligner options >> 
my $opts; 

open XMAP, $inputs{xmap} or die $!;
while (my $line = <XMAP>) {
	chomp $line;
	if ($line =~ m=tools/RefAligner=) { 
		$opts = $line;
		last;
	} else { next; }
}
close XMAP; 

$opts =~ s/^.+RefAligner\s//;
$opts =~ s/-i\s[^\s]+\s//g; 
$opts =~ s/-o\s[^\s]+\s//g; 
$opts =~ s/-ref\s[^\s]+\s//g; 
$opts =~ s/-maxthreads\s[^\s]+\s//g; 
$opts =~ s/-maxmem\s[^\s]+\s//g; 
$opts =~ s/-output-veto-filter\s[^\s]+\s//g; 
$opts =~ s/-readparameters\s[^\s]+\s//g; 
$opts =~ s/-stderr\s+//g; 
$opts =~ s/-stdout\s+//g; 
# print $opts, "\n"; 


# << get unaligned cmap >> 
my $unmappedfile = $inputs{outdir}."/".$oriprefix."_unmapped";
# my $veto = q/-output-veto-filter '(_intervals.txt|.err|.maprate|[].map)$'/;
my $ofilter = q/-output-filter '(.[cx]map)$'/;
my $cmd = "~/tools/RefAligner -ref ".$inputs{ref}." -i ".$inputs{oricmap}." -o ".$inputs{outdir}."/".$oriprefix."  -stdout ".$ofilter." -readparameters ".$inputs{errbin}." -BestRef 1 -f -unmapped ".$unmappedfile; 
## ~~~~~~~~~~~~~
# $cmd = $cmd." -maxmem ".$mem."  -maxthreads ".$cpuCount; 
# $cmd = $cmd." ".$opts; 
# Time with xmap options
# real	0m45.296s
# user	21m39.457s
# sys	2m34.604s
## ~~~~~~~~~~~~~

## =========
# $cmd = $cmd." -maxmem ".$mem."  -maxthreads ".$cpuCount; 
# $cmd = $cmd." -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 -hashMultiMatch 100  ";
## Time with hashing
# real	0m18.771s
# user	1m50.554s
# sys	1m54.615s
## =========

## ----- 
# $cmd = $cmd." -maxmem ".$mem."  -maxthreads ".$cpuCount; 
## Time with this command 
# real	9m13.562s
# user	197m14.658s
# sys	12m9.891s
## ------ 
$cmd = $cmd." -maxmem ".$mem."  -maxthreads ".$cpuCount; 
$cmd = $cmd." ".$opts; 

print "Running command: \n\n $cmd \n\n";
system($cmd); 


# << merge unmapped with fsiteRepaired cmap >> 
my $finalcmap = $inputs{outdir}."/FSITEREPAIRED_FINAL"; 
$cmd = "~/tools/RefAligner -merge -i ".$unmappedfile.".cmap  -i ".$inputs{oricmap}." -o ".$finalcmap." -stdout -f "; 
print "Running command: \n\n $cmd \n\n";
system($cmd); 


# << calculate cmap statistics >> 
# usage: calc_cmap_stats.pl <CMAP_File>
my $script = glob("~/scripts/HybridScaffold/scripts/calc_cmap_stats.pl");
if (-e $script) {
	my $tmpcmap = $inputs{oricmap};
	print "\nOriginal consensus CMAP $tmpcmap stats:\n";
	my $cmd = "perl $script $tmpcmap"; system($cmd);
	print "\n";

	$tmpcmap = $unmappedfile.".cmap";
	print "\nUnmapped query CMAP $tmpcmap stats:\n";
	$cmd = "perl $script $tmpcmap"; system($cmd);
	
	$tmpcmap = $inputs{fsrcmap};
	print "\nInput fsiteRepaired CMAP $tmpcmap stats:\n";
	$cmd = "perl $script $tmpcmap"; system($cmd);

	$tmpcmap = $finalcmap.".cmap";
	print "\nFinal merged fragile-site-repaired CMAP $tmpcmap stats:\n";
	$cmd = "perl $script $tmpcmap"; system($cmd);
}
else { print "perl script $script not found\n"; }


