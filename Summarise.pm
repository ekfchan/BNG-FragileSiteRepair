#!/usr/bin/perl

package Summarise;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

BEGIN { 
	require Exporter; 
	our $VERSION     = 1.00;
	our @ISA         = qw(Exporter);
	our @EXPORT      = ();
	our @EXPORT_OK   = qw(GetCmapIds CalcCmapStats GetFsiteStats);
}

sub GetCmapIds {
	my $cmapfile = shift; 
	my @cmapIds; 
	open (CMAP, "<", $cmapfile) or die "ERROR: Could not open $cmapfile: $!\n";
	while (my $line = <CMAP>) {
		chomp($line);
		if( $line !~ /^#/ ) {
			my @s = split(/\t/,$line); 
			my $id = $s[0];
			$id =~ s/^\s+|\s+$//g;
			push @cmapIds, $id; 
		}
	}
	close CMAP;
	@cmapIds = unique(@cmapIds);

	# return( \@cmapIds ); 
	return( @cmapIds ); 
}


sub CalcCmapStats {
	my $cmapfile = shift;
	# usage: calc_cmap_stats.pl <CMAP_File>
	# my $script = glob("$ENV{'HOME'}/scripts/HybridScaffold/scripts/calc_cmap_stats.pl");
	my $script = $ENV{'HOME'}."/scripts/HybridScaffold/scripts/calc_cmap_stats.pl";
	if (-e $script) {
		my $cmd = "perl $script $cmapfile";
		print "\n---------------\n"; 
		system($cmd);
		# if ($? == -1) { print "failed to execute command: $!\n\t$cmd"; }
		# elsif ($? & 127) {
		# 	printf "child died with signal %d, %s coredump\n",
		# 	($? & 127),  ($? & 128) ? 'with' : 'without';
		# }
		# else { printf "child exited with value %d\n", $? >> 8; }
		print "\n---------------\n"; 
	}
	else {
		print "perl script calc_cmap_stats.pl not found at $script\n"; }
}


sub GetFsiteStats {
	my $bedfile = shift @_; 

	my %mapCount;
	my %typeCount; 
	my $fragileSites = 0;
	open (FSITES, "<", $bedfile) or die "ERROR: Could not read $bedfile: $!\n";
	while (my $line = <FSITES>) {
		chomp($line);
		if ( $line !~ /^#/) {
			$fragileSites++; 
			#CMapId	Start	End	Type	Score	Strand	ThickStart	ThickEnd	ItemRgba	Sequence
			my @s = split(/\t/,$line);
			$mapCount{$s[0]}++;
			if ($s[3] =~ m/Type/i) { $typeCount{$s[3]}++; }
		}

	}
	close FSITES; 
	print "\nFound a total of $fragileSites in $bedfile :\n";
	print "\tBy Map ID:\n"; 
	foreach my $map (sort keys %mapCount) { print "\t\tMap ID $map: $mapCount{$map}\n"; }
	print "\tBy fragile site Type:\n"; 
	foreach my $type (sort keys %typeCount) { print "\t\t$type: $typeCount{$type}\n"; }
	print "\n"; 
}


sub unique {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

1; 

