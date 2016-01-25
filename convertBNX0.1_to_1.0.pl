#!/usr/bin/perl

# Designed to work with 1 color BNX 0.1

use strict;
use warnings;

open FILE, "$ARGV[0]" or die $!;

my $outname = "$ARGV[0]_temp";

#first convert BNX 0.1 to CMAP
#RefAligner -i all0918_stretch.bnx -o all0918_stretch -merge -minsites 0
my $cmd = "~/tools/RefAligner -i $ARGV[0] -o $outname -merge -minsites 0 -f -stdout -stderr";
system($cmd);

#second convert CMAP 0.1 to BNX 1.0
#RefAligner -i all0918_stretch.cmap -o all0918_stretch_122314 -bnx -bnxversion 1 -merge -minsites 0
$cmd = "~/tools/RefAligner -i $outname.cmap -o $outname -bnx -bnxversion 1 -merge -minsites 0 -f -stdout -stderr";
system($cmd);

#third clean up intermediary CMAP
$cmd = "rm $outname.cmap";
system($cmd);

close FILE;

#fourth fix new BNX 1.0 header and add in fake QX lines
open FILE, "$outname.bnx" or die $!;
my $outname2 = "$ARGV[0]_bnx1.0.bnx";
open OUT, ">$outname2" or die $!;
chomp(my @bnxin = <FILE>);
close FILE;
my @bnxout;

#fix header
for (my $i=0; $i<8; $i++) {
	print OUT "$bnxin[$i]\n";
}
print OUT "# Min Molecule Length (Kb):	0\n";
print OUT "# Label SNR Filter Type:	Static\n";
print OUT "# Min Label SNR:	0.000\n";
print OUT "# Software Version:	\n";
print OUT "$bnxin[8]\n";
#rh SourceFolder	InstrumentSerial	Time	NanoChannelPixelsPerScan	StretchFactor	BasesPerPixel	NumberofScans	ChipId	FlowCell	SNRFilterType	MinMoleculeLength	MinLabelSNR	RunId
print OUT "# Run Data	/DummyData_2014-07-17_18_38/Detect Molecules	ALPHAUNIT09	7/17/2014 6:38:51 PM	68819821	0.85	500	1	20249,11843,07/17/2014,840014289	1	Static	0	0	1\n";
for (my $i=10; $i<17; $i++) {
	print OUT "$bnxin[$i]\n";
}
print OUT "# Quality Score QX11: SNR for channel 1\n";
print OUT "# Quality Score QX12: Ave Intensity for channel 1\n";

#fix body
for (my $i=17; $i<scalar(@bnxin); $i+=2) {
	my $backbone = $bnxin[$i];
	my $labels = $bnxin[$i+1];
	print OUT "$backbone\n";
	print OUT "$labels\n";
	my @lab = split("\t",$labels);
	my $sites = scalar(@lab) - 2;
	my $qx11 = "QX11\t";
	my $qx12 = "QX12\t";
	for (my $i=0; $i<($sites-1); $i++) {
		$qx11 = $qx11."10\t";
		$qx12 = $qx12."10\t";
	}
	$qx11 = $qx11."10";
	$qx12 = $qx12."10";
	print OUT "$qx11\n";
	print OUT "$qx12\n";
}

#clean up
$cmd = "rm $outname.bnx";
system($cmd);
	
		
		
		
		
		
		
		

