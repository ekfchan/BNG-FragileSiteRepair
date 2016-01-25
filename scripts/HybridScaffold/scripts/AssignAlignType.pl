# $Id: AssignAlignType.pl 4268 2015-11-12 01:40:05Z apang $

#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
# This script sould sit one level above the additional Perl modules directory.
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path2 = abs_path($script_path . "/perl5");
	unshift @INC, $module_path2;
	my $lib4;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.016000 && $] <= 5.017000) {
		$module_path2 = $module_path2."/5.10.1";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.018000 && $] <= 5.019000) {
		$module_path2 = $module_path2."/5.18.2";
		$lib4 = $module_path2."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X and 5.16.X and 5.18.X\n"; 
		exit; }
	unshift @INC, $module_path2;
	unshift @INC, $lib4;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

print "\nInfo: Running the command: $0 @ARGV\n";

use BNG::Utility;
my $xmap_in = $ARGV[0];
my $ngs_cmap_in = $ARGV[1];
my $bng_cmap_in = $ARGV[2];
my $sticky_xmap_fn = $ARGV[3];
my $filtered_ngs_fn = $ARGV[4];
my $filtered_bng_fn = $ARGV[5];
my $T_cutoff = $ARGV[6];
my $max_overhang = $ARGV[7];
my $orig_ngs_cmap_in = $ARGV[8];
my $orig_bng_cmap_in = $ARGV[9];
my $break_point_file = $ARGV[10];

my $xmap = readXMap($xmap_in);
my ($ngs_cmap, undef, undef) = readCMap($ngs_cmap_in);
my ($bn_cmap, undef, undef) = readCMap($bng_cmap_in);
my ($orig_ngs_cmap, undef, undef) = readCMap($orig_ngs_cmap_in);
my ($orig_bn_cmap, undef, undef) = readCMap($orig_bng_cmap_in);
my @xmap_data_name = @{$xmap->{dataName}};

## we shift 20 bp to find out only labels extending outside alignment
my @refLeftLabCnt = ();
my @refRightLabCnt = ();
my @qryLeftLabCnt = ();
my @qryRightLabCnt = ();

for (my $i=0; $i < $xmap->{totalHits}; $i++) {
	my ($numL, $numR) = countOverangLables( $ngs_cmap, $xmap->{hits}->{RefContigID}->[$i],
					($xmap->{hits}->{RefStartPos}->[$i]) - 20,
					($xmap->{hits}->{RefEndPos}->[$i]) + 20       );
	push(@refLeftLabCnt, {numLabel => $numL, bkpt => $xmap->{hits}->{RefStartPos}->[$i]});
	push(@refRightLabCnt, {numLabel => $numR, bkpt => $xmap->{hits}->{RefEndPos}->[$i]});

	if (($xmap->{hits}->{Orientation}->[$i]) eq "+") {
		($numL, $numR) = countOverangLables( $bn_cmap, $xmap->{hits}->{QryContigID}->[$i],
					($xmap->{hits}->{QryStartPos}->[$i]) - 20,
					($xmap->{hits}->{QryEndPos}->[$i]) + 20      );
	} else {
		($numR, $numL) = countOverangLables( $bn_cmap, $xmap->{hits}->{QryContigID}->[$i],
					($xmap->{hits}->{QryEndPos}->[$i]) - 20,
					($xmap->{hits}->{QryStartPos}->[$i]) + 20      );
	}
	push(@qryLeftLabCnt, {numLabel => $numL, bkpt => $xmap->{hits}->{QryStartPos}->[$i]});
	push(@qryRightLabCnt, {numLabel => $numR, bkpt => $xmap->{hits}->{QryEndPos}->[$i]});
} # for i

# create sticky xmap:
my $stickyXMap={};
$stickyXMap->{headers} = $xmap->{headers};
$stickyXMap->{dataName} = $xmap->{dataName};
$stickyXMap->{dataType} = $xmap->{dataType};
$stickyXMap->{MAPFileVersion} = $xmap->{MAPFileVersion};
my $numc = @xmap_data_name;
for (my $i=0; $i<$numc; $i++) {
	$stickyXMap->{hits}->{$xmap_data_name[$i]} = [];
} # for i
my $tg = 0;
my $stickyNGSContigs={};
my $stickyBNContigs={};
my @breakPointInfo = ();
for (my $i=0; $i < $xmap->{totalHits}; $i++) {
	if ((($xmap->{hits}->{Confidence}->[$i]) >= $T_cutoff) && (
			($refLeftLabCnt[$i]{numLabel} > $max_overhang &&
			$qryLeftLabCnt[$i]{numLabel} > $max_overhang) ||
			($refRightLabCnt[$i]{numLabel} > $max_overhang &&
			$qryRightLabCnt[$i]{numLabel} > $max_overhang) ) ) {
		# significant overhang
		$tg++;
		
		for (my $m=0; $m<$numc; $m++) {
			my $a = $stickyXMap->{hits}->{$xmap_data_name[$m]};
			my $b = $xmap->{hits}->{$xmap_data_name[$m]};
			my $c = $b->[$i];
			push(@$a, $c);
			$stickyXMap->{hits}->{$xmap_data_name[$m]}=$a;
		} # for m
		$stickyNGSContigs->{ $xmap->{hits}->{RefContigID}->[$i] } = 1;
		$stickyBNContigs->{ $xmap->{hits}->{QryContigID}->[$i] } = 1;

		my ($leftRefBkpt, $leftQryBkpt) = ($refLeftLabCnt[$i]{numLabel} > $max_overhang && $qryLeftLabCnt[$i]{numLabel} > $max_overhang) ? ($refLeftLabCnt[$i]{bkpt}, $qryLeftLabCnt[$i]{bkpt}) : (-1, -1);
		my ($rightRefBkpt, $rightQryBkpt) = ($refRightLabCnt[$i]{numLabel} > $max_overhang && $qryRightLabCnt[$i]{numLabel} > $max_overhang) ? ($refRightLabCnt[$i]{bkpt}, $qryRightLabCnt[$i]{bkpt}) : (-1, -1);
		push(@breakPointInfo, {xmapId => $xmap->{hits}->{XmapEntryID}->[$i], refId => $xmap->{hits}->{RefContigID}->[$i], qryId => $xmap->{hits}->{QryContigID}->[$i], 
					leftRefBkpt => $leftRefBkpt, leftQryBkpt => $leftQryBkpt, 
					rightRefBkpt => $rightRefBkpt, rightQryBkpt => $rightQryBkpt, 
					alignOrientation => $xmap->{hits}->{Orientation}->[$i]});
	} # if xmap->{hits}
} # for i


$stickyXMap->{totalHits} = $tg;
#writeXMapFile($sticky_xmap_fn, $stickyXMap, "noheader" );

my $conflictingNgsCmap = getSubsetCMap($orig_ngs_cmap, $stickyNGSContigs, 'include');
#writeCMapFile($conflictingNgsCmap, $conflicting_ngs_fn, 1);
my $conflictingBngCmap = getSubsetCMap($orig_bn_cmap, $stickyBNContigs, 'include');
#writeCMapFile($conflictingBngCmap, $conflicting_bng_fn, 1);
writeXMapFile($sticky_xmap_fn, $stickyXMap, 1, $ngs_cmap_in, $bng_cmap_in);     # header is 1 (or 0), and the _r cmap is conflicting_ngs_fn, and the _q cmap is conflicting_bng_fn

my $filtered_ngs_cmap = getSubsetCMap($orig_ngs_cmap, $stickyNGSContigs, 'exclude');
writeCMapFile($filtered_ngs_cmap, $filtered_ngs_fn, 1);
my $filtered_bionano_cmap = getSubsetCMap($orig_bn_cmap, $stickyBNContigs, 'exclude');
writeCMapFile($filtered_bionano_cmap, $filtered_bng_fn, 1);

printBkptFile($break_point_file, \@breakPointInfo);
print "Info: $0 finished successfully.\n\n";
exit 0;

sub printBkptFile	{
	my ($file, $breakPointInfoRef) = @_;
	my $shiftAmount = 30;
	open(OUT, ">$file") or die "ERROR: AssignAlignType printBkptFile: cannot write to $file: $!\n";
	print OUT "xMapId\trefQry\trefId\tleftRefBkpt\trightRefBkpt\talignmentOrientation\trefQry\tqryId\tleftQryBkpt\trightQryBkpt\talignmentOrientation\n";
	for (my $i = 0; $i < scalar(@{$breakPointInfoRef}); $i += 1)	{
		# print in one row
		# if the breakpoint is the left breakpoint, then lower $shftAmount bp of the breakpoint coordinates
		# if the breakpoint is the right breakpoint, then raise $shiftAmount bp of the breakpoint coordinates (reason: avoid cutting at  the aligned label in subsequence cut conflict step)
		# $shiftAmount bp is okay since there must be at least $max_overhang unaligned labels beyond the aligned region
		my $bRef = $breakPointInfoRef->[$i];
		my ($toPrintLeftRefBkpt, $toPrintRightRefBkpt) = (($bRef->{leftRefBkpt}==-1) ? ($bRef->{leftRefBkpt}) : ($bRef->{leftRefBkpt} - $shiftAmount), ($bRef->{rightRefBkpt}==-1) ? ($bRef->{rightRefBkpt}) : ($bRef->{rightRefBkpt} + $shiftAmount));
		my ($toPrintLeftQryBkpt, $toPrintRightQryBkpt) = (($bRef->{leftQryBkpt}==-1) ? ($bRef->{leftQryBkpt}) : ($bRef->{leftQryBkpt} - $shiftAmount), ($bRef->{rightQryBkpt}==-1) ? ($bRef->{rightQryBkpt}) : ($bRef->{rightQryBkpt} + $shiftAmount));

		print OUT join("\t", $bRef->{xmapId},	"ref", $bRef->{refId}, $toPrintLeftRefBkpt, $toPrintRightRefBkpt, $bRef->{alignOrientation})."\t";
		print OUT join("\t", 			"qry", $bRef->{qryId}, $toPrintLeftQryBkpt, $toPrintRightQryBkpt, $bRef->{alignOrientation})."\n";
	} # for i
	close OUT;
} # printBkptFile 


