#usage: perl split_xmap_standalone.pl [xmap] [_q.cmap] [_r.cmap] [_contig prefix] [output folder]
#example: perl split_xmap_standalone.pl C24_BspQI_to_EXP_REFINEFINAL1.xmap C24_BspQI_to_EXP_REFINEFINAL1_r.cmap C24_BspQI_to_EXP_REFINEFINAL1_q.cmap C24_BspQI_to_EXP_REFINEFINAL1_contig ./contigs

use strict;
use warnings;
use Cwd qw(getcwd cwd abs_path realpath);
use Config::Simple;
use IPC::System::Simple qw(system systemx capture capturex);
use File::Path qw(mkpath rmtree);
use File::Slurp;
use File::Copy qw(copy move);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use XML::Simple;
use DateTime;
use DateTime::Format::Human::Duration;
use File::Basename;
use threads;
use IO::Select;
use IPC::Open3;
use File::Spec;

#my $x_map="C24_BspQI_to_EXP_REFINEFINAL1.xmap";
#my $r_map="C24_BspQI_to_EXP_REFINEFINAL1_r.cmap";
#my $q_map="C24_BspQI_to_EXP_REFINEFINAL1_q.cmap";
#my $string="C24_BspQI_to_EXP_REFINEFINAL1_contig";
#my $path="./contigs_lion4";

my $x_map=$ARGV[0];
my $r_map=$ARGV[2];
my $q_map=$ARGV[1];
my $string=$ARGV[3];
my $path=$ARGV[4];

split_xmap($x_map,$r_map,$q_map,$string,$path);

sub split_xmap {
	my ($xmap, $rcmap, $qcmap, $contig_prefix, $molPath) = @_;
	 	
	my @xmap_in = read_file($xmap);	
	my @rcmap_in = read_file($rcmap);
	my @qcmap_in = read_file($qcmap);
	
	my %rcmap;
	my $rcmap_header="";
	my %qcmap;
	my $qcmap_header="";
	
	my @xmapheader;
	
	print "\nSplitting $xmap...\n";
	
	for (my $i=0; $i <= 4; $i++) {
		push (@xmapheader, $xmap_in[$i]); }
		
	eval { mkpath($molPath) };
	if ($@) {
		print "ERROR: Couldn't create $molPath: $@"; }
		
	
	foreach my $rline (@rcmap_in) {
	if ($rline =~ /^#/ ) {
		$rcmap_header = $rcmap_header . $rline; }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$Position,$StdDev,$Coverage,$Occurrence,$GmeanSNR,$lnSNRsd) = split("\t",$rline);
		if ($CMapId) {
			push @{ $rcmap{$CMapId} }, $rline; }}}
			
	foreach my $qline (@qcmap_in) {
	if ($qline =~ /^#/ ) {
		$qcmap_header = $qcmap_header . $qline;  }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$PositionStdDev,$Coverage,$Occurrence) = split("\t",$qline);
		if ($CMapId) {
			push @{ $qcmap{$CMapId} }, $qline; }}}
	
	
	my $prevContigID = 0;
	foreach my $xline (@xmap_in) {
		if ($xline !~ /^#/ ) {			
			my ($XmapEntryID, $QryContigID, $RefContigID, $QryStartPos, $QryEndPos, $RefStartPos, $RefEndPos, $Orientation, $Confidence, $HitEnum) = split("\t",$xline);  
			if ($RefContigID != $prevContigID) {
				
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				my $rcmapName = $contig_prefix.$RefContigID."_r.cmap";
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				
				open (XMAP, ">>$molPath/$xmapName") or die "ERROR: $!\n";
				open (RCMAP, ">>$molPath/$rcmapName") or die "ERROR: $!\n";
				open (QCMAP, ">>$molPath/$qcmapName") or die "ERROR: $!\n";
				
				foreach my $s (@xmapheader) {
					$s =~ s/^\s+//;
					print XMAP $s; }
				print XMAP "# Reference Maps From:	".$rcmapName."\n";
				print XMAP "# Query Maps From:	".$qcmapName."\n";
				print XMAP "#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum\n";
				print XMAP "#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string\n"; 
				print XMAP $xline;
				

				
				# RCMAP output
				print RCMAP $rcmap_header; 
				if($rcmap{$RefContigID}){
					print RCMAP join("", @{ $rcmap{$RefContigID} }); }
				
				#QCMAP output
				print QCMAP $qcmap_header;
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }
				
				$prevContigID = $RefContigID;
					
				close QCMAP;
				close RCMAP;
				close XMAP;
				}				
				
			elsif ($RefContigID == $prevContigID) {
				
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				
				open (XMAP, ">>$molPath/$xmapName") or die "ERROR: $!\n";
				open (QCMAP, ">>$molPath/$qcmapName") or die "ERROR: $!\n";				

				print XMAP $xline; 
				
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }

				close QCMAP;
				close XMAP;				
				
				} } 
				
				

				
				
				}
	
	print "Splitting $xmap complete.\n";
	exit 0;
	
}
