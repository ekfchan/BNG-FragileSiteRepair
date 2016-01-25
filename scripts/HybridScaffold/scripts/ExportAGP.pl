# $Id: ExportAGP.pl 4317 2015-12-02 23:21:46Z apang $
#!/usr/bin/perl -w

########################################################################
# File: ExportAGP.pl                                                   #
# Date: 9/23/2015                                                      #
# Purpose: Convert xmap output from HybridScaffold to agp format       #
#                                                                      #
# Author: Jian Wang                                                    #
# Email : jwang@bionano.com                                            #
# Affiliation: Research Department, BioNano Genomics Inc.              #
## Usage:                                                              #
#  
########################################################################

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
	my $lib2; my $lib4;
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
use Getopt::Std;
use Tie::File;

#global varaibles

our $DEBUG = 0;
our $TYPE_GAP=0;
our $TYPE_CONTIG=1;
our $MAX_MEM = 2000000;
our $path_sep="/";


#enzyme cut site sequences
our %enzyme = (
	"BSPQI" => "GCTCTTC",
	"BBVCI" => "CCTCAGC",
	"BSMI"  => "GAATGC",
	"BSRDI" => "GCAATG",
	"BSECI" => "ATCGAT",
	"BSSSI" => "CACGAG"
	# You can add more enzymes here ...
	
);

sub initCli;
sub sortXMap;
sub updateNGSFiles;
sub printUnUsedNGS;
sub cleanUp;

###################Main entry point of the export script########################

our $INPUTS = initCli() || die "";
our $sorted_xmap_file = sortXMap($INPUTS->{xmap_file}, $INPUTS->{out_dir});
our $xmap = readXMap($sorted_xmap_file);

if($INPUTS->{cut_coord_file}){
	print "Detect cut-coord file, updating fast and fasta name map file\n";
	my ($fasta_file, $ngs_namemap_file) = updateNGSFiles($INPUTS->{fasta_file}, 
														 $INPUTS->{ngs_namemap_file}, 
										                 $INPUTS->{cut_coord_file}, $INPUTS->{out_dir}); 
	$INPUTS->{fasta_file} = $fasta_file;
	$INPUTS->{ngs_namemap_file} = $ngs_namemap_file;										                 
}	

print("Processing fasta related files: ".$INPUTS->{fasta_file}."\t".$INPUTS->{ngs_namemap_file});
our $ngsMap = getNGSMap($INPUTS->{ngs_namemap_file});
our $seqMap = readFasta($INPUTS->{fasta_file}, $INPUTS->{out_dir});


our ($hybridCmap, $numcontig, $contigLength) = readCMap($INPUTS->{hybrid_cmap_file});



processAlign($xmap, $INPUTS->{gap_out_file}, 
					$INPUTS->{agp_out_file}, 
					$INPUTS->{begin_end_file},
					$INPUTS->{paddingGapLen},
					$ngsMap);
printFasta($seqMap, $xmap, $INPUTS->{fasta_out_file}, $hybridCmap, $INPUTS->{cut_site}, $INPUTS->{paddingGap});
printUnUsedNGS($INPUTS->{agp_out_file}, $ngsMap, $INPUTS->{begin_end_file}, $INPUTS->{Unused_fasta_out}, $seqMap);
cleanUp($INPUTS);

###################Utility functions ############################################
#get the file name from a path string
sub getFileName{
	my $path = $_[0];
	my $sep = $_[1];
	if(not defined $sep){
		$sep = '/';
	}
	return substr($path, rindex($path, $sep)+1);
}

#strip file extension from file path
sub stripExtension{
	my $file = $_[0];
	my $sep = $_[1];
	if(not defined $sep){
		$sep = '.';
	}
	my $endInd = rindex($file, $sep);
	if($endInd < 0){
		$endInd = length($file);
	}
	return substr($file, 0, $endInd);
}

#obtain the basename from a path
sub getFileBasename{
	my $basename = `basename $_[0]`;
	chomp $basename;
	return $basename;
}


sub round{
	return int($_[0] + 0.5);
}

#creates array by repeating some values n times
sub rep{
	my $unit = $_[0];
	my $len = $_[1];
	my @arry=();
	for(my $i = 0; $i < $len; $i++){
		push(@arry, $unit);
	}
	return(\@arry);
}

###############End of Utility functions ###########################################################



#############Function for pre-processing inputs for the exporter ###################################

#Print help usage message
sub Usage{
	my $usage="Usage:                                                       
     ExportAGP.pl <-h> <-i xmap_file> <-o output_file>                   
       -h    : This help message                                        
       -i    : Input xmap file
       -m    : fasta to cmap name mapping file (i.e. generated from fa2cmap.pl)                                          
       -o    : Output directory   
       -s    : Input fasta file from NGS data
       -c    : CMap file for the hybrid contigs
       -e    : Names of the enzymes or the sequence pattern of the cut-site
       -r    : The output file when hybrid-scaffold was ran with -M option (i.e cut ngs contig/bionano map to resolve conflict)
       -p    : Length of gap to be inserted between overlapping NGF contigs in a super-scaffold contig\n";
    print $usage;
}



#sort xmap by refID and increasing refStartPos and decreasing refEndPos
#this simpifly checking if one contigs is completely embedded in another contigs
sub sortXMap{
	my $xmapFile = $_[0];
	my $out_dir = $_[1];
	my $headerCnt = `grep -c '^#' $xmapFile`;
	#print "header count: $headerCnt\n";
	chomp $headerCnt;
	$headerCnt = $headerCnt; 
	my $sorted_xmap = "$out_dir/".getFileBasename($xmapFile)."_sorted.xmap";
	print $sorted_xmap."\n";
	`head -n $headerCnt $xmapFile > $sorted_xmap`;  #export header to sorted xmap 
     $headerCnt = $headerCnt+1; #set to appropriate param for tail
	`tail -n+$headerCnt $xmapFile | sort -k3,3n -k6,6n -k7,7nr >> $sorted_xmap`;
	return($sorted_xmap);
}


#Parse command line options and extract input parameters
#Parameters are stored in the Hash %INPUTS
sub initCli(){
	my $opt_string = "hi:o:s:c:e:m:r:p:";
	my %options=();
	my %INPUTS =();
	if(!getopts("$opt_string", \%options)){
		print("ERROR: Invalid parameter(s)! Must have -i and -o options.\n");
		Usage();
		exit(1);
	}

	if(!defined $options{"i"} || !defined $options{"o"} || !defined $options{"m"}){
		print("ERROR: Missing parameter(s)! Try -h for more information.\n");
		Usage();
		exit(1);
	}
	
	$INPUTS{xmap_file} = $options{"i"};
	$INPUTS{out_dir} = $options{"o"};
	$INPUTS{agp_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".agp";
	$INPUTS{begin_end_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file})))."_trimHeadTailGap.coord";
	$INPUTS{fasta_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".fasta";
	$INPUTS{Unused_fasta_out} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file})))."_NOT_SCAFFOLDED.fasta"; #fasta file for unused ngs contigs
	$INPUTS{gap_out_file} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".gap";
	#$auxiliary_file = stripExtension($out_file)."_trimmedcontig.txt";
	print "output agp file: $INPUTS{agp_out_file}\n";
	
	$INPUTS{ngs_namemap_file} = $options{"m"};
	
	if(defined $options{"r"}){
		$INPUTS{cut_coord_file} = $options{"r"};
	}
	
	if(defined $options{"p"}){
		$INPUTS{paddingGapLen} = $options{"p"};
		if($INPUTS{paddingGapLen} < 1){
			warn("Padding gap length must be at least 1, it is now $INPUTS{paddingGapLen}, using default value 500 instead");
			$INPUTS{paddingGapLen} = 500;
		}
	}else{
		$INPUTS{paddingGapLen} = 500;
	}
	print "Padding gap length between overlapping ngs contig: $INPUTS{paddingGapLen}\n";
	my $gap="";
	for(my $c=0; $c < $INPUTS{paddingGapLen}; $c++){
		$gap = $gap."N";
	}
	$INPUTS{paddingGap}=$gap;
	
	
	if(defined $options{"s"} || defined $options{"c"} || defined $options{"e"}){
		if(defined $options{"s"} && defined $options{"c"} && defined $options{"e"}){		
			$INPUTS{fasta_file} = $options{"s"};
			$INPUTS{hybrid_cmap_file} = $options{"c"};
			$INPUTS{cut_site} = $options{"e"};
			$INPUTS{fasta_out} = $INPUTS{out_dir}.$path_sep.(stripExtension(getFileName($INPUTS{xmap_file}))).".fasta";
			#if we can find the name of the enzyme in table we used the known cut-site sequence
			#otherwise we assume the cut_site is input instead
			#cut_sites is in multi-color format: enzyme1 label1 enzyme2 label2
			my @tokens = split(/\s+/, $INPUTS{cut_site});
			#print "argument values for -e: $cut_site\n";
			my $all_cap_name = uc $tokens[0];  #only use first enzyme for now 
			if($enzyme{$all_cap_name}){
				print "Found enzyme $INPUTS{cut_site} with cut site: ".($enzyme{$all_cap_name})."\n";
				$INPUTS{cut_site} = $enzyme{$all_cap_name};
			}else{
				print "Cannot find enzyme $INPUTS{cut_site}, using input as cut-site sequence\n";
			}
		}else{
			print("ERROR: Missing -s or -c or -e options! Try -h for more information.\n");
			print("-s: ".$options{"s"}."\n");
			print("-c: ".$options{"c"}."\n");
			print("-e: ".$options{"e"}."\n");
			Usage();
			return 0;
		}
	}
	if(defined $options{"x"}){
	 
	}
	return \%INPUTS;
}

#parse the ngs mapping file to create name mapping between contig id in cmap and xmap and its ngs name and other stats
sub getNGSMap{
	my $name_map_file = $_[0];
	my $cut_coord_file = $_[1];
	my %name_map=();
	
	open(NAMEMAP, $name_map_file) || dieLog "ERROR: cannot open sequene name map file $name_map_file";
	while(my $line = <NAMEMAP>){
		if($line =~ /^#/ || $line !~ /^\d/){
			next;
		}
		chomp $line;
		my @tokens = split(/\t/, $line);
		#we store the ngs stat as 1) name; 2) length; 3) a flag whether it is used scaffolding; 4) the begin coord 
		#of this contig in case it is split by the pipline into multiple smaller contigs
		$name_map{$tokens[0]} = [$tokens[1], $tokens[2], 0, 1];
	}
	if(defined $cut_coord_file && length($cut_coord_file) > 0){
		print "cut_coord_file is: $cut_coord_file\n";
		open(CUT_COORD, $cut_coord_file) || dieLog "ERROR: cannot open cut coordinate files $cut_coord_file\n";
		while(my $line = <CUT_COORD>){
			if($line =~ /^oldId/){ #skipping header
				next;
			}
			chomp $line;
			my @tokens = split(/\t/, $line);
			if(exists $name_map{$tokens[0]}){
				my @ngsstat = @{$name_map{$tokens[0]}};
				my $newlength = $tokens[2] - $tokens[1] + 1; 
				$name_map{$tokens[3]} = [$ngsstat[0], $newlength, 0, $tokens[1]+1];
			}else{
				warn "Warning map with Id: $tokens[0] is not recognized to be a valide ID, check $cut_coord_file or $name_map_file\n";
			}
		}
	}
	if($DEBUG){
		print "Total ngs named parsed: ".(keys %name_map)."\n";
		print "printing out ngs map file\n";
		my $count=1;
		foreach my $key(keys %name_map){
			print $key."\t";
			print join "\t", @{$name_map{$key}};
			print "\n";
			$count = $count + 1;
		}
	}	
	return \%name_map;
}

#when hybrid-scaffold cut ngs contigs to resolve conflicts
#we need to generate new fasta files and ngs name map files for 
#the cutted ngs contigs
sub updateNGSFiles{
	my $fasta_file = $_[0];
	my $ngs_namemap_file = $_[1];
	my $cut_coord_file = $_[2];
	my $out_dir = $_[3];
	
	my $cutted_fasta_file = "$out_dir/".getFileBasename($fasta_file).".cutted.fasta";
	my $cutted_ngsNameMap_file = "$out_dir/".getFileBasename($ngs_namemap_file).".cutted.txt";
	
	my  %ngs_map = %{getNGSMap($ngs_namemap_file, $cut_coord_file)};
	my $fasta_accessor = readFasta($fasta_file, $out_dir);
	
	

	open(my $cutted_fasta, ">".$cutted_fasta_file) || dieLog "ERROR: cannot open file for writing: $cutted_fasta_file";
	open(my $cutted_ngs_nameMap, ">".$cutted_ngsNameMap_file) || dieLog "ERROR: cannot open file for writing: $cutted_ngsNameMap_file";

	print $cutted_ngs_nameMap "CompntId\tCompntName\tCompntLength\n";
	foreach my $key(keys %ngs_map){
		my @ngscontig = @{$ngs_map{$key}};
		my $fasta_seq = $fasta_accessor->($ngscontig[0]);
		if($ngscontig[1] == length($fasta_seq)){
			print $cutted_fasta ">$ngscontig[0]\n";
			print $cutted_fasta $fasta_seq."\n";
			print $cutted_ngs_nameMap $key."\t".$ngscontig[0]."\t".$ngscontig[1]."\n";
			#print 			
		}else{
			my $new_ngs_name = $ngscontig[0]."_subseq_".$ngscontig[3].":".($ngscontig[3] + $ngscontig[1]-1);
			print $cutted_fasta ">$new_ngs_name\n";
			print $cutted_fasta substr($fasta_seq, $ngscontig[3] - 1, $ngscontig[1])."\n";
			print $cutted_ngs_nameMap $key."\t".$new_ngs_name."\t".$ngscontig[1]."\n";
		}	
	}
	print "DONE updating NGS files\n";
	close($cutted_fasta);
	close($cutted_ngs_nameMap);
	return ($cutted_fasta_file, $cutted_ngsNameMap_file);
}



#Read in a fasta file, return a function pointer which can be used 
#to access the fasta sequence by NGS name
sub readFasta{
	#open FASTA, $_[0] || die "cannot open fasta file";
	my $fasta_file = $_[0];
	my $out_dir = $_[1];
	
	my ($fastaMap, $tmp_file) = processFastaFile($fasta_file, $out_dir);
	tie my @fasta_arry, 'Tie::File',  $tmp_file,  or dieLog "ERROR: cannot process fasta file";
	my %fasta_map = %{$fastaMap};
	
	my $accesser = sub{
		my $ID = $_[0];
		if(exists $fasta_map{$ID}){
			my $ind = $fasta_map{$ID};
			return $fasta_arry[$ind];
		}else{
			print "Warning: cannot find sequence for ID: $ID\n";
			return "";
		}
	};	
	return $accesser;
}

#This function remove extra white space in a fasta file, so the sequence remains in one line
#It also creates an index table of seq name to the which line in the file
#corresponding to that sequence
sub processFastaFile{
	my $tmp_file = getFileBasename($_[0]);
	my $out_dir = $_[1];
	$tmp_file =  "$out_dir/$tmp_file\.tmp.fa";
	open FASTA, $_[0] || dieLog "ERROR: cannot open fasta file";
	open TMP, ">".$tmp_file || dieLog "ERROR: cannot access file system for temporary file";
	my %fastaTable=();
	my $curr="";
	my $currHeader="";
	my $numLines=1;
	while(<FASTA>){
		my $line = $_;
		chomp $line;
		if($line =~ /^>/){
			if(length($curr) > 0){
				$fastaTable{substr($currHeader,1)}=$numLines;
				print TMP $currHeader."\n";
				print TMP $curr."\n";
				$numLines=$numLines+2;
			}else{
				warn("Sequence  ".$currHeader." has zero length...\n")
			}
			$currHeader = $line;
			$curr="";
		}else{
			$curr=$curr.$line;
		}
	}
	#taking care of the last seq contig
	$fastaTable{substr($currHeader,1)}=$numLines;	
	print TMP $currHeader."\n";
	print TMP $curr."\n";
	
	close(FASTA);
	close(TMP);
	return (\%fastaTable, $tmp_file)
}




##########################Functions for processing the hybrid-scaffold to NGS contigs alignment (i.e. the xmap file) ##########################

#This function read-in the alignments btw hybrid-scaffold contigs and the NGS contigs and compute gap length
#bewteen each consecutive pair of scaffolded NGS contigs and output this information into a file

sub processAlign{
	my $xmap = $_[0]; 
	my $gapLen_out_file = $_[1];
	my %alignMap = %{getAlignMap($xmap)};
	my $agp_out_file = $_[2];
	my $begin_end_file = $_[3];
	my $paddingGapLen = $_[4];
	my $ngs_map = $_[5];
	$xmap = mapNGSStat($xmap, $ngs_map);
	$xmap = computeGapLengths($xmap);
	printGapFile($xmap, $gapLen_out_file);
	printAGPFile($xmap, $agp_out_file, $begin_end_file, $paddingGapLen, $ngs_map);
}


#This function maps the QryContigID in xmap to the ngs stats such as name and length 
#from the original ngs contigs 
sub mapNGSStat{
	my $xmap = $_[0];
	my $ngs_map = $_[1];
	$xmap->{hits}->{NGSName} = rep("", $xmap->{totalHits});
	$xmap->{hits}->{NGSLen} = rep(0, $xmap->{totalHits});
	$xmap->{hits}->{NGSStart} = rep(1, $xmap->{totalHits});
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my $ID = $xmap->{hits}->{QryContigID}[$i];
		if(exists $ngs_map->{$ID}){
			#print "NSG name: ".$ngs_map->{$ID}->[0]."\n"; 
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{$ID}->[0];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{$ID}->[1];
			$xmap->{hits}->{NGSStart}[$i] = $ngs_map->{$ID}->[3];
		}else{
			print "Warning: cannot find ngs name for QryContig ID: ".$ID.". Using stats from qry in Xmap instead.\n";
			#if we cannot find ngs for some reasons use stat from original xmap
			$xmap->{hits}->{NGSName}[$i] =$ngs_map->{hits}->{QryContigID}[$i];
			$xmap->{hits}->{NGSLen}[$i] =$ngs_map->{hits}->{QryLen}[$i];
		} 
	}
	return $xmap;
}


#This function group alignments that has the same RefContigID
#It return a map of row indices of all alignment has the same RefContigID
sub getAlignMap{
	my $xmap = $_[0];	
	my %alignMap = ();
	print "Total alignments: ".($xmap->{totalHits})."\n";
	for(my $i = 0; $i < $xmap->{totalHits}; $i++){
		my $hybridContigID = $xmap->{hits}->{RefContigID}[$i];
		#print "hybrid: $hybridContigID\n";
		my @alignmentInds=();
		if(exists $alignMap{$hybridContigID}){
			@alignmentInds = @{$alignMap{$hybridContigID}};
		}
		#print "i: $i\n";
		push(@alignmentInds, $i);
		#print "arry: @alignments\n";
		$alignMap{$hybridContigID}=\@alignmentInds;
	}
	return \%alignMap;
}


#check direction
sub isForwardDir{
	return $_[0] eq '+';
}	

#adjust the ref start/end positions of alignments taking into account
#un-aligned sequence and direction of the query contig
sub adjustRefStart{
	my($qryStart, $qryEnd, $refStart, $qryLen, $direction) = @_;
	return $refStart - beginGap($qryStart, $qryEnd, $qryLen, $direction);	
}

sub beginGap{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;
	
	if(isForwardDir($direction)){
		return ($qryStart - 1);
	}else{
		return ($qryLen - $qryStart);
	}
}

sub adjustRefEnd{
	my($qryStart, $qryEnd, $refEnd, $qryLen, $direction) = @_;
	#print "$qryStart\t$qryEnd\t$refEnd\t$qryLen\t$direction\n";
	return $refEnd + EndGap($qryStart, $qryEnd, $qryLen, $direction);
}

sub EndGap{
	my($qryStart, $qryEnd, $qryLen, $direction) = @_;

	#print "$qryStart\t$qryEnd\t$refEnd\t$qryLen\t$direction\n";
	if(isForwardDir($direction)){
		return ($qryLen - $qryEnd);
	}else{
		return ($qryEnd - 1);
	}
}

#compute gaps for a particular hybrid-scaffolded contig
#given the alignments of all ngs contigs to the hybrid
sub computeGapLengths{
	my $xmap = $_[0];
	#the xmap data table will be used as the main table to store data
	#we added three-more columns to store the gap length information
	$xmap->{hits}->{GapLen} = rep(0, $xmap->{totalHits});
	$xmap->{hits}->{AdjustedGapLen} = rep(0, $xmap->{totalHits}); 
	$xmap->{hits}->{IsEmbedded} = rep(0, $xmap->{totalHits});
	
	my $alignGapLen = 0;
	my $adjustedGapLen = 0;
	for(my $i = 0; $i < $xmap->{totalHits}-1; $i++){
		my $refID1 = $xmap->{hits}->{RefContigID}[$i];
		my $refID2 = $xmap->{hits}->{RefContigID}[$i+1];
		if($refID1 ne $refID2){
			$alignGapLen = 0;
			$adjustedGapLen = 0;
		}else{
			my $qryStartPos1 = $xmap->{hits}->{QryStartPos}[$i];
			my $qryStartPos2 = $xmap->{hits}->{QryStartPos}[$i+1];
			my $qryEndPos1 = $xmap->{hits}->{QryEndPos}[$i];
			my $qryEndPos2 = $xmap->{hits}->{QryEndPos}[$i+1];
		
			my $refStartPos1 = $xmap->{hits}->{RefStartPos}[$i];
			my $refStartPos2 = $xmap->{hits}->{RefStartPos}[$i+1];
			my $refEndPos1 = $xmap->{hits}->{RefEndPos}[$i];
			my $refEndPos2 = $xmap->{hits}->{RefEndPos}[$i+1];

			my $qryLen1 = $xmap->{hits}->{QryLen}[$i];
			my $qryLen2 = $xmap->{hits}->{QryLen}[$i+1];
		
			my $direction1 = $xmap->{hits}->{Orientation}[$i];
			my $direction2 = $xmap->{hits}->{Orientation}[$i+1];
			
			my $qryID1 = $xmap->{hits}->{QryContigID}[$i];
			my $qryID2 = $xmap->{hits}->{QryContigID}[$i+1];
						
			$alignGapLen = round($refStartPos2) - round($refEndPos1);  #11/12/2015 remove plus one from gap length computation
			my $adjustedRefEnd1 = adjustRefEnd($qryStartPos1, $qryEndPos1, $refEndPos1, $qryLen1, $direction1);
			my $adjustedRefStart2=adjustRefStart($qryStartPos2, $qryEndPos2, $refStartPos2, $qryLen2, $direction2);
			if($DEBUG){
				print "$qryStartPos1\t$qryEndPos1\t$qryLen1\t$refStartPos1\t$refEndPos1\t$direction1\t$qryID1\t$refID1\t$adjustedRefEnd1\n";
				print "$qryStartPos2\t$qryEndPos2\t$qryLen2\t$refStartPos2\t$refEndPos2\t$direction2\t$qryID2\t$refID2\t$adjustedRefStart2\n";
			}
			
			$adjustedGapLen = $adjustedRefStart2 - $adjustedRefEnd1;
			if($DEBUG){print("Gap is ".$alignGapLen."\t".$adjustedGapLen."\n");}
			#detecting contigs that is completely covered by the previous contigs
			if(-1*$adjustedGapLen  > $qryLen2){
				$xmap->{hits}->{IsEmbedded}[$i+1] = 1;
			}
		}
		$xmap->{hits}->{GapLen}[$i] = $alignGapLen;
		$xmap->{hits}->{AdjustedGapLen}[$i] = $adjustedGapLen;
		
	}
	return $xmap;
}



###########################Functions for outputting agp/fasta or its auxilliary files###########################################


#print the gap informaton from alignments to a file
sub printGapFile{
	my $xmap = $_[0];
	my $out_put = $_[1];
	print "Gap length output file $out_put\n";
	#print "Hello\n";
	open(my $fh, ">".$out_put);
	print $fh "#NGSId1\tNGSId2\tSuperScaffoldId\tXmapGapLength\tAdjustedGapLength\tNGSLength1\tNGSLength2\n";
	my $numAlign = $xmap->{totalHits};
	for(my $i = 0; $i < $numAlign; $i++){
		if($i < $numAlign - 1 
			&& $xmap->{hits}->{RefContigID}[$i] == $xmap->{hits}->{RefContigID}[$i+1]){
			print $fh $xmap->{hits}->{NGSName}[$i]."\t".
				$xmap->{hits}->{NGSName}[$i+1]."\t".
				$xmap->{hits}->{RefContigID}[$i]."\t".
				$xmap->{hits}->{GapLen}[$i]."\t".
				$xmap->{hits}->{AdjustedGapLen}[$i]."\t".
				$xmap->{hits}->{QryLen}[$i]."\t".
				$xmap->{hits}->{QryLen}[$i+1]."\n";
		}else{
			print $fh $xmap->{hits}->{NGSName}[$i]."\t".
				$xmap->{hits}->{NGSName}[$i]."\t".
				$xmap->{hits}->{RefContigID}[$i]."\t".
				$xmap->{hits}->{GapLen}[$i]."\t".
				$xmap->{hits}->{AdjustedGapLen}[$i]."\t".
				$xmap->{hits}->{QryLen}[$i]."\t".
				$xmap->{hits}->{QryLen}[$i]."\n";
				
		}
	}
}

#print header of AGP file
sub printAGPHeader{
my $fh = $_[0];
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $agp_header = "##agp-version\t2.0\n".
				"# Organism:  \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
             "# Obj_Name\tObj_Start\tObj_End\tPartNum\tCompnt_Type\tCompntId_GapLength\tCompntStart_GapType\tCompntEnd_Linkage\tOrientation_LinkageEvidence\n";
	print $fh $agp_header;
}


sub printAGPAuxHeader{
	my $fh = $_[0];
	my $header= "##agp-version\t2.0\n".
				"# Organism:   \n".
				"# Platform:     \n".
				"# Model:        \n".
				"# Enzyme(s):    \n".
				"# BioSample:    \n".
				"# BioProject:   \n".
				"Obj_Id\tHeadTrimmedLength\tTailTrimmedLength\n";
	print $fh $header;
}



#print the agp files
sub printAGPFile{
	my ($xmap, $out_put, $aux_out_file, $dummyGapLen, $ngs_map)= @_;
	print "AGP output file $out_put\n";
	open(my $fh, ">".$out_put);
	printAGPHeader($fh);
	open(my $aux_fh, ">".$aux_out_file);
	printAGPAuxHeader($aux_fh);
	my $numAlign = $xmap->{totalHits};
	my $currRefStart = 1;
	my $currRefEnd = 1;
	my $Ind = 1;
	
	#outputing ngs contigs
	my $currInd = 0;
	my $nextInd = 1;	
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			#marking the ngs contig as embedded
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = -1;
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			#marking the ngs contig as embedded
			$ngs_map->{$xmap->{hits}->{QryContigID}[$nextInd]}->[2] = -1;
			$nextInd = $nextInd + 1;
			next;
		}
					
		#printing begin gap for first contig
		#note first contig cannot be embedded so this condition is valid
		if($currInd == 0){
			my $begin_gap = $xmap->{hits}->{RefStartPos}[$currInd]-1;
			print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t$begin_gap\t";
		}
		#printing ngs sequence
		#print "Index $i\n";
		$currRefEnd = $currRefStart + round($xmap->{hits}->{NGSLen}[$currInd]) - 1;
		print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
				"$currRefStart\t".	
				$currRefEnd."\t".
				"$Ind\tW\t".
				$xmap->{hits}->{NGSName}[$currInd]."\t".
				#"1\t".$xmap->{hits}->{QryLen}[$currInd]."\t".
				$xmap->{hits}->{NGSStart}[$currInd]."\t".($xmap->{hits}->{NGSStart}[$currInd]+$xmap->{hits}->{NGSLen}[$currInd]-1)."\t".
				$xmap->{hits}->{Orientation}[$currInd]."\n";
		$currRefStart = $currRefEnd + 1;
		#marking in the ngs table that this contig was scaffolded
		if(exists $ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}){
			$ngs_map->{$xmap->{hits}->{QryContigID}[$currInd]}->[2] = 1;
		}
		
		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			my $gapLen = $xmap->{hits}->{GapLen}[$currInd];
			if($gapLen < 0){
				$gapLen = $dummyGapLen;
			}
			$currRefEnd = $currRefStart + $gapLen - 1;
			if($gapLen > 0){
				$Ind++;
				print $fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\t".
					"$currRefStart\t".	
					"$currRefEnd\t".
					"$Ind\tN\t".$gapLen."\tscaffold\tyes\tmap\n";
			}
			$currRefStart = $currRefEnd + 1;		
		}else{								
			#handle beginning and end gap sequence
			#end gap for the previous hybrid contig
			my $end_gap = $xmap->{hits}->{RefLen}[$currInd] - $xmap->{hits}->{RefEndPos}[$currInd];
			#print "printing END gap: ".$xmap->{hits}->{RefLen}[$currInd]."\t". $xmap->{hits}->{RefEndPos}[$currInd]."\n";
			print $aux_fh round($end_gap)."\n"; 
			#begin gap for the current hybrid contig
			if($nextInd < $numAlign){ 
				my $begin_gap = $xmap->{hits}->{RefStartPos}[$nextInd]-1;
				print $aux_fh "Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\t".round($begin_gap)."\t";
			}
			$currRefStart = 1;
			$currRefEnd = 0;
			$Ind = 0;
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		$Ind++;
	}
	close($fh);
	close($aux_fh);
}


#convert the interval objects to fasta file format
sub printFasta{
	my ($fasta_map, $xmap, $out_file, $hybridCmap, $cut_site, $paddingGap) = @_;
	my $numAlign = $xmap->{totalHits};
	my $isBegin = 1;
	open(my $fh, ">".$out_file) || die "cannot open output fasta file $out_file\n";
	
	my $currInd = 0;
	my $nextInd = 1;	
	while($currInd < $numAlign){
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			$currInd = $currInd + 1;
			next;
		}
		
		if($xmap->{hits}->{IsEmbedded}[$nextInd]){
			$nextInd = $nextInd + 1;
			next;
		}

		my $numAlign = $xmap->{totalHits};
		#if contig is embbed inside another contig, skip this contig
		if($xmap->{hits}->{IsEmbedded}[$currInd]){
			next;
		}
		#handling the first line
		if($isBegin){
			print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$currInd]."\n";
			$isBegin=0;
		}
		
		#printing ngs sequence
		print "NGSname $xmap->{hits}->{NGSName}[$currInd]\n";
		my $ngsseq = $fasta_map->($xmap->{hits}->{NGSName}[$currInd]);
		if(1 | $DEBUG){print "NGS seq length ".length($ngsseq)." expected: ".$xmap->{hits}->{NGSLen}[$currInd]."\n"};
		if(isForwardDir($xmap->{hits}->{Orientation}[$currInd])){
			print $fh $ngsseq;
		}else{
			print $fh reverseComplement($ngsseq);
		}

		#printing gap		
		if($nextInd < $numAlign && $xmap->{hits}->{RefContigID}[$currInd] == $xmap->{hits}->{RefContigID}[$nextInd]){
			my $gapBegin = round($xmap->{hits}->{RefEndPos}[$currInd]);
			my $gapEnd = round($xmap->{hits}->{RefStartPos}[$nextInd])-1; #gap end at one b.p.before next start
			my $RefContigID = $xmap->{hits}->{RefContigID}[$currInd];
			my $gapSeq = $paddingGap;
			if($gapBegin < $gapEnd){
				$gapSeq = getSeqFromCmap([$gapBegin, $gapEnd], $hybridCmap, $RefContigID, $cut_site);
			}
			if(1 | $DEBUG){print "Gap seq length: ".length($gapSeq)." expected ".$xmap->{hits}->{GapLen}[$currInd]."\n"};
			print $fh $gapSeq;
		}else{								
			print $fh "\n";	
			if($nextInd < $numAlign){
				print $fh ">Super-Scaffold_".$xmap->{hits}->{RefContigID}[$nextInd]."\n";
			}
		}
		$currInd = $nextInd;
		$nextInd = $nextInd + 1;
		
	}
	close($fh);
}


#printing out a sequence accoriding to Cmap coordinate
#This is essentially a N-sequence with location of cut sites
#specified by the cut-site sequences of the enzyme/enzymes

sub getSeqFromCmap{
	my @interval = @{$_[0]}; #begin end coordinate of interval
	if($interval[0] >= $interval[1]){
		print STDERR "Begin index is not smaller than end Index in interval: \[$interval[0], $interval[1]\]\n";
		return "";
	}
	my $cmap = $_[1]; #cmap pointer
	my $contigID = $_[2]; #contigID of the reference contig
	#print "hybrid map ID: ".$interval[$refIDInd]; 
	my $contig =$cmap->{contigs}->{$contigID};
	my $cut_site = $_[3];
	my $seq="";
	(my $begin, my $end) = getCoordinatesSubset(\@interval, $contig);

	#print "hybrid interval $interval[0]\t$interval[1]\n";
	my $coord = $interval[0];
	for(my $i = $begin; $i < $end; $i++){
		#print "end pos: ".($contig->{$Position}[$i])."\n";
		for(; $coord < $contig->{Position}[$i]; $coord++){
			$seq = $seq."N";
		}
		$seq=$seq.$cut_site;
		$coord = $coord+length($cut_site);
	}
	for(; $coord <= $interval[1]; $coord++){
			$seq = $seq."N";
	}
	if($DEBUG){print "length of gap: ".(length($seq))."\n"};
	return $seq;
}

#return begin and end index (right exclusive) in a cmap contig within 
#a specified coordinate range note that to avoid scanning 
#the whole cmap data all the time, one can "batch lookup" 
#the coordinates in successive calls of this function on
#many non-overlapping intervals in increasing coordinate order
#the method will remember the last coordinate left-off using
#a static/state variable

{
	my $prevInd = 0;
	my $prevContig;
	sub getCoordinatesSubset{
		my @interval = @{$_[0]};
		my $contig = $_[1];
		my $begin = 0;
		my $end = 0;
		#resetting the remembered index
		if($DEBUG){print "previous stopping at $prevInd\t$prevContig\t$contig\t".($prevContig == $contig)."\n"};
		if(!(defined $prevInd) || !(defined $prevContig) || !($prevContig == $contig)){
			#print "resetting\n";
			$prevInd = 0;
		}
		#user can reset the index and decide to scan from a specified index
		if($#_ == 2){
			$prevInd = $_[2];
		}
		
		my @positions = @{$contig->{Position}};
		for(my $i = $prevInd; $i < $#positions; $i++){
			if($DEBUG){print "$interval[0]\t$interval[1]\t".$positions[$i]."\n"};
			if($positions[$i] > $interval[0] && $begin <= 0){
				$begin = $i;
			}
			
			if($positions[$i] > $interval[1]){
				$end = $i;
				$prevInd = $end-1;
				last;
			}
			
			if($i == $#positions && $begin > 0){
				$end = $i+1;
				$prevInd = $end;
			}
		}
		$prevContig = $contig;
		if($begin > $end){
			$begin = 0;
			$end = 0;
		}
		if($DEBUG){print "Interval $begin $end\n";}
		return ($begin, $end);		
	}
}


sub reverseComplement{
	my $str = $_[0];
	#print "forward: ".substr($str, length($str)-50, 50)."\n";
	my @match_forward = $str =~ /GCTCTTC/g;
	$str =~ tr/ATCGatcg/TAGCtagc/;
	$str = reverse($str);
	#my @match_reverse = $str =~ /GAAGAGC/g;
	#if($#match_forward != $#match_reverse){
	#	print "foward and reverse matches are not the same: $#match_forward\t$#match_reverse\n"
	#}
	#print "rcmp: ".substr(reverse($str), 0, 50)."\n";
	return $str;
}



#printing out the un-used ngs object as "singleton" to agp file as well as fasta
#file
#@TODO when only part of the ngs contig is not used only print out 
#part of the contigs
sub printUnUsedNGS(){
	my $agp_out_file = $_[0];
	my %ngs_map = %{$_[1]};
	my $aux_out_file = $_[2];
	my $fasta_out=$_[3]; 
	my $seqMap = $_[4];
	open(my $fh, ">>".$agp_out_file) || die "Cannot open agp output file";
	open(my $aux_fh, ">>".$aux_out_file) || die "Cannot open agp auxilliary file";
	open(my $fasta_fh, ">>".$fasta_out) || die "Cannot open fasta output file";
	
	my $unUsedCount=0;
	foreach my $key(keys %ngs_map){
		my @ngs = @{$ngs_map{$key}};
		my $ngsStart = $ngs[3];
		my $ngsEnd = $ngsStart + $ngs[1] - 1;
		#print "key is: $key\t$ngs[2]\n";
		if($ngs[2] <= 0){
			print $fh "$ngs[0]_obj\t1\t$ngs[1]\t1\tW\t$ngs[0]\t1\t$ngs[1]\t\+\n";
			print $aux_fh "$ngs[0]_obj\t0\t0\n";
			#if($ngs[2] < 0){
				#print $fasta_fh ">".$ngs[0]."_overlap\n";	
			#}else{
			if(length($seqMap->($ngs[0])) > 0){				
				print $fasta_fh ">".$ngs[0]."\n";
				print $fasta_fh $seqMap->($ngs[0])."\n";
			}else{
				warn("Sequence $ngs[0] has zero length, skipping it in output fasta\n");
			}
			$unUsedCount++;
		}
	}
	close($fh);
	close($aux_fh);
	close($fasta_fh);
	#print "Unused NGS: $unUsedCount\n";
}

sub cleanUp(){
	my %INPUTS = %{$_[0]};
	my $dir = $INPUTS{out_dir};
	opendir(DIR, $dir);
	while(my $file = readdir(DIR)){
		if($file =~ /tmp\.fa/){
			print "Removing temp file: ". $file."\n";
			unlink $dir."/".$file;
		}
	}
	close(DIR);
}

