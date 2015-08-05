#!/usr/bin/perl

package stitchFsites;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

BEGIN { 
	require Exporter; 
	our $VERSION     = 1.00;
	our @ISA         = qw(Exporter);
	our @EXPORT      = ();
	our @EXPORT_OK   = qw(mergeContigs doalignment appendStitched getXmapOpts);
}

sub mergeContigs {
	# usage: mergeContigs(\%q_cmap{$firstQryContigID}, \%q_cmap{$secondQryContigID}, $firstOrientation, $secondOrientation, abs($secondRefStartPos - $firstRefEndPos))
	my ($firstContigRef, $secondContigRef, $firstOrientation, $secondOrientation, $positionOffset, $verbose) = @_; #$firstCOntigRef and $secondContigRef are referencess to hashes with keys: ContigLength, NumSites, and Positions of each label site 
print $verbose; 
	if ( $verbose ) { print "\t\tOffset: $positionOffset\n"; }
	
	my $firstContigSites = $firstContigRef->{'NumSites'};	# number of sites in 1st contig
	my $firstContigLength = $firstContigRef->{'ContigLength'};	# length of 1st contig
	my @firstPositions = (keys %$firstContigRef);	#all label positions in 1st contig
	@firstPositions = grep ! /ContigLength/, @firstPositions;
	@firstPositions = grep ! /NumSites/, @firstPositions;
	@firstPositions =  sort {$a <=> $b} @firstPositions;	#sort them
	my $firstContigEnd = $firstPositions[$#firstPositions - 1];	#actual last label position (not lenght of contig) 

	my $secondContigSites = $secondContigRef->{'NumSites'}; #number of sites in 2nd contig
	my $secondContigLength = $secondContigRef->{'ContigLength'}; #length of 2nd contig
	my @secondPositions = (keys %$secondContigRef);
	@secondPositions = grep ! /ContigLength/, @secondPositions;
	@secondPositions = grep ! /NumSites/, @secondPositions;
	@secondPositions =  sort {$a <=> $b} @secondPositions;
	my $secondContigEnd = $secondPositions[$#secondPositions - 1]; #actual last label position (not contig length)

	my $mergedSites = $firstContigSites + $secondContigSites;
	if( $positionOffset==0 ) { $mergedSites = $mergedSites-1; } #if contig pair start/end at same label, do not double count the overlapping label
	my %newQCmap; 
	$newQCmap{'NumSites'} = $mergedSites; 

	if ( $firstOrientation eq '+' && $secondOrientation eq '+' ) {
		my $mergedLength = $firstContigEnd + $positionOffset + ($secondContigLength - $secondPositions[0]); 
		$newQCmap{'ContigLength'} = $mergedLength; 
		print "\t\tMerged Contig: $mergedSites sites, $mergedLength bp\n"; 

		my $cursite = 1; 
		for ( my $i=0; $i < $#firstPositions; $i++ ) { #excludes dude label where SiteID > NumSites
			my $pos = $firstPositions[$i]; 
			$newQCmap{$pos}{'SiteID'} = $cursite++;
			$newQCmap{$pos}{'LabelChannel'} = $firstContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$pos}{'TheRest'} = $firstContigRef->{$pos}{'TheRest'}; 
			if( $verbose ) { print "\t\t\t\t1\t",sprintf("%.1f",$pos),"\t$newQCmap{$pos}{'SiteID'}\t$newQCmap{$pos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $firstContigRef->{$pos}{'SiteID'},"\n"; }
		}
		for ( my $i=0; $i <= $#secondPositions; $i++ ) { #include and reuse dude label where SiteID > NumSites
			if( $positionOffset==0 && $i==0 ) { next; }	#if contig pair start/end at same label, exclude start label of second contig
			my $pos = $secondPositions[$i]; 
			my $newpos = $firstContigEnd + $positionOffset + ($pos - $secondPositions[0]); 
			$newQCmap{$newpos}{'LabelChannel'} = $secondContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $secondContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $secondContigRef->{$pos}{'SiteID'},"\n"; }
		}
	}
	elsif ( $firstOrientation eq '+' && $secondOrientation eq '-' ) {
		my $mergedLength = $firstContigEnd + $positionOffset + $secondContigEnd;
		$newQCmap{'ContigLength'} = $mergedLength; 
		print "\t\tMerged Contig: $mergedSites sites, $mergedLength bp\n"; 

		my $cursite=1; 
		for ( my $i=0; $i < $#firstPositions; $i++ ) { #excludes dude label where SiteID > NumSites
			my $pos = $firstPositions[$i]; 
			$newQCmap{$pos}{'SiteID'} = $cursite++;
			$newQCmap{$pos}{'LabelChannel'} = $firstContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$pos}{'TheRest'} = $firstContigRef->{$pos}{'TheRest'}; 
			if( $verbose ) { print "\t\t\t\t1\t",sprintf("%.1f",$pos),"\t$newQCmap{$pos}{'SiteID'}\t$newQCmap{$pos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $firstContigRef->{$pos}{'SiteID'},"\n"; }
		}
		for ( my $i=($#secondPositions-1); $i >= 0; $i-- ) { #excludes dude label where SiteID > NumSites
			if( $positionOffset==0 && $i==($#secondPositions-1) ) { next; }	#if contig pair start/end at same label, exclude start (i.e. reversed end) label of second contig
			my $pos = $secondPositions[$i]; 
			my $newpos = $firstContigEnd + $positionOffset + ($secondContigEnd-$pos); 
			$newQCmap{$newpos}{'LabelChannel'} = $secondContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $secondContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $secondContigRef->{$pos}{'SiteID'},"\n"; }
		}
		# append end of contig (dud label site) 
		$newQCmap{$mergedLength}{"SiteID"} = $cursite++;
		$newQCmap{$mergedLength}{"LabelChannel"} = "0";
		$newQCmap{$mergedLength}{"TheRest"} = $secondContigRef->{$secondPositions[$#secondPositions]}{'TheRest'};
		if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$mergedLength),"\t$newQCmap{$mergedLength}{'SiteID'}\t$newQCmap{$mergedLength}{'TheRest'}\t<--\t",sprintf("%.1f",$secondPositions[0]),"\t", $secondContigRef->{$secondPositions[0]}{'SiteID'},"\n"; }
	}
	elsif ( $firstOrientation eq '-' && $secondOrientation eq '+' ) {
		my $mergedLength = ($firstContigLength - $firstPositions[0]) + $positionOffset + ($secondContigLength - $secondPositions[0]); 
		$newQCmap{'ContigLength'} = $mergedLength; 
		print "\t\tMerged Contig: $mergedSites sites, $mergedLength bp\n"; 
		
		my $cursite=1; 
		for ( my $i=($#firstPositions-1); $i >=0 ; $i-- ) { #exclude dude label where SiteID > NumSites
			my $pos = $firstPositions[$i];  
			my $newpos = $firstContigLength - $pos;
			$newQCmap{$newpos}{'LabelChannel'} = $firstContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $firstContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t1\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $firstContigRef->{$pos}{'SiteID'},"\n"; }
		}
		for ( my $i=0; $i <= $#secondPositions; $i++ ) { #include and reuse dude label where SiteID > NumSites
			if( $positionOffset==0 && $i==0 ) { next; }	#if contig pair start/end at same label, exclude start label of second contig
			my $pos = $secondPositions[$i]; 
			my $newpos = ($firstContigLength - $firstPositions[0]) + $positionOffset + ($pos - $secondPositions[0]); 
			$newQCmap{$newpos}{'LabelChannel'} = $secondContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $secondContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $secondContigRef->{$pos}{'SiteID'},"\n"; }
		}
	}
	elsif ( $firstOrientation eq '-' && $secondOrientation eq '-' ) {
		my $mergedLength = ($firstContigLength - $firstPositions[0]) + $positionOffset + $secondContigEnd; 
		$newQCmap{'ContigLength'} = $mergedLength; 
		print "\t\tMerged Contig: $mergedSites sites, $mergedLength bp\n"; 
		
		my $cursite=1; 
		for ( my $i=($#firstPositions-1); $i >=0 ; $i-- ) { #exclude dude label where SiteID > NumSites
			my $pos = $firstPositions[$i]; 
			my $newpos = $firstContigLength - $pos; 
			$newQCmap{$newpos}{'LabelChannel'} = $firstContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $firstContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t1\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $firstContigRef->{$pos}{'SiteID'},"\n"; }
		}
		for ( my $i=($#secondPositions-1); $i >= 0; $i-- ) { #excludes dude label where SiteID > NumSites
			if( $positionOffset==0 && $i==($#secondPositions-1) ) { next; }	#if contig pair start/end at same label, exclude start (i.e. reversed end) label of second contig
			my $pos = $secondPositions[$i]; 
			my $newpos = ($firstContigEnd - $firstPositions[0]) + $positionOffset + ($secondContigEnd-$pos); 
			$newQCmap{$newpos}{'LabelChannel'} = $secondContigRef->{$pos}{'LabelChannel'}; 
			$newQCmap{$newpos}{'TheRest'} = $secondContigRef->{$pos}{'TheRest'}; 
			$newQCmap{$newpos}{'SiteID'} = $cursite++;
			if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$newpos),"\t$newQCmap{$newpos}{'SiteID'}\t$newQCmap{$newpos}{'TheRest'}\t<--\t",sprintf("%.1f",$pos),"\t", $secondContigRef->{$pos}{'SiteID'},"\n"; }
		}
		# append end of contig (dud label site) 
		$newQCmap{$mergedLength}{"SiteID"} = $cursite++;
		$newQCmap{$mergedLength}{"LabelChannel"} = "0";
		$newQCmap{$mergedLength}{"TheRest"} = $secondContigRef->{$secondPositions[$#secondPositions]}{'TheRest'};
		if( $verbose ) { print "\t\t\t\t2\t",sprintf("%.1f",$mergedLength),"\t$newQCmap{$mergedLength}{'SiteID'}\t$newQCmap{$mergedLength}{'TheRest'}\t<--\t",sprintf("%.1f",$secondPositions[0]),"\t", $secondContigRef->{$secondPositions[0]}{'SiteID'},"\n"; }
	}

	return( \%newQCmap ); 
}


sub doalignment {

	my ($xmap, $rcmap, $qcmap, $out, $errbin, $mem, $cpuCount, $verbose) = @_; 
	if( $verbose ) {
		print "XMAP: $xmap\n";
		print "Reference CMAP: $rcmap\n";
		print "Query CMAP: $qcmap\n";
		print "OUTPUT prefix: $out\n";
		print "ERRBIN: $errbin\n";
	}

	# get RefAligner options
	my $opts = getXmapOpts($xmap); 
	# my $opts; 
	# open (XMAP, "<", $xmap) or die $!;
	# while (my $line = <XMAP>) {
	# 	chomp $line;
	# 	if ($line =~ m=tools/RefAligner=) { 
	# 		$opts = $line;
	# 		last;
	# 	} else { next; }
	# }
	# close XMAP; 
	# $opts =~ s/^.+RefAligner\s//;
	# $opts =~ s/-i\s[^\s]+\s//g; 
	# $opts =~ s/-o\s[^\s]+\s//g; 
	# $opts =~ s/-ref\s[^\s]+\s//g; 
	# $opts =~ s/-maxthreads\s[^\s]+\s//g; 
	# $opts =~ s/-maxmem\s[^\s]+\s//g; 
	# $opts =~ s/-output-veto-filter\s[^\s]+\s//g; 
	# $opts =~ s/-output-filter\s[^\s]+\s//g; 
	# $opts =~ s/-readparameters\s[^\s]+\s//g; 
	# $opts =~ s/-stderr\s+//g; 
	# $opts =~ s/-stdout\s+//g; 

	# RefAligner filters 
	my $ofilter = q/-output-filter '(.[cx]map)$'/;
	# my $ofilter = q/-output-filter '(q.cmap)$'/;
	my $veto = q/-output-veto-filter '(_intervals.txt|.err|.maprate|[a-z|A-Z].map)$'/;

	# RefAligner command
	my $cmd = "~/tools/RefAligner -ref ".$rcmap." -i ".$qcmap." -o ".$out."  -stdout ".$ofilter." -readparameters ".$errbin." -BestRef 1 -f "; 
	$cmd = $cmd." -maxmem ".$mem."  -maxthreads ".$cpuCount; 
	$cmd = $cmd." ".$opts; 

	if( $verbose ) { print "\nRunning command: \n\n $cmd \n\n"; }
	system($cmd); 

	if( $verbose ) {
		if ($? == -1) {
	        print "failed to execute: $!\n";
	    }
	    elsif ($? & 127) {
	        printf "child died with signal %d, %s coredump\n",
	            ($? & 127),  ($? & 128) ? 'with' : 'without';
	    }
	    else {
	        printf "child exited with value %d\n", $? >> 8;
	    }
	}
}

sub printstitched {
	# example usage: printstitched($outbed, $RefContigID, int($firstRefEndPos-1), int($secondRefStartPos+1), $fsiteFound);
	my ($outbed, $RefContigID, $leftPos, $rightPos, $fsiteline) = @_; 
	if (!-e $outbed) { 
		#print header 
		open ( OUT, ">", $outbed) or die "ERROR: Cannot open $outbed for writing! $!\n";
		print OUT "#CMapId\tStart\tEnd\tType\tScore\tStrand\tThickStart\tThickEnd\tItemRgba\tSequence\n"; 
	}
	else {
		# open file for appending
		open (OUT, ">>", $outbed) or die "ERROR: Cannot open $outbed for appending! $!\n";
	}
	#append stitch info
	my @s = split(/\t/,$fsiteline); 
	print OUT "$RefContigID\t$leftPos\t$rightPos\t", join("\t",@s[3..($#s)]), "\n"; 
	close OUT; 
}

sub appendStitched {
	# example usage: appendStitched($outbed, $RefContigID, int($firstRefEndPos-1), int($secondRefStartPos+1), $fsiteFound); 
	my ($outbed, $RefContigID, $leftPos, $rightPos, $fsiteline) = @_; 

	open (OUT, ">>", $outbed) or die "ERROR: Cannot open $outbed for writing! $!\n";
	my @s = split('\t',$fsiteline); 
	print OUT "$RefContigID\t$leftPos\t$rightPos\t", join("\t",@s[3..($#s)]), "\n"; 
	close OUT; 
}

sub getXmapOpts {
	#ussage: getXmapOpts(xmap_file)
	#returns:  string of options to RefAligner, with dataset specific opts stripped
	my $xmap = shift; 

	my $opts; 
	open (XMAP, "<", $xmap) or die $!;
	while (my $line = <XMAP>) {
		chomp $line;
		if ($line =~ m=tools/RefAligner=) { 
			$opts = $line;
			last;
		} else { next; }
	}
	close XMAP; 
	$opts =~ s/^.+RefAligner\s//;
	$opts =~ s/-i[f]?\s\S+(\s|$)//g;	#input maps
	$opts =~ s/-o\s\S+(\s|$)//g;	#output prefix
	$opts =~ s/-ref\s\S+(\s|$)//g;	#input reference cmap 
	$opts =~ s/-maxthreads\s\S+(\s|$)//g;	#maximum threads to use
	$opts =~ s/-maxmem\s\S+(\s|$)//g;	#max memory in Gbytes
	$opts =~ s/-output-veto-filter\s\S+(\s|$)//g;	#outputs to veto
	$opts =~ s/-output-filter\s\S+(\s|$)//g;	#outputs to keep 
	$opts =~ s/-readparameters\s\S+(\s|$)//g;	# error (noise) parameters
	$opts =~ s/-stderr(\s|$)//g;	#redirect stdout to file.stdout
	$opts =~ s/-stdout(\s|$)//g;	#redirect stderr to file.stdout, else to file.stdout

	return( $opts ); 
}

1; 
