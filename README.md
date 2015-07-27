#BNG-FragileSiteRepair
Scripts for repairing (stitching) fragile sites of BioNano genome maps. 

```
Usage: perl fragileSiteRepair.pl [--bed reference_fsites.bed] [--fasta reference.fasta] --xmap input.xmap --qcmap input_q.cmap --rcmap input_r.cmap --errbin input.errbin --output output_folder --maxlab max_label_gap_tolerence --maxfill max_basepairs_to_fill_between_contigs --wobble fragile_site_wobble_in_bp --force 
```
``` 
OPTIONS
	--bed reference_fsites.bed
			Optional. BED file of fragile sites predicted from reference.fasta. 
			If provided, calcFragileSites will not be performed and --fasta will be ignored. 
			File need not be in same dir as scripts or cwd. 
	--fasta reference.fasta
			Optional. FASTA of reference genome. If --bed not provided, fragile sites 
			will be calculated based on this FASTA. 
			File need not be in same dir as scripts or cwd. 
	--xmap input.xmap
			Required. Alignment file from Pipeline from which fragile site stitching is 
			to be performed. 
			File need not be in same dir as scripts or cwd. 
	--qcmap input_q.cmap
			Required. Query cmap corresponding to --xmap
			File need not be in same dir as scripts or cwd. 
	--rcmap input_r.cmap
			Required. Reference cmap corresponding to --xmap
			File need not be in same dir as scripts or cwd. 
	--errbin input.errbin
			Required. Error parameters corresponding to --xmap.  
			File need not be in same dir as scripts or cwd. 
			**TODO:  May be possible to make this optional in the future**
	--output output_folder
			Required.  Directory to store outputs. 
			Dir need not be in same as scripts or in cwd. 
	--maxlab integer
			Optional. The maximum label allowed between contigs. Default = 0.
	--maxfill integer
			Optional. The maximum bases to fill between contigs. Default = 10000. 
	--wobble integer
			Optional. The basepair tolerance allowed between fragile sites. Default = 250. 
	--force 
			Optional. To force overwrite of output folder. Deafult off. 
```

###Installing dependencies:
1. Install BioNano Genomics de-novo assembly Linux pipeline in user's home directory (~/scripts and ~/tools) according to [instructions](http://www.bnxinstall.com/training/docs/IrysViewSoftwareInstallationGuide.pdf)  
2. login as root [or use sudo]
3. `yum install perl-CPAN.noarch perl-Config-Tiny.noarch perl-Data-Dumper.x86_64 perl-ExtUtils-CBuilder.noarch perl-ExtUtils-ParseXS.noarch perl-HTML-Parser.x86_64 perl-IPC-Run.noarch perl-Module-Build.noarch perl-Module-CoreList.noarch perl-Module-Install.noarch perl-Module-Metadata.noarch perl-Sys-MemInfo.x86_64 perl-TimeDate.noarch perl-XML-Parser.x86_64 perl-XML-Simple.noarch perl-devel.x86_64 perl-local-lib.noarch perl-threads.x86_64 perl-threads-shared.x86_64 perl-App-cpanminus.noarch perl-CPAN-Meta.noarch perl-Config-Simple.noarch perl-DateTime.x86_64 perl-ExtUtils-Install.noarch perl-ExtUtils-MakeMaker.noarch perl-ExtUtils-Manifest.noarch perl-File-Copy-Recursive.noarch perl-File-Path.noarch perl-File-Slurp.noarch perl-IPC-Run3.noarch perl-Sys-CPU.x86_64 perl-Time-HiRes.x86_64 perl-YAML.noarch`
4. `cpanm install -f IPC::System::Simple Sys::Info Bio::Seq Parallel::ForkManager DateTime::Format::Human::Duration`
