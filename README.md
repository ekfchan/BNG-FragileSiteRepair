# BNG-FragileSiteRepair
Scripts for repairing (stitching) fragile sites of BioNano genome maps. 

Usage: perl fragileSiteRepair.pl --fasta [reference.fasta] --xmap [input.xmap] --qcmap [input_q.cmap] --rcmap [input_r.cmap] --errbin [input.errbin] --output [output folder] --maxlab [max_label_gap_tolerence] --maxfill [max basepairs to fill between contigs] --wobble [fragile site wobble in bp] --force [overwrite output folder]


Installing dependencies:
1) login as root [or use sudo]
2) yum install perl-CPAN.noarch perl-Config-Tiny.noarch perl-Data-Dumper.x86_64 perl-ExtUtils-CBuilder.noarch perl-ExtUtils-ParseXS.noarch perl-HTML-Parser.x86_64 perl-IPC-Run.noarch perl-Module-Build.noarch perl-Module-CoreList.noarch perl-Module-Install.noarch perl-Module-Metadata.noarch perl-Sys-MemInfo.x86_64 perl-TimeDate.noarch perl-XML-Parser.x86_64 perl-XML-Simple.noarch perl-devel.x86_64 perl-local-lib.noarch perl-threads.x86_64 perl-threads-shared.x86_64 perl-App-cpanminus.noarch perl-CPAN-Meta.noarch perl-Config-Simple.noarch perl-DateTime.x86_64 perl-ExtUtils-Install.noarch perl-ExtUtils-MakeMaker.noarch perl-ExtUtils-Manifest.noarch perl-File-Copy-Recursive.noarch perl-File-Path.noarch perl-File-Slurp.noarch perl-IPC-Run3.noarch perl-Sys-CPU.x86_64 perl-Time-HiRes.x86_64 perl-YAML.noarch
3) cpanm install -f IPC::System::Simple Sys::Info Bio::Seq Parallel::ForkManager DateTime::Format::Human::Duration
