# BNG-FragileSiteRepair
Scripts for repairing (stitching) across fragile sites breaks in BioNano genome maps.

```
Usage: perl fragileSiteRepair.pl [OPTIONS] --fasta <ref.fa> --cmap <assembly.cmap> --output <output folder>
```
``` 
=== REQUIRED ===
--fasta <.fasta> : input reference FASTA file to be used for fragile site prediction and alignment
--cmap <.cmap> : input BioNano assembly CMAP to perform fragileSiteRepair on
--output <path/to/output> : output folder

=== KEY OPTIONS ===
--bnx <.bnx> : input NMX used for assembly. Requires --err to perform single molecule alignments. Default: OFF
--err <.err> : molecules ERR file from assembly alignmol/merge OR autoNoise1. Requires --bnx to perform single molecule alignments. Default: OFF
--bam <.bam> : BAM file of NGS scaffolds/contigs/reads aligned to --fasta. Used to supplement scoring of stitchPositions. Default: OFF

--break : Flag to enable breaking of maps at stitchPositions with score less than threshold value. Requires --bnx and --err. Default: OFF
--runSV : Flag to run SV module on final cmap. Default: OFF
--force : Flag to overwrite --output folder. Default: OFF
--aggressive : Flag to calculate TypeIII and TypeIV fragile sites in addition. Default: OFF

=== OTHER OPTIONS ===
--enzyme <nickase sequence> : Nickase enzyme for in-silico digestion and fragile site prediction. Default: GCTCTTC
--ngsBuffer <basepairs> : Number of basepairs that a single NGS alignment must extend past the ends of a stitchPosition to supplement score. Default: 500
--threshold <score> : Minimum stitchPoisitons score below which to break maps. Default: 1.0
--maxlab <label count> : Maximum number of reference labels to allow between adjacent maps. Default: 1
--maxfill <basepairs> : Maximum number of basepairs to allow between adjacent maps. Default: 30000
--wobble <basepairs> : Maximum number of basepairs to allow the fragile site to vary. Default: 30000
--seq <basepairs> : Number of basepairs of reference sequence +/- fragile site to output into BED file. Default: OFF
--optArgs <optArguments.xml> : optArguments.xml to use for alignment. Default: ~/scripts/optArguments_human.xml
--endoutlier <pvalue> : endoutlier penalty for single molecule alignments. Default: 1e-3

--n <CPU cores> : Maximum number of CPU cores/threads to use. Default: nproc
--j <number jobs> : Maximum number of parallel jobs. Default: nproc/6
