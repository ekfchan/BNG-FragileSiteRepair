# BNG-FragileSiteRepair
Scripts for repairing (stitching) across fragile sites breaks in BioNano genome maps.

```
Usage: perl fragileSiteRepair.pl --fasta <reference.fasta> --cmap <assembly.cmap> --output <output folder> [--bnx <input.bnx to run alignmol>] [--err <original assembly alignmol_merge.err or autonoise1.err>] [--enzyme <sequence of enzyme to use =GCTCTTC>] [--bam <.bam alignment of NGS reads/contigs to input ref FASTA>] [--ngsBuffer <bp +/- stitchPosition to require NGS alignment =500>] [--threshold <minimum score threshold below which to cut stitched maps =1>] [--maxlab <max_label_gap_tolerence =1>] [--maxfill <max basepairs to fill between contigs =30000>] [--wobble <fragile site wobble in bp =30000>] [--seq <sequence bp +/- fragile site to print in final BED =off>] [--force <overwrite output folder =off>] [--n <number of CPU cores =nproc>] [--j <number of parallel jobs =nproc/6>] [--optArgs <optArguments.xml =optArguments_human.xml>] [--runSV <enable run SV detection =off>]  [--break <break maps at stitch locations with score below threshold =off>] [--endoutlier <Pvalue penalty for endoutlier =1e-3>] [--aggressive <calculate TypeIII and TypeIV fragile sites =off>]
```
``` 
