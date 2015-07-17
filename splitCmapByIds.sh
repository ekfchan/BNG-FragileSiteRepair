#!/usr/bin/bash

# Script to split merged (by chromosome) fragile site repaired cmap into individual query contigs. 

# Usage: splitCmapByIds.sh <cmap> <dir> <prefix>
if [ $# -ne 3 ]; then 
    echo "Usage: splitCmapByIds.sh <cmap> <dir> <prefix>"
    exit 1
fi
mergedCmap=$1
cmapDir=$2
cmapPrefix=$3

# Create directory <dir> for storing results
# Exit with error if <dir> already exists
if [ -d "$cmapDir" ]; then
 	echo "$cmapDir already exists"
 	exit 1
else 
	mkdir $cmapDir
fi

# Get list of ids to merge
cat $1 | cut -f 1 | uniq | sed "/#/d" | sort > "$cmapDir"/cmapIdsToSplit.txt

# Extract each contig cmap using RefAligner -merge, via qsub
while read i; do
	# Create qsub script for each contig
	echo '#!/usr/bin/bash' > "$cmapDir"/contig"$i".sh
	echo 'export SGE_ROOT=/opt/sge' >> "$cmapDir"/contig"$i".sh
	echo "~/tools/RefAligner -merge -selectid "$i" -i "$mergedCmap" -o "$cmapDir"/"$cmapPrefix"_contig"$i"" >> "$cmapDir"/contig"$i".sh
	# qsub on Xeon Phi; note the specified queue phihost
	qsub -e /dev/null -o /dev/null -q phihost -cwd -N contig"$i" "$cmapDir"/contig"$i".sh
	rm "$cmapDir"/contig"$i".sh
done < "$cmapDir"/cmapIdsToSplit.txt 