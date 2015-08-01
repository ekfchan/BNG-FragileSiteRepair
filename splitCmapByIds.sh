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
host=$(hostname)
phihost="phihost"

echo "HOSTNAME: "$host""
echo ""

# Create directory <dir> for storing results
# Exit with error if <dir> already exists
if [ -d "$cmapDir" ]; then
 	echo "$cmapDir already exists and will overwritten"
 	rm -rf $cmapDir
	mkdir $cmapDir
else 
	mkdir $cmapDir
fi
echo "Splitting "$1" into individual maps..."

# Get list of ids to merge
cat $1 | cut -f 1 | uniq | sed "/#/d" | sort > "$cmapDir"/cmapIdsToSplit.txt

# Extract each contig cmap using RefAligner -merge, via qsub
while read i; do
	# Create qsub script for each contig
	echo '#!/usr/bin/bash' > "$cmapDir"/contig"$i".sh
	echo 'export SGE_ROOT=/opt/sge' >> "$cmapDir"/contig"$i".sh
	echo 'export O64_OMP_SET_AFFINITY=false' >> "$cmapDir"/contig"$i".sh
	echo "~/tools/RefAligner -merge -maxthreads 2 -selectid "$i" -i "$mergedCmap" -o "$cmapDir"/"$cmapPrefix"_contig"$i"" >> "$cmapDir"/contig"$i".sh
	# qsub on Xeon Phi; note the specified queue phihost
	if [[ $host =~ $phihost ]]; then
		qsub -e /dev/null -o /dev/null -pe openmp 2 -q phihost -cwd -N contig"$i" "$cmapDir"/contig"$i".sh > /dev/null
	else
		/bin/bash "$cmapDir"/contig"$i".sh > /dev/null &
	fi
	rm "$cmapDir"/contig"$i".sh
done < "$cmapDir"/cmapIdsToSplit.txt 
echo "All jobs submitted"