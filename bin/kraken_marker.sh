#!/bin/bash
set -e

Query=$1
N_threads=$2
DIR=$3
Simcutoff=$4
Lencuoff=$5
nowt=$6

##### taxonomy annotation based on KRAKEN  #################
Query=$1
out=${Query}_KRAKEN
KRAKEN_DB=`grep "KRAKEN_DB" ${DIR}/ARGpore_CONFIG | head -1 | sed 's/KRAKEN_DB=//' | sed 's/"//g'`

echo "
Start KRAKEN @ `date +"%Y-%m-%d %T"`"

if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists. previous results are overwrited"
	rm -rf $out
	mkdir -p $out

fi
echo "
KRAKEN annotating"
${DIR}/bin/kraken/kraken --db ${KRAKEN_DB} --fasta-input --threads $N_threads --output ${out}/combined.krak $Query

echo "
Finish KRAKEN @ `date +"%Y-%m-%d %T"`
"

###############################################################################
# taxonomy based on lastal 1D.fa or 2D.fa agaisnt MetaPhlan markers.fasta
################################################################################
Query=$1
out=${Query}_marker

if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists. previous results are overwrited"
	rm -rf $out
	mkdir -p $out

fi
echo "
Start markergene @ `date +"%Y-%m-%d %T"`"
echo "
Last $Query agaisnt MetaPhlan 2.0 markergene database with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"

${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -b 1 -q 2 -P $N_threads -f BlastTab ${DIR}/database/markers.lastindex $Query > /tmp/argpore_${nowt}_${Query}_tmp.blast3

echo "finish last agaisnt markergene database"
echo "parsing markergene results"
grep -v "#" /tmp/argpore_${nowt}_${Query}_tmp.blast3 > /tmp/argpore_${nowt}_${Query}_tmp.blast3.modified

ruby ${DIR}/bin/BlastTab.addlen.rb -s -f ${DIR}/database/markers.fasta < /tmp/argpore_${nowt}_${Query}_tmp.blast3.modified > /tmp/argpore_${nowt}_${Query}_tmp.blast4

ruby ${DIR}/bin/BlastTab.addlen.rb -f $Query < /tmp/argpore_${nowt}_${Query}_tmp.blast4 > ${out}/${Query}_marker.last

echo "
Finish markergene @ `date +"%Y-%m-%d %T"`"

## merge taxonomy annotation by KRAKEN， taxator-tk and MetaPhlan2.0 markergene ########
echo "
merging KRAKEN,taxator-tk and markergene annoation in R"
Query=$1
Simcutoff=$4
Lencuoff=$5

# out1=${Query}_taxator-tk
out2=${Query}_KRAKEN
out3=${Query}_marker

# Rscript ${DIR}/bin/combine.kraken.taxator.R ${DIR}/database/2020-06-16_lineage.tab ${out2}/combined.krak ${out1}/contig_taxonomy.tab.modified ${out3}/${Query}_marker.last $DIR/database/taxa.info.RData $Simcutoff $Lencuoff  ${Query}_taxa.tab

Rscript ${DIR}/bin/combine.kraken.marker.R \
	${DIR}/database/2020-06-16_lineage.tab \
	${out2}/combined.krak \
	${out3}/${Query}_marker.last \
	$DIR/database/taxa.info.RData \
	$Simcutoff \
	$Lencuoff  \
	${Query}_taxa.tab