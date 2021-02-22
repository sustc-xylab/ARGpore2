#!/bin/bash
set -e

Query=$1
N_threads=$2
DIR=$3
Simcutoff=$4
Lencuoff=$5
nowt=$6

##### taxonomy annotation based on Centrifuge #################
Query=$1
out=${Query}_Centrifuge
CEN_DB="${DIR}/database/p_compressed+h+v"

echo "
Start Centrifuge @ `date +"%Y-%m-%d %T"`"

if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists. previous results are overwrited"
	rm -rf $out
	mkdir -p $out

fi
echo "
Centrifuge annotating"
${DIR}/bin/centrifuge/centrifuge -x ${CEN_DB} -f ${Query} -p ${N_threads} -S ${out}/${Query}_centrifuge.report.txt --report-file ${out}/${Query}_centrifuge.report.summary.tsv

echo "
Finish Centrifuge @ `date +"%Y-%m-%d %T"`
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

${DIR}/bin/BlastTab.addlen.sh \
		${DIR}/database/markers.fasta \
		${Query} \
		/tmp/argpore_${nowt}_${Query}_tmp.blast3.modified \
		$DIR \
		${out}/${Query}_marker.last

echo "
Finish markergene @ `date +"%Y-%m-%d %T"`"

## merge taxonomy annotation by KRAKEN and MetaPhlan2.0 markergene ########
echo "
merging Centrifuge and markergene annoation in R"
Query=$1
Simcutoff=$4
Lencuoff=$5

out2=${Query}_Centrifuge
out3=${Query}_marker

Rscript ${DIR}/bin/combine.centrifuge.marker.R \
	${DIR}/database/2020-06-16_lineage.tab \
	${out2}/${Query}_centrifuge.report.txt \
	${N_threads} \
	${out3}/${Query}_marker.last \
	$DIR/database/taxa.info.RData \
	$Simcutoff \
	$Lencuoff  \
	${Query}_taxa.tab