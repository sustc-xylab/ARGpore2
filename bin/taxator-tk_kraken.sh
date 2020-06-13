#!/bin/bash
set -e

# taxonomy assignment of contig by taxator-tk ########################
# taxator-tk is homology based on blastn search agaisnt nt ####
Query=$1
N_threads=$2
DIR=$3
Simcutoff=$4
Lencuoff=$5
nowt=$6

out=${Query}_taxator-tk

BLASTDB=`grep "BLASTDB" ${DIR}/ARGpore_CONFIG | head -1 | sed 's/BLASTDB=//'`
TAXDUMP=`grep "TAXDUMP" ${DIR}/ARGpore_CONFIG | head -1 | sed 's/TAXDUMP=//' | sed 's/"//g'`
TAXDUMP2=`echo $TAXDUMP | sed 's/\/nodes.dmp//'`


echo "
Start taxator-tk @ `date +"%Y-%m-%d %T"`"
export TAXATORTK_TAXONOMY_NCBI=$TAXDUMP2

if [ ! -d $out ]; then
        mkdir $out;
else
        echo "Warning: $out already exists. previous results are overwrited"
		rm -rf $out
		mkdir -p $out
fi


echo "
aligning $Query to ${BLASTDB} database with MEGABLAST. This is the slowest step, please stay patient.We will check the progress of megablast every one minute for you:" 

rm -f tmp.megablast.jobs
$DIR/bin/ncbi-blast-2.9.0+/bin/blastn -task megablast -num_threads $N_threads\
 -db ${BLASTDB}\
 -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length'\
 -query ${Query} > ${out}/megablast_out.raw.tab &
PID=$!
echo "$PID" >> tmp.megablast.jobs

echo ""
touch ${out}/megablast_out.raw.tab
bash $DIR/bin/blastab_monitor.sh ${out}/megablast_out.raw.tab $Query tmp.megablast.jobs
rm -f tmp.megablast.jobs

echo "removing unnecessary lines that lead to bad tax IDs (without a proper rank)"

python2.7 ${DIR}/bin/taxator-tk/prune_blast_hits.py ${TAXDUMP} ${out}/megablast_out.raw.tab > ${out}/megablast_out.pruned.tab

cat ${out}/megablast_out.pruned.tab | cut -f1,2,3,4,5,7,8,9,10,11,12 > ${out}/megablast_out.tab

echo "making mapping file"
cat ${out}/megablast_out.pruned.tab | cut -f 5,6 > ${out}/mapping.tax


echo "pulling out classifications with taxator with megan-lca algorithm"
cat ${out}/megablast_out.tab | ${DIR}/bin/taxator-tk/taxator -a megan-lca -t 0.3 -e 0.00001 -g ${out}/mapping.tax -p $N_threads > ${out}/predictions.gff3

echo "binning and consolidating classifications for each contig"
sort -k1,1 ${out}/predictions.gff3 | ${DIR}/bin/taxator-tk/binner -n classification -i genus:0.6 > ${out}/binned_predictions.txt
mv binning.log ${out}

echo "pulling out full taxonomy path with taxknife"
cat ${out}/binned_predictions.txt | ${DIR}/bin/taxator-tk/taxknife -f 2 --mode annotate -s path | grep -v "Could not" | cut -f1,2 > ${out}/contig_taxonomy.tab

# modify for merging with KRAKEN
cat ${out}/contig_taxonomy.tab| sed 's/;/\t/g' > ${out}/contig_taxonomy.tab.modified

echo "
Finish taxator-tk @ `date +"%Y-%m-%d %T"`"


# ##### taxonomy annotation based on KRAKEN  #################
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
Last $Query agaisnt MetaPhlan 2.0 markergene database with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"

${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -b 1 -q 2 -P $N_threads -f BlastTab ${DIR}/database/markers.lastindex $Query > /tmp/argpore_${nowt}_${Query}_tmp.blast3

echo "finish last agaisnt markergene database"
echo "parsing markergene results"
grep -v "#" /tmp/argpore_${nowt}_${Query}_tmp.blast3 > /tmp/argpore_${nowt}_${Query}_tmp.blast3.modified

ruby ${DIR}/bin/BlastTab.addlen.rb -s -f ${DIR}/database/markers.fasta < /tmp/argpore_${nowt}_${Query}_tmp.blast3.modified > /tmp/argpore_${nowt}_${Query}_tmp.blast4

ruby ${DIR}/bin/BlastTab.addlen.rb -f $Query < /tmp/argpore_${nowt}_${Query}_tmp.blast4 > ${out}/${Query}_marker.last


## merge taxonomy annotation by KRAKEN， taxator-tk and MetaPhlan2.0 markergene ########
echo "
merging KRAKEN,taxator-tk and markergene annoation in R"
Query=$1
Simcutoff=$4
Lencuoff=$5

out1=${Query}_taxator-tk
out2=${Query}_KRAKEN
out3=${Query}_marker

Rscript ${DIR}/bin/combine.kraken.taxator.R ${DIR}/database/lineages-2019-02-20.csv ${out2}/combined.krak ${out1}/contig_taxonomy.tab.modified ${out3}/${Query}_marker.last $DIR/database/taxa.info.RData $Simcutoff $Lencuoff  ${Query}_taxa.tab
