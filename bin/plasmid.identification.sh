#!/bin/bash
set -e

# bash ${DIR}/scr/plasmid.identification.sh ${Input_fa2} ${DIR} $N_threads $Query $gff
Query=${1}
DIR=${2}
N_threads=$3
nowt=$5
CONDA_BASE=$(conda info --base)

out=${Query}_Plasmid

Plasdb=${DIR}/database/PLSDB_2020_03_04.fna.lastindex
Plasdb_fasta=${DIR}/database/PLSDB_2020_03_04.fna
Plasdb_fastaname=${DIR}/database/PLSDB_2020_03_04.fna.name

if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists. previous results are overwrited"
	rm -rf $out
	mkdir -p $out

fi

########### plasmid identification by Plasflow ####################
echo "predict plasmid using plasflow with probability threshold of 0.95"
# need this line to use conda in a script
source ${CONDA_BASE}/etc/profile.d/conda.sh

conda activate plasflow

Query=$1
out=${Query}_Plasmid
plasflow_cutoff="0.95"
PATH=${CONDA_BASE}/envs/plasflow/bin:$PATH


n=`grep ">" ${Query} | wc -l`
if [ $n -gt 10000 ]; then
	n2=`echo $(($n/10000))`
	echo "split $Query into $n2 pices to run PlasFlow"
	${DIR}/bin/FastA.split.pl ${Query} /tmp/${Query}_${nowt}.split $n2
	
	rm -f tmp.plasflow_${nowt}.jobs
	find /tmp -name "${Query}_${nowt}.split*.fa" | while read line
	do 
		PlasFlow.py --input ${line} --output ${line}_plasflow --threshold ${plasflow_cutoff} &
		PID=$!
		echo "$PID" >> tmp.plasflow_${nowt}.jobs
	done
	
	bash $DIR/bin/monitor.bgpid.sh tmp.plasflow_${nowt}.jobs
	rm -f tmp.plasflow_${nowt}.jobs
	
	cat /tmp/${Query}_${nowt}.split*.fa_plasflow_plasmids.fasta > ${out}/${Query}_plasflow_plasmids.fasta
	cat /tmp/${Query}_${nowt}.split*.fa_plasflow_chromosomes.fasta > ${out}/${Query}_plasflow_chromosomes.fasta
	cat /tmp/${Query}_${nowt}.split*.fa_plasflow_unclassified.fasta > ${out}/${Query}_plasflow_unclassified.fasta
	rm -f /tmp/${Query}_${nowt}.split*
else
	PlasFlow.py --input ${Query} --output ${out}/${Query}_plasflow --threshold ${plasflow_cutoff} --batch_size 2000
fi

grep ">" ${out}/${Query}_plasflow_plasmids.fasta | sed 's/>//' | cut -d " " -f 1 > ${out}/${Query}_plasflow.tab

conda deactivate


if [ -s ${out}/${Query}_plasflow.tab ]; then
	###### Last agaisnt refseq_plasmid database############################
	echo "
	further filter plasmids by last agaisnt PLSDB database"
	Query=$1
	Simcutoff=0.8
	
	echo "last agaisnt PLSDB database using similarity cutoff $Simcutoff length cutoff 0.7"
	${DIR}/bin/last-983/scripts/parallel-fasta "${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -P ${N_threads} -f BlastTab ${Plasdb}" < ${Query} > ${out}/${Query}_last.plasmid.tab 

	grep -v "#" ${out}/${Query}_last.plasmid.tab  > /tmp/${nowt}_${Query}_tmp.modified

	${DIR}/bin/BlastTab.addlen.rb -s -f ${Plasdb_fasta} < /tmp/${nowt}_${Query}_tmp.modified > /tmp/${nowt}_${Query}_tmp.modified2

	${DIR}/bin/BlastTab.addlen.rb -f ${Query} < /tmp/${nowt}_${Query}_tmp.modified2 > ${out}/${Query}_last.plasmid.tab_wlength 

	###### merge last and plasflow results in R
	echo "merge plasflow and PLSDB result"
	Rscript ${DIR}/bin/parse.plasmid.last.R ${out}/${Query}_plasflow.tab ${out}/${Query}_last.plasmid.tab_wlength $Plasdb_fastaname $Simcutoff ${Query}_plasmid.like.tab ${out}/${Query}_plasmid.like_last.hit
	
else
	echo "
	No Plasmid identified"
fi

