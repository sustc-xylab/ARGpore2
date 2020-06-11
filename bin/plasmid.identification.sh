#!/bin/bash
set -e

# bash ${DIR}/scr/plasmid.identification.sh ${Input_fa2} ${DIR} $N_threads $Query $gff
Query=${1}
DIR=${2}
N_threads=$3
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
echo "predict plasmid using plasflow with cutoff 0.95"
# need this line to use conda in a script
source ${CONDA_BASE}/etc/profile.d/conda.sh

conda activate plasflow

Query=$1
out=${Query}_Plasmid
plasflow_cutoff="0.95"
PATH=/home/sustc-xy/miniconda2/envs/plasflow/bin:$PATH

PlasFlow.py --input ${Query} --output ${out}/${Query}_plasflow --threshold ${plasflow_cutoff}

grep ">" ${out}/${Query}_plasflow_plasmids.fasta | sed 's/>//' | cut -d " " -f 1 > ${out}/${Query}_plasflow.tab

conda deactivate


if [ -s ${out}/${Query}_plasflow.tab ]; then
	###### Last agaisnt refseq_plasmid database############################
	echo "
	further filter plasmids by last agaisnt PLSDB database"
	Query=$1
	Simcutoff=0.8
	
	echo "last agaisnt PLSDB database using similarity cutoff $Simcutoff"
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

