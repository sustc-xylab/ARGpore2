#!/bin/bash
set -e

# bash ${DIR}/scr/plasmid.identification.sh ${Input_fa2} ${DIR} $N_threads $Query $gff
Query=${1}
DIR=${2}
N_threads=$3
nowt=$5

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


###### Last agaisnt PLSDB database############################
echo "
Identify plasmids by last agaisnt PLSDB database"
Query=$1
Simcutoff=0.7
	
echo "last agaisnt PLSDB database using similarity cutoff $Simcutoff length cutoff 0.7"
${DIR}/bin/last-983/scripts/parallel-fasta "${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -P ${N_threads} -f BlastTab ${Plasdb}" < ${Query} > ${out}/${Query}_last.plasmid.tab 
	
if [ -s ${out}/${Query}_last.plasmid.tab ]; then

	grep -v "#" ${out}/${Query}_last.plasmid.tab  > /tmp/${nowt}_${Query}_tmp.modified

	${DIR}/bin/BlastTab.addlen.sh \
		${Plasdb_fasta} \
		${Query} \
		/tmp/${nowt}_${Query}_tmp.modified \
		$DIR \
		${out}/${Query}_last.plasmid.tab_wlength
		
	
	###### merge last and plasflow results in R
	# echo "merge plasflow and PLSDB result"
		Rscript ${DIR}/bin/parse.plasmid.last.R \
		${out}/${Query}_last.plasmid.tab_wlength \
		$Plasdb_fastaname \
		$Simcutoff \
		${Query}_plasmid.like.tab \
		${out}/${Query}_plasmid.like_last.hit \
		$N_threads
	
	
else
	echo "
	No Plasmid identified"
fi

rm -f ${out}/${Query}_last.plasmid.tab_wlength # to save disk space