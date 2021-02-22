#!/bin/bash
set -e

Query=$1
N_threads=$2
DIR=$3
Simcutoff=$4
Lencuoff=$5
out=${Query}_sarg

if [ ! -d $out ]; then
        mkdir $out;
else
        echo "Warning: $out already exists. previous results are overwrited"
		rm -rf $out
		mkdir -p $out
fi

##########################
echo "
searching $Query agaisnt SARG-nt with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"


${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -P $N_threads -f BlastTab ${DIR}/database/SARG_20170328_5020.ffn $Query > /tmp/argpore_${nowt}_${Query}_tmp.blast

echo "parsing SARG-nt last alignment"
grep -v "#" /tmp/argpore_${nowt}_${Query}_tmp.blast > /tmp/argpore_${nowt}_${Query}_tmp.blast.modified

${DIR}/bin/BlastTab.addlen.sh \
		${DIR}/database/SARG_20170328_5020.ffn \
		${Query} \
		/tmp/argpore_${nowt}_${Query}_tmp.blast.modified \
		$DIR \
		${out}/${Query}_sarg.last


############################
echo "
searching $Query agaisnt ESCG database with similarity cutoff $Simcutoff and alignment length cutoff $Lencuoff using $N_threads threads"


${DIR}/bin/last-983/src/lastal -s 2 -T 0 -Q 0 -a 1 -P $N_threads -f BlastTab ${DIR}/database/ESCG.fna $Query > /tmp/argpore_${nowt}_${Query}_tmp.blast

echo "parsing ESCG last alignment"
grep -v "#" /tmp/argpore_${nowt}_${Query}_tmp.blast > /tmp/argpore_${nowt}_${Query}_tmp.blast.modified

${DIR}/bin/BlastTab.addlen.sh \
		${DIR}/database/ESCG.fna \
		${Query} \
		/tmp/argpore_${nowt}_${Query}_tmp.blast.modified \
		$DIR \
		${out}/${Query}_escg.last


echo "parsing ESCG last alignment"