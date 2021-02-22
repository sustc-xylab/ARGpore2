#!/bin/bash
set -e

db_fasta=$1
Query=$2
blast_tab=$3
DIR=$4
output=$5

cut -f 2 ${blast_tab}| sort -u > tmp_${Query}
fgrep -f tmp_${Query} ${db_fasta}.length > tmp_${Query}_db.length



$DIR/bin/fastaNameLengh.pl $Query > ${Query}.length
grep -v "#" $blast_tab | cut -f 1 | sort -u > tmp_${Query}_2
fgrep -f tmp_${Query}_2  ${Query}.length > tmp_${Query}_query.length

Rscript $DIR/bin/merge.length.R \
	${blast_tab} \
	tmp_${Query}_query.length \
	tmp_${Query}_db.length \
	${output}

rm -f tmp_${Query}
rm -f tmp_${Query}_2
rm -f tmp_${Query}_db.length
rm -f tmp_${Query}_query.length
rm -f ${Query}.length
