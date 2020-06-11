# ARGpore2
Author: Yu XIA 

Email: shuixia100@gmail.com

version 2.0

ARGpore2 is designed to do **ARG identification of nanopore 1D/2D reads (fasta format)** and then carried out **taxonomy and plasmid annotation** of ARG-containing long reads.

NOTE: You may still use ARGpore2 for your metagenomic-assembled contig/scaffold.

Please read this manual carefully to avoid unnecessary errors for ARGpore implementation.

## Installation 
### Pre-requisites for ARGpore 
	ruby 2.3.1p112
	python2.7
	GNU parallel
	conda
	R and library: plyr, data.table

### Setup ARGpore2
	
	git clone https://github.com/sustc-xylab/ARGpore2.git
	
	cd ARGpore2
	
	bash setup.sh


The ARGpore_CONFIG contains the PATH for database required for ARGpore, this file should always be stored in the same directory with ARGpore.sh.
 
Before runing ARGpore, users should modify ARGpore_CONFIG with their specific database PATH


### Database to prepare before ARGpore run 

Download comprehensive database usually takes quite long time. please stay patient :)
 

	** NCBI nt database for BLAST+ and nodes.dmp of corresponding taxonomy ** 
		
		1. download NCBI preformatted nt database with:
			
			perl  $PATH_to_ARGpore/bin/ncbi-blast-2.9.0+/bin/update_blastdb.pl --passive --decompress nt
			
		specify the name of the nt database in ARGpore_CONFIG
		
		2. download NCBI taxonomy with：
			
			wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		
		Place the files names.dmp and nodes.dmp in a folder and specify its path in ARGpore_CONFIG. keep in mind that the taxonomy files are modified on a regular basis.  
		
	** Default database of KRAKEN **
		
		To create the standard Kraken database, you can use the following command:
		
			$PATH_to_ARGpore/bin/kraken/kraken-build --standard --threads 24 --db $DBNAME
		
		Replace "$DBNAME" with your preferred database name/location. 
		
		specify the path of the KRAKEN database in ARGpore_CONFIG
	
## Using ARGpore 
Once download ARGpore package, all needed analysis is wrapped up in one executable named **ARGpore.sh**. Please **use bash instead of sh** to initiate ARGpore.

	bash $PATH_to_ARGpore/argpore.sh -f test.fa -t 60 > ARGpore.log


	
#### Output files 
All output files of ARGpore are stored in a directory named $INPUT_FASTA_ARGpore_nowtime, main output files include:
	
	input_taxa.tab	phylogenetic assignment of nanopore reads by combining taxator-tk, KRAKEN and MetaPhlan2 marker gene 
	input_arg.w.taxa.tab	ARGs-containing nanopore reads with valid taxonomy assignment and plausible plasmid identification
	input_plasmid.like	plasible plasmids identified by combining Plasflow and last against PLSDB database
	input_circular.tsv	circular nanopore reads identified by ccontigs
	

## *Citation:*

If you use ARGpore2 in your nanopore dataset analysis please cite:

Xia, Y., Li, A.-D., Deng, Y., Jiang, X.-T., Li, L.-G., and Zhang, T. (**2017**) MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. *Front Microbiol* 8: 2105.

##### Tools included in ARGpore2 should be also cited, these tools includes: 

last, blast+, taxator-tk, kraken, MetaPhlan2, ccontigs, plasflow, conda, GNU parallel, ruby, R, python




