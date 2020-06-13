# ARGpore2

**ARG identification** from nanopore 1D/2D reads (fasta format) and then carried out **taxonomy and plasmid annotation** of nanopore reads.

Please read below instructions carefully to avoid unnecessary errors.

## Installation 
### Pre-requisites for ARGpore 
	ruby 2.3.1p112
	python2.7
	GNU parallel
	conda
	R and library: plyr, data.table, doParallel,foreach

### Setup ARGpore2
	
	git clone https://github.com/sustc-xylab/ARGpore2.git
	
	cd ARGpore2
	
	bash setup.sh	

The setup.sh will install ccontigs, blast+,  plasflow, kraken and then download SARG-nt,ESCG, MetaPhlan2 Markergene and PLSDB database for you. It will take quite some time to finish, please stay patient :)


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
		
The ARGpore_CONFIG contains the PATH for database required for ARGpore, this file should always be stored in the same directory with ARGpore.sh. Before runing ARGpore, users should **modify ARGpore_CONFIG with their specific database PATH**

## Using ARGpore2 
Once download ARGpore2 package, all needed analysis is wrapped up in one executable named **ARGpore.sh**. Please **use bash instead of sh** to initiate argpore.sh.

	bash $PATH_to_ARGpore2/argpore.sh -f test.fa -t 60 > ARGpore.log


	
#### Output files 
All output files of ARGpore are stored in a directory named $INPUT_FASTA_ARGpore_nowtime, main output files include:
	
	input_arg.tab	ARG quntification (copy per cell)
	input_arg.w.taxa.tab	ARGs-containing nanopore reads with taxonomy assignment and plausible plasmid identification
	input_circular.tab	circular nanopore reads identified
	input_plasmid.like.tab	plasmid-like nanopore reads identified
	input_taxa.tab	taxonomy assignment of all nanopore reads

plasmid-like nanopore reads are identified by firstly using plasflow to identify plasmids (probability threshold of 0.95), then plasflow-plasmids are further filtered by last against PLSDB (only hit showing alignment with > 80% similarity over 70% of its lenth to a known plasmid in PLSDB is considered as valid hit). **NOTICE**: This method cannot fully distinguish plasmids from chromosome, as a result, it only report plasmid-like nanopore reads in **input_plasmid.like.tab**. If such a plasmid-like nanopre read also showed circular nature as indicated in **Input_circular.tab**, it is more likely to be a real plasmid. 

Taxonomy annotation of nanopore reads were derived by combining results of taxator-tk, KRAKEN and MetaPhlan2 markergene database. If case of inconsistent annotations among these tools, to maximize classification ratio, ARGpore2 combines results with priority as kraken > taxator-tk > markergene. 

## *Citation:*

If you use ARGpore2 in your nanopore dataset analysis please cite:

Xia, Y., Li, A.-D., Deng, Y., Jiang, X.-T., Li, L.-G., and Zhang, T. (**2017**) MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. *Front Microbiol* 8: 2105.

##### Tools included in ARGpore2 should be also cited, these tools includes: 

last, blast+, taxator-tk, kraken, MetaPhlan2, ccontigs, plasflow, conda, GNU parallel, ruby, R, python




