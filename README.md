# ARGpore2

**ARGs identification** from nanopore 1D/2D reads

ARGpore2 is a easy-to-use bioinformatics pipeline which codifies the current best practice to identify antibiotic resistance genes (ARGs) and its host populations from nanopore reads (longer than 1kb in fasta format).

Please read below instructions carefully to avoid unnecessary errors.

## Installation 
### Pre-requisites for ARGpore 
	
	python2.7	### sudo apt install python2.7
	GNU parallel	### sudo apt install parallel
	git lfs	        ### sudo apt install git-lfs
	R and library: plyr, data.table, doParallel, foreach 
	

### Setup ARGpore2
	
	git clone https://github.com/sustc-xylab/ARGpore2.git
	
	cd ARGpore2
	
	bash ./setup.sh	

The setup.sh will install **blast+, Centrifuge** and then download **bacteria+archaea+virus database for Centrifuge, MetaPhlan2 Markergene and PLSDB database** for you. It will take at least **4 hour** to finish, please stay patient :)



## Using ARGpore2 
Once installed ARGpore2 package, all needed analysis is wrapped up in one executable named **argpore.sh**. Please **use bash instead of sh** to initiate argpore.sh.

**NOTICE:**
	To avoid cross-writing of intermediate files, each ARGpore run should has a independent working directory. To improve annotation accuracy, your input fasta should be **longer than 1kb** 
	
	

	mkdir -p demo
	cp test.fa demo 
	cd demo 
	bash $PATH_to_ARGpore2/argpore.sh -f test.fa -t 60 > ARGpore.log


	
#### Output files 
All output files of ARGpore are stored in a folder named $INPUT_FASTA_ARGpore_nowtime in the working directory.

Main output files include:
	
	input_arg.tab		ARG quntification (copy per cell)
	input_arg.w.taxa.tab	ARGs-containing nanopore reads with taxonomy assignment and plausible plasmid identification
	input_plasmid.like.tab	plasmid-like nanopore reads identified
	input_taxa.tab		taxonomy assignment of all nanopore reads

plasmid-like nanopore reads are identified by last query against PLSDB (only hit showing alignment with > 70% similarity over 70% of its lenth to a known plasmid in PLSDB is considered as valid plasmid hit). **NOTICE**: This method cannot fully distinguish plasmids from chromosome, as a result, it only reports plasmid-like nanopore reads in **input_plasmid.like.tab**. If such a plasmid-like nanopre read also showed circular nature,it is more likely to be a real plasmid. You may use ccontigs (https://github.com/Microbiology/ccontigs.git) to check the circular nature of nanopore reads. Althernatively, users may want to use Plasflow (https://github.com/smaegol/PlasFlow) to double-confirm these plasmid-like nanopore reads by comparing their kmer freqeucy to that of known plasmids. 

Taxonomy annotation of nanopore reads were derived by combining results of Centrifuge and MetaPhlan2 markergene database. If case of inconsistent annotations among these tools, to maximize classification ratio, ARGpore2 combines results with priority as Centrifuge > markergene. 

## *Citation:*

If you use ARGpore2 in your nanopore dataset analysis please cite:

Xia, Y., Li, A.-D., Deng, Y., Jiang, X.-T., Li, L.-G., and Zhang, T. (**2017**) MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. *Front Microbiol* 8: 2105.

##### Tools included in ARGpore2 should be also cited, these tools includes: 

last, blast+, Centrifuge, MetaPhlan2, PLSDB, GNU parallel, R, python


