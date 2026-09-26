# ARGpore2

**ARGs identification** from nanopore 1D/2D reads

ARGpore2 is a easy-to-use bioinformatics pipeline which codifies the current best practice to identify antibiotic resistance genes (ARGs) and its host populations from nanopore reads (longer than 1kb in fasta format).

Please read below instructions carefully to avoid unnecessary errors.

## Installation

Program installation and reference-data preparation are separate. Reference data stay on their original sites (git LFS, Figshare, S3). This repository does not republish those files. The versions used here are those of the 2022 analysis (PLSDB `2020_03_04`, SARG `20211207`, Centrifuge `p_compressed+h+v`). A failed download is not replaced by a later release.

After installation the command is `argpore2` (the shell script `argpore.sh` remains the analysis engine).

R cannot be solved into the same conda environment as BioPerl and centrifuge, so two environments are required.

```bash
git clone https://github.com/sustc-xylab/ARGpore2.git
cd ARGpore2
conda env create -f environment-tools.yml
conda env create -f environment-r.yml
conda activate argpore2
pip install -e .
argpore2 setup-tools
argpore2 setup-db
argpore2 check
```

`setup-tools` and `setup-db` look for conda environments `argpore2` and `argpore2-r` (`ARGPORE_CONDA_ENV`, `ARGPORE_R_ENV`). `setup-db` can be repeated. Source archives are kept. HTTP 403/404 or a truncated file is an error: the message names the file and the path under `database/` where a manual copy should be placed. Already completed files are left in place.

| Item | Channel | Version used |
| --- | --- | --- |
| MetaPhlAn2 markers, SARG, ESCG, lineage | git LFS / files in `database/` | as archived in this repository |
| PLSDB | Figshare file `21961095` | `2020_03_04` |
| Centrifuge index | `genome-idx.s3.amazonaws.com` | `p_compressed+h+v` |

If Figshare returns 403, download version `2020_03_04` from the [PLSDB previous-versions list](https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download) and place `PLSDB.zip` or `PLSDB_2020_03_04.fna` in `database/`, then run `argpore2 setup-db` again.

`setup.sh` still mixes program downloads and reference data. Use `argpore2 setup-tools` and `argpore2 setup-db` instead.

## Using ARGpore2

```bash
argpore2 -f test.fa -t 60 -o demo
```

Accepted FASTA suffixes are `.fa`, `.fasta`, `.fas`, `.fna`, and `.ffn`. `-o` receives a hard link to the input when the filesystem allows it, and a copy otherwise. Products are written only under `-o`. Reads longer than 1 kb are recommended.

`bash argpore.sh -f test.fa -t 60 -o demo` calls the same analysis.


	
#### Output files

With `-o`, results are written under that directory, in a folder named `$INPUT_FASTA_ARGpore2_<timestamp>`.

Main output files include:
	
	input_arg.tab		ARG quntification (No. of reads)
	input_arg.w.taxa.tab	ARGs-containing nanopore reads with taxonomy assignment and plausible plasmid identification
	input_plasmid.like.tab	plasmid-like nanopore reads identified
	input_taxa.tab		taxonomy assignment of all nanopore reads

plasmid-like nanopore reads are identified by last query against PLSDB (only hit showing alignment with > 70% similarity over 70% of its lenth to a known plasmid in PLSDB is considered as valid plasmid hit). **NOTICE**: This method cannot fully distinguish plasmids from chromosome, as a result, it only reports plasmid-like nanopore reads in **input_plasmid.like.tab**. If such a plasmid-like nanopre read also showed circular nature,it is more likely to be a real plasmid. You may use ccontigs (https://github.com/Microbiology/ccontigs.git) to check the circular nature of nanopore reads. Althernatively, users may want to use Plasflow (https://github.com/smaegol/PlasFlow) to double-confirm these plasmid-like nanopore reads by comparing their kmer freqeucy to that of known plasmids. 

Taxonomy annotation of nanopore reads were derived by combining results of Centrifuge and MetaPhlan2 markergene database. If case of inconsistent annotations among these tools, to maximize classification ratio, ARGpore2 combines results with priority as Centrifuge > markergene. 

## *Citation:*

If you use ARGpore2 in your nanopore dataset analysis please cite:

Ziqi Wu, You Che, Chenyuan Dang, Miao Zhang, Xuyang Zhang, Yuhong Sun, Xiang Li, Tong Zhang, Yu Xia*. 2022. Nanopore-based Long-read Metagenomics Uncover the Resistome Intrusion by Antibiotic Resistant Bacteria from Treated Wastewater in Receiving Water Body. Water Research. 226:119282. (https://doi.org/10.1016/j.watres.2022.119282)


##### Tools included in ARGpore2 should be also cited, these tools includes: 

last, blast+, Centrifuge, MetaPhlan2, PLSDB, GNU parallel, R, python


