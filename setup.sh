SCRIPT=`realpath $0`
DIR=`dirname $SCRIPT`

cd ${DIR}/bin

# seqkit
wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_amd64.tar.gz
tar -zxvf seqkit_linux_amd64.tar.gz
rm seqkit_linux_amd64.tar.gz

# blast+ 2.9.0
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
tar -zvxf ncbi-blast-2.9.0+-x64-linux.tar.gz
rm -f ncbi-blast-2.9.0+-x64-linux.tar.gz

# kraken
git clone https://github.com/DerrickWood/kraken.git
cd kraken
./install_kraken.sh ./
cd ..

# plasflow 
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name plasflow python=3.5
source activate plasflow
conda install -c jjhelmus tensorflow=0.10.0rc0
conda install plasflow -c smaegol
source deactivate

# ccontigs #####
git clone https://github.com/Microbiology/ccontigs.git
wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz
tar -xvzf julia-1.4.2-linux-x86_64.tar.gz
rm -f julia-1.4.2-linux-x86_64.tar.gz


cd ${DIR}/database
# MetaPhlan2.0 markergene database ##############################
echo "
Downloading Metaphlan2 markergene database from git lfs"

git lfs install
git lfs pull

tar jxvf markers.fasta.tar.xz
echo "
Building lastdb for Metaphlan2 markergene database"
${DIR}/bin/last-983/src/lastdb -Q 0 markers.lastindex markers.fasta -P 10

rm -f markers.fasta.tar.xz

###### SARG-nt database ################
echo "
Building lastdb for SARG-nt database"
tar jxvf SARG_20170328_5020.ffn.tar.xz
${DIR}/bin/last-983/src/lastdb -Q 0 SARG_20170328_5020.ffn SARG_20170328_5020.ffn -P 10

rm -f SARG_20170328_5020.ffn.tar.xz

########### ESCG database ##########################################
echo "
Building lastdb for ESCG database"
tar jxvf ESCG.fna.tar.xz
${DIR}/bin/last-983/src/lastdb -Q 0 ESCG.fna ESCG.fna -P 10
rm -f ESCG.fna.tar.xz


########## lineage database ###################################################################
echo "
Downloading lineage information for NCBI taxonomy"
wget https://gitlab.com/zyxue/ncbitax2lin-lineages/blob/master/lineages-2019-02-20.csv.gz
	
gunzip lineages-2019-02-20.csv.gz

######### PLSDB database #################################
echo "
Downloading PLSDB"
wget https://ndownloader.figshare.com/files/21961095 --output-document 'PLSDB.zip'
unzip PLSDB.zip
mv data/pls/PLSDB_2020_03_04 .

$DIR/bin/ncbi-blast-2.9.0+/bin/blastdbcmd -db PLSDB_2020_03_04/plsdb.fna -dbtype nucl -entry all -outfmt "%f" -out PLSDB_2020_03_04.fna

grep ">" PLSDB_2020_03_04.fna | sed 's/>//' | sed 's/\ /\t/' > PLSDB_2020_03_04.fna.name

rm -rf data
rm PLSDB.zip


# make last database using PLSDB database
echo "
Building lastdb for PLSDB, this step is quite slow, please stay patient :)"
${DIR}/bin/last-983/src/lastdb -Q 0 PLSDB_2020_03_04.fna.lastindex PLSDB_2020_03_04.fna -P 50


