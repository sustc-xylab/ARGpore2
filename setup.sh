SCRIPT=`realpath $0`
DIR=`dirname $SCRIPT`

cd bin

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
cd ../database

# MetaPhlan2.0 markergene database ##############################
echo "downloading MetaPhlan 2 markergene database "
wget --header 'Host: public.bn.files.1drv.com' --user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://onedrive.live.com/' --header 'DNT: 1' --header 'Cookie: WLSRDSecAuth=FABSARQL3KgEDBNbW84gMYrDN0fBab7xkQNmAAAEgAAACIDtreS4OWgkEAEnReiBtSdW5c1dTktkpYAMYhNOSIaCJf52VQrio3sR4F5XoLxtKMbFTLfb6xPIqMuuOE5p2/EuMQsMJ4A6FBDSKXYmVM6lf1M/ijKAUSzN3%2Beuh8dd39i5%2BM3qA%2B8jSK4ylHaN2LHL1lpwjZvltYYKcJCi/RK19txt3jKcc7gZCOtgXs5zr714uPXUdPLNDlzct/Mcl2WaRtUNKxcRjVqt9srDxCCgMWC8wCrlr5AuDeapElflkBU0MYiwIfkN6wnnpmNaNwMtUocZN5SG2tV8AniF7zB3esdHTjD9G5Vhfw/DJc9U2iXmBfOOiOJYRs5MRTE6dqw%2Bc6uSFUxTZ9Xid2SPmmtQiOVvBBw9tUJ9axQA1XLONcYfon8P77J5RNUJ/9g8uAQ%3D; ANON=; NAP=' --header 'Upgrade-Insecure-Requests: 1' 'https://public.bn.files.1drv.com/y4mzqkZIN1mC4HB0ADZnApUNVmLXTUsoF7TH-MRwaG3X_OJpEzJReAxyjbDpz3P7KGP-ufyH9pWOzEaqDqH1sdjf08iiOzbKlk-zNg7rFwcJR58uV7-4Q_l-oja5jMwWn1fog2QKWTX44w75OenYKdN4kxFeHbUfil4cvSXyhBTcQydvfit0ij0RGLrHkc3UZSI4Kanp4Adi4V1Vwy5YTDgWjbLqC9RLl9WZHPbazEuGHo?access_token=EwAAA61DBAAUzl/nWKUlBg14ZGcybuC4/OHFdfEAAeXQh5gHi4ymZUa1kxH7weXYfxZSJJ3fiJsHHpiakCtJaMbrQanlyaLNsgA7Xm3Wim%2b4%2bns9ovBzLaGxXHi0yLeL4/UcdoffK7tn6VCRNBlXHky7ccvNFmrYQJxCz/RQl20/NOPt0hPOUGlQR0qlYBmqlYdMnIVZvyRcx9iC1HaTsNcALCppKeG5G1CnFl6UEsy9GD4WCQDMQ5oafXJDNbi8BNKwv%2bmWjLtZwHBSgBaCW9Mt7P786r1jrCdnTwia84r6k15XbLSxjSAw24YeM555dd/G/0mWotDz6nDu0FwHxUbQPfzhs3p9gdK5MxbJvx8lABbl52hyf0VS3bbgPhUDZgAACIP3vSHbjttj0AGp2XUny54cJJplaZbXj3ayXHpjZqmfrJ%2bM/eOo3fv3SUgb4V1CZIAXAgBcU6tngKXw1I8cI223pb66nJKnlHmvQ05Fnp5M2LX24wMJfJDV9gxS/sUZFys0mh4gKtwIA5Pu1S62fLrqp0iWrhO80EeKJTN3zGnYTL0hMpsd6x2Au1ZjIszYH0eLHa4uk79i31NHgLEhxMXA85qNgUovIvQ5YRbAXcL42nNvi/x1mpCa41iLMmbkGgKuCKjctseaBcca9eRb%2b%2bDEyAoc%2bvc6iduFMzO2NUFDSw4U7j0BtFiLQ3bhZUMS8fqUFxWBlh0Yiz/Q2vQdkXQZ%2bW5CSaS1CIDMh%2b0HwAJfVy9i5SEQNWGxVmjI/MEH2ZEwosGwXRwSgu5/bukAyoJdp7fuIV8I3fniHuSY6LXfg%2b4DM9nFd%2bPvj699aBlObveZX7QX7ky2gycV138qVgxfiwa8yyuTddqzeHpCDpE0RsML3ktF3rIbo7XfvNWTH3gIlmyysnCU8ctUrgXdsZnuMigq1q8uriSDqRHDopZALc8jc5s2tn%2b03t22gtqSuYBW4PH%2bo9TnB3zNBEVeRadBi279od84zREZp1NT5Ee9P%2b8X2agiVmTtWAAC' --output-document 'markers.fasta.zip'

unzip markers.fasta.zip
rm -f markers.fasta
${DIR}/bin/last-983/src/lastdb -Q 0 markers.lastindex markers.fasta -P 10

###### SARG-nt database ################
echo "downloading SARG-nt database "
wget --header 'Host: public.bn.files.1drv.com' --user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://onedrive.live.com/' --header 'DNT: 1' --header 'Cookie: WLSRDSecAuth=FABSARQL3KgEDBNbW84gMYrDN0fBab7xkQNmAAAEgAAACIDtreS4OWgkEAEnReiBtSdW5c1dTktkpYAMYhNOSIaCJf52VQrio3sR4F5XoLxtKMbFTLfb6xPIqMuuOE5p2/EuMQsMJ4A6FBDSKXYmVM6lf1M/ijKAUSzN3%2Beuh8dd39i5%2BM3qA%2B8jSK4ylHaN2LHL1lpwjZvltYYKcJCi/RK19txt3jKcc7gZCOtgXs5zr714uPXUdPLNDlzct/Mcl2WaRtUNKxcRjVqt9srDxCCgMWC8wCrlr5AuDeapElflkBU0MYiwIfkN6wnnpmNaNwMtUocZN5SG2tV8AniF7zB3esdHTjD9G5Vhfw/DJc9U2iXmBfOOiOJYRs5MRTE6dqw%2Bc6uSFUxTZ9Xid2SPmmtQiOVvBBw9tUJ9axQA1XLONcYfon8P77J5RNUJ/9g8uAQ%3D; ANON=; NAP=' --header 'Upgrade-Insecure-Requests: 1' 'https://public.bn.files.1drv.com/y4mPntow_6n32sfP1xdCn0KxNyNoX5lVrjuDizTOBkHhWyaZtaCNEEYQrX0e4PbuqmCyfrsJpMU5HX3qAV1AY1how26KCWpCCpgFjiXBCC1SDZDyUIc2yrkbZw67wWuN677PVRR5Je-DONxtpJewdvTl40Bv3rRGplENo1G4VwoE9Hv2ZpMAvLAv-OFsWsCRq7Mn7YQLXzuLpITosxZ9cX8FnhD2V46hd5jFf1tXf8LGUI?access_token=EwAAA61DBAAUzl/nWKUlBg14ZGcybuC4/OHFdfEAAeXQh5gHi4ymZUa1kxH7weXYfxZSJJ3fiJsHHpiakCtJaMbrQanlyaLNsgA7Xm3Wim%2b4%2bns9ovBzLaGxXHi0yLeL4/UcdoffK7tn6VCRNBlXHky7ccvNFmrYQJxCz/RQl20/NOPt0hPOUGlQR0qlYBmqlYdMnIVZvyRcx9iC1HaTsNcALCppKeG5G1CnFl6UEsy9GD4WCQDMQ5oafXJDNbi8BNKwv%2bmWjLtZwHBSgBaCW9Mt7P786r1jrCdnTwia84r6k15XbLSxjSAw24YeM555dd/G/0mWotDz6nDu0FwHxUbQPfzhs3p9gdK5MxbJvx8lABbl52hyf0VS3bbgPhUDZgAACIP3vSHbjttj0AGp2XUny54cJJplaZbXj3ayXHpjZqmfrJ%2bM/eOo3fv3SUgb4V1CZIAXAgBcU6tngKXw1I8cI223pb66nJKnlHmvQ05Fnp5M2LX24wMJfJDV9gxS/sUZFys0mh4gKtwIA5Pu1S62fLrqp0iWrhO80EeKJTN3zGnYTL0hMpsd6x2Au1ZjIszYH0eLHa4uk79i31NHgLEhxMXA85qNgUovIvQ5YRbAXcL42nNvi/x1mpCa41iLMmbkGgKuCKjctseaBcca9eRb%2b%2bDEyAoc%2bvc6iduFMzO2NUFDSw4U7j0BtFiLQ3bhZUMS8fqUFxWBlh0Yiz/Q2vQdkXQZ%2bW5CSaS1CIDMh%2b0HwAJfVy9i5SEQNWGxVmjI/MEH2ZEwosGwXRwSgu5/bukAyoJdp7fuIV8I3fniHuSY6LXfg%2b4DM9nFd%2bPvj699aBlObveZX7QX7ky2gycV138qVgxfiwa8yyuTddqzeHpCDpE0RsML3ktF3rIbo7XfvNWTH3gIlmyysnCU8ctUrgXdsZnuMigq1q8uriSDqRHDopZALc8jc5s2tn%2b03t22gtqSuYBW4PH%2bo9TnB3zNBEVeRadBi279od84zREZp1NT5Ee9P%2b8X2agiVmTtWAAC' --output-document 'SARG_20170328_5020.ffn.zip'

unzip SARG_20170328_5020.ffn.zip
rm -f SARG_20170328_5020.ffn.zip
${DIR}/bin/last-983/src/lastdb -Q 0 SARG_20170328_5020.ffn SARG_20170328_5020.ffn -P 10

########### ESCG database ##########################################
echo "downloading ESCG database "
wget --header 'Host: public.bn.files.1drv.com' --user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:77.0) Gecko/20100101 Firefox/77.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://onedrive.live.com/' --header 'DNT: 1' --header 'Cookie: WLSRDSecAuth=FABSARQL3KgEDBNbW84gMYrDN0fBab7xkQNmAAAEgAAACIDtreS4OWgkEAEnReiBtSdW5c1dTktkpYAMYhNOSIaCJf52VQrio3sR4F5XoLxtKMbFTLfb6xPIqMuuOE5p2/EuMQsMJ4A6FBDSKXYmVM6lf1M/ijKAUSzN3%2Beuh8dd39i5%2BM3qA%2B8jSK4ylHaN2LHL1lpwjZvltYYKcJCi/RK19txt3jKcc7gZCOtgXs5zr714uPXUdPLNDlzct/Mcl2WaRtUNKxcRjVqt9srDxCCgMWC8wCrlr5AuDeapElflkBU0MYiwIfkN6wnnpmNaNwMtUocZN5SG2tV8AniF7zB3esdHTjD9G5Vhfw/DJc9U2iXmBfOOiOJYRs5MRTE6dqw%2Bc6uSFUxTZ9Xid2SPmmtQiOVvBBw9tUJ9axQA1XLONcYfon8P77J5RNUJ/9g8uAQ%3D; ANON=; NAP=' --header 'Upgrade-Insecure-Requests: 1' 'https://public.bn.files.1drv.com/y4med-XcolaRQyVDyPOzBYyAu0VRVv64YZ_x_WZK6AKxcXx6LymtS8AeGY6y3_blfF2hF-idUHnCO-MsNC8Q8VcnE_txbYY-zr-ycBionKj3IKNgmIFVw2ip6WHlKGVHMQpqvN5nvDJx0fyWgaOuKvrfvF4PQdHt8T4g1TD4exv6B_ZpGWcYFN8yVaoTtXzk2eAvGyr28L5hInpEkolZ9B5pMmyVva9MZgOeN5mIGKcZ_Y?access_token=EwAAA61DBAAUzl/nWKUlBg14ZGcybuC4/OHFdfEAAeXQh5gHi4ymZUa1kxH7weXYfxZSJJ3fiJsHHpiakCtJaMbrQanlyaLNsgA7Xm3Wim%2b4%2bns9ovBzLaGxXHi0yLeL4/UcdoffK7tn6VCRNBlXHky7ccvNFmrYQJxCz/RQl20/NOPt0hPOUGlQR0qlYBmqlYdMnIVZvyRcx9iC1HaTsNcALCppKeG5G1CnFl6UEsy9GD4WCQDMQ5oafXJDNbi8BNKwv%2bmWjLtZwHBSgBaCW9Mt7P786r1jrCdnTwia84r6k15XbLSxjSAw24YeM555dd/G/0mWotDz6nDu0FwHxUbQPfzhs3p9gdK5MxbJvx8lABbl52hyf0VS3bbgPhUDZgAACIP3vSHbjttj0AGp2XUny54cJJplaZbXj3ayXHpjZqmfrJ%2bM/eOo3fv3SUgb4V1CZIAXAgBcU6tngKXw1I8cI223pb66nJKnlHmvQ05Fnp5M2LX24wMJfJDV9gxS/sUZFys0mh4gKtwIA5Pu1S62fLrqp0iWrhO80EeKJTN3zGnYTL0hMpsd6x2Au1ZjIszYH0eLHa4uk79i31NHgLEhxMXA85qNgUovIvQ5YRbAXcL42nNvi/x1mpCa41iLMmbkGgKuCKjctseaBcca9eRb%2b%2bDEyAoc%2bvc6iduFMzO2NUFDSw4U7j0BtFiLQ3bhZUMS8fqUFxWBlh0Yiz/Q2vQdkXQZ%2bW5CSaS1CIDMh%2b0HwAJfVy9i5SEQNWGxVmjI/MEH2ZEwosGwXRwSgu5/bukAyoJdp7fuIV8I3fniHuSY6LXfg%2b4DM9nFd%2bPvj699aBlObveZX7QX7ky2gycV138qVgxfiwa8yyuTddqzeHpCDpE0RsML3ktF3rIbo7XfvNWTH3gIlmyysnCU8ctUrgXdsZnuMigq1q8uriSDqRHDopZALc8jc5s2tn%2b03t22gtqSuYBW4PH%2bo9TnB3zNBEVeRadBi279od84zREZp1NT5Ee9P%2b8X2agiVmTtWAAC' --output-document 'ESCG.fna.zip'

unzip ESCG.fna.zip
rm -f ESCG.fna.zip
${DIR}/bin/last-983/src/lastdb -Q 0 ESCG.fna ESCG.fna -P 10


########## lineage database ###################################################################
wget https://gitlab.com/zyxue/ncbitax2lin-lineages/blob/master/lineages-2019-02-20.csv.gz
	
gunzip lineages-2019-02-20.csv.gz

######### PLSDB database #################################
wget https://ndownloader.figshare.com/files/21961095 --output-document 'PLSDB.zip'
unzip PLSDB.zip
mv data/pls/PLSDB_2020_03_04 .

$DIR/bin/ncbi-blast-2.9.0+/bin/blastdbcmd -db PLSDB_2020_03_04/plsdb.fna -dbtype nucl -entry all -outfmt "%f" -out PLSDB_2020_03_04.fna

grep ">" PLSDB_2020_03_04.fna | sed 's/>//' | sed 's/\ /\t/' > PLSDB_2020_03_04.fna.name

rm -rf data
rm PLSDB.zip


# make last database using PLSDB database
echo "Building lastdb for PLSDB, this step is quite slow, please stay patient"
${DIR}/bin/last-983/src/lastdb -Q 0 PLSDB_2020_03_04.fna.lastindex PLSDB_2020_03_04.fna -P 50


