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
wget --header 'Host: public.bn.files.1drv.com' --user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:74.0) Gecko/20100101 Firefox/74.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://onedrive.live.com/' --header 'DNT: 1' --header 'Cookie: WLSRDSecAuth=FABSARQL3KgEDBNbW84gMYrDN0fBab7xkQNmAAAEgAAACO/iQIMcitU4EAF06T5zGrTYbps3KiwkqrZjWpB5AYKIOUsoEXyT9CVefR8%2BMkPlONY6lkUVH%2B/X9JFhiwqabigrSkdX1Q9T70%2BtOhBZ4Y/J%2BFJNVlYo6dI4b98O2mjZYYk0%2BdouPnVFPYGldVXwK2MDjYASgeT4IcItJv8D9CQfY3kDnY/9gW2zEc3dtNA8/15udCfVj5LlagS8BkfkoAcTWlT%2BwSfwv9afSgLGQ/3inJjmOybI3UHFYS3/FV%2BhO0851sMXmLzOrDfd24yhCfgUzC2JQhFq19AUFwqVnRjKFrMHk5Cq9FQq6POC1MnVucvD2vw7mjuTNCagkamNHK1P61phdzV%2B6VRhkQFGbAA2AVEdyQaPidrLIxQA6DZwQwLde0WMv/p3ZuAFivLZhrs%3D; ANON=; NAP=' --header 'Upgrade-Insecure-Requests: 1' 'https://public.bn.files.1drv.com/y4mFHIf37q1KlQ0s7iW_tsqGlo4Bw8DV6GBfkBAEVmgQ_jhAIRNMeMDDmLXGeWtMkgTomfUdcSdvGXn1jt7C1-GerwDXZaEBdTP4VanKFa4GoqMlFu78YqSLgj9_Zg1Fkm28ytKd3sf-Xny-QG3UtHq2fCOUjNR0JwJbi2nMkMuO3tdgykbwWVWrabLZ7YndqRc8-WmiJDG3ZKd5uuwMB7uTRC2o0fu0TB6TTsxgQJtWYY?access_token=EwAIA61DBAAUzl/nWKUlBg14ZGcybuC4/OHFdfEAAbXcGHYVo8JfJGzevmRQmIgEZFfu4CvpMMlgDIvm8PKtOSDszPlVLf1gdtDvPZejFg5CDjzEwKY%2bg%2bKJVbm0g54ipThl%2bKdAB0aQ9DbQQPpTv1a2migU9w7rqGxesZz%2bzhi9UitrAYHXRfwZH9%2b3DM8TI7X7lCbOLIE%2bYPWm2Rg00yze2hUti%2bc2RcfjOsKzUudikpANy4suwIaqy9S7xv71DYdfIuB02lAh/a50YeZ2VxM/Gs67N5%2bZ%2bL4Re3Y2pTpuJE7bS%2bl2SxZVNs8rcN75F2FO7Re05WuNXRY7QUXnjnYj7bF5ozkyNfnoA0SdJeaYDYciV6C5nv5/3wWrD5QDZgAACIoFPknMfk/z2AGic%2bo54lqT6R6wBJ0eQpdOd6pTUNClEg8zQR0jkZNKv8Zu5H6rnd%2b7UxuT5xnEY12J6nasW70hyaktMLZboe1Z8XJrS/Cz5yrvYGGhtkSGc4%2bz6T/gDOcg9G6aOXRTpSk8hwif9NfEgLiUyCZq2Egd%2bvEpW9/y9D4A74oC4Jm1WtwH%2bwopS3VHjhwJgIGDNHE4UCz5AY22udfsbUf83Eyjf8nLNxgq%2bCczfiHa2jrqealY00c1YAVLnHfZV0K%2bomsjdUqAO2OzI/VlL3VTHp5UwzDqGs9fDiSa7PdGZvpfgRRwet2e3p8L1YVZTOhhTHFh19ZFZEariMilVpQ%2b0OsC90vOn8H3ceqFEdMVLkmrySlKkNpROi0pLxsix1/xG8PYa2JQ8PO0Qp9FC0mMgodHWeqw1GMyyyu%2bxqH3n9xWOpECjRWYCUpsDgvheI6afWE%2bTxdca7me7QaZJ9jZIJQyey2J7dNFbt2ifQ%2ba6189zBkaIK%2bCZzhC0sE4ZeEc5Bmb0mLnKpk/%2bpVKtWb/xomt4r8MYP89eaR2BcuEhr8JBmvdkK3tglfTxXnTPoGDswPrYcB17/X07hjPJhU5rya29cwIZVp43Ity55RTY38F2cd09YNObcH9AAI%3d' --output-document 'markers.fasta'

# build last database using  markergene database
${DIR}/bin/last-983/src/lastdb -Q 0 markers.lastindex markers.fasta -P 10

###### SARG-nt database ################
wget --header 'Host: public.bn.files.1drv.com' --user-agent 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:74.0) Gecko/20100101 Firefox/74.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://onedrive.live.com/' --header 'DNT: 1' --header 'Cookie: WLSRDSecAuth=FABSARQL3KgEDBNbW84gMYrDN0fBab7xkQNmAAAEgAAACIl8lLNQ3jzNEAHhRW/iihZBXqBQvjrEEsGo/g60Z%2BhpBIL6eVFC1CNfjeI9%2B4f3L6J2kS7L9QRg5%2BQ7nK8MUbG7WQ6bt7sZwFylYz0CGq1pyoraK/erkIKkYFVCsH0okb6bksZeKhJdl8JgDXIHPSY6r7SdIini33CHBJ80kMkiwXCzbAd2yUBep379szBs2J67y8nBBxmLs9esiM9o0LOq47wgd5Io4s2%2BAxkNdx8wPG2EhTI9yKJ0%2BMKeetgi%2B5iypCoy9Bn%2BKrbdRDMNleT%2B8R6TrbXLEN/iTZ1CL5jDeeQLOFcNOq/9X2LzgW7O2YG%2ByvQJTphpCgNwDkOpqYgLJsigepDOidwmubO7ga%2BtZQlUPbqtpvVxLxQA2MZjLSfsD%2B9zN08DwhqGP6xaDBE%3D; ANON=; NAP=' --header 'Upgrade-Insecure-Requests: 1' 'https://public.bn.files.1drv.com/y4mlFg3yfgbx-zYZLFpHDmAgnpagTn1MfJtjZPmcEOwMNmBp29FaeiYwsM1gBlMi5McNqNbKcaPzADjbRusRyrqv7IDjuyDoNuFJ0s8LQwuKtnmjjlcgTha1aSWBZpzLbw_wjCUZeuIu_Q8yRwh_CQec9giipP4Nrz8PJl-zOE9YSEezfKe1qafG51DfWcPy57nd0xZSbVhURH5-FCby-Nj1O51Py53KFNG_JjowIBZsT4?access_token=EwAAA61DBAAUzl/nWKUlBg14ZGcybuC4/OHFdfEAARgkX8E/yNql6mvuImZNHVIfipNHo6ctyZmmuvRQz6CS1SS6NvS1UEyBC%2bnH/mmAvR9SSBnsCLBbwhKaBseJ2LmzIN%2b5/YmAn3oeQuqcGdlI47aAsHGrF5Gl2Ervj5FvGna/Ei/3sexF/wCa8DlYMeClcyGjc2c0wHw9PsEDBjGUD7jtCn%2b0b2HyegOjB0gDgBOt8Fr4lcQEsyHe/o%2b41f/OQ%2bFq%2beXQ8Q8nmgK3usPJGNwA7CIlN7B3oJvjmZCs/IUf1dsHQKf%2bJPVdYoymVno65Xm%2bB13sHaasUNXOM5S%2b8lBta7wwzovgrUV%2bn0WgzfbJsausUrNboPqfgw1RtYYDZgAACL7p5tP0llzX0AEuus4gmCXHyidXcRRO5%2brj%2boHzRvr3sW8RujaIGmS4ehpENLbCXUX1%2b2UFg4LrmDOgqNNLIfyRNrWjddFIVVg2PbyyxsQorP7KM8EtCi0ELlYF3BeaVQe%2b8fYqhnqte2Q%2bBBUyprRSMJ5VBZOpNs5yVfAr%2b3NaObp24US8/B6ibT%2b4h5xYisbXn9vMEpXv%2bwmBtxkETFqKR1rtUlUUjvnHat4xNH9c%2b7x6ZwxrlhK/0oVSJjVuPuOqO%2b4GJ454MvXUI3yj/8F3Q/yP1JqepdDLGzmkeiVfx8v5Fj29dDigdJZXLxSmLbSQg6%2b7auKuctGcv4MdfvEy3jG%2b12R0ijPLNg6VNXsff0R0IN4MdCqTuOe3CDR5w9ZJzxCGqOA7Q1Fl98AS3J/bmja1hvyy/uO3b2sJ4fYhvADLedqp50HcXVccZRY9oNq39mu1ycOxWmv3XnsPjVW3nTksYstdgSumJmkwUeS%2bvUK9xRttVf61vXvnjz06cA5LrOYrlAM8OP5wRWL7XSeiz1ZUGbEN219EWP7boJf3U2EvKBHl4xRCJ22eCA%2bD1hpCtgZApgzsl1t%2bCTbbjL9EzWbTMWWbZ9E1AWxXbKHL/H%2bA7157O0D81AAC' --output-document 'SARG_20170328_5020.ffn'

${DIR}/bin/last-983/src/lastdb -Q 0 SARG_20170328_5020.ffn SARG_20170328_5020.ffn -P 10

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


