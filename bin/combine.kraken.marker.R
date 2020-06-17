options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(data.table)


lineage<-fread(args[1])
colnames(lineage)<-c("tax_id","kingdom","phylum","class","order","family","genus","species")

krak<-fread(args[2],header=F)
krak<-krak[,2:3]
colnames(krak)<-c("contig","tax_id")

kraken<-merge(krak,lineage, by="tax_id")
kraken<-kraken[,-1]
kraken<-as.data.frame(kraken)
for(i in 2:ncol(kraken)){
	kraken[,i]<-sub("__","",kraken[,i])
}


# read in the 2D.fa last marker gene result
taxa<-fread(args[3],header=F)
colnames(taxa)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
load(args[4])

# filtering hit based on similarity & alignment length of the marker gene length
lookat<-which(taxa$similarity> args[5] & taxa$align.lenth/taxa$s.len>args[6])
taxa<-taxa[lookat,]
taxa<-merge(taxa,taxa.info,by="subject")
taxa<-arrange(taxa,query,desc(bitscore))

# for each nanopore read only the taxa with the highest bitscore is kept。
taxa<-taxa[!duplicated(taxa$query),] # default deduplicated 就是取第一条
colnames(taxa)[2]<-c("contig")
taxa$kindom<-sub("k__","",taxa$kindom)
colnames(taxa)[15]<-c("kingdom")
taxa$phylum<-sub("p__","",taxa$phylum)
taxa$class<-sub("c__","",taxa$class)
taxa$order<-sub("o__","",taxa$order)
taxa$family<-sub("f__","",taxa$family)
taxa$genus<-sub("g__","",taxa$genus)
taxa$species<-sub("s__","",taxa$species)
taxa<-taxa[,c("contig","kingdom","phylum","class","order","family","genus","species")]
taxa[is.na(taxa)]<-""


rank<-c("kingdom","phylum","class","order","family","genus","species")
minl.k<-ncol(kraken)-1
minl.t<-ncol(taxa)-1
if(minl.k<7){for(j in (minl.k+1):7){kraken[,rank[j]]<-""}} # subset rank in case none of the ont read is assigned to species
if(minl.t<7){for(j in (minl.t+1):7){taxa[,rank[j]]<-""}} 

merge.taxa<-function(k,t){
	# k as the basis to merge in t, taxa from t is merged in only when t's above level taxa is equal to k

	m<-merge(k,t, by="contig", all=T)
	m[is.na(m)]<-""
	m<-as.data.frame(m)
	for(i in 1:length(rank)){
		#target level
		x<-paste(rank[i],"x",sep=".") # corespond to k 
		y<-paste(rank[i],"y",sep=".") # corespond to t
		lookat<-which(m[,x]=="")
		if(i==1) {
		m[,rank[i]]<-m[,x]
		m[,rank[i]][lookat]<-m[,y][lookat]
		} 
		else {
			# previous level
			x2<-paste(rank[i-1],"x",sep=".") # corespond to k
			y2<-paste(rank[i-1],"y",sep=".") # corespond to t
			lookat2<-which(m[,x2]==m[,y2])
			lookat3<-intersect(lookat,lookat2)
			if(length(lookat3)>0) {
			m[,rank[i]]<-m[,x]
			m[,rank[i]][lookat3]<-m[,y][lookat3]}
		}
	}
	# m<-m[,1:8]
	m<-m[,c("contig",rank)]
	colnames(m)<-c("contig",rank)
	return(m)
}

# use which kraken as basis to merge in markergene
# merge in taxa from marker only when marker's above level taxa is equal to kraken
m<-merge.taxa(data.frame(kraken),data.frame(taxa))

write.table(m,file=args[7],row.name=F,col.name=T, quote=F, sep="\t")

#### print out the unclassified ratio########
m2<-m
c<-vector()
for(i in 1:length(rank)){
	n<-length(which(m2[,rank[i]]==""))
	c[i]<-paste(round(n/nrow(m2)*100,2),"%",sep="")
	}
t<-matrix(c,ncol=length(rank))
colnames(t)<-rank
rownames(t)<-c("unclassified ratio")
classified<-paste(args[7],"_unclassified.ratio.tab",sep="")
write.table(t,file=classified,row.name=T,col.name=T, quote=F, sep="\t")

