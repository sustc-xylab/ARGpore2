options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(data.table)
library(foreach)
library(doParallel)


lineage<-fread(args[1])
colnames(lineage)<-c("tax_id","kingdom","phylum","class","order","family","genus","species")

cen<-fread(args[2])
#cen<-cen[which(cen$hitLength>150),] # only kept hitlenght > 150
colnames(cen)[3]<-c("tax_id")
cen<-merge(cen,lineage, by="tax_id")
cen<-as.data.frame(cen)
for(i in 9:ncol(cen)){
	cen[,i]<-sub("__","",cen[,i])
}

n_threads<-as.numeric(args[3])

fasta.name<-fread(args[4])
colnames(fasta.name)<-c("contig","length")
fasta.name<-data.frame(fasta.name)
# fasta.name<-fasta.name[,1]
fasta.name$kingdom<-""
fasta.name$phylum<-""
fasta.name$class<-""
fasta.name$order<-""
fasta.name$family<-""
fasta.name$genus<-""
fasta.name$species<-""
fasta.name<-fasta.name[,c("contig","kingdom","phylum","class","order","family","genus","species")]

# --- taxa voting for centrifuge taxa : cen2---------
cat("taxa voting for centrifuge result\n")

cen<-arrange(cen, readID, score)

cen2.nv<-cen[cen$numMatches<=1,]
cen2.v<-cen[cen$numMatches>1,]

if(nrow(cen2.v)>=1){
tmp.readID<-unique(cen2.v$readID)
rank.lst<-list()
cl<-makeCluster(n_threads, outfile="")
registerDoParallel(cl)
pb <- txtProgressBar(min=0, max= length(tmp.readID), style = 3) # progress bar
  
rank.lst<-foreach(i =1:length(tmp.readID),
                    .packages = c("plyr")) %dopar% {
                      
					  setTxtProgressBar(pb,i)
                      # for(i in 1:length(tmp.readID)){
                      x<-cen2.v[cen2.v$readID==tmp.readID[i],]
                      tmprank.lst<-list()
                      
                      # voting cutoff of 50% for superkingdom, phylum, class 
                      tmprank<-c("kingdom","phylum","class")
					  cutoff.spc<-0.5
                      for(j in 1:length(tmprank)){
                        tmpd<-data.frame(table(x[,tmprank[j]]))
                        tmpd$Freq.per<-tmpd$Freq/sum(tmpd$Freq)
                        tmpd<-arrange(tmpd, desc(Freq))
                        tmpd<-tmpd[tmpd$Freq.per>cutoff.spc,]
                        if(nrow(tmpd)==0){tmprank.lst[tmprank[j]]<-NA}
                        else{tmprank.lst[tmprank[j]]<-as.character(tmpd$Var1[1])}
                        
                      }
                      
                      # voting cutoff of 30% for order,family,genus,species
                      tmprank<-c("order","family","genus","species")
					  cutoff.ofgss<-0.3
                      for(j in 1:length(tmprank)){
                        tmpd<-data.frame(table(x[,tmprank[j]]))
                        tmpd$Freq.per<-tmpd$Freq/sum(tmpd$Freq)
                        tmpd<-arrange(tmpd, desc(Freq))
                        tmpd<-tmpd[tmpd$Freq.per>cutoff.ofgss,] 
                        if(nrow(tmpd)==0){tmprank.lst[tmprank[j]]<-NA}
                        else{tmprank.lst[tmprank[j]]<-as.character(tmpd$Var1[1])}
                      }
                      
                      # if the top level is NA then the following level should all be NA
                      result<-unlist(tmprank.lst)
                      for(j in 2:length(result)){
                        if(is.na(result[j-1])){result[j:length(result)]<-NA
                        break}
                      }
                      
                      return(result)
                      # rank.lst[[i]]<-result
                      #cat(paste("finish ",i,"\n", sep=" "))
                    }
stopCluster(cl)
stopImplicitCluster()
cat("\nfinish voting\n")
  
# merge voting result with those don't need to do voting to get cen2
tmp.lst2<-lapply(rank.lst, function(x) data.frame(t(matrix(x))))
tmp<-rbind.fill(tmp.lst2)
rm(tmp.lst2 ) # rm this big list
tmp2<-data.frame(readID=unique(cen2.v$readID))
cen2.v<-cbind(tmp2,tmp)
rank<-c("kingdom","phylum","class","order","family","genus","species")
colnames(cen2.v)<-c("readID",rank)
cen2.nv<-cen2.nv[,c("readID",rank)]
cen2<-rbind(cen2.nv,cen2.v)
} else {
rank<-c("kingdom","phylum","class","order","family","genus","species")
cen2.nv<-cen2.nv[,c("readID",rank)]
cen2<-cen2.nv}


colnames(cen2)<-c("contig",rank)


# read in the 2D.fa last marker gene result #####
taxa<-fread(args[5],header=F)
colnames(taxa)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
load(args[6])

# filtering hit based on similarity & alignment length of the marker gene length
lookat<-which(taxa$similarity> args[7] & taxa$align.lenth/taxa$s.len>args[8])
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
minl.c<-ncol(cen2)-1
minl.t<-ncol(taxa)-1
if(minl.c<7){for(j in (minl.c+1):7){cen2[,rank[j]]<-""}} # subset rank in case none of the ont read is assigned to species
if(minl.t<7){for(j in (minl.t+1):7){taxa[,rank[j]]<-""}} 

merge.taxa<-function(k,t){
	# k as the basis for taxa of t to merge into, taxa from t is merged in only when t's above level taxa is equal to k

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
			lookat2<-which(m[,x2]==m[,y2] & all(m[,x2]=="",m[,y2]==""))
			lookat3<-intersect(lookat,lookat2)
			m[,rank[i]]<-m[,x]
			if(length(lookat3)>0) {
			
			m[,rank[i]][lookat3]<-m[,y][lookat3]}
		}
	}
	# m<-m[,1:8]
	m<-m[,c("contig",rank)]
	colnames(m)<-c("contig",rank)
	return(m)
}

# use which centrifuge as basis to merge in markergene
# merge in taxa from marker only when marker's above level taxa is equal to centrifuge
m<-merge.taxa(data.frame(cen2),data.frame(taxa))

m<-rbind(m,fasta.name[which(!fasta.name$contig %in% m$contig),])


write.table(m,file=args[9],row.name=F,col.name=T, quote=F, sep="\t")

# #### print out the unclassified ratio########
# m2<-m
# c<-vector()
# for(i in 1:length(rank)){
	# n<-length(which(m2[,rank[i]]==""))
	# c[i]<-paste(round(n/nrow(m2)*100,2),"%",sep="")
	# }
# t<-matrix(c,ncol=length(rank))
# colnames(t)<-rank
# rownames(t)<-c("unclassified ratio")
# classified<-paste(args[9],"_unclassified.ratio.tab",sep="")
# write.table(t,file=classified,row.name=T,col.name=T, quote=F, sep="\t")

