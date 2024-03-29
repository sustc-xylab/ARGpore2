options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(data.table)
# library(parallel)
library(foreach)
library(doParallel)

filter.plasmid<-function(df,S=0.7,nL.plasmid=0.7,nL.chimera=0.5){
   #S=similarity cutoff to filter initial last alignment
   #S=0.7
   #nL.plsmid= the contig is consider as plasmid if its hit to plasmid is longer than $nL.plasmid of its length 
   #nL.plasmid=0.7
   #nL.chimera= the contig is consider a chimera like plasmid if its hit to plasmid is shorter than $nL.chimer of its length
   #nL.chimera=0.5
   
    library(foreach)
	library(doParallel)
	library(plyr)
   colnames(df)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
   # df$query<-sapply(strsplit(df$query,"-"),"[[",1)
   
   # filtering hit based on similarity & alignment length of the Plasmid length
   lookat<-which(df$similarity>=S*100)
   df<-df[lookat,]
   
   # merge in plasmid name
   df<-merge(df,pname,by="subject")
   
   #filter overlap regions
   if(nrow(df)>0){
      #### 去除同一个region hit到多个plasmid的情况，防止对ARG的重复count ####
      tmp.lst<-split(df,df$query)
      
      # for one region exactly hit to multiple ARG, only the best hit (the one with highest bitscore also the first hit ) was kept
      tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.start),])
      tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.end),])
      # lapply(tmp.lst,nrow)
      
      # if one region hit to multiple ARG, then if the hited region is overlaped > 50% alignment length with the first hit (the hit with highest bitscore) then it will be removed, otherwise it will be kept.  ##
      
	  
	  # start_time <- Sys.time()
	    cl<-makeCluster(no_threads, outfile="")
		registerDoParallel(cl)
		pb <- txtProgressBar(min=0, max= length(tmp.lst), style = 3) # progress bar
		
		
		tmp.lst2<-list()
		tmp.lst2<-foreach(g =1:length(tmp.lst),
                    .packages = c("plyr")) %dopar% {
		# for(g in 1:length(test.lst)) {
		# tmp.lst2<-mclapply(tmp.lst,function(x){
         setTxtProgressBar(pb,g)
		 x<-tmp.lst[[g]]
         # 保证q.start<q.end 对于反向的情况就把q.start,q.end互换
         lookat<-which(apply(x,1,function(y) as.numeric(y[7])> as.numeric(y[8])))
         tmp<-x[lookat,7]
         tmp2<-x[lookat,8]
         x[lookat,7]<-tmp2
         x[lookat,8]<-tmp
         # 对于多于两行的x,做clustering,然后每个cluster filtering
         if(nrow(x)>=2){
            x<-x[order(x$bitscore,decreasing=T),] # line with highest bitscore as first line
            tmp6<-list()
            tmp5<-vector() # store the line overlaped more than 80% with the first line, these lines should be deleted
            for(i in 1:(nrow(x)-1)){
               tmp5<-vector()
               for(j in (i+1):nrow(x)){
                  tmp.start<-max(x[i,]$q.start,x[j,]$q.start)
                  tmp.end<-min(x[i,]$q.end,x[j,]$q.end)
                  overlap<-tmp.end-tmp.start
                  # if和第一条没有overlap， then overlap should be <=0, and this line should be kept for another loop
                  # overlap with first line for more than 50% of alignment length
                  if(overlap>x[j,]$align.lenth*0.5) {tmp5[j-1]<-j}
               }
               tmp6[[i]]<-tmp5
            }
            tmp6<-unique(unlist(tmp6))
            tmp6<-tmp6[!is.na(tmp6)]
            x<-x[-tmp6,]
         }
		 return(x)
      }
	  stopCluster(cl)
      # end_time <- Sys.time()
	  # end_time - start_time
      result<-rbind.fill(tmp.lst2)
      
   }  else { 
      cat("Warning!: NO Plasmid identified\n")
      result<-df
   }
   
   # since the refseq plasmid database contain two parts: 1) CDS of MGE 2) plasmid sequence, we need to further filter them with different cutoff for CDS and plasmid sequence
   # firstly seperate by CDS and plasmid sequence
   result.cds<-result[grep("complete CDS",result$plasmid.name),]
   result.plasmid<-result[grep("complete CDS",result$plasmid.name, invert = T),]
   
   
   # --- for hitted to plasmid sequence  ----
   # cases we can confirm this nanopore read is a plasmid, only when almost the whole nanopore read (df$align.length/df$q.len>0.8)is hitted to the plasmid
   lookat2<-which(result.plasmid$align.lenth/result.plasmid$q.len>nL.plasmid) 
   result.plasmid2<-result.plasmid[lookat2,]
   
   # half plasmid and half not, then probably the nanopore reads is a chimera should be removed
   lookat3<-which(result.plasmid$align.lenth/result.plasmid$q.len<nL.chimera) 
   result.chimera<-result.plasmid[lookat3,]
   
   # --- for hitted to CDS of MGE sequence  ----
   # the above filtering is enough so
   result.cds<-result.cds
   
   result.lst<-list(cds=result.cds, plasmid=result.plasmid2, chimera=result.chimera)
   
   
   return(result.lst)
}


# fread(args[1],header=F,fill=T)->plasflow
fread(args[1],header=F,fill=T)->plasmid
fread(args[2],header=F)->pname #name of plasmid
colnames(pname)<-c("subject","plasmid.name")
no_threads<-as.numeric(args[6])
plasmid.f<-filter.plasmid(plasmid,
                          S=as.numeric(args[3]),
                          nL.plasmid=0.7,
                          nL.chimera=0.5)


# plasmid.c<-intersect(plasmid.f$plasmid$query,plasflow$V1)
plasmid.c<-plasmid.f$plasmid$query
plasmid.d<-plasmid.f$plasmid[which(plasmid.f$plasmid$query %in% plasmid.c),]

if(length(plasmid.c)>0){
	
	write.table(plasmid.d[,c("query","plasmid.name")],file=args[4],sep="\t",quote=F, row.names=F)
	
	write.table(plasmid.d,file=args[5],sep="\t",quote=F, row.names=F)
} else {
	cat("	No plasmid hit against PLSDB database\n")
}



