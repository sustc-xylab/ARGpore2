#!/usr/bin/env Rscript
system("echo \n")
#### read in system arguements ######
options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
# args<-c(
# "test.fa.uniq_sarg/test.fa.uniq_sarg.last",
# "test.fa.uniq_sarg/test.fa.uniq_escg.last",
# "20",
# "../database/structure.RData",
# "test.fa.uniq_taxa.tab",
# "75",
# "0.9",
# "test.fa.uniq_plasmid.like.tab",
# "test.fa.uniq_arg.w.taxa.tab",
# "test.fa.uniq_arg.tab",
# "test.fa.uniq_arg.summary.tab")

library(plyr)
library(data.table)
library(foreach)
library(doParallel)

# read in new SNP-free SARG ----
arg<-fread(args[1],header=F)
colnames(arg)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")

escg<-fread(args[2],header = F)
colnames(escg)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")

no_threads<-as.numeric(args[3])

load(args[4])
overlap.db<-l

taxa<-fread(args[5])
colnames(taxa)[1]<-c("query")

simcutoff=as.numeric(args[6])
lencutoff=as.numeric(args[7])

plasmid<-fread(args[8])
plasmid$source<-c("plasmid.like")


########################
# SARG filtering
########################
arg.filter<-function(arg.coliform,simcutoff=70,lencutoff=0.9){
	colnames(arg.coliform)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
	# arg.coliform$query<-sapply(strsplit(arg.coliform$query,"-"),"[[",1)

	# filtering hit based on similarity & alignment length of the ARG length
	lookat<-which(arg.coliform$similarity>simcutoff & arg.coliform$align.lenth/arg.coliform$s.len>lencutoff)
	arg.coliform2<-arg.coliform[lookat,]

	if(nrow(arg.coliform2)>0){
	  #### remove the case where the same region on nanopore read hit to multiple ARGs ####
	  tmp.lst<-split(arg.coliform2,arg.coliform2$query)
	  
	  # for one region exactly hit to multiple ARG, only the best hit (the one with highest bitscore also the first hit ) was kept
	  tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.start),])
	  tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.end),])
	  # lapply(tmp.lst,nrow)
	  
	  # if one region hit to multiple ARG, then if the hited region is overlaped > 50% alignment length with the first hit (the hit with highest bitscore) then it will be removed, otherwise it will be kept.  ##
	  
	  for(g in 1:length(tmp.lst)) {
		x<-tmp.lst[[g]]
		# makesure q.start<q.end, flip q.end and q.start if not satisfy this standard
		lookat<-which(apply(x,1,function(y) as.numeric(y[7])> as.numeric(y[8])))
		tmp<-x[lookat,7]
		tmp2<-x[lookat,8]
		x[lookat,7]<-tmp2
		x[lookat,8]<-tmp
		# if nrow(x) > 2 then need to do clustering and then filter each cluster 
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
			  # if no overlap with the first line， then overlap should be <=0, and this line should be kept for another loop
			  # overlap with first line for more than 50% of alignment length
			  if(overlap>x[j,]$align.lenth*0.5) {tmp5[j-1]<-j}
			}
			tmp6[[i]]<-tmp5
		  }
		  tmp6<-unique(unlist(tmp6))
		  tmp6<-tmp6[!is.na(tmp6)]
		  if(length(tmp6>0)){tmp.lst[[g]]<-x[-tmp6,]}
		  else {tmp.lst[[g]]<-x}
		}
	  }

	  
	  arg.colifom4<-rbind.fill(tmp.lst)
	  arg.colifom4$acc<-sapply(strsplit(as.character(arg.colifom4$subject),":"),"[[",1)

	  
	  # merge ARG annotation with ARDB type #
	  tmp<-merge(arg.colifom4,overlap.db,by="acc")
	  arg.colifom4<-tmp
	  
	}  else { 
	  cat("Warning!: NO ARG identified\nbelow items will be empty；\narg.tab\narg.w.taxa.tab\n")
	  arg.colifom4<-arg.coliform2
	}
	
	return(arg.colifom4)


}


arg.f<-arg.filter(arg,
				simcutoff=simcutoff,
				lencutoff=lencutoff
)
arg.f<-arg.f[,c("query","subtype","type","q.start","q.end")]

##############################
# ESCG filtering
##############################
escg.filter<-function(df,S=70,L=0.9,no_threads=1){
  library(foreach)
  library(doParallel)
  library(plyr)
  
  colnames(df)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")

  # filtering hit based on similarity & alignment length of the ESCG length
  lookat<-which(df$similarity>S & df$align.lenth/df$s.len>L)
  df<-df[lookat,]
  
  # df$query<-sapply(strsplit(as.character(df$query),"-"),"[[",1)
  df$ko<-sapply(strsplit(as.character(df$subject),"::"),"[[",2)
  
  #filter overlap
  cat("\nfiltering overlap escg hits\n")
  if(nrow(df)==0) { 
    cat("Warning!: NO ESCG identified")
    df.filtered<-df
  } else   {
    
    tmpdf<-data.frame(table(df$query))
    tmp.query1<-tmpdf$Var1[tmpdf$Freq>1] # queries need filtering
    tmp.query2<-tmpdf$Var1[tmpdf$Freq==1] # queries don't need filtering
    
    # if one region hit to multiple ESCG, then if the hited region is overlaped
    # > 50% alignment length with the first hit (the hit with highest bitscore)
    # then it will be removed, otherwise it will be kept.  ##
    cl<-makeCluster(no_threads, outfile="")
    registerDoParallel(cl)
    pb <- txtProgressBar(min=0, max= length(tmp.query1), style = 3) # progress bar
    
    
    result<-list()
    result<-foreach(g =1:length(tmp.query1),
                    .packages = c("plyr","doParallel")) %dopar% {
                      
                      # for(g in 1:length(tmp.query1)) {
                      
                      x<-df[which(df$query==tmp.query1[g]),]
                      setTxtProgressBar(pb,g)
                      
                      # for one region exactly hit to multiple ESCG, only the best hit (the one
                      # with highest bitscore also the first hit ) was kept
                      x$s_e<-paste(x$q.start,x$q.end,sep="_")
                      x<-x[!duplicated(x$s_e),] 
                      x<-x[,1:(ncol(x)-1)] # remove this s_e column 
                      
                      
                      # make sure the nubmer of q.start is smaller than q.end 
                      lookat<-which(apply(x,1,function(y) as.numeric(y[7])> as.numeric(y[8])))
                      tmp<-x[lookat,7]
                      tmp2<-x[lookat,8]
                      x[lookat,7]<-tmp2
                      x[lookat,8]<-tmp
                      
                      if(nrow(x)>=2){
                        x<-x[order(x$bitscore,decreasing=T),] # line with highest bitscore as first line
                        tmp6<-list() # store the line overlaped more than 50% with the first line, these lines should be deleted
                        tmp5<-vector() 
                        for(i in 1:(nrow(x)-1)){
                          tmp5<-vector()
                          for(j in (i+1):nrow(x)){
                            tmp.start<-max(x[i,]$q.start,x[j,]$q.start)
                            tmp.end<-min(x[i,]$q.end,x[j,]$q.end)
                            overlap<-tmp.end-tmp.start
                            # ifåç¬¬ä¸æ¡æ²¡æoverlapï¼ then overlap should be <=0, and this line should be kept for another loop
                            # overlap with first line for more than 50% of alignment length
                            if(overlap>x[j,]$align.lenth*0.5) {tmp5[j-1]<-j}
                          }
                          tmp6[[i]]<-tmp5
                        }
                        tmp6<-unique(unlist(tmp6))
                        tmp6<-tmp6[!is.na(tmp6)]
                        if(length(tmp6)>0){x2<-x[-tmp6,]} else {x2<-x}
                        
                        # one reads shall not kit to multiple ESCG, therefore only the best hit one is kept
                        x2<-x2[!duplicated(x2$ko),] 
                      } else {x2<-x}
                      
                      return(x2)
                    }
    
    stopCluster(cl)
    tmp1<-rbind.fill(result) # after filtering 
    tmp2<-df[which(df$query %in% tmp.query2),] # those without filtering
    df.filtered<-rbind(tmp1,tmp2)
    
  } 
  # calclate no_c based on df.filtered
  # cat("\nCalculating cell number based on KO annotation \n")
  tmp.ko<-table(df.filtered$ko)
  # no_c<-mean(tmp.ko) 
  result2<-list(df.filtered=df.filtered, no_c=tmp.ko)
  return(result2)
  
}


escg.f<-escg.filter(escg,
				S=simcutoff,
				L=lencutoff,
				no_threads=no_threads
)

# cell number estimation qual to ESCG showing the largest copy 
NO.c<-max(escg.f$no_c)

#cat("done calculating cell number\n")
###############
# combine ARG profile with taxa profile
# arg.summary: the final results
################

# get ARG profile of nanopore query with taxa classification
# results in arg.summary
tmp<-arg.f[,c("query","subtype","type")]
write.table(tmp,file=args[10],quote=F,row.names = F,sep="\t")

if(nrow(arg.f)>0){
	arg.w.taxa<-merge(arg.f,taxa,by="query",all.x=T)
	arg.w.taxa<-merge(arg.w.taxa,plasmid,by="query",all.x=T)
	arg.w.taxa[is.na(arg.w.taxa)]<-""
	cat("done here\n")
	
	
	arg.c<-aggregate( query~subtype+type,arg.w.taxa,length)
	colnames(arg.c)[3]<-c("No.reads")
	
	# --  write out ----
	write.table(arg.w.taxa,file=args[9], quote=F,row.names = F,sep="\t")
	write.table(arg.c,file=args[11], quote=F,row.names = F,sep="\t")
} else { 
  cat("Warning: NO ARG identified\nOnly taxa annotations were generated \n")
  
}
