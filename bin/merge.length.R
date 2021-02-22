options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(data.table)
library(dplyr)


blast_tab<-fread(args[1],header=F,sep="\t")
query.len<-fread(args[2],header=F,sep="\t")
db.len<-fread(args[3],header=F,sep="\t")

colnames(blast_tab)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore")
colnames(query.len)<-c("query","q.len")
colnames(db.len)<-c("subject","s.len")

df<-full_join(blast_tab,query.len,by="query")
df2<-full_join(df,db.len,by="subject")
df2<-arrange(df2,query,desc(bitscore))
df2<-df2[,c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")]

write.table(df2,file=args[4],row.names=F,col.names=F,quote=F,sep="\t")
