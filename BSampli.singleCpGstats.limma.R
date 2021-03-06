#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))


options(stringsAsFactors=FALSE,na.rm=TRUE)
require("limma")
require("car")

###read in sample sheet

spath<-commandArgs(trailingOnly=TRUE)[2]
sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)


############read in methylation tracks; use (reference)-sorted bed files #############################################
#sorted bed file was used at the beginning to extract CGs from intervals of interest + 0-padding was used -> all files should have same row ordering (although possibly different from the reference)
###deal with non-0-padded tables from methylDackel and BisSNP

mpath<-commandArgs(trailingOnly=TRUE)[3]
mdir<-dir(mpath,pattern="*CpG.filt2.bed",full.names=TRUE)
mshort<-gsub("_prin.CpG.filt2.bed","",basename(mdir))

cC<-c(rep("NULL",3),"numeric",rep("NULL",3),"character")

require(data.table)

mlist<-vector("list",length(mdir))
for(i in seq_along(mdir)){
    tabi<-fread(mdir[i],colClasses=cC,sep="\t",header=TRUE)
    colnames(tabi)[colnames(tabi) %in% "Beta"]<-mshort[i]
    mlist[[i]]<-tabi
}

limdat<-Reduce(function(...) merge(..., all=T,by="ms",sort=FALSE), mlist)
limdat<-limdat[,c(1,match(sampleInfo$SampleID,colnames(limdat))),with=FALSE]


limdat.LG<-limdat 
if(unique(grepl("MethylDackel",mdir))){
    limdat.LG[,2:ncol(limdat.LG)]<-limdat[,2:ncol(limdat),with=FALSE]/100
}
if("Merge" %in% colnames(sampleInfo)){
    if(unique(grepl("methylCtools",mdir))){ colnames(limdat.LG)[2:ncol(limdat.LG)] <- gsub(".CG.call","",colnames(limdat.LG)[2:ncol(limdat.LG)])}
    colnames(limdat.LG)[2:ncol(limdat.LG)]<-sampleInfo$PlottingID[match(colnames(limdat.LG[,2:ncol(limdat.LG)]),sampleInfo$Merge)]}else{colnames(limdat.LG)[2:ncol(limdat.LG)]<-sampleInfo$PlottingID[match(colnames(limdat.LG[,2:ncol(limdat.LG)]),sampleInfo$SampleID)]}

limdat.LG.CC<-limdat.LG[complete.cases(limdat.LG),] 


limdat.LG.CC.logit<-logit(limdat.LG.CC[,2:ncol(limdat.LG.CC),with=FALSE],percents=FALSE,adjust=0.025)
rownames(limdat.LG.CC.logit)<-limdat.LG.CC$ms


require("FactoMineR")
x1<-PCA(limdat.LG.CC[,-1,with=FALSE],graph=FALSE)

pdf("limdat.LG.CC.PCA.pdf",paper="a4",bg="white")
plot.PCA(x1,choix="var")
dev.off()

########################## prepare density plots per group ##################################################

require("ggplot2")
require("reshape2")
require(dplyr)

#calculate and save row means
limdat.LG.CC.L<-melt(limdat.LG.CC,id.vars="ms",value.name="Beta",variable.name="SampleID")
limdat.LG.CC.L$SampleID<-as.character(limdat.LG.CC.L$SampleID)
limdat.LG.CC.L$Group<-sampleInfo$Group[match(limdat.LG.CC.L$SampleID,sampleInfo$PlottingID)]
limdat.LG.CC.Means<-data.table(summarize(group_by(limdat.LG.CC.L,ms,Group),Beta.Mean=mean(Beta)))

head(limdat.LG.CC.Means)

ggplot(data=limdat.LG.CC.Means,aes(x=Beta.Mean))+geom_density(aes(group=Group,colour=Group,fill=Group),alpha=0.3)+ggtitle("All CG sites")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
ggsave(filename="Beta.MeanXgroup.all.dens.png")

save(limdat.LG,file="limdat.LG.RData")
save(limdat.LG.CC.Means,limdat.LG,file="singleCpG.RData")



