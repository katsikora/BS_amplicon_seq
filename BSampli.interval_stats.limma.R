#run in R-3.3.1
#set working directory
wdir<-commandArgs(trailingOnly=TRUE)[1]
#system(paste0('mkdir -p ',wdir)) #for debugging
setwd(wdir)
message(sprintf("working directory is %s",getwd()))

options(stringsAsFactors=FALSE,na.rm=TRUE)

bedF<-commandArgs(trailingOnly=TRUE)[2]
message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))

bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
cnpool<-c("CHROM","START","END")#message(sprintf("processing %s",bedF))
bedshort<-gsub(".bed","",basename(bedF))
bedtab<-read.table(bedF,header=FALSE,sep="\t",as.is=TRUE,quote="")
colnames(bedtab)[1:3]<-cnpool#[1:ncol(bedtab)]
if(!unique(grepl("STRAND",colnames(bedtab)))){bedtab$STRAND<-"*"}
if(!unique(grepl("Name",colnames(bedtab)))){bedtab$Name<-paste(bedtab$CHROM,bedtab$START,bedtab$END,sep="_")}

#load single CpG data, count NAs per row/CpG
sCpGF<-commandArgs(trailingOnly=TRUE)[3]
message(sprintf("loading %s",sCpGF))

require(data.table)
load(sCpGF)

noNA<-apply(limdat.LG,1,function(X)sum(is.na(X)))
NAf<-ifelse(noNA>1,1,0)
limdat.LG$NAf<-NAf

require(GenomicRanges)

bedGR<-GRanges(seqnames=bedtab$CHROM,ranges=IRanges(start=bedtab$START,end=bedtab$END,names=bedtab$Name),strand=bedtab$STRAND)
limdatGR<-GRanges(seqnames=gsub("_.+","",limdat.LG$ms),ranges=IRanges(start=as.numeric(gsub(".+_","",limdat.LG$ms)),width=2,names=limdat.LG$ms),strand="*")
mtch <- as.data.frame(findOverlaps(query=limdatGR, subject=bedGR))

limdat.LG.inCGI<-limdat.LG[mtch$queryHits,]
limdat.LG.inCGI$IntID<-names(ranges(bedGR))[mtch$subjectHits]

NA.inCGI<-with(limdat.LG.inCGI,ave(NAf,factor(IntID),FUN=sum))
limdat.LG.inCGI$NA.inCGI<-NA.inCGI

to.m<-limdat.LG.inCGI[,c("IntID", "NA.inCGI"),with=FALSE]

CGI.map<-unique(to.m)
bedtab$N.CG.NA<-CGI.map$NA.inCGI[match(bedtab$Name,CGI.map$IntID)]

N.CG.tot<-with(limdat.LG.inCGI,ave(NAf,IntID,FUN=length))
bedtab$N.CG.tot<-N.CG.tot[match(bedtab$Name,limdat.LG.inCGI$IntID)]

bedtab$CGI.NAf<-ifelse(bedtab$N.CG.NA>(0.8*bedtab$N.CG.tot),NA,1)
bedtab.CC<-bedtab[complete.cases(bedtab),]

CGI.limdat<-as.data.frame(apply(limdat.LG.inCGI[,2:(ncol(limdat.LG.inCGI)-3)],2,function(X){ave(X,factor(limdat.LG.inCGI$IntID),FUN=function(X)mean(X,na.rm=TRUE))}),stringsAsFactors=FALSE)

CGI.limdat$IntID<-limdat.LG.inCGI$IntID
CGI.limdat<-unique(CGI.limdat)
rownames(CGI.limdat)<-CGI.limdat$IntID
CGI.limdat<-CGI.limdat[,-grep("IntID",colnames(CGI.limdat))]

CGI.limdat.CC<-CGI.limdat[bedtab.CC$Name,]



####for differential interval methylation
### limma + ebayes + BH p value adjustment

require("limma")
require("car")
require("FactoMineR")
require("reshape2")
require("ggplot2")
require("dplyr")

CGI.limdat.CC.logit<-logit(CGI.limdat.CC,percents=FALSE,adjust=0.025)
x1<-PCA(CGI.limdat.CC,graph=FALSE)

pdf(paste0(bedshort,".CGI.limdat.CC.PCA.pdf"),paper="a4",bg="white") 
plot.PCA(x1,choix="var")
dev.off()

#calculate row means
spath<-commandArgs(trailingOnly=TRUE)[4]
sampleInfo<-read.table(spath,header=TRUE,sep="\t",as.is=TRUE)
#calculate and save row means
CGI.limdat.CC$IntID<-rownames(CGI.limdat.CC)
CGI.limdat.CC.L<-melt(CGI.limdat.CC,id.vars="IntID",value.name="Beta",variable.name="SampleID")
CGI.limdat.CC.L$Group<-sampleInfo$Group[match(CGI.limdat.CC.L$SampleID,sampleInfo$PlottingID)]
CGI.limdat.CC.Means<-data.table(summarize(group_by(CGI.limdat.CC.L,IntID,Group),Beta.Mean=mean(Beta)))

##violin plots
ggplot(data=CGI.limdat.CC.Means)+geom_violin(aes(x=Group,y=Beta.Mean,fill=Group))+geom_boxplot(aes(x=Group,y=Beta.Mean),width=0.1)+ggtitle("All CG sites")+
    theme(text = element_text(size=16),axis.text = element_text(size=12),axis.title = element_text(size=14))+xlab("Mean methylation ratio")
ggsave(paste0(bedshort,".Beta.MeanXgroup.all.violin.png"))

save(bedtab,limdat.LG.inCGI,CGI.limdat.CC,CGI.limdat.CC.Means,file=paste0(bedshort,".aggCpG.RData"))

#############################################

if(!is.null(sampleInfo$Timepoint)){
    tu<-unique(sampleInfo$Timepoint)
    CGI.limdat.CC.Means$Timepoint<-sampleInfo$Timepoint[match(CGI.limdat.CC.Means$SampleID,sampleInfo$SampleID)]}else{
    sampleInfo$Timepoint<-""
    CGI.limdat.CC.Means$Timepoint<-sampleInfo$Timepoint[match(CGI.limdat.CC.Means$SampleID,sampleInfo$SampleID)]}

for(i in seq_along(tu)){
    sampleInfo.tu<-sampleInfo[sampleInfo$Timepoint %in% tu[i],]
    CGI.limdat.CC.tu<-CGI.limdat.CC[,colnames(CGI.limdat.CC) %in% sampleInfo.tu$SampleID]
    x1<-PCA(CGI.limdat.CC.tu,graph=FALSE)

    pdf(paste0(tu[i],".CGI.limdat.CC.PCA.pdf"),paper="a4",bg="white") 
    plot.PCA(x1,choix="var")
    dev.off()

  
    CGI.limdat.CC.Means.tu<-CGI.limdat.CC.Means[CGI.limdat.CC.Means$Timepoint %in% tu[i],]
    CGI.limdat.CC.tu$IntID<-rownames(CGI.limdat.CC.tu)
    CGI.limdat.CC.tu.L<-melt(CGI.limdat.CC.tu,id.vars="IntID",value.name="Beta.Mean",variable.name="SampleID")
    CGI.limdat.CC.tu.L$Group<-sampleInfo.tu$Group[match(CGI.limdat.CC.tu.L$SampleID,sampleInfo.tu$SampleID)]
    CGI.limdat.CC.tu.L$Group<-factor(CGI.limdat.CC.tu.L$Group,levels=c("WT","MUT"))
    CGI.limdat.CC.tu.L$IntID<-factor(CGI.limdat.CC.tu.L$IntID,levels=c("krt8","rarga","Vasa","Ntla"))

###dot per interval per sample
    ggplot(data=CGI.limdat.CC.tu.L,aes(x=Group,y=Beta.Mean))+geom_jitter(size=4,width=0.1,height=0.00005,alpha=0.6,aes(colour=Group))+facet_grid(.~IntID)+ggtitle(paste0("Mean methylation rate per region in ",tu[i]))+theme(text = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=14),axis.text.x=element_text(size=12,angle=90))+xlab("")+ylab("Mean methylation rate")+ylim(0,1)+scale_colour_manual(values=c("grey28","red"))
    ggsave(paste0(tu[i],".BetaMeanXInt.perGroup.dotplot.png"),height=10,width=16)

    write.table(CGI.limdat.CC.tu.L,file=paste0("CGI.limdat.CC.",tu[i],".L.xls"),sep="\t",quote=FALSE,row.names=FALSE)


####barplot per interval per group
    plotdat2<-summarySE(CGI.limdat.CC.tu.L, measurevar="Beta.Mean", groupvars=c("Group","IntID"))
    ggplot(plotdat2, aes(x=Group, y=Beta.Mean)) + 
    geom_bar(stat="identity",position="identity",aes(fill=Group))+
    geom_errorbar(aes(ymin=Beta.Mean-sd, ymax=Beta.Mean+sd),width=.2)+facet_grid(.~IntID)+
   theme(text = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=14),axis.text.x=element_text(size=12,angle=90))+xlab("")+ylab("Mean methylation rate")+ylim(0,1)+scale_fill_manual(values=c("grey28","red"))+ggtitle(paste0("Mean methylation rate per region in ",tu[i]))
    ggsave(paste0(tu[i],".BetaMeanXInt.perGroup.barplot.png"),height=8,width=12)

    write.table(plotdat2,file=paste0("plotdat2.",tu[i],".L.xls"),sep="\t",quote=FALSE,row.names=FALSE)

    print(paste0(i,"_processed"))

###stats WT vs MUT
    plotdat<-CGI.limdat.CC.tu.L
    plotdat$logit25<-logit(p=plotdat$Beta.Mean,percents=FALSE,adjust=0.025)
    GSv<-unique(as.character(plotdat$IntID))
    res<-data.frame(GSv,stringsAsFactors=FALSE)
    res$MUT<-NA
    res$WT<-NA
    res$logit.MUT<-NA
    res$logit.WT<-NA
    res$pval<-NA
    plotdat<-plotdat
    for(k in seq_along(GSv)){
        td<-plotdat[plotdat$IntID %in% GSv[k],]

        res$MUT[k]<-mean(td$Beta.Mean[td$Group %in% "MUT"])
        res$WT[k]<-mean(td$Beta.Mean[td$Group %in% "WT"])
        res$logit.MUT[k]<-mean(td$logit25[td$Group %in% "MUT"])
        res$logit.WT[k]<-mean(td$logit25[td$Group %in% "WT"])

        z<-t.test(x=td$logit25[td$Group %in% "MUT"],y=td$logit25[td$Group %in% "WT"],var.equal=TRUE)
        res$pval[k]<-z$p.value

    }

    res<-res[order(res$pval),]

    write.table(res,file=paste0("interval.mean.",tu[i],".ttest.xls"),row.names=FALSE,quote=FALSE,sep="\t")


}

