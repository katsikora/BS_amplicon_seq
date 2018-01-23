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

###dot per interval per sample
ggplot(data=CGI.limdat.CC.Means,aes(x=Group,y=Beta.Mean))+geom_jitter(size=4,width=0.1,height=0.00005,alpha=0.6,aes(colour=Group))+facet_grid(IntID~.)+ggtitle("Mean methylation rate per region")+
   theme(text = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=14),axis.text.x=element_text(size=12,angle=90))+xlab("")+ylab("Mean methylation rate")+ylim(0,1)+scale_colour_manual(values=c("blue","darkgreen","grey28","red","turquoise4","orange4"))
ggsave(paste0(bedshort,".BetaMeanXInt.perGroup.dotplot.png"),height=10,width=16)

##barplots mean +/- stdev
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    datac$ci.sd <- datac$sd * ciMult

    return(datac)
}


plotdat2<-summarySE(CGI.limdat.CC.Means, measurevar="Beta.Mean", groupvars=c("Group","IntID"))

ggplot(plotdat2, aes(x=Group, y=Beta.Mean)) + 
    geom_errorbar(aes(ymin=Beta.Mean-sd, ymax=Beta.Mean+sd),width=.2) +
    geom_bar(stat="identity",position="identity",aes(fill=Group))+facet_grid(.~IntID)+
   theme(text = element_text(size=16),axis.text.y = element_text(size=12),axis.title = element_text(size=14),axis.text.x=element_text(size=12,angle=90))+xlab("")+ylab("Mean methylation rate")+ylim(0,1)
ggsave(paste0(bedshort,".BetaMeanXInt.perGroup.barplot.png"),height=8,width=12)

