import os
import sys
import time
import argparse

#to-do: reduce the imports, take too long
import zipfile
import re
import subprocess
import pandas
import numpy
import commands
import math
import time
from collections import OrderedDict
import logging
from operator import is_not
from functools import partial
import statistics
import pysam
import fnmatch

sys.path.insert(0, "/data/boehm/group/pipelines/ruffus")
from ruffus import *
from ruffus.proxy_logger import *
from ruffus.combinatorics import *
from ruffus.drmaa_wrapper import run_job, run_job_using_drmaa, error_drmaa_job

from PIL import Image
import string

parser = argparse.ArgumentParser(prog="methPipe", version="0.1.0", description="runs amplicon CpG methylation pipeline", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-ri", "--readIn", dest="readdir", action="store", default=False, help="input read folder")
parser.add_argument("-g", "--ref", dest="refpath", action="store", default=False, help="path to indexed reference genome")
parser.add_argument("-cg", "--cref", dest="convrefpath", action="store", default=False, help="path to converted and indexed reference genome")
parser.add_argument("-il", "--intList", dest="intList", action="store", default=None, help="target interval file")
parser.add_argument("-bl", "--blackList", dest="blackList", action="store", default=None, help="SNP black list")
parser.add_argument("-si", "--sampleInfo", dest="sampleInfo", action="store", default=None, help="sample sheet")
parser.add_argument("-wd", "--wdir", dest="wdir", action="store", default=False, help="output folder")
parser.add_argument("-vb", "--verbose", dest="verbose", action="store_true", help="more detailed log")
parser.add_argument("-bs", "--batchSize", dest="bsize", action="store",metavar="INT",type=int,default=10, help="number of samples to process in parallel")
parser.add_argument("-nt", "--numThr", dest="nthreads", action="store",metavar="INT",type=int,default=8, help="number of threads to use per sample")
parser.add_argument('--mbias', dest="mbias_ignore",action="store",default="auto",help="number of nucleotides with mbias to ignore during methylation extraction")
parser.add_argument("--touchOnly", dest="touchOnly", action="store_true", help="only touch files")
parser.add_argument("--target_tasks", dest="target_tasks", action="store",default=[], help="target tasks")
parser.add_argument("--forcedtorun_tasks", dest="forcedtorun_tasks", action="store",default=[], help="forced to run tasks")

args = parser.parse_args()


#setup central working directory
wdir=args.wdir
if not os.path.exists(wdir):
    os.makedirs(wdir)
os.chdir(wdir)
    
#setup logging
logger = logging.getLogger(__name__)
fhandler = logging.FileHandler(filename=os.path.join(wdir,'pipeline.log'), mode='a')
fformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fhandler.setFormatter(fformatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)

logger.debug(subprocess.check_output('echo $DRMAA_LIBRARY_PATH',shell=True))


import drmaa

#initiate 1 drmaa session for the whole pipeline
mySession=drmaa.Session()
mySession.initialize()


logger.info('Working directory is : ' + wdir)
logger.info(args)

readdir=args.readdir
#os.chdir(readdir)

libset2 = []
# Walk through directory
for dName, sdName, fList in os.walk(readdir):
    for fileName in fList:
        if fnmatch.fnmatch(fileName, "*fastq.gz"): # Match search string
            libset2.append(os.path.join(dName, fileName))

libset2_R1=filter(lambda x:'_R1.fastq.gz' in x, libset2)
libset2_R1.sort()
libset2_R2=filter(lambda x:'_R2.fastq.gz' in x, libset2)
libset2_R2.sort()
read_root=[ re.sub('_R1.fastq.gz','',R1f) for R1f in libset2_R1 ]
INfiles=list(zip(libset2_R1,libset2_R2))
    
logger.debug(INfiles) 

######################## paths to unconverted reference genome #######################################
if os.path.isfile(args.refpath):
    refG=args.refpath
    refGpath=os.path.dirname(refG)
    
elif os.path.isfile(os.path.join('/data/repository/organisms',(args.refpath+'_ensembl'),'genome_fasta/genome.fa')):
    refG=os.path.join('/data/repository/organisms',(args.refpath+'_ensembl'),'genome_fasta/genome.fa')
    refGpath=os.path.dirname(refG)
    
else:
    logger.info('Reference genome not recognized. Please check spelling and/or path.')
        
######################## paths to converted reference genome #######################################

if os.path.isdir(str(args.convrefpath)):
    #crefG=args.convrefpath
    #crefGpath=os.path.dirname(crefG)
    crefGpath=args.convrefpath
elif os.path.isdir(os.path.join('/data/repository/organisms',(str(args.refpath)+'_ensembl'),'BismarkIndex')):
    
    crefG=os.path.join('/data/repository/organisms',(str(args.refpath)+'_ensembl'),'BismarkIndex') #legacy from the WGBS pipeline
    crefGpath=os.path.join('/data/repository/organisms',(str(args.refpath)+'_ensembl'),'BismarkIndex')
elif not args.convrefpath:
    print('Converted reference genome not specified, will be generated by the pipeline.')
    args.convRef=True
else:
    logger.error('Converted reference genome not recognized. Please check spelling and/or path.')    
    
logger.debug('Reference genome is '+refG)
#logger.debug('Converted reference genome is '+crefG)
logger.debug('Convert genome? ' + 'False')#legacy from the WGBS pipeline

###poz files
auxdir=os.path.join(wdir,'aux_files')
pozF=os.path.join(auxdir,re.sub('.fa','.poz.gz',os.path.basename(refG)))
pozFsub=os.path.join(os.path.dirname(pozF),re.sub('poz.gz',re.sub('.bed','.poz.gz',os.path.basename(args.intList)),pozF))
logger.debug('Interval index file will be: '+pozFsub)

##################PATHS TO EXECUTABLES###############################################################
FQCpath='/package/FastQC-0.11.3'
cutpath='/package/cutadapt-1.9.1/bin'
mCTpath='/data/manke/repository/scripts/DNA_methylation/methylCtools'
tabpath='/package/tabix-1.2.1'
bwapath='/package/bwa-0.7.4/bin'
bismpath='/data/manke/repository/scripts/DNA_methylation/Bismark-master'
BTpath='/package/bowtie2-2.2.8/'
bmethpath='/data/manke/repository/scripts/DNA_methylation/bwa-meth-master_2016'
sampath='/package/samtools-1.3/bin'
Picpath='/package/picard-tools-1.136'
GATKpath='/package/GenomeAnalysisTK-3.5'
POMpath='/package/MethylDackel-0.3.0/bin'
bedpath='/package/bedtools2-2.25.0/bin'
BisSNPpath='/data/manke/repository/scripts/DNA_methylation/BisSNP-0.82.2.jar'
metipath='/data/manke/repository/scripts/DNA_methylation/metilene_v0.2-6'
Rpath='/package/R-3.3.1/bin'
prinpath='/data/boehm/sikora/tools/prinseq-lite-0.20.4'

######################################################################################################

#############################PIPELINE#################################################################################
#TRIM READS############################################################################################################
fqcout=os.path.join(wdir,'fastqc_cut')
cutout=os.path.join(wdir,'reads_cut')

        
@mkdir(cutout,os.path.join(cutout,'logs'))
@transform(input=INfiles,filter=suffix('_R1.fastq.gz'),output=['_R1.fastq.gz','_R2.fastq.gz'],output_dir=cutout) 
def adapter_trim(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    ii3 = os.path.join(fqcout,re.sub('_R1.fastq.gz','.R12.ct.txt',os.path.basename(ii1)))
    
    oo1 = output_files[0]
    oo2 = output_files[1]
    #hard-code threshold values
    ct1=2
    ct2=2
    
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(cutpath,'cutadapt') + ' -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 30  -n 5 -u ' + str(ct1) + ' -U ' + str(ct2) + '  -o ' + oo1 + ' -p ' + oo2 + ' ' + ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(cutout,"logs","%s.cut_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.cut_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = bshcmd,
                                      job_name          = 'cut_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("Cut_reads error: %s" % err)
            raise
        else:
            logger.info('Adapter trimming complete')

@transform(adapter_trim,suffix('_R1.fastq.gz'),['_prin_1.fastq.gz','_prin_2.fastq.gz'],output_dir=cutout)
def trim_BQ(infiles,outfiles):
    ii1=infiles[0]
    ii2=infiles[1]
    oo1=outfiles[0]
    oo2=outfiles[1]
    read_root=re.sub('_R1.fastq.gz','',os.path.basename(ii1))
    uzcmd1='zcat -v '+ ii1 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1)))
    uzcmd2='zcat -v '+ ii2 + ' > ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2)))
    bshcmd='perl '+ os.path.join(prinpath,'prinseq-lite.pl') + ' -fastq ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' -fastq2 ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' -out_good ' + os.path.join(cutout,re.sub('_1.fastq.gz','',os.path.basename(oo1))) +' -trim_qual_right 24 -trim_qual_type mean -trim_qual_window 6 -trim_qual_step 3 -min_len 50 -ns_max_p 10 -out_bad null'
    zcmd1='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' > ' + oo1
    zcmd2='gzip -c '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2))) + ' > ' + oo2
    clcmd='rm -v '+ os.path.join(cutout,re.sub('.gz','',os.path.basename(ii1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(ii2))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo1))) + ' ' + os.path.join(cutout,re.sub('.gz','',os.path.basename(oo2)))
    cmd_all=';'.join([uzcmd1,uzcmd2,bshcmd,zcmd1,zcmd2,clcmd])
    logger.info(cmd_all)           
    with open(os.path.join(cutout,"logs","%s.BQtrim_reads.out" % read_root),'w+') as stdoutF, open(os.path.join(cutout,"logs","%s.BQtrim_reads.err" % read_root),'w+') as stderrF:
        try:
        
            stdout_res, stderr_res  = run_job(cmd_str   = cmd_all,
                                      job_name          = 'BQtrim_reads',
                                      logger            = logger,
                                      drmaa_session     = mySession,
                                      run_locally       = False,
                                      working_directory = os.getcwd(),
                                      job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        except error_drmaa_job as err:
            logger.error("BQtrim_reads error: %s" % err)
            raise
        else:
            logger.info('Base quality trimming complete')    
        

@mkdir(fqcout,os.path.join(fqcout,'logs'))            
@transform(trim_BQ,suffix('_1.fastq.gz'),output=['_1.zip','_1.html','_2.zip','_2.html'],output_dir=fqcout)
def postTrim_fqc(input_files,output_files):
    ii1 = input_files[0]
    ii2 = input_files[1]
    read_root=re.sub('_prin_1.fastq.gz','',os.path.basename(ii1))
    bshcmd=os.path.join(FQCpath,'fastqc ')+' --outdir ' + fqcout + ' -t 8 '+ ii1 + ' ' + ii2
    logger.info(bshcmd)       
    with open(os.path.join(fqcout,"logs","%s.post_fqc.out" % read_root),'w+') as stdoutF, open(os.path.join(fqcout,"logs","%s.post_fqc.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = bshcmd,
                                          job_name          = 'post_fqc',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=8')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Post_trim_fastqc error: %s" % err)
            raise
        else:
            logger.info('Post trim fastqc complete')           

#CONVERT AND MAP READS############################################################################################
bamoutO=os.path.join(wdir,'bams_Bismark')
################################################################################################################### 

@mkdir(bamoutO,os.path.join(bamoutO,'logs'))
@transform(trim_BQ,suffix('_1.fastq.gz'),'_pe.bam',output_dir=bamoutO)
def map_reads(input_files,output_file):
    ii1 = input_files[0]
    ii2 = input_files[1]
    oo = output_file
    read_root=re.sub('_1.fastq.gz','',os.path.basename(ii1))
    mapcmd= os.path.join(bismpath,'bismark') + ' -p ' + str(args.nthreads) + ' --non_directional --dovetail --temp_dir /data/extended --path_to_bowtie /package/bowtie2-2.2.8 --output_dir ' + bamoutO + ' --basename ' + read_root + ' --genome_folder ' + crefGpath + ' -1 ' + ii1 + ' -2 ' + ii2
    logger.info(mapcmd)
    with open(os.path.join(bamoutO,"logs","%s.readmap.out.log" % read_root),'w') as stdoutF, open(os.path.join(bamoutO,"logs","%s.readmap.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mapcmd,
                                          job_name          = 'BSmap',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus='+str(args.nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("Map_reads error: %s" % err)
            raise
        else:
            logger.info('Mapping complete')
    return



################## BAM POSTPROCESSING #########################################################
#bam postprocessing; bwa-meth is sorted and indexed, Bismark and methylCtools require these steps

@transform(map_reads,suffix('_pe.bam'),'.sorted.bam',output_dir=bamoutO)
def sort_bam(input_file,output_file):
    ii = input_file
    oo = output_file
    read_root=re.sub('_pe.bam','',os.path.basename(ii))
    cmd=os.path.join(sampath,'samtools') + ' sort -T ' + os.path.join('/data/extended',read_root) + ' -m 6G -@ ' + str(args.nthreads) + ' -o ' + oo + ' ' + ii
    logger.info(cmd)
    with open(os.path.join(bamoutO,"logs","%s.bamsort.out" % read_root),'w+') as stdoutF, open(os.path.join(bamoutO,"logs","%s.bamsort.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd,
                                          job_name          = 'bamsort',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mincpus=' + str(args.nthreads))
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bam sorting error: %s" % err)
            raise
        else:
            logger.info('Bam sorting complete')
            
            
            
@transform(sort_bam,suffix('.sorted.bam'),'.sorted.bam.bai',output_dir=bamoutO) 
def index_bam(input_file,output_file):
    ii = input_file
    oo = output_file
    read_root=re.sub('.sorted.bam','',os.path.basename(ii))
    cmd=os.path.join(sampath,'samtools') + ' index ' + ii
    logger.info(cmd)
    with open(os.path.join(bamoutO,"logs","%s.bam_index.out" % read_root),'w+') as stdoutF, open(os.path.join(bamoutO,"logs","%s.bam_index.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd,
                                          job_name          = 'bam_index',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Bam indexing error: %s" % err)
            raise
        else:
            logger.info('Bam indexing complete')
    
@follows(index_bam)    
@transform(sort_bam,suffix('.sorted.bam'),'.RGi.bam') #for compatibility with GATK; depth of coverage
def add_ReadGroupInfo(input_file,output_file):
    ii=input_file
    oo=output_file
    read_root=re.sub('.sorted.bam','',os.path.basename(ii))
    read_str=subprocess.check_output(os.path.join(sampath,'samtools') + ' view ' + ii + '| head -n 1',shell=True)
    logger.debug(read_str)
    PL=read_str.split(":")[0]
    PU=read_str.split(":")[2]
    cmd_addRGI='java -Xmx20g -Djava.io.tmpdir=/data/extended -jar ' +  os.path.join(Picpath,'picard.jar') + ' AddOrReplaceReadGroups I=' + ii + ' O=' + oo + ' SORT_ORDER=coordinate CREATE_INDEX=true RGPL=' + PL + ' RGSM=' + read_root + ' RGLB=' + read_root + ' RGPU=' + PU + ' VALIDATION_STRINGENCY=SILENT; ln -f -s '+ re.sub('RGi.bam','RGi.bai',oo) + ' ' + re.sub('RGi.bam','RGi.bam.bai',oo)
    logger.info(cmd_addRGI)
    with open(os.path.join(bamoutO,"logs","%s.add_RGi.out" % read_root),'w+') as stdoutF, open(os.path.join(bamoutO,"logs","%s.add_RGi.err" % read_root),'w+') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str           = cmd_addRGI,
                                          job_name          = 'addRGi',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except error_drmaa_job as err:
            logger.error("Adding read group info error: %s" % err)
            raise
        else:
            logger.info('Adding read group info complete')
            
######################################################################################################################  
#RUN VARIOUS COVERAGE AND QUALITY METRICS#######################################################################
metout=os.path.join(wdir,'QC_metrics')

#methylation bias
@mkdir(metout,os.path.join(metout,'logs'))
@transform(add_ReadGroupInfo,suffix('.RGi.bam'),'.Mbias.txt',output_dir=metout)
def calc_Mbias(input_file,output_file):
    ii=input_file
    oo=output_file
    oos=re.sub('.txt','',oo)
    read_root=re.sub('.bam','',os.path.basename(ii))
    Mb_cmd=os.path.join(POMpath,'MethylDackel') + ' mbias --txt --keepDupes -@ '+ str(args.nthreads) + ' ' + refG + ' ' + ii +' ' + oos +' > ' + oo #+ '.txt'
    logger.info(Mb_cmd)
    with open(os.path.join(metout,"logs","%s.mbias.out" % read_root),'w') as stdoutF, open(os.path.join(metout,"logs","%s.mbias.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = Mb_cmd,
                                          job_name          = 'mbias',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem-per-cpu=10000 --mincpus='+str(args.nthreads))
            
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))


        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("Methylation bias error: %s" % err)
            raise
        else:
            logger.info('Methylation bias calculation complete')

#flagstat -> mapping rate

@mkdir(metout,os.path.join(metout,'logs'))
@transform(add_ReadGroupInfo,suffix('.RGi.bam'),'.flagstat',output_dir=metout)
def get_flagstat(input_file,output_file):
    ii=input_file
    oo=output_file
    read_root=re.sub('.RGi.bam','',os.path.basename(ii))
    cmd=os.path.join(sampath,'samtools') + ' flagstat ' + ii +' > ' + oo 
    logger.info(cmd)
    with open(os.path.join(metout,"logs","%s.flagstat.out" % read_root),'w') as stdoutF, open(os.path.join(metout,"logs","%s.flagstat.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd,
                                          job_name          = 'fstat',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo ' )
            
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))


        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("Flagstat error: %s" % err)
            raise
        else:
            logger.info('Flagstat calculation complete')



#GATK depth of coverage
if ( args.intList ): 
    int_dest=re.sub('.bed','.mean.doc.sample_summary',os.path.basename(args.intList))
    @mkdir(metout,os.path.join(metout,'logs'))
    @transform(add_ReadGroupInfo,suffix('RGi.bam'),int_dest,output_dir=metout)#
    def depth_of_cov(input_file,output_file):
        ii=input_file
        oos=output_file 
        oos2=oos.replace('.sample_summary', '') 
        read_root=re.sub('.RGi.bam','',os.path.basename(ii))
        #OUTlist2=oos2[2:]
        cmd_all='java -Xmx50g -Djava.io.tmpdir=/data/extended -jar '+ os.path.join(GATKpath,'GenomeAnalysisTK.jar')+' -R '+ refG + ' -T DepthOfCoverage -o ' + oos2 + ' -I ' + ii + ' -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50  -omitBaseOutput -mmq 10 --partitionType sample -L ' + args.intList 
        logger.info(cmd_all)
        with open(os.path.join(metout,"logs","%s.depth_cov.out.log" % read_root),'w') as stdoutF, open(os.path.join(metout,"logs","%s.depth_cov.err.log" % read_root),'w') as stderrF:
            try:
                stdout_res, stderr_res  = run_job(cmd_str       = cmd_all,
                                          job_name          = 'depth_cov',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo --mem=50000')
                stdoutF.write("".join(stdout_res))
                stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
            except Exception as err:
                logger.error("Depth of coverage error: %s" % err)
                raise
            else:
                logger.info('Depth of coverage calculation complete')
    
    
#conversion rate: currently from fastq, implement phiX control!
@mkdir(metout,os.path.join(metout,'logs'))
@transform(trim_BQ,suffix('_1.fastq.gz'),'.conv.rate.txt',output_dir=metout) #it's ok to have both reads summarized in 1 file
def conv_rate(input_files,output_file):
    ii1=input_files[0]
    ii1sub=re.sub('_1.fastq.gz','',ii1)
    oo=output_file
    read_root=os.path.basename(ii1sub)
    CR_cmd='/data/boehm/group/pipelines/BS_amplicon_seq/v0.1.0/conversionRate_prin_KS.sh '+ ii1sub + ' ' + oo
    logger.info(CR_cmd)
    with open(os.path.join(metout,"logs","%s.conv_rate.out.log" % read_root),'w') as stdoutF, open(os.path.join(metout,"logs","%s.conv_rate.err.log" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = CR_cmd,
                                          job_name          = 'conv_rate',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("Conversion rate error: %s" % err)
            raise
        else:
            logger.info('Conversion rate calculation complete')
        
################################################################################################################## 
###final QC report: only if started from reads or bam files
if not os.path.exists(metout):
        os.makedirs(metout)
os.chdir(metout)
@merge(input=[depth_of_cov,conv_rate,calc_Mbias,get_flagstat],output='QC_report.pdf')
def produce_report(input_file,output_file):
    ii=input_file
    oo=output_file
    cmd=os.path.join(Rpath,'Rscript -e "rmarkdown::render(\'"/data/boehm/group/pipelines/BS_amplicon_seq/v0.1.0/BSampli_QC_report_template.Rmd"\',params=list(QCdir=\'"' + metout +'"\' ),output_file =\'"'+ os.path.join(metout,'QC_report.pdf"\'')+')"')
    logger.info(cmd)
    logger.debug(subprocess.check_output(cmd,shell=True))
    logger.info('Generating QC report complete')
    
#####EXTRACT METHYLATION COUNTS  ########################################################################
mextout=os.path.join(wdir,'methXT')
auxdir=os.path.join(wdir,'aux_files') 
@mkdir(auxdir)
@transform(refG,suffix('.fa'),['.poz.P.txt','.poz.M.txt'],output_dir=auxdir)
def get_poz(input_file,output_files):
    ii=input_file
    oo1=output_files[0]
    oo2=output_files[1]
    logger.info('Preparing an index of CG positions')
    cmd_fapos=os.path.join(mCTpath,'methylCtools') + ' fapos ' + refG + ' - | ' + os.path.join(tabpath,'bgzip') + ' > ' + pozF
    idx_fapos=os.path.join(tabpath,'tabix') + ' -s 1 -b 2 -e 2 ' + pozF
    #awk column selection adapted to ensembl genome stored under organisms
    prep_poz='gzip -dc ' + pozF + '| grep \'+\'  | awk \'{print $1, $6, $6+1, $7, $8, $9}\' - | tr " " "\t" > ' + re.sub('gz','P.txt',pozF) + ' ; gzip -dc ' + pozF + ' | grep \'-\'  | awk \'{print $1, $2, $2+1, $3, $4, $5}\' - | tr " " "\t" > ' + re.sub('gz','M.txt',pozF)
    sel_int=os.path.join(tabpath,'tabix ') + pozF + ' -R ' + args.intList + ' | ' + os.path.join(tabpath,'bgzip') + ' > ' + pozFsub 
    re_idx=os.path.join(tabpath,'tabix') + ' -s 1 -b 2 -e 2 ' + pozFsub
    cmd_all=';'.join([cmd_fapos,idx_fapos,prep_poz,sel_int,re_idx])
    logger.info(cmd_all)
    try:
        logger.debug(subprocess.check_call(cmd_all,shell=True))
    except Exception as poze:
        logger.error("CpG index preparation error: %s" % poze)
        raise IOError
    else:
        logger.info('CpG index preparation complete')

    
@follows('calc_Mbias','get_poz',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))
@transform(add_ReadGroupInfo,suffix('.RGi.bam'),'.CG.call.gz',output_dir=mextout)
def methyl_extract(input_file,output_file):
    ii = input_file
    oo = output_file
    read_root=re.sub('.RGi.bam','',os.path.basename(ii))
    if len(args.mbias_ignore) < 3:
        mCT_cmd=os.path.join(mCTpath,'methylCtools') + ' bcall --mapQ 0 --snv --trimPE --skipend ' + str(args.mbias_ignore) + ' --genomic --zero ' + pozFsub + ' ' + ii + ' - | ' +  os.path.join(tabpath,'bgzip') + ' > ' + oo 
    elif args.mbias_ignore=="auto":
        OT=pandas.read_table(os.path.join(metout,"logs",read_root)+'.mbias.err',sep=' ',header=None)[4]
        OT=OT.str.split(pat=',').tolist()
        OB=pandas.read_table(os.path.join(metout,"logs",read_root)+'.mbias.err',sep=' ',header=None)[6]
        OB=OB.str.split(pat=',').tolist()
        auto_ignore=max([OT[0][0],OT[0][2],OB[0][0],OB[0][2]])
        auto_ignore_corr=min(10,auto_ignore)
        mCT_cmd=os.path.join(mCTpath,'methylCtools') + ' bcall --mapQ 0 --snv --trimPE --skipend ' + auto_ignore_corr + ' --genomic --zero ' + pozFsub + ' ' + ii + ' - | ' +  os.path.join(tabpath,'bgzip') + ' > ' + oo 
    logger.info(mCT_cmd)
    with open(os.path.join(mextout,"logs","%s.mCT_extract.out" % read_root),'w') as stdoutF, open(os.path.join(mextout,"logs","%s.mCT_extract.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = mCT_cmd,
                                          job_name          = 'mCT_extract',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("Methylation extraction error: %s" % err)
            raise
        else:
            logger.info('Methylation extraction complete')

@follows('get_poz',mkdir(mextout),mkdir(os.path.join(mextout,'logs')))    
@transform(methyl_extract,suffix('.CG.call.gz'),'.CpG.filt2.bed',output_dir=mextout) 
def CpG_filt(input_file,output_file):
    ii = input_file
    oo = output_file 
    read_root=re.sub('.CG.call.gz','',os.path.basename(ii))
    gz_cmd='gzip -dc '+ ii + ' > '+ re.sub('.gz','',ii)
    filt_cmd=os.path.join(Rpath,'Rscript') + ' --no-save --no-restore /data/boehm/group/pipelines/BS_amplicon_seq/v0.1.0/BSampli.mCT.filt.R ' + mextout + ' ' + re.sub('.gz','',ii) + ' ' + pozFsub
    clean_cmd='rm -v ' + re.sub('.gz','',ii)
    cmd_all=';'.join([gz_cmd,filt_cmd,clean_cmd])
    logger.info(cmd_all)
    with open(os.path.join(mextout,"logs","%s.CpG_filt.out" % read_root),'w') as stdoutF, open(os.path.join(mextout,"logs","%s.CpG_filt.err" % read_root),'w') as stderrF:
        try:
            stdout_res, stderr_res  = run_job(cmd_str       = cmd_all,
                                          job_name          = 'CpG_filt',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
            stdoutF.write("".join(stdout_res))
            stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
        except Exception as err:
            logger.error("CpG filtering error: %s" % err)
            raise
        else:
            logger.info('CpG filtering complete') 
            
####RUN SINGLE CYTOSINE STATISTICS IN R VIA LIMMA#########################################################
if args.sampleInfo :
    CpGstat_out=os.path.join(wdir,'singleCpG_stats_limma')
    if not os.path.exists(CpGstat_out) or not os.path.exists(os.path.join(CpGstat_out,'logs')):
        os.makedirs(CpGstat_out)
        os.makedirs(os.path.join(CpGstat_out,'logs'))
    os.chdir(CpGstat_out)
    @merge(input=CpG_filt,output=[os.path.join(CpGstat_out,'singleCpG.RData'),os.path.join(CpGstat_out,'limdat.LG.RData')])
    def CpG_stats(input_files,output_file):
        ii = os.path.dirname(input_files[0])
        oo = output_file
        Rstat_cmd=os.path.join(Rpath,'Rscript') + ' --no-save --no-restore /data/boehm/group/pipelines/BS_amplicon_seq/v0.1.0/BSampli.singleCpGstats.limma.R ' + CpGstat_out + ' ' + args.sampleInfo + ' ' + ii
        logger.info(Rstat_cmd)
        with open(os.path.join(CpGstat_out,"logs","singleCpG_stats.out" ),'w') as stdoutF, open(os.path.join(CpGstat_out,"logs","singleCpG_stats.err"),'w') as stderrF:
            try:
                stdout_res, stderr_res  = run_job(cmd_str       = Rstat_cmd,
                                          job_name          = 'sCpG',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
                stdoutF.write("".join(stdout_res))
                stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
            except Exception as err:
                logger.error("Single CpG stats error: %s" % err)
                raise
            else:
                logger.info('Single CpG stats calculation complete')

####RUN INTERVAL AGGREGATE STATISTICS IN R ####################################################
if args.sampleInfo:
    int_dest=re.sub('.bed','.aggCpG.RData',os.path.basename(args.intList)) 
    intStat_out=os.path.join(wdir,'aggregate_stats_limma')
    @mkdir(intStat_out,os.path.join(intStat_out,'logs'))
    @transform(CpG_stats,suffix('singleCpG.RData'),int_dest,output_dir=intStat_out)
    def intAgg_stats(input_files,output_files):
        ii = os.path.join(CpGstat_out,input_files[1])
        oo = output_files
        Rcmd=os.path.join(Rpath,'Rscript') + ' --no-save --no-restore /data/boehm/group/pipelines/BS_amplicon_seq/v0.1.0/BSampli.interval_stats.limma.R ' + intStat_out + ' ' + args.intList + ' ' + ii + ' ' + args.sampleInfo 
        logger.info(Rcmd)
        with open(os.path.join(intStat_out,"logs","interval_stats.out" ),'w') as stdoutF, open(os.path.join(intStat_out,"logs","interval_stats.err"),'w') as stderrF:
            try:
                stdout_res, stderr_res  = run_job(cmd_str       = Rcmd,
                                          job_name          = 'agg_stats',
                                          logger            = logger,
                                          drmaa_session     = mySession,
                                          run_locally       = False,
                                          working_directory = os.getcwd(),
                                          job_other_options = '-p bioinfo')
                stdoutF.write("".join(stdout_res))
                stderrF.write("".join(stderr_res))

        # relay all the stdout, stderr, drmaa output to diagnose failures
            except Exception as err:
                logger.error("Interval stats error: %s" % err)
                raise
            else:
                logger.info('Interval stats calculation complete')
        

#############################################################################################

#####main

if __name__ == '__main__':
    with open(os.path.join(wdir,"pipelineGraph.png"),'w') as pipeGraph:
        pipeline_printout_graph(stream=pipeGraph,output_format='png',pipeline_name='WGBS',target_tasks=args.target_tasks)
    with open (os.path.join(wdir,"pipelinePrint.txt"),'w') as pipePrint:
        pipeline_printout(verbose_abbreviated_path=0,output_stream=pipePrint,target_tasks=args.target_tasks)    

    pipeline_run(touch_files_only=args.touchOnly,multiprocess=args.bsize,target_tasks=args.target_tasks,forcedtorun_tasks=args.forcedtorun_tasks,logger=logger)
    mySession.exit()           
