######################################################gRNAinfo ########################################
################################## Generating analysis files for input gRNAs ##########################
#######################################################################################################

# To run the script: set the working directory to the location of the data files that are used
# for the analysis. Next, define the locations of interest. After adjusting these parameters, 
# run Code to generate output files (5 in total). 

# Load libraries
library(rtracklayer)
library(Gviz)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggsci)
library(genomation)
library(chromoMap)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(reshape2)
library(dplyr)
library(ggplot2)
library(plyr)

# Set working directory to folder where input files are located
setwd("")

########################### Read in the gRNA locations for the analysis ####################
# read in .txt file with gRNA location information
gRNAloc<-read.table("gRNAloc.txt")
# assign columns with the chromosome and the location of the first base of the gRNA
# for analysis of only a small list of locations, a vector with locations and chromosome numbers can also be inserted here
loclist<-as.vector(gRNAloc$V4)
chrlist<-as.vector(gRNAloc$V3)

# Define the window around the gRNA location for which the information is returned (default: 2000)
# eg 2000 means a window of 2 kb left and 2 kb right to the location of interest (4 kb window)
range <- 2000



#######################################################################################################
############################## import of input datafiles ##############################################
#######################################################################################################

# import DNAseI data 
DNaseI<-import("ENCFF128BPC.bigwig")
# import DAMID data
DAMID<-import("RPE_LMNB1-5kb-combined_AD.bed")
    values(DAMID)<-DataFrame(score="LAD")
    # ADD ranges to remove gaps
    DAMgap<-gaps(DAMID)
    values(DAMgap)<-DataFrame(score="iLAD")
    DAMID_full<-c(DAMgap,DAMID)
# import ChromHMM data
chromHMM<-import("Rep1_10_dense.bed")
# import RNA seq data
RNAseq<-read.delim("RPE_RNA_seq.txt")
    RNAseq$Mean_expression<-as.numeric(sub(",", ".", RNAseq$Mean_expression, fixed = TRUE))
    chromo<-RNAseq$chromosome_name
    chromo<-paste("chr",chromo)
    chromo<-gsub(" ","",chromo,fixed=TRUE)
    start<-RNAseq$start_position
    w<-RNAseq$end_position-RNAseq$start_position
    expr<-RNAseq$Mean_expression
    gene<-RNAseq$external_gene_id
    RNAseq$Mean_expression<-as.numeric(RNAseq$Mean_expression)
    generange<-GRanges(seqnames = chromo,ranges = IRanges(start,width=w),strand = "*",expression=expr,gene=gene)
    gapsrange<-gaps(generange)
    rangeRNA<-c(generange,gapsrange)
    RNAseq<-sort(rangeRNA)
# import H3K9me3
H3K9me3<-import("RPE1_H3K9me3_ChIPSeq_GSM3105086.bw")
    
    
    
############################################### START ANALYSIS #################################################

################################################################################################################
############################## 1.DNaseI signal analysis at gRNA locations ######################################
################################################################################################################

# Function for analysis of DNAseI data
DNaseI_env<- function(location,chrom,range) {
  w<-(range*2)+1
  start<-location-range
  gr1<-GRanges(seqnames = chrom,ranges = IRanges(start,width=w),strand = "*",name=11,score=1,itemRgb="0000")
  hitlist<-DNaseI[queryHits(findOverlaps(DNaseI,gr1)),]
  sc<-as.numeric(hitlist$score)
  s<-summary(sc)
  max<-as.numeric(s[6])
  mean<-as.numeric(s[4])
  table<-cbind(location,max,mean)
  print(table)
}  

# Generate output file
DNaseIlist<-mapply(DNaseI_env,loclist,chrlist,range,SIMPLIFY = FALSE)
DNaseIlist2 = do.call(rbind,DNaseIlist)
DNaseIlist2<-as.data.frame(DNaseIlist2)
# Save output file
write.table(DNaseIlist2,"1.DNaseI_hypersensitivity.txt",sep="\t",row.names=FALSE)

################################################################################################################
############################## 2.DAMID signal analysis at gRNA locations #######################################
################################################################################################################

# Function for analysis of DAMID data
DAMID_env<- function(location,chrom,range) {
  w<-(range*2)+1
  start<-location - range
  gr1<-GRanges(seqnames = chrom,ranges = IRanges(start,width=w),strand = "*",score="nodata")
  findOverlaps(DAMID_full,gr1)
  hitlist<-DAMID_full[queryHits(findOverlaps(DAMID_full,gr1)),]
  combi<-(disjoin(c(gr1@ranges,hitlist@ranges)))
  combi<-combi[-NROW(combi),]
  combi<-combi[-1,]
  score<-as.character(hitlist$score)
  v1<-width(combi)
  proportion<-v1/sum(v1)
  table<-cbind(location,score,proportion)
  print(table)
}

# Generate output file
DAMIDlist<-mapply(DAMID_env,loclist,chrlist,range,SIMPLIFY = FALSE)
DAMIDlist = do.call(rbind,DAMIDlist)
DAMIDlist<-as.data.frame(DAMIDlist)
# Save output file
write.table(DAMIDlist,"2.DAMID.txt",sep="\t",row.names=FALSE)

###############################################################################################################
############################## 3.ChromHMM signal analysis at cut sites ########################################
###############################################################################################################

# Function for analysis of ChromHMM data
chromHMM_env<- function(location,chrom,range) {
  w<-(range*2)+1
  start<-location-range
  gr1<-GRanges(seqnames = chrom,ranges = IRanges(start,width=w),strand = "*",name=11,score=1,itemRgb="0000")
  findOverlaps(chromHMM,gr1)
  hitlist<-chromHMM[queryHits(findOverlaps(chromHMM,gr1)),]
  combi<-(disjoin(c(gr1@ranges,hitlist@ranges)))
  combi<-combi[-NROW(combi),]
  combi<-combi[-1,]
  v1<-width(combi)
  ChromHMMdata<-as.numeric(hitlist$name)
  proportion<-v1/sum(v1)
  t1<-cbind(location,ChromHMMdata,proportion)
  print(t1)
}

# Generate output file
chromHMM_list<-mapply(chromHMM_env,loclist,chrlist,range,SIMPLIFY = FALSE)
chromHMM_list = do.call(rbind,chromHMM_list)
chromHMM_list<-as.data.frame(chromHMM_list)
# Save output file
write.table(chromHMM_list,"3.chromHMM_score.txt",sep="\t",row.names=FALSE)



###############################################################################################################
############################## 4.RNAseq signal analysis at gRNA locations #####################################
###############################################################################################################

# Function for analysis of RNAseq data
RNAseq_env<- function(location,chrom,range) {
  w<-(range*2)+1
  start<-location-range
  gr1<-GRanges(seqnames = chrom,ranges = IRanges(start,width=w),strand = "*",expression="nvt",gene="nvt")
  findOverlaps(RNAseq,gr1)
  hitlist<-RNAseq[queryHits(findOverlaps(RNAseq,gr1)),]
  start<-start(hitlist)
  end<-end(hitlist)
  seqnames<-as.character(seqnames(hitlist))
  expression<-hitlist$expression
  gene<-hitlist$gene
  gene<-as.character(levels(gene))[gene]
  ref<-as.character(paste(chrom,location,sep=":"))
  t1<-cbind(ref,location,chrom,start,end,seqnames,expression,gene)
  print(t1)
}

# Generate output file
RNAlist<-mapply(RNAseq_env,loclist,chrlist,range,SIMPLIFY = FALSE)
RNAlist = do.call(rbind,RNAlist)
RNAlist<-as.data.frame(RNAlist)
# Save output file
write.table(RNAlist,"4.RNAseq.txt",sep="\t",row.names=FALSE)



###############################################################################################################
############################## 5.H3K9me3 signal analysis at gRNA locations ####################################
###############################################################################################################

# Function for analysis of H3K9me3 data
H3K9me3_env<- function(location,chrom,range) {
  w<-(range*2)+1
  start<-location-range
  gr1<-GRanges(seqnames = chrom,ranges = IRanges(start,width=w),strand = "*",score="nodata")
  findOverlaps(H3K9me3,gr1)
  hitlist<-H3K9me3[queryHits(findOverlaps(H3K9me3,gr1)),]
  meanpeak<-mean(hitlist$score)
  maximum<-max(hitlist$score)
  lengthlist<-length(hitlist)
  weightav<-(sum(hitlist$score * width(hitlist)))/(sum(width(hitlist)))
  t1<-cbind(location,maximum,lengthlist,meanpeak,weightav)
  print(t1)
}

# Generate output file
H3K9me3list<-mapply(H3K9me3_env,loclist,chrlist,range,SIMPLIFY = FALSE)
H3K9me3list = do.call(rbind,H3K9me3list)
H3K9me3list<-as.data.frame(H3K9me3list)
# Save output file
write.table(H3K9me3list,"5.H3K9me3_weightmean.txt",sep="\t",row.names=FALSE)
