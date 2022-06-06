######################################################gRNAfilter ########################################
################################## Filter data files for selection of gRNAs ##########################

# After generation of data files for all locations of interest for each environment feature,
# this gRNAfilter script is used for filtering of the gRNAs. For filtering, multiple parameters need to 
# be adjusted throughout the script, depending on the characteristics of interest. 

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
library(chromoMap)

# read in .txt file with gRNA location information
gRNAloc<-read.table("gRNAloc.txt")
# assign columns with the chromosome and the location of the first base of the gRNA
# for analysis of only a small list of locations, a vector with locations and chromosome numbers can also be inserted here
chrlist<-as.vector(gRNAloc$V4)
loclist<-as.vector(gRNAloc$V3)

# load the datasets from the single_cutter_gRNAinfo script
# DNAseI
DNaseI_file<-read.delim("1.DNaseI_hypersensitivity.txt", header = TRUE, sep = "\t")
# DAMID
DAMID_file<-read.delim("2.DAMID.txt", header = TRUE, sep = "\t")
# ChromHMM
ChromHMM_file<-read.delim("3.chromHMM_score.txt", header = TRUE, sep = "\t")
# RNA-seq
RNAseq_file<-read.delim("4.RNAseq.txt", header = TRUE, sep = "\t")
# H3K9me3
H3K9me3_file<-read.delim("5.H3K9me3_weightmean.txt", header = TRUE, sep = "\t")



#####################################################################################################################
############################## ChromHMM signal analysis at cut sites #############################################
####################################################################################################################

##### Filtering step 1: selection on chromHMM state
# First filter the data for locations with a specific state. (default is state 6) 
select_state <- ChromHMM_file %>% filter(ChromHMM_file[,2] == 6)
# Then filter the % of this state that is required to be present (default is 0.8, so 40% of the window)
state_perc <- select_state %>% filter(select_state[,3] > 0.3)

# After filtering, the percentage of the chosen state for each location is plotted for visualisation
ggplot(state_perc,aes(x=as.character(location),y=state_perc[,3]))+
  geom_bar(stat="identity",fill="blue")+
  theme(axis.text.x = element_text(angle=70)) +
  ggtitle("title") +
  xlab("Location") +
  ylab("State quantification")+
  scale_y_continuous(limits =c(0,1))+
  coord_flip()



#####################################################################################################################
############################## DNaseI signal analysis at cut sites ##################################################
####################################################################################################################

##### Filter step 2: selection on DNaseI peak maximum (column 2) or mean (column 3)
# merge the data from DNaseI with this subset of guides that were filtered based on ChromHMM score
ChromHMM_sel<-merge(x=DNaseI_file,y=state_perc,by="location",all.y=TRUE)
# set the desired threshold for maximal DnaseI peak 2 kb around the break site (default: filter on peak maximum, peak height = 2)
ChromHMM_DNAseI_sel<- ChromHMM_sel %>% filter(ChromHMM_sel[,2] > 2)

# plot the DNAseI peaks max
ggplot(ChromHMM_DNAseI_sel,aes(y=max,x=as.character(location)))+
  geom_bar(stat="identity",fill="purple")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ggtitle("DNaseI") +
  scale_x_discrete(name="Location") +
  scale_y_continuous(name="DNaseI value",limits =c(0,21))

# After filtering, retrieve a new list with the gRNA locations that were filtered
gRNAloc_filtered<-unique(ChromHMM_DNAseI_sel$location)
gRNAselection<-gRNAloc[gRNAloc$V4 %in% gRNAloc_filtered,]

loclist<-gRNAselection$V4
loclist<-as.vector(loclist)

chrlist<-gRNAselection$V3
chrlist<-as.vector(chrlist)

# optional: save the current filtered list as .txt, so the different filter steps can be traced back later
write.table(gRNAselection,"gRNAselection_chromHMM_DNaseI.txt",sep="\t",row.names=FALSE)



#####################################################################################################################
############################## DAMID signal analysis at cut sites ##################################################
####################################################################################################################

##### filter step 3: DAMID
# The DAMID file has already been loaded.
# Use the DAMIDlist file for filtering (default: 0.5, meaning 50% in an iLAD region)
DAMIDlist2 <- DAMIDlist %>% filter(DAMIDlist[,3] > 0.5)
DAMIDlist3<-unique(DAMIDlist2$location)

# Generate a table with the DAMID data for the gRNA selection
Filteredlist<-gRNAselection[gRNAselection$V4 %in% DAMIDlist3,]
# Generate a table with the ChromHMM data for the selection
Filteredlist2<-ChromHMM_file[ChromHMM_file$location %in% DAMIDlist3,]

# optional: save the current filtered list as .txt, so the different filter steps can be traced back later
write.table(Filteredlist,"gRNAselection_DAMID.txt",sep="\t",row.names=FALSE)



#####################################################################################################################
############################## RNA-seq signal analysis at cut sites ##################################################
####################################################################################################################

# for the analysis we use the RNAseq_file 

# based on filtering preferences, remove either genic or non-genic
RNAseq_sel<- RNAseq_file %>% filter(expression == "NA")

# for genic locations, breaks in introns were manually selected by checking RefSeq data.

# plot the results
ggplot(RNAseq_sel, aes(x=as.character(location), y=expression)) + 
  geom_bar(stat = "identity",fill="darkgrey") +
  ggtitle("Expression levels of genes at break sites")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  scale_x_discrete(name="Location") +
  scale_y_continuous(name="Gene expression")



#####################################################################################################################
############################## H3K9me3 signal analysis at cut sites ##################################################
####################################################################################################################

# filter the selection of gRNAs based on highest H3K9me3 from the H3K9me3_file

# Set filtering of the weighted peak height in the windown around the break (default: 5)
H3K9me3sel <- H3K9me3_file %>% filter(H3K9me3_file[,5] >1)

# plot results for weighted average
ggplot(H3K9me3sel,aes(y=weightav,x=as.character(location)))+
  geom_bar(stat="identity",fill="purple")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ggtitle("H3K9me3 weighted average") +
  scale_x_discrete(name="gRNA location") +
  scale_y_continuous(name="weighted average H3K9me3 normalized fragment coverage within window")

# plot for maximums
ggplot(H3K9me3sel,aes(y=maximum,x=as.character(location)))+
  geom_bar(stat="identity",fill="purple")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  ggtitle("H3K9me3 maximal peak height") +
  scale_x_discrete(name="gRNA location") +
  scale_y_continuous(name="maximal H3K9me3 peak height within window")

H3K9me3sel
# After filtering, retrieve a new list with the gRNA locations that were filtered
gRNAloc_H3K9me3_filtered<-unique(H3K9me3sel$location)
gRNAselection_H3K9me3<-gRNAloc[gRNAloc$V4 %in% gRNAloc_H3K9me3_filtered,]

# Generate the loclist and chrlist vectors for the filtered selection. With new assignments of these parameters, the gRNAinfo script can be
# executed again, to get information on all variables for the filtered gRNA selection.
loclist<-gRNAselection_H3K9me3$V4
loclist<-as.vector(loclist)

chrlist<-gRNAselection_H3K9me3$V3
chrlist<-as.vector(chrlist)



#####################################################################################################################
############################## Check distribution of gRNAs with a chromomap plot #######################
####################################################################################################################

# for visualisation of the location of the filtered gRNAs, chromomap is used.

# load the chromosome data
genome = BSgenome.Hsapiens.UCSC.hg19
si=seqinfo(genome)
si=si[paste0('chr',c(1:22,"X","Y"))]
sid<-as.data.frame(si)
write.csv(sid,file="chromosomes")

# create 2 files: 1 with chromosome lengths and 1 with the break locations which should be visualised
# Currently, a chromosome length file for Homo Sapiens is used. As an example, a .txt file is included with a number of locations.
# The input gRNA locations should be located in the "chromomap_annotation.txt" file. 
chromoMap("chromomap_chromosome.txt","chromomap_annotation.txt",
          chr_color=c("grey"),       
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("orange","yellow","red","blue","purple")),
          legend=T, lg_x=100,lg_y=250
)