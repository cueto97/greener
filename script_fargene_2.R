#extract fargene data
library(data.table)
library("vegan")
library(BiodiversityR)
library(microeco)
library("tibble")
library("dplyr")
library(phyloseq); library(tidyverse)
library("mirlyn")
library("BiocManager")
library(data.table)
library(file2meco)
library(rgl)
library(writexl)
library(ggpubr)
library(ggplot2)

isEmpty <- function(x) { #This function checks if a data frame is empty or not
  return(length(x)==0)
}

# Modificar el path y ejecutar con Rscript desde l�nea de comando
setwd("C:/Users/psanchez/Documents/R_fargene")
getwd()
path<-("C:/Users/psanchez/Documents/loop_fargene/")
getwd()
files_l<-list.files(path)#listamos cada una de las carpetas del fargene (n =84)

results<-data.frame(matrix(ncol = 3))

colnames(results)<-c("Name","Contigs","ARGs")
    
for(i in 1:length(files_l))
{
  
  results[i,1]<-files_l[i]
  
  var_cont<-0
  
  if(file.exists(paste0(path,files_l[i],"/predictedGenes/retrieved-contigs.fasta")) & file.size(paste0(path,files_l[i],"/predictedGenes/retrieved-contigs.fasta"))>0)
  {
    path_temp<-read.table(paste0(path,files_l[i],"/predictedGenes/retrieved-contigs.fasta"))
    
    for(o in 1:dim(path_temp)[1])
    {
      if(grepl(">",path_temp[o,]))
      {
        var_cont<-var_cont+1
      }
    }
  }
  
  results[i,2]<-var_cont
  
  var_cont<-0
  
  if(file.exists(paste0(path,files_l[i],"/predictedGenes/predicted-orfs.fasta")) & file.size(paste0(path,files_l[i],"/predictedGenes/retrieved-contigs.fasta"))>0)
  {
    path_temp<-read.delim(paste0(path,files_l[i],"/predictedGenes/predicted-orfs.fasta"),header = F)
    
    for(o in 1:dim(path_temp)[1])
    {
      if(grepl(">",path_temp[o,]))
      {
        var_cont<-var_cont+1
      }
    }
    
    results[i,3]<-var_cont
  }
}

#yp remove or replace data under certain value
#data<-replace(df$Marks, df$Marks<0, 0) 

library(stringr)
results$Sample <- sapply(strsplit((results$Name), '\\_'), "[", 4) # extract sampleID
results$model <- sapply(strsplit((results$Name), '\\_'), "[", 1) # extract sampleID
results$mode2 <- sapply(strsplit((results$Name), '\\_'), "[", 2) # extract sampleID

library(writexl)
getwd()
write_xlsx(results,"results_arg_full.xlsx")
#MODIFICAR EN EXCEL PARA A�ADIR LAS VARIABLES A LA METADATA
write.table(results, sep = "\t",quote = F,row.names = F)
results_2 <- read.table('results_arg_full.txt', header = T, row.names = 1, sep = '\t')#dataset at OTU level provided by Martin
results_2<- subset(results_2, model !="qnr")
#delete.na <- function(results_2, n=0) {
#  results_2[rowSums(is.na(results_2)) <= n,]
#}
#delete.na(results_2)
#NORMALIZAE THE DATA TO GIVE NUMBERS IN 20M READS
#results_2$ARGs[results_2$Sample=="GREEN48"]<- results_2$ARGs[results_2$Sample=="GREEN48"]*0.4


results_2$full <- interaction(results_2$Type, results_2$Location, results_2$Replicate)
results_2$full2 <- interaction(results_2$Type, results_2$Location)
results_c1 <- subset(results_2, Location =="C1")
results_c2<- subset(results_2, Location =="C2")
results_c3 <- subset(results_2, Location =="C3")
results_mm <- subset(results_2, Location =="MM")
ggplot(results_c1, aes(model, ARGs, fill = model)) +
  geom_col(width = 1) +
  facet_grid(full2~Replicate, scales = "free_x") +
  #facet_wrap(~full, scales = "free_x", drop = TRUE)+
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),axis.text.x = element_blank(), axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0,0.1)) + #scale_fill_taylor(palette="lover") +
  labs(x = NULL, y = 'N� features identified', fill = 'ARGs models')
dev.off()

ggplot(results_c2, aes(model, ARGs, fill = model)) +
  geom_col(width = 1) +
  facet_grid(full2~Replicate, scales = "free_x") +
  #facet_wrap(~full, scales = "free_x", drop = TRUE)+
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),axis.text.x = element_blank(), axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  #scale_y_continuous(expand = c(0,0.1)) + #scale_fill_taylor(palette="lover") +
  labs(x = NULL, y = 'N� features identified', fill = 'ARGs models')
dev.off()

ggplot(results_c3, aes(model, ARGs, fill = model)) +
  geom_col(width = 1) +
  facet_grid(full2~Replicate, scales = "free_x") +
  #facet_wrap(~full, scales = "free_x", drop = TRUE)+
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),axis.text.x = element_blank(), axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  #scale_y_continuous(expand = c(0,0.01)) + #scale_fill_taylor(palette="lover") +
  labs(x = NULL, y = 'N� features identified', fill = 'ARGs models')
dev.off()

ggplot(results_mm, aes(model, ARGs, fill = model)) +
  geom_col(width = 1) +
  facet_grid(full2~Replicate, scales = "free_x") +
  #facet_wrap(~full, scales = "free_x", drop = TRUE)+
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),axis.text.x = element_blank(), axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0,0.01)) + #scale_fill_taylor(palette="lover") +
  labs(x = NULL, y = 'N� features identified', fill = 'ARGs models')
dev.off()

#subir  metadata y cambiar eje x

#o<-ggplot(results,aes(Name,Contigs)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave(o,filename="genes_GREENER.pdf", device = "pdf")

#data from blast(resfinder)
getwd()
setwd("C:/Users/psanchez/Documents/R_fargene/resfinder")
#create a.list.of.objects (try to do)
GREEN48_RB <- read.table('results.blast_48.txt', header = FALSE, sep = '\t')#dataset at OTU level provided by Martin
GREEN48_RB$sample <- "GREEN48"
GREEN48_RB$Type <- "Biofilm"
GREEN48_RB$Location <- "C1"
#crear columna ge los genes
GREEN48_RB$gene <- sapply(strsplit(as.character(GREEN48_RB$V2), '\\_'), `[`, 1)
GREEN48_RB$gene2 <- sapply(strsplit(as.character(GREEN48_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN48_RB_filter <- GREEN48_RB[GREEN48_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN48_RB_filter$gene2)
GREEN48_RB_filter$count <- countsDF[match(GREEN48_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN48_RB_filter <- GREEN48_RB_filter[GREEN48_RB_filter$count> 80,]#I romove 5.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)



#GREEN51
GREEN51_RB <- read.table('results.blast_51.txt', header = FALSE, sep = '\t')#dataset at OTU level provided by Martin
GREEN51_RB$sample <- "GREEN51"
GREEN51_RB$Type <- "Biomass"
GREEN51_RB$Location <- "C1"
#crear columna ge los genes
GREEN51_RB$gene <- sapply(strsplit(as.character(GREEN51_RB$V2), '\\_'), `[`, 1)
GREEN51_RB$gene2 <- sapply(strsplit(as.character(GREEN51_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN51_RB_filter <- GREEN51_RB[GREEN51_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN51_RB_filter$gene2)
GREEN51_RB_filter$count <- countsDF[match(GREEN51_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN51_RB_filter <- GREEN51_RB_filter[GREEN51_RB_filter$count> 59,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)


#GREEN71
GREEN71_RB <- read.table('results.blast_71.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN71_RB$sample <- "GREEN71"
GREEN71_RB$Type <- "Biomass"
GREEN71_RB$Location <- "MM"
#crear columna ge los genes
GREEN71_RB$gene <- sapply(strsplit(as.character(GREEN71_RB$V2), '\\_'), `[`, 1)
GREEN71_RB$gene2 <- sapply(strsplit(as.character(GREEN71_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN71_RB_filter <- GREEN71_RB[GREEN71_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN71_RB_filter$gene2)
GREEN71_RB_filter$count <- countsDF[match(GREEN71_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN71_RB_filter <- GREEN71_RB_filter[GREEN71_RB_filter$count> 20,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)

#GREEN54
GREEN54_RB <- read.table('results.blast_54.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN54_RB$sample <- "GREEN54"
GREEN54_RB$Type <- "Biofilm"
GREEN54_RB$Location <- "C2"
#crear columna ge los genes
GREEN54_RB$gene <- sapply(strsplit(as.character(GREEN54_RB$V2), '\\_'), `[`, 1)
GREEN54_RB$gene2 <- sapply(strsplit(as.character(GREEN54_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN54_RB_filter <- GREEN54_RB[GREEN54_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN54_RB_filter$gene2)
GREEN54_RB_filter$count <- countsDF[match(GREEN54_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN54_RB_filter <- GREEN54_RB_filter[GREEN54_RB_filter$count> 69,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)

#GREEN57
GREEN57_RB <- read.table('results.blast_57.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN57_RB$sample <- "GREEN57"
GREEN57_RB$Type <- "Biomass"
GREEN57_RB$Location <- "C2"
#crear columna ge los genes
GREEN57_RB$gene <- sapply(strsplit(as.character(GREEN57_RB$V2), '\\_'), `[`, 1)
GREEN57_RB$gene2 <- sapply(strsplit(as.character(GREEN57_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN57_RB_filter <- GREEN57_RB[GREEN57_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN57_RB_filter$gene2)
GREEN57_RB_filter$count <- countsDF[match(GREEN57_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN57_RB_filter <- GREEN57_RB_filter[GREEN57_RB_filter$count> 73,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)

#GREEN61
GREEN61_RB <- read.table('results.blast_61.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN61_RB$sample <- "GREEN61"
GREEN61_RB$Type <- "Biofilm"
GREEN61_RB$Location <- "C3"
#crear columna ge los genes
GREEN61_RB$gene <- sapply(strsplit(as.character(GREEN61_RB$V2), '\\_'), `[`, 1)
GREEN61_RB$gene2 <- sapply(strsplit(as.character(GREEN61_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN61_RB_filter <- GREEN61_RB[GREEN61_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN61_RB_filter$gene2)
GREEN61_RB_filter$count <- countsDF[match(GREEN61_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN61_RB_filter <- GREEN61_RB_filter[GREEN61_RB_filter$count> 17,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)

#GREEN63
GREEN63_RB <- read.table('results.blast_63.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN63_RB$sample <- "GREEN63"
GREEN63_RB$Type <- "Biomass"
GREEN63_RB$Location <- "C3"
#crear columna ge los genes
GREEN63_RB$gene <- sapply(strsplit(as.character(GREEN63_RB$V2), '\\_'), `[`, 1)
GREEN63_RB$gene2 <- sapply(strsplit(as.character(GREEN63_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN63_RB_filter <- GREEN63_RB[GREEN63_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN63_RB_filter$gene2)
GREEN63_RB_filter$count <- countsDF[match(GREEN63_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN63_RB_filter <- GREEN63_RB_filter[GREEN63_RB_filter$count> 55,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)

#GREEN66
GREEN66_RB <- read.table('results.blast_66.txt', header = FALSE, sep = '\t', na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)#dataset at OTU level provided by Martin
GREEN66_RB$sample <- "GREEN66"
GREEN66_RB$Type <- "Biofilm"
GREEN66_RB$Location <- "MM"
#crear columna ge los genes
GREEN66_RB$gene <- sapply(strsplit(as.character(GREEN66_RB$V2), '\\_'), `[`, 1)
GREEN66_RB$gene2 <- sapply(strsplit(as.character(GREEN66_RB$gene), '\\('), `[`, 1)
#filtrar por evalue <1x10^-10
GREEN66_RB_filter <- GREEN66_RB[GREEN66_RB$V11< .000000000000001,]
#contar cuanto aparece cada grupo de genes
countsDF <- table(GREEN66_RB_filter$gene2)
GREEN66_RB_filter$count <- countsDF[match(GREEN66_RB_filter$gene2,  names(countsDF))]
#eliminar los que no est�n en el 10 (o 1)% de abundancia
GREEN66_RB_filter <- GREEN66_RB_filter[GREEN66_RB_filter$count> 14,]#I romove 4.000 features
#NORMALIZAR A 10M READS (BASADOS EN LOS READS VISTOS EN EL FASTQC)


#JUNTAR TODAS LAS TABLAS
results_blast <- rbind(GREEN71_RB_filter, GREEN51_RB_filter,GREEN48_RB_filter, GREEN54_RB_filter, GREEN57_RB_filter, GREEN61_RB_filter, GREEN63_RB_filter, GREEN66_RB_filter)
results_blast$counts.10M <- "1"
results_blast$counts.10M<- as.numeric(results_blast$counts.10M)
results_blast <- results_blast[complete.cases(results_blast),] #aqu� eliminamos los missing cases del dataframe

results_blast<- subset(results_blast, gene2 !="lsa")
results_blast<- subset(results_blast, gene2 !="lnu")
results_blast<- subset(results_blast, gene2 !="ole")
results_blast<- subset(results_blast, gene2 !="lnu")
results_blast<- subset(results_blast, gene2 !="msr")
results_blast<- subset(results_blast, gene2 !="dfrA15")
results_blast<- subset(results_blast, gene2 !="cmlB1")
results_blast<- subset(results_blast, gene2 !="cmlA1")
results_blast<- subset(results_blast, gene2 !="blaOXA-551")
#results_blast<- subset(results_blast, gene2 !="blaOXA-347")
results_blast<- subset(results_blast, gene2 !="blaPAO")
results_blast<- subset(results_blast, gene2 !="TOprJ4")
results_blast<- subset(results_blast, gene2 !="mph")
results_blast<- subset(results_blast, gene2 !="mcr-5.2")
results_blast<- subset(results_blast, gene2 !="qepA4")
#final plot
ggplot(data= results_blast, aes(x =gene2, y = counts.10M, fill=gene2))+
  facet_grid(Location~Type)+ #, labeller = labeller(temperature = temp.labs)) +
  #facet_wrap(type~location) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "ARGs genes", y = "N� Counts",  fill = 'ARGs genes')+
  scale_y_continuous(expand = c(0,0.1)) + #scale_fill_taylor(palette="lover") +
  geom_bar(stat ="identity")


GREEN51_RB <- read.table('results.blast_48.txt', header = FALSE, sep = '\t')#dataset at OTU level provided by Martin
GREEN51_RB$sample <- "GREEN48"



ggplot(data= GREEN48_RB, aes(sample, gene2))+ geom_bar(stat = "identity", position = "fill")
#IMPORT FULL BLAST TABLE


RESFINDER_blast<- read.table('results.blast.txt', header = FALSE, sep = '\t',  na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)
RESFINDER_blast$gene <- sapply(strsplit(as.character(RESFINDER_blast$V2), '\\_'), `[`, 1)
RESFINDER_blast$gene2 <- sapply(strsplit(as.character(RESFINDER_blast$gene), '\\('), `[`, 1)
RESFINDER_blast <- RESFINDER_blast[RESFINDER_blast$V11< .000000000000001,]#I romove 2.000 genes

countsDF <- table(RESFINDER_blast$gene2)
#to understand what are the most abundant genes
RESFINDER_blast$count <- countsDF[match(RESFINDER_blast$gene2,  names(countsDF))]
RESFINDER_blast <- RESFINDER_blast[order(-RESFINDER_blast$count),]
RESFINDER_blast <- RESFINDER_blast[RESFINDER_blast$count> 500,]#I romove 2.000 genes
ggplot(data= RESFINDER_blast, aes(V13, gene2))+ geom_bar(stat = "identity", position = "fill")
RESFINDER_blast$unique_counts <- "1"
RESFINDER_blast$unique_counts<- as.numeric(RESFINDER_blast$unique_counts)
library(tayloRswift)
ggplot(data= RESFINDER_blast, aes(x =gene2, y = unique_counts, fill=gene2))+
  facet_wrap(V13~.) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(expand = c(0,0.1)) +  scale_fill_taylor(palette="lover") +
  geom_bar(stat ="identity")

genus <- genus[order(-genus$Abundance),]
genus$Taxonomy <- factor(genus$Taxonomy, levels = unique(genus$Taxonomy))
paste(unique(genus$Taxonomy)[1:20],collapse = "','")#generamos una listilla con los phylum con comillas, para facilitar copia yh pega
genus$taxonomy.plot <- as.character(genus$Taxonomy)
top20 <- c('Bacteria|Proteobacteria|Gammaproteobacteria|Betaproteobacteriales|Burkholderiaceae|Burkholderia-Caballeronia-Paraburkholderia','Bacteria|Proteobacteria|Gammaproteobacteria|Vibrionales|Vibrionaceae|Enterovibrio','Bacteria|Bacteroidetes|Bacteroidia|Flavobacteriales|Flavobacteriaceae|Aureimarina','Bacteria|Proteobacteria|Gammaproteobacteria|Pseudomonadales|Moraxellaceae|Psychrobacter','Bacteria|Epsilonbacteraeota|Campylobacteria|Campylobacterales|Thiovulaceae|Sulfurimonas','Bacteria|Bacteroidetes|Bacteroidia|Flavobacteriales|Flavobacteriaceae|NS3a marine group','Bacteria|Proteobacteria|Gammaproteobacteria|Pseudomonadales|Moraxellaceae|Acinetobacter','Bacteria|Proteobacteria|Alphaproteobacteria|Rhodobacterales|Rhodobacteraceae|Thalassobius','Bacteria|Proteobacteria|Gammaproteobacteria|Alteromonadales|Pseudoalteromonadaceae|Pseudoalteromonas','Unassigned|||||','Bacteria|Proteobacteria|Gammaproteobacteria|Alteromonadales|Alteromonadaceae|Glaciecola','Bacteria|Proteobacteria|Alphaproteobacteria|Rhodobacterales|Rhodobacteraceae|Nautella','Bacteria|Bacteroidetes|Bacteroidia|Flavobacteriales|Cryomorphaceae|uncultured','Bacteria|Proteobacteria|Gammaproteobacteria|Oceanospirillales|Saccharospirillaceae|Oleibacter','Bacteria|Bacteroidetes|Bacteroidia|Flavobacteriales|Flavobacteriaceae|Polaribacter 4','Archaea|Euryarchaeota|Halobacteria|Halobacteriales|Halococcaceae|Halococcus','Bacteria|Bacteroidetes|Bacteroidia|Flavobacteriales|Flavobacteriaceae|Polaribacter','Bacteria|Bacteroidetes|Bacteroidia|Chitinophagales|Saprospiraceae|uncultured','Bacteria|Proteobacteria|Alphaproteobacteria|Rhodobacterales|Rhodobacteraceae|Tropicimonas','Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Pantoea')


genus.top <- genus[genus$Taxonomy%in%top20,]
genus.low <- genus[!(genus$Taxonomy%in%top20),]
genus.low$taxonomy.plot <- c('Other')

genus <- rbind(genus.top, genus.low)

temp.labs <- c("Temp. 1 (24�C)", "Temp. 2 (29�C)", "Temp. 3 (33�C)")
names(temp.labs) <- c("T1", "T2", "T3")


genus$taxonomy.plot <- factor(genus$taxonomy.plot, levels = c(sort(top20), 'Other'))
genus$location <- factor(genus$location, levels = c("Fish_Gut", "Fish_Gills", "Fish_Skin", "Water_Tank", "Water_Inlet"))
pdf("genus-total.pdf", height = 6, width =  12)
ggplot(genus, aes(month, Abundance, fill = taxonomy.plot)) +
  geom_col(width = 1) +
  facet_grid(temperature~location, labeller = labeller(temperature = temp.labs)) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA), 
        strip.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA)) +
  scale_y_continuous(expand = c(0,0.1)) +  scale_fill_taylor(palette="lover") +
  labs(x = 'Sampling month', y = 'Relative abundance (%)', fill = 'Genus')

ggplot(data=diamonds, aes(x=clarity)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)

ggplot(data= RESFINDER_blast, aes(x = gene2))+ geom_bar()

library(ggplot2)
#repasar ggplot package


GREEN51_RB <- read.table('results.blast_51.txt', header = T, row.names = 1, sep = '\t')#dataset at OTU level provided by Martin
GREEN71_RB <- read.table('results.blast_71.txt', header = T, row.names = 1, sep = '\t')#dataset at OTU level provided by Martin
#ad a header to each table

