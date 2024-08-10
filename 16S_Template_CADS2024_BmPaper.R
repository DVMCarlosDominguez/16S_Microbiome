#============================================================================================================#
#Microbiome Analysis

# DADA2 ----
setwd("path/to/raw_requences")
getwd()
path <- setwd("path/to/raw_requences")

# Forward and reverse fastq filenames have format: SAMPLENAME_F_filt.fastqsanger and SAMPLENAME_R_filt.fastqsanger
fnFs <- sort(list.files(path, pattern="_F.fastqsanger", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R.fastqsanger", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastqsanger
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
View(sample.names)

## INSPECT READ QUALITY PROFILES ----

plotQualityProfile(fnFs[1:20])
plotQualityProfile(fnRs[1:20])

##FILTER AND TRIM ----

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastqsanger"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastqsanger"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     maxN=0, maxEE=c(2,2), truncQ=25, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

save(out, file = "out.RData")
View(out)

##LEARN THE ERROR RATES ----
errF <- learnErrors(filtFs, multithread=TRUE)
save(errF, file = "errF.RData")
errR <- learnErrors(filtRs, multithread=TRUE)
save(errR, file = "errR.RData")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ = TRUE)

##DEREPLICATE -> UNIQUE SEQUENCES ----

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##SAMPLE INFERENCE ----
dadaFsPooled <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
save(dadaFsPooled, file = "dadaFsPooled.RData")
dadaFsPooled

dadaRsPooled <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)
save(dadaRsPooled, file = "dadaRsPooled.RData")
dadaRsPooled

##MERGE PAIRED READS ----
mergers <- mergePairs(dadaFsPooled, derepFs, dadaRsPooled, derepRs, verbose=TRUE)

##Frequency table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Explorar la distribucion de los largos de secuencias
table(nchar(getSequences(seqtab)))

##REMOVING CHIMERA ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) 

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_FASTA.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.csv(asv_tab, "ASVs_counts")

##TAXONOMY ----
taxa <- assignTaxonomy(seqtab.nochim, "/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa[is.na(taxa)] <- "unclassified"
write.csv(taxa, file = "ASVs_taxonomy")

# PHYLOSEQ ----
DADA2 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  tax_table(taxa))

# Change row names
taxa_names(DADA2) #dna strings

# Save all the dna strings as a "DNAStringSet"
dna <- Biostrings::DNAStringSet(taxa_names(DADA2))
dna

#Match dna stringsets with the taxa_table
names(dna) <- taxa_names(DADA2)

#Add dna strings to the phyliseq object
DADA2 <- merge_phyloseq(DADA2, dna)
DADA2

#Rename using ASV_
taxa_names(DADA2) <- paste0("ASV", seq(ntaxa(DADA2)))

DADA2_sub <- DADA2

## FILTER ----
#Now we will filter out Eukaryotes, Archaea, chloroplasts and mitochondria, because 
#we only intended to amplify bacterial sequences.
length(which(DADA2_sub@tax_table == "Bacteria"))
length(which(DADA2_sub@tax_table == "Archaea"))
length(which(DADA2_sub@tax_table == "Eukaryota"))
length(which(DADA2_sub@tax_table == "unclassified"))
length(which(DADA2_sub@tax_table == "Chloroplast"))
length(which(DADA2_sub@tax_table == "Mitochondria"))

erie <- DADA2_sub %>%
    subset_taxa(
        Kingdom == "Bacteria" &
            Order != "Chloroplast" &
            Family  != "Mitochondria"
    )

erie

#Checkind seqs after filtering
length(which(erie@tax_table == "Bacteria"))
length(which(erie@tax_table == "Archaea"))
length(which(erie@tax_table == "Eukaryota"))
length(which(erie@tax_table == "unclassified"))
length(which(erie@tax_table == "Chloroplast"))
length(which(erie@tax_table == "Mitochondria"))

#Eliminar "unclassified Phylum"
erie <- erie %>%
    subset_taxa(Phylum != "unclassified")

erie

##DECONTAM ----

library(decontam)

metadata_decom <- map
rownames(metadata_decom) <- metadata_decom$SampleID
psdec <- erie
sample_data(psdec) <- metadata_decom

contamdf.freq <- isContaminant(psdec, method="frequency", conc="quant_reading")
contamdf.freq
table(contamdf.freq$contaminant)

sample_data(psdec)$is.neg <- sample_data(psdec)$SampleType == "Control"
contamdf.prev <- isContaminant(psdec, method="prevalence", neg="is.neg")
contamdf.prev
table(contamdf.prev$contaminant)

pop_taxa = function(physeq, badTaxa){
    allTaxa = taxa_names(physeq)
    allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
    return(prune_taxa(allTaxa, physeq))
}

badTaxa = c("ASV89", "ASV280", "ASV300", "ASV443", "ASV488", "ASV656", 
            "ASV670", "ASV703", "ASV722", "ASV788", "ASV835", "ASV856", 
            "ASV881", "ASV936", "ASV948", "ASV1009", "ASV1028", "ASV1036", 
            "ASV1070") 

eriefilt = pop_taxa(erie, badTaxa)
eriefilt
erie

rank_names(erie)
phylum_table <- table(tax_table(erie)[, "Phylum"], exclude = NULL) # Create table, number of features for each phyla
phylum_table

#As a first analysis, we will look at the distribution of read counts
#Data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(erie))
sample_sum_df

##ANALYSIS ----
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
    geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
    ggtitle("Distribution of sample sequencing depth") + 
    xlab("Read counts") +
    theme(axis.title.y = element_blank())

# mean, max and min of sample read counts
min(sample_sums(erie))
mean(sample_sums(erie))
max(sample_sums(erie))

###Abundance table----
erie_phylum <- erie %>%
    tax_glom(taxrank = "Species") %>%                     # agglomerate at phylum level - modificar esto depende del nivel taxonÃ³mico
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt() %>%                                         # Melt to long format
    filter(Abundance > 0) %>%                         # Filter out low abundance taxa modificar segun la abundancia deseada
    arrange(Phylum)

#Relative abundance barplot, by phylum
p_abundancia<-ggplot(erie_phylum,aes(x=erie_phylum$Sample, y=erie_phylum$Abundance, fill=erie_phylum$Phylum)) +  geom_bar(stat = "identity")
p_abundancia + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + 
    labs(x = "", y = "") +       
    theme(axis.text.y = element_text(size = 8)) +
    scale_fill_discrete(name = "") 

### Diversity analysis ----
diversidad_alfa <- estimate_richness(erie, measures=c("Observed", "Simpson", "Shannon"))
diversidad_alfa
Alfa<- plot_richness(erie, measures =c("Observed", "Simpson", "Shannon"))
Alfa + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +#Cambiar el ángulo del texto en el eje x
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = "") +
    theme(axis.title.y = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8))

#

#Eliminar las muestras control
NoControl = subset_samples(erie, Type != "Control")
NoControl

###Core microbiome----
#Transfor to relative abundance
Abundancia.relativa <- transform(SinControl, "compositional")

# Taxa with over 50% prevance at X% relative abundance
Gen16S.core <- core(Abundancia.relativa, 
                    detection = 0.02/100, prevalence = 50/100, include.lowest = FALSE) #Todo lo que aparez en mínimo en 80% de los animales

Gen16S.core
Core_16S <- plot_heatmap(Gen16S.core, taxa.label = "Family")
Core_16S

Core_16S +  scale_x_discrete(labels= c("Mucus"= "Bm042")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +#Cambiar el ángulo del texto en el eje x
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = "") +
    theme(axis.title.y = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8))
