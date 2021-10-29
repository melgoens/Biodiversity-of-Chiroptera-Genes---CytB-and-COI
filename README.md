# Biodiversity-of-Chiroptera-Genes---CytB-and-COI
####Assignment 2----

#Main Question: Is the COI gene more effective than the CytB gene when it comes to identifying and clustering Chiropteran species (bats)?

###Loading the Packages----

library(rentrez)
library(tidyverse)
library(stringi)
library(ape)
library(RSQLite)
library(BiocManager)
library(Biostrings)
library(muscle)
library(DECIPHER)
library(kmer)
library(vegan)
library(ggdendro)
library(factoextra)
library(cluster)

####Entrez Functions - Gathering NCBI Data----
entrez_dbs() #From this, I want to find the nucleotide sequences so I will then investigate terms within the nuccore database
entrez_db_searchable(db = "nuccore")

ls_Chiroptera_results <- entrez_search(db = "nuccore", term = "Chiroptera[ORGN]")
ls_Chiroptera_results #2774301 hits
class(ls_Chiroptera_results)

#retmax is set to 20 by default so let's see how many hits we actually get
ls_Chiroptera_results$count #2774301 hits
#I will change the retmax value from the default 20 to 5000 to accommodate for more ids
ls_Chiroptera_results_max <- entrez_search(db = "nuccore", term = "Chiroptera", retmax = 5000)
length(ls_Chiroptera_results_max$ids) #checking
class(ls_Chiroptera_results_max) #checking

rm(ls_Chiroptera_results, ls_Chiroptera_results_max) #cleaning up

source("Entrez_Functions.R")
getwd() #making sure I'm in the right working directory (My Bioinformatics Assignment 2 folder where Entrez.Functions.R is saved)
Fetch_Chiroptera_CytB <- FetchFastaFiles(searchTerm = "Chiroptera[ORGN] AND 600:1000[SLEN] AND (cytochromeb[Gene] OR cytochrome b[Gene] OR cyt b[Gene] OR cytb[Gene])", seqsPerFile = 1000, fastaFileName = "CytB_fetch.fasta") #fetching the files with the ids and sequences

write(Fetch_Chiroptera_CytB, "CytB_fetch.fasta", sep = "\n") #writing the files to my folder
stringSetCytB <- readDNAStringSet("CytB_fetch.fasta0.fasta") #keeping one of the files with 1000 sequences in it for use in this study and writing the sequences so I can use them instead of id numbers
class(stringSetCytB) #checking
head(names(stringSetCytB)) #checking an example of what titles I have
head(stringSetCytB) #checking an example of what sequences I have

dfCytB <- data.frame(CytB_Title = names(stringSetCytB), CytB_Sequence = paste(stringSetCytB)) #converting the string set sequences and appropriate names to a data frame
View(dfCytB)

dfCytB$Species_Name <- word(dfCytB$CytB_Title, 2L, 3L) #adding a species name column
dfCytB <- dfCytB[, c("CytB_Title", "Species_Name", "CytB_Sequence")] #naming the titles appropriately
View(dfCytB)


#Again for COI
Fetch_Chiroptera_COI <- FetchFastaFiles(searchTerm = "Chiroptera[ORGN] AND 600:1000[SLEN] AND (COI[Gene])", seqsPerFile = 1000, fastaFileName = "COI_fetch.fasta") #Keeping both genes the same at 1000 sequences for easier comparison and shorter run times

write(Fetch_Chiroptera_COI, "COI_fetch.fasta", sep = "\n")
stringSetCOI <- readDNAStringSet("COI_fetch.fasta0.fasta")
class(stringSetCOI)
head(names(stringSetCOI))

dfCOI <- data.frame(COI_Title = names(stringSetCOI), COI_Sequence = paste(stringSetCOI))
View(dfCOI)

dfCOI$Species_Name <- word(dfCOI$COI_Title, 2L, 3L)
dfCOI <- dfCOI[, c("COI_Title", "Species_Name", "COI_Sequence")]
View(dfCOI)


###Subsetting and Merging the Two Genes for Comparison----

#Filtering the CytB data
length(unique(dfCytB$Species_Name)) #checking - 84 unique species

dfCytB_Subset <- dfCytB %>% #subsetting by one CytB sequence per species
  group_by(Species_Name) %>%
  sample_n(1)

all.equal(length(unique(dfCytB$Species_Name)), nrow(dfCytB_Subset)) #making sure there are the same number of rows and columns

#Filtering the COIB data
length(unique(dfCOI$Species_Name)) #checking - 153 unique species

dfCOI_Subset <- dfCOI %>% #subsetting by one COI sequence per species
  group_by(Species_Name) %>% 
  sample_n(1)

all.equal(length(unique(dfCOI$Species_Name)), nrow(dfCOI_Subset)) #making sure there are the same number of rows and columns

#Merge
dfAllSeqs <- merge(dfCytB_Subset, dfCOI_Subset, by = "Species_Name", all = TRUE)
View(dfAllSeqs) #merging so that each species has a COI gene and a CytB gene associated with it

dfOverlap <- merge(dfCytB_Subset, dfCOI_Subset, by = "Species_Name", all = FALSE)
View(dfOverlap) #84 common sequences

rm(dfAllSeqs)

###Filtering the CytB Sequences----
missing.data <- 0.02
length.var <- 40
chosen.model <- "TN93"
clustering.threshold <- 0.1
clustering.method <- "single"

dfCytB <- dfCytB %>% #creating a new column called Sequences2 and filtering those sequences to rid gaps and N's
  select(CytB_Sequence) %>%
  filter(!is.na(CytB_Sequence)) %>%
  mutate(CytB_Sequence2 = str_remove_all(CytB_Sequence, "^N+|N+$|-")) %>%
  filter(str_count(CytB_Sequence2, "N") <= (missing.data * str_count(CytB_Sequence2))) %>%
  filter(str_count(CytB_Sequence2) >= median(str_count(CytB_Sequence2)) - length.var & str_count(CytB_Sequence2) <= median(str_count(CytB_Sequence2)) + length.var)

###CytB Histogram of Sequence lengths----
class(dfCytB$CytB_Sequence2)
CytBSeq.Len <- nchar(dfCytB$CytB_Sequence2)#getting sequence lengths
length(CytBSeq.Len)
class(CytBSeq.Len)#checking for integer data
hist(CytBSeq.Len, main = "Histogram of CytB Sequence Lengths", xlab = "Sequence lengths (bp)", ylab = "Frequency", col = "blue") #making a histogram of the data

#CytB Distance Matrix
dfCytB <- as.data.frame(dfCytB) #setting it as a data frame
StringsetCytB.dend <- DNAStringSet(dfCytB$CytB_Sequence2)
names(StringsetCytB.dend) <- dfCytB$processid

dfCytB.alignment <- DNAStringSet(muscle::muscle(StringsetCytB.dend)) #aligning the sequences after viewing them and aligning them seperately in separate software
BrowseSeqs(dfCytB.alignment) #looking at the sequences in another window to make sure all looks well
dfCytB.alignment <- as.DNAbin(dfCytB.alignment) #converting format to a DNAbin

distanceMatrix.pairwise.CytB <- dist.dna(x = dfCytB.alignment, model = "raw", as.matrix = TRUE, pairwise.deletion = TRUE) #creating a distance matrix to be used in dendrogram

###CytB Dendrogram----
clusters.CytB.dend <- IdClusters(distanceMatrix.pairwise.CytB,
                           method = clustering.method,
                           cutoff = clustering.threshold,
                           showPlot = TRUE,
                           type = "dendrogram",
                           verbose = TRUE)


class(clusters.CytB.dend)
ggdendrogram(clusters.CytB.dend, rotate = FALSE, theme_dendro = FALSE) +
  labs(title = "CytB Dendrogram", x = "Species Names", y = "Divergence")

cut <- c(seq(0, 0.2, by = 0.005))
number.clusters <- vector("numeric", length(cut))
class(number.clusters)
number.clusters

for (i in 1:length(cut)) {
  clusters.temp <- IdClusters(distanceMatrix.pairwise.CytB, method = "complete", cutoff = cut[i], type = "clusters")
  number.clusters[i] <- length(unique(clusters.temp[, 1]))
} #creating the clusters in the dendrogram

DivergenceCytB <-plot(cut, number.clusters, xlab = "Clustering Threshold for Complete Linkage", ylab = "Number of Clusters", main = "CytB Sequence Divergence vs. Cluster Number in Chiroptera", col = "green") #dendrogram


###Filtering the COI Sequences----
dfCOI <- dfCOI %>%
  select(COI_Sequence) %>%
  filter(!is.na(COI_Sequence)) %>%
  mutate(COI_Sequence2 = str_remove_all(COI_Sequence, "^N+|N+$|-")) %>%
  filter(str_count(COI_Sequence2, "N") <= (missing.data * str_count(COI_Sequence2))) %>%
  filter(str_count(COI_Sequence2) >= median(str_count(COI_Sequence2)) - length.var & str_count(COI_Sequence2) <= median(str_count(COI_Sequence2)) + length.var)

###COI Histogram of Sequence lengths----
COISeq.Len <- nchar(dfCOI$COI_Sequence2)
length(COISeq.Len)
  class(COISeq.Len)
hist(COISeq.Len, main = "Histogram of COI Sequence Lengths", xlab = "Sequence lengths (bp)", ylab = "Frequency", col = "blue")

#Distance Matrix
dfCOI <- as.data.frame(dfCOI)
StringsetCOI.dend <- DNAStringSet(dfCOI$COI_Sequence2)
names(StringsetCOI.dend) <- dfCOI$processid

dfCOI.alignment <- DNAStringSet(muscle::muscle(StringsetCOI.dend))
BrowseSeqs(dfCOI.alignment)

dfCOI.alignment <- as.DNAbin(dfCOI.alignment)

distanceMatrix.pairwise.COI <- dist.dna(x = dfCOI.alignment, model = "raw", as.matrix = TRUE, pairwise.deletion = TRUE)
class(distanceMatrix.pairwise.COI) #checking

###COI Dendrogram----
clusters.COI.dend <- IdClusters(distanceMatrix.pairwise.COI,
                           method = clustering.method,
                           cutoff = clustering.threshold,
                           showPlot = TRUE,
                           type = "dendrogram",
                           verbose = TRUE)

class(clusters.COI.dend)
ggdendrogram(clusters.COI.dend, rotate = FALSE, theme_dendro = FALSE) +
  labs(title = "COI Dendrogram", x = "Species Names", y = "Divergence")
  

cut <- c(seq(0, 0.2, by = 0.005))
cut

number.clusters <- vector("numeric", length(cut))
class(number.clusters)
number.clusters

for (i in 1:length(cut)) {
  clusters.temp <- IdClusters(distanceMatrix.pairwise.COI, method = "complete", cutoff = cut[i], type = "clusters")
  number.clusters[i] <- length(unique(clusters.temp[, 1]))
}

DivergenceCOI <- plot(cut, number.clusters, xlab = "Clustering Threshold for Complete Linkage", ylab = "Number of Clusters", main = "COI Sequence Divergence vs. Cluster Number in Chiroptera", col = "green")

all.equal(DivergenceCOI, DivergenceCytB) #checking the difference - same divergence/cluster pattern


###Silhouette Indexes----
#CytB gene
clusters.CytB <- IdClusters(distanceMatrix.pairwise.CytB,
                            method = clustering.method,
                            cutoff = clustering.threshold,
                            showPlot = TRUE,
                            type = "both",
                            verbose = TRUE) #changing the "type" to both as opposed to "dendrogram" from earlier

class(clusters.CytB) #this is a list of a list and it needs to be a data frame

dfcluster.CytB <- clusters.CytB[[1]] #extracting from the list into a data frame
dfcluster.CytB

clustersCytB.vect <- dfcluster.CytB$cluster #naming
names(clustersCytB.vect) <- rownames(dfcluster.CytB)
silhouette.CytB <- silhouette(clustersCytB.vect, distanceMatrix.pairwise.CytB) #making the silhouette plot
plot(silhouette.CytB, do.n.k = TRUE, title = title(main = "example", xlab = "examplex"))

fviz_silhouette(silhouette.CytB, label = FALSE, print.summary = TRUE) + labs(title = "CytB Silhouette Plot") #adding colour

#COI gene
clusters.COI <- IdClusters(distanceMatrix.pairwise.COI,
                           method = clustering.method,
                           cutoff = clustering.threshold,
                           showPlot = TRUE,
                           type = "both",
                           verbose = TRUE)

class(clusters.COI)

dfcluster.COI <- clusters.COI[[1]]
dfcluster.COI
clustersCOI.vect <- dfcluster.COI$cluster
names(clustersCOI.vect) <- rownames(dfcluster.COI)
silhouette.COI <- silhouette(clustersCOI.vect, distanceMatrix.pairwise.COI)
plot(silhouette.COI)

fviz_silhouette(silhouette.COI, label = FALSE, print.summary = TRUE) + labs(title = "COI Silhouette Plot")


###Comparing Dendrograms----

cluster(dfCytB.alignment, k = 5, residues = NULL, gap = "-")
cluster(dfCOI.alignment, k = 5, residues = NULL, gap = "-")

#how many clusters do I have on the dendrograms?
max(dfcluster.CytB$cluster) #28 clusters
max(dfcluster.COI$cluster) #97 clusters
