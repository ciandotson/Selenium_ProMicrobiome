setwd("~/sel_test2")

#### FastQC on the reads ####
system('mkdir QC')
system('mkdir ./QC/raw_qc')
system(paste0("fastqc --noextract ./reads/raw/*fastq.gz -o ./QC/raw_qc"))

system("rm ./QC/raw_qc/*.zip")

#### Primer trimming using cutadapt ####
raw.fp <- "./reads/raw"
# change the path to your raw reads #

raw_for.fp <- sort(list.files(raw.fp, pattern = "_R1_001.fastq.gz", full.names = TRUE))
raw_rev.fp <- sort(list.files(raw.fp, pattern = "_R2_001.fastq.gz", full.names = TRUE))

sample.names <- strsplit(basename(raw_for.fp), "_L001_R1_001.fastq.gz")
sample.names <- as.character(sample.names)
# change the patterns to match those of your reads #


library(Biostrings); packageVersion('Biostrings')

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
for.primer <- "GTGCCAGCMGCCGCGGTAA"
rev.primer <- "GGACTACHVGGGTWTCTAAT"

for.ori <- allOrients(for.primer)
rev.ori <- allOrients(rev.primer)

library(ShortRead); packageVersion('ShortRead')

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}  

system('mkdir ./reads/pretrim')
pre_for.fp <- file.path('./reads/pretrim', paste0(sample.names, '_pretrim_R1.fastq.gz'))
pre_rev.fp <- file.path('./reads/pretrim', paste0(sample.names, '_pretrim_R2.fastq.gz'))

library(dada2); packageVersion('dada2')
prefilt.track <- filterAndTrim(raw_for.fp, pre_for.fp, raw_rev.fp, pre_rev.fp, minLen = 75)

for.revprim <- dada2::rc(for.primer)
rev.revprim <- dada2::rc(rev.primer)

system('mkdir ./reads/ptrimmed')
pt_for.fp <- file.path("./reads/ptrimmed", paste0(sample.names, "_ptrimmed_R1.fastq.gz"))
pt_rev.fp <- file.path("./reads/ptrimmed", paste0(sample.names, "_ptrimmed_R2.fastq.gz"))

for(i in seq_along(pre_for.fp)){
  system(paste0('cutadapt -g ', for.primer, ' -a ', rev.revprim, ' -G ', rev.primer, ' -A ', for.revprim, ' --pair-filter=any -m 50:50 -o ', pt_for.fp[i], ' -p ', pt_rev.fp[i], ' ', pre_for.fp[i], ' ', pre_rev.fp[i]))
}

# reads with primers trimmed will be found in the ./reads/ptrimmed/ file created within your working directory # 

#### dada2 Implementation ####
system('mkdir ./reads/postfilt')
post_for.fp <- file.path('./reads/postfilt', paste0(sample.names, '_filt_R1.fastq.gz'))
post_rev.fp <- file.path('./reads/postfilt', paste0(sample.names, '_filt_R2.fastq.gz'))

postfilt.track <- filterAndTrim(pt_for.fp, post_for.fp, pt_rev.fp, post_rev.fp, truncLen=c(230,230),
                                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                compress=TRUE, multithread=TRUE)

# postfilt.track is a nifty object that tells you number of reads that went into this filtering step and passed the filtering criteria # 

# Learning error rates # 
for.er <- learnErrors(post_for.fp, multithread=TRUE, verbose = TRUE)
rev.er <- learnErrors(post_rev.fp, multithread=TRUE, verbose = TRUE)

# Plotting Errors #
plotErrors(for.er, nominalQ=TRUE)
plotErrors(rev.er, nominalQ=TRUE)

# Constructing dada-class objects #
for.dada <- dada(post_for.fp, err=for.er, multithread=TRUE)
rev.dada <- dada(post_rev.fp, err=rev.er, multithread=TRUE)

# merging denoised forward and reversed reads #
remerged <- mergePairs(for.dada, post_for.fp, rev.dada, post_rev.fp, verbose=TRUE)

# Constructing Sequence Table #
all.st <- makeSequenceTable(remerged)
dim(all.st)
table(nchar(getSequences(all.st)))

# removing chimeras #
all_nochim.st <- removeBimeraDenovo(all.st, method="consensus", multithread=TRUE, verbose=TRUE)
dim(all_nochim.st)

# determine ratio of non-chimeras to all reads #
sum(all_nochim.st)/sum(all.st)
all_nochim.st <- t(all_nochim.st)

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
final_track <- cbind(prefilt.track[,1], prefilt.track[,2], postfilt.track[,2], sapply(for.dada, getN), sapply(rev.dada, getN), sapply(remerged, getN), rowSums(all_nochim.st))
colnames(final_track) <- c("pre-cutadapt", "post-cutadapt", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(final_track) <- sample.names
final_track <- as.data.frame(final_track)

# Taxonomy Assignment #
rdp.taxa <- assignTaxonomy(rownames(all_nochim.st), "./reference/rdp_19_toGenus_trainset.fa.gz")

rdp.taxa <- as.matrix(rdp.taxa)
silva.taxa

# Metadata table #
sel.met <- read.csv2('./Metadata.csv', sep = ',')
rownames(sel.met) <- sel.met$Sample
colnames(all_nochim.st) <- rownames(sel.met)

#### Phyloseq Object Construction and Filtering ####
raw_sel.ps <- phyloseq(otu_table(all_nochim.st, taxa_are_rows = TRUE),
                   sample_data(sel.met),
                   tax_table(rdp.taxa))

# Filter out contaminating reads and unused samples # 
raw_sel.ps <- subset_taxa(raw_sel.ps, Kingdom != 'Eukaryota')
raw_sel.ps <- subset_taxa(raw_sel.ps, Class != 'Chloroplast')
raw_sel.ps <- subset_samples(raw_sel.ps, sample_sums(raw_sel.ps) > 1000)
raw_sel.ps <- subset_taxa(raw_sel.ps, taxa_sums(raw_sel.ps) > 0)

# Adding DNA sequences to phyloseq object #
raw_sel.dna <- DNAStringSet(taxa_names(raw_sel.ps))
names(raw_sel.dna) <- taxa_names(raw_sel.ps)
raw_sel.ps <- merge_phyloseq(raw_sel.ps, raw_sel.dna)
taxa_names(raw_sel.ps) <- paste0("ASV", seq(ntaxa(raw_sel.ps)))

# Breaking down phyloseq object into data.frames and refseq objects #
raw_sel.tax <- as.data.frame(tax_table(raw_sel.ps))
raw_sel.met <- as(sample_data(sel.ps), 'data.frame')
raw_sel.otu <- as.data.frame(otu_table(raw_sel.ps))
raw_sel.fra <- cbind(raw_sel.tax, raw_sel.otu)
raw_sel.dna <- refseq(raw_sel.ps)

#### Cross-Validation of Reads Using blast ####
library(rBLAST)

# Create local blast database from the 16S rRNA database using rBLAST #
blast.tar <- blast_db_get("16S_ribosomal_RNA.tar.gz", baseURL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/', check_update = TRUE)
untar(blast.tar, exdir = './reference/16S_database')
list.files('./reference/16S_database')
blast.db <- blast(db = './reference/16S_database/16S_ribosomal_RNA')

# Performs the blast for each read and returns the best hit # 
sel.hits <- matrix(nrow = nrow(sel.tax), ncol = 12)
sel.hits <- as.data.frame(sel.hits) 
hold <- c()
for(i in 1:length(sel.dna)){
  hold <- predict(blast.db, raw_sel.dna[i])
  sel.hits[i,] <- hold[1,]
  raw_sel.tax$Best_Hit[i] <- hold[1, 2]
}

# Filter out reads that do not correspond to a NCBI entry #
library(dplyr); packageVersion('dplyr')
filt_sel.tax <- filter(raw_sel.tax, !is.na(raw_sel.tax$Best_Hit))

# Output the resulting NCBI entry names to a list #
if(!dir.exists("./blast_hits")){
  dir.create('./blast_hits')
}
write.table(filt_sel.tax$Best_Hit, './blast_hits/blast_hits.txt')

# calls a python script that assigns taxonomies based on the NCBI entry ID # 
system('python3 ./rRNA_blast.py')

# Read in the output from the python script and make new taxonomy table "ncbi_fin.tax" #
ncbi.taxa <- read.csv2('./blast_hits/ncbi_hits.csv', header = FALSE, fill = TRUE)
ncbi.int <- strsplit(as.character(ncbi.taxa$V1), ",")
ncbi_fin.tax <- do.call(rbind, lapply(ncbi.int, function(x) { length(x) <- max(sapply(ncbi.int, length)); x }))
ncbi_fin.tax <- as.data.frame(ncbi_fin.tax, stringsAsFactors = FALSE)
rownames(ncbi.taxa) <- rownames(filt_sel.tax)
colnames(ncbi_fin.tax) <- c('Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'hold')
for(i in 1:nrow(ncbi_fin.tax)){
  if(!is.na(ncbi_fin.tax$hold[i])){
    ncbi_fin.tax$Genus[i] <- ncbi_fin.tax$hold[i]
  }
}

ncbi_fin.tax <- ncbi_fin.tax[,1:7]

# Reorganize colnames of rdp taxonomy to denote they are the rdp taxonomy #
for(i in 1:ncol(filt_sel.tax)){
  colnames(filt_sel.tax)[i] <- paste0('rdp_', colnames(filt_sel.tax)[i])
}
filt_sel.tax <- cbind(filt_sel.tax, ncbi_fin.tax)

# Reorganize rownames of the ASV table and dna refseq object to match those of the new taxonomy table #
filt_sel.otu <- filter(raw_sel.otu, rownames(raw_sel.otu) %in% rownames(filt_sel.tax))
seldna.df <- as.data.frame(raw_sel.dna)
seldna.df <- filter(seldna.df, rownames(seldna.df) %in% rownames(filt_sel.tax))
filt_sel.dna <- DNAStringSet(seldna.df$x)
names(filt_sel.dna) <- rownames(filt_sel.tax)
filt_sel.tax <- as.matrix(filt_sel.tax)
filt_sel.met <- as(sample_data(raw_sel.ps), 'data.frame')

filt_sel.ps <- phyloseq(otu_table(filt_sel.otu, taxa_are_rows = TRUE),
                   sample_data(filt_sel.met),
                   tax_table(filt_sel.tax),
                   refseq(filt_sel.dna))

taxa_names(filt_sel.ps) <- refseq(filt_sel.ps)
taxa_names(sel.ps)

#### Tree Construction ####
# Produce a multiple sequence alignment of all reads that have a total sum >= 1,000 #
library(msa); packageVersion('msa')
sel.msa <- msa(filt_sel.dna, method = "Muscle", type = "dna", order = "input", verbose = TRUE)

# Construct a basic tree using "JC69" as a simple model, then use a gtr model to mkae an optimized tree # 
library(phangorn); packageVersion('phangorn')
sel.phang <- as.phyDat(sel.msa, type="DNA", names = taxa_names(sel.ps))
sel.dm <- dist.ml(sel.phang)
sel.tree <- NJ(sel.dm) # Note, tip order != sequence order
sel.tfit = pml(sel.tree, data=sel.phang)
sel.gtrfit <- update(sel.tfit, k=4, inv=0.2)
sel.gtrfit <- optim.pml(sel.gtrfit, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))

# Make phyloseq object that incldues the phylogenetic tree #
sel.otu <- as.data.frame(otu_table(sel.ps))
sel.tax <- as.matrix(tax_table(sel.ps))
sel.fme <- as(sample_data(sel.ps), 'data.frame')
sel.dna <- refseq(sel.ps)

sel.tax <- as.matrix(filt_sel.tax)
sel.ps <- phyloseq(otu_table(filt_sel.otu, taxa_are_rows = TRUE),
                   tax_table(filt_sel.tax),
                   sample_data(filt_sel.met),
                   refseq(filt_sel.dna),
                   phy_tree(sel.gtrfit$tree))


sel.met <- as(sample_data(sel.ps),'data.frame')

# Rename the ASVs such that it is in order of total abundance across all remaining samples and includes the most specific taxonomic classification in parentheses #
taxa_names(sel.ps) <- paste0('ASV', seq(ntaxa(sel.ps)))
sel.tax <- as.data.frame(tax_table(sel.ps))
for(i in 1:nrow(sel.tax)){
  if(!is.na(sel.tax$Genus[i])){
    taxa_names(sel.ps)[i] = paste0(taxa_names(sel.ps)[i], '(', sel.tax$Genus[i], ')')
  }else if(!is.na(sel.tax$Family[i])){
    taxa_names(sel.ps)[i] = paste0(taxa_names(sel.ps)[i], '(', sel.tax$Family[i], ')')
  }else if(!is.na(sel.tax$Order[i])){
    taxa_names(sel.ps)[i] = paste0(taxa_names(sel.ps)[i], '(', sel.tax$Order[i], ')')
  }else{
    taxa_names(sel.ps)[i] = paste0(taxa_names(sel.ps)[i], '(NA)')
  }
}

# Create a phyloseq object that just contains the single nodule sample #
sel_no.ps <- subset_samples(sel.ps, Type == 'Nodule')
sel_no.ps <- subset_taxa(sel_no.ps, taxa_sums(sel_no.ps) > 0)

# New ohyloseq object without nodule sample #
sel.ps <- subset_samples(sel.ps, Type != 'Nodule')
sel.ps <- subset_taxa(sel.ps, taxa_sums(sel.ps) > 0)

# Create a new variable that joins all values for the three test variables #
sel.met <- as(sample_data(sel.ps), 'data.frame')
for(i in 1:nrow(sel.met)){
  sel.met$Tres[i] <- paste0(substr(sel.met$Treatment[i], 1,1), substr(sel.met$Type[i],1,2), substr(sel.met$Species[i],4,4))
}
sample_data(sel.ps) <- sel.met

#### Beta Diversity Measurements and Visualizations ####

# Construct a Weighted Unifrac distance matrix and PCoA ordination #
sel.ps <- subset_taxa(sel.ps, taxa_sums(sel.ps) > 0)

sel_prop.ps <- transform_sample_counts(sel.ps, function(x) x/sum(x))

set.seed(248)
sel.wuni <- phyloseq::distance(sel_prop.ps, method = 'wunifrac')
sel.pcoa <- phyloseq::ordinate(sel_prop.ps, 'PCoA', distance = sel.wuni)


library(vegan); packageVersion('vegan')
# Perform an NMDS analysis using the weighted Unifrac distance matrix, with the PCoA ordination as the starting ordination # 
sel.nmds <- metaMDS(sel.wuni, 
                    k = 5, try = 20, trymax = 1000, maxit = 999,
                    model = 'global', 
                    autotransform = FALSE, previous.best = sel.pcoa$vectors[,1:5])

# Save the loading scores for all axes and make a distance matrix from these scores #
sel_nmds.scores <- scores(sel.nmds, display = 'sites')
sel_nmds.dist <- dist(sel_nmds.scores)

# Fit a linear model using the vectorized weighted Unifrac ditsance matrix as a predictor of the vectorized loading score distance matrix to fine total R^2 of the model # 
sel_nmds.ffit <- lm(as.vector(sel_nmds.dist) ~ as.vector(sel.wuni))
summary(sel_nmds.ffit)
sel_nmds.totr2 <- summary(sel_nmds.ffit)$r.squared

# Axes Variance Calculation #
# Fit linear models as before expect to preidtc the distance matrix of each individual axis #
sel_nmds.dist1 <- dist(sel_nmds.scores[,1])
sel_nmds.fit1 <- lm(as.vector(sel_nmds.dist1)~as.vector(sel.wuni))
sel_nmds.r1 <- summary(sel_nmds.fit1)$r.squared

sel_nmds.dist2 <- dist(sel_nmds.scores[,2])
sel_nmds.fit2 <- lm(as.vector(sel_nmds.dist2)~as.vector(sel.wuni))
sel_nmds.r2 <- summary(sel_nmds.fit2)$r.squared

sel_nmds.dist3 <- dist(sel_nmds.scores[,3])
sel_nmds.fit3 <- lm(as.vector(sel_nmds.dist3)~as.vector(sel.wuni))
sel_nmds.r3 <- summary(sel_nmds.fit3)$r.squared

sel_nmds.dist4 <- dist(sel_nmds.scores[,4])
sel_nmds.fit4 <- lm(as.vector(sel_nmds.dist4)~as.vector(sel.wuni))
sel_nmds.r4 <- summary(sel_nmds.fit4)$r.squared

sel_nmds.dist5 <- dist(sel_nmds.scores[,5])
sel_nmds.fit5 <- lm(as.vector(sel_nmds.dist5)~as.vector(sel.wuni))
sel_nmds.r5 <- summary(sel_nmds.fit5)$r.squared

# Take the sum the R^2 value from each axis #
sel_nmds.comb <- sel_nmds.r1 + sel_nmds.r2 + sel_nmds.r3 + sel_nmds.r4 + sel_nmds.r5

# Divide each axis R^2 by the total of all axes and then multiply by the variation explained by the whole model
sel_nmds.axisr <- c()
for(i in 1:ncol(sel_nmds.scores)){
  sel_nmds.axisr[i] <- (get(paste0('sel_nmds.r', i)) / sel_nmds.comb) * sel_nmds.totr2 
}

### Calculating Variance Components ###
library(lme4); packageVersion('lme4')

# Construct a data.frame that has sample info and their loading scores #
sel.met <- as(sample_data(sel.ps), 'data.frame')
sel.met$Plants <- factor(sel.met$Species, levels = c('A. crotalariae', 'A. lentiginosus'))
sel.met$Comps <- factor(sel.met$Type, c('Bulk Soil', 'Rhizosphere', 'Root Endosphere', 'Leaf Endosphere'))
sel.met$Treats <- factor(sel.met$Treatment, levels = c(0, 20, 80, 100))
sample_data(sel.ps) <- sel.met
sel_nmds.load <- cbind(sel.met, sel_nmds.scores)

# Test the mixed linear model on the first NMDS axis #
sel_nmds.vfit1 <- lmer(NMDS1 ~  (1|Plants) + (1|Comps) + (1|Treats) +(1|Plants:Comps) + (1|Plants:Treats) + (1|Comps:Treats) + (1|Comps:Plants:Treats), data = sel_nmds.load, REML = TRUE)
summary(sel_nmds.vfit1)
sel_nmds.vca1 <- as.data.frame(VarCorr(sel_nmds.vfit1))

# Using Loop to do each NMDS axis #
sel_nmds.vca <- matrix(nrow = 8, ncol = ncol(sel_nmds.scores))
hold <- c()
for(i in 1:ncol(sel_nmds.scores)){
  hold <- lmer(sel_nmds.scores[,i] ~  (1|Species) + (1|Type) + (1|Treatment) +(1|Species:Type) + (1|Species:Treatment) + (1|Type:Treatment) + (1|Type:Species:Treatment), data = sel_nmds.load, REML = TRUE)
  hold <- as.data.frame(VarCorr(hold))
  sel_nmds.vca[1,i] <- hold[1,4]
  sel_nmds.vca[2,i] <- hold[2,4]
  sel_nmds.vca[3,i] <- hold[3,4]
  sel_nmds.vca[4,i] <- hold[4,4]
  sel_nmds.vca[5,i] <- hold[5,4]
  sel_nmds.vca[6,i] <- hold[6,4]
  sel_nmds.vca[7,i] <- hold[7,4]
  sel_nmds.vca[8,i] <- hold[8,4]
}

# Save the variance components to their assigned variable/variable interaction and their NMDS loading axis #
rownames(sel_nmds.vca) <- c('Compartment x Plant x Treatment', 'Plant x Treatment', 'Plant x Compartment', 'Treatment x Compartment', 'Plant', 'Compartment', 'Treatment', 'Residuals')
colnames(sel_nmds.vca) <- colnames(sel_nmds.scores)

# Calculate the total variance of each variance component#
sel_nmds.vtot <- colSums(sel_nmds.vca)

# Weight each variance component by the amount of variation each axis explains #
sel_nmds.wvca <- matrix(nrow = nrow(sel_nmds.vca), ncol = length(sel_nmds.axisr))
for(i in 1:length(sel_nmds.axisr)){
  for(j in 1:nrow(sel_nmds.vca)){
    sel_nmds.wvca[j,i] <- sel_nmds.vca[j,i]*sel_nmds.axisr[i] 
  }
}
# Take the total variance explained by each predictor and take the sum of those values #
rownames(sel_nmds.wvca) <- rownames(sel_nmds.vca); colnames(sel_nmds.wvca) <- colnames(sel_nmds.vca)
sel_nmds.tvca <- rowSums(sel_nmds.wvca)
sel_nmds.tvca <- as.data.frame(sel_nmds.tvca)
sel_nmds.ptot <- colSums(sel_nmds.tvca)

# Take the variance explained by each predictor and divide by the total vraiance explained and multiply by 100% #
sel_nmds.pvca <- matrix(nrow = nrow(sel_nmds.tvca), ncol = 1)
for(i in 1:nrow(sel_nmds.tvca)){
  sel_nmds.pvca[i,1] <- sel_nmds.tvca[i,1] / sel_nmds.ptot * 100
}

# Save the variation explained percentages of each predictor/predictor interaction #
rownames(sel_nmds.pvca) <- rownames(sel_nmds.vca); colnames(sel_nmds.pvca) <- 'Variance Explained'
sel_nmds.pvca

# Perform PermANOVA using all samples #
sel.adon <- adonis2(sel.wuni~Plants*Comps*Treats, sel.met, permutations = 9999)
sel.adon_by <- adonis2(sel.wuni~Plants*Comps*Treats, sel.met, permutations = 9999, by = 'terms')
sel.adon
sel.adon_by

# Perform PermDISP using all samples # 
sel.bdis <- betadisper(sel.wuni, group = sel.met$Tres)
anova(sel.bdis)
TukeyHSD(sel.bdis)

# Make a new variable that combines plant species and compartment #
for(i in 1:nrow(sel.met)){
  sel.met$PC[i] <- paste0(sel.met$Species[i], "; ", sel.met$Type[i]) 
}
sel.met$PC_fact <- factor(sel.met$PC, levels = c("A. crotalariae; Bulk Soil", "A. crotalariae; Rhizosphere", "A. lentiginosus; Rhizosphere", "A. crotalariae; Root Endosphere", "A. lentiginosus; Root Endosphere"))

# Make an object with the data to be plotted #
sel_nmds.load <- cbind(sel.met, sel_nmds.scores)
library(ggplot2);packageVersion('ggplot2')
# plot the NMDS ordination of all samples that will be used to make a patchwork plot #
sel_nmds.plot <- ggplot(sel_nmds.load, aes(NMDS1, NMDS2, color = PC_fact, shape = Treats)) +
  geom_point(size = 8) +
  scale_color_manual(name = expression(italic('Astragalus')~'Species; Compartment'), labels = c('Both Species; Bulk Soil',  italic('A. crotalariae;')~' Rhizosphere', italic('A. lentiginosus;')~' Rhizosphere', italic('A. crotalariae;')~' Root Endosphere', italic('A. lentiginosus;')~' Root Endosphere') , values = c("#cc79a7","#D55e00","#F0e442", "#0072b2","#009e73")) +
  scale_shape_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"),expression("20" * mu * "M"), expression("80" * mu * "M"), expression("100" * mu * "M")), values = c(16,15,18,17)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.6724)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.1457)")) +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid = element_blank(),
        legend.text = element_text(size = 12, family = "Garamond"),
        legend.title = element_text(size = 16, family = "Garamond"),
        axis.title = element_text(size=20, family = "Garamond"),
        axis.text = element_text(color = "black", size = 8, family = "Garamond")) +
    annotate('text', x = 0.015, y = -0.15,
             label = expression("(PERMANOVA) F"["14,44"] ~ "= 9.7161, P < 0.001;  (PERMDISP) F"["14,44"] ~  "= 3.1477, P = 0.0095"),
             size = 4, family = 'Garamond') +
  coord_cartesian(ylim = c(-0.145, 0.08)) +
  labs(tag = 'A.')

sel_nmds.plot

### Per Species Breakdown ###
# Make a phyloseq object that only includes A. crotalariae samples #
crot.ps <- subset_samples(sel.ps, Species == 'A. crotalariae')
crot.ps <- subset_samples(crot.ps, Type != 'Bulk Soil')
crot.ps <- subset_taxa(crot.ps, taxa_sums(crot.ps) > 0)

# Perform same analysis with this subset of the data as seen previously in lines 281-440 for lines 448-867#  
crot_prop.ps <- transform_sample_counts(crot.ps, function(x) x/sum(x))
set.seed(248)
crot.wuni <- phyloseq::distance(crot_prop.ps, method = 'wunifrac')
crot.pcoa <- phyloseq::ordinate(crot_prop.ps, 'PCoA', distance = crot.wuni)


# 3 axes is the smallest number that explains >99% of the variation #
crot.nmds <- metaMDS(crot.wuni, 
                     k = 3, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = crot.pcoa$vectors[,1:3])

crot_nmds.scores <- scores(crot.nmds, display = 'sites')
crot_nmds.dist <- dist(crot_nmds.scores)

### Calculating the Variance of the Principal Components (or NMDS axes in our case) ###
## Total Model Variance Calculation ##
crot_nmds.ffit <- lm(as.vector(crot_nmds.dist) ~ as.vector(crot.wuni))
crot_nmds.totr2 <- summary(crot_nmds.ffit)$r.squared

# Axes Variance Calculation #
crot_nmds.dist1 <- dist(crot_nmds.scores[,1])
crot_nmds.fit1 <- lm(as.vector(crot_nmds.dist1)~as.vector(crot.wuni))
crot_nmds.r1 <- summary(crot_nmds.fit1)$r.squared

crot_nmds.dist2 <- dist(crot_nmds.scores[,2])
crot_nmds.fit2 <- lm(as.vector(crot_nmds.dist2)~as.vector(crot.wuni))
crot_nmds.r2 <- summary(crot_nmds.fit2)$r.squared

crot_nmds.dist3 <- dist(crot_nmds.scores[,3])
crot_nmds.fit3 <- lm(as.vector(crot_nmds.dist3)~as.vector(crot.wuni))
crot_nmds.r3 <- summary(crot_nmds.fit3)$r.squared

crot_nmds.comb <- crot_nmds.r1 + crot_nmds.r2 + crot_nmds.r3
crot_nmds.axisr <- c()
for(i in 1:ncol(crot_nmds.scores)){
  crot_nmds.axisr[i] <- (get(paste0('crot_nmds.r', i)) / crot_nmds.comb) * crot_nmds.totr2 
}

### Calculating Variance Components ###

## Constructing a data.frame that has sample info and their loading scores ##
crot.met <- as(sample_data(crot.ps), 'data.frame')
crot_nmds.load <- cbind(crot.met, crot_nmds.scores)

## Testing the mixed linear model on the first NMDS axis
crot_nmds.vfit1 <- lmer(NMDS1 ~ (1|Comps) + (1|Treats)  + (1|Comps:Treats), data = crot_nmds.load)
summary(crot_nmds.vfit1)
crot_nmds.vca1 <- as.data.frame(VarCorr(crot_nmds.vfit1))
## Using Loop to do each NMDS axis ##
crot_nmds.vca <- matrix(nrow = 4, ncol = ncol(crot_nmds.scores))
hold <- c()
for(i in 1:ncol(crot_nmds.scores)){
  hold <- lmer(crot_nmds.scores[,i] ~ (1|Treats) + (1|Comps) + (1|Comps:Treats), data = crot_nmds.load)
  hold <- as.data.frame(VarCorr(hold))
  crot_nmds.vca[1,i] <- hold[1,4]
  crot_nmds.vca[2,i] <- hold[2,4]
  crot_nmds.vca[3,i] <- hold[3,4]
  crot_nmds.vca[4,i] <- hold[4,4]
}

rownames(crot_nmds.vca) <- c('Treatment x Compartment', 'Treatment', 'Compartment', 'Residuals')
colnames(crot_nmds.vca) <- colnames(crot_nmds.scores)
crot_nmds.vtot <- colSums(crot_nmds.vca)

crot_nmds.wvca <- matrix(nrow = nrow(crot_nmds.vca), ncol = length(crot_nmds.axisr))
for(i in 1:length(crot_nmds.axisr)){
  for(j in 1:nrow(crot_nmds.vca)){
    crot_nmds.wvca[j,i] <- crot_nmds.vca[j,i]*crot_nmds.axisr[i] 
  }
}

rownames(crot_nmds.wvca) <- rownames(crot_nmds.vca); colnames(crot_nmds.wvca) <- colnames(crot_nmds.vca)
crot_nmds.tvca <- rowSums(crot_nmds.wvca)
crot_nmds.tvca <- as.data.frame(crot_nmds.tvca)
crot_nmds.ptot <- colSums(crot_nmds.tvca)

crot_nmds.pvca <- matrix(nrow = nrow(crot_nmds.tvca), ncol = 1)
for(i in 1:nrow(crot_nmds.tvca)){
  crot_nmds.pvca[i,1] <- crot_nmds.tvca[i,1] / crot_nmds.ptot * 100
}

rownames(crot_nmds.pvca) <- rownames(crot_nmds.vca); colnames(crot_nmds.pvca) <- 'Variance Explained'
crot_nmds.pvca

## permanova ##
crot.adon <- adonis2(crot.wuni~Treats*Comps, crot.met, permutations = 9999)
crot.adon

crot.bdis <- betadisper(crot.wuni, group = crot.met$Tres, type = 'median')
anova(crot.bdis)
TukeyHSD(crot.bdis)

crot_nmds.plot <- ggplot(crot_nmds.load, aes(NMDS1, NMDS3, color = Comps, shape = Treatment)) +
  geom_point(size = 8) +
  scale_color_manual(name = expression('Compartment'), labels = c(italic('A. crotalariae;')~' Rhizosphere', italic('A. crotalariae;')~' Root Endosphere'), values = c("#D55e00", "#0072b2")) +
  scale_shape_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"),expression("20" * mu * "M"), expression("100" * mu * "M")), values = c(16,15,17)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.9073)")) +
  ylab(expression("NMDS3 (R"^2 ~ "= 0.0510)")) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        legend.text = element_text(size = 12, family = "Garamond"),
        legend.title = element_text(size = 16, family = "Garamond"),
        axis.title = element_text(size=20, family = "Garamond"),
        axis.text = element_text(color = "black", size = 8, family = "Garamond")) +
  annotate('text', x = 0.135, y = -0.26,
           label = expression("(PERMANOVA) F"["5,19"] ~ "= 6.4897, P < 0.001;  (PERMDISP) F"["5,19"] ~  "= 3.1267, P = 0.0581"),
           size = 4, family = 'Garamond') +
  coord_cartesian(ylim = c(-0.25,0.125)) +
  labs(tag = 'B.')

crot_nmds.plot

## Rhizosphere Comparison ##
crot_rhiz.ps <- subset_samples(crot.ps, Type == 'Rhizosphere')
crot_rhiz.ps <- subset_taxa(crot_rhiz.ps, taxa_sums(crot_rhiz.ps) > 0)
crot_rhiz.met <- as(sample_data(crot_rhiz.ps), 'data.frame')

crot_rhiz_prop.ps <- transform_sample_counts(crot_rhiz.ps, function(x) x/sum(x))
crot_rhiz_wuni <- phyloseq::distance(crot_rhiz_prop.ps, 'wunifrac', "samples")

crot_rhiz.met <- as(sample_data(crot_rhiz_prop.ps), 'data.frame')
crot_rhiz.adon <- adonis2(crot_rhiz_wuni~Treatment, crot_rhiz.met, permutations = 9999)
crot_rhiz.adon

crot_rhiz.bdis <- betadisper(crot_rhiz_wuni, group = crot_rhiz.met$Treats)
anova(crot_rhiz.bdis)
## Root Comparison ##
crot_root.ps <- subset_samples(crot.ps, Type == 'Root Endosphere')
crot_root.ps <- subset_taxa(crot_root.ps, taxa_sums(crot_root.ps) > 0)
crot_root.met <- as(sample_data(crot_root.ps), 'data.frame')

crot_root_prop.ps <- transform_sample_counts(crot_root.ps, function(x) x/sum(x))
crot_root_wuni <- phyloseq::distance(crot_root_prop.ps, 'wunifrac', "samples")
crot_root.pcoa <- phyloseq::ordinate(crot_root_prop.ps, method = "PCoA", distance = crot_root_wuni)

crot_root.met <- as(sample_data(crot_root_prop.ps), 'data.frame')
crot_root.adon <- adonis2(crot_root_wuni~Treatment, crot_root.met, permutations = 999)
crot_root.adon

crot_root.bdis <- betadisper(crot_root_wuni, group = crot_root.met$Treats)
anova(crot_root.bdis)

plot_ordination(crot_root_prop.ps, crot_root.pcoa, color = 'Treats')
  
## A. lentiginosus ##
lent.ps <- subset_samples(sel.ps, Species == 'A. lentiginosus')
lent.ps <- subset_taxa(lent.ps, taxa_sums(lent.ps) > 0)

lent_prop.ps <- transform_sample_counts(lent.ps, function(x) x/sum(x))

set.seed(248)
lent.wuni <- phyloseq::distance(lent_prop.ps, method = 'wunifrac')
lent.pcoa <- phyloseq::ordinate(lent_prop.ps, 'PCoA', distance = lent.wuni)


# 4 axes is the smallest number that explains >99% of the variation #
lent.nmds <- metaMDS(lent.wuni, 
                     k = 4, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = lent.pcoa$vectors[,1:4])

lent_nmds.scores <- scores(lent.nmds, display = 'sites')
lent_nmds.dist <- dist(lent_nmds.scores)

### Calculating the Variance of the Principal Components (or NMDS axes in our case) ###
## Total Model Variance Calculation ##
lent_nmds.ffit <- lm(as.vector(lent_nmds.dist) ~ as.vector(lent.wuni))
lent_nmds.totr2 <- summary(lent_nmds.ffit)$r.squared

# Axes Variance Calculation #
lent_nmds.dist1 <- dist(lent_nmds.scores[,1])
lent_nmds.fit1 <- lm(as.vector(lent_nmds.dist1)~as.vector(lent.wuni))
lent_nmds.r1 <- summary(lent_nmds.fit1)$r.squared

lent_nmds.dist2 <- dist(lent_nmds.scores[,2])
lent_nmds.fit2 <- lm(as.vector(lent_nmds.dist2)~as.vector(lent.wuni))
lent_nmds.r2 <- summary(lent_nmds.fit2)$r.squared

lent_nmds.dist3 <- dist(lent_nmds.scores[,3])
lent_nmds.fit3 <- lm(as.vector(lent_nmds.dist3)~as.vector(lent.wuni))
lent_nmds.r3 <- summary(lent_nmds.fit3)$r.squared

lent_nmds.dist4 <- dist(lent_nmds.scores[,4])
lent_nmds.fit4 <- lm(as.vector(lent_nmds.dist4)~as.vector(lent.wuni))
lent_nmds.r4 <- summary(lent_nmds.fit4)$r.squared

lent_nmds.comb <- lent_nmds.r1 + lent_nmds.r2 + lent_nmds.r3 + lent_nmds.r4
lent_nmds.axisr <- c()
for(i in 1:ncol(lent_nmds.scores)){
  lent_nmds.axisr[i] <- (get(paste0('lent_nmds.r', i)) / lent_nmds.comb) * lent_nmds.totr2 
}

### Calculating Variance Components ###

## Constructing a data.frame that has sample info and their loading scores ##
lent.met <- as(sample_data(lent.ps), 'data.frame')
lent_nmds.load <- cbind(lent.met, lent_nmds.scores)

## Testing the mixed linear model on the first NMDS axis
lent_nmds.vfit1 <- lmer(NMDS1 ~ (1|Treats) + (1|Comps) + (1|Comps:Treats), data = lent_nmds.load)
summary(lent_nmds.vfit1)
lent_nmds.vca1 <- as.data.frame(VarCorr(lent_nmds.vfit1))
## Using Loop to do each NMDS axis ##
lent_nmds.vca <- matrix(nrow = 4, ncol = ncol(lent_nmds.scores))
hold <- c()
for(i in 1:ncol(lent_nmds.scores)){
  hold <- lmer(lent_nmds.scores[,i] ~ (1|Treats) + (1|Comps) + (1|Comps:Treats), data = lent_nmds.load)
  hold <- as.data.frame(VarCorr(hold))
  lent_nmds.vca[1,i] <- hold[1,4]
  lent_nmds.vca[2,i] <- hold[2,4]
  lent_nmds.vca[3,i] <- hold[3,4]
  lent_nmds.vca[4,i] <- hold[4,4]
}

rownames(lent_nmds.vca) <- c('Treatment x Compartment', 'Treatment', 'Compartment', 'Residuals')
colnames(lent_nmds.vca) <- colnames(lent_nmds.scores)
lent_nmds.vtot <- colSums(lent_nmds.vca)

lent_nmds.wvca <- matrix(nrow = nrow(lent_nmds.vca), ncol = length(lent_nmds.axisr))
for(i in 1:length(lent_nmds.axisr)){
  for(j in 1:nrow(lent_nmds.vca)){
    lent_nmds.wvca[j,i] <- lent_nmds.vca[j,i]*lent_nmds.axisr[i] 
  }
}

rownames(lent_nmds.wvca) <- rownames(lent_nmds.vca); colnames(lent_nmds.wvca) <- colnames(lent_nmds.vca)
lent_nmds.tvca <- rowSums(lent_nmds.wvca)
lent_nmds.tvca <- as.data.frame(lent_nmds.tvca)
lent_nmds.ptot <- colSums(lent_nmds.tvca)

lent_nmds.pvca <- matrix(nrow = nrow(lent_nmds.tvca), ncol = 1)
for(i in 1:nrow(lent_nmds.tvca)){
  lent_nmds.pvca[i,1] <- lent_nmds.tvca[i,1] / lent_nmds.ptot * 100
}

rownames(lent_nmds.pvca) <- rownames(lent_nmds.vca); colnames(lent_nmds.pvca) <- 'Variance Explained'
lent_nmds.pvca

## permanova ##
lent.adon <- adonis2(lent.wuni~Treats*Comps, lent.met, permutations = 9999)
lent.adon

lent.bdis <- betadisper(lent.wuni, group = lent.met$Tres, type = 'median')
anova(lent.bdis)
TukeyHSD(lent.bdis)

lent_nmds.plot <- ggplot(lent_nmds.load, aes(NMDS1, NMDS2, color = Comps, shape = Treats)) +
  geom_point(size = 10) +
  scale_color_manual(name = expression('Compartment'), labels = c(italic('A. lentiginosus;')~' Rhizosphere', italic('A. lentiginosus;')~' Root Endosphere') , values = c("#F0e442", "#009e73")) +
  scale_shape_manual(name = expression("NaSeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M")), values = c(16,15,18)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.8159)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.0979)")) +
  theme_bw() +
  theme(legend.position= 'none',
        panel.grid = element_blank(),
        legend.text = element_text(size = 12, family = "Garamond"),
        legend.title = element_text(size = 16, family = "Garamond"),
        axis.title = element_text(size=20, family = "Garamond"),
        axis.text = element_text(color = "black", size = 8, family = "Garamond")) +
  annotate('text', x = -0.045, y = -0.26,
           label = expression("(PERMANOVA) F"["5,15"] ~ "= 5.1392, P < 0.001;  (PERMDISP) F"["5,15"] ~  "= 6.2152, P = 0.0095"),
           size = 4, family = 'Garamond') +
  coord_cartesian(ylim = c(-0.25, 0.12))+
  labs(tag = "C.")
lent_nmds.plot

## Rhizosphere Comparison ##
lent_rhiz.ps <- subset_samples(lent.ps, Type == 'Rhizosphere')
lent_rhiz.ps <- subset_taxa(lent_rhiz.ps, taxa_sums(lent_rhiz.ps) > 0)

lent_rhiz_prop.ps <- transform_sample_counts(lent_rhiz.ps, function(x) x/sum(x))
lent_rhiz_wuni <- phyloseq::distance(lent_rhiz_prop.ps, 'wunifrac', "samples")

lent_rhiz.met <- as(sample_data(lent_rhiz_prop.ps), 'data.frame')
lent_rhiz.adon <- adonis2(lent_rhiz_wuni~Treatment, lent_rhiz.met, permutations = 9999)
lent_rhiz.adon

lent_rhiz.bdis <- betadisper(lent_rhiz_wuni, group = lent_rhiz.met$Treats)
anova(lent_rhiz.bdis)

## Root Comparison ##
lent_root.ps <- subset_samples(lent.ps, Type == 'Root Endosphere')
lent_root.ps <- subset_taxa(lent_root.ps, taxa_sums(lent_root.ps) > 0)

lent_root_prop.ps <- transform_sample_counts(lent_root.ps, function(x) x/sum(x))
lent_root.wuni <- phyloseq::distance(lent_root_prop.ps, 'wunifrac', "samples")
lent_root.pcoa <- phyloseq::ordinate(lent_root_prop.ps, method = "PCoA")

lent_root.met <- as(sample_data(lent_root_prop.ps), 'data.frame')
lent_root.adon <- adonis2(lent_root.wuni~Treatment, lent_root.met, permutations = 9999)
lent_root.adon

lent_root.bdisp <- betadisper(lent_root.wuni, group = lent_root.met$Treats)
anova(lent_root.bdisp)
TukeyHSD(lent_root.bdisp)

lent_root.nmds <- metaMDS(lent_root.wuni, 
                     k = 2, try = 20, trymax = 1000, maxit = 999,
                     model = 'global', 
                     autotransform = FALSE, previous.best = lent.pcoa$vectors[,1:2])

lent_root_nmds.scores <- scores(lent_root.nmds, display = 'sites')
lent_root_nmds.dist <- dist(lent_root_nmds.scores)

lent_root_nmds.ffit <- lm(as.vector(lent_root_nmds.dist) ~ as.vector(lent_root.wuni))
lent_root_nmds.totr2 <- summary(lent_root_nmds.ffit)$r.squared

lent_root_nmds.dist1 <- dist(lent_root_nmds.scores[,1])
lent_root_nmds.fit1 <- lm(as.vector(lent_root_nmds.dist1)~as.vector(lent_root.wuni))
lent_root_nmds.r1 <- summary(lent_root_nmds.fit1)$r.squared

lent_root_nmds.dist2 <- dist(lent_root_nmds.scores[,2])
lent_root_nmds.fit2 <- lm(as.vector(lent_root_nmds.dist2)~as.vector(lent_root.wuni))
lent_root_nmds.r2 <- summary(lent_root_nmds.fit2)$r.squared

lent_root_nmds.comb <- lent_root_nmds.r1 + lent_root_nmds.r2
lent_root_nmds.axisr <- c()
for(i in 1:ncol(lent_root_nmds.scores)){
  lent_root_nmds.axisr[i] <- (get(paste0('lent_root_nmds.r', i)) / lent_root_nmds.comb) * lent_root_nmds.totr2 
}

lent_root_nmds.load <- cbind(lent_root.met, lent_root_nmds.scores)

lent_root.adon <- adonis2(lent_root.wuni~Treats, data = lent_root.met, permutations = 9999)
lent_root.adon

lent_root_nmds.plot <- ggplot(lent_root_nmds.load, aes(NMDS1, NMDS2, color = Comps, shape = Treats)) +
  geom_point(size = 10) +
  scale_color_manual(name = expression('Compartment'), labels = c(italic('A. lentiginosus;')~' Root Endosphere') , values = c("#009e73")) +
  scale_shape_manual(name = expression("NaSeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M")), values = c(16,15,18)) +
  xlab(expression("NMDS1 (R"^2 ~ "= 0.8788)")) +
  ylab(expression("NMDS2 (R"^2 ~ "= 0.0693)")) +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        legend.text = element_text(size = 12, family = "Garamond"),
        legend.title = element_text(size = 16, family = "Garamond"),
        axis.title = element_text(size=20, family = "Garamond"),
        axis.text = element_text(color = "black", size = 8, family = "Garamond"),) +
  annotate('text', x = -0.045, y = -0.21,
           label = expression("(PERMANOVA) F"["2,8"] ~ "= 4.206, P = 0.0341;  (PERMDISP) F"["2,8"] ~  "= 8.2494, P = 0.0313"),
           size = 4, family = 'Garamond') +
  coord_cartesian(ylim = c(-0.2,0.37)) +
  labs(tag = 'D.')
lent_root_nmds.plot

# F(2,8) = 4.2096; P = 0.01 

## Root No to high Comparison ##
lent_root_nl.ps <- subset_samples(lent_root.ps, Treatment != '20')
lent_root_nl.ps <- subset_taxa(lent_root_nl.ps, taxa_sums(lent_root_nl.ps) > 0)

lent_root_nl_prop.ps <- transform_sample_counts(lent_root_nl.ps, function(x) x/sum(x))
lent_root_nl.wuni <- phyloseq::distance(lent_root_nl_prop.ps, 'wunifrac', "samples")
lent_root_nl.pcoa <- phyloseq::ordinate(lent_root_nl_prop.ps, method = "NMDS")

lent_root_nl.met <- as(sample_data(lent_root_nl_prop.ps), 'data.frame')
lent_root_nl.adon <- adonis2(lent_root_nl.wuni~Treatment, lent_root_nl.met, permutations = 9999)
lent_root_nl.adon

lent_root_nl.bdisp <- betadisper(lent_root_nl.wuni, group = lent_root_nl.met$Treats)
anova(lent_root_nl.bdisp)
TukeyHSD(lent_root_nl.bdisp)

plot_ordination(lent_root_nl_prop.ps, lent_root_nl.pcoa, color = "Treats")

## Root No to low Comparison ##
lent_root_nh.ps <- subset_samples(lent_root.ps, Treatment != '80')
lent_root_nh.ps <- subset_taxa(lent_root_nh.ps, taxa_sums(lent_root_nh.ps) > 0)

lent_root_nh_prop.ps <- transform_sample_counts(lent_root_nh.ps, function(x) x/sum(x))
lent_root_nh.wuni <- phyloseq::distance(lent_root_nh_prop.ps, 'wunifrac', "samples")
lent_root_nh.pcoa <- phyloseq::ordinate(lent_root_nh_prop.ps, method = "NMDS")

lent_root_nh.met <- as(sample_data(lent_root_nh_prop.ps), 'data.frame')
lent_root_nh.adon <- adonis2(lent_root_nh.wuni~Treatment, lent_root_nh.met, permutations = 9999)
lent_root_nh.adon

lent_root_nh.bdisp <- betadisper(lent_root_nh.wuni, group = lent_root_nh.met$Treats)
anova(lent_root_nh.bdisp)
TukeyHSD(lent_root_nh.bdisp)

plot_ordination(lent_root_nh_prop.ps, lent_root_nh.pcoa, color = "Treats")

## Root No to low Comparison ##
lent_root_nn.ps <- subset_samples(lent_root.ps, Treatment != '0')
lent_root_nn.ps <- subset_taxa(lent_root_nn.ps, taxa_sums(lent_root_nn.ps) > 0)

lent_root_nn_prop.ps <- transform_sample_counts(lent_root_nn.ps, function(x) x/sum(x))
lent_root_nn.wuni <- phyloseq::distance(lent_root_nn_prop.ps, 'wunifrac', "samples")
lent_root_nn.pcoa <- phyloseq::ordinate(lent_root_nn_prop.ps, method = "NMDS")

lent_root_nn.met <- as(sample_data(lent_root_nn_prop.ps), 'data.frame')
lent_root_nn.adon <- adonis2(lent_root_nn.wuni~Treatment, lent_root_nn.met, permutations = 9999)
lent_root_nn.adon

lent_root_nn.bdisp <- betadisper(lent_root_nn.wuni, group = lent_root_nn.met$Treats)
anova(lent_root_nn.bdisp)
TukeyHSD(lent_root_nn.bdisp)

plot_ordination(lent_root_nn_prop.ps, lent_root_nn.pcoa, color = "Treats")


soil.ps <- subset_samples(sel.ps, Type == 'Bulk Soil')
soil.ps <- subset_taxa(soil.ps, taxa_sums(sel.ps) > 0)

soil_prop.ps <- transform_sample_counts(soil.ps, function(x) x/sum(x))
soil_prop.ps <- subset_samples(soil_prop.ps, Sample != 'CD_0C5S')
soil.wuni <- phyloseq::distance(soil_prop.ps, 'wunifrac', "samples")
soil.nmds <- phyloseq::ordinate(soil_prop.ps, method = "NMDS")

soil.met <- as(sample_data(soil_prop.ps), 'data.frame')
soil.adon <- adonis2(soil.wuni~Treatment, soil.met, permutations = 9999)
soil.adon

soil.bdisp <- betadisper(soil.wuni, group = soil.met$Treats)
anova(soil.bdisp)
TukeyHSD(soil.bdisp)

## Multiple test corrections ##
# Save the p-values from all PermANOVA tests # 
adon.pval <- cbind(sel.adon$`Pr(>F)`[1],
                   crot.adon$`Pr(>F)`[1],
                   crot_rhiz.adon$`Pr(>F)`[1],
                   crot_root.adon$`Pr(>F)`[1],
                   lent.adon$`Pr(>F)`[1],
                   lent_rhiz.adon$`Pr(>F)`[1],
                   lent_root.adon$`Pr(>F)`[1],
                   lent_root_nl.adon$`Pr(>F)`[1],
                   lent_root_nh.adon$`Pr(>F)`[1],
                   lent_root_nn.adon$`Pr(>F)`[1],
                   soil.adon$`Pr(>F)`[1])

# Perform test corrections #
colnames(adon.pval) <- seq(ncol(adon.pval))
adon.fdr <- p.adjust(adon.pval, method = "fdr")
adon.adj <- rbind(adon.pval, adon.fdr)

# Save the p-values from all PermDISP tests # 
bdis.pval <- cbind(anova(sel.bdis)$`Pr(>F)`[1],
                   anova(crot.bdis)$`Pr(>F)`[1],
                   anova(crot_rhiz.bdis)$`Pr(>F)`[1],
                   anova(crot_root.bdis)$`Pr(>F)`[1],
                   anova(lent.bdis)$`Pr(>F)`[1],
                   anova(lent_rhiz.bdis)$`Pr(>F)`[1],
                   anova(lent_root.bdisp)$`Pr(>F)`[1],
                   anova(lent_root_nl.bdisp)$`Pr(>F)`[1],
                   anova(lent_root_nh.bdisp)$`Pr(>F)`[1],
                   anova(lent_root_nn.bdisp)$`Pr(>F)`[1],
                   anova(soil.bdisp)$`Pr(>F)`[1])

# Perform test corrections #
colnames(bdis.pval) <- seq(ncol(bdis.pval))
bdis.fdr <- p.adjust(bdis.pval, method = 'fdr') 
bdis.adj <- rbind(bdis.pval, bdis.fdr)

# Create a plot that contains the NMDS ordinations of plots in which the PermANOVA test was significant # 
library(patchwork)
(sel_nmds.plot | crot_nmds.plot) /
(lent_nmds.plot | lent_root_nmds.plot) +
  plot_layout(guides = 'collect') &
  theme(plot.tag = element_text(size = 20, family = 'Liberation Sans', face = 'bold', color = 'black'))

#### Alpha Diversity Measurements and Visualizations ####
# Produce a data.frame that has caluclated alpha diversity metrics for each sample #
sel.rich <- estimate_richness(sel.ps)
sel.rich <- as.data.frame(sel.rich)
sel.rich <- cbind(sel.met, sel.rich)

# Calculate Shannon Evenness by dividing the Shannon Diversity value by the natual log of the number of observed species #
for(i in 1:nrow(sel.rich)){
  sel.rich$ShaEvn[i] <- sel.rich$Shannon[i]/log(sel.rich$Chao1[i], base = 2.718) 
}

# Perform ANOVAs for each alpha diversity metric across all samples #
sel_cha.aov <- aov(Chao1~Plants*Treats*Comps, sel.rich, )
summary(sel_cha.aov)
sel_evn.aov <- aov(ShaEvn~Plants*Treats*Comps, sel.rich)
summary(sel_evn.aov)
sel_sha.aov <- aov(Shannon~Plants*Treats*Comps, sel.rich)
summary(sel_sha.aov)

# Subset Data to include only A. lentiginosus samples and test for significance #
lent.rich <- filter(sel.rich, Plants == 'A. lentiginosus')
lent_cha.aov <- aov(Chao1~Treats*Comps, lent.rich)
summary(lent_cha.aov)
lent_evn.aov <- aov(ShaEvn~Treats*Comps, lent.rich)
summary(lent_evn.aov)
lent_sha.aov <- aov(Shannon~Treats*Comps, lent.rich)
summary(lent_sha.aov)

# Subset Data to include only A. crotalariae samples and test for significance #
crot.rich <- filter(sel.rich, Plants == 'A. crotalariae')
crot.rich <- filter(crot.rich, Comps != 'Bulk Soil')
crot_cha.aov <- aov(Chao1~Treats*Comps, crot.rich)
summary(crot_cha.aov)
crot_evn.aov <- aov(ShaEvn~Treats*Comps, crot.rich)
summary(crot_evn.aov)
crot_sha.aov <- aov(Shannon~Treats*Comps, crot.rich)
summary(crot_sha.aov)

# Subset by each compartment and plant species group (5 groups) and test significance of treatment within these groups #
## Shannon ##
crot_rhiz.rich <- filter(crot.rich, Comps == "Rhizosphere")
crot_rhiz_sha.aov <- aov(Shannon ~ Treats, data = crot_rhiz.rich)
summary(crot_rhiz_sha.aov)

crot_rhiz_sha.hsd <- TukeyHSD(crot_rhiz_sha.aov)
crot_rhiz_sha.hsd

library(multcompView)

crot_rhiz_sha.let <- multcompLetters4(crot_rhiz_sha.aov, crot_rhiz_sha.hsd)
crot_rhiz_sha.let

crot_root.rich <- filter(crot.rich, Comps == "Root Endosphere")
crot_root_sha.aov <- aov(Shannon ~ Treats, data = crot_root.rich)
summary(crot_root_sha.aov)

crot_root_sha.hsd <- TukeyHSD(crot_root_sha.aov)
crot_root_sha.hsd

crot_root_sha.let <- multcompLetters4(crot_root_sha.aov, crot_root_sha.hsd)
crot_root_sha.let

bulk.rich <- filter(sel.rich, Comps == 'Bulk Soil')
bulk_sha.aov <- aov(Shannon~Treats, bulk.rich)
summary(bulk_sha.aov)

bulk_sha.hsd <- TukeyHSD(bulk_sha.aov)
bulk_sha.hsd

bulk_sha.let <- multcompLetters4(bulk_sha.aov, bulk_sha.hsd)
bulk_sha.let

lent_rhiz.rich <- filter(lent.rich, Comps == "Rhizosphere")
lent_rhiz_sha.aov <- aov(Shannon ~ Treats, data = lent_rhiz.rich)
summary(lent_rhiz_sha.aov)

lent_rhiz_sha.hsd <- TukeyHSD(lent_rhiz_sha.aov)
lent_rhiz_sha.hsd

lent_rhiz_sha.let <- multcompLetters4(lent_rhiz_sha.aov, lent_rhiz_sha.hsd)
lent_rhiz_sha.let

lent_root.rich <- filter(lent.rich, Comps == "Root Endosphere")
lent_root_sha.aov <- aov(Shannon ~ Treats, data = lent_root.rich)
summary(lent_root_sha.aov)

lent_root_sha.hsd <- TukeyHSD(lent_root_sha.aov)
lent_root_sha.hsd

lent_root_sha.let <- multcompLetters4(lent_root_sha.aov, lent_root_sha.hsd)
lent_root_sha.let

## Chao1 ##
crot_rhiz_cha.aov <- aov(Chao1 ~ Treats, data = crot_rhiz.rich)
summary(crot_rhiz_cha.aov)

crot_rhiz_cha.hsd <- TukeyHSD(crot_rhiz_cha.aov)
crot_rhiz_cha.hsd

crot_rhiz_cha.let <- multcompLetters4(crot_rhiz_cha.aov, crot_rhiz_cha.hsd)
crot_rhiz_cha.let

crot_root_cha.aov <- aov(Chao1 ~ Treats, data = crot_root.rich)
summary(crot_root_cha.aov)

crot_root_cha.hsd <- TukeyHSD(crot_root_cha.aov)
crot_root_cha.hsd

crot_root_cha.let <- multcompLetters4(crot_root_cha.aov, crot_root_cha.hsd)
crot_root_cha.let

bulk_cha.aov <- aov(Chao1~Treats, bulk.rich)
summary(bulk_cha.aov)

bulk_cha.hsd <- TukeyHSD(bulk_cha.aov)
bulk_cha.hsd

bulk_cha.let <- multcompLetters4(bulk_cha.aov, bulk_cha.hsd)
bulk_cha.let

lent_rhiz_cha.aov <- aov(Chao1 ~ Treats, data = lent_rhiz.rich)
summary(lent_rhiz_cha.aov)

lent_rhiz_cha.hsd <- TukeyHSD(lent_rhiz_cha.aov)
lent_rhiz_cha.hsd

lent_rhiz_cha.let <- multcompLetters4(lent_rhiz_cha.aov, lent_rhiz_cha.hsd)
lent_rhiz_cha.let

lent_root_cha.aov <- aov(Chao1 ~ Treats, data = lent_root.rich)
summary(lent_root_cha.aov)

lent_root_cha.hsd <- TukeyHSD(lent_root_cha.aov)
lent_root_cha.hsd

lent_root_cha.let <- multcompLetters4(lent_root_cha.aov, lent_root_cha.hsd)
lent_root_cha.let

## Evenness ##
crot_rhiz_evn.aov <- aov(ShaEvn ~ Treats, data = crot_rhiz.rich)
summary(crot_rhiz_evn.aov)

crot_rhiz_evn.hsd <- TukeyHSD(crot_rhiz_evn.aov)
crot_rhiz_evn.hsd

crot_rhiz_evn.let <- multcompLetters4(crot_rhiz_evn.aov, crot_rhiz_evn.hsd)
crot_rhiz_evn.let

crot_root_evn.aov <- aov(ShaEvn ~ Treats, data = crot_root.rich)
summary(crot_root_evn.aov)

crot_root_evn.hsd <- TukeyHSD(crot_root_evn.aov)
crot_root_evn.hsd

crot_root_evn.let <- multcompLetters4(crot_root_evn.aov, crot_root_evn.hsd)
crot_root_evn.let

bulk_evn.aov <- aov(ShaEvn~Treats, bulk.rich)
summary(bulk_evn.aov)

bulk_evn.hsd <- TukeyHSD(bulk_evn.aov)
bulk_evn.hsd

bulk_evn.let <- multcompLetters4(bulk_evn.aov, bulk_evn.hsd)
bulk_evn.let

lent_rhiz_evn.aov <- aov(ShaEvn ~ Treats, data = lent_rhiz.rich)
summary(lent_rhiz_evn.aov)

lent_rhiz_evn.hsd <- TukeyHSD(lent_rhiz_evn.aov)
lent_rhiz_evn.hsd

lent_rhiz_evn.let <- multcompLetters4(lent_rhiz_evn.aov, lent_rhiz_evn.hsd)
lent_rhiz_evn.let

lent_root_evn.aov <- aov(ShaEvn ~ Treats, data = lent_root.rich)
summary(lent_root_evn.aov)

lent_root_evn.hsd <- TukeyHSD(lent_root_evn.aov)
lent_root_evn.hsd

lent_root_evn.let <- multcompLetters4(lent_root_evn.aov, lent_root_evn.hsd)
lent_root_evn.let

# Create lists that contain the letters as denoted from the TukeyHSD for each comparison #
sha.let <- c(bulk_sha.let$Treats$Letters,
             crot_rhiz_sha.let$Treats$Letters,
             rev(lent_rhiz_sha.let$Treats$Letters),
             crot_root_sha.let$Treats$Letters,
             rev(lent_root_sha.let$Treats$Letters))
evn.let <- c(bulk_evn.let$Treats$Letters,
            crot_rhiz_evn.let$Treats$Letters,
            rev(lent_rhiz_evn.let$Treats$Letters),
            crot_root_evn.let$Treats$Letters,
            rev(lent_root_evn.let$Treats$Letters))
cha.let <- c(bulk_cha.let$Treats$Letters,
             crot_rhiz_cha.let$Treats$Letters,
             rev(lent_rhiz_cha.let$Treats$Letters),
             crot_root_cha.let$Treats$Letters,
             rev(lent_root_cha.let$Treats$Letters))

# Create a data.frame that has the mean and standard deviation of each unique group by compartment, plant species, and treatment for each alpha diversity metric #
sel_rich.mnsd <- sel.rich %>%
  group_by(Plants, Treats, Comps) %>%
  summarize(
    sha.mean = mean(Shannon),
    evn.mean = mean(ShaEvn),
    cha.mean = mean(Chao1),
    sha.sd = sd(Shannon),
    evn.sd = sd(ShaEvn),
    cha.sd = sd(Chao1),
    .groups = "drop") # Prevent grouping in the result
sel_rich.mnsd <- as.data.frame(sel_rich.mnsd)
sel_rich.mnsd <- arrange(sel_rich.mnsd, Comps)

# Create a variable that joins plant species and compartment #
for(i in 1:nrow(sel_rich.mnsd)){
  sel_rich.mnsd$PC[i] <- paste0(as.character(sel_rich.mnsd$Plants[i]), "\n", sel_rich.mnsd$Comps[i]) 
}
sel_rich.mnsd$PC <- gsub("A. crotalariae\nBulk Soil", "Both Species\nBulk Soil", sel_rich.mnsd$PC)
sel_rich.mnsd$PC <- factor(sel_rich.mnsd$PC, levels = c("A. crotalariae\nRhizosphere", "A. crotalariae\nRoot Endosphere", "Both Species\nBulk Soil", "A. lentiginosus\nRhizosphere", "A. lentiginosus\nRoot Endosphere"))
sel_rich.mnsd <- cbind(sel_rich.mnsd, sha.let, evn.let, cha.let)

# Plot the results #
library(ggprism)
sha.plot <- ggplot(sel_rich.mnsd, aes(x = PC, y = sha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = sha.mean - sha.sd, ymax = sha.mean + sha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = sha.let, y = sha.mean + sha.sd + 0.1), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Diversity (H)') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration     "), labels = c(expression("0"* mu * "M     "), expression("20" * mu * "M     "), expression("80" * mu * "M     "), expression("100" * mu * "M     ")), values = c('white', "gray", "#808080", "#4D4D4D")) +
  scale_color_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration     "), labels = c(expression("0"* mu * "M     "), expression("20" * mu * "M     "), expression("80" * mu * "M     "), expression("100" * mu * "M     ")), values = c('black', 'black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,7.5),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .)) +
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 24, family = "Liberation Sans", face = "bold"),
        legend.title = element_text(size = 20, face = "plain"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  labs(tag = 'C.')
sha.plot

sel_rich.mnsd$evn.let[4] <- 'ab'
sel_rich.mnsd$evn.let[5] <- 'a'

evn.plot <- ggplot(sel_rich.mnsd, aes(x = PC, y = evn.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = evn.mean - evn.sd, ymax = evn.mean + evn.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = evn.let, y = evn.mean + evn.sd + 0.01), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Shannon Evenness (H/log(S))') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M"), expression("100" * mu * "M")), values = c('white', "gray", "#808080", "#4D4D4D")) +
  scale_color_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M"), expression("100" * mu * "M")), values = c('black', 'black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,1),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .)) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 20, face = "plain"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  labs(tags = 'B.')
evn.plot

cha.plot <- ggplot(sel_rich.mnsd, aes(x = PC, y = cha.mean, fill = Treats, color = Treats)) +
  geom_bar(stat = 'summary', position = 'dodge', width = 0.7) +
  geom_errorbar(aes(ymin = cha.mean - cha.sd, ymax = cha.mean + cha.sd), show.legend = FALSE, position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = cha.let, y = cha.mean + cha.sd + 10), show.legend = FALSE, position = position_dodge(width = 0.7), vjust = 0, size = 4) +
  ylab('Observed ASV Richness (S) ') +
  xlab("") +
  theme_prism() +
  scale_fill_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M"), expression("100" * mu * "M")), values = c('white', "gray", "#808080", "#4D4D4D")) +
  scale_color_manual(name = expression("Na"[2]* "SeO"[4] ~ "Treatment Concentration"), labels = c(expression("0"* mu * "M"), expression("20" * mu * "M"), expression("80" * mu * "M"), expression("100" * mu * "M")), values = c('black', 'black', 'black', 'black')) +
  scale_y_continuous(limits = c(0,710),expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(~ .)) +
  theme(legend.position = 'none',
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 20, face = "plain"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  labs(tag = 'A.')
cha.plot

# Join a larger plot using patchwork #
(cha.plot) /
(evn.plot) /
(sha.plot) +
  plot_layout(guides = 'keep') &
  theme(plot.tag = element_text(size = 20, face = 'bold', color = 'black'))

#### Stacked Histograms ####
library(microbiome)
library(microbiomeutilities)

# Find the most abundant genera within each sample group #
bulk_glom.ps <- tax_glom(bulk.ps, taxrank = 'Genus')
top_bulk.ps <- aggregate_top_taxa2(bulk_glom.ps, top = 19, level = "Genus")
top_bulk.name <- names(sort(taxa_sums(top_bulk.ps), decreasing = TRUE))


crot_rhiz.ps
crot_rhiz_glom.ps <- tax_glom(crot_rhiz.ps, taxrank = "Genus")
top_crot_rhiz.ps <- aggregate_top_taxa2(crot_rhiz_glom.ps, top = 19, level = "Genus")
top_crot_rhiz.name <- names(sort(taxa_sums(top_crot_rhiz.ps), decreasing = TRUE))

crot_root_glom.ps <- tax_glom(crot_root.ps, taxrank = "Genus")
top_crot_root.ps <- aggregate_top_taxa2(crot_root_glom.ps, top = 19, level = "Genus")
top_crot_root.name <- names(sort(taxa_sums(top_crot_root.ps), decreasing = TRUE))
top_crot_root.name <- top_crot_root.name[-3]
top_crot_root.name <- c("Other", top_crot_root.name)

lent_rhiz_glom.ps <- tax_glom(lent_rhiz.ps, taxrank = "Genus")
top_lent_rhiz.ps <- aggregate_top_taxa2(lent_rhiz_glom.ps, top = 19, level = "Genus")
top_lent_rhiz.name <- names(sort(taxa_sums(top_lent_rhiz.ps), decreasing = TRUE))

lent_root.ps
lent_root_glom.ps <- tax_glom(lent_root.ps, taxrank = "Genus")
top_lent_root.ps <- aggregate_top_taxa2(lent_root_glom.ps, top = 19, level = "Genus")

top_lent_root.name <- names(sort(taxa_sums(top_lent_root.ps), decreasing = TRUE))
top_lent_root.name <- top_lent_root.name[-2]
top_lent_root.name <- c("Other", top_lent_root.name)

# Select the most abudant and prevlaent genera across all samples #
all.name <- c(top_bulk.name, top_crot_rhiz.name, top_crot_root.name, top_lent_rhiz.name, top_lent_root.name)

final.name <- c()
all_tab.name <- table(all.name)
for(i in 1:length(all_tab.name)){
  if(all_tab.name[i] >= 3){
    final.name = c(final.name, names(all_tab.name[i]))
  }
}
final.name <- final.name[-8]
final.name <- c("Other", final.name)

other_ps <- function(phy.obj, name.list){
  # function that takes all of the taxa that are not in name list and names them other and produces a new phyloseq object #
  other.ps <- prune_taxa(!as.data.frame(tax_table(phy.obj))$Genus %in% name.list, phy.obj)
  keep.ps <- prune_taxa(as.data.frame(tax_table(phy.obj))$Genus %in% name.list, phy.obj)
  other.tax <- as.data.frame(tax_table(other.ps))
  other.tax$Genus <- "Other"
  final.otu <- rbind(as.data.frame(otu_table(keep.ps)), as.data.frame(otu_table(other.ps)))
  final.tax <- rbind(as.data.frame(tax_table(keep.ps)), other.tax)
  final.tax <- as.matrix(final.tax)
  final.met <- as(sample_data(keep.ps), 'data.frame')
  ot.ps <- phyloseq(otu_table(final.otu, taxa_are_rows = TRUE),
                    tax_table(final.tax),
                    sample_data(final.met))
  return(ot.ps)
}

# Make a color-blind friendly palette #
prev.colr <- c(
  "#D4D4D4", "#4477AA", "#CC6677", "#117733", "#DDCC77",
  "#88CCEE", "#AA4499", "#44AA99", "#882255",
  "#999933", "salmon", "#332288", "#117744"
)

# make a phyloseq object that specifically plots for the stacked histogram for each group #
bulk_hg.ps <- other_ps(bulk_glom.ps, final.name)

bulk_hg.df <- psmelt(bulk_hg.ps)
bulk_hg.df$genera <- factor(bulk_hg.df$Genus, levels = final.name)
bulk_hg.df$order <- factor(bulk_hg.df$Treatment, levels = c('0', '20', '100'))

bulk_hg.colr <- sel.colr[final.name,]

# plot the phyloseq object #
bulk_hg.plot <- ggplot(bulk_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = prev.colr) +
  scale_x_discrete(labels = c(expression("0"* mu * "M"), expression("20"* mu * "M"), expression("100"* mu * "M"))) +
  scale_y_continuous(sec.axis = dup_axis(name = "Bulk Soil")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'none') +
  labs(tag = "A.")

crot_rhiz_hg.ps <- other_ps(crot_rhiz_glom.ps, final.name)

crot_rhiz_hg.df <- psmelt(crot_rhiz_hg.ps)
crot_rhiz_hg.df$genera <- factor(crot_rhiz_hg.df$Genus, levels = final.name)
crot_rhiz_hg.df$order <- factor(crot_rhiz_hg.df$Treatment, levels = c('0', '20', '100'))

crot_rhiz_hg.colr <- sel.colr[final.name,]

crot_rhiz_hg.plot <- ggplot(crot_rhiz_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = prev.colr) +
  scale_x_discrete(labels = c(expression("0"* mu * "M"), expression("20"* mu * "M"), expression("100"* mu * "M"))) +
  scale_y_continuous(sec.axis = dup_axis(name = "Rhizosphere")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'none') +
  labs(tag = "B.")
crot_rhiz_hg.plot

# crot_root hg with selected genera #
crot_root_hg.ps <- other_ps(crot_root_glom.ps, final.name)

crot_root_hg.df <- psmelt(crot_root_hg.ps)
crot_root_hg.df$genera <- factor(crot_root_hg.df$Genus, levels = final.name)
crot_root_hg.df$order <- factor(crot_root_hg.df$Treatment, levels = c('0', '20', '100'))

crot_root_hg.colr <- sel.colr[final.name,]

crot_root_hg.plot <- ggplot(crot_root_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('A. crotalariae') +
  ylab('') +
  scale_fill_manual(values = prev.colr) +
  scale_x_discrete(labels = c(expression("0"* mu * "M"), expression("20"* mu * "M"), expression("100"* mu * "M"))) +
  scale_y_continuous(sec.axis = dup_axis(name = "Root Endosphere")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 20, face = 'bold.italic', family = 'Liberation Sans'),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'none') +
  labs(tag = "D.")
crot_root_hg.plot

# lent_rhiz hg with selected genera #
lent_rhiz_hg.ps <- other_ps(lent_rhiz_glom.ps, final.name)

lent_rhiz_hg.df <- psmelt(lent_rhiz_hg.ps)
lent_rhiz_hg.df$genera <- factor(lent_rhiz_hg.df$Genus, levels = final.name)
lent_rhiz_hg.df$order <- factor(lent_rhiz_hg.df$Treatment, levels = c('0', '20', '80'))

lent_rhiz_hg.colr <- sel.colr[final.name,]

lent_rhiz_hg.plot <- ggplot(lent_rhiz_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('') +
  scale_fill_manual(values = prev.colr) +
  scale_x_discrete(labels = c(expression("0"* mu * "M"), expression("20"* mu * "M"), expression("80"* mu * "M"))) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Rhizosphere')) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'none') +
  labs(tag = "C.")
lent_rhiz_hg.plot

# lent_root hg with selected genera #
lent_root_hg.ps <- other_ps(lent_root_glom.ps, final.name)

lent_root_hg.df <- psmelt(lent_root_hg.ps)
lent_root_hg.df$genera <- factor(lent_root_hg.df$Genus, levels = final.name)
lent_root_hg.df$order <- factor(lent_root_hg.df$Treatment, levels = c('0', '20', '80'))

lent_root_hg.colr <- sel.colr[final.name,]

lent_root_hg.plot <- ggplot(lent_root_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('A. lentiginosus') +
  ylab('') +
  scale_fill_manual(values = prev.colr) +
  scale_x_discrete(labels = c(expression("0"* mu * "M"), expression("20"* mu * "M"), expression("80"* mu * "M"))) +
  scale_y_continuous(sec.axis = dup_axis(name = "Root Endosphere")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 18),
        axis.title = element_text(size = 20, face = 'bold.italic', family = 'Liberation Sans'),
        strip.text = element_text(size =18),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, family = "Liberation Sans", face = 'bold', angle = -90),
        legend.position = 'none') +
  labs(tag = "E.")
lent_root_hg.plot

hg.legend <- ggplot(lent_root_hg.df, aes(x = order, y = Abundance, fill = genera,)) +
  geom_bar(stat='identity', position = 'fill') +
  scale_fill_manual(values = prev.colr) +
  theme_prism() +
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size = 20, family = "Liberation Sans", face = "bold"))

hg.legend <- ggpubr::get_legend(hg.legend)

# Create a plot that joins all of the stacked histograms into one image #
(bulk_hg.plot | hg.legend) /
(crot_rhiz_hg.plot | lent_rhiz_hg.plot) /
(crot_root_hg.plot | lent_root_hg.plot) +
  plot_layout(guides = 'keep') & 
  theme(plot.tag = element_text(size = 20, face = 'bold', color = 'black'))

#### Differential Abundance ####
library(Maaslin2)
system(paste0('mkdir ./maas_results'))

# Create a phyloseq object that contains all samples and all taxa with reads greater than or equal to 1000 #
sel_maas.ps <- subset_taxa(sel.ps, taxa_sums(sel.ps) > 1000)
sel_maas.otu <- as.data.frame(otu_table(sel_maas.ps))
sel_maas.met <- as(sample_data(sel_maas.ps), 'data.frame')

# Perform a differential abudanc analysis using root endosphere samples of A. lentiginosus with 0 uM treatment as the baseline #
lent_root_css.maas <- Maaslin2(input_data = sel_maas.otu,
                                       input_metadata = sel_maas.met,
                                       output = "./maas_results/lent_root_css.maas",
                                       fixed_effects = c("Tres"),
                                       analysis_method = "LM",
                                       normalization = "CSS",
                                       transform = "NONE",
                                       min_prevalence = 0.05,
                                       correction = "BH",
                                       max_significance = 0.05,
                                       reference = c("Tres,0Rol"))
# Save the significant results and then ONLY the results that compare to other A. lentiginosus root endosphere samples #
lent_root_css.res <- lent_root_css.maas$results
lent_root_css.dares <- c()
for(i in 1:nrow(lent_root_css.res)){
  if(lent_root_css.res$value[i] == "8Rol" | lent_root_css.res$value[i] == "2Rol"){
    if(lent_root_css.res$qval[i] < 0.05){
      lent_root_css.dares <- rbind(lent_root_css.dares, lent_root_css.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

lent_root.res <- dplyr::arrange(lent_root_css.dares, value, coef)

# Same as lines 1447-1470 but for the rhizosphere samples of A. lentiginosus #
lent_rhiz_css.maas <- Maaslin2(input_data = sel_maas.otu,
                               input_metadata = sel_maas.met,
                               output = "./maas_results/lent_rhiz_css.maas",
                               fixed_effects = c("Tres"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.05,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Tres,0Rhl"))

lent_rhiz_css.res <- lent_rhiz_css.maas$results
lent_rhiz_css.dares <- c()
for(i in 1:nrow(lent_rhiz_css.res)){
  if(lent_rhiz_css.res$value[i] == "8Rhl" | lent_rhiz_css.res$value[i] == "2Rhl"){
    if(lent_rhiz_css.res$qval[i] < 0.05){
      lent_rhiz_css.dares <- rbind(lent_rhiz_css.dares, lent_rhiz_css.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

lent_rhiz.res <- dplyr::arrange(lent_rhiz_css.dares, value, coef)

# Same as lines 1447-1470 but for the root endosphere samples of A. crotalariae #
crot_root_css.maas <- Maaslin2(input_data = sel_maas.otu,
                               input_metadata = sel_maas.met,
                               output = "./maas_results/crot_root_css.maas",
                               fixed_effects = c("Tres"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.05,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Tres,0Roc"))
crot_root_css.res <- crot_root_css.maas$results
crot_root_css.dares <- c()
for(i in 1:nrow(crot_root_css.res)){
  if(crot_root_css.res$value[i] == "1Roc" | crot_root_css.res$value[i] == "2Roc"){
    if(crot_root_css.res$qval[i] < 0.05){
      crot_root_css.dares <- rbind(crot_root_css.dares, crot_root_css.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

crot_root.res <- dplyr::arrange(crot_root_css.dares, value, coef)

# Same as lines 1447-1470 but for the rhizosphere samples of A. crotalariae #
crot_rhiz_css.maas <- Maaslin2(input_data = sel_maas.otu,
                               input_metadata = sel_maas.met,
                               output = "./maas_results/crot_rhiz_css.maas",
                               fixed_effects = c("Tres"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.05,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Tres,0Rhc"))
crot_rhiz_css.res <- crot_rhiz_css.maas$results
crot_rhiz_css.dares <- c()
for(i in 1:nrow(crot_rhiz_css.res)){
  if(crot_rhiz_css.res$value[i] == "1Rhc" | crot_rhiz_css.res$value[i] == "2Rhc"){
    if(crot_rhiz_css.res$qval[i] < 0.05){
      crot_rhiz_css.dares <- rbind(crot_rhiz_css.dares, crot_rhiz_css.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

crot_rhiz.res <- dplyr::arrange(crot_rhiz_css.dares, value, coef)

# Same as lines 1447-1470 but for the Bulk Soil Samples #
bulk_css.maas <- Maaslin2(input_data = sel_maas.otu,
                               input_metadata = sel_maas.met,
                               output = "./maas_results/bulk_css.maas",
                               fixed_effects = c("Tres"),
                               analysis_method = "LM",
                               normalization = "CSS",
                               transform = "NONE",
                               min_prevalence = 0.05,
                               correction = "BH",
                               max_significance = 0.05,
                               reference = c("Tres,0Buc"))
bulk_css.res <- bulk_css.maas$results
bulk_css.dares <- c()
for(i in 1:nrow(bulk_css.res)){
  if(bulk_css.res$value[i] == "1Buc" | bulk_css.res$value[i] == "2Buc"){
    if(bulk_css.res$qval[i] < 0.05){
      bulk_css.dares <- rbind(bulk_css.dares, bulk_css.res[i, c("feature", "metadata", "value", "coef", "stderr", "pval", "name", "qval", "N", "N.not.zero")]) 
    }
  }
}

bulk.res <- dplyr::arrange(bulk_css.dares, value, coef)

# Save all of the differential abundance results to a csv file #
all.res <- rbind(bulk.res, crot_rhiz.res, crot_root.res,lent_rhiz.res, lent_root.res)
fwrite(all.res, file = "./new_figs/all_res.csv")

#### ASV 1 DA ####
# lent #
# Make a phyloseq object to plot the root endosphere relative abundances #
plant_maas.ps <- subset_samples(sel_maas.ps, Type != "Bulk Soil")
plant_maas.ps <- transform_sample_counts(plant_maas.ps, function(x) x/sum(x))
root_maas.ps <- subset_samples(plant_maas.ps, Type == 'Root Endosphere')
root_maas.met <- as(sample_data(root_maas.ps), 'data.frame')
root_maas.otu <- as.data.frame(otu_table(root_maas.ps))
root_maas.fra <- cbind(root_maas.met, t(root_maas.otu))

# Add annotations for the significance levels as determined from Maaslin2 #
lent_root_nl.anno <- '***'
lent_root_nh.anno <- '***'
crot_root_nl.anno <- 'ns'
crot_root_nh.anno <- 'ns'

# Produce the plot #
library(ggtext)
root_da.plot <- ggplot(root_maas.fra, aes(x = Plants, y = `ASV1(Allomesorhizobium)`)) +
  geom_boxplot(position = 'dodge', aes(fill = `Treats`)) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1, by = 0.2), sec.axis = dup_axis(name = "Root Endosphere")) +
  scale_fill_manual(name = "NaSeO<sub>4</sub> Concentration", values = c("0" = "white", "20" = "gray", "80" = "#808080", "100" = "#4D4D4D"), labels = c("0&mu;M", "20&mu;M", "80&mu;M", "100&mu;M")) +
  scale_x_discrete(labels = c('A. crotalariae', 'A. lentiginosus')) +
  ylab('Relative Abundance of ASV1(Allomesorhizobium) Reads') +
  xlab('') +
  theme_prism() +
  theme(axis.text.x = element_text(color = 'black', size = 32, face = 'bold.italic'),
        axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.y.left = element_text(size = 24),
        axis.ticks.length.y.right = unit(0, 'cm'),
        axis.text.y.right = element_blank(),
        axis.title.y.right = element_text(size = 24, angle = -90),
        legend.text = ggtext::element_markdown(size = 28, family= 'Liberation Sans', face = 'bold'),
        legend.title = element_blank(),
        legend.key.spacing.x = unit(2,'cm'),
        legend.position = 'bottom') + 
  geom_signif(annotation = formatC(lent_root_nh.anno, digits = 3),
              y_position = 1.05, xmin = 1.75, xmax = 2,
              tip_length = c(0.01, 0.01),
              textsize = 12,
              vjust = 0.65) +
  geom_signif(annotation = formatC(lent_root_nl.anno, digits = 3),
              y_position = 1.1, xmin = 1.75, xmax = 2.25,
              tip_length = c(0.01, 0.01),
              textsize = 12,
              vjust = 0.65) +
  geom_signif(annotation = formatC(crot_root_nh.anno, digits = 3),
              y_position = 0.93, xmin = 0.75, xmax = 1,
              tip_length = c(0.01, 0.01),
              textsize = 12) +
  geom_signif(annotation = formatC(crot_root_nl.anno, digits = 3),
              y_position = 1, xmin = 0.75, xmax = 1.25,
              tip_length = c(0.01, 0.01),
              textsize = 12) 
root_da.plot

#### Selenium Quantification ####
# create a data.frame that has all of the selenium quantification data #
sel_levs <- matrix(nrow = 6, ncol = 5)
colnames(sel_levs) <- c("Sample Source", "Concentration", "Sample Mass", "Selenium Level", "Percent DW")
sel_levs <- as.data.frame(sel_levs)
sel_levs[,1] <- c("A. lentiginosus", "A. lentiginosus", "A. crotalariae", "A. crotalariae", "Bulk Soil", "Bulk Soil")
sel_levs[,2] <- c("0", "80", "0", "100", "0", "100")
sel_levs[1,3] <- 0.0751
sel_levs[2,3] <- 0.1561 
sel_levs[3,3] <- 0.2815
sel_levs[4,3] <- 0.2984
sel_levs[1,4] <- 13
sel_levs[2,4] <- 480
sel_levs[3,4] <- 3.6
sel_levs[4,4] <- 490
sel_levs[5,4] <- 2.0 
sel_levs[6,4] <- 2.0
sel_levs$group <- factor(sel_levs$Concentration, levels = c(0, 80, 100))
sel_levs$Plant <- factor(sel_levs$`Sample Source`, levels = c("Bulk Soil", "A. lentiginosus", "A. crotalariae"))

# Plot the data #
ggplot(data = sel_levs, aes(x = Plant, y = `Selenium Level`, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(), aes(fill = group), color = "black") +
  ylab(expression(bold('Leaf Se Concentration (mg ' ~kg^-1~'DW)'))) +
  xlab("") +
  theme_prism() +
  scale_fill_manual(name = "Na<sub>2</sub>SeO<sub>4</sub> Treatment Concentration", labels = c("0&mu;M", "80&mu;M", "100&mu;M"), values = c('white', 'gray', 'black')) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_y_continuous(limits = c(0,500),expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c('Bulk&nbsp;Soil', "<i>A. lentiginosus<i>", "<i>A. crotalariae")) +
  theme(axis.text.y = element_text(color = "black", size = 30, face = 'bold'),
        axis.text.x = ggtext::element_markdown(size = 30, face = 'bold', color = 'black'),
        legend.text = ggtext::element_markdown(size = 26, face = 'bold'),
        legend.title = ggtext::element_markdown(size = 26, face = 'bold'),
        legend.position = 'inside',
        legend.position.inside = c(0.20, 0.8),
        legend.background = element_rect(color = 'black'),
        axis.title.y = element_text(color = 'black', face = 'bold', size = 32)) +
  coord_cartesian(ylim = c(0, 500))

#### Final Touches ####
# Save the phyloseq Object as two data.frames #
sel.otu <- as.data.frame(otu_table(sel.ps))
sel.tax <- as.data.frame(tax_table(sel.ps))
sel.dna <- as.data.frame(refseq(sel.ps))
sel.fra <- cbind(sel.dna, sel.tax, sel.otu)
sel.met <- as(sample_data(sel.ps), 'data.frame')
colnames(sel.fra)[1] <- "DNA Sequence"
colnames(sel.fra)[8] <- "NCBI Best Hit"
sel.fra <- cbind(rownames(sel.fra), sel.fra)
colnames(sel.fra)[1] <- "ASV"


fwrite(sel.fra, file = "~/sel_test2/new_figs/sel_taxXasv.csv", )
fwrite(sel.met, file = "~/sel_test2/new_figs/sel_met.csv")
