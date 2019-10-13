###############################################
# Genomic Ratio Set (Methylation data)
###############################################

# Loading GRSet
library('minfi')
load("/Lacie_CRW10023/HELIX_preproc/methylation/
Final_data/methylome_subcohort_ComBatSlide_6cells_v3.Rdata")
gset <- methylome_subcohort_ComBatSlide_6cells
pheno <- pData(gset)

# ---- Phenotype
# Transcriptome ID list loading
tran <- read.table("/Lacie_CRW10023/HELIX_preproc/gene_expression/
Final_data/ID_list_transcriptome_v2.txt", header = T)
head(tran)

# Period 1B filtering
tran <- tran[tran$Period != "1B",]
nrow(tran)
table(tran$Period)

# Methylome ID list loading
met <- read.table("/Lacie_CRW10023/HELIX_preproc/methylation/
Final_data/ID_list_methylome_v3.txt", header = T)
head(met)

# Period IB filtering
met <- met[met$Period != "1B",]
nrow(met)
table(met$Period)

# Transcriptome and methylome IDs merging
merg <- merge(tran, met, by = "SampleID")
dim(merg) # 933 11

# Caucasians only filtering
pheno <- pheno[pheno$h_ethnicity_cauc == "yes",]
dim(pheno)
table(pheno$h_ethnicity_cauc)

# Sample IDs filtering on GRset
filter <- pheno$SampleID %in% merg$SampleID
pheno <- pheno[filter,]
dim(pheno)

ids <- rownames(pheno)
save(ids, file = "ids.RData")

# Adding the IDs to the GRset
gset <- gset[,ids]

# ---- Annotation
# Loading annotation
rd <- getAnnotation(gset)
rd <- rd[,-c(5:14, 16, 20:21)]

# Adding relative position group columns
rd$TSS200 <- ifelse(grepl("TSS200", rd$UCSC_RefGene_Group), T, F)
rd$TSS1500 <- ifelse(grepl("TSS1500", rd$UCSC_RefGene_Group), T, F)
rd$UTR5 <- ifelse(grepl("5'UTR", rd$UCSC_RefGene_Group), T, F)
rd$FirstExon <- ifelse(grepl("1stExon", rd$UCSC_RefGene_Group), T, F)
rd$Body <- ifelse(grepl("Body", rd$UCSC_RefGene_Group), T, F)
rd$UTR3 <- ifelse(grepl("3'UTR", rd$UCSC_RefGene_Group), T, F)

# Adding dhs and crom15
crom15 <- read.csv("crom15.csv", header = T, row.names = 1)
crom15 <- crom15[,c(1,7:21)]
dhs <- read.csv("dhs.csv", header = T, row.names = 1)
dhs <- dhs[,c(1,7:19)]

states <- merge(crom15, dhs, by = "HT12v4.ArrayAddress", sort = F)
rd <- merge(rd, states, by.x = "row.names", by.y = "HT12v4.ArrayAddress", sort = F)

# Converting genes column into lists
a <- strsplit(rd$UCSC_RefGene_Name, ";")
a[a == "character(0)"] <- ""
rd$UCSC_RefGene_Name <- a

# ---- Saving
rowData(gset) <- rd
save(gset, file = "gset_sm_1.RData")


###############################################
# Summarized Experiment (Gene expression data)
###############################################

library(Biobase)
library(SummarizedExperiment)

setwd("/Lacie_CRW10023/HELIX_analyses/expr_met_SM")
load("anno_expr_sm_2.RData")
names(anno)
dim(anno)

load("/Lacie_CRW10023/HELIX_preproc/gene_expression/Final_data/
transcriptome_subcohort_f1_residuals3_v2.Rdata")
eset <- transcriptome_subcohort_f1_residuals3
feat <- featureData(eset)
eanno <- pData(feat)
names(eanno)
dim(eanno)

# Filter IDs
load("ids.RData")
eset <- eset[,ids] #832

# Filter Transcript Clusters
# eanno: original expression set annotation
# anno: filtered annotation (chr gl, chr M and haplotypes)
eset <- eset[rownames(anno),] #59130

# Adding annotation
pData(feat) <- anno
featureData(eset) <- feat

# Making Summarized Experiment
se <- makeSummarizedExperimentFromExpressionSet(eset)
save(se, file = "summ_exp_sm_1.RData")