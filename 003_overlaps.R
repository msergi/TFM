#####################################################################
# FINDOVERLAPS

library(minfi)

# load("summ_exp_sm_2.RData")
load("summ_exp_sm_res2.RData")
# se

# load("gset_sm_1.RData")
load("gset_sm_res1.RData")
# gset

# Find CpG - TC pairs 500kb away
fo <- findOverlaps(rowRanges(gset_res)[,0] + 5e5, rowRanges(se_res)[,0])
# 13,671,858 CpG - TC pairs

# Prepare some data related to CpGs and TCs to use it later
froms <- from(fo)
tos <- to(fo) 
cpg <- rownames(gset_res)
tc <- rownames(se_res)
froms2 <- cpg[froms]
tos2 <- tc[tos]

overlaps <- data.frame("CpG" = froms2, "TC" = tos2)
overlaps <- as.matrix(overlaps)

gr_met <- as.data.frame(rowRanges(gset_res)[,0])
gr_se <- as.data.frame(rowRanges(se_res)[,0])

overlaps_d <- merge(overlaps, gr_se[,c(1,2,5)], by.x = "TC", by.y = "row.names")
x <- merge(overlaps_d, gr_met[,c(1,2)], by.x = "CpG", by.y = "row.names")
overlaps_d$seqnames.x <- as.character(overlaps_d$seqnames.x)
overlaps_d$seqnames.y <- as.character(overlaps_d$seqnames.y)
overlaps_d$strand <- as.character(overlaps_d$strand)

names(overlaps_d)[3:7] <- c("chrTC", "TSS_Affy", "strTC", "chrCpG", "posCpG")
overlaps_d$dif <- overlaps_d$TSS_Affy - overlaps_d$posCpG
summary(overlaps_d$dif)
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -500000.0 -235681.0     -49.0      56.9  236375.0  500000.0

save(overlaps, file = "overlaps.RData")
save(overlaps_d, file = "overlaps_d.RData")