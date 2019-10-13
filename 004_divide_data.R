setwd('/home/isglobal.lan/smari/data/smari')
library('minfi')

load('gset_sm_1.RData') # gset
load('summ_exp_sm_res2.RData') # se_res
stopifnot(colnames(gset) == colnames(se_res))
print("Data loaded")

ma <- assays(gset)$Beta # methylation assay
ea <- assays(se_res)$exprs # expression assay


mr <- rowRanges(gset)[,0] # methylation ranges
er <- rowRanges(se_res)[,0] # expression ranges
pheno <- colData(gset)[,c(4,7,10,40:45)]
pheno <- as.data.frame(pheno)
stopifnot(rownames(pheno) == colnames(gset))

# Automatic subsetting + saving
setwd('/home/isglobal.lan/smari/data/smari/data_res2')
save(pheno, file = 'pheno.RData')

count = 0
for (x in 1:22) {
  chr <- paste('chr', x, sep = '')
  mrsub <- mr[seqnames(mr) == chr]
  ersub <- er[seqnames(er) == chr]
  fo <- findOverlaps(mrsub + 5e5, ersub)
  
  cpgnames <- names(mrsub)
  tcnames <- names(ersub)
  cpg <- cpgnames[from(fo)]
  tc <- tcnames[to(fo)]
  
  overlaps <- cbind(cpg, tc)
  masub <- ma[cpgnames,]
  easub <- ea[tcnames,]
  save(overlaps, file = paste('overlaps', x, '.RData', sep = ''))
  save(masub, file = paste('masub', x, '.RData', sep = ''))
  save(easub, file = paste('easub', x, '.RData', sep = ''))
  print(paste(x, nrow(overlaps), sep = ': '))
  count <- count + nrow(overlaps)
}

print(count)

# [1] "1: 1404423"
# [1] "2: 797401"
# [1] "3: 527904"
# [1] "4: 370049"
# [1] "5: 552642"
# [1] "6: 1484057"
# [1] "7: 691998"
# [1] "8: 459977"
# [1] "9: 310614"
# [1] "10: 533132"
# [1] "11: 996292"                                                                        
# [1] "12: 711989"
# [1] "13: 243912"
# [1] "14: 451613"
# [1] "15: 442527"
# [1] "16: 722951"
# [1] "17: 1030021"
# [1] "18: 81119"
# [1] "19: 1074534"
# [1] "20: 328452"
# [1] "21: 148177"
# [1] "22: 308074"

# [1] 13671858