library("minfi")
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/summ_exp_sm_res2.RData") # se_res
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/gset_sm_1.RData") # gset
tcanno <- rowData(se_res)
load("/home/isglobal.lan/smari/data/smari/gset_sm_1.RData")
cpganno <- rowData(gset)
cpganno <- as.data.frame(cpganno[, c(5,9,12,23:28)])
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/overlaps_d.RData")
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/fil.RData") # extra filtering added later
overlaps_d$dif <- ifelse(overlaps_d$strTC == "+", overlaps_d$dif, -overlaps_d$dif)


setwd("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM")
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/data_res2/all_output.RData")
res2 <- all_output
rm(all_output)
res2 <- res2[!(res2$TC %in% fil), ]
res2$FDR <- p.adjust(res2$p.value, method = "fdr")
res2 <- merge(res2, tcanno[,c(1,24)], by.x = "TC",
              by.y = "transcript_cluster_id", sort = F)
save(res2, file = "res2.RData")

res2.data <- merge(res2, overlaps_d[, c(1,2,6,8)], by.x = c("CpG", "TC"),
                   by.y = c("CpG", "TC"), sort = F)
res2.data <- merge(res2.data, cpganno, by.x = "CpG", by.y = "Name", sort = F)
res2.data <- merge(res2.data, tcanno[,c(1,15)], by.x = "TC",
                   by.y = "transcript_cluster_id", sort = F)
save(res2.data, file = "res2.data.RData")




### PLOTS
library(ggplot2)
library(cowplot)
library(Cairo)
setwd("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/plots")
bonf <- 0.05/nrow(res2.data)

# width - p-value
png("a2pos.png", type = "cairo", width = 1024, height = 768)
ggplot(res2.data[res2.data$Estimate >= 0, ],
       aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5, color = "seagreen") +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ggtitle("Width & P-value, res2, effect+")
dev.off()

png("a2neg.png", type = "cairo", width = 1024, height = 768)
ggplot(res2.data[res2.data$Estimate < 0, ],
       aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5, color = "red") +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ggtitle("Width & P-value, res2, effect-")
dev.off()

png("a3pos.png", type = "cairo", width = 1024, height = 768)
ggplot(res3.data[res3.data$Estimate >= 0, ],
       aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5, color = "seagreen") +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ggtitle("Width & P-value, res3, effect+")
dev.off()

png("a3neg.png", type = "cairo", width = 1024, height = 768)
ggplot(res2.data[res3.data$Estimate < 0, ],
       aes(x = dif, y = -log10(p.value))) +
  geom_point(alpha = 1/10, size = 1.5, color = "red") +
  theme_classic(base_size = 20) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ggtitle("Width & P-value, res3, effect-")
dev.off()

# width - estimate (+ density width)
bonf <- 0.05/nrow(res2.data)
res2.data$sig <- ifelse(res2.data$p.value < bonf, "sig", "nosig")
res3.data$sig <- ifelse(res3.data$p.value < bonf, "sig", "nosig")

s2 <- ggplot(res2.data, aes(x = dif, y = Estimate, color = sig)) +
  geom_point(alpha = 1/50, size = 1.5) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
s2z <- ggplot(res2.data, aes(x = dif, y = Estimate, color = sig)) +
  geom_point(alpha = 1/50, size = 1.5) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt")) +
  ylim(-15, 15)
s2x <- ggplot(res2.data, aes(x = dif, color = sig)) + 
  geom_density() + theme_void() +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(10, 0, 0, 10), "pt"))
b2 <- plot_grid(s2x, s2, s2z, ncol = 1, align = "v",
                rel_heights = c(1.5, 3, 3))
ggsave("b2.png", b2, w = 16, h = 16, type = "cairo-png")

s3 <- ggplot(res3.data, aes(x = dif, y = Estimate, color = sig)) +
  geom_point(alpha = 1/50, size = 1.5) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
s3z <- ggplot(res3.data, aes(x = dif, y = Estimate, color = sig)) +
  geom_point(alpha = 1/50, size = 1.5) +
  theme_classic(base_size = 22) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt")) +
  ylim(-15, 15)
s3x <- ggplot(res3.data, aes(x = dif, color = sig)) + 
  geom_density() + theme_void() +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(10, 0, 0, 10), "pt"))
b3 <- plot_grid(s3x, s3, s3z, ncol = 1, align = "v",
                rel_heights = c(1.5, 3, 3))
ggsave("b3.png", b3, w = 16, h = 16, type = "cairo-png")


summary(res1.data$Estimate)
cor(abs(res1.data$Estimate), abs(res1.data$dif))

summary(res2.data$Estimate)
cor(abs(res2.data$Estimate), abs(res2.data$dif))

summary(res3.data$Estimate)
cor(abs(res3.data$Estimate), abs(res3.data$dif))

# VennDiagram
bonf <- 0.05/nrow(res1)
library(limma)
sig1 <- res1[res1$p.value < bonf, ] # 15861
sig2 <- res2[res2$p.value < bonf, ] # 15403
sig3 <- res3[res3$p.value < bonf, ] # 21241
set1 <- paste(sig1$CpG, sig1$TC, sep = ":")
set2 <- paste(sig2$CpG, sig2$TC, sep = ":")
set3 <- paste(sig3$CpG, sig3$TC, sep = ":")
universe <- sort(unique(c(set1, set2, set3)))
Counts <- matrix(0, nrow = length(universe), ncol = 3)
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% set1
  Counts[i,2] <- universe[i] %in% set2
  Counts[i,3] <- universe[i] %in% set3
}
colnames(Counts) <- c("set1","set2","set3")

cols <- c("Red", "Green", "Blue")
png("venn.png", type = "cairo")
vennDiagram(vennCounts(Counts), circle.col = cols)
dev.off()

# estimate - p-value
res2$cr <- NULL
res2$cr <- "75-100"
res2$cr[res2$CallRate < 75] <- "50-75"
res2$cr[res2$CallRate < 50] <- "25-50"
res2$cr[res2$CallRate < 25] <- "1-25"
l <- c("1-25", "25-50", "50-75", "75-100")
res2$cr <- factor(res2$cr, levels = l)

png("c2_1.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res2, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("blue", "grey", "grey", "grey")) +
  scale_size_manual(values = c(2, 1, 1, 1)) +
  scale_shape_manual(values = c(16, 3, 3, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res2, callrate 1-25")
dev.off()

png("c2_2.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res2, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "blue", "grey", "grey")) +
  scale_size_manual(values = c(1, 2, 1, 1)) +
  scale_shape_manual(values = c(3, 16, 3, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res2, callrate 25-50")
dev.off()

png("c2_3.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res2, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "grey", "blue", "grey")) +
  scale_size_manual(values = c(1, 1, 2, 1)) +
  scale_shape_manual(values = c(3, 3, 16, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res2, callrate 50-75")
dev.off()

png("c2_4.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res2, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "grey", "grey", "blue")) +
  scale_size_manual(values = c(1, 1, 1, 2)) +
  scale_shape_manual(values = c(3, 3, 3, 16)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res2, callrate 75-100")
dev.off()

res3$cr <- NULL
res3$cr <- "75-100"
res3$cr[res3$CallRate < 75] <- "50-75"
res3$cr[res3$CallRate < 50] <- "25-50"
res3$cr[res3$CallRate < 25] <- "1-25"
l <- c("1-25", "25-50", "50-75", "75-100")
res3$cr <- factor(res3$cr, levels = l)

png("c3_1.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res3, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("blue", "grey", "grey", "grey")) +
  scale_size_manual(values = c(2, 1, 1, 1)) +
  scale_shape_manual(values = c(16, 3, 3, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res3, callrate 1-25")
dev.off()

png("c3_2.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res3, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "blue", "grey", "grey")) +
  scale_size_manual(values = c(1, 2, 1, 1)) +
  scale_shape_manual(values = c(3, 16, 3, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res3, callrate 25-50")
dev.off()

png("c3_3.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res3, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "grey", "blue", "grey")) +
  scale_size_manual(values = c(1, 1, 2, 1)) +
  scale_shape_manual(values = c(3, 3, 16, 3)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res3, callrate 50-75")
dev.off()

png("c3_4.png", type = "cairo", width = 1024, height = 768)
ggplot(data = res3, aes(x = Estimate/100, y = -log10(p.value),
                        color = cr, shape = cr)) +
  geom_point(alpha = 1/4, aes(size = cr)) +
  scale_color_manual(values = c("grey", "grey", "grey", "blue")) +
  scale_size_manual(values = c(1, 1, 1, 2)) +
  scale_shape_manual(values = c(3, 3, 3, 16)) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  xlim(-1, 1) + ylim(0, 150) +
  ggtitle("Est & P-value, res3, callrate 75-100")
dev.off()


# width - significant p.values
res2.data$sigeff <- ifelse(res2.data$p.value < bonf,
                           ifelse(res2.data$Estimate >= 0, "Pos", "Neg"),
                           "NoSig")
res2.data$sigeff <- as.factor(res2.data$sigeff)
#  Neg    NoSig      Pos
# 9477 13600479     5926

res3.data$sigeff <- ifelse(res3.data$p.value < bonf,
                           ifelse(res3.data$Estimate >= 0, "Pos", "Neg"),
                           "NoSig")
res3.data$sigeff <- as.factor(res3.data$sigeff)
#   Neg    NoSig      Pos
# 13887 13594641     7354

png("d2.png", type = "cairo", width = 1024, height = 768)
ggplot(res2.data, aes(x = dif,  color = sigeff)) +
  geom_density() +
  theme_classic(base_size = 22) +
  scale_color_manual(values=c("red", "black", "green"))
dev.off()

png("d3.png", type = "cairo", width = 1024, height = 768)
ggplot(res3.data, aes(x = dif,  color = sigeff)) +
  geom_density() +
  theme_classic(base_size = 22) +
  scale_color_manual(values=c("red", "black", "green"))
dev.off()


# sign p.values/chromosome
library(reshape2)
chr <- sapply(1:22, function(x) paste("chr", x, sep = ""))
count.chr <- function(x, d){
  subd <- d[d$chrCpG == x, ]
  a <- nrow(subd[subd$FDR < 0.05, ])
  b <- nrow(subd[subd$FDR >= 0.05, ])
  print(x)
  return(c(a/(a+b)*100, b/(a+b)*100))
}

data2 <- melt(t(sapply(chr, count.chr, d = res2.data)))
data3 <- melt(t(sapply(chr, count.chr, d = res3.data)))

png("e2.png", type = "cairo", width = 1024, height = 768)
ggplot(data2, aes(x = Var1, y = round(value, 2), fill = as.factor(Var2))) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = round(value,2)), vjust = 1.5,
            color = "black", size = 4, position = "stack") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(90,100))
dev.off()

png("e3.png", type = "cairo", width = 1024, height = 768)
ggplot(data3, aes(x = Var1, y = round(value, 2), fill = as.factor(Var2))) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = round(value,2)), vjust = 1.5,
            color = "black", size = 4, position = "stack") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(90,100))
dev.off()


# p-values between results
cor(-log10(res1.data$p.value), -log10(res2.data$p.value)) # 0.9976679
cor(-log10(res1.data$p.value), -log10(res3.data$p.value)) # 0.8591621
cor(-log10(res2.data$p.value), -log10(res3.data$p.value)) # 0.8601841

d <- data.frame("d1" = -log10(res1.data$p.value),
                "d2" = -log10(res2.data$p.value),
                "d3" = -log10(res3.data$p.value))

png("f_1_2.png", type = "cairo", height = 720, width = 720)
ggplot(d, aes(x = d1, y = d2)) +
  geom_point(alpha = 1/200, size = 1.5) + 
  theme_classic(base_size = 22)
dev.off()

png("f_1_3.png", type = "cairo", height = 720, width = 720)
ggplot(d, aes(x = d1, y = d3)) +
  geom_point(alpha = 1/200, size = 1.5) +
  theme_classic(base_size = 22)
dev.off()

png("f_2_3.png", type = "cairo", height = 720, width = 720)
ggplot(d, aes(x = d2, y = d3)) +
  geom_point(alpha = 1/200, size = 1.5) +
  theme_classic(base_size = 22)
dev.off()

sig12 <- res1.data$p.value < bonf | res2.data$p.value < bonf
sig13 <- res1.data$p.value < bonf | res3.data$p.value < bonf
sig23 <- res2.data$p.value < bonf | res3.data$p.value < bonf

png("f_1_2s.png", type = "cairo", height = 720, width = 720)
ggplot(d[sig12,], aes(x = d1, y = d2)) +
  geom_point(alpha = 1/200, size = 1.5) + 
  theme_classic(base_size = 22) +
  ggtitle("P-values comparison, sig only")
dev.off()

png("f_1_3s.png", type = "cairo", height = 720, width = 720)
ggplot(d[sig13, ], aes(x = d1, y = d3)) +
  geom_point(alpha = 1/200, size = 1.5) +
  theme_classic(base_size = 22) +
  ggtitle("P-values comparison, sig only")
dev.off()

png("f_2_3s.png", type = "cairo", height = 720, width = 720)
ggplot(d[sig23,], aes(x = d2, y = d3)) +
  geom_point(alpha = 1/200, size = 1.5) +
  theme_classic(base_size = 22) +
  ggtitle("P-values comparison, sig only")
dev.off()

# estimates between results
cor(res1.data$Estimate, res2.data$Estimate) # 0.9934421
cor(res1.data$Estimate, res3.data$Estimate) # 0.9301162
cor(res2.data$Estimate, res3.data$Estimate) # 0.9379744

e <- data.frame("d1" = res1$Estimate,
                "d2" = res2$Estimate,
                "d3" = res3$Estimate)

png("g_1_2.png", type = "cairo", height = 720, width = 720)
ggplot(e, aes(x = d1, y = d2)) +
  geom_point(alpha = 1/20, size = 1.5) +
  theme_classic(base_size = 22) +
  xlim(-40, 40) + ylim(-40, 40)
dev.off()

png("g_1_3.png", type = "cairo", height = 720, width = 720)
ggplot(e, aes(x = d1, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) +
  theme_classic(base_size = 22) +
  xlim(-40, 40) + ylim(-40, 40)
dev.off()

png("g_2_3.png", type = "cairo", height = 720, width = 720)
ggplot(e, aes(x = d2, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  theme_classic(base_size = 22) +
  xlim(-40, 40) + ylim(-40, 40)
dev.off()

sig12 <- res1.data$p.value < bonf | res2.data$p.value < bonf
sig13 <- res1.data$p.value < bonf | res3.data$p.value < bonf
sig23 <- res2.data$p.value < bonf | res3.data$p.value < bonf
e$s12 <- ifelse(res1.data$p.value < bonf & res2.data$p.value < bonf,
                "Both", ifelse(res1.data$p.value < bonf, "res1", "res2"))
e$s13 <- ifelse(res1.data$p.value < bonf & res3.data$p.value < bonf,
                "Both", ifelse(res1.data$p.value < bonf, "res1", "res3"))
e$s23 <- ifelse(res2.data$p.value < bonf & res3.data$p.value < bonf,
                "Both", ifelse(res2.data$p.value < bonf, "res2", "res3"))

png("g_1_2s.png", type = "cairo", height = 720, width = 720)
ggplot(e[sig12,], aes(x = d1, y = d2, color = s12)) +
  geom_point(alpha = 1/20, size = 1.5) +
  theme_classic(base_size = 22) +
  xlim(-40, 40) + ylim(-40, 40)
dev.off()

png("g_1_3s.png", type = "cairo", height = 720, width = 720)
ggplot(e[sig13,], aes(x = d1, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) +
  theme_classic(base_size = 22)+
  xlim(-40, 40) + ylim(-40, 40)
dev.off()

png("g_2_3s.png", type = "cairo", height = 720, width = 720)
ggplot(e[sig23,], aes(x = d2, y = d3)) +
  geom_point(alpha = 1/20, size = 1.5) + 
  theme_classic(base_size = 22) +
  xlim(-40, 40) + ylim(-40, 40)
dev.off()


### CpG vs TC (res1 and res2)
library(minfi)
library(reshape2)
a <- res1[res1$p.value < bonf, c(1,2,3,5)]
b <- res2[res2$p.value < bonf, c(1,2,3,5)]
sig <- merge(a, b, by = c("CpG", "TC"), sort = FALSE)
sig <- sig[sample(nrow(sig), 8),]
sig$pair <- paste(sig$CpG, sig$TC, sep = ":")

a <- res1[res1$p.value > bonf, c(1,2,3,5)]
b <- res2[res2$p.value > bonf, c(1,2,3,5)]
nosig <- merge(a, b, by = c("CpG", "TC"), sort = FALSE)
nosig <- nosig[sample(nrow(nosig), 8),]
nosig$pair <- paste(nosig$CpG, nosig$TC, sep = ":")

load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/summ_exp_sm_res1.RData") # se_res
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/gset_sm_res1.RData") # gset_res
r1exp <- assays(se_res)$expr
r1met <- assays(gset_res)$Beta
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/summ_exp_sm_res2.RData") # se_res
load("/home/isglobal.lan/smari/data/WS_HELIX/HELIX_analyses/expr_met_SM/gset_sm_1.RData") # gset
r2exp <- assays(se_res)$expr
r2met <- assays(gset)$Beta
rm(se_res, gset, gset_res)

# Sig Res 1
e1 <- melt(r1exp[sig$TC, ])
m1 <- melt(r1met[sig$CpG, ])
d1 <- cbind(m1, e1)
colnames(d1) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")
d1$pair <- paste(d1$CpG, d1$TC, sep = ":")

# No Sig Res 1
e2 <- melt(r1exp[nosig$TC, ])
m2 <- melt(r1met[nosig$CpG, ])
d2 <- cbind(m2, e2)
colnames(d2) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")
d2$pair <- paste(d2$CpG, d2$TC, sep = ":")

# Sig Res 2
e3 <- melt(r2exp[sig$TC, ])
m3 <- melt(r2met[sig$CpG, ])
d3 <- cbind(m3, e3)
colnames(d3) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")
d3$pair <- paste(d3$CpG, d3$TC, sep = ":")

# No Sig Res 2
e4 <- melt(r2exp[nosig$TC, ])
m4 <- melt(r2met[nosig$CpG, ])
d4 <- cbind(m4, e4)
colnames(d4) <- c("CpG", "met.sample", "met", "TC", "exp.sample", "exp")
d4$pair <- paste(d4$CpG, d4$TC, sep = ":")

## PLOTTING
spair <- unique(d1$pair) # = unique(d3$pair)

myplot <- function(i) {
  print(
  ggplot(d1[d1$pair == spair[i],], aes(x = met, y = exp)) +
    geom_point(alpha = 0.5) +
    theme_classic(base_size = 28) + ggtitle(spair[i]) +
    geom_smooth(method = 'lm') + 
    labs(subtitle = paste(sig[sig$pair == spair[i], 4],
                          sig[sig$pair == spair[i], 6],
                          sep = "\n"))
  )
  print(
  ggplot(d3[d3$pair == spair[i],], aes(x = met, y = exp)) +
    geom_point(alpha = 0.5) +
    theme_classic(base_size = 28) + ggtitle(spair[i]) +
    geom_smooth(method = 'lm') + 
    labs(subtitle = paste(sig[sig$pair == spair[i], 5],
                          sig[sig$pair == spair[i], 3],
                          sep = "\n"))
  )
}
png("h.png",type = "cairo", width = 720, height = (720/2)*8)
par(mfrow = c(8,2))
for (i in 1:8) {myplot(i)}
dev.off()

# TCs call rate vs TCs min p-value
pv <- tapply(res1$p.value, res1$TC, min)
pv <- as.data.frame(pv)
t <- res1[,c(1,9)]
t <- t[!(duplicated(t$TC)), ]
df1 <- merge(pv, t, by.x = "row.names", by.y = "TC")
est <- tapply(abs(res1$Estimate), res1$TC, max)
est <- as.data.frame(est)
df1 <- merge(df1, est, by.x = "Row.names", by.y = "row.names")

pv <- tapply(res2$p.value, res2$TC, min)
pv <- as.data.frame(pv)
t <- res2[,c(1,9)]
t <- t[!(duplicated(t$TC)), ]
df2 <- merge(pv, t, by.x = "row.names", by.y = "TC")
est <- tapply(abs(res2$Estimate), res2$TC, max)
est <- as.data.frame(est)
df2 <- merge(df2, est, by.x = "Row.names", by.y = "row.names")

pv <- tapply(res3$p.value, res3$TC, min)
pv <- as.data.frame(pv)
t <- res1[,c(1,9)]
t <- t[!(duplicated(t$TC)), ]
df3 <- merge(pv, t, by.x = "row.names", by.y = "TC")
est <- tapply(abs(res3$Estimate), res3$TC, max)
est <- as.data.frame(est)
df3 <- merge(df3, est, by.x = "Row.names", by.y = "row.names")

cor(-log10(df1$pv), df1$CallRate) # 0.1516837
cor(-log10(df2$pv), df2$CallRate) # 0.150663
cor(-log10(df3$pv), df3$CallRate) # 0.1755687

png("i_1.png", type = "cairo", width = 720, height = 720)
i1 <- ggplot(df1, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i1z <- ggplot(df1, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ylim(0, 10) +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i1x <- ggplot(df1, aes(x = CallRate)) + 
  geom_density() + theme_void() +
  theme(plot.margin = unit(c(10, 0, 0, 10), "pt"))
plot_grid(i1x, i1, i1z, ncol = 1, align = "v", rel_heights = c(1.5, 3, 3))
dev.off()

png("i_2.png", type = "cairo", width = 720, height = 720)
i2 <- ggplot(df2, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i2z <- ggplot(df2, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ylim(0, 10) +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i2x <- ggplot(df2, aes(x = CallRate)) + 
  geom_density() + theme_void() +
  theme(plot.margin = unit(c(10, 0, 0, 10), "pt"))
plot_grid(i2x, i2, i2z, ncol = 1, align = "v", rel_heights = c(1.5, 3, 3))
dev.off()

png("i_3.png", type = "cairo", width = 720, height = 720)
i3 <- ggplot(df3, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i3z <- ggplot(df3, aes(x = CallRate, y = -log10(pv))) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22) +
  geom_hline(yintercept = -log10(bonf), color = "blue") +
  ylim(0, 10) +
  theme(plot.margin = unit(c(0, 0, 10, 10), "pt"))
i3x <- ggplot(df3, aes(x = CallRate)) + 
  geom_density() + theme_void() +
  theme(plot.margin = unit(c(10, 0, 0, 10), "pt"))
plot_grid(i3x, i3, i3z, ncol = 1, align = "v", rel_heights = c(1.5, 3, 3))
dev.off()


# TCs call rate vs TCs max abs effect

cor(df1$est, df1$CallRate) # 0.1521714
cor(df2$est, df2$CallRate) # 0.1546927
cor(df3$est, df3$CallRate) # 0.1665243

png("j_1.png", type = "cairo", width = 720, height = 720)
j1p <- ggplot(df1[df1$pv < bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
j1n <- ggplot(df1[df1$pv >= bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
plot_grid(j1p, j1n, ncol = 1, align = "v", rel_heights = c(3, 3))
dev.off()

png("j_2.png", type = "cairo", width = 720, height = 720)
j2p <- ggplot(df2[df2$pv < bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
j2n <- ggplot(df2[df2$pv >= bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
plot_grid(j2p, j2n, ncol = 1, align = "v", rel_heights = c(3, 3))
dev.off()

png("j_3.png", type = "cairo", width = 720, height = 720)
j3p <- ggplot(df3[df3$pv < bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
j3n <- ggplot(df3[df3$pv >= bonf, ], aes(x = CallRate, y = est)) +
  geom_point(alpha = 0.1) +
  theme_classic(base_size = 22)
plot_grid(j3p, j3n, ncol = 1, align = "v", rel_heights = c(3, 3))
dev.off()