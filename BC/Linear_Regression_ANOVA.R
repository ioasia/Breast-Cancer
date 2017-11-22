 ### Manual Linear Regression and ANOVA

## Define ESR status with mclust
# install.packages('mclust')
# library(mclust)
# estr_clust <- Mclust(mrna["ESR1",])
## plot(density(mrna["ESR1",]))

## Define ESR status with logistic regression
estr_status <- meta_data[tumours, "ER"]
id_To_Predict <- which(estr_status == "NA")
id_Train <- which(estr_status != "NA")
estr_status <- ifelse(meta_data[tumours, "ER"] == "neg", 0, 1)
p <- prot["ESR1",id_Train]

model <- glm(formula = estr_status[id_Train] ~ p, family = binomial)
# summary(model)
predict.glm(model, newdata = data.frame(p=prot["ESR1",id_To_Predict]) ,type = "response")
estr_status[id_To_Predict] <- 1

## Use residual dataset based on PCA
compute.residual.matrix <- function(data, type) {
  if (type=="PCA") {
    PCA_data = prcomp(t(data), scale. = FALSE)

  # Estimate correlation of subtype with 1st principal component
  # cor(PCA_data$x[,1],as.numeric(setSubtypes))

  # Plot 1st and 2nd Principal Component
  # PCAPlots(data)

  model <- lm(t(data) ~ PCA_data$x[,1])
  residual_matrix <- t(model$residuals)
  } 
  
  # Use residual dataset based on mRNA ESR1 expression
  else if (type =="ESR" && all(data==mrna)) {
  model <- lm(t(data) ~ mrna["ESR1",])
  residual_matrix <- t(model$residuals)
  }
  
  # Use residual dataset based on Protein ESR1 expression
  else if (type =="ESR" && all(data==prot)) {
    model <- lm(t(data) ~ prot["ESR1",])
    residual_matrix <- t(model$residuals)
  }
  
  # Use residual dataset based on  ESR1 status
  else if (type =="ESR_binary") {
    model <- lm(t(data) ~ estr_status)
    residual_matrix <- t(model$residuals)
  }
  
  # Adjust Protein expression for mRNA transcript
  else if (type=="mrna_adjust" && all(data==prot)) {
    residual_matrix <- t(sapply(genes, function(i) {
      model <- lm(data[i,] ~ mrna[i,])
      res <- model$residuals
    }))
  } else(residual_matrix=NA)
  return(residual_matrix)
}

## ANOVA without permutation 
do.anova.init <- function(data,gene) {
  anova_init <- aov(data[gene, ] ~ factor(cna_bin[gene,]))
  anova_init <- summary(anova_init)
  anova_init_pval <- anova_init[[1]][[5]][1]
  anova_init_pval
}

# ANOVA on ranked values
do.anova.rank <- function(data,gene) {
  mod <- lm(rank(data[gene, ]) ~ factor(cna_bin[gene, ]))
  anova_rank <- anova(mod)
  anova_rank$`Pr(>F)`[[1]]
}

do.anova.adj <- function(data,gene) {
  mod <- lm(rank(data[gene, ]) ~  rank(mrna_res[gene, ]) + factor(cna_bin[gene, ]))
  anova_rank <- anova(mod)
  anova_rank$`Pr(>F)`[[2]]
}

## Anova with permutation
do.anova.permutation <- function (data,gene) {
  cna_status <- revalue(factor(cna_bin[gene,]), c("-1"="Loss","0"="Neutral","1"="Gain"))
  anova_perm <- aovp(data[gene, ] ~ cna_status)
  # sceffe_test <- ScheffeTest(anova_perm)
  # n <- names(sceffe_test$cna_status[,4])
  # n <- n[!n %in% "Gain-Loss"]
  # n <- if (length(n)==0) {NA} else {n}
  anova_perm <- summary(anova_perm)
  anova_perm_pval <- anova_perm[[1]][[5]][1]
  return(anova_perm_pval)
}

## Kruskal-Wallis test
do.kruskal.wall.test <- function (data,gene) {
  test = kruskal.test(data[gene, ] ~ factor(cna_bin[gene,])) 
  pval = test$p.value
  return(pval)
}
  

# Plot residuals
# anova <- aov(mrna["ERBB2", ] ~ factor(cna_bin["ERBB2",]))
# pdf("qqplot_ERBB2_mRNA_CNA.pdf", width=10, height=6, paper="special")
# qqnorm(resid(anova), main="Normal Q-Q Plot ERBB2")
# qqline(resid(anova), col="red")
# dev.off()

find.significant.genes <- function(i, testName, data, genes, method_adj) {
  tab <- data.frame(do.call(rbind,lapply(genes, i, data=data)))
  colnames(tab) <- "pval"
  rownames(tab) <- genes
  pval_adj <- p.adjust(tab$pval, method=method_adj)
  res <- data.frame(genes, tab$pval, pval_adj)
  res <- res[order(res$pval_adj),]
  colnames(res) = c("genes",paste('pval', testName, sep="_"), paste('qval', testName, sep="_"))
  return(res)
}

calculate.wilcox.pval.fc <- function(omicsType, binary_CNA, levelName, level) {
  testTable <- data.frame(do.call(rbind,lapply(genes, difExprByCNstatus, omicsType = omicsType, binary_CNA = binary_CNA, level = level)))
  testTable$pvalWilcoxonAdjusted = p.adjust(testTable[,1], method = "BH")
  rownames(testTable) <- genes
  colnames(testTable) = c(paste("pvalWilcoxon",levelName,sep="_"),paste("FoldChange",levelName,sep="_"),
                          paste("Effect_Size", levelName, sep="_"),
                          paste("adjustedPvalueWilcoxon",levelName,sep="_"))
  return(testTable)
}

# Correlation function

molecular.species.correlation <- function(species1,species2, cor_method) {
  cordata = corAndPvalue(t(species1),t(species2), method=cor_method, use= "all.obs")
  a = melt(cordata$cor)
  b = melt(cordata$p)
  cordata = data.frame(CNA=a[,1], Gene=a[,2], correlation=a[,3], pvalue=b[,3])
  cordata$adj.pvalue = p.adjust(cordata$pvalue, method = "fdr")
  # adjcordata = cordata[cordata$adj.pvalue<=0.05,]
  # identAdjcordata = adjcordata[adjcordata$CNA==adjcordata$mRNA,]
  # identAdjcordata = identAdjcordata[order(identAdjcordata$adj.pvalue),]
  cordataCis = cordata[cordata$CNA==cordata$Gene,]; rownames(cordataCis) = NULL
  cordataCis = cordataCis[!is.na(cordataCis$correlation),]
  rownames(cordataCis) = cordataCis$CNA
  return(cordataCis)
}


make.summary.table <- function(data, molecular_view) {
  anova_sig <- find.significant.genes(do.anova.init, "ANOVA", data, genes,"fdr")
  anova_perm_sig <- find.significant.genes(do.anova.permutation, "ANOVA_Perm", data, genes, "fdr")
  kruskal_sig <- find.significant.genes(do.kruskal.wall.test, "Kruskal", data, genes, "fdr")

  wilcox_gains <- calculate.wilcox.pval.fc(data,cna_bin,"Gains", 1)
  wilcox_losses <- calculate.wilcox.pval.fc(data,cna_bin,"Losses", -1)

  wilcox <- cbind(wilcox_gains, wilcox_losses)
  wilcox$genes <- rownames(wilcox)

  pearson <- molecular.species.correlation(cna_log, data, 'pearson')
  spearman <- molecular.species.correlation(cna_log, data, 'spearman')
  
  cor_table <- merge(pearson,spearman, by=c("CNA","Gene"), all=TRUE)
  cor_table <- cor_table[, -1]
  colnames(cor_table) <- c("genes", "Pearson_correlation", "Pearson_pval", "Pearson_adj_pval",
                         "Spearman_correlation", "Spearman_pval", "Spearman_adj_pval")
  
  res_all <- merge(merge(merge(merge(anova_sig, anova_perm_sig, by="genes"), kruskal_sig, by="genes"),wilcox,by="genes"),
                 cor_table, by="genes")

  # res_all <- res_all[order(res_all$qval_ANOVA),]
  rownames(res_all) = NULL
  colnames(res_all) <- c( "genes" , paste(molecular_view, colnames(res_all)[-1], sep="_"))
  return(res_all)
}
#--------------------------------------------------------------------------------------------------------
## ANOVA on ranked values as an omnibus test, preceeding pairwise Wilcoxon test
mrna_res <- compute.residual.matrix (mrna,"ESR_binary")
prot_res <- compute.residual.matrix (prot,"ESR_binary")

# Correlation of 1st PC with Estrogen status
estr_stat_1st_PC_mRNA <- cor(prcomp(t(mrna))$x[,1],estr_status)
estr_stat_1st_PC_Protein<- cor(prcomp(t(prot))$x[,1],estr_status)

mrna_pvals <- sapply(genes, do.anova.rank, data = mrna_res)
prot_pvals <- sapply(genes, do.anova.rank, data = prot_res)

mrna_pvals_adj <- p.adjust(mrna_pvals, method = "BH")
prot_pvals_adj <- p.adjust(prot_pvals, method = "BH")

perc_cna_mrna <- round(sum(mrna_pvals_adj<=0.1)/nrow(mrna) * 100,2)
perc_cna_prot <- round(sum(prot_pvals_adj<=0.1)/nrow(prot) * 100,2)

# sum(mrna_pvals_adj<=0.1 & prot_pvals_adj<=0.1)/ sum(mrna_pvals_adj<=0.1) * 100
# sum(mrna_pvals_adj<=0.1 & prot_pvals_adj<=0.1)/ sum(prot_pvals_adj<=0.1) * 100

pval_dat <- data.frame(gene = genes,mrna = -log10(mrna_pvals_adj), prot = -log10(prot_pvals_adj),
                       sig = ifelse(rownames(mrna) %in% intCl754, "yes", "no"))

scatter <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
  geom_rect(aes(xmin=4, xmax=5.5, ymin=2, ymax=3.5), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
  geom_rect(aes(xmin=4, xmax=5.5, ymin=0, ymax=1), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
  geom_rect(aes(xmin=0, xmax=1, ymin=2, ymax=3.5), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
  
  geom_point() + 
  theme_bw() + 
  scale_y_continuous(name =  "Protein (-log10)", expand = c(0.05, 0.05), breaks = seq(0,4,1),
                     labels = seq(0,4,1)) +
  scale_x_continuous(name = "mRNA (-log10)", expand = c(0.05, 0.05), breaks = seq(0,6,1),
                     labels = seq(0,6,1)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = alpha('grey',0.50),yes ='blue')) +
  # ggtitle("mRNA/Protein ANOVA adjusted P values") +
  guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 3))) + 
  geom_vline(xintercept=1, colour = "black", linetype=2) + 
  geom_hline(yintercept=1, colour = "black", linetype=2) + 
  geom_text(aes(x=1, label="mRNA", y=3.5), colour=alpha("black",0.02), angle=90, vjust = -1, size = 6) +
  geom_text(aes(y=1, label="Protein", x= 5), colour=alpha("black",0.02), vjust = -1, size = 6) +
  geom_text(aes(x=1, label="threshold", y=3.5), colour=alpha("black",0.02), angle=90, vjust = 1, size = 6) +
  geom_text(aes(y=1, label="threshold", x=5), colour=alpha("black",0.02), vjust = 1, size = 6) + 
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text.x=element_text(size=15, hjust = 1),
        axis.text.y=element_text(size=15, hjust = 1),
        legend.key.height = unit(1, "cm"),
        legend.key.width= unit(1, "cm"),
        legend.key.size= unit(1, "cm")) 

density_x <- ggplot(data = pval_dat, aes(x = mrna, colour = sig, fill = sig)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  scale_y_continuous(name =  "Protein", expand = c(0.05, 0.05)) +
  scale_x_continuous(name = "mRNA", expand = c(0.05, 0.05)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) + 
  scale_fill_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
  guides(colour = FALSE, fill=FALSE) + 
  ggtitle("mRNA/Protein ANOVA adjusted P values") +
  theme(plot.margin      = unit(c(1, 0, 0, 1.4), "cm"),
        plot.title = element_text(size = 26, face = "bold", hjust=0.5),
        axis.ticks       = element_blank(), 
        panel.background = element_blank(), 
        panel.grid       = element_blank(),
        panel.border = element_blank(),
        axis.text.x      = element_blank(), 
        axis.text.y      = element_blank(),           
        axis.title.x     = element_blank(), 
        axis.title.y     = element_blank()) 


density_y <- ggplot(data = pval_dat, aes(x = prot, colour = sig, fill = sig)) +
  geom_density(alpha = .5) +
  theme_bw() + 
  scale_y_continuous(name =  "Protein", expand = c(0.05, 0.05)) +
  scale_x_continuous(name = "mRNA", expand = c(0.05, 0.05)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) + 
  scale_fill_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
  guides(colour = FALSE, fill=FALSE) + 
  theme(plot.margin      = unit(c(0, 0, 1.45, 0), "cm"),
        axis.ticks       = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(),
        panel.grid       = element_blank(),
        axis.text.x      = element_blank(), 
        axis.text.y      = element_blank(),           
        axis.title.x     = element_blank(), 
        axis.title.y     = element_blank()) + 
  coord_flip()

# empty <- ggplot() +
#   geom_point(aes(1,1), colour="white") +
#   theme(plot.margin      = unit(c(0, 0, 0, 0), "cm"),
#         axis.ticks       = element_blank(), 
#         panel.background = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.text.x      = element_blank(), 
#         axis.text.y      = element_blank(),           
#         axis.title.x     = element_blank(), 
#         axis.title.y     = element_blank())

legend_scatter <- get.legend(scatter)

scatter <- scatter + theme(legend.position="none")

# pdf("IntClust_Anova_Plot.pdf", width=20, height = 10)
grid.arrange(density_x, legend_scatter, scatter, density_y, 
                                ncol = 2, 
                                nrow = 2, 
                                widths  = c(4, 1), 
                                heights = c(1, 4))
# dev.off()

ur <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
  geom_point() + 
  theme_bw() + 
  scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(2,3.5)) +
  scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(4,5.5)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
  ggtitle("mRNA/Protein significant") + 
  geom_label_repel(aes(label = gene)) + 
  # guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) + 
  guides(colour=FALSE) + 
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text.x=element_text(size=7, hjust = 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width= unit(0.5, "cm")) 

# pdf("ur_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
plot(ur)
# dev.off() 

dr <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
  geom_point() + 
  theme_bw() + 
  scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(0,1)) +
  scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(4,5.5)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
  ggtitle("mRNA only significant") + 
  geom_label_repel(aes(label = gene)) + 
  # guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) + 
  guides(colour=FALSE) + 
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text.x=element_text(size=7, hjust = 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width= unit(0.5, "cm")) 

# pdf("dr_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
plot(dr)
# dev.off() 

ul <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
  geom_point() + 
  theme_bw() + 
  scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(2,3.5)) +
  scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(0,1)) + 
  scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +  ggtitle("Protein only significant") + 
  geom_label_repel(aes(label = gene)) + 
  guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) +
  # guides(colour=FALSE) + 
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text.x=element_text(size=7, hjust = 1),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width= unit(0.5, "cm")) 

# pdf("ul_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
plot(ul)
# dev.off() 

# pdf("IntCl_quartiles.pdf", width = 30, height = 10, paper = 'special')
grid.arrange(ur, dr, ul,
             ncol = 3, 
             nrow = 1)
# dev.off()

## Venn plots

total_genes <- nrow(mrna)
total_IntClGenes <- length(intCl754)
overlap_IntClGenes <- length(intersect(rownames(mrna),intCl754))

write.table(pval_dat, "venn.txt", quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = "\t")

# pdf("Total_study_IntClust_Venn.pdf", width=20, height=15, paper='special')
draw.pairwise.venn(
  area1 = total_genes,
  area2 = total_IntClGenes,
  cross.area = overlap_IntClGenes,
  fill = c("black", "blue"),
  label.col = 'white',
  category =  c("Total study genes", "IntClust genes"),
  cat.pos = c(-45, 45),
  euler.d = TRUE,
  cat.dist = c(0.03,0.03),
  # sep.dist = 0.03,
  offset = 1,
  cat.cex=2,
  ext.text = FALSE,
  cex=1.5,
  rotation.degree = 45,
  lwd = 2,
  print.mode = c('raw', 'percent')
)
# dev.off()

mrna_sig_IntCl <- sum(pval_dat$mrna > 1 & pval_dat$sig == 'yes')
prot_sig_IntCl <- sum(pval_dat$prot > 1 & pval_dat$sig == 'yes')
both_sig_IntCl <- sum(pval_dat$mrna > 1 & pval_dat$sig == 'yes' & pval_dat$prot > 1)

# pdf("mRNA_Protein_IntClust.pdf", width=20, height=15, paper='special')
draw.pairwise.venn(
  area1 = mrna_sig_IntCl,
  area2 = prot_sig_IntCl,
  cross.area = both_sig_IntCl,
  fill = c('palegreen', 'steelblue1'),
  category =  c("mRNA significant", " Protein significant"),
  cat.pos = c(-45, 45),
  euler.d = TRUE,
  cat.dist = c(0.03,0.03),
  # sep.dist = 0.03,
  offset = 1,
  cat.cex=2,
  cex = 1.5,
  ext.text = FALSE,
  rotation.degree = 45,
  lwd = 2,
  print.mode = c('raw', 'percent')
)
# dev.off()

## Boxplots for mRNA buffered genes
genes_both <- c("HEATR6","GGA3", "SUPT6H", "COQ9")
genes_mrna <- c("PSMB3","JUP","CDC16","SNF8")
genes_prot <- c("TLK2", "NUDT5", "IPO9","STAU2")

genesToBoxplot = genes_both
boxplot_data <- list('mRNA'=mrna_res,"Protein"=prot_res,"CNA"=cna_log)

# pdf(paste0("mRNA_Protein_Boxplots_mRNA_Protein_IntClust",".pdf"), width=25, height=10, paper="special")
par(mfrow=c(2,4), oma = c(4, 6, 4, 9))
for(i in c("mRNA","Protein")) {
  for (j in genesToBoxplot) {
    per.gene.boxplot(boxplot_data[[i]],cna_bin, j, i)
  }
}
mtext(paste0('Significant CNA-mRNA/Protein associations'), outer = TRUE, cex = 3)
mtext("Copy Number Status", side = 1, outer = TRUE, line = 2, cex=2)
mtext("    Protein ratio                       mRNA expression value", 
      side = 2, outer = TRUE, line = 2, cex=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend = molecular_subtype,
       title = "Subtype", cex = 2, pch = 16, col = PAM50cols, xpd=TRUE)
# dev.off()
#--------------------------------------------------------------------------------------------------------
## ANOVA ranked 

# test.results <- function(test, thres) {
#   mrna <- find.significant.genes(test,"anova",mrna, genes, "BH")
#   prot <- find.significant.genes(test,"anova",prot, genes, "BH")
#   combine <- merge(mrna, prot, by = "genes")
#   combine$mrna_sig <- combine$qval_anova.x<= thres
#   combine$prot_sig <- combine$qval_anova.y<= thres
#   combine$both_sig <- combine$mrna_sig & combine$prot_sig
#   return(combine)
# }
# res_anova_init <- test.results(do.anova.init,0.05) 
# 
# perc_mrna <- sum(res_anova_init$mrna_sig)/nrow(mrna)
# perc_prot <- sum(res_anova_init$prot_sig)/nrow(prot)
# perc_mrna_prot <- sum(res_anova_init$both_sig)/sum(res_anova_init$prot_sig)
# 
# testName = "anova_init"
# 
# ## Venn plot 1
# pdf(paste0("mRNA_Sig_Overlap","_",testName,".pdf"), width=20, height=15, paper='special')
# draw.pairwise.venn(
#   area1 = nrow(mrna),
#   area2 = sum(res_anova_init$mrna_sig),
#   cross.area =  sum(res_anova_init$mrna_sig),
#   fill = c('grey', 'green'),
#   category =  c("Total genes", "mRNA significant"),
#   cat.pos = c(0, 0),
#   euler.d = TRUE,
#   cat.dist = c(-0.03,-0.03),
#   # sep.dist = 0.03,
#   offset = 1,
#   cat.cex=2,
#   cex = 1.5,
#   ext.text = FALSE,
#   rotation.degree = 45,
#   lwd = 2,
#   print.mode = c('raw', 'percent')
# )
# dev.off()
# 
# 
# ## Venn plot 2
# pdf(paste0("mRNA_Sig_Prot_Sig_Overlap","_", testName,".pdf"), width=20, height=15, paper='special')
# draw.pairwise.venn(
#   area1 =   sum(res_anova_init$mrna_sig),
#   area2 = sum(res_anova_init$prot_sig),
#   cross.area =  sum(res_anova_init$both_sig),
#   fill = c('green', 'blue'),
#   category =  c("mRNA significant", "Protein significant"),
#   cat.pos = c(-45, 45),
#   euler.d = TRUE,
#   cat.dist = c(0.02,0.02),
#   # sep.dist = 0.03,
#   offset = 1,
#   cat.cex=2,
#   cex = 1.5,
#   ext.text = FALSE,
#   rotation.degree = 45,
#   lwd = 2,
#   print.mode = c('raw', 'percent')
# )
# 
# dev.off()

# ## Calculate pairwise comparisons significance level
# calculate.wilcox.pval.fc <- function(omicsType, binary_CNA, levelName, level) {
#   testTable <- data.frame(do.call(rbind,lapply(genes, difExprByCNstatus, omicsType = omicsType, binary_CNA = binary_CNA, level = level)))
#   testTable$pvalWilcoxonAdjusted = p.adjust(testTable[,1], method = "BH")
#   rownames(testTable) <- genes
#   colnames(testTable) = c(paste("pvalWilcoxon",levelName,sep="_"),paste("FoldChange",levelName,sep="_"),
#                           paste("Effect_Size", levelName, sep="_"),
#                           paste("adjustedPvalueWilcoxon",levelName,sep="_"))
#   return(testTable)
# }
#----------------------------------------------------------------------------------------------------------------------------------------------