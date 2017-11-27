### Correlation analysis

## Cis 
cis.cor <- function(gene, data) {
  cor = corAndPvalue(cna_log[gene, ], data[gene, ], method = "pearson")
  return(c(cor$cor, cor$p))
}

cis_mrna <- sapply(genes, cis.cor, data = mrna)
rownames(cis_mrna) <- c("cor", "pvalue")

cis_prot <- sapply(genes, cis.cor, data = prot)
rownames(cis_prot) <-  c("cor", "pvalue")

attenuate <- cis_mrna["cor", ] - cis_prot["cor",] 
# ----------------------------------------------------------------------------------------------------

# install.packages('e1071')
# library(e1071)

# gaussian_means <- mixture_gaussian$parameters$mean
# gaussian_var <- mixture_gaussian$parameters$variance$sigmasq
# 
# g1 <- data.frame(value = rnorm(n = 1000, mean = gaussian_means[1], sd = gaussian_var[1]), stat = "0")
# g2 <- data.frame(value = rnorm(n = 1000, mean = gaussian_means[2], sd = gaussian_var[2]), stat = "1")
# 
# g <- rbind(g1, g2)
# 
# svm_tune <- tune(svm, train.x=g$value, train.y=g$stat, 
#                  kernel="radial", ranges=list(cost=10^(-1:3), gamma=c(.5,1,2,3)))
# svm_model <- svm(stat ~., data=g, kernel="radial", cost=0.1, gamma=3)
# summary(svm_model)
# 
# table(pred,g$stat)


# 
# install.packages("mixtools")
# library(mixtools)

# library(stats) # for bw.SJ() and density()
# install.packages("modeest")
# library(modeest) # for mlv()
# # Assume data for sample i is in D
# b<-bw.SJ(attenuate) # Find bandwith with Shafer-Jones method
# d<-density(attenuate, bw=b, kernel="gaussian") # Gaussian kernel density estimate
# mode<- d$x[which.max(d$y)]
# model<-normalmixEM(attenuate,mu=c(mode,mode), k=2)
# sigma = min(model$sigma)

## Cluster according to density function
# install.packages("mclust")
# library(mclust)

# mixture_gaussian <- Mclust(attenuate)
# attenuate_group <- names(attenuate)[mixture_gaussian$classification==2]

# ----------------------------------------------------------------------------------------------------
## Classify according to attenuation level
attenuate_group <- sapply(attenuate, function(i) {
  if (i>0.2) {
    group = "High (>0.2)"
  } else {
    group = "Low"
    }
  if (i>0.3) {
    group = "High (>0.3)"
  }
  if (i>0.4) {
    group = "High (>0.4)"
  }
  if (i>0.5) {
    group = "High (>0.5)"
  } 
  return(group)
})


## Create data frame
attenuate <- data.frame(gene = I(names(attenuate)), attenuation = attenuate, group = attenuate_group)

# Load ubiquitination data
ubiq_table <- read.xlsx2("ubiquitin_matrix.xls", sheetIndex = 2, colClasses = "character", stringsAsFactors = FALSE)
ubiq_table <- ubiq_table[, !grepl("X", colnames(ubiq_table))]
# str(ubiq_table)
cols <- c("gene_symbol", "bortezomib_0hrs_median", "bortezomib_2hrs_median", "bortezomib_4hrs_median", "bortezomib_8hrs_median")
ubiq_table <- ubiq_table[, cols]
ubiq_val <- apply(ubiq_table[,-1], 2, as.numeric)
ubiq_table <- data.frame(gene = I(ubiq_table[, "gene_symbol"]), ubiq_val)

# Collapse to gene symbols
ubiq_table <- as.data.frame(ubiq_table %>% group_by(gene) %>% summarise_each(funs(median(.,na.rm=TRUE))))

# Remove "-" gene symbol
ubiq_table <- ubiq_table[-1,]

# Correct gene symbols mistakes
ubiq_table$gene[1:6] <- c("MARCH6", "SEPT2", "SEPT7", "SEPT9", "SEPT10", "SEPT11")
rownames(ubiq_table) <- NULL
# toReplace <- c("40243", "40423", "40428", "40430", "40431", "40432")
# toAdd <- c("MARCH2","SEPT2","SEPT7","SEPT9", "SEPT10", "SEPT11")
# 
# for (i in 1:length(toReplace))  {
#   idx = which(ubiq_table$gene == toReplace[i])
#   ubiq_table$gene[idx] = toAdd[i]
# }
# 

alias <- sapply(strsplit(ubiq_table$gene,";"), length)

new_genes <- unlist(strsplit(ubiq_table$gene,";"))

ubiq_table <- data.frame(apply(ubiq_table, 2,  function(i) {
   rep(i,alias)}), stringsAsFactors = FALSE)

ubiq_table <- data.frame(gene = new_genes, ubiq_table[,-1], stringsAsFactors = FALSE)


# ubiq_table$gene <- new_genes

# Merge attenuate with ubiquitination data
merge_attenuation_ubiq <- merge(attenuate, ubiq_table, by = "gene")
merge_attenuation_ubiq$group <- factor(merge_attenuation_ubiq$group, levels = c("Low", "High (>0.2)", "High (>0.3)", 
                                                                                "High (>0.4)", "High (>0.5)"))
merge_attenuation_ubiq <- melt(merge_attenuation_ubiq, id.vars = c("gene", "attenuation", "group"))
merge_attenuation_ubiq$value <- as.numeric(as.character(merge_attenuation_ubiq$value)) 

colfunct <- colorRampPalette(c("grey", "red"))
colours <- colfunct(5)

# Remove NAs
merge_attenuation_ubiq <- merge_attenuation_ubiq[!is.na(merge_attenuation_ubiq$value), ]

# Find mean, sum and label position for each group 

boxplot.top <- function(value) {
  thres <- quantile(value,0.75) + 1.5*IQR(value)
  idx1 <- which(value - thres <= 0)
  val <- value[idx1]
  idx2 <- which.min(abs(val - thres))
  new_val <- val[idx2]
  return(new_val)
  }

per_col_sum <- as.data.frame(merge_attenuation_ubiq %>% 
  dplyr::group_by(variable, group) %>%
  dplyr::summarise(n = n(), m = mean(value,na.rm = TRUE),
                   q3 = boxplot.top(value)))
  

# Kruskal-Wallis P-value
pval_KW <- sapply(levels(merge_attenuation_ubiq$variable) , function(i) {
  k <- kruskal.test(value ~ group, merge_attenuation_ubiq[merge_attenuation_ubiq$variable==i,]) 
  pval <- k$p.value
  })

pval_Wilc <- sapply(levels(merge_attenuation_ubiq$variable), function(i) {
  w <- wilcox.test(merge_attenuation_ubiq$value[merge_attenuation_ubiq$variable==i & 
                                                  merge_attenuation_ubiq$group=="High (>0.5)"],
            merge_attenuation_ubiq$value[merge_attenuation_ubiq$variable==i & 
                                           merge_attenuation_ubiq$group=="Low"])
  pval <- w$p.value
  })

pval_dat <- data.frame(pval_KW, pval_Wilc)

p <- ggplot(data = merge_attenuation_ubiq, aes(x = variable, y = value, fill = group)) + 
  stat_boxplot(geom ='errorbar', linetype = 2, position = position_dodge(width=0.9)) + 
  geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = TRUE,
               position = position_dodge(width=0.9)) +
  geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2)) +
  scale_fill_manual(values = colours) +
  scale_x_discrete(name = "Hours", labels = c("0","2","4","8")) +
  scale_y_continuous(name = "Ubiquitination sites (log2 FC)", breaks = seq(-10,10,2), 
                     labels = seq(-10,10,2)) + 
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("n=",n)), size = 3, 
           position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("mean=",round(m,2))), size = 3, 
            position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  annotate(geom = "text", x = 0.7, y = -5.5, label = "Kruskal-Wallis \nP-value", size = 3, fontface =2) +
  annotate(geom = "text", x = 0.7, y = -7, label = "Wilcoxon \nP-value \n(Low vs High (>0.5)", size = 3, fontface =2) +
  annotate(geom = "text", x = 1:4, y = -6, label = paste0("P=",round(pval_dat$pval_KW,5)), size = 5) +
  annotate(geom = "text", x = 1:4, y = -7, label = paste0("P=",round(pval_dat$pval_Wilc,5)), size = 5) +
  theme_bw() + 
  ggtitle("Proteasome Inhibition (Bortezomid)") + 
  guides(fill = guide_legend(title = "Protein attenuation", override.aes = list(fill = colours, alpha = 1,
                                                                                size = 10))) +
  
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text=element_text(size=10, hjust = 1))

pdf("Protein_ubiquitination_residual.pdf" , width = 20, height = 7, paper = "special")
plot(p)
dev.off()
# -----------------------------------------------------------------------------------------------------------------------------
### Protein stability evaluation
stability_table <- read.xlsx2("protein_stability.xls", sheetIndex = 1, colClasses = "character", stringsAsFactors = FALSE)
# str(stability_table)
new_genes <- unlist(strsplit(stability_table$Gene.Names, ";"))

alias <- sapply(strsplit(as.character(stability_table$Gene.Names),";"), length)
stability_table <- apply(stability_table, 2,  function(i) {
  rep(i,alias)})
stability_table <- data.frame(MGI.symbol = new_genes, stability_table, stringsAsFactors = FALSE)
# Convert human to mouse gene names

# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
# mouse <-  useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = genes , mart = ensembl_75, attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)

matched_genes <- lapply(1:nrow(genesV2), function(i) {
  idx <- match(genesV2$HGNC.symbol[i],genes)
  attenuate[idx,]
})

matched_genes <- do.call(rbind, matched_genes)
attenuate_mouse <- cbind(genesV2, matched_genes)
attenuate_mouse <- attenuate_mouse[,-1]

merged_stability_table <- merge(attenuate_mouse, stability_table, by = "MGI.symbol")


thres_mrna_stable = quantile(as.numeric(merged_stability_table$mRNA.half.life.average..h.), 2/3 ,na.rm = TRUE)
thres_mrna_unstable = quantile(as.numeric(merged_stability_table$mRNA.half.life.average..h.), 1/3 ,na.rm = TRUE)

thres_prot_stable = quantile(as.numeric(merged_stability_table$Protein.half.life.average..h.), 2/3, na.rm = TRUE)
thres_prot_unstable = quantile(as.numeric(merged_stability_table$Protein.half.life.average..h.), 1/3, na.rm = TRUE)

merged_stability_table$stat_mrna <- ifelse(as.numeric(merged_stability_table$mRNA.half.life.average..h.)>thres_mrna_stable, 
                                           "mrna_stable", ifelse(as.numeric(merged_stability_table$mRNA.half.life.average..h.)<thres_mrna_unstable,
                                           "mrna_unstable","NA"))

merged_stability_table$stat_protein <- ifelse(as.numeric(merged_stability_table$Protein.half.life.average..h.)>thres_prot_stable, 
                                             "prot_stable", ifelse(as.numeric(merged_stability_table$Protein.half.life.average..h.)<thres_prot_unstable,
                                                                 "prot_unstable","NA"))

merged_stability_table$stat_both <- factor(paste0(merged_stability_table$stat_mrna,"-",merged_stability_table$stat_protein), 
                                           levels = c("mrna_stable-prot_stable", "mrna_stable-prot_unstable", "mrna_unstable-prot_stable",
                                                      "mrna_unstable-prot_unstable"))
merged_stability_table_intCl <- merged_stability_table
merged_stability_table <- merged_stability_table[!is.na(merged_stability_table$stat_both), ]
merge_attenuation_ubiq$group <- factor(merge_attenuation_ubiq$group, levels = c("Low", "High (>0.2)", "High (>0.3)", 
                                                                                "High (>0.4)", "High (>0.5)"))

per_col_sum <- as.data.frame(merged_stability_table %>% 
                               dplyr::group_by(stat_both) %>%
                               dplyr::summarise(n = n(), m = mean(attenuation,na.rm = TRUE),
                                                q3 = boxplot.top(attenuation)))


# Kruskal-Wallis P-value
k <- kruskal.test(attenuation ~ stat_both, merged_stability_table) 
pval_KW <- k$p.value


w <- wilcox.test(merged_stability_table$attenuation[merged_stability_table$stat_both=="mrna_unstable-prot_unstable"],
                   merged_stability_table$attenuation[merged_stability_table$stat_both=="mrna_stable-prot_stable"])
pval_Wilc <- w$p.value


install.packages("wesanderson")
library(wesanderson)

p <- ggplot(merged_stability_table, aes(x = stat_both, y = attenuation, fill = stat_both)) + 
  stat_boxplot(geom ='errorbar', linetype = 2) + 
  geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = TRUE) + 
  geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge()) +
  geom_text(data = per_col_sum, aes(x = stat_both, y = q3, label = paste0("n=",n)), size = 3, 
            position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = stat_both, y = q3, label = paste0("mean=",round(m,2))), size = 3, 
            position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  scale_fill_manual(values = wes_palette("GrandBudapest"),  label = c("stable mRNA\nstable Protein","stable mRNA\nunstable Protein",
                                                                      "unstable mRNA\nstable Protein","unstable mRNA\nunstable Protein")) +
  scale_x_discrete(name = "", labels = c("stable mRNA/stable Protein \nI","stable mRNA/unstable Protein \nII",
                                         "unstable mRNA/stable Protein \nIII","unstable mRNA/unstable Protein \nIV")) +
  scale_y_continuous(name = "Attenuation") + 
  theme_bw() + 
  ggtitle("mRNA/Protein stability") + 
  guides(fill = FALSE) +
  annotate(geom = "text", x = 4, y = -0.5, label = paste0("Kruskal Wallis P-value=",round(pval_KW,6)), size = 3) +
  annotate(geom = "text", x = 4, y = -0.6, label = paste0("Wilcoxon P-value (I vs IV)=",round(pval_Wilc,5)), size = 3) +
  theme(axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        plot.title = element_text(size = 28, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=25) , legend.text=element_text(size=20),
        axis.text=element_text(size=10, hjust = 0.5))

pdf("Protein_stability_initial.pdf" , width = 10, height = 7, paper = "special")
plot(p)
dev.off()
# --------------------------------------------------------------------------------------------------
# ## Apply ubiquitination/stability data to IntClust genes
# 
# intCl_both_significant <- as.character(pval_dat$gene[pval_dat$mrna>1 & pval_dat$prot>1])
# intCl_mrna_significant <- as.character(pval_dat$gene[pval_dat$mrna>1 & pval_dat$prot<1])
# intCl_prot_significant <- as.character(pval_dat$gene[pval_dat$mrna<1 & pval_dat$prot>1])
# 
# 
# # intersect(ubiq_table[,1], intCl_both_significant)
# # intersect(ubiq_table[,1], intCl_mrna_significant)
# # intersect(ubiq_table[,1], intCl_prot_significant)
# 
# ubiq_table$intCl <- NA
# ubiq_table$intCl[ubiq_table$gene %in% intCl_both_significant] <- 'both'
# ubiq_table$intCl[ubiq_table$gene %in% intCl_mrna_significant] <- 'mrna'
# ubiq_table$intCl[ubiq_table$gene %in% intCl_prot_significant] <- 'protein'
# 
# ubiq_table$intCl <- factor(ubiq_table$intCl, levels = c("mrna","protein","both"))
# ubiq_table <- melt(ubiq_table, id.vars =  c("gene", "intCl"))
# ubiq_table$value <- as.numeric(as.character(ubiq_table$value)) 
# 
# 
# str(ubiq_table)
# 
# 
# ubiq_table <- ubiq_table[!(is.na(ubiq_table$intCl) | is.na(ubiq_table$value)),]
# 
# per_col_sum <- as.data.frame(ubiq_table %>% 
#                                dplyr::group_by(variable,intCl) %>%
#                                dplyr::summarise(n = n(), m = mean(value,na.rm = TRUE),
#                                                 q3 = boxplot.top(value)))
# 
# # Kruskal-Wallis P-value
# pval_KW <- sapply(levels(ubiq_table$variable) , function(i) {
#   k <- kruskal.test(value ~ intCl, ubiq_table[ubiq_table$variable==i,]) 
#   pval <- k$p.value
# })
# 
# pval_Wilc <- sapply(levels(ubiq_table$variable), function(i) {
#   w <- wilcox.test(ubiq_table$value[ubiq_table$variable==i & 
#                                                   ubiq_table$intCl=="mrna"],
#                   ubiq_table$value[ubiq_table$variable==i & 
#                                      ubiq_table$intCl=="both"])
#   pval <- w$p.value
# })
# 
# pval_dat <- data.frame(pval_KW, pval_Wilc)
# 
# 
# ggplot(ubiq_table, aes(x = variable, y = value, fill = intCl)) + 
#   stat_boxplot(geom ='errorbar', linetype = 2) + 
#   geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = FALSE) + 
#   geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge()) +
#   scale_fill_manual(labels = c("mRNA","Protein","Both"), values = c("mrna"="green","protein"="blue", "both" = "black")) + 
#   scale_x_discrete(name = "Hours", labels = c("0","2","4","8")) +
#   scale_y_continuous(name = "Ubiquitination sites (log2 FC)", breaks = seq(-6,6,2), 
#                      labels = seq(-6,6,2)) + 
#   geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("n=",n)), size = 3, 
#             position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
#   geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("mean=",round(m,2))), size = 3, 
#             position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
#   annotate(geom = "text", x = 1:4, y = -6, label = paste0("P=",round(pval_dat$pval_KW,5)), size = 5) +
#   annotate(geom = "text", x = 1:4, y = -7, label = paste0("P=",round(pval_dat$pval_Wilc,5)), size = 5) +
#   theme_bw() + 
#   ggtitle("Proteasome Inhibition (Bortezomid)") + 
#   guides(fill = guide_legend(title = "Level of \nsignificance", 
#                              override.aes = list(fill = c("green","blue","black"), alpha = 1,
#                                                                                 size = 10))) +
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         panel.grid.major.x  = element_blank(),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text=element_text(size=10, hjust = 1))
# 
# # Stability
# 
# merged_stability_table_intCl$intCl <- NA
# merged_stability_table_intCl$intCl[merged_stability_table_intCl$gene %in% intCl_both_significant] <- 'both'
# merged_stability_table_intCl$intCl[merged_stability_table_intCl$gene %in% intCl_mrna_significant] <- 'mrna'
# merged_stability_table_intCl$intCl[merged_stability_table_intCl$gene %in% intCl_prot_significant] <- 'protein'
# merged_stability_table_intCl$intCl <- factor(merged_stability_table_intCl$intCl, levels = c("mrna","protein","both"))
# 
# merged_stability_table_intCl <- merged_stability_table_intCl[!is.na(merged_stability_table_intCl$intCl), ]
# 
# merged_stability_table_intCl$stab <- as.factor(ifelse(merged_stability_table_intCl$intCl=="both" & merged_stability_table_intCl$stat_mrna=='mrna_stable' & merged_stability_table_intCl$stat_protein=='prot_stable',"stable",
#                                       ifelse(merged_stability_table_intCl$intCl=="mrna" & merged_stability_table_intCl$stat_mrna=='mrna_stable', "stable",
#                                              ifelse(merged_stability_table_intCl$intCl=="protein" & merged_stability_table_intCl$stat_protein=='prot_stable', "stable",
#                                              'unstable'))))
# 
# merged_stability_table_intCl$Protein.half.life.average..h. <- as.numeric(merged_stability_table_intCl$Protein.half.life.average..h.)
# merged_stability_table_intCl$mRNA.half.life.average..h. <- as.numeric(merged_stability_table_intCl$mRNA.half.life.average..h.)
# 
# 
# ggplot(merged_stability_table_intCl, aes(x = intCl, y = as.numeric(translation.rate.constant..ksp..average..molecules..mRNA.h..), fill = intCl), na.rm = TRUE) + 
#   stat_boxplot(geom ='errorbar', linetype = 2) + 
#   geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = TRUE) + 
#   geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge(jitter.width = 0.1)) + 
#   scale_fill_manual(labels = c("stable","unstable"), values = wes_palette("GrandBudapest")) + 
#   scale_x_discrete(name = "", labels = c("mRNA","Protein",
#                                          "Both")) +
#   scale_y_continuous(name = "Attenuation", limits = c(0,50)) + 
#   theme_bw() + 
#   ggtitle("mRNA/Protein stability") + 
#   guides(fill = guide_legend(title = "Stability", 
#                              override.aes = list(fill = c("green","blue"), alpha = 1,
#                                                  size = 10))) +
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         panel.grid.major.x  = element_blank(),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text=element_text(size=10, hjust = 0.5))

# --------------------------------------------------------------------------------------------------
## Trans correlations
# a <- corAndPvalue(t(cna_log), t(mrna_res))
# 
# pval <- a$p
# pval.adj <- p.adjust(pval, method="BH")
# 
# b <- matrix(pval.adj, ncol = ncol(a$cor))
# dimnames(b) <- list(genes,genes)
# 
# ma <- melt(a$cor)
# mb <- melt(b)
# mc <- ma[mb$value<0.01, ]
# 
# aug_mc <- augment_biomart(mc, values = mc$X1)
# aug_mc2 <- augment_biomart(mc, values = mc$X2)
# 
# mc3 <- data.frame(cna = aug_mc$hgnc_symbol, pos_cna = aug_mc$pos, gene = aug_mc2$hgnc_symbol, pos_gene = aug_mc2$pos)
# ggplot(mc3, aes(x = pos_cna, y = pos_gene)) + 
#   geom_point() + 
#   scale_x_continuous(name="Chromosome (cna)", breaks= (chrlengths$start +chrlengths$end)/2,
#                      labels=chrlengths$chr, limits=c(-1e8,chrlengths$end[23]), expand = c(0,0)) +
#   scale_y_continuous(name= ("Chromosome (rna)"), breaks= (chrlengths$start +chrlengths$end)/2,
#                      labels=chrlengths$chr, limits=c(-1e8,chrlengths$end[23]), expand = c(0,0)) + 
#   geom_vline(xintercept = c(chrlengths$start[1], chrlengths$end),
#              linetype = 1,
#              color= alpha("black", 0.05)) + 
#   geom_hline(yintercept = c(chrlengths$start[1], chrlengths$end),
#              linetype = 1,
#              color= alpha("black", 0.05)) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(), axis.ticks.x = element_line(colour = "black"),axis.ticks.y = element_line(colour = "black"), axis.ticks.length=unit(0.1,"cm")) 
#   


