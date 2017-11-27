### Correlation analysis

# CNA-mRNA
cis_cna_mrna <- sapply(genes, cis.cor, data1 = cna_log, data2 = mrna)

# CNA - Protein
cis_cna_prot <- sapply(genes, cis.cor, data1 = cna_log, data2 = prot)

# 1st metric - mRNA - Protein correlation
cis_mrna_prot <- t(sapply(genes, cis.cor, data1 = mrna, data2 = prot))
cis_mrna_prot <- as.data.frame(cis_mrna_prot)
cis_mrna_prot$gene <- rownames(cis_mrna_prot)

# # 2nd  metric - Attenuation level
# attenuate <- cis_cna_mrna["cor", ] - cis_cna_prot["cor",]
# attenuate <- data.frame(gene = names(attenuate), val = attenuate)
# 
# # 3rd anti-metric - Semi-partial correlation
# install.packages('ppcor')
# library(ppcor)

# attenuate <- t(mclapply(1:nrow(cna_log), function(i) spcor.test(as.numeric(cna_log[i, ]), as.numeric(prot[i, ]), as.numeric(mrna[i, ]))))
# attenuate <- as.data.frame(apply(attenuate, 2, unlist))
# rownames(attenuate) <- genes
# attenuate$adj.pvalue <- p.adjust(attenuate$p.value)
# genes[which(attenuate$adj.pvalue < 0.1)]

# # Plot an example
# 
# gene <- 'ASH1L'
# 
# per.gene.scatterplot(data1 = cna_log, 
#                      data2 = mrna, 
#                      gene = gene,
#                      molecular_view = 'mRNA',
#                      cor_method = 'spearman')
# 
# per.gene.scatterplot(data1 = cna_log, 
#                      data2 = prot, 
#                      gene = gene,
#                      molecular_view = 'Protein',
#                      cor_method = 'spearman')


# Load ubiquitination data
ubiq_table <- read.xlsx2(paste0(filefolder, "ubiquitin_matrix.xls"), sheetIndex = 2, colClasses = "character", stringsAsFactors = FALSE)
ubiq_table <- ubiq_table[, !grepl("X", colnames(ubiq_table))]
# str(ubiq_table)
cols <- c("gene_symbol", "bortezomib_0hrs_median", "bortezomib_2hrs_median", "bortezomib_4hrs_median", "bortezomib_8hrs_median")
ubiq_table <- ubiq_table[, cols]
ubiq_val <- suppressWarnings(apply(ubiq_table[,-1], 2, as.numeric))
ubiq_table <- data.frame(gene = I(ubiq_table[, "gene_symbol"]), ubiq_val)

# Collapse to gene symbols
ubiq_table <- as.data.frame(ubiq_table %>% group_by(gene) %>% summarise_all(funs(median(.,na.rm=TRUE))))

# Remove "-" gene symbol
# ubiq_table <- ubiq_table[-1,]

# Correct gene symbols mistakes
# ubiq_table$gene[1:6] <- c("MARCH6", "SEPT2", "SEPT7", "SEPT9", "SEPT10", "SEPT11")
# rownames(ubiq_table) <- NULL
toReplace <- c("40243", "40423", "40428", "40430", "40431", "40432")
toAdd <- c("MARCH2","SEPT2","SEPT7","SEPT9", "SEPT10", "SEPT11")

for (i in 1:length(toReplace))  {
  idx = which(ubiq_table$gene == toReplace[i])
  ubiq_table$gene[idx] = toAdd[i]
}

# Alias genes
alias <- sapply(strsplit(ubiq_table$gene,";"), length)
new_genes <- unlist(strsplit(ubiq_table$gene,";"))
ubiq_table <- apply(ubiq_table, 2,  function(i) {
   rep(i,alias)})

ubiq_table <- as.data.frame(ubiq_table, stringsAsFactors = FALSE)

# Remove proteins with unknown gene symbol
ubiq_table <- ubiq_table[!grepl('^-', ubiq_table$gene), ] 
ubiq_table[,-1] <- apply(ubiq_table[,-1], 2, as.numeric)

merge_attenuation_ubiq <- merge(cis_mrna_prot, ubiq_table, by = "gene")

groups <- sapply(colnames(merge_attenuation_ubiq)[-c(1:4)], function(i) {
  new_col <- ifelse(abs(merge_attenuation_ubiq[,i]) > 1, 'High', 'Low')
  return(new_col)
  }
  )

colnames(groups) <- paste0('group_', substr(colnames(groups), 12,15))

merge_attenuation_ubiq <- cbind.data.frame(merge_attenuation_ubiq[, c("gene", "cor")], groups)                                                                                
merge_attenuation_ubiq <- melt(merge_attenuation_ubiq, id.vars = c("gene", "cor"))
# merge_attenuation_ubiq$value <- as.numeric(as.character(merge_attenuation_ubiq$value)) 

colfunct <- colorRampPalette(c("grey", "red"))
colours <- colfunct(2)


# Find mean, sum and label position for each group 
boxplot.top <- function(value) {
  thres <- quantile(value,0.75) + 1.5*IQR(value)
  idx1 <- which(value - thres <= 0)
  val <- value[idx1]
  idx2 <- which.min(abs(val - thres))
  new_val <- val[idx2]
  return(new_val)
}

  
per_col_sum <- as.data.frame(na.omit(merge_attenuation_ubiq) %>% 
                               dplyr::group_by(variable, value) %>%
                               dplyr::summarise(n = n(), m = median(cor, na.rm = TRUE),
                                                q3 = boxplot.top(cor)))


pval_Wilc <- sapply(levels(merge_attenuation_ubiq$variable), function(i) {
  w <- wilcox.test(merge_attenuation_ubiq$cor[merge_attenuation_ubiq$variable == i &
                                                  merge_attenuation_ubiq$value == "High"],
            merge_attenuation_ubiq$cor[merge_attenuation_ubiq$variable == i &
                                           merge_attenuation_ubiq$value == "Low"])
  pval <- w$p.value
  })


merge_attenuation_ubiq$value <- factor(merge_attenuation_ubiq$value, levels = c('High', 'Low'), ordered = TRUE)

p <-  ggplot(data = na.omit(merge_attenuation_ubiq), aes(x = variable, y = cor, fill = value)) + 
  stat_boxplot(geom ='errorbar', linetype = 2, position = position_dodge(width=0.9)) +
  geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = TRUE,
               position = position_dodge(width=0.9)) +
  geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2),
              size = 4) +
  scale_fill_manual(values = c('black', 'grey'), labels = c("High (|log2FC| > 1)", "Low (|log2FC| < 1)")) +
  scale_x_discrete(name = "Time (hours)", labels = c(2,4,8)) +
  scale_y_continuous(name = "mRNA - Protein Spearman correlation", breaks = seq(-0.4, 1, 0.2)) + 
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("n=",n)), size = 10, 
            position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("median=",round(m,2))), size = 10, 
            position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  geom_signif(y_position=rep(1.1, 3), xmin=c(0.8, 1.8, 2.8), xmax=c(1.2, 2.2, 3.2),
              annotation= paste0("P=",round(pval_Wilc,5)), tip_length = 0.005, textsize = 10) + 
  # annotate(geom = "text", x = 0.7, y = -5.5, label = "Kruskal-Wallis \nP-value", size = 3, fontface =2) +
  # annotate(geom = "text", x = 0.7, y = -7, label = "Wilcoxon \nP-value \n(Low vs High (>0.5)", size = 3, fontface =2) +
  # annotate(geom = "text", x = 1:3, y = -6, label = paste0("P=",round(pval_Wilc,5)), size = 5) +
  # annotate(geom = "text", x = 1:4, y = -7, label = paste0("P=",round(pval_dat$pval_Wilc,5)), size = 5) +
  theme_bw() + 
  ggtitle("Proteasome Inhibition (Bortezomib)") + 
  guides(fill = guide_legend(title = "Differential \nUbiquitination", override.aes = list(alpha = 1,
                                                                                size = 10))) +
  
  theme(axis.title.x = element_text(size=35),
        axis.title.y = element_text(size=35),
        plot.title = element_text(size = 40, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=30) , legend.text=element_text(size=28),
        legend.position = "bottom",
        axis.text=element_text(size=20))

# pdf("Protein_ubiquitination_Spearman.pdf" , width = 20, height = 15, paper = "special")
plot(p)
# dev.off()
# -----------------------------------------------------------------------------------------------------------------------------