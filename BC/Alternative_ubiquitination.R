## Attenution plot

# CNA - mRNA correlation
cis_mrna <- sapply(genes, cis.cor, data1 = cna_log, data2 = mrna)
rownames(cis_mrna) <- c("cor", "pvalue")

# CNA - Protein correlation
cis_prot <- sapply(genes, cis.cor,data1 = cna_log, data2 = prot)
rownames(cis_prot) <-  c("cor", "pvalue")

# Attenuation
attenuate <- cis_mrna["cor", ] - cis_prot["cor",] 

# Cluster according to density function

# Model attenuation as as mixture of Gaussians
set.seed(123)
mixture_gaussian <- Mclust(attenuate, G = 2)
group = as.factor(mixture_gaussian$classification)

# Set negative attenuation to 0
group[attenuate < 0 ] <- 1

# Plot
figureData <- data.frame('mRNA' = cis_mrna[1,], 'Protein' = cis_prot[1,], 'group' = group, 'gene' = names(attenuate))

p <- ggplot(data = figureData, mapping = aes(x = mRNA, y= Protein, colour = group)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_x_continuous(name = 'Copy-number ~ Transcriptomics \n(Spearman)',limits = c(-0.6, 1)) + 
  scale_y_continuous(name = 'Copy-number ~ Proteomics \n(Spearman)', limits = c(-0.6, 1)) +
  scale_colour_manual(name = "Protein \nattenuation", values = c('1' = alpha('grey', 0.3), '2' = alpha('orangered1', 0.3)), breaks = c(2,1), labels = c('High', 'Low')) + 
  geom_density2d(aes(group = group), colour = 'white', size = 0.5) + 
  geom_vline(xintercept = 0, colour = alpha('grey', 0.3)) + 
  geom_hline(yintercept = 0, colour = alpha('grey', 0.3))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 15)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
  

density_x <- ggplot(data = figureData, mapping = aes(x = mRNA, colour = group, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_hline(yintercept = 0, colour = 'white') +
  theme_bw() +
  scale_y_continuous(name =  "Protein") +
  scale_x_continuous(name = "mRNA") +
  scale_colour_manual(name = "Protein \nattenuation", 
                      values = c('1' = alpha('grey', 1), '2' = alpha('orangered1', 1)), 
                      labels = c('Low', 'High')) + 
  scale_fill_manual(name = "Protein \nattenuation", 
                    values = c('1' = alpha('grey', 0.3), '2' = alpha('orangered1', 0.3)), 
                    labels = c('Low', 'High')) +
  guides(colour = FALSE, fill=FALSE) +
  theme(plot.margin      = unit(c(0, 0, 0, 2.2), "cm"),
        plot.title = element_text(size = 26, face = "bold", hjust=0.5),
        axis.ticks       = element_blank(),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        panel.border = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.title.x     = element_blank(),
        axis.title.y     = element_blank())


density_y <- ggplot(data = figureData, mapping = aes(x = Protein, colour = group, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_hline(yintercept = 0, colour = 'white') +
  theme_bw() +
  scale_y_continuous(name =  "Protein") +
  scale_x_continuous(name = "mRNA") +
  scale_colour_manual(name = "Protein \nattenuation", 
                      values = c('1' = alpha('grey', 1), '2' = alpha('orangered1', 1)), 
                      labels = c('Low', 'High')) + 
  scale_fill_manual(name = "Protein \nattenuation", 
                    values = c('1' = alpha('grey', 0.3), '2' = alpha('orangered1', 0.3)), 
                    labels = c('Low', 'High')) +
  guides(colour = FALSE, fill=FALSE) +
  theme(plot.margin      = unit(c(0.3, 0, 1.8, 0), "cm"),
        plot.title = element_text(size = 26, face = "bold"),
        axis.ticks       = element_blank(),
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        panel.border = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.title.x     = element_blank(),
        axis.title.y     = element_blank()) + 
  coord_flip()



legend_p <- get.legend(p)

p <- p + theme(legend.position="none")

# pdf("Spearman_correlations_Attenuation.pdf", width=10, height = 10)
grid.arrange(density_x, legend_p, p, density_y,
                                ncol = 2,
                                nrow = 2,
                                widths  = c(4, 1),
                                heights = c(1, 4))
# dev.off()

## Export table
# attenuationGroupTable <- figureData[ , c('gene', 'group')]
# attenuationGroupTable$group <- revalue(attenuationGroupTable$group, replace = c('1' = 'Low', '2' = 'High'))
# colnames(attenuationGroupTable) <- c('Gene symbol', 'Attenuation group')
# write.table(x = attenuationGroupTable, file = 'Protein_attenuation_group.csv', row.names = FALSE, sep = ',')


# Load ubiquitination data
ubiq_table <- read.xlsx2(paste0(filefolder, "ubiquitin_matrix.xls"), sheetIndex = 2, colClasses = "character", stringsAsFactors = FALSE)
ubiq_table <- ubiq_table[, !grepl("X", colnames(ubiq_table))]
cols <- c("gene_symbol", "bortezomib_0hrs_median", "bortezomib_2hrs_median", "bortezomib_4hrs_median", "bortezomib_8hrs_median")
ubiq_table <- ubiq_table[, cols]
ubiq_val <- suppressWarnings(apply(ubiq_table[,-1], 2, as.numeric))
ubiq_table <- data.frame(gene = I(ubiq_table[, "gene_symbol"]), ubiq_val)

# Collapse to gene symbols
ubiq_table <- as.data.frame(ubiq_table %>% group_by(gene) %>% summarise_all(funs(median(.,na.rm=TRUE))))

# Remove "-" gene symbol
ubiq_table <- ubiq_table[-1,]

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

merge_attenuation_ubiq <- merge(figureData, ubiq_table, by = "gene")

merge_attenuation_ubiq <- melt(merge_attenuation_ubiq, id.vars = c("gene", "mRNA", "Protein", 'group'))
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

merge_attenuation_ubiq$group <- revalue(merge_attenuation_ubiq$group, replace = c('1' = 'Low', '2' = 'High'))

per_col_sum <- as.data.frame(na.omit(merge_attenuation_ubiq) %>% 
                               dplyr::group_by(variable, group) %>%
                               dplyr::summarise(n = n(), m = median(value, na.rm = TRUE),
                                                q3 = boxplot.top(value)))


pval_Wilc <- sapply(levels(merge_attenuation_ubiq$variable), function(i) {
  w <- wilcox.test(merge_attenuation_ubiq$value[merge_attenuation_ubiq$variable == i &
                                                merge_attenuation_ubiq$group == "High"],
                   merge_attenuation_ubiq$value[merge_attenuation_ubiq$variable == i &
                                                merge_attenuation_ubiq$group == "Low"])
  pval <- w$p.value
})



p <-  ggplot(data = na.omit(merge_attenuation_ubiq), aes(x = variable, y = value, fill = group)) + 
  stat_boxplot(geom ='errorbar', linetype = 2, position = position_dodge(width=0.9)) +
  geom_boxplot(outlier.colour = NA, outlier.size = NA, coef = 0, show.legend = FALSE, notch = TRUE,
               position = position_dodge(width=0.9)) +
  geom_jitter(pch=21, colour = "white", alpha = 0.3, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2),
              size = 4) +
  scale_fill_manual(name = 'Protein \nattenuation', values = c('grey', 'orangered1')) +
  scale_x_discrete(name = "Time (hours)", labels = c(0,2,4,8)) +
  scale_y_continuous(name = "Ubiquitination (log2FC)", breaks = seq(-6, 10, 2)) + 
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("n=",n)), size = 7, 
            position=position_dodge(width=0.9), vjust=-1, hjust = 0.5) +
  geom_text(data = per_col_sum, aes(x = variable, y = q3, label = paste0("median=",round(m,2))), size = 7, 
            position=position_dodge(width=0.9),  vjust=-3, hjust = 0.5) + 
  geom_signif(y_position=rep(7, 4), xmin=c(0.75, 1.75, 2.75, 3.75), xmax=c(1.25, 2.25, 3.25, 4.25),
              annotation= paste0("P=",round(pval_Wilc ,5)), tip_length = 0.005, textsize = 10) + 
  # annotate(geom = "text", x = 0.7, y = -5.5, label = "Kruskal-Wallis \nP-value", size = 3, fontface =2) +
  # annotate(geom = "text", x = 0.7, y = -7, label = "Wilcoxon \nP-value \n(Low vs High (>0.5)", size = 3, fontface =2) +
  # annotate(geom = "text", x = 1:3, y = -6, label = paste0("P=",round(pval_Wilc,5)), size = 5) +
  # annotate(geom = "text", x = 1:4, y = -7, label = paste0("P=",round(pval_dat$pval_Wilc,5)), size = 5) +
  theme_bw() + 
  ggtitle("Proteasome Inhibition (Bortezomib)") + 
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 10))) +
  
  theme(axis.title.x = element_text(size=35),
        axis.title.y = element_text(size=35),
        plot.title = element_text(size = 40, face = "bold", hjust=0.5),
        panel.grid.major.x  = element_blank(),
        legend.title=element_text(size=30) , legend.text=element_text(size=28),
        legend.position = "bottom",
        axis.text=element_text(size=20))

# pdf("UbiquitinationFC_Protein_attenuation.pdf" , width = 20, height = 15, paper = "special")
plot(p)
# dev.off()
# ------------------------------------------------------------------------------------------------------------
