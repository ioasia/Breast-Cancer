### Linear Regression, ANOVA

# ANOVA ranked test
mrna_pvals <- sapply(genes, do.anova.rank, data = mrna_res)
prot_pvals <- sapply(genes, do.anova.rank, data = prot_res)

mrna_pvals_adj <- p.adjust(mrna_pvals, method = "BH")
prot_pvals_adj <- p.adjust(prot_pvals, method = "BH")

# mRNA -Protein Pearson correlation
cis_mrna_prot <- t(sapply(genes, cis.cor, data1 = mrna_res, data2 = prot_res))
cis_mrna_prot <- as.data.frame(cis_mrna_prot)
cis_mrna_prot$gene <- rownames(cis_mrna_prot)

# Merge data
anovaQvalueCor <- data.frame(intClgene = cis_mrna_prot$gene %in% intCl754, 
                             qval = mrna_pvals_adj,
                             cor = cis_mrna_prot$cor)


p <- ggplot(data = anovaQvalueCor, mapping = aes(x = -log10(qval), y = cor, colour = intClgene)) + 
  geom_point() + 
  scale_x_continuous(name = 'q-value CNA-mRNA (-log10)', expand = c(0,0), limits = c(0,6), labels = math_format()) + 
  scale_y_continuous(name = 'mRNA - Protein correlation \n(Spearman)', limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.5)) +
  scale_colour_manual(name = '', labels = c("", "Curtis et al., Inclust classifier genes"),
                      values = c('TRUE' =  'red', 'FALSE' = alpha('grey', 0.3)))+
  
  theme(legend.position = 'top') + 
  guides(colour = guide_legend(override.aes = list(colour = c(NA, 'red'))))
  
 
plot(p)
## Percentage of significant results
# perc_cna_mrna <- round(sum(mrna_pvals_adj<=0.1)/nrow(mrna) * 100,2)
# perc_cna_prot <- round(sum(prot_pvals_adj<=0.1)/nrow(prot) * 100,2)
# sum(mrna_pvals_adj<=0.1 & prot_pvals_adj<=0.1)/ sum(mrna_pvals_adj<=0.1) * 100
# sum(mrna_pvals_adj<=0.1 & prot_pvals_adj<=0.1)/ sum(prot_pvals_adj<=0.1) * 100

# # Plot q-values 
# pval_dat <- data.frame(gene = genes,mrna = -log10(mrna_pvals_adj), prot = -log10(prot_pvals_adj),
#                        sig = ifelse(rownames(mrna) %in% intCl754, "yes", "no"))
# 
# scatter <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
#   geom_rect(aes(xmin=4, xmax=5.5, ymin=2, ymax=3.5), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
#   geom_rect(aes(xmin=4, xmax=5.5, ymin=0, ymax=1), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
#   geom_rect(aes(xmin=0, xmax=1, ymin=2, ymax=3.5), color="grey", fill= alpha("aliceblue",0.3), linetype = 2) + 
#   
#   geom_point() + 
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein (-log10)", expand = c(0.05, 0.05), breaks = seq(0,4,1),
#                      labels = seq(0,4,1)) +
#   scale_x_continuous(name = "mRNA (-log10)", expand = c(0.05, 0.05), breaks = seq(0,6,1),
#                      labels = seq(0,6,1)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = alpha('grey',0.50),yes ='blue')) +
#   # ggtitle("mRNA/Protein ANOVA adjusted P values") +
#   guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 3))) + 
#   geom_vline(xintercept=1, colour = "black", linetype=2) + 
#   geom_hline(yintercept=1, colour = "black", linetype=2) + 
#   geom_text(aes(x=1, label="mRNA", y=3.5), colour=alpha("black",0.02), angle=90, vjust = -1, size = 6) +
#   geom_text(aes(y=1, label="Protein", x= 5), colour=alpha("black",0.02), vjust = -1, size = 6) +
#   geom_text(aes(x=1, label="threshold", y=3.5), colour=alpha("black",0.02), angle=90, vjust = 1, size = 6) +
#   geom_text(aes(y=1, label="threshold", x=5), colour=alpha("black",0.02), vjust = 1, size = 6) + 
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text.x=element_text(size=15, hjust = 1),
#         axis.text.y=element_text(size=15, hjust = 1),
#         legend.key.height = unit(1, "cm"),
#         legend.key.width= unit(1, "cm"),
#         legend.key.size= unit(1, "cm")) 
# 
# density_x <- ggplot(data = pval_dat, aes(x = mrna, colour = sig, fill = sig)) + 
#   geom_density(alpha = 0.5) + 
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein", expand = c(0.05, 0.05)) +
#   scale_x_continuous(name = "mRNA", expand = c(0.05, 0.05)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) + 
#   scale_fill_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
#   guides(colour = FALSE, fill=FALSE) + 
#   ggtitle("mRNA/Protein ANOVA adjusted P values") +
#   theme(plot.margin      = unit(c(1, 0, 0, 1.4), "cm"),
#         plot.title = element_text(size = 26, face = "bold", hjust=0.5),
#         axis.ticks       = element_blank(), 
#         panel.background = element_blank(), 
#         panel.grid       = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x      = element_blank(), 
#         axis.text.y      = element_blank(),           
#         axis.title.x     = element_blank(), 
#         axis.title.y     = element_blank()) 
# 
# 
# density_y <- ggplot(data = pval_dat, aes(x = prot, colour = sig, fill = sig)) +
#   geom_density(alpha = .5) +
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein", expand = c(0.05, 0.05)) +
#   scale_x_continuous(name = "mRNA", expand = c(0.05, 0.05)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) + 
#   scale_fill_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
#   guides(colour = FALSE, fill=FALSE) + 
#   theme(plot.margin      = unit(c(0, 0, 1.45, 0), "cm"),
#         axis.ticks       = element_blank(), 
#         panel.background = element_blank(), 
#         panel.border = element_blank(),
#         panel.grid       = element_blank(),
#         axis.text.x      = element_blank(), 
#         axis.text.y      = element_blank(),           
#         axis.title.x     = element_blank(), 
#         axis.title.y     = element_blank()) + 
#   coord_flip()
# 
# # empty <- ggplot() +
# #   geom_point(aes(1,1), colour="white") +
# #   theme(plot.margin      = unit(c(0, 0, 0, 0), "cm"),
# #         axis.ticks       = element_blank(), 
# #         panel.background = element_blank(), 
# #         panel.grid.minor = element_blank(),
# #         axis.text.x      = element_blank(), 
# #         axis.text.y      = element_blank(),           
# #         axis.title.x     = element_blank(), 
# #         axis.title.y     = element_blank())
# 
# legend_scatter <- get.legend(scatter)
# 
# scatter <- scatter + theme(legend.position="none")
# 
# # pdf("IntClust_Anova_Plot.pdf", width=20, height = 10)
# grid.arrange(density_x, legend_scatter, scatter, density_y, 
#                                 ncol = 2, 
#                                 nrow = 2, 
#                                 widths  = c(4, 1), 
#                                 heights = c(1, 4))
# # dev.off()
# 
# ur <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(2,3.5)) +
#   scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(4,5.5)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
#   ggtitle("mRNA/Protein significant") + 
#   geom_label_repel(aes(label = gene)) + 
#   # guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) + 
#   guides(colour=FALSE) + 
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         panel.grid.major.x  = element_blank(),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text.x=element_text(size=7, hjust = 1),
#         legend.key.height = unit(0.5, "cm"),
#         legend.key.width= unit(0.5, "cm")) 
# 
# # pdf("ur_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
# plot(ur)
# # dev.off() 
# 
# dr <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(0,1)) +
#   scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(4,5.5)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +
#   ggtitle("mRNA only significant") + 
#   geom_label_repel(aes(label = gene)) + 
#   # guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) + 
#   guides(colour=FALSE) + 
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         panel.grid.major.x  = element_blank(),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text.x=element_text(size=7, hjust = 1),
#         legend.key.height = unit(0.5, "cm"),
#         legend.key.width= unit(0.5, "cm")) 
# 
# # pdf("dr_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
# plot(dr)
# # dev.off() 
# 
# ul <- ggplot(data = pval_dat, aes(x = mrna, y = prot, colour = sig)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_y_continuous(name =  "Protein(-log10)", expand = c(0.05, 0.05), limits = c(2,3.5)) +
#   scale_x_continuous(name = "mRNA(-log10)", expand = c(0.05, 0.05), limits = c(0,1)) + 
#   scale_color_manual(name = '', labels = c("","IntClust"), values=c("no" = 'black',yes ='blue')) +  ggtitle("Protein only significant") + 
#   geom_label_repel(aes(label = gene)) + 
#   guides(colour = guide_legend(override.aes = list(fill=c(NA, 'blue'), colour=c(NA, 'blue'), size = 1))) +
#   # guides(colour=FALSE) + 
#   theme(axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         plot.title = element_text(size = 28, face = "bold", hjust=0.5),
#         panel.grid.major.x  = element_blank(),
#         legend.title=element_text(size=25) , legend.text=element_text(size=20),
#         axis.text.x=element_text(size=7, hjust = 1),
#         legend.key.height = unit(0.5, "cm"),
#         legend.key.width= unit(0.5, "cm")) 
# 
# # pdf("ul_IntCl_quartiles.pdf", width = 10, height = 10, paper = 'special')
# plot(ul)
# # dev.off() 
# 
# # pdf("IntCl_quartiles.pdf", width = 30, height = 10, paper = 'special')
# grid.arrange(ur, dr, ul,
#              ncol = 3, 
#              nrow = 1)
# # dev.off()


## Venn plot
# total_genes <- nrow(mrna)
# total_IntClGenes <- length(intCl754)
# overlap_IntClGenes <- length(intersect(rownames(mrna),intCl754))
# 
# write.table(pval_dat, "venn.txt", quote = FALSE, row.names = FALSE, 
#             col.names = FALSE, sep = "\t")
# 
# # pdf("Total_study_IntClust_Venn.pdf", width=20, height=15, paper='special')
# draw.pairwise.venn(
#   area1 = total_genes,
#   area2 = total_IntClGenes,
#   cross.area = overlap_IntClGenes,
#   fill = c("black", "blue"),
#   label.col = 'white',
#   category =  c("Total study genes", "IntClust genes"),
#   cat.pos = c(-45, 45),
#   euler.d = TRUE,
#   cat.dist = c(0.03,0.03),
#   # sep.dist = 0.03,
#   offset = 1,
#   cat.cex=2,
#   ext.text = FALSE,
#   cex=1.5,
#   rotation.degree = 45,
#   lwd = 2,
#   print.mode = c('raw', 'percent')
# )
# # dev.off()
# 
# mrna_sig_IntCl <- sum(pval_dat$mrna > 1 & pval_dat$sig == 'yes')
# prot_sig_IntCl <- sum(pval_dat$prot > 1 & pval_dat$sig == 'yes')
# both_sig_IntCl <- sum(pval_dat$mrna > 1 & pval_dat$sig == 'yes' & pval_dat$prot > 1)
# 
# # pdf("mRNA_Protein_IntClust.pdf", width=20, height=15, paper='special')
# draw.pairwise.venn(
#   area1 = mrna_sig_IntCl,
#   area2 = prot_sig_IntCl,
#   cross.area = both_sig_IntCl,
#   fill = c('palegreen', 'steelblue1'),
#   category =  c("mRNA significant", " Protein significant"),
#   cat.pos = c(-45, 45),
#   euler.d = TRUE,
#   cat.dist = c(0.03,0.03),
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
#----------------------------------------------------------------------------------------------------------------------------------------------