### Plot CNA landscape with superimposed significance score 
melted_cna_Binary_Samples_Augmented = melt(cna_Binary_Samples_Augmented,id = c("pos","chr","gene"), 
                                           measure.vars = c("gains","losses"))

deviseScore <- function(data) {
  # ranked_Pval_Adjusted <- rank(-data$pvalAdjusted)/nrow(data)
  # 
  # ranked_FoldChange <- rank(abs(data$FoldChange))/nrow(data)
  # 
  # score <- (ranked_Pval_Adjusted + ranked_FoldChange)/2
  # 
  # score_new <- sign(data$FoldChange) * score
  score_new <- -log10(data$pvalAdjusted) * data$FoldChange
  
  return(score_new)
} 


manhat.matrix <- function(data1, data2, molecular_view, level) {
  
  tab <- returnVolcanoMatrix(data1, data2, cna_bin, level = level)
  
  tab_specific <- tab[[molecular_view]]
  
  augment_tab_specific <- augment_biomart(tab_specific)
  
  augment_tab_specific$colour <- ifelse(augment_tab_specific$bothPvalFoldChangeSignificant,"both",ifelse(
    augment_tab_specific$mrna_pvalFoldChangeSignificantOnly,"mrna",ifelse(augment_tab_specific$prot_pvalFoldChangeSignificantOnly,
                                                                          "protein","none")))
  augment_tab_specific$score <- deviseScore(augment_tab_specific)
  return(augment_tab_specific)
}

molecular_view <- "Protein"
level <- 1
levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications")) 

noSamples <- ncol(mrna)
manhat_Significant_Score <- manhat.matrix(mrna, prot, molecular_view, level)
# write.xlsx(manhat_Significant_Score, file = paste0("FC_pval_", molecular_view,"_", levelName,".xlsx"), row.names = FALSE, col.names = TRUE)


dotplot.manhattan <- function(toPlotCNA, molecular_view, level) {
  p <- ggplot() + 
    geom_bar(data = toPlotCNA, aes(x = pos,y = 2*value/noSamples, col = variable), stat = "identity", position = "identity") + 
    geom_point(data = manhat_Significant_Score, aes(x = pos ,y = score, fill = colour, colour=colour)) +
    scale_colour_manual(labels = c("","gains","losses","","",""), values = c("gains"="indianred1", "losses"="steelblue1", "both"="black", "protein"="blue","mrna"="green","none"=alpha("gray",0.3)), name="Effect") +
    scale_fill_manual(breaks = c("both","protein","mrna","none"), labels = c("Both mRNA/Protein","Protein only","mRNA only","none"),  values = c("both"="black", "protein"="blue","mrna"="green","none"=alpha("gray",1)), name="Significance level") +
    scale_x_continuous(name="Chromosomes", breaks= (chrlengths$start + chrlengths$end)/2,
                       labels=chrlengths$chr, limits=c(-1e8,chrlengths$end[23]), expand = c(0,0)) +
    scale_y_continuous(name = paste0(molecular_view," Significance Score \n( -log2FC * log10(q-value) )"), 
                       sec.axis = sec_axis(~ . * 0.5, name = "Samples with Copy Number Aberrations (%)", breaks = seq(from=-1,to=1,by=0.2),
                                           labels = abspos(seq(-100,100,by=20)))) +
    # scale_y_continuous(name = 'Copy Number Frequency (%)', breaks = seq(-1,1, by = 0.2), limits = c(-1,1),
    # labels = abspos(seq(-100,100,by=20)), expand=c(0,0) ) +
    geom_vline(xintercept = c(chrlengths$start, chrlengths$end), lty = 2, lwd = 0.5, colour = alpha("black",0.3)) +
    # geom_text(aes(label = chrlengths$chr,x = (chrlengths$start+chrlengths$end)/2, y = -0.95), size = 4) +
    # facet_grid(. ~ chr, scales="free", space = "free_x" ) +
    #geom_segment(aes(x = (chrlengths$start[1] + chrlengths$end[1])/2, y = -7.5, xend = (chrlengths$start[23]+chrlengths$end[23])/2, yend = -7.5),
    # colour = "black", size = 0.1) +
    # geom_segment(aes(x = -1e8, y = -5, xend = -1e8, yend = 5), 
    # colour = "black", size = 0.1) +
    guides(fill = guide_legend(order=1, override.aes = list(colour=c("black","blue","green","gray"), size=4,shape=rep(21,4))), colour = guide_legend(order=2, override.aes = list(size=2, colour=c(NA,"indianred1","steelblue1",NA,NA,NA), fill = c(NA,"indianred1","steelblue1",NA,NA,NA)))) +
    ggtitle(paste0('Significantly changed mRNA/protein abundances ', "across the copy number aberration landscape \n",
                   levelName," vs neutrals comparison (Protein Perspective"))   +
    theme(axis.title.x = element_text(size=15,face="bold"),
          axis.title.y = element_text(size=15),
          plot.title = element_text(size = 20, face = "bold")) 
  # theme_minimal() +
  # theme(panel.grid = element_blank(), axis.ticks.x = element_line(colour = "black"),axis.ticks.y = element_line(colour = "black"), axis.ticks.length=unit(0.1,"cm"),
  # legend.position = "none", panel.grid.major.x = element_blank(),plot.title = element_text(hjust = 0.5)) 
  # pdf(paste(molecular_view,levelName,"CN_score.pdf",sep="_"),
  # width = 25, height = 10, paper = "special")
plot(p)
# dev.off()
}

dotplot.manhattan(melted_cna_Binary_Samples_Augmented,molecular_view, level)

#----------------------------------------------------------------------------------------------------------------------------
## Zoom to chromosome regions

# Read band position
chromosome_Bands_Table <- read.table(paste0(filefolder,"cytoBand_hg19.txt"))
chromosome_Bands_Table$V1 <- as.character(chromosome_Bands_Table$V1)
chromosome_Bands_Table$V1[chromosome_Bands_Table$V1=='chrX'] = 'chr23'

# Extract position of significant CNA/mRNA/Protein correlation
idSig <- which(manhat_Significant_Score$bothPvalFoldChangeSignificant==TRUE)
chromosome_Bands_Significants <- paste(paste0('chr',manhat_Significant_Score$chromosome_name[idSig]),manhat_Significant_Score$band[idSig], sep="_")

chromosome_Bands_Id <- paste(chromosome_Bands_Table$V1, chromosome_Bands_Table$V4, sep="_")
# all(chromosome_Bands_Significants %in% chromosome_Bands_Id)
matched_Chromosome_Bands_Table <- unique(chromosome_Bands_Table[match(chromosome_Bands_Significants,chromosome_Bands_Id), ])
rownames(matched_Chromosome_Bands_Table) <- unique(chromosome_Bands_Significants)
colnames(matched_Chromosome_Bands_Table) <- c('chromosome','start','stop','band','stain')

manhat_Significant_Score$chromosome_Bands <- paste(paste0('chr',manhat_Significant_Score$chromosome_name),manhat_Significant_Score$band, sep="_")

# Set region
for(i in 1:nrow(matched_Chromosome_Bands_Table)) {
  genesToZoom <- manhat_Significant_Score$hgnc_symbol[manhat_Significant_Score$chromosome_Bands %in% rownames(matched_Chromosome_Bands_Table)[i]] 
  
  # Find the probe position that corresponds to the chromosome band
  chromosome_Bands_Significants_Ranges <- GRanges(seqnames=Rle(matched_Chromosome_Bands_Table$chromosome),
                                                  ranges=IRanges(matched_Chromosome_Bands_Table$start, matched_Chromosome_Bands_Table$stop),
                                                  strand=rep("+",nrow(matched_Chromosome_Bands_Table)))
  chromosome_Bands_Ranges <- GRanges(seqnames=Rle(paste0('chr',melted_cna_Binary_Samples_Augmented$chr)),
                                     ranges=IRanges(as.integer(melted_cna_Binary_Samples_Augmented$pos-chrlengths$start[as.integer(melted_cna_Binary_Samples_Augmented$chr)]), width=1),
                                     strand=rep("+",nrow(melted_cna_Binary_Samples_Augmented)))
  
  overlap_Regions <- findOverlaps(chromosome_Bands_Significants_Ranges[i],chromosome_Bands_Ranges)
  zoomed_melted_cna <- melted_cna_Binary_Samples_Augmented[overlap_Regions@to,]
  zoomed_manhat <- manhat_Significant_Score[manhat_Significant_Score$hgnc_symbol %in% genesToZoom,]
  zoomed_manhat$colour <- factor(zoomed_manhat$colour, levels = c("both","protein","mrna","none"))
  
  p <- ggplot() + 
    geom_bar(data = zoomed_melted_cna, aes(x = pos-chrlengths$start[as.integer(zoomed_melted_cna$chr)], y = value/noSamples, colour = variable), stat = "identity", position = "identity") +
    geom_line(data = zoomed_melted_cna, aes(x = pos-chrlengths$start[as.integer(zoomed_melted_cna$chr)], y = value/noSamples, colour = variable), stat = "identity", position = "identity") +
    geom_area(data = zoomed_melted_cna, aes(x = pos-chrlengths$start[as.integer(zoomed_melted_cna$chr)], y = value/noSamples, fill = variable), stat="identity", alpha = 0.2) + 
    geom_point(data = zoomed_manhat, aes(x = pos-chrlengths$start[as.integer(zoomed_manhat$chromosome_name)] ,y = score)  , size = 2) +
    geom_hline(yintercept = 0, lty = 1, lwd = 0.5, colour = alpha("black",0.2)) +
    geom_label_repel(data = zoomed_manhat, aes(x=pos-chrlengths$start[as.integer(zoomed_manhat$chromosome_name)], y=score,
                                               label = hgnc_symbol, fill = colour), col = "white", 
                     point.padding = unit(0.5, "lines"),
                     segment.color = 'grey50',
                     size = 6, fontface = "bold", force=2, na.rm=TRUE) +
    scale_colour_manual(values = c("gains"="indianred1", "losses"="steelblue1"), name="Copy Number Status") +
    scale_fill_manual(breaks = c("both","protein","mrna","none", "gains","losses"), labels = c("Both mRNA/Protein","Protein only","mRNA only","none", "",""),  values = c("both"="black", "protein"="blue","mrna"="green","none"=alpha("gray",1),"gains"="indianred1", "losses"="steelblue1"), name="Significance level", drop=FALSE) +
    scale_x_continuous(name= paste(matched_Chromosome_Bands_Table$chromosome[i],paste(matched_Chromosome_Bands_Table$start[i],
                                                                                      matched_Chromosome_Bands_Table$stop[i],sep="-"), sep = ": "),
                       # breaks = seq(from = matched_Chromosome_Bands_Table$start[i], to = matched_Chromosome_Bands_Table$stop[i], by= 1e6), 
                       limits = c(matched_Chromosome_Bands_Table$start[i],matched_Chromosome_Bands_Table$stop[i]), expand = c(0,0)) +
    scale_y_continuous(name = expression(atop(bold('Significance Score'), italic(paste("(",-log[10],"(q-value) * ",log[2],"FC)", sep="")))),
                       sec.axis = sec_axis(~ . * 1, name = "Samples with Copy Number Aberrations (%)", breaks = seq(from=-1,to=1,by=0.2),
                                           labels = abspos(seq(-100,100,by=20)))) +
    guides(fill = guide_legend(order=1, override.aes = list(colour=NA, fill=c("black","blue","green",alpha("grey",1),NA,NA))), colour = guide_legend(order=2, override.aes = list(size=2, fill=c("indianred1","steelblue1")))) +
    ggtitle(paste("Genomic Band:",rownames(matched_Chromosome_Bands_Table)[i], '\n', paste0('(', molecular_view, ' Level',')'))) +
    # scale_y_continuous(name = 'Copy Number Frequency (%)', breaks = seq(-1,1, by = 0.2), limits = c(-1,1),
    # labels = abspos(seq(-100,100,by=20)), expand=c(0,0) ) +
    # geom_text(aes(label = chrlengths$chr,x = (chrlengths$start+chrlengths$end)/2, y = -0.95), size = 4) +
    # facet_grid(. ~ chr, scales="free", space = "free_x" ) +
    # geom_segment(aes(x = matched_Chromosome_Bands_Table$start[10], y = -7.5, xend = matched_Chromosome_Bands_Table$stop[10], yend = -7.5),
    # colour = "black", size = 0.1) +
    # geom_segment(aes(x = matched_Chromosome_Bands_Table$start[10], y = -5, xend = matched_Chromosome_Bands_Table$start[10], yend = 5), 
    # colour = "black", size = 0.1) +
    # theme_minimal() +
    # theme(panel.grid = element_blank(), axis.ticks.x = element_line(colour = "black"),axis.ticks.y = element_line(colour = "black"), axis.ticks.length=unit(0.1,"cm"),
    # legend.position = "none", panel.grid.major.x = element_blank(),plot.title = element_text(hjust = 0.5)) 
  theme(axis.title.x = element_text(size=15,face="bold"),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 20, face = "bold"))  
  # tiff(filename = paste0(molecular_view, "_", levelName,"_", "Copy_Number_Band_Specific_score","_",rownames(matched_Chromosome_Bands_Table)[i],'.tiff'), res=300,
  # width = 30, height = 10, units = "in", pointsize = 12)
  # plot(p)
  # dev.off()
}
#----------------------------------------------------------------------------------------------------------------------------
# Scatterplot of Significance score

## Show IntCl genes
Nik_Zainal <- read.table(paste0(filefolder, "DriverGenes_Nik_Zainal.txt"), sep="\t")
Nik_Zainal <-  unique(as.character(Nik_Zainal$V1))

level = 1

levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications"))

significance_Score_ScatterPlot <- function(genesToPlot, level) {
  manhat_Significant_Score_mRNA <- manhat.matrix(mrna, prot, molecular_view ='mRNA', level)
  manhat_Significant_Score_Protein <- manhat.matrix(mrna, prot, molecular_view ='Protein', level)
  mRNA_Protein_Score <- data.frame(mRNA=manhat_Significant_Score_mRNA$score, Protein=manhat_Significant_Score_Protein$score) 
  rownames(mRNA_Protein_Score) = manhat_Significant_Score_mRNA$hgnc_symbol
  
  # Set threshold
  thres_mRNA <- quantile(defSigFoldChangePvalue(mrna, cna_bin,level)$FoldChange,c(0.05,0.95), na.rm=TRUE)
  
  thres_Protein <- quantile(defSigFoldChangePvalue(prot, cna_bin,level)$FoldChange,c(0.05,0.95), na.rm=TRUE)
  
  genesToColour <- rownames(mRNA_Protein_Score)  %in% genesToPlot
  
  # Number of IntCLuster genes in the comparison level
  # print(sum(genesToColour))
  
  thresSig <- if(levelName=="Gains" || levelName=="Amplifications") {
    mRNA_Protein_Score$mRNA >= thres_mRNA[2]  & mRNA_Protein_Score$Protein >= thres_Protein[2]} else {
      mRNA_Protein_Score$mRNA <= thres_mRNA[1]  & mRNA_Protein_Score$Protein <= thres_Protein[1]
    }
  if (length(genesToPlot)>1) {
    genesToLabel = thresSig & genesToColour
  } else {
    genesToLabel = thresSig 
    genesToColour = thresSig
  }
  
  p <- ggplot(data = mRNA_Protein_Score, aes(x = mRNA, y = Protein)) + 
    geom_hline(yintercept = if(level==1) {c(0, thres_Protein[2])} else {c(0, thres_Protein[1])}, col = alpha("black",0.5), lty = c(1,2)) + 
    geom_vline(xintercept = if(level==1) {c(0, thres_mRNA[2])} else {c(0, thres_mRNA[1])}, col = alpha("black",0.5), lty = c(1,2)) + 
    annotate("text",x = ifelse(level==1, thres_mRNA[2], thres_mRNA[1]), y=-3.5, 
             label="threshold q-value<=0.1 \ntop 5% log2FC", angle=90, size=2.5) + 
    annotate("text",x = 6, y=ifelse(level==1, thres_Protein[2], thres_Protein[1]), 
             label="threshold q-value<=0.1 \ntop 5% log2FC", size=2.5) +
    annotate("text",x = 4, y=-4, 
             label = paste0("% of labeled IntCluster genes", " = ", paste0(100*round(sum(genesToLabel,na.rm = TRUE)/sum(genesToColour),2),"%"), " (", sum(genesToLabel,na.rm = TRUE),"/", sum(genesToColour),")")) +
    geom_point() + 
    geom_point(data = mRNA_Protein_Score[genesToColour,], aes(x = mRNA, y = Protein, col='red')) + 
    geom_label_repel(data = mRNA_Protein_Score[genesToLabel,], aes(x = mRNA, y = Protein, col='red', label = rownames(mRNA_Protein_Score[genesToLabel,]))) +
    ggtitle(paste("mRNA - Protein Significance score", "\n", paste(levelName, "vs", "Neutrals"))) + 
    xlab(expression(atop(bold('mRNA'), italic(paste("(",-log[10],"(q-value) * ",log[2],"FC)", sep=""))))) +
    ylab(expression(atop(bold('Protein'), italic(paste("(",-log[10],"(q-value) * ",log[2],"FC)", sep=""))))) +
    guides(colour=FALSE) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size=20),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15))
  return(p)
}

# pdf(paste0(levelName," SpecialGenes_mRNA_Protein_Score_Scatterplot",".pdf"), width=15, height=9, paper='special')
significance_Score_ScatterPlot(intCl754, 1)
# dev.off()
