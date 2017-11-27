### Plot CNA landscape with superimposed significance score 
melted_cna_Binary_Samples_Augmented = melt(cna_Binary_Samples_Augmented,
                                           id = c("pos","chr","gene"), 
                                           measure.vars = c("gains","losses"))


noSamples <- ncol(prot)
level = 1
levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications"))
molecular_view = 'Protein'

manhat_Significant_Score <- manhat.matrix(data1 = mrna, 
                                          data2 = prot, 
                                          molecular_view = molecular_view, 
                                          level = level)

# Plot 
dotplot.manhattan(toPlotCNA = melted_cna_Binary_Samples_Augmented, 
                  molecular_view = molecular_view, 
                  level = level)


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

manhat_Significant_Score$chromosome_Bands <- paste(paste0('chr',manhat_Significant_Score$chromosome_name),
                                                   manhat_Significant_Score$band, sep="_")

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


## Scatterplot of mRNA - Protein Significance score

# pdf(paste0(levelName," SpecialGenes_mRNA_Protein_Score_Scatterplot",".pdf"), width=15, height=9, paper='special')
significance_Score_ScatterPlot(genesToPlot = intCl754, 
                               level = level)
# dev.off()
