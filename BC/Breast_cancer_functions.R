### Useful functions and data to load

# Separate legend plotting
get.legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Axis labels as absolute values
abspos <- function(x) {
  format(abs(x))
}

augment_biomart <- function(data, ensembl = ensembl_75, filter ='hgnc_symbol',values = rownames(data)) {
  # Find chromosomal position of gene symbols
  annot = getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position','band'),
                filters = filter, values = values, mart = ensembl)
  annot = annot[!annot[[filter]]=="",]
  annot$chromosome_name[annot$chromosome_name=="X"] = "23"
  chrnames = as.character(1:23)
  annot = annot[annot$chromosome_name %in% chrnames,]
  
  # Match rownames(genes) to annnotation 
  idgene1 = match(values,annot[[filter]])
  idgene1 = idgene1[!is.na(idgene1)]
  annot = annot[idgene1, ]
  
  # Match annotation to rownames(genes)
  idgene2 = match(annot[[filter]],values)
  data = data[idgene2, ]
  pos <- as.integer(annot$start_position) + chrlengths$start[as.integer(annot$chromosome_name)]
  cbind(annot,pos,data)
}


difExprByCNstatus <- function(omicsType,binary_CNA,gene,level) {
  pval = NA
  foldChange = NA
  z_statistic = NA
  toTest = omicsType[gene,binary_CNA[gene,]==level]
  Neutral = omicsType[gene,binary_CNA[gene,]==0]
  # Neutral = omicsType[gene,binary_CNA[gene,]!=level]
  ## No of CN gains to consider
  if (length(toTest) > 2) {
    # res = t.test(toTest,Neutral, var.equal = FALSE)
    res <- wilcox.test(toTest,Neutral)
    z_stat_factor <- factor(c(rep('toTest',length(toTest)),rep("Neutral",length(Neutral))))
    z_stat_values <- c(toTest,Neutral)
    z_statistic <- wilcox_test(z_stat_values~z_stat_factor,distribution="exact")
    z_statistic <- abs(z_statistic@statistic@teststatistic)/sqrt(length(c(toTest,Neutral)))
    pval <- res$p.value
    foldChange = log2(mean(2^toTest)/mean(2^Neutral))
  }
  c(pval,foldChange,z_statistic)
}

subtypEnrich <- function(binary_CNA,gene,level,subtype) {
  hypToTest = sum(binary_CNA[gene,]==level & setSubtypes == subtype)
  hypTotal = ncol(binary_CNA)
  hypSubtype = sum(setSubtypes == subtype)
  hypToPick = sum(binary_CNA[gene,]==level)
  pvalHyp = 1 - phyper(hypToTest-1, hypSubtype, hypTotal-hypSubtype, hypToPick, lower.tail=TRUE)
  ## Fisher test to match hypergeometric distribution
  # fishRest = sum(binary_CNA[gene,]!=level & setSubtypes != subtype)
  # toTestMatrix <- matrix(c(hypToTest, hypToPick-hypToTest, hypSubtype-hypToTest, fishRest),nrow=2)
  # pval = fisher.test(toTestMatrix, alternative="greater")$p.value
  pvalHyp
}


setThreshold <- function(foldchange,level,quantiles) {
  #  Ensure positive (negative) thresholds for gains (losses)
  
  #  if (level==1) {
  #    thres <- quantile(foldchange,quantiles[2])
  #    foldchange>=thres
  #  } else if(level==-1) {
  #    thres <- quantile(foldchange,quantiles[1])
  #    foldchange<=thres
  # }
  
  # Set threshold regardless of the sign
  thres <- quantile(foldchange,quantiles, na.rm = TRUE)
  idx = foldchange<=thres[1] | foldchange>=thres[2]
  idx
}


defSigFoldChangePvalue <- function(omicsType,binary_CNA,level) {
  testTable <- data.frame(do.call(rbind,lapply(genes, difExprByCNstatus, omicsType=omicsType, binary_CNA=binary_CNA, level=level)))
  colnames(testTable) = c("pvalTest","FoldChange", "Effect_Size")
  rownames(testTable) <- genes
  # testTable = testTable[!is.na(testTable$pvalTest),]
  testTable$pvalAdjusted = p.adjust(testTable$pvalTest, method = "BH")
  testTable$foldChangeSignificant =  setThreshold(testTable$FoldChange,level,c(0.05,0.95))
  testTable$pvalAdjustedSignificant =  testTable$pvalAdjusted <= 0.1
  testTable$PvalueFoldChangeSignificant = testTable$foldChangeSignificant & testTable$pvalAdjustedSignificant
  colnames(testTable) =  c("pvalTest","FoldChange", "Effect_Size", "pvalAdjusted","foldChangeSignificant","pvalAdjustedSignificant","pvalueFoldChangeSignificant")
  testTable
}


returnVolcanoMatrix <- function(omicsType1,omicsType2,binary_CNA,level) {
  # Prepares a matrix to be displayed as a volcano plot
  table1 <- defSigFoldChangePvalue(omicsType1,binary_CNA,level)
  table2 <- defSigFoldChangePvalue(omicsType2,binary_CNA,level)
  table3 <- data.frame(do.call(rbind,lapply(rownames(omicsType1), function(x) {
    sapply(as.character(molecular_subtype), subtypEnrich, gene = x, binary_CNA=binary_CNA, level=level)
  })))
  colnames(table3) = c(paste(as.character(molecular_subtype),"HyperTest.pvalue", sep="_"))
  rownames(table3) = rownames(omicsType1)
  table3 = table3[rownames(table3) %in% rownames(table2),]
  
  # sum(sapply(1:nrow(table3), function(i) {
  #   sum(table3[i,] <= 0.05)
  #   })>1)
  
  table3$enrichedType = sapply(1:nrow(table3), function(i) {
    ifelse(sum(table3[i,]<=0.05)>0, as.character(molecular_subtype[which(table3[i,]<=0.05)]), NA)
  })
  table3$enrichedType <- factor(table3$enrichedType, levels = molecular_subtype)
  # levels(table3$enrichedType) <- c("Basal","Her2","LumB","LumA","Normal")
  
  combinedTable1 = data.frame(table1,table3)
  combinedTable2 = data.frame(table2,table3)
  bothPvalFoldChangeSignificant <- combinedTable1$pvalueFoldChangeSignificant & combinedTable2$pvalueFoldChangeSignificant
  mrna_pvalFoldChangeSignificantOnly <-  combinedTable1$pvalueFoldChangeSignificant & !bothPvalFoldChangeSignificant
  prot_pvalFoldChangeSignificantOnly <-  combinedTable2$pvalueFoldChangeSignificant & !bothPvalFoldChangeSignificant
  enrichedBothPvalFoldChangeSignificant <- bothPvalFoldChangeSignificant & table3$enrichedType %in% molecular_subtype
  combinedTable1 <- cbind(combinedTable1,mrna_pvalFoldChangeSignificantOnly,prot_pvalFoldChangeSignificantOnly,bothPvalFoldChangeSignificant,enrichedBothPvalFoldChangeSignificant)
  combinedTable2 <- cbind(combinedTable2,mrna_pvalFoldChangeSignificantOnly,prot_pvalFoldChangeSignificantOnly,bothPvalFoldChangeSignificant,enrichedBothPvalFoldChangeSignificant)
  list(mRNA=combinedTable1, Protein=combinedTable2)
}

# write.table(rownames(dat[dat$BothPvalFoldChangeSignificant==TRUE,]),"BothGenes.txt", quote=FALSE, row.names = FALSE,col.names = FALSE)



toPlotVolcano <- function(omicsType1, omicsType2, binary_CNA, level, molecular_view) {
  dat <- returnVolcanoMatrix(omicsType1, omicsType2, binary_CNA, level)
  dat <- dat[[molecular_view]]
  xaxisneg = min(dat$FoldChange, na.rm=TRUE)
  xaxispos = max(dat$FoldChange, na.rm=TRUE)
  xaxis = ifelse(xaxispos>abs(xaxisneg),xaxispos,abs(xaxisneg))
  imag <- ggplot(data = dat, aes(x=FoldChange, y=-log10(pvalAdjusted))) + 
    # geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf),fill= alpha("khaki1",0.03)) +
    # geom_rect(data=NULL,aes(xmin=0,xmax=-Inf,ymin=-Inf,ymax=Inf),fill= alpha("lightsalmon",0.03)) +
    geom_point(size=1,na.rm=TRUE) +   
    geom_point(data = dat[dat$mrna_pvalFoldChangeSignificantOnly==TRUE,], aes(x=FoldChange, y=-log10(pvalAdjusted), col = "mRNA only"), size = 4) + 
    geom_point(data = dat[dat$prot_pvalFoldChangeSignificantOnly==TRUE,], aes(x=FoldChange, y=-log10(pvalAdjusted), col = "Protein only"), size = 4) +
    geom_point(data = dat[dat$bothPvalFoldChangeSignificant==TRUE,], aes(x=FoldChange, y=-log10(pvalAdjusted), col = "mRNA/Protein"), size = 4) + 
    scale_x_continuous(name = paste0(levelName," vs Neutrals (log2 Fold Change)"),breaks=seq(round(-xaxis),round(xaxis),1), lim = c(-xaxis-0.1,xaxis+0.1), expand = c(0,0)) +
    scale_y_continuous(name = "-log10(qValue)", expand=c(0,0)) +
    geom_vline(xintercept = if(molecular_view=="mRNA") {
      quantile(defSigFoldChangePvalue(omicsType1, binary_CNA,level)$FoldChange,c(0.05,0.95), na.rm=TRUE)} else {
        quantile(defSigFoldChangePvalue(omicsType2, binary_CNA,level)$FoldChange,c(0.05,0.95), na.rm=TRUE)
      }, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 1, lwd = 0.2) +
    geom_hline(yintercept = 1, linetype = 2) + 
    geom_label_repel(data = dat[dat$enrichedBothPvalFoldChangeSignificant==TRUE,], aes(x=FoldChange, y=-log10(pvalAdjusted),
                                                                                       label = rownames(dat[dat$enrichedBothPvalFoldChangeSignificant==TRUE,]), 
                                                                                       fill = unlist(dat$enrichedType[dat$enrichedBothPvalFoldChangeSignificant==TRUE])), col = "white", size = 6, fontface = "bold", force=2, na.rm=TRUE,
                     segment.color = 'grey50') +
    scale_color_manual(name = "Level of significance", breaks = c("mRNA/Protein","mRNA only", "Protein only"),  values = c("mRNA/Protein"="orange", "mRNA only"="green","Protein only" = "blue")) +
    scale_fill_manual(name = "Subtype Enrichment", breaks =c("Basal","Her2","LumA","LumB","Normal"), values = c("Basal" = "#E31A1C","Her2"="#FB9A99","LumA"="#1F78B4","LumB"="#A6CEE3", "Normal"="#33A02C"), drop=FALSE) +
    guides(fill = guide_legend(order=1,override.aes = list(colour=NA)),colour = guide_legend(order=2,override.aes = list(size=10))) +
    theme_cowplot() +
    theme(text = element_text(size=25),plot.title = element_text(size=30),legend.key.size = unit(3, 'lines')) +
    ggtitle(paste0(levelName," Vs Neutral Volcano plot \n(",molecular_view," Level)"))
  imag
}

# -----------------------------------------------------------------------------------------------------------

## Set variables
omicsType1 <- mrna
omicsType2 <- prot
binary_CNA <- cna_bin
level <- 1
molecular_view <- "Protein"

levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications"))  
p <- toPlotVolcano(omicsType1, omicsType2, binary_CNA, level, molecular_view)

# pdf(paste0(levelName,".",molecular_view, ".","Volcano_plot",".pdf"), width=20, height=15, paper='special')
plot(p)
# multiplot(imag[[1]], imag[[2]], cols=2)
dev.off()

#--------------------------------------------------------------------------------------------------




