### BC functions


# Useful function to generate array data (Bengt Sennbland)
convert2Array=function(x){
  if(is.null(dim(x))){  # for vectors
    return(array(unlist(x)))
  }
  else{ # for lists
    return(array(unlist(x), dim=dim(x), dimnames=sapply(dimnames(x),make.names,allow_=FALSE)))
  }
}


# Wrapper for reading excel files 
readArrayFromExcel = function(file, cols.to.remove=c(), row.name.col=c(), sheet=1){
  library(readxl)
  ret=as.data.frame(read_excel(file,sheet=sheet, col_names = TRUE))
  # Set column row.name.col as row names and convert to array
  if(!is.null(row.name.col)){
    row.names(ret)=ret[,row.name.col]
  }
  ret=ret[,!names(ret) %in% c(row.name.col,cols.to.remove)]
  ret=convert2Array(ret)
  return(ret)
}


readListFromExcel = function(file, sheet=1) {
  library(readxl)
  raw = read_excel(file,sheet=sheet, col_names = TRUE)
  ret = list()
  for(key in names(raw)){
    tmp=as.vector(raw[key])
    ret[[key]]=as.data.frame(tmp)[!is.na(tmp)]
  }
  return(ret)
}


# Create gene-centric matrix for CNA data
gene.centric.cna <- function(fileToLoad) {
  cna = read.table(fileToLoad, header=TRUE, sep="\t") 
  col = colnames(cna)
  tumours = intersect(colnames(mrna), col)
  cna = cna[, c(col[1:6], tumours)]
  annot = getBM(attributes=c('efg_agilent_sureprint_g3_ge_8x60k','hgnc_symbol','chromosome_name','start_position','end_position','band'),
                filters ='efg_agilent_sureprint_g3_ge_8x60k',values = cna$probes, mart = ensembl_75)
  annot = annot[!annot$hgnc_symbol=="", ]
  annot$chromosome_name[annot$chromosome_name=="X"] = "23"
  chrnames = as.character(1:23)
  annot = annot[annot$chromosome_name %in% chrnames,]
  annot = annot[!duplicated(annot$efg_agilent_sureprint_g3_ge_8x60k),]
  idprobe = match(annot$efg_agilent_sureprint_g3_ge_8x60k, cna$probes)
  cna = cna[idprobe,]
  # print(paste("NAs =", sum(is.na(cna))))
  
  # Median of absolute cna probes
  cna = data.frame(GeneName = annot$hgnc_symbol, cna[, -c(1:6)])
  cna <- as.data.frame(cna %>% group_by(GeneName) %>% summarise_all(funs(median)))
  rownames(cna) = cna$GeneName; cna$GeneName = NULL
  cna = as.matrix(cna)
  return(cna)
}


# Estimate threshold for ploidy calling ("Gain", "Loss") for the samples 
# and for the whole cohort respectively
binarize.copy.number <- function(cna) {
  if ("gene" %in% colnames(cna)) {
    cna = cna[,-c(1:5)]
  }
  thres_gain =  0.6 + ploidy$Ploidy[match(colnames(cna), ploidy$SampleID)]
  # thres_hypgain = 4.6 + ploidy$Ploidy[match(tumours,ploidy$SampleID)]
  thres_loss = -0.6 + ploidy$Ploidy[match(colnames(cna), ploidy$SampleID)]
  
  ## Define threshold without ploidy adjustment
  # thres_gain = 3
  # thres_loss = 1
  
  res = apply(cna, 1, function(i) {
    ifelse(i >= thres_gain, 1, ifelse(i <= thres_loss, -1, 0))
  })
  res = t(res)
  return(res)
}

# Count copy number status 
count.gain.loss <- function(x) {
  gains <- apply(x,1, function(i) sum(i==1))
  losses <- apply(x,1, function(i) sum(i==-1))
  list(gains = gains, losses = losses)
}


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


# Add genomic location
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


# Differential expressionfunction
difExprByCNstatus <- function(omicsType, binary_CNA, gene, level) {
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


# Subtype enrichment function
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


# THreshold function
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


# Determine significance function
defSigFoldChangePvalue <- function(omicsType, binary_CNA, level) {
  testTable <- data.frame(do.call(rbind,lapply(genes, difExprByCNstatus, omicsType=omicsType, binary_CNA=binary_CNA, level=level)))
  colnames(testTable) = c("pvalTest","FoldChange", "Effect_Size")
  rownames(testTable) <- genes
  # testTable = testTable[!is.na(testTable$pvalTest),]
  testTable$pvalAdjusted = p.adjust(testTable$pvalTest, method = "BH")
  testTable$foldChangeSignificant =  setThreshold(testTable$FoldChange,level,c(0.05, 0.95))
  testTable$pvalAdjustedSignificant =  testTable$pvalAdjusted <= 0.1
  testTable$PvalueFoldChangeSignificant = testTable$foldChangeSignificant & testTable$pvalAdjustedSignificant
  colnames(testTable) =  c("pvalTest","FoldChange", "Effect_Size", "pvalAdjusted","foldChangeSignificant","pvalAdjustedSignificant","pvalueFoldChangeSignificant")
  testTable
}


# Generate volcano input data
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


# Plot function
toPlotVolcano <- function(omicsType1, omicsType2, binary_CNA, level, molecular_view) {
  dat <- returnVolcanoMatrix(omicsType1, omicsType2, binary_CNA, level)
  dat <- dat[[molecular_view]]
  xaxisneg = min(dat$FoldChange, na.rm=TRUE)
  xaxispos = max(dat$FoldChange, na.rm=TRUE)
  xaxis = ifelse(xaxispos>abs(xaxisneg), xaxispos,abs(xaxisneg))
  levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses", "Amplifications"))
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

# Manhattan plot
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

# Cis correlations function
cis.cor <- function(gene, data1, data2) {
  cor = corAndPvalue(data1[gene, ], data2[gene, ], method = "spearman")
  return(c('cor' = cor$cor, 'pvalue' = cor$p))
}


### Boxplot/Scatterplot per individual gene

## Boxplot
per.gene.boxplot <- function(data, binary_CNA, gene, molecular_view) {
  factor_new = revalue(as.factor(binary_CNA[gene,]), c("-1"="Loss", "0"="Neutral", "1"="Gain"))
  boxplot(data[gene, ] ~ factor_new,
          outline = TRUE,
          main = gene,
          xlab ="",
          ylab ="",
          cex.main=2,cex.axis=2)
  
  beeswarm(data[gene,] ~ factor_new,
           pwcol = colour_subtype, labels = colnames(data), pch = 16, cex = 4, add = TRUE)
  # plot(data[gene,] ~ factor_new,
  #        pwcol = colour_subtype, labels = colnames(data), pch = 16, cex = 4, add = TRUE)
  # legend("topleft", legend = molecular_subtype,
  # title = "Subtype", cex = 1.5, pch = 16, col = PAM50cols)
}


# Scatterplot
per.gene.scatterplot <- function(data1, data2, gene, molecular_view, cor_method) {
  plot(data1[gene,], data2[gene,], pch = 16, cex = 2.5, col = colour_subtype,
         main = ifelse(all(data2==mrna),gene,""),
         xlab = ifelse(all(data1==cna_log), "log2(cna)", "log2 (mRNA)"),
         ylab = paste0("log2","(", molecular_view,")"), cex.lab = 1.5, cex.main = 2)
    # legend("topleft", legend = molecular_subtype,
    #        title = "Subtype", pch = 16, cex = 2, col = PAM50cols,2)
    abline(lm(data2[gene,]~data1[gene,]))
    # text(x = floor(max(data1[gene,])), y = ceiling(min(data2[gene,])), cex=2, labels = paste0("rho=",round(cor(data1[gene,],data2[gene,],
    # method=cor_method),2)))
    legend("bottomright", legend = paste0("rho=",round(cor(data1[gene,],data2[gene,],method=cor_method),2)),
           title = "", cex = 2, pch = 16, col = NA, xpd=TRUE,bty="n")
  }


# Final table
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


generate.final.table <- function(level) {
  levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications")) 
  dat_mRNA <- manhat.matrix(data1 = omicsType1, data2 = omicsType2, molecular_view="mRNA", level)
  dat_Protein <- manhat.matrix(data1 = omicsType1, data2 = omicsType2, molecular_view="Protein", level)
  
  id <- colnames(dat_mRNA) %in% c("pvalTest","FoldChange","Effect_Size","pvalAdjusted", "foldChangeSignificant", 
                                  "pvalAdjustedSignificant","pvalueFoldChangeSignificant", "score")
  colnames(dat_mRNA)[id] <- paste("mRNA",colnames(dat_mRNA)[id],sep="_")
  
  colnames(dat_Protein)[id] <- paste("Protein",colnames(dat_Protein)[id],sep="_")
  
  dat_all <- merge(dat_mRNA, dat_Protein, by=colnames(dat_mRNA)[!id])
  colnames(dat_all)[-c(1:6)] <- paste(colnames(dat_all)[-c(1:6)], levelName, sep="_")
  colnames(dat_all)[1] = "genes"
  return(dat_all)
  # write.xlsx(dat_all, "dat_all.xlsx", row.names = FALSE, col.names = TRUE)
}

