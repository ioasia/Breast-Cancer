### Setup libraries
install.packages("easypackages")
library(easypackages)

source("http://bioconductor.org/biocLite.R")
toInstallBio <- c("impute", "preprocessCore", "GO.db", "biomaRt","GenomicRanges","CGHregions")

source("http://bioconductor.org/biocLite.R")
biocLite(toInstallBio)

toInstall=c("R6","stringi","reshape", "WGCNA","readxl",'rJava', 'xlsxjars',"xlsx","matrixStats","data.table", "beeswarm","dplyr", "plyr", "ggthemes", "perm",
            "cowplot","ggrepel","tibble", "reshape2","ggplot2","doParallel","ggrepel", "lmPerm","DescTools",
            'lsr',"coin","gsubfn","gridExtra","VennDiagram", 'gsubfn', 'venneuler', 'ggsignif', 'scales')
packages(toInstall)

libraries(c(toInstall,toInstallBio))
rm(toInstall,toInstallBio)


# Specify filefolder
filefolder <- "/Users/johnsiavelis/Downloads/Breast Cancer/data/"


# Protein data
prot = readArrayFromExcel(paste0(filefolder, "diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx"), 
                               row.name.col="GeneSymbol")
prot = log2(prot)


# Correct typos
typos=c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newcolnames=colnames(prot)
for(typo in typos){
  newcolnames=gsub(typo,paste(typo,"T1",sep=""), newcolnames)
}
colnames(prot) = newcolnames


# mRNA data
mrna = read.table(paste0(filefolder, "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne_updated.txt"), 
            'lsr',"coin","gsubfn","gridExtra","VennDiagram", 'gsubfn', 'venneuler')
packages(toInstall)

libraries(c(toInstall,toInstallBio))

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

# PROTEIN DATA
##############
# Read the expression data extracted from the original preoteome data from 
# Henrik Johansson, Lehtiö group

filefolder <- "/Users/johnsiavelis/Downloads/Breast Cancer/data/"
prot_data = readArrayFromExcel(paste0(filefolder, "diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx"), 
                               row.name.col="GeneSymbol")
prot_data = apply(prot_data, c(1,2), log2)
head(prot_data, 5)
# Correct 
typos=c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newcolnames=colnames(prot_data)
for(typo in typos){
  print(paste(typo,"T1",sep=""))
  newcolnames=gsub(typo,paste(typo,"T1",sep=""), newcolnames)
}
colnames(prot_data) = newcolnames
colnames(prot_data) =  gsub("_",".",colnames(prot_data))
colnames(prot_data) =  gsub("-",".",colnames(prot_data))

## Brief exploratory data analysis
# summary(prot_data)
# hist(prot_data, breaks=1000)
# boxplot(prot_data)

# MRNA ADATA
############
# Read the original RNA data from Kristine Kleivi in Oslo
rna_data = read.table(paste0(filefolder, "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne_updated.txt"),                      header=TRUE, 
                      sep="\t", 
                      stringsAsFactors=FALSE)

# Correct erroneous value jan.27 with 1.3527
mrna$OSL2U.0383T1[mrna$ProbeUID.SampleArray.gDetrendedSignal == "12629"] = 1.3527
mrna$OSL2U.0383T1 = as.double(mrna$OSL2U.0383T1)

# Remove duplicate
mrna$OSL2U.0219T1 = NULL

# Remove superfluous columns
mrna=mrna[,!names(mrna) %in% c("ProbeUID.SampleArray.gDetrendedSignal","ProbeName","SystematicName","chrom","seq_beg","seq_end","accessions")]

# Replace the rows of individual probes mapping to the same genename by their median 
# row (for each column) and finally convert to array
mrna <- as.data.frame(mrna %>% group_by(GeneName) %>% summarise_all(funs(median)))
rownames(mrna) = mrna$GeneName
mrna$GeneName = NULL
mrna = as.matrix(mrna)

# Metadata
# Correct erroneous value jan.27 with 1.3527 as per Miriam email
# Notice that this error caused the col to be coerced to factor and we need to
# convert back to double
rna_data$OSL2U.0383T1[rna_data$ProbeUID.SampleArray.gDetrendedSignal=="12629"] = 1.3527
rna_data$OSL2U.0383T1=as.double(rna_data$OSL2U.0383T1)

# #Correct "sep" date values to "SEP" gene names
# idx = grep("sep",rna_data$GeneName)
# to_rep_genes = rna_data$GeneName[idx]
# to_rep_genes = gsub("\\.","", to_rep_genes); to_rep_genes = gsub("p0","p",to_rep_genes);to_rep_genes = gsub("p","pt",to_rep_genes); to_rep_genes = toupper(to_rep_genes)
# rna_data$GeneName[idx] = to_rep_genes

# remove duplicate according to Kristines email
rna_data$OSL2U.0219T1 = NULL

# remove superfluous columns
rna_data=rna_data[,!names(rna_data) %in% c("ProbeUID.SampleArray.gDetrendedSignal","ProbeName","SystematicName","chrom","seq_beg","seq_end","accessions")]
# replace the rows of individual probes mapping to the same genename by their median 
# row (for each column) and finally convert to array
#For the whole cohort: colnames(rna_data)[1] = 'GeneName'
rna_data <- as.data.frame(rna_data %>% group_by(GeneName) %>% summarise_all(funs(median)))
rownames(rna_data) = rna_data$GeneName; rna_data$GeneName = NULL
rna_data = as.matrix(rna_data)
# rna_data=convert2Array(mymean(rna_data))

# META DATA
###########
# Read in the tumour meta data, i.e., PAM50 classification of the tumours
meta_tumours <- readArrayFromExcel(paste0(filefolder, "clin_info_connect_proteomics.xlsx"),
                                   row.name.col= "OSL2_full_name")
# Correct 
typos=c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newrownames=rownames(meta_tumours)
for(typo in typos){
  newrownames=gsub(typo,paste(typo,"T1",sep="/"), newrownames)
}
rownames(meta_tumours) = newrownames


meta_data <- as.data.frame(readArrayFromExcel(paste0(filefolder,"clin_info_connect_proteomics.xlsx")))
meta_data$OSL2.full.name <- newrownames
meta_data$OSL2.full.name <- gsubfn(".", list("-" = ".", "_" = ".", "/" = ""), meta_data$OSL2.full.name)
rownames(meta_data) <- meta_data$OSL2.full.name


# Read in the gene meta data, i.e., a list of PAM50/IntCl genes
pam50 = readListFromExcel(paste0(filefolder,"PAM50_and_Polyak_Cell_Snapshot-gene_lists.xlsx"))
# intCl1000 <- read.table(paste0(filefolder, "intCl1000.txt"))
# intCl1000 <- as.character(unique(intCl1000$V1))
intCl754 <- unique(read.table(paste0(filefolder,"short_IntCl_genes.txt")))
intCl754 <- as.character(intCl754$V1)
Nik_Zainal <- read.table(paste0(filefolder, "DriverGenes_Nik_Zainal.txt"), sep="\t")
Nik_Zainal <-  unique(as.character(Nik_Zainal$V1))

meta_genes = readListFromExcel(paste0(filefolder,"Genelists-summary_v2.xlsx"))
meta_genes["PAM50"] = pam50
meta_genes["IntCl754"] <- list('intCl754'=intCl754)
newnames = names(meta_genes)
newnames=gsub(" ","_", newnames)
names(meta_genes) = newnames

# Ensembl
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# Copy number data
cna_log <- gene.centric.cna(fileToLoad = paste0(filefolder, 'Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt'))
cna_tot <- gene.centric.cna(fileToLoad = paste0(filefolder, 'Gene_probe_centric.totCN.OSL2_n331_updated.txt'))

# Common samples
tumours = intersect(intersect(colnames(prot),
                              colnames(mrna)),
                    colnames(cna_log))

# Absolute Copy Number file  
cna_Absolute_number_Cohort = read.table(paste0(filefolder, "Gene_probe_centric.totCN.OSL2_n331_updated.txt"), header=TRUE, sep="\t")
cna_Absolute_number_Samples = cna_Absolute_number_Cohort[ ,c("probes","gene","chr","start","stop",tumours)]
cna_tot_probes = cna_Absolute_number_Samples[ ,tumours]

# Load ploidy file
ploidy = read.table(paste0(filefolder, "OSL2_n331_ploidy_t.perc.txt"), sep = "\t", header = TRUE)


head(pam50, 5)
intCl1000 <- read.table(paste0(filefolder, "intCl1000.txt"))
intCl1000 <- as.character(unique(intCl1000$V1))
intCl754 <- unique(read.table(paste0(filefolder,"short_IntCl_genes.txt")))
intCl754 <- as.character(intCl754$V1)

meta_genes = readListFromExcel(paste0(filefolder,"Genelists-summary_v2.xlsx"))
meta_genes["PAM50"] = pam50
meta_genes["IntCl1000"] <- list('intCl1000'=intCl1000)
newnames = names(meta_genes)
newnames=gsub(" ","_", newnames)
names(meta_genes) = newnames
names(meta_genes)

# CNA data
##########
# read in the CNA data,
filename = paste0(filefolder, "Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt")
cna=read.table(filename, header=TRUE, sep="\t") 
# This data contains dupliate lines -- remove them!
cna=cna[!duplicated(cna,MARGIN=1),]
#cnv$chr[cnv$chr==23] = "X"
colnames(cna)[colnames(cna)=="gene"]="GeneName"
cna = cna[, names(cna) %in% c("probes", "GeneName", "SystematicName", "chr", "start", "stop", colnames(prot_data))]

#Replace -Inf values to the row minimum
cna = cna[-which(apply(cna=="-Inf",1,sum)>=1),]

# OUTPUT
########
# Bundle protein and rna data arrays together with some info in a list
exprdata=list(name="Tumour expressionData",
              readme=paste(
                "This list comprise proteome expression data and mRNA",
                "expression data for 45 breast cancer tumours as 2 separate",
                "arrays indexed by tumours (columns) and gene symbol (rows),",
                "as well as GW cnv data for the same tumours. It is created",
                "by the script formatData.R from the original data files",
                "diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx",
                "(protein data from Henrik Johansson, Lehtiö group),",
                "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne.txt",
                "(mRNA data, from Kristine Kleivi, Oslo group),",
                "PAM50_and_Polyak_Cell_Snapshot-gene_lists.xlsx and",
                "Genelists-summary.xlsx (Meta data from Henrik Johansson), and",
                "140219_OSLO2.n325.CNA.byGeneExpressionProbePositons.HKMV_updated.Sample.Ids_MRA_SJ_41samples.txt",
                "(CNV data, from Miriam Ragle Aure)"), 
              protexpr=prot_data, rnaexpr=rna_data, metatumour=meta_tumours,
              metagene=meta_genes, cna=cna) 
# Save the list for later access
save(exprdata,file="tumourExpressionData")
#-----------------------------------------------------------------------------------------------------------------------------
# Common samples
tumours = intersect(intersect(colnames(exprdata$protexpr),
                              colnames(exprdata$rnaexpr)),
                    colnames(exprdata$cna))

# Ensembl to use 
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


# Create gene-centric matrix
gene.centric.cna <- function(fileToLoad) {
  cna = read.table(fileToLoad, header=TRUE, sep="\t") 
  col = colnames(cna)
  cna = cna[,c(col[1:6],tumours)]
  annot = getBM(attributes=c('efg_agilent_sureprint_g3_ge_8x60k','hgnc_symbol','chromosome_name','start_position','end_position','band'),
                filters ='efg_agilent_sureprint_g3_ge_8x60k',values = cna$probes, mart = ensembl_75)
  annot = annot[!annot$hgnc_symbol=="",]
  annot$chromosome_name[annot$chromosome_name=="X"] = "23"
  chrnames = as.character(1:23)
  annot = annot[annot$chromosome_name %in% chrnames,]
  annot = annot[!duplicated(annot$efg_agilent_sureprint_g3_ge_8x60k),]
  idprobe = match(annot$efg_agilent_sureprint_g3_ge_8x60k, cna$probes)
  cna = cna[idprobe,]
  # print(paste("NAs =", sum(is.na(cna))))
  
  # Median of absolute cna probes
  cna = data.frame(GeneName = annot$hgnc_symbol, cna[,-c(1:6)])
  cna <- as.data.frame(cna %>% group_by(GeneName) %>% summarise_all(funs(median)))
  rownames(cna) = cna$GeneName; cna$GeneName = NULL
  cna = as.matrix(cna)
  return(cna)
}

cna_log <- gene.centric.cna(paste0(filefolder, 'Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt'))
cna_tot <- gene.centric.cna(paste0(filefolder, 'Gene_probe_centric.totCN.OSL2_n331_updated.txt'))

# Absolute copy number file  
cna_Absolute_number_Cohort = read.table(paste0(filefolder, "Gene_probe_centric.totCN.OSL2_n331_updated.txt"), header=TRUE, sep="\t")

cna_Absolute_number_Samples = cna_Absolute_number_Cohort[ ,c("probes","gene","chr","start","stop",tumours)]

cna_tot_probes = cna_Absolute_number_Samples[ ,tumours]
# Load ploidy file
ploidy = read.table(paste0(filefolder, "OSL2_n331_ploidy_t.perc.txt"), sep = "\t", header = TRUE)

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

# Binary CNA data
cna_Binary_Samples <- binarize.copy.number(cna_Absolute_number_Samples)
rownames(cna_Binary_Samples) <- cna_Absolute_number_Samples$probes

count.gain.loss <- function(x) {
  gains <- apply(x,1, function(i) sum(i==1))
  losses <- apply(x,1, function(i) sum(i==-1))
  list(gains = gains, losses = losses)
}

gains <- count.gain.loss(cna_Binary_Samples)$gains
losses <- -count.gain.loss(cna_Binary_Samples)$losses

cna_Binary_Samples_Augmented <- cbind(cna_Absolute_number_Samples[,c("probes","gene","chr","start","stop")], gains=gains, losses=losses)


## Get genomic position for probes
# ChrLengths
setClass("num.with.commas")
setAs("character", "num.with.commas", 
      function(from) as.numeric(gsub(",", "", from) ) )

chrlengths=read.table(file=paste0(filefolder, "GRCh37_chromLengths_NCBI.txt"),header=TRUE, sep="\t", 
                      colClasses=c("character","num.with.commas",
                                   "character","character"))[1:24,1:2]
# Cumulative lengths
chrlengths$end=cumsum(chrlengths$total.length)
chrlengths$start=0
chrlengths$start[2:24]=chrlengths$end[1:23]

# Get genomic position for probes
source("chrLength.R")

cna_Binary_Samples_Augmented$pos <- as.integer(cna_Binary_Samples_Augmented$start) + chrlengths$start[as.integer(cna_Binary_Samples_Augmented$chr)]

# PAM50 subtypes
typ = read.table(paste0(filefolder, "Pam50subtype.txt"), header = TRUE, sep = "\t")
typ$Id = gsub("-",".",typ$Id)
setSubtypes = typ$Category[match(tumours,typ$Id)]
names(setSubtypes) <- tumours
molecular_subtype = unique(setSubtypes)


# Common genes (include cna_log also)
rownames(prot)[grepl("\\.", rownames(prot))] = gsub("\\.","-",rownames(prot)[grepl("\\.", rownames(prot))])
genes = intersect(intersect(intersect(rownames(prot),
                                      rownames(mrna)),
                            rownames(cna_tot)),
                  rownames(cna_log))


# Final matrices
cna_tot <- cna_tot[genes, tumours]
cna_log <- cna_log[genes, tumours]


## Impute cna_log -Infs with the minimum values
id_row <- which(apply(cna_log,1,min)=="-Inf")
cna_log[cna_log=="-Inf"] = Inf
cna_log[cna_log=="Inf"] = apply(cna_log[id_row,],1,min)

cna_bin <- binarize.copy.number(cna_tot)
mrna <- mrna[genes,tumours]
prot <- prot[genes,tumours]

# Common genes (include cna_log also)
rownames(exprdata$protexpr)[grepl("\\.", rownames(exprdata$protexpr))] = gsub("\\.","-",rownames(exprdata$protexpr)[grepl("\\.", rownames(exprdata$protexpr))])
genes = intersect(intersect(intersect(rownames(exprdata$protexpr),
                                      rownames(exprdata$rnaexpr)),rownames(cna_tot)),
                  rownames(cna_log))

# Find annotation mismatches
# setdiff(rownames(exprdata$protexpr),rownames(exprdata$rnaexpr))

## Load IntCl genes table
# intCl_Table <- read.xlsx("intCl_Table.xls", sheetName = "ANOVA_ALL_subtypes")


# Final matrices
cna_tot <- cna_tot[genes,tumours]

cna_log <- cna_log[genes,tumours]
# Inf values
# id <- which(apply(cna_log,1,min)=="-Inf")


cna_bin <- binarize.copy.number(cna_tot)

mrna <- exprdata$rnaexpr[genes,tumours]
prot <- exprdata$protexpr[genes,tumours]

# PAM50 colours
PAM50cols <- c('Basal' = '#E31A1C', 'Her2'='#FB9A99','LumA'='#1F78B4','LumB'='#A6CEE3','Normal'='#33A02C')

colour_subtype <- as.character(revalue(setSubtypes, PAM50cols))

# iCluster colours
IntCLcols <- c("1" = '#FF5500',"2" = '#00EE76', "3" = '#CD3278', "4" = '#00C5CD', "5" = '#8B0000',
               "6" = '#FFFF40',"7" = '#0000CD', "8" = '#FFAA00', "9" = '#EE82EE', "10" = '#7D26CD')

setICluster <- meta_data$iCluster[match(tumours, meta_data$OSL2.full.name)]
names(setICluster) <- tumours
colour_IntCl <- as.character(revalue(setICluster, IntCLcols))


## Define ESR status with logistic regression
estr_status <- meta_data[tumours, "ER"]
id_To_Predict <- which(estr_status == "NA")
id_Train <- which(estr_status != "NA")
estr_status <- ifelse(meta_data[tumours, "ER"] == "neg", 0, 1)
p <- prot["ESR1",id_Train]

model <- suppressWarnings(glm(formula = estr_status[id_Train] ~ p, family = binomial))
# summary(model)
predict.glm(model, newdata = data.frame(p=prot["ESR1",id_To_Predict]) ,type = "response")
estr_status[id_To_Predict] <- 1


# Estrogen receptor residual matrices
mrna_res <- compute.residual.matrix (mrna,"ESR_binary")
prot_res <- compute.residual.matrix (prot,"ESR_binary")

# meta_data$iCluster
## Impute cna_lod -Infs with the minimum values
id_row <- which(apply(cna_log,1,min)=="-Inf")
cna_log[cna_log=="-Inf"] = Inf
cna_log[cna_log=="Inf"] = apply(cna_log[id_row,],1,min)
#----------------------------------------------------------------------------------------------------
