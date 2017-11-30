## Data

# Specify filefolder
filefolder <- "/Users/johnsiavelis/Downloads/Breast Cancer/data/"


# Ensembl annotation
ensembl_75 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                     host="feb2014.archive.ensembl.org", 
                     path="/biomart/martservice", 
                     dataset="hsapiens_gene_ensembl")


# mRNA
mrna <- read.table(paste0(filefolder, "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne_updated.txt"),                      header=TRUE, 
                  sep="\t", 
                  stringsAsFactors=FALSE)

# Correct erroneous value jan.27 with 1.3527
mrna$OSL2U.0383T1[mrna$ProbeUID.SampleArray.gDetrendedSignal == "12629"] <- 1.3527
mrna$OSL2U.0383T1 <- as.double(mrna$OSL2U.0383T1)

# Remove duplicate
mrna$OSL2U.0219T1 <- NULL

# Remove superfluous columns
mrna<-mrna[,!names(mrna) %in% c("ProbeUID.SampleArray.gDetrendedSignal","ProbeName","SystematicName","chrom","seq_beg","seq_end","accessions")]

# Replace the rows of individual probes mapping to the same genename by their median 
# row (for each column) and convert to array
mrna <- as.data.frame(mrna %>% group_by(GeneName) %>% summarise_all(funs(median)))
rownames(mrna) <- mrna$GeneName
mrna$GeneName <- NULL
mrna <- as.matrix(mrna)


# Protein data
prot <- readArrayFromExcel(paste0(filefolder, "diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx"), 
                               row.name.col="GeneSymbol")
prot <- log2(prot)

# Correct typos
typos <- c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newcolnames <- colnames(prot)
for(typo in typos){
  print(paste(typo,"T1",sep=""))
  newcolnames<-gsub(typo,paste(typo,"T1",sep=""), newcolnames)
}
colnames(prot) <- newcolnames
colnames(prot) <-  gsub("_",".",colnames(prot))
colnames(prot) <-  gsub("-",".",colnames(prot))


# Copy number data
cna_log <- gene.centric.cna(fileToLoad = paste0(filefolder, 'Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt'))
cna_tot <- gene.centric.cna(fileToLoad = paste0(filefolder, 'Gene_probe_centric.totCN.OSL2_n331_updated.txt'))


# Metadata
meta_data <- as.data.frame(readArrayFromExcel(paste0(filefolder,"clin_info_connect_proteomics.xlsx")))
meta_data$OSL2.full.name <- newcolnames
rownames(meta_data) <- meta_data$OSL2.full.name


# Read in the gene meta data, i.e., a list of PAM50/IntCl genes
pam50 <- readListFromExcel(paste0(filefolder,"PAM50_and_Polyak_Cell_Snapshot-gene_lists.xlsx"))
intCl754 <- unique(read.table(paste0(filefolder,"short_IntCl_genes.txt")))
intCl754 <- as.character(intCl754$V1)
Nik_Zainal <- read.table(paste0(filefolder, "DriverGenes_Nik_Zainal.txt"), sep="\t")
Nik_Zainal <-  unique(as.character(Nik_Zainal$V1))


# Load ploidy file
ploidy <- read.table(paste0(filefolder, "OSL2_n331_ploidy_t.perc.txt"), sep = "\t", header = TRUE)


# Binary CNA data
# Common samples
tumours <- intersect(intersect(colnames(prot),
                              colnames(mrna)),
                    colnames(cna_log))

# Absolute Copy Number file  
cna_Absolute_number_Cohort <- read.table(paste0(filefolder, "Gene_probe_centric.totCN.OSL2_n331_updated.txt"), header=TRUE, sep="\t")
cna_Absolute_number_Samples <- cna_Absolute_number_Cohort[ ,c("probes","gene","chr","start","stop",tumours)]
cna_tot_probes <- cna_Absolute_number_Samples[ ,tumours]


cna_Binary_Samples <- binarize.copy.number(cna_Absolute_number_Samples)
rownames(cna_Binary_Samples) <- cna_Absolute_number_Samples$probes


gains <- count.gain.loss(cna_Binary_Samples)$gains
losses <- -count.gain.loss(cna_Binary_Samples)$losses

cna_Binary_Samples_Augmented <- cbind(cna_Absolute_number_Samples[,c("probes","gene","chr","start","stop")], gains=gains, losses=losses)


## Get genomic position for probes
# ChrLengths
setClass("num.with.commas")
setAs("character", "num.with.commas", 
      function(from) as.numeric(gsub(",", "", from) ) )

chrlengths<-read.table(file=paste0(filefolder, "GRCh37_chromLengths_NCBI.txt"),header=TRUE, sep="\t", 
                      colClasses=c("character","num.with.commas",
                                   "character","character"))[1:24,1:2]
# Cumulative lengths
chrlengths$end<-cumsum(chrlengths$total.length)
chrlengths$start<-0
chrlengths$start[2:24]<-chrlengths$end[1:23]

cna_Binary_Samples_Augmented$pos <- as.integer(cna_Binary_Samples_Augmented$start) + chrlengths$start[as.integer(cna_Binary_Samples_Augmented$chr)]

# PAM50 subtypes
typ <- read.table(paste0(filefolder, "Pam50subtype.txt"), header = TRUE, sep = "\t")
typ$Id <- gsub("-",".",typ$Id)
setSubtypes <- typ$Category[match(tumours,typ$Id)]
names(setSubtypes) <- tumours
molecular_subtype <- unique(setSubtypes)


# Common genes (include cna_log also)
rownames(prot)[grepl("\\.", rownames(prot))] <- gsub("\\.","-",rownames(prot)[grepl("\\.", rownames(prot))])
genes <- intersect(intersect(intersect(rownames(prot),
                                      rownames(mrna)),
                            rownames(cna_tot)),
                  rownames(cna_log))


# Final matrices
cna_tot <- cna_tot[genes, tumours]
cna_log <- cna_log[genes, tumours]


## Impute cna_log -Infs with the minimum values
id_row <- which(apply(cna_log,1,min)=="-Inf")
cna_log[cna_log=="-Inf"] <- Inf
cna_log[cna_log=="Inf"] <- apply(cna_log[id_row,],1,min)

cna_bin <- binarize.copy.number(cna_tot)
mrna <- mrna[genes,tumours]
prot <- prot[genes,tumours]


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
#----------------------------------------------------------------------------------------------------
