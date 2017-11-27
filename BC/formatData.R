### Setup libraries
install.packages("easypackages")
library(easypackages)

source("http://bioconductor.org/biocLite.R")
toInstallBio <- c("impute", "preprocessCore", "GO.db", "biomaRt","GenomicRanges","CGHregions")
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
                      header=TRUE, 
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


# Binary CNA data
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

chrlengths=read.table(file=paste0(filefolder, "GRCh37_chromLengths_NCBI.txt"),header=TRUE, sep="\t", 
                      colClasses=c("character","num.with.commas",
                                   "character","character"))[1:24,1:2]
# Cumulative lengths
chrlengths$end=cumsum(chrlengths$total.length)
chrlengths$start=0
chrlengths$start[2:24]=chrlengths$end[1:23]


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