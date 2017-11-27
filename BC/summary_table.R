### Create summary table
omicsType1 <- mrna
omicsType2 <- prot
binary_CNA <- cna_bin

data_all <- lapply(c(1,-1), function(i) {
  
  generate.final.table(level = i)})

dat_all_Gains <- data_all[[1]]
dat_all_Losses <- data_all[[2]]
#-------------------------------------------------------------------------------------------------------
common_cols <- intersect(colnames(dat_all_Gains), colnames(dat_all_Losses))

dat_all <- merge(dat_all_Gains, dat_all_Losses, by = common_cols)

res_all_mrna <- make.summary.table(omicsType1, "mRNA")
res_all_prot <- make.summary.table(omicsType2, "Protein")

res_all <- merge(res_all_mrna, res_all_prot, by = "genes")

colsToRemove <- grepl("Gains", colnames(res_all)) | grepl("Losses", colnames(res_all))
res_all <- res_all[,!colsToRemove]

# write.xlsx(res_all, "res_all.xlsx", row.names = FALSE, col.names = TRUE)

final_table <- merge(dat_all,res_all,by="genes")
# write.table(final_table, "final_table_res_ESR.txt", row.names = FALSE, col.names = TRUE, sep="\t")
