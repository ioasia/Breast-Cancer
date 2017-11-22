### Create summary table

mrna_res <- compute.residual.matrix(mrna,"ESR")
prot_res <- compute.residual.matrix(prot,"ESR")

omicsType1 <- mrna
omicsType2 <- prot
binary_CNA <- cna_bin
level = -1
levelName <- ifelse(level==1, "Gains", ifelse(level==-1,"Losses","Amplifications")) 


generate.final.table <- function(level) {
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

dat_all_Gains <- generate.final.table(level)
# dat_all_Losses <- generate.final.table(level)
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
write.table(final_table, "final_table_res_ESR.txt", row.names = FALSE, col.names = TRUE, sep="\t")
