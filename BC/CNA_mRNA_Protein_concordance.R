### mRNA, protein expression abundacne relative to ploidy

meltedData <-  lapply(list('Ploidy' = cna_tot, 'CNA' = cna_log, 'mRNA' = mrna, 'Protein' = prot), function(i) {
  res <- melt(i)
  res <- res$value
})

meltedData <- as.data.frame(do.call(cbind, meltedData))

meltedData <- melt(meltedData, id.vars = 'Ploidy')
meltedData$Ploidy <- as.factor(meltedData$Ploidy)

ggplot(data = meltedData, mapping = aes(x = Ploidy, y = value, fill = variable)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(name = 'log2 values', limits = c(-4,5)) + 
  guides(fill = guide_legend(title = '')) + 
  theme(text = element_text(size = 15))
# ---------------------------------------------------------------------------