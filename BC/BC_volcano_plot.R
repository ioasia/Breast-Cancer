## Volcano plots

p1 <- toPlotVolcano(omicsType1 = mrna, 
                   omicsType2 = prot, 
                   binary_CNA = cna_bin, 
                   level = 1, 
                   molecular_view = 'Protein')

p2 <- toPlotVolcano(omicsType1 = mrna, 
                    omicsType2 = prot, 
                    binary_CNA = cna_bin, 
                    level = -1, 
                    molecular_view = 'Protein')

# pdf(paste0(levelName,".",molecular_view, ".","Volcano_plot",".pdf"), width=20, height=15, paper='special')
plot(p1)
plot(p2)
# multiplot(p1, p2, cols=2)
# dev.off()
#--------------------------------------------------------------------------------------------------




