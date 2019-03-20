#' Volcan Plot
#'
#' Draws a volcanplot showing diffentially expressed genes. X axis is log2 fold change, y axis is adjusted p.value
#'
#' @param df Dataframe of results
#' @param logFC_thresh Log fold change threshold for highlighting genes
#' @param padj_thresh Adjusted p. value threshold for highlighting genes
#' @param contrast_name Title of the plot
#'
#' @return a ggplot object
#'
#' @export

VolcanoPlot = function(df, logFC_thresh, padj_thresh, contrast_name){
  df = as.data.frame(df)
  mydata<-df%>%mutate(threshold = ifelse(log2FoldChange >= logFC_thresh & padj<padj_thresh,"A", ifelse(log2FoldChange<=-logFC_thresh &padj<padj_thresh, "B", "C")))
  mydata$Gene = as.character(rownames(df))
  lim_logFC = round(max(abs(min(mydata$log2FoldChange)), abs(max(mydata$log2FoldChange)))*1.1)
  lim_padj = round(-log10(min(mydata$padj))*1.1)
  g = ggplot(data=mydata, aes(x=log2FoldChange, y=-log10(padj),label=Gene, colour=threshold)) +
    geom_point(alpha=0.3, size=2) +
    theme(legend.position = "none") +
    xlim(c(-lim_logFC, lim_logFC)) + ylim(c(0, lim_padj)) +
    xlab("log2 fold change") + ylab("-log10 padj") + scale_colour_manual(values = c("A"= "darkgreen", "B"="red",  "C"= "black")) + geom_text_repel(
      data          = subset(mydata, abs(log2FoldChange) > logFC_thresh & padj<padj_thresh),
      segment.size  = 0.2,
      segment.color = "grey50") + geom_hline(yintercept = -log10(padj_thresh), colour='blue') + geom_vline(xintercept = c(logFC_thresh, -logFC_thresh), colour='blue') + ggtitle(paste(contrast_name,'\nVolcano plot for abs(LogFC) >', logFC_thresh,'and padj <',padj_thresh))
  return(g)
  #ggsave(g, filename = paste0('volcano_plot_LogFC_',logFC_thresh, '_padj_',gsub(pattern = '\\.','',padj_thresh), '_',gsub(pattern = ' ', replacement = '_', contrast_name), '.pdf'))
}
