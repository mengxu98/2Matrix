if (!require("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}
if (!require("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", version = "1.32")
}
if (!require("umap", quietly = TRUE)) {
  install.packages("umap")
}
library(GEOquery)
library(limma)

##### GSE184097
geo_set <- getGEO("GSE184097")
if (length(geo_set) > 1) idx <- grep(geo_set[[1]]@annotation, attr(geo_set, "names")) else idx <- 1
geo_set <- geo_set[[idx]]
geo_set_info <- geo_set@featureData@data
expression_matrix <- as.data.frame(geo_set@assayData$exprs)
expression_matrix$ID <- rownames(expression_matrix)
expression_matrix <- merge(expression_matrix,
                           geo_set_info,
                           by = c("ID"),
                           all.x = TRUE)
rownames(expression_matrix) <- expression_matrix$circRNA
expression_matrix_GSE184097 <- expression_matrix[, 2:9]



##### GSE182192
geo_set <- getGEO("GSE182192")
if (length(geo_set) > 1) idx <- grep(geo_set[[1]]@annotation, attr(geo_set, "names")) else idx <- 1
geo_set <- geo_set[[idx]]
geo_set_info <- geo_set@featureData@data
expression_matrix <- as.data.frame(geo_set@assayData$exprs)
expression_matrix$ID <- rownames(expression_matrix)
expression_matrix <- merge(expression_matrix,
                           geo_set_info,
                           by = c("ID"),
                           all.x = TRUE)
expression_matrix <- expression_matrix[!duplicated(expression_matrix$`GENE SYMBOL`), ]
rownames(expression_matrix) <- expression_matrix$`GENE SYMBOL`
expression_matrix_GSE182192 <- expression_matrix[, 2:7]

#####
de_function <- function(data,
                        num,
                        GEO) {
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!require("ggrepel", quietly = TRUE)) {
    install.packages("ggrepel")
  }
  library(limma)
  library(Biobase)
  list <- c(rep("Treat", num), rep("CK", num)) %>% factor(., levels = c("CK", "Treat"), ordered = F)
  
  list <- model.matrix(~factor(list)+0)
  colnames(list) <- c("CK", "Treat")
  expression_data <- data
  df.fit <- lmFit(expression_data, list)
  df.matrix <- makeContrasts(Treat - CK , levels = list)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit, n = Inf, adjust = "fdr")
  
  foldChange <- 1
  padjThreshold <- 0.05
  tempOutput$reg <- "not DE"
  tempOutput$reg[which(tempOutput$logFC >= 0.5)] <- "min-up-regulated"
  tempOutput$reg[which(tempOutput$logFC <= -0.5)] <- "min-down-regulated"
  tempOutput$reg[which(tempOutput$logFC >= foldChange & tempOutput$adj.P.Val <= padjThreshold)] <- "up-regulated"
  tempOutput$reg[which(tempOutput$logFC <= -foldChange & tempOutput$adj.P.Val <= padjThreshold)] <- "down-regulated"
  
  dataPlot <- as.data.frame(tempOutput)
  dataPlot$change <- ifelse(dataPlot$P.Value < padjThreshold & abs(dataPlot$logFC) >= foldChange,
                            ifelse(dataPlot$logFC > foldChange,
                                   "Up (log2FC > 1, P-value < 0.05)",
                                   "Down (log2FC < -1, P-value < 0.05)"),
                            "no sig")
  
  dataPlot$gene <- rownames(dataPlot)
  dataPlot$label <- ifelse(dataPlot$P.Value < padjThreshold & abs(dataPlot$logFC) >= 1,
                           as.character(dataPlot$gene), "")
  
  library(ggplot2)
  library(ggrepel)
  p <- ggplot(dataPlot, aes(x = logFC, y = -log10(P.Value), colour = change)) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c("#006699", "#d2dae2", "#990033")) +
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(padjThreshold), lty = 4, col = "black", lwd = 0.8) +
    labs(x = "logFC", y = "-Log10(P-value)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank()) +
    geom_text_repel(data = dataPlot, 
                    aes(x = logFC, y = -log10(P.Value),label = label),
                    size = 3,
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),
                    segment.color = "black",
                    show.legend = FALSE)
  
  write.csv(tempOutput,
            file = paste0("DE_results_", GEO, ".csv"))
  
  ggsave(p, file = paste0(GEO, ".png"), height = 4.5, width = 7)
}

de_function(expression_matrix_GSE182192, 3, "GSE182192")
de_function(expression_matrix_GSE184097, 4, "GSE184097")
