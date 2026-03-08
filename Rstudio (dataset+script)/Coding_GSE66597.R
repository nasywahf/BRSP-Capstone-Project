#Analisis Ekspresi Gen sel astrocyte terinfeksi H5N1 
#Dataset: GSE66597 (Infeksi H5N1 vs Normal)
#Platform: Microarray (Agilent GPL6480)

###############################
# PART A. PERSIAPAN
###############################

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GEOquery","limma"), ask=FALSE, update=FALSE)

install.packages(c("pheatmap","ggplot2","dplyr"))

if (!requireNamespace("umap", quietly=TRUE)){
  install.packages("umap")
}

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(umap)

###############################
# PART B. DOWNLOAD DATA GEO
###############################

gset <- getGEO("GSE66597", GSEMatrix=TRUE, AnnotGPL=TRUE)[[1]]

###############################
# PART C. PREPROCESSING
###############################

ex <- exprs(gset)

qx <- as.numeric(quantile(ex,
                          c(0,0.25,0.5,0.75,0.99,1),
                          na.rm=TRUE))

LogTransform <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)

#Jika True
if (LogTransform){
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

exprs(gset) <- ex

###############################
# PART D. METADATA SAMPEL
###############################

pdata <- pData(gset)

head(pdata[,c("title","source_name_ch1")])
#Digunakan untuk memisahkan data berdasarkan perlakuan infeksi (Control vs infeksius)
group_info <- ifelse(
  grepl("control", pdata$source_name_ch1, ignore.case=TRUE),
  "control",
  "infection"
)

gset$group <- factor(group_info)

table(gset$group)

###############################
# PART E. MEMBERSIHKAN DATA
###############################

ex <- exprs(gset)

ex <- ex[rowSums(is.na(ex)) == 0,]

ex <- ex[apply(ex,1,var) != 0,]

###############################
# PART F. DESIGN MATRIX
###############################

groups <- factor(gset$group)

design <- model.matrix(~0 + groups)

colnames(design) <- levels(groups)

grup_infeksi <- levels(groups)[2]
grup_normal <- levels(groups)[1]

contrast_formula <- paste(grup_infeksi,"-",grup_normal)

print(paste("Kontras:",contrast_formula))

###############################
# PART G. ANALISIS LIMMA
###############################

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust="fdr",
  sort.by="B",
  number=Inf
)

head(topTableResults)

###############################
# PART H. ANOTASI GEN
###############################

feature_data <- fData(gset)

colnames(feature_data)

symbol_col <- grep("symbol",
                   colnames(feature_data),
                   ignore.case=TRUE,
                   value=TRUE)[1]

title_col <- grep("title",
                  colnames(feature_data),
                  ignore.case=TRUE,
                  value=TRUE)[1]

gene_annotation <- feature_data[,c("ID",symbol_col,title_col)]

colnames(gene_annotation) <- c("PROBEID",
                               "Gene.symbol",
                               "Gene.title")

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by="PROBEID",
  all.x=TRUE
)

head(topTableResults)

###############################
# PART I. BOXPLOT
###############################

group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col=group_colors,
  las=2,
  outline=FALSE,
  main="Boxplot Ekspresi Distribusi Nilai Ekspresi per Sampel",
  ylab="Expression (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

###############################
# PART J. DENSITY PLOT
###############################

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() + 
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

###############################
# PART K. UMAP
###############################
#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#agar plot konsisten (opsional)
set.seed(123)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1=umap_result$layout[,1],
  UMAP2=umap_result$layout[,2],
  Group=gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )
###############################
# PART L. VOLCANO PLOT
###############################


volcano_data <- data.frame(
  logFC=topTableResults$logFC,
  adj.P.Val=topTableResults$adj.P.Val,
  Gene=topTableResults$Gene.symbol
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG H5N1")

###############################
# PART M. HEATMAP
###############################
#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),]

top50 <- head(topTableResults,50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID,]

gene_label <- ifelse(
  is.na(top50$Gene.symbol) |
    top50$Gene.symbol=="",
  top50$PROBEID,
  top50$Gene.symbol
)

rownames(mat_heatmap) <- gene_label

annotation_col <- data.frame(
  Group=gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

###############################
# PART N. SIMPAN HASIL
###############################

write.csv(topTableResults,
          "Hasil_GSE66597_DEG.csv")

write.csv(top50,
          "Top50_GSE66597.csv")

message("Analisis selesai.")
