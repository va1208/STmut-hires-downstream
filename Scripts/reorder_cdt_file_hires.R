#!/usr/bin/env Rscript
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

## load the library for the input argument ###
library(optparse)
option_list <- list(make_option(c("-i", "--cdt_file"), type = "character", help = "STmut CDT file"),
                    make_option(c("-l", "--barlegend"), type = "character", help="Barcode legend obtained from STmut"),
                    make_option(c("-g", "--cluster_info"), type = "character", help = "Cluster info from cloupe file"),
                    make_option(c("-a", "--annotate"), type = "character", help = "annotation file containa two columns-cluster and annotations"),
                    make_option(c("-c", "--order_cluster"), type = "character", help = "txt file with order of cluster, want to arrange, no header"),
                    make_option(c("-o", "--output"), type = "character", help = "outfile name with prefix"))
opt_parser <- OptionParser(option_list = option_list, description = "To reoder the STmut cdt file and create the cdt and heatmap depends on your ordering file supplied")
opt <- parse_args(opt_parser)

args <- c(opt$cdt_file, opt$barlegend, opt$cluster_info, opt$annotate, opt$order_cluster, opt$output)

if (is.null(opt$cdt_file) | is.null(opt$barlegend) | is.null(opt$annotate) | is.null(opt$order_cluster)) {
  print_help(opt_parser)
  stop("All inputs must be supplied", call. = FALSE)
}

### this script will arange the cnv cdt file to make the heatmap ####
## import libraries ###
packages <- c("data.table", "ggplot2", "dplyr")

for(pkg in packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies=TRUE, repos = "https://cloud.r-project.org")
  } else {
    library(pkg, character.only = TRUE)
  }
}
library(ComplexHeatmap)
library(circlize)
print("packages loaded")
gc()
### function to extract the cluster number form the cluster column
extract_cluster_num <- function(x) {
  as.integer(sub(".*?(\\d+)\\s*$", "\\1", x))
}

### set working directory ###
#setwd("~/Library/CloudStorage/OneDrive-UCSF/vignesh/STmut_run/KA03/bins16_out/")
### import the cdt file ###
#cdt <- fread("CNVs_OrganizedByGEcluster_UMIcount.cdt")
cdt <- fread(args[1])
print(cdt[1:5,1:5])
#spot_included <- colnames(cdt)[-c(1,2)]
print("imported cdt file")
### import the barcode legned file ###
legend <- fread(args[2])
print(colnames(legend))
setnames(legend, old = "main_barcode", new = "barcode")
#legend <- fread("BarcodeLegend.csv")
#head(legend)
print("imported barcode legend file")
## import the annotation file ###
annotate <- fread(args[4])
print(colnames(annotate))
annotate$cluster <- gsub(" ", "", gsub("cluster", "Cluster", annotate$cluster))
print("imported annotation file")

### import the cluster info ###
cluster_info <- fread(args[3])
colnames(cluster_info) <- c("barcode", "cluster")
cluster_info$cluster <- gsub(" ", "", cluster_info$cluster)

#### merged cluster and annotation file ####
cluster_info[annotate, annotate := i.annotate, on = .(cluster)]
print(head(cluster_info))
### merge cluster annotation and legend file ###
legend <- left_join(legend, cluster_info, by = "barcode")
print(head(legend))
#legend[annotate, annotate := i.annotate, on = .(cluster)] ### combine the annotation and barlegend file ###
stopifnot(!is.na(legend$annotate))

### get the column names for spot
annot_col  <- c(1L,2L)
annot_name <- names(cdt)[annot_col]
spot_cols  <- setdiff(names(cdt), annot_name)
print(length(spot_cols))
#tie meta to CDT columns, preserving original intra-cluster order
meta <- data.table(barcode = spot_cols)[legend, on="barcode"]
meta$cluster <- extract_cluster_num(meta$cluster) 
print(head(meta))
#### cluster_order ##
cl <- fread(args[5], header = F)
colnames(cl) <- c("cluster", "annotate")
cl$cluster <- extract_cluster_num(cl$cluster)
cl_order <- cl$cluster
tu_order <- cl$annotate
tu_levels <- unique(tu_order)
cl$cluster <- as.factor(cl$cluster)
#cl_order <- c(4,5,6,10,12,13,14,17,8,3,7,9,16,19,20,21,22,23)
remainder <- setdiff(unique(meta$cluster), cl_order)

#### add the levels to the meta to remain the order same ###
final_cluster_levels <- c(cl_order, sort(remainder)) 
meta[, cluster_f := factor(cluster, levels = final_cluster_levels)] ## add levels, so it wont change
meta$cluster <- as.factor(meta$cluster)
meta[cl, annotate2 := i.annotate, on =.(cluster)]
meta[, annotate_f := factor(annotate2, levels = tu_levels)]

### get the meta and expand by barcodes ###
#print(head(meta))
#barcode <- fread(args[2])
meta_long <- meta[, .(barcode2 = unlist(strsplit(merged_barcodes, ","))), by = .(barcode, annotate_f)] ### convert the data into long format ### so we can combine the FDR, CNV score 
meta_long[, barcode2 := trimws(barcode2)]
meta_long <- subset(meta_long, barcode %in% spot_cols)
meta_long <- meta_long[,c("barcode2", "annotate_f")]
setnames(meta_long, old = "barcode2", new = "barcode")
print(head(meta_long))
fwrite(meta_long, paste0(args[6],"_tumorclones_projection.csv"), row.names = F, quote = F, sep = ",")
print("saved the tumor clones projection file")
## subset the meta file and keep the columns in cdt only
meta <- subset(meta, barcode %in% spot_cols)
print(dim(meta))
print(table(meta$cluster_f, useNA = "always"))
print(table(meta$annotate_f, useNA = "always"))
meta[, orig_idx := match(barcode, spot_cols)] ### create index on the meta file #####
setorder(meta, cluster_f, orig_idx) ### set the oder based on the index file ####
ordered_spots <- meta$barcode ### get the spot from the meta and can be used to arrange the cdt file #####
first_two <- names(cdt)[annot_col] ### get the first two columnes (NAME and CLID)
### arrange the based on the order ###
out <- cbind(cdt[,..first_two], cdt[, ..ordered_spots]) ## arrange the cdt file ####
setnames(out, c(annot_name, ordered_spots))
print("ordering completed")
#fwrite(out, paste0(args[6], "_reordered.cdt"), row.names = F, quote = F, sep = "\t")
#print("saved the cdt file and use jave tree to create the plot")

### to generate heatmap #####
out_to_plot <- out[,..ordered_spots] ### export the data wanted to plot
out_to_plot <- as.matrix(out_to_plot) ### convert into matrix 
out_to_plot <- apply(out_to_plot, 2, as.numeric)
print(dim(out_to_plot))
print("data extrcated")

chr_num <- tstrsplit(out[[2]], ":", fixed = TRUE)[[1]] ### get the chr number from the data ###
chr_fac <- factor(as.character(chr_num), levels = as.character(1:23)) 
chr_values <- c("#000000", rep(c("#CBCBCB", "#000000"), 11)) ## add the colors 
names(chr_values) <- as.character(1:23)

### get the cluster factor and its corresponding annotations
clus_labes <- meta$cluster_f
clus_levels <- levels(clus_labes)
print(length(clus_levels))
#tn_labels <- ifelse(tolower(meta$annotate) %in% c("tumor","t"), "T", "N")
tn_labels <- tolower(meta$annotate_f)
cluster_colors <- rep(c("#666666", "#D3D3D3"), length.out = length(clus_levels))
print(length(cluster_colors))
names(cluster_colors) <- clus_levels
#tn_colors <- c(rep("#E41A1C", length(tu_levels)-1), "#66A61E")
#tn_colors <- c(rep("#E41A1C", length(grep("tumor", tu_levels))), rep("#66A61E", length(grep("normal", tu_levels))))
tn_colors <- ifelse(grepl("tumor", tu_levels), "#E41A1C", "#66A61E")
#names(tn_colors) <- c("T", "N")
names(tn_colors) <- tu_levels
print(cluster_colors)
print(tn_colors)
### cluster annotation#####
## chr annotation to get the strips ###
left_ha <- rowAnnotation(
  ChrStrip = anno_simple(as.vector(chr_num), col = chr_values, border = FALSE),
  width = unit(8, "mm"), show_annotation_name = FALSE, annotation_label = FALSE)

#### columns annotation for cluster labels and tumor labels ####
top_ha <- HeatmapAnnotation(
  #TN = factor(tn_labels, levels = c("T","N")),
  col = list(Cluster = cluster_colors, TN = tn_colors),
  TN = factor(tn_labels, levels = tu_levels),
  Cluster = factor(clus_labes, levels = clus_levels),
  show_legend = FALSE, show_annotation_name = FALSE)

### create the heatmap using the annotation ###
hm_col <- colorRamp2(c(-1, 0, 1), c("#37C8FF", "#FFFFFF", "#A50026")) ## give a color the values inside the matrix ####
col_split <- factor(tn_labels, levels = tu_levels)  ### split by tumor 
#col_split2 <- factor(clus_labes, levels = clus_levels) ### split by cluster ###
print("heatmap_starting")
print(class(out_to_plot))
print(sum(is.na(out_to_plot)))
print(which(colSums(is.na(out_to_plot)) > 0))
#pdf(paste0(args[5], "_reordered_heatmap.pdf"), width = 8, height = 10, onefile = F, paper = "a4", useDingbats = TRUE)
#pdf("KA03_CNV_reordered_heatmap.pdf", width = 8, height = 10, onefile = F, paper = "a4", useDingbats = TRUE)
ht <- Heatmap(
  out_to_plot,
  show_heatmap_legend = FALSE,
  top_annotation = top_ha,
  left_annotation = left_ha,
  column_split = col_split,
  column_gap = unit(1, "mm"),
  cluster_rows = FALSE,
  border = FALSE,
  rect_gp = gpar(col = NA),
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = hm_col, use_raster = FALSE)
#draw(ht)
#dev.off()
#print("saved heatmap in pdf format")
png(paste0(args[6], "_reordered_heatmap.png"), w=5000, h=5000, units = "px", res = 290)
#png("KA03_CNV_reordered_heatmap.png", w = 5000, h = 5000, units = "px", res = 400)
draw(ht)
dev.off()
#pdf(paste0(args[6], "_reordered_heatmap.pdf"))
#draw(ht)
#dev.off()
print("saved heatmap in png and pdf format")
