args <- commandArgs(trailingOnly = TRUE)
print(args)

### import packages ###
library(data.table)
library(ggplot2)

## define the args ###
##args[1] - cdt file from the STmut output
##args[2] - clustering file used to run STmut ###
##args[3] - output folder ##


### import cdt ###
cdt <- fread(args[1], header = TRUE)
cdt <- as.data.frame(cdt)

### import clusters ###
clusters <- fread(args[2], header = TRUE)
colnames(clusters) <- c("Barcode", "Cluster")
clusters$Cluster <- gsub(" ", "", clusters$Cluster)

### import annotation ###
#annotate <- fread("annotate_v3.csv")
#annotate$cluster <- gsub("cluster ", "Cluster", annotate$cluster)
#setnames(annotate, "cluster", "Cluster")


## combine cluster and annotation ####
#clusters[annotate, annotate := i.annotate, on = "Cluster"]
#clusters <- subset(clusters, annotate == "tumor")

### find common bins ###
common_bins <- intersect(colnames(cdt), clusters$Barcode)
length(common_bins)
### subset the cdt file ####
cdt <- cdt[, common_bins]

### match with cluster file 
clusters <- as.data.frame(clusters)
clusters <- clusters[,c("Barcode", "Cluster")]
clusters <- clusters[match(common_bins, clusters$Barcode),]

### calculate the average cnv for each cluster ####
cluster_profile <- sapply(
    unique(clusters$Cluster),
    function(cl) {
        cells <- clusters$Barcode[clusters$Cluster == cl]
        rowMeans(cdt[ ,cells, drop = FALSE])
    }
)
colnames(cluster_profile) <- gsub(" ", "", colnames(cluster_profile))

cluster_matrix <- t(cluster_profile)

d <- dist(cluster_matrix, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

### prepare the output folder ###
if (!dir.exists(args[3])) {
  dir.create(args[3], recursive = TRUE)
}

setwd(args[3])
#### draw the dentogram ## and save it in the output folder ###
pdf("cluster_dendrogram.pdf", width = 8, height = 6)
plot(hc, xlab = "", sub = "")
dev.off()
message("Dendrogram saved in the output folder.")

### save the matrix, distance and hc object in the output folder ###
files <- list(cluster_matrix = cluster_matrix, distance_matrix = d, hc = hc)
saveRDS(files, file = "cluster_analysis_objects.RDS")
message("Cluster analysis objects saved in the output folder."
