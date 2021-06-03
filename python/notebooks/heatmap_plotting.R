library(RColorBrewer)
library(ggplot2)
library(gplots)
library(data.tree)
library(data.table)
library(htmlwidgets)
library(webshot)
library(DiagrammeR)
library(stringr)
library(heatmap3)


counts <- read.table("counts", header=F, sep=";")
ccs <- read.table("ccs", header=F, sep=";")
diffs <- read.table("ratios", header=F, sep=";")
############# heatmaps for inferred_counts
############# 
#############

image_qual <- 200

CNs_heatmap <- as.matrix(counts)
CNs_heatmap[CNs_heatmap>3.9] <- 4

x_labels <- c("1")
x_ticks <- c(1)
for (i in 1:(ncol(CNs_heatmap)-1)) {

  if (FALSE)
  {
    x_labels <- append(x_labels, colnames(CNs_heatmap)[i+1])
    x_ticks <- append(x_ticks, i+1)
  }
  else {
    x_labels <- append(x_labels, NA)

    
  }
}

x_labels[which(x_labels==23)] <-"X"
x_labels[which(x_labels==24)] <-"Y"


### outside clustering for faster heatmap
row_hc = hclust(dist(CNs_heatmap))

col_breaks <- seq(-0.5, 4.5, by = 1)
my_palette <- c("blue1", "#8FB7D9", "white", "#FF8870","red")
tiff("./heatmap_CNs.tiff", units="in", width=13, height=8, res=image_qual)
heatmap3(CNs_heatmap, col=my_palette, breaks=col_breaks, Colv = NA, Rowv = as.dendrogram(row_hc), scale = "none", labCol = x_labels, cexCol = 1.5, lasCol=1, add.expr = abline(v=x_ticks, col="black", lwd=0.2, lty=3))


dev.off()


############# heatmaps for corrected_counts


CCs_heatmap <- as.matrix(ccs)

#### for correct color plotting need to change all values of 4 and above 
CCs_heatmap[CCs_heatmap>3.9] <- 3.9




col_breaks_2 <- seq(0, 4, by = 0.2)
n_col <- length(col_breaks_2)-1

my_palette_2 <- c(hcl.colors(n_col/2, palette = "Blues"), 
  hcl.colors(n_col/2, palette = "Reds", rev=-1))


tiff("./heatmap_CCs.tiff", units="in", width=13, height=8, res=image_qual)

heatmap3(CCs_heatmap, col=my_palette_2, breaks=col_breaks_2, Colv = NA, Rowv = as.dendrogram(row_hc), scale = "none", labCol = x_labels,  cexCol = 1.5, lasCol=1, add.expr = abline(v=x_ticks, col="black", lwd=0.2, lty=3))

dev.off()


hist(unlist(ccs[ ccs < 11]))


hist(unlist(diffs))
