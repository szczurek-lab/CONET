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


READS_DIR <- "data/"
RESULTS_DIR <- "./"

all_reads <- read.table(paste0(READS_DIR, "SA501X3F_filtered_corrected_counts.csv"), sep =",", header = T, stringsAsFactors = F)

all_reads$chr[all_reads$chr == "X"] <- "23"
all_reads$chr[all_reads$chr == "Y"] <- "24"
all_reads$chr <- as.numeric(all_reads$chr)


TREE_FILE_NAME <- paste(RESULTS_DIR, "inferred_tree", sep="")
ATTACHMENT_FILE_NAME <- paste(RESULTS_DIR, "inferred_attachment", sep="")

create_edges_matrix <- function() {
  edges <- read.table(TREE_FILE_NAME, sep="-", stringsAsFactors = F)
  edges$hash <- as.character(edges$V1)
  edges$child_hash <- as.character(edges$V2)
  attachment <- read.table(ATTACHMENT_FILE_NAME, sep =";", stringsAsFactors = F)
  attachment$node <- paste("(", as.character(attachment$V2), ",", as.character(attachment$V3), ")", sep ="")
  EDGES <- dim(edges)[1]
  edges$start1 <- rep(0, EDGES)
  edges$start2 <- rep(0, EDGES)
  edges$start3 <- rep(0, EDGES)
  edges$start4 <- rep(0, EDGES)
  edges$chromosome1 <- rep(0, EDGES)
  edges$chromosome2 <- rep(0, EDGES)
  edges$cells <- rep(0, EDGES)  
  return (edges)
}

attachment <- read.table(ATTACHMENT_FILE_NAME, sep=";", stringsAsFactors = F)
attachment$node <- paste("(", as.character(attachment$V2), ",", as.character(attachment$V3), ")", sep ="")
edges <- create_edges_matrix()
EDGES <- dim(edges)[1]

get_start_of_node <- function(chromosome, start) {
  reads_local <- all_reads[all_reads$chr == chromosome & all_reads$start > start, ]
  return (min(reads_local$start))
}

get_end_of_node <- function(chromosome, start) {
  reads_local <- all_reads[all_reads$chr == chromosome & all_reads$start < start, ]
  return (max(reads_local$start))
}

get_number_of_attached_cells <- function(node) {
  return (sum(attachment$node == node))
}



fill_edges_matrix <- function(edges) {
  for (i in 1:EDGES) {
    label <- as.character(edges$child_hash[i])
    event1 <- strsplit(strsplit(label, ",")[[1]][1], "\\(")[[1]][2]
    event2 <- strsplit(strsplit(label, ",")[[1]][2], "\\)")[[1]][1]
    chromosome <- as.integer(strsplit(event1, "_")[[1]][1])
    X1 <- strsplit(event1, "_")[[1]][2]
    X2 <- strsplit(event2, "_")[[1]][2]
    start1 <- as.integer(X1)
    start2 <- get_end_of_node(chromosome, as.integer(X2))
    edges$start3[i] <- start1
    edges$start4[i] <- start2
    edges$chromosome2[i] <- chromosome
    edges$cells[i] <- get_number_of_attached_cells(label)
    
    if (i != 1 &  as.character(edges$hash[i]) != "(0,0)" ) {
      label <- as.character(edges$hash[i])
      event1 <- strsplit(strsplit(label, ",")[[1]][1], "\\(")[[1]][2]
      event2 <- strsplit(strsplit(label, ",")[[1]][2], "\\)")[[1]][1]
      chromosome <- as.integer(strsplit(event1, "_")[[1]][1])
      X1 <- strsplit(event1, "_")[[1]][2]
      X2 <- strsplit(event2, "_")[[1]][2]
      start1 <- as.integer(X1)
      start2 <- get_end_of_node(chromosome, as.integer(X2))
      edges$start1[i] <- start1
      edges$start2[i] <- start2
      edges$chromosome1[i] <- chromosome
  
    }
  }
  return (edges)
}

edges <- fill_edges_matrix(edges)

buildTree <- function(root, edges) {
  for (i in 1:length(edges$hash)[1]) {
    if (edges$hash[i] == root$name) {
      child <- root$AddChild(edges$child_hash[i])
      buildTree(child, edges)
    }
  }
}


GetNodeLabel <- function(node, dict_for_genes=genes_location) {
  if (as.character(node$name) == "(0,0)") {
    return (paste(node$name, " " , as.character(get_number_of_attached_cells(node$name)), sep = ""))
  }
  
  edge <- edges[edges$child_hash == node$name, ]
  return (paste(
    paste(rownames(edges[edges$child_hash == node$name, ]), edge$chromosome2, " [",edge$start3, ",", edge$start4, "] " ),
    paste("number of cells:", as.character(get_number_of_attached_cells(node$name))),
    sep = "\n"))
}


root <- Node$new(edges$hash[1])
buildTree(root, edges = edges)
SetNodeStyle(root, fontname = 'helvetica', fontsize=45,  label = GetNodeLabel)


###
### TREE PLOTTING
###
tree_plt = plot(root)
saveWidget(tree_plt, "tree.html")




all_non_root_nodes <- unique(edges$child_hash)
leafs <- all_non_root_nodes[!(all_non_root_nodes %in%edges$hash)]
post_order <- Traverse(root, traversal = "post-order")
post_order_nodes <- c()
for (i in 1:length(all_non_root_nodes)) {
  post_order_nodes <- append(post_order_nodes, post_order[[i]]$name)
}

#### list of all non-root nodes with list of cells attached
#### NA if no cells attached
#### position on attachment list corresponds to position on all_non_root_nodes list

cell_attachment <- list()
for (n in all_non_root_nodes) {
  cell_id_from_att <- attachment$V1[attachment$node == n]
  if (length(cell_id_from_att)!=0) {
    
    cell_id_from_names <- cell_id_from_att + 1
    names_of_cells <- list(names(all_reads)[cell_id_from_names+ 4])
    cell_attachment <- append(cell_attachment, names_of_cells)
  }
  else {
    cell_attachment <- append(cell_attachment, list(NA))
  }
}
names(cell_attachment) <- all_non_root_nodes

#### list of all non-root nodes by child hash with list of CUMULATIVE cells
#### i.e. cells attached to this node and all it's children
#### i.e. all cells that underwent the event described in node
#### position on attachment list corresponds to position on all_non_root_nodes list

### create empty list
node_and_ancestors_cells_att <- list()
checked_nodes <- c()

### iterate over leafs
for (l in leafs) {
  leaf = l
  checked_nodes <- append(checked_nodes, leaf)
  node_and_ancestors_cells_att <- append(node_and_ancestors_cells_att, cell_attachment[leaf])
  leaf_parent <- edges[edges$child_hash==leaf, "hash"]
  while (leaf_parent != "(0,0)") {
    new_cells <- list(na.omit(unique(unlist(c(cell_attachment[leaf_parent], node_and_ancestors_cells_att[leaf])))))
    names(new_cells) <- leaf_parent

    if (leaf_parent %in% names(node_and_ancestors_cells_att)) {
      more_new_cells <- list(na.omit(unique(unlist(c(new_cells, node_and_ancestors_cells_att[leaf_parent])))))
      node_and_ancestors_cells_att[leaf_parent] <- more_new_cells
    }
    else {
      node_and_ancestors_cells_att <- append(node_and_ancestors_cells_att, new_cells)
    }
    checked_nodes <- append(checked_nodes, leaf_parent)
    leaf <- leaf_parent
    leaf_parent <- edges[edges$child_hash==leaf, "hash"]

  }
}

####################################################################################
####################################################################################

# list of pre ordered nodes by hash

pre_order <- Traverse(root, traversal = "pre-order")
pre_order_nodes <- c()
for (i in 1:length(all_non_root_nodes)) {
  pre_order_nodes <- append(pre_order_nodes, pre_order[[i]]$name)
}
pre_order_nodes <- append(pre_order_nodes, post_order_nodes[!(post_order_nodes%in%pre_order_nodes)])



### create data frame of bins with their an(i) - history of vertices/events
bins_an <- all_reads
bins_an[,5:ncol(bins_an)] <- root$name


### iterate over nodes
for (node in pre_order_nodes[2:length(pre_order_nodes)]) {
  
  start_chromosome <- edges[edges$child_hash==node, "chromosome2"]
  start <- edges[edges$child_hash==node, "start3"]
  end <- edges[edges$child_hash==node, "start4"]
  start_row <- which(all_reads$chr == start_chromosome & all_reads$start == start)
  end_row <- which(all_reads$chr == start_chromosome & all_reads$start == end)
  names_of_cells <- node_and_ancestors_cells_att[[node]]
  
  bins_an[start_row:end_row, names_of_cells] <- sapply(bins_an[start_row:end_row, names_of_cells], function(x) paste(x, node))
}


clusters <- unique(unlist(bins_an[,5:ncol(bins_an)]))
no_of_clusters <- length(clusters[clusters!=root$name])

no_bins_changed_by_tree <- length(which(bins_an[,5:ncol(bins_an)]!=root$name))


### create data frame of inferred CNs with neutral CN
CNs <- as.matrix(all_reads)
CNs[,5:ncol(CNs)] <- 2

all_reads_m <- as.matrix(all_reads)
max_CN <- 10



########################## CNs inference procedure ######################
  
for (c in clusters[clusters!=root$name]) {
  filter <- which(bins_an==c)
  counts <- all_reads_m[filter]
  inffered_CN <- min(round(mean(counts), digits = 0), max_CN)

  CNs[filter] <- inffered_CN
}



CNs <- as.data.frame(CNs)
colnames(CNs) <- colnames(all_reads)
rownames(CNs) <- rownames(all_reads)



############# heatmaps for inferred_counts
############# 
#############

image_qual <- 600

CNs_heatmap <- as.matrix(transpose(CNs[,5:ncol(CNs)]))

rownames(CNs_heatmap) <- colnames(CNs)[5:ncol(CNs)]
colnames(CNs_heatmap) <- CNs$chr

x_labels <- c("1")
x_ticks <- c(1)
for (i in 1:(ncol(CNs_heatmap)-1)) {

  if (colnames(CNs_heatmap)[i] != colnames(CNs_heatmap)[i+1])
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


CCs_heatmap <- as.matrix(transpose(all_reads[,5:ncol(all_reads)]))

#### for correct color plotting need to change all values of 4 and above 
CCs_heatmap[CCs_heatmap>3.9] <- 3.9

rownames(CCs_heatmap) <- colnames(all_reads)[5:ncol(all_reads)]
colnames(CCs_heatmap) <- all_reads$chr


col_breaks_2 <- seq(0, 4, by = 0.2)
n_col <- length(col_breaks_2)-1

my_palette_2 <- c(hcl.colors(n_col/2, palette = "Blues"), 
  hcl.colors(n_col/2, palette = "Reds", rev=-1))


tiff("./heatmap_CCs.tiff", units="in", width=13, height=8, res=image_qual)

heatmap3(CCs_heatmap, col=my_palette_2, breaks=col_breaks_2, Colv = NA, Rowv = as.dendrogram(row_hc), scale = "none", labCol = x_labels,  cexCol = 1.5, lasCol=1, add.expr = abline(v=x_ticks, col="black", lwd=0.2, lty=3))

dev.off()



