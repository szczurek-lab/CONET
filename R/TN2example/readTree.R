# install and load packages
library(RColorBrewer)
library(ggplot2)
library(data.tree)
library(data.table)
library(htmlwidgets)
library(webshot)
library(DiagrammeR)
library(stringr)
library(DescTools)
library(PerformanceAnalytics)
library(heatmap3)
library(cowplot)

############################################################################################
#######################    START of user input    #######################
############################################################################################

################## set directory and import CONET functions ###################
path <- getwd()
################## import CONET functions ###################
source(paste0(path, "/functionsCONET.R"))


#######################    cancer genes data    #######################
genes_location <- readRDS(paste0(path, "/cosmic_cancer_genes.RDS"))
available_cancer_types <- readRDS(paste0(path, "/available_cancer_types.RDS"))

# print available cancer types
sort(available_cancer_types)

#### USER can set cancer type of cancer from the ones printed above
# list of character cancer types
# or "no" - all cosmic cancer genes will be printed
cancer_type_input <- c("breast", "breast carcinoma", "lobular breast", "luminal A breast", "secretory breast")


#######################  read CONET data  #######################

## read in per-bin data input
## here we read in a subsample of 50 cells from TN2 breast cancer sample sequenced using ACT
corrected_counts_file <- "/TN2_corrected_counts_with_indices_50cells.csv"
all_reads <- read.table(paste0(path, corrected_counts_file), sep =",", header = T, stringsAsFactors = F)

## give bin width - will be considered when checking which genes overlap with CN events
bin_width_input <- 220000

## set basal ploidy
neutral_CN <- 3

## set maximum inferred integer CN
## if higher than 2*neutral_CN the colors in heatmaps may not reflect higher CNs
maximum_inferred_CN <- 2*neutral_CN

## set id of the results iteration (is usefull when running CONET multiple times)
## can be left empty ""
VERSION <- "iteration_1"

# set model output file names
# in the example we have output from our CONET inferred for TN2 sample
# please remember that the tree looks differently from publication 
# because we only use subsample for illustration purposes
TREE_FILE_NAME <- "inferred_tree"
ATTACHMENT_FILE_NAME <- "inferred_attachment"

############################################################################################
#######################    END of user input    #######################
############################################################################################



# create attachment matrix
attachment <- create_attachment_matrix(attachment_file = paste0(path, "/", ATTACHMENT_FILE_NAME))
attachment[1:10,]

# create edges matrix
edges <- create_edges_matrix(tree_file=paste0(path, "/", TREE_FILE_NAME), attachment_file=paste0(path, "/", ATTACHMENT_FILE_NAME))
# this one may take some time for the whole sample
edges <- fill_edges_matrix(edges)
edges[1:10,]

# build tree object 
tree <- build_tree(root_name=edges$hash[1], edges_matrix=edges)

# tree plotting - a little time consuming
# saving the tree plot in a chosen version as a pdf
# in node label we have:
# chromosome strat breakpoint locus and end breakpoint locus
# breast cancer genes if they overlap with above genomic locations
# number of cells attached to a given node
plot_tree(tree_object = tree, output_file=paste0(path, "/tree_", VERSION, "_with_genes.pdf"))


##### save tree to Newick format
# change separator between start_breakpoint and end_breakpoint
sep <- "__"
edgesNewick <- edges
for (i in 1:nrow(edgesNewick)) {
  edgesNewick$child_hash[i] <- gsub(",", sep, edgesNewick$child_hash[i])
  edgesNewick$hash[i] <- gsub(",", sep, edgesNewick$hash[i])
}

# create data.tree tree object
rootNewick <- Node$new(edgesNewick$hash[1])
buildTree(rootNewick, edges = edgesNewick)

# switch to Newick
newick <- ToNewick(rootNewick, heightAttribute = NULL)
write(newick, paste0(TREE_FILE_NAME, "_Newick_", VERSION))


####################################################################################
####################################################################################

#### infer CN matrix and additional statistics
#### the most time consuming function
results <- CN_inference(tree_object=tree, cc_input=all_reads, edges_matrix=edges, attachment=attachment, only_CNs_for_heatmaps=F, main_CN = neutral_CN, max_CN = maximum_inferred_CN)


# save inferred CN matrix
CNs <- results[[1]]
CNs[1:10,1:10]
write.csv(CNs, paste0(path,"/inferred_CN_matrix_", VERSION, ".csv"))

# save calculated quality measures (described in Additional File 1: Section S6.1)
# is calculated only if results were calculated with only_CNs_for_heatmaps=F
stats <- results[[3]]
write.csv(stats, paste0(path,"/CONET_quality_measures_", VERSION, ".csv"))

############# heatmaps for inferred_counts
############# 
############# important set resolution


#############
############# please note that heatmap color scales work automatically only if
############# neutral_CN in {2, 3, 4} and maximum_inferred_CN <- 2*neutral_CN


image_qual <- 600

CNs_heatmap_input <- prepare_CNs_heatmap_input(CNs_data = CNs, main_CN = neutral_CN) 
CNs_heatmap <- CNs_heatmap_input[[1]]

### outside clustering for faster heatmap
row_hc = hclust(dist(CNs_heatmap))

tiff(paste0(path, "/heatmap_CNs_", VERSION,".tiff"), units="in", width=13, height=8, res=image_qual)
heatmap3(CNs_heatmap, col=CNs_heatmap_input[[5]], breaks=CNs_heatmap_input[[4]], Colv = NA, Rowv = as.dendrogram(row_hc), scale = "none", labCol = CNs_heatmap_input[[2]], cexCol = 1.5, lasCol=1, add.expr = abline(v=CNs_heatmap_input[[3]], col="black", lwd=0.2, lty=3),
  legendfun=function() showLegend(legend=c(0:(2*neutral_CN)),col=CNs_heatmap_input[[5]]))
dev.off()

############# heatmaps for corrected_counts

CCs_heatmap_input <- prepare_CCs_heatmap_input(CCs_data = all_reads, main_CN = neutral_CN)
CCs_heatmap <- CCs_heatmap_input[[1]]

  #!!!! laeve the same clustering as for CNs

tiff(paste0(path, "/heatmap_CCs_", VERSION, ".tiff"), units="in", width=13, height=8, res=image_qual)

heatmap3(CCs_heatmap, col=CCs_heatmap_input[[3]], breaks=CCs_heatmap_input[[2]], Colv = NA, Rowv = as.dendrogram(row_hc), scale = "none", labCol = CNs_heatmap_input[[2]],  cexCol = 1.5, lasCol=1, add.expr = abline(v=CNs_heatmap_input[[3]], col="black", lwd=0.2, lty=3))
dev.off()



