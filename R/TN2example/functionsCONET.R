


#' function getting genes in node
#' @param label_child_hash name of the node from CONET edges$child_hash matrix
#' @param genes_dict dictionary with genes names, coordinates and cancer types: cosmic_cancer_genes.csv
#' @param bin_width average no of bases in one corrected counts bin, if the bins differ much in width, the function will not work properly
#' @param cancer_type the cancer of the data set - only genes characteristic to this cancer type will be printed in the tree nodes
get_gene_in_node <- function(label_child_hash, genes_dict, bin_width, cancer_type="no") {
  genes_in_node <- c()
  sub_genes <- genes_dict[which(genes_dict$chr==edges[edges$child_hash==label_child_hash,"chromosome2"]),]
  
  if (nrow(sub_genes)==0){
    return(genes_in_node) 
  }
  
  for (j in 1:nrow(sub_genes)) {
    if (
      as.integer(sub_genes[j,"start"]) > edges[edges$child_hash==label_child_hash,"start3"] && as.integer(sub_genes[j,"start"]) < edges[edges$child_hash==label_child_hash,"start4"]+bin_width
      |
        as.integer(sub_genes[j,"end"]) > edges[edges$child_hash==label_child_hash,"start3"] && as.integer(sub_genes[j,"end"]) < edges[edges$child_hash==label_child_hash,"start4"]+bin_width
    ) {
      
      if (cancer_type=="no") {
        genes_in_node <- append(genes_in_node, rownames(sub_genes[j,]))
      }
      
      else if (sum(cancer_type %in% strsplit(sub_genes[j,"cancer_type"], ", ")[[1]])) {
        genes_in_node <- append(genes_in_node, rownames(sub_genes[j,]))
      }
      
    }
    
    
  }
  return(genes_in_node)
  
}


###### function creating edges matrix
create_edges_matrix <- function(tree_file, attachment_file) {
  edges <- read.table(tree_file, sep="-", stringsAsFactors = F)
  edges$hash <- as.character(edges$V1)
  edges$child_hash <- as.character(edges$V2)
  attachment <- read.table(attachment_file, sep =";", stringsAsFactors = F)
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


##### node functions
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



#### function filling the edges matrix with tree data
fill_edges_matrix <- function(edges) {
  for (i in 1:nrow(edges)) {
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



##### functions for tree printing
buildTree <- function(root, edges) {
  for (i in 1:length(edges$hash)[1]) {
    #print(i)
    if (edges$hash[i] == root$name) {
      child <- root$AddChild(edges$child_hash[i])
      #print(root)
      buildTree(child, edges)
    }
  }
}


GetNodeLabel <- function(node, genes_dict=genes_location, cancer_type=cancer_type_input, bin_width=bin_width_input) {
  if (as.character(node$name) == "(0,0)") {
    return (paste(node$name, " " , as.character(get_number_of_attached_cells(node$name)), sep = ""))
  }
  
  edge <- edges[edges$child_hash == node$name, ]
  end_for_label <- edge$start4 + bin_width
  return (paste(
    paste(edge$chromosome2, " [",edge$start3, ",", end_for_label, "] " ),
    paste(get_gene_in_node(node$name, genes_dict = genes_dict, bin_width=bin_width, cancer_type = cancer_type), collapse=" "),
    paste("number of cells:", as.character(get_number_of_attached_cells(node$name))),
    sep = "\n"))
}



GetEdgeLabel <- function(node) {
  confidence <- max(edge_confidence[as.character(edge_confidence$V2) == node$name, 3])
  return (as.character(confidence))
}



##### reading attachment info
create_attachment_matrix <- function(attachment_file) {

attachment <- read.table(attachment_file, sep =";", stringsAsFactors = F)
attachment$node <- paste("(", as.character(attachment$V3), ",", as.character(attachment$V4), ")", sep ="")
return(attachment)
}

##### building a tree object
build_tree <- function(root_name, edges_matrix) {
  root <- Node$new(root_name)
  buildTree(root, edges = edges_matrix)
  return(root)
}

##### plot and save tree
plot_tree <- function(tree_object, output_file) {
  
  # fontname = 'helvetica', fontsize=45, label = GetNodeLabel,
  SetNodeStyle(tree_object,  label = GetNodeLabel, fontcolor = "black", fontname = 'helvetica', fontsize=45)
  
  #style = "filled", fillcolor = GetNodeColor
  #SetEdgeStyle(root, fontname = 'helvetica', fontsize=45, label = GetEdgeLabel)
  
  #RYSOWANIE DRZEWA
  # plotting and saving the plot
  x = plot(tree_object)
  saveWidget(x, "temp.html")
  webshot("temp.html", output_file)
  #, vwidth = 1591, vheight = 1451
}


#### function creating attachment dictionary
#### list of all non-root nodes with list of cells attached
#### NA if no cells attached
#### position on attachment list corresponds to position on all_non_root_nodes list

create_node_cell_dict <- function(list_of_nodes, cell_node_attachment) {
  
  cell_attachment <- list()
  for (n in list_of_nodes) {
    names_of_cells <- list(cell_node_attachment$V1[which(cell_node_attachment$node == n)])
    if (length(names_of_cells)!=0) {
      cell_attachment <- append(cell_attachment, names_of_cells)
    }
    else {
      cell_attachment <- append(cell_attachment, list(NA))
    }
  }
  
  names(cell_attachment) <- list_of_nodes
  
  return(cell_attachment)
}


#################


#### function prepearing list of all non-root nodes by child hash with list of CUMULATIVE cells
#### i.e. cells attached to this node and all it's children
#### i.e. all cells that underwent the event described in a node
#### position on attachment list corresponds to position on all_non_root_nodes list

create_node_ancestor_cell_dict <- function(node_cell_dict, list_of_leafs, edges_matrix){
  ### create empty list
  node_and_ancestors_cells_att <- list()
  checked_nodes <- c()
  
  ### iterate over leafs
  for (l in list_of_leafs) {
    leaf = l
    checked_nodes <- append(checked_nodes, leaf)
    node_and_ancestors_cells_att <- append(node_and_ancestors_cells_att, node_cell_dict[leaf])
    leaf_parent <- edges_matrix[edges_matrix$child_hash==leaf, "hash"]
    while (leaf_parent != "(0,0)") {
      new_cells <- list(na.omit(unique(unlist(c(node_cell_dict[leaf_parent], node_and_ancestors_cells_att[leaf])))))
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
      leaf_parent <- edges_matrix[edges_matrix$child_hash==leaf, "hash"]
      
    }
  }
  return(node_and_ancestors_cells_att)
}


####################################################################################
####################################################################################

##### creating list of pre ordered nodes by hash
create_pre_order_nodes_list <- function(tree_object, list_of_nodes, list_of_leafs) {
  
  pre_order <- Traverse(tree_object, traversal = "pre-order")
  pre_order_nodes <- c()
  for (i in 1:length(list_of_nodes)) {
    pre_order_nodes <- append(pre_order_nodes, pre_order[[i]]$name)
  }
  pre_order_nodes <- append(pre_order_nodes, list_of_leafs[!(list_of_leafs%in%pre_order_nodes)])
  
  return(pre_order_nodes)
  
}

##### list of in ordered nodes by hash
create_in_order_nodes_list <- function(tree_object, list_of_nodes) {
  
  in_order <- Traverse(tree_object, traversal = "level")
  
  in_order_nodes <- c()
  for (i in 1:length(list_of_nodes)+1) {
    in_order_nodes <- append(in_order_nodes, in_order[[i]]$name)
  }
  
  return(in_order_nodes)
}


CN_inference <- function(tree_object, cc_input, edges_matrix, attachment, only_CNs_for_heatmaps=F, main_CN, max_CN, use_median=T) {
  
  # list of all non-root nodes and genes in them
  all_non_root_nodes <- unique(edges_matrix$child_hash)
  
  # list of leafs
  leafs <- all_non_root_nodes[!(all_non_root_nodes %in%edges_matrix$hash)]
  
  # dictionary of node -- attached cells
  cell_attachment <- create_node_cell_dict(list_of_nodes = all_non_root_nodes, cell_node_attachment = attachment)

  
  #### list of all non-root nodes by child hash with list of CUMULATIVE cells
  #### i.e. cells attached to this node and all it's children
  #### i.e. all cells that underwent the event described in node
  #### position on attachment list corresponds to position on all_non_root_nodes list
  
  ### create empty list
  att_dict <- create_node_ancestor_cell_dict(node_cell_dict = cell_attachment, list_of_leafs = leafs, edges_matrix = edges_matrix)

  
  # list of pre ordered nodes by hash
  pre_order <- Traverse(tree_object, traversal = "pre-order")
  pre_order_nodes <- c()
  for (i in 1:length(all_non_root_nodes)) {
    pre_order_nodes <- append(pre_order_nodes, pre_order[[i]]$name)
  }
  pre_order_nodes <- append(pre_order_nodes, leafs[!(leafs%in%pre_order_nodes)])
  
  ### create data frame of bins with their an(i) - history of vertices/events
  bins_an <- cc_input
  bins_an[,6:ncol(bins_an)] <- tree_object$name
  
  ### iterate over nodes
  for (node in pre_order_nodes[2:length(pre_order_nodes)]) {
    
    start_chromosome <- edges_matrix[edges_matrix$child_hash==node, "chromosome2"]
    start <- edges_matrix[edges_matrix$child_hash==node, "start3"]
    end <- edges_matrix[edges_matrix$child_hash==node, "start4"]
    
    start_row <- which(cc_input$chr == start_chromosome & cc_input$start == start)
    end_row <- which(cc_input$chr == start_chromosome & cc_input$start == end)
    
    names_of_cells <- att_dict[[node]]
    
    bins_an[start_row:end_row, names_of_cells] <- sapply(bins_an[start_row:end_row, names_of_cells], function(x) paste(x, node))
    
  }
  
  
  clusters <- unique(unlist(bins_an[,6:ncol(bins_an)]))
  length(clusters)
  no_of_clusters <- length(clusters[clusters!=tree_object$name])

  
  no_bins_changed_by_tree <- length(which(bins_an[,6:ncol(bins_an)]!=tree_object$name))
  
  
  ### create data frame of inferred CNs with the size of counts.csv and neutral CN
  CNs <- as.matrix(cc_input)
  CNs[,6:ncol(CNs)] <- main_CN
  
  ### create data frame of inferred CNs with the size of counts.csv and neutral CN
  ### BUT for checking when we infer CN == main_CN (will be then changing CN to -1 and coloring green)
  CNs2 <- as.matrix(cc_input)
  CNs2[,6:ncol(CNs2)] <- main_CN

  all_reads_m <- as.matrix(cc_input)


  GI <- c() 
  Ent <- c()
  support <- c()
  cluster_CNs <- c()
  EDGES <- nrow(edges_matrix)
  
  
  if (only_CNs_for_heatmaps) {
    
    for (c in clusters[clusters!=tree_object$name]) {
      
      filter <- which(bins_an==c)
      counts <- all_reads_m[filter]
      
      ### infer CN
      if (use_median){
        inffered_CN <- min(round(Median(counts), digits = 0), max_CN)
      }
      else {
        inffered_CN <- min(round(mean(counts), digits = 0), max_CN)
      }
      
      CNs[filter] <- inffered_CN
      if (inffered_CN==main_CN) {
        CNs2[filter] <- -1
      }
      else { 
        CNs2[filter] <- inffered_CN
      }
    }
  }
  
  else {
    
    c <-  tree_object$name
    filter <- which(bins_an==c)
    counts <- all_reads_m[filter]
    rounded_counts_hist <- hist(counts, breaks = c(seq(-0.5, (max_CN-0.5), 1), 50), plot=F)
    res_hist <- rounded_counts_hist$counts
    ENT_2 <- Entropy(res_hist, base = exp(1))/log(length(res_hist))
    GI_2 <- Gini(counts)
    
    #### infer CN
    if (use_median){
      inffered_CN <- min(round(Median(counts), digits = 0), max_CN)
    }
    else {
      inffered_CN <- min(round(mean(counts), digits = 0), max_CN)
    }
    
    support_2 <- res_hist[3]/length(counts)
    
    for (c in clusters[clusters!=tree_object$name]) {
      
      filter <- which(bins_an==c)
      counts <- all_reads_m[filter]
      no_of_bins <- length(counts)
      rounded_counts_hist <- hist(counts, breaks = c(seq(-0.5, (max_CN-0.5), 1), 50), plot=F)
      res_hist <- rounded_counts_hist$counts
      if (no_of_bins<2)
      {e <- NaN}
      else{e <- Entropy(res_hist, base = exp(1))/log(length(res_hist))}
      
      ### infer CN
      if (use_median){
        inffered_CN <- min(round(Median(counts), digits = 0), max_CN)
      }
      else {
        inffered_CN <- min(round(mean(counts), digits = 0), max_CN)
      }
      
      cluster_CNs <- append(cluster_CNs, inffered_CN)
      s <- res_hist[inffered_CN+1]/length(counts)
      support <- append(support, s)
      GI <- append(GI, Gini(counts))
      Ent <- append(Ent, e)
      
      CNs[filter] <- inffered_CN
      if (inffered_CN==main_CN) {
        CNs2[filter] <- -1
      }
      else { 
        CNs2[filter] <- inffered_CN
      }
    }
    
    
    avg_GI_inside <- mean(na.omit(GI))
    avg_Ent_inside <- mean(na.omit(Ent))
    avg_cluster_support <- mean(support[which(support>0.01)])
    perc_clusters_with_at_least_0_7_support <- length(which(support>=0.7))/no_of_clusters
    f_bins_changed_by_tree <- which(bins_an[,6:ncol(bins_an)]!=tree_object$name)
    f_bins_far_from_2 <- c(which(all_reads[,6:ncol(all_reads)]<(main_CN-0.5)), which(all_reads[,6:ncol(all_reads)]>(main_CN+0.5)))
    TPCNR <- length(which(f_bins_far_from_2 %in% f_bins_changed_by_tree))/length(f_bins_far_from_2)
    f_bins_close_2 <- which(all_reads[,6:ncol(all_reads)]>(main_CN-0.5) & all_reads[,6:ncol(all_reads)]<(main_CN+0.5))
    FPCNR <- length(which(f_bins_close_2 %in% f_bins_changed_by_tree))/length(f_bins_close_2)
    
    CNs <- as.matrix(CNs)
    CCs_MSE <- all_reads_m[,6:ncol(all_reads_m)]
    CCs_MSE[CCs_MSE>(max_CN+0.5)] <- (max_CN+0.5)
    MSE_CNs <- (sum((CCs_MSE-CNs[,6:ncol(CNs)])^2)/length(CCs_MSE))^0.5

    perc_clusters_CN_2 <- length(which(cluster_CNs==main_CN))/no_of_clusters
    
    stats <- c(no_of_clusters, EDGES, avg_GI_inside, avg_Ent_inside, avg_cluster_support, perc_clusters_with_at_least_0_7_support, GI_2, ENT_2, TPCNR, FPCNR, MSE_CNs, perc_clusters_CN_2)
    names(stats) <- c("no of bin clusters", "no of edges", "avg GI inside", "avg Ent inside", "avg_bin_cluster_support", "perc_bin_clusters_with_at_least_0_7_support", "avg GI outside", "avg Ent outside", "TP-CNR", "FP-CNR", "CNCC-RMSE", "perc_bin_clusters_neutral_CN")
    
  }

  
  CNs <- as.data.frame(CNs)
  colnames(CNs) <- colnames(all_reads)
  rownames(CNs) <- rownames(all_reads)
  CNs2 <- as.data.frame(CNs2)
  colnames(CNs2) <- colnames(all_reads)
  rownames(CNs2) <- rownames(all_reads)
  
  if (only_CNs_for_heatmaps) {
    return(list(CNs, CNs2))
  }
  else {
    return(list(CNs, CNs2, stats))
  }
  
}

################## heatmaps ####################
################## 

# preparing CN input
prepare_CNs_heatmap_input <- function(CNs_data, main_CN) {
  CNs_heatmap <- as.matrix(transpose(CNs_data[,6:ncol(CNs_data)]))
  
  rownames(CNs_heatmap) <- colnames(CNs_data)[6:ncol(CNs_data)]
  colnames(CNs_heatmap) <- CNs_data$chr
  
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
  
  ### color palette
  ### 
  if (main_CN==2) {
    col_breaks <- seq(-0.5, 4.5, by = 1)
    my_palette <- c("#273871", "#8FB7D9", "white", "#FF8870","red4")
  }
  else if (main_CN==3) {
    col_breaks <- seq(-0.5, 6.5, by = 1)
    my_palette <- c("#273871", "steelblue", "lightskyblue2", "white", "#FF8870", "#FF0000", "red4")
  }
  
  else if (main_CN==4) {
    col_breaks <- seq(-0.5, 8.5, by = 1)
    my_palette <- c("#273871","royalblue4", "steelblue", "lightskyblue2", "white", "#FF8870", "#FF0000", "firebrick", "red4")
  }


  return(list(CNs_heatmap, x_labels, x_ticks, col_breaks, my_palette))
}



prepare_CCs_heatmap_input <- function(CCs_data, main_CN) {
  CCs_heatmap <- as.matrix(transpose(CCs_data[,6:ncol(CCs_data)]))
  max_CN <- 2*main_CN

  CCs_heatmap[CCs_heatmap>(max_CN-0.1)] <- (max_CN-0.1)
  
  rownames(CCs_heatmap) <- colnames(CCs_data)[6:ncol(CCs_data)]
  colnames(CCs_heatmap) <- CCs_data$chr
  
  col_breaks_2 <- seq(0, max_CN, by = 0.2)
  n_col <- length(col_breaks_2)-1
  
  my_palette_2 <- c(hcl.colors(n_col/2, palette = "Blues"), 
    hcl.colors(n_col/2, palette = "Reds", rev=-1))
  
  return(list(CCs_heatmap, col_breaks_2, my_palette_2))
}
