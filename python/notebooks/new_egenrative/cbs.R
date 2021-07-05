library(DNAcopy)
library(aCGH)


get_smoothed_cell <- function(cell_, output) {
  result <- rep(0, length(cell_$chr))
  for (i in 1:length(result)) {
    start <- cell_$start[i]
    chr <- cell_$chr[i]
    region <- output[output$chrom == chr & output$loc.start <= start & output$loc.end >= start, ]
    if (dim(region)[1] != 1) {
      region <- cell_$cell[i]
    } else {
      region <- max(region$seg.mean)
    }
    result[i] <- region
  }
  return (result)
}

get_cell_ratios <- function(cell) {
  return (log2(cell))
}




get_brkp_matrix <- function(counts) {
  res <- counts
  res[] <- 0
  for (i in 1:dim(counts)[2]) {
    for (j in 1:dim(counts)[1]) {
      if (( j== 1 && counts[j, i] != 2) ){
        res[j,i] <- 1
      } else if (j > 1 && (counts[j,i] != counts[j-1,i])) {
        res[j,i] <- 1
      }
    }
  }
  return (res)
}



IDs <- c()
for (cell in c(200, 1000)) {
  for (i in 0:4) {
    id2 <- paste(paste("NEW_NOISE", cell, sep = ""), "40" , i, sep = "_")
    id <- paste(paste("NEW_NOISE", cell, sep = ""), "20" , i, sep = "_")
    IDs <- c(IDs, id, id2)
  }
}

file_conn <- file("cbs_results", open = "wt")
writeLines("ID;MSE;SD;FP;FN;AvergaeEvents;Singletons", con=file_conn)



PIPELINE <- function(ID, save) {
  all_counts <- t(read.csv(paste("models/corrected_counts_", ID,sep=""), sep = ";", header = F))
  real_counts <- t(read.csv(paste("models/counts_", ID,sep=""), sep = ";", header = F))
  
  get_cell_counts <- function(i) {
    print(i)
    cell <- get_cell_ratios(all_counts[,i])
    
    cell_ <- data.frame(rep(10, dim(all_counts)[1]), cell, 15000 * 1:(dim(all_counts)[1]))
    colnames(cell_) <- c("chr", "cell", "start")
    
    
    
    CNA.object <- CNA(cbind(cell_$cell),
                      cell_$chr, cell_$start,
                      data.type="logratio",sampleid="c05296")
    
    
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    
    
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=0, alpha = 0.01, undo.prune = 0.05)
    
    #plot(segment.smoothed.CNA.object)
    smoothed_cell <- get_smoothed_cell(cell_, segment.smoothed.CNA.object$output)
    vecObs <- cell_$cell
    vecPred<- smoothed_cell
    
    # Run merge function
    merge.obj <- mergeLevels(vecObs,vecPred, verbose = 0)
    
    # Examine optimum threshold
    merge.obj$sq
    
    
    counts <- round(2^merge.obj$vecMerged, digits=0)
    
    return (counts)
    
  }
  
  inferred_counts <- all_counts
  for ( i in 1:dim(all_counts)[2]) {
    inferred_counts[, i] <- get_cell_counts(i)
  }
  
  MSE <- mean((inferred_counts - real_counts)^2)
  
  
  real_brp <- get_brkp_matrix(real_counts)
  inf_brp <- get_brkp_matrix(inferred_counts)
  
  SD <- sum(abs(real_brp - inf_brp)) / dim(all_counts)[2]
  
  FP <- sum(inf_brp == 1 & real_brp == 0) / sum(inf_brp == 1)
  
  FN <- sum(inf_brp == 0 & real_brp == 1) / sum(real_brp == 1)

  
  result <- matrix(0L, nrow = 1000000, ncol = 4)
  

  extract_events <- function(cell, counts, result, result_index) {
    i <- 1
    while (i < dim(counts)[1]) {
      count <- counts[i, cell]
      if (count != 2) {
        start <- i
        i <- i +1
        while (i < dim(counts)[1] && counts[i, cell] == count) {
          i <- i +1
        }
        end <- i - 1
        result[result_index, ] <- c(cell, count, start, end)
        result_index <- result_index + 1
      } else {
        i <- i +1
      }
      
    }
    return ( list(result, result_index))
  }
  
  
  
  ind <- 1
  for (i in 1:dim(inferred_counts)[2]) {
    res <- extract_events(i, inferred_counts, result, ind)
    ind <- res[[2]]
    result <- res[[1]]
  }
  
  
  ind <- ind-1
  EN <- print(ind/260)
  
  result <- data.frame(result[1:ind, ])
  
  result$hash <- paste(result$X2, result$X3, result$X4, sep="_")
  res <- 0
  for (i in 1:ind) {
    hash <- result$hash[i]
    events <- length(which(result$hash == hash))
    if (events == 1) {
      res <- res + 1
    }
  }
  
  SINGLETONS <- res
  writeLines(paste(ID,MSE, SD, FP, FN, EN, SINGLETONS, sep = ";"), con=file_conn)
  if (save) {
    write.table(inferred_counts, "CBS_counts", row.names = FALSE, col.names = FALSE)
  }
}




for (id in IDs) {
  PIPELINE(id, FALSE)
}
close(file_conn)
