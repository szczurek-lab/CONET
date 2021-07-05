library(HMMcopy)



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
  for (i in 0:49) {
    id <- paste(paste("NEW_NOISE", cell, sep = ""), "40" , i, sep = "_")
    id3 <- paste(paste("NEW_NOISE", cell, sep = ""), "20" , i, sep = "_")
    
    id2 <- paste(paste("NO_NOISE", cell, sep = ""), "40" , i, sep = "_")
    id4 <- paste(paste("NO_NOISE", cell, sep = ""), "20" , i, sep = "_")
    IDs <- c(IDs, id, id3, id, id3, id2, id4, id2, id4)
  }
}
file_conn <- file("hmm_results2", open = "wt")
writeLines("ID;MSE;SD;FP;FN;AvergaeEvents;Singletons", con=file_conn)

PIPELINE <- function(ID, save) {
  all_counts <- t(read.csv(paste("models/corrected_counts_", ID,sep=""), sep = ";", header = F))
  real_counts <- t(read.csv(paste("models/counts_", ID,sep=""), sep = ";", header = F))
  
  
  
  
  create_cell <- function(cell) {
    dummy <- rep(1, dim(all_counts)[1])
    res <- data.frame(dummy, (1:dim(all_counts)[1]), (1:(dim(all_counts)[1])), dummy, dummy,
                      dummy, rep(FALSE, dim(all_counts)[1]), rep(TRUE, dim(all_counts)[1]),
                      cell, cell, cell)
    names(res) <- c("chr" ,    "start"  , "end"  ,   "reads"  , "gc"   ,   "map"   ,  
                    "valid" ,  "ideal" ,  "cor.gc" , "cor.map" ,"copy"  )
    
    return (res)
  }
  
  get_cell_counts <- function(i) {
    res <- HMMsegment(create_cell(all_counts[, i]), getparam=F,maxiter = 50000)[1]$state -1 
    return (res)
    
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
    write.table(inferred_counts, "HMM_counts", row.names = FALSE, col.names = FALSE)
  }
  
  
  brkp <- get_brkp_matrix(inferred_counts)
  
  brkp2 <- rowSums(brkp)
  
  x <- unique(unlist(which(brkp2 >0) -1))
  diffs <- read.table(paste("models/diffs_", ID, sep = ""), sep = ";")
  for (r in 1:dim(brkp)[1]) {
    if (quantile(abs(diffs[,r]), 0.7) > 1) {
      x <- c(x, r-1)
    }
  }
  x <- sort(unique(x))
  # 
  write.table(x, paste("models/indices_", ID, sep = ""), row.names= FALSE, col.names = FALSE)
  
  r <- read.table(paste("models/brkp_", ID, sep = ""), sep = ";")
  x <- which(colSums(r) > 0) - 1
  
  write.table(x, paste("models/real_indices_", ID, sep = ""), row.names= FALSE, col.names = FALSE)
}





for (id in IDs) {
  PIPELINE(id, TRUE)
}
close(file_conn)
