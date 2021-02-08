#!/usr/bin/env Rscript

###------------------------------------------------------------------------------
### 
### This script predicts tsRNA classification types
### 
###-------------------------------------------------------------------------------

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

### Check if the correct number of command line arguments were provide. If not, return an error.
if (length(args)==0) {
  stop("Error: Not enough command line arguments provided. Input file and output file names required.")
} 

##### Define functions
### Compare two boolean vectors. Plus one point if both vectors have same TRUE/FALSE switch at each position
vector.comparison <- function(new.bool.vector, template.bool.vector){
  compared.vectors <- ifelse(new.bool.vector == template.bool.vector, 1, 0) # + 1 if the same (e.g. TRUE and TRUE), otherwise 0
  return(sum(compared.vectors)) # Return sum of score 
}
### Compare boolean vector to all tsRNA profiles
tsRNA.classification <- function(bool.vector){
  tiRNA5.score <- vector.comparison(bool.vector, tiRNA5)
  tiRNA3.score <- vector.comparison(bool.vector, tiRNA3)
  itRNA.score <- vector.comparison(bool.vector, itRNA)
  tRF5.score <- vector.comparison(bool.vector, tRF5)  
  tRF3.score <- vector.comparison(bool.vector, tRF3)
  tsRNA.scores <- data.frame(t(data.frame(tiRNA5.score,
                                          tiRNA3.score, 
                                          itRNA.score, 
                                          tRF5.score, 
                                          tRF3.score)))
  names.vector <- c("tiRNA5", "tiRNA3", "itRNA", "tRF5", "tRF3")
  tsRNA.scores <- cbind(names.vector, tsRNA.scores)
  colnames(tsRNA.scores) <- c("tsRNA.type", "score")
  tsRNA.scores$percentage.score <- (tsRNA.scores$score/70) * 100
  tsRNA.scores <- arrange(tsRNA.scores, desc(score))
  return(tsRNA.scores)
}
##### Define function end

### Define profiles for different types of tsRNAs.
### TRUE at a nucleotide position indicates that the expected 
### coverage should be above average. FALSE is vice versa
tiRNA5 <- c(rep("TRUE", 33), rep("FALSE", 37))
tiRNA3 <- c(rep("FALSE", 33), rep("TRUE", 37))
itRNA <- c(rep("FALSE", 15), rep("TRUE", 35), rep("FALSE", 20))
tRF5 <- c(rep("TRUE", 15), rep("FALSE", 55))
tRF3 <- c(rep("FALSE", 50), rep("TRUE", 20))

#### Input file
input <- read.table(args[1])

df <- split( input , f = input$V1 )  # Split dataframe based on column 1 elements

results.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("feature",
                                                                 "tsRNA.type",
                                                                 "tsRNA.classification.scores")) 

for(subset in df) {
  feature <- as.character(subset[1,1])
  subset.length <- nrow(subset)
  subset.cov.mean <- mean(subset$V3)
  if (subset.cov.mean < 5) { # If coverage low, don't calculate most likely tsRNA type
    tsRNA.type <- "Undetermined"
    classification.results = "NA"
  } else {
    subset <- subset %>%
      mutate("HighCov" = ifelse(subset$V3 > subset.cov.mean, "TRUE", "FALSE"))
    subset.highcov.bool <- subset$HighCov
    number.of.rows.to.remove <- subset.length - 70
    rows.to.remove <- sample(1:subset.length, number.of.rows.to.remove) # Randomly choose n numbers between 1 and the length of the tRNA
    new.subset <- subset[-c(rows.to.remove),] # Remove the randomly selected rows from the tRNA dataframe to create DF 70 rows long
    rownames(new.subset) <- seq(1:nrow(new.subset))
    ### Determine best tsRNA type using functions from above
    classification.df <- tsRNA.classification(new.subset$HighCov)
    ### Choose most suitable tsRNA type(s)
    
    if (classification.df[[3]][1] > 90) { # If the top hit scored over 90%
      tsRNA.type <- as.character(classification.df[[1]][1])
    } else if (classification.df[[3]][1] > 80) { # Best match is between 65% - 80%
      if (classification.df[[3]][2] > 80) { # If the second hit is also over 80%
        tsRNA.type <- paste0(as.character(classification.df[[1]][1]), "-or-", as.character(classification.df[[1]][2]))
      } else { # Only one tsRNA type above 80%
        tsRNA.type <- as.character(classification.df[[1]][1])
      } 
    } else if (classification.df[[3]][1] > 65) { # Best match is between 65% - 80%
      tsRNA.type <- paste0("Possibly-", as.character(classification.df[[1]][1]))
    } else {
      tsRNA.type <- "Undetermined"
    }
    is.num <- sapply(classification.df, is.numeric) # Convert to numeric
    classification.df[is.num] <- lapply(classification.df[is.num], round, 2) # Round down to 2 decimal places
    classification.df <- lapply(classification.df, as.vector) # Convert to vector
    classification.results <- toString(classification.df) # Convert to string
  }
  results.df[nrow(results.df) + 1,] = list(feature,
                                           tsRNA.type,
                                           classification.results)
}

write.table(results.df, 
            file = args[2],
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

