library(gplots)

args = commandArgs(trailingOnly=TRUE)

input1 <- try(read.table(args[1]), silent = TRUE)
input2 <- try(read.table(args[2]), silent = TRUE)
input3 <- try(read.table(args[3]), silent = TRUE)

if(class(input1)=="try-error"){
  input1 <- data.frame()
  DE <- NA
} else {
  DE <- as.vector(input1[,1])
}

if(class(input2)=="try-error"){
  input2 <- data.frame()
  HD <- NA
} else {
  HD <- as.vector(input2[,1])
}

if(class(input2)=="try-error"){
  input2 <- data.frame()
  PC <- NA
} else {
  PC <- as.vector(input3[,1])
}

pdf(file = paste0(args[4], ".VennDiagram.pdf"), 
    width = 8, 
    height = 8)
venn(list(Differentially.Expressed=input1$V1, 
          High.Distribution=input2$V1,
          Potentially.Cleaved=input3$V1))
dev.off()

intersect.DE.HD <- Reduce(intersect, list(DE, HD))
intersect.PC.HD <- Reduce(intersect, list(PC, HD))
intersect.DE.PC <- Reduce(intersect, list(DE, PC))
intersect.all <- Reduce(intersect, list(DE, HD, PC))

write.table(intersect.all, 
            file = paste0(args[4], ".intersect.all.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersect.DE.HD, 
            file = paste0(args[4], ".intersect.DiffExpr-vs-HighDistr.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersect.DE.PC, 
            file = paste0(args[4], ".intersect.DiffExpr-vs-PotCleav.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersect.PC.HD, 
            file = paste0(args[4], ".intersect.PotCleav-vs-HighDistr.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
