#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("VennDiagram")
library("gplots")

# Borrowed this code from https://www.r-bloggers.com/working-with-venn-diagrams/

input1 <- try(read.table(args[1]), silent = TRUE)
input2 <- try(read.table(args[2]), silent = TRUE)
input3 <- try(read.table(args[3]), silent = TRUE)
geneLists <- c("DE" = input1, #DESeq2 
               "DS" = input2, #Distribution score
               "CS" = input3) #Cleavage score

print(venn.diagram(geneLists, 
             filename = paste0(args[4], "_VennDiagram.pdf"), 
             #fill=c("darkmagenta", "darkblue", "red"), 
             fill=c("#3e4574", "#00a9ff", "#ff0c3e"),
             alpha=c(0.5,0.5,0.5), 
             cex = 2, 
             cat.fontface=4,
             #cat.pos = c(340,20,0),
             #cat.col = c("#3e4574", "#00a9ff", "#ff0c3e"),
             cat.cex = 1.5,
             category.names=c("DESeq2", "Distribution", "Cleavage"),
             main.cex = 1.5,
             main="Features identified by the three tsRNAsearch_DE methods"))


a <- venn(geneLists, show.plot=FALSE)

# You can inspect the contents of this object with the str() function
#str(a)

# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
intersections <- lapply(inters, head) 

write.table(intersections$`DE.V1:DS.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.all.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$`DE.V1:DS.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Distribution.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$`DE.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Cleavage.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$`DS.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.Distribution_Cleavage.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$DE.V1, 
            file = paste0(args[4], ".intersect.DESeq2-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$DS.V1, 
            file = paste0(args[4], ".intersect.Distribution-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(intersections$CS.V1, 
            file = paste0(args[4], ".intersect.Cleavage-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
