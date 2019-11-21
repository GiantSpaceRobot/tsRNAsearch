#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(VennDiagram)
library(gplots)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # Suppress VennDiagram messages/log file generation

# Borrowed this code from https://www.r-bloggers.com/working-with-venn-diagrams/

input1 <- try(read.table(args[1]), silent = TRUE)
if (inherits(input1, 'try-error')){ 
  input1 <- data.frame(feature=character())
} 

input2 <- try(read.table(args[2]), silent = TRUE)
if (inherits(input2, 'try-error')){ 
  input2 <- data.frame(feature=character())
} 

input3 <- try(read.table(args[3]), silent = TRUE)
if (inherits(input3, 'try-error')){ 
  input3 <- data.frame(feature=character())
} 

input4 <- try(read.table(args[4], sep = "\t"), silent = TRUE)
if (inherits(input4, 'try-error')){ 
  input4 <- data.frame(feature=character())
}
input4 <- data.frame(gsub(' .*', '' , input4$V1))
colnames(input4) <- c("V1")

geneLists <- c("DE" = input1, #DESeq2 
               "DS" = input2, #Distribution score
               "CS" = input3, #Cleavage score
               "FM" = input4) #Fisher's method

my.venn <- venn.diagram(geneLists, 
                   filename = NULL,
                   #fill=c("darkmagenta", "darkblue", "red"), 
                   fill=c("#3e4574", "#00a9ff", "#ff0c3e", "#ffc60c"),
                   alpha=c(0.5,0.5,0.5,0.5), 
                   cex = 2, 
                   cat.fontface=4,
                   imagetype = "png",
                   cat.cex = 1.5,
                   category.names=c("DESeq2", "Distribution", "Cleavage", "Fisher"),
                   main.cex = 1.5,
                   main="Features identified by the four tsRNAsearch methods")

png(file = paste0(args[4], "_VennDiagram.png"))
grid.newpage() # Create new grid
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc"))) # Resize this grid (80%)
grid.draw(my.venn) # Draw venn diagram in new smaller plot
dev.off()

a <- venn(geneLists, show.plot=FALSE)

# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB

write.table(inters$`DE.V1:DS.V1:CS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.all.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:DS.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Distribution_Cleavage.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DS.V1:CS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.Distribution_Cleavage_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:CS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Cleavage_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:DS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Distribution_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:DS.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Distribution.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Cleavage.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DE.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.DESeq2_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DS.V1:CS.V1`, 
            file = paste0(args[4], ".intersect.Distribution_Cleavage.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`DS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.Distribution_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$`CS.V1:FM.V1`, 
            file = paste0(args[4], ".intersect.Cleavage_Fisher.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$DE.V1, 
            file = paste0(args[4], ".intersect.DESeq2-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$DS.V1, 
            file = paste0(args[4], ".intersect.Distribution-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$CS.V1, 
            file = paste0(args[4], ".intersect.Cleavage-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
write.table(inters$FM.V1, 
            file = paste0(args[4], ".intersect.Fisher-only.txt"),
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
