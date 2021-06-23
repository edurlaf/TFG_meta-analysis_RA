# Script to generate a "pyramid" plot
# comparing the percentages of enriched GO terms assinged
# to each category of Biological Process GOslims

library(dplyr)
library(ggplot2)


#####################################################
# Set the following variables for the appropriate comparisons/files
#####################################################
# Comparison
comparison <- "(Mujer_AR-Mujer_C)-(Hombre_AR-Hombre_C)"

# Grupos
Female <- "Mujer_AR"
Male <- "Hombre_AR"

# Read in first comparison files

# df1 <- read.csv("analyses/infected-uninfected/P0.05_C1.infected-UP.subset.GOseq.enriched.flattened.FDR_1.0.BP.GOslims.csv")
# 
# # Read in second comparison file
# df2 <- read.csv("analyses/infected-uninfected/P0.05_C1.uninfected-UP.subset.GOseq.enriched.flattened.FDR_1.0.BP.GOslims.csv")


######################################################
# CHANGES BELOW HERE ARE PROBABLY NOT NECESSARY
######################################################

# GOslim categories
# ontologies <- c("BP", "CC", "MF")
ontology <- "BP"

# Remove generic "biological_process" category
# df1 <- df1[df1$GOslim != "GO:0008150",]
maledf <- male_slimdf[rownames(male_slimdf) != "GO:0008150",]
# df2 <- df2[df2$GOslim != "GO:0008150",]
femaledf <- female_slimdf[rownames(female_slimdf) != "GO:0008150",]

# Select columns
# df1 <- df1 %>% select(Term, Percent)
malesel <- maledf[,c("Term", "Count")]
# df2 <- df2 %>% select(Term, Percent)
femalesel <- femaledf[,c("Term", "Count")]

# Create treatment column and assign term to all rows
# df1$treatment <- treatment_01
malesel$group <- as.factor(Male)
# df2$treatment <- treatment_02
femalesel$group <- as.factor(Female)

# Concatenate dataframes
# df3 <- rbind(df1, df2)
df3 <- rbind(malesel, femalesel)
df3sel <- df3[df3$Count>0,]
df3sel <- df3sel[order(df3sel$Count, decreasing = T),]
# Filename for plot
pyramid <- paste("GOslims", "BP", "sel", "png", sep = ".")
pyramid_path <- paste(pyramid, sep = "/")
pyramid_dest <- file.path("./", pyramid_path)

# "Open" PNG file for saving subsequent plot
png(pyramid_dest, width = 600, height = 1200, units = "px", pointsize = 12)

# Create "pyramid" plot
prueba <- ggplot(df3sel, aes(x = reorder(Term, Count), fill = group, 
                y = ifelse(test = group == Male, 
                           yes = -Count, 
                           no = Count))) + 
  geom_bar(stat = "identity") + 
  scale_y_continuous(labels = abs, limits = max(df3sel$Count) * c(-1,1)) + 
  labs(title = "Counts of GO terms assigned to BP GOslims", x = "GOslim", y = "Counts of GO terms in GOslim") + 
  scale_x_discrete(expand = c(-1,0))
prueba + coord_flip()
# Close PNG file
dev.off()
