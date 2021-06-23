# Cargar fichero input (resultados MA genes ordenados por logFC de mayor a menor)
load("MA_genes_results.rda")
MA_genes <- MAgenes

# Perform GSEA analysis ------------------------------------------------------
library(biomaRt)
# Fist we get the annotation
# listMarts() # to see which database options are present
ensembl <-  useMart("ENSEMBL_MART_ENSEMBL")
# listDatasets(ensembl) # function to see which datasets are present in ensembl
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
# listFilters(ensembl) # check which filters are available
# listAttributes(ensembl) # check attributes are available to select
annot <- getBM(attributes = c("entrezgene_id","go_id"),
               filters = "entrezgene_id",
               values = rownames(MAgenes),
               mart = ensembl)

# annot <- select(org.Hs.eg.db, keys = rownames(fit3), 
#                 columns = c("ENTREZID", "GO"), keytype = "ENTREZID")

#colnames(annot) <- c("ENTREZID", "GO")
annot <- annot[annot[,"go_id"] != "",]

# Función doGSEA: biological process
load("doGSEA.R")
GOresults <- doGSEA(fit = as.data.frame(MAgenes), annot = annot, ontology = "bp",
                    propagate = TRUE, contrasts = 1)


save(GOresults, file = "Global_enrichment.rda")
bpsig <- GOresults[GOresults$padj < 0.05,]

# Añadir columna name con el nombre de los GOs
bpsig[,"name"] <- mdgsa::getGOnames(rownames(bpsig))

# Reordenar
bpsig <- bpsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
bpsig <- bpsig[order(bpsig$lor),]

save(bpsig, file = "Global_enrichment_sig.rda")
write.csv2(bpsig, file = 'Global_enrichment_sig.csv')

# prueba slim
library(GSEABase)
library(tidyverse)
maleup <- bpsig[bpsig$lor<0,]
femaleup <- bpsig[bpsig$lor>0,]
myIds <- rownames(maleup)
myCollection <- GOCollection(myIds)
#fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
slim <- getOBOCollection("./goslim_generic.obo")
# slim <- getOBOCollection(fl)

# Sets ontology category to "Biological Process"
slimdf <- goSlim(myCollection, slim, ontology = "BP")

# This should match the ontology used above.
# E.g. GOBPOFFSPRING or GOCCOFFSPRING
GO.db::GOBPOFFSPRING


mappedIds <-
  function(df, collection, OFFSPRING) {
    map <- as.list(OFFSPRING[rownames(df)])
    mapped <- lapply(map, intersect, ids(collection))
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L))
    df
  }

male_slimdf <- mappedIds(slimdf, myCollection, GOBPOFFSPRING)

save(male_slimdf, file = "maleGOslim.rda")

# para female
myIds <- rownames(femaleup)
myCollection <- GOCollection(myIds)
#fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
slim <- getOBOCollection("./goslim_generic.obo")
# slim <- getOBOCollection(fl)

# Sets ontology category to "Biological Process"
fem_slimdf <- goSlim(myCollection, slim, ontology = "BP")

# This should match the ontology used above.
# E.g. GOBPOFFSPRING or GOCCOFFSPRING
GO.db::GOBPOFFSPRING


mappedIds <-
  function(df, collection, OFFSPRING) {
    map <- as.list(OFFSPRING[rownames(df)])
    mapped <- lapply(map, intersect, ids(collection))
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L))
    df
  }

female_slimdf <- mappedIds(fem_slimdf, myCollection, GOBPOFFSPRING)

save(female_slimdf, file = "femaleGOslim.rda")


write.csv2(slimdf, file = "gsea_goslim.csv", quote = FALSE)

slimmer <- slimdf[slimdf$Count>0,]
write.csv2(slimmer, file = "gsea_goslimmer.csv", quote = FALSE)

# Funcion doGSEA: KEGG pathway
# grep("Kegg", listAttributes(ensembl)[,1], ignore.case = T)
# listAttributes(ensembl)[70,]
# annot2 <- getBM(attributes = c("entrezgene_id","kegg_enzyme"),
#                filters = "entrezgene_id",
#                values = rownames(fit3),
#                mart = ensembl)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
annot2 <- select(org.Hs.eg.db, keys = rownames(MAgenes), 
                 columns = c("ENTREZID", "PATH"), keytype = "ENTREZID")
# colnames(annot) <- c("ENTREZID", "GO")
annot2 <- annot2[complete.cases(annot2),]
KEGGresults <- doGSEA(fit = MAgenes, annot = annot2, ontology = "path",
                    propagate = TRUE, contrasts = 1)

save(KEGGresults, file = "Global_enrichment_KEGGresults.rda")

keggsig <- KEGGresults[KEGGresults$padj<0.05,]

# Añadir columna name con el nombre de los KEGGs
keggsig[,"name"] <- mdgsa::getKEGGnames(rownames(keggsig))

# Reordenar
keggsig <- keggsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
keggsig <- keggsig[order(keggsig$lor),]
save(keggsig, file = "Globalenrichment_KEGGresults_sig.rda")
write.csv2(keggsig, file = 'Globalenrichment_Keggsig.csv')
