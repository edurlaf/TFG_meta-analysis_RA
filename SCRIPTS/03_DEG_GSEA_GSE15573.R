library(GEOquery)
library(edgeR)
load(file = "../GITHUB/GSE_AR_annotated.RData")
data <- exprs(GSE[["GSE15573"]])
variable <- pData(GSE[["GSE15573"]])$Diagnosis_Sex
contrasts <- c("(RA_Female - Control_Female)",  "(RA_Male - Control_Male)", "(RA_Female - Control_Female) - (RA_Male - Control_Male)")
batch = FALSE
paired = FALSE
trend = FALSE

# Crear matriz de diseño
design <- stats::model.matrix( ~ 0 + variable)
colnames(design) <- c(levels(variable))

# Make contrast matrix, depending on the number of contrasts
contrastMatrix <- limma::makeContrasts(
  contrasts = contrasts,
  levels = design
)

# Fit linear model for each gene
fit <- limma::lmFit(data, design)

fit2 <- limma::contrasts.fit(fit, contrastMatrix)
fit3 <- limma::eBayes(fit2, trend = trend)

# p-valor ajustado por metodo Benjamini-Hochberg
fit3$p.adj <- apply (fit3$p.value, 2, p.adjust, method = "BH")

# Escribir archivo con los resultados del ajuste
write.fit(fit3, results = NULL, file = "DEGfile_15573.txt", digits = NULL,
          method = "separate", quote = FALSE, sep = "\t", row.names = TRUE)

# Contraste interés (Coeficiente 3)
top.table <- topTable(fit3, sort.by = "P", confint=TRUE, n = Inf, coef = 3)
length(which(top.table$adj.P.Val < 0.05))

# Calcular SE de las 2 maneras:
top.table[, "SE"] <- (top.table[, "CI.R"] - top.table[, "CI.L"])/ 3.92 #STEVE
SE.coef <- sqrt(fit3$s2.post) * fit3$stdev.unscaled #GORDON
head(SE.coef)
SE.coef <- SE.coef[rownames(top.table),]
top.table[, "SE.coef"] <- SE.coef[,3]
head(top.table)

# Contraste 2
top.table2 <- topTable(fit3, sort.by = "P", confint=TRUE, n = Inf, coef = 2)
length(which(top.table2$adj.P.Val < 0.05))

# Contraste 1
top.table1 <- topTable(fit3, sort.by = "P", confint=TRUE, n = Inf, coef = 1)
head(top.table1)
length(which(top.table1$adj.P.Val < 0.05))

# Significativos
sig1 <- top.table1[which(top.table1$adj.P.Val < 0.05),]
length(which(sig1$logFC < 0))
length(which(sig1$logFC >= 0))

# Generar RData como input para el MA
top.table <- tibble::rownames_to_column(top.table, "EntrezID")
save(top.table, file = "GSE15573inputMA.RData")


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
               values = rownames(fit),
               mart = ensembl)

# annot <- select(org.Hs.eg.db, keys = rownames(fit3), 
#                 columns = c("ENTREZID", "GO"), keytype = "ENTREZID")

#colnames(annot) <- c("ENTREZID", "GO")
annot <- annot[annot[,"go_id"] != "",]

# Función doGSEA: biological process
load("doGSEA.R")
GOresults <- doGSEA(fit = fit3, annot = annot, ontology = "bp",
                    propagate = TRUE, contrasts = 3)
save(GOresults, file = "GSE15573_Goresultsbp.rda")
bp <- GOresults[[3]]
bpsig <- bp[bp$padj<0.05,]
save(bpsig, file = "GSE15573_GOresultsbp_sig.rda")
write.csv2(bpsig, file = 'GSE15573bpsig.csv')

# Funcion doGSEA: KEGG pathway
# grep("Kegg", listAttributes(ensembl)[,1], ignore.case = T)
# listAttributes(ensembl)[70,]
# annot2 <- getBM(attributes = c("entrezgene_id","kegg_enzyme"),
#                filters = "entrezgene_id",
#                values = rownames(fit3),
#                mart = ensembl)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
annot2bis <- select(org.Hs.eg.db, keys = rownames(fit), 
                 columns = c("ENTREZID", "PATH"), keytype = "ENTREZID")
# colnames(annot) <- c("ENTREZID", "GO")
annot2bis <- annot2bis[complete.cases(annot2bis),]
GOresults <- doGSEA(fit = fit3, annot = annot2bis, ontology = "path",
                    propagate = TRUE, contrasts = 3)

save(GOresults, file = "GSE15573_KEGGresults.rda")
kegg <- GOresults[[3]]
keggsig <- kegg[kegg$padj<0.05,]
save(keggsig, file = "GSE15573_KEGGresults_sig.rda")
write.csv2(keggsig, file = 'GSE15573Keggsig.csv')
