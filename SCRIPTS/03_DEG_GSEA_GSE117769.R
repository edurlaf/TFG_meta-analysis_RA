library(GEOquery)
library(edgeR)
load(file = "../GITHUB/GSE_AR_annotated.RData")
data <- exprs(GSE[["GSE117769"]])
variable <- pData(GSE[["GSE117769"]])$Diagnosis_Sex
contrasts <- c("(RA_Female - Control_Female)",  "(RA_Male - Control_Male)", "(RA_Female - Control_Female) - (RA_Male - Control_Male)")
batch = FALSE
paired = FALSE
trend = TRUE

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

# p-value ajustado por método Benjamini-Hochberg
fit3$p.adj <- apply (fit3$p.value, 2, p.adjust, method = "BH")

# Genes diferenciales
summary(decideTests(fit3, method="separate"))

# Escribir archivo con los resultados del ajuste
write.fit(fit3, results = NULL, file = "DEGfile_117769.txt", digits = NULL,
          method = "separate", quote = FALSE, sep = "\t", row.names = TRUE)

# Contraste de interes (coef 3)
top.table <- topTable(fit3, sort.by = "P", confint=TRUE, n = Inf, coef = 3)
head(top.table, 20)
dim(top.table)
length(which(top.table$adj.P.Val < 0.05))

# Calcular el error estándar (SE) de las 2 maneras:
top.table[, "SE"] <- (top.table[, "CI.R"] - top.table[, "CI.L"])/ 3.92 #STEVE
SE.coef <- sqrt(fit3$s2.post) * fit3$stdev.unscaled #GORDON
head(SE.coef)
SE.coef <- SE.coef[rownames(top.table),]
top.table[, "SE.coef"] <- SE.coef[,3]

# # significativos del contraste de interes (coef 3)
# sig3 <- top.table[which(top.table$adj.P.Val < 0.05),]
# entrezsig3 <- rownames(sig3)
# fdata <- fData(GSE[["GSE93272"]])
# sig3$Genename <- fdata[rownames(sig3),"GENENAME"]
# sig <- sig3[,-c(2,3,6)]
# length(which(sig$logFC >= 0))
# length(which(sig$logFC < 0))

# contraste 1
top.table1 <- topTable(fit3, sort.by = "P", n = Inf, coef = 1)
length(which(top.table1$adj.P.Val < 0.05))
# genes significativos
sig1 <- top.table1[which(top.table1$adj.P.Val < 0.05),]
length(which(sig1$logFC < 0))

# contraste 2
top.table2 <- topTable(fit3, sort.by = "P", n = Inf, coef = 2)
length(which(top.table2$adj.P.Val < 0.05))
# genes significativos
sig2 <- top.table2[which(top.table2$adj.P.Val < 0.05),]
length(which(sig2$logFC < 0))

# Generar RData como input para el MA
top.table <- tibble::rownames_to_column(top.table, "EntrezID")
save(top.table, file = "GSE117769_inputMA.RData")
save(top.table1, file = "inputMA1_GSE117769.RData")
save(top.table2, file = "inputMA2_GSE117769.RData")

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
               values = rownames(fit3),
               mart = ensembl)

# annot <- select(org.Hs.eg.db, keys = rownames(fit3), 
#                 columns = c("ENTREZID", "GO"), keytype = "ENTREZID")

#colnames(annot) <- c("ENTREZID", "GO")
annot <- annot[annot[,"go_id"] != "",]

# Función doGSEA: biological process
load("doGSEA.R")
GOresults <- doGSEA(fit = fit3, annot = annot, ontology = "bp",
                    propagate = TRUE, contrasts = 3)
save(GOresults, file = "GSE117769_Goresultsbp.rda")
bp <- GOresults[[3]]
bpsig <- bp[bp$padj<0.05,]

# Añadir columna name con el nombre de los GOs
bpsig[,"name"] <- mdgsa::getGOnames(rownames(bpsig))

# Reordenar
bpsig <- bpsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
bpsig <- bpsig[order(bpsig$lor),]

# Guardar
save(bpsig, file = "GSE117769_GOresultsbp_sig_name.rda")
write.csv2(bpsig, file = 'GSE117769bpsig_name.csv')

# Funcion doGSEA: KEGG pathway
# grep("Kegg", listAttributes(ensembl)[,1], ignore.case = T)
# listAttributes(ensembl)[70,]
# annot2 <- getBM(attributes = c("entrezgene_id","kegg_enzyme"),
#                filters = "entrezgene_id",
#                values = rownames(fit3),
#                mart = ensembl)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
annot2 <- select(org.Hs.eg.db, keys = rownames(fit), 
                 columns = c("ENTREZID", "PATH"), keytype = "ENTREZID")
# colnames(annot) <- c("ENTREZID", "GO")
annot2 <- annot2[complete.cases(annot2),]
KEGGresults <- doGSEA(fit = fit3, annot = annot2, ontology = "path",
                    propagate = TRUE, contrasts = 3)

save(KEGGresults, file = "GSE117769_KEGGresults.rda")
kegg <- KEGGresults[[3]]
keggsig <- kegg[kegg$padj<0.05,]

# Añadir columna name con el nombre de los KEGGs
keggsig[,"name"] <- mdgsa::getKEGGnames(rownames(keggsig))

# Reordenar
keggsig <- keggsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
keggsig <- keggsig[order(keggsig$lor),]

# Guardar
save(keggsig, file = "GSE117769_KEGGresults_sig_name.rda")
write.csv2(keggsig, file = 'GSE117769keggsig_name.csv')
