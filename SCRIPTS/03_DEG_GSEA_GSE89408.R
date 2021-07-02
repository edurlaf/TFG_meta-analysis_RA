library(GEOquery)
library(edgeR)
load(file = "../GITHUB/GSE_AR_annotated.RData")
data <- exprs(GSE[["GSE89408"]])
variable <- pData(GSE[["GSE89408"]])$Diagnosis_Sex
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

# Aplicar voom
vdata <- limma::voom(data, design, plot = T)

# Fit linear model for each gene
fit <- limma::lmFit(vdata, design)
fit2 <- limma::contrasts.fit(fit, contrastMatrix)
fit3 <- limma::eBayes(fit2, trend = trend)

# p-valor ajustado por metodo Benjamini-Hochberg
fit3$p.adj <- apply (fit3$p.value, 2, p.adjust, method = "BH")

# Escribir archivo con los resultados del ajuste
write.fit(fit3, results = NULL, file = "DEGfile_89408.txt", digits = NULL,
          method = "separate", quote = FALSE, sep = "\t", row.names = TRUE)

# Genes diferenciales
summary(decideTests(fit3, adjust.method = "fdr"))

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

# significativos del contraste de interes (coef 3)
sig3 <- top.table[which(top.table$adj.P.Val < 0.05),]
entrezsig3 <- rownames(sig3)
fdata <- fData(GSE[["GSE89408"]])
sig3$Genesymbol <- fdata[rownames(sig3),"hgnc_symbol"]
sig3$Genename <- fdata[rownames(sig3),"description"]
sig <- sig3[,-c(2,3,6)]

# Escribir excel con genes significativos del contraste de interes
write.csv2(sig, file = 'GSE89408sig.csv')

# Contraste 2
top.table2 <- topTable(fit3, sort.by = "P", n = Inf, coef = 2)
length(which(top.table2$adj.P.Val < 0.05))

# Contraste 1
top.table1 <- topTable(fit3, sort.by = "P", n = Inf, coef = 1)
head(top.table1, 20)
length(which(top.table1$adj.P.Val < 0.05))

# Generar RData como input para el MA
top.table <- tibble::rownames_to_column(top.table, "EntrezID")
save(top.table, file = "GSE89408_inputMA.RData")
save(top.table1, file = "inputMA1_GSE89408.RData")
save(top.table2, file = "inputMA2_GSE89408.RData")

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

annot <- annot[annot[,"go_id"] != "",]

# Función doGSEA: biological process
load("doGSEA.R")
GOresults <- doGSEA(fit = fit3, annot = annot, ontology = "bp",
                    propagate = TRUE, contrasts = 3)
save(GOresults, file = "GSE89408_Goresultsbp.rda")
bp <- GOresults[[3]]
bpsig <- bp[bp$padj < 0.05,]

# Añadir columna name con el nombre de los GOs
bpsig[,"name"] <- mdgsa::getGOnames(rownames(bpsig))

# Reordenar
bpsig <- bpsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
bpsig <- bpsig[order(bpsig$lor),]

save(bpsig, file = "GSE89408_GOresultsbp_sig_name.rda")
write.csv2(bpsig, file = 'GSE89408bpsig_name.csv')

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

save(KEGGresults, file = "GSE89408_KEGGresults.rda")
kegg <- KEGGresults[[3]]
keggsig <- kegg[kegg$padj<0.05,]

# Añadir columna name con el nombre de los KEGGs
keggsig[,"name"] <- mdgsa::getKEGGnames(rownames(keggsig))

# Reordenar
keggsig <- keggsig[c("name", "lor", "pval", "padj", "sd", "t", "N", "conv")]
keggsig <- keggsig[order(keggsig$lor),]

# Guardar
save(keggsig, file = "GSE89408_KEGGresults_sig_name.rda")
write.csv2(keggsig, file = 'GSE89408Keggsig_name.csv')
