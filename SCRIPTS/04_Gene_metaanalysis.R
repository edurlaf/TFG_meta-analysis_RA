######################################
### metaanalysis for 5 GEO studies: GSE93272, GSE89408, GSE110169, GSE117769, GSE15573
### 05/2021
######################################


# STEP 1. Preparing input for meta-analysis: LOR and SE matrix
# ===============================================================

## starting
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
R.version.string 

## clean the working space
rm (list = ls ())

## load libraries
library(Biobase)
library(metafor)

## load input

### previously for each study, we need results from limma in a file .RData
### including this structure (dataframe): 

# EntrezID logFC      CI.L      CI.R  AveExpr        t      P.Value adj.P.Val
# 728489 0.4044977 0.2172210 0.5917743 9.211481 4.473528 0.0001809354 0.8971414
# 283989 0.3552017 0.1867096 0.5236937 7.963093 4.366299 0.0002357869 0.8971414
# 149095 0.2438707 0.1236510 0.3640904 7.126072 4.201477 0.0003541950 0.8971414
 
# EntrezID      B         SE      coef   se.coef   se2.coef
# 728489 -0.2008863 0.09554933 0.4044977 0.2120544 0.09042029
# 283989 -0.3626524 0.08596534 0.3552017 0.1907844 0.08135075
# 149095 -0.6133722 0.06133656 0.2438707 0.1361253 0.05804404


ficheros <- dir(pattern= "inputMA.RData")
ficheros 

# we search a list including all unique ID genes for all studies
genes<-NULL
for (fi in ficheros){
  load (fi)
  genes <- c(genes, rownames(top.table))
}

length (genes) #84208
genes <- unique (genes)
length (genes) #23062
genes <- sort (genes)


### generating matrix with all logFC for all studies
mat.logFC <- matrix (NA, nrow = length (genes), ncol = length (ficheros))
rownames (mat.logFC) <- genes
colnames (mat.logFC) <- strsplit (ficheros, "inputMA.RData")
head (mat.logFC)

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, "inputMA.RData"))
  logFC <- top.table$logFC
  names (logFC) <- top.table$EntrezID
  mat.logFC[, co] <- logFC[rownames(mat.logFC)] 
}

head (mat.logFC)
tail(mat.logFC)
table (is.na(mat.logFC))
dim (mat.logFC)

# select genes included at least in 2 or more studies
mat.logFC.NA <- is.na(mat.logFC) #donde hay NA pondrá TRUE
head(mat.logFC.NA)
sum.NA <-  apply(mat.logFC.NA, 1, sum) # suma por filas. Porque True = 1 y F = 0. Va a sacar en cuantos estudios NO está
table(sum.NA)
min.sum.NA <- sum.NA < ncol(mat.logFC) - 1
table(min.sum.NA) #F:3087, T:19975

# filter by min.sum.NA. Se queda con los que esté en al menos 2 estudios
mat.logFC <- mat.logFC[min.sum.NA == T, ]
dim(mat.logFC)

### generating matrix with all SE for all studies  
mat.SE <- matrix (NA, nrow = length (genes), ncol = length (ficheros))
rownames (mat.SE) <- genes
colnames (mat.SE) <- strsplit (ficheros, "inputMA.RData")
head (mat.SE)


# SE de la opción Gordon
for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, "inputMA.RData"))
  SE <- top.table$SE.coef
  names (SE) <- top.table$EntrezID
  mat.SE[, co] <- SE[rownames(mat.SE)] # se asegura de que estén en el mismo orden que la matriz en la que los va a meter
}

head (mat.SE)
tail(mat.SE)
table (is.na(mat.SE))
dim (mat.SE)

# filter by min.sum.NA
mat.SE <- mat.SE[min.sum.NA == T, ] #se queda solo con los que estan en al menos 2 estudios
dim(mat.SE)


# STEP 2. Meta-analysis for genes
# ===============================================================

# suppose between-study variance is non-zero.
# there are different methods to estimate this variance:
# DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
# Now we have logFC and SE  (not VARIANCE), so:
# yi -> logFC   sei -> SE
# result.lor <- rma(yi = mat.logFC[1, ], 
#                   sei = mat.SE[1, ],   #pay attention, not vi (varianze)   
#                   method = "DL") # DerSimonian-Laird.


# explore the function to do the meta-analysis
?rma # Function to fit the meta-analytic fixed- and random/mixed-effects models with or without moderators via linear (mixed-effects) models. 

MA <- lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SE[x, ],
                             method = "DL")})
# me han salido más de 50 warnings. Son los NAs

# MA <- lapply(1:length(rownames(mat.logFC)),
#              function(x){rma(yi = mat.logFC[x, ],
#                              sei = mat.SE[x, ],
#                              method = "FE")})

names (MA) <- rownames(mat.logFC)
class (MA)
length(MA)
head (MA)
MA[[1]]

# Random-Effects Model (k = 3; tau^2 estimator: DL)
# 
# tau^2 (estimated amount of total heterogeneity): 0.0086 (SE = 0.0232)
# tau (square root of estimated tau^2 value):      0.0927
# I^2 (total heterogeneity / total variability):   37.29%
# H^2 (total variability / sampling variability):  1.59
# 
# Test for Heterogeneity:
#   Q(df = 2) = 3.1892, p-val = 0.2030
# 
# Model Results:
#   
#   estimate      se    zval    pval    ci.lb   ci.ub 
# 0.0555  0.0863  0.6431  0.5202  -0.1137  0.2247    
# 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#result.logFC$pval      #p-value about logFC = 0
#result.logFC$ci.lb     #IC down
#result.logFC$ci.ub     #IC up
#result.logFC$b         #estimation of combined logFC

#data.frame including all detailed results:
result_meta <- as.data.frame(do.call("rbind",
                                     lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub, 
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                              })))

colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                           "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")

p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr") #coge los pvalue
p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
result_meta <- round(cbind(result_meta, p.adjust.fdr, p.adjust.BY), 3) #añade los nuevos pvalue como nuevas columna, redondeando a 3 decimales
head(result_meta)
#preparar para input GSEA
result_meta <- result_meta[order(result_meta$logFC, decreasing = T),]
save(result_meta, file = "../MA/result_meta_gen.rda")
class(result_meta)
MAgenes <- data.matrix(result_meta, rownames.force = NA)
save(MAgenes, file = "MA_genes_results.rda")


# significant genes
corte = 0.05
table(result_meta[, "pvalue"] < corte) #F: 18829 T: 1155
table(result_meta[, "p.adjust.fdr"] < corte) # F: 19973 T: 5
table(result_meta[, "p.adjust.BY"] < corte) # F: 19975


# add number of studies where the gene is evaluated
n.studies <-  ncol(mat.logFC) - sum.NA 
table(n.studies) #    1    2    3    4    5 
                  # 3087 2079 4412 3693 9791 
n.studies <- n.studies [rownames(mat.logFC)] #reordena para que tenga mismo orden que mat.logFC
length(n.studies)
result_meta[, "n.studies"]  <- n.studies #añade la columna con el numero de estudios en los que se ha evaluado
head(result_meta)
summary(result_meta$p.adjust.fdr)

# corte = 0.05
sig.genes.df = result_meta[result_meta$p.adjust.fdr < corte,] 
dim(sig.genes.df)
symbol <- getSYMBOL(rownames(sig.genes.df), data='org.Hs.eg')
sig.genes.df$Symbol <- symbol




# corte = 0.1
sig.genes01.df = result_meta[result_meta$p.adjust.fdr < 0.1,] 
dim(sig.genes01.df)

sig.genes01.df <- sig.genes01.df[order(sig.genes01.df$logFC, decreasing = T),]
library(annotate)
library(org.Hs.eg.db)
symbol <- getSYMBOL(rownames(sig.genes01.df), data='org.Hs.eg')
sig.genes01.df$Symbol <- symbol


# STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS
# ===============================================================

#add 4 new variables about influence & sensitivity analysis:  

for (i in rownames(sig.genes.df)){
  #define studies for each function (not NA)
  estudios <- colnames(mat.logFC)[!mat.logFC.NA[i,]] # saca los nombres de estudios con valor, que no tengan NA
  
  #influence info 1: 
  #number of studies where the sign of the logOR is the same  of the global logOR:
  sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios))) #compara logFC con b ¿?. Mira el signo a ver si en los estudios individuales coincide con el global del MA. Suma cuantos coinciden si es positivo habrá sobreexpresión
  
  #influence info 2: how many studies could be influencers?
  inf <- influence(MA[[i]])
  res <- paste(estudios[inf$is.infl], collapse = ",")  
  sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)
  
  #sensivity analysis # como de robusta es mi conclusión (si dejo un estudio fuera me sale lo mismo?)
  l1 <-as.data.frame(leave1out(MA[[i]])) # me hace todas las combinaciones dejando un estudio fuera
  rownames(l1) <- estudios
  
  #1. p.value about differences between all estimates from leave one out
  #   and global estimate) # si hay dif significativa en alguna de esas combinaciones. pvalue peque, resultado muy diferente respecto del resto> conclusion no robusta
  sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                         mu=as.numeric(MA[[i]]$b))$p.value
  #2. number of  studies where pvalue > 0.05 
  # (we hope p-values < 0.05, significant estimates) 
  res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")  
  sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}


## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:

#1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
table(sig.genes.df$infl.same.sign.logFC)

#2. INFLUENCE STUDIES. How many functions including influence studies?
table(sig.genes.df$infl.nstudies=="non")

#3. SENSITIVITY. In global, are there many functions with differences in the estimate?
table(sig.genes.df$sensi.global < 0.05)

#4. SENSITIVITY.  How many functions including changes in the significance about 
# its new estimate  after leave1out? 
table(sig.genes.df$sensi.specific == "all.p.values < 0.05")


#save final results:
cat ("ID\t", file = "sig.genes.df.txt")
write.table(sig.genes.df, file = "sig.genes.df.txt", sep ="\t", quote = F, 
            append = TRUE, row.names = T)



# STEP 4. Visualization of significant genes
# ===============================================================

#select significant functions to visualize:
sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.05,]
symbol <- getSYMBOL(rownames(sig.results), data='org.Hs.eg')
sig.results$Symbol <- symbol

sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.2,]
sig.results
dim(sig.results)

setwd("../MA/prueba")
selMethod <- "DL"

for (i in 1:nrow(sig.results)){
  mygenes <- rownames(sig.results)[i]
  symbols <- sig.results[mygenes, "Symbol"]
  res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")
  
  #FOREST PLOT
  png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 , 
       height = 960, res = 200) 
  forest(res, 
         xlab="logFC", cex=0.7,
         mlab=bquote(paste("RE Model for All Studies (", I^2, " = ",.(formatC(res$I2, digits=1, format="f")), "%)")),
         col = "blue", 
         main = paste("\n", mygenes, "(", symbols, ")", sep=" "))    
  text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
# paste(selMethod, " Model for All Studies (", I^2, " = ",.(formatC(res$I2, digits=1, format="f")),"%)", sep = " ")#   text(-3, 0, pos=4, cex=0.75, bquote(
# paste("RE Model for All Studies (Q = ",
# .(formatC(res$QE, digits=2, format="f")), 
# ", df = ", .(res$k - res$p),", p = ", 
# .(formatC(res$QEp, digits=2, format="f")),
#  "; ", I^2, " = ",
# .(formatC(res$I2, digits=1, format="f")), 
#  "%)", "; ", tau^2 == 
# .(formatC(res$tau2, digits=2, format="f")))))        
  
  dev.off()
  
  #FUNNEL PLOT
  png (filename = paste("FUNNEL_", mygenes,".png", sep =""), width = 960 , 
       height = 960, res = 200) 
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance", 
         back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
  par(mfrow=c(1,1))
  dev.off()
  
  #INFLUENCE PLOTS 
  # That shows various diagnostic measures
  png (filename = paste("INFLUENCE_", mygenes,".png", sep =""), width = 960 , 
       height = 960, res = 200) ##CAMBIAR
  inf <- influence(res)
  #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
  plot(inf)
  dev.off()
  
}


 
### EXIT
warnings ()
sessionInfo ()
q ("no")
   

