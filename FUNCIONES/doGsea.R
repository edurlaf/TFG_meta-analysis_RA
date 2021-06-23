#' Gene Set Enrichment Analysis
#'
#'\code{doGSEA} performs Gene Set Enrichment Analysis from a limma fit object
#'
#'
#'@param fit MArrayLM limma object
#'@param annot data-frame with ENTREZIDs in the first column and GO IDs in the second
#'@param ontology string indicating whether you want to use biological
#'process ("bp"), cellular component ("cc") or molecular function ("mf"). If
#'GSEA is going to be performed on KEGG pathways, this must be ("path").
#'@param propagate boolean, uses mdgsa::propagateGO() or not
#'@param minBlockSize minimum block size kept by annotFilter
#'@param maxBlockSize maximum block size kept by annotFilter
#'@param contrasts number of contrasts
#'@return results from analysis filtered over 0.05 adjusted p-value
#'@import limma
#'@import mdgsa
#'@import biomaRt
#'@export

doGSEA <- function(fit, annot, ontology, propagate = FALSE, minBlockSize = 10,
                   maxBlockSize = 500, contrasts = 1){

  # Necesitamos que los ENTREZID sean caracteres, si son enteros la N de uvGsa
  # será errónea
  annot[,1] <- as.character(annot[,1])
  annot <- mdgsa::annotMat2list(annot)
  if (ontology != "path"){
    # Clasificamos los GO en las tres ontologías y elegimos la que queremos
    annot <- mdgsa::splitOntologies(annot, na.rm = TRUE, verbose = TRUE)
    ifelse(ontology == "bp", annotBp <- annot$bp,
           ifelse(ontology == "cc", annotBp <- annot$cc,
                  ifelse(ontology == "mf", annotBp <- annot$mf, annotBp <- annot)))
  } else {
    # En el caso de los pathways, no habrá que separar nada
    annotBp <- annot
  }

  if (propagate == TRUE) {
    # Propagación de la anotación. En el caso de los pathways, no propagará nada
    annotBp <- mdgsa::propagateGO(annotBp, verbose = TRUE)
  }
  # Filtrado de las anotaciones muy específicas (10) o muy genéricas (500)
  annot <- mdgsa::annotFilter(annot = annotBp, minBlockSize = minBlockSize, maxBlockSize = maxBlockSize)


  if (contrasts == 3) {
    rindex1 <- mdgsa::pval2index(pval = fit$p.value[,1], sign = fit$t[,1])
    rindex1 <- mdgsa::indexTransform(index = rindex1, method = "normalize")

    rindex2 <- mdgsa::pval2index(pval = fit$p.value[,2], sign = fit$t[,2])
    rindex2 <- mdgsa::indexTransform(index = rindex2, method = "normalize")

    rindex3 <- mdgsa::pval2index(pval = fit$p.value[,3], sign = fit$t[,3])
    rindex3 <- mdgsa::indexTransform(index = rindex3, method = "normalize")

    goResults <- mdgsa::uvGsa(index = rindex1, annot = annot, fulltable = TRUE)
    goRes1 <- goResults#[goResults$padj < 0.05,]
    goResults <- mdgsa::uvGsa(index = rindex2, annot = annot, fulltable = TRUE)
    goRes2 <- goResults#[goResults$padj < 0.05,]
    goResults <- mdgsa::uvGsa(index = rindex3, annot = annot, fulltable = TRUE)
    goRes3 <- goResults#[goResults$padj < 0.05,]

    return(list("GORes1" = goRes1, "GORes2" = goRes2, "GORes3" = goRes3))
  } else {
    rindex0 <- mdgsa::pval2index(pval = fit$p.value, sign = fit$t)
    rindex <- mdgsa::indexTransform(index = rindex0, method = "normalize")
    goResults <- mdgsa::uvGsa(index = rindex, annot = annot, fulltable = TRUE, p.adjust.method = "BH")
    goRes <- goResults
    return(goRes)
  }

}
