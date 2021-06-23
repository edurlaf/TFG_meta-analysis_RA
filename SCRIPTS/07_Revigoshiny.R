if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install("rrvgo")

library(rrvgo)
# go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
go_analysis <- bpsig[,"lor"]
names(go_analysis) <- rownames(bpsig)
simMatrix <- calculateSimMatrix(names(go_analysis),
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                go_analysis,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
scatterPlot(simMatrix, reducedTerms)
go_abs <- abs(go_analysis)
reducedTerms <- reduceSimMatrix(simMatrix,
                                go_abs,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
scatterPlot(simMatrix, reducedTerms)


rrvgo::shiny_rrvgo()
install.packages("heatmaply")
