
Todo el código empleado se presenta a ccontinuación, este y el resto de archivos generados durante las práctica pueden ser encontrado en el siguiente repositorio: https://github.com/corgazp/ADO-PEC2 

````
## ----setup, include=FALSE------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, comment = NULL)
# Guardamos en el vector packages los paquetes que vamos a usar
packages_bioconductor <- c("Biobase","oligo", "arrayQualityMetrics", "genefilter","limma","rae230a.db",
"pd.rae230a","xtable", "annotate", "GOstats")
packages<-c("dplyr", "gplots")

# Comprobamos si algún paquete no está instalado para instalarlo
installed_packages_bioconductor <- packages_bioconductor %in% rownames(installed.packages())
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages_bioconductor == FALSE)) {
 BiocManager::install(packages[!installed_packages])
}
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Cargamos los paquetes, usamos invisible para no mostrar output
invisible(lapply(packages_bioconductor, library, character.only = TRUE))
invisible(lapply(packages, library, character.only = TRUE))


## ------------------------------------------------------------------------------------------------------------------------------
targets_df<- read.csv("targets.csv", header=TRUE, sep = ";")
sampleNames <- as.character(targets_df$shortName)
sampleColor <- as.character(targets_df$colors)
targets <- AnnotatedDataFrame(targets_df)
CELfiles <- targets_df$fileName
rawData <- read.celfiles(CELfiles, phenoData=targets)


## ------------------------------------------------------------------------------------------------------------------------------
#BOXPLOT
boxplot(rawData, which="all",las=1, main="Distribución de intensidad de datos RAW",cex.axis=0.6, 
col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE,
formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale)
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000)
       ,ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="rawData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)


## ------------------------------------------------------------------------------------------------------------------------------
rerun <- FALSE
if(rerun){
  arrayQualityMetrics(rawData,  reporttitle="QC_RawData", force=TRUE)
}
eset_rma<-rma(rawData) 
write.exprs(eset_rma, "normalized.txt")


## ------------------------------------------------------------------------------------------------------------------------------
annotation(eset_rma)<-"rae230a.db"
filtered <- nsFilter(eset_rma, var.func=IQR, var.cutoff=0.75, var.filter=TRUE, 
filterByQuantile=TRUE, require.entrez = TRUE)
filtered_eset<-filtered$eset
filtered_data <- exprs(filtered_eset)
colnames(filtered_data) <- pData(filtered$eset)$ShortName


## ------------------------------------------------------------------------------------------------------------------------------
treat<-pData(filtered_eset)$group
lev<- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design) <- levels(lev)
rownames(design) <- sampleNames
print(design)


## ------------------------------------------------------------------------------------------------------------------------------
cont.matrix1 <- makeContrasts( 
        control.vs.ethanol = control-ethanol,
        levels = design)
comparisonName <- "Efecto del alcohol"
print(cont.matrix1)


## ------------------------------------------------------------------------------------------------------------------------------
fit1 <- lmFit(filtered_data, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1 <- eBayes(fit.main1)
topTab <-  topTable (fit.main1, number=nrow(fit.main1), 
coef="control.vs.ethanol", adjust="fdr",lfc=2, p.value=0.05)


## ------------------------------------------------------------------------------------------------------------------------------
topTab


## ------------------------------------------------------------------------------------------------------------------------------
keytypes(rae230a.db)
anotaciones<- AnnotationDbi::select (rae230a.db, 
keys=rownames(filtered_data), columns=c("ENTREZID", "SYMBOL"))

## ------------------------------------------------------------------------------------------------------------------------------
topTabAnotada <- topTab %>%  
  mutate(PROBEID=rownames(topTab)) %>%
  left_join(anotaciones) %>% 
  arrange(P.Value) %>%
  select(7,8,9, 1:6)

head(topTabAnotada)

## ------------------------------------------------------------------------------------------------------------------------------
write.csv2(topTabAnotada, file="Genes seleccionados.csv")
print(xtable(topTab,align="lllllll"),type="html",
html.table.attributes="",file="Genes seleccionados.html")


## ------------------------------------------------------------------------------------------------------------------------------
genenames <- AnnotationDbi::select(rae230a.db, 
                    rownames(fit.main1), c("SYMBOL"))
volcanoplot(fit.main1, highlight=10, names=genenames, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))
abline(v = c(-3, 3))


pdf("Volcanos.pdf")
volcanoplot(fit.main1, highlight = 10, names = genenames, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()

## ------------------------------------------------------------------------------------------------------------------------------
plot(1:30)
selectedRows <- rownames(filtered_data) %in% rownames(topTab)
selectedData <- filtered_data[selectedRows,]

#HEATMAP PLOT
my_palette <- colorRampPalette(c("yellow", "blue"))(n = 299)
library()
heatmap.2(selectedData,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap Control.vs.Ethanol FC>=2",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=c(rep("yellow",3),rep("blue",3)),
          tracecol=NULL,
          srtCol=30)
pdf("Heatmap.pdf")
dev.off()


## ------------------------------------------------------------------------------------------------------------------------------
probesUniverse <- rownames(filtered_data)
entrezUniverse<- AnnotationDbi::select(rae230a.db, probesUniverse, "ENTREZID")$ENTREZID

topProbes <-   rownames(selectedData)
entrezTop<- AnnotationDbi::select(rae230a.db, topProbes, "ENTREZID")$ENTREZID

# Eliminamos posibles duplicados

topGenes <- entrezTop[!duplicated(entrezTop)]
entrezUniverse <- entrezUniverse[!duplicated(entrezUniverse)]


## ------------------------------------------------------------------------------------------------------------------------------
GOparams = new("GOHyperGParams",
    geneIds=topGenes, universeGeneIds=entrezUniverse,
    annotation="rae230a.db", ontology="BP",
    pvalueCutoff=0.01)


## ----runORAnalysis-------------------------------------------------------------------------------------------------------------
GOhyper = hyperGTest(GOparams)


## ----summarizeORAesults--------------------------------------------------------------------------------------------------------
head(summary(GOhyper))
dim(summary(GOhyper))


## ----ORAreport-----------------------------------------------------------------------------------------------------------------
# Creamos un informe html con los resultados
GOfilename ="GOResults.html"
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))


## ------------------------------------------------------------------------------------------------------------------------------

````
