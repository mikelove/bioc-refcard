# Bioconductor quick reference card

## install

   source("http://bioconductor.org/biocLite.R")
   biocLite()
   biocLite(c("package1","package2"))

more information at http://www.bioconductor.org/install/

## eSet

   library(Biobase)
   exprs(e)
   pData(e)
   fData(e)

## Annotation

   library(biomaRt)
   ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
   entrezmap <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), 
                      filters = "ensembl_gene_id", 
                      values = some.ensembl.genes, 
                      mart = ensembl)

## Genomic ranges

   library(GenomicRanges)
   z <- GRanges("chr1",IRanges(1000001,1001000))
   start(z)
   end(z)
   width(z)
   flank(z, both=TRUE, width=100)

## Sequencing data

   library(Rsamtools)
   which <- GRanges("chr1",IRanges(1000001,1001000))
   what <- c("rname","strand","pos","qwidth","seq")
   param <- ScanBamParam(which=which, what=what)
   reads <- scanBam(bamfile, param=param)

## RNA-seq

   library(DESeq)
   cds <- newCountDataSet(counts, condition)
   cds <- estimateSizeFactors(cds)
   cds <- estimateDispersions(cds)
   res <- nbinomTest(cds, "A", "B")
   # or with glm
   fit1 <- fitNbinomGLMs(cds, count ~ group + treatment)
   fit0 <- fitNbinomGLMs(cds, count ~ group)
   pvals <- nbinomGLMTest(fit1, fit0)

   library(edgeR)
   y <- DGEList(counts=counts,group=group)
   y <- estimateCommonDisp(y)
   y <- estimateTagwiseDisp(y)
   et <- exactTest(y)
   topTags(et)
   # or with glm
   design <- model.matrix(~group)
   y <- estimateGLMCommonDisp(y,design)
   y <- estimateGLMTrendedDisp(y,design)
   y <- estimateGLMTagwiseDisp(y,design)
   fit <- glmFit(y,design)
   lrt <- glmLRT(fit,coef=2)
   topTags(lrt)

## Microarray

   library(affy)
   library(limma)
   phenoData <- read.AnnotatedDataFrame("sample-description.csv")
   eset <- justRMA("/celfile-directory", phenoData=phenoData)
   design <- model.matrix(~ Disease, pData(eset))
   fit <- lmFit(eset, design)
   efit <- eBayes(fit)
   topTable(efit, coef=2)

