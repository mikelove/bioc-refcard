# Bioconductor quick reference card

## Install

    source("http://bioconductor.org/biocLite.R")
    biocLite()
    biocLite(c("package1","package2"))

more information at http://www.bioconductor.org/install/

## help within R
   
    ?functionName
    ?"eSet-class"         # classes need the '-class' on the end
    vignette("topic")
    browseVignettes(package="package") # show vignettes for the package
    functionName                 # prints source code
    getMethod("method","class")  # prints source code for methods
    showMethods(classes="class") # show all methods for class

## Annotations

    # get a transcript database
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    
    # or alternatively build a transcript database from biomart
    txdb <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    saveDb(txdb,file="txdb.RData")
    loadDb("txdb.RData")
   
    # get GRanges or GRangesList of genomic features
    tx <- transcripts(txdb)
    exons <- exons(txdb)
    exonsByGenes <- exonsBy(txdb, by="gene")

    # map from one annotation to another
    library(biomaRt)
    ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    entrezmap <- getBM(mart = ensembl,
    	      	       attributes = c("ensembl_gene_id", "entrezgene"), 
    	               filters = "ensembl_gene_id", 
                       values = some.ensembl.genes)

## GenomicRanges

    library(GenomicRanges)
    z <- GRanges("chr1",IRanges(1000001,1001000))
    start(z)
    end(z)
    width(z)
    flank(z, both=TRUE, width=100)
    x %over% y                      # logical vector of overlaps
    fo <- findOverlaps(x,y)         # returns a Hits object
    queryHits(fo)                   # which in x
    subjectHits(fo)                 # which in y 
    xsub <- keepSeqlevels(x, seqs)  # subsets x based on the seqlevels seqs

## SummarizedExperiment

    library(GenomicRanges)
    library(Rsamtools)
    fls <- list.files(pattern="*.bam$")
    bamlst <- BamFileList(fls)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    tx <- exonsBy(txdb, by="gene")
    library(parallel)
    options(mc.cores=4)                      # summarizeOverlaps uses mclapply if parallel loaded
    txhits <- summarizeOverlaps(tx, bamlst)  # lots of options in the man page
                                             # singleEnd, ignore.strand, fragments, etc.

    # operations on SummarizedExperiments
    assay(txhits)
    colData(txhits)
    rowData(txhits)

## Sequencing/sequence data

    library(Rsamtools)
    which <- GRanges("chr1",IRanges(1000001,1001000))
    what <- c("rname","strand","pos","qwidth","seq")
    param <- ScanBamParam(which=which, what=what)
    reads <- scanBam(bamfile, param=param)
    dnastringset <- scanFa(fastaFile, param=granges)
    # DNAStringSet is defined in the Biostrings package

## RNA-Seq analysis

    library(DESeq)
    cds <- newCountDataSet(counts, condition)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    res <- nbinomTest(cds, "A", "B")
    # glm
    fit1 <- fitNbinomGLMs(cds, count ~ group + treatment)
    fit0 <- fitNbinomGLMs(cds, count ~ group)
    pvals <- nbinomGLMTest(fit1, fit0)

    library(edgeR)
    y <- DGEList(counts=counts,group=group)
    y <- calcNormFactors(y)
    # exact test
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    topTags(et)
    # glm
    design <- model.matrix(~group)
    y <- estimateGLMCommonDisp(y,design)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    topTags(lrt)

## Expression set

    library(Biobase)
    data(sample.ExpressionSet)
    e <- sample.ExpressionSet
    exprs(e)
    pData(e)
    fData(e)

## Get GEO dataset

    library(GEOquery)
    e <- getGEO("GSE9514")

## Microarray analysis

    library(affy)
    library(limma)
    phenoData <- read.AnnotatedDataFrame("sample-description.csv")
    eset <- justRMA("/celfile-directory", phenoData=phenoData)
    design <- model.matrix(~ Disease, pData(eset))
    fit <- lmFit(eset, design)
    efit <- eBayes(fit)
    topTable(efit, coef=2)

## Gene ontology set enrichment

    library(org.Hs.eg.db)
    library(GOstats)
    params <- new("GOHyperGParams", geneIds = entrezlist, 
                  universeGeneIds = NULL, annotation ="org.Hs.eg.db", 
                  ontology="BP", pvalueCutoff=.001, conditional=FALSE, 
                  testDirection="over")
    hgOver <- hyperGTest(params)
    hgSummary <- summary(hgOver)

