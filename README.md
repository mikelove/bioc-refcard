# Bioconductor quick reference card

## Install

```
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("package1","package2"))
```

more information at http://www.bioconductor.org/install/

## help within R

```
?functionName
?"eSet-class" # classes need the '-class' on the end
vignette("topic")
browseVignettes(package="package") # show vignettes for the package
functionName # prints source code
getMethod("method","class")  # prints source code for methods
showMethods(classes="class") # show all methods for class
?"method,class-method" # for S4 objects, where the second 'method' here is the word method
# e.g.
?"plotMA,data.frame-method" # from library(geneplotter)
?"method.class" # for S3 objects
# e.g.
?"plot.lm"
```

## Annotations

```
# using one of the *.db annotation packges
library(pkg.db)
library(AnnotationDbi)
columns(pkg.db)
keytypes(pkg.db)
head(keys(pkg.db, keytype="PROBEID"))
# generates warning for 1:many mappings
res <- select(pkg.db, keys=k,
  columns=c("ENTREZID","ENSEMBL","SYMBOL"),
  keytype="PROBEID")
idx <- match(k, res$PROBEID)
res[idx,]
```

```
# get a transcript database, which stores exon, trancript, and gene information
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# or build a txdb from GTF file (e.g. downloadable from Ensembl FTP site)
txdb <- makeTranscriptDbFromGFF("file.GTF", format="gtf")

# or build a txdb from Biomart (however, not as easy to reproduce later)
txdb <- makeTranscriptDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# saving and loading
saveDb(txdb,file="txdb.RData")
loadDb("txdb.RData")

# extracting information from txdb into: 
#  - "Genomic Ranges" = a vector of ranges, each with sequence name, start, end, strand
#   e.g. exons or transcripts
#  - "Genomic Ranges List" = a list of Genomic Ranges
#  e.g. exons/transcripts split into a list by gene
tx <- transcripts(txdb)
exons <- exons(txdb)
exonsByGenes <- exonsBy(txdb, by="gene")
```

```
# map from one annotation to another
library(biomaRt)
m <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map <- getBM(mart = m,
  attributes = c("ensembl_gene_id", "entrezgene"),
  filters = "ensembl_gene_id", 
  values = some.ensembl.genes)
```

## GenomicRanges

```
library(GenomicRanges)
z <- GRanges("chr1",IRanges(1000001,1001000))
start(z)
end(z)
width(z)
flank(z, both=TRUE, width=100)
x %over% y  # logical vector of overlaps
fo <- findOverlaps(x,y) # returns a Hits object
queryHits(fo)   # which in x
subjectHits(fo) # which in y 
xsub <- keepSeqlevels(x, seqs)  # subsets x based on the seqlevels seqs
```

## SummarizedExperiment

```
library(GenomicRanges)
library(Rsamtools)
fls <- list.files(pattern="*.bam$")
bamlst <- BamFileList(fls)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx <- exonsBy(txdb, by="gene")
library(parallel)
options(mc.cores=4)  # summarizeOverlaps uses mclapply if parallel loaded
txhits <- summarizeOverlaps(tx, bamlst)  # lots of options in the man page
 # singleEnd, ignore.strand, fragments, etc.


# operations on SummarizedExperiments
assay(txhits)
colData(txhits)
rowData(txhits)
```

## Sequencing/sequence data

```
library(Rsamtools)
which <- GRanges("chr1",IRanges(1000001,1001000))
what <- c("rname","strand","pos","qwidth","seq")
param <- ScanBamParam(which=which, what=what)
reads <- scanBam(bamfile, param=param)
dnastringset <- scanFa(fastaFile, param=granges)
# DNAStringSet is defined in the Biostrings package
```

## RNA-Seq analysis

```
library(DESeq2)
dds <- DESeqDataSet(counts, DataFrame(condition), ~ condition)
	dds <- DESeq(dds)
	res <- results(dds)
```

```
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
```

```
library(limma)
design <- model.matrix(~ group)
dgel <- DGEList(exprs(e))
dgel <- calcNormFactors(dgel)
v <- voom(dgel,design,plot=FALSE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
topTable(fit,coef=2)
```

## Expression set

```
library(Biobase)
data(sample.ExpressionSet)
e <- sample.ExpressionSet
exprs(e)
pData(e)
fData(e)
```

## Get GEO dataset

```
library(GEOquery)
e <- getGEO("GSE9514")
```

## Microarray analysis

```
library(affy)
library(limma)
phenoData <- read.AnnotatedDataFrame("sample-description.csv")
eset <- justRMA("/celfile-directory", phenoData=phenoData)
design <- model.matrix(~ Disease, pData(eset))
fit <- lmFit(eset, design)
efit <- eBayes(fit)
topTable(efit, coef=2)
```

