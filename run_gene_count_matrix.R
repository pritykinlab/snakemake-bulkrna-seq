library(rtracklayer)
library(SummarizedExperiment)
library(DESeq2)
library(Rsubread)
library(data.table)
library(tidyverse)
library(pheatmap)
library(yaml)

config <- yaml.load_file("config.yaml")

gtf <- rtracklayer::import(config$annotation)
gencode <- keepStandardChromosomes(gtf, pruning.mode = "coarse")
saveRDS(gencode, "processed_data/gencode.rds")
# input
bamdir <- paste(config$results, 'bam-RNA/', sep="/")
bamfiles <- list.files(bamdir, pattern = "*.bam$")
bamfiles.full <- paste0(bamdir, "/", bamfiles)

# genes and exons
print("Prepare collection of all exons for genes in GENCODE")
genes <- levels(as.factor(gencode$gene_name))
exons <- gencode[(gencode$gene_type == 'protein_coding') & (gencode$type == 'exon')]
exons.genes <- disjoin(split(exons, factor(exons$gene_name, levels = genes)))
print("prepared")

# count reads in exons
print("Count overlaps of genes (all exons) with reads")
print("Use files:")
print(bamfiles.full)
exons.table <- as.data.frame(exons) %>%
                   dplyr::select(GeneID = gene_name, Chr = seqnames,
                                 Start = start, End = end, Strand = strand)
gene.counts <- featureCounts(bamfiles.full, annot.ext = exons.table,
                             isPairedEnd = TRUE, requireBothEndsMapped = TRUE,
                             minOverlap = 80,
                             countChimericFragments = FALSE,
                             nthreads = 15)
print("count complete, save the result")
saveRDS(gene.counts, "processed_data/rna-counts.raw.rds")
