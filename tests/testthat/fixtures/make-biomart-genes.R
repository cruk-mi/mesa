# Regenerate biomart-genes.csv, the cached Ensembl lookups that
# localMockBiomart() (helper-biomart.R) serves to plotGeneHeatmap() so the
# tests run offline. Needs internet. Run from the package root:
#   Rscript tests/testthat/fixtures/make-biomart-genes.R
# Each query matches one made by plotGeneHeatmap() in the tests: the archive
# host decides the genome build (jul2022 = GRCh38, feb2014 = GRCh37,
# jul2023 = GRCm39). exampleMouse stores a nov2020 (GRCm38) Mart, but that
# archive is retired, so its row takes the jul2023 gene and mm10
# coordinates from TxDb.Mmusculus.UCSC.mm10.knownGene.

queries <- data.frame(
    host = c("jul2022", "feb2014", "feb2014", "jul2023"),
    dataset = c("hsapiens_gene_ensembl", "hsapiens_gene_ensembl",
        "hsapiens_gene_ensembl", "mmusculus_gene_ensembl"),
    filter = c("hgnc_symbol", "hgnc_symbol", "hgnc_symbol", "mgi_symbol"),
    gene = c("HOXA10", "JAM2", "EWSR1", "Fbxl18")
)

attributes <- c("hgnc_symbol", "description", "chromosome_name",
    "start_position", "end_position", "strand", "ensembl_gene_id")

rows <- lapply(seq_len(nrow(queries)), function(i) {
    q <- queries[i, ]
    mart <- biomaRt::useMart("ensembl", dataset = q$dataset,
        host = paste0("https://", q$host, ".archive.ensembl.org"))
    res <- biomaRt::getBM(attributes = attributes, filters = q$filter,
        values = q$gene, mart = mart)
    stopifnot(nrow(res) >= 1)
    cbind(host = q$host, filter = q$filter, value = q$gene, res)
})

genes <- do.call(rbind, rows)

# nov2020 (GRCm38) row: jul2023 gene with mm10 coordinates.
mouse <- genes[genes$host == "jul2023", ]
entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = "Fbxl18",
    keytype = "SYMBOL", column = "ENTREZID")
mm10 <- GenomicFeatures::genes(
    TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
    filter = list(gene_id = entrez))
mouse$host <- "nov2020"
mouse$chromosome_name <- sub("^chr", "", as.character(
    GenomicRanges::seqnames(mm10)))
mouse$start_position <- GenomicRanges::start(mm10)
mouse$end_position <- GenomicRanges::end(mm10)
mouse$strand <- ifelse(as.character(GenomicRanges::strand(mm10)) == "+",
    1L, -1L)
genes <- rbind(genes, mouse)

utils::write.csv(genes,
    "tests/testthat/fixtures/biomart-genes.csv", row.names = FALSE)
