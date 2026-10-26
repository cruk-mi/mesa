# Offline stand-ins for biomaRt::useMart() and biomaRt::getBM(), so the
# plotGeneHeatmap() tests need no network (#123). Lookups are served from
# fixtures/biomart-genes.csv, which make-biomart-genes.R in that folder
# regenerates.
#
# plotGeneHeatmap() retries getBM() and then reports any failure as "Could not
# retrieve gene information from biomart ... check your internet connection".
# In a mocked test that message means the gene or archive host is missing from
# the CSV, not that the network is down.

.cachedBiomartGenes <- function() {
    utils::read.csv(testthat::test_path("fixtures", "biomart-genes.csv"),
        stringsAsFactors = FALSE)
}

localMockBiomart <- function(env = parent.frame()) {
    genes <- .cachedBiomartGenes()
    # plotGeneHeatmap() probes the default archive before connecting (#183)
    testthat::local_mocked_bindings(
        ensemblReachable = function(host, ...) TRUE,
        .env = env
    )
    testthat::local_mocked_bindings(
        useMart = function(biomart, dataset, host, ...) {
            methods::new("Mart", biomart = biomart, dataset = dataset,
                host = host)
        },
        getBM = function(attributes, filters, values, mart, ...) {
            atHost <- vapply(genes$host, grepl, logical(1), x = mart@host)
            keep <- atHost & genes$filter == filters & genes$value %in% values
            hits <- genes[keep, ]
            if (nrow(hits) == 0) {
                stop("No cached gene for ", values, " at ", mart@host)
            }
            hits[, attributes]
        },
        .package = "biomaRt",
        .env = env
    )
}
